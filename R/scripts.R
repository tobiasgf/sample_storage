library(pairwiseAdonis)
library(rstatix)
library(here)
library(tidyverse)

# some functions used in the data processing

itsx_on_dada2 <- function(dada2_table, asv_file_name = "asv_seqs_from_dada2.fasta", itsout_name = "itsxout", org="F"){
  require("here")
  require("dada2")
  require("seqinr")
  asv_file_namex <- here::here(asv_file_name)
  its2_file_namex <- here::here(paste0(itsout_name,".ITS2.fasta"))
  uniquesToFasta(getUniques(dada2_table), fout=asv_file_namex, ids=paste0("Seq", seq(length(getUniques(dada2_table)))))
  
  #perform ITSx
  commandX <- paste0('ITSx -i "',asv_file_namex,'" --complement F -t "', org, '" --cpu 8 --preserve T -o ',itsout_name)
  system(commandX)
  
  #read in ITSx data
  its_seq <- read.fasta(its2_file_namex, as.string = T, seqonly = F)
  its_seq <- unlist(its_seq)
  
  #restrict result to accepted fungal ITS2 sequences and transpose
  positions <- as.numeric(gsub("Seq","",names(its_seq)))
  
  #aggregate sequences that are now identical
  pos_agg <- aggregate(x = t(dada2_table)[positions,], by = list(its_seq), FUN = "sum")
  rownames(pos_agg) <- pos_agg[,1]
  pos_agg <- t(pos_agg[,-1])
  return(pos_agg)
}


sumSequenceTables <- function(table1, table2, ..., orderBy = "abundance") {
 # Combine passed tables into a list
 tables <- list(table1, table2)
 tables <- c(tables, list(...))
 # Validate tables
 if(!(all(sapply(tables, dada2:::is.sequence.table)))) {
  stop("At least two valid sequence tables, and no invalid objects, are expected.")
 }
 sample.names <- rownames(tables[[1]])
 for(i in seq(2, length(tables))) {
  sample.names <- c(sample.names, rownames(tables[[i]]))
 }
 seqs <- unique(c(sapply(tables, colnames), recursive=TRUE))
 sams <- unique(sample.names)
 # Make merged table
 rval <- matrix(0L, nrow=length(sams), ncol=length(seqs))
 rownames(rval) <- sams
 colnames(rval) <- seqs
 for(tab in tables) {
  rval[rownames(tab), colnames(tab)] <- rval[rownames(tab), colnames(tab)] + tab
 }
 # Order columns
 if(!is.null(orderBy)) {
  if(orderBy == "abundance") {
   rval <- rval[,order(colSums(rval), decreasing=TRUE),drop=FALSE]
  } else if(orderBy == "nsamples") {
   rval <- rval[,order(colSums(rval>0), decreasing=TRUE),drop=FALSE]
  }
 }
 rval
}

#### function for extracting mean and calculating sd ####
data_summary <- function(data, varname, groupnames){
 require(plyr)
 summary_func <- function(x, col){
  c(mean = mean(x[[col]], na.rm=TRUE),
    sd = sd(x[[col]], na.rm=TRUE),
    sem = sd(x[[col]], na.rm=TRUE)/sqrt(length(x[[col]])))
 }
 data_sum<-ddply(data, groupnames, .fun=summary_func,
                 varname)
 data_sum <- rename(data_sum, c("mean" = varname))
 return(data_sum)
}

#### script for extracting richness measures etc ####
prepare_data <- function(community_table = table, sam_info = info, use_partial = FALSE){
  require(vegan)
  require(tidyr)
  
  rawtab <- community_table
 
 #synchronize content and sequence of table/list
  sampleinfo1 <- sam_info[sam_info$sample %in% rownames(rawtab),]
  rawtab <- rawtab[sampleinfo1$sample,]

  # richness of abundant taxa
  rawtab_ab <- rawtab
  rawtab_ab <- rawtab_ab[,colSums(rawtab_ab)>0]
  rawtab_ab <- sweep(rawtab_ab,2,colSums(rawtab_ab),"/")
  rawtab_ab[rawtab_ab<0.01] <- 0
  ab_richness <- rowSums(rawtab_ab>0)
  
  #richness & seq depth
  raw_richness <- rowSums(rawtab>0)
  #seq_depth <- rowSums(rawtab)
 
  #rarefy and normalize
  rr_tab <- rrarefy(rawtab, quantile(rowSums(rawtab))[2])
  #rr_tab[rr_tab<0] <- 0
  rr_tab_norm <- sweep(rr_tab, 1, rowSums(rr_tab), '/')
  
  #specaccum time zero samples
  rr_tab_norm_spa <- rr_tab[sampleinfo1$harvest == 0,]
  rr_tab_norm_spa <- rr_tab_norm_spa[,colSums(rr_tab_norm_spa)>0]
  #rr_tab_norm_spa[rr_tab_norm_spa>0] <- 1
  s1 <- vector()
  s2 <- vector()
  trx <- vector()
  for(i in 1:100){
    samp <- sample(seq(dim(rr_tab_norm_spa)[1]),3)
    A <- colSums(rr_tab_norm_spa[samp,])
    #A[A>1] <- 1
    B <- colSums(rr_tab_norm_spa[-samp,])
    #B[B>1] <- 1
    #s2[i] <- sum(A == 1)
    #s1[i] <- sum(B==1 & A==1)
    s2[i] <- sum(A > 0)
    s1[i] <- sum(B>0 & A>0)
    trx[i] <- sum(A[A>0 & B == 0])/sum(A)
  }
  
  new_data <- list(mean_new_otus = mean(s2-s1), sd_new_otus = sd(s2-s1), percent_new_otus = round(mean(round((s2-s1)*100/s2,2)),1), mean_percent_new_reads = round(mean(trx*100),1), sdn_percent_new_reads = round(sd(trx*100),1))
  
  
 #rarefied richness and shannon
  rr_richness <- rowSums(rr_tab>0)
  rr_norm_shannon <- diversity(rr_tab_norm)
  #rr_norm_simp <- diversity(rr_tab_norm, "simp")
  rr_norm_even <- rr_norm_shannon/log(specnumber(rr_tab_norm))
 
 #distance metrics
 dis_mat <- as.matrix(vegdist(decostand(rr_tab, method = "hellinger"), method = "bray"))
 mds <- metaMDS(dis_mat, distance = "bray", k=2, try = 500, trymax = 4000)
 
 #####what is this used for??
 # sampleinfo2 <- sampleinfo1 #%>% filter(batch %in% c("original","lts"))
 # dis_mat2 <- dis_mat[sampleinfo2$sample,sampleinfo2$sample]
 # dis_mat3 <- as.data.frame(dis_mat2[sampleinfo2$sample[sampleinfo2$harvest == 0],]) #sampleinfo2$harvest != 0
 # distances <- gather(dis_mat3, sample, dissimilarity)
 # dis_data <- left_join(distances, sampleinfo1, by = "sample")
 
 strata_seq <- rep("discard",107)
 strata_seq[sampleinfo1$harvest ==0] <- "time0"
 cen_dist <- dist_to_centroids(dis_mat, strata_seq)
 cen_dist2 <- cen_dist$CentroidDistance[cen_dist$CentroidGroup == "time0"]
 
 mdsp <- mds$points[,1:2]
 names(mdsp) <- c("NMDS1","NMDS2")
 
 nmds_tab <- as.data.frame(mds$points)
 
 sampleinfox <- sampleinfo1
 sampleinfox$group <- "nul"
 sampleinfox$group <- paste0(sampleinfox$exposure, "_",sampleinfox$temperature, "_", sampleinfox$harvest)
 sampleinfox$group[sampleinfox$harvest == 0] <- "reference"
 sampleinfox$group <- factor(sampleinfox$group , ordered = FALSE )
 sampleinfox <- within(sampleinfox, group <- relevel(group, ref = "reference"))
 pwa <- pairwise.adonis(dis_mat, sampleinfox$group, reduce="reference", sim.method = "bray")

 #combine output
 metrics <- data.frame(rr_richness, rr_norm_even, dist_t0= cen_dist2, ab_richness)
 metadata <- cbind(sampleinfo1, metrics, mdsp)
 sampleinfo3 <- metadata %>% filter(harvest == 0)
 sampleinfo3x <- metadata %>% filter(harvest != 0)
 sampleinfo3a <- sampleinfo3 %>% mutate(temperature = 20, exposure = "open")
 sampleinfo3b <- sampleinfo3 %>% mutate(temperature = 5, exposure = "open")
 sampleinfo3c <- sampleinfo3 %>% mutate(temperature = 40, exposure = "closed")
 sampleinfo3d <- sampleinfo3 %>% mutate(temperature = 20, exposure = "closed")
 sampleinfo3e <- sampleinfo3 %>% mutate(temperature = 10, exposure = "closed")
 sampleinfo3f <- sampleinfo3 %>% mutate(temperature = 5, exposure = "closed")
 sampleinfo3g <- sampleinfo3 %>% mutate(temperature = 0, exposure = "closed")
 
 metadata_full <- rbind(sampleinfo3x, sampleinfo3a, sampleinfo3b, sampleinfo3c, sampleinfo3d, sampleinfo3e, sampleinfo3f, sampleinfo3g)
 
 if(use_partial){
   metadata_full <- metadata
 }
 
 richness_statistics <- metadata_full %>% group_by(exposure, temperature) %>% t_test(rr_richness ~ harvest, ref.group = "0", p.adjust.method = "bonferroni")
 evenness_statistics <- metadata_full %>% group_by(exposure, temperature) %>% t_test(rr_norm_even ~ harvest, ref.group = "0", p.adjust.method = "bonferroni")
 distance_statistics <- metadata_full %>% group_by(exposure, temperature) %>% t_test(dist_t0 ~ harvest, ref.group = "0", p.adjust.method = "bonferroni")
 ab_richness_statistics <- metadata_full %>% group_by(exposure, temperature) %>% t_test(ab_richness ~ harvest, ref.group = "0", p.adjust.method = "bonferroni")
 
 result <- list(metadata = metadata, metadata_full = metadata_full, rare_tab = rr_tab, nmds = mds, pair_adonis = pwa, richness_statistics = richness_statistics, evenness_statistics = evenness_statistics, distance_statistics = distance_statistics, ab_richness_statistics = ab_richness_statistics, new_otu_stat = new_data)
 return(result)
}



#### plot main stats related to time - line version ####
overall_line_graphs <- function(prepared_data = xxx, markername = "xxx"){
  prepared_data$metadata <- prepared_data$metadata_full
  df1 <- data_summary(prepared_data$metadata,varname="rr_richness", 
                      groupnames=c("temperature", "harvest", "exposure"))
  
  df1$fTemperature <- factor(df1$temperature)
  df1$fHarvest <- factor(df1$harvest)
  
  p1 <- ggplot(df1, aes(x=fHarvest, y=rr_richness, group=interaction(fTemperature,exposure), color=fTemperature, shape = exposure)) + 
    geom_errorbar(aes(ymin=rr_richness-sem, ymax=rr_richness+sem), width=.1, size = 0.3, 
                  position=position_dodge(0.05)) + scale_shape_manual(values=c(17,1)) +
    geom_line(position=position_dodge(0.05), size = 0.3) + geom_point(position=position_dodge(0.05))+
    scale_color_brewer(type = 'div', palette = 'Spectral', direction = -1, "Temperature") + theme_minimal() + theme_bw() + labs(shape = "Exposure") +
    theme(panel.grid.minor = element_blank(),panel.grid.major = element_blank(),text=element_text(family="Helvetica", face="plain", size=10),
          axis.text = element_text(colour="black"),
          panel.border = element_blank(),
          axis.line = element_line(colour = 'black', size = 0.25),
          strip.background = element_rect(colour="white", fill="white"))+ ylab("OTU richness") + xlab("Storage time (days)") + 
    theme(legend.position = "bottom") + 
    guides(colour = guide_legend(nrow = 1))
  
  df1 <- data_summary(prepared_data$metadata,varname="ab_richness", 
                      groupnames=c("temperature", "harvest", "exposure"))
  
  df1$fTemperature <- factor(df1$temperature)
  df1$fHarvest <- factor(df1$harvest)
  
  p1b <- ggplot(df1, aes(x=fHarvest, y=ab_richness, group=interaction(fTemperature,exposure), color=fTemperature, shape = exposure)) + 
    geom_errorbar(aes(ymin=ab_richness-sem, ymax=ab_richness+sem), width=.1, size = 0.3, 
                  position=position_dodge(0.05)) + scale_shape_manual(values=c(17,1)) +
    geom_line(position=position_dodge(0.05), size = 0.3) + geom_point(position=position_dodge(0.05))+
    scale_color_brewer(type = 'div', palette = 'Spectral', direction = -1, "Temperature") + theme_minimal() + theme_bw() + labs(shape = "Exposure") +
    theme(panel.grid.minor = element_blank(),panel.grid.major = element_blank(),text=element_text(family="Helvetica", face="plain", size=10),
          axis.text = element_text(colour="black"),
          panel.border = element_blank(),
          axis.line = element_line(colour = 'black', size = 0.25),
          strip.background = element_rect(colour="white", fill="white"))+ theme(legend.position = "none") + ylab("Dominant OTU richness") + xlab("Storage time (days)") #+ facet_wrap(.~exposure)
  
  df1 <- data_summary(prepared_data$metadata, varname = "rr_norm_even", 
                      groupnames = c("temperature", "harvest", "exposure"))
  df1$fTemperature <- factor(df1$temperature)
  df1$fHarvest <- factor(df1$harvest)
  
  p4 <- ggplot(df1, aes(x=fHarvest, y=rr_norm_even, group=interaction(fTemperature,exposure), color=fTemperature, shape = exposure)) + 
    geom_errorbar(aes(ymin=rr_norm_even-sem, ymax=rr_norm_even+sem), width=.1, size = 0.3, 
                  position=position_dodge(0.05)) + scale_shape_manual(values=c(17,1)) +
    geom_line(position=position_dodge(0.05), size = 0.3) + geom_point(position=position_dodge(0.05))+
    scale_color_brewer(type = 'div', palette = 'Spectral', direction = -1) + theme_minimal() + theme_bw() + 
    theme(panel.grid.minor = element_blank(),panel.grid.major = element_blank(),text=element_text(family="Helvetica", face="plain", size=10),
          axis.text = element_text(colour="black"),
          panel.border = element_blank(),
          axis.line = element_line(colour = 'black', size = 0.25),
          strip.background = element_rect(colour="white", fill="white")) + ylim(0, 1) + theme(legend.position = "none") + ylab("Pilou's evenness H'/ln(S)") + xlab("Storage time (days)") #+ facet_wrap(.~exposure)
  
  df1 <- data_summary(prepared_data$metadata, varname = "dist_t0", 
                      groupnames = c("temperature", "harvest", "exposure"))
  df1$fTemperature <- factor(df1$temperature)
  df1$fHarvest <- factor(df1$harvest)
  
  p42 <- ggplot(df1, aes(x=fHarvest, y=dist_t0, group=interaction(fTemperature,exposure), color=fTemperature, shape = exposure)) + 
    geom_errorbar(aes(ymin=dist_t0-sem, ymax=dist_t0+sem), width=.1, size = 0.3, 
                  position=position_dodge(0.05)) + scale_shape_manual(values=c(17,1)) +
    geom_line(position=position_dodge(0.05), size = 0.3) + geom_point(position=position_dodge(0.05))+
    scale_color_brewer(type = 'div', palette = 'Spectral', direction = -1) + theme_minimal() + theme_bw() + 
    theme(panel.grid.minor = element_blank(),panel.grid.major = element_blank(),text=element_text(family="Helvetica", face="plain", size=10),
          axis.text = element_text(colour="black"),
          panel.border = element_blank(),
          axis.line = element_line(colour = 'black', size = 0.25),
          strip.background = element_rect(colour="white", fill="white")) + ylim(0, 1) + theme(legend.position = "none") + ylab("BC dissimilarity to t0") + xlab("Storage time (days)") #+ facet_wrap(.~exposure)
  
  legend <- g_legend(p1)
  p1 <- p1 +  theme(legend.position = "none")
  plots <- list(p1, p1b, p4, p42, legend)
  return(plots)
}





#function for plotting taxonomic composition per harvest time. Gives a stacked bar chart and line graphs of change per taxon.
taxonomic_turnover <- function(prepared_data, taxonomy = taxa, cutoff = 0, markername = "bac", temp = NULL, exp = NULL, harv = NULL, useallzero=T, tax_level = 3, num_tax = 20, legend_title = "taxon", plotting = F){
  my_samp <- prepared_data$metadata
  if(useallzero){
    my_samp <- my_samp %>% filter((temperature %in% temp) & (exposure %in% exp) & (harvest %in% harv) | harvest == 0)
  } else {
    my_samp <- my_samp %>% filter((temperature %in% temp) & (exposure %in% exp) & (harvest %in% harv))
  }
  
  my_norm <- prepared_data$norm_table
  my_tab <- prepared_data$rare_tab
  #colnames(my_tab) <- NULL
  my_tab[prepared_data$norm_table<cutoff] <- 0
  my_tab <- my_tab[my_samp$sample,]
  
  #collapse harvest times
  harv_col <- aggregate(x = my_tab, by = list(harvest = my_samp$harvest), FUN = "sum")
  rownames(harv_col) <- harv_col[,1]
  t_harv_col <- t(harv_col)[-1,]
  
  #collapse at taxon level
  harv_tax <- aggregate(x = t_harv_col, by = list(taxon = taxonomy[,tax_level]), FUN = "sum")
  
  #get relative abundances
  harv_tax_rel <- cbind(taxon=harv_tax[,1],sweep(harv_tax[,-1], 2, colSums(harv_tax[,-1]), '/'), stringsAsFactors = F)
  harv_tax_rel_full <- harv_tax_rel
  
  #assign other to low abundant taxa for both tables
  harv_tax_rel$taxon[apply(harv_tax_rel[, -1], 1, max) < sort(apply(harv_tax_rel[, -1], 1, max), decreasing = T)[num_tax]] <- "Other"
  harv_tax$taxon <- harv_tax_rel$taxon
  
  #collapse "other"
  harv_tax <- aggregate(. ~ taxon, data = harv_tax, FUN = sum)
  harv_tax_rel <- aggregate(. ~ taxon, data = harv_tax_rel, FUN = sum)
  
  #make a table with relative to t0
  harv_tax_rel_t0 <- cbind(taxon=harv_tax_rel[,1],sweep(harv_tax_rel[,-1], 1, (0.001+harv_tax_rel[,2]), '/'))
  harv_tax_rel_full_t0 <- cbind(taxon=harv_tax_rel_full[,1],sweep(0.000000000000001+harv_tax_rel_full[,-1], 1, (0.000000000000001+harv_tax_rel_full[,2]), '/'))
  harv_tax_rel_full_t0[,c(2,3)] <- log10(harv_tax_rel_full_t0[,c(2,3)]) #### den her
  
  #gather data
  harv_tax_long <- gather(as.data.frame(harv_tax), harvest, value, -taxon)
  harv_tax_rel_long <- gather(as.data.frame(harv_tax_rel), harvest, rel_value, -taxon)
  harv_tax_rel_t0_long  <- gather(as.data.frame(harv_tax_rel_t0), harvest, rel_t0, -taxon)
 
  harv_tax_rel_full_long <- gather(as.data.frame(harv_tax_rel_full), harvest, rel_value, -taxon)
  harv_tax_rel_full_t0_long <- gather(as.data.frame(harv_tax_rel_full_t0), harvest, rel_t0, -taxon)
  
  #combine
  harv_tax_long$rel_value <- harv_tax_rel_long$rel_value
  harv_tax_long$rel_t0 <- harv_tax_rel_t0_long$rel_t0
  
  #make harvet to factor
  harv_tax_long$harvest <- factor(harv_tax_long$harvest, levels = as.character(harv))
  harv_tax_long$temperature <- temp
  harv_tax_long$exposure <- exp
  harv_tax_long$taxonomic_level <- tax_level
  harv_tax_long$class <- paste0(harv_tax_long$temperature,"°C/",harv_tax_long$exp)
  
  harv_tax_rel_full_t0_long$harvest <- factor(harv_tax_rel_full_t0_long$harvest, levels = as.character(harv))
  harv_tax_rel_full_t0_long$temperature <- temp
  harv_tax_rel_full_t0_long$exposure <- exp
  harv_tax_rel_full_t0_long$taxonomic_level <- tax_level
  harv_tax_rel_full_t0_long$class <- paste0(harv_tax_rel_full_t0_long$temperature,"°C/",harv_tax_rel_full_t0_long$exp)
  
  custom_final28 = c("#771155", "#AA4488", "#EA6CC0", "#CC99BB", "#114477", "#4477AA","#1E78D2", "#77AADD", "#117777", "#44AAAA", "#3FE4E4", "#77CCCC", "#117744","#44AA77", "#1ED278", "#88CCAA", "#771122", "#AA4455", "#D21E2C","#DD7788","#777711", "#AAAA44", "#D2D21E", "#DDDD77","#774411", "#AA7744", "#D2781E", "#DDAA77")
  
  if(plotting){
  custom_final28 = brewer.pal(n = 12, name = "Paired")
  p1 <- ggplot(harv_tax_long,aes(harvest, y=rel_value, fill = factor(taxon))) + geom_bar(position = "fill", stat = "identity") + xlab("Storage time (days)") + ylab("Proportion of total reads") +
    theme_bw() + 
    scale_fill_manual(values= custom_final28, legend_title) +
    theme(panel.grid.minor = element_blank(),panel.grid.major = element_blank(),text=element_text(family="Helvetica", face="plain", size=10),
          axis.text = element_text(colour="black"),
          panel.border = element_blank(),
          axis.line = element_line(colour = 'black', size = 0.25),
          strip.background = element_rect(colour="white", fill="white"))
  
  p2 <- ggplot(harv_tax_long,aes(harvest, y=rel_t0, group = taxon, col = factor(taxon))) + geom_line() + xlab("Storage time (days)") + ylab("Relative change in read abundance") +
    theme_bw() + 
    scale_colour_manual(values= custom_final28, legend_title) +
    theme(panel.grid.minor = element_blank(),panel.grid.major = element_blank(),text=element_text(family="Helvetica", face="plain", size=10),
          axis.text = element_text(colour="black"),
          panel.border = element_blank(),
          axis.line = element_line(colour = 'black', size = 0.25),
          strip.background = element_rect(colour="white", fill="white")) 
  
  } else {
    p1 <- 0
    p2 <- 0
  }
  result <- list(tax_comp_bar = p1,  relative_tax_change = p2, data = harv_tax_long, data_full = harv_tax_rel_full_t0_long)
  return(result)
}


taxonomic_turnover_wrapper <- function(prepared_data=lll, taxonomy = lll, markername = "bact", taxlevels = c("Kingdom","Phylum","Class","Order","Family","Genus"), num_tax = 11){
  cutoffV <- c(0,0,0,0,0,0,0)
  tempV <- c(0,5,5,10,20,20,40)
  expV <- c("closed","closed","open","closed","closed","open","closed")
  short <- c(0,1,7,28)
  long <- c(0,1,7,28,60,120,240,480)
  harvV <- list(short, short, short, short, long, long, short, short)
  useallzeroV <- rep(T, times = 7)
  exclude_timezeroV = rep(T, times = 7)
  inc_legendV <- c(F,F,F,F,T,F,F)
  result <- data.frame(taxon = NA, harvest = NA, value = NA, rel_value = NA, rel_t0 = NA, temperature = NA, exposure = NA, taxonomic_level = NA, class= NA, number = NA, taxon_name= NA)
  result_full <- data.frame(taxon = NA, harvest = NA, rel_t0 = NA, temperature = NA, exposure = NA, taxonomic_level = NA, class= NA, number = NA, taxon_name= NA)
  i=1
  for(t_lev in 1:length(taxlevels)){
    for (x in 1:7){
      current_level <- which(colnames(taxonomy) == taxlevels[t_lev])
      #tp <- taxonomic_turnover(prepared_data, taxonomy = taxonomy, cutoff = cutoffV[x], markername = markername, temp = tempV[x], exp = expV[x], harv = harvV[[x]], useallzero = useallzeroV[x], tax_level = current_level, num_tax = num_tax, legend_title = taxlevels[t_lev], plotting = F)
      ###
      cutoff = cutoffV[x]
      temp = tempV[x]
      exp = expV[x]
      harv = harvV[[x]]
      useallzero = useallzeroV[x]
      tax_level = current_level
      legend_title = taxlevels[t_lev]
      
      my_samp <- prepared_data$metadata
      if(useallzero){
        my_samp <- my_samp %>% filter((temperature %in% temp) & (exposure %in% exp) & (harvest %in% harv) | harvest == 0)
      } else {
        my_samp <- my_samp %>% filter((temperature %in% temp) & (exposure %in% exp) & (harvest %in% harv))
      }
      
      #my_norm <- prepared_data$norm_table
      my_tab <- prepared_data$rare_tab
      #colnames(my_tab) <- NULL
      my_tab[prepared_data$norm_table<cutoff] <- 0
      my_tab <- my_tab[my_samp$sample,]
      
      #collapse harvest times
      harv_col <- aggregate(x = my_tab, by = list(harvest = my_samp$harvest), FUN = "sum")
      rownames(harv_col) <- harv_col[,1]
      t_harv_col <- t(harv_col)[-1,]
      
      #collapse at taxon level
      harv_tax <- aggregate(x = t_harv_col, by = list(taxon = taxonomy[,tax_level]), FUN = "sum")
      
      #get relative abundances
      harv_tax_rel <- cbind(taxon=harv_tax[,1],sweep(harv_tax[,-1], 2, colSums(harv_tax[,-1]), '/'), stringsAsFactors = F)
      harv_tax_rel_full <- harv_tax_rel
      
      #assign other to low abundant taxa for both tables
      harv_tax_rel$taxon[apply(harv_tax_rel[, -1], 1, max) < sort(apply(harv_tax_rel[, -1], 1, max), decreasing = T)[num_tax]] <- "Other"
      harv_tax$taxon <- harv_tax_rel$taxon
      
      #collapse "other"
      harv_tax <- aggregate(. ~ taxon, data = harv_tax, FUN = sum)
      harv_tax_rel <- aggregate(. ~ taxon, data = harv_tax_rel, FUN = sum)
      
      #make a table with relative to t0
      harv_tax_rel_t0 <- cbind(taxon=harv_tax_rel[,1],sweep(harv_tax_rel[,-1], 1, (0.001+harv_tax_rel[,2]), '/'))
      harv_tax_rel_full_t0 <- cbind(taxon=harv_tax_rel_full[,1],sweep(0.000000000000001+harv_tax_rel_full[,-1], 1, (0.000000000000001+harv_tax_rel_full[,2]), '/'))
      harv_tax_rel_full_t0[,c(2,3)] <- log10(harv_tax_rel_full_t0[,c(2,3)]) #### den her
      
      #gather data
      harv_tax_long <- gather(as.data.frame(harv_tax), harvest, value, -taxon)
      harv_tax_rel_long <- gather(as.data.frame(harv_tax_rel), harvest, rel_value, -taxon)
      harv_tax_rel_t0_long  <- gather(as.data.frame(harv_tax_rel_t0), harvest, rel_t0, -taxon)
      
      harv_tax_rel_full_long <- gather(as.data.frame(harv_tax_rel_full), harvest, rel_value, -taxon)
      harv_tax_rel_full_t0_long <- gather(as.data.frame(harv_tax_rel_full_t0), harvest, rel_t0, -taxon)
      
      #combine
      harv_tax_long$rel_value <- harv_tax_rel_long$rel_value
      harv_tax_long$rel_t0 <- harv_tax_rel_t0_long$rel_t0
      
      #make harvet to factor
      harv_tax_long$harvest <- factor(harv_tax_long$harvest, levels = as.character(harv))
      harv_tax_long$temperature <- temp
      harv_tax_long$exposure <- exp
      harv_tax_long$taxonomic_level <- tax_level
      harv_tax_long$class <- paste0(harv_tax_long$temperature,"°C/",harv_tax_long$exp)
      
      harv_tax_rel_full_t0_long$harvest <- factor(harv_tax_rel_full_t0_long$harvest, levels = as.character(harv))
      harv_tax_rel_full_t0_long$temperature <- temp
      harv_tax_rel_full_t0_long$exposure <- exp
      harv_tax_rel_full_t0_long$taxonomic_level <- tax_level
      harv_tax_rel_full_t0_long$class <- paste0(harv_tax_rel_full_t0_long$temperature,"°C/",harv_tax_rel_full_t0_long$exp)
      
      tp <- list(data = harv_tax_long, data_full = harv_tax_rel_full_t0_long)
      
      ###
      tp$data$number <- i
      tp$data$taxon_name <- taxlevels[t_lev]
      result <- rbind(result, tp$data)
      
      tp$data_full$number <- i
      tp$data_full$taxon_name <- taxlevels[t_lev]
      i <- i+1
      result_full <- rbind(result_full, tp$data_full)

    }
  }
  result <- result[-1,]
  result_full <- result_full[-1,]
  result$harvest <- factor(result$harvest, levels = as.character(long))
  result_full$harvest <- factor(result_full$harvest, levels = as.character(long))
  res <- list(result, result_full)
  return(res)
}

#### Function for plotting above function across full dataset ####
taxonomic_turnover_plotter <- function(tax_data = yy, tax_levels_used = c("Kingdom","Phylum","Class","Order","Family","Genus"), legend_columns = 2){
  plots <- list()
  for (i in 1:length(tax_levels_used)){
    tax_lev_sub <- tax_data[[1]] %>% filter(taxon_name == tax_levels_used[i])
    colourCount = length(unique(tax_lev_sub$taxon))
    getPalette = colorRampPalette(brewer.pal(8, "Set2"))
    plots[[i]] <-  tax_lev_sub %>% ggplot(aes(harvest, y=rel_value, fill = factor(taxon))) + 
      geom_bar(position = "fill", stat = "identity") + 
      xlab("Storage time (days)") + 
      ylab("Proportion of total reads") +
      theme_bw() + 
      facet_grid(. ~ class,scale="free_x", space = "free_x") + 
      scale_fill_manual(values= getPalette(colourCount), tax_levels_used[i]) +
      guides(fill=guide_legend(ncol=legend_columns)) + 
      theme(panel.grid.minor = element_blank(),panel.grid.major = element_blank(),text=element_text(family="Helvetica", face="plain", size=7),
            axis.text = element_text(colour="black"),
            panel.border = element_blank(),
            axis.line = element_line(colour = 'black', size = 0.25),
            strip.background = element_rect(colour="white", fill="white"), legend.key.size = unit(0.4, "cm"))
  }
  return(plots)
}

#delete?
#### plot nmds results ####
plotmds <- function(data=prepared_data){
  p1 <- data$metadata %>% 
  ggplot(aes(MDS1,MDS2, shape = exposure, col = factor(temperature))) + 
  geom_point() + 
  facet_wrap(.~harvest, nrow = 2) + 
  scale_shape_manual(values=c(17,1)) +
  guides(color=guide_legend(title="Temperature"), shape = guide_legend(title = "Exposure")) + 
    scale_color_brewer(palette="Set1") +
  xlab("MDS1") + 
  ylab("MDS2") +
  theme_bw() + 
  theme(panel.grid.minor = element_blank(),panel.grid.major = element_blank(),text=element_text(size=10),
        axis.text = element_text(colour="black"),
        #panel.border = element_blank(),
        axis.line = element_line(colour = 'black', size = 0.25),
        strip.background = element_rect(colour="white", fill="white"))
  
  p1
}

#### plot nmds results version 2####
plotmds2 <- function(data=prepared_data){
  data$metadata$temperature <- paste0(data$metadata$temperature,"°C")
  data$metadata$temperature <- factor(data$metadata$temperature, levels = c("0°C","5°C","10°C","20°C","40°C"))
  p1 <- data$metadata %>% 
    ggplot(aes(MDS1,MDS2, shape = exposure, col = as.factor(harvest))) + 
    geom_point(alpha = 0.8, size = 3) + 
    facet_wrap(.~temperature, nrow = 1) + 
    scale_shape_manual(values=c(16,1)) +
    guides(color=guide_legend(title="Storage time", nrow = 1), shape = guide_legend(title = "Exposure", nrow = 1)) + 
    scale_color_brewer(type = 'div', palette = 'Spectral', direction = -1) +
    xlab("NMDS1") + 
    ylab("NMDS2") +
    theme_bw() +
    theme(panel.grid.minor = element_blank(),panel.grid.major = element_blank(),text=element_text(size=10),
          axis.text = element_text(colour="black"),
          #panel.border = element_blank(),
          axis.line = element_line(colour = 'black', size = 0.25),
          strip.background = element_rect(colour="white", fill="white"), legend.position = "bottom")
  
  p1
}

###
plot_model_mds <- function(nmdsx=nmdsx, habitats_to_test = habitats_to_test, mdata = mdata ){
  plots <- list()
  i <- 1
  for(hab_no in 1:length(habitats_to_test)){
    nm_data <- as.data.frame(nmdsx$points)
    habitat <- c(ifelse(mdata[,habitats_to_test[hab_no]] == 1, habitats_to_test[hab_no], paste0("NOT ",habitats_to_test[hab_no])))
    hn <- names(table(habitat))
    hn <- hn[order(nchar(hn))]
    bw_info <- data.frame(sample = seq(1,130), no = seq(1,130), exposure = c(rep("Ref",80),"Ref - same site",rep("Ref",49)), temperature = habitat, harvest = rep("Ref",130), batch = "Ref")
    nm_data2 <- cbind(rbind(bw_info, sampleinfo),nm_data)
    nm_data2$temperature[nchar(nm_data2$temperature)<4] <- paste0(nm_data2$temperature[nchar(nm_data2$temperature)<4],"°C")
    nm_data2$temperature <- factor(nm_data2$temperature, levels = c("0°C","5°C","10°C","20°C","40°C", hn))
    
    p1 <- nm_data2 %>% 
      ggplot(aes(MDS1,MDS2, shape = exposure, col = temperature)) + 
      geom_point(alpha = 0.9, size = 2) + 
      #facet_wrap(.~temperature, nrow = 1) + 
      scale_shape_manual(values=c(16,1,8, 9)) +
      guides(color=guide_legend(title="", nrow = 1), shape = guide_legend(title = ""), nrow = 2) + 
      scale_color_manual(values=c("#2B83BA", "#ABDDA4", "#FFFFBF", "#FDAE61", "#D7191C","#000000","#AAAAAA")) + 
      #scale_color_brewer(type = 'div', palette = 'Spectral', direction = -1) +
      xlab("MDS1") + 
      ylab("MDS2") +
      theme_bw() +
      theme(panel.grid.minor = element_blank(),panel.grid.major = element_blank(),text=element_text(size=10),
            axis.text = element_text(colour="black"),
            #panel.border = element_blank(),
            axis.line = element_line(colour = 'black', size = 0.25),
            strip.background = element_rect(colour="white", fill="white"), legend.position = "bottom", legend.box="vertical", legend.margin=margin())
    
    plots[[i]] <- p1
    i <- i+1
    
    p1 <- nm_data2 %>% 
      ggplot(aes(MDS3,MDS4, shape = exposure, col = temperature)) + 
      geom_point(alpha = 0.9, size = 2) + 
      #facet_wrap(.~temperature, nrow = 1) + 
      scale_shape_manual(values=c(16,1,8, 9)) +
      guides(color=guide_legend(title="", nrow = 1), shape = guide_legend(title = ""), nrow = 2) + 
      scale_color_manual(values=c("#2B83BA", "#ABDDA4", "#FFFFBF", "#FDAE61", "#D7191C","#000000","#AAAAAA")) + 
      #scale_color_brewer(type = 'div', palette = 'Spectral', direction = -1) +
      xlab("MDS3") + 
      ylab("MDS4") +
      theme_bw() +
      theme(panel.grid.minor = element_blank(),panel.grid.major = element_blank(),text=element_text(size=10),
            axis.text = element_text(colour="black"),
            #panel.border = element_blank(),
            axis.line = element_line(colour = 'black', size = 0.25),
            strip.background = element_rect(colour="white", fill="white"), legend.position = "bottom", legend.box="vertical", legend.margin=margin())
    
    plots[[i]] <- p1
    i <- i+1
  }
  return(plots)
}

#### plot relative distance facet plot
plot_ratio_facet <- function(dist_obj = vdm, prep = prep, sampleinfo = sampleinfo){
  require(tidyverse)
  #
  mmb2 <- dist_obj[grepl("SE", rownames(dist_obj)), !grepl("SE", colnames(dist_obj))]
  ns081 <- mmb2[,"SN081"]
  mmb3 <- mmb2[,colnames(mmb2) != "SN081"]
  #minx <- apply(mmb3, 1, FUN=min)
  
  rest <- gather(as.data.frame(mmb3))
  rest$site <- rownames(mmb3)
  ddd <- data.frame(sample=names(ns081), rest = rest$value, ref_id = rest$key)
  sampleinfo$ns081 <- ns081
  #sampleinfo$minx <- minx
  ggg <- right_join(sampleinfo, ddd, by = "sample")
  
  prep2 <- prep$metadata[,c("sample","dist_t0")]
  
  lll <- left_join(ggg,prep2, by = "sample")
  lll$class <- paste0(lll$temperature,"°C/",lll$exposure)
  lll$class <- factor(lll$class, levels = c("0°C/closed", "5°C/closed", "5°C/open",  "10°C/closed", "20°C/closed", "20°C/open", "40°C/closed"))
  
  miny <- min(c(min(lll$rest/lll$dist_t0),0.95))
  p1 <- ggplot(lll,(aes(factor(harvest),rest/dist_t0, col = factor(harvest), shape = exposure, group = exposure))) +
    geom_jitter(size=1, alpha = 0.8) + 
    facet_grid(.~class, scales = "free", space = "free_x") + 
    geom_abline(intercept = 1, slope = 0, size = 0.3, linetype = 2) + ylim(0.9,5.2) + #scale_y_continuous(breaks=c(1,2,3,4,5)) +
    scale_shape_manual(values=c(17,1)) +
    scale_color_brewer(type = 'div', palette = 'Spectral', direction = -1) + 
    guides(color=guide_legend(title="Storage time", nrow = 1), shape = guide_legend(title = "Exposure", nrow = 1)) +
    theme_minimal() + theme_bw() + 
    theme(panel.grid.minor = element_blank(), panel.grid.major = element_blank(),text=element_text(family="Helvetica", face="plain", size=10),
          axis.text = element_text(colour="black"),
          #panel.border = element_blank(),
          axis.line = element_line(colour = 'black', size = 0.25),
          strip.background = element_rect(colour="white", fill="white"), legend.position = "bottom") + xlab("Storage time (days)") + ylab("Dissim. ratio (other site / t0)")
  
  theme(legend.position = "bottom") + 
    guides(colour = guide_legend(nrow = 1))
  
  miny <- min(c(min(lll$rest/lll$ns081),0.95))
  maxy <- max(c(min(lll$rest/lll$ns081),2.22))
  p2 <- ggplot(lll,(aes(factor(harvest),rest/ns081, col = factor(harvest), shape = exposure, group = exposure))) +
    geom_jitter(size=1, alpha = 0.8) + 
    facet_grid(.~class, scales = "free", space = "free_x") + 
    geom_abline(intercept = 1, slope = 0, size = 0.3, linetype = 2) + ylim(0.9,2.4) +# scale_y_continuous(breaks=c(1,2)) + 
    scale_shape_manual(values=c(17,1)) +
    scale_color_brewer(type = 'div', palette = 'Spectral', direction = -1) + theme_minimal() + theme_bw() + 
    guides(color=guide_legend(title="Storage time", nrow = 1), shape = guide_legend(title = "Exposure", nrow = 1)) +
    theme(panel.grid.minor = element_blank(), 
          panel.grid.major = element_blank(),
          text=element_text(family="Helvetica", face="plain", size=10),
          axis.text = element_text(colour="black"),
          #panel.border = element_blank(),
          axis.line = element_line(colour = 'black', size = 0.25),
          strip.background = element_rect(colour="white", fill="white"),
          legend.position = "bottom") + xlab("Storage time (days)") + 
    ylab("Dissim. ratio (other site / same)")
  
  p3 <- ggplot(lll,(aes(factor(harvest),rest, col = factor(harvest), shape = exposure, group = exposure))) +
    geom_jitter(size=1, alpha = 0.8) + 
    facet_grid(.~class, scales = "free", space = "free_x") + 
    ylim(0.5,1.1) +# scale_y_continuous(breaks=c(1,2)) + 
    scale_shape_manual(values=c(17,1)) +
    scale_color_brewer(type = 'div', palette = 'Spectral', direction = -1) + theme_minimal() + theme_bw() + 
    guides(color=guide_legend(title="Storage time", nrow = 1), shape = guide_legend(title = "Exposure", nrow = 1)) +
    theme(panel.grid.minor = element_blank(), 
          panel.grid.major = element_blank(),
          text=element_text(family="Helvetica", face="plain", size=10),
          axis.text = element_text(colour="black"),
          #panel.border = element_blank(),
          axis.line = element_line(colour = 'black', size = 0.25),
          strip.background = element_rect(colour="white", fill="white"),
          legend.position = "bottom") + xlab("Storage time (days)") + 
    ylab("BC dissim. to other sites")
  
  result <- list(p1, p2, p3)
  return(result)
  
}


#### NOT USED FROM HERE ----

#### prepare data with reference data included ####
prepare_data_ref <- function(community_table = table, sam_info = info, data_set = "unnamed", bw_between = NA){
  require(vegan)
  require(tidyr)
  
  rawtab <- community_table
  
  #synchronize content and sequence of table/list
  target <- tibble(sample = rownames(rawtab))
  sampleinfo1 <- sam_info %>% inner_join(target, by = "sample")
  rawtab <- rawtab[sampleinfo1$sample,]
  
  #richness & seq depth
  raw_richness <- rowSums(rawtab>0)
  seq_depth <- rowSums(rawtab)
  
  #rarefy and normalize
  rr_tab <- rrarefy(rawtab, quantile(rowSums(rawtab))[2])
  rr_tab[rr_tab<0] <- 0
  rr_tab_norm <- sweep(rr_tab, 1, rowSums(rr_tab), '/')
  
  #rarefied richness and shannon
  rr_richness <- rowSums(rr_tab>0)
  rr_norm_shannon <- diversity(rr_tab_norm)
  
  #distance metrics
  dis_mat <- as.matrix(vegdist(decostand(rr_tab_norm, method = "hellinger"), method = "bray"))
  sim_mat <- 1-dis_mat
  ddd <- dist_groups(sim_mat, sampleinfo1$dataset)
  ddd$Label <- as.character(ddd$Label)
  ddd <- ddd[!ddd$Label == "Within treatment",]
  ddd$Label[ddd$Label == "Between SN081 and treatment"] <- "Correct site"
  ddd$Label[grepl("and treatment",ddd$Label)] <- "Incorrect site"
  ddd$Label[grepl("Within",ddd$Label)] <- "Within site"
  ddd$Label[grepl("Between",ddd$Label)] <- "Between site"
  eee <- ddd %>% filter(Label %in% c("Correct site","Incorrect site")) %>% group_by(Item1, Label) %>% arrange(-Distance) %>% slice(1) %>% select(Item1,Label,Distance) %>% spread(key=Label, value = Distance)
  names(eee) <- c("sample","max_sim_correct","max_sim_incorrect")
  dis_data <- left_join(ddd, sampleinfo1, by = c("Item1"="sample"))
  dis_data <- left_join(dis_data, eee, by = c("Item1"="sample"))
  ###
  
  between_site <- ddd$Distance[ddd$Label == "Between site"]
  within_site <- ddd$Distance[ddd$Label == "Within site"]
  if(!is.na(bw_between)){
    between_site <- 1-as.vector(bw_between)
  }
  
  emp_w <-  ecdf(within_site)
  emp_b <- ecdf(between_site)
  dens_w <- CDF(density(within_site))
  dens_b <- CDF(density(between_site))
  
  dis_data$emp_w <- emp_w(dis_data$Distance)
  dis_data$emp_b <- 1-emp_b(dis_data$Distance)
  dis_data$dens_w <- dens_w(dis_data$Distance)
  dis_data$dens_b <- 1-dens_b(dis_data$Distance)
  
  
  mds <- metaMDS(dis_mat, distance = "bray", k=2, try = 500, trymax = 4000)
  mdsp <- mds$points[,1:2]
  names(mdsp) <- c("NMDS1","NMDS2")
  
  #combine output
  metadata <- cbind(sampleinfo1,mdsp)
  result <- list(metadata=metadata, matrix = dis_mat, distance_data = dis_data, emp_w=emp_w, emp_b=emp_b, dens_w=dens_w, dens_b=dens_b, within_site=within_site, between_site=between_site)
  return(result)
}

plot_model_outputs <- function(prepared_data = NULL){
  d2 <- prepared_data$distance_data
  d2$harvest <- factor(d2$harvest, levels = as.character(c(0,1,7,28,60,120,240,480)))
  d2$class <- paste0(d2$temperature,"°C/",d2$exposure)
  d2$class <- factor(d2$class, levels = c("0°C/closed", "5°C/closed", "5°C/open",  "10°C/closed", "20°C/closed", "20°C/open", "40°C/closed"))
  d2 <- d2 %>% filter(Label %in% c("Correct site","Incorrect site"))
  d2$prop <- log10((d2$emp_w+0.00001)/(d2$emp_b+0.00001))
  d2$prop_d <- log10((d2$dens_w+0.00001)/(d2$dens_b+0.00001))
  
sim_plot <- ggplot(d2,aes(factor(harvest), Distance, col = Label)) + geom_jitter(width = 0.2, size = 1, alpha=0.5, shape = 21) + 
  facet_grid(.~class, scales = "free", space = "free") +
  theme_bw() + scale_color_brewer(palette = "Set2", name = "Similarity with") + 
  theme(panel.grid.minor = element_blank(),panel.grid.major = element_blank(),text=element_text(size=10),
        axis.text = element_text(colour="black"),
        panel.border = element_blank(),
        axis.line = element_line(colour = 'black', size = 0.25),
        strip.background = element_rect(colour="white", fill="white")) + xlab("Storage time (days)") + ylab("Similarity (1-BC)")

#ggplot(d2,aes(factor(harvest), emp_w, col = Label)) + geom_jitter(width = 0.2, size = 1, alpha=0.5, shape = 21) + 
#  facet_grid(.~class, scales = "free", space = "free") +
#  theme_bw() + scale_color_brewer(palette = "Set2", name = "Compared with sample from") + 
#  theme(panel.grid.minor = element_blank(),panel.grid.major = element_blank(),text=element_text(size=10),
#        axis.text = element_text(colour="black"),
#        panel.border = element_blank(),
#        axis.line = element_line(colour = 'black', size = 0.25),
#        strip.background = element_rect(colour="white", fill="white")) + xlab("Storage time (days)") + ylab("Probability")

#ggplot(d2,aes(factor(harvest), dens_w, col = Label)) +  geom_jitter(width = 0.2, size = 1, alpha=0.5, shape = 21) + 
#  facet_grid(.~class, scales = "free", space = "free") +
#  theme_bw() + scale_color_brewer(palette = "Set2", name = "Compared with sample from") + 
#  theme(panel.grid.minor = element_blank(),panel.grid.major = element_blank(),text=element_text(size=10),
#        axis.text = element_text(colour="black"),
#        panel.border = element_blank(),
#        axis.line = element_line(colour = 'black', size = 0.25),
#        strip.background = element_rect(colour="white", fill="white")) + xlab("Storage time (days)") + ylab("Probability")

ratio_plot <- ggplot(d2,aes(factor(harvest), prop, col = Label)) +  geom_jitter(width = 0.2, size = 1, alpha=0.5, shape = 21) + 
  facet_grid(.~class, scales = "free", space = "free") +
  theme_bw() + scale_color_brewer(palette = "Set2", name = "Compared with sample from") + 
  theme(panel.grid.minor = element_blank(),panel.grid.major = element_blank(),text=element_text(size=10),
        axis.text = element_text(colour="black"),
        panel.border = element_blank(),
        axis.line = element_line(colour = 'black', size = 0.25),
        strip.background = element_rect(colour="white", fill="white")) + xlab("Storage time (days)") + ylab("log10 probability ratio ")

ratio_plot_d <- ggplot(d2,aes(factor(harvest), prop_d, col = Label)) +  geom_jitter(width = 0.2, size = 1, alpha=0.5, shape = 21) + 
  facet_grid(.~class, scales = "free", space = "free") +
  theme_bw() + scale_color_brewer(palette = "Set2", name = "Compared with sample from") + 
  theme(panel.grid.minor = element_blank(),panel.grid.major = element_blank(),text=element_text(size=10),
        axis.text = element_text(colour="black"),
        panel.border = element_blank(),
        axis.line = element_line(colour = 'black', size = 0.25),
        strip.background = element_rect(colour="white", fill="white")) + xlab("Storage time (days)") + ylab("log10 probability ratio ")

ratio_plot_f <- ggplot(d2,aes(factor(harvest), max_sim_correct/max_sim_incorrect)) + geom_point(size = 1, alpha=1, shape = 21) + 
  facet_grid(.~class, scales = "free", space = "free") +
  theme_bw() + scale_shape_manual(name = "Ratio") + 
  theme(panel.grid.minor = element_blank(),panel.grid.major = element_blank(),text=element_text(size=10),
        axis.text = element_text(colour="black"),
        panel.border = element_blank(),
        axis.line = element_line(colour = 'black', size = 0.25),
        strip.background = element_rect(colour="white", fill="white")) + xlab("Storage time (days)") + ylab("ratio H0/H1")

#model_data_plot <- prepared_data$distance_data %>% filter(!Label %in% c("Correct site","Incorrect site")) %>% 
#  ggplot(aes(Distance, col = Label, fill = Label)) + 
#  geom_histogram(aes(y=..density..),binwidth=.005, alpha=1, position="dodge", size=.1, col = "black") + 
#  geom_density(alpha=0, show.legend=FALSE) + 
#  scale_fill_brewer(palette = "Set2", name = "Pairwise comparison", direction = -1) + 
#  scale_colour_brewer(palette = "Set2", direction = -1) + 
#  theme_bw() +
#  theme(panel.grid.minor = element_blank(),panel.grid.major = element_blank(),text=element_text(size=10),
#       axis.text = element_text(colour="black"),
#        panel.border = element_blank(),
#        axis.line = element_line(colour = 'black', size = 0.25),
#       strip.background = element_rect(colour="white", fill="white")) + xlab("Similarity (1-BC)") + ylab("Density ")

#version on the reference model distribition
modeldata <- data.frame(value = c(prepared_data$within_site, prepared_data$between_site), type = c(rep("within site", length(prepared_data$within_site)), rep("between site",length(prepared_data$between_site))))

model_data_plot <- ggplot(modeldata, aes(value, col = type, fill = type)) + 
  geom_histogram(aes(y=..density..),binwidth=.005, alpha=1, position="dodge", size=.1, col = "black") + 
  geom_density(alpha=0, show.legend=FALSE) + 
  scale_fill_brewer(palette = "Set2", name = "Pairwise comparison", direction = -1) + 
  scale_colour_brewer(palette = "Set2", direction = -1) + 
  theme_bw() +
  theme(panel.grid.minor = element_blank(),panel.grid.major = element_blank(),text=element_text(size=10),
        axis.text = element_text(colour="black"),
        panel.border = element_blank(),
        axis.line = element_line(colour = 'black', size = 0.25),
        strip.background = element_rect(colour="white", fill="white")) + xlab("Similarity (1-BC)") + ylab("Density ")


data2 <- prepared_data$metadata

every_facet_data = subset(data2, dataset != "treatment")
every_facet_data$temperature <- "reference"
every_facet_data$exposure <- as.character(every_facet_data$exposure)
every_facet_data$exposure[every_facet_data$dataset == "SN081"] <- "Ref. same site"
every_facet_data$exposure[every_facet_data$dataset != "SN081"] <- "Ref. other site"
every_facet_data$harvest <- 0
individual_facet_data = subset(data2, dataset == "treatment")
individual_facet_data$facet = individual_facet_data$temperature
individual_facet_data$exposure <- as.character(individual_facet_data$exposure)
every_facet_data = merge(every_facet_data, data.frame(facet = unique(individual_facet_data$facet)), all = T)
plot_data = rbind(every_facet_data, individual_facet_data)

plot_data$facet <- paste0(plot_data$facet,"°C")
plot_data$facet <- factor(plot_data$facet, levels = c("0°C","5°C","10°C","20°C","40°C"))

nmds_ref_plot <- plot_data %>%
  ggplot(aes(MDS1,MDS2, shape = exposure, col = as.factor(harvest))) + 
  geom_point(alpha = 0.7, size = 2) + 
  facet_wrap(.~facet, nrow = 1) + 
  scale_shape_manual(values=c(16,1,0,7),breaks=c("Ref. same site", "Ref. other site", "closed", "open"),labels=c("Ref. same site", "Ref. other site", "Exp. closed", "Exp. open")) +
  guides(color=guide_legend(title="Storage time"), shape = guide_legend(title = "Sample")) + 
  scale_color_brewer(type = 'div', palette = 'Spectral', direction = -1) +
  xlab("MDS1") + 
  ylab("MDS2") +
  theme_bw() + 
  theme(panel.grid.minor = element_blank(),panel.grid.major = element_blank(),text=element_text(size=10),
        axis.text = element_text(colour="black"),
        #panel.border = element_blank(),
        axis.line = element_line(colour = 'black', size = 0.25),
        strip.background = element_rect(colour="white", fill="white"))

nmds_ref_plot_one_row <- nmds_ref_plot + guides(color=guide_legend(title="Storage time", nrow = 1), shape = guide_legend(title = "Sample", nrow = 1)) + theme(legend.position = "bottom")

result <- list(similarity_plot = sim_plot, ratio_plot = ratio_plot, model_data_plot = model_data_plot, nmds_ref_plot = nmds_ref_plot, ratio_plot_d=ratio_plot_d, ratio_plot_f=ratio_plot_f, nmds_ref_plot_one_row=nmds_ref_plot_one_row)
return(result)
}


plot_new_models <- function(sampleinfo = sampleinfo, knn_data = knn_data){
  dist_data <- cbind(sampleinfo, knn_data$distance_to_centroid)
  mean_dist_ref <- mean(knn_data$distance_to_centroid_references$CentroidDistance)
  sd_dist_ref <- sd(knn_data$distance_to_centroid_references$CentroidDistance)
  
  df1 <- data_summary(dist_data,varname="CentroidDistance", 
                      groupnames=c("temperature", "harvest", "exposure"))
  
  df1$fTemperature <- factor(df1$temperature)
  df1$fHarvest <- factor(df1$harvest)
  
  p1 <- ggplot(df1, aes(x=fHarvest, y=CentroidDistance, group=interaction(fTemperature,exposure), color=fTemperature, shape = exposure)) + 
    geom_errorbar(aes(ymin=CentroidDistance-sem, ymax=CentroidDistance+sem), width=.1, size = 0.3, 
                  position=position_dodge(0.05)) + scale_shape_manual(values=c(17,1)) +
    geom_line(position=position_dodge(0.05), size = 0.3) + geom_point(position=position_dodge(0.05))+
    scale_color_brewer(type = 'div', palette = 'Spectral', direction = -1, "Temperature") + theme_minimal() + theme_bw() + labs(shape = "Exposure") +
    theme(panel.grid.minor = element_blank(),panel.grid.major = element_blank(),text=element_text(family="Helvetica", face="plain", size=10),
          axis.text = element_text(colour="black"),
          panel.border = element_blank(),
          axis.line = element_line(colour = 'black', size = 0.25),
          strip.background = element_rect(colour="white", fill="white"))+ ylab("BC dissim. to habitat cluster centroid") + xlab("Storage time (days)") +
    geom_hline(yintercept = mean_dist_ref+0*sd_dist_ref, linetype = "longdash", size = 0.3) + geom_hline(yintercept = mean_dist_ref+1*sd_dist_ref, linetype = "dashed", size = 0.3) + geom_hline(yintercept = mean_dist_ref+2*sd_dist_ref, linetype = "dotted", size = 0.3 ) + ylim(0.3, 0.9) + theme(legend.position = "bottom") + guides(colour = guide_legend(nrow = 1))

  
  sampleinfo$classification_success <- knn_data$biostratum_knn[,knn_data$correct_cluster]
  
  df1 <- data_summary(sampleinfo,varname="classification_success", 
                      groupnames=c("temperature", "harvest", "exposure"))
  
  df1$fTemperature <- factor(df1$temperature)
  df1$fHarvest <- factor(df1$harvest)
  
  p2 <- df1 %>% ggplot(aes(as.factor(harvest),as.factor(temperature), fill = classification_success, label = round(classification_success,2))) + geom_tile(color="black", size = 0.2) + 
    scale_fill_gradientn(colours = c("red", "white", "green"), limits=c(0,1), name = "Mean probability (correct habitat cluster)") + theme_classic() + ylab("Temperature °C") + xlab("Storage time (days)") + facet_grid(exposure ~ ., scales = "free", space = "free_y")  + geom_text(size = 2) + theme(legend.position = "bottom") + guides(colour = guide_legend(nrow = 1)) +
    theme(panel.grid.minor = element_blank(),
          panel.grid.major = element_blank(),
          text=element_text(family="Helvetica", face="plain", size=10),
          axis.text = element_text(colour="black"),
          panel.border = element_blank(),
          axis.line = element_line(colour = 'black', size = 0.25),
          #strip.background = element_rect(colour="white", fill="white")
    )

  sampleinfo$classification_success <- knn_data$biostratum_knn[,"LateDryRich"]
  
  df1 <- data_summary(sampleinfo,varname="classification_success", 
                      groupnames=c("temperature", "harvest", "exposure"))
  
  df1$fTemperature <- factor(df1$temperature)
  df1$fHarvest <- factor(df1$harvest)
  df1[,4:6][is.na(df1[,4:6])] <- 0
  p3 <- df1 %>% ggplot(aes(as.factor(harvest),as.factor(temperature), fill = classification_success, label = round(classification_success,2))) + geom_tile(color="black", size = 0.2) + 
    scale_fill_gradientn(colours = c("red", "white", "green"), limits=c(0,1), name = "Mean probability (correct habitat cluster)") + theme_classic() + ylab("Temperature °C") + xlab("Storage time (days)") + facet_grid(exposure ~ ., scales = "free", space = "free_y")  + geom_text(size = 2) + theme(legend.position = "bottom") + guides(colour = guide_legend(nrow = 1))+
    theme(panel.grid.minor = element_blank(),
          panel.grid.major = element_blank(),
          text=element_text(family="Helvetica", face="plain", size=10),
          axis.text = element_text(colour="black"),
          panel.border = element_blank(),
          axis.line = element_line(colour = 'black', size = 0.25),
          #strip.background = element_rect(colour="white", fill="white")
    )
  
  result <- list(p1,p2,p3)
  return(result)
}


plot_knn_site <- function(knn_data=NULL, sampleinfo = NULL){
  sampleinfo$classification_success <- knn_data

df1 <- data_summary(sampleinfo,varname="classification_success", 
                    groupnames=c("temperature", "harvest", "exposure"))

df1$fTemperature <- factor(df1$temperature)
df1$fHarvest <- factor(df1$harvest)

p1 <- df1 %>% ggplot(aes(as.factor(harvest),as.factor(temperature), fill = classification_success, label = round(classification_success,2))) + geom_tile(color="black", size = 0.2) + 
  scale_fill_gradientn(colours = c("red", "white", "green"), limits=c(0,1), name = "Mean probability (correct locality)") + 
  theme_classic() + 
  ylab("Temperature °C") + 
  xlab("Storage time (days)") + 
  facet_grid(exposure ~ ., scales = "free", space = "free_y")  + 
  geom_text(size = 2) + 
  theme(legend.position = "bottom") + 
  guides(colour = guide_legend(nrow = 1)) + 
  theme(panel.grid.minor = element_blank(),
        panel.grid.major = element_blank(),
        text=element_text(family="Helvetica", face="plain", size=10),
        axis.text = element_text(colour="black"),
        panel.border = element_blank(),
        axis.line = element_line(colour = 'black', size = 0.25),
        #strip.background = element_rect(colour="white", fill="white")
        )


return(p1)
}

#### rarefaction and richness correlation plots ----
rare_rich_calc <- function(dataset=bactX, dataname = "My_data"){
 out <- rarecurve(dataset, step = 100, sample = 100, label = FALSE)
 names(out) <- paste("species", 1:107, sep = "")
 protox <- mapply(FUN = function(x, y) {
  mydf <- as.data.frame(x)
  colnames(mydf) <- "value"
  mydf$species <- y
  mydf$subsample <- attr(x, "Subsample")
  mydf
 }, x = out, y = as.list(names(out)), SIMPLIFY = FALSE)
 
 xy <- do.call(rbind, protox)
 rownames(xy) <- NULL  # pretty
 xy$organism <- dataname
 
 rarefy_level <- quantile(rowSums(dataset))[2]
 raretab <- rrarefy(dataset, rarefy_level)
 rr_richness <- rowSums(raretab>0)
 raw_richness <- rowSums(dataset>0)
 richtab <- data.frame(raw_richness = raw_richness, resampled_richness = rr_richness)
 result <- list(rarefaction_curve = xy, rarelevel = rarefy_level, richness = richtab, dataset = dataname)
 return(result)
}

plot_rrare <- function(dataset = NA, dataname = ""){
 p1 <- ggplot(dataset$rarefaction_curve, aes(x=subsample, y=value, group = species)) + 
  geom_line(size = 0.3) +
  theme_minimal() + theme_bw() +
  theme(panel.grid.minor = element_blank(),panel.grid.major = element_blank(),text=element_text(family="Helvetica", face="plain", size=10),
        axis.text = element_text(colour="black"),
        panel.border = element_blank(),
        axis.line = element_line(colour = 'black', size = 0.25),
        strip.background = element_rect(colour="white", fill="white"))+ ylab("ASVs") + xlab("Sequencing depth") + ggtitle(dataname) + geom_vline(xintercept = dataset$rarelevel, linetype = 2)
 return(p1)
}

plot_rare_rich <- function(dataset = NA, dataname = ""){
 p1 <- ggplot(dataset$richness, aes(x=raw_richness, y=resampled_richness)) + 
  geom_point(size = 0.3) +
  theme_minimal() + theme_bw() +
  theme(panel.grid.minor = element_blank(),panel.grid.major = element_blank(),text=element_text(family="Helvetica", face="plain", size=10),
        axis.text = element_text(colour="black"),
        panel.border = element_blank(),
        axis.line = element_line(colour = 'black', size = 0.25),
        strip.background = element_rect(colour="white", fill="white"))+ ylab("Resamp ASV richness") + xlab("Raw ASV richness") + ggtitle(dataname)
 return(p1)
}
