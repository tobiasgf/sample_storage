#### preparation of basic datasets ----

# this script expects that all the reads have been processed with dada2 and that asv tables are available in the folder "asv_tables"

#### loading of libraries ####
library(tidyverse)
library(here)
source(here::here("R","scripts.R"))
library(usedist)
library(pairwiseAdonis)
library(rstatix)
library(BiodiversityR)
library(vegan)
library(lemon)
library(cowplot)
library(dada2)
library(RColorBrewer)
library(ggpubr)
library(caret)
library(data.table)

#### Fungi data preparation ----
# reading the dada2 otutables from the three MiSeqruns (run1: short storage, run2: long storage, run3: reruns of failed samples)
org <- readRDS(here::here("asv_tables","fun1_dada2_otutable.rds"))
lts <- readRDS(here::here("asv_tables","fun2_dada2_otutable.rds"))
rerun <- readRDS(here::here("asv_tables","fun3_dada2_otutable.rds"))
sampleinfo <- read_csv(here::here("basic_data","sample_data.csv"))

#combine tables (storage only)
tab1 <- mergeSequenceTables(org,lts, rerun, repeats = "sum")
saveRDS(tab1, here::here("asv_tables","fungi_combined_otutable.rds"))

#quick check of clustering
tab1 <- readRDS(here::here("asv_tables","fungi_combined_otutable.rds"))
rare_tab <- rrarefy(tab1, quantile(rowSums(tab1))[2])
plot(hclust(vegdist(rare_tab^0.25)),cex = 1)

#removing controls and including reruns samples insted of failed samples
rownames(tab1)[which(rownames(tab1) == "SE058b")] <- "SE058"
rownames(tab1)[which(rownames(tab1) == "SE060b")] <- "SE060"
rownames(tab1)[which(rownames(tab1) == "SE059b")] <- "SE059"
rownames(tab1)[which(rownames(tab1) == "SE082b")] <- "SE082"
rownames(tab1)[which(rownames(tab1) == "SE073o")] <- "SE073"
rownames(tab1)[which(rownames(tab1) == "SE063b")] <- "SE063"
rownames(tab1)[which(rownames(tab1) == "SE075b")] <- "SE075"
rownames(tab1)[which(rownames(tab1) == "SE079b")] <- "SE079"
rownames(tab1)[which(rownames(tab1) == "SE074b")] <- "SE074"
rownames(tab1)[which(rownames(tab1) == "SE081b")] <- "SE081"
rownames(tab1)[which(rownames(tab1) == "SE083b")] <- "SE083"

tab1 <- tab1[rownames(tab1) %in% c('SE001','SE002','SE003','SE004','SE005','SE006','SE007','SE008','SE009','SE010','SE011','SE012','SE013','SE014','SE015','SE016','SE017','SE018','SE019','SE020','SE021','SE022','SE023','SE024','SE025','SE026','SE027','SE028','SE029','SE030','SE031','SE032','SE033','SE034','SE035','SE036','SE037','SE038','SE039','SE040','SE041','SE042','SE043','SE044','SE045','SE046','SE047','SE048','SE049','SE050','SE051','SE052','SE053','SE054','SE055','SE056','SE057','SE058','SE059','SE060','SE061','SE062','SE063','SE064','SE065','SE066','SE067','SE068','SE069','SE070','SE071','SE072','SE073','SE074','SE075','SE076','SE077','SE078','SE079','SE080','SE081','SE082','SE083','SE084','SE088','SE089','SE090','SE091','SE092','SE093','SE094','SE095','SE096','SE097','SE098','SE099','SE100','SE101','SE102','SE103','SE104','SE105','SE106','SE107','SE108','SE109','SE110','SE111'),]

#checking for right labeling and completeness
org_set <- as.character(sampleinfo$sample[sampleinfo$batch %in% c("lts","original")])
tab_set <- rownames(tab1)
setdiff(org_set,tab_set) # S109 is missing from data (we knew that)
setdiff(tab_set, org_set)

tab1 <- tab1[order(as.numeric(substr(rownames(tab1),3,5))),]

rare_tab <- rrarefy(tab1, quantile(rowSums(tab1))[2])
plot(hclust(vegdist(rare_tab^0.25)),cex = .8)

saveRDS(tab1, here::here("asv_tables","fungi_combined_adjusted.rds"))

tab1 <- readRDS(here::here("asv_tables","fungi_combined_adjusted.rds"))

#trimming to ITS2 using custom script and assigning taxonomy using UNITE (via dada2)
tab_fmi <- itsx_on_dada2(tab1, asv_file_name = "asv_seqs_from_dada2.fasta", itsout_name = "fun_itsxout", org="T")
taxa <- assignTaxonomy(tab_fmi, here::here("ressources","sh_general_release_dynamic_s_10.05.2021.fasta"), multithread=TRUE) #needs to be downloaded from unite.ut.ee
taxa <- gsub("[kpcofgs]__","",taxa)
tab_fmi <- tab_fmi[,!is.na(taxa[,1])]
taxa <- taxa[!is.na(taxa[,1]),]
saveRDS(taxa,here::here("taxonomic_assignments","fun_taxa.rds"))
saveRDS(tab_fmi, here::here("asv_tables","fungi_combined_adjusted_its2_fungi.rds"))

#### Fungi Biodiversity calculations ----
sampleinfo <- read_csv(here::here("basic_data","sample_data.csv"))
tab2 <- readRDS(here::here("asv_tables","fungi_combined_adjusted_its2_fungi.rds"))
prep <- prepare_data(tab2, sampleinfo, use_partial = F)
saveRDS(prep, here::here("processed_data","fungi_prepared.rds"))

#### fungi taxonomy calculations ----
taxa <- readRDS(here::here("taxonomic_assignments","fun_taxa.rds"))
prep <- readRDS(here::here("processed_data","fungi_prepared.rds"))
fun_tax_profile <- taxonomic_turnover_wrapper(prepared_data=prep, taxonomy = taxa, markername = "Fungi", taxlevels = c("Kingdom","Phylum","Class","Order","Family","Genus"),num_tax = 10)
saveRDS(fun_tax_profile,here::here("processed_data","fungi_taxon_profiles.rds"))

prep <- readRDS(here::here("processed_data","fungi_prepared.rds"))
taxa <- readRDS(here::here("taxonomic_assignments","fun_taxa.rds"))
tax_turnover <- fun_tax_profile

#### fungi line graph ----
line_graphs <- overall_line_graphs(prep)
saveRDS(line_graphs, here::here("partial_plots","fungi_line_graph_single.rds"))

#### NMDS plot ----
mdsplot <- plotmds2(prep)
saveRDS(mdsplot, here::here("partial_plots","fungi_nmds.rds"))

#### fungi taxonomic plots
tax_turnover <- readRDS(here::here("processed_data","fungi_taxon_profiles.rds"))
tax_turnover[[1]]$class <- factor(tax_turnover[[1]]$class, levels = c("0°C/closed", "5°C/closed", "5°C/open",  "10°C/closed", "20°C/closed", "20°C/open", "40°C/closed"))
tax_turnover[[1]]$taxon <- substr(tax_turnover[[1]]$taxon,1,19)
tax_turnover[[2]]$class <- factor(tax_turnover[[2]]$class, levels = c("0°C/closed", "5°C/closed", "5°C/open",  "10°C/closed", "20°C/closed", "20°C/open", "40°C/closed"))
tax_turnover[[2]]$taxon <- substr(tax_turnover[[2]]$taxon,1,19)

#Design plots per taxon level and combine
tax_levels_used = c("Phylum","Class","Order","Family","Genus")
plotz <- taxonomic_turnover_plotter(tax_data = tax_turnover, tax_levels_used = tax_levels_used, legend_columns = 2)
saveRDS(plotz, here::here("partial_plots","fungi_raw_tax_plots.rds"))
tax_turn_plot <- plot_grid(plotlist=plotz, ncol=1, align = "v")
saveRDS(tax_turn_plot, here::here("partial_plots","fungi_combined_tax_plots.rds"))

#### supervised classification
#loading and combining tables from this study and the biowide reference study
tab2 <- readRDS(here::here("asv_tables","fungi_combined_adjusted.rds"))

# reading in the reference data and standardizing and correcting some sample names
# this reference data has pcr replicates that need to be combined.
bw <- readRDS(here::here("asv_tables_referencedata","bw_fungi_300bp_seqtab.nochim_Both.RDS"))
sampleinfo_bw <- read.table(here::here("asv_tables_referencedata","fun_bw_sample_info.txt"), header = T, stringsAsFactors = F)
names(sampleinfo_bw) <- c("sample","site", "replicate")
sampleinfo_bw <- sampleinfo_bw[!(grepl("Bla",sampleinfo_bw$site) | grepl("neg",sampleinfo_bw$site)),]
target <- tibble(sample = rownames(bw))
sampleinfo_bw2 <- sampleinfo_bw %>% inner_join(target, by = "sample")
bw <- bw[sampleinfo_bw2$sample,]
bw_col <- aggregate(x = bw, by = list(sampleinfo_bw2$site), FUN = "sum")
rownames(bw_col) <- bw_col[,1]
bw_col <- bw_col[order(as.numeric(substr(rownames(bw_col),3,5))),]
bw_col <- as.matrix(bw_col[,-1])
tab_m <- mergeSequenceTables(bw_col,tab2, repeats = "sum")

#merge with 
saveRDS(tab_m, here::here("asv_tables","fungi_bw_and_storage_raw.rds"))
tab_m <- readRDS(here::here("asv_tables","fungi_bw_and_storage_raw.rds"))

#restrict to ITS2 region
tab_mi <- itsx_on_dada2(tab_m) # itsx needs to be installed to run this wrapping script
saveRDS(tab_mi, here::here("asv_tables","fungi_bw_storage_itsx.rds"))

#extract sequences and map against UNITE to exclude non-fungi
uniquesToFasta(getUniques(tab_mi), fout=here::here("asv_tables","fungi_bw_otus.fasta"), ids=paste0("Seq", seq(length(getUniques(tab_mi)))))
otus <- here::here("taxonomic_assignments","fungi_bw_otus.fasta")
sh_file <- here::here("ressources", "sh_general_release_dynamic_10.05.2021.fasta") #needs to be downloaded from unite.ut.ee
commandX <- paste0('vsearch --usearch_global  "',otus,'" --dbmask none --qmask none --query_cov .9 --notrunclabels --userfields query --maxaccepts 1 --maxrejects 32 --maxhits 1 --db "', sh_file, '" --id 0.9 --iddef 0 --threads 8 --userout fungi_hits.hits')
system(commandX) # vsearch needs to be installed
fun_hits <- read_csv(here::here("taxonomic_assignments","fungi_hits.hits"), col_names = F)
otunumber <- sort(as.numeric(gsub("Seq","",fun_hits$X1)))
tab_mif <- tab_mi[,otunumber]
saveRDS(tab_mif, here::here("asv_tables","fungi_bw_storage_itsx_ingroup_only.rds"))

#read classifications and sampleinfo
tab_mif <- readRDS(here::here("asv_tables","fungi_bw_storage_itsx_ingroup_only.rds"))
sampleinfo <- read_csv(here::here("basic_data","sample_data.csv"))
sampleinfo <- sampleinfo[sampleinfo$sample != "SE109",]
sampleinfo$sample == rownames(tab_mif[-(1:130),]) # test

#resample and make dissimilarity matrix
tab_mif <- tab_mif[,colSums(tab_mif[1:130,])>10] # removing otus with low readcounts to facilitate faster analyses (no significant impact on results, has been tested)

tab3 <- rrarefy(tab_mif, quantile(rowSums(tab_mif))[2])
tab3_tr <- decostand(tab3, method = "hellinger")
vdm <- as.matrix(vegan::vegdist(tab3_tr, method = "bray"))

#### distance ratio plots ----
prep <- readRDS(here::here("processed_data","fungi_prepared.rds"))
ratio_plots <- plot_ratio_facet(dist_obj = vdm, prep = prep, sampleinfo = sampleinfo)
saveRDS(ratio_plots, here::here("partial_plots","fungi_ratio_plot.rds"))

#### KNN site matching
t0_references <- which((sampleinfo$harvest == 0 & sampleinfo$temperature == 5) | (sampleinfo$harvest == 0 &  sampleinfo$temperature == 0)) 
correct <- vector()
vdm_x <- vdm[131:237,c(1:80,82:130,130+t0_references)]
classifier <- c(rep(0,129), rep(1, length(t0_references)))
for(i in 1:nrow(sampleinfo)){
  correct[i] <- sum(classifier[order(vdm_x[i,], decreasing=F)[1:7]])/7
}
saveRDS(correct,here::here("processed_data","fungi_site_classification.rds"))

#### KNN Habitat classification ----
#KNN assigning of stored samples into pre-defined habitat classes
correct <- list()
knncount <- 9
classifier <- readRDS(here::here("processed_data","habitat_classifications.rds"))[c(1:80,82:130)]
k_clust <- length(table(classifier))

#classify the storage samples
for(i in 1:nrow(sampleinfo)){
  correct[[i]] <- as.list(table(classifier[order(vdm[130+i,c(1:80,82:130)], decreasing=F)[1:knncount]])/knncount)
}

correct_sn081 <- as.list(table(classifier[order(vdm[81,c(1:80,82:130)], decreasing=F)[1:knncount]])/knncount)
as.data.frame(correct_sn081)

#LateDryPoor LateDryRich LateWetRich
#  0.5555556   0.3333333   0.1111111

#put the data into a full table
biostratum_knn_result <- as.data.frame(rbindlist(correct, fill=TRUE))
biostratum_knn_result <- biostratum_knn_result[,order(as.numeric(names(biostratum_knn_result)))]
biostratum_knn_result[1:3,]

correct_cluster <- names(biostratum_knn_result)[which(biostratum_knn_result[1,] == max(biostratum_knn_result[1,], na.rm = T))]
#correct_cluster <- "LateDryRich"

# calculate distabce to data defined cluster centroid
dis_mat <- vdm[c(1:80,82:237), c(1:80,82:237)]
strata_seq <- c(classifier[1:129], rep("storage",107))
cen_dist <- dist_to_centroids(dis_mat, strata_seq)

distance_to_centroid <- cen_dist[grepl("SE",cen_dist$Item) & cen_dist$CentroidGroup == correct_cluster,]
distance_to_centroid_references <- cen_dist[cen_dist$Item %in% rownames(vdm)[which(classifier == correct_cluster)] & cen_dist$CentroidGroup == correct_cluster,]

result <- list(correct_cluster = correct_cluster, biostratum_knn = biostratum_knn_result, distance_to_centroid = distance_to_centroid, distance_to_centroid_references = distance_to_centroid_references)

saveRDS(result, here::here("processed_data","fungi_habitat_classification.rds"))

