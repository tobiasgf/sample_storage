#### preparation of basic datasets ####

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

#### Protist data preparation ----
# reading the dada2 otutables from the two MiSeqruns (run1: short storage, run2: long storage)
org <- readRDS(here::here("asv_tables","pro1_dada2_otutable.rds"))
lts <- readRDS(here::here("asv_tables","pro2_dada2_otutable.rds"))

sampleinfo <- read_csv(here::here("basic_data","sample_data.csv"))

#combine tables (storage only)
tab1 <- mergeSequenceTables(org,lts, repeats = "sum")
saveRDS(tab1, here::here("asv_tables","pro_combined_otutable.rds"))

tab1 <- readRDS(here::here("asv_tables","pro_combined_otutable.rds"))

#assign with DADA2. SKIPPED - takes forever
#taxa_pro <- assignTaxonomy(tab1, "~/tax/pr2_version_4.13.0_18S_dada2.fasta", multithread=TRUE)

#assign taxomnomy using pr2 database (download from here: https://github.com/pr2database/pr2database/releases)
uniquesToFasta(getUniques(tab1), fout=here::here("taxonomic_assignments","protist_otus.fasta"), ids=paste0("Seq", seq(length(getUniques(tab1)))))
# blastn -db ~/PATHTODIR/pr2_version_4.13.0_18S_dada2.fasta -max_target_seqs 10 -num_threads 50 -outfmt "6 std qlen qcovs sgi sseq ssciname staxid" -out pro_blasthits.txt -qcov_hsp_perc 80 -perc_identity 50 -query protist_otus.fasta

#quick check of clustering
rare_tab <- rrarefy(tab1, quantile(rowSums(tab1))[2])
plot(hclust(vegdist(rare_tab^0.25)),cex = 1)

tab1 <- tab1[rownames(tab1) %in% c('SE001','SE002','SE003','SE004','SE005','SE006','SE007','SE008','SE009','SE010','SE011','SE012','SE013','SE014','SE015','SE016','SE017','SE018','SE019','SE020','SE021','SE022','SE023','SE024','SE025','SE026','SE027','SE028','SE029','SE030','SE031','SE032','SE033','SE034','SE035','SE036','SE037','SE038','SE039','SE040','SE041','SE042','SE043','SE044','SE045','SE046','SE047','SE048','SE049','SE050','SE051','SE052','SE053','SE054','SE055','SE056','SE057','SE058','SE059','SE060','SE061','SE062','SE063','SE064','SE065','SE066','SE067','SE068','SE069','SE070','SE071','SE072','SE073','SE074','SE075','SE076','SE077','SE078','SE079','SE080','SE081','SE082','SE083','SE084','SE088','SE089','SE090','SE091','SE092','SE093','SE094','SE095','SE096','SE097','SE098','SE099','SE100','SE101','SE102','SE103','SE104','SE105','SE106','SE107','SE108','SE109','SE110','SE111'),]

#checking for right labeling and completeness
org_set <- as.character(sampleinfo$sample[sampleinfo$batch %in% c("lts","original")])
tab_set <- rownames(tab1)
setdiff(org_set,tab_set) # S109 is missing from data (we knew that)
setdiff(tab_set, org_set)
rare_tab <- rrarefy(tab1, quantile(rowSums(tab1))[2])
plot(hclust(vegdist(rare_tab^0.25)),cex = .8)

saveRDS(tab1, here::here("asv_tables","pro_combined_adjusted.rds"))

#taxonomy
pro_tax <- read.csv(here::here("taxonomic_assignments","pro_blasthits.txt"),sep='\t',header=F,as.is=TRUE)
pro_tax <- pro_tax[,1:3]
names(pro_tax) <- c("otu_id","SH_header","pident")

#### adjust annotation ####
#split taxon string into elements
pro_tax$kingdom <- str_split_fixed(pro_tax$SH_header, ";", 8)[,1]
pro_tax$supergroup <- str_split_fixed(pro_tax$SH_header, ";", 8)[,2]
pro_tax$division <- str_split_fixed(pro_tax$SH_header, ";", 8)[,3]
pro_tax$class <- str_split_fixed(pro_tax$SH_header, ";", 8)[,4]
pro_tax$order <- str_split_fixed(pro_tax$SH_header, ";", 8)[,5]
pro_tax$family <- str_split_fixed(pro_tax$SH_header, ";", 8)[,6]
pro_tax$genus <- str_split_fixed(pro_tax$SH_header, ";", 8)[,7]
#pro_tax$species <- str_split_fixed(pro_tax$SH_header, ";", 8)[,8]

pro_tax2 <- pro_tax %>% add_count(otu_id, name="total_n") %>% #calculate scores per taxonomic level
 # add_count(otu_id, species, name = "sp_n") %>% mutate(sp_ratio = sp_n/total_n) %>% 
  add_count(otu_id, genus, name = "gen_n") %>% mutate(gen_ratio = gen_n/total_n) %>%
  add_count(otu_id, family, name = "fam_n") %>% mutate(fam_ratio = fam_n/total_n) %>%
  add_count(otu_id, order, name = "ord_n") %>% mutate(ord_ratio = ord_n/total_n) %>%
  add_count(otu_id, class, name = "cla_n") %>% mutate(cla_ratio = cla_n/total_n) %>%
  add_count(otu_id, division, name = "div_n") %>% mutate(div_ratio = div_n/total_n) %>%
  add_count(otu_id, supergroup, name = "sup_n") %>% mutate(sup_ratio = sup_n/total_n) %>%
  add_count(otu_id, kingdom, name = "kin_n") %>% mutate(kin_ratio = kin_n/total_n) %>%
  group_by(otu_id) %>% slice_max(gen_ratio, n = 1, with_ties = FALSE) %>%  # select best scoring species per OTU
  #mutate(best_species_match=species) %>%
  #mutate(species=replace(species, sp_ratio < 0.5, "ambiguous")) %>%  # evaluate level of ambiguity per taxon (scores under 0.5 are ambuguous)
  mutate(genus=replace(genus, gen_ratio < 0.5, "ambiguous")) %>%
  mutate(family=replace(family, fam_ratio < 0.5, "ambiguous")) %>%
  mutate(order=replace(order, ord_ratio < 0.5, "ambiguous")) %>%
  mutate(class=replace(class, cla_ratio < 0.5, "ambiguous")) %>% 
  mutate(division=replace(division, div_ratio < 0.5, "ambiguous")) %>%
  mutate(supergroup=replace(supergroup, sup_ratio < 0.5, "ambiguous")) %>%
  mutate(kingdom=replace(kingdom, kin_ratio < 0.5, "ambiguous")) #%>%
  #mutate(species=replace(species, pident < 97, "imprecise")) # evaluate %match. Species level accepted above 97% 

otunumber <- as.numeric(gsub("Seq","",pro_tax2$otu_id))
table_dim <- dim(tab1)[2]
dada2_format_tax <- matrix("unidenfified", nrow = table_dim, ncol = 7)
dada2_format_tax[otunumber,1] <- pro_tax2$kingdom
dada2_format_tax[otunumber,2] <- pro_tax2$supergroup
dada2_format_tax[otunumber,3] <- pro_tax2$division
dada2_format_tax[otunumber,4] <- pro_tax2$class
dada2_format_tax[otunumber,5] <- pro_tax2$order
dada2_format_tax[otunumber,6] <- pro_tax2$family
dada2_format_tax[otunumber,7] <- pro_tax2$genus
#dada2_format_tax[otunumber,8] <- pro_tax2$species

colnames(dada2_format_tax) <- c("kingdom", "supergroup","division", "class", "order", "family", "genus")

saveRDS(dada2_format_tax,here::here("taxonomic_assignments","pro_taxa.rds"))

#### Protozoa Biodiversity calculations ----
sampleinfo <- read_csv(here::here("basic_data","sample_data.csv"))
tab2 <- readRDS(here::here("asv_tables","pro_combined_adjusted.rds"))
prep <- prepare_data(tab2, sampleinfo)
saveRDS(prep, here::here("processed_data","protist_prepared.rds"))

#### protozoa taxonomy calculations ----
taxa <- readRDS(here::here("taxonomic_assignments","pro_taxa.rds"))
prep <- readRDS(here::here("processed_data","protist_prepared.rds"))
pro_tax_profile <- taxonomic_turnover_wrapper(prepared_data=prep, taxonomy = taxa, markername = "Protozoa", taxlevels = c("kingdom","supergroup","division","class","order","family","genus"),num_tax = 11)

saveRDS(pro_tax_profile,here::here("data","protist_taxon_profiles.rds"))

prep <- readRDS(here::here("processed_data","protist_prepared.rds"))
taxa <- readRDS(here::here("taxonomic_assignments","pro_taxa.rds"))
tax_turnover <- readRDS(here::here("processed_data","protist_taxon_profiles.rds"))

#### protozoa line graph ----
line_graphs <- overall_line_graphs(prep)
line_graphs_compound <- plot_grid(plotlist=line_graphs, nrow=1)
saveRDS(line_graphs, here::here("partial_plots","protist_line_graph_single.rds"))
saveRDS(line_graphs_compound, here::here("partial_plots","protist_line_graph_compound.rds"))

#### NMDS plot ----
mdsplot <- plotmds2(prep)
saveRDS(mdsplot, here::here("partial_plots","protist_nmds.rds"))

#### taxonomic plots
tax_turnover <- readRDS(here::here("processed_data","protist_taxon_profiles.rds"))
tax_turnover[[1]]$class <- factor(tax_turnover[[1]]$class, levels = c("0°C/closed", "5°C/closed", "5°C/open",  "10°C/closed", "20°C/closed", "20°C/open", "40°C/closed"))
tax_turnover[[1]]$taxon <- substr(tax_turnover[[1]]$taxon,1,19)
tax_turnover[[2]]$class <- factor(tax_turnover[[2]]$class, levels = c("0°C/closed", "5°C/closed", "5°C/open",  "10°C/closed", "20°C/closed", "20°C/open", "40°C/closed"))
tax_turnover[[2]]$taxon <- substr(tax_turnover[[2]]$taxon,1,19)

#Design plots per taxon level and combine
tax_levels_used = c("supergroup","division","class","order","family","genus")
plotz <- taxonomic_turnover_plotter(tax_data = tax_turnover, tax_levels_used = tax_levels_used)
saveRDS(plotz, here::here("partial_plots","protist_raw_tax_plots.rds"))
tax_turn_plot <- plot_grid(plotlist=plotz, ncol=1, align = "v")
saveRDS(tax_turn_plot, here::here("partial_plots","protist_combined_tax_plots.rds"))

####  supervised classification ----
#loading and combining tables from this study and the biowide reference study
tab2 <- readRDS(here::here("asv_tables","pro_combined_adjusted.rds"))

# reading in the reference data and standardizing and correcting some sample names
bw <- readRDS(here::here("asv_tables_referencedata","bw_pro_seqtab.nochim_Both.rds"))
sampleinfo_bw <- read.table(here::here("asv_tables_referencedata","BW_PZ_sampleinfo.txt"), header = T, stringsAsFactors = F)
names(sampleinfo_bw) <- c("sample","site", "replicate")
sampleinfo_bw <- sampleinfo_bw[!(grepl("Bla",sampleinfo_bw$site) | grepl("Neg",sampleinfo_bw$site)),]
target <- tibble(sample = rownames(bw))
sampleinfo_bw2 <- sampleinfo_bw %>% inner_join(target, by = "sample")
bw <- bw[sampleinfo_bw2$sample,]
bw_col <- aggregate(x = bw, by = list(sampleinfo_bw2$site), FUN = "sum")
rownames(bw_col) <- bw_col[,1]
bw_col <- bw_col[order(as.numeric(substr(rownames(bw_col),3,5))),]
bw_col <- as.matrix(bw_col[,-1])

rownames(bw_col) <-  c("NV001", "NV002", "NV003", "NV004", "NV005", "NV006", "NV007", "NV008", "NV009", "NT010", "NT011", "NT012", "NT013", "NT014", "NT015", "NT016", "NT017", "NT018", "NH019", "NH020", "NH021", "NH022", "NH023", "NH024", "NH025", "NH026", "VU027", "VU028", "VU029", "VU030", "VU031", "VU032", "VU033", "VU034", "VU035", "VO036", "VO037", "VO038", "VO039", "VO040", "VO041", "VO042", "VO043", "VD044", "VD045", "VD046", "VD047", "VD048", "VD049", "VD050", "VD051", "VD052", "EM053", "EM054", "EM055", "EM056", "EM057", "EM058", "EM059", "EM060", "EM061", "ES062", "ES063", "ES064", "ES065", "ES066", "ES067", "ES068", "ES069", "ES070", "EV071", "EV072", "EV073", "EV074", "EV075", "EV076", "EV077", "EV078", "SN079", "SN080", "SN081", "SN082", "SN083", "SN084", "SN085", "SN086", "SN087", "SV088", "SV089", "SV090", "SV091", "SV092", "SV093", "SV094", "SV095", "SM096", "SM097", "SM098", "SM099", "SM100", "SM101", "SM102", "SM103", "SM104", "FF105", "FF106", "FF107", "FF108", "FF109", "FF110", "FF111", "FF112", "FL113", "FL114", "FL115", "FL116", "FL117", "FL118", "FL119", "FL120", "FM121", "FM122", "FM123", "FM124", "FM125", "FM126", "FM127", "FM128", "FM129", "FM130")

tab_m <- mergeSequenceTables(bw_col,tab2, repeats = "sum")
saveRDS(tab_m, here::here("sequence_data","pro_storage_and_biowide_otutable.rds"))

tab_m <- readRDS(here::here("asv_tables","pro_storage_and_biowide_otutable.rds"))
sampleinfo <- read.csv(here::here("basic_data","sampleinfo.txt"), sep = "\t")
sampleinfo <- sampleinfo[c(1:105,107:108),]
sampleinfo$sample == rownames(tab_m[-(1:130),]) # test

#resample and make dissimilarity matrix
tab_m <- tab_m[,colSums(tab_m[1:130,])>10]# removing otus with low readcounts to facilitate faster analyses (no significant impact on results, has been tested)

tab3 <- rrarefy(tab_m, quantile(rowSums(tab_m))[2])
tab3_tr <- decostand(tab3, method = "hellinger")
vdm <- as.matrix(vegan::vegdist(tab3_tr, method = "bray"))

#distance ratio plots
prep <- readRDS(here::here("processed_data","protist_prepared.rds"))
ratio_plots <- plot_ratio_facet(dist_obj = vdm, prep = prep, sampleinfo = sampleinfo)
saveRDS(ratio_plots, here::here("partial_plots","protist_ratio_plot.rds"))

#### KNN site matching
t0_references <- which((sampleinfo$harvest == 0 & sampleinfo$temperature == 5) | (sampleinfo$harvest == 0 &  sampleinfo$temperature == 0)) 
correct <- vector()
vdm_x <- vdm[131:237,c(1:80,82:130,130+t0_references)]
classifier <- c(rep(0,129), rep(1, length(t0_references)))
for(i in 1:nrow(sampleinfo)){
  correct[i] <- sum(classifier[order(vdm_x[i,], decreasing=F)[1:7]])/7
}
saveRDS(correct,here::here("processed_data","protist_site_classification.rds"))


#### KNN Biostratum
#KNN assigning of stored samples into pre-defined clusters
correct <- list()
knncount <- 9

classifier <- readRDS(here::here("processed_data","habitat_classifications.rds"))[c(1:80,82:130)] # supervised habitat classes of reference data
k_clust <- length(table(classifier))

#classify the storage samples
for(i in 1:nrow(sampleinfo)){
  correct[[i]] <- as.list(table(classifier[order(vdm[130+i,c(1:80,82:130)], decreasing=F)[1:knncount]])/knncount)
}

correct_sn081 <- as.list(table(classifier[order(vdm[81,c(1:80,82:130)], decreasing=F)[1:knncount]])/knncount)
as.data.frame(correct_sn081)

#LateDryPoor LateDryRich
#1   0.7777778   0.2222222

#put the data into a full table
biostratum_knn_result <- as.data.frame(rbindlist(correct, fill=TRUE))
biostratum_knn_result <- biostratum_knn_result[,order(as.numeric(names(biostratum_knn_result)))]
biostratum_knn_result[1:3,]

correct_cluster <- names(biostratum_knn_result)[which(biostratum_knn_result[1,] == max(biostratum_knn_result[1,], na.rm = T))]
#correct_cluster <- "LateDryRich"

# calculate distance to data defined cluster centroid
dis_mat <- vdm[c(1:80,82:237), c(1:80,82:237)]
strata_seq <- c(classifier[1:129], rep("storage",107))
cen_dist <- dist_to_centroids(dis_mat, strata_seq)

distance_to_centroid <- cen_dist[grepl("SE",cen_dist$Item) & cen_dist$CentroidGroup == correct_cluster,]
distance_to_centroid_references <- cen_dist[cen_dist$Item %in% rownames(vdm)[which(classifier == correct_cluster)] & cen_dist$CentroidGroup == correct_cluster,]

result <- list(correct_cluster = correct_cluster, biostratum_knn = biostratum_knn_result, distance_to_centroid = distance_to_centroid, distance_to_centroid_references = distance_to_centroid_references)

saveRDS(result, here::here("processed_data","protists_habitat_classification.rds"))
