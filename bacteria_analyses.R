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

#reading in basic metadata
sampleinfo <- read.csv(here::here("basic_data","sampleinfo.txt"), sep = "\t", stringsAsFactors = F)
sampleinfo <- sampleinfo[1:108,]
write_csv(sampleinfo, here::here("basic_data","sample_data.csv"))

#### Bacteria data preparation ----
# reading the dada2 otutables from the three MiSeqruns (run1: short storage, run2: long storage, run3: reruns of failed samples)
org <- readRDS(here::here("asv_tables","bac1_dada2_otutable.rds"))
lts <- readRDS(here::here("asv_tables","bac2_dada2_otutable.rds"))
rerun <- readRDS(here::here("asv_tables","bac3_dada2_otutable.rds"))
# combine tables and save
tab1 <- mergeSequenceTables(org,lts,rerun, repeats = "sum")
saveRDS(tab1, here::here("asv_tables","bac_combined_otutable.rds"))

# assign taxomnomy and save
tab1 <- readRDS(here::here("asv_tables","bac_combined_otutable.rds"))
taxa <- assignTaxonomy(tab1, here::here("ressources","silva_nr_v132_train_set.fa.gz"), multithread=TRUE) #CHECK
saveRDS(taxa,here::here("taxonomic_assignments","bac_taxa.rds"))

# quick check of clustering
plot(hclust(vegdist(tab1)),cex = .5)

#focussing on study samples and adding one sample from rerun (missing from first run)
tab2 <- tab1[grepl("SE",rownames(tab1)) & !grepl("SE008b",rownames(tab1)),]
rownames(tab2)[which(rownames(tab2) == "SE008c")] <- "SE008"
plot(hclust(vegdist(tab2)),cex = .8) # looks ok

#reorder samples according to sequence in metadata
tab2 <- tab2[order(as.numeric(substr(rownames(tab2),3,5))),]
saveRDS(tab2, here::here("asv_tables","bac_combined_adjusted.rds"))

#### Bacteria Biodiversity calculations ----
sample_data <- read_csv(here::here("basic_data","sample_data.csv"))
tab2 <- readRDS(here::here("asv_tables","bac_combined_adjusted.rds"))
prep <- prepare_data(tab2, sample_data, use_partial = F)
saveRDS(prep, here::here("processed_data","bacteria_prepared.rds"))

#### Bacteria taxonomy calculations ----
taxa <- readRDS(here::here("taxonomic_assignments","bac_taxa.rds"))
prep <- readRDS(here::here("processed_data","bacteria_prepared.rds"))
bac_tax_profile <- taxonomic_turnover_wrapper(prepared_data=prep, taxonomy = taxa, markername = "Bacteria", taxlevels = c("Kingdom","Phylum","Class","Order","Family","Genus"),num_tax = 11)
saveRDS(bac_tax_profile,here::here("processed_data","bacteria_taxon_profiles.rds"))

prep <- readRDS(here::here("processed_data","bacteria_prepared.rds"))
taxa <- readRDS(here::here("taxonomic_assignments","bac_taxa.rds"))
tax_turnover <- readRDS(here::here("processed_data","bacteria_taxon_profiles.rds"))

#### Bacteria line graph ----
line_graphs <- overall_line_graphs(prep)
#line_graphs_compound <- plot_grid(plotlist=line_graphs, nrow=1)
saveRDS(line_graphs, here::here("partial_plots","bacteria_line_graph_single.rds"))

#### NMDS plot ----
mdsplot <- plotmds2(prep)
saveRDS(mdsplot, here::here("partial_plots","bacteria_nmds.rds"))

#### taxonomic plots
tax_turnover <- readRDS(here::here("processed_data","bacteria_taxon_profiles.rds"))
tax_turnover[[1]]$class <- factor(tax_turnover[[1]]$class, levels = c("0°C/closed", "5°C/closed", "5°C/open",  "10°C/closed", "20°C/closed", "20°C/open", "40°C/closed"))
tax_turnover[[1]]$taxon <- substr(tax_turnover[[1]]$taxon,1,19)
tax_turnover[[2]]$class <- factor(tax_turnover[[2]]$class, levels = c("0°C/closed", "5°C/closed", "5°C/open",  "10°C/closed", "20°C/closed", "20°C/open", "40°C/closed"))
tax_turnover[[2]]$taxon <- substr(tax_turnover[[2]]$taxon,1,19)

#Design plots per taxon level and combine
tax_levels_used = c("Kingdom","Phylum","Class","Order","Family","Genus")
plotz <- taxonomic_turnover_plotter(tax_data = tax_turnover, tax_levels_used = tax_levels_used, legend_columns = 2)
saveRDS(plotz, here::here("partial_plots","bacteria_raw_tax_plots.rds"))
tax_turn_plot <- plot_grid(plotlist=plotz, ncol=1, align = "v")
saveRDS(tax_turn_plot, here::here("partial_plots","bacteria_combined_tax_plots.rds"))


####  supervised classification ----
#loading and combining tables from this study and the biowide reference study
tab2 <- readRDS(here::here("asv_tables","bac_combined_adjusted.rds"))

# reading in the reference data and standardizing and correcting some sample names
bw_col <- readRDS(here::here("asv_tables_referencedata","bw_bac_seqtab.nochim_Both_newrun.RDS"))
rownames(bw_col) <- gsub("URVC_","",rownames(bw_col))
bw_col <- bw_col[order(as.numeric(substr(rownames(bw_col),3,5))),]
rownames(bw_col)[c(27,28,65,66,67,68,69,70,71,72,73,74,75,76,77,78,79,80,81,82,83,84)] <- c("VU027", "VU028", "ES065", "ES066", "ES067", "ES068", "ES069", "ES070", "EV071", "EV072", "EV073", "EV074", "EV075", "EV076", "EV077", "EV078", "SN079", "SN080", "SN081", "SN082", "SN083", "SN084") 

tab_m <- mergeSequenceTables(bw_col,tab2, repeats = "sum")
saveRDS(tab_m, here::here("asv_tables","bac_storage_and_biowide_otutable.rds"))

tab_m <- readRDS(here::here("asv_tables","bac_storage_and_biowide_otutable.rds"))
sampleinfo <- read.csv(here::here("basic_data","sampleinfo.txt"), sep = "\t")
sampleinfo <- sampleinfo[c(1:105,107:108),]

sampleinfo$sample == rownames(tab_m[-(1:130),]) # test if names are aligned

tab_m <- tab_m[,colSums(tab_m[1:130,])>10] # removing otus with low readcounts to facilitate faster analyses (no significant impact on results, has been tested)

#resample (rarefy) and make dissimilarity matrix
tab3 <- rrarefy(tab_m, quantile(rowSums(tab_m))[2])
tab3_tr <- decostand(tab3, method = "hellinger")
vdm <- as.matrix(vegan::vegdist(tab3_tr, method = "bray"))

#distance ratio plots
prep <- readRDS(here::here("processed_data","bacteria_prepared.rds"))
ratio_plots <- plot_ratio_facet(dist_obj = vdm, prep = prep, sampleinfo = sampleinfo)
saveRDS(ratio_plots, here::here("partial_plots","bacteria_ratio_plot.rds"))

#### KNN site matching
t0_references <- which((sampleinfo$harvest == 0 & sampleinfo$temperature == 5) | (sampleinfo$harvest == 0 &  sampleinfo$temperature == 0)) 
correct <- vector()
vdm_x <- vdm[131:237,c(1:80,82:130,130+t0_references)]
classifier <- c(rep(0,129), rep(1, length(t0_references)))
for(i in 1:nrow(sampleinfo)){
  correct[i] <- sum(classifier[order(vdm_x[i,], decreasing=F)[1:7]])/7
}
saveRDS(correct,here::here("processed_data","bacteria_site_classification.rds"))


#### KNN Habitat classification ----
#KNN assigning of stored samples into pre-defined habitat classes
correct <- list()
knncount <- 9

classifier <- readRDS(here::here("processed_data","habitat_classifications.rds"))[c(1:80,82:130)] # supervised habitat classes of reference data
k_clust <- length(table(classifier)) # number of classes

#classify the storage samples
for(i in 1:nrow(sampleinfo)){
  correct[[i]] <- as.list(table(classifier[order(vdm[130+i,c(1:80,82:130)], decreasing=F)[1:knncount]])/knncount)
}

correct_sn081 <- as.list(table(classifier[order(vdm[81,c(1:80,82:130)], decreasing=F)[1:knncount]])/knncount)
as.data.frame(correct_sn081)

#EarlyDryPoor LateDryPoor LateDryRich
#1    0.1111111   0.7777778   0.1111111

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

saveRDS(result, here::here("processed_data","bacteria_habitat_classification.rds"))
