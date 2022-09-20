# This script produces the supervised classification of the reference data.

library(vegan)
library(readxl)
library(tidyverse)
library(vegan)

# Here we read in the full species observation data from the project Biowide with the 130 sites that are used as a reference dataset for the study here.

bw_data <- as.data.frame(read_xlsx(here::here("ressources","SamletArtsdata.xlsx")))
bw_data$sitex <- paste0(bw_data$cluster, str_pad(bw_data$site_nr, 3, pad = "0"))
bw_data <- bw_data[,c("latin","sitex")]
bw_data$present <- 1
bw_wide <- bw_data %>% spread(sitex, present, fill = 0)
bw_wide2 <- bw_wide[,-1]
bw_wide3 <- as.data.frame(t(as.matrix(bw_wide2)))
bw_wide4 <- bw_wide3[order(as.numeric(substr(rownames(bw_wide3),3,5))),]

#We do an nmds ordination in 6 dimensions
biota_tobias.nms <- metaMDS(bw_wide4, k=6, try=100, trymax = 200)

#read in the expert based classifications of the 130 sites, and use this to make a supervised classification of the dataset.
bibi <- read.delim(here::here("ressources","expert_classification.txt"))
bibi2 <- cbind(bibi,cbind(scores(biota_tobias.nms)$sites[,1:4]))
library(MASS)
biota.qda <- qda(expert_classification~NMDS1+NMDS2+NMDS3+NMDS4, data=bibi2)

#save the predicted habitat classes as a vector (this is the habitat types used for the habitat predictions in the supervised classification of the stored samples)

saveRDS(as.vector(predict(biota.qda)$class), here::here("processed_data","habitat_classifications.rds"))
