library(flextable)
library(officer)
library(tidyverse)
library(cowplot)
library(ggpubr)
library(lemon)
source(here::here("R","scripts.R"))
library(readxl)
library(plyr)
library(here)

#### co2 plot ----
sampleinfo <- read.csv(here::here("basic_data","sampleinfo.txt"), sep = "\t", stringsAsFactors = F)

co2data2 <- readxl::read_xlsx(here::here("basic_data","co2_data_100621.xlsx"), col_names = T, range = "R2C1:R62C11")

df1 <- data_summary(co2data2, varname="ppm4", 
                    groupnames=c("temperature", "harvest","exposure"))

df1$fTemperature <- factor(df1$temperature)
df1$fHarvest <- factor(df1$harvest)
co_df <- df1

p1 <- ggplot(df1, aes(x=fHarvest, y=ppm4, group=fTemperature, color=fTemperature)) + 
 geom_errorbar(aes(ymin=ppm4-sem, ymax=ppm4+sem), width=.1, size = 0.3, 
               position=position_dodge(0.05)) + scale_shape_manual(values=c(17,1)) +
 geom_line(position=position_dodge(0.05), size = 0.3) + geom_point(position=position_dodge(0.05))+
 scale_color_brewer(type = 'div', palette = 'Spectral', direction = -1, "Temperature") + theme_minimal() + theme_bw() + 
 theme(panel.grid.minor = element_blank(),panel.grid.major = element_blank(),text=element_text(family="Helvetica", face="plain", size=10),
       axis.text = element_text(colour="black"),
       panel.border = element_blank(),
       axis.line = element_line(colour = 'black', size = 0.25),
       strip.background = element_rect(colour="white", fill="white")) + xlab("Storage time (days)") + ylab(bquote('CO'[2]~'(ppm/hour)'))

saveRDS(p1, here::here("partial_plots","co2_plot.rds"))

#### DNA plot ----
sampleinfo <- read.csv(here::here("basic_data","sampleinfo.txt"), sep = "\t", stringsAsFactors = F)
samples1 <- readxl::read_xlsx(here::here("basic_data","Qubit and dilution of SE Soil Samples.xlsx"), range = "A3:B87", col_names = T)
samples2 <- readxl::read_xlsx(here::here("basic_data","Qubit and dilution of SE Soil Samples.xlsx"), range = "A94:B117", col_names = F)
names(samples2) <- names(samples1)
dna <- rbind(samples1,samples2)
sample_dna <- cbind(sampleinfo[1:108,], dna)

sample_dna[106,"Qubit"] <- NA
sample_dna$Qubit <- as.numeric(sample_dna$Qubit)
df1 <- data_summary(sample_dna,varname="Qubit", 
                    groupnames=c("temperature", "harvest", "exposure"))

df1$fTemperature <- factor(df1$temperature)
df1$fHarvest <- factor(df1$harvest)

p1 <- ggplot(df1, aes(x=fHarvest, y=Qubit, group=interaction(fTemperature,exposure), color=fTemperature, shape = exposure)) + 
 geom_errorbar(aes(ymin=Qubit-sem, ymax=Qubit+sem), width=.1, size = 0.3, 
               position=position_dodge(0.05)) + scale_shape_manual(values=c(17,1)) +
 geom_line(position=position_dodge(0.05), size = 0.3) + geom_point(position=position_dodge(0.05))+
 scale_color_brewer(type = 'div', palette = 'Spectral', direction = -1, "Temperature") + theme_minimal() + theme_bw() + labs(shape = "Exposure") +
 theme(panel.grid.minor = element_blank(),panel.grid.major = element_blank(),text=element_text(family="Helvetica", face="plain", size=10),
       axis.text = element_text(colour="black"),
       panel.border = element_blank(),
       axis.line = element_line(colour = 'black', size = 0.25),
       strip.background = element_rect(colour="white", fill="white"))+ ylab("DNA conc (ng/µl)") + xlab("Storage time (days)")

saveRDS(p1, here::here("partial_plots","dna_concentration.rds"))

#### Supplementary table 4 - co2 and dna concentration ----
co_df <- co_df %>% select(temperature, harvest, exposure, ppm4, sem) %>% mutate(ppm4 = round(ppm4, digits = 1), sem = round(sem, digits = 1))
dna_df <- df1 %>%select(temperature, harvest, exposure, Qubit, sem) %>% mutate(Qubit = round(Qubit, digits = 1), sem = round(sem, digits = 1))
tot_tab <- left_join(dna_df,co_df, by = c("temperature" = "temperature", "harvest" = "harvest", "exposure" = "exposure")) %>% arrange(temperature, exposure, harvest)

write.table(tot_tab, here::here("tables","dna_co2.txt"), sep = "\t", row.names = F, col.names = T, quote = F)


#### Figure 2 line graphs ----
p1b <- readRDS(here::here("partial_plots","bacteria_line_graph_single.rds"))
p1f <- readRDS(here::here("partial_plots","fungi_line_graph_single.rds"))
p1p <- readRDS(here::here("partial_plots","protist_line_graph_single.rds"))
p1e <- readRDS(here::here("partial_plots","eukaryote_line_graph_single.rds"))

dna <- readRDS(here::here("partial_plots","dna_concentration.rds"))
co2 <- readRDS(here::here("partial_plots","co2_plot.rds"))
dna <- dna + theme(legend.position = "none")
co2 <- co2 + theme(legend.position = "none")
bottom_row <- plot_grid(NULL, dna, co2, NULL, ncol=4, rel_widths=c(1,2,2,1), labels = c("","e","f",""),label_size = 10)

line_list <- list(p1b[[1]], p1b[[3]], p1b[[4]], p1f[[1]], p1f[[3]], p1f[[4]], p1p[[1]],p1p[[3]],p1p[[4]],p1e[[1]],p1e[[3]],p1e[[4]])
line_compound <- plot_grid(plotlist=line_list, rel_widths=c(2,2,2),ncol=3, labels = c("a","","","b","","","c","","","d","",""), label_size = 10)

compound_compound <- plot_grid(plotlist=list(line_compound, bottom_row, p1f[[5]]), ncol=1, rel_heights = c(20,5,1))

ggsave(here::here("final_plots/a_from_r","Fig2.pdf"), compound_compound, width = 25, height = 30, units = "cm")

#### Figure 3 NMDS ----
p2b <- readRDS(here::here("partial_plots","bacteria_nmds.rds"))
p2f <- readRDS(here::here("partial_plots","fungi_nmds.rds"))
p2p <- readRDS(here::here("partial_plots","protist_nmds.rds"))
p2e <- readRDS(here::here("partial_plots","eukaryote_nmds.rds"))

legend <- g_legend(p2b)
p2b <- p2b + theme(legend.position = "none")
p2f <- p2f + theme(legend.position = "none")
p2p <- p2p + theme(legend.position = "none")
p2e <- p2e + theme(legend.position = "none")

nmds_list <- list(p2b, p2f, p2p, p2e, legend) 

compound_compound <- plot_grid(plotlist=nmds_list, ncol=1, rel_heights = c(5,5,5,5,1), labels = c("a","b","c","d",""))
ggsave(here::here("final_plots/a_from_r","Fig3.pdf"), compound_compound, width = 25, height = 26, units = "cm")

#### Figure 4 combined taxonomic class_plot ----
tax_p_b <- readRDS(here::here("partial_plots","bacteria_raw_tax_plots.rds"))
tax_p_f <- readRDS(here::here("partial_plots","fungi_raw_tax_plots.rds"))
tax_p_p <- readRDS(here::here("partial_plots","protist_raw_tax_plots.rds"))
tax_p_e <- readRDS(here::here("partial_plots","eukaryote_raw_tax_plots.rds"))

compound_compound <- plot_grid(plotlist=list(tax_p_b[[3]], tax_p_f[[2]], tax_p_p[[3]], tax_p_e[[3]]), ncol=1, labels = c("a","b","c","d"), align = "v") #+ theme(plot.margin = margin(0.1,0.1,0.5,0.1, "cm"))

ggsave(here::here("final_plots/a_from_r","Fig4.pdf"), compound_compound, width = 22, height = 22, units = "cm")


#### Figure 5 - site classification KNN ----
sampleinfo <- read.csv(here::here("basic_data","sampleinfo.txt"), sep = "\t", stringsAsFactors = F)
sampleinfo <- sampleinfo[c(1:105,107:108),]

site_b <- readRDS(here::here("processed_data","bacteria_site_classification.rds"))
bac_plot <- plot_knn_site(site_b, sampleinfo)
site_f <- readRDS(here::here("processed_data","fungi_site_classification.rds"))
fun_plot <- plot_knn_site(site_f, sampleinfo)
site_p <- readRDS(here::here("processed_data","protist_site_classification.rds"))
pro_plot <- plot_knn_site(site_p, sampleinfo)
site_e <- readRDS(here::here("processed_data","eukaryote_site_classification.rds"))
tar_plot <- plot_knn_site(site_e, sampleinfo)

legend <- g_legend(bac_plot)
bac_plot <- bac_plot +  theme(legend.position = "none")
tar_plot <- tar_plot +  theme(legend.position = "none")
fun_plot <- fun_plot +  theme(legend.position = "none")
pro_plot <- pro_plot +  theme(legend.position = "none")

line_list <- list(bac_plot, fun_plot, pro_plot, tar_plot)
site_knn_compound <- plot_grid(plotlist=line_list, rel_widths=c(2,2,2,2),nrow=1, labels = c("a","b","c","d"), label_size = 10)

site_knn_compound_legend <- plot_grid(plotlist=list(site_knn_compound, legend),ncol=1, rel_heights = c(5,1))

ggsave(here::here("final_plots/a_from_r","Fig5.pdf"), site_knn_compound_legend, width = 26, height = 8, units = "cm")


#### Figure 6 ratio_plots  ----
rp_b <- readRDS(here::here("partial_plots","bacteria_ratio_plot.rds"))
rp_f <- readRDS(here::here("partial_plots","fungi_ratio_plot.rds"))
rp_p <- readRDS(here::here("partial_plots","protist_ratio_plot.rds"))
rp_e <- readRDS(here::here("partial_plots","eukaryote_ratio_plot.rds"))

legendx <- g_legend(rp_p[[1]])
rp_b <- rp_b[[1]] + theme(legend.position = "none")
rp_f <- rp_f[[1]] + theme(legend.position = "none")
rp_p <- rp_p[[1]] + theme(legend.position = "none")
rp_e <- rp_e[[1]] + theme(legend.position = "none")

rp_list <- list(rp_b, rp_f, rp_p, rp_e)
compound_compound <- plot_grid(plotlist=rp_list, ncol=1, labels = c("a","b","c","d"))+ theme(plot.margin = margin(0.1,0.1,0.5,0.1, "cm"))
compound_compound2 <- ggarrange(compound_compound, legendx, ncol=1, heights = c(30,1)) + theme(plot.margin = margin(0.1,0.1,0.5,0.1, "cm"))
ggsave(here::here("final_plots/a_from_r","Fig6.pdf"), compound_compound2, width = 24, height = 28, units = "cm")

#### Figure 7 - KNN Habitat classification ----
sampleinfo <- read.csv(here::here("basic_data","sampleinfo.txt"), sep = "\t", stringsAsFactors = F)
sampleinfo <- sampleinfo[c(1:105,107:108),]
knn_data <- readRDS(here::here("processed_data","bacteria_habitat_classification.rds"))
bac_plot <- plot_new_models(sampleinfo, knn_data)
knn_data <- readRDS(here::here("processed_data","fungi_habitat_classification.rds"))
fun_plot <- plot_new_models(sampleinfo, knn_data)
knn_data <- readRDS(here::here("processed_data","protist_habitat_classification.rds"))
pro_plot <- plot_new_models(sampleinfo, knn_data)
knn_data <- readRDS(here::here("processed_data","eukaryote_habitat_classification.rds"))
tar_plot <- plot_new_models(sampleinfo, knn_data)

legend <- g_legend(bac_plot[[2]] + theme(legend.position = "bottom")) #+ guides(fill=guide_legend(nrow=1, reverse = T)))
bac_plot[[2]] <- bac_plot[[2]] +  theme(legend.position = "none")
tar_plot[[2]] <- tar_plot[[2]] +  theme(legend.position = "none")
fun_plot[[2]] <- fun_plot[[2]] +  theme(legend.position = "none")
pro_plot[[2]] <- pro_plot[[2]] +  theme(legend.position = "none")

line_list <- list(bac_plot[[2]], fun_plot[[2]], pro_plot[[2]], tar_plot[[2]])
hab_cluster_compound <- plot_grid(plotlist=line_list, rel_widths=c(2,2,2,2),nrow=1, labels = c("a","b","c","d"), label_size = 10) #+ theme(plot.margin = margin(0.1,0.1,0.5,0.1, "cm"))

hab_cluster_legend <- plot_grid(plotlist=list(hab_cluster_compound, legend),ncol=1, rel_heights = c(5,1)) #+ theme(plot.margin = margin(0.1,0.1,0.5,0.1, "cm"))

# Dist to habitat cluster centroid ---
sampleinfo <- read.csv(here::here("basic_data","sampleinfo.txt"), sep = "\t", stringsAsFactors = F)
sampleinfo <- sampleinfo[c(1:105,107:108),]
knn_data <- readRDS(here::here("processed_data","bacteria_habitat_classification.rds"))
bac_plot <- plot_new_models(sampleinfo, knn_data)
knn_data <- readRDS(here::here("processed_data","fungi_habitat_classification.rds"))
fun_plot <- plot_new_models(sampleinfo, knn_data)
knn_data <- readRDS(here::here("processed_data","protist_habitat_classification.rds"))
pro_plot <- plot_new_models(sampleinfo, knn_data)
knn_data <- readRDS(here::here("processed_data","eukaryote_habitat_classification.rds"))
tar_plot <- plot_new_models(sampleinfo, knn_data)

legend <- g_legend(bac_plot[[1]])
bac_plot[[1]] <- bac_plot[[1]] +  theme(legend.position = "none")
tar_plot[[1]] <- tar_plot[[1]] +  theme(legend.position = "none")
fun_plot[[1]] <- fun_plot[[1]] +  theme(legend.position = "none")
pro_plot[[1]] <- pro_plot[[1]] +  theme(legend.position = "none")

line_list <- list(bac_plot[[1]], fun_plot[[1]], pro_plot[[1]], tar_plot[[1]])
hab_change_compound <- plot_grid(plotlist=line_list, rel_widths=c(2,2,2,2),nrow=1, labels = c("e","f","g","h"), label_size = 10)

hab_change_compound_legend <- plot_grid(plotlist=list(hab_change_compound, legend),ncol=1, rel_heights = c(5,1))

#Combine the three plots
knn_plot <- plot_grid(hab_cluster_legend, hab_change_compound_legend, ncol=1, rel_heights = c(4,6))

ggsave(here::here("final_plots/a_from_r","Fig7.pdf"), knn_plot, width = 26, height = 20, units = "cm")


#### Supplementary Figures S1-S4 taxonomic plots ----
tax_b <- readRDS(here::here("partial_plots","bacteria_combined_tax_plots.rds"))
ggsave(here::here("final_plots/a_from_r","Fig_S1.pdf"), tax_b, width = 25, height = 32, units = "cm")

tax_f <- readRDS(here::here("partial_plots","fungi_combined_tax_plots.rds"))
ggsave(here::here("final_plots/a_from_r","Fig_S2.pdf"), tax_f, width = 25, height = 32, units = "cm")

tax_p <- readRDS(here::here("partial_plots","protozoa_combined_tax_plots.rds"))
ggsave(here::here("final_plots/a_from_r","Fig_S3.pdf"), tax_p, width = 25, height = 32, units = "cm")

tax_e <- readRDS(here::here("partial_plots","tar_combined_tax_plots.rds"))
ggsave(here::here("final_plots/a_from_r","Fig_S4.pdf"), tax_e, width = 25, height = 32, units = "cm")


#### Supplementary Figure S5 - dist other plots same plot ----
rp_b <- readRDS(here::here("partial_plots","bacteria_ratio_plot.rds"))
rp_f <- readRDS(here::here("partial_plots","fungi_ratio_plot.rds"))
rp_p <- readRDS(here::here("partial_plots","protist_ratio_plot.rds"))
rp_e <- readRDS(here::here("partial_plots","eukaryote_ratio_plot.rds"))

legendx <- g_legend(rp_p[[1]])
rp_b <- rp_b[[3]] + theme(legend.position = "none")
rp_f <- rp_f[[3]] + theme(legend.position = "none")
rp_p <- rp_p[[3]] + theme(legend.position = "none")
rp_e <- rp_e[[3]] + theme(legend.position = "none")

rp_list <- list(rp_b, rp_f, rp_p, rp_e)
compound_compound <- plot_grid(plotlist=rp_list, ncol=1, labels = c("a","b","c","d"))+ theme(plot.margin = margin(0.1,0.1,0.5,0.1, "cm"))
compound_compound2 <- ggarrange(compound_compound, legendx, ncol=1, heights = c(20,1)) + theme(plot.margin = margin(0.1,0.1,0.5,0.1, "cm"))

ggsave(here::here("final_plots/a_from_r","Fig_S5.pdf"), compound_compound2, width = 20, height = 30, units = "cm")


#### Supplementary figure 6 - KNN Habitat classification second best class ----
sampleinfo <- read.csv(here::here("basic_data","sampleinfo.txt"), sep = "\t", stringsAsFactors = F)
sampleinfo <- sampleinfo[c(1:105,107:108),]
knn_data <- readRDS(here::here("processed_data","bacteria_habitat_classification.rds"))
bac_plot <- plot_new_models(sampleinfo, knn_data)
knn_data <- readRDS(here::here("processed_data","fungi_habitat_classification.rds"))
fun_plot <- plot_new_models(sampleinfo, knn_data)
knn_data <- readRDS(here::here("processed_data","protist_habitat_classification.rds"))
pro_plot <- plot_new_models(sampleinfo, knn_data)
knn_data <- readRDS(here::here("processed_data","eukaryote_habitat_classification.rds"))
tar_plot <- plot_new_models(sampleinfo, knn_data)

legend <- g_legend(bac_plot[[3]] + theme(legend.position = "bottom")) #+ guides(fill=guide_legend(nrow=1, reverse = T)))
bac_plot[[3]] <- bac_plot[[3]] +  theme(legend.position = "none")
tar_plot[[3]] <- tar_plot[[3]] +  theme(legend.position = "none")
fun_plot[[3]] <- fun_plot[[3]] +  theme(legend.position = "none")
pro_plot[[3]] <- pro_plot[[3]] +  theme(legend.position = "none")

line_list <- list(bac_plot[[3]], fun_plot[[3]], pro_plot[[3]], tar_plot[[3]])
hab_cluster_compound <- plot_grid(plotlist=line_list, rel_widths=c(2,2,2,2),nrow=1, labels = c("a","b","c","d"), label_size = 10) #+ theme(plot.margin = margin(0.1,0.1,0.5,0.1, "cm"))

hab_cluster_legend <- plot_grid(plotlist=list(hab_cluster_compound, legend),ncol=1, rel_heights = c(5,1)) #+ theme(plot.margin = margin(0.1,0.1,0.5,0.1, "cm"))

ggsave(here::here("final_plots/a_from_r","Fig_S6.pdf"), hab_cluster_compound, width = 26, height = 7, units = "cm")





#### Supplementary Figure 7 - rarefaction and richness correlation plots ----
bact <- readRDS(here::here("asv_tables","bac_combined_adjusted.rds"))
funt <- readRDS(here::here("asv_tables","fungi_combined_adjusted_its2_fungi.rds"))
tart <- readRDS(here::here("asv_tables","euk_combined_adjusted.rds"))
prot <- readRDS(here::here("asv_tables","pro_combined_adjusted.rds"))

bactx <- rare_rich_calc(bact, dataname = "Bacteria")
funtx <- rare_rich_calc(funt, dataname = "Fungi")
tartx <- rare_rich_calc(tart, dataname = "Eukaryotes")
protx <- rare_rich_calc(prot, dataname = "Protists")

bactx1 <- plot_rrare(bactx, dataname = "Bacteria")
bactx2 <- plot_rare_rich(bactx, dataname = "Bacteria")

funtx1 <- plot_rrare(funtx, dataname = "Fungi")
funtx2 <- plot_rare_rich(funtx, dataname = "Fungi")

protx1 <- plot_rrare(protx, dataname = "Protists")
protx2 <- plot_rare_rich(protx, dataname = "Protists")

tartx1 <- plot_rrare(tartx, dataname = "Eukaryotes")
tartx2 <- plot_rare_rich(tartx, dataname = "Eukaryotes")

ddd <- plot_rrare(bactx$richness, dataname = "Bacteria")
ddd2 <- plot_rare_rich(bactx, dataname = "Bacteria")

plotlistx <- list(bactx1, bactx2, funtx1, funtx2, protx1, protx2, tartx1, tartx2)
rareplot <- plot_grid(plotlist = plotlistx, ncol = 2)
ggsave(here::here("final_plots/a_from_r","Fig_S7.pdf"), rareplot, width = 15, height = 20, units = "cm")


#### check stress of nmds plots ----
b_prep <- readRDS(here::here("processed_data","bacteria_prepared.rds"))
b_prep$nmds$stress # [1] 0.07462213

f_prep <- readRDS(here::here("processed_data","fungi_prepared.rds"))
f_prep$nmds$stress # [1] 0.08163422

p_prep <- readRDS(here::here("processed_data","protist_prepared.rds"))
p_prep$nmds$stress # [1] [1] 0.119655

e_prep <- readRDS(here::here("data","eukaryote_prepared.rds"))
e_prep$nmds$stress # [1] 0.110721

#### Supplementary table 2 - changes in biodiversity measures  ----
# table of statistics
t1b <- readRDS(here::here("processed_data","bacteria_prepared.rds"))
t1f <- readRDS(here::here("processed_data","fungi_prepared.rds"))
t1p <- readRDS(here::here("processed_data","protist_prepared.rds"))
t1e <- readRDS(here::here("processed_data","eukaryote_prepared.rds"))
listx <- list(t1b, t1f, t1p, t1e)

colx <- 1
df <- as.data.frame(matrix(NA, nrow=29, ncol = 16))
df_p <- matrix(NA, nrow=29, ncol = 16)
order <- c(1,4,2,3)
for(i in 1:4){
 for (o in 1:4){
  df[,colx] <- paste0(Hmisc::format.pval(listx[[o]][[order[i]+5]]$p.adj,nsmall = 3, digits=3),"\n", listx[[o]][[order[i]+5]]$p.adj.signif)
  df_p[,colx] <- as.numeric(listx[[o]][[order[i]+5]]$p.adj,1)
  colx <- colx+1
 }
}

dfe <- t1b$richness_statistics[,c(1,2,5)]
names(dfe) <- c("exp","temp","days")
dft <- cbind(dfe,df)

ft <- flextable(dft) %>% 
 merge_v(j = "exp") %>% 
 merge_v(j = "temp") %>% 
 add_header_row(colwidths = c(3, 4,4,4,4), values = c("Storage conditions", "Richness","Ricness (dominant species)","Evenness","Distance to time 0")) %>% 
 valign(valign = "top") %>% 
 theme_box() %>% 
 fontsize(part = "all", size = 5) %>% 
 autofit()

read_docx() %>% 
 body_add_flextable(value = ft) %>% 
 print(target = here::here("tables","Supplementary_table2.docx"))


#### Supplementary Table 3 - permanova  ----
# table of statistics
t1b <- readRDS(here::here("processed_data","bacteria_prepared.rds"))
t1f <- readRDS(here::here("processed_data","fungi_prepared.rds"))
t1p <- readRDS(here::here("processed_data","protist_prepared.rds"))
t1e <- readRDS(here::here("processed_data","eukaryote_prepared.rds"))
listx <- list(t1b, t1f, t1p, t1e)

dfv <- data.frame(matrix(unlist(lapply(listx, function(x) paste0(x$pair_adonis$p.adjusted, " ", x$pair_adonis$sig))), ncol = 4))
dfv <- cbind(str_split_fixed(gsub("reference vs ","",t1b$pair_adonis[,1]),"_",3), dfv )
names(dfv) <- c("Exposure","Temperature","Storage time","Bacteria","Fungi","Protozoa","Eukaryotes")

dfv <- dfv %>% arrange(Exposure, as.numeric(Temperature), as.numeric(`Storage time`))

ft <- flextable(dfv) %>% 
 valign(valign = "top") %>% 
 theme_box() %>% 
 fontsize(part = "all", size = 10) %>% 
 autofit()

read_docx() %>% 
 body_add_flextable(value = ft) %>% 
 print(target = here::here("tables","supplementary_table3_permanova.docx"))


