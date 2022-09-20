# It is expected that all demultiplexed reads are be placed in directory "raw_reads" in ten sub-directories ("bac1", "bac2", ...), each having subdirectories called "DADA2_SS" (with paired reads in true orientation) and "DADA2_AS" (with paired reads in reverse orientation)

# you need to set the correct path to your folder

library(here)
library(dada2)

source(here::here("R","dada2_processing.r"))

setwd("~/THE_PATH/bac1")
dada2_tgf(trunc_reads = TRUE, cut_r1 = 240, cut_r2 = 210, minimum_reads = 10, plotting = FALSE, filter_r1 = TRUE, filter_r2 = TRUE, run_dada2_r1 = TRUE, run_dada2_r2 = TRUE, merge_antisense = TRUE, get_stats = TRUE, use_sickle = F)
tab <- readRDS(here::here("raw_reads/bac1","seqtab.nochim_Both_RDS"))
saveRDS(tab, here::here("asv_tables","bac1_dada2_otutable.rds"))

setwd("~/THE_PATH/bac2")
dada2_tgf(trunc_reads = TRUE, cut_r1 = 240, cut_r2 = 210, minimum_reads = 10, plotting = FALSE, filter_r1 = TRUE, filter_r2 = TRUE, run_dada2_r1 = TRUE, run_dada2_r2 = TRUE, merge_antisense = TRUE, get_stats = TRUE, use_sickle = F)
tab <- readRDS(here::here("raw_reads/bac2","seqtab.nochim_Both_RDS"))
saveRDS(tab, here::here("asv_tables","bac2_dada2_otutable.rds"))

setwd("~/THE_PATH/bac3")
dada2_tgf(trunc_reads = TRUE, cut_r1 = 240, cut_r2 = 210, minimum_reads = 10, plotting = FALSE, filter_r1 = TRUE, filter_r2 = TRUE, run_dada2_r1 = TRUE, run_dada2_r2 = TRUE, merge_antisense = TRUE, get_stats = TRUE, use_sickle = F)
tab <- readRDS(here::here("raw_reads/bac3","seqtab.nochim_Both_RDS"))
saveRDS(tab, here::here("asv_tables","bac3_dada2_otutable.rds"))

setwd("~/THE_PATH/euk1")
dada2_tgf(trunc_reads = TRUE, cut_r1 = 230, cut_r2 = 190, minimum_reads = 10, plotting = FALSE, filter_r1 = TRUE, filter_r2 = TRUE, run_dada2_r1 = TRUE, run_dada2_r2 = TRUE, merge_antisense = TRUE, get_stats = TRUE, use_sickle = F)
tab <- readRDS(here::here("raw_reads/euk1","seqtab.nochim_Both_RDS"))
saveRDS(tab, here::here("asv_tables","euk1_dada2_otutable.rds"))

setwd("~/THE_PATH/euk2")
dada2_tgf(trunc_reads = TRUE, cut_r1 = 230, cut_r2 = 190, minimum_reads = 10, plotting = FALSE, filter_r1 = TRUE, filter_r2 = TRUE, run_dada2_r1 = TRUE, run_dada2_r2 = TRUE, merge_antisense = TRUE, get_stats = TRUE, use_sickle = F)
tab <- readRDS(here::here("raw_reads/euk1","seqtab.nochim_Both_RDS"))
saveRDS(tab, here::here("asv_tables","euk1_dada2_otutable.rds"))

setwd("~/THE_PATH/pro1")
dada2_tgf(trunc_reads = TRUE, cut_r1 = 190, cut_r2 = 170, minimum_reads = 10, plotting = FALSE, filter_r1 = TRUE, filter_r2 = TRUE, run_dada2_r1 = TRUE, run_dada2_r2 = TRUE, merge_antisense = TRUE, get_stats = TRUE, use_sickle = F)
tab <- readRDS(here::here("raw_reads/pro1","seqtab.nochim_Both_RDS"))
saveRDS(tab, here::here("asv_tables","pro1_dada2_otutable.rds"))

setwd("~/THE_PATH/pro1")
dada2_tgf(trunc_reads = TRUE, cut_r1 = 190, cut_r2 = 170, minimum_reads = 10, plotting = FALSE, filter_r1 = TRUE, filter_r2 = TRUE, run_dada2_r1 = TRUE, run_dada2_r2 = TRUE, merge_antisense = TRUE, get_stats = TRUE, use_sickle = F)
tab <- readRDS(here::here("raw_reads/pro2","seqtab.nochim_Both_RDS"))
saveRDS(tab, here::here("asv_tables","pro2_dada2_otutable.rds"))

setwd("~/THE_PATH/fun1")
dada2_tgf(trunc_reads = FALSE, cut_r1 = NA, cut_r2 = NA, minimum_reads = NA, plotting = FALSE, filter_r1 = NA, filter_r2 = NA, run_dada2_r1 = NA, run_dada2_r2 = NA, merge_antisense = NA, get_stats = NA, use_sickle = TRUE)
tab <- readRDS(here::here("raw_reads/fun1","seqtab.nochim_Both_RDS"))
saveRDS(tab, here::here("asv_tables","fun1_dada2_otutable.rds"))

setwd("~/THE_PATH/fun2")
dada2_tgf(trunc_reads = FALSE, cut_r1 = NA, cut_r2 = NA, minimum_reads = NA, plotting = FALSE, filter_r1 = NA, filter_r2 = NA, run_dada2_r1 = NA, run_dada2_r2 = NA, merge_antisense = NA, get_stats = NA, use_sickle = TRUE)
tab <- readRDS(here::here("raw_reads/fun2","seqtab.nochim_Both_RDS"))
saveRDS(tab, here::here("asv_tables","fun2_dada2_otutable.rds"))

setwd("~/THE_PATH/fun3")
dada2_tgf(trunc_reads = FALSE, cut_r1 = NA, cut_r2 = NA, minimum_reads = NA, plotting = FALSE, filter_r1 = NA, filter_r2 = NA, run_dada2_r1 = NA, run_dada2_r2 = NA, merge_antisense = NA, get_stats = NA, use_sickle = TRUE)
tab <- readRDS(here::here("raw_reads/fun3","seqtab.nochim_Both_RDS"))
saveRDS(tab, here::here("asv_tables","fun3_dada2_otutable.rds"))

