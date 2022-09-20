#### prepare reads for submission to ENA ----

# This file is just for internal documentation
# this was done after the demultiplexing

# true/reverse tag added to all reads (i.e. name all paired reads uniquely, preferably so they can be sorted into directories relevant for separate analyses, e.g.)
# In all DADA2_SS directories: find . -type f | awk -F_ '{system("mv "$0" "$1"_true_"$2"_"$3)}'
# In all DADA2_AS directories: find . -type f | awk -F_ '{system("mv "$0" "$1"_reverse_"$2"_"$3)}'
# filenams are then e.g. :
# bac1_reverse_Blank2_R1.fastq
# bac1_reverse_Blank2_R2.fastq
# bac1_true_Blank2_R1.fastq
# bac1_true_Blank2_R2.fastq

# move all R1 reads to one directory
# gzip all reads: gz *
# move all R2 reads to one directory
# gzip all reads: gz *
# move controls etc to a separate dictory: mv */*Neg* controls
# move blanks etc to a separate dictory: mv */*lank* controls

# get md5 hash for all files:  md5 * > ../R1_md5.txt
# get md5 hash for all files:  md5 * > ../R2_md5.txt

#### upload reads to ENA with fta ----

# lftp webin2.ebi.ac.uk -u Webin-xxxxx PASSWORD #ENA Webin username and password

#### Register project ----
# Go to: https://www.ebi.ac.uk/ena/submit/webin/
# press "Register Study"
# add relevant information
# note down the Accession number: PRJEB56039	alt: ERP140961

#### Register samples ----
# Press "Register Samples"
# select "Download Spreadsheet to register samples"
# select "Environmental Checklists"
# select "GSC MIxS soil"
# tick off these optional fields "sample weight for DNA extraction",	"storage conditions (fresh/frozen/other)",	"soil horizon",	"soil type"
# press "next"
# press "Download tsv template"
# save to directory "ENA_submission"
# prepare sample information for submission

sampleinfo <- read.csv(here::here("basic_data","sampleinfo.txt"), sep = "\t", stringsAsFactors = F)
sampleinfo$batch[sampleinfo$batch == "new"] <- "new_pcr"
df1 <- data.frame(
 tax_id	= "410658",
 scientific_name = "soil metagenome",
 sample_alias = sampleinfo$sample,
 sample_title = sampleinfo$sample,
 sample_description = paste0(sampleinfo$exposure,"_temp",sampleinfo$temperature,"_time",sampleinfo$harvest,ifelse(sampleinfo$batch %in% c("lts","original"),"",paste0("_",sampleinfo$batch))),
 "project name" = "PRJEB56039", check.names = F,
 "sequencing method" = "Illumina",
 "collection date" = "2017-08-31",
 "geographic location (country and/or sea)" = "Denmark",
 "geographic location (latitude)" = 55.95689,
 "geographic location (longitude)" = 12.27209,
 "geographic location (region and locality)" = "Denmark, Sjælland, Strødam",
 depth = 0,
 "broad-scale environmental context" = "tempreate biome",
 "local environmental context"	= "temperate woodland biome",
 "environmental medium"	= "beech forest soil",
 elevation	= 20,
 "sample weight for DNA extraction"	= 0.25,
 "storage conditions (fresh/frozen/other)" = "frozen",
 "soil horizon"	= "A Horizon",
 "soil type" = "Luvisol")

write_tsv(df1, here::here("ENA_submission","sample_info2.tsv"))

# copy columns into downloaded template and upload

#### submit samples ----

# Press "Submit Reads"
# press "Download spreadsheet template for Read submission"
# select "Submit paired reads using two Fastq files"
# press "next" and "Download TSV Template"

# Prepare data for template
R1_md5 <- read_tsv(here::here("raw_reads/all_reads","R1_md5.txt"), col_names = F)
R2_md5 <- read_tsv(here::here("raw_reads/all_reads","R2_md5.txt"), col_names = F)

df2 <- data.frame(
 sample = str_split_fixed(R1_md5$X1,"_",4)[,3],
 study = "PRJEB56039",
 instrument_model = "Illumina MiSeq",
 library_name = NA,
 library_source = "METAGENOMIC",
 library_selection = "PCR",
 library_strategy = "AMPLICON",
 library_layout = "PAIRED",
 forward_file_name = stringr::str_extract(string = R1_md5$X1, pattern = "(?<=\\().*(?=\\))"),
 forward_file_md5 = gsub(".*= ","",R1_md5$X1),
 reverse_file_name = stringr::str_extract(string = R2_md5$X1, pattern = "(?<=\\().*(?=\\))"),
 reverse_file_md5 = gsub(".*= ","",R2_md5$X1))

df1x <- df1 %>% select(sample_alias, sample_description)
df2x <- left_join(df2, df1x, by = c("sample" = "sample_alias"))
df2x$library_name <- paste0(gsub("_S.*","",df2x$forward_file_name),"_",df2x$sample_description)
df2x2 <- df2x %>% select(-sample_description)
write_tsv(df2x2, here::here("ENA_submission","read_info.tsv"))
# copy columns into downloaded template and upload

#read in the accession numbers and merge with readinfo for supplementary table.
read_info <- read_csv(here::here("ENA_submission","run-files-2022-09-20T11_28_36.csv"))
read_info2 <- read_info %>% select(fileName, id)
df2x3 <- df2x2 %>% left_join(read_info2, by = c("forward_file_name" = "fileName")) %>% rename(forward_file_accession = id) %>% left_join(read_info2, by = c("reverse_file_name" = "fileName")) %>% rename(reverse_file_accession = id)

write_tsv(df2x3, here::here("tables","supplementary_table5_ENA_read_info.tsv"))
