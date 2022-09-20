# Sample Storage - analyses
___
This repository (R project) contains all data and scripts necessary to run the analyses and produce the figures from the study Treated like dirt: Robust forensic and ecological inferences from soil eDNA after challenging sample storage (Frøslev et al. (2022) - published in Environmental DNA).
The intention is that all data is available to jump into the analyses at any given step (as intermediate data is available). All you need to do is to download the R project including data from this GitHub repository (https://github.com/tobiasgf/sample_storage).
If you wish to start from raw un-demultiplexed data, you also need to download the raw Illumina data from Dryad (LINK)
If you want to start from sample-wise demultiplexed paired reads, you need to download these from ENA/SRA (https://www.ebi.ac.uk/ena/browser/view/PRJEB56039).
Most of the analyses are carried out in R (with standard packages), but there are also some command line steps. To be able tu run all analyses, you need some additional tools installed (ITSx, cutadapt, snasltx, sickle)

*An overview and description of the directories and scripts in this repository*  

Directories|contents  
--- | ---   
asv_tables|contains all the 10 tables from the primary dada2 analyses (see above), as well as the downstream versions, where datasets are combined   
Asv_tables_referencedata|contais the four asv_tables from the reference data (and some information on replicates for two of them)
Basic_data|contains data for co2 measurements, dna concentration, and some basic sample meta data
Final_plots|contains plots from r as well as the publication ready plots as pdf (with some editing done in Adobe Illustrator)
Partial_plots|contains marker wise plots (that are later combined into final compound plots)
Processed_data|contains marker-wise stats and results from different analyses
R|contains some scripts for statistical analyses and plot making
Raw_reads|the directory for the raw fastq reads and demultiplexed reads (both needs downloading), as well as datafile necessary for the demultiplexing
Ressources|contains some external ressources for the habitat classification as well as reference datasets for the taxonomic annotation
Tables|output tables for publication
Taxonomic_assignments| contains info and tables related to the taxonomic annotation of reads
ENA_submission|files associated with the preparation and upload of reads to ENA

The following scripts are used in the listed order to do the analyses, statistics and plotting of the study  

script|function
--- | ---
Analyses_sample_storage.Rproj|the r-project file
reproducing_analyses_start here.R|some infor on initial steps, if you wish to download either raw un-demultiplexed reads (from dryad) an demultiplex these, OR downolad demultiplexed sample-wise paired end reads from ENA/SRA
primary_dada2_analyses|this script is used to construct the initial 10 asv tables with dada2
construct_habitat_classes_for_supervised_classification.R|used to perform a supervised classification of the reference data from the Biowide study
bacteria_analyses.R|do all the analyses related to the bacterial data
fungi_analyses.R|do all the analyses related to the fungi data
protist_analyses.R|do all the analyses related to the protist data
eukaryote_analyses.R| do all the analyses related to the eukaryote data
produce_plots_and_tables.R| produce plots and tables
prepare_and_submit_reads_to_ENA.R| a documentation of the demultiplexing of the reads for ENA submission, and a bit of further information


Link to the ENA project containing demultiplexed samples:
prepare_and_submit_reads_to_ENA.R

Link to Dryad containing the un-demultiplexed paired read files:
