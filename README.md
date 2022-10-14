# Sample Storage - analyses
___
This repository (R project) contains all data and scripts necessary to run the analyses and produce the figures from the study Treated like dirt: Robust forensic and ecological inferences from soil eDNA after challenging sample storage (Frøslev et al. (2022) - published in Environmental DNA).
The intention is that all data is available to jump into the analyses at any given step (as intermediate data is available). All you need to do is to download the R project including data from this GitHub repository (https://github.com/tobiasgf/sample_storage).
If you wish to start from raw un-demultiplexed data, you also need to download the raw Illumina data from Dryad (https://doi.org/10.5061/dryad.k0p2ngfbs)
If you want to start from sample-wise demultiplexed paired reads, you need to download these from ENA/SRA (https://www.ebi.ac.uk/ena/browser/view/PRJEB56039).
Most of the analyses are carried out in R (with standard packages), but there are also some command line steps. To be able to run all analyses, you need some additional tools installed (ITSx, cutadapt, blastx, sickle)

*An overview and description of the directories and scripts in this repository*  

Directories|contents  
--- | ---   
asv_tables|contains all the 10 tables from the primary dada2 analyses (see above), as well as the downstream versions, where datasets are combined   
asv_tables_referencedata|contais the four asv_tables from the reference data (and some information on replicates for two of them)
basic_data|contains data for co2 measurements, dna concentration, and some basic sample meta data
final_plots|contains plots from r as well as the publication ready plots as pdf (with some editing done in Adobe Illustrator)
partial_plots|contains marker wise plots (that are later combined into final compound plots)
processed_data|contains marker-wise stats and results from different analyses
R|contains some scripts for statistical analyses and plot making
raw_reads|the directory for the raw fastq reads and demultiplexed reads (both needs downloading), as well as datafile necessary for the demultiplexing
ressources|contains some external ressources for the habitat classification as well as reference datasets for the taxonomic annotation
tables|output tables for publication
taxonomic_assignments| contains info and tables related to the taxonomic annotation of reads
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
https://www.ebi.ac.uk/ena/browser/view/PRJEB56039


Link to Dryad containing the un-demultiplexed paired read files:
https://doi.org/10.5061/dryad.k0p2ngfbs

