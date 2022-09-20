# To reproduce some or all the analyses of this study:
#    NB: Be sure to read the information in the supplementary material of the publication
#
# The intention is that all data is available to jump into the analyses at any given step (as intermediate data is available)
# All you need to do is to download the R project including data from GitHub (https://github.com/tobiasgf/sample_storage).
# If you want to start from raw un-demultiplexed data, you also need to download the raw Illumina data from Dryad (see step 2A)
# If you want to start from sample wise demultiplexed paired reads (step 2B), you need to download these from ENA/SRA 
# 
# 2A. Demultiplex the raw fastq libraries (or download the demultiplexed reads from ENA, see 2B)
#    To demultiplex using the scripts,  dependencies are needed in addition to the r-libraries listed in the scripts:
#      Cutadapt (https://cutadapt.readthedocs.io/)
#      you also need the taglists and batchfile (in folder)
#    
#    Place the 9 paired R1 and R2 files (18 in total) in the same directory as the 10 taglists and batchfileDADA2.txt 
#      and the bash script demultiplex_miseq_data_paired_end.sh
#    
#    Run the script: bash demultiplex_miseq_data_paired_end.sh &>logfile.txt
#       All the demultiplexed files will now be in the directories
#       DADA2_SS (read pairs with true orientation) and
#       DADA2_AS (read pairs with reverse orientation) (see supplementary material for details)
#
#    For each of the 10 sequence pools:
#      Make a separate directory with two sub-directories named DADA2_SS and DADA2_AS
#      Place the relevant "true" oriented reads pairs in DADA2_SS
#      Place the relevant "reverse" oriented reads pairs in DADA2_AS
#      This can be done with the scripts:
#      cd DADA2_SS
#      find . -type f | awk -F_ '{system("mkdir -p ../"$1"/DADA2_SS ;mv "$0" ../"$1"/DADA2_SS/"$2"_"$3)}'
#      cd ../DADA2_AS
#      find . -type f | awk -F_ '{system("mkdir -p ../"$1"/DADA2_AS ;mv "$0" ../"$1"/DADA2_AS/"$2"_"$3)}'

# OR!!!!
# 2B. (alternative to 2A above) Download the demultiplexed read pairs from ENA/SRA
#   https://www.ebi.ac.uk/ena/browser/view/PRJEB56039
#   (or use the accessions from supplementary table supplementary_table_ENA_read_info.tsv)
#   
#   These sample pairs also need to be sorted into a directory system as mentioned above.
#   This can be done by placing all the ~1800 reads in one directory and the run this:
#   find . -type f | awk -F_ '{system("mkdir -p "$1"/"$2";mv "$0" "$1"/"$2"/"$3"_"$4)}'
#   rename all sub-directories called "true" to "DADA2_SS", and all sub-directories called "reverse" to "DADA2_AS", to be able to use the dada2 processing scripts supplied here.

# 3. Analyse the paired reads with DADA2
#     This can be done using the customized script supplied here, and expects all demultiplexed reads to be placed in raw_reads in ten sub-directories (bac1, bac2, ...), each having subdirectories called DADA2_SS (with paired reads in true orientation) and DADA2_AS (with paired reads in reverse orientation).
#    use primary_dada2_analyses.R