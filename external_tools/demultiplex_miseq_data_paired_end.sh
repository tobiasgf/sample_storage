### make directories
#DADA2_SS contains reads where forward primer was located in the R1 read
#DADA2_AS contains reads where forward primer was located in the R2 read

mkdir -p DADA2_SS
mkdir -p DADA2_AS

## Read from the batchlist file, which has 
## per line, the read1 filename, read2 filename 
## tag filename, the forward primer, the reverse primer, 
## and the minimum length 


while read INPUT_R1 INPUT_R2 TAGS PRIMER_F PRIMER_R MIN_LENGTH ; do
  # Define binaries, temporary files and output files
  CUTADAPT="$(which cutadapt) --discard-untrimmed --minimum-length ${MIN_LENGTH}"
  CUTADAPT2="$(which cutadapt) --minimum-length ${MIN_LENGTH}"

  REV_PRIMER_F="$(echo $PRIMER_F | rev | tr ATUGCYRKMBDHVN TAACGRYMKVHDBN)"
  REV_PRIMER_R="$(echo $PRIMER_R | rev | tr ATUGCYRKMBDHVN TAACGRYMKVHDBN)"
 
  while read TAG_NAME TAG_SEQ RTAG_SEQ; do

    FTFP="$TAG_SEQ$PRIMER_F"
    RTRP="$RTAG_SEQ$PRIMER_R"

    MINTAGOVERLAP=${#FTFP}
    if [ ${#FTFP} -gt ${#RTRP} ]; then
        MINTAGOVERLAP=${#RTRP}
    fi 

    ## straight strand (forward primer in R1)
    LOG="DADA2_SS/${TAG_NAME}.log"
    FINAL1_FASTQ="DADA2_SS/${TAG_NAME}_R1.fastq"
    FINAL2_FASTQ="DADA2_AS/${TAG_NAME}_R2.fastq"

    # Trim tags, forward & reverse primers
    ${CUTADAPT} -O $MINTAGOVERLAP -e 0 -g ${FTFP} -G ${RTRP} -o T1.fastq -p T2.fastq $INPUT_R1 $INPUT_R2 --cores 0
    #${CUTADAPT} -O $MINPRIMEROVERLAP -g ^${PRIMER_F} -G ^${PRIMER_R} -o $T3_FASTQ -p $T4_FASTQ $T1_FASTQ $T2_FASTQ     
    ${CUTADAPT2} -a ${REV_PRIMER_R} -A ${REV_PRIMER_F} -o $FINAL1_FASTQ -p $FINAL2_FASTQ T1.fastq T2.fastq --cores 0

    ## reverse order reads (forward primer in R2)
    LOG="DADA2_AS/${TAG_NAME}.log"
    FINAL1_FASTQ="DADA2_AS/${TAG_NAME}_R1.fastq"
    FINAL2_FASTQ="DADA2_AS/${TAG_NAME}_R2.fastq"

    # Trim tags, forward & reverse primers
    ${CUTADAPT} -O $MINTAGOVERLAP -e 0 -g ${RTRP} -G ${FTFP} -o T1.fastq -p T2.fastq $INPUT_R1 $INPUT_R2 --cores 0
    #${CUTADAPT} -O $MINPRIMEROVERLAP -g ^${PRIMER_R}  -G ^${PRIMER_F} -o $T3_FASTQ -p $T4_FASTQ $T1_FASTQ $T2_FASTQ 
    ${CUTADAPT2} -a ${REV_PRIMER_F} -A ${REV_PRIMER_R} -o $FINAL1_FASTQ -p $FINAL2_FASTQ T1.fastq T2.fastq --cores 0

    echo $TAG_NAME $TAG_SEQ $RTAG_SEQ

  done < <(cat ${TAGS})

  echo $INPUT_R1 $INPUT_R2 $TAGS $PRIMER_F $PRIMER_R $MIN_LENGTH

done < batchfileDADA2.txt
