dada2_tgf <- function(trunc_reads = TRUE, cut_r1 = 230, cut_r2 = 180, minimum_reads = 10, plotting = FALSE, filter_r1 = TRUE, filter_r2 = TRUE, run_dada2_r1 = TRUE, run_dada2_r2 = TRUE, merge_antisense = TRUE, get_stats = TRUE, use_sickle = F){
 message("Running dada2 on demultiplexed fastq files in two directories (DADA2_SS and DADA2_AS)")
 message(paste("Discarding samples with less than", minimum_reads))
 if(trunc_reads){
  message(paste("Truncating R1 reads to", cut_r1, "nucleotides"))
  message(paste("Truncating R2 reads to", cut_r2, "nucleotides"))
 } else {
  message("Not truncating reads")
 }
 
 require(dada2)
 require(readr)
 require(dplyr)
 require(digest)
 library(data.table)
 
 start.time <- Sys.time()
 
 if(!use_sickle){
  print("not using sickle")
  
  pwd <- getwd()
  main_path <- pwd
  
  if(!file_test("-d", file.path(main_path, "DADA2_SS"))) stop("DADA2_SS directory not in current directory!")
  if(!file_test("-d", file.path(main_path, "DADA2_AS"))) stop("DADA2_AS directory not in current directory!")
  
  
  
  #filtering of the Sense-reads:
  path <- file.path(main_path, "DADA2_SS")
  
  if(filter_r1){
   fns <- list.files(path)
   fastqs <- fns[grepl("fastq", fns)]
   fastqs <- sort(fastqs)
   
   fnFs <- fastqs[grepl("_R1\\.fastq", fastqs)]
   fnRs <- fastqs[grepl("_R2\\.fastq", fastqs)]
   sample.names <- sapply(strsplit(fnFs, "_R1\\.fastq"), `[`, 1) # Extract sample names, assuming filenames have format: SAMPLENAME_XXX.fastq
   
   fnFs <- file.path(path, fnFs)
   fnRs <- file.path(path, fnRs)
   filt_path <- file.path(path, "filtered") # Place filtered files in filtered/ subdirectory
   
   if(!file_test("-d", filt_path)) dir.create(filt_path)
   filtFs <- file.path(filt_path, paste0(sample.names, "_F_filt.fastq.gz"))
   filtRs <- file.path(filt_path, paste0(sample.names, "_R_filt.fastq.gz"))
   
   indexX <- vector()
   for(i in seq_along(fnFs)) {
    if (file.info(fnFs[i])$size == 0 | file.info(fnRs[i])$size == 0) {
     print(paste(fnFs[i], "empty."))
     indexX <- c(indexX,i)
    }
   }
   
   if (length(indexX)>0) {
    fnFs <- fnFs[-indexX]
    fnRs <- fnRs[-indexX]
    filtFs <- filtFs[-indexX]
    filtRs <- filtRs[-indexX]
   }
   
   if(trunc_reads){
    print(paste0("filtering... truncating R1 to ",cut_r1," and R2 to ",cut_r2))
    SSout <- filterAndTrim(fnFs, filtFs, fnRs, filtRs,minLen=10,matchID=TRUE, maxN=0, maxEE=c(2,2),truncQ=2,verbose=TRUE,truncLen=c(cut_r1,cut_r2),compress=TRUE, multithread=TRUE)
   } else {
    print(paste0("filtering without truncation"))
    SSout <- filterAndTrim(fnFs, filtFs, fnRs, filtRs,minLen=10,matchID=TRUE, maxN=0, maxEE=c(2,2),truncQ=2,verbose=TRUE,compress=TRUE, multithread=TRUE)
   }
   
   message("taking a look at how many reads passed the filter")
   SSout
   
   #removing samples with too few reads (set in the variable minimum_reads in the start)
   low_read_files <- SSout[,"reads.out"] < 10 & SSout[,"reads.out"] != 0 # Or other cutoff
   print(paste0("Removing these files as they have fewer than ",minimum_reads," reads:\n"))
   print(filtFs[low_read_files])
   print(filtRs[low_read_files])
   file.remove(filtFs[low_read_files])
   file.remove(filtRs[low_read_files])
   
  }
  
  if(run_dada2_r1){
   
   #Processing the set of files containing the forward primer in the R1 reads (the sense reads):
   filt_path <- file.path(main_path, "DADA2_SS/filtered") 
   fns <- list.files(filt_path)
   fastqs <- fns[grepl(".fastq.gz", fns)]
   fastqs <- sort(fastqs)
   fnFs <- fastqs[grepl("_F_", fastqs)]
   fnRs <- fastqs[grepl("_R_", fastqs)]
   fSSsample.names <- sapply(strsplit(fnFs, "_F_"), `[`, 1)
   rSSsample.names <- sapply(strsplit(fnRs, "_R_"), `[`, 1)
   filtFs <- file.path(filt_path, fnFs)
   filtRs <- file.path(filt_path, fnRs)
   
   errF <- learnErrors(filtFs, multithread=TRUE)
   errR <- learnErrors(filtRs, multithread=TRUE)
   plotErrors(errF, nominalQ=TRUE) 
   plotErrors(errR, nominalQ=TRUE)
   
   derepFs <- derepFastq(filtFs, verbose=TRUE)
   derepRs <- derepFastq(filtRs, verbose=TRUE)
   names(derepFs) <- fSSsample.names
   names(derepRs) <- rSSsample.names
   
   # ---Sample Inference---
   
   # We are now ready to apply the core sample inference algorithm to the dereplicated data.
   SSdadaFs <- dada(derepFs, err=errF, multithread = TRUE)
   SSdadaRs <- dada(derepRs, err=errR, multithread = TRUE)
   SSmergers <- mergePairs(SSdadaFs, derepFs, SSdadaRs, derepRs, verbose=TRUE,minOverlap = 5)
   # ---Construct sequence table---
   
   # We can now construct an amplicon sequence variant table (ASV) table, a higher-resolution version of the OTU table produced by traditional methods.
   seqtab_SS <- makeSequenceTable(SSmergers[names(SSmergers)]) # The sequence table is a matrix with rows corresponding to (and named by) the samples, and columns corresponding to (and named by) the sequence variants. 
   seqtab.nochim_SS <- removeBimeraDenovo(seqtab_SS, verbose=TRUE) # The core dada method corrects substitution and indel errors, but chimeras remain. Fortunately, the accuracy of sequence variants after denoising makes identifying chimeric ASVs simpler than when dealing with fuzzy OTUs. Chimeric sequences are identified if they can be exactly reconstructed by combining a left-segment and a right-segment from two more abundant "parent" sequences.
   
   stSS <- file.path(main_path,"seqtab_SS_RDS")
   stnsSS <- file.path(main_path,"seqtab.nochim_SS_RDS")
   saveRDS(seqtab_SS,stSS)
   saveRDS(seqtab.nochim_SS,stnsSS)
  }
  
  if(filter_r2){
   #filtering of the antiSense-reads:
   path <- file.path(main_path, "DADA2_AS") 
   fns <- list.files(path)
   fastqs <- fns[grepl("fastq", fns)]
   fastqs <- sort(fastqs)
   fnFs <- fastqs[grepl("_R2\\.fastq", fastqs)] # Reverse direction compared to above for the "sense reads". In practice the reads are here complement reversed to be in the same orientation as the "sense" reads.
   fnRs <- fastqs[grepl("_R1\\.fastq", fastqs)] # See above
   sample.names <- sapply(strsplit(fnFs, "_R2\\.fastq"), `[`, 1)
   fnFs <- file.path(path, fnFs)
   fnRs <- file.path(path, fnRs)
   filt_path <- file.path(path, "filtered")
   
   if(!file_test("-d", filt_path)) dir.create(filt_path)
   filtFs <- file.path(filt_path, paste0(sample.names, "_F_filt.fastq.gz"))
   filtRs <- file.path(filt_path, paste0(sample.names, "_R_filt.fastq.gz"))
   
   indexX <- vector()
   for(i in seq_along(fnFs)) {
    if (file.info(fnFs[i])$size == 0 | file.info(fnRs[i])$size == 0) {
     print(paste(fnFs[i], "empty."))
     indexX <- c(indexX,i)
    }
   }
   
   if (length(indexX)>0) {
    fnFs <- fnFs[-indexX]
    fnRs <- fnRs[-indexX]
    filtFs <- filtFs[-indexX]
    filtRs <- filtRs[-indexX]
   }
   
   if(trunc_reads){
    print(paste0("filtering... truncating R1 to ",cut_r1," and R2 to ",cut_r2))
    ASout <- filterAndTrim(fnFs, filtFs, fnRs, filtRs,minLen=10,matchID=TRUE, maxN=0, maxEE=c(2,2),truncQ=2,verbose=TRUE,truncLen=c(cut_r2,cut_r1),compress=TRUE, multithread=TRUE)
   } else {
    print(paste0("filtering without truncation"))
    ASout <- filterAndTrim(fnFs, filtFs, fnRs, filtRs,minLen=10,matchID=TRUE, maxN=0, maxEE=c(2,2),truncQ=2,verbose=TRUE,compress=TRUE, multithread=TRUE)
   }
   
   #taking a look at how many reads passed the filter
   ASout
   
   #removing samples with too few reads (set in the variable minimum_reads in the start)
   low_read_files <- ASout[,"reads.out"] < 10 & ASout[,"reads.out"] != 0 # Or other cutoff
   print(paste0("Removing these files as they have fewer than ",minimum_reads," reads:\n"))
   print(filtFs[low_read_files])
   print(filtRs[low_read_files])
   file.remove(filtFs[low_read_files])
   file.remove(filtRs[low_read_files])
  }
  
  if(run_dada2_r2){
   #Then DADA2 processing of "the antisense" reads (same as for SS):
   filt_path <- file.path(main_path, "DADA2_AS/filtered") 
   fns <- list.files(filt_path)
   fastqs <- fns[grepl(".fastq.gz", fns)]
   fastqs <- sort(fastqs)
   fnFs <- fastqs[grepl("_F_", fastqs)]
   fnRs <- fastqs[grepl("_R_", fastqs)]
   fASsample.names <- sapply(strsplit(fnFs, "_F_"), `[`, 1)
   rASsample.names <- sapply(strsplit(fnRs, "_R_"), `[`, 1)
   filtFs <- file.path(filt_path, fnFs)
   filtRs <- file.path(filt_path, fnRs)
   errF <- learnErrors(filtFs, multithread=TRUE)
   errR <- learnErrors(filtRs, multithread=TRUE)
   derepFs <- derepFastq(filtFs, verbose=TRUE)
   derepRs <- derepFastq(filtRs, verbose=TRUE)
   names(derepFs) <- fASsample.names
   names(derepRs) <- rASsample.names
   ASdadaFs <- dada(derepFs, err=errF, multithread = TRUE)
   ASdadaRs <- dada(derepRs, err=errR, multithread = TRUE)
   ASmergers <- mergePairs(ASdadaFs, derepFs, ASdadaRs, derepRs, verbose=TRUE,minOverlap = 5)
   seqtab_AS <- makeSequenceTable(ASmergers[names(ASmergers)])
   seqtab.nochim_AS <- removeBimeraDenovo(seqtab_AS, verbose=TRUE)
   stAS <- file.path(main_path,"seqtab_AS_RDS")
   stnsAS <- file.path(main_path,"seqtab.nochim_AS_RDS")
   saveRDS(seqtab_AS,stAS)
   saveRDS(seqtab.nochim_AS,stnsAS)
   
  }
  
  if(merge_antisense){
   stAS <- file.path(main_path,"seqtab_AS_RDS")
   stnsAS <- file.path(main_path,"seqtab.nochim_AS_RDS")
   stSS <- file.path(main_path,"seqtab_SS_RDS")
   stnsSS <- file.path(main_path,"seqtab.nochim_SS_RDS")
   seqtab.nochim_AS <- readRDS(stnsAS)
   seqtab.nochim_SS <- readRDS(stnsSS)
   seqtab_AS <- readRDS(stAS)
   seqtab_SS <- readRDS(stSS)
   total_sumtable <- mergeSequenceTables(seqtab_SS,seqtab_AS, repeats = "sum")
   total_nochim_sumtable <- mergeSequenceTables(seqtab.nochim_SS,seqtab.nochim_AS, repeats = "sum")
   stBoth <- file.path(main_path,"seqtab_Both_RDS")
   stnsBoth <- file.path(main_path,"seqtab.nochim_Both_RDS")
   saveRDS(total_sumtable,stBoth)
   saveRDS(total_nochim_sumtable,stnsBoth)
   
  }
  
  if(get_stats){
   #Get statistics
   getN <- function(x) sum(getUniques(x))
   SSout <- as.data.frame(SSout)
   SSout$sample <-  gsub("_R.\\.fastq","",row.names(SSout))
   #SSout <- as.data.frame(SSout)
   r1 <- data.frame(sample=names(SSdadaFs), denoisedF=sapply(SSdadaFs, getN),stringsAsFactors=FALSE)
   r2 <- data.frame(sample=names(SSdadaRs), denoisedR=sapply(SSdadaRs, getN),stringsAsFactors=FALSE)
   r3 <- data.frame(sample=names(SSmergers), merged=sapply(SSmergers, getN),stringsAsFactors=FALSE)
   r4 <- data.frame(sample=row.names(seqtab.nochim_SS), nonchim=rowSums(seqtab.nochim_SS),stringsAsFactors=FALSE)
   
   SStrack <- left_join(SSout, r1, by='sample') %>% left_join(., r2, by='sample') %>% left_join(., r3, by='sample') %>% left_join(., r4, by='sample')
   row.names(SStrack) <- SStrack$sample
   SStrack <- SStrack[,-3]
   
   #saveRDS(SStrack, "read_track_SS_RDS")
   #write.table(SStrack,paste0(pwd, "read_track_SS_.txt"),col.names = NA, quote=FALSE,sep="\t")
   
   ASout <- as.data.frame(ASout)
   ASout$sample <-  gsub("_R.\\.fastq","",row.names(ASout))
   #ASout <- as.data.frame(ASout)
   r1 <- data.frame(sample=names(ASdadaFs), denoisedF=sapply(ASdadaFs, getN),stringsAsFactors=FALSE)
   r2 <- data.frame(sample=names(ASdadaRs), denoisedR=sapply(ASdadaRs, getN),stringsAsFactors=FALSE)
   r3 <- data.frame(sample=names(ASmergers), merged=sapply(ASmergers, getN),stringsAsFactors=FALSE)
   r4 <- data.frame(sample=row.names(seqtab.nochim_AS), nonchim=rowSums(seqtab.nochim_AS),stringsAsFactors=FALSE)
   
   AStrack <- left_join(ASout, r1, by='sample') %>% left_join(., r2, by='sample') %>% left_join(., r3, by='sample') %>% left_join(., r4, by='sample')
   row.names(AStrack) <- AStrack$sample
   AStrack <- AStrack[,-3]
   
   #saveRDS(AStrack, "read_track_AS_RDS")
   #write.table(AStrack,paste0(pwd,"read_track_AS_",appendix,".txt"),col.names = NA, quote=FALSE,sep="\t")
   
   total_stat <- data.frame(merged_reads=rowSums(total_sumtable),merged_nochim=rowSums(total_nochim_sumtable))
   #saveRDS(total_stat, "total_stat_RDS")
   #write.table(total_stat,paste0(pwd,"total_stat_",appendix,".txt"),col.names = NA, quote=FALSE,sep="\t")
   
   otu_stat <- data.frame(ss=ncol(seqtab_SS), ss_nochim=ncol(seqtab.nochim_SS), as=ncol(seqtab_AS), as_nochim=ncol(seqtab.nochim_AS), total=ncol(total_sumtable), total_nochim=ncol(total_nochim_sumtable))
   #write.table(otu_stat,paste0(pwd,"otu_stat_",appendix,".txt"),col.names = NA, quote=FALSE,sep="\t")
  } else {
   SStrack <- NULL
   AStrack <- NULL
   total_stat <- NULL
   otu_stat <- NULL
  }
  
  end.time <- Sys.time()
  time.taken <- end.time - start.time
  settings <- list(trunc_reads = trunc_reads, cut_r1 = cut_r1, cut_r2 = cut_r2, minimum_reads = minimum_reads, plotting = plotting)
  result <- list(dada2_table_raw = total_sumtable, dada2_table_nochimeras = total_nochim_sumtable, SStrack = SStrack, AStrack = AStrack, total_stat = total_stat, otu_stat = otu_stat, time.taken = time.taken, settings = settings)
  return(result)
 }
 
 if(use_sickle){
  #Filtering the "Sense"-reads.
  path <- file.path(main_path, "DADA2_SS")
  fns <- list.files(path)
  fastqs <- fns[grepl("fastq$", fns)]
  fastqs <- sort(fastqs)
  fnFs <- fastqs[grepl("_R1.", fastqs)]
  fnRs <- fastqs[grepl("_R2.", fastqs)]
  sample.names <- sapply(strsplit(fnFs, "_R1"), `[`, 1)
  fnFs <- file.path(path, fnFs)
  fnRs <- file.path(path, fnRs)
  match_path <- file.path(path, "filtered")
  if(!file_test("-d", match_path)) dir.create(match_path)
  filtFs <- file.path(match_path, paste0(sample.names, "_F_filtered.fastq"))
  filtRs <- file.path(match_path, paste0(sample.names, "_R_filtered.fastq"))
  
  for(i in seq_along(fnFs)) {
   if (file.info(fnFs[i])$size == 0 | file.info(fnRs[i])$size == 0) {
    print(paste(fnFs[i], "empty.")) } else {
     print(paste("processing", fnFs[i]))
     commandX <- paste0("/Users/cdr273/bin/sickle se -l 100 -q 28 -x -t sanger -f '",fnFs[i],"' -o '",filtFs[i],"'")
     commandY <- paste0("/Users/cdr273/bin/sickle se -l 100 -q 28 -x -t sanger -f '",fnRs[i],"' -o '",filtRs[i],"'")
     print(commandX)
     system(commandX)
     print(commandY)
     system(commandY)
    }
  }
  
  
  #Filtering the "Anti-Sense"-reads.
  path <- file.path(main_path, "DADA2_AS")
  fns <- list.files(path)
  fastqs <- fns[grepl("fastq$", fns)]
  fastqs <- sort(fastqs)
  fnFs <- fastqs[grepl("_R2.", fastqs)] # Reverse direction compared to above for the "sense reads". In practice the reads are here complement reversed to be in the same orientation as the "sense" reads.
  fnRs <- fastqs[grepl("_R1.", fastqs)] # See above
  sample.names <- sapply(strsplit(fnFs, "_R2"), `[`, 1)
  fnFs <- file.path(path, fnFs)
  fnRs <- file.path(path, fnRs)
  match_path <- file.path(path, "filtered")
  if(!file_test("-d", match_path)) dir.create(match_path)
  filtFs <- file.path(match_path, paste0(sample.names, "_F_filtered.fastq"))
  filtRs <- file.path(match_path, paste0(sample.names, "_R_filtered.fastq"))
  
  for(i in seq_along(fnFs)) {
   if (file.info(fnFs[i])$size == 0 | file.info(fnRs[i])$size == 0) {
    print(paste(fnFs[i], "empty.")) } else {
     print(paste("processing", fnFs[i]))
     commandX <- paste0("/Users/cdr273/bin/sickle se -l 100 -q 28 -x -t sanger -f '",fnFs[i],"' -o '",filtFs[i],"'")
     commandY <- paste0("/Users/cdr273/bin/sickle se -l 100 -q 28 -x -t sanger -f '",fnRs[i],"' -o '",filtRs[i],"'")
     print(commandX)
     system(commandX)
     print(commandY)
     system(commandY)
    }
  }
  
  
  #Matching the "Sense"-reads.
  path <- file.path(main_path, "DADA2_SS/filtered")
  fns <- list.files(path)
  fastqs <- fns[grepl("fastq$", fns)]
  fastqs <- sort(fastqs)
  fnFs <- fastqs[grepl("_F_", fastqs)]
  fnRs <- fastqs[grepl("_R_", fastqs)]
  sample.names <- sapply(strsplit(fnFs, "_F_"), `[`, 1)
  fnFs <- file.path(path, fnFs)
  fnRs <- file.path(path, fnRs)
  match_path <- file.path(path, "matched")
  if(!file_test("-d", match_path)) dir.create(match_path)
  filtFs <- file.path(match_path, paste0(sample.names, "_F_matched.fastq.gz"))
  filtRs <- file.path(match_path, paste0(sample.names, "_R_matched.fastq.gz"))
  for(i in seq_along(fnFs)) {
   if (file.info(fnFs[i])$size == 0 | file.info(fnRs[i])$size == 0) {
    print(paste(fnFs[i], "empty.")) } else {
     print(paste("processing", fnFs[i]))
     fastqPairedFilter(c(fnFs[i], fnRs[i]), c(filtFs[i], filtRs[i]),
                       matchIDs=TRUE)
    }
  }
  
  #Matching the "Anti-Sense"-reads.
  path <- file.path(main_path, "DADA2_AS/filtered")
  fns <- list.files(path)
  fastqs <- fns[grepl("fastq$", fns)]
  fastqs <- sort(fastqs)
  fnFs <- fastqs[grepl("_F_", fastqs)] # Reverse direction compared to above for the "sense reads". In practice the reads are here complement reversed to be in the same orientation as the "sense" reads.
  fnRs <- fastqs[grepl("_R_", fastqs)] # See above
  sample.names <- sapply(strsplit(fnFs, "_F_"), `[`, 1)
  fnFs <- file.path(path, fnFs)
  fnRs <- file.path(path, fnRs)
  match_path <- file.path(path, "matched")
  if(!file_test("-d", match_path)) dir.create(match_path)
  filtFs <- file.path(match_path, paste0(sample.names, "_F_matched.fastq.gz"))
  filtRs <- file.path(match_path, paste0(sample.names, "_R_matched.fastq.gz"))
  for(i in seq_along(fnFs)) {
   if (file.info(fnFs[i])$size == 0 | file.info(fnRs[i])$size == 0) {
    print(paste(fnFs[i], "empty.")) } else {
     print(paste("processing", fnFs[i]))
     fastqPairedFilter(c(fnFs[i], fnRs[i]), c(filtFs[i], filtRs[i]),
                       matchIDs=TRUE)
    }
  }
  
  
  #Processing the set of files containing the forward primer in the R1 reads (the sense reads):
  filt_path <- file.path(main_path, "DADA2_SS/filtered/matched") 
  fns <- list.files(filt_path)
  fastqs <- fns[grepl(".fastq.gz", fns)]
  fastqs <- sort(fastqs)
  fnFs <- fastqs[grepl("_F_", fastqs)]
  fnRs <- fastqs[grepl("_R_", fastqs)]
  SSsample.names <- sapply(strsplit(fnFs, "_F_"), `[`, 1)
  filtFs <- file.path(filt_path, fnFs)
  filtRs <- file.path(filt_path, fnRs)
  errF <- learnErrors(filtFs, multithread=TRUE)
  errR <- learnErrors(filtRs, multithread=TRUE)
  derepFs <- derepFastq(filtFs, verbose=TRUE)
  derepRs <- derepFastq(filtRs, verbose=TRUE)
  names(derepFs) <- SSsample.names
  names(derepRs) <- SSsample.names
  SSdadaFs <- dada(derepFs, err=errF, multithread = TRUE)
  SSdadaRs <- dada(derepRs, err=errR, multithread = TRUE)
  SSmergers <- mergePairs(SSdadaFs, derepFs, SSdadaRs, derepRs, verbose=TRUE,minOverlap = 5)
  seqtab_SS <- makeSequenceTable(SSmergers[names(SSmergers)])
  seqtab.nochim_SS <- removeBimeraDenovo(seqtab_SS, verbose=TRUE)
  stSS <- file.path(main_path,"seqtab_SS_RDS")
  stnsSS <- file.path(main_path,"seqtab.nochim_SS_RDS")
  saveRDS(seqtab_SS,stSS)
  saveRDS(seqtab.nochim_SS,stnsSS)
  
  #Then DADA2 processing of "the antisense" reads:
  filt_path <- file.path(main_path, "DADA2_AS/filtered/matched") 
  fns <- list.files(filt_path)
  fastqs <- fns[grepl(".fastq.gz", fns)]
  fastqs <- sort(fastqs)
  fnFs <- fastqs[grepl("_F_", fastqs)]
  fnRs <- fastqs[grepl("_R_", fastqs)]
  ASsample.names <- sapply(strsplit(fnFs, "_F_"), `[`, 1)
  filtFs <- file.path(filt_path, fnFs)
  filtRs <- file.path(filt_path, fnRs)
  errF <- learnErrors(filtFs, multithread=TRUE)
  errR <- learnErrors(filtRs, multithread=TRUE)
  derepFs <- derepFastq(filtFs, verbose=TRUE)
  derepRs <- derepFastq(filtRs, verbose=TRUE)
  names(derepFs) <- ASsample.names
  names(derepRs) <- ASsample.names
  ASdadaFs <- dada(derepFs, err=errF, multithread = TRUE)
  ASdadaRs <- dada(derepRs, err=errR, multithread = TRUE)
  ASmergers <- mergePairs(ASdadaFs, derepFs, ASdadaRs, derepRs, verbose=TRUE,minOverlap = 5)
  seqtab_AS <- makeSequenceTable(ASmergers[names(ASmergers)])
  seqtab.nochim_AS <- removeBimeraDenovo(seqtab_AS, verbose=TRUE)
  stAS <- file.path(main_path,"seqtab_AS_RDS")
  stnsAS <- file.path(main_path,"seqtab.nochim_AS_RDS")
  saveRDS(seqtab_AS,stAS)
  saveRDS(seqtab.nochim_AS,stnsAS)
  
  #Define a function for combining two or more tables, collapsing samples with similar names:  
  sumSequenceTables <- function(table1, table2, ..., orderBy = "abundance") {
   # Combine passed tables into a list
   tables <- list(table1, table2)
   tables <- c(tables, list(...))
   # Validate tables
   if(!(all(sapply(tables, dada2:::is.sequence.table)))) {
    stop("At least two valid sequence tables, and no invalid objects, are expected.")
   }
   sample.names <- rownames(tables[[1]])
   for(i in seq(2, length(tables))) {
    sample.names <- c(sample.names, rownames(tables[[i]]))
   }
   seqs <- unique(c(sapply(tables, colnames), recursive=TRUE))
   sams <- unique(sample.names)
   # Make merged table
   rval <- matrix(0L, nrow=length(sams), ncol=length(seqs))
   rownames(rval) <- sams
   colnames(rval) <- seqs
   for(tab in tables) {
    rval[rownames(tab), colnames(tab)] <- rval[rownames(tab), colnames(tab)] + tab
   }
   # Order columns
   if(!is.null(orderBy)) {
    if(orderBy == "abundance") {
     rval <- rval[,order(colSums(rval), decreasing=TRUE),drop=FALSE]
    } else if(orderBy == "nsamples") {
     rval <- rval[,order(colSums(rval>0), decreasing=TRUE),drop=FALSE]
    }
   }
   rval
  }
  
  stAS <- file.path(main_path,"seqtab_AS_RDS")
  stnsAS <- file.path(main_path,"seqtab.nochim_AS_RDS")
  stSS <- file.path(main_path,"seqtab_SS_RDS")
  stnsSS <- file.path(main_path,"seqtab.nochim_SS_RDS")
  seqtab.nochim_AS <- readRDS(stnsAS)
  seqtab.nochim_SS <- readRDS(stnsSS)
  seqtab_AS <- readRDS(stAS)
  seqtab_SS <- readRDS(stSS)
  total_sumtable <- mergeSequenceTables(seqtab_SS,seqtab_AS, repeats = "sum")
  total_nochim_sumtable <- mergeSequenceTables(seqtab.nochim_SS,seqtab.nochim_AS, repeats = "sum")
  stBoth <- file.path(main_path,"seqtab_Both_RDS")
  stnsBoth <- file.path(main_path,"seqtab.nochim_Both_RDS")
  saveRDS(total_sumtable,stBoth)
  saveRDS(total_nochim_sumtable,stnsBoth)
 }
 
 
 
}

