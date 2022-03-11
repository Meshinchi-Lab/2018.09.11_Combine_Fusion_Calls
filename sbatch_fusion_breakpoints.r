#! /usr/bin/env Rscript

library(methods)
library(stringr)
library(magrittr)

library(dplyr)
library(tidyr)
library(tibble)
library(readr)
suppressPackageStartupMessages(library(Biostrings))

# setwd(file.path(SCRATCH,"jlsmith3/Fusion_Breakpoints/"))
setwd(file.path(PROJHOME,"2018.09.11_Combine_Fusion_Calls"))
source(file.path(SCRIPTS,"RNAseq_Analysis/Analysis/Fusions_Analysis/Breakpoint_Analysis_Function.r"))


#define files 
inputs <- dir(getwd(), pattern="^[A-G]_.+\\.RDS$")
print(inputs)

if(length(inputs) == 0){
  print("Missing Inputs. Check regular expression.")
  quit(save="no")
}else if(length(inputs) > 7){
  print("Too Many Inputs. Check regular expression.")
  quit(save="no")
}

#Read in the input files 
in_files <- lapply(inputs, readRDS)
names(in_files) <- inputs

#Define the Sample
Sample <- Sys.getenv("SAMPLE_ID")

#Run the workflow
if(Sample == ""){
  print("Sample_ID not found. Retry")
  quit(save="no")
}else{
  print(Sample)
  #run the function 
  breakpoints <- examine_breakpoint_evidence(Sample_ID=Sample,
                                             fusion_data=in_files[[1]],
                                             file_manifest=in_files[[2]],
                                             mutseq=in_files[[3]],
                                             fusion_theoretical=in_files[[4]],
                                             Seqlens=in_files[[5]],
                                             cicero.contig=in_files[[6]],
                                             TA.contigs=in_files[[7]])
  
  try(traceback(), silent=T) #for printing the error if one occurs
  saveRDS(breakpoints, file.path(paste0(Sample,"_breakpoint_sequence_evidence.RDS")))
}

