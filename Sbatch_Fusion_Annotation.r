#!/app/easybuild/software/R/3.5.3-foss-2016b/bin/Rscript
#Jenny Smith 
#10/8/18
#Purpose: to split the cleaned and combined fusion data into parts for annotation of all fusions. Its too memory intenstive to be done without paralellization. 



#print R version to ensure correct environment was used. 
print(version)

library(magrittr)
library(dplyr)
library(tibble)
library(tidyr)
library(methods)
source("~/scripts/RNAseq_Analysis/Fusions_Analysis/Annotate_Fusions_Function.r")

#Set the working directory 
setwd('/fh/fast/meshinchi_s/workingDir/TARGET/AML_TARGET/RNA/mRNAseq/analysis/2018.09.11_Combine_Fusion_Calls/')
options(stringsAsFactors = FALSE)


#Read in the files 
combined.clean.annot <- read.csv("TARGET_AML_0531_1031_TransAbyss_STAR_TargetedAlignment_Combined_AnnotatedBlacklist.csv")
dim(combined.clean.annot)

dbs <-  readRDS("Fusion_Reference/Fusion_ReferenceFiles_10.8.18.RDS")

gene_aliases <- read.csv("~/RNA_seq_Analysis/0000.00.02_Reference_GeneInfo/Hs_GeneInfo.Symb.Alias.csv",
                         header = TRUE, stringsAsFactors = FALSE)
gene_aliases$alias <- toupper(gene_aliases$alias)
head(gene_aliases)


#RANGE is a variable that will be exported during the sbatch command 
#The exact code from command line: 
# > ranges=( "1:8434" "8435:16867"  "16868:25299" "25300:33732" "33733:42165" "42166:50598" "50599:59030" "59031:67463" "67464:75896")
# > for i in $(echo ${ranges[*]}) ; do   sbatch --mail-type=FAIL,END --mail-user=jlsmith3@fredhutch.org --export=RANGE=$i ~/scripts/RNAseq_Analysis/Fusions_Analysis/Sbatch_Fusion_Annotation.r ; done
rows.range <- Sys.getenv('RANGE') 
forEval <- paste("c(", rows.range, ")")
range <- eval(parse(text= forEval))
print(rows.range)


#Define a filename and the path
dir.create(path = "/fh/scratch/delete90/meshinchi_s/jlsmith3/Annotate_Fusions/", recursive = TRUE)
path="/fh/scratch/delete90/meshinchi_s/jlsmith3/Annotate_Fusions/"
filename <- paste0("TARGET_AML_0531_1031_TransAbyss_STAR_TargetedAlignment_Combined_Annotated_",rows.range,"_.csv")

print(dir(path = path))
print("Starting Annotation")

#Annotate the subset of the whole fusions file. 
temp <- combined.clean.annot %>%
  #select a range of rows. 
  slice(range) %>%
  
  #Annotate Breakpoints in other Databases. 
  rowwise() %>%
  mutate(Present_inFusionCancer_byFusionName=AnnotFusions(c1 = All_Fusions_Called,
                                                          df2.annot = dbs$db_fusionCancer,
                                                          alias.df = gene_aliases)) %>%
  
  mutate(Present_inMitelman_byFusionName=AnnotFusions(c1 = All_Fusions_Called,
                                                      df2.annot = dbs$db_mitelman,
                                                      alias.df = gene_aliases)) %>%
  
  mutate(Present_inTicdb_byFusionName=AnnotFusions(c1 = All_Fusions_Called,
                                                   df2.annot = dbs$db_ticdb,
                                                   alias.df = gene_aliases)) %>%
  
  mutate(Present_inCOSMIC_byFusionName=AnnotFusions(c1 = All_Fusions_Called,
                                                    df2.annot = dbs$db_cosmic,
                                                    alias.df = gene_aliases)) %>%
  
  mutate(Present_inTumorFusion_byFusionName=AnnotFusions(c1 = All_Fusions_Called,
                                                         df2.annot = dbs$db_tumorFusion,
                                                         alias.df = gene_aliases)) %>%
  
  #Other Annotation Columns
  mutate(CancerType_inFusionCancer_byFusionName=AnnotFusions(c1 = All_Fusions_Called,
                                                             df2.annot = dbs$db_fusionCancer,
                                                             alias.df = gene_aliases,
                                                             col2 = "Cancer_type")) %>%
  
  mutate(Morphology_inMitelman_byFusionName=AnnotFusions(c1 = All_Fusions_Called,
                                                         df2.annot = dbs$db_mitelman,
                                                         alias.df = gene_aliases,
                                                         col2="Morphology")) %>%
  
  mutate(FAB_inMitelman_byFusionName=AnnotFusions(c1 = All_Fusions_Called,
                                                  df2.annot = dbs$db_mitelman,
                                                  alias.df = gene_aliases,
                                                  col2="FAB_Type")) %>%
  
  mutate(CancerType_inCOSMIC_byFusionName=AnnotFusions(c1 = All_Fusions_Called,
                                                       df2.annot = dbs$db_cosmic,
                                                       alias.df = gene_aliases,
                                                       col2="Primary.histology")) %>%
  
  mutate(CancerType_inTumorFusion_byFusionName=AnnotFusions(c1 = All_Fusions_Called,
                                                            df2.annot = dbs$db_tumorFusion,
                                                            alias.df = gene_aliases,
                                                            col2="Cancer")) %>%
  ungroup()




dim(temp) 


#write the results to a file 
write.csv(temp,
          file = paste0(path,filename),
          row.names = FALSE )











