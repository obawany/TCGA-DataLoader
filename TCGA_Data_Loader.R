# Preamble
################################################################################
# Title : TCGA Data Loader
# INPUT : TCGA website datasets of your choosing
# OUTPUT: TCGA datasets
# Date  : 20-09-2017 
# TCGA_Data_Loader.R finds datasets of interest from TCGA, downloads them into R
# and saves them into a big dataframe.
################################################################################

## Load libraries---------------------------------------------------------------
library(TCGAbiolinks)
library(SummarizedExperiment)
library(devtools)  # <- This needs to be installed 
devtools::install_github("BioinformaticsFMRP/TCGAbiolinks") # <- patches for api

## Directory--------------------------------------------------------------------
setwd("/Users/obawany/Desktop/Winter 2018/Honours Research")

i=17

# Projects
project_names <- c("TCGA-LIHC", "TCGA-COAD", "TCGA-PAAD", "TARGET-WT", 
                  "TCGA-UCS",  "TCGA-UVM",  "TCGA-STAD", "TCGA-LUSC", "TCGA-UCEC", 
                   "TCGA-CHOL", "TCGA-OV",   "TCGA-ACC",  "TCGA-MESO", "TARGET-AML", 
                   "TCGA-LGG",  "TCGA-ESCA", 
                   "TCGA-SARC", "TCGA-THCA", "TCGA-READ", "TCGA-SKCM",
                   "TCGA-GBM",  "TCGA-KIRC", "TCGA-BRCA", "TCGA-PCPG", "TCGA-CESC", 
                   "TCGA-LUAD", "TCGA-BLCA", "TCGA-KIRP", "TCGA-PRAD", "TCGA-KICH", 
                   "TCGA-TGCT", "TCGA-LAML", "TCGA-THYM", "TCGA-DLBC")



## RNA-Seq Gene expression------------------------------------------------------
# Downloading GDC data
query_RNA_Seq <- GDCquery(project = project_names[i],
                          data.category = "Transcriptome Profiling",
                          data.type = "Gene Expression Quantification", 
                          workflow.type = "HTSeq - Counts")

# Preparing GDC data
GDCdownload(query_RNA_Seq)

# Saving GDC data into .Rda
exp_filename = paste(project_names[i],"RNA_HTSeq_Counts.rda",sep="_")
data_RNA_Seq <- GDCprepare(query_RNA_Seq, save = TRUE, save.filename = exp_filename)

# Saving .Rda into .txt
exp_filename = paste(project_names[i],"RNA_HTSeq_Counts.txt",sep="_")
test <- assay(data_RNA_Seq)
write.csv(test, file=exp_filename)

# Clear memory-intensive variables
rm(test, query_RNA_Seq, data_RNA_Seq, exp_filename)

## miRNA-Seq Gene expression----------------------------------------------------
# Downloading GDC data
query_miRNA_Seq <- GDCquery(project = project_names[i],
                          data.category = "Transcriptome Profiling",
                          data.type = "miRNA Expression Quantification")

# Preparing GDC data
GDCdownload(query_miRNA_Seq)

# Saving GDC data into .Rda
exp_filename = paste(project_names[i],"miRNA_Quantification.rda",sep="_")
data_miRNA_Seq <- GDCprepare(query_miRNA_Seq, save = TRUE, save.filename = exp_filename)

# Saving .Rda into .txt
exp_filename = paste(project_names[i],"miRNA_Expression_Quantification.txt",sep="_")
test <- as.data.frame(data_miRNA_Seq)
write.csv(test, file=exp_filename)

# Clear memory-intensive variables
rm(test, query_miRNA_Seq, data_miRNA_Seq, exp_filename)

