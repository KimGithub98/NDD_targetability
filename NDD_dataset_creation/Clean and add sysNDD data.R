# ------------------------------------------------------------------------------
# Title: clean sysNDD database files
# Author: KN Wijnant
# Date: 09-06-2022
# Purpose: To deduplicate and prepare for integration with other NDD databases 
# Inputs: sysNDD.xlsx, mimGenes, mimDiseaseTitles, mimDiseaseGenes files
# Outputs: EPI/MR_release_OMIM+ cleaned files
# Dependencies: romim, stringr, dplyr, readxl
# Notes: R version 4.4.1 (2024-06-14 ucrt)
# ------------------------------------------------------------------------------


library(romim)
library(stringr)
library(dplyr)
library(readxl)

sysNDD <- read_excel("sysNDD.xlsx", )
sysNDD <- sysNDD[,c(4,5,6,9,13)]
colnames(sysNDD) <- c("HGNC_symbol", "OMIM_Disease", "Disease_name", "Inheritance", "Category")

#get OMIM files
mimGenes <- read.delim("mimGenes.txt", header=TRUE, comment.char="#")
mimGenes <- mimGenes[,c(1,4)]
colnames(mimGenes) <- c("OMIM_GeneID", "HGNC_symbol")
mimDiseaseTitles <- read.delim("mimDiseaseTitles.txt", header=TRUE, comment.char="#")
mimDiseaseTitles <- mimDiseaseTitles[,c(2,3)]
colnames(mimDiseaseTitles) <- c("OMIM_Disease", "Disease_name")
mimDiseaseTitles$OMIM_Disease <- as.character(mimDiseaseTitles$OMIM_Disease)
mimDiseaseGenes <- read.delim("mimDiseaseGenes.txt", comment.char="#")

translate_inheritance <- function(inheritance_words){
  inheritance_words = inheritance_words[1]
  if (inheritance_words == "Autosomal dominant"){
    inheritance_abbr = "AD"
  }
  else if (inheritance_words == "Autosomal recessive"){
    inheritance_abbr = "AR" 
  }
  else if (inheritance_words == "X-linked"){
    inheritance_abbr = "XL"
  }
  else if (inheritance_words == "Other"){
    inheritance_abbr = "NA"
  }
  else {
    inheritance_abbr = "NA"
  }
}

sysNDD <- left_join(sysNDD, mimGenes, by = "HGNC_symbol")

#extract only OMIM numbers and remove "OMIM:"
sysNDD[!grepl("OMIM", sysNDD$OMIM_Disease, fixed = TRUE), "OMIM_Disease"] <- rep(NA, sum(!grepl("OMIM", sysNDD$OMIM_Disease, fixed = TRUE)))
sysNDD$OMIM_Disease <- sapply(sysNDD$OMIM_Disease, function(x) str_extract(x, "[0-9]{6}"))
sysNDD$Inheritance <- sapply(sysNDD$Inheritance, translate_inheritance)

write.csv(sysNDD, "cleaned_sysNDD.csv", row.names = TRUE)
