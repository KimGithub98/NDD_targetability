# ------------------------------------------------------------------------------
# Title: integrate cleaned data from HPO, sysNDD, DDG2P and the EPI & MR RadboudUMC gene panels
# Author: KN Wijnant
# Date: 14-11-2023
# Purpose: To deduplicate and prepare for integration with other NDD databases 
# Inputs: cleaned and manually assessed DDG2P files, cleaned and manually assessed MR & EPI gene panel release files, cleaned sysNDD file, cleaned HPO file
# Outputs: NDD_dataset, file for manual assessment
# Dependencies: readxl, plyr, dplyr, stringr, romim, data.table, tidyverse
# Notes: R version 4.4.1 (2024-06-14 ucrt)
# ------------------------------------------------------------------------------

library(readxl)
library(plyr)
library(dplyr)
library(stringr)
library(romim)
library(data.table)
library(tidyverse)
outersect <- function(x, y) {
  sort(setdiff(x, y))
}

options(scipen = 999)

#fetch datasets
DDG2P_ID_genelist <- read.csv("DDG2P_ID_genelist+_manualadded_+wrongOMIM.csv", row.names=1, sep=",")
DDG2P_ID_genelist$OMIM_GeneID <- as.character(DDG2P_ID_genelist$OMIM_GeneID)
DDG2P_ID_genelist$OMIM_Disease <- as.numeric(DDG2P_ID_genelist$OMIM_Disease)
MR_release <- read.csv("MR_release_OMIM+ (manual adaptations).csv", row.names=1, sep=";")
EPI_release <- read.csv("EPI_release_OMIM+ (manual adaptations).csv", row.names = 1, sep = ";")
colnames(MR_release)[5] <- "OMIM_GeneID"
MR_release$OMIM_GeneID <- as.character(MR_release$OMIM_GeneID)
MR_release$OMIM_Disease <- as.numeric(MR_release$OMIM_Disease)
colnames(EPI_release)[5] <- "OMIM_GeneID"
EPI_release$OMIM_GeneID <- as.character(EPI_release$OMIM_GeneID)
EPI_release$OMIM_Disease <- as.numeric(EPI_release$OMIM_Disease)
data_HP_0001249 <-read.csv("data_HP_0001249.csv", row.names=1)
data_HP_0001250 <-read.csv("data_HP_0001250.csv", row.names=1)
data_HP_0012758 <-read.csv("data_HP_0012758.csv", row.names=1)
sysNDD <- read.csv("cleaned_sysNDD.csv", row.names=1)

#prepare all dfs for merge
DDG2P_ID_genelist <- data.frame(DDG2P_ID_genelist$HGNC_symbol, DDG2P_ID_genelist$OMIM_GeneID, DDG2P_ID_genelist$OMIM_Disease, DDG2P_ID_genelist$Disease_name , DDG2P_ID_genelist$Confidence, DDG2P_ID_genelist$Inheritance, rep("DDG2P_BrainCog",nrow(DDG2P_ID_genelist)))
colnames(DDG2P_ID_genelist) <- c("HGNC_symbol", "OMIM_GeneID", "OMIM_Disease", "Disease_name", "Confidence_DDG2P", "Inheritance", "sources_DDG2P")

colnames(sysNDD)[5] <- "Confidence_sysNDD"
sysNDD[,"Source_sysNDD"] <- rep("sysNDD",nrow(sysNDD))

EPI_release["Source_EPI"] <- "EPIpanelRadboudumc"
MR_release["Source_MR"] <- "IDpanelRadboudumc"

data_HP_0001249 <- data_HP_0001249[,-1]
data_HP_0001250 <- data_HP_0001250[,-1]
data_HP_0012758 <- data_HP_0012758[,-1]
data_HP_0001249["Source_HP_0001249"] <- "HPOmim_ID"
data_HP_0001250["Source_HP_0001250"] <- "HPOmim_EPI"
data_HP_0012758["Source_HP_0012758"] <- "HPOmim_NDD"

#make combined df
df_list <- list(DDG2P_ID_genelist, sysNDD, EPI_release, MR_release, data_HP_0001249, data_HP_0001250, data_HP_0012758)
Gen_dis_treat_df <- df_list %>% reduce(full_join, by = c('OMIM_Disease', "HGNC_symbol", "Inheritance"))

#combine disease_names, gene_omim and sources
OMIM_GeneID_cols <- colnames(Gen_dis_treat_df)[grepl("OMIM_GeneID", colnames(Gen_dis_treat_df))]
Disease_name_cols <- colnames(Gen_dis_treat_df)[grepl("Disease_name", colnames(Gen_dis_treat_df))]
Source_cols <- colnames(Gen_dis_treat_df)[grepl("ource", colnames(Gen_dis_treat_df))]
for (rowi in rownames(Gen_dis_treat_df)){
  row <- Gen_dis_treat_df[rowi,]
  GeneIDs <- c(row[OMIM_GeneID_cols])
  GeneID <- unique(as.character(GeneIDs))[!(unique(as.character(GeneIDs)) == "NA" | is.na(unique(as.character(GeneIDs))))]
  if (length(GeneID)>1){
    print("multiple GeneIDs") 
    break}
  else if (length(GeneID)==0){
    GeneID = NA
  }
  Disease_names <- unlist(c(row[Disease_name_cols]))
  Disease_name <- names(sort(table(Disease_names), decreasing = TRUE)[1])
  if (is.null(Disease_name)){
    Disease_name = NA
  }
  Sources <- unlist(c(row[Source_cols]))
  Sources <- paste(Sources[!is.na(Sources)], collapse = ";")
  comments <- unique(row$comments_PMID.x,row$comments_PMID.y)[!(unique(row$comments_PMID.x,row$comments_PMID.y) == "-" | is.na(unique(row$comments_PMID.x,row$comments_PMID.y)))]
  if (length(comments) == 0){
    comments = NA
  }
  Gen_dis_treat_df[rowi, c("OMIM_GeneID_combined", "Disease_name_combined", "Sources", "comments_PMID")] <- c(GeneID, Disease_name, Sources, comments)
}
drop <- c(OMIM_GeneID_cols, Disease_name_cols, Source_cols, "comments_PMID.y", "comments_PMID.x")
Gen_dis_treat_df <- Gen_dis_treat_df[,!(names(Gen_dis_treat_df) %in% drop)]
colnames(Gen_dis_treat_df)[8] <- "Disease_name"
colnames(Gen_dis_treat_df)[7] <- "OMIM_GeneID_combined"

#remove duplicated entries (execute 2x to delete all duplications)
Gen_dis_treat_df_dup <- Gen_dis_treat_df[Gen_dis_treat_df$HGNC_symbol %in% unique(Gen_dis_treat_df$HGNC_symbol[duplicated(Gen_dis_treat_df$HGNC_symbol)]) & Gen_dis_treat_df$OMIM_Disease %in% unique(Gen_dis_treat_df$OMIM_Disease[duplicated(Gen_dis_treat_df$OMIM_Disease)]),]
num_remove <- c()
manual_asses <- data.frame()
i = 0
for (gene in unique(Gen_dis_treat_df_dup$HGNC_symbol)){
  duplicates <- Gen_dis_treat_df[Gen_dis_treat_df$HGNC_symbol == gene,]
  if (nrow(duplicates) < 2){next}
  for (OMIM_Disease in unique(duplicates$OMIM_Disease)[!is.na(unique(duplicates$OMIM_Disease))]){
    duplicates2 <- duplicates[duplicates$OMIM_Disease == OMIM_Disease & !is.na(duplicates$OMIM_Disease),]
    if (nrow(duplicates2) <2){next}
    Confidence_DDG2P <- paste(unique(duplicates2$Confidence_DDG2P)[!is.na(unique(duplicates2$Confidence_DDG2P))], collapse = ";")
    Confidence_sysNDD <- paste(unique(duplicates2$Confidence_sysNDD)[!is.na(unique(duplicates2$Confidence_sysNDD))], collapse = ";")
    Linked_OMIM_Disease <- paste(unique(duplicates2$Linked_OMIM_Disease)[!is.na(unique(duplicates2$Linked_OMIM_Disease))], collapse = ";")
    combined_inheritance <- c()
    nsource <- lapply(str_split(duplicates2$Sources, ";"), length)
    dup <- duplicated(duplicates2$Inheritance)
    for (rowd in 1:nrow(duplicates2)){
      nsou <- nsource[rowd]
      if (dup[rowd]){
        first_dup <- which(duplicates2$Inheritance == duplicates2$Inheritance[rowd])[1]
        new_nsource <- length(unique(unlist(str_split(duplicates2$Sources[rowd], ";")), unlist(str_split(duplicates2$Sources[first_dup], ";"))))
        combined_inheritance[first_dup] <- gsub("[0-9]", new_nsource, combined_inheritance[first_dup])
      }
      else {
        inh <- duplicates2$Inheritance[rowd]
        combined_inheritance <- c(combined_inheritance, paste(inh, "(", nsou, ")", sep = ""))
      }
    }
    combined_inheritance <- paste(combined_inheritance, collapse = ";", sep = "")
    diseasename <- duplicates2$Disease_name[order(unlist(nsource), decreasing = T)[1]]
    Sources <- paste(unique(unlist(str_split(duplicates2$Sources, ";"))), collapse = ";", sep = "")
    Gen_dis_treat_df[rownames(duplicates2)[1],"Confidence_DDG2P"] <- Confidence_DDG2P
    Gen_dis_treat_df[rownames(duplicates2)[1],"Inheritance"] <- combined_inheritance
    Gen_dis_treat_df[rownames(duplicates2)[1],"Confidence_sysNDD"] <- Confidence_sysNDD
    Gen_dis_treat_df[rownames(duplicates2)[1],"Linked_OMIM_Disease"] <- Linked_OMIM_Disease
    Gen_dis_treat_df[rownames(duplicates2)[1],"Disease_name"] <- diseasename
    Gen_dis_treat_df[rownames(duplicates2)[1],"Sources"] <- Sources
    num_remove <- c(num_remove, rownames(duplicates2)[-1])
  }
}
Gen_dis_treat_df <- Gen_dis_treat_df[-c(as.numeric(num_remove)),]
#write.csv(Gen_dis_treat_df, file = "NDD_dataset_removed_dup.csv", row.names = TRUE)

#remove duplicate rows with unknown OMIM_disease
Gen_dis_treat_df_dup <- Gen_dis_treat_df[Gen_dis_treat_df$HGNC_symbol %in% unique(Gen_dis_treat_df$HGNC_symbol[duplicated(Gen_dis_treat_df$HGNC_symbol)]) & Gen_dis_treat_df$HGNC_symbol %in% Gen_dis_treat_df$HGNC_symbol[is.na(Gen_dis_treat_df$OMIM_Disease)],]

#create file for for manual assessment
manual_asses <- data.frame()
for (gene in unique(Gen_dis_treat_df_dup$HGNC_symbol)){
  duplicates <- Gen_dis_treat_df[Gen_dis_treat_df$HGNC_symbol == gene,]
  if (nrow(duplicates) < 2){next}
  for (rowi in rownames(duplicates)[is.na(duplicates$OMIM_Disease)]){
    row = duplicates[rowi,]
    duplicates_min = duplicates[rownames(duplicates) != rowi,]
    INH = row$Inheritance
    if (is.na(INH)){next}
    SOURCE = row$Sources
    for (row_min in rownames(duplicates_min)){
      row_min = duplicates_min[row_min,]
      assess <- grepl(INH, row_min$Inheritance) & !grepl(SOURCE, row_min$Sources)
      if (assess){
        manual_asses<- rbind(manual_asses, row_min)
        manual_asses<- rbind(manual_asses, row)
      }
    }
  }
}
#write.csv(manual_asses, file = "manual_assess_duplications.csv", row.names = TRUE)


