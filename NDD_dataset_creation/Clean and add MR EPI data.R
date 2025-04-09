# ------------------------------------------------------------------------------
# Title: clean MR EPI gene panels RadboudUMC
# Author: KN Wijnant
# Date: 26-5-2022
# Purpose: To deduplicate and prepare for integration with other NDD databases 
# Inputs: MR/EPI_release RadboudUMC gene panel files, mimGenes, mimDiseaseTitles, mimDiseaseGenes files
# Outputs: EPI/MR_release_OMIM+ cleaned files
# Dependencies: romim, stringr, dplyr, plyr, readxl
# Notes: R version 4.4.1 (2024-06-14 ucrt)
# ------------------------------------------------------------------------------


library(readxl)
library(plyr)
library(dplyr)
library(stringr)
library(romim)
outersect <- function(x, y) {
  sort(setdiff(x, y))
}

get_inheritance <- function(my_xml){
  my_inheritance_node <- getNodeSet(my_xml, path = "/omim/entryList/entry/phenotypeMapList/phenotypeMap/phenotypeInheritance")
  if (length(my_inheritance_node) == 0){
    my_inheritance_node <- getNodeSet(my_xml, path = "/omim/entryList/entry/geneMap/phenotypeMapList/phenotypeMap/phenotypeInheritance")
  }
  xmlSApply(my_inheritance_node, xmlValue)
}

get_PMID <- function(my_xml){
  my_inheritance_node <- getNodeSet(my_xml, path = "/omim/entryList/entry/referenceList/reference/pubmedID")
  xmlSApply(my_inheritance_node, xmlValue)
}

#import and change colnames
MR_release <- read.delim2("MR_release.txt")
EPI_release <- read.delim2("EPI_release.txt")
MR_release <- MR_release[,c(2,3,4,5)]
EPI_release <- EPI_release[,c(2,3,4,5)]
colnames(MR_release) <- c("HGNC_symbol", "OMIM_Disease", "Inheritance", "comments_PMID")
colnames(EPI_release) <- c("HGNC_symbol", "OMIM_Disease", "Inheritance", "comments_PMID")

#get OMIM files
mimGenes <- read.delim("mimGenes.txt", header=TRUE, comment.char="#")
mimGenes <- mimGenes[,c(1,4)]
colnames(mimGenes) <- c("OMIM_gene", "HGNC_symbol")
mimDiseaseTitles <- read.delim("mimDiseaseTitles.txt", header=TRUE, comment.char="#")
mimDiseaseTitles <- mimDiseaseTitles[,c(2,3)]
colnames(mimDiseaseTitles) <- c("OMIM_Disease", "Disease_name")
mimDiseaseTitles$OMIM_Disease <- as.character(mimDiseaseTitles$OMIM_Disease)
mimDiseaseGenes <- read.delim("mimDiseaseGenes.txt", comment.char="#")

#make row for each disease (MR_release, EPI_release) and add OMIM data
for (release in list(MR_release, EPI_release)){
  release <- left_join(release, mimGenes, by = "HGNC_symbol") #add OMIM_Genes
  if (release[1,1] == MR_release[1,1]){
    dataf = "MR_release"
    MR_release = release
    }
  if (release[1,1] == EPI_release[1,1]){
    dataf = "EPI_release"
    EPI_release = release
  }
  new_rows <- NULL
  for (rowi in rownames(release)){ #look if double and if that is the case split into seperate rows
    Disease_ID <- release[rowi,"OMIM_Disease"]
    if (Disease_ID == "-" | Disease_ID == ""){next}
    if (grepl(";", Disease_ID, fixed = TRUE)){
      row <- release[rowi,]
      if (dataf == "MR_release"){MR_release <- MR_release[-c(which(MR_release$OMIM_Disease == Disease_ID)),]} #remove row
      if (dataf == "EPI_release"){EPI_release <- EPI_release[-c(which(EPI_release$OMIM_Disease == Disease_ID)),]} #remove row
      unique_DiD <- strsplit(as.character(row$OMIM_Disease), split=";")
      for (DiD in unique_DiD[[1]]){
        new_row <- data.frame(as.character(row$HGNC_symbol), DiD, as.character(row$Inheritance), NA, as.character(row$OMIM_gene)) 
        colnames(new_row) <- colnames(release)
        new_rows <- rbind(new_rows, new_row)
      }
      colnames(new_rows) <- colnames(release)
    }
  }
  if (dataf == "MR_release"){MR_release <- rbind(MR_release, new_rows)} #create new rows
  if (dataf == "EPI_release"){EPI_release <- rbind(EPI_release, new_rows)} #create new rows
  if (release[1,1] == MR_release[1,1]){ #add mimDisease Titles by OMIM_disease
    MR_release <- left_join(MR_release, mimDiseaseTitles, by = "OMIM_Disease")
  }
  if (release[1,1] == EPI_release[1,1]){
    EPI_release <- left_join(EPI_release, mimDiseaseTitles, by = "OMIM_Disease")
  }
}

#write.csv(MR_release, "MR_release_OMIM.csv", row.names = TRUE)
#write.csv(EPI_release, "EPI_release_OMIM.csv", row.names = TRUE)

#function rewrite inheritance
translate_inheritance <- function(inheritance_words){
  inheritance_words = inheritance_words[1]
  if (inheritance_words == "Autosomal dominant"){
    inheritance_abbr = "AD"
  }
  else if (inheritance_words == "Autosomal recessive"){
    inheritance_abbr = "AR" 
  }
  else if (inheritance_words == "X-linked recessive"){
    inheritance_abbr = "XLR"
  }
  else if (inheritance_words == "X-linked"){
    inheritance_abbr = "XL"
  }
  else if (inheritance_words == "Autosomal dominant; Somatic mosaicism"){
    inheritance_abbr = "AD"
  }
  else if (inheritance_words == "X-linked dominant"){
    inheritance_abbr = "XLD"
  }
  else if (inheritance_words == "X-linked dominant; X-linked recessive"){
    inheritance_abbr = "XLD,XLR"
  }
  else {
    inheritance_abbr = NA
  }
}

#add OMIM data for MR
#add missing omim_disease ID data to MR and EPI
set_key("") #request and insert OMIM key
#PAY attention that you can only send 5000 requests a day (time is set at 7am NY time)
for (row in rownames(MR_release)[1:nrow(MR_release)]){
  skip_to_next = FALSE
  if (MR_release[row, "OMIM_Disease"] == "-" | MR_release[row, "OMIM_Disease"] == ""){
    gene <- MR_release[row, "HGNC_symbol"]
    my_list <- tryCatch(gene_to_omim(gene, show_query = FALSE), error = function(e) { skip_to_next <<- TRUE})
    if (skip_to_next == TRUE | is.null(my_list)){next}
    if (grepl("\\d{8}", MR_release[row, "comments_PMID"])){
      my_list_omim <- sapply(my_list, get_omim, geneMap = TRUE, referenceList = TRUE)
      my_list_titles <- sapply(my_list_omim, get_title) #get disease names
      diseaseID = names(my_list_titles) #get OMIM_ID
      PMID <- str_extract(MR_release[row, "comments_PMID"], "\\d{8}")
      my_list_PMID <- sapply(my_list_omim, get_PMID)
      if (length(names(my_list_PMID)) > 1){
        print(gene)
        which_list <- sapply(my_list_PMID, intersect, PMID)
        DiseaseID <- unique(names(which_list)[sapply(which_list, length) > 0])
        Inheritance <- sapply(my_list_omim, get_inheritance)
        Inheritance <- sapply(Inheritance, translate_inheritance)
        diseases_include <- which(DiseaseID == my_list & Inheritance == MR_release[row, "Inheritance"])
        MR_release[row, c("Disease_name", "OMIM_Disease")] <- c(paste(my_list_titles[diseases_include], collapse='|' ), paste(diseaseID[diseases_include], collapse="|"))
      }
      else if (PMID %in% my_list_PMID){
        if (all(class(my_list_PMID) == c("matrix", "array"))){
          DiseaseID <- colnames(my_list_PMID)
        }
        else if (class(my_list_PMID) == "character"){
          DiseaseID <- names(my_list_PMID)
        }
        MR_release[row, c("Disease_name", "OMIM_Disease")] <- c(my_list_titles[DiseaseID], DiseaseID)
      }
      else {message(paste(MR_release[row,], collapse='|' ))}
    }
    closeAllConnections()
  }
}

# Update Inheritance data for MR
for (row in rownames(MR_release)[1:nrow(MR_release)]){
  if ("XL" %in% MR_release[row, "Inheritance"] | is.na(MR_release[row, "Inheritance"])){
    if (grepl("\\d", MR_release[row, "OMIM_Disease"])){
      OMIM_Disease <- MR_release[row, "OMIM_Disease"]
      my_list_omim <- get_omim(OMIM_Disease, geneMap = TRUE)
      my_list_inheritance1 <- get_inheritance(my_list_omim) #get inheritance #this function is changed!!!!!
      if (length(my_list_inheritance1) > 1 & length(unique(my_list_inheritance)) == 1){
        my_list_inheritance1 = my_list_inheritance1[1]
      }
      my_list_inheritance <- translate_inheritance(my_list_inheritance1)
      if (!is.na(my_list_inheritance)){
        MR_release[row, "Inheritance"] <- my_list_inheritance
        message(paste(row, my_list_inheritance1, my_list_inheritance))
      }
    }
  }
  else if (as.numeric(row) >= 1416){
    if (grepl("\\d", MR_release[row, "OMIM_Disease"])){
      OMIM_Disease <- MR_release[row, "OMIM_Disease"]
      my_list_omim <- get_omim(OMIM_Disease, geneMap = TRUE)
      my_list_inheritance1 <- get_inheritance(my_list_omim) #get inheritance #this function is changed!!!!!
      if (length(my_list_inheritance1) > 1 & length(unique(my_list_inheritance)) == 1){
        my_list_inheritance1 = my_list_inheritance1[1]
      }
      my_list_inheritance <- translate_inheritance(my_list_inheritance1)
      if (!is.na(my_list_inheritance)){
        MR_release[row, "Inheritance"] <- my_list_inheritance
        message(paste(row, my_list_inheritance1, my_list_inheritance))
      }
    }
  }
}
  
#write.csv(MR_release, "MR_release_OMIM+.csv", row.names = TRUE)

#add OMIM data for MR
#add missing omim_disease ID data to MR and EPI
set_key("") #request and insert OMIM key
#PAY attention that you can only send 5000 requests a day (time is set at 7am NY time)
for (row in rownames(EPI_release)[1:nrow(EPI_release)]){
  skip_to_next = FALSE
  if (EPI_release[row, "OMIM_Disease"] == "-" | EPI_release[row, "OMIM_Disease"] == ""){
    gene <- EPI_release[row, "HGNC_symbol"]
    my_list <- tryCatch(gene_to_omim(gene, show_query = FALSE), error = function(e) { skip_to_next <<- TRUE})
    if (skip_to_next == TRUE | is.null(my_list)){next}
    if (grepl("\\d{8}", EPI_release[row, "comments_PMID"])){
      my_list_omim <- sapply(my_list, get_omim, geneMap = TRUE, referenceList = TRUE)
      my_list_titles <- sapply(my_list_omim, get_title) #get disease names
      diseaseID = names(my_list_titles) #get OMIM_ID
      PMID <- str_extract(EPI_release[row, "comments_PMID"], "\\d{8}")
      my_list_PMID <- sapply(my_list_omim, get_PMID)
      if (length(names(my_list_PMID)) > 1){
        print(gene)
        which_list <- sapply(my_list_PMID, intersect, PMID)
        DiseaseID <- unique(names(which_list)[sapply(which_list, length) > 0])
        Inheritance <- sapply(my_list_omim, get_inheritance)
        Inheritance <- sapply(Inheritance, translate_inheritance)
        diseases_include <- which(DiseaseID == my_list & Inheritance == EPI_release[row, "Inheritance"])
        EPI_release[row, c("Disease_name", "OMIM_Disease")] <- c(paste(my_list_titles[diseases_include], collapse='|' ), paste(diseaseID[diseases_include], collapse="|"))
      }
      else if (PMID %in% my_list_PMID){
        if (all(class(my_list_PMID) == c("matrix", "array"))){
          DiseaseID <- colnames(my_list_PMID)
        }
        else if (class(my_list_PMID) == "character"){
          DiseaseID <- names(my_list_PMID)
        }
        EPI_release[row, c("Disease_name", "OMIM_Disease")] <- c(my_list_titles[DiseaseID], DiseaseID)
      }
      else {message(paste(EPI_release[row,], collapse='|' ))}
    }
    closeAllConnections()
  }
}

# Update Inheritance data for EPI
for (row in rownames(EPI_release)[1:nrow(EPI_release)]){
  if ("XL" %in% EPI_release[row, "Inheritance"] | is.na(EPI_release[row, "Inheritance"])){
    if (grepl("\\d", EPI_release[row, "OMIM_Disease"])){
      OMIM_Disease <- EPI_release[row, "OMIM_Disease"]
      my_list_omim <- get_omim(OMIM_Disease, geneMap = TRUE)
      my_list_inheritance1 <- get_inheritance(my_list_omim) #get inheritance #this function is changed!!!!!
      if (length(my_list_inheritance1) > 1 & length(unique(my_list_inheritance)) == 1){
        my_list_inheritance1 = my_list_inheritance1[1]
      }
      my_list_inheritance <- translate_inheritance(my_list_inheritance1)
      if (!is.na(my_list_inheritance)){
        EPI_release[row, "Inheritance"] <- my_list_inheritance
        message(paste(row, my_list_inheritance1, my_list_inheritance))
      }
    }
  }
  else if (as.numeric(row) >= 369){
    if (grepl("\\d", EPI_release[row, "OMIM_Disease"])){
      OMIM_Disease <- EPI_release[row, "OMIM_Disease"]
      my_list_omim <- get_omim(OMIM_Disease, geneMap = TRUE)
      my_list_inheritance1 <- get_inheritance(my_list_omim) #get inheritance #this function is changed!!!!!
      if (length(my_list_inheritance1) > 1 & length(unique(my_list_inheritance)) == 1){
        my_list_inheritance1 = my_list_inheritance1[1]
      }
      my_list_inheritance <- translate_inheritance(my_list_inheritance1)
      if (!is.na(my_list_inheritance)){
        EPI_release[row, "Inheritance"] <- my_list_inheritance
        message(paste(row, my_list_inheritance1, my_list_inheritance))
      }
    }
  }
}
#write.csv(EPI_release, "EPI_release_OMIM+.csv", row.names = TRUE)

#Identify not matching OMIM IDs in EPI/MR riles for Manual assessment
i = 0
for (rowi in rownames(EPI_release)){ #change to MR_release for MR files
  row = EPI_release[rowi,] #change to MR_release for MR files
  if (!grepl("\\d{6}", row$OMIM_gene)){next}
  if (!grepl("\\d{6}", row$OMIM_Disease) | row$OMIM_Disease == ""){next}
  if (nchar(row$OMIM_Disease)>6 ){
    Omim_diseases = str_split(row$OMIM_Disease, ";")[[1]]
    for (disease in Omim_diseases){
      genes_OMIM = unique(mimDiseaseGenes$MIM.Number[disease == mimDiseaseGenes$OMIM_Disease & !is.na(mimDiseaseGenes$OMIM_Disease)])
      if (length(genes_OMIM) == 0){next}
      B <- sapply(row$OMIM_gene, grepl, genes_OMIM)
      if (!any(B)){
        print(row)
        print(disease)
        i = i+1
      }
    }
  }
  else{
    genes_OMIM = unique(mimDiseaseGenes$MIM.Number[row$OMIM_Disease == mimDiseaseGenes$OMIM_Disease & !is.na(mimDiseaseGenes$OMIM_Disease)])
    if (length(genes_OMIM) == 0){next}
    B <- sapply(row$OMIM_gene, grepl, genes_OMIM)
    if (!any(B)){
      print(row)
      i = i+1
    }
  }
}

