# ------------------------------------------------------------------------------
# Title: clean DDG2P data
# Author: KN Wijnant
# Date: 14-11-2023
# Purpose: To deduplicate and prepare for integration with other NDD databases 
# Inputs: DDG2P database csv, mimDiseaseTitles file
# Outputs: DDG2P_ID_genelist, manual assessment list & combined manual assessment and DDG2P_ID_genelist
# Dependencies: romim, stringi, stringr
# Notes: R version 4.4.1 (2024-06-14 ucrt)
# ------------------------------------------------------------------------------


library(romim)
library(stringi)
library(stringr)

get_PMID <- function(my_xml){
  my_inheritance_node <- getNodeSet(my_xml, path = "/omim/entryList/entry/referenceList/reference/pubmedID")
  xmlSApply(my_inheritance_node, xmlValue)
}

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
  else if (inheritance_words == "Autosomal dominant; Autosomal recessive"){
    inheritance_abbr = "AD,AR"
  }
  else {
    inheritance_abbr = "NA"
  }
}

DDG2P_genelist <- read.csv("DDG2P_24_5_2022.csv")
DDG2P_genelist <- DDG2P_genelist[,c(1,2,3,4,5,6,7,9,10)]
colnames(DDG2P_genelist) <- c("HGNC_symbol", "OMIM_GeneID", "Disease_name", "OMIM_Disease", "Confidence", "Inheritance", "mutational_consequence", "Organ_system","PMID")

#filter DDG2P for IDs and only select Brain/Cognition
DDG2P_ID_genelist <- DDG2P_genelist[grepl("Brain/Cognition", DDG2P_genelist$Organ_system, fixed = TRUE),]
#DDG2P_ID_genelist$Inheritance == "biallelic_autosomal"
DDG2P_ID_genelist <- data.frame(lapply(DDG2P_ID_genelist, function(x){gsub("biallelic_autosomal", "AR", x)}))
DDG2P_ID_genelist <- data.frame(lapply(DDG2P_ID_genelist, function(x){gsub("monoallelic_autosomal", "AD", x)}))
DDG2P_ID_genelist <- data.frame(lapply(DDG2P_ID_genelist, function(x){gsub("monoallelic_X_het", "XLD", x)}))
DDG2P_ID_genelist <- data.frame(lapply(DDG2P_ID_genelist, function(x){gsub("monoallelic_X_hem", "XLR", x)}))
DDG2P_ID_genelist <- data.frame(lapply(DDG2P_ID_genelist, function(x){gsub("mitochondrial", "mtDNA", x)}))

#replace No disease mim for NA
DDG2P_ID_genelist[DDG2P_ID_genelist[,"OMIM_Disease"] == "No disease mim" & !is.na(DDG2P_ID_genelist[,"OMIM_Disease"]), "OMIM_Disease"] <- rep(NA, length(DDG2P_ID_genelist[DDG2P_ID_genelist[,"OMIM_Disease"] == "No disease mim" & !is.na(DDG2P_ID_genelist[,"OMIM_Disease"]), "OMIM_Disease"]))

#add missing OMIM_Disease data
set_key() #request and insert OMIM key
#PAY attention that you can only send 5000 requests a day (time is set at 7am NY time)
for (row in rownames(DDG2P_ID_genelist)[is.na(DDG2P_ID_genelist$OMIM_Disease)]){
  skip_to_next = FALSE
  gene <- DDG2P_ID_genelist[row, "HGNC_symbol"]
  my_list <- tryCatch(gene_to_omim(gene, show_query = FALSE), error = function(e) { skip_to_next <<- TRUE})
  if (skip_to_next == TRUE | is.null(my_list)){next}
  my_list_omim <- sapply(my_list, get_omim, geneMap = TRUE, referenceList = TRUE)
  my_list_titles <- sapply(my_list_omim, get_title) #get disease names
  diseaseID = names(my_list_titles) #get OMIM_ID
  PMID <- str_split(DDG2P_ID_genelist[row,"PMID"], pattern = ";")[[1]]
  my_list_PMID <- lapply(my_list_omim, get_PMID)
  Inheritance <- sapply(my_list_omim, get_inheritance)
  Inheritance <- sapply(Inheritance, translate_inheritance)
  if (length(names(my_list_PMID)) > 1){
      which_list <- sapply(my_list_PMID, intersect, PMID)
      if (!any(sapply("\\d{6}", grepl, which_list))){next}
      which_list <- which_list[sapply(which_list, length) == max(sapply(which_list, length))] #pick disease with most overlapping references
      DiseaseID <- unique(names(which_list)[sapply(which_list, length) > 0])
      if (length(which_list) >1 & length(DiseaseID) == 1){which_list <- which_list[1]}
      else if (length(which_list) >1 & length(DiseaseID) > 1){
        DiseaseID = DiseaseID[which(lapply(DDG2P_ID_genelist[row, "Inheritance"], grepl, Inheritance)[[1]])]
        if (length(DiseaseID) > 1){message(message(paste("Warning: needs manual evaluation", row, gene, DiseaseID, Inheritance[diseases_include], DDG2P_ID_genelist[row, "Inheritance"], sep = " ")))
          next}
      }
      diseases_include <- which(DiseaseID == my_list)[1]
      if (Inheritance[diseases_include] != DDG2P_ID_genelist[row, "Inheritance"]){
        if (Inheritance[diseases_include] == "AD,AR" & DDG2P_ID_genelist[row, "Inheritance"] == "AR" | DDG2P_ID_genelist[row, "Inheritance"] == "AD"){
          DDG2P_ID_genelist[row, c("Disease_name", "OMIM_Disease")] <- c(paste(my_list_titles[diseases_include], collapse='|' ), paste(diseaseID[diseases_include], collapse="|"))
        }
        else if (Inheritance[diseases_include] == "XL" & DDG2P_ID_genelist[row, "Inheritance"] == "XLD" | DDG2P_ID_genelist[row, "Inheritance"] == "XLR"){
          DDG2P_ID_genelist[row, c("Disease_name", "OMIM_Disease")] <- c(paste(my_list_titles[diseases_include], collapse='|' ), paste(diseaseID[diseases_include], collapse="|"))
        }
        else if (Inheritance[diseases_include] == "XLD,XLR" & DDG2P_ID_genelist[row, "Inheritance"] == "XLD" | DDG2P_ID_genelist[row, "Inheritance"] == "XLR"){
          DDG2P_ID_genelist[row, c("Disease_name", "OMIM_Disease")] <- c(paste(my_list_titles[diseases_include], collapse='|' ), paste(diseaseID[diseases_include], collapse="|"))
        }
        else {message(paste("Warning: Inheritance not equal for", row, gene, DiseaseID, Inheritance[diseases_include], DDG2P_ID_genelist[row, "Inheritance"], sep = " "))}
      }
      else {DDG2P_ID_genelist[row, c("Disease_name", "OMIM_Disease")] <- c(paste(my_list_titles[diseases_include], collapse='|' ), paste(diseaseID[diseases_include], collapse="|"))}
    }
    else if (any(PMID %in% unlist(my_list_PMID))){
      if (all(class(my_list_PMID) == c("matrix", "array"))){
        DiseaseID <- colnames(my_list_PMID)
      }
      else if (class(my_list_PMID) == "list"){
        DiseaseID <- names(my_list_PMID)
      }
      diseases_include <- which(DiseaseID == my_list)[1]
      if (Inheritance != DDG2P_ID_genelist[row, "Inheritance"]){
        if (Inheritance[diseases_include] == "AD,AR" & DDG2P_ID_genelist[row, "Inheritance"] == "AR" | DDG2P_ID_genelist[row, "Inheritance"] == "AD"){
          DDG2P_ID_genelist[row, c("Disease_name", "OMIM_Disease")] <- c(paste(my_list_titles[diseases_include], collapse='|' ), paste(diseaseID[diseases_include], collapse="|"))
        }
        else if (Inheritance[diseases_include] == "XL" & DDG2P_ID_genelist[row, "Inheritance"] == "XLD" | DDG2P_ID_genelist[row, "Inheritance"] == "XLR"){
          DDG2P_ID_genelist[row, c("Disease_name", "OMIM_Disease")] <- c(paste(my_list_titles[diseases_include], collapse='|' ), paste(diseaseID[diseases_include], collapse="|"))
        }
        else if (Inheritance[diseases_include] == "XLD,XLR" & DDG2P_ID_genelist[row, "Inheritance"] == "XLD" | DDG2P_ID_genelist[row, "Inheritance"] == "XLR"){
          DDG2P_ID_genelist[row, c("Disease_name", "OMIM_Disease")] <- c(paste(my_list_titles[diseases_include], collapse='|' ), paste(diseaseID[diseases_include], collapse="|"))
        }
        else if (Inheritance[diseases_include] == "NA"){
          DDG2P_ID_genelist[row, c("Disease_name", "OMIM_Disease")] <- c(paste(my_list_titles[diseases_include], collapse='|' ), paste(diseaseID[diseases_include], collapse="|"))
        }
        else {message(paste("Warning: Inheritance not equal for", row, gene, DiseaseID, Inheritance[diseases_include], DDG2P_ID_genelist[row, "Inheritance"], sep = " "))}
      }
      else {DDG2P_ID_genelist[row, c("Disease_name", "OMIM_Disease")] <- c(paste(my_list_titles[diseases_include], collapse='|' ), paste(diseaseID[diseases_include], collapse="|"))}
    }
}
write.csv(DDG2P_ID_genelist, "DDG2P_ID_genelist.csv", row.names = TRUE)

manual_assesment_db <- data.frame()
#create list for manual evaluation (OMIM disease IDs that could not be matched)
for (rowi in rownames(DDG2P_ID_genelist)[is.na(DDG2P_ID_genelist$OMIM_Disease)]){
  row <- DDG2P_ID_genelist[rowi,]
  skip_to_next = FALSE
  gene <- DDG2P_ID_genelist[rowi, "HGNC_symbol"]
  my_list <- tryCatch(gene_to_omim(gene, show_query = FALSE), error = function(e) { skip_to_next <<- TRUE})
  if (skip_to_next == TRUE | is.null(my_list)){next}
  else {manual_assesment_db <- rbind(manual_assesment_db, row)}
}

write.csv(manual_assesment_db, "manual_assessment_DDG2PR2.csv", row.names = TRUE)

#add correct MIM titles
mimDiseaseTitles <- read.delim("OMIM_database/mimDiseaseTitles.txt", header=TRUE, comment.char="#")
mimDiseaseTitles <- mimDiseaseTitles[,c(2,3)]
colnames(mimDiseaseTitles) <- c("OMIM_Disease", "Disease_name")
mimDiseaseTitles$OMIM_Disease <- as.character(mimDiseaseTitles$OMIM_Disease)

DDG2P_ID_genelist <- read.csv("DDG2P_ID_genelist_manual_assess.csv", row.names = 1)
DDG2P_ID_genelist_manual <- read.csv("manual_assessment_DDG2PR_done.csv", sep=";")
for (row in DDG2P_ID_genelist_manual$Column1[!is.na(DDG2P_ID_genelist_manual$OMIM_Disease)]){
  rowM <- DDG2P_ID_genelist_manual[DDG2P_ID_genelist_manual$Column1 == row,]
  rowD <- DDG2P_ID_genelist[as.numeric(row),]
  if (nchar(rowM$OMIM_Disease) > 6){
    DDG2P_ID_genelist <- DDG2P_ID_genelist[-c(which(rownames(DDG2P_ID_genelist) == row)),]
    for (disease in str_split(rowM$OMIM_Disease, pattern = ",")[[1]]){
      rowD[c("Disease_name", "OMIM_Disease")] <- c(mimDiseaseTitles$Disease_name[mimDiseaseTitles$OMIM_Disease == disease], disease)
      print(rowD)
      DDG2P_ID_genelist <- rbind(DDG2P_ID_genelist, rowD)
    }
  }
  else {DDG2P_ID_genelist[row, c("Disease_name", "OMIM_Disease")] <- c(mimDiseaseTitles$Disease_name[mimDiseaseTitles$OMIM_Disease == rowM$OMIM_Disease], rowM$OMIM_Disease)}
}
#write.csv(DDG2P_ID_genelist, "DDG2P_ID_genelist_manualadded.csv", row.names = TRUE)
