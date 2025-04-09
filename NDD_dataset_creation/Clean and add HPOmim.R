# ------------------------------------------------------------------------------
# Title: clean HPO data
# Author: KN Wijnant
# Date: 9-6-2022
# Purpose: To deduplicate and prepare for integration with other NDD databases 
# Inputs: genes/diseases_for_HP_<0012758, 0001250, 0001249>, mimGenes, mimDiseaseTitles, mimDiseaseGenes files
# Outputs: genes_for_HP file
# Dependencies: romim, stringi, stringr, dplyr, DescTools, readxl
# Notes: R version 4.4.1 (2024-06-14 ucrt)
# ------------------------------------------------------------------------------

library(romim)
library(stringi)
library(stringr)
library(dplyr)
library(DescTools)
library(readxl)

genes_for_HP_0012758 <- read_excel("genes_for_HP_0012758.xlsx")
diseases_for_HP_0012758 <- read_excel("diseases_for_HP_0012758.xlsx")
genes_for_HP_0001250 <- read_excel("genes_for_HP_0001250.xlsx")
diseases_for_HP_0001250 <- read_excel("diseases_for_HP_0001250.xlsx")
genes_for_HP_0001249 <- read_excel("genes_for_HP_0001249.xlsx")
diseases_for_HP_0001249 <- read_excel("diseases_for_HP_0001249.xlsx")

#get OMIM files
mimGenes <- read.delim("mimGenes.txt", header=TRUE, comment.char="#")
mimGenes <- mimGenes[,c(1,4)]
colnames(mimGenes) <- c("OMIM_gene", "HGNC_symbol")
mimDiseaseTitles <- read.delim("mimDiseaseTitles.txt", header=TRUE, comment.char="#")
mimDiseaseTitles <- mimDiseaseTitles[,c(2,3)]
colnames(mimDiseaseTitles) <- c("OMIM_Disease", "Disease_name")
mimDiseaseTitles$OMIM_Disease <- as.character(mimDiseaseTitles$OMIM_Disease)
mimDiseaseGenes <- read.delim("mimDiseaseGenes.txt", comment.char="#")

#translate inheritance
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

#additions to romim package
get_gene_symbols <- function(my_xml){
  my_symbol_node <- getNodeSet(my_xml, path = "/omim/entryList/entry/phenotypeMapList/phenotypeMap/geneSymbols")
  if (length(my_symbol_node) == 0){
    my_symbol_node <- getNodeSet(my_xml, path = "/omim/entryList/entry/geneMap/phenotypeMapList/phenotypeMap/geneSymbols")
  }
  xmlSApply(my_symbol_node, xmlValue)
}

get_approved_symbols <- function(my_xml){
  my_symbol_node <- getNodeSet(my_xml, path = "/omim/entryList/entry/phenotypeMapList/phenotypeMap/approvedGeneSymbols")
  if (length(my_symbol_node) == 0){
    my_symbol_node <- getNodeSet(my_xml, path = "/omim/entryList/entry/geneMap/phenotypeMapList/phenotypeMap/approvedGeneSymbols")
  }
  xmlSApply(my_symbol_node, xmlValue)
}

get_Gene_OMIM_ID <- function(my_xml){
  my_symbol_node <- getNodeSet(my_xml, path = "/omim/entryList/entry/phenotypeMapList/phenotypeMap/mimNumber")
  if (length(my_symbol_node) == 0){
    my_symbol_node <- getNodeSet(my_xml, path = "/omim/entryList/entry/geneMap/phenotypeMapList/phenotypeMap/mimNumber")
  }
  xmlSApply(my_symbol_node, xmlValue)
}

#make dataset per HP number
#"HP_0012758", "HP_0001249", "HP_0001250"
HP_nums <- c("HP_0012758", "HP_0001249", "HP_0001250")
for (HP_num in HP_nums){
  print(HP_num)
  if (HP_num == "HP_0012758"){
    genes_for_HP <- genes_for_HP_0012758
    diseases_for_HP <- diseases_for_HP_0012758
  }
  else if (HP_num == "HP_0001249"){
    genes_for_HP <- genes_for_HP_0001249
    diseases_for_HP <- diseases_for_HP_0001249
  }
  else if (HP_num == "HP_0001250"){
    genes_for_HP <- genes_for_HP_0001250
    diseases_for_HP <- diseases_for_HP_0001250
  }
  #make a row for every disease
  x <- strsplit(as.character(genes_for_HP$DISEASE_IDS), ",", fixed = T)
  genes_for_HP <- cbind(genes_for_HP[rep(1:nrow(genes_for_HP), lengths(x)), 1:2], content = unlist(x))
  
  #extract only OMIM numbers and remove "OMIM:"
  genes_for_HP <- genes_for_HP[grepl("OMIM", genes_for_HP$content, fixed = TRUE),]
  diseases_for_HP <- diseases_for_HP[grepl("OMIM", diseases_for_HP$DISEASE_ID, fixed = TRUE),]
  genes_for_HP$content <- sapply(genes_for_HP$content, function(x) gsub("OMIM:", "", x))
  diseases_for_HP$DISEASE_ID <- sapply(diseases_for_HP$DISEASE_ID, function(x) gsub("OMIM:", "", x))
  colnames(genes_for_HP) <- c("Entrez_ID", "HGNC_symbol", "OMIM_Disease")
  
  #only rows with HP_0012758 phenotype
  genes_for_HP <- genes_for_HP[genes_for_HP$OMIM_Disease %in% diseases_for_HP$DISEASE_ID,]
  
  #add mimGenes
  genes_for_HP <- left_join(genes_for_HP, mimGenes, by = "HGNC_symbol")
  #add disease titles
  genes_for_HP <- left_join(genes_for_HP, mimDiseaseTitles, by = "OMIM_Disease")
  
  #add missing OMIM_Disease data
  set_key() #request and insert OMIM key
  #PAY attention that you can only send 5000 requests a day (time is set at 7am NY time)
  for (row in rownames(genes_for_HP)){
    #skip_to_next = FALSE
    gene <- genes_for_HP[row, "HGNC_symbol"]
    disease <- genes_for_HP[row, "OMIM_Disease"]
    #my_list <- tryCatch(gene_to_omim(gene, show_query = FALSE), error = function(e) { skip_to_next <<- TRUE})
    #if (skip_to_next == TRUE | is.null(my_list)){next}
    my_list_omim <- get_omim(disease, geneMap = TRUE)
    Inheritance <- get_inheritance(my_list_omim)
    Inheritance <- translate_inheritance(Inheritance)
    genes_for_HP[row, "Inheritance"] <- Inheritance
    if (is.na(genes_for_HP[row, "OMIM_gene"])){
      alternative_symbols <- get_gene_symbols(my_list_omim)
      if (length(alternative_symbols) == 0){next}
      symbol_in_alternative <- gene %in% alternative_symbols
      if (symbol_in_alternative == TRUE){
        genes_for_HP[row, "HGNC_symbol"] <- get_approved_symbols(my_list_omim)
        genes_for_HP[row, "OMIM_GeneID"] <- get_Gene_OMIM_ID(my_list_omim)
      }
    }
  }
  write.csv(genes_for_HP, paste("data_", HP_num, ".csv", sep = ""), row.names = TRUE)
}


