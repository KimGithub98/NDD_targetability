# ------------------------------------------------------------------------------
# Title: Add external links and prot length information to NDD datase
# Author: KN Wijnant
# Date: 17-2-2023
# Purpose: To add external links and prot length information to NDD datase
# Inputs: NDD_dataset file with source scores & ClinVar links
# Outputs: NDD_dataset file with external links
# Dependencies: stringr, biomaRt
# Notes: R version 4.4.1 (2024-06-14 ucrt)
# ------------------------------------------------------------------------------


library(stringr)
library(biomaRt)
mart <- useMart('ENSEMBL_MART_ENSEMBL')
mart <- useDataset('hsapiens_gene_ensembl', mart)

HGNC_information_externa_links <- read.delim("HGNC_information_externa_links.txt")
NDD_database <- read.csv("NDD_dataset.csv", row.names=1)

HGNC_addtoNDD <- HGNC_information_externa_links[,c("Approved.symbol", "Ensembl.ID.supplied.by.Ensembl.", "MANE.Select.Ensembl.transcript.ID..supplied.by.NCBI.", "UniProt.ID.supplied.by.UniProt.", "RefSeq.IDs")]
colnames(HGNC_addtoNDD) <- c("HGNC_symbol", "EnsemblID", "MANE_EnsemblID", "Uniprot_ID", "RefseqIDs")

#change old HGNC symbols to new symbols
not_in_HGNC <- NDD_database$HGNC_symbol[!NDD_database$HGNC_symbol %in% HGNC_addtoNDD$HGNC_symbol]
old_sym <- HGNC_information_externa_links$Previous.symbols[HGNC_information_externa_links$Previous.symbols %in% not_in_HGNC]
new_sym <- HGNC_information_externa_links$Approved.symbol[HGNC_information_externa_links$Previous.symbols %in% not_in_HGNC]
for (i in 1:length(old_sym)){
  NDD_database$HGNC_symbol[NDD_database$HGNC_symbol == old_sym[i]] <- new_sym[i]
}

#add external links to NDD database
NDD_join <- left_join(NDD_database, HGNC_addtoNDD, by = "HGNC_symbol")

#retrieve coding sequence and add length to NDD db
ENSTs <- str_extract(NDD_join$MANE_EnsemblID, "ENST\\d{11}")
ENSTs_valid <- ENSTs[nchar(ENSTs) == 15]
ENST_coding_length <- sapply(ENSTs, function(x) (nchar(getSequence(id = x, 
                                                                   type = "ensembl_transcript_id", 
                                                                   seqType = "coding", 
                                                                   mart = mart)[[1]])/3)-1)
ENST_coding_length[lapply(names(ENST_coding_length), function (x) is.null(x)|is.na(x)) == TRUE] <- NA
names(ENST_coding_length[lapply(names(ENST_coding_length), function (x) is.null(x)|is.na(x)) == TRUE]) <- "NO_ESTN"
NDD_join <- cbind(NDD_join, unlist(ENST_coding_length))
colnames(NDD_join)[17] <- "ENST_coding_length"

#choose unique UNIprot ID
double_UNIprot <- NDD_join[lapply(str_split(NDD_join$Uniprot_ID, ", "), length) > 1,]
#biomaRt does not have an option for this, do this manually

#calculate protein length
for (rowi in 1:length(NDD_database$HGNC_symbol)){
  UNI_ID <- NDD_database$Uniprot_ID[rowi]
  if (length(UNI_ID) > 1 | any(is.na(UNI_ID)) | any(nchar(UNI_ID) == 0)){next}
  rel_json <- drawProteins::get_features(UNI_ID)
  rel_data <- drawProteins::feature_to_dataframe(rel_json)
  length = rel_data$end[rel_data$type == "CHAIN"][1]
  length_ENS = NDD_database$ENST_coding_length[rowi]
  if (is.na(length_ENS)){next}
  if (length - length_ENS != 0){
    print(length - length_ENS)
  }
  NDD_database$UNIprot_length[rowi] <- length
}
NDD_database <- cbind(rownames(NDD_database), NDD_database)
colnames(NDD_database)[1] <- "Gene_Disease_UI"
#write.csv(NDD_database, "NDD_dataset_externaldb_links_UNIprotLength.csv")
