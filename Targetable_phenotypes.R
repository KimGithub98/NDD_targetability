# ------------------------------------------------------------------------------
# Title: get targetable phenotypic features and create upset plots
# Author: KN Wijnant
# Date: 21-1-25
# Purpose: To create upset plots of targetable phenotypic features 
# Inputs: NDD dataset , HPO phenotype_to_genes file
# Outputs: NDD dataset with phenotypic feature targetability
# Dependencies: ComplexUpset, stringr, ggplot2
# Notes: R version 4.4.1 (2024-06-14 ucrt)
# ------------------------------------------------------------------------------

library(stringr)
library(ComplexUpset)
library(ggplot2)

#Add phenotype
NDD_database_targetable <- read.csv("*/NDD_database_targetable.csv", row.names=1)
HPO_mim_phenotypes <- read.delim("*/phenotype_to_genes.txt", header=FALSE, comment.char="#")
colnames(HPO_mim_phenotypes) <- c("HPO_term", "Name", "entrez_gene_id", "entrez_gene_symbol", "Additional_Info", "Source", "disease_ID")

extract_OMIM <- function(disease_ID){
  if (grepl("OMIM:(\\d+)", disease_ID))
    OMIM_ID <- gsub("OMIM:(\\d+)", "\\1", disease_ID)
  else (OMIM_ID = NA)
  return(OMIM_ID)
}
HPO_mim_phenotypes$OMIM_Disease <- unlist(lapply(HPO_mim_phenotypes$disease_ID, function (x) extract_OMIM(x)))


# HP:0002376, HP:0001268 HP:0002505 deterioration
# HP:0007367 neurodegeneration
# HP:0001250 seizures
# HP:0003811 neonathal lethality


#get OMIM for deterioration
HPO_deterioration <- HPO_mim_phenotypes[HPO_mim_phenotypes$HPO_term == "HP:0002376" | HPO_mim_phenotypes$HPO_term == "HP:0001268"| HPO_mim_phenotypes$HPO_term == "HP:0002505",]
OMIM_deterioration <- intersect(NDD_database_targetable$OMIM_Disease, HPO_deterioration$OMIM_Disease)
OMIM_deterioration <- OMIM_deterioration[!is.na(OMIM_deterioration)]

#get OMIM for degeneration
HPO_degeneration <- HPO_mim_phenotypes[HPO_mim_phenotypes$HPO_term == "HP:0007367",]
OMIM_degeneration <- intersect(NDD_database_targetable$OMIM_Disease, HPO_degeneration$OMIM_Disease)
OMIM_degeneration <- OMIM_degeneration[!is.na(OMIM_degeneration)]

#get OMIM for seizures
HPO_seizures <- HPO_mim_phenotypes[HPO_mim_phenotypes$HPO_term == "HP:0001250",]
OMIM_seizures <- intersect(NDD_database_targetable$OMIM_Disease, HPO_seizures$OMIM_Disease)
OMIM_seizures <- OMIM_seizures[!is.na(OMIM_seizures)]

#get OMIM for movement
HPO_movement <- HPO_mim_phenotypes[HPO_mim_phenotypes$HPO_term == "HP:0100022",]
OMIM_movement <- intersect(NDD_database_targetable$OMIM_Disease, HPO_movement$OMIM_Disease)
OMIM_movement <- OMIM_movement[!is.na(OMIM_movement)]

#get OMIM for neonathal lethality
HPO_neo_leth <- HPO_mim_phenotypes[HPO_mim_phenotypes$HPO_term == "HP:0003811",]
OMIM_neo_leth <- intersect(NDD_database_targetable$OMIM_Disease, HPO_neo_leth$OMIM_Disease)
OMIM_neo_leth <- OMIM_neo_leth[!is.na(OMIM_neo_leth)]
OMIM_not_neo_leth <- setdiff(NDD_database_targetable$OMIM_Disease, HPO_neo_leth$OMIM_Disease)

HPO_targetable <- unique(c(OMIM_deterioration, OMIM_degeneration, OMIM_seizures, OMIM_movement))
HPO_targetable <- HPO_targetable[!HPO_targetable %in% OMIM_neo_leth]

#add to NDD_database_targetable
NDD_database_targetable$HPO_deterioration <- NDD_database_targetable$OMIM_Disease %in% OMIM_deterioration
NDD_database_targetable$HPO_degeneration <- NDD_database_targetable$OMIM_Disease %in% OMIM_degeneration
NDD_database_targetable$HPO_seizures <- NDD_database_targetable$OMIM_Disease %in% OMIM_seizures
NDD_database_targetable$HPO_movement <- NDD_database_targetable$OMIM_Disease %in% OMIM_movement
NDD_database_targetable$HPO_neo_lethal <- NDD_database_targetable$OMIM_Disease %in% OMIM_neo_leth
NDD_database_targetable$HPO_targetable <- NDD_database_targetable$OMIM_Disease %in% HPO_targetable
sum(!NDD_database_targetable$HPO_targetable[NDD_database_targetable$source_scores >= 3])

#figure target ability
library(ComplexUpset)
binary_matrix <- NDD_database_targetable[NDD_database_targetable$source_scores >= 3,c("HPO_deterioration", "HPO_degeneration", "HPO_seizures", "HPO_movement", "HPO_neo_lethal")]
binary_matrix$No_HPO_Sets <- apply(binary_matrix, 1, function(row) all(row == FALSE))
# Convert to factor for the upset plot
binary_matrix$No_HPO_Sets <- factor(binary_matrix$No_HPO_Sets, levels = c(TRUE, FALSE))
intersections <- list(
  "HPO_seizures",
  c("HPO_seizures", "HPO_movement"),
  c("HPO_seizures", "HPO_deterioration"),
  c("HPO_seizures", "HPO_degeneration"),
  c("HPO_seizures", "HPO_movement", "HPO_deterioration"),
  c("HPO_seizures", "HPO_movement", "HPO_degeneration"),
  c("HPO_seizures", "HPO_deterioration", "HPO_degeneration"),
  c("HPO_seizures", "HPO_movement", "HPO_deterioration", "HPO_degeneration"),
  "HPO_movement",
  c("HPO_movement", "HPO_deterioration"),
  c("HPO_movement", "HPO_degeneration"),
  c("HPO_movement", "HPO_deterioration", "HPO_degeneration"),
  "HPO_deterioration",
  c("HPO_deterioration", "HPO_degeneration"),
  "HPO_degeneration",
  "HPO_neo_lethal",
  c("HPO_neo_lethal", "HPO_seizures"),
  c("HPO_neo_lethal", "HPO_movement"),
  c("HPO_neo_lethal", "HPO_deterioration"),
  c("HPO_neo_lethal", "HPO_degeneration"),
  c("HPO_neo_lethal", "HPO_seizures", "HPO_movement"),
  c("HPO_neo_lethal", "HPO_seizures", "HPO_deterioration"),
  c("HPO_neo_lethal", "HPO_seizures", "HPO_degeneration"),
  c("HPO_neo_lethal", "HPO_movement", "HPO_deterioration"),
  c("HPO_neo_lethal", "HPO_movement", "HPO_degeneration"),
  c("HPO_neo_lethal", "HPO_deterioration", "HPO_degeneration"),
  c("HPO_neo_lethal", "HPO_seizures", "HPO_movement", "HPO_deterioration"),
  c("HPO_neo_lethal", "HPO_seizures", "HPO_movement", "HPO_degeneration"),
  c("HPO_neo_lethal", "HPO_seizures", "HPO_deterioration", "HPO_degeneration"),
  c("HPO_neo_lethal", "HPO_movement", "HPO_deterioration", "HPO_degeneration"),
  c("HPO_neo_lethal", "HPO_seizures", "HPO_movement", "HPO_deterioration", "HPO_degeneration"),
  "No_HPO_Sets"
)
upset_phenotypes <- upset(
  binary_matrix,
  intersect = c("HPO_neo_lethal", "HPO_degeneration",  "HPO_deterioration", "HPO_movement","HPO_seizures" ,"No_HPO_Sets"),
  intersections = intersections,
  name = "HPO Intersections",
  width_ratio = 0.2, 
  mode = "distinct", sort_intersections = FALSE, sort_sets = FALSE, set_sizes=FALSE
)
