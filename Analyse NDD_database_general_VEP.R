# ------------------------------------------------------------------------------
# Title: General analysis of NDD dataset
# Author: KN Wijnant
# Date: 20-1-25
# Purpose: For the general analysis of the NDD dataset generated for the systematic evaluation of targetability of AONs with code to re-generate figures
# Inputs: NDD dataset & clinvar variants (supplementary files 1&2)
# Outputs: Figures used to create Figure 2 & supplementary figure 1
# Dependencies: RColorBrewer, readxl, viridis, stringr, dplyr, functions in 'Create_sunburst_plot_var_data_VEP.R'
# Notes: R version 4.4.1 (2024-06-14 ucrt)
# ------------------------------------------------------------------------------

library(RColorBrewer)
library(readxl)
library(viridis)
library(stringr)
library(dplyr)
source("*/Create_Sunburst_plot_var_data_VEP.R")

#analyse NDD database
NDD_database <- read.csv("*/NDD_dataset.csv", row.names=1)

clinvar_variants_NDDonly <- read.csv("*/clinvar_variants.csv", row.names=1, stringsAsFactors = FALSE)


#amount of entries
#analyse NDD list
##Identify number of NDDs with a source score of 3 of higher
sum(NDD_database$source_scores >= 3) #number of NDDs with source score of >3 = 1996
length(unique(NDD_database$HGNC.ID[NDD_database$source_scores >= 3])) #number of unique genes = 1773
length(NDD_database$OMIM_Disease[NDD_database$source_scores >= 3][!is.na(NDD_database$OMIM_Disease[NDD_database$source_scores >= 3])]) #final number of NDDs = 1885


#SOURCES
sources_sep <- lapply(NDD_database$Sources, str_split, ";|:|,")
source_df <- data.frame()
for (source in unique(unlist(sources_sep))){
  in_source <- grepl(source, sources_sep)
  for (score in 1:11){
    is_score <- score == NDD_database$source_scores
    n_disease <- sum(is_score & in_source)
    source_df <- rbind(source_df, c(source, score, n_disease))
  }
}
colnames(source_df) <- c("source", "source_score", "number_of_diseases")
source_df$source_score <- factor(source_df$source_score,levels = c("1", "2", "3", "4", "5", "6", "7", "8", "9", "10", "11"))
source_score_plot <- ggplot(source_df, aes(fill=source, y=as.numeric(number_of_diseases)+5, x=source_score)) + 
  geom_bar(width = 0.7, stat = "identity", position = "dodge") +
  scale_fill_manual(values=c("#003f5c", "#374c80", "#7a5195", "#bc5090", "#ef5675", "#ff764a", "#ffa600")) +
  labs(x = "Source score", y = "Number of entries") +
  scale_y_continuous(limits = c(0,950), expand = expansion(mult = c(0, .1))) +
  theme_bw()
#Venn diagram
#create named list
library(eulerr)
sorted_sources <- sapply(NDD_database$Sources, function(x) paste(sort(unlist(strsplit(x, ";|:|,"))), collapse = ";"))
sources_sep <- lapply(NDD_database$Sources, function(x) sort(unlist(strsplit(x, ";|:|,"))))
source_list <- c()
for (combination in unique(sorted_sources)){
  #print(combination)
  all_sources <- sort(unlist(strsplit(combination, ";|:|,")))
  filter_sorted_sources_all <- rep(0, length(sorted_sources))
  n_comb <- sum(grepl(combination, sorted_sources))
  name <- paste(all_sources, collapse = "&")
  add_source <- setNames(n_comb, name)
  source_list <- c(source_list, add_source)
}
fit <- euler(source_list)
sources_ven <- plot(fit, fill = c("gold2", "coral1", "coral2", "coral3", "cadetblue1", "olivedrab3", "cadetblue3"), quantities = TRUE)


NDD_database = NDD_database[NDD_database$source_scores >= 3 & !is.na(NDD_database$OMIM_Disease),] #apply source score selection & omit NDDs without OMIM annotation
#INHERITANCE
NDD_database_inheritance_class <- NDD_database
inheritance_class <- function(inheritance = inheritance){
  if (grepl("X", inheritance)){
    if (grepl("XLD", inheritance) & grepl("XLR", inheritance)){
      inheritance_class <- "XL"
    }else if (grepl("XLD", inheritance)){
      inheritance_class <- "XLD"
    }else if (grepl("XLR", inheritance)){
      inheritance_class <- "XLR"
    }else{
      inheritance_class <- "XL"
    }
  }else if (grepl("IMP", inheritance)){
    inheritance_class <- "Imprinting"
  }else if (grepl("mtDNA", inheritance)){
    inheritance_class <- "Mitochondrial"
  }else{
    if (grepl("AD", inheritance) & grepl("AR", inheritance)){
      inheritance_class <- "AD, AR"
    }else if (grepl("AD", inheritance)){
      inheritance_class <- "AD"
    }else if (grepl("AR", inheritance)){
      inheritance_class <- "AR"
    }else if (is.na(inheritance)){
      inheritance_class <- "Unknown"
    }else{
      inheritance_class <- "Other"
    }
  }
}
#merge Others plot
NDD_database_inheritance_class$Inheritance_class <- unlist(lapply(NDD_database$Inheritance, inheritance_class))
NDD_database_inheritance_class$Inheritance_class[NDD_database_inheritance_class$Inheritance_class %in% c("Imprinting", "Mitochondrial", "Unknown")] <- "Other"
NDD_data_inheritance <- NDD_database_inheritance_class %>% group_by(Inheritance_class) %>% summarise(Freq = n(), Perc = (n()/nrow(.))*100)
colnames(NDD_data_inheritance) <- c("L1", "Frequency", "Percentage")
sunburst_NDD_inheritance <- create_sunburst(Level_dataframe = NDD_data_inheritance, levels = c("L1"), label = "Levels_count")
#Others pie-chart
NDD_database_inheritance_class$Inheritance_class <- unlist(lapply(NDD_database$Inheritance, inheritance_class))
NDD_database_inheritance_class_others <- NDD_database_inheritance_class[NDD_database_inheritance_class$Inheritance_class %in% c("Imprinting", "Mitochondrial", "Unknown"),]
NDD_data_inheritance_others <- NDD_database_inheritance_class_others %>% group_by(Inheritance_class) %>% summarise(Freq = n(), Perc = (n()/nrow(.))*100)
colnames(NDD_data_inheritance_others) <- c("L1", "Frequency", "Percentage")
sunburst_NDD_inheritance_others <- create_sunburst(Level_dataframe = NDD_data_inheritance_others, levels = c("L1"), label = "Levels_count")


clinvar_variants_NDDonly$consequence_VEP[is.na(clinvar_variants_NDDonly$consequence_VEP)] <- "Not defined"
#CLINVAR figures
#Submission number
clinvar_sep <- unlist(lapply(NDD_database$clinvar_mut_IDs, str_split, ";"))
Clinvar_source_score_select <- clinvar_variants_NDDonly$ID %in% clinvar_sep
sum(Clinvar_source_score_select) #number of unique Clinvar Entries
clinvar_variants_NDDonly_rep <- clinvar_variants_NDDonly[Clinvar_source_score_select,]
sum(as.numeric(clinvar_variants_NDDonly_rep$num_submit), na.rm = TRUE) #number of Clinvar Variants
table(clinvar_variants_NDDonly_rep$consequence_VEP)
clinvar_variants_NDDonly_rep$num_submit[is.na(clinvar_variants_NDDonly_rep$num_submit)] <- 1
clinvar_variants_NDDonly_rep <- clinvar_variants_NDDonly_rep[rep(seq_len(nrow(clinvar_variants_NDDonly_rep)), clinvar_variants_NDDonly_rep$num_submit), ]
clinvar_variants_for_sunburst <- clinvar_variants_NDDonly_rep
clinvar_variants_for_sunburst$consequence_VEP[clinvar_variants_for_sunburst$consequence_VEP %in% c("intergenic_variant", "transcript_ablation", "protein_altering_variant", "coding_sequence_variant", "5_prime_UTR_variant", "stop_lost", "synonymous_variant", "3_prime_UTR_variant", "start_lost", "Not defined", "upstream_gene_variant", "inframe_deletion","intron_variant", "inframe_insertion", "downstream_gene_variant")] <- "Others"
clinvar_variants_for_sunburst$consequence_VEP[clinvar_variants_for_sunburst$consequence_VEP %in% c("splice_polypyrimidine_tract_variant", "splice_donor_region_variant", "splice_donor_5th_base_variant", "splice_region_variant", "splice_acceptor_variant", "splice_donor_variant")] <- "splice_region_variant"
clinvar_variants_for_sunburst$NMD[clinvar_variants_for_sunburst$NMD == "-"] <- NA
var_data_sub_class <- clinvar_variants_for_sunburst %>% group_by(consequence_VEP, NMD) %>% summarise(Freq = n(), Perc = (n()/nrow(.))*100)
colnames(var_data_sub_class) <- c("L1", "L2", "Frequency", "Percentage")
sunburst_sub_class <- create_sunburst(Level_dataframe = var_data_sub_class, levels = c("L1", "L2"), label = "Levels_count")
#others plot
clinvar_variants_for_sunburst <- clinvar_variants_NDDonly_rep
clinvar_variants_for_sunburst$consequence_VEP[clinvar_variants_for_sunburst$consequence_VEP %in% c("intergenic_variant", "transcript_ablation", "protein_altering_variant", "coding_sequence_variant")] <- "Others"
clinvar_variants_for_sunburst <- clinvar_variants_for_sunburst[clinvar_variants_for_sunburst$consequence_VEP %in% c("intergenic_variant", "transcript_ablation", "protein_altering_variant", "coding_sequence_variant", "5_prime_UTR_variant", "stop_lost", "synonymous_variant", "3_prime_UTR_variant", "start_lost", "Not defined", "upstream_gene_variant", "inframe_deletion","intron_variant", "inframe_insertion", "Others", "downstream_gene_variant"),]
clinvar_variants_for_sunburst$NMD[clinvar_variants_for_sunburst$NMD == "-"] <- NA
var_data_sub_class <- clinvar_variants_for_sunburst %>% group_by(consequence_VEP, NMD) %>% summarise(Freq = n(), Perc = (n()/nrow(.))*100)
colnames(var_data_sub_class) <- c("L1", "L2", "Frequency", "Percentage")
sunburst_sub_class_others <- create_sunburst(Level_dataframe = var_data_sub_class, levels = c("L1"), label = "Levels_count")
#Unique variants
clinvar_sep <- unlist(lapply(NDD_database$clinvar_mut_IDs, str_split, ";"))
table(clinvar_variants_NDDonly$class)
clinvar_variants_for_sunburst <- clinvar_variants_NDDonly[Clinvar_source_score_select,]
clinvar_variants_for_sunburst$consequence_VEP[clinvar_variants_for_sunburst$consequence_VEP %in% c("intergenic_variant", "transcript_ablation", "protein_altering_variant", "coding_sequence_variant", "5_prime_UTR_variant", "stop_lost", "synonymous_variant", "3_prime_UTR_variant", "start_lost", "Not defined", "upstream_gene_variant", "inframe_deletion","intron_variant", "inframe_insertion", "downstream_gene_variant")] <- "Others"
clinvar_variants_for_sunburst$consequence_VEP[clinvar_variants_for_sunburst$consequence_VEP %in% c("splice_polypyrimidine_tract_variant", "splice_donor_region_variant", "splice_donor_5th_base_variant", "splice_region_variant", "splice_acceptor_variant", "splice_donor_variant")] <- "splice_region_variant"
clinvar_variants_for_sunburst$NMD[clinvar_variants_for_sunburst$NMD == "-"] <- NA
var_data_sub_class <- clinvar_variants_for_sunburst %>% group_by(consequence_VEP, NMD) %>% summarise(Freq = n(), Perc = (n()/nrow(.))*100)
colnames(var_data_sub_class) <- c("L1", "L2", "Frequency", "Percentage")
sunburst_sub_class <- create_sunburst(Level_dataframe = var_data_sub_class, levels = c("L1", "L2"), label = "Levels_count")
#others plot
clinvar_variants_for_sunburst <- clinvar_variants_NDDonly[Clinvar_source_score_select,]
clinvar_variants_for_sunburst$consequence_VEP[clinvar_variants_for_sunburst$consequence_VEP %in% c("intergenic_variant", "transcript_ablation", "protein_altering_variant", "coding_sequence_variant")] <- "Others"
clinvar_variants_for_sunburst <- clinvar_variants_for_sunburst[clinvar_variants_for_sunburst$consequence_VEP %in% c("intergenic_variant", "transcript_ablation", "protein_altering_variant", "coding_sequence_variant", "5_prime_UTR_variant", "stop_lost", "synonymous_variant", "3_prime_UTR_variant", "start_lost", "Not defined", "upstream_gene_variant", "inframe_deletion","intron_variant", "inframe_insertion", "Others", "downstream_gene_variant"),]
clinvar_variants_for_sunburst$NMD[clinvar_variants_for_sunburst$NMD == "-"] <- NA
var_data_sub_class <- clinvar_variants_for_sunburst %>% group_by(consequence_VEP, NMD) %>% summarise(Freq = n(), Perc = (n()/nrow(.))*100)
colnames(var_data_sub_class) <- c("L1", "L2", "Frequency", "Percentage")
sunburst_sub_class_others <- create_sunburst(Level_dataframe = var_data_sub_class, levels = c("L1"), label = "Levels_count")
