# ------------------------------------------------------------------------------
# Title: Add clinVar data to NDD_dataset
# Author: KN Wijnant
# Date: 27-2-2024
# Purpose: To link clinVar data and VEP predition to the NDD dataset
# Inputs: NDD_dataset file with source scores
# Outputs: clinvar_variants_NDDonly_VEP_predition, NDD_dataset_linked_clinvar
# Dependencies: readxl, plyr, dplyr, stringr, romim, data.table, tidyverse
# Notes: R version 4.4.1 (2024-06-14 ucrt)
# ------------------------------------------------------------------------------


#filter Clinvar variants
NDD_database <- read.csv("NDD_dataset.csv", row.names=1)

#filter clinvar_variant_summary for only likely pathogenic variants in NDD genes
variant_summary <- read.delim("variant_summary.txt", row.names=NULL)
variant_summary_a <- variant_summary[!((variant_summary$GeneID == "-1"| variant_summary$HGNC_ID == "-") & abs(variant_summary$Stop - variant_summary$Start) > 50),]
##only HG38
variant_summary_a <- variant_summary_a[variant_summary_a$Assembly == "GRCh38",]
##only genes in NDD_database or multiple genes & only likely pathogenic
variant_summary_b <- variant_summary_a[(variant_summary_a$HGNC_ID %in% NDD_database$HGNC.ID | variant_summary_a$GeneSymbol %in% NDD_database$HGNC_symbol) | (variant_summary_a$HGNC_ID == "-" & variant_summary_a$GeneSymbol != "-"),]
CLNSGN_select <- c("Conflicting interpretations of pathogenicity", "Likely pathogenic", "Pathogenic", "Pathogenic/Likely pathogenic")
variant_summary_b <- variant_summary_b[variant_summary_b$ClinicalSignificance %in% CLNSGN_select & variant_summary_b$ClinSigSimple == "1",]
for (rowi in 1:nrow(variant_summary_b)){
  row <- variant_summary_b[rowi,]
  HGNC_symbols <- unlist(str_split(row$GeneSymbol, ";"))
  if (any(HGNC_symbols %in% NDD_database$HGNC_symbol)){
    variant_summary_b[rowi,"GeneSymbol"] <- paste(HGNC_symbols[HGNC_symbols %in% NDD_database$HGNC_symbol], collapse = ";")
    if (length(HGNC_symbols[HGNC_symbols %in% NDD_database$HGNC_symbol]) > 1){
      variant_summary_b <- variant_summary_b[-rowi,]
    } else {variant_summary_b[rowi,"GeneSymbol"] <- HGNC_symbols[HGNC_symbols %in% NDD_database$HGNC_symbol]}
  } else {variant_summary_b <- variant_summary_b[-rowi,]}
}
#write.csv(variant_summary_b, "clinvar_variants2_NDDonly.csv", row.names = FALSE)

#link OMIM IDs to ClinVar variants and filter
clinvar_variants_NDDonly <- read.csv("clinvar_variants2_NDDonly.csv")
clinvar_variants_MIM <- clinvar_variants_NDDonly
clinvar_variants_MIM$OMIM_num <- str_count(clinvar_variants_NDDonly$CLNDISDB, "OMIM")
clinvar_variants_MIM$best_OMIM_ID <- NA

for (clinvar_ID in clinvar_variants_MIM$ID){
  row <- clinvar_variants_MIM[clinvar_variants_MIM$ID == clinvar_ID,] 
  if (is.na(row$best_OMIM_ID) | nchar(row$best_OMIM_ID) == 0){
    if (row$OMIM_num == 1){OMIM_ID <- str_extract(str_extract(row$CLNDISDB, "OMIM:\\d+"), "\\d+")
    } else if (row$OMIM_num >= 1){
      OMIM_ID_df <- get_MIMdf_for_multOMIMID(clinvar_ID)
      if (all(is.na(OMIM_ID_df))){OMIM_ID = NA}
      else if (all(is.na(OMIM_ID_df$MIM_ID))){OMIM_ID = NA}
      else{
        OMIM_ID_df <- OMIM_ID_df[!grepl("NA", OMIM_ID_df$MIM_ID),]
        if (all(is.na(OMIM_ID_df))){OMIM_ID = NA} else {OMIM_ID <- paste(OMIM_ID_df$MIM_ID[OMIM_ID_df$num_submit == max(OMIM_ID_df$num_submit)], collapse = ";")}
      }
    }
    else if (row$OMIM_num == 0){OMIM_ID = NA}
    clinvar_variants_MIM$best_OMIM_ID[clinvar_variants_MIM$ID == clinvar_ID] <- OMIM_ID
  }
}

#filter out non-NDD clinvar variants
clinvar_variants_MIM_select <- clinvar_variants_MIM
clinvar_variants_MIM_select <- clinvar_variants_MIM_select[!is.na(clinvar_variants_MIM_select$best_OMIM_ID),]
clinvar_OMIM_disease <- lapply(clinvar_variants_MIM_select$best_OMIM_ID, str_split, ";") 
disease_OMIM_in_NDD_database <- lapply(clinvar_OMIM_disease, function (x) any(unlist(lapply(unlist(x), function(y) y %in% as.character(NDD_database$OMIM_Disease)))))
clinvar_variants_MIM_select$OMIM_in_NDD_database <- unlist(disease_OMIM_in_NDD_database)
clinvar_variants_MIM_select <- clinvar_variants_MIM_select[clinvar_variants_MIM_select$OMIM_in_NDD_database,]
clinvar_variants_MIM_select <- subset(clinvar_variants_MIM_select, select = -c(X, OMIM_num, OMIM_in_NDD_database))
clinvar_variants_MIM_select$class[clinvar_variants_MIM_select$class == "Other"] <- "Others"
#write.csv(clinvar_variants_MIM_select, "clinvar_variants2_NDDonly_OMIM_select.csv", row.names = TRUE)

#add classes
var_names <- clinvar_variants_NDDonly$variation_name
writeLines(var_names, "HGVS_for_VEP_list.txt")
clinvar_variants_NDDonly <- read.csv("clinvar_variants2_NDDonly_OMIM_select.csv")

VEP_output <- read.delim("VEP_output.txt") #can be obtained from VEP (clinvar_variants_NDDonly variants as input)
clinvar_variants_NDDonly$variation_name <- gsub(" ", "", clinvar_variants_NDDonly$variation_name)
VEP_clinvar_df <- left_join(clinvar_variants_NDDonly, VEP_output, by = c("variation_name" = "Uploaded_variation"))

#checks if any NA Consequences and try to obtain the consequence of these
RS_NA <- VEP_clinvar_df$RS...dbSNP.[is.na(VEP_clinvar_df$Consequence)]
RS_NA <- RS_NA[RS_NA != "-1"]
RS_NA <- paste("rs", RS_NA, sep = "")
writeLines(RS_NA, "RSNA_for_VEP_list.txt") #can be obtained from VEP (clinvar_variants_NDDonly variants as input)
RSNA_VEP_output <- read.delim("RSNA_VEP_output.txt")
clinvar_variants_NDDonly$RS...dbSNP. <- paste("rs", clinvar_variants_NDDonly$RS...dbSNP., sep = "")
for (rs in unique(RSNA_VEP_output$X.Uploaded_variation)){
  #print(rs)
  input <- clinvar_variants_NDDonly$variation_name[rs == clinvar_variants_NDDonly$RS...dbSNP.]
  RSNA_VEP_output$X.Uploaded_variation[RSNA_VEP_output$X.Uploaded_variation == rs] <- input
}
colnames(RSNA_VEP_output) <- colnames(VEP_output)

#combine RSNA_VEP_output & VEP_output files
VEP_output_combined <- rbind(VEP_output, RSNA_VEP_output)
VEP_clinvar_df <- left_join(clinvar_variants_NDDonly, VEP_output_combined, by = c("variation_name" = "Uploaded_variation"))
class_dict <- c("frameshift_variant" = "Frameshift",
                "stop_gained" = "Nonsense",
                "splice_acceptor_variant" = "Splice/Intron",
                "splice_donor_variant" = "Splice/Intron",
                "intron_variant" = "Splice/Intron",
                "splice_region_variant" = "Splice/Intron",
                "splice_donor_5th_base_variant" = "Splice/Intron",
                "splice_polypyrimidine_tract_variant" = "Splice/Intron",
                "splice_donor_region_variant" = "Splice/Intron",
                "missense_variant" = "Missense",
                "inframe_deletion" = "Inframe Indel",  
                "inframe_insertion" = "Inframe Indel",
                "synonymous_variant" = "Synonymous",
                "start_lost" = "Start lost",            
                "5_prime_UTR_variant" = "5'UTR",
                "stop_lost" = "Stop lost",
                "3_prime_UTR_variant" = "3'UTR", 
                "downstream_gene_variant" = "downstream_gene_variant",
                "upstream_gene_variant" = "upstream_gene_variant",
                "protein_altering_variant" = "protein_altering_variant",
                "coding_sequence_variant" = "coding_sequence_variant",
                "transcript_ablation" = "transcript_ablation",
                "intergenic_variant" = "intergenic_variant"
                )
severity_hierarchy <- c("transcript_ablation", "splice_acceptor_variant", "splice_donor_variant", 
                        "stop_gained", "frameshift_variant", "stop_lost", "start_lost",
                        "transcript_amplification", "feature_elongation", "feature_truncation", "inframe_insertion", "inframe_deletion",
                        "missense_variant", "protein_altering_variant", "splice_donor_5th_base_variant" ,"splice_region_variant", "splice_donor_region_variant", "splice_polypyrimidine_tract_variant",
                        "incomplete_terminal_codon_variant", "start_retained_variant", "stop_retained_variant",
                        "synonymous_variant", "coding_sequence_variant", "mature_miRNA_variant",
                        "5_prime_UTR_variant", "3_prime_UTR_variant", "non_coding_transcript_exon_variant",
                        "intron_variant", "NMD_transcript_variant", "non_coding_transcript_variant", "coding_transcript_variant",
                        "upstream_gene_variant", "downstream_gene_variant", "TFBS_ablation",
                        "TFBS_amplification", "TF_binding_site_variant", "regulatory_region_ablation",
                        "regulatory_region_amplification", "regulatory_region_variant",
                        "intergenic_variant", "sequence_variant")

#selects most severe consequence
most_severe_consequence <- function(consequences) {
  split_consequences <- strsplit(consequences, ",")[[1]]
  max_severity_index <- which.min(match(split_consequences, severity_hierarchy))
  return(split_consequences[max_severity_index])
}

VEP_clinvar_df$consequence_VEP <- NA
for (rowi in 1:nrow(VEP_clinvar_df)){
  row <- VEP_clinvar_df[rowi,]
  if (!is.na(row$Consequence)){
    if(is.na(class_dict[row$Consequence])){
      #print(row$Consequence)
      if (grepl(",", row$Consequence)){
        if (grepl("splice", row$Consequence)){
          SpliceAI_cols <- colnames(row)[grepl("SpliceAI_pred_DS", colnames(row))]
          if (any(row[,SpliceAI_cols] > 0.2)){
            VEP_clinvar_df[rowi,"consequence_VEP"] <- "splice_region_variant"
          }else(VEP_clinvar_df[rowi,"consequence_VEP"]<- most_severe_consequence(row$Consequence))
        }else(VEP_clinvar_df[rowi,"consequence_VEP"]<- most_severe_consequence(row$Consequence))
      }else(VEP_clinvar_df[rowi,"consequence_VEP"]<- most_severe_consequence(row$Consequence))
    }else(VEP_clinvar_df[rowi,"consequence_VEP"] <- row$Consequence)
  }
}

#correct SNV deletions/
head(VEP_clinvar_df)
VEP_clinvar_df$var_length <- VEP_clinvar_df$Stop - VEP_clinvar_df$Start
VEP_SNVs <- VEP_clinvar_df[VEP_clinvar_df$var_length >= 50,]
VEP_clinvar_df$consequence_VEP[VEP_clinvar_df$var_length >= 50] <- "SNV"
#write.csv(VEP_clinvar_df, "clinvar_variants_NDDonly_VEP_predition.csv")

#link clinvar & NDD dataset
clinvar_per_NDD <- sapply(NDD_database$OMIM_Disease[!is.na(NDD_database$OMIM_Disease)], function(x) paste(VEP_clinvar_df$ID[grepl(x, VEP_clinvar_df$best_OMIM_ID) & !is.na(VEP_clinvar_df$best_OMIM_ID)], collapse = ";"))
NDD_database$clinvar_mut_IDs[!is.na(NDD_database$OMIM_Disease)] <- clinvar_per_NDD
#write.csv(NDD_database, "NDD_dataset_linked_clinvars.csv", row.names= TRUE)


