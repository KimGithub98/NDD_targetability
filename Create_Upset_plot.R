# ------------------------------------------------------------------------------
# Title: create NDD & patient upset plots
# Author: KN Wijnant
# Date: 21-1-25
# Purpose: To create upset plots on patient and NDD level and to calculate how many patients one can treat with a single AON
# Inputs: NDD dataset & clinvar dataset (with targetability data)
# Outputs: Figures for upset plot
# Dependencies: ComplexUpset, stringr, ggplot2
# Notes: R version 4.4.1 (2024-06-14 ucrt)
# ------------------------------------------------------------------------------

library(ComplexUpset)
library(stringr)
library(ggplot2)
NDD_database_targetable <- read.csv("*/NDD_database_targetable.csv", row.names=1)
NDD_database_targetable$HI <- NDD_database_targetable$HI*1
head(NDD_database_targetable)
NDD_database_targetable[,c("lncRNA", "HPO_deterioration", "HPO_degeneration", "HPO_seizures", "HPO_movement","HPO_targetable", "AD")] <- NDD_database_targetable[,c("lncRNA", "HPO_deterioration", "HPO_degeneration", "HPO_seizures","HPO_movement", "HPO_targetable", "AD")]*1
NDDs_source_score_select <- NDD_database_targetable$source_scores >= 3 & !is.na(NDD_database_targetable$OMIM_Disease)
NDD_database_targetable$perc_LoF_0.2[is.na(NDD_database_targetable$perc_LoF_0.2)] <- FALSE
boolean_cols <- sapply(NDD_database_targetable, is.logical)
NDD_database_targetable[boolean_cols] <- NDD_database_targetable[boolean_cols] * 1

#patient impact analysis
Clinvar_targetable <- read.csv("*/clinvar_variants_targetable.csv", row.names = 1)
Clinvar_for_upset <- data.frame()
for (rowi in 1:nrow(Clinvar_targetable)){
  row <- Clinvar_targetable[rowi,]
  Disease_targetablility <- c(0,0,0,0,0,0)
  OMIM_id <- str_split(row$best_OMIM_ID, ";")[[1]]
  if (length(OMIM_id) >= 1){
    #splice_correct
    splice <- row$splice_correct == "deep_intronic" & !is.na(row$splice_correct)
    #rec_splice_correct
    splice_rec <- row$splice_correct == "recurrent_deep_intronic" & !is.na(row$splice_correct)
    #fs_restore
    ES_RF <- row$exon_skip
    #inframe
    ES_VE <- row$IF_exon_skip > 0 & row$IF_exon_skip < 10 & !is.na(row$IF_exon_skip)
    #disease treatability
    for (OMIM_id_sep in OMIM_id){
      if (OMIM_id_sep %in% NDD_database_targetable$OMIM_Disease){
        Disease_targetablility_entry <- NDD_database_targetable[NDD_database_targetable$OMIM_Disease == OMIM_id_sep & !is.na(NDD_database_targetable$OMIM_Disease) & NDD_database_targetable$HGNC_symbol == row$HGNC_symbol,c(35, 36, 37, 39:41, 46, 52)][1,]
        if (all(is.na(Disease_targetablility_entry))){Disease_targetablility_entry <- data.frame(0,0,0,0,0,0,0,0)}
        if (nrow(Disease_targetablility_entry) != 0){
          Disease_targetablility <- Disease_targetablility + Disease_targetablility_entry
        }
      }
    }
    Disease_targetablility <- Disease_targetablility/length(OMIM_id)
    Disease_targetablility[Disease_targetablility != 0] = 1
    Disease_targetablility <- Disease_targetablility == 1
    #LoF
    LoF <- (row$consequence_VEP == "stop_gained" | row$consequence_VEP == "frameshift_variant" | 
              row$consequence_VEP == "splice_acceptor_variant" | row$consequence_VEP ==  "splice_donor_variant") & row$NMD != "NMD_escaping_variant"
  }
  Clinvar_for_upset <- rbind(Clinvar_for_upset, c(row$AlleleID, row$HGNC_symbol, row$num_submit, 
                                                  row$variation_name, row$class, row$NMD, row$best_OMIM_ID,
                                                  splice, splice_rec, ES_RF, ES_VE, Disease_targetablility, LoF))
}
NDD_OMIM_source_score_select <- NDD_database_targetable$OMIM_Disease[NDDs_source_score_select]
Clinvar_source_score_select <- c()
for (OMIM in Clinvar_targetable$best_OMIM_ID){
  OMIM_ids <- str_split(OMIM, ";")[[1]]
  ClinVar_select<- any(OMIM_ids %in% NDD_OMIM_source_score_select)
  Clinvar_source_score_select <- c(Clinvar_source_score_select, ClinVar_select)
}
colnames(Clinvar_for_upset) <- c("AlleleID", "HGNC_symbol", "num_submit", "variation_name", "class", "is_NMD_escape", "best_OMIM_ID", "Deep_splice_variants", "rec_deep_splice_variants", "fs_restore", "inframe_targetable", "pLI", "AD", "perc_LoF", "TANGO_all", "uORF", "ARE_3UTR", "lncRNA", "HPO_targetable", "LoF")
Clinvar_for_upset$Deep_splice_variants[Clinvar_for_upset$rec_deep_splice_variants == TRUE] = TRUE
Clinvar_for_upset$is_NMD_escape <- Clinvar_for_upset$is_NMD_escape == "NMD_escaping_variant"
Clinvar_for_upset$perc_LoF[is.na(Clinvar_for_upset$perc_LoF)] <- FALSE
ClinVar_ID_source_score_select <- Clinvar_targetable$AlleleID[Clinvar_source_score_select]
ClinVar_ID_source_score_select
#create upset plots
intersections <- list("Deep_splice_variants", c("Deep_splice_variants", "rec_deep_splice_variants"), c("Deep_splice_variants", "HPO_targetable"), c("Deep_splice_variants", "rec_deep_splice_variants", "HPO_targetable"),
                      "fs_restore", c("fs_restore", "HPO_targetable"),
                      "inframe_targetable", c("inframe_targetable", "HPO_targetable"),
                      "TANGO_all", c("TANGO_all", "AD", "pLI"), c("TANGO_all", "AD", "pLI", "LoF"), c("TANGO_all", "AD", "pLI", "LoF", "HPO_targetable"),
                      "uORF", c("uORF", "AD", "pLI"), c("uORF", "AD", "pLI", "LoF"), c("uORF", "AD", "pLI", "LoF", "HPO_targetable"),
                      "ARE_3UTR", c("ARE_3UTR", "AD", "pLI"), c("ARE_3UTR", "AD", "pLI", "LoF"), c("ARE_3UTR", "AD", "pLI", "LoF", "HPO_targetable"),
                      "lncRNA", c("lncRNA", "AD", "pLI"), c("lncRNA", "AD", "pLI", "LoF"), c("lncRNA", "AD", "pLI", "LoF", "HPO_targetable"))
sets <- c("HPO_targetable", "LoF", "pLI","AD", "lncRNA", "ARE_3UTR", "uORF", "TANGO_all", "inframe_targetable", "fs_restore", "rec_deep_splice_variants", "Deep_splice_variants")
Clinvar_for_upset_select <- Clinvar_for_upset[Clinvar_source_score_select,]
Clinvar_for_upset_select[,8:20] <- sapply(Clinvar_for_upset_select[,8:20], \(x) + as.logical(x))
upset_all_ClinVar <- upset(data = Clinvar_for_upset_select[1:50,], intersect = sets, intersections = intersections,
                           mode = "intersect", sort_intersections = FALSE, sort_sets = FALSE, set_sizes=FALSE) #only 1:50 to avoid crashing upset function

#create barplot to combine with upset plot with all the data
Clinvar_for_upset_select$zero <- 0
intersections_n <- c()
for (intersection in intersections){
  intersect_n <- sum(as.numeric(Clinvar_for_upset_select$num_submit[rowSums(Clinvar_for_upset_select[,c(intersection, "zero")]) == length(intersection)]))
  intersections_n <- c(intersections_n, intersect_n)
}
int_as_string <- unlist(lapply(intersections, function(x) paste(x, collapse = "-")))
for_barplot_df <- data.frame(int_as_string, intersections_n)
for_barplot_df$int_as_string <- factor(for_barplot_df$int_as_string, levels = for_barplot_df$int_as_string)
barplot_all_ClinVar <- ggplot(data = for_barplot_df, aes(x = int_as_string, y = intersections_n)) +
  geom_bar(stat = "identity") + 
  theme(axis.text.x = element_text(angle = 90))
for_barplot_df$percent <- for_barplot_df$intersections_n/sum(as.numeric(Clinvar_for_upset_select$num_submit[]))*100 #barplot can be combined with upset plot figure to create a full upset plot without crashing of upset function for a high number of combinations


#get disease impact
NDD_database_targetable_select <- NDD_database_targetable[NDDs_source_score_select & !is.na(NDD_database_targetable$OMIM_Disease),]
colnames(NDD_database_targetable_select)[colnames(NDD_database_targetable_select) == "HI"] <- "pLI"
intersections <- list("Deep_splice_variants", c("Deep_splice_variants", "rec_deep_splice_variants"), c("Deep_splice_variants", "HPO_targetable"), c("Deep_splice_variants", "rec_deep_splice_variants", "HPO_targetable"),
                      "fs_restore", c("fs_restore", "HPO_targetable"),
                      "inframe_targetable", c("inframe_targetable", "HPO_targetable"),
                      "TANGO_all", c("TANGO_all", "AD", "pLI"), c("TANGO_all", "AD", "pLI", "perc_LoF_0.2"), c("TANGO_all", "AD", "pLI", "perc_LoF_0.2", "HPO_targetable"),
                      "uORF", c("uORF", "AD", "pLI"), c("uORF", "AD", "pLI", "perc_LoF_0.2"), c("uORF", "AD", "pLI", "perc_LoF_0.2", "HPO_targetable"),
                      "ARE_3UTR", c("ARE_3UTR", "AD", "pLI"), c("ARE_3UTR", "AD", "pLI", "perc_LoF_0.2"), c("ARE_3UTR", "AD", "pLI", "perc_LoF_0.2", "HPO_targetable"),
                      "lncRNA", c("lncRNA", "AD", "pLI"), c("lncRNA", "AD", "pLI", "perc_LoF_0.2"), c("lncRNA", "AD", "pLI", "perc_LoF_0.2", "HPO_targetable"))

sets <- c("HPO_targetable", "AD", "perc_LoF_0.2", "pLI", "lncRNA", "ARE_3UTR", "uORF", "TANGO_all", "inframe_targetable", "fs_restore", "rec_deep_splice_variants", "Deep_splice_variants")
#create upset plot
upset_select_NDDs <- upset(data = NDD_database_targetable_select, intersect = sets, intersections = intersections,
                           mode = "intersect", sort_intersections = FALSE, sort_sets = FALSE)

#get the number of overall treatable NDDs (union of all strategies)
select_GUI <- NDD_database_targetable$Gene_Disease_UI[NDDs_source_score_select & !is.na(NDD_database_targetable$OMIM_Disease)]
intersections <- list("Deep_splice_variants", c("Deep_splice_variants", "rec_deep_splice_variants"), c("Deep_splice_variants", "HPO_targetable"), c("Deep_splice_variants", "rec_deep_splice_variants", "HPO_targetable"),
                      "fs_restore", c("fs_restore", "HPO_targetable"),
                      "inframe_targetable", c("inframe_targetable", "HPO_targetable"),
                      "TANGO_all", c("TANGO_all", "AD", "pLI"), c("TANGO_all", "AD", "pLI", "perc_LoF_0.2"), c("TANGO_all", "AD", "pLI", "perc_LoF_0.2", "HPO_targetable"),
                      "uORF", c("uORF", "AD", "pLI"), c("uORF", "AD", "pLI", "perc_LoF_0.2"), c("uORF", "AD", "pLI", "perc_LoF_0.2", "HPO_targetable"),
                      "ARE_3UTR", c("ARE_3UTR", "AD", "pLI"), c("ARE_3UTR", "AD", "pLI", "perc_LoF_0.2"), c("ARE_3UTR", "AD", "pLI", "perc_LoF_0.2", "HPO_targetable"),
                      "lncRNA", c("lncRNA", "AD", "pLI"), c("lncRNA", "AD", "pLI", "perc_LoF_0.2"), c("lncRNA", "AD", "pLI", "perc_LoF_0.2", "HPO_targetable"))
colnames(NDD_database_targetable)[colnames(NDD_database_targetable) == "HI"] <- "pLI"
GUI_splice <-NDD_database_targetable$Gene_Disease_UI[rowSums(NDD_database_targetable[,c("Deep_splice_variants", "HPO_targetable")]) == 2]
GUI_fsrestore <- NDD_database_targetable$Gene_Disease_UI[rowSums(NDD_database_targetable[,c("fs_restore", "HPO_targetable")]) == 2]
GUI_inframe <- NDD_database_targetable$Gene_Disease_UI[rowSums(NDD_database_targetable[,c("inframe_targetable", "HPO_targetable")]) == 2]
GUI_TANGO <- NDD_database_targetable$Gene_Disease_UI[rowSums(NDD_database_targetable[,c("TANGO_all", "AD", "pLI", "perc_LoF_0.2", "HPO_targetable")]) == 5]
GUI_uORF <- NDD_database_targetable$Gene_Disease_UI[rowSums(NDD_database_targetable[,c("uORF", "AD", "pLI", "perc_LoF_0.2", "HPO_targetable")]) == 5]
GUI_ARE <- NDD_database_targetable$Gene_Disease_UI[rowSums(NDD_database_targetable[,c("ARE_3UTR", "AD", "pLI", "perc_LoF_0.2", "HPO_targetable")]) == 5]
GUI_lncRNA <- NDD_database_targetable$Gene_Disease_UI[rowSums(NDD_database_targetable[,c("lncRNA", "AD", "pLI", "perc_LoF_0.2", "HPO_targetable")]) == 5]
treatable_GUI <- unique(c(GUI_splice, GUI_fsrestore, GUI_inframe, GUI_TANGO, GUI_uORF, GUI_ARE, GUI_lncRNA))
treatable_GUI_sc3 <- intersect(treatable_GUI, select_GUI)
length(treatable_GUI) #number of targetable NDDs
length(treatable_GUI_sc3) #number of targetable NDDs with source score of higher then 3


#get the number of overall treatable Clinvar submittions (union of all strategies)
GUI_splice <- Clinvar_for_upset_select$AlleleID[rowSums(Clinvar_for_upset_select[,c("Deep_splice_variants", "HPO_targetable")]) == 2]
GUI_recurrent <- Clinvar_for_upset_select$AlleleID[rowSums(Clinvar_for_upset_select[,c("rec_deep_splice_variants", "HPO_targetable")]) == 2]
GUI_fsrestore <- Clinvar_for_upset_select$AlleleID[rowSums(Clinvar_for_upset_select[,c("fs_restore", "HPO_targetable")]) == 2]
GUI_inframe <- Clinvar_for_upset_select$AlleleID[rowSums(Clinvar_for_upset_select[,c("inframe_targetable", "HPO_targetable")]) == 2]
GUI_TANGO <- Clinvar_for_upset_select$AlleleID[rowSums(Clinvar_for_upset_select[,c("TANGO_all", "AD", "HPO_targetable", "LoF")]) == 4]
GUI_uORF <- Clinvar_for_upset_select$AlleleID[rowSums(Clinvar_for_upset_select[,c("uORF", "AD", "HPO_targetable", "LoF")]) == 4]
GUI_ARE <- Clinvar_for_upset_select$AlleleID[rowSums(Clinvar_for_upset_select[,c("ARE_3UTR", "AD", "HPO_targetable", "LoF")]) == 4]
GUI_lncRNA <- Clinvar_for_upset_select$AlleleID[rowSums(Clinvar_for_upset_select[,c("lncRNA", "AD", "HPO_targetable", "LoF")]) == 4]
treatable_ID <- unique(c(GUI_splice, GUI_fsrestore, GUI_inframe, GUI_TANGO, GUI_uORF, GUI_ARE, GUI_lncRNA))
sum(as.numeric(Clinvar_for_upset_select$num_submit[Clinvar_for_upset_select$AlleleID %in% treatable_ID])) #number of targetable ClinVar submittions
sum(as.numeric(Clinvar_for_upset_select$num_submit)) #total number of targetable ClinVar submittions

#number of patients treatable with a single AON
Clinvar_targetable_select <- Clinvar_targetable[Clinvar_source_score_select,]
tot <- sum(as.numeric(Clinvar_for_upset_select$num_submit[Clinvar_for_upset_select$AlleleID %in% GUI_splice]))
aons <- length(as.numeric(Clinvar_for_upset_select$num_submit[Clinvar_for_upset_select$AlleleID %in% GUI_splice]))
100/(100*aons/tot) #number of patients per AON (splice correction)
tot <- sum(as.numeric(Clinvar_for_upset_select$num_submit[Clinvar_for_upset_select$AlleleID %in% GUI_recurrent]))
aons <- length(as.numeric(Clinvar_for_upset_select$num_submit[Clinvar_for_upset_select$AlleleID %in% GUI_recurrent]))
100/(100*aons/tot) #number of patients per AON (splice correction recurrent)
tot <- sum(as.numeric(Clinvar_for_upset_select$num_submit[Clinvar_for_upset_select$AlleleID %in% GUI_fsrestore]))
aons <- c()
for (allele_ID in GUI_fsrestore){
  row <- Clinvar_targetable_select[Clinvar_targetable_select$AlleleID == allele_ID,]
  exon_inf <- unlist(str_split(row$exons_inf, ";"))
  exon_inf <- lapply(exon_inf, function(x) str_split(x, "-"))
  up_bound_i_start <- min(which(as.numeric(unlist(exon_inf))[order(as.numeric(unlist(exon_inf)), decreasing = T)] < row$Start))
  if (!is.finite(up_bound_i_start)){up_bound_i_start = length(unlist(exon_inf))+1}
  low_bound_i_start <- up_bound_i_start -1
  up_bound_i_stop <- min(which(as.numeric(unlist(exon_inf))[order(as.numeric(unlist(exon_inf)), decreasing = T)] < row$Stop))
  if (!is.finite(up_bound_i_start)){up_bound_i_start = 1}
  low_bound_i_stop <- up_bound_i_stop -1
  if (up_bound_i_start/2 == round(up_bound_i_start/2)| up_bound_i_stop/2 == round(up_bound_i_stop/2) | up_bound_i_start == up_bound_i_stop){in_exon = TRUE}else{
    in_exon = FALSE
    if (low_bound_i_stop == 0){low_bound_i_stop = 1}
    if (up_bound_i_start == 0){up_bound_i_start = 1}
    del_exon_start <- as.numeric(as.numeric(unlist(exon_inf))[order(as.numeric(unlist(exon_inf)), decreasing = T)][low_bound_i_start])
    del_exon_stop <- as.numeric(as.numeric(unlist(exon_inf))[order(as.numeric(unlist(exon_inf)), decreasing = T)][up_bound_i_stop])
    exon_df <- sapply(exon_inf, "[[", 1)
    if (ncol(exon_df) > 1){
      intron_df <- rbind(exon_df[2,1:length(exon_inf)-1], exon_df[1,2:length(exon_inf)])
      intron_length <- -abs(as.numeric(intron_df[1,]) - as.numeric(intron_df[2,]))
      remove_intron <- sum(intron_length[(dplyr::between(as.numeric(intron_df[1,]), row$Start, row$Stop) + dplyr::between(as.numeric(intron_df[2,]), row$Start, row$Stop)) > 1])
    }else {remove_intron = 0}
    in_frame <- round((abs(del_exon_start- del_exon_stop)-abs(remove_intron )+1)/3) == (abs(del_exon_start- del_exon_stop)-abs(remove_intron )+1)/3
    if (!in_frame){
      phase <- round(round((abs(del_exon_start- del_exon_stop)+1)/3) - (abs(del_exon_stop- del_exon_start)+1)/3, digits = 2)
      exon_df <- rbind(exon_df, abs(as.numeric(exon_df[1,]) - as.numeric(exon_df[2,]))+1)
      exon_df <- rbind(exon_df, round(as.numeric(exon_df[3,])/3) == as.numeric(exon_df[3,])/3)
      exon_df <- rbind(exon_df, round(as.numeric(exon_df[3,])/3) - as.numeric(exon_df[3,])/3)
      exon_num <- which(row$Start < as.numeric(exon_df[1,]) & row$Stop > as.numeric(exon_df[2,]))
      if (ncol(exon_df)<4){next}
      if (any(!exon_num %in% 2:(ncol(exon_df)-1))){next}
      pre_exon <- min(exon_num) - 1
      fol_exon <- max(exon_num) + 1
      if(pre_exon %in% 2:(ncol(exon_df)-1)){
        pre_exon_f <- round(as.numeric(exon_df[5,min(exon_num) - 1]), digits = 2)
        if(phase+pre_exon_f == 0){
          nt_del<- sum(as.numeric(exon_df[3,c(min(exon_num) - 1, exon_num)]))
          aons <- c(aons, paste(c(row$HGNC_symbol, pre_exon), collapse = "_"))
        }
      }
      if(fol_exon %in% 2:(ncol(exon_df)-1)){
        fol_exon_f <- round(as.numeric(exon_df[5,max(exon_num) + 1]), digits= 2)
        if(phase+fol_exon_f == 0){
          nt_del<- sum(as.numeric(exon_df[3,c(exon_num, max(exon_num) + 1)]))
          aons <- c(aons, paste(c(row$HGNC_symbol, fol_exon), collapse = "_"))
        }
      }
    }
  }
}
table(aons)

## targetable frameshift splice mutations
skips_exon_or_intron_retention <- function(allele_ID){
  row <- Clinvar_targetable_select[Clinvar_targetable_select$AlleleID == allele_ID,]
  exon_inf <- row$exons_inf
  exon_inf <- unlist(str_split(exon_inf, ";"))
  exon_inf <- lapply(exon_inf, function(x) str_split(x, "-"))
  exon_df <- data.frame(sapply(exon_inf, "[[", 1))
  exon_df <- mutate_all(exon_df, function (x) as.numeric(x))
  if (ncol(exon_df)<4){return(-1)}
  distance_df <- abs(exon_df - row$POS)
  if(any(distance_df < 3)){
    closest_exon_index <- which(distance_df[1,] == min(distance_df) | distance_df[2,] == min(distance_df))
    if (closest_exon_index == 1 | closest_exon_index == ncol(exon_df)){return(-1)}
    closest_exon_num <- exon_df[,distance_df[1,] == min(distance_df) | distance_df[2,] == min(distance_df)]
    in_frame <- round((abs(closest_exon_num[2]- closest_exon_num[1])+1)/3) == (abs(closest_exon_num[2]- closest_exon_num[1])+1)/3
    if (!in_frame){
      perc_del = -1
      exon_df <- rbind(exon_df, abs(as.numeric(exon_df[1,]) - as.numeric(exon_df[2,]))+1)
      tot_len <- row$cds_length
      exon_df <- rbind(exon_df, round(as.numeric(exon_df[3,])/3) == as.numeric(exon_df[3,])/3)
      exon_df <- rbind(exon_df, round(as.numeric(exon_df[3,])/3) - as.numeric(exon_df[3,])/3)
      phase <- round(exon_df[5,closest_exon_index], digits = 2)
      pre_exon <- min(closest_exon_index) - 1
      fol_exon <- max(closest_exon_index) + 1
      if(pre_exon %in% 2:(ncol(exon_df)-1)){
        pre_exon_f <- round(as.numeric(exon_df[5,pre_exon]), digits = 2)
        if(phase+pre_exon_f == 0){
          nt_del<- sum(as.numeric(exon_df[3,c(pre_exon, closest_exon_index)]))
          aon <- paste(c(row$HGNC_symbol, pre_exon), collapse = "_")
        }
      }
      if(fol_exon %in% 2:(ncol(exon_df)-1)){
        fol_exon_f <- round(as.numeric(exon_df[5,fol_exon]), digits= 2)
        if(phase+fol_exon_f == 0){
          nt_del<- sum(as.numeric(exon_df[3,c(closest_exon_index, fol_exon)]))
          aon <- paste(c(row$HGNC_symbol, fol_exon), collapse = "_")
        }
      }
      return(aon)
    }
    return(NA)
  }
}
aons_splice <- lapply(GUI_fsrestore, skips_exon_or_intron_retention)
aons <- length(unique(c(aons_splice, aons)))
100/(100*aons/tot) #number of patients per AON (ES-RF targeting)
tot <- sum(as.numeric(Clinvar_for_upset_select$num_submit[Clinvar_for_upset_select$AlleleID %in% GUI_inframe]))
in_inframe_exon <- function(allele_ID){
  del_perc = -1
  row <- Clinvar_targetable_select[Clinvar_targetable_select$AlleleID == allele_ID,]
  row <- row[1,]
  exon_inf <- row$exons_inf
  if (is.na(exon_inf)){return(NA)}
  position <- row$POS
  exon_num <- which(unlist(lapply(str_split(unlist(str_split(exon_inf, ";")), "-"), function(x) dplyr::between(as.numeric(row$POS), as.numeric(x[1]), as.numeric(x[2])) | dplyr::between(as.numeric(row$POS), as.numeric(x[2]), as.numeric(x[1])))))
  return(paste(c(row$HGNC_symbol, exon_num), collapse  = "_"))
}
aons <- length(unique(lapply(GUI_inframe,in_inframe_exon)))
100/(100*aons/tot) #number of patients per AON (ES-VE targeting)
tot <- sum(as.numeric(Clinvar_for_upset_select$num_submit[Clinvar_for_upset_select$AlleleID %in% GUI_TANGO]))
aons <- length(unique(Clinvar_for_upset_select$HGNC_symbol[Clinvar_for_upset_select$AlleleID %in% GUI_TANGO]))
100/(100*aons/tot) #number of patients per AON (TANGO targeting)
tot <- sum(as.numeric(Clinvar_for_upset_select$num_submit[Clinvar_for_upset_select$AlleleID %in% GUI_uORF]))
aons <- length(unique(Clinvar_for_upset_select$HGNC_symbol[Clinvar_for_upset_select$AlleleID %in% GUI_uORF]))
100/(100*aons/tot) #number of patients per AON (uORF targeting)
tot <- sum(as.numeric(Clinvar_for_upset_select$num_submit[Clinvar_for_upset_select$AlleleID %in% GUI_ARE]))
aons <- length(unique(Clinvar_for_upset_select$HGNC_symbol[Clinvar_for_upset_select$AlleleID %in% GUI_ARE]))
100/(100*aons/tot) #number of patients per AON (ARE targeting)
tot <- sum(as.numeric(Clinvar_for_upset_select$num_submit[Clinvar_for_upset_select$AlleleID %in% GUI_lncRNA]))
aons <- length(unique(Clinvar_for_upset_select$HGNC_symbol[Clinvar_for_upset_select$AlleleID %in% GUI_lncRNA]))
100/(100*aons/tot) #number of patients per AON (lncRNA targeting)
