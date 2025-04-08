# ------------------------------------------------------------------------------
# Title: Targetability analysis of NDD dataset
# Author: KN Wijnant
# Date: 20-1-25
# Purpose: For the targetability analysis of the NDD dataset generated for the systematic evaluation of targetability of AONs with code to re-generate figures & re-assess the targetability of NDDs/ClinVar variants of 7 AON strategies
# Inputs: NDD dataset & clinvar variants (supplementary files 1&2) & TANGO target list (supplementary file 2, Lim et al 2020) & uORF data (Ji et al 2015) & ARE data (ARED-plus, Bakheet et al 2018)
# Outputs: Figures used to create Figure 3/4 & supplementary table 1&2
# Dependencies: RColorBrewer, readxl, viridis, stringr, dplyr, svglite, ggplot2, biomaRt, functions in 'Create_sunburst_plot_var_data_VEP.R'
# Notes: R version 4.4.1 (2024-06-14 ucrt)
# ------------------------------------------------------------------------------

source("*/Create_Sunburst_plot_var_data_VEP.R")
library(stringr)
library(ggplot2)
library(svglite)
library(dplyr)
library(RColorBrewer)
library(readxl)
library(viridis)
library(biomaRt)

#analyse Database targetable
NDD_database <- read.csv("*/NDD_database.csv", row.names=1)
clinvar_variants_NDDonly <- read.csv("*/clinvar_variants.csv", row.names=1, stringsAsFactors = FALSE)
Sup2_TANGO_article <- read_excel("*/Sup2.xlsx")

#add inheritance class
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
NDD_database_inheritance_class$Inheritance_class <- unlist(lapply(NDD_database$Inheritance, inheritance_class))
NDD_database <- NDD_database_inheritance_class

#remove not defined clinvar variants
clinvar_variants_NDDonly <- clinvar_variants_NDDonly[!is.na(clinvar_variants_NDDonly$consequence_VEP),]

#create files targetable
NDD_database_targetable <- NDD_database_inheritance_class
clinvar_variants_NDDonly_targetability <- clinvar_variants_NDDonly

#NDDs and ClinVar entries for source_scores of 3 and higher
NDDs_source_score_select <- NDD_database$source_scores >= 3 & !is.na(NDD_database$OMIM_Disease)
NDD_OMIM_source_score_select <- NDD_database$OMIM_Disease[NDDs_source_score_select]
Clinvar_source_score_select <- c()
for (OMIM in clinvar_variants_NDDonly$best_OMIM_ID){
  OMIM_ids <- str_split(OMIM, ";")[[1]]
  ClinVar_select<- any(OMIM_ids %in% NDD_OMIM_source_score_select)
  Clinvar_source_score_select <- c(Clinvar_source_score_select, ClinVar_select)
}
ClinVar_ID_source_score_select <- clinvar_variants_NDDonly$AlleleID[Clinvar_source_score_select]

#SPLICE SITE RESTAURATION
#targetable splice mutations and hotspot splice
distance_closest_canonical <- function(exon_inf, location){
  splice_sites <- as.numeric(unlist(str_split(unlist(str_split(exon_inf, ";")), "-")))
  distance_closest_canonical_site <- abs(splice_sites - location)[order(abs(splice_sites - location))][1]
  return(distance_closest_canonical_site)
}
clinvar_splice_variants <- clinvar_variants_NDDonly[clinvar_variants_NDDonly$consequence_VEP == "intron_variant",]
##50 basepairs from canonical site
location_intron <- mapply(distance_closest_canonical, clinvar_splice_variants$exons_inf, clinvar_splice_variants$POS)
clinvar_splice_variants <- cbind(clinvar_splice_variants, location_intron)
clinvar_splice_variants$location_intron <- abs(as.numeric(clinvar_splice_variants$location_intron))
clinvar_splice_variants$deep_intronic <- clinvar_splice_variants$location_intron >= 50 & clinvar_splice_variants$location_intron < 620060 & clinvar_splice_variants$consequence_VEP == "intron_variant" & !is.na(clinvar_splice_variants$consequence_VEP)
deep_intronic <- clinvar_splice_variants$ID[!is.na(clinvar_splice_variants$ID[clinvar_splice_variants$deep_intronic])]
clinvar_variants_NDDonly_targetability$splice_correct[clinvar_variants_NDDonly_targetability$ID %in% deep_intronic] <- "deep_intronic"
##identify recurrent targetable splice variants
recurrent_targetable_splice <- clinvar_splice_variants[clinvar_splice_variants$deep_intronic & clinvar_splice_variants$num_submit >= 3 & !is.na(clinvar_splice_variants$deep_intronic) & !is.na(clinvar_splice_variants$num_submit),]
clinvar_variants_NDDonly_targetability$splice_correct[clinvar_variants_NDDonly_targetability$ID %in% recurrent_targetable_splice$ID] <- "recurrent_deep_intronic"
#make figures
## deep vs canonical piechard
clinvar_variants_for_sunburst <- clinvar_variants_NDDonly[Clinvar_source_score_select,]
clinvar_variants_for_sunburst <- clinvar_variants_for_sunburst[clinvar_variants_for_sunburst$consequence_VEP == "intron_variant" & !is.na(clinvar_variants_for_sunburst$consequence_VEP),]
clinvar_variants_for_sunburst$splice_correct[clinvar_variants_for_sunburst$ID %in% clinvar_splice_variants$ID[clinvar_splice_variants$deep_intronic]] <- "deep_intronic"
clinvar_variants_for_sunburst$recurrent[clinvar_variants_for_sunburst$ID %in% recurrent_targetable_splice$ID] <- "recurrent_deep_intronic"
var_data_splice <- clinvar_variants_for_sunburst %>% group_by(consequence_VEP, splice_correct, recurrent) %>% summarise(Freq = n(), Perc = (n()/nrow(.))*100)
colnames(var_data_splice) <- c("L1", "L2", "L3", "Frequency", "Percentage")
sunburst_rec_splice <- create_sunburst(Level_dataframe = var_data_splice, levels = c("L1", "L2", "L3"), label = "Levels_count")
#figure recurent splice variants
recurrent_targetable_splice <- recurrent_targetable_splice[recurrent_targetable_splice$AlleleID %in% ClinVar_ID_source_score_select,]
recurrent_targetable_splice <- recurrent_targetable_splice[recurrent_targetable_splice$num_submit >= 3,]
recurrent_targetable_splice <- recurrent_targetable_splice %>% group_by(HGNC_symbol) %>% mutate(FreqSum = sum(num_submit)) %>% arrange(-FreqSum, -num_submit)
recurrent_targetable_splice$HGNC_symbol <- factor(recurrent_targetable_splice$HGNC_symbol, levels = unique(recurrent_targetable_splice$HGNC_symbol))
recurrent_targetable_splice$protein_change <- factor(recurrent_targetable_splice$protein_change, levels = unique(recurrent_targetable_splice$protein_change))
recurrent_targetable_splice$color <- NA
recurrent_targetable_splice$labelx <- NA
recurrent_targetable_splice$labelx[recurrent_targetable_splice$num_submit > 5] <- as.character(recurrent_targetable_splice$protein_change[recurrent_targetable_splice$num_submit > 5])
for (gene in recurrent_targetable_splice$HGNC_symbol){
  n <- length(recurrent_targetable_splice$ID[recurrent_targetable_splice$HGNC_symbol == gene])
  recurrent_targetable_splice$color[recurrent_targetable_splice$HGNC_symbol == gene] <- viridis(n, option = "turbo")
}
recurrent_targetable_splice$labels <- unlist(lapply(recurrent_targetable_splice$variation_name, function(x) str_remove(x, "NM_\\d*.\\d")))
recurrent_splice_plot<- ggplot(data = recurrent_targetable_splice[recurrent_targetable_splice$num_submit >= 5,], aes(x = reorder(labels, -num_submit), y = as.numeric(num_submit))) +
  geom_col(position = "stack") + 
  theme_classic() +
  theme(axis.text.x = element_text(hjust = 1, angle = 45, size = 8),
        strip.placement = 'outside',
        strip.background = element_rect(fill='white', color = 'white', size = 0),
        strip.text.x = ggplot2::element_text(size = 9, angle=90),
        panel.grid.major.x = element_blank(),
        legend.position = "none",
        panel.spacing.x = unit(0, "cm")
  ) +
  scale_y_continuous(expand = c(0,0))

#HI prediction
gnomad_metrics <- read.delim("*/gnomad.v2.1.1.lof_metrics.by_gene.txt")
HI_genes <- gnomad_metrics$gene_id[gnomad_metrics$pLI > 0.9]
NDD_database_targetable$HI <- NA
NDD_database_targetable$HI <- NDD_database_targetable$EnsemblID %in% HI_genes
NDD_database_targetable$AD <- NDD_database_targetable$Inheritance_class == "AD" & !is.na(NDD_database_targetable$Inheritance_class)
AD_NDD <- NDD_database_targetable$OMIM_Disease[NDD_database_targetable$AD]
AD_NDD <- AD_NDD[!is.na(AD_NDD)]
for (AD_disease in AD_NDD){
  clinvar_disease <- clinvar_variants_NDDonly[clinvar_variants_NDDonly$best_OMIM_ID == AD_disease,]
  Perc_LoF <- 0
  LoF_count <- 0
  if (nrow(clinvar_disease) > 1){
    n_LoF <- sum(clinvar_disease$num_submit[(clinvar_disease$consequence_VEP == "splice_acceptor_variant" | 
                                                 clinvar_disease$consequence_VEP == "stop_gained" | 
                                                 clinvar_disease$consequence_VEP == "splice_donor_variant" |
                                                 clinvar_disease$consequence_VEP == "frameshift_variant") &
                                                 !is.na(clinvar_disease$consequence_VEP) & 
                                                 clinvar_disease$NMD != "NMD_escaping_variant"])
    uni_var<- sum(clinvar_disease$num_submit)
    perc_LoF <- n_LoF/uni_var
  }
  NDD_database_targetable$perc_LoF[NDD_database_targetable$OMIM_Disease == AD_disease & !is.na(NDD_database_targetable$Inheritance_class)] <- perc_LoF
}
NDD_database_targetable$perc_LoF_0.2[NDD_database_targetable$perc_LoF > 0.2] <- TRUE
NDD_database_targetable_filt <- NDD_database_targetable[NDDs_source_score_select,]
NDD_database_targetable_filt$HI[!NDD_database_targetable_filt$AD] <- NA
NDD_database_targetable_filt$perc_LoF_0.2[!NDD_database_targetable_filt$AD] <- NA
NDD_database_targetable_filt$HI[!NDD_database_targetable_filt$HI] <- NA
NDD_database_targetable_filt$perc_LoF_0.2[!NDD_database_targetable_filt$HI & !is.na(NDD_database_targetable_filt$HI)] <- NA
NDD_data_HI <- NDD_database_targetable_filt %>% group_by(Inheritance_class, HI, perc_LoF_0.2) %>% summarise(Freq = n(), Perc = (n()/nrow(.))*100)
colnames(NDD_data_HI) <- c("L1", "L2", "L3", "Frequency", "Percentage")
sunburst_NDD_HI <- create_sunburst(Level_dataframe = NDD_data_HI, levels = c("L1", "L2", "L3"), label = "Levels_count")

#coding length for exon skipping
biolist <- as.data.frame(listMarts())
ensembl = useMart("ensembl")
ensemblist <- as.data.frame(listDatasets(ensembl))
ensembl = useDataset("hsapiens_gene_ensembl", mart = ensembl)
filters = listFilters(ensembl)
attributes = listAttributes(ensembl)
cds_length <- getBM(attributes = c("ensembl_transcript_id_version", "cds_length"), filters = c("ensembl_transcript_id_version"), values = NDD_database$MANE_EnsemblID, mart = ensembl)
colnames(cds_length)[1] <- "MANE_EnsemblID"
NDD_database <- left_join(NDD_database, cds_length, by = c("MANE_EnsemblID"))
clinvar_variants_NDDonly <- left_join(clinvar_variants_NDDonly, cds_length, by = c("MANE_EnsemblID"))

#skipping for reading frame restoration
clinvar_variants_NDDonly_targetability$exon_skip <- NA
i_oof_full_exon_SNVs <- c()
for (ID in clinvar_variants_NDDonly$ID[clinvar_variants_NDDonly$CLNVC == "Deletion"]){
  row <- clinvar_variants_NDDonly[clinvar_variants_NDDonly$ID == ID,]
  if (nrow(row)> 1){
    print(row)
    row = row[1,]
  }
  if (row$CLNVC == "Deletion")
    affected_INTRON <- str_extract(row$INTRON, "\\d*-\\d*|\\d*")
    affected_INTRON <- as.numeric(unlist(str_split(affected_INTRON, "-")))
    if (length(affected_INTRON) > 1){
      affected_EXON <- str_extract(row$EXON, "\\d*-\\d*|\\d*")
      affected_EXON <- as.numeric(unlist(str_split(affected_EXON, "-")))
      if (length(affected_EXON) == 1){
        affected_EXON <- c(affected_EXON, affected_EXON)
      }
      whole_exon_del <- affected_INTRON[1] < affected_EXON[1] & affected_INTRON[2] == affected_EXON[2]
      exon_inf <- unlist(str_split(row$exons_inf, ";"))
      #if (all(is.na(exon_inf))){next}
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
          i_oof_full_exon_SNVs <- c(i_oof_full_exon_SNVs,row$best_OMIM_ID)
          phase <- round(round((abs(del_exon_start- del_exon_stop)+1)/3) - (abs(del_exon_stop- del_exon_start)+1)/3, digits = 2)
          exon_df <- rbind(exon_df, abs(as.numeric(exon_df[1,]) - as.numeric(exon_df[2,]))+1)
          tot_len <- row$cds_length
          exon_df <- rbind(exon_df, round(as.numeric(exon_df[3,])/3) == as.numeric(exon_df[3,])/3)
          exon_df <- rbind(exon_df, round(as.numeric(exon_df[3,])/3) - as.numeric(exon_df[3,])/3)
          exon_num <- which(row$Start < as.numeric(exon_df[1,]) & row$Stop > as.numeric(exon_df[2,]))
          if (ncol(exon_df)<4){next}
          if (any(!exon_num %in% 2:(ncol(exon_df)-1))){next}
          pre_exon <- min(exon_num) - 1
          fol_exon <- max(exon_num) + 1
          if(pre_exon %in% 2:(ncol(exon_df)-1)){
            pre_exon <- round(as.numeric(exon_df[5,min(exon_num) - 1]), digits = 2)
            if(phase+pre_exon == 0){
              nt_del<- sum(as.numeric(exon_df[3,c(min(exon_num) - 1, exon_num)]))
              perc_del <- (nt_del/tot_len)*100
              clinvar_variants_NDDonly_targetability$exon_skip[clinvar_variants_NDDonly_targetability$ID == ID] <- perc_del
            }
          }
          if(fol_exon %in% 2:(ncol(exon_df)-1)){
            fol_exon <- round(as.numeric(exon_df[5,max(exon_num) + 1]), digits= 2)
            if(phase+fol_exon == 0){
              nt_del<- sum(as.numeric(exon_df[3,c(exon_num, max(exon_num) + 1)]))
              perc_del <- (nt_del/tot_len)*100
              if (is.na(clinvar_variants_NDDonly_targetability$exon_skip[clinvar_variants_NDDonly_targetability$ID == ID])|clinvar_variants_NDDonly_targetability$exon_skip[clinvar_variants_NDDonly_targetability$ID == ID]>perc_del){
                clinvar_variants_NDDonly_targetability$exon_skip[clinvar_variants_NDDonly_targetability$ID == ID] <- perc_del
              }
            }
          }
        }
      }
    }
  }
}
## targetable frameshift splice mutations
skips_exon_or_intron_retention <- function(clinvar_ID){
  row <- clinvar_variants_NDDonly[clinvar_variants_NDDonly$ID == clinvar_ID,] 
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
        pre_exon <- round(as.numeric(exon_df[5,pre_exon]), digits = 2)
        if(phase+pre_exon == 0){
          nt_del<- sum(as.numeric(exon_df[3,c(pre_exon, closest_exon_index)]))
          perc_del <- (nt_del/tot_len)*100
        }
      }
      if(fol_exon %in% 2:(ncol(exon_df)-1)){
        fol_exon <- round(as.numeric(exon_df[5,fol_exon]), digits= 2)
        if(phase+fol_exon == 0){
          nt_del<- sum(as.numeric(exon_df[3,c(closest_exon_index, fol_exon)]))
          perc_del <- (nt_del/tot_len)*100
        }
      }
      return(perc_del)
    }
  }
  return(-1)
}
canonical_splice <- clinvar_variants_NDDonly_targetability$ID[clinvar_variants_NDDonly_targetability$consequence_VEP == "splice_donor_variant" | clinvar_variants_NDDonly_targetability$consequence_VEP == "splice_acceptor_variant" & !is.na(clinvar_variants_NDDonly_targetability$consequence_VEP)]
is_in_frame <- unlist(mapply(skips_exon_or_intron_retention, canonical_splice))
clinvar_variants_NDDonly_targetability$is_in_frame <- NA
clinvar_variants_NDDonly_targetability$is_in_frame[clinvar_variants_NDDonly_targetability$consequence_VEP == "splice_donor_variant" | clinvar_variants_NDDonly_targetability$consequence_VEP == "splice_acceptor_variant" & !is.na(clinvar_variants_NDDonly_targetability$consequence_VEP)] <- is_in_frame
out_frame_splice <- clinvar_variants_NDDonly_targetability$ID[clinvar_variants_NDDonly_targetability$is_in_frame > 0 & !is.na(clinvar_variants_NDDonly_targetability$is_in_frame)]
out_frame_splice_perc <- clinvar_variants_NDDonly_targetability$ID[clinvar_variants_NDDonly_targetability$is_in_frame > 0 & clinvar_variants_NDDonly_targetability$is_in_frame < 10 & !is.na(clinvar_variants_NDDonly_targetability$is_in_frame)]
out_frame_splice_genes <- clinvar_variants_NDDonly_targetability$best_OMIM_ID[clinvar_variants_NDDonly_targetability$is_in_frame > 0 & clinvar_variants_NDDonly_targetability$is_in_frame < 10& !is.na(clinvar_variants_NDDonly_targetability$is_in_frame)]
#visualize
targetable_variants <- clinvar_variants_NDDonly_targetability$ID[!is.na(clinvar_variants_NDDonly_targetability$exon_skip)& clinvar_variants_NDDonly_targetability$exon_skip > 0 ]
targetable_genes <- clinvar_variants_NDDonly_targetability$best_OMIM_ID[!is.na(clinvar_variants_NDDonly_targetability$exon_skip)& clinvar_variants_NDDonly_targetability$exon_skip > 0& Clinvar_source_score_select]
targetable_variants_10perc <- clinvar_variants_NDDonly_targetability$ID[!is.na(clinvar_variants_NDDonly_targetability$exon_skip) & clinvar_variants_NDDonly_targetability$exon_skip < 10 & clinvar_variants_NDDonly_targetability$exon_skip > 0]
fs_restore_genes <- clinvar_variants_NDDonly_targetability$best_OMIM_ID[!is.na(clinvar_variants_NDDonly_targetability$exon_skip) & clinvar_variants_NDDonly_targetability$exon_skip < 10 & clinvar_variants_NDDonly_targetability$exon_skip > 0]
clinvar_variants_for_sunburst <- clinvar_variants_NDDonly[Clinvar_source_score_select,]
clinvar_variants_for_sunburst$exon_skip[clinvar_variants_for_sunburst$ID %in% targetable_variants | clinvar_variants_for_sunburst$ID %in% out_frame_splice] <- "restore_FS"
clinvar_variants_for_sunburst$exon_skip_perc[clinvar_variants_for_sunburst$ID %in% targetable_variants_10perc | clinvar_variants_for_sunburst$ID %in% out_frame_splice_perc] <- "restore_FS_10perc"
length(clinvar_variants_for_sunburst$ID[clinvar_variants_for_sunburst$best_OMIM_ID == "310200" & clinvar_variants_for_sunburst$exon_skip_perc == "restore_FS_10perc" & !is.na(clinvar_variants_for_sunburst$exon_skip_perc)])
length(clinvar_variants_for_sunburst$ID[clinvar_variants_for_sunburst$best_OMIM_ID == "310200" ])
clinvar_variants_for_sunburst$consequence_VEP[clinvar_variants_for_sunburst$consequence_VEP %in% c("intergenic_variant", "transcript_ablation", "protein_altering_variant", "coding_sequence_variant", "5_prime_UTR_variant", "stop_lost", "synonymous_variant", "3_prime_UTR_variant", "start_lost", "Not defined", "upstream_gene_variant", "inframe_deletion","intron_variant", "inframe_insertion", "downstream_gene_variant")] <- "Others"
clinvar_variants_for_sunburst$consequence_VEP[clinvar_variants_for_sunburst$consequence_VEP %in% c("splice_polypyrimidine_tract_variant", "splice_donor_region_variant", "splice_donor_5th_base_variant", "splice_region_variant", "splice_acceptor_variant", "splice_donor_variant")] <- "splice_region_variant"
var_data_restore_fs <- clinvar_variants_for_sunburst %>% group_by(consequence_VEP, exon_skip, exon_skip_perc) %>% summarise(Freq = n(), Perc = (n()/nrow(.))*100)
colnames(var_data_restore_fs) <- c("L1", "L2", "L3", "Frequency", "Percentage")
sunburst_frame_shift_restore <- create_sunburst(Level_dataframe = var_data_restore_fs, levels = c("L1", "L2", "L3"), label = "Levels_count")


## fs/stop in in-frame exon
in_inframe_exon <- function(clinvar_id){
  del_perc = -1
  row <- clinvar_variants_NDDonly[clinvar_variants_NDDonly$ID == clinvar_id,]
  row <- row[1,]
  exon_inf <- row$exons_inf
  if (is.na(exon_inf)){return(NA)}
  position <- row$POS
  exon_num <- which(unlist(lapply(str_split(unlist(str_split(exon_inf, ";")), "-"), function(x) between(as.numeric(row$POS), as.numeric(x[1]), as.numeric(x[2])) | between(as.numeric(row$POS), as.numeric(x[2]), as.numeric(x[1])))))
  if (length(exon_num) == 0){return(NA)}
  exon <- as.numeric(unlist(str_split(unlist(str_split(exon_inf, ";"))[exon_num], "-")))
  nt_del <- abs(exon[1] - exon[2])
  #if (all(is.na(exon_inf))){next}
  exon_df <- data.frame(sapply(lapply(unlist(str_split(exon_inf, ";")), function(x) str_split(x, "-")), "[[", 1))
  exon_df <- rbind(exon_df, abs(as.numeric(exon_df[1,]) - as.numeric(exon_df[2,]))+1)
  exon_df <- mutate_all(exon_df, function (x) as.numeric(x))
  nt_tot <- row$cds_length
  del_perc <- (nt_del/nt_tot)*100
  if (exon_num == 1 | exon_num == ncol(exon_df)){
    del_perc = -1
  } else if (round((abs(exon[2]- exon[1])+1)/3) == (abs(exon[2]- exon[1])+1)/3){
    return(del_perc)
  } else{del_perc = -1}
  return(del_perc)
}
trunc_ID <- clinvar_variants_NDDonly$ID[clinvar_variants_NDDonly$consequence_VEP == "stop_gained" | clinvar_variants_NDDonly$consequence_VEP == "frameshift_variant" & !is.na(clinvar_variants_NDDonly$consequence_VEP) & !clinvar_variants_NDDonly$NMD == "NMD_escaping_variant"]
skippable_trunc_ID <- lapply(trunc_ID, in_inframe_exon)
clinvar_variants_NDDonly_targetability$IF_exon_skip <- NA
clinvar_variants_NDDonly_targetability$IF_exon_skip[clinvar_variants_NDDonly$consequence_VEP == "stop_gained" | clinvar_variants_NDDonly$consequence_VEP == "frameshift_variant" & !is.na(clinvar_variants_NDDonly$consequence_VEP) & !clinvar_variants_NDDonly$NMD == "NMD_escaping_variant"] <- unlist(skippable_trunc_ID)
inframe_targetable_variants <- clinvar_variants_NDDonly_targetability$ID[clinvar_variants_NDDonly_targetability$IF_exon_skip > 0 & !is.na(clinvar_variants_NDDonly_targetability$IF_exon_skip)]
inframe_targetable_genes <- clinvar_variants_NDDonly_targetability$best_OMIM_ID[clinvar_variants_NDDonly_targetability$IF_exon_skip > 0 & !is.na(clinvar_variants_NDDonly_targetability$IF_exon_skip) & Clinvar_source_score_select]
inframe_targetable_variants_10_perc <- clinvar_variants_NDDonly_targetability$ID[clinvar_variants_NDDonly_targetability$IF_exon_skip > 0 & clinvar_variants_NDDonly_targetability$IF_exon_skip < 10 & !is.na(clinvar_variants_NDDonly_targetability$IF_exon_skip)]
inframe_targetable_genes_10_perc <- clinvar_variants_NDDonly_targetability$best_OMIM_ID[clinvar_variants_NDDonly_targetability$IF_exon_skip > 0 & clinvar_variants_NDDonly_targetability$IF_exon_skip < 10 & !is.na(clinvar_variants_NDDonly_targetability$IF_exon_skip) & Clinvar_source_score_select]
#create sunburst plot
clinvar_variants_for_sunburst <- clinvar_variants_NDDonly[Clinvar_source_score_select,]
clinvar_variants_for_sunburst$consequence_VEP[clinvar_variants_for_sunburst$consequence_VEP %in% c("intergenic_variant", "transcript_ablation", "protein_altering_variant", "coding_sequence_variant", "5_prime_UTR_variant", "stop_lost", "synonymous_variant", "3_prime_UTR_variant", "start_lost", "Not defined", "upstream_gene_variant", "inframe_deletion","intron_variant", "inframe_insertion", "downstream_gene_variant")] <- "Others"
clinvar_variants_for_sunburst$consequence_VEP[clinvar_variants_for_sunburst$consequence_VEP %in% c("splice_polypyrimidine_tract_variant", "splice_donor_region_variant", "splice_donor_5th_base_variant", "splice_region_variant", "splice_acceptor_variant", "splice_donor_variant")] <- "splice_region_variant"
clinvar_variants_for_sunburst$exon_skip[clinvar_variants_for_sunburst$ID %in% inframe_targetable_variants] <- "skip_inframe"
clinvar_variants_for_sunburst$exon_skip_perc[clinvar_variants_for_sunburst$ID %in% inframe_targetable_variants_10_perc] <- "skip_inframe_10perc"
var_data_inframe_skip <- clinvar_variants_for_sunburst %>% group_by(consequence_VEP, exon_skip, exon_skip_perc) %>% summarise(Freq = n(), Perc = (n()/nrow(.))*100)
colnames(var_data_inframe_skip) <- c("L1", "L2", "L3", "Frequency", "Percentage")
sunburst_inframe_skip <- create_sunburst(Level_dataframe = var_data_inframe_skip, levels = c("L1", "L2", "L3"), label = "Levels_count")

#TANGO NDDs
amenable_genes_TANGO <- unique(intersect(Sup2_TANGO_article$Gene, NDD_database$HGNC_symbol))
amenable_genes_TANGO_AD <- intersect(amenable_genes_TANGO, NDD_database$HGNC_symbol)
NDD_database_targetable$TANGO_all <- NDD_database_targetable$HGNC_symbol %in% amenable_genes_TANGO
NDD_database_targetable$Inheritance_class[NDD_database_targetable$Inheritance_class %in% c("Imprinting", "Mitochondrial", "Unknown")] <- "Other"
#visualize TANGO
NDD_database_targetable$TANGO_all[!NDD_database_targetable$TANGO_all] <- NA
NDD_data_TANGO <- NDD_database_targetable[NDDs_source_score_select,]
NDD_data_TANGO <- NDD_data_TANGO %>% group_by(Inheritance_class, TANGO_all) %>% summarise(Freq = n(), Perc = (n()/nrow(.))*100)
colnames(NDD_data_TANGO) <- c("L1", "L2", "Frequency", "Percentage")
sunburst_NDD_TANGO <- create_sunburst(Level_dataframe = NDD_data_TANGO, levels = c("L1", "L2"), label = "Levels_count")

#uORF
uORF <- read_excel("*/uORF.xlsx", col_names = FALSE, sheet = "uORF")
ouORF <- read_excel("*/uORF.xlsx", col_names = FALSE, sheet = "overlapping.uORF")
colnames(uORF)[1:3] <- c("HGNC_symbol", "loc", "chr")
colnames(ouORF)[1:3] <- c("HGNC_symbol", "loc", "chr")
genes_uORF <- unique(intersect(NDD_database$HGNC_symbol, uORF$HGNC_symbol), intersect(NDD_database$HGNC_symbol, ouORF$HGNC_symbol))
NDD_database_targetable$uORF <- NDD_database_targetable$HGNC_symbol %in% genes_uORF
#visualize uORF
NDD_database_targetable$uORF[!NDD_database_targetable$uORF] <- NA
NDD_data_uORF <- NDD_database_targetable[NDDs_source_score_select,]
NDD_data_uORF <- NDD_data_uORF %>% group_by(Inheritance_class, uORF) %>% summarise(Freq = n(), Perc = (n()/nrow(.))*100)
colnames(NDD_data_uORF) <- c("L1", "L2", "Frequency", "Percentage")
sunburst_NDD_uORF <- create_sunburst(Level_dataframe = NDD_data_uORF, levels = c("L1", "L2"), label = "Levels_count")

#ARED
Complete_3_UTR_Data <- read_excel("*/Complete 3 UTR Data.xls")
genes_ARE <- unique(intersect(NDD_database$HGNC_symbol, Complete_3_UTR_Data$GeneName))
NDD_database_targetable$ARE_3UTR <- NDD_database_targetable$HGNC_symbol %in% genes_ARE
#visualize ARED
NDD_database_targetable$ARE_3UTR[!NDD_database_targetable$ARE_3UTR] <- NA
NDD_data_ARE_3UTR <- NDD_database_targetable[NDDs_source_score_select,]
NDD_data_ARE_3UTR <- NDD_data_ARE_3UTR %>% group_by(Inheritance_class, ARE_3UTR) %>% summarise(Freq = n(), Perc = (n()/nrow(.))*100)
colnames(NDD_data_ARE_3UTR) <- c("L1", "L2", "Frequency", "Percentage")
sunburst_NDD_ARE_3UTR <- create_sunburst(Level_dataframe = NDD_data_ARE_3UTR, levels = c("L1", "L2"), label = "Levels_count")

#ORFs combined
NDD_data_ARE_uORF <- NDD_database_targetable[NDDs_source_score_select,]
NDD_data_ARE_uORF <- NDD_data_ARE_uORF %>% group_by(Inheritance_class, ARE_3UTR, uORF) %>% summarise(Freq = n(), Perc = (n()/nrow(.))*100)
colnames(NDD_data_ARE_uORF) <- c("L1", "L2","L3", "Frequency", "Percentage")
sunburst_NDD_data_ARE_uORF <- create_sunburst(Level_dataframe = NDD_data_ARE_uORF, levels = c("L1", "L2", "L3"), label = "Levels_count")

#lncRNAs
NDD_data_lncRNAs <- NDD_database_targetable
NDD_data_lncRNAs$lncRNA_in_TAD <- NA
NDD_data_lncRNAs$lncRNA_in_TAD[nchar(NDD_data_lncRNAs$lncRNAs) > 0] <-  "Yes"
NDD_data_lncRNAs$lncRNA_in_TAD[nchar(NDD_data_lncRNAs$lncRNAs) == 0] <-  NA
NDD_data_lncRNAs$is_conf_lncRNAs <- NA
NDD_data_lncRNAs$is_conf_lncRNAs[nchar(NDD_data_lncRNAs$conf_lncRNAs) > 0 & !is.na(NDD_data_lncRNAs$conf_lncRNAs)] <-  "Yes"
NDD_data_lncRNAs$is_conf_lncRNAs[is.na(NDD_data_lncRNAs$conf_lncRNAs)] <-  NA
NDD_data_lncRNAs <- NDD_data_lncRNAs[NDDs_source_score_select,]
NDD_data_lncRNAs <- NDD_data_lncRNAs %>% group_by(Inheritance_class, lncRNA_in_TAD, is_conf_lncRNAs) %>% summarise(Freq = n(), Perc = (n()/nrow(.))*100)
colnames(NDD_data_lncRNAs) <- c("L1", "L2","L3","Frequency", "Percentage")
sunburst_NDD_data_lncRNAs <- create_sunburst(Level_dataframe = NDD_data_lncRNAs, levels = c("L1", "L2", "L3"), label = "Levels_count")

#add targetable variants to NDD database
NDD_splice_variants <- clinvar_variants_NDDonly_targetability$best_OMIM_ID[(clinvar_variants_NDDonly_targetability$splice_correct == "recurrent_deep_intronic" | clinvar_variants_NDDonly_targetability$splice_correct == "deep_intronic") & !is.na(clinvar_variants_NDDonly_targetability$splice_correct)]
NDD_splice_variants <- unlist(lapply(NDD_splice_variants, function(x) str_split(x, ";")))
NDD_splice_variants <- NDD_splice_variants[!is.na(NDD_splice_variants)]
NDD_database_targetable$Deep_splice_variants <- NDD_database_targetable$OMIM_Disease %in% NDD_splice_variants
NDD_rec_splice_variants <- clinvar_variants_NDDonly_targetability$best_OMIM_ID[clinvar_variants_NDDonly_targetability$splice_correct == "recurrent_deep_intronic" & !is.na(clinvar_variants_NDDonly_targetability$splice_correct)]
NDD_rec_splice_variants <- unlist(lapply(NDD_rec_splice_variants, function(x) str_split(x, ";")))
NDD_rec_splice_variants <- NDD_rec_splice_variants[!is.na(NDD_rec_splice_variants)]
NDD_database_targetable$rec_deep_splice_variants <- NDD_database_targetable$OMIM_Disease %in% NDD_rec_splice_variants

inframe_targetable_genes <- unique(clinvar_variants_NDDonly_targetability$best_OMIM_ID[clinvar_variants_NDDonly_targetability$IF_exon_skip > 0 & clinvar_variants_NDDonly_targetability$IF_exon_skip < 10 & !is.na(clinvar_variants_NDDonly_targetability$IF_exon_skip)])
inframe_targetable_genes <- unique(unlist(lapply(inframe_targetable_genes, function(x) str_split(x, ";"))))
inframe_targetable_genes <- inframe_targetable_genes[!is.na(inframe_targetable_genes)]
NDD_database_targetable$inframe_targetable <- NDD_database_targetable$OMIM_Disease %in% inframe_targetable_genes

fs_restore_genes_indel <- unique(clinvar_variants_NDDonly_targetability$best_OMIM_ID[!is.na(clinvar_variants_NDDonly_targetability$exon_skip) & clinvar_variants_NDDonly_targetability$exon_skip > 0 & clinvar_variants_NDDonly_targetability$exon_skip < 10])
fs_restore_genes_splice <- unique(clinvar_variants_NDDonly_targetability$best_OMIM_ID[clinvar_variants_NDDonly_targetability$is_in_frame > 0 & clinvar_variants_NDDonly_targetability$is_in_frame < 10 & !is.na(clinvar_variants_NDDonly_targetability$is_in_frame)])
fs_restore_genes <- unique(c(fs_restore_genes_indel, fs_restore_genes_splice))
fs_restore_genes <- unlist(lapply(fs_restore_genes, function(x) str_split(x, ";")))
fs_restore_genes <- fs_restore_genes[!is.na(fs_restore_genes)]
fs_restore_splice <- unique(clinvar_variants_NDDonly_targetability$ID[clinvar_variants_NDDonly_targetability$is_in_frame > 0 & clinvar_variants_NDDonly_targetability$is_in_frame < 10 & !is.na(clinvar_variants_NDDonly_targetability$is_in_frame)])
fs_restore_indel <- unique(clinvar_variants_NDDonly_targetability$ID[!is.na(clinvar_variants_NDDonly_targetability$exon_skip) & clinvar_variants_NDDonly_targetability$exon_skip > 0 & clinvar_variants_NDDonly_targetability$exon_skip < 10])
fs_restore_vars <- unique(c(fs_restore_splice, fs_restore_indel))
clinvar_variants_NDDonly_targetability$exon_skip <- clinvar_variants_NDDonly_targetability$ID %in% fs_restore_vars
NDD_database_targetable$fs_restore <- NDD_database_targetable$OMIM_Disease %in% fs_restore_genes

#add TANGO, HI UTRs, lncRNAs
NDD_database_targetable$HI <- NDD_database_targetable$EnsemblID %in% HI_genes
NDD_database_targetable$TANGO_all <- NDD_database_targetable$HGNC_symbol %in% amenable_genes_TANGO
NDD_database_targetable$uORF <- NDD_database_targetable$HGNC_symbol %in% genes_uORF
NDD_database_targetable$ARE_3UTR <- NDD_database_targetable$HGNC_symbol %in% genes_ARE
NDD_database_targetable$lncRNA <- !is.na(NDD_database_targetable$conf_lncRNAs)

clinvar_variants_NDDonly_targetability$'Source_score_>3'<- clinvar_variants_NDDonly_targetability$AlleleID %in% ClinVar_ID_source_score_select
write.csv(NDD_database_targetable, "*/NDD_database_targetable.csv", row.names = TRUE)
write.csv(clinvar_variants_NDDonly_targetability, "*/clinvar_variants_targetable.csv", row.names = TRUE)

