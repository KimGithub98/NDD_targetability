# ------------------------------------------------------------------------------
# Title: clean DDG2P data
# Author: KN Wijnant
# Date: 14-11-2023
# Purpose: To deduplicate and prepare for integration with other NDD databases 
# Inputs: NDD_dataset (manually assessed)
# Outputs: NDD_dataset with source scores
# Dependencies: UpSetR
# Notes: R version 4.4.1 (2024-06-14 ucrt)
# ------------------------------------------------------------------------------

library(UpSetR)

options(scipen = 999)

NDD_database <- read.csv("NDD_dataset_manually_assessed.csv", row.names=1) #manually assessed NDD dataset file

#visualize overlap databases
listx <- list(
  sysNDD_Def = rownames(NDD_database)[grepl("Definitive", NDD_database$Confidence_sysNDD) & grepl("sysNDD", NDD_database$Sources)], 
  sysNDD_Lim = rownames(NDD_database)[grepl("Limited", NDD_database$Confidence_sysNDD) & grepl("sysNDD", NDD_database$Sources)], 
  DDG2P_def = rownames(NDD_database)[grepl("definitive", NDD_database$Confidence_DDG2P) & grepl("DDG2P", NDD_database$Sources)],
  DDG2P_lim = rownames(NDD_database)[grepl("limited", NDD_database$Confidence_DDG2P) & grepl("DDG2P", NDD_database$Sources)],
  DDG2P_other = rownames(NDD_database)[grepl("moderate", NDD_database$Confidence_DDG2P) | grepl("strong", NDD_database$Confidence_DDG2P)| grepl("both", NDD_database$Confidence_DDG2P) & grepl("DDG2P", NDD_database$Sources)],
  HPOmim_ID =  rownames(NDD_database)[grepl("HPOmim_ID", NDD_database$Sources)],
  HPOmim_NDD =  rownames(NDD_database)[grepl("HPOmim_NDD", NDD_database$Sources)],
  HPOmim_EPI =  rownames(NDD_database)[grepl("HPOmim_EPI", NDD_database$Sources)],
  IDpanel =  rownames(NDD_database)[grepl("IDpanel", NDD_database$Sources)],
  EPIpanel =  rownames(NDD_database)[grepl("EPIpanel", NDD_database$Sources)]
)

#figure with confidence categories
boolean_df <- t(data.frame(lapply(rownames(NDD_database), function (y) unlist(lapply(listx, function(x) any(grepl(paste("^", y, "$", sep = ""), x)))))))
rownames(boolean_df) <- rownames(NDD_database)
bit_df <- data.frame(1*boolean_df)
upsetplot_all <- upset(bit_df, sets = colnames(bit_df), sets.bar.color = "#56B4E9", order.by = "freq", keep.order = TRUE)
upsetplot_all_order_score <- upset(bit_df, sets = colnames(bit_df), sets.bar.color = "#56B4E9", order.by = "degree", keep.order = TRUE, nintersects = 100, decreasing = TRUE)

#figurer with sysNDD, DDG2P total
bit_df$sysNDD <- bit_df$sysNDD_Def + bit_df$sysNDD_Lim
bit_df$DDG2P <- bit_df$DDG2P_def + bit_df$DDG2P_lim + bit_df$DDG2P_other
bit_df[bit_df == 2] <- 1
upsetplot_sum <- upset(bit_df, sets = c("sysNDD","DDG2P","HPOmim_ID","HPOmim_NDD","HPOmim_EPI", "IDpanel", "EPIpanel"), sets.bar.color = "#56B4E9", order.by = "freq", empty.intersections = "on", keep.order = TRUE)
upsetplot_sum_all_order_score <- upset(bit_df, sets = c("sysNDD","DDG2P","HPOmim_ID","HPOmim_NDD","HPOmim_EPI", "IDpanel", "EPIpanel"), sets.bar.color = "#56B4E9", order.by = "degree", keep.order = TRUE, nintersects = 60, decreasing = TRUE)


#add scores based on Sources checklist file
bit_df <- data.frame(1*boolean_df)
bit_df$sysNDD_Def <- bit_df$sysNDD_Def*2
bit_df$DDG2P_def <- bit_df$DDG2P_def*2
bit_df$IDpanel <- bit_df$IDpanel*2
bit_df$EPIpanel <- bit_df$EPIpanel*2

scores <- rowSums(bit_df)

NDD_database$source_scores <- scores

write.csv(NDD_dataset, file = "NDD_dataset_source_scores.csv", row.names = TRUE)
