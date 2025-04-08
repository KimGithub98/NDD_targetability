# ------------------------------------------------------------------------------
# Title: targetable lncRNA identification
# Author: KN Wijnant
# Date: 20-1-25
# Purpose: For the identication of targetable lncRNAs (NDD genes with correlated lncRNAs in their TAD domains)
# Inputs: NDD dataset & TADMap TAD regions & libd_stemcell_timecourse_rseGene_n157.rda
# Outputs: Figures used to create Figure 3/4 & supplementary table 1&2
# Dependencies: RColorBrewer, readxl, viridis, stringr, dplyr, svglite, ggplot2, biomaRt, functions in 'Create_sunburst_plot_var_data_VEP.R'
# Notes: R version 4.4.1 (2024-06-14 ucrt)
# ------------------------------------------------------------------------------

library(SummarizedExperiment)
library(biomaRt)
library(dplyr)
library(stringr)
library(tidyr)
library(data.table)

#NDD database
NDD_database <- read.csv("*/NDD_database.csv", row.names=1)
NDD_database <- NDD_database[,-c(28:30)]
TADmap_bed <- read.table("*/TADMap TAD regions.txt", quote="\"", comment.char="")
colnames(TADmap_bed) = c("chrom", "TAD_start", "TAD_end")

biolist <- as.data.frame(listMarts())
ensembl = useMart("ensembl")
ensemblist <- as.data.frame(listDatasets(ensembl))
ensembl = useDataset("hsapiens_gene_ensembl", mart = ensembl)
filters = listFilters(ensembl)
attributes = listAttributes(ensembl)

t2g <- getBM(c("ensembl_gene_id", "ensembl_gene_id_version","external_gene_name", "chromosome_name","transcript_start","transcript_end", "gene_biotype"), mart = ensembl)

#Gtex data
GTEx_Analysis_2017.06.05_v8_RNASeQCv1.1.9_gene_median_tpm.gct <- read.delim("C:/Users/Z978196/OneDrive - Radboudumc/z978196/Documenten/Overview (O1)/genes_lists/raw/lncRNAs/GTEx_Analysis_2017-06-05_v8_RNASeQCv1.1.9_gene_median_tpm.gct", header=FALSE, comment.char="#")
GTEx_TPM <- GTEx_Analysis_2017.06.05_v8_RNASeQCv1.1.9_gene_median_tpm.gct
GTEx_TPM <- GTEx_TPM[-1,]
colnames(GTEx_TPM) <- GTEx_TPM[1,]
GTEx_TPM <- GTEx_TPM[-1,]
GTEx_name <- GTEx_TPM$Name
GTEx_TPM$Name <- unlist(lapply(GTEx_TPM$Name, function (x) str_extract(x, "\\D*\\d*")))


#lncRNAs within TAD
genes_ENS = getBM(attributes = c("hgnc_symbol","ensembl_gene_id", "chromosome_name", "start_position", "end_position","gene_biotype"), mart = ensembl)
lncRNAs_ENS = getBM(attributes = c("hgnc_symbol","ensembl_gene_id", "chromosome_name", "start_position", "end_position","gene_biotype"), filters = c("biotype"), values = list(biotype="LncRNA"), mart = ensembl)
genes_ENS <- genes_ENS[nchar(lncRNAs_ENS$chromosome_name) < 3,]
lncRNAs_ENS <- lncRNAs_ENS[nchar(lncRNAs_ENS$chromosome_name) < 3,]
colnames(genes_ENS)[2] <- "EnsemblID"

for (gene in unique(NDD_database$HGNC_symbol)){
  row <- NDD_database[NDD_database$HGNC_symbol == gene,]
  row <- row[1,]
  gene_range <- unique(c(row$start_position, row$end_position))
  TADmap_bed_chrom <- TADmap_bed[str_extract(TADmap_bed$chrom, "\\d+|X") == row$chromosome_name,]
  TAD_start_pos <- which(data.table::between(row$start_position, TADmap_bed_chrom$TAD_start, TADmap_bed_chrom$TAD_end))
  TAD_end_pos <- which(data.table::between(row$end_position, TADmap_bed_chrom$TAD_start, TADmap_bed_chrom$TAD_end))
  if (length(TAD_start_pos) == 0){
    if (length(TAD_end_pos) == 0){
      search_range = c(-1,0)
    }else{
      TAD <- TADmap_bed_chrom[TAD_end_pos,]
      search_range <- c(TAD$TAD_start, TAD$TAD_end)
    }
  } else{
    if (length(TAD_end_pos) == 0){search_range = c(-1,0)}
    else if (TAD_start_pos != TAD_end_pos){
      TAD <- TADmap_bed_chrom[TAD_start_pos,]
      TAD2 <- TADmap_bed_chrom[TAD_end_pos,]
      search_range <- range(c(TAD$TAD_start, TAD$TAD_end, TAD2$TAD_start, TAD2$TAD_end))
      }
    TAD <- TADmap_bed_chrom[TAD_start_pos,]
    search_range <- c(TAD$TAD_start, TAD$TAD_end)
  }
  close_lncRNAs <- lncRNAs_ENS$ensembl_gene_id[(dplyr::between(lncRNAs_ENS$start_position, search_range[1], search_range[2]) | dplyr::between(lncRNAs_ENS$end_position, search_range[1], search_range[2])) & lncRNAs_ENS$chromosome_name == unique(row$chromosome_name)]
  test_close_lncRNAs <- unique(close_lncRNAs)[nchar(unique(close_lncRNAs))>0]
  NDD_database$lncRNAs[NDD_database$HGNC_symbol == gene] <- paste(test_close_lncRNAs, collapse = ";")
}

cor_lncRNAs_gtex <- function(gene_name, lncRNA_name){
  #x <- GTEx_TPM[GTEx_TPM$Name == gene_name,][1,-c(1:2)]
  #y <- GTEx_TPM[GTEx_TPM$Name == lncRNA_name,][1,-c(1:2)]
  x <- GTEx_TPM[GTEx_TPM$Name == gene_name,][1,10:22] #only brain
  y <- GTEx_TPM[GTEx_TPM$Name == lncRNA_name,][1,10:22] #only brain
  if (all(is.na(y))| all(is.na(x))){
    p_value = NA
  }else{
    cor_inf <- cor.test(as.numeric(x), as.numeric(y))
    p_value <- cor_inf$p.value
  }
  return(p_value)
}

#test correlation and adjust according to bonferroni
for (gene in unique(NDD_database$HGNC_symbol)){
  row <- NDD_database[NDD_database$HGNC_symbol == gene,]
  test_close_lncRNAs <- unlist(str_split(unique(row$lncRNAs), ";"))
  p_values_GTX <- lapply(test_close_lncRNAs, function (x) cor_lncRNAs_gtex(row$EnsemblID, x))
  p_values <- cbind(p_values_GTX)
  rownames(p_values) <- test_close_lncRNAs
  NDD_database[NDD_database$HGNC_symbol == gene, c(colnames(p_values))] <- unlist(apply(p_values, 2, function (x) paste(x, collapse = ";")))
}

#p.adjust from bonferroni to BH
NDD_database_p_seperated <- separate_rows(NDD_database,c("lncRNAs", "p_values_GTX"),sep = ";")
p_values_GTX <- NDD_database_p_seperated$p_values_GTX[NDD_database_p_seperated$p_values_GTX != "NA" & !is.na(NDD_database_p_seperated$p_values_GTX)]
p_values_GTX_new <- p.adjust(NDD_database_p_seperated$p_values_GTX[NDD_database_p_seperated$p_values_GTX != "NA" & !is.na(NDD_database_p_seperated$p_values_GTX)], method = "BH")
NDD_database_p_seperated$p_adj_GTX <- NA
NDD_database_p_seperated$p_adj_GTX[NDD_database_p_seperated$p_values_GTX != "NA" & !is.na(NDD_database_p_seperated$p_values_GTX)] <- p_values_GTX_new
NDD_database_p_seperated$conf_lncRNAs <- NA
NDD_database_p_seperated$trend_lncRNAs <- NA
confident <- as.numeric(NDD_database_p_seperated$p_adj_GTX) < 0.05 & !is.na(as.numeric(NDD_database_p_seperated$p_adj_GTX))
trend <- as.numeric(NDD_database_p_seperated$p_adj_GTX) < 0.1 & !is.na(as.numeric(NDD_database_p_seperated$p_adj_GTX))
NDD_database_p_seperated$conf_lncRNAs[confident] <- NDD_database_p_seperated$lncRNAs[confident]
NDD_database_p_seperated$trend_lncRNAs[trend] <- NDD_database_p_seperated$lncRNAs[trend]

NDD_database_lncRNAs_new <- NDD_database_p_seperated %>% ###merge conf_lncRNAs
  group_by(!!!syms(colnames(NDD_database_p_seperated)[c(1:23, 30)])) %>% 
  summarise(lncRNAs = paste(lncRNAs, collapse = ";"), p_values_GTX = paste(p_values_GTX, collapse = ";"), p_adj_GTX=paste(p_adj_GTX, collapse = ";"), conf_lncRNAs = paste(conf_lncRNAs[!is.na(conf_lncRNAs)], collapse = ";"), trend_lncRNAs = paste(trend_lncRNAs[!is.na(trend_lncRNAs)], collapse = ";"))
NDD_database_lncRNAs_new$conf_lncRNAs[nchar(NDD_database_lncRNAs_new$conf_lncRNAs) == 0] <- NA
NDD_database_lncRNAs_new$trend_lncRNAs[nchar(NDD_database_lncRNAs_new$trend_lncRNAs) == 0] <- NA
