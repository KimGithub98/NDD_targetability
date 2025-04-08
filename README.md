# NDD targetability
Code used in "systematic analysis of genetic and phenotypic characteristics reveals antisense oligonucleotide therapy potential for one-third of neurodevelopmental disorders" (Wijnant et al. 2025)

This repository contains all code to reanalyse the targetability of the NDD dataset published in (Wijnant et al. 2025) https://www.biorxiv.org/content/10.1101/2025.03.20.644369v1
it contains in order:
1. find_targetable_lncRNAs_clinvar2.R: Identication of targetable lncRNAs (NDD genes with correlated lncRNAs in their TAD domains)
2. Analyse NDD_database_general_VEP.R: General analysis of the NDD dataset generated for the systematic evaluation of targetability of AONs with code to re-generate figures
3. Analyse NDD_database_targetable_VEP.R: Targetability analysis of the NDD dataset generated for the systematic evaluation of targetability of AONs with code to re-generate figures & re-assess the targetability of NDDs/ClinVar variants of 7 AON strategies
4. Targetable_phenotypes.R: Create upset plots of targetable phenotypic features
5. Create_Upset_plot.R: Create upset plots on patient and NDD level and to calculate how many patients one can treat with a single AON

The NDD dataset analysed in this work (as used in find_targetable_lncRNAs_clinvar2.R is listed in the files folder). This NDD dataset is a combined dataset from different sources that were integrated and deduplicated using a combination of code and manual curation. Thus, the standalone code is not sufficient to regenerate created datasets, this code is therefore not included in this repository. All original datasets can be retrieved from the public resources listed in (Wijnant et al. 2025) https://www.biorxiv.org/content/10.1101/2025.03.20.644369v1. 

Please cite: [link & citation available upon publication] or https://www.biorxiv.org/content/10.1101/2025.03.20.644369v1
