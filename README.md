# NDD_targetability
Code used in "systematic analysis of genetic and phenotypic characteristics reveals antisense oligonucleotide therapy potential for one-third of neurodevelopmental disorders" (Wijnant et al. 2025)

This repository contains all code to reanalyse the targetability of the NDD dataset published in (Wijnant et al. 2025) https://www.biorxiv.org/content/10.1101/2025.03.20.644369v1
it contains:
- Analyse NDD_database_general_VEP.R: General analysis of the NDD dataset generated for the systematic evaluation of targetability of AONs with code to re-generate figures
- find_targetable_lncRNAs_clinvar2.R: Identication of targetable lncRNAs (NDD genes with correlated lncRNAs in their TAD domains)
- Analyse NDD_database_targetable_VEP.R: Targetability analysis of the NDD dataset generated for the systematic evaluation of targetability of AONs with code to re-generate figures & re-assess the targetability of NDDs/ClinVar variants of 7 AON strategies
- Create_Upset_plot.R: Create upset plots on patient and NDD level and to calculate how many patients one can treat with a single AON
- Targetable_phenotypes.R: Create upset plots of targetable phenotypic features

All data generated or analyzed during this study are included in the published article [and its supplementary information files]. All original datasets can be retrieved from the public resources listed in the work. Datasets from different sources were integrated and deduplicated using a combination of code and manual curation. Thus, the standalone code is not sufficient to regenerate created datasets, this code is therefore not included in this repository.

Please cite: [link & citation available upon publication] or https://www.biorxiv.org/content/10.1101/2025.03.20.644369v1
