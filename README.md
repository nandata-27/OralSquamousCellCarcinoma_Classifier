# OralSquamousCellCarcinoma_Classifier

# Gene Expression Analysis and Classification in Oral Squamous Cell Carcinoma (OSCC)

This project investigates molecular differences in Oral Squamous Cell Carcinoma (OSCC) using RNA-seq gene expression data. It combines differential expression analysis, machine learning, and statistical testing to identify biomarkers and classify samples across disease stages.

## Project Overview

- Used a publicly available RNA-seq dataset (GSE227919) containing samples from healthy controls, premalignant lesions (PML), and OSCC cases.
- Performed preprocessing, quality control, normalization, and outlier filtering.
- Conducted unsupervised and supervised machine learning to identify top-ranked genes that distinguish between disease states.
- Built classifiers to predict disease class and evaluated their performance.

## Tools and Methods

*R (DESeq2)*
- Differential gene expression analysis  
- Volcano plots and heatmaps  

*MATLAB*
- Data integration, normalization, and log-transformation  
- Feature selection (fsrftest, ReliefF)  
- Dimensionality reduction (t-SNE, hierarchical clustering)  
- Classification using Linear Discriminant and Bagged Trees (accuracy: 73.3%)

*SPSS*
- Statistical tests on selected genes (Shapiro–Wilk, Mann–Whitney U, Kruskal–Wallis)  
- Verified group-level significance of key biomarkers

# Key Results

- Genes such as HDGF, FBLIM1, and LARGE1 showed consistent upregulation in OSCC and PML samples.
- 7 of the top 10 genes selected via ReliefF were significantly different across groups.
- Classifiers built in MATLAB using selected features achieved ~73% weighted F1-score in cross-validation.

## Notes
This analysis was conducted as part of a graduate project for PATH 828 (Queen’s University, 2025). Data were sourced from the Gene Expression Omnibus (GSE227919).
