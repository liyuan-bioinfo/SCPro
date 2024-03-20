

- Datasets and Methods of uSCPro
- update: 20240319
- author: Yuan
  
```
# ------------------------------------------------
##     Bioinformastic analysis of datasets
# ------------------------------------------------
```

dataset_01: 3 cell types, each has 4 replications, including Acinar, Lymph and Tumor. 
- Quantified Analysis: PCA; Corr; ANOVA; 
- Stat Analysis: ANOVA; 
- Functional Analysis: GO-BP; 

dataset_02: 2 cell types, each has 4 replications, including PCC and Stroma. 
- Quantified Analysis: PCA; Corr; 
- Stat Analysis: Student's test;

## dataset_03: 6 cell types, each has 3 replications, including [Acinar, PanIN and PDAC]; CAF, [IT and LN]. 
- Quantified Analysis: PCA; 
- Stat Analysis: ANOVA[ct3]; Student's test[ct2];
- Functional Analysis: GO-BP; 

## dataset_04: 14 cell types / 4 cell lineages, each has 3~4 replications, total 69 samples, including PCC, CAF, LYM and MYE. 
### 4 cell lineages
- Quantified Analysis: PCA; Pearson Corr.
- Stat Analysis: LIMMA[ct4]
- Functional Analysis: GO-BP; 

### 14 celltypes
- Quantified Analysis: Plasma Membrane Annotation; Molecular Function Annotation; 
- Stat Analysis: LIMMA[ct14]; 
- Futher Analysis: score for identify novel celltypes

## dataset_05: 8 cell types, each has 3~4 replications, total 69 samples, including PCC, CAF, LYM and MYE. 
- Quantified Analysis: PCA; 
- Stat Analysis: Student's test;
- Functional Analysis: GSEA; 

## dataset_06, dataset_04 and dataset_05
- Deconvolution Analysis: Tangram; 


```
# ------------------------------------------------
##             Quantification analysis
# ------------------------------------------------
```
## Basic workflow
- Filter: at least two valid in one group
- Imputation: normal distribution[downshift=1.8;width=0.3]
- PCA: not filter / quantified protein
- corr: using data after impute. [Note]. When calculating Pearson correlation, the comparison of a variable with itself was replaced with the maximum value of comparisons between pair-wise groups.
- pheatmap: z-score of log2LFQ intensity. [Note]. All heatmap visualizations were performed on log2-normalized and scaled data with a standard deviation of 1 and a mean of 0. Values exceeding the 99th percentile or falling below the 1st percentile were replaced with the values of the 99th and 1st percentiles, respectively, prior to visualization.

```
# ------------------------------------------------
##                  Stat analysis
# ------------------------------------------------
```
### Two groups
- Significance: two-tail Student's test, Benjamini–Hochberg correction for multiple hypothesis testing.
- Fold change: mean

### Three groups
- Significance: one-way ANOVA, followed by Tukey's post-hoc test and  Fisher’s Method to determine pvalues. Benjamini–Hochberg correction for multiple hypothesis testing.
- Fold change: One vs the Rest, mean of each group

### more than three groups
- Significance: LIMMA, followed by Fisher’s Method to determine pvalues. Benjamini–Hochberg correction for multiple hypothesis testing.
- Fold change: One vs the Rest. sum of each group

```
# ------------------------------------------------
##                Functional analysis
# ------------------------------------------------
```
### GO enrichment
- Significance: adj.pvalue[FDR] < 0.05

### GSEA enrichment
- Rank: -log10pvalue * log2FC
- Significance: NES > 0 and pvalue<0.05

