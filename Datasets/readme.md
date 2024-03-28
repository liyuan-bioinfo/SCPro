# Datasets of spatial proteomics and cell-type proteomics
sp_dataset_01: 3 regions (n = 4), including Acinar, Lymph and Tumor.
sp_dataset_02: 2 regions (n = 4), including PCC and Stroma.
sp_dataset_03: 6 regions (n = 3), including Acinar, PanIN, PDAC, CAF, IT and LN.
ct_dataset_01: 14 cell types (n = 4 for apCAF, n = 5 for the other 13 cell types), including PCC, CAF, iCAF, myCAF, apCAF, T4, T8, Treg, B, MYE, NEU, MO, MAC and DC.
ct_dataset_02: 8 cell types (n = 3) of subtype of Treg, including CD4+CD25+Klrg1+ Treg (also named Klrg1+ Treg) and CD4+CD25+Klrg1- Treg (also named Klrg1- Treg).

# Quantification and Statistical analysis
- For spatial proteomics dataset 1 and dataset 2, we filtered the quantified protein groups for at least three valid values in at least one region. Missing values were imputed by drawing from a normal distribution (width = 0.3; downshift = 1.8) from sample 's proteome abundance. The same imputation method was applied for all spatial proteomics datasets. Significance was calculated by one-way ANOVA for dataset 1 followed by permutation-based FDR for multiple hypothesis testing (FDR < 0.05). Two-tailed Students’s t-test was performed for dataset 2, average proteins’ abundance of each region was used to calculate difference.
- For spatial proteomics dataset 3, we filtered the quantified protein groups for at least three valid values in at least one region. Significance of EpCAM+ cells (n = Acinar, PanIN and PDAC) was calculated by one-way ANOVA and median proteins’ abundance of one vs the rest cells were used to calculated difference. Significance of IT and LN (n = 3) was calculated by two-tailed Students’s t-test, median proteins’ abundance of each region was used to calculate difference.
- For cell-type proteomics dataset 1 and dataset 2, we filtered the quantified proteins for at least two valid values in at least one cell type. Missing values were imputed by 0.1* minimum intensity of each protein. Significance and difference were calculated by LIMMA package (version 3.48.0) in R with p-value < 0.05 and fold change > 2.
