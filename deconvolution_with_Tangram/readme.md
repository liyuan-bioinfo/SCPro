# Tangram for Integrated Analysis of Cell-Type Proteomics and Spatial-Proteomics Data

This repository contains code and documentation for using Tangram, a machine learning library, to perform integrated analysis of cell-type proteomics and spatial-proteomics data.

## Overview
Spatial proteomics involves mapping proteins to specific locations within cells or tissues, while cell-type proteomics involves analyzing protein expression levels within specific cell types. Integrated analysis of these two types of data can provide insights into the heterogeneity of cell types and the spatial organization of proteins within cells.

We performed deconvolution of the raw mouse pancreatic cancer spatial data using the Tangram (version 1.0.4) on log2-normalized mouse pancreatic cancer cell-type proteomics data on top 1000 proteins, and other parameters were set to default.

Tangram is a machine learning library that can be used for integrated analysis of multi-omics data, including spatial proteomics and cell-type proteomics data. It provides a number of tools for model building and evaluation, including cross-validation, hyperparameter tuning, and feature selection.

This repository provides code and documentation for using Tangram to perform integrated analysis of cell-type proteomics and spatial-proteomics data. It includes scripts for preprocessing and normalizing the data, building and evaluating models, and making predictions on new data.

Installation
To use Tangram for integrated analysis of cell-type proteomics and spatial-proteomics data, you will need to install the following dependencies:

- Python 3.x
- NumPy
- Pandas
- Scikit-learn
- Matplotlib
You can install these dependencies using pip:

<code>pip install numpy pandas scikit-learn matplotlib</code>
You will also need to install Tangram. You can install Tangram using pip:

<code> pip install tangram </code>
## Usage
To perform integrated analysis of cell-type proteomics and spatial-proteomics data using Tangram, you can follow these steps:

Preprocess and normalize your data. This might involve aligning spatial proteomics data with imaging data, removing batch effects, and normalizing protein expression levels.
Use Tangram to build a model that integrates both types of data. This might involve using spatial proteomics data to identify protein expression patterns within specific cell types, and then using cell-type proteomics data to identify proteins that are differentially expressed between cell types.
Evaluate your model using cross-validation, hyperparameter tuning, and feature selection.
Use your model to make predictions on new data, or to identify novel biomarkers or pathways that are associated with specific cell types or disease states.
This repository includes scripts for performing each of these steps, as well as documentation on how to use Tangram for integrated analysis of cell-type proteomics and spatial-proteomics data.