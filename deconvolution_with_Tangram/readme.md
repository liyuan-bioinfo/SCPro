# Integrated Analysis of Cell-Type Proteomics and Spatial-Proteomics Data with Tangram

This repository contains code and documentation for using Tangram, a deep-learning, to perform integrated analysis of cell-type proteomics and spatial-proteomics data.

## Overview
Spatial proteomics involves mapping proteins to specific locations within cells or tissues, while cell-type proteomics involves analyzing protein expression levels within specific cell types. Integrated analysis of these two types of data can provide insights into the heterogeneity of cell types and the spatial organization of proteins within cells.

We performed deconvolution of the raw mouse pancreatic cancer spatial data using the Tangram (version 1.0.4) on log2-normalized mouse pancreatic cancer cell-type proteomics data on top 1000 proteins, and other parameters were set to default.

This repository provides code and documentation for using Tangram to perform integrated analysis of cell-type proteomics and spatial-proteomics data. It includes scripts for preprocessing and normalizing the data, building and evaluating models, and making predictions on new data.

Installation
To use Tangram for integrated analysis of cell-type proteomics and spatial-proteomics data, you will need to install the following dependencies:

- Python 3.x
- NumPy
- Pandas
- Matplotlib
You can install these dependencies using pip:
<code>pip install numpy pandas scikit-learn matplotlib</code>

You will also need to install Tangram. You can install Tangram using pip:
<code> pip install tangram </code>
