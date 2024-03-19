#'@time 202403
#'@author Yuan
#'@desc Integration cell-type proteomics and spatial proteomics with Tangram

import scanpy as sc
import tangram as tg
import pandas as pd
import sys
import os
import numpy as np
import pandas as pd

import warnings
warnings.filterwarnings("ignore")
os.chdir("") # change this to YOUR dir


def run_Tangram(sc_adata, sp_adata, output_file_name, celltype_key="celltype"):
    ## Find DEG for sc
    sc.pp.log1p(sc_adata)
    sc.tl.rank_genes_groups(sc_adata, groupby=celltype_key, use_raw=False)

    markers_df = pd.DataFrame(sc_adata.uns["rank_genes_groups"]["names"]).iloc[0:1000, :] # top 1000 sig. proteins use for deconvolution

    genes_sc = np.unique(markers_df.melt().value.values)
    genes_st = sp_adata.var_names.values
    genes = list(set(genes_sc).intersection(set(genes_st)))
    tg.pp_adatas(sc_adata, sp_adata, genes=genes)

    ad_map = tg.map_cells_to_space(
                    sc_adata,
                    sp_adata,
                    mode='clusters',
                    cluster_label=celltype_key)

    tg.project_cell_annotations(ad_map, sp_adata, annotation=celltype_key)

    celltype_density = sp_adata.obsm['tangram_ct_pred']
    celltype_density = (celltype_density.T/celltype_density.sum(axis=1)).T
    
    celltype_density.to_csv(output_file_name)

    return(celltype_density)

# ----------------------------------------------------------------------------
# run Tangram
# Predicted proportion of 4 celllineage for the spatial-data (PCC, CAF, MYE and LYM)
# ----------------------------------------------------------------------------

sp_adata_path = "raw/Tangram_deconvolution/PDAC2023_Spatial_20240319_5844.h5ad" # spatial-proteomics data
sc_adata_path = "raw/Tangram_deconvolution/PDAC2023_ct14_20240319_5900.h5ad" #cell-type proteomics data
save_path = "raw/Tangram_deconvolution/ct14_cellperc_tangram.csv"

sp_adata = sc.read_h5ad(sp_adata_path)
sc_adata = sc.read_h5ad(sc_adata_path)

intersect = np.intersect1d(sc_adata.var_names, sp_adata.var_names)
sc_adata = sc_adata[:, intersect].copy()
sp_adata = sp_adata[:, intersect].copy()

# The subtypes of CAFs usual have immune cell characteristics, so we do not consider them here.
sc_adata = sc_adata[~sc_adata.obs["celltype"].str.contains('^iCAF$', case=False)].copy()
sc_adata = sc_adata[~sc_adata.obs["celltype"].str.contains('^myCAF$', case=False)].copy()
sc_adata = sc_adata[~sc_adata.obs["celltype"].str.contains('^apCAF$', case=False)].copy()
run_Tangram(sc_adata=sc_adata, sp_adata=sp_adata, output_file_name=save_path, celltype_key="celltype")

# ----------------------------------------------------------------------------
# run Tangram
# Predicted proportion of 8 Treg-related celltypes for the spatial-data
# ----------------------------------------------------------------------------
sp_adata_path = "raw/Tangram_deconvolution/PDAC2023_Spatial_20240319_5844.h5ad"
sc_adata_path = "raw/Tangram_deconvolution/PDAC2023_ct8_20240319_3824.h5ad"
save_path = "raw/Tangram_deconvolution/ct8_cellperc_tangram.csv"

sp_adata = sc.read_h5ad(sp_adata_path)
sc_adata = sc.read_h5ad(sc_adata_path)
intersect = np.intersect1d(sc_adata.var_names, sp_adata.var_names)
sc_adata = sc_adata[:, intersect].copy()
sp_adata = sp_adata[:, intersect].copy()
run_Tangram(sc_adata=sc_adata, sp_adata=sp_adata, output_file_name=save_path, celltype_key="celltype")
