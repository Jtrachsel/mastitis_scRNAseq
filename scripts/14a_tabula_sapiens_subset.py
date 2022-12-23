import numpy as np
import pandas as pd
import anndata as ad
from scipy.sparse import csr_matrix
print(ad.__version__)

tabula=ad.read('./reference_mapping_data/TabulaSapiens.h5ad')
# tabula
# tabula.obs
# tabula.obs.organ_tissue

TISSUES=tabula.obs['organ_tissue'].unique()
for tissue in TISSUES:
        print(tissue)

# create some boolean 'series' to subset the anndata object
cells_blood=tabula.obs['organ_tissue'] == 'Blood'
cells_10X=tabula.obs['method'] == '10X'

# keep only cells from Blood and sequenced on the 10X platform
tabula_10X_blood = tabula[cells_blood & cells_10X,]

# remove raw counts because we only are interested in the counts with 
# ambient RNA removed
# tabula_10X_blood.layers.keys()
del tabula_10X_blood.layers["raw_counts"]

# remove some other data to help ease the conversion to seurat object
del tabula_10X_blood.uns
del tabula_10X_blood.obsm
del tabula_10X_blood.obsp



tabula_10X_blood.write('reference_mapping_data/tabula_10X_blood.h5ad', compression="gzip")

