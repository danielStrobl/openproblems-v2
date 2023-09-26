from typing import Any, Callable, Dict, Tuple
import scanpy as sc
import pandas as pd
import anndata as ad
import utils


def optical_screen(csv_path):
    csv = pd.read_csv(csv_path, index_col=0)
    layers = utils.pivot(csv, 'gene', 'Cell line', ['IL1b_defect', 'TNFa_defect'])
    adata = ad.AnnData(X=next(iter(layers.values())), layers=layers)
    adata.uns['layers'] = list(layers.keys())
    adata.write_h5ad('optical_screen.h5ad', compression="gzip")

