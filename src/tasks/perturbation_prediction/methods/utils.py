import pandas as pd
import anndata as ad
import numpy as np
from sklearn.preprocessing import OneHotEncoder
from scipy.stats import spearmanr
def convert_to_df(adata):
    return pd.concat([adata.to_df(), adata.obs], axis=1)

def X_cats_to_emb(
    data_cat_df: pd.DataFrame,
    gene_emb: pd.DataFrame = None,
    cell_line_emb: pd.DataFrame = None,
    drug_emb: pd.DataFrame = None,
    cell_line_col: str = None,
    gene_col: str = None,
    drug_col: str = None,
    phenotype_col: str = None,
    verbose = True
):
    """
    Converts dataframe with X in categorical from to X in embedded form.

    Parameters
    ----------
    data_cat_df - pd.DataFrame with observations as rows, and X as categorical
        in columns (e.g. gene column specifying gene name perturbed per row, or
        drug column specifying the drug applied)
    gene_emb - pd.DataFrame with gene name as index and embedding dimensions as
        columns, with one embedding vector per gene. Only needed if gene part of
        data is present.
    drug_emb - pd.DataFrame with drug name as index and embedding dimensions as
        columns, with one embedding vector per gene. Only needed if drug part of
        data is present.
    cell_line_emb - as gene_emb, but for cell lines
    cell_line_col - name of the cell line column in data_cat_df specifying cell
        line per observation. Only needed if cell line in X and cell_line_emb
        provided.
    gene_col - as cell_line_col but for genes
    drug_col - as cell_line_col but for drugs
    phenotype_col - as cell_line_col but for phenotypes. No embedding needed,
        phenotype will be one-hot encoded
    verbose - whether to print extra info

    Returns:
    np.array with concatenated observation-specific embeddings, one row per
    observation
    """
    dims_to_convert = dict()
    # check if cell line should be included
    if isinstance(cell_line_col, str):
        if cell_line_col not in data_cat_df.columns:
            raise ValueError(f"cell_line_col {cell_line_col} not in your data_cat_df.")
        dims_to_convert['cell_line'] = (cell_line_col)
        if verbose:
            print("Using cell line embedding as part of output X embedding.")
    # check if perturbed gene should be included
    if isinstance(gene_col, str):
        if gene_col not in data_cat_df.columns:
            raise ValueError(f"gene_col {gene_col} not in your data_cat_df.")
        dims_to_convert['gene'] = gene_col
        if isinstance(drug_col, str):
            raise ValueError(f"You cannot use both the gene and drug column for embedding.")
        if verbose:
            print("Using gene embedding as part of output X embedding.")
    # check if drug should be included:
    elif isinstance(drug_col, str):
        if drug_col not in data_cat_df.columns:
            raise ValueError(f"drug_col {drug_col} not in your data_cat_df.")
        dims_to_convert['drug'] = drug_col
        if verbose:
            print("Using drug embedding as part of output X embedding.")
    # check if phenotype should be included:
    if isinstance(phenotype_col, str):
        if phenotype_col not in data_cat_df.columns:
            raise ValueError(f"phenotype_col {phenotype_col} not in your data_cat_df.")
        dims_to_convert['phenotype'] = phenotype_col
        # convert to one-hot encoded array
        oh_enc = OneHotEncoder()
        oh_enc.fit(data_cat_df[phenotype_col].values[:,np.newaxis])
        phenotype_one_hot = oh_enc.transform(
            data_cat_df[phenotype_col].values[:,np.newaxis]
        ).toarray()
        if verbose:
            print("Using one-hot encoded phenotype as part of output X embedding.")
        
    # generate embedding:
    embs = list()
    for dim, dim_col_name in dims_to_convert.items():
        if dim == "cell_line":
            embs.append(cell_line_emb.loc[data_cat_df[dim_col_name],:])
        elif dim == "gene":
            embs.append(gene_emb.loc[data_cat_df[dim_col_name],:])
        elif dim == "drug":
            embs.append(drug_emb.loc[data_cat_df[dim_col_name],:])
        elif dim == "phenotype":
            embs.append(phenotype_one_hot)
    X_emb = np.concatenate(tuple(embs),axis=1)
    return X_emb
