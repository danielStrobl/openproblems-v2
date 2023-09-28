import pandas as pd
import anndata as ad

def convert_to_adata(df):
    label_cols = list(set(df.columns) - set(['value']))
    return ad.AnnData(X=df[['value']], obs=df[label_cols])

def convert_to_df(adata):
    return pd.concat([adata.to_df(), adata.obs], axis=1)
