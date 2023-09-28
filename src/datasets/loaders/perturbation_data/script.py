import pandas as pd
import anndata as ad
import os

## VIASH START
# The following code has been auto-generated by Viash.
par = {
  'output': r'resources_test/common/pancreas/raw.h5ad'
}
meta = {
  'functionality_name': r'optical_screen',
  'resources_dir': r'/private/tmp/viash_inject_optical_screen6586960705516935516',
  'executable': r'/private/tmp/viash_inject_optical_screen6586960705516935516/optical_screen',
  'config': r'/private/tmp/viash_inject_optical_screen6586960705516935516/.config.vsh.yaml',
  'temp_dir': r'/var/folders/jl/1_41b62j2bl1bqykfj7q9z0h6c58h5/T/',
  'cpus': int(r'123'),
  'memory_b': int(r'123'),
  'memory_kb': int(r'123'),
  'memory_mb': int(r'123'),
  'memory_gb': int(r'123'),
  'memory_tb': int(r'123'),
  'memory_pb': int(r'123')
}

## VIASH END

def convert_to_adata(df):
    label_cols = list(set(df.columns) - set(['value']))
    return ad.AnnData(X=df[['value']], obs=df[label_cols])

def optical_screen(csv_path):
    csv = pd.read_csv(csv_path, index_col=0)
    adata = convert_to_adata(csv)
    print(adata)
    print(par)
    adata.write_h5ad(par["output"], compression="gzip")

print(os.listdir(meta['resources_dir']))
optical_screen(meta['resources_dir']+"/optical_secondary_processed.csv")
