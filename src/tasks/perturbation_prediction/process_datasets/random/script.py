import sys
import anndata as ad
import numpy as np

## VIASH START
par = {
    'input': 'resources_test/common/pancreas/dataset.h5ad',
    'frac': 0.7,
    'train_output': 'train_output.h5ad',
    'test_output': 'test_output.h5ad'
}
meta = {}
## VIASH END

# Remove this after upgrading to Viash 0.7.5
sys.dont_write_bytecode = True

adata_raw = ad.read(par['input'])
frac=0.7


def subsample(l, perc):
    return np.random.choice(l, size=int(len(l)*perc), replace=False)


var_split = subsample(adata_raw.var_names, frac)
obs_split = subsample(adata_raw.obs_names, frac)

adata_train = adata_raw[obs_split, var_split]
adata_test = adata_raw[list(set(adata_raw.obs_names)-set(obs_split)), list(set(adata_raw.var_names)-set(var_split))]

print(adata_train)
print(adata_test)


print('Writing adatas to file', flush=True)
adata_train.write(par['train_output'], compression='gzip')
adata_test.write(par['test_output'], compression='gzip')
