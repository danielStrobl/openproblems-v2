import sys
import anndata as ad
import numpy as np
from utils import *

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
def subsample(l, perc):
    return np.random.choice(l, size=int(len(l)*perc), replace=False)

adata = ad.read(par['input'])
df = convert_to_df(adata)
frac=0.7

train_indices = subsample(df.index, frac)
train = df[df.index.isin(train_indices)]
test = df[~df.index.isin(train_indices)]

adata_train = convert_to_adata(train)
adata_test = convert_to_adata(test)
adata_test.X = np.zeros(adata_test.shape)  # set test data to zero so user cannot see

print(adata_train)
print(adata_test)

print('Writing adatas to file', flush=True)
adata_train.write(par['train_output'], compression='gzip')
adata_test.write(par['test_output'], compression='gzip')
