## VIASH START
par = {
    'adata': './src/batch_integration/resources/datasets_pancreas.h5ad',
    'output': './src/batch_integration/resources/pancreas_bbknn.h5ad',
    'hvg': True,
    'scaling': True,
    'debug': True
}
## VIASH END

print('Importing libraries')
from pprint import pprint
import scanpy as sc
from scib.integration import scanorama

if par['debug']:
    pprint(par)

adata_file = par['adata']
output = par['output']
hvg = par['hvg']
scaling = par['scaling']

print('Read adata')
adata = sc.read(adata_file)

if hvg:
    print('Select HVGs')
    adata = adata[:, adata.var['highly_variable']]

if scaling:
    print('Scale')
    adata.X = adata.layers['logcounts_scaled']
else:
    adata.X = adata.layers['logcounts']

print('Integrate')
adata.X = scanorama(adata, batch='batch').X.todense()

print('Postprocess data')
sc.pp.pca(
    adata,
    n_comps=50,
    use_highly_variable=True,
    svd_solver='arpack',
    return_info=True
)
sc.pp.neighbors(adata, use_rep='X_pca')

print('Save HDF5')
adata.uns['hvg'] = hvg
adata.uns['scaled'] = scaling

adata.write(output, compression='gzip')
