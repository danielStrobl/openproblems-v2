import yaml
import pandas as pd
import anndata as ad
from sklearn.linear_model import Lasso

## VIASH START
par = {
    'gene_emb':'gene_emb.csv',
    'cl_emb':'cl_emb.csv',
    'train':'A549splittrain.h5ad',
    'test':'A549splittest.h5ad'  # must not contain real values
}
meta = {
    'functionality_name': 'foo',
    'config': 'bar'
}
## VIASH END

print('Read input', flush=True)
input = ad.read_h5ad(par['train'])

train = convert_to_df(input)

# load embeddings
url = 'https://raw.githubusercontent.com/danielStrobl/perturbation_embeddings/main/CCLE_expression_gene_emb.csv'
gene_emb = pd.read_csv(url, index_col=0)
url = 'https://raw.githubusercontent.com/danielStrobl/perturbation_embeddings/main/CCLE_expression_cell_line_emb.csv'
cl_emb = pd.read_csv(url, index_col=0)

X = X_cats_to_emb(
    train,
    pd.read_csv(gene_emb, index_col=0),
    pd.read_csv(cl_emb, index_col=0),
    cell_line_col='cell_line',
    gene_col='gene',
    phenotype_col='phenotype'
)
X_test = X_cats_to_emb(
    test,
    pd.read_csv(gene_emb, index_col=0),
    pd.read_csv(cl_emb, index_col=0),
    cell_line_col='cell_line',
    gene_col='gene',
    phenotype_col='phenotype'
)
y = adata.X.toarray()


reg = Lasso()
reg.fit(X, y)
pred = reg.predict(X_test)

# convert predictions back to anndata
adata_test = test.copy()
adata_test.X = pred

print("Store outputs", flush=True)
adata_test.write_h5ad(par['output'], compression='gzip')
