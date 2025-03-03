
cat(">> Loading dependencies\n")
library(anndata, warn.conflicts = FALSE)
library(Matrix, warn.conflicts = FALSE)
library(ALRA, warn.conflicts = FALSE)

## VIASH START
par <- list(
  input_train = "resources_test/denoising/pancreas/train.h5ad",
  # input_train = "resources_test/common/pancreas/dataset.h5ad",
  output = "output.h5ad"
)
## VIASH END

cat(">> Load input data\n")
adata <- read_h5ad(par$input_train)

counts <- t(adata$layers[["counts"]])

cat(">> Run ALRA\n")
# alra doesn't work with sparce matrices
out <- alra(as.matrix(counts))

cat(">> Store output\n")
adata$layers[["denoised"]] <- as(t(out$A_norm_rank_k_cor_sc), "CsparseMatrix")
adata$uns[["method_id"]] <- meta[["functionality_name"]]

cat(">> Write output to file\n")
adata$write_h5ad(par$output, compression = "gzip")
