from os import path
import subprocess
import pandas as pd
import numpy as np

np.random.seed(42)

metric = 'nmi'
metric_file = metric + '.tsv'

print(">> Running script")
out = subprocess.check_output([
    "./" + metric,
    "--adata", 'mnn.h5ad',
    '--hvgs', '2000',
    "--output", metric_file
]).decode("utf-8")

print(">> Checking whether file exists")
assert path.exists(metric_file)

print(">> Check that score makes sense")
result = pd.read_table(metric_file)
assert result.shape == (1, 5)
score = result.loc[0, 'value']
print(score)

assert 0 < score < 1
assert score == 0.4871368591999889

print(">> All tests passed successfully")
