def pivot(csv, var, obs, pheno):
    layers = {}
    for i in pheno:
        layers[i] = csv.pivot(index=obs, columns=var, values=i)
    return layers