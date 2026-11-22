# scripts/10_enrichment.py -- run from the PROJECT ROOT:  .venv\Scripts\python.exe scripts\10_enrichment.py
# Over-representation of the signature genes in MSigDB gene sets (hypergeometric test, BH adjustment), using the
# genes that passed the expression filter as the background. Same test as clusterProfiler::enricher (not installed here).
# Needs data/genesets/msigdb_sets_used.csv, exported from R (see the R block in the instructions).
import numpy as np, pandas as pd
from scipy.stats import hypergeom

bg = set(pd.read_csv("data/processed/filtered_gene_background.csv", keep_default_na=False).iloc[:, 0])
ms = pd.read_csv("data/genesets/msigdb_sets_used.csv", keep_default_na=False)
ms = ms[ms.gene_symbol.isin(bg)].drop_duplicates()
sets = {k: set(v) for k, v in ms.groupby("gs_name")["gene_symbol"]}
sets = {k: v for k, v in sets.items() if 10 <= len(v) <= 500}          # same size limits as the guide
print("gene sets tested:", len(sets), "| background genes:", len(bg))

def bh(p):
    p = np.asarray(p, float); n = len(p); order = np.argsort(p)
    adj = np.empty(n); prev = 1.0
    for rank, i in zip(range(n, 0, -1), order[::-1]):
        prev = min(prev, p[i] * n / rank); adj[i] = prev
    return adj

N = len(bg)
for target in ["log2_day", "peg_c"]:
    s = pd.read_csv(f"results/signature_{target}_A.csv", keep_default_na=False, na_values=[""])
    sig = set(s[s["in_signature"].astype(str) == "True"].gene) & bg
    if len(sig) < 5:
        print(target, ": fewer than 5 signature genes, skipped"); continue
    n = len(sig); rows = []
    for name, gs in sets.items():
        K = len(gs); k = len(sig & gs)
        p = hypergeom.sf(k - 1, N, K, n)
        rows.append(dict(ID=name, Count=k, GeneRatio=f"{k}/{n}", BgRatio=f"{K}/{N}", pvalue=p))
    out = pd.DataFrame(rows)
    out["p.adjust"] = bh(out["pvalue"].values)
    out = out.sort_values("p.adjust")
    out.to_csv(f"results/enrichment_{target}.csv", index=False)
    print("\n", target, "- signature genes:", n)
    print(out.head(10)[["ID", "Count", "p.adjust"]].to_string(index=False))
