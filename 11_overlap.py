# scripts/11_overlap.py -- run from the PROJECT ROOT:  .venv\Scripts\python.exe scripts\11_overlap.py
# Fisher exact test (one-sided): does each signature overlap the frozen panels / disease genes more than chance?
# Outputs results/Table_overlap.csv. Version B signatures come from inside the frozen panel, so their overlap with
# panel gene sets is partly circular (flagged in the 'note' column).
import glob, os, numpy as np, pandas as pd
from scipy.stats import fisher_exact

bg = set(pd.read_csv("data/processed/filtered_gene_background.csv", keep_default_na=False).iloc[:, 0])
sets = {}
for f in sorted(glob.glob("data/genesets/*.txt")):
    name = os.path.basename(f).rsplit("_", 1)[0]
    sets[name] = {g.strip() for g in open(f) if g.strip()}

def bh(p):                                              # Benjamini-Hochberg adjusted p-values
    p = np.asarray(p, float); n = len(p); order = np.argsort(p)
    adj = np.empty(n); prev = 1.0
    for rank, i in zip(range(n, 0, -1), order[::-1]):
        prev = min(prev, p[i] * n / rank); adj[i] = prev
    return adj

rows = []
for version in ["A", "B"]:
    for target in ["log2_day", "peg_c"]:
        s = pd.read_csv(f"results/signature_{target}_{version}.csv", keep_default_na=False, na_values=[""])
        sig = set(s[s["in_signature"].astype(str) == "True"].gene) & bg
        if len(sig) == 0:
            continue
        for name, gs in sets.items():
            gs = gs & bg
            a = len(sig & gs); b = len(sig - gs); c = len(gs - sig); d = len(bg) - a - b - c
            p = fisher_exact([[a, b], [c, d]], alternative="greater")[1]
            note = "partly circular (B comes from the panel)" if version == "B" else ""
            rows.append(dict(version=version, target=target, geneset=name, signature_size=len(sig),
                             overlap=a, set_size=len(gs), p=p, note=note))
df = pd.DataFrame(rows)
df["FDR"] = float("nan")
primary = (df["version"] == "A").values            # FDR over the primary (non-circular) Version A tests only
df.loc[primary, "FDR"] = bh(df.loc[primary, "p"].values)
df.to_csv("results/Table_overlap.csv", index=False)
pd.set_option("display.width", 220)
print(df.round(4).to_string(index=False))
