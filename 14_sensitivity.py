# scripts/14_sensitivity.py -- run from the PROJECT ROOT:  .venv\Scripts\python.exe scripts\14_sensitivity.py
# POST HOC, exploratory (decided after viewing 10_enrichment output; logged in docs/code_changes.md 2026-11-18).
# Culture-day model with ribosomal/translation genes removed. Writes results/sensitivity_no_translation.csv only.
import numpy as np, pandas as pd
from mlcore import *

X, meta = load("primary")
y = meta["log2_day"].values
enr = pd.read_csv("results/enrichment_log2_day.csv").sort_values("p.adjust").head(10)
sets = pd.read_csv("data/genesets/msigdb_sets_used.csv")
drop = set(sets[sets["gs_name"].isin(enr["ID"])]["gene_symbol"])
drop |= {g for g in X.columns if pd.Series([g]).str.match("^(RPL|RPS|MRPL|MRPS)[0-9]").iloc[0]}
drop = [g for g in X.columns if g in drop]
Xs = X.drop(columns=drop)
print("genes before:", X.shape[1], "| removed:", len(drop), "| after:", Xs.shape[1])
null_pred = np.array([np.delete(y, i).mean() for i in range(len(y))])   # leave-one-out training mean
null_mae = float(np.mean(np.abs(y - null_pred)))
rows = []
for model in ["enet", "ridge"]:
    for label, Xm in [("all genes", X), ("translation genes removed", Xs)]:
        rho, mae = rho_mae(y, loo_predict(model, "A", Xm.values, y))
        rows.append(dict(model=model, input=label, n_genes=Xm.shape[1], rho=rho, mae=mae, null_mae=null_mae,
                         mae_gain_vs_null_pct=100 * (null_mae - mae) / null_mae))
        print(rows[-1], flush=True)
out = pd.DataFrame(rows)
out.to_csv("results/sensitivity_no_translation.csv", index=False)
print(out.round(4).to_string(index=False))
