# scripts/09_composition.py -- run from the PROJECT ROOT:  .venv\Scripts\python.exe scripts\09_composition.py
# Composition vs regulation. Does the signature score still follow day / crosslinking after the
# endothelial:pericyte marker ratio is removed? Outputs results/composition_v2.csv and figures/Fig5_markers.pdf
# NOTE: the signature was selected using these same 11 samples, so rho_score_vs_target is optimistic (circular);
# use it descriptively. The informative comparison is rho_ratio_removed and the model with/without marker genes.
# Patch 2026-11-18: added MAE (with/without marker genes) to the output; figure gets group lines and panel labels.
# Output renamed to composition_v2.csv because results/composition.csv could not be overwritten (file locked).
import numpy as np, pandas as pd
from scipy.stats import spearmanr
from mlcore import *

X_all, meta = load("primary")
Z = (X_all - X_all.mean()) / X_all.std()
endo = [g for g in ["PECAM1", "CDH5", "KDR"] if g in Z.columns]
peri = [g for g in ["PDGFRB", "CSPG4", "RGS5"] if g in Z.columns]
print("Endothelial markers used:", endo, "| Pericyte markers used:", peri)
endo_score, peri_score = Z[endo].mean(axis=1), Z[peri].mean(axis=1)
ratio = endo_score - peri_score                      # endothelial minus pericyte score


def resid_on(score, x):                              # score with the cell ratio removed (ordinary least squares)
    A = np.c_[np.ones(len(x)), np.asarray(x, float)]
    coef = np.linalg.lstsq(A, np.asarray(score, float), rcond=None)[0]
    return np.asarray(score, float) - A @ coef


rows = []
for target in ["log2_day", "peg_c"]:
    sig = pd.read_csv(f"results/signature_{target}_A.csv", keep_default_na=False, na_values=[""])
    sig = sig[sig["in_signature"].astype(str) == "True"]
    sig = sig[~sig.gene.isin(MARKERS)]                # avoid circular reasoning
    if len(sig) < 3:
        rows.append(dict(target=target, sig_size=len(sig), note="fewer than 3 signature genes; skipped"))
        print(target, ": fewer than 3 signature genes, skipped"); continue
    up, down = list(sig[sig.direction > 0].gene), list(sig[sig.direction < 0].gene)
    score = (Z[up].mean(axis=1) if up else 0) - (Z[down].mean(axis=1) if down else 0)
    score = pd.Series(score, index=Z.index)
    resid = resid_on(score.values, ratio.values)
    r1, p1 = spearmanr(meta[target], score)
    r2, p2 = spearmanr(meta[target], resid)
    r3, p3 = spearmanr(score, ratio)
    Xn = X_all.drop(columns=[g for g in MARKERS if g in X_all.columns])   # model without marker genes
    y = meta[target].values
    r_no, mae_no = rho_mae(y, loo_predict("enet", "A", Xn.values, y))
    r_with, mae_with = rho_mae(y, loo_predict("enet", "A", X_all.values, y))
    rows.append(dict(target=target, sig_size=len(sig), rho_score_vs_target=r1, rho_ratio_removed=r2,
                     rho_score_vs_ratio=r3, p_score_vs_target=p1, p_ratio_removed=p2, p_score_vs_ratio=p3,
                     enet_rho_with_markers=r_with, enet_rho_without_markers=r_no,
                     enet_mae_with_markers=mae_with, enet_mae_without_markers=mae_no))
    print(rows[-1], flush=True)
df = pd.DataFrame(rows)
df.to_csv("results/composition_v2.csv", index=False)
pd.set_option("display.width", 200)
print(df.T)

# ---- Figure 5: marker z-scores against day and crosslinking % (descriptive; one sample per point, n = 11)
try:
    import matplotlib; matplotlib.use("Agg"); import matplotlib.pyplot as plt
    fig, ax = plt.subplots(1, 3, figsize=(11, 3.6))
    mk = {35: "o", 50: "^", 60: "s"}
    for k, (ttl, v) in enumerate([("Endothelial markers", endo_score), ("Pericyte markers", peri_score),
                                  ("Endothelial minus pericyte", ratio)]):
        for peg, m in mk.items():
            idx = meta.peg_pct.values == peg
            d, s = meta.day.values[idx], v.values[idx]
            o = np.argsort(d)
            ax[k].plot(d[o], s[o], color="0.8", lw=1, zorder=1)          # faint line joining each group
            ax[k].scatter(d, s, marker=m, s=45, label=f"{peg}%", zorder=2)
        ax[k].set_xlabel("Culture day"); ax[k].set_ylabel("Mean z-score"); ax[k].set_title(ttl, fontsize=10)
        ax[k].text(-0.18, 1.05, "ABC"[k], transform=ax[k].transAxes, fontsize=12, fontweight="bold")
    ax[0].legend(title="Crosslinking", fontsize=8)
    fig.tight_layout(); fig.savefig("figures/Fig5_markers.pdf"); plt.close(fig)
    print("saved figures/Fig5_markers.pdf")
except Exception as e:
    print("Figure 5 not made:", type(e).__name__, str(e)[:200])
