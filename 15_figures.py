# scripts/15_figures.py -- run from the PROJECT ROOT:  .venv\Scripts\python.exe scripts\15_figures.py
# Figures 3 and 4 from existing results files only (read-only on results; writes figures/Fig3_ML.pdf, figures/Fig4_signature.pdf)
import numpy as np, pandas as pd
from scipy.stats import spearmanr
import matplotlib; matplotlib.use("Agg"); import matplotlib.pyplot as plt

pd.set_option("display.width", 200)
t1 = pd.read_csv("results/Table1_performance.csv", keep_default_na=False)
print("Table1 columns:", list(t1.columns))
for c in ["mae", "null_mae"]:
    t1[c] = pd.to_numeric(t1[c], errors="coerce")

def panel_scatter(ax, target, ver, model, xlabel, ylabel, title):
    p = pd.read_csv(f"results/pred_{target}_{ver}_{model}.csv")
    ax.scatter(p["observed"], p["predicted"], s=40, zorder=2)
    lo, hi = min(p["observed"].min(), p["predicted"].min()), max(p["observed"].max(), p["predicted"].max())
    ax.plot([lo, hi], [lo, hi], color="0.6", lw=1, ls="--", zorder=1)
    rho = spearmanr(p["observed"], p["predicted"])[0]
    mae = float(np.mean(np.abs(p["observed"] - p["predicted"])))
    ax.text(0.04, 0.95, f"MAE = {mae:.3f}\nrho = {rho:.2f}", transform=ax.transAxes, va="top", fontsize=8)
    ax.set_xlabel(xlabel); ax.set_ylabel(ylabel); ax.set_title(title, fontsize=10)

# ---------------- Figure 3
fig, ax = plt.subplots(1, 3, figsize=(12, 3.8))
panel_scatter(ax[0], "log2_day", "A", "enet", "Observed log2(culture day)", "Predicted (leave-one-sample-out)",
              "Culture day, elastic net (Version A)")
labels, specs = [], [("ridge", "A"), ("enet", "A"), ("ridge", "B"), ("enet", "B")]
targets = [("log2_day", "Culture day"), ("peg_c", "Crosslinking %")]
w = 0.16
for gi, (tg, tl) in enumerate(targets):
    sub = t1[t1["target"] == tg]
    null_mae = sub["null_mae"].dropna().iloc[0]
    vals = [null_mae] + [float(sub[(sub["model"] == m) & (sub["version"] == v)]["mae"].iloc[0]) for m, v in specs]
    names = ["null"] + [f"{m} {v}" for m, v in specs]
    for k, (val, nm) in enumerate(zip(vals, names)):
        ax[1].bar(gi + (k - 2) * w, val, width=w, color=f"C{k}", label=nm if gi == 0 else None)
ax[1].set_xticks([0, 1]); ax[1].set_xticklabels([t[1] for t in targets])
ax[1].set_ylabel("MAE (leave-one-sample-out)")
ax[1].set_title("Error versus null model", fontsize=10); ax[1].set_ylim(0, 1.3)
ax[1].legend(fontsize=7, ncol=3, loc="upper center", title="Model, version", title_fontsize=7)
panel_scatter(ax[2], "peg_c", "A", "enet", "Observed crosslinking (centered %)", "Predicted (leave-one-sample-out)",
              "Crosslinking, elastic net (Version A): not detectable")
for k, a in enumerate(ax):
    a.text(-0.16, 1.06, "ABC"[k], transform=a.transAxes, fontsize=12, fontweight="bold")
fig.tight_layout(); fig.savefig("figures/Fig3_ML.pdf"); plt.close(fig)
print("saved figures/Fig3_ML.pdf")

# ---------------- Figure 4
enr = pd.read_csv("results/enrichment_log2_day.csv")
enr = enr.sort_values("p.adjust").head(10).iloc[::-1]
def short(s):
    for pre in ["REACTOME_", "GOBP_", "KEGG_", "HALLMARK_"]:
        if s.startswith(pre): s = s[len(pre):]
    s = s.replace("_", " ").lower()
    return s if len(s) <= 48 else s[:45] + "..."
sig = pd.read_csv("results/signature_log2_day_A.csv", keep_default_na=False)
sig = sig[sig["in_signature"].astype(str) == "True"].copy()
sig["ribo"] = sig["gene"].str.match("^(RPL|RPS|MRPL|MRPS)[0-9]")
sig["direction"] = pd.to_numeric(sig["direction"])
fig, ax = plt.subplots(1, 2, figsize=(11, 4), gridspec_kw={"width_ratios": [2.2, 1]})
ax[0].barh([short(s) for s in enr["ID"]], -np.log10(enr["p.adjust"]), color="C0")
for i, c in enumerate(enr["Count"]):
    ax[0].text(0.2, i, f"n = {c}", va="center", color="white", fontsize=7)
ax[0].set_xlabel("-log10 adjusted p (BH)"); ax[0].set_title("Top enriched gene sets, culture-day signature", fontsize=10)
ax[0].tick_params(axis="y", labelsize=7)
up, dn = sig[sig["direction"] > 0], sig[sig["direction"] < 0]
ax[1].bar(["Up with day", "Down with day"], [(~up["ribo"]).sum(), (~dn["ribo"]).sum()], label="Other genes", color="C0")
ax[1].bar(["Up with day", "Down with day"], [up["ribo"].sum(), dn["ribo"].sum()],
          bottom=[(~up["ribo"]).sum(), (~dn["ribo"]).sum()], label="Ribosomal (RPL/RPS/MRPL/MRPS)", color="C1")
ax[1].set_ylabel("Signature genes"); ax[1].set_title(f"Signature composition (n = {len(sig)})", fontsize=10)
ax[1].set_ylim(0, 115); ax[1].legend(fontsize=7, loc="upper center")
for k, a in enumerate(ax):
    a.text(-0.12 if k == 1 else -0.5, 1.06, "AB"[k], transform=a.transAxes, fontsize=12, fontweight="bold")
fig.tight_layout(); fig.savefig("figures/Fig4_signature.pdf"); plt.close(fig)
print("saved figures/Fig4_signature.pdf")
