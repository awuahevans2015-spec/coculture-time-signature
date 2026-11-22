# scripts/23_figure6.py -- run from the PROJECT ROOT. Reads results/hahn_scores.csv, hahn_random_auc.csv, Table_validation_Hahn.csv
import numpy as np, pandas as pd, matplotlib; matplotlib.use("Agg"); import matplotlib.pyplot as plt
from scipy.stats import rankdata
sc = pd.read_csv("results/hahn_scores.csv"); ra = pd.read_csv("results/hahn_random_auc.csv"); tb = pd.read_csv("results/Table_validation_Hahn.csv")
def auc(a, b):
    r = rankdata(np.r_[a, b]); n1, n0 = len(a), len(b); return (r[:n1].sum() - n1 * (n1 + 1) / 2) / (n1 * n0)
def boot(a, b, n=2000):
    rng = np.random.default_rng(42)
    v = [auc(rng.choice(a, len(a)), rng.choice(b, len(b))) for _ in range(n)]
    return np.percentile(v, [2.5, 97.5])
fig, ax = plt.subplots(2, 2, figsize=(11, 8.5)); C = {"Control": "#8c8c8c", "HFpEF": "#c0504d", "HFrEF": "#4f81bd"}
def box(a, name, groups, title):
    d = sc[sc["analysis"] == name]; pos = 0
    for g in groups:
        s = d[d["group"] == g]["score"].values; pos += 1
        a.boxplot(s, positions=[pos], widths=.5, showfliers=False, medianprops=dict(color="k"))
        a.scatter(np.random.default_rng(1).normal(pos, .06, len(s)), s, s=14, color=C[g], alpha=.75, zorder=3)
    a.set_xticks(range(1, len(groups) + 1)); a.set_xticklabels(groups); a.set_title(title, fontsize=10); a.set_ylabel("Culture-time score")
box(ax[0, 0], "PRIMARY HFpEF vs control RVS", ["Control", "HFpEF"], "A  HFpEF vs control (RVS), primary")
a = ax[0, 1]; pos = 0
for nm, tis in [("secondary HFrEF vs control RVS", "RVS"), ("secondary HFrEF vs control LVS", "LVS")]:
    d = sc[sc["analysis"] == nm]
    for g in ["Control", "HFrEF"]:
        pos += 1; s = d[d["group"] == g]["score"].values
        a.boxplot(s, positions=[pos], widths=.5, showfliers=False, medianprops=dict(color="k"))
        a.scatter(np.random.default_rng(1).normal(pos, .06, len(s)), s, s=14, color=C[g], alpha=.75, zorder=3)
a.set_xticks([1, 2, 3, 4]); a.set_xticklabels(["Control\nRVS", "HFrEF\nRVS", "Control\nLVS", "HFrEF\nLVS"]); a.set_title("B  HFrEF vs control (secondary)", fontsize=10)
p = tb.iloc[0]; r = ra["PRIMARY HFpEF vs control RVS"]
ax[1, 0].hist(r, bins=40, color="#bbbbbb"); ax[1, 0].axvline(p["AUC"], color="#c0504d", lw=2, label=f"Observed AUC {p['AUC']:.3f}")
ax[1, 0].axvline(.5, color="k", ls="--", lw=1, label="No separation (0.5)"); ax[1, 0].legend(fontsize=8)
ax[1, 0].set_title(f"C  1,000 random gene sets (two-sided p = {p['random_p_twosided']:.2f})", fontsize=10); ax[1, 0].set_xlabel("AUC, failing vs control"); ax[1, 0].set_ylabel("Random sets")
a = ax[1, 1]
for i, row in tb.iterrows():
    d = sc[sc["analysis"] == row["analysis"]]
    lo, hi = boot(d[d["group"] != "Control"]["score"].values, d[d["group"] == "Control"]["score"].values)
    a.errorbar(row["AUC"], -i, xerr=[[row["AUC"] - lo], [hi - row["AUC"]]], fmt="o", color="#c0504d" if i == 0 else "#4f81bd", capsize=3)
a.axvline(.5, color="k", ls="--", lw=1); a.set_yticks([-i for i in range(len(tb))])
a.set_yticklabels([n.replace("secondary ", "").replace("sens ", "S: ").replace("PRIMARY ", "") for n in tb["analysis"]], fontsize=7)
a.set_xlabel("AUC (95% bootstrap CI); expected < 0.5 (locked)"); a.set_title("D  Primary, secondary and exploratory sensitivities", fontsize=10)
for x in ax.ravel(): x.spines[["top", "right"]].set_visible(False)
plt.tight_layout(); plt.savefig("figures/Fig6_validation.pdf"); plt.savefig("figures/Fig6_validation.png", dpi=300); print("saved")
