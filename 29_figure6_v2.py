import numpy as np, pandas as pd, matplotlib; matplotlib.use("Agg"); import matplotlib.pyplot as plt
from scipy.stats import rankdata
hs = pd.read_csv("results/hahn_scores.csv"); gs = pd.read_csv("results/gse57338_scores.csv")
hr = pd.read_csv("results/hahn_random_auc.csv"); gr = pd.read_csv("results/gse57338_random_auc.csv")
ht = pd.read_csv("results/Table_validation_Hahn.csv"); gt = pd.read_csv("results/Table_validation_GSE57338.csv")
COL = {"Control": "#8c8c8c", "HFpEF": "#c0504d", "HFrEF": "#4f81bd", "Failing": "#2e7d32"}
def auc(a, b):
    r = rankdata(np.r_[a, b]); n1, n0 = len(a), len(b); return (r[:n1].sum() - n1 * (n1 + 1) / 2) / (n1 * n0)
def ci(a, b, n=2000):
    rng = np.random.default_rng(42); return np.percentile([auc(rng.choice(a, len(a)), rng.choice(b, len(b))) for _ in range(n)], [2.5, 97.5])
def vals(df, g): return df[df["group"] == g]["score"].values
def box(a, items, title, note):
    for i, (lab, s, g) in enumerate(items, 1):
        a.boxplot(s, positions=[i], widths=.5, showfliers=False, medianprops=dict(color="k"))
        a.scatter(np.random.default_rng(1).normal(i, .06, len(s)), s, s=10, color=COL[g], alpha=.7, zorder=3)
    a.set_xticks(range(1, len(items) + 1)); a.set_xticklabels([x[0] for x in items], fontsize=8); a.set_title(title, fontsize=10, loc="left", fontweight="bold")
    lo, hi = a.get_ylim(); a.set_ylim(lo, hi + 0.3 * (hi - lo))
    a.text(.02, .98, note, transform=a.transAxes, va="top", fontsize=7.5); a.set_ylabel("Culture-time score")
P = "PRIMARY HFpEF vs control RVS"; R = "secondary HFrEF vs control RVS"; L = "secondary HFrEF vs control LVS"
H = lambda n: hs[hs["analysis"] == n]; row = lambda n: ht[ht["analysis"] == n].iloc[0]
fig, ax = plt.subplots(2, 3, figsize=(14, 8.5))
box(ax[0, 0], [("Control", vals(H(P), "Control"), "Control"), ("HFpEF", vals(H(P), "HFpEF"), "HFpEF")], "A  Hahn: HFpEF vs control (RVS)", f"AUC {row(P)['AUC']:.2f}, p = {row(P)['MWU_p']:.1e}")
box(ax[0, 1], [("Control\nRVS", vals(H(R), "Control"), "Control"), ("HFrEF\nRVS", vals(H(R), "HFrEF"), "HFrEF"), ("Control\nLVS", vals(H(L), "Control"), "Control"), ("HFrEF\nLVS", vals(H(L), "HFrEF"), "HFrEF")],
    "B  Hahn: HFrEF vs control", f"RVS AUC {row(R)['AUC']:.2f}, p = {row(R)['MWU_p']:.2f}\nLVS AUC {row(L)['AUC']:.2f}, p = {row(L)['MWU_p']:.1e}")
box(ax[0, 2], [("Non-failing", vals(gs, "Control"), "Control"), ("Failing", vals(gs, "Failing"), "Failing")], "C  GSE57338 (LV)", f"AUC {gt['AUC'][0]:.2f}, p = {gt['MWU_p'][0]:.1e}")
for a, r, o, t, p in [(ax[1, 0], hr[P], row(P)["AUC"], "D  Hahn primary vs 1,000 random gene sets", row(P)["random_p_twosided"]),
                      (ax[1, 1], gr["GSE57338"], gt["AUC"][0], "E  GSE57338 vs 1,000 random gene sets", gt["random_p_twosided"][0])]:
    a.hist(r, bins=40, color="#bbbbbb"); a.axvline(o, color="#c0504d", lw=2, label=f"Observed AUC {o:.3f}"); a.axvline(.5, color="k", ls="--", lw=1, label="No separation")
    a.legend(fontsize=7); a.set_title(f"{t} (two-sided p = {p:.2f})", fontsize=9, loc="left", fontweight="bold"); a.set_xlabel("AUC, failing vs control"); a.set_ylabel("Random sets")
a = ax[1, 2]; labs = []; k = 0
for _, r in ht.iterrows():
    d = H(r["analysis"]); lo, hi = ci(d[d["group"] != "Control"]["score"].values, vals(d, "Control"))
    a.errorbar(r["AUC"], -k, xerr=[[r["AUC"] - lo], [hi - r["AUC"]]], fmt="o", color="#c0504d" if k == 0 else "#4f81bd", capsize=3)
    labs.append(r["analysis"].replace("secondary ", "").replace("sens ", "S: ").replace("PRIMARY ", "")); k += 1
lo, hi = ci(vals(gs, "Failing"), vals(gs, "Control")); a.errorbar(gt["AUC"][0], -k, xerr=[[gt["AUC"][0] - lo], [hi - gt["AUC"][0]]], fmt="o", color="#2e7d32", capsize=3); labs.append("GSE57338 failing vs non-failing")
a.axvline(.5, color="k", ls="--", lw=1); a.set_yticks([-i for i in range(len(labs))]); a.set_yticklabels(labs, fontsize=7)
a.set_xlabel("AUC (95% bootstrap CI)\nlocked expectation: AUC < 0.5", fontsize=8); a.set_title("F  All analyses", fontsize=10, loc="left", fontweight="bold")
for x in ax.ravel(): x.spines[["top", "right"]].set_visible(False)
plt.tight_layout(); plt.savefig("figures/Fig6_validation_v2.pdf"); plt.savefig("figures/Fig6_validation_v2.png", dpi=300); print("saved")
