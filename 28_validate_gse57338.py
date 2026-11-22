# scripts/28_validate_gse57338.py -- run from the PROJECT ROOT. Locked: docs/expected_direction_GSE57338_2026-11-21.txt
# AUC = P(score failing > score non-failing). H0 (original): AUC < 0.5. H1 (replicate Hahn): AUC > 0.5. Exploratory.
import numpy as np, pandas as pd
from scipy.stats import rankdata, mannwhitneyu, chi2
from sklearn.linear_model import LogisticRegression
from sklearn.metrics import log_loss
sig = pd.read_csv("results/signature_log2_day_A.csv", keep_default_na=False)
sig = sig[sig["in_signature"].astype(str) == "True"].copy()
sig["sign"] = sig["direction"].map(lambda x: 1 if float(x) > 0 else -1)
UP = sig.loc[sig["sign"] > 0, "gene"].tolist(); DOWN = sig.loc[sig["sign"] < 0, "gene"].tolist(); ALLSIG = set(UP + DOWN)
E = pd.read_csv("data/processed/gse57338_expr_symbols.csv", index_col=0)
meta = pd.read_csv("data/processed/gse57338_meta.csv"); meta["group"] = meta["failing"].map({"yes": "Failing", "no": "Control"})
assert list(E.columns) == list(meta["sample"]) and meta["group"].notna().all()
Z = ((E.T - E.T.mean()) / E.T.std(ddof=0))
u = [g for g in UP if g in Z.columns]; d = [g for g in DOWN if g in Z.columns]
sc = Z[u].mean(axis=1) - Z[d].mean(axis=1)
y = (meta["group"] == "Failing").astype(int).values
def auc(a, b):
    r = rankdata(np.r_[a, b]); n1, n0 = len(a), len(b); return (r[:n1].sum() - n1 * (n1 + 1) / 2) / (n1 * n0)
a = auc(sc.values[y == 1], sc.values[y == 0])
p = mannwhitneyu(sc.values[y == 1], sc.values[y == 0], alternative="two-sided").pvalue
cov = pd.get_dummies(meta[["Age", "Sex"]], columns=["Sex"], drop_first=True).astype(float)
C = ((cov - cov.mean()) / (cov.std(ddof=0) + 1e-9)).values; s = ((sc - sc.mean()) / sc.std(ddof=0)).values; X1 = np.column_stack([s, C])
m0 = LogisticRegression(penalty=None, max_iter=5000).fit(C, y); m1 = LogisticRegression(penalty=None, max_iter=5000).fit(X1, y)
ll0 = -log_loss(y, m0.predict_proba(C)[:, 1], normalize=False); ll1 = -log_loss(y, m1.predict_proba(X1)[:, 1], normalize=False)
OR = float(np.exp(m1.coef_[0][0])); plr = float(chi2.sf(max(2 * (ll1 - ll0), 0), 1))
pool = [g for g in Z.columns if g not in ALLSIG]; Zv = Z[pool].values; rng = np.random.default_rng(42); ra = []
for _ in range(1000):
    idx = rng.choice(len(pool), len(u) + len(d), replace=False)
    z = Zv[:, idx[:len(u)]].mean(axis=1) - Zv[:, idx[len(u):]].mean(axis=1); ra.append(auc(z[y == 1], z[y == 0]))
ra = np.array(ra)
row = dict(analysis="GSE57338 failing vs non-failing", tissue="LV", n_case=int(y.sum()), n_ctrl=int((1 - y).sum()), n_up=len(u), n_down=len(d),
           AUC=round(a, 3), H0_lower_met=bool(a < 0.5), H1_higher_met=bool(a > 0.5), MWU_p=p, adj_OR_per_SD=round(OR, 3), adj_LR_p=plr,
           random_p_onesided_lower=(1 + (ra <= a).sum()) / 1001, random_p_onesided_higher=(1 + (ra >= a).sum()) / 1001,
           random_p_twosided=(1 + (abs(ra - .5) >= abs(a - .5)).sum()) / 1001)
out = pd.DataFrame([row]); out.to_csv("results/Table_validation_GSE57338.csv", index=False)
pd.DataFrame({"sample": sc.index, "group": meta["group"].values, "score": sc.values}).to_csv("results/gse57338_scores.csv", index=False)
pd.DataFrame({"GSE57338": ra}).to_csv("results/gse57338_random_auc.csv", index=False)
# BH over the declared secondary family (Mann-Whitney p): HFrEF RVS, HFrEF LVS, GSE57338
h = pd.read_csv("results/Table_validation_Hahn.csv"); h = h[h["analysis"].str.startswith("secondary")][["analysis", "AUC", "MWU_p"]]
fam = pd.concat([h, out[["analysis", "AUC", "MWU_p"]]], ignore_index=True)
pv = fam["MWU_p"].values; mm = len(pv); o = np.argsort(pv)
adj = np.empty(mm); adj[o] = np.minimum(np.minimum.accumulate((pv[o] * mm / np.arange(1, mm + 1))[::-1])[::-1], 1)
fam["BH_adj_p"] = adj; fam.to_csv("results/Table_secondary_BH.csv", index=False)
pd.set_option("display.width", 250); pd.set_option("display.max_columns", 30)
print(out.to_string(index=False)); print(fam.to_string(index=False)); print("done")
