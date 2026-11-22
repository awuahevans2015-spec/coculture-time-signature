# scripts/22_validate_hahn.py -- run from the PROJECT ROOT:  .venv\Scripts\python.exe scripts\22_validate_hahn.py
# Locked design: docs/expected_direction_Hahn2021_2026-11-21.txt. Score = mean z(up genes) - mean z(down genes), Version A day signature.
# Expected: LOWER in failing hearts, so AUC = P(score_failing > score_control) expected < 0.5. Exploratory, hypothesis-generating.
import glob, re, numpy as np, pandas as pd
from scipy.stats import rankdata, mannwhitneyu, chi2
from sklearn.linear_model import LogisticRegression
from sklearn.metrics import log_loss

sig = pd.read_csv("results/signature_log2_day_A.csv", keep_default_na=False)
sig = sig[sig["in_signature"].astype(str) == "True"].copy()
def sgn(x):
    s = str(x).strip().lower()
    try: return 1 if float(s) > 0 else -1
    except ValueError: pass
    return 1 if s.startswith(("up", "+", "pos")) else (-1 if s.startswith(("down", "-", "neg")) else 0)
sig["sign"] = sig["direction"].map(sgn)
print("direction values:", sig["direction"].value_counts().to_dict(), "| unparsed:", int((sig["sign"] == 0).sum()))
assert (sig["sign"] != 0).all(), "direction not parsed - send me this output"
UP = sig.loc[sig["sign"] > 0, "gene"].tolist(); DOWN = sig.loc[sig["sign"] < 0, "gene"].tolist(); ALLSIG = set(UP + DOWN)

def lab(d):
    s = str(d).lower()
    if "hfpef" in s: return "HFpEF"
    if "hfref" in s: return "HFrEF"
    if any(k in s for k in ("control", "donor", "ctrl", "normal")): return "Control"
    return "other"
def zmat(V, samples):
    X = V[samples].T
    return (X - X.mean()) / X.std(ddof=0)
def auc(case, ctrl):
    r = rankdata(np.r_[case, ctrl]); n1, n0 = len(case), len(ctrl)
    return (r[:n1].sum() - n1 * (n1 + 1) / 2) / (n1 * n0)
def adj_lr(sc, y, m):
    cov = pd.get_dummies(m[["Age", "Sex"]], columns=["Sex"], drop_first=True).astype(float)
    C = ((cov - cov.mean()) / (cov.std(ddof=0) + 1e-9)).values
    s = ((sc - sc.mean()) / sc.std(ddof=0)).values
    X1 = np.column_stack([s, C])
    m0 = LogisticRegression(penalty=None, max_iter=5000).fit(C, y)
    m1 = LogisticRegression(penalty=None, max_iter=5000).fit(X1, y)
    ll0 = -log_loss(y, m0.predict_proba(C)[:, 1], normalize=False)
    ll1 = -log_loss(y, m1.predict_proba(X1)[:, 1], normalize=False)
    return float(np.exp(m1.coef_[0][0])), float(chi2.sf(max(2 * (ll1 - ll0), 0), 1))

RAND = {}; SCORES = []
def run(V, meta, case, ctrl, up, down, name, tissue, excl=()):
    m = meta[meta["group"].isin([case, ctrl]) & ~meta["sample"].isin(excl)].copy()
    Z = zmat(V, list(m["sample"]))
    u = [g for g in up if g in Z.columns]; d = [g for g in down if g in Z.columns]
    sc = Z[u].mean(axis=1) - Z[d].mean(axis=1)
    m = m.set_index("sample"); sc = sc.loc[m.index]
    y = (m["group"] == case).astype(int).values
    a = auc(sc.values[y == 1], sc.values[y == 0])
    p = mannwhitneyu(sc.values[y == 1], sc.values[y == 0], alternative="two-sided").pvalue
    OR, plr = adj_lr(sc, y, m)
    pool = [g for g in Z.columns if g not in ALLSIG]; Zv = Z[pool].values
    rng = np.random.default_rng(42); ra = []
    for _ in range(1000):
        idx = rng.choice(len(pool), len(u) + len(d), replace=False)
        s = Zv[:, idx[:len(u)]].mean(axis=1) - Zv[:, idx[len(u):]].mean(axis=1)
        ra.append(auc(s[y == 1], s[y == 0]))
    ra = np.array(ra); RAND[name] = ra
    SCORES.append(pd.DataFrame({"analysis": name, "tissue": tissue, "sample": sc.index, "group": m["group"].values, "score": sc.values}))
    return dict(analysis=name, tissue=tissue, case=case, n_case=int(y.sum()), n_ctrl=int((1 - y).sum()), n_up=len(u), n_down=len(d),
                AUC=round(a, 3), direction_as_expected=bool(a < 0.5), MWU_p=p, adj_OR_per_SD=round(OR, 3), adj_LR_p=plr,
                random_p_onesided=(1 + (ra <= a).sum()) / 1001, random_p_twosided=(1 + (abs(ra - .5) >= abs(a - .5)).sum()) / 1001)

meta = {}; V = {}
for t in ["RVS", "LVS"]:
    meta[t] = pd.read_csv(f"data/processed/hahn2021_meta_{t}.csv"); meta[t]["sample"] = meta[t]["sample"].astype(str)
    meta[t]["group"] = meta[t]["disease"].map(lab)
    print(t, "disease labels:", meta[t]["disease"].value_counts().to_dict(), "->", meta[t]["group"].value_counts().to_dict())
    assert (meta[t]["group"] != "other").all(), "unmapped disease label - send me this output"
    V[t] = pd.read_csv(f"data/processed/hahn2021_vst_{t}.csv", index_col=0); V[t].columns = [str(c) for c in V[t].columns]

rows = []
rows.append(run(V["RVS"], meta["RVS"], "HFpEF", "Control", UP, DOWN, "PRIMARY HFpEF vs control RVS", "RVS"))
rows.append(run(V["RVS"], meta["RVS"], "HFrEF", "Control", UP, DOWN, "secondary HFrEF vs control RVS", "RVS"))
rows.append(run(V["LVS"], meta["LVS"], "HFrEF", "Control", UP, DOWN, "secondary HFrEF vs control LVS", "LVS"))
# Sensitivities (post hoc, labelled exploratory; primary contrast only)
rows.append(run(V["RVS"], meta["RVS"], "HFpEF", "Control", UP, DOWN, "sens A exclude flagged samples", "RVS", excl=["7709", "7546", "7715", "7559", "7718"]))
f = glob.glob("data/hahn2021_v1.0.1/**/ExcludeRBCgenes.xlsx", recursive=True)
if f:
    cells = {str(v).strip() for df in pd.read_excel(f[0], sheet_name=None).values() for v in df.to_numpy().ravel() if pd.notna(v)}
    rbc = ALLSIG & cells; print("sens B: RBC list cells", len(cells), "| signature genes removed:", len(rbc))
    rows.append(run(V["RVS"], meta["RVS"], "HFpEF", "Control", [g for g in UP if g not in rbc], [g for g in DOWN if g not in rbc], "sens B remove RBC genes", "RVS"))
else: print("sens B skipped: ExcludeRBCgenes.xlsx not found")
en = pd.read_csv("results/enrichment_log2_day.csv", keep_default_na=False); en["p.adjust"] = pd.to_numeric(en["p.adjust"], errors="coerce")
top = en.sort_values("p.adjust")["ID"].head(10).tolist()
gs = pd.read_csv("data/genesets/msigdb_sets_used.csv")
setg = set(gs.loc[gs["gs_name"].isin(top), "gene_symbol"])
ribo = {g for g in ALLSIG if re.match(r"^(MRPL|MRPS|RPL|RPS)\d|^RPLP\d|^RPSA$", g)}
rem = ribo | (ALLSIG & setg); print("sens C: ribosomal", len(ribo), "| top-10-set genes in signature", len(ALLSIG & setg), "| total removed", len(rem))
rows.append(run(V["RVS"], meta["RVS"], "HFpEF", "Control", [g for g in UP if g not in rem], [g for g in DOWN if g not in rem], "sens C remove ribosomal and top-10-set genes", "RVS"))

out = pd.DataFrame(rows); out.to_csv("results/Table_validation_Hahn.csv", index=False)
pd.concat(SCORES).to_csv("results/hahn_scores.csv", index=False)
pd.DataFrame({k: v for k, v in RAND.items()}).to_csv("results/hahn_random_auc.csv", index=False)
pd.set_option("display.width", 250); pd.set_option("display.max_columns", 30)
print(out.drop(columns=["case"]).to_string(index=False))
print("done")
