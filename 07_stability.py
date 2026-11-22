# scripts/07_stability.py -- run from the PROJECT ROOT:  .venv\Scripts\python.exe scripts\07_stability.py --jobs 7
# Stability selection: refit the elastic net on 200 random subsets (8 of 11 samples); keep genes picked in >= 70% of
# fits with the same sign every time. Outputs results/signature_TARGET_VERSION.csv for both targets and versions.
import argparse, numpy as np, pandas as pd
from joblib import Parallel, delayed
from mlcore import *

FRAC, THRESH = 8, 0.70


def one_fit(r, X, y, version):
    idx = np.random.default_rng(SEED + 100 + r).choice(len(y), size=FRAC, replace=False)
    m = make_model("enet", version)
    m.set_params(enet__cv=4)                       # smaller inner CV for 8 samples
    m.fit(X[idx], y[idx])
    coef = m.named_steps["enet"].coef_
    cols = m.named_steps["var"].idx_ if version == "A" else np.arange(X.shape[1])
    full = np.zeros(X.shape[1]); full[cols] = coef
    return full


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--nrep", type=int, default=200)
    ap.add_argument("--jobs", type=int, default=-1)
    a = ap.parse_args()
    X_all, meta = load("primary")
    panel = load_panel(X_all, "core")
    for target in ["log2_day", "peg_c"]:
        for version in ["A", "B"]:
            X, m, y = select(X_all, meta, target, version, panel)
            Xv = X.values
            C = np.array(Parallel(n_jobs=a.jobs)(delayed(one_fit)(r, Xv, y, version) for r in range(a.nrep)))
            freq = (C != 0).mean(axis=0)
            sign_ok = np.array([np.unique(np.sign(C[:, g][C[:, g] != 0])).size <= 1 for g in range(C.shape[1])])
            direction = np.sign(np.where(C != 0, C, 0).sum(axis=0))
            df = pd.DataFrame({"gene": X.columns, "freq": freq, "direction": direction, "sign_consistent": sign_ok})
            df["in_signature"] = (df.freq >= THRESH) & df.sign_consistent
            df = df.sort_values("freq", ascending=False)
            df.to_csv(f"results/signature_{target}_{version}.csv", index=False)
            print(target, version, "signature size:", int(df.in_signature.sum()),
                  "| genes ever selected:", int((freq > 0).sum()), flush=True)


if __name__ == "__main__":
    main()
