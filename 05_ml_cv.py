# scripts/05_ml_cv.py -- run from the PROJECT ROOT (Raw_Dat2), e.g.
#   .venv\Scripts\python.exe scripts\05_ml_cv.py --target log2_day --version A --model enet --nperm 20 --jobs 4
import argparse, json, numpy as np, pandas as pd
from joblib import Parallel, delayed
from mlcore import *


def one_perm(i, model, version, Xv, y):
    yp = np.random.default_rng(SEED + 1 + i).permutation(y)
    return rho_mae(yp, loo_predict(model, version, Xv, yp))[0]


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--target", required=True)     # log2_day or peg_c
    ap.add_argument("--version", required=True)    # A or B
    ap.add_argument("--model", required=True)      # enet, ridge, rf, null
    ap.add_argument("--nperm", type=int, default=1000)
    ap.add_argument("--jobs", type=int, default=-1)
    a = ap.parse_args()

    X_all, meta = load("primary")
    panel = load_panel(X_all, "core")
    X, m, y = select(X_all, meta, a.target, a.version, panel)
    Xv = X.values
    pred = loo_predict(a.model, a.version, Xv, y)
    rho, mae = rho_mae(y, pred)
    out = dict(target=a.target, version=a.version, model=a.model, n=len(y), n_genes=int(Xv.shape[1]),
               rho=rho, mae=mae, nperm=0, p_perm=None)
    tag = f"{a.target}_{a.version}_{a.model}"
    if a.model != "null" and a.nperm > 0:
        null = np.array(Parallel(n_jobs=a.jobs)(
            delayed(one_perm)(i, a.model, a.version, Xv, y) for i in range(a.nperm)))
        out.update(nperm=a.nperm, p_perm=float((1 + np.sum(null >= rho)) / (1 + a.nperm)))
        np.save(f"results/null_{tag}.npy", null)
    json.dump(out, open(f"results/perf_{tag}.json", "w"))
    pd.DataFrame({"sample": X.index, "observed": y, "predicted": pred}).to_csv(f"results/pred_{tag}.csv", index=False)
    print(out)


if __name__ == "__main__":
    main()
