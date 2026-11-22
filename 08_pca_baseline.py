# scripts/08_pca_baseline.py -- run from the PROJECT ROOT:  .venv\Scripts\python.exe scripts\08_pca_baseline.py --jobs 7
# Baseline: top 3000 variable genes -> scale -> PCA (2-5 components, tuned by 5-fold CV) -> ridge,
# with the same leave-one-sample-out validation and 1000 permutations. Outputs results/Table_PCA_baseline.csv
# The PCA step is a small SVD-based transformer (same projections as sklearn's PCA, up to sign), used so the
# script does not need sklearn.decomposition (some compiled scikit-learn files are blocked on this computer).
import argparse, numpy as np, pandas as pd
from joblib import Parallel, delayed
from sklearn.base import BaseEstimator, TransformerMixin
from sklearn.model_selection import GridSearchCV, KFold, LeaveOneOut, cross_val_predict
from sklearn.pipeline import Pipeline
from sklearn.preprocessing import StandardScaler
from sklearn.linear_model import RidgeCV
from mlcore import *


class SvdPCA(TransformerMixin, BaseEstimator):
    def __init__(self, n_components=2):
        self.n_components = n_components

    def fit(self, X, y=None):
        self.mean_ = X.mean(axis=0)
        _, _, vt = np.linalg.svd(X - self.mean_, full_matrices=False)
        self.components_ = vt[: self.n_components]
        return self

    def transform(self, X):
        return (X - self.mean_) @ self.components_.T


def pca_model():
    pipe = Pipeline([("var", TopVar(3000)), ("scale", StandardScaler()), ("pca", SvdPCA()),
                     ("reg", RidgeCV(alphas=np.logspace(-2, 3, 20)))])
    return GridSearchCV(pipe, {"pca__n_components": [2, 3, 4, 5]}, cv=KFold(5, shuffle=True, random_state=SEED))


def one_perm(i, Xv, y):
    yp = np.random.default_rng(SEED + 1 + i).permutation(y)
    return rho_mae(yp, cross_val_predict(pca_model(), Xv, yp, cv=LeaveOneOut()))[0]


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--nperm", type=int, default=1000)
    ap.add_argument("--jobs", type=int, default=-1)
    a = ap.parse_args()
    X_all, meta = load("primary")
    rows = []
    for target in ["log2_day", "peg_c"]:
        X, m, y = select(X_all, meta, target, "A", [])
        Xv = X.values
        pred = cross_val_predict(pca_model(), Xv, y, cv=LeaveOneOut())
        rho, mae = rho_mae(y, pred)
        null = np.array(Parallel(n_jobs=a.jobs)(delayed(one_perm)(i, Xv, y) for i in range(a.nperm)))
        p = (1 + np.sum(null >= rho)) / (1 + len(null))
        rows.append(dict(target=target, model="PCA+ridge", rho=rho, mae=mae, nperm=a.nperm, p_perm=p))
        print(rows[-1], flush=True)
    pd.DataFrame(rows).to_csv("results/Table_PCA_baseline.csv", index=False)
    print(pd.DataFrame(rows))


if __name__ == "__main__":
    main()
