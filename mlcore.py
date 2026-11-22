# scripts/mlcore.py -- shared machine-learning code (guide Step 5a)
# Changes from the guide:
#   ElasticNetCV(n_alphas=50) -> ElasticNetCV(alphas=50)   (n_alphas was removed in scikit-learn 1.9)
#   class TopVar(BaseEstimator, TransformerMixin) -> class TopVar(TransformerMixin, BaseEstimator)
#   RandomForestRegressor is imported only when model "rf" is requested
import glob, numpy as np, pandas as pd
from scipy.stats import spearmanr
from sklearn.base import BaseEstimator, TransformerMixin
from sklearn.pipeline import Pipeline
from sklearn.preprocessing import StandardScaler
from sklearn.linear_model import ElasticNetCV, RidgeCV
from sklearn.dummy import DummyRegressor
from sklearn.model_selection import LeaveOneOut, LeaveOneGroupOut, cross_val_predict

SEED = 42
MARKERS = ["PECAM1", "CDH5", "KDR", "PDGFRB", "CSPG4", "RGS5"]


class TopVar(TransformerMixin, BaseEstimator):
    """Keep the k most variable genes, learned from the TRAINING fold only."""
    def __init__(self, k=3000):
        self.k = k

    def fit(self, X, y=None):
        self.idx_ = np.argsort(X.var(axis=0))[::-1][: self.k]
        return self

    def transform(self, X):
        return X[:, self.idx_]


def load(kind="primary"):
    X = pd.read_csv(f"data/processed/vst_{kind}.csv", index_col=0).T   # samples x genes
    meta = pd.read_csv(f"data/processed/meta_{kind}.csv", index_col=0).loc[X.index]
    return X, meta


def load_panel(X, which="core"):
    f = sorted(glob.glob(f"data/genesets/MVD_panel_{which}_*.txt"))[-1]   # latest frozen file
    genes = [g.strip() for g in open(f) if g.strip()]
    return [g for g in genes if g in X.columns]


def select(X, meta, target, version, panel, k=3000):
    ok = meta[target].notna().values          # drop samples with no value (e.g. Matrigel PEG)
    Xs = X.loc[ok] if version == "A" else X.loc[ok, panel]
    return Xs, meta.loc[ok], meta.loc[ok, target].values


def make_model(name, version, k=3000):
    steps = []
    if version == "A":
        steps.append(("var", TopVar(k)))
    if name == "enet":
        est = ElasticNetCV(l1_ratio=[0.1, 0.5, 0.9], alphas=50, cv=5, max_iter=20000, random_state=SEED)
    elif name == "ridge":
        est = RidgeCV(alphas=np.logspace(-2, 5, 30))
    elif name == "rf":
        # imported only when needed: importing sklearn.ensemble loads compiled files
        # that an application-control policy may block on some computers
        from sklearn.ensemble import RandomForestRegressor
        est = RandomForestRegressor(n_estimators=500, random_state=SEED, n_jobs=1)
    elif name == "null":
        est = DummyRegressor(strategy="mean")
    else:
        raise ValueError(name)
    steps += [("scale", StandardScaler()), (name, est)]
    return Pipeline(steps)


def loo_predict(name, version, X, y, k=3000, groups=None):
    cv = LeaveOneOut() if groups is None else LeaveOneGroupOut()
    return cross_val_predict(make_model(name, version, k), X, y, cv=cv, groups=groups)


def rho_mae(y, pred):
    rho = spearmanr(y, pred)[0] if np.std(pred) > 0 else 0.0
    return float(np.nan_to_num(rho)), float(np.mean(np.abs(y - pred)))
