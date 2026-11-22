Code changes vs the guide (scikit-learn 1.9.1, Python 3.13.14):
ElasticNetCV(n_alphas=50) -> ElasticNetCV(alphas=50) because n_alphas was removed (same 50-value grid).
TopVar base-class order changed to (TransformerMixin, BaseEstimator) (no change in behavior).
Scripts run from the project root; script 05 restructured with a main() guard for parallel runs on Windows.
RandomForestRegressor import made lazy in mlcore.py because an application-control policy blocked a compiled scikit-learn file (sklearn.ensemble histogram module) on this computer. Random forest status: undecided; to be recorded before any results.
11_overlap: FDR computed over Version A tests only (Version B signatures come from the panel, so their panel overlap is circular). Decided before viewing Table_overlap.csv.

## Hahn 2021 preparation (scripts 16-19)
- Excel count sheets converted to CSV (16); layout and header checked (17, 18).
- Symbol map = geneid and external_gene_name columns of DEG_allresultsRVS.csv only (13,266 genes); no group result was read.
- Labelled samples only (RVS 95, LVS 50); all-zero genes dropped; duplicate symbols keep the highest mean count.
- Signature coverage 131/161 in both tissues (above the 50% locked rule); the 30 missing genes were not in the authors' tested-gene table.
