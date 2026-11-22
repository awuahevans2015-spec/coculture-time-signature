# scripts/21_vst_hahn.py -- run from the PROJECT ROOT:  .venv\Scripts\python.exe scripts\21_vst_hahn.py
# Variance-stabilizing transform within tissue, blind to group (intercept-only design), PyDESeq2 0.6.0.
# Deviation from plan: R DESeq2 could not load (Application Control blocked rlang.dll). No group information is used.
import pandas as pd
from pydeseq2.dds import DeseqDataSet
for t in ["RVS", "LVS"]:
    cnt = pd.read_csv(f"data/processed/hahn2021_counts_{t}_symbols.csv", index_col=0)
    cnt.columns = [str(c) for c in cnt.columns]
    meta = pd.read_csv(f"data/processed/hahn2021_meta_{t}.csv"); meta["sample"] = meta["sample"].astype(str)
    assert list(cnt.columns) == list(meta["sample"]), "sample order mismatch"
    counts = cnt.T                                            # samples x genes
    md = pd.DataFrame({"dummy": "all"}, index=counts.index)   # no group column is passed
    dds = DeseqDataSet(counts=counts, metadata=md, design="~1", quiet=True)
    dds.vst(use_design=False)
    v = pd.DataFrame(dds.layers["vst_counts"], index=counts.index, columns=counts.columns).T
    v.to_csv(f"data/processed/hahn2021_vst_{t}.csv")
    print(t, "VST done:", v.shape[0], "genes x", v.shape[1], "samples | min/max:", round(v.values.min(), 2), round(v.values.max(), 2))
print("done")
