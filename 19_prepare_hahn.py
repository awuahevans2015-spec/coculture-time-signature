# scripts/19_prepare_hahn.py -- run from the PROJECT ROOT:  .venv\Scripts\python.exe scripts\19_prepare_hahn.py
# Builds symbol-labelled count tables and sample metadata for Hahn 2021 (RVS and LVS). No group comparison is run.
# Gene map: ONLY the geneid and external_gene_name columns of DEG_allresultsRVS.csv are read.
import numpy as np, pandas as pd
p = "data/hahn2021_v1.0.1/HFpEF_2020-v1.0.1/baderzone-HFpEF_2020-71fc3a4/Data/"
mp = pd.read_csv(p + "DEG_allresultsRVS.csv", usecols=["geneid", "external_gene_name"])
mp = mp.dropna(); mp = mp[mp["external_gene_name"].astype(str).str.len() > 0].drop_duplicates("geneid")
samples = pd.read_csv(p + "SamplesTable_ALLpatients.csv"); samples["samplenames"] = samples["samplenames"].astype(str)
clin = pd.read_csv("data/hahn2021_v1.0.1/ClinicalCharacteristicsRNAseqSamples.csv")
clin["SequenceRunId"] = clin["SequenceRunId"].astype(str)
sig = pd.read_csv("results/signature_log2_day_A.csv", keep_default_na=False)
sig_genes = sig[sig["in_signature"].astype(str) == "True"]["gene"].tolist()
print("signature genes:", len(sig_genes))
for t in ["RVS", "LVS"]:
    c = pd.read_csv(f"data/processed/hahn2021_RawReads_{t}.csv")
    c.columns = ["geneid"] + [str(x) for x in c.columns[1:]]
    keep = [k for k in samples.loc[samples["tissue"] == t, "samplenames"] if k in c.columns]
    c = c[["geneid"] + keep].merge(mp, on="geneid", how="inner")
    vals = c[keep].to_numpy(float)
    print(t, "| labelled samples used:", len(keep), "| genes with a symbol:", len(c),
          "| all counts integer-valued:", bool(np.allclose(vals, np.round(vals))))
    c = c[vals.sum(axis=1) > 0].copy()                                   # drop genes with zero counts in all samples
    c["mean_count"] = c[keep].mean(axis=1)
    c = c.sort_values("mean_count", ascending=False).drop_duplicates("external_gene_name")   # duplicated symbol: keep highest mean
    c = c.set_index("external_gene_name")[keep]
    present = [g for g in sig_genes if g in c.index]
    print(t, "| genes after filters:", len(c), "| signature genes present:", len(present), "of", len(sig_genes),
          "| missing (first 10):", [g for g in sig_genes if g not in c.index][:10])
    c.round().astype(int).to_csv(f"data/processed/hahn2021_counts_{t}_symbols.csv")
    m = samples[samples["tissue"] == t][["samplenames", "disease", "tissue"]].rename(columns={"samplenames": "sample"})
    m = m[m["sample"].isin(keep)].merge(clin[["SequenceRunId", "Age", "Sex"]], left_on="sample", right_on="SequenceRunId", how="left")
    m = m.drop(columns="SequenceRunId"); m.to_csv(f"data/processed/hahn2021_meta_{t}.csv", index=False)
    print(t, "| metadata rows:", len(m), "| missing Age:", int(m["Age"].isna().sum()), "| missing Sex:", int(m["Sex"].isna().sum()))
print("done")
