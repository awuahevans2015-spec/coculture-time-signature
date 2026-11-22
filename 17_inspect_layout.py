# scripts/17_inspect_layout.py -- run from the PROJECT ROOT:  .venv\Scripts\python.exe scripts\17_inspect_layout.py
# Layout check only: group counts, ID matching and gene labels. Prints no expression values.
import pandas as pd
p = "data/hahn2021_v1.0.1/HFpEF_2020-v1.0.1/baderzone-HFpEF_2020-71fc3a4/Data/"
s = pd.read_csv(p + "SamplesTable_ALLpatients.csv"); s["samplenames"] = s["samplenames"].astype(str)
print(s.head(3).to_string()); print(pd.crosstab(s["disease"], s["tissue"]))
for t in ["RVS", "LVS"]:
    f = f"data/processed/hahn2021_RawReads_{t}.csv"
    cols = [str(c) for c in pd.read_csv(f, nrows=0).columns][1:]
    ids = set(s.loc[s["tissue"] == t, "samplenames"])
    print(t, "count columns:", len(cols), "| matched to samples table:", len(set(cols) & ids),
          "| columns not in table:", sorted(set(cols) - set(s["samplenames"]))[:10],
          "| in table but no counts:", sorted(ids - set(cols))[:10])
    g = pd.read_csv(f, usecols=[0]); print(t, "rows:", len(g), "| first gene labels:", g.iloc[:5, 0].tolist())
cl = pd.read_csv("data/hahn2021_v1.0.1/ClinicalCharacteristicsRNAseqSamples.csv")
print("clinical table:", cl.shape, list(cl.columns)[:15])
