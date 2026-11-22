# scripts/16_xlsx_to_csv.py -- run from the PROJECT ROOT:  .venv\Scripts\python.exe scripts\16_xlsx_to_csv.py
# Format conversion only (xlsx -> csv) for the Hahn 2021 count files; no values are changed.
import pandas as pd
src = "data/hahn2021_v1.0.1/HFpEF_2020-v1.0.1/baderzone-HFpEF_2020-71fc3a4/Data/"
for tissue in ["RVS", "LVS"]:
    df = pd.read_excel(src + f"RawReads_{tissue}_Ensembl_geneName.xlsx", sheet_name="Reads")
    print(tissue, "shape:", df.shape, "| first columns:", list(df.columns[:8]), flush=True)
    df.to_csv(f"data/processed/hahn2021_RawReads_{tissue}.csv", index=False)
print("done")
