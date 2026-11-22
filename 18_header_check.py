# scripts/18_header_check.py -- run from the PROJECT ROOT:  .venv\Scripts\python.exe scripts\18_header_check.py
# Prints only the column NAMES of the authors' result tables (no values), to find the Ensembl ID and gene-name columns.
import pandas as pd
p = "data/hahn2021_v1.0.1/HFpEF_2020-v1.0.1/baderzone-HFpEF_2020-71fc3a4/Data/"
for name in ["DEG_allresultsRVS.csv", "DEG_allresults_HFpEFgroups.csv"]:
    cols = list(pd.read_csv(p + name, nrows=0).columns)
    print(name, "| n columns:", len(cols)); print(cols)
