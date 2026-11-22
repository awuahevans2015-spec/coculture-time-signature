# scripts/06b_loco_summary.py -- run from the PROJECT ROOT:  .venv\Scripts\python.exe scripts\06b_loco_summary.py
# Prints results/Table_LOCO.csv (leave-one-condition-out). Read-only; changes nothing.
import pandas as pd
pd.set_option("display.width", 220); pd.set_option("display.max_columns", 50); pd.set_option("display.max_rows", 200)
df = pd.read_csv("results/Table_LOCO.csv", keep_default_na=False)   # keep the model name "null" as text, not NaN
print("columns:", list(df.columns))
print("rows:", len(df))
print(df.round(4).to_string(index=False))
