# scripts/05b_collect.py -- run from the PROJECT ROOT:  .venv\Scripts\python.exe scripts\05b_collect.py
# Merges all results/perf_*.json into results/Table1_performance.csv
import glob, json, pandas as pd

rows = [json.load(open(f)) for f in glob.glob("results/perf_*.json")]
df = pd.DataFrame(rows)
null = df[df.model == "null"].groupby("target")["mae"].first().rename("null_mae")
df = df.merge(null, on="target")
df["mae_gain_vs_null_%"] = 100 * (df["null_mae"] - df["mae"]) / df["null_mae"]
# the guide's rule: useful only if rho > 0 AND MAE below the null MAE AND permutation p < 0.05
df["meets_rule"] = (df["model"] != "null") & (df["rho"] > 0) & (df["mae"] < df["null_mae"]) & (df["p_perm"] < 0.05)
df.loc[df.model == "null", "rho"] = float("nan")           # the null model's rho is not reported
df = df.sort_values(["target", "version", "model"])
df.to_csv("results/Table1_performance.csv", index=False)
pd.set_option("display.width", 200)
print(df[["target", "version", "model", "n_genes", "rho", "mae", "null_mae",
          "mae_gain_vs_null_%", "p_perm", "meets_rule"]].round(4).to_string(index=False))
print("rows:", len(df))
