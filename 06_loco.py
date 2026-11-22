# scripts/06_loco.py -- run from the PROJECT ROOT:  .venv\Scripts\python.exe scripts\06_loco.py
# Leave-one-condition-out: hold out a whole day (or a whole crosslinking level), train on the rest, predict it.
# Random forest is not included (blocked on this computer); see docs/panel_decisions.md.
import pandas as pd
from mlcore import *

MODELS = ["enet", "ridge", "null"]
X_all, meta = load("primary")
panel = load_panel(X_all, "core")
rows = []
for target in ["log2_day", "peg_c"]:
    for gname, gcol in [("by_day", "day"), ("by_peg", "peg_pct")]:
        for version in ["A", "B"]:
            for model in MODELS:
                X, m, y = select(X_all, meta, target, version, panel)
                groups = m[gcol].values
                pred = loo_predict(model, version, X.values, y, groups=groups)
                rho, mae = rho_mae(y, pred)
                rows.append(dict(target=target, hold_out=gname, version=version, model=model,
                                 rho=(float("nan") if model == "null" else rho), mae=mae))
                for g in sorted(set(groups)):                  # error for each held-out group
                    idx = groups == g
                    rows.append(dict(target=target, hold_out=f"{gname}={g}", version=version, model=model,
                                     rho=float("nan"), mae=float(abs(y[idx] - pred[idx]).mean())))
                print("done", target, gname, version, model, flush=True)
out = pd.DataFrame(rows)
out.to_csv("results/Table_LOCO.csv", index=False)
pd.set_option("display.width", 200)
print(out.query("hold_out in ['by_day','by_peg']").round(4).to_string(index=False))
