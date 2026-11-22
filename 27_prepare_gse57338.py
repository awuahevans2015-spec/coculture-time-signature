# scripts/27_prepare_gse57338.py -- run from the PROJECT ROOT. Builds symbol-labelled expression + metadata. No group comparison.
import gzip, io, numpy as np, pandas as pd
p = "data/human/GSE57338/GSE57338_series_matrix.txt.gz"
chars = {}; geo = None; tbl = []; inside = False
with gzip.open(p, "rt", encoding="utf-8", errors="replace") as f:
    for line in f:
        if line.startswith("!series_matrix_table_begin"): inside = True; continue
        if line.startswith("!series_matrix_table_end"): break
        if inside: tbl.append(line)
        elif line.startswith("!Sample_geo_accession"): geo = [x.strip('"') for x in line.rstrip("\n").split("\t")[1:]]
        elif line.startswith("!Sample_characteristics_ch1"):
            r = [x.strip('"') for x in line.rstrip("\n").split("\t")[1:]]
            chars[r[0].split(":")[0].strip()] = [x.split(":", 1)[1].strip() for x in r]
X = pd.read_csv(io.StringIO("".join(tbl)), sep="\t", index_col=0); X.index = X.index.astype(str)
assert list(X.columns) == geo, "sample order mismatch"
meta = pd.DataFrame({"sample": geo, "failing": chars["heart failure"], "disease_status": chars["disease status"],
                     "Sex": chars["gender"], "Age": pd.to_numeric(chars["age"], errors="coerce")})
v = X.values.astype(float)
print("samples:", X.shape[1], "| probes:", X.shape[0], "| NaN values:", int(np.isnan(v).sum()),
      "| expression min/median/max:", round(np.nanmin(v), 2), round(np.nanmedian(v), 2), round(np.nanmax(v), 2))
lines = open("data/human/GSE57338/GPL11532-32230.txt", encoding="utf-8", errors="replace").read().splitlines()
s = next(i for i, l in enumerate(lines) if not l.startswith("#"))
A = pd.read_csv(io.StringIO("\n".join(lines[s:])), sep="\t", dtype=str, quoting=3)
print("platform columns:", list(A.columns)[:12])
def syms(c):
    if not isinstance(c, str) or c.strip() in ("", "---"): return set()
    out = set()
    for e in c.split(" /// "):
        q = e.split(" // ")
        if len(q) > 1 and q[1].strip() not in ("", "---"): out.add(q[1].strip())
    return out
A["S"] = A["gene_assignment"].map(syms); A["n"] = A["S"].map(len)
print("transcript clusters:", len(A), "| one symbol:", int((A.n == 1).sum()), "| ambiguous (>1):", int((A.n > 1).sum()), "| none:", int((A.n == 0).sum()))
mp = A[A.n == 1].copy(); mp["symbol"] = mp["S"].map(lambda z: next(iter(z))); mp = mp[["ID", "symbol"]]; mp["ID"] = mp["ID"].astype(str)
Xr = X.copy(); Xr["ID"] = Xr.index; E = Xr.merge(mp, on="ID")
E["mean_expr"] = E[geo].mean(axis=1)
E = E.sort_values("mean_expr", ascending=False).drop_duplicates("symbol").set_index("symbol")[geo]
sig = pd.read_csv("results/signature_log2_day_A.csv", keep_default_na=False)
sg = sig[sig["in_signature"].astype(str) == "True"]["gene"].tolist()
pres = [g for g in sg if g in E.index]
print("genes after mapping:", len(E), "| signature genes present:", len(pres), "of", len(sg), "| missing (first 10):", [g for g in sg if g not in E.index][:10])
print("missing Age:", int(meta["Age"].isna().sum()), "| missing Sex:", int(meta["Sex"].isna().sum()))
E.to_csv("data/processed/gse57338_expr_symbols.csv"); meta.to_csv("data/processed/gse57338_meta.csv", index=False)
print("done")
