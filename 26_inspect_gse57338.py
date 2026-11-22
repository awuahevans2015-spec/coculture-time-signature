# scripts/26_inspect_gse57338.py -- layout check only. No expression values and no group comparison are printed.
import gzip, collections
p = "data/human/GSE57338/GSE57338_series_matrix.txt.gz"
meta = collections.OrderedDict(); inside = False; n = 0; header = None; ids = []
with gzip.open(p, "rt", encoding="utf-8", errors="replace") as f:
    for line in f:
        if line.startswith("!series_matrix_table_begin"): inside = True; continue
        if line.startswith("!series_matrix_table_end"): break
        if inside:
            n += 1; parts = line.rstrip("\n").split("\t")
            if n == 1: header = parts
            elif n <= 4: ids.append(parts[0])
        elif line.startswith("!Sample_"):
            k, *v = line.rstrip("\n").split("\t"); meta.setdefault(k, []).append([x.strip('"') for x in v])
print("samples in table:", len(header) - 1, "| table rows (probes):", n - 1)
print("first sample ids:", header[1:4], "| first probe ids:", ids[:3])
for k, rows in meta.items():
    if k in ("!Sample_title", "!Sample_characteristics_ch1", "!Sample_geo_accession"):
        for r in rows:
            c = collections.Counter(r)
            print(k, "| unique:", len(c), "|", dict(c.most_common(6)) if len(c) <= 8 else list(c.items())[:4])
g = "data/human/GSE57338/GPL11532-32230.txt"
lines = open(g, encoding="utf-8", errors="replace").read().splitlines()
print("platform file lines:", len(lines))
for l in lines:
    if l.startswith("#"): print(l[:160])
b = [i for i, l in enumerate(lines) if l.startswith("!platform_table_begin")]
s = b[0] + 1 if b else 0
print("platform header:", lines[s][:300])
for l in lines[s + 1:s + 4]: print(l[:300])
print("done")
