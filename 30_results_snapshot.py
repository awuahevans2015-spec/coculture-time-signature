import glob, os, pandas as pd
skip = ("pred_", "hahn_scores", "gse57338_scores", "hahn_random", "gse57338_random")
out = []
for f in sorted(glob.glob("results/*.csv")):
    b = os.path.basename(f)
    if b.startswith(skip): continue
    try: d = pd.read_csv(f, keep_default_na=False)
    except Exception as e: out.append(f"## {b}: unreadable ({e})"); continue
    if b.startswith("signature_") and "in_signature" in d.columns:
        out.append(f"## {b}: {len(d)} genes, in_signature = {int((d['in_signature'].astype(str) == 'True').sum())}"); continue
    out.append(f"## {b}  shape={d.shape}"); out.append(d.head(40).to_string(index=False, max_colwidth=40))
open("docs/results_snapshot.txt", "w", encoding="utf-8").write("\n".join(out))
print(len(out), "blocks written to docs/results_snapshot.txt")
