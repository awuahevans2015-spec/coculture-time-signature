import hashlib, glob, datetime
files = sorted(glob.glob("data/processed/hahn2021_*.csv"))
with open("docs/checksums_hahn_processed.txt", "w") as o:
    o.write(f"# SHA-256, written {datetime.date.today()} before any group comparison\n")
    for f in files:
        h = hashlib.sha256(open(f, "rb").read()).hexdigest()
        o.write(f"{h}  {f}\n"); print(h[:16], f)
with open("docs/code_changes.md", "a") as o:
    o.write("\n## Hahn 2021 preparation (scripts 16-19)\n"
            "- Excel count sheets converted to CSV (16); layout and header checked (17, 18).\n"
            "- Symbol map = geneid and external_gene_name columns of DEG_allresultsRVS.csv only (13,266 genes); no group result was read.\n"
            "- Labelled samples only (RVS 95, LVS 50); all-zero genes dropped; duplicate symbols keep the highest mean count.\n"
            "- Signature coverage 131/161 in both tissues (above the 50% locked rule); the 30 missing genes were not in the authors' tested-gene table.\n")
print("logged")
