import hashlib, glob, datetime
files = sorted(glob.glob("data/human/GSE57338/*"))
with open("docs/checksums_GSE57338_raw.txt", "w") as o:
    o.write(f"# SHA-256 of downloaded GSE57338 files, {datetime.date.today()}, written before the files were opened\n")
    for f in files:
        h = hashlib.sha256(open(f, "rb").read()).hexdigest()
        o.write(f"{h}  {f}\n"); print(h[:16], f)
print("done")
