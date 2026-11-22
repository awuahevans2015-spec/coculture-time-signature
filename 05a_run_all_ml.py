# scripts/05a_run_all_ml.py -- runs every model x feature-set x target combination from guide Step 5b.
# Run from the PROJECT ROOT:  .venv\Scripts\python.exe scripts\05a_run_all_ml.py --nperm 1000 --jobs 4
# Safe to stop and restart: finished combinations are skipped. Progress is logged to results/run_log.txt.
import argparse, json, os, subprocess, sys, time, datetime

ap = argparse.ArgumentParser()
ap.add_argument("--nperm", type=int, default=1000)
ap.add_argument("--jobs", type=int, default=4)
ap.add_argument("--models", default="null,ridge,enet")   # random forest (rf) not included: blocked on this computer
a = ap.parse_args()

HERE = os.path.dirname(os.path.abspath(__file__))
script = os.path.join(HERE, "05_ml_cv.py")
os.makedirs("results", exist_ok=True)

def log(msg):
    line = f"{datetime.datetime.now():%Y-%m-%d %H:%M:%S}  {msg}"
    print(line, flush=True)
    with open("results/run_log.txt", "a") as fh:
        fh.write(line + "\n")

def done(tag, model):
    f = f"results/perf_{tag}.json"
    if not os.path.exists(f):
        return False
    d = json.load(open(f))
    return True if model == "null" else d.get("nperm") == a.nperm

# cheap models first, then feature set B (frozen panel), then A (top-variable genes)
for model in a.models.split(","):
    for version in ["B", "A"]:
        for target in ["log2_day", "peg_c"]:
            tag = f"{target}_{version}_{model}"
            if done(tag, model):
                log(f"skip   {tag} (already finished)")
                continue
            log(f"start  {tag}")
            t0 = time.time()
            r = subprocess.run([sys.executable, script, "--target", target, "--version", version,
                                "--model", model, "--nperm", str(a.nperm), "--jobs", str(a.jobs)])
            mins = (time.time() - t0) / 60
            log(f"{'done ' if r.returncode == 0 else 'FAILED'} {tag}  ({mins:.1f} min)")
log("all requested combinations processed")
