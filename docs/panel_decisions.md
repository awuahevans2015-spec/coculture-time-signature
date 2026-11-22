Open Targets panel decisions, written 2026-11-14 before any modeling
Release: UI v26.9.1, API 26.9.0.
Heart failure: MONDO_0005252, associated targets table (5,544 targets).
HFpEF: search returned 6 GWAS studies and 0 disease entries; no separate HFpEF disease page exists. Not used as a separate source. HFpEF GWAS studies are mapped by Open Targets to heart failure.
Coronary microvascular dysfunction: searched in Open Targets (UI v26.9.1); 'No search result found'. Not used as a source.
MSigDB set 'GOBP_ANGIOGENESIS' not found in this msigdbr release (generic term absent). Replaced by GOBP_SPROUTING_ANGIOGENESIS (closest by name; narrower). Chosen by name before the panel was frozen, without reference to any results.
Cardiac fibrosis: closest available entry used = interstitial cardiac fibrosis (HP_0031329), 171 targets. Decided before any modeling.
Final panel (built 2026-11-14): core 1212 genes; with markers 1218; pathway genes 1363.
ML run plan, recorded 2026-11-15 before the full runs:
Models: null, ridge, elastic net; feature sets A (top 3000 variable genes inside folds) and B (frozen panel); targets log2_day and peg_c; leave-one-sample-out; 1000 permutations; seed 42.
Random forest: not run on this computer (application-control policy blocked a compiled scikit-learn file). Handling to be decided before any results: run on another computer if available, otherwise report as not performed.
