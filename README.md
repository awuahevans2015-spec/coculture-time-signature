# Co-culture culture-time signature (GEO GSE93511)

Analysis code for the manuscript:

"Machine-learning identification of a culture-time signature in 3D PEG-RGD
endothelial-pericyte co-culture, with limited predictability of crosslinking and
exploratory cross-dataset analysis in human heart failure"

Author: Evans Awuah, University of Bridgeport

## What this repository contains
- Scripts are in the top level of the repository, numbered in the order they were run
  (04 to 30). mlcore.py is a shared module used by the machine-learning scripts.
- results/  Result tables produced by the scripts.
- figures/  Figure files produced by the scripts.
- docs/     Provenance notes, session information, pre-specified expectation files
            and checksums.

## Data (not included; download from the original sources)
- GEO GSE93511: co-culture RNA-seq (35%, 50% and 60% PEG, days 1, 2, 3 and 5; n = 11).
- GEO GSE57338: human heart failure myocardium.
- Hahn et al. 2021 (HFpEF) data, Zenodo: [add the DOI; see docs/data_source.md]

Place the downloaded data in a data/ folder in the project root before running the scripts.

## How to run
Run each script from the project root, in numerical order. The random seed is 42
throughout. Scripts 04 and 04a are R scripts; the others are Python.

## Environment
- Python 3.13.17 (scikit-learn 1.9.1, NumPy 2.5.4, pandas 3.1.1, SciPy 1.18.1,
  matplotlib 3.11.2, PyDESeq2 0.6.0)
- R 4.6.1 (DESeq2 1.52.0, msigdbr 26.1.1)
- Full package list: docs/pip_freeze.txt and docs/sessionInfo*.txt

## Notes
- The human cross-check is exploratory and hypothesis-generating. The expected direction
  was recorded before the expression data were opened, and was not met.
- The GSE57338 expectation file was locked after the Hahn results had been viewed.

## Licence
MIT
