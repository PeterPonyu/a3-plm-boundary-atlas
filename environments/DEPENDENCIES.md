# Dependency layers and recorded environments

No dependency was installed for this export.

## Default evidence verifier

Python >=3.10, standard library only. pytest is needed for the optional test suite.
No GPU, torch, fair-esm, NumPy, R, LaTeX or network is required by the default
checksum/row-score verification command.

## Completed inference provenance (not the export runtime)

`row-scores/provenance.json` records torch **2.12.1**, fair-esm **2.0.0** and
scikit-learn **1.6.1**, CUDA, fp32, TF32 false and an NVIDIA GeForce RTX 5090 Laptop
GPU. The original handoff names environment `cofold5090`. Python/NumPy/SciPy
versions for that historical inference environment are not recorded there and
are not guessed. Weights are not bundled.

## Export verification runtime observed locally on 2026-10-02

Python 3.13.5; pytest 8.4.2; NumPy 2.2.6; SciPy 1.16.3; scikit-learn 1.8.0;
pandas 2.3.3; PyMuPDF 1.27.1; Pillow 12.1.1; matplotlib 3.10.8; RDKit 2026.3.3.
These are recorded observations, not a claim that new scientific results were
fitted in this environment or that this is a universal lockfile.

R 4.3.3; ggplot2 4.0.3; patchwork 1.3.2; dplyr 1.2.1; tidyr 1.3.2; digest 0.6.39;
systemfonts 1.3.2; ragg 1.5.2; scales 1.4.0; base grid. Arial regular and bold are
required by the current rendering source and are not included as font files.

## Optional analysis/build modules

- Current clinical analysis: NumPy, SciPy, scikit-learn plus exported Python
  fold/calibrator modules. Preserved code is not run by the default verifier.
- Current figures: the R packages above, Python/PyMuPDF normalizer, matplotlib
  native-vector GA, pandas/RDKit/Pillow binding illustrations.
- Documents: latexmk, pdfLaTeX/TeX packages declared in Article/SI source,
  Poppler `pdfinfo`/`pdffonts`/`pdftoppm`. Optional coordinate rerender: PyMOL.
- Full scoring: torch, fair-esm, scikit-learn, official weights, matching upstream
  clinical/reference inputs and recorded N6c source dependencies.

Use an existing environment satisfying the selected tier; no install or download
is automatically attempted. An environment observation is not a dependency
resolution guarantee. The minimal default remains independently portable.
