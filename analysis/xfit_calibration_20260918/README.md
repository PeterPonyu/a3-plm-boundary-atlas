# A3-XFIT-CAL-20260918

New 5-fold cross-fitted calibration. **Not** a restoration of the retired
historical full-20 Platt pair (`0.2443 → 0.0207` / 91%).

```bash
python3 run_xfit_calibration.py
python3 run_antibody_interface_robustness.py
python3 run_xfit_round2.py
```

Round 2 (`A3-XFIT-CAL-R2-20260918`) reuses the same seed/folds. It is not a
restoration of the retired historical pair.

Inputs are copied under `inputs/` with SHA-256 recorded in `outputs/xfit_calibration.json`.
The historical STOP file `manuscript/figures/plotdata/F7_20gene_perbin_STATUS.txt` is left unchanged.
