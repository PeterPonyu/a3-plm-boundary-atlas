# Replacement identity-grouped nested refit

This directory contains the replacement uncertainty estimate for the 20-gene cross-fitted calibration audit. It supersedes the *absence* of a corrected interval in the active author-review record; it does not restore the withdrawn historical interval.

- Protocol: `A3-XFIT-NESTED-IDENTITY-GROUPED-CANDIDATE`
- Input: the locked 4,677-row manifest and unchanged five-fold OOF maps
- Bootstrap: 200 draws, seed `20260919`; copies of each original `example_id` stay in one fold
- Calibrators: z-score sigmoid, Platt, and isotonic are refit inside every draw
- Uncalibrated ECE: mean `0.248735`, 95% CI `[0.237547, 0.260993]`
- Pooled Platt ECE: mean `0.020427`, 95% CI `[0.013016, 0.029036]`
- Pooled isotonic ECE: mean `0.013056`, 95% CI `[0.006587, 0.020948]`

The estimate is conditional on this panel and its locked OOF maps. It is not an independent replication, an external validation result, or a recovery of the withdrawn historical `0.0207`/`91%` claim. The historical STOP file is unchanged; its SHA-256 is recorded in `round3_summary.json`.
