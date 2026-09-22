# Oxeiptosis CRC model reproduction

This package reproduces Cox models from derived analysis frames and fixed random-score matrices. It does not download or reconstruct every original expression matrix.

## Requirements and commands

R 4.6.1 and survival 3.8-6 were used for verification. Python's standard library is sufficient for checksum verification. Run the following commands from the reproducibility directory after preserving its directory structure:

```text
python verify_inputs.py
Rscript reproduce.R smoke
Rscript reproduce.R full run_verification
```

The manifest covers 15 released input and expected-output files. It applies to the exact bytes of the Zenodo archive; Git checkouts that alter line endings may not match. Keep input and expected files unchanged. Use a distinct result-directory argument when preserving previous runs.

The smoke mode fits 45 overall-survival models. Full mode additionally fits 300 recurrence diagnostics and 6,000 matched random models. Full mode was run from the packaged inputs on 2026-09-22; both output model CSV files matched the expected CSV files byte for byte. Thirty recurrence models carry convergence warnings and are diagnostic only. The 6,000 background models carry no Cox convergence warning.

## Included resources

- inputs: derived model frames, fixed random memberships and score matrices.
- expected: reference outputs and R session information.
- reproduce.R: model-fitting entry point.
- verify_inputs.py and manifest.json: checksums for the released inputs and expected outputs.

The package does not include audit_source or provenance directories. See the manuscript and supporting methods for score construction, population definitions and interpretation limits. This numerical reproduction is not independent human source assessment or functional calibration of the score.

OX4 averages KEAP1, PGAM5, AIFM1 and OTUD1. Constituent standardisation uses ddof=0; composite standardisation uses ddof=1. GSE39582 uses 566 expression profiles as the reference population and 561 patients in the OS model; do not silently rescale using the smaller fitted subset.

Version-specific Zenodo DOI: https://doi.org/10.5281/zenodo.22882797 (published and verified)
Companion repository: https://github.com/wangjinglin0905-gif/oxeiptosis-crc-reproducibility
