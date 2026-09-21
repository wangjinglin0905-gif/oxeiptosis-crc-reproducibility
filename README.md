# Oxeiptosis CRC analysis — public reproducibility package (v04)

This draft release contains the reproducibility inputs, expected outputs, analysis code, figures, and supplementary tables for the colorectal-cancer OX4 transcript-score analysis.

## Scope

The package reproduces the declared association analysis in three public cohorts (TCGA-COADREAD, GSE39582, and GSE17538). It is a computational reproducibility package, not a claim that the OX4 score measures functional oxeiptosis activity. The analysis uses public secondary data and does not contain newly collected human specimens.

The primary OX4 score is the equal-weight transcript score of **KEAP1, PGAM5, AIFM1, and OTUD1**, with the standardization rules documented in `supplementary_methods_v04.md` and `reproducibility/README.md`. The package also includes the 6,000 matched random-background models and the 300 recurrence diagnostics; the latter are diagnostic outputs and are not treated as valid inferential evidence when convergence warnings are present.

## Reproduce

From the `reproducibility` directory:

```text
python verify_inputs.py
Rscript reproduce.R smoke
Rscript reproduce.R full
```

The package was validated with R 4.6.1 and survival 3.8-6. `smoke` fits the 45 primary OS models. `full` additionally runs the recurrence diagnostics and 6,000 fixed random backgrounds. Outputs are written to a user-selected `run_*` directory and do not overwrite the frozen inputs or expected results.

## Data provenance

The expression and clinical data are reused from public accessions and are represented here by the derived model-input tables needed to reproduce the reported analyses. The original public datasets remain governed by their respective repositories and accession terms. Source accession and reconstruction notes are given in the supplementary methods and reference files.

## Contents

- `reproducibility/inputs/`: frozen model inputs, endpoint frames, random registries, and random score matrices.
- `reproducibility/expected/`: frozen expected numerical outputs and session information.
- `reproducibility/reproduce.R`: portable reproduction script.
- `reproducibility/verify_inputs.py`: SHA256 input verification.
- `figures/`: editable SVG, vector PDF, 300 dpi PNG, and 600 dpi RGB LZW TIFF exports.
- `supplementary_tables/`: the 27 supplementary CSV tables.
- `supplementary_methods_v04.md` and `.docx`: methods and interpretation boundaries.
- `references_v04.md`: references used by the revised manuscript.

## Interpretation boundary

The results support reproducibility of this particular score and its association analysis. They do not establish a universal oxeiptosis activity assay, causal mechanism, treatment response, or absence of biological relevance. Observational associations should be interpreted with the stated cohort, endpoint, covariate, and zero-time assumptions.

## Release status

Version 0.4.0 is archived at Zenodo: https://doi.org/10.5281/zenodo.22865462. The companion public repository is https://github.com/wangjinglin0905-gif/oxeiptosis-crc-reproducibility. The Zenodo record is the citable archive; the GitHub repository is the working code-and-file mirror.

## Authors

- Xia Zhang — Department of Oncology, Guihang Guiyang Hospital (ORCID: https://orcid.org/0009-0003-1978-4976)
- Jinglin Wang — Department of Digestive Disease, Guizhou Provincial People's Hospital (ORCID: https://orcid.org/0000-0002-7389-1022)

Code is licensed under MIT; original documentation, derived tables, model-input tables, and figures are licensed under CC BY 4.0. See LICENSE_SCOPE.md.
