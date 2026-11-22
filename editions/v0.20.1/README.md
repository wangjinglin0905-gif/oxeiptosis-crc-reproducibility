[Download the complete 0.20.1 archive](https://github.com/wangjinglin0905-gif/oxeiptosis-crc-reproducibility/releases/download/v0.20.1/OX4_CRC_reproducibility_0.20.1.zip). The ZIP contains the source panel, all tables, scripts and figures referenced below. Run the instructions from its unpacked root.

# OX4 in colorectal cancer

Reproducibility materials for the tumour–stroma comparison and count-split analysis, version 0.20.1. Jinglin Wang, Pan Sun, Wei Zhu, Yan Zhang, Yao Tang and Juan Yang. Juan Yang is the corresponding author. The initial bulk archive was released by Xia Zhang and Jinglin Wang (https://doi.org/10.5281/zenodo.22882797). The earlier six-creator spatial edition is https://doi.org/10.5281/zenodo.24330126. Reuse these version-specific citations for their respective contributions.

This edition includes the annotated spatial map and corrected reference-overlap figure. See [figure inputs, attribution and rebuild instructions](FIGURES.md). The previous numerical baseline is [version 0.20.0](https://doi.org/10.5281/zenodo.24421516). Cite version 0.20.1 at https://doi.org/10.5281/zenodo.24432899.

## Contents and verification

The tables directory contains results S1 through S43. The two pathology cohorts are Valdeolivas (6 patients, 12 sections) and GSE294385 (5 patients, 8 sections). GSE280315 and GSE226997 are technical assessments. See methods/Archive_index.md for the baseline identifiers. Historical figure and document labels in the methods identify the source edition.

From the unpacked top-level directory, run these commands with Python 3.11+ and NumPy and pandas installed:

    python checks/verify_frozen_morphology.py
    python checks/verify_external.py
    python checks/verify_direction.py

For the bulk models, change into reproducibility/ and run python verify_inputs.py. Follow its README to refit the 45 or 6,345 models in R. The historical_spatial/ directory contains the original analysis routines for the earlier spatial modules as source documentation; those require their original reference directories, package layouts and inputs. The archived full earlier workflows and purchase-free source data locations are linked above and in the methods.

## Recalculate the external morphology endpoint

Obtain the eight original GEO archives and the pinned annotation RDS using the URLs and checksums in tables/S40_input_provenance.json. Place them in inputs/ at the top of this package. Create run/ first. From the top level run:

    Rscript external/export_annotations.R inputs/visium_meta_after_qc.rds run/source_primary_labels.csv
    python external/prepare_inputs.py
    python external/replicate_endpoint.py

Copy external/analysis_contract.md to run/ before running the last command. Python dependencies are NumPy, pandas, SciPy, openpyxl, requests and threadpoolctl. Preparation also fetches the published cohort table, GEO file catalogue and pinned source annotation for identity checks. Override OX4_INPUTS and OX4_OUTPUT for alternative directories and adjust the R paths accordingly. The distributed table check works entirely offline.

Public source matrices and tissue images remain with their originating repositories. Original analysis code is available under MIT and original derived results under CC BY 4.0. Source data and upstream contributed materials retain their original terms and attributions.
