# OX4 in colorectal cancer

Version 0.17.0 combines three bulk survival cohorts with spatial count-split analyses and a patient-level comparison of source-annotated tumour and stroma. Six patients contributed paired sections to the morphology analysis. All six patient averages were higher in tumour, with a mean Cliff's delta of 0.239. Median OX4 split correlations across sections ranged from 0.015 to 0.159.

The six authors of this version are Jinglin Wang, Pan Sun, Wei Zhu, Yan Zhang, Yao Tang and Juan Yang. The earlier bulk release retains its original creators and is available at https://doi.org/10.5281/zenodo.22882797.

- [Manuscript and figures](editions/v0.17.0/)
- [Supplement index](editions/v0.17.0/supporting_information/Supplement_index.md)
- [Spatial scripts, contracts and frozen-table verification](editions/v0.17.0/supporting_information/S1_Code_spatial.zip)
- [Bulk reproduction inputs, expected results and instructions](reproducibility/README.md)
- [Citation](CITATION.cff) and [license scope](LICENSE_SCOPE.md)

Run the bulk check from the reproducibility directory with `python verify_inputs.py`. The R entry point is `Rscript reproduce.R smoke` (or `full` with an output directory as documented there). For the spatial check, unzip Code S1 and run `python S1_Code_spatial_extension/verify_frozen_morphology.py` with NumPy and pandas installed. This verifies patient aggregation, bootstrap intervals, exact tests, matched-set percentiles and split summaries from the supplied tables. The historical spatial scripts provide source-level provenance and require local path configuration.

The version tag fixes the public snapshot: https://github.com/wangjinglin0905-gif/oxeiptosis-crc-reproducibility/tree/v0.17.0. Source expression matrices and tissue images remain in their original archives.
