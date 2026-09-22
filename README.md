# Oxeiptosis CRC analysis: reproducibility package and manuscript materials

Version 0.8.0 accompanies the revised PLOS ONE manuscript materials dated 2026-09-22. It preserves the frozen analytical results and updates scientific wording, figure formatting, supporting-file numbering and reproduction instructions. The manuscript is an author-review version, not a peer-reviewed publication.

## Current materials

- [Manuscript](manuscript/manuscript_PLOS_v08.docx) and [readable text](manuscript/manuscript_PLOS_v08.md).
- [Cover letter](manuscript/cover_letter_PLOS_ONE_v08.docx).
- [Figures](figures/): four 600 dpi RGB LZW TIFF files.
- [Figure sources](figure_sources/): SVG, PDF, PNG and the generation script.
- [Supporting information](supporting_information/): S1 Text, S2 Text and S1–S27 Table.
- [Supplementary methods](supplementary_methods_v08.md).
- [Old-to-new supplement mapping](SUPPLEMENT_FILE_MAP.csv).

The two creators are Xia Zhang (Department of Oncology, Guihang Guiyang Hospital; ORCID 0009-0003-1978-4976) and Jinglin Wang (Department of Digestive Disease, Guizhou Provincial People's Hospital; ORCID 0000-0002-7389-1022).

## Reproduce the numerical models

From the reproducibility directory:

```text
python verify_inputs.py
Rscript reproduce.R smoke
Rscript reproduce.R full run_verification
```

R 4.6.1 and survival 3.8-6 were used. Full mode was executed from the packaged derived inputs on 2026-09-22: 45 OS models, 300 recurrence diagnostics and 6,000 fixed-member matched-background models. Both model-result CSV files matched the archived expected files byte for byte. Thirty recurrence diagnostics retain convergence warnings; none of the 6,000 random-background fits had a Cox convergence warning. The separate comparison between Python and R background estimates has a maximum absolute log-HR difference of approximately 1.97e-9.

The SHA256 manifest verifies 15 input/expected files against the exact previously released bytes. The previously omitted manifest is now included. Text conversion is disabled in .gitattributes to preserve those bytes. This is a derived-input model-reproduction package; it does not download or rebuild every original expression matrix. See [the reproduction instructions](reproducibility/README.md).

## Interpretation and data provenance

OX4 is an equal-weight transcript score of KEAP1, PGAM5, AIFM1 and OTUD1. Constituent standardisation uses population SD and composite standardisation uses sample SD. The GSE39582 expression reference comprises 566 profiles; its OS analysis contains 561 patients. The other primary OS cohorts contain 367 TCGA-COADREAD and 232 GSE17538 patients. Original public matrices remain at their source repositories.

The study evaluates the construction, survival associations and transcript context of this particular score. It does not provide functional calibration of oxeiptosis activity. Literature checks are targeted and AI assisted, without completed independent human ratings. AI use, uncertainty and exploratory analyses are disclosed in the manuscript. Institutional and submission declarations remain subject to author confirmation.

## Versions and citation

The previous v0.4.0 archive remains available at https://doi.org/10.5281/zenodo.22865462. That DOI identifies the earlier release, not these updated v0.8.0 manuscript and figure files. The corresponding legacy document/table/figure copies are retained in archive/v0.4.0; shared model inputs stay in reproducibility.

Version 0.8.0 is published at https://doi.org/10.5281/zenodo.22882797 (record: https://zenodo.org/records/22882797). The public ZIP was downloaded and verified against the release package: 132 files, 130 manifest entries, and matching archive SHA256 c215a45802f492ea11bcd25fafad87633f75f07f8179e095a7479121a768487a. The immutable Zenodo ZIP matches Git commit d591af2de8ef5797e819095df5432897cea0d5b7. Later Git commits may update publication-status documentation and repository checksums; the manuscript, figures, supporting files and analytical inputs remain those of the archived snapshot. Historical draft-status wording inside the frozen ZIP records the state when it was prepared. See PUBLICATION_STATUS.md for the completed verification.

Code is MIT licensed; original documentation, manuscripts, derived tables and figures are CC BY 4.0. Third-party data and references retain their source terms. See LICENSE_SCOPE.md.
