# Supplementary methods and results for the Oxeiptosis transcript score analysis

This supplement accompanies manuscript_revised_v04. Tables are supplied as UTF-8 CSV files in supplementary_tables. The CSV files, their definitions below and the R reproduction package are the operative supplement. Historical outcomes not eligible for inference are retained with explicit labels rather than silently dropped.

## Population and endpoint definitions

S1a reports characteristics of the 367 TCGA, 561 GSE39582 and 232 GSE17538 overall-survival populations. S1b records the historical denominator, positive-time sensitivity and complete-stage sensitivity. TCGA uses primary-tumour patient identifiers and OS time in days. GSE39582 uses os.event and os.delay; GSE17538 uses source OS status and duration in months. Time units are not pooled. Positive-time exclusion in GSE39582 removes six patients, including four deaths. Clinical case-count agreement does not establish the precise meaning of a zero interval.

Gene-wise standardization uses ddof 0; the final composite uses ddof 1. The TCGA reference population is already restricted to the 367 eligible patients, while the GSE39582 expression reference contains 566 profiles before survival-field restriction. The same reference is used for all candidate definitions in a cohort. Recalculating gene-wise standardization in another population can alter the relative contributions of the genes; it is not necessarily a single affine transformation of the composite. Rescaling an already formed composite alone is an affine transformation. Analyses here preserve the frozen convention rather than interchange these operations.

## Score and reference registry

S2a is the corrected registry. OX6 contains six genes, including CUL3 and EZH2; the historical registry that labelled it as four genes is superseded. S2b contains one record per declared gene, set and cohort, with flags for whether it is measured and whether it remains after the fixed OX6 exclusion. S2c records the actual reference and fitted sample sizes, fitted SD and reconstruction discrepancy. S2d provides the primary GPL570 probe mapping. The larger reference-set parser permits a probe mapped to more than one valid symbol to contribute to each mapped gene; that retained convention can introduce additional dependence and should be considered when interpreting reference correlations.

The fixed exclusion set is KEAP1, PGAM5, AIFM1, OTUD1, CUL3 and EZH2. It is not definition-specific. SINGH_NFE2L2_TARGETS is a C2 CGP list with 14 human genes, not a C3 list or a direct measurement of NRF2 activity. Resource members and their snapshots are authoritative for this computation; no claim is made that they exhaust the current ontology or that no other oxeiptosis list exists. The GO response-to-oxidative-stress and KEGG-derived files use the retained mapped members rather than a later silently refreshed database release.

## Primary survival results and model verification

S3a contains all 45 candidate-definition by cohort by model OS estimates from the current independent R computation. S3b contains OX4 joint PH diagnostics. The frozen M2 model uses the fixed six-gene reference exclusion. The separately named OX5 member-clean sensitivity in the historical analysis is not substituted into this 45-row table. Nominal exploratory findings are retained but no best-performing model is selected from them. The three OX4 M1 rows are the principal clinical estimates; repeated correlated models are not independent validations.

The nine OX4 fits reproduce the historical estimates to reporting precision. Across 45 observed fits, the maximum absolute differences between R and the original lifelines table were 0.0000362 in HR and 0.0002811 in P. The current revision uses R values throughout its principal numerical tables. This agreement concerns computational implementation, not independent author review, independent patients or validation of the biological interpretation.

## Transcriptional context and finite reference pairs

S4a contains raw, subcohort-adjusted and conditional rank correlations. Adjustment for a subcohort uses indicator variables, while conditioning on a transcript reference uses its ranks. Pearson correlation is then computed between ranked-variable residuals; these residuals are not reranked. S4b reports diagonal variance fractions and the sum of signed covariance terms. The equality of diagonal fractions follows equal weighting and equal constituent variances in the standardization population. Covariance can make their sum exceed one; a negative aggregate cross-term is not a negative biological importance value.

S4c records all 108 pair rows and the agreement with the upstream pair file. S4d gives descriptive statistics separately for each cohort. Correlations use the expression reference populations, including all 566 GSE39582 profiles. Rows share sets and patient-level data, so no pooled independent-row P value or row-resampling bootstrap CI is presented. The 49 low-overlap/high-correlation rows are 4, 15 and 30 across the three cohorts and are not 49 zero-overlap pairs. There are 18 zero-overlap rows in total. S4e is the in-silico two-member coverage-loss diagnostic; its squared rank correlation must not be labelled variance explained in the original score.

## Clinical sensitivity models

S5a preserves extended clinical models with provenance from the WorkBuddy clinical-adjustment analysis. Its fields are not equivalent across cohorts. GSE39582 uses adjuvant chemotherapy Y/N, dMMR/pMMR and proximal/distal tumour location. TCGA uses any postoperative treatment, loss of mismatch-repair protein expression by IHC and Colon/Rectum. IHC is not the genomic MSI classification used in S7. Complete-case restriction differs by added covariate. In the GSE39582 499-patient joint complete-case subset, M1 gives HR 0.989 and the extended model 0.987; the full-sample M1 is 0.939. This comparison isolates case restriction from the small additional shift attributable to the recorded terms. It does not rule out unmeasured confounding. S5b preserves independently verified OS time and stage sensitivities from the preceding audit; these stable checks were reused, not claimed as newly rerun here.

## Recurrence risk-set diagnostics

S6a contains all 300 historical recurrence specifications independently refitted for diagnosis. Thirty fits retain convergence warnings and must not be used for inference. The TCGA construction uses new tumour events and is explicitly excluded from recurrence interpretation. In GSE39582, variants range from the old recorded set to positive-time restriction, removal of stage IV, removal of Unknown stage, stage II–III restriction and five-year censoring. V3 retains stage II–III zero times; V4 also requires positive time; V5 additionally censors after 60 months. S6b contains the 30 GSE39582/GSE17538 V5 rows, with all candidate definitions and M0–M2 models.

Marisa et al. provide a source-informed basis for stage II–III and a five-year window in GSE39582. Excluding recorded zero intervals remains an explicit assumption because stage II–III zero records persist. Applying this same specification in GSE17538 is a harmonized sensitivity, not proof of the original disease-free risk set. Stage IV at diagnosis cannot universally be equated with absence of a later disease-free state; patient-level treatment and resection information would be needed. A cohort's recurrence/no-recurrence counts do not settle that question.

The corrected primary-score rows below replace the OX5 values previously mislabelled as OX4. The unchanged fact that both P values exceed 0.05 does not make the mislabelling immaterial.

| Cohort | Patients | Events | OX4 HR (95% CI) | Nominal P |
| --- | --- | --- | --- | --- |
| GSE39582 | 453 | 130 | 0.976 (0.821–1.160) | 0.783 |
| GSE17538 | 145 | 31 | 1.393 (0.914–2.123) | 0.123 |

No recurrence matched-null analysis is restored. It was not rerun on these reconstructed risk sets. Its absence is a not-performed analysis, not a negative finding.

## KEAP1 and MSI sensitivity records

S7a records the available case-level KEAP1 alteration routes. Any mutation, nonsynonymous mutation or deep deletion, and the GDC nonsynonymous route identify 9, 6 and 4 patients in the frozen population. These groups are not equivalent to experimentally confirmed KEAP1 loss of function. The strict group has mean standardized OX4 1.151 versus -0.019 in patients without a detected specified alteration, difference -1.170 when expressed as comparison minus altered; Cohen d is -1.181 and the Mann–Whitney nominal P is approximately 0.0036. Linear rescaling changes the difference but not Cohen d or the rank test. The minimum-group-size gate was not met. It is a feasibility rule, not a formal power calculation, and does not erase the observed direction.

S7b contains three MSI-model versions. The operative row is unknown_and_clean: OX4 HR 0.936937, CI 0.729792–1.202878, P 0.609368, 367 patients and 85 deaths with 10 parameters. Historical submitted and unknown_fixed rows are retained for traceability, not mixed into the operative estimate. The GDC MSIsensor2 label used a 20% threshold and case-level aggregation; 49 patients were labelled MSI-H. It is not automatically interchangeable with a clinical IHC or PCR classification.

## Published member-list uses and targeted readout inventory

S8a separates the five shared seed-member lists from their downstream objects. S8b records supplementary-file cell anchors. NeoSAC explicitly calls its analysis an oxeiptosis score; its description does not establish an executable algorithm identity with the equally weighted OX5 score. Du's terminal index and Ma's broader PCD-derived three-gene model are distinct from that seed set. The two lncRNA papers derive other variables from the seed genes. Two papers cite Zou et al.; the same genealogy is not assumed for all five.

S9a replaces the historical ordinal and five-axis evidence tables. It lists 20 records, source access, publication class, inspected readouts, source anchors and interpretation limits. Its classes describe publications and do not claim deduplicated independent experiments. Keyword absence is not accepted as absence of a measurement. In particular, Kang Fig 3C and Pallichankandy Fig 5–7 contradict the previous absence-of-abundance labels. Modification, total protein, transcript, localization and viability measurements must be distinguished at assay level. The current source audit does not calculate a field-wide proportion of studies with paired abundance and modification data. Two misassigned abstracts remain excluded. Primary sources outside this 20-record set, including OTUD1 and cohort reports, are kept separately.

## Matched backgrounds and platform sensitivity

S10a reports the six matched-background summaries and S10b all 6,000 R refits. MSD denotes the mean/variability/status strata. The second design, historically labelled M2_msd_coexpr in files, is a random-set construction rule and must not be confused with the survival M2 model. All random outcomes reported here use survival M1. Candidate-to-OX4 correlation is not within-draw correlation. No convergence warning was produced in these R refits, with at most five iterations and maximum stored-versus-R absolute log-HR discrepancy 1.98e-9. This new complete check supersedes the earlier 120-fit verification sample for the matched M1 background only.

S11a/b retain the independently checked GSE71187 platform sensitivity from v03. The source has 189 records, 99 tumour expression profiles and 52 eligible survival records with 22 deaths. The score is standardized using the 99-profile expression reference. AIRE is unmeasured, so OX5 and its NFE2L2-excluded variant cannot be evaluated without changing the construct. The OX4 M0 HR is 1.109 (0.717–1.716); the M1 estimate is exploratory because of the small number of events. The two-member coverage simulation in S4e is a different analysis and is not used to claim cross-platform assay equivalence.

## References and reproducibility

S12a verifies every main-manuscript reference against current DOI/PMID metadata and records support scope. Source access and support were inspected separately from identifier matching. Comment-in records are distinguished from correction/retraction notices; absence of such a flag in one metadata response is not an absolute certification of publication status. The 22-entry baseline bibliography is retained in the audit's separate verification record, including the four source records without a fresh exact Europe PMC match in the initial refresh; Smith was subsequently resolved by PMID. Those not needed for the revised argument remain in the source inventory rather than being used to support new claims.

The code package reruns R models from supplied analysis frames and random score matrices. Raw-matrix reconstruction scripts and manifests document the local source dependencies. A complete fresh public-data download is not claimed to have been reproduced on a clean machine. User and software confirmation, source-data access arrangements and author declarations remain separate submission requirements.

## File index

| file | rows |
| --- | --- |
| S10a_matched_null_summary.csv | 6 |
| S10b_matched_null_all6000_R.csv | 6000 |
| S11a_GSE71187_R_recheck.csv | 6 |
| S11b_GSE71187_probe_coverage.csv | 8 |
| S12a_main_reference_verification.csv | 17 |
| S1a_baseline.csv | 3 |
| S1b_OS_population_variants.csv | 9 |
| S2a_score_registry_corrected.csv | 6 |
| S2b_reference_member_universe.csv | 10818 |
| S2c_reference_population_and_reconstruction.csv | 15 |
| S2d_primary_probe_mapping.csv | 7 |
| S3a_OS_all45_R.csv | 45 |
| S3b_OX4_PH_diagnostics.csv | 9 |
| S4a_transcript_and_subcohort_correlations.csv | 9 |
| S4b_variance_decomposition.csv | 15 |
| S4c_reference_pair_all108.csv | 108 |
| S4d_reference_pairs_by_cohort.csv | 3 |
| S4e_coverage_rank_diagnostic.csv | 3 |
| S5a_extended_clinical_adjustment_source.csv | 95 |
| S5b_OS_time_and_stage_sensitivities.csv | 27 |
| S6a_recurrence_diagnostic_all300.csv | 300 |
| S6b_fiveyear_sensitivity_all30.csv | 30 |
| S7a_KEAP1_case_classification.csv | 14 |
| S7b_MSI_model_versions.csv | 6 |
| S8a_published_seed_uses.csv | 5 |
| S8b_supplement_cell_trace.csv | 985 |
| S9a_targeted_source_inventory20.csv | 20 |
