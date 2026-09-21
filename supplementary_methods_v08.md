# Supplementary methods and results

This supplement describes the analysis populations, score definitions, model verification and exploratory diagnostics. S1–S27 Tables contain the corresponding machine-readable results. The primary estimands are the three cohort-specific overall-survival associations for OX4 under model M1.

## Analysis populations and endpoints

S1 Table describes 367 TCGA-COADREAD, 561 GSE39582 and 232 GSE17538 patients in the overall-survival analyses. S2 Table records the full analysis population, positive-time restriction and complete-stage restriction. TCGA patients were linked by 12-character identifiers, with one primary-tumour aliquot per patient and preference for sample vial 01A. Its processed HiSeqV2 expression and clinical matrices were obtained through UCSC Xena. Overall-survival time was recorded in days. GSE39582 used os.event and os.delay; GSE17538 used the recorded overall-survival status and duration in months. Cohorts were analysed separately, without pooling expression values or survival time units.

GSE39582 contained 566 eligible tumour expression profiles and 561 records with the required survival and clinical fields. Six records had zero recorded follow-up, including four deaths. They were retained in the primary model and removed in a positive-time sensitivity analysis. Their temporal interpretation could not be resolved from the available metadata. GSE17538 contributed 232 eligible human tumour profiles. Its GSE17536 component was not counted as another independent cohort. Missing stage was represented by an Unknown category; other clinical model covariates were required to be observed.

## Score definitions and standardisation

OX4 consists of KEAP1, PGAM5, AIFM1 and OTUD1. OX3 contains KEAP1, PGAM5 and AIFM1. OX5 contains KEAP1, PGAM5, AIFM1, NFE2L2 and AIRE; OX5 without NFE2L2 contains the other four members. OX6 contains the OX4 members together with CUL3 and EZH2. S3 Table records these definitions and diagnostic constructions. Each member received the same positive weight.

Within each expression reference population, each gene was centred and divided by its population SD (ddof = 0). The gene-wise standardised values were averaged, and the composite was centred and divided by its sample SD (ddof = 1). Reference populations contained 367 TCGA-COADREAD, 566 GSE39582 and 232 GSE17538 profiles. Consequently, GSE39582 survival-subset composite SDs differ slightly from one. Recalculating gene standardisation in a smaller population can change relative gene contributions; rescaling an already constructed composite is a different operation.

Array probes were mapped using the archived platform annotation. Multiple probe values for the same gene were averaged. All four OX4 members were required to be measurable; their seven GPL570 probes mapped unambiguously in the retained annotation (S6 Table). For larger reference sets, a probe assigned to multiple valid gene symbols could contribute to each mapped gene. This convention can add dependence between the constructed reference scores.

## Transcript references and model definitions

S4 Table lists each reference member by cohort with measurement and exclusion flags. S5 Table records reference-population and model-population sizes, fitted score SDs and reconstruction differences. The reference resources include FerrDb, KEGG, Gene Ontology and MitoCarta-derived lists. The mitochondrial inventory contains 1,136 MitoCarta3.0 human genes before measurement filtering. The 14-member SINGH_NFE2L2_TARGETS list belongs to MSigDB C2 chemical and genetic perturbations and does not include NFE2L2 itself.

The fixed exclusion set comprises KEAP1, PGAM5, AIFM1, OTUD1, CUL3 and EZH2. It was removed from the competing transcript references for primary exploratory M2 models and reference-pair comparisons. A separately labelled OX5 sensitivity also removed NFE2L2 where relevant. Retained resource-member lists, rather than a subsequent database update, define the computations.

Cox models used Efron handling of ties. M0 included only the score; M1 added continuous age, sex and categorical stage; M2 further added eight transcript references. These references were ferroptosis, apoptosis, necroptosis, parthanatos, mitophagy, oxidative-stress response, proliferation and mitochondrial transcripts. Stage I was the reference where present, and unused levels were dropped. HRs are per reference-population composite SD. Confidence intervals are 95% model-based intervals; two-sided P values are nominal, without multiplicity adjustment. All exploratory specifications are reported without selection by significance.

## Numerical verification

S7 Table contains all 45 candidate-definition by cohort by model overall-survival fits from the separate R implementation. S8 Table contains global rank-transformed cox.zph diagnostics for the nine OX4 fits. R 4.6.1 and survival 3.8-6 were used for verification. Across the 45 fits, the maximum absolute differences between R and the Python lifelines output were 0.0000362 for HR and 0.0002811 for P. Principal survival tables use the R results. This comparison evaluates numerical implementation; it does not add independent patients or independent human assessment.

## Transcript context and reference-set pairs

S9 Table reports raw, subcohort-adjusted and conditional rank correlations. Ranked variables were linearly regressed on the specified adjustment variables, and Pearson correlation was calculated between their residuals. Residuals were not reranked. Subcohort adjustment used Moffitt/Vanderbilt indicators in GSE17538 and discovery/validation indicators in GSE39582.

S10 Table decomposes equal-weight score variance into diagonal terms and signed covariance terms. Equal diagonal contributions follow from equal weights and constituent variances in the standardisation population. Their sum can exceed one when aggregate covariance is negative. These algebraic quantities do not rank genes by biological importance.

S11 Table contains 108 cohort-by-pair rows: 36 unordered pairs among nine references in each cohort. S12 Table summarises the pairs within cohorts. These comparisons used 367, 566 and 232 expression profiles. Four, 15 and 30 pairs, respectively, had Jaccard overlap below 0.05 and Spearman correlation at least 0.50. Six pairs in each cohort had zero overlap. Pairs share reference sets and patient data, so summaries are descriptive; no pooled independent-row test or row-bootstrap CI was used. S13 Table reports a two-member KEAP1–AIFM1 coverage-loss simulation. Its squared rank correlation is a descriptive measure of rank agreement, not variance explained in the original expression scale.

## Clinical and recurrence sensitivities

S14 Table contains extended clinical models. GSE39582 covariates included adjuvant chemotherapy, mismatch-repair status and proximal/distal location. TCGA used postoperative treatment, mismatch-repair protein expression by immunohistochemistry and colon/rectum location. These covariates are not identical across cohorts. In the GSE39582 joint complete-case subset of 499 patients, M1 gave HR 0.989 and the extended model gave HR 0.987, compared with 0.939 in the full M1 population. S15 Table provides overall-survival time and stage sensitivities.

S16 Table contains 300 recurrence risk-set specifications. Thirty fits carry convergence warnings and are unsuitable for inferential interpretation. TCGA specifications use new tumour events and are excluded from recurrence interpretation. GSE39582 specifications progressively restrict recorded time and stage. V3 retains stage II–III zero-time records; V4 requires positive time; V5 also censors follow-up at 60 months. S17 Table contains the 30 V5 rows for GSE39582 and GSE17538.

The GSE39582 cohort report supports a stage II–III and five-year analysis window. Zero-time exclusion remains an assumption. Applying the same specification to GSE17538 creates a harmonised sensitivity because patient-level treatment and disease-free-state information were incomplete. The OX4 M1 estimates were:

| Cohort | Patients | Events | HR (95% CI) | Nominal P |
| --- | --- | --- | --- | --- |
| GSE39582 | 453 | 130 | 0.976 (0.821–1.160) | 0.783 |
| GSE17538 | 145 | 31 | 1.393 (0.914–2.123) | 0.123 |

Matched-background analyses were not performed for these recurrence risk sets.

## KEAP1 and MSI sensitivities

S18 Table records alternative KEAP1 classification routes. Any mutation, nonsynonymous mutation or deep deletion, and the GDC nonsynonymous route identified 9, 6 and 4 patients, respectively. These calls were not experimentally confirmed protein inactivation. In the six-patient group, mean standardised OX4 was 1.151 versus −0.019 in the comparison group. The comparison-minus-altered difference was −1.170; Cohen d was −1.181 and the nominal Mann–Whitney P was approximately 0.0036. The minimum-group-size feasibility criterion of ten was not met. This criterion was not a formal power calculation.

S19 Table retains alternative MSI model specifications. The analysed row is labelled unknown_and_clean: OX4 HR 0.936937, 95% CI 0.729792–1.202878 and P 0.609368, with 367 patients, 85 deaths and ten parameters. The other rows document alternative missing-stage and reference-overlap handling. The GDC MSIsensor2 classification used a 20% threshold and case-level aggregation; 49 patients were labelled MSI-H. This genomic label is distinct from clinical immunohistochemistry or PCR classifications.

## Published lists and targeted source audit

S20 Table distinguishes five shared seed-member lists from their downstream analysis objects; S21 Table provides supplementary-file locations. The endometrial and gastric studies derived lncRNA signatures. Ma and colleagues constructed a three-gene model from a broader programmed-cell-death pool, and Du and colleagues developed a broader machine-learning index. NeoSAC reported an oxeiptosis process score, but available method details did not establish algorithmic identity with the equal-weight OX5 score. Two reports attributed their seed list to Zou and colleagues; the same provenance was not established for all five.

S22 Table documents 20 targeted records, source access, publication type, inspected molecular readouts and source locations. Ten records had full-text access, eight had abstracts only and two lacked usable material. Five were reviews and one was commentary. No field-wide evidence proportion was calculated. Protein modification, protein abundance, transcript abundance, localisation and viability were assessed as distinct readouts. Abstract-only access did not justify coding a measurement as absent. Kang and colleagues measured protein abundance and phosphorylation in Fig 3C; Pallichankandy and colleagues combined protein measurements and functional perturbations in Figs 5–7. In the methyl-cantharidimide report, displayed abundance measurements and a discussion statement about PGAM5 reduction differed; this unresolved discrepancy is retained in the source record.

## Matched backgrounds and platform sensitivity

S23 Table contains six background summaries; S24 Table contains 6,000 R refits. Each design sampled 1,000 four-gene sets per cohort, with seed 20260911. MSD matched expressed status, mean-expression decile and variability decile. The second design additionally constrained mean absolute candidate-to-OX4 Spearman correlation to within 0.1 of the target value. The archive labels these designs M1_msd and M2_msd_coexpr; these names denote sampling designs, not the survival-model adjustment sets. All background outcome models used survival M1.

Original OX4 genes were excluded from the candidate pool and sampling was without replacement within each set. The second design did not constrain internal pairwise correlations among the four newly drawn genes. Percentiles report the proportion of random absolute log HRs below the observed value. Upper-tail fractions use (1 + number at least as large as observed)/(1 + number of draws). These are comparisons with the specified empirical backgrounds, not equivalence tests. All 6,000 R fits completed without Cox convergence warnings; the maximum absolute log-HR discrepancy was approximately 1.97 × 10⁻⁹.

S25 and S26 Tables describe GSE71187 platform sensitivity. The source contains 189 records, 99 tumour profiles and 52 eligible survival records with 22 deaths. Scores were standardised in the 99-profile expression population. AIRE was unmeasured, preventing complete OX5 and OX5-without-NFE2L2 construction. OX4 M0 gave HR 1.109 (95% CI 0.717–1.716). Its M1 result was exploratory given the limited event count. The two-member simulation in S13 Table is a separate diagnostic.

## Reference verification and reproduction resources

S27 Table covers the 17 main references, their DOI/PMID identities, publication-notice metadata and claim-support scope. Bibliographic matches were distinguished from support for a specific statement. Linked commentaries were not classified as corrections or retractions. Metadata checks cannot certify the absence of every unindexed publication notice.

Zenodo v0.8.0 is assigned DOI 10.5281/zenodo.22882797 for the release containing derived model inputs, random-score matrices, expected outputs, R code, tables, publication figures and manuscript materials. The GitHub companion provides the directory-preserving working copy. The package starts from derived analysis inputs; downloading and rebuilding every original expression matrix is outside its executable scope. Figure files supplied with this manuscript were regenerated from the same frozen result tables for PLOS ONE formatting. Supplementary table values remain those of the frozen analysis, with the reference-verification record updated for this submission review.

S2 Text supplies the SHA256 manifest for the 15 released input and expected-output files; the same manifest is included as reproducibility/manifest.json in v0.8.0. After extracting the archive, run python verify_inputs.py from its reproducibility directory; a successful check reports 15 checked files and no mismatches. The manifest applies to the exact released bytes, including line endings.

In the reproduction directory, Rscript reproduce.R smoke fits the 45 overall-survival models. Rscript reproduce.R full additionally fits 300 recurrence diagnostics and 6,000 background models. A second argument selects a new output directory. The full mode was executed from the packaged inputs during submission preparation; both model-result CSV files matched the archived expected files byte for byte. The executable package comprises inputs, expected, reproduce.R and verify_inputs.py. It does not include a complete workflow for downloading or rebuilding the source expression matrices.
