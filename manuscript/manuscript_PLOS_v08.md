Construction and interpretation of an oxeiptosis-associated transcript score in colorectal cancer: a reproducible public-data analysis

Short title: Oxeiptosis-associated transcript score in colorectal cancer

Xia Zhang1, Jinglin Wang2,*

1 Department of Oncology, Guihang Guiyang Hospital, Guiyang, China

2 Department of Digestive Disease, Guizhou Provincial People’s Hospital, Guiyang, China

Corresponding author: Jinglin Wang; wangjinglin@gz5055.com; ORCID 0000-0002-7389-1022

Xia Zhang ORCID: 0009-0003-1978-4976

Article type: Research Article

Abstract

Background: Oxeiptosis-associated gene lists are increasingly used in cancer transcriptomics, but the survival associations and transcriptional context of the resulting scores require evaluation.

Methods: We evaluated an equally weighted KEAP1–PGAM5–AIFM1–OTUD1 score in TCGA-COADREAD, GSE39582 and GSE17538. Overall-survival analyses included 367, 561 and 232 patients, respectively, with adjustment for age, sex and stage. We examined four alternative member lists, two matched random-set designs, transcriptional reference scores and platform coverage. Scores were reconstructed from expression matrices, and 6,000 matched random-set models were refitted using a separate R implementation. A targeted source audit distinguished published gene lists from their downstream scoring procedures.

Results: Adjusted hazard ratios were 0.939 (95% confidence interval 0.753–1.171), 0.939 (0.809–1.091) and 1.175 (0.935–1.476). The observed absolute log hazard ratios lay at percentiles 40.1, 44.5 and 45.5 of the expression-and-variability-matched backgrounds. The score correlated with mitochondrial, proliferation and NRF2-related transcript references; these correlations attenuated after subcohort adjustment in GSE17538. Low member overlap between reference sets did not ensure low score correlation. Published five-gene lists supported heterogeneous downstream analyses. The largest absolute log-hazard-ratio difference across the 6,000 R refits was approximately 1.97 × 10⁻⁹.

Conclusions: This analysis provides a reproducible assessment of score construction, clinical association and transcriptional context across three colorectal cancer cohorts. The candidate score lacked a consistent adjusted survival association and shared variation with broader transcript programmes. These findings support explicit reporting of score definitions and comparator structure, with functional calibration required before interpreting the score as oxeiptosis activity.

Keywords oxeiptosis; colorectal cancer; gene expression; gene-set scoring; overall survival; reproducibility

Introduction

Oxeiptosis is a reactive-oxygen-species-responsive, caspase-independent cell-death pathway involving KEAP1, PGAM5 and AIFM1. Its original characterisation combined molecular perturbation, AIFM1 phosphorylation and cellular phenotypes [1]. Subsequent work examined OTUD1 in KEAP1-associated stress responses and cell death [2]. In colorectal cancer models, sanguinarine treatment implicated the KEAP1–PGAM5–AIFM1 axis through protein measurements, genetic perturbation and cell-death assays [3]. These studies motivate the investigation of pathway-associated genes in tumour datasets. Whether their steady-state transcript abundances capture the stress-dependent death process remains a separate measurement question.

Small gene lists enter public-data analyses in several ways. They can identify mechanistic participants, seed the selection of correlated transcripts, or define an enrichment score. Oxeiptosis-associated long noncoding RNA (lncRNA) studies in endometrial and gastric cancer used pathway genes to derive other prognostic variables [4, 5]. A five-gene list also appeared in a treatment-response analysis in the NeoSAC trial and in broader programmed-cell-death analyses [6–8]. Reuse of the same members therefore does not establish that these studies measured the same quantity.

Small transcript scores may share variation with broader expression programmes. Random-signature analyses in breast cancer showed that many unrelated gene sets could associate with outcome, partly through common transcriptional structure [9]. This motivates explicit background comparisons when evaluating pathway-associated scores. Such comparisons place an observed clinical association in context, while analyses of score composition address what the expression construct captures.

We evaluated an investigator-assembled, equally weighted four-gene score in three colorectal cancer cohorts. The primary question was whether its overall-survival association was consistent after adjustment for age, sex and stage. We then examined member selection, matched random-set backgrounds, platform coverage and correlated transcript references. A targeted source audit established how published oxeiptosis-related lists had been used. The study combines clinical benchmarking with construction diagnostics and a public model-reproduction package, including 108 cohort-specific reference-pair comparisons and 6,000 background models.

Methods

Study design and analysis populations

This retrospective secondary analysis used public expression matrices and accompanying clinical annotations from TCGA-COADREAD, GSE39582 and GSE17538. The GEO studies and their patient populations are described by Marisa and colleagues and Smith and colleagues [10, 11]. The processed TCGA HiSeqV2 matrix and its clinical matrix were obtained through the UCSC Xena TCGA hub. The present analysis used the archived clinical fields described in Supplementary Methods, not the harmonised TCGA Clinical Data Resource endpoints [12]. Original processed expression values were retained; no cross-cohort batch correction or pooled expression model was applied.

For TCGA, one primary-tumour aliquot per patient was selected, preferring sample vial 01A when available. Patients were linked by the 12-character TCGA identifier. The frozen analysis contained 367 patients with positive overall-survival time, an event indicator and required clinical data. GSE39582 contributed 566 tumour expression profiles and 561 patients with analysable overall survival. GSE17538 contributed 232 eligible human tumour profiles after excluding non-target samples from the SuperSeries. GSE17536 is a component of GSE17538 and was not treated as an additional validation cohort. Patient identifiers, inclusion flags and cohort-specific denominators are retained in the analysis package.

The primary overall-survival analysis retained the full GSE39582 denominator, including six records with zero recorded follow-up, four of which were deaths. Excluding these records was examined as a sensitivity analysis. Their exact temporal interpretation cannot be established from the available public metadata. Missing stage was represented by a separate Unknown category. This coding preserves those patients but does not eliminate bias from missingness. Other covariates were required to be observed. Table 1 presents baseline characteristics; S2 Table records population and endpoint rules.

Table 1. Baseline characteristics of the overall-survival analysis populations.

| Cohort | n | Deaths | Age, years / Median (IQR) | Male, n | Stage I / II / III / IV / Unknown |
| --- | --- | --- | --- | --- | --- |
| TCGA-COADREAD | 367 | 85 | 66 (55.5–75) | 202 | 56 / 132 / 111 / 51 / 17 |
| GSE39582 | 561 | 191 | 68 (59–76) | 308 | 31 / 262 / 204 / 60 / 4 |
| GSE17538 | 232 | 93 | 65.5 (56–74) | 122 | 28 / 72 / 76 / 56 / 0 |

n, patients; IQR, interquartile range. Stage counts are shown in the order indicated; Unknown denotes missing stage.

Candidate definitions and score construction

The primary score, OX4, contained KEAP1, PGAM5, AIFM1 and OTUD1. OX3 contained the first three genes. OX5 contained KEAP1, PGAM5, AIFM1, NFE2L2 and AIRE; OX5 without NFE2L2 omitted that member. OX6 expanded OX4 with CUL3 and EZH2. We applied the same equal-weight procedure to all five definitions to isolate sensitivity to member selection. The positive weights define a transcript average; they do not encode the direction or strength of each member’s mechanistic contribution.

Array probes were mapped to gene symbols with the archived platform annotation, and available probe values for each gene were averaged. The seven GPL570 probes mapping to the four primary members were unambiguous in that annotation. The handling of multi-symbol probes in larger reference sets follows the recorded mapping code and is a limitation of those reference constructions. All four primary members were required to be measurable before a score was computed. Probe mappings and exact lists are provided in S3 and S6 Tables.

Within each reference population, each gene was centred and divided by its population standard deviation (SD). The standardised values were averaged, and the composite was centred and divided by its sample SD. Reference populations comprised 367 TCGA patients, all 566 eligible GSE39582 tumour profiles, and 232 GSE17538 tumour profiles. Hazard ratios (HRs) are reported per one reference-population composite SD. In the GSE39582 survival subset, score SDs were approximately 0.996–1.004. These unsupervised transformations were independent of outcome and were performed separately within each cohort.

Survival models and sensitivity analyses

Cox proportional-hazards models used Efron handling of tied event times. M0 included the score alone. The primary model, M1, included the score, continuous age, sex and categorical stage. The R formula was Surv(OS_time, OS_event) ~ OX4 + age + male + stage4. Stage I was the reference where present, and unused levels were dropped. M2 additionally included eight transcript references: ferroptosis, apoptosis, necroptosis, parthanatos, mitophagy, oxidative-stress response, proliferation and mitochondrial transcripts. M2 was an exploratory model for conditional association.

The reference lists came from the archived FerrDb, KEGG, Gene Ontology and MitoCarta-derived files. The nine-set correlation panel additionally included the 14-member SINGH_NFE2L2_TARGETS list. This list belongs to the MSigDB C2 chemical and genetic perturbations collection and derives from an NRF2-knockdown study in lung cancer [13]. It does not contain NFE2L2 itself. The mitochondrial list corresponds to the 1,136-member MitoCarta3.0 human inventory before measurement filtering [14]. It is used here as a transcript reference, not as a measure of mitochondrial number, mass or function. Resource identifiers, measured members and frozen file hashes are supplied in S4 and S5 Tables.

To reduce direct member overlap, the primary exploratory M2 and reference-pair analyses removed a fixed six-gene list, KEAP1, PGAM5, AIFM1, OTUD1, CUL3 and EZH2, from the reference sets. For OX5, a separately labelled sensitivity also excluded NFE2L2 from the relevant competing references. The SINGH comparison did not require such removal. Models with correlated reference scores and few events per fitted parameter were interpreted cautiously. Parameter counts and model warnings accompany the estimates in S7 and S8 Tables.

Proportional hazards were assessed for OX4 using the global rank-transformed cox.zph test. Confidence intervals (CIs) were 95% model-based intervals. P values were two-sided, nominal and unadjusted for multiple comparisons. The primary clinical estimands were the three cohort-specific OX4 M1 associations. Alternative definitions, reference adjustment, clinical and missing-stage sensitivities, altered-gene groups and recurrence analyses were exploratory and reported without significance-based selection. The study evaluated continuous-score associations and did not develop a prognostic cutpoint or trained prediction model.

Matched backgrounds and transcriptional diagnostics

Two matched random-set designs each contained 1,000 four-gene sets per cohort, with seed 20260911. The first matched each OX4 member to a candidate in the same expressed-status, mean-expression-decile and variability-decile stratum. Expressed status was defined using the cohort-specific 25th percentile of gene means. Original OX4 members were excluded from the candidate pool, and genes were sampled without replacement within a set. The second design additionally constrained each candidate's mean absolute Spearman correlation with the original OX4 members to within 0.1 of its target. This constraint matches candidate-to-OX4 correlation and does not match the internal pairwise correlations of the newly drawn set. Gene length and other unrecorded attributes were not matched.

Each random score underwent the same standardisation and M1 adjustment as OX4. We report the fraction of random absolute log hazard ratios below the observed absolute log hazard ratio. The associated upper-tail fraction was calculated as one plus the number at least as large as the observed value, divided by one plus the number of draws. These quantities describe the specified empirical backgrounds; they are not equivalence tests or tests of cell-death specificity. Membership was kept fixed during verification. The matched analyses were added during exploratory development; no external preregistration was performed.

Spearman correlations quantified associations between OX4 and transcript references. Partial rank correlations were calculated as Pearson correlations between residuals of ranked variables after linear adjustment. Subcohort labels entered the adjustment as indicator variables. These analyses quantify conditional association within the recorded cohort structure.

For the nine reference sets, all 36 unordered set pairs were evaluated separately in each expression reference population. Jaccard overlap was calculated from measured members after the fixed six-gene exclusion; score correlation used the corresponding expression profiles. The 108 cohort-by-pair rows share genes, reference scores and patients within cohorts and repeat the same reference pairs across cohorts. We therefore present them descriptively as cohort-faceted summaries rather than as pooled independent observations. A two-member KEAP1–AIFM1 score was also compared with OX4 as an in-silico coverage-loss diagnostic. Squared Spearman correlation is reported as a descriptive overlap measure, not as explained variance of the original expression scores.

Source audit and computational verification

The literature component was a targeted audit of 20 records with oxeiptosis in the title, supplemented by primary mechanistic, cohort and gene-set sources. Full texts, figure legends and supplements were inspected where available. Abstract-only records were not assigned negative measurement findings. For five published uses of the common five-gene list, we recorded member identities, the downstream analysis object and the described scoring procedure. This targeted audit contextualised score definitions and was not designed as a systematic review.

Candidate and random scores were reconstructed from archived expression matrices using Python. A separate R implementation refitted 45 overall-survival models, 300 recurrence specifications and 6,000 matched random M1 models, using R 4.6.1 and survival 3.8-6. The original candidate-model estimates used Python lifelines, whereas the random-background pipeline used a custom Python Cox routine. R estimates are reported in the principal survival tables. Convergence warnings were retained, and recurrence fits were classified as risk-set diagnostics. The public package reproduces models from supplied analysis frames and random-score matrices.

Ethics statement

This secondary analysis used de-identified, publicly available data. It involved no new participant recruitment, participant contact or collection of human specimens. No ethics application was submitted for the present analysis.

Use of AI-assisted tools

WorkBuddy and OpenAI Codex assisted with code drafting and debugging, file organisation, source retrieval, computational checks and manuscript drafting and language editing. Their use affected analysis scripts, supporting files and manuscript wording. Numerical checks compared results with frozen tables and a separate R implementation. Source checks compared bibliographic records and available primary texts with the claims they supported. The authors are responsible for reviewing these checks, the scientific interpretation and the final submitted content.

Results

Overall-survival estimates varied in direction across cohorts

The primary analyses included 1,160 patients and 369 deaths, analysed separately by cohort. Under M1, the OX4 HRs were 0.939 (95% CI 0.753–1.171), 0.939 (0.809–1.091) and 1.175 (0.935–1.476) for TCGA-COADREAD, GSE39582 and GSE17538, respectively (Table 2 and Fig 1). All three intervals crossed one. The point estimates differed in direction, and their uncertainty remained compatible with potentially relevant protective or adverse associations.

Fig 1. Overall-survival associations of OX4. Points show cohort-specific Cox estimates, with bars indicating 95% model-based CIs. M0 includes the score only; M1 additionally adjusts for age, sex and stage; M2 adds eight transcript references after the fixed six-gene exclusion. HRs are per reference-population score SD. The dashed line marks HR 1. Counts indicate patients/deaths. The models use the same patients within each cohort. Nominal P values and parameter counts are in S7 Table; global proportional-hazards diagnostics are in S8 Table.

Table 2. OX4 overall-survival associations.

| Cohort | Model | HR (95% CI) | Nominal P |
| --- | --- | --- | --- |
| TCGA-COADREAD | M0 | 0.867 (0.705–1.065) | 0.174 |
| TCGA-COADREAD | M1 | 0.939 (0.753–1.171) | 0.574 |
| TCGA-COADREAD | M2 | 0.941 (0.704–1.259) | 0.684 |
| GSE39582 | M0 | 0.875 (0.759–1.009) | 0.066 |
| GSE39582 | M1 | 0.939 (0.809–1.091) | 0.414 |
| GSE39582 | M2 | 1.089 (0.892–1.329) | 0.404 |
| GSE17538 | M0 | 1.073 (0.871–1.322) | 0.508 |
| GSE17538 | M1 | 1.175 (0.935–1.476) | 0.166 |
| GSE17538 | M2 | 0.856 (0.505–1.452) | 0.564 |

HR, hazard ratio; CI, confidence interval. M0 is unadjusted; M1 adjusts for age, sex and stage; M2 additionally includes eight transcript references. Estimates are per reference-population score standard deviation.

The primary OX4 estimates changed under the reference-adjusted M2 models to 0.941, 1.089 and 0.856, respectively. M2 included 15 fitted coefficients in TCGA and GSE39582 and 14 in GSE17538, corresponding to 5.7, 12.7 and 6.6 deaths per coefficient. Its correlated predictors and limited events constrain interpretation, particularly in the first and third cohorts. Joint proportional-hazards tests for the nine OX4 cohort-by-model fits did not reject at the nominal 0.05 level; diagnostics and parameter counts are reported in S7 and S8 Tables.

Alternative member lists also lacked a consistent M1 association across cohorts. Selected unadjusted or M2 estimates had intervals excluding one, including OX5 and OX5 without NFE2L2 in GSE17538 M2. These exploratory findings are reported alongside all 45 models in S7 Table.

Observed effects were not unusually large in the matched backgrounds

The OX4 absolute log hazard ratios under M1 were 0.0633, 0.0625 and 0.1612 in TCGA, GSE39582 and GSE17538. They lay at percentiles 40.1, 44.5 and 45.5 of the expression-and-variability-matched distributions. Under the candidate-to-OX4 correlation constraint, the percentiles were 41.5, 35.3 and 59.8 (Table 3 and Fig 2). Thus, neither design placed the observed estimate near the upper tail in these cohorts.

Fig 2. Matched-background distributions for the OX4 M1 association. Empirical cumulative distributions use 1,000 four-gene sets per design and cohort. The dashed vertical marker denotes OX4. MSD matches expressed status, mean-expression decile and variability decile. MSD + OX4 correlation adds a candidate-to-OX4 correlation constraint. Panels a–c show TCGA-COADREAD, GSE39582 and GSE17538. Percentiles compare the observed effect with the specified backgrounds; the second design does not match within-draw co-expression.

Table 3. Observed adjusted effects relative to matched backgrounds.

| Cohort | Observed absolute log HR | MSD percentile | MSD + OX4 correlation percentile |
| --- | --- | --- | --- |
| TCGA-COADREAD | 0.0633 | 40.1 | 41.5 |
| GSE39582 | 0.0625 | 44.5 | 35.3 |
| GSE17538 | 0.1612 | 45.5 | 59.8 |

HR, hazard ratio. MSD matches expressed status, mean-expression decile and variability decile. The second design adds candidate-to-OX4 correlation matching. Each background contains 1,000 sets per cohort.

All 6,000 matched models were refitted without a Cox convergence warning in R. The largest absolute difference from the stored random-model log HR was approximately 1.97 × 10⁻⁹, indicating numerical agreement between implementations.

The score shared variation with several transcriptional references

OX4 correlated positively with the mitochondrial transcript reference, proliferation reference and SINGH-NRF2 target score in all three cohorts (Fig 3). The respective correlations were 0.386, 0.305 and 0.323 in TCGA; 0.566, 0.520 and 0.473 in GSE39582; and 0.797, 0.793 and 0.704 in GSE17538. The mitochondrial reference excluded the directly shared AIFM1 and PGAM5 members. These associations therefore persisted beyond their direct contribution to both scores.

Fig 3. Transcript-reference correlations and subcohort composition. Panel a shows Spearman correlations of OX4 with mitochondrial, proliferation and SINGH-NRF2 references in the overall-survival populations. Panel b compares raw and subcohort-adjusted rank correlations in GSE17538, using the recorded Moffitt/Vanderbilt label. The mitochondrial reference excludes the directly shared AIFM1 and PGAM5 members. All references are transcript-based scores.

Adjustment for Moffitt versus Vanderbilt labels in GSE17538 reduced mitochondrial, proliferation and NRF2 correlations to 0.587, 0.582 and 0.484, respectively. In GSE39582, adjustment for the discovery-versus-validation label had little effect. Conditional OX4–mitochondrial and OX4–NRF2 rank correlations remained positive after adjustment for the other reference. Subcohort comparisons used subsets of the existing cohorts.

For the equal-weight score, each standardised member had the same diagonal variance contribution in the standardisation population. Covariance terms altered the total score variance. In TCGA, diagonal shares were approximately 27.217% per member and signed cross-terms summed to −8.87%. These quantities describe the variance of a correlated average. Small differences in the GSE39582 fitted subset followed its restriction from the expression reference population (S10 Table).

Low member overlap did not ensure low reference-score correlation. Among 36 pairs per cohort, 4, 15 and 30 had Jaccard overlap below 0.05 and Spearman correlation at least 0.50 in TCGA, GSE39582 and GSE17538, respectively. Each cohort contained six zero-overlap pairs, with maximum correlations of 0.405, 0.698 and 0.878. Within-cohort descriptive correlations between Jaccard overlap and score correlation were 0.523, 0.501 and 0.645 (Fig 4). No pair had Jaccard overlap of 0.20 or greater.

Fig 4. Member overlap and reference-score correlation by cohort. Each point represents one of 36 unordered pairs among nine transcript references. Jaccard overlap uses measured members after the fixed six-gene exclusion. Correlations use 367, 566 and 232 expression profiles in panels a–c. Dashed lines mark Jaccard 0.05 and Spearman correlation 0.50. Pairs share reference sets and patients and are shown descriptively. No pair had Jaccard overlap of 0.20 or greater.

An in-silico KEAP1–AIFM1 score correlated with OX4 at 0.770, 0.795 and 0.770 across the three expression populations. Removing PGAM5 and OTUD1 therefore changed patient ranks (S13 Table).

Clinical and endpoint sensitivity analyses

Complete-stage restriction left the OX4 M1 point estimates close to their primary values. Excluding the six zero-time overall-survival records in GSE39582 gave HR 0.950 (95% CI 0.816–1.106). In that cohort, extended adjustment on 499 complete cases gave HR 0.987 (0.838–1.163), while M1 fitted to the same 499 patients already gave HR 0.989. The change from the full-cohort estimate was therefore largely associated with case restriction in this comparison. TCGA used different treatment and pathology variables, and the corresponding complete-case model contained only 33 deaths; those exploratory results are documented separately in S14 and S15 Tables.

In GSE39582, the stage II–III, positive-time, five-year-censored recurrence sensitivity gave OX4 HR 0.976 (95% CI 0.821–1.160; 453 patients, 130 events). Applying the same restriction to GSE17538 gave HR 1.393 (0.914–2.123; 145 patients, 31 events). The GSE17538 restriction approximates a harmonised risk set because treatment and disease-free-state information were incomplete. S16 and S17 Tables report the diagnostic specifications and sensitivities. The TCGA new-tumour-event field was excluded from recurrence interpretation.

Six TCGA patients met the specified nonsynonymous-KEAP1-mutation or deep-deletion definition and had higher mean OX4 values than the comparison group. This exploratory comparison did not meet the minimum-group-size feasibility criterion. Alteration calls alone did not establish functional protein inactivation. Adding genomic microsatellite instability (MSI) status to the age-, sex-, stage-, proliferation- and mitochondrial-reference-adjusted model gave OX4 HR 0.937 (95% CI 0.730–1.203; S18 and S19 Tables).

Source verification distinguished member lists from measured constructs

Five inspected publications used the same five seed symbols after NRF2 was harmonised to NFE2L2. Their downstream uses differed. The endometrial and gastric studies derived lncRNA prognostic signatures; the Ma study derived a three-gene model from a broader programmed-cell-death pool; and the Du study used a broader cell-death analysis and a machine-learning index. NeoSAC explicitly reported an oxeiptosis process score for a treatment-response comparison, but the available description did not provide enough detail to reproduce every scoring parameter. None of these source descriptions was substituted for the equal-weight OX5 implementation in this study. Two publications attributed the seed set to the same earlier source [15]; provenance for the other three was not established from the inspected material (S20 and S21 Tables).

Of the 20 targeted records, ten had accessible full texts, eight had abstracts only and two lacked usable abstracts or full texts. Five records were reviews and one was commentary. The two records without usable source material were excluded from measurement-content assessment.

The inspected mechanistic studies measured different levels of the pathway. The sanguinarine colorectal cancer and vitiligo studies included protein abundance and phosphorylation measurements alongside functional assays [3, 16]. In the methyl-cantharidimide liver-cancer study, phosphorylation changed while the displayed measurements showed little change in total AIFM1 and PGAM5 abundance [17]. Its discussion also described reduced PGAM5 protein, an unresolved discrepancy recorded in S22 Table. The targeted sources provide context for molecular readouts but no patient-level calibration of OX4 in the three analysed cohorts.

Discussion

This study integrates survival benchmarking, transcript-context analysis and source tracing for a defined four-gene score. Across three colorectal cancer cohorts, clinically adjusted OX4 estimates differed in direction and were not unusually large within the specified matched backgrounds. The results characterise the behaviour of this expression construct while identifying two practical requirements for its interpretation: explicit score construction and explicit comparator structure.

Holding the scoring algorithm constant clarified how member choice changes the measurement. OX4 averages standardised transcript abundances with positive weights. Protein interactions, phosphorylation and the timing of oxidative stress are outside that definition. Substituting members changes the construct even when each candidate is described as oxeiptosis-associated. Likewise, the published lncRNA models and enrichment procedures represent distinct uses of seed lists. Direct comparisons require both the member list and the scoring algorithm.

OX4 shared variation with mitochondrial, proliferation and NRF2-related transcript references. The attenuation after GSE17538 subcohort adjustment illustrates how cohort composition can contribute to a strong transcript-level association. Institutional labels combine technical and biological differences, so the adjustment cannot distinguish their origins. Removing shared genes addresses direct overlap but leaves broader expression covariance and cell-composition effects. Accordingly, partial correlations and M2 coefficients describe conditional associations rather than causal sources of OX4 variation.

Reference-set overlap and score correlation captured different properties of the comparator panel. Some zero-overlap pairs remained strongly correlated, even though overlap and correlation were positively associated within each cohort. Both quantities should therefore be reported with their measurement universe and cohort. Because pairs share sets and patients, the 108 comparisons provide descriptive context rather than 108 independent tests.

Matched backgrounds placed the adjusted OX4 effect near the middle of distributions from four-gene sets with similar expression and variability. Adding a candidate-to-OX4 correlation constraint produced the same broad pattern. This extends the comparison logic of random-signature studies [9] to an oxeiptosis-associated candidate in colorectal cancer. The inference remains conditional on the matching rules: it is neither an equivalence result nor a test of pathway function. A calibrated score or a stress-dependent association could behave differently.

The combined analyses identify reporting practices that readers can use to evaluate other pathway-associated scores. These include the exact members and algorithm, standardisation population, platform coverage, comparator overlap and correlation, and outcome and risk-set definitions. The archived model inputs and verification scripts make these comparisons reproducible. These are practical recommendations from a case study, rather than a consensus standard. The next measurement question is whether a specified transcript construct tracks a matched functional phenotype under defined conditions.

Limitations

The study used retrospective processed bulk-expression data with incomplete treatment, molecular and cell-composition information. The primary cohorts differ in clinical composition, assay and follow-up, and residual confounding is plausible. We did not estimate tumour purity or test a cell-type-specific model. Missingness and zero follow-up records require assumptions, particularly for recurrence endpoints. M2 and some sensitivity models have limited events relative to their parameter counts, and model-based intervals may understate uncertainty from analytic choices.

Candidate and reference lists were investigator specified, and some sensitivity analyses were added after initial results without external preregistration. Background comparisons cover finite gene universes and matching choices. Within-cohort standardisation does not provide a locked assay for individual patients, and survival association does not establish functional calibration. Confidence intervals remained compatible with potentially relevant effects. The targeted literature audit was access limited and AI assisted, without completed independent human ratings. The small fourth-platform analysis was exploratory, and the two-member simulation did not reproduce an empirical cross-platform comparison.

Conclusions

Benchmarking across three colorectal cancer cohorts showed that the equally weighted KEAP1–PGAM5–AIFM1–OTUD1 score had imprecise adjusted survival associations and shared variation with broader transcript programmes. Its observed effect sizes were not unusual in the two specified matched backgrounds. The study provides a reproducible example of assessing clinical association together with score construction and comparator structure. Functional calibration remains necessary before treating this transcript score as a measurement of oxeiptosis activity.

Data availability

The original datasets are available from GEO (GSE39582, GSE17538 and GSE71187) and the UCSC Xena TCGA hub (TCGA-COADREAD). Derived model inputs, fixed random-set memberships and scores, R code, results and the original figure exports are archived at https://doi.org/10.5281/zenodo.22865462. The companion repository is https://github.com/wangjinglin0905-gif/oxeiptosis-crc-reproducibility. These resources support reproduction of the reported models from derived inputs; complete original expression matrices remain at the source repositories. Supplementary methods describe the analysis populations, transformations and scope of the reproduction package. S2 Text supplies the checksum manifest used with the archived inputs.

Acknowledgments

The authors thank the investigators and participants of TCGA-COADREAD, GSE39582, GSE17538 and GSE71187, whose publicly deposited data made this analysis possible, and the maintainers of the UCSC Xena, MSigDB, FerrDb, KEGG, Gene Ontology and MitoCarta resources used in the reference constructions.

References

1. Holze C, Michaudel C, Mackowiak C, Haas DA, Benda C, Hubel P, et al. Oxeiptosis, a ROS-induced caspase-independent apoptosis-like cell-death pathway. Nat Immunol. 2018;19(2):130-140. doi: 10.1038/s41590-017-0013-y.

2. Oikawa D, Gi M, Kosako H, Shimizu K, Takahashi H, Shiota M, et al. OTUD1 deubiquitinase regulates NF-κB- and KEAP1-mediated inflammatory responses and reactive oxygen species-associated cell death pathways. Cell Death Dis. 2022;13(8):694. doi: 10.1038/s41419-022-05145-5.

3. Pallichankandy S, Thayyullathil F, Cheratta AR, Subburayan K, Alakkal A, Sultana M, et al. Targeting oxeiptosis-mediated tumor suppression: a novel approach to treat colorectal cancers by sanguinarine. Cell Death Discov. 2023;9(1):94. doi: 10.1038/s41420-023-01376-3.

4. Niu L, Wu Z. Identification and validation of oxeiptosis-associated lncRNAs and prognosis-related signature genes to predict the immune status in uterine corpus endometrial carcinoma. Aging (Albany NY). 2023;15(10):4236-4252. doi: 10.18632/aging.204726.

5. Wen L, Xu K, Huang M, Pan Q. Identification of oxeiptosis-associated lncRNAs and prognosis-related signature to predict the immune status in gastric cancer. Medicine (Baltimore). 2024;103(7):e37189. doi: 10.1097/MD.0000000000037189.

6. Shen G, Liu Z, Wang M, Zhao Y, Liu X, Hou Y, et al. Neoadjuvant apatinib addition to sintilimab and carboplatin-taxane based chemotherapy in patients with early triple-negative breast cancer: the phase 2 NeoSAC trial. Signal Transduct Target Ther. 2025;10(1):41. doi: 10.1038/s41392-025-02137-7.

7. Du JO, Huang QK, Sun XC, Hu SK, Lin TS. Machine learning-driven transcriptomic and single-cell profiling of programed cell death patterns in colon cancer. Sci Prog. 2026;109(1):368504261432995. doi: 10.1177/00368504261432995.

8. Ma JY, Wang YX, Zhao ZY, Xiong ZY, Zhang ZL, Cai J, et al. Identification of key programmed cell death genes for predicting prognosis and treatment sensitivity in colorectal cancer. Front Oncol. 2024;14:1483987. doi: 10.3389/fonc.2024.1483987.

9. Venet D, Dumont JE, Detours V. Most random gene expression signatures are significantly associated with breast cancer outcome. PLoS Comput Biol. 2011;7(10):e1002240. doi: 10.1371/journal.pcbi.1002240.

10. Marisa L, de Reyniès A, Duval A, Selves J, Gaub MP, Vescovo L, et al. Gene expression classification of colon cancer into molecular subtypes: characterization, validation, and prognostic value. PLoS Med. 2013;10(5):e1001453. doi: 10.1371/journal.pmed.1001453.

11. Smith JJ, Deane NG, Wu F, Merchant NB, Zhang B, Jiang A, et al. Experimentally derived metastasis gene expression profile predicts recurrence and death in patients with colon cancer. Gastroenterology. 2010;138(3):958-968. doi: 10.1053/j.gastro.2009.11.005.

12. Liu J, Lichtenberg T, Hoadley KA, Poisson LM, Lazar AJ, Cherniack AD, et al. An Integrated TCGA Pan-Cancer Clinical Data Resource to Drive High-Quality Survival Outcome Analytics. Cell. 2018;173(2):400-416.e11. doi: 10.1016/j.cell.2018.02.052.

13. Singh A, Boldin-Adamsky S, Thimmulappa RK, Rath SK, Ashush H, Coulter J, et al. RNAi-mediated silencing of nuclear factor erythroid-2-related factor 2 gene expression in non-small cell lung cancer inhibits tumor growth and increases efficacy of chemotherapy. Cancer Res. 2008;68(19):7975-7984. doi: 10.1158/0008-5472.CAN-08-1401.

14. Rath S, Sharma R, Gupta R, Ast T, Chan C, Durham TJ, et al. MitoCarta3.0: an updated mitochondrial proteome now with sub-organelle localization and pathway annotations. Nucleic Acids Res. 2021;49(D1):D1541-D1547. doi: 10.1093/nar/gkaa1011.

15. Zou Y, Xie J, Zheng S, Liu W, Tang Y, Tian W, et al. Leveraging diverse cell-death patterns to predict the prognosis and drug sensitivity of triple-negative breast cancer patients after surgery. Int J Surg. 2022;107:106936. doi: 10.1016/j.ijsu.2022.106936.

16. Kang P, Chen J, Zhang W, Guo N, Yi X, Cui T, et al. Oxeiptosis: a novel pathway of melanocytes death in response to oxidative stress in vitiligo. Cell Death Discov. 2022;8(1):70. doi: 10.1038/s41420-022-00863-3.

17. Li YD, Chen X, Dong XD, Ding Z, Chen ZS, Fang S. Methyl-cantharidimide suppresses cyclin-dependent kinase 1 (CDK1) and induces oxeiptosis in liver cancer. J Exp Clin Cancer Res. 2026;45(1):90. doi: 10.1186/s13046-026-03669-8.

Supporting information

S1 Text. Supplementary methods and results, including population definitions, score construction, model diagnostics, source assessment and reproduction instructions.

S2 Text. SHA256 manifest for the 15 archived reproduction inputs and expected outputs. Save as manifest.json in the reproduction directory; instructions are in S1 Text.

S1 Table. Baseline characteristics of the 367 TCGA-COADREAD, 561 GSE39582 and 232 GSE17538 overall-survival populations.

S2 Table. Population and endpoint variants: full denominator, positive-time sensitivity and complete-stage sensitivity.

S3 Table. Score registry for the candidate gene-set definitions.

S4 Table. Reference member universe: one record per declared gene, set and cohort with measurement and exclusion flags.

S5 Table. Reference and fitted sample sizes, fitted standard deviations and reconstruction discrepancies.

S6 Table. Primary GPL570 probe mapping.

S7 Table. All 45 candidate-definition by cohort by model overall-survival estimates from the separate R implementation.

S8 Table. OX4 joint proportional-hazards diagnostics.

S9 Table. Raw, subcohort-adjusted and conditional rank correlations.

S10 Table. Diagonal variance fractions and signed covariance-term sums for the equal-weight score.

S11 Table. All 108 reference-pair overlap and correlation rows.

S12 Table. Reference-pair descriptive statistics by cohort.

S13 Table. In-silico two-member coverage-loss diagnostic.

S14 Table. Extended clinical adjustment models with cohort-specific covariates.

S15 Table. R-verified overall-survival time and stage sensitivities.

S16 Table. All 300 recurrence risk-set specifications refitted in R for diagnostic assessment.

S17 Table. Thirty GSE39582/GSE17538 five-year-censored sensitivity rows.

S18 Table. Case-level KEAP1 alteration routes and group definitions.

S19 Table. MSI-adjusted model versions.

S20 Table. Published uses of the shared five-gene seed list and their downstream analysis objects.

S21 Table. Supplementary-file cell anchors for the published-use audit.

S22 Table. Targeted 20-record source inventory: access, publication class, inspected readouts, source anchors and interpretation limits.

S23 Table. Six matched-background summaries.

S24 Table. All 6,000 matched random-set models refitted in R.

S25 Table. GSE71187 platform-sensitivity recheck results.

S26 Table. GSE71187 probe coverage.

S27 Table. Bibliographic identity, publication-notice checks and claim-support scope for the 17 main-manuscript references.