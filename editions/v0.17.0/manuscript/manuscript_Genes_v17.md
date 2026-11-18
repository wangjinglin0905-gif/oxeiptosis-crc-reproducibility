Compartment association and count-split stability of an oxeiptosis-associated transcript score in colorectal cancer

Jinglin Wang¹˒², Pan Sun¹, Wei Zhu¹, Yan Zhang¹, Yao Tang¹ and Juan Yang¹,*

¹ Department of Geriatrics, The Affiliated Hospital of Guizhou Medical University, Guiyang, China

² Department of Digestive Disease, Guizhou Provincial People's Hospital, Guiyang, China

* Correspondence: Juan Yang, MD, PhD (yj63yj63@163.com). ORCID: 0009-0001-2081-6950.

Article

Abstract

Background/Objectives: Does a pathway-associated expression score distinguish tumour from stroma when its local ranks are sensitive to sampling depth? We examined this question for OX4, the mean of four standardised transcripts (KEAP1, PGAM5, AIFM1 and OTUD1), in colorectal cancer.

Methods: Clinical association analyses used TCGA-COADREAD, GSE39582 and GSE17538 (367, 561 and 232 patients, respectively). Sequential spatial analyses examined five GSE280315 sections, four GSE226997 patients and a separate six-patient dataset with two sections per patient and source-provided hematoxylin-and-eosin annotations. We compared patient-level tumour–stroma effects, conditional count-split correlations, epithelial and stromal references, and expression-matched gene sets.

Results: Adjusted overall-survival hazard ratios were 0.939 (95% confidence interval 0.753–1.171), 0.939 (0.809–1.091) and 1.175 (0.935–1.476), with effects near the centres of matched backgrounds. In the morphology dataset, OX4 was higher in tumour in all six patient averages (mean Cliff's delta 0.239; patient-bootstrap 95% interval 0.139–0.325). Epithelial and stromal references showed expected opposing directions. OX4 section-wise median count-split correlations were 0.015–0.159, extending the low absolute stability observed in the earlier spatial analyses. Absolute compartment effects occupied the 28th–96th percentiles of expression-and-detection-matched backgrounds; the mean effect decreased to 0.114 after conditioning on total counts.

Conclusions: OX4 distinguished tumour from stroma modestly at the patient level, even though individual spots ranked inconsistently after count splitting. Compartment effects and local rank stability answer distinct questions about a transcript score. Functional interpretation requires a separate measure of cell death.

Keywords: oxeiptosis; colorectal cancer; gene-set scoring; spatial transcriptomics; histological compartments; count-split stability; matched controls

1. Introduction

Holze and colleagues identified oxeiptosis by tracing how cells respond to oxidative stress. Their experiments linked KEAP1, PGAM5 and AIFM1 to caspase-independent death through perturbations, phosphorylation measurements and cellular phenotypes [1]. Subsequent studies linked OTUD1 to KEAP1-mediated stress responses [2] and implicated the KEAP1–PGAM5–AIFM1 axis in sanguinarine-treated colorectal cancer models [3]. In archived tumours, steady-state RNA abundance offers a different measurement. What can an average of these transcripts tell us about the tissue?

Small gene lists enter public-data analyses in several ways. They can identify mechanistic participants, seed the selection of correlated transcripts, or define an enrichment score. Oxeiptosis-associated long noncoding RNA (lncRNA) studies in endometrial and gastric cancer used pathway genes to derive other prognostic variables [4, 5]. A five-gene list also appeared in a treatment-response analysis in the NeoSAC trial and in broader programmed-cell-death analyses [6, 7, 8]. Reuse of the same members therefore does not establish that these studies measured the same quantity.

Small transcript scores may share variation with broader expression programmes. Random-signature analyses in breast cancer showed that many unrelated gene sets could associate with outcome, partly through common transcriptional structure [9]. This motivates explicit background comparisons when evaluating pathway-associated scores. Such comparisons place an observed clinical association in context, while analyses of score composition address what the expression construct captures.

We first asked whether an equally weighted four-gene score showed a consistent survival association across colorectal cancer cohorts, relative to matched gene sets. Spatial data then allowed us to ask how well individual tissue locations retained their ranks after count splitting. Finally, we examined whether tumour and stroma differed in score even when local ranks were unstable. For this comparison we used source-provided hematoxylin-and-eosin (H&E) annotations, epithelial and stromal reference scores, and expression-and-detection-matched gene sets.

2. Materials and Methods

2.1. Study design and analysis populations

We analysed expression and clinical records from The Cancer Genome Atlas colon and rectal adenocarcinoma cohort (TCGA-COADREAD) and two Gene Expression Omnibus (GEO) series. Marisa et al. described the GSE39582 cohort [10] and Smith et al. described the cohorts represented in GSE17538 [11]. For TCGA, we used the UCSC Xena HiSeqV2 matrix and the accompanying archived clinical fields listed in Text S1. These fields differ from the endpoints in the harmonised TCGA Clinical Data Resource [12]. We retained each cohort's processed expression scale and fitted models separately.

For TCGA, one primary-tumour aliquot per patient was selected, preferring sample vial 01A when available. Patients were linked by the 12-character TCGA identifier. The frozen analysis contained 367 patients with positive overall-survival time, an event indicator and required clinical data. GSE39582 contributed 566 tumour expression profiles and 561 patients with analysable overall survival. GSE17538 contributed 232 eligible human tumour profiles after excluding non-target samples from the SuperSeries. GSE17536 is a component of GSE17538 and was not treated as an additional validation cohort. Patient identifiers, inclusion flags and cohort-specific denominators are retained in the analysis package.

The primary overall-survival analysis retained the full GSE39582 denominator, including six records with zero recorded follow-up, four of which were deaths. Excluding these records was examined as a sensitivity analysis. Their exact temporal interpretation cannot be established from the available public metadata. Missing stage was represented by a separate Unknown category. This coding preserves those patients but does not eliminate bias from missingness. Other covariates were required to be observed. Table 1 presents baseline characteristics, with cohort-level counts also given in Table S1; Table S2 records population and endpoint rules.

Table 1. Baseline characteristics of the overall-survival analysis populations.

| Cohort | n | Deaths | Age, years; Median (IQR) | Male, n | Stage I / II / III / IV / Unknown |
| --- | --- | --- | --- | --- | --- |
| TCGA-COADREAD | 367 | 85 | 66 (55.5–75) | 202 | 56 / 132 / 111 / 51 / 17 |
| GSE39582 | 561 | 191 | 68 (59–76) | 308 | 31 / 262 / 204 / 60 / 4 |
| GSE17538 | 232 | 93 | 65.5 (56–74) | 122 | 28 / 72 / 76 / 56 / 0 |

n, patients; IQR, interquartile range. Stage counts are shown in the order indicated; Unknown denotes missing stage.

2.2. Candidate definitions and score construction

The primary score, OX4, contained KEAP1, PGAM5, AIFM1 and OTUD1. OX3 contained the first three genes. OX5 contained KEAP1, PGAM5, AIFM1, NFE2L2 and AIRE; OX5 without NFE2L2 omitted that member. OX6 expanded OX4 with CUL3 and EZH2. We applied the same equal-weight procedure to all five definitions to isolate sensitivity to member selection. The positive weights define a transcript average; they do not encode the direction or strength of each member’s mechanistic contribution.

Array probes were mapped to gene symbols with the archived platform annotation, and available probe values for each gene were averaged. The seven GPL570 probes mapping to the four primary members were unambiguous in that annotation. The handling of multi-symbol probes in larger reference sets follows the recorded mapping code and is a limitation of those reference constructions. All four primary members were required to be measurable before a score was computed. Probe mappings and exact lists are provided in Tables S3 and S6.

Within each reference population, each gene was centred and divided by its population standard deviation (SD). The standardised values were averaged, and the composite was centred and divided by its sample SD. Reference populations comprised 367 TCGA patients, all 566 eligible GSE39582 tumour profiles, and 232 GSE17538 tumour profiles. Hazard ratios (HRs) are reported per one reference-population composite SD. In the GSE39582 survival subset, score SDs were approximately 0.996–1.004. These unsupervised transformations were independent of outcome and were performed separately within each cohort.

2.3. Survival models and sensitivity analyses

We fitted Cox proportional hazards models with Efron’s method for tied event times. M0 included the score alone. The primary model, M1, included the score, continuous age, sex and categorical stage. The R formula was Surv(OS_time, OS_event) ~ OX4 + age + male + stage4; OS_time and OS_event denote overall-survival time and event status, respectively. Stage I was the reference where present, and unused levels were dropped. M2 additionally included eight transcript references: ferroptosis, apoptosis, necroptosis, parthanatos, mitophagy, oxidative-stress response, proliferation and mitochondrial transcripts. M2 was an exploratory model for conditional association.

The reference lists came from the archived FerrDb, KEGG, Gene Ontology and MitoCarta-derived files. The nine-set correlation panel additionally included the 14-member SINGH_NFE2L2_TARGETS list. This list belongs to the MSigDB C2 chemical and genetic perturbations collection and derives from an NRF2-knockdown study in lung cancer [13]. It does not contain NFE2L2 itself. The mitochondrial list corresponds to the 1,136-member MitoCarta3.0 human inventory before measurement filtering [14]. It is used here as a transcript reference, not as a measure of mitochondrial number, mass or function. Measured reference members and exclusion flags are supplied in Table S4; reference-population sizes, fitted score SDs and reconstruction checks are in Table S5. Text S2 provides checksums for the released reproduction inputs and expected outputs.

To reduce direct member overlap, the primary exploratory M2 and reference-pair analyses removed a fixed six-gene list, KEAP1, PGAM5, AIFM1, OTUD1, CUL3 and EZH2, from the reference sets. For OX5, a separately labelled sensitivity also excluded NFE2L2 from the relevant competing references. The SINGH comparison did not require such removal. Models with correlated reference scores and few events per fitted parameter were interpreted cautiously. Parameter counts and model warnings accompany the estimates in Tables S7 and S8.

Proportional hazards were assessed for OX4 using the global rank-transformed cox.zph test. Confidence intervals (CIs) were 95% model-based intervals. P values in the bulk analyses were two-sided, nominal and unadjusted for multiple comparisons. The primary clinical estimands were the three cohort-specific OX4 M1 associations. Alternative definitions, reference adjustment, clinical and missing-stage sensitivities, altered-gene groups and recurrence analyses were exploratory and reported without significance-based selection. The study evaluated continuous-score associations and did not develop a prognostic cutpoint or trained prediction model.

2.4. Matched backgrounds and transcriptional diagnostics

Two matched random-set designs each contained 1,000 four-gene sets per cohort, with seed 20260911. The first matched each OX4 member to a candidate in the same expressed-status, mean-expression-decile and variability-decile stratum. Expressed status was defined using the cohort-specific 25th percentile of gene means. Original OX4 members were excluded from the candidate pool, and genes were sampled without replacement within a set. The second design additionally constrained each candidate's mean absolute Spearman correlation with the original OX4 members to within 0.1 of its target. This constraint matches candidate-to-OX4 correlation and does not match the internal pairwise correlations of the newly drawn set. Gene length and other unrecorded attributes were not matched.

Each random score underwent the same standardisation and M1 adjustment as OX4. We report the fraction of random absolute log hazard ratios below the observed absolute log hazard ratio. The associated upper-tail fraction was calculated as one plus the number at least as large as the observed value, divided by one plus the number of draws. These quantities describe the specified empirical backgrounds; they are not equivalence tests or tests of cell-death specificity. Membership was kept fixed during verification. The matched analyses were added during exploratory development; no external preregistration was performed.

Spearman correlations quantified associations between OX4 and transcript references. Partial rank correlations were calculated as Pearson correlations between residuals of ranked variables after linear adjustment. Subcohort labels entered the adjustment as indicator variables. These analyses quantify conditional association within the recorded cohort structure.

For the nine reference sets, all 36 unordered set pairs were evaluated separately in each expression reference population. Jaccard overlap was calculated from measured members after the fixed six-gene exclusion; score correlation used the corresponding expression profiles. The 108 cohort-by-pair rows share genes, reference scores and patients within cohorts and repeat the same reference pairs across cohorts. We therefore present them descriptively as cohort-faceted summaries rather than as pooled independent observations. A two-member KEAP1–AIFM1 score was also compared with OX4 as an in-silico coverage-loss diagnostic. Squared Spearman correlation is reported as a descriptive overlap measure, not as explained variance of the original expression scores.

2.5. Source audit and computational verification

The literature component was a targeted audit of 20 records with oxeiptosis in the title, supplemented by primary mechanistic, cohort and gene-set sources. Full texts, figure legends and supplements were inspected where available. Abstract-only records were not assigned negative measurement findings. For five published uses of the common five-gene list, we recorded member identities, the downstream analysis object and the described scoring procedure. This targeted audit contextualised score definitions and was not designed as a systematic review. Table S27 links each reference to its role in this study.

Candidate and random scores were reconstructed from archived expression matrices using Python. A separate R implementation refitted 45 overall-survival models, 300 recurrence specifications and 6,000 matched random M1 models, using R 4.6.1 and survival 3.8-6. The original candidate-model estimates used Python lifelines, whereas the random-background pipeline used a custom Python Cox routine. R estimates are reported in the principal survival tables. Convergence warnings were retained, and recurrence fits were classified as risk-set diagnostics. The public reproduction entry point covers these 6,345 fits from supplied analysis frames and fixed random-score matrices; other exploratory results are documented in the supporting tables.

2.6. Exploratory spatial extension and input provenance

The spatial extension was conducted after completion of the bulk analysis. It used five GSE280315 Visium HD sections (P1CRC, P2CRC, P5CRC, P3NAT and P5NAT) from the published colorectal study [15] and four GSE226997 patients [16]. The five HD sections represent three tumour and two normal-adjacent sections; their labels indicate four donors, including a tumour/normal-adjacent pair. Sections, regions and technical splits were not treated as independent patients. GSE226997 expression tables, coordinates and H&E images were obtained from STimage-1K4M [17], fixed at repository commit b2073c6e8bb7e6d8dfd15c0db6134ba4f2a8a84f. All four complete filtered_feature_bc_matrix.h5 files were recovered from bounded prefixes of the original GEO archives. We compared every count entry, feature position and barcode with the derivative tables, documented duplicate-symbol suffixes, and verified the depth denominator. We did not reprocess FASTQ files or download and validate the complete archives.

For GSE280315, 8-μm bins were summed into square regions of 32, 64, 128 and 256 μm, with 64 μm designated as the primary scale before the stability analyses. A region required at least 50% of its grid bins to be represented and at least 500 total counts. GSE226997 was analysed at native Visium spots with at least 500 counts across the features represented in each matrix. The platforms and sampling geometries were analysed separately. Zero member counts were retained; joint detection of all four members was a reported property, not an inclusion criterion.

Spatial OX4 was defined as the mean of four gene-wise z scores after log1p transformation of counts per 10,000 (CP10K), using each section's eligible units as the reference and sample SD for each gene. This spatial analogue was not a transfer of the bulk-score calibration. For each technical split, member counts and counts for the remaining represented genes were allocated between two halves by binomial sampling with probability 0.5. Each half used its own total count for normalisation and the full-data gene means and SDs. We reported Spearman correlations between the two half-scores and, as a sensitivity analysis, correlations between rank residuals after conditioning on both half-depth ranks. Working stability thresholds were a median correlation of at least 0.7 and a 10th percentile of at least 0.6; these were operational decision rules, not established biological cutoffs. Depth associations and joint detection were reported separately. The initial GSE280315 screen additionally required at least 50 regions, joint detection of at least 80% and absolute score–depth correlation no greater than 0.5 for its composite feasibility gate (Text S3).

The GSE280315 scale pilot used 30 splits per section and scale. The matched comparison used 100 unique four-gene sets per section, each evaluated with 10 splits, alongside OX4 and a stromal construct. Candidates matched each OX4 member's mean CP10K within a factor of two and detection rate within 0.05. Original members, reference members, duplicate symbols and MT-/RPL/RPS-prefixed genes were excluded. Up to 30 nearest candidates per member were retained; each pool required at least 10 candidates. The background percentile was the fraction of matched sets with a lower median split correlation than OX4, not a P value. Matching did not control all covariance or spatial structure. GSE226997 used 30 splits per complete construct and did not repeat the matched-set analysis. Seeds, full matching rules and provenance are in Text S3 and Tables S28–S31. The complete epithelial construct was unavailable in GSE280315 because KRT18 was not represented, and in three GSE226997 samples because KRT19 was not represented; it was not shortened to fit the measured features.

The predefined stromal construct comprised COL1A1, COL1A2, COL3A1 and DCN; the epithelial construct comprised EPCAM, KRT8, KRT18 and KRT19. These investigator-defined marker averages served as technical references, not mechanism-positive controls or independently validated histological labels. Incomplete constructs were not scored. Image–coordinate overlays were checked for gross alignment. No new human pathology annotation was performed for these two datasets. The published GSE226997 compartment definitions combined expression clustering and H&E morphology [16] and were not used as expression-independent ground truth. Technical-split distributions in these two datasets were reported descriptively, without patient-level significance tests or pooled cross-platform effect estimates.

2.7. Morphology-defined compartment analysis

A subsequent analysis used the spatial colorectal dataset of Valdeolivas et al. [18], fixed at Zenodo record 7760264. The deposited H&E-based spot annotations provided a morphological reference independent of the OX4 score. Six tumour-bearing patients contributed two sections each; a seventh, non-neoplastic tissue source lacked the target tumour/stromal contrast and was excluded from the primary comparison. Source annotations distinguished tumour, fibroblastic stromal categories with differing immune content, mixed tissue and other compartments. The stromal comparison is therefore a morphological compartment contrast rather than a pure fibroblast comparison. No new pathologist rating, blinding procedure or inter-rater assessment was performed by our team.

The analysis contract was fixed after input eligibility checks and before calculating scores in this dataset. The broader study and its sequential extensions were exploratory and were not externally preregistered. Spots required in-tissue status and at least 500 total Gene Expression counts; source-excluded and blank annotations were then removed. All remaining spots, including mixed and other tissue, defined section-specific normalisation and technical stability. Only frozen tumour and stromal categories entered the primary compartment contrast. All four members of each construct were present; zero counts were retained. Gene-wise log1p(CP10K) z scores were averaged with equal positive weights and scaled by the full-section composite sample SD.

The primary compartment effect was Cliff's delta for tumour versus stroma, computed as twice the rank area under the curve (AUC) minus one with ties retained. A positive effect indicates higher tumour scores; this AUC is an effect-size representation, not the performance of a trained classifier. The two section effects were averaged equally within each patient. Uncertainty was described by 10,000 patient bootstrap resamples and exploratory two-sided exact sign-flip tests across the 64 sign combinations of six patient effects, with Benjamini–Hochberg (BH) correction across OX4 and the two reference constructs. The reference constructs checked the expected direction of each compartment contrast. These tests assume sign symmetry under the null; no spot-level significance tests were used. Patient and paired-section direction agreement and leave-one-patient means were also reported.

Each section included 100 unique expression-and-detection-matched four-gene sets and 30 conditional count splits for all 103 constructs. Matching bounds followed the earlier spatial comparison, with section-specific candidate selection and the distance rule detailed in Text S4. Percentiles used the fraction of matched values less than or equal to the OX4 value, for absolute compartment effects and median split correlations. Matched sets were not pooled across patients because their members differed by section. Split correlations were calculated across all eligible spots and separately within tumour and stroma. All constructs shared the same count split; both halves used the full-data standardisation reference. No spatial smoothing or additional scale selection was performed.

Depth sensitivity regressed the full score on log1p(total counts) within each section and repeated the compartment contrast on residuals. This differs from conditioning technical half-score correlations on half-depth ranks in the preceding analyses. A second sensitivity restricted stroma to the source fibroblastic immune-content (IC)-low/medium classes, requiring at least 20 such spots and 30 tumour spots per section; both sections had to be eligible for a patient summary. Depth conditioning and restricted stromal labels were sensitivity analyses rather than causal adjustments or new independent cohorts. Source checks, seeds, software versions and all patient/section results are supplied in Text S4 and Tables S32–S37.

2.8. Data access dates

Based on retained project records, the research access dates were reconstructed as 9 September 2026 for TCGA-COADREAD and 11 September 2026 for GSE39582. GSE17538 and GSE71187 were downloaded for this study on 11 September 2026.

The spatial extension accessed GSE280315 and GSE226997 resources for research on 17 November 2026. Acquisition records distinguish the STimage derivative files from the subsequently retrieved original GEO processed matrices. The Zenodo 7760264 source data and annotations were accessed on 17 November 2026; the compartment analysis was performed on 18 November 2026.

2.9. Use of AI-assisted tools

WorkBuddy and OpenAI Codex assisted with code development, source retrieval, numerical checks and manuscript preparation. Codex also assisted with the spatial analysis and plotting code for Figures 1–5. Kimi assisted with manuscript review and English editing. Numerical checks used frozen result tables and, for the bulk survival models, a separate R implementation. Source checks used bibliographic records, primary texts and supplements where available.

3. Results

3.1. Overall-survival estimates varied in direction across cohorts

The primary analyses included 1,160 patients and 369 deaths, analysed separately by cohort. Under M1, the OX4 HRs were 0.939 (95% CI 0.753–1.171), 0.939 (0.809–1.091) and 1.175 (0.935–1.476) for TCGA-COADREAD, GSE39582 and GSE17538, respectively (Table 2 and Figure 1). All three intervals crossed one. The point estimates differed in direction, and their uncertainty remained compatible with potentially relevant protective or adverse associations.

Figure 1. Overall-survival associations of OX4. Points show cohort-specific Cox estimates, with bars indicating 95% model-based CIs. M0 includes the score only; M1 additionally adjusts for age, sex and stage; M2 adds eight transcript references after the fixed six-gene exclusion. HRs are per reference-population score SD. The dashed line marks HR 1. Counts indicate patients/deaths. The models use the same patients within each cohort. Nominal P values and parameter counts are in Table S7; global proportional-hazards diagnostics are in Table S8.

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

The primary OX4 estimates changed under the reference-adjusted M2 models to 0.941, 1.089 and 0.856, respectively. M2 included 15 fitted coefficients in TCGA and GSE39582 and 14 in GSE17538, corresponding to 5.7, 12.7 and 6.6 deaths per coefficient. Its correlated predictors and limited events constrain interpretation, particularly in the first and third cohorts. Joint proportional-hazards tests for the nine OX4 cohort-by-model fits did not reject at the nominal 0.05 level; diagnostics and parameter counts are reported in Tables S7 and S8.

Alternative member lists also lacked a consistent M1 association across cohorts. Selected unadjusted or M2 estimates had intervals excluding one, including OX5 and OX5 without NFE2L2 in GSE17538 M2. These exploratory findings are reported alongside all 45 models in Table S7.

3.2. Adjusted effects lay within the central ranges of the matched backgrounds

The OX4 absolute log hazard ratios under M1 were 0.0633, 0.0625 and 0.1612 in TCGA, GSE39582 and GSE17538. They lay at percentiles 40.1, 44.5 and 45.5 of the expression-and-variability-matched distributions. Under the candidate-to-OX4 correlation constraint, the percentiles were 41.5, 35.3 and 59.8 (Table 3 and Figure 2). Both designs therefore placed the observed estimate near the middle of their distributions (Table S23).

Figure 2. Matched-background distributions for the OX4 M1 association. Empirical cumulative distributions use 1,000 four-gene sets per design and cohort. The dashed vertical marker denotes OX4. MSD matches expressed status, mean-expression decile and variability decile. MSD + OX4 correlation adds a candidate-to-OX4 correlation constraint. Panels a–c show TCGA-COADREAD, GSE39582 and GSE17538. Percentiles compare the observed effect with the specified backgrounds; the second design does not match within-draw co-expression.

Table 3. Observed adjusted effects relative to matched backgrounds.

| Cohort | Observed absolute log HR | MSD percentile | MSD + OX4 correlation percentile |
| --- | --- | --- | --- |
| TCGA-COADREAD | 0.0633 | 40.1 | 41.5 |
| GSE39582 | 0.0625 | 44.5 | 35.3 |
| GSE17538 | 0.1612 | 45.5 | 59.8 |

HR, hazard ratio. MSD matches expressed status, mean-expression decile and variability decile. The second design adds candidate-to-OX4 correlation matching. Each background contains 1,000 sets per cohort.

All 6,000 matched models were refitted without a Cox convergence warning in R. The largest absolute difference from the stored random-model log HR was approximately 1.97 × 10⁻⁹, indicating numerical agreement between implementations (Table S24).

3.3. The score shared variation with several transcriptional references

OX4 correlated positively with the mitochondrial transcript reference, proliferation reference and SINGH-NRF2 target score in all three cohorts (Figure 3). The respective correlations were 0.386, 0.305 and 0.323 in TCGA; 0.566, 0.520 and 0.473 in GSE39582; and 0.797, 0.793 and 0.704 in GSE17538. The mitochondrial reference excluded the directly shared AIFM1 and PGAM5 members. These associations therefore persisted beyond their direct contribution to both scores.

Figure 3. Transcript-reference correlations and subcohort composition. Panel a shows Spearman correlations of OX4 with mitochondrial, proliferation and SINGH-NRF2 references in the overall-survival populations. Panel b compares raw and subcohort-adjusted rank correlations in GSE17538, using the recorded Moffitt/Vanderbilt label. The mitochondrial reference excludes the directly shared AIFM1 and PGAM5 members. All references are transcript-based scores.

Adjustment for Moffitt versus Vanderbilt labels in GSE17538 reduced mitochondrial, proliferation and NRF2 correlations to 0.587, 0.582 and 0.484, respectively. In GSE39582, adjustment for the discovery-versus-validation label had little effect. Conditional OX4–mitochondrial and OX4–NRF2 rank correlations remained positive after adjustment for the other reference (Table S9). Subcohort comparisons used subsets of the existing cohorts.

For the equal-weight score, each standardised member had the same diagonal variance contribution in the standardisation population. Covariance terms altered the total score variance. In TCGA, diagonal shares were approximately 27.217% per member and signed cross-terms summed to −8.87%. These quantities describe the variance of a correlated average. Small differences in the GSE39582 fitted subset followed its restriction from the expression reference population (Table S10).

Low member overlap did not ensure low reference-score correlation. Among 36 pairs per cohort, 4, 15 and 30 had Jaccard overlap below 0.05 and Spearman correlation at least 0.50 in TCGA, GSE39582 and GSE17538, respectively. Each cohort contained six zero-overlap pairs, with maximum correlations of 0.405, 0.698 and 0.878. Within-cohort descriptive correlations between Jaccard overlap and score correlation were 0.523, 0.501 and 0.645 (Figure 4; Tables S11 and S12). No pair had Jaccard overlap of 0.20 or greater.

Figure 4. Member overlap and reference-score correlation by cohort. Each point represents one of 36 unordered pairs among nine transcript references. Jaccard overlap uses measured members after the fixed six-gene exclusion. Correlations use 367, 566 and 232 expression profiles in panels a–c. Dashed lines mark Jaccard 0.05 and Spearman correlation 0.50. Pairs share reference sets and patients and are shown descriptively. No pair had Jaccard overlap of 0.20 or greater.

An in-silico KEAP1–AIFM1 score correlated with OX4 at 0.770, 0.795 and 0.770 across the three expression populations. Removing PGAM5 and OTUD1 therefore changed patient ranks (Table S13).

3.4. Clinical and endpoint sensitivity analyses

Complete-stage restriction left the OX4 M1 point estimates close to their primary values. Excluding the six zero-time overall-survival records in GSE39582 gave HR 0.950 (95% CI 0.816–1.106). In that cohort, extended adjustment on 499 complete cases gave HR 0.987 (0.838–1.163), while M1 fitted to the same 499 patients already gave HR 0.989. The change from the full-cohort estimate was therefore largely associated with case restriction in this comparison. TCGA used different treatment and pathology variables, and the corresponding complete-case model contained only 33 deaths; those exploratory results are documented separately in Tables S14 and S15.

In GSE39582, the stage II–III, positive-time, five-year-censored recurrence sensitivity gave OX4 HR 0.976 (95% CI 0.821–1.160; 453 patients, 130 events). Applying the same restriction to GSE17538 gave HR 1.393 (0.914–2.123; 145 patients, 31 events). The GSE17538 restriction approximates a harmonised risk set because treatment and disease-free-state information were incomplete. Tables S16 and S17 report the diagnostic specifications and sensitivities. The TCGA new-tumour-event field was excluded from recurrence interpretation.

Among the 367 analysed TCGA patients, nine had a recorded KEAP1 alteration in Table S18; six of these met the specified nonsynonymous-mutation or deep-deletion definition and had higher mean OX4 values than the comparison group. This exploratory comparison did not meet the minimum-group-size feasibility criterion. Alteration calls alone did not establish functional protein inactivation. Adding genomic microsatellite instability (MSI) status to the age-, sex-, stage-, proliferation- and mitochondrial-reference-adjusted model gave OX4 HR 0.937 (95% CI 0.730–1.203; Tables S18 and S19).

3.5. Source verification distinguished member lists from measured constructs

Five inspected publications used the same five seed symbols after NRF2 was harmonised to NFE2L2. Their downstream uses differed. The endometrial and gastric studies derived lncRNA prognostic signatures; the Ma study derived a three-gene model from a broader programmed-cell-death pool; and the Du study used a broader cell-death analysis and a machine-learning index. NeoSAC explicitly reported an oxeiptosis process score for a treatment-response comparison, but the available description did not provide enough detail to reproduce every scoring parameter. None of these source descriptions was substituted for the equal-weight OX5 implementation in this study. Two publications attributed the seed set to the same earlier source [19]; provenance for the other three was not established from the inspected material (Tables S20 and S21).

Of the 20 targeted records, ten had accessible full texts, eight had abstracts only and two lacked usable abstracts or full texts. Five records were reviews and one was commentary. The two records without usable source material were excluded from measurement-content assessment.

The inspected mechanistic studies measured different levels of the pathway. The sanguinarine colorectal cancer and vitiligo studies included protein abundance and phosphorylation measurements alongside functional assays [3, 20]. In the methyl-cantharidimide liver-cancer study, phosphorylation changed while the displayed measurements showed little change in total AIFM1 and PGAM5 abundance [21]. Its discussion also described reduced PGAM5 protein, an unresolved discrepancy recorded in Table S22. The targeted sources provide context for molecular readouts but no patient-level calibration of OX4 in the three analysed cohorts.

3.6. Spatial OX4 had limited absolute count-split stability without an additional deficit against matched sets

In GSE280315, the 30-split pilot yielded median OX4 correlations of 0.358, 0.501, 0.253, 0.358 and 0.285 at 64 μm for P1CRC, P2CRC, P5CRC, P3NAT and P5NAT, respectively. All were below the working median threshold. None of the 20 section-by-scale combinations met every composite feasibility condition. P2CRC reached a median correlation of 0.822 at 256 μm, but its score–depth correlation was 0.729 (Table S28).

In the separately seeded matched comparison, OX4 ranked at the 66th–95th percentiles of the matched backgrounds and exceeded the matched-set median in every section (Table 4). After conditioning on half-depth ranks, OX4 ranked at the 31st–90th percentiles, with no section below the 10th percentile. The predefined criterion for an additional deficit, at least three of five sections below the matched 10th percentile, was therefore not met. Stromal median correlations were 0.955–0.964 in tumour sections and 0.649 and 0.845 in the two normal-adjacent sections. The five supplied matrices lacked the exact KRT18 feature symbol, so the complete epithelial construct was not calculated. These findings separate limited absolute stability from evidence of an OX4-specific deficit.

Table 4. Spatial count-split stability relative to matched backgrounds in GSE280315.

| Section | OX4 median rho | Matched-set median rho | OX4 percentile | Stromal median rho |
| --- | ---: | ---: | ---: | ---: |
| P1CRC | 0.361 | 0.244 | 95 | 0.964 |
| P2CRC | 0.496 | 0.445 | 66 | 0.962 |
| P3NAT | 0.357 | 0.328 | 87 | 0.649 |
| P5CRC | 0.250 | 0.192 | 84 | 0.955 |
| P5NAT | 0.285 | 0.258 | 74 | 0.845 |

The primary 64-μm analysis used 10 splits per construct and 100 matched sets per section. Percentiles are descriptive background ranks. The initial 30-split pilot used a separate random seed; its results are retained in Table S28. CRC, colorectal cancer; NAT, normal-adjacent tissue.

3.7. Original-file verification supported the GSE226997 technical replication

All 294,437,824 count entries agreed between the four GEO-deposited filtered HDF5 matrices and their STimage tables. Removing the STimage sample prefix gave identical barcode order. The total-count denominator also agreed for every spot. Duplicate symbols received suffixes without count aggregation; these affected neither OX4 nor the predefined reference members. Three original matrices contained 17,944 features and lacked KRT19; GSM7089857 contained 19,144 features and included KRT19. Thus, the epithelial coverage gap preceded the STimage representation (Table S31).

Of 16,096 native spots, 15,502 met the depth criterion. OX4 median split correlations were 0.181–0.328, with joint member detection in 0.86%–6.16% of eligible spots (Table 5). Conditioning on half-depth ranks reduced the median correlations to 0.023–0.091; full-score correlations with depth were 0.388–0.534. Stromal median correlations were 0.778–0.937. The only complete epithelial construct, in GSM7089857, had median correlation 0.384. All available results, including the low epithelial stability, were retained (Table S30).

Table 5. Native-spot technical stability in GSE226997.

| Sample | Eligible spots | OX4 joint detection, % | OX4 median rho | OX4 P10–P90 | Stromal median rho |
| --- | ---: | ---: | ---: | ---: | ---: |
| GSM7089855 | 4142 | 5.87 | 0.206 | 0.193–0.209 | 0.937 |
| GSM7089856 | 3588 | 6.16 | 0.322 | 0.312–0.327 | 0.894 |
| GSM7089857 | 4508 | 1.46 | 0.181 | 0.168–0.194 | 0.856 |
| GSM7089858 | 3264 | 0.86 | 0.328 | 0.319–0.343 | 0.778 |

Each complete construct used 30 splits. Joint detection was calculated among eligible spots and was not an inclusion filter. P10–P90 describes technical repetitions, not a patient-level confidence interval.

All four image–coordinate overlays passed the gross alignment inspection, but no independent histological compartment evaluation was performed. These analyses establish correspondence of the processed expression inputs and characterize technical splitting behaviour. They do not establish biological reproducibility between sections or spatial localization of oxeiptosis activity.

3.8. Morphological compartment association coexisted with limited local count-split rank stability

The independent morphology analysis included 16,225 eligible spots from 12 sections of six patients, with 4,027 tumour and 5,221 stromal spots in the primary comparison. OX4 was higher in tumour in the two-section averages of all six patients (mean Cliff's delta 0.239; patient-bootstrap 95% interval 0.139–0.325; Table 6 and Figure 5a). Five patients had same-direction section effects; A551763 had effects of −0.049 and 0.092, giving a near-zero patient average of 0.022. Leave-one-patient mean effects ranged from 0.204 to 0.282. Epithelial and stromal references showed the expected opposite directions in all sections, with patient mean effects of 0.553 and −0.769, respectively. All three exact two-sided sign-flip P values and three-score BH q values were 0.03125, the smallest attainable two-sided value for six nonzero, same-direction effects. The small patient count limits inference despite the large number of spots.

Table 6. Morphological compartment effects and conditional count-split stability.

| Construct | Mean patient delta | Patient-bootstrap 95% interval | Patient direction | Section median split rho range |
| --- | ---: | --- | --- | --- |
| OX4 | 0.239 | 0.139 to 0.325 | 6/6 higher in tumour | 0.015–0.159 |
| Epithelial | 0.553 | 0.352 to 0.751 | 6/6 higher in tumour | 0.243–0.886 |
| Stromal | -0.769 | -0.897 to -0.631 | 6/6 higher in stroma | 0.491–0.971 |

Patient effects average two section-level Cliff's deltas with equal weights. Positive delta denotes higher tumour scores. OX4 direction was assessed exploratorily; only the epithelial/stromal reference directions were expected. Bootstrap intervals resample six patients and do not incorporate uncertainty from all analytic choices. Split-correlation ranges describe the 12 section medians of 30 technical splits, not 12 independent patients.

OX4 section-wise median split correlations ranged from 0.015 to 0.159, compared with 0.243–0.886 for the epithelial and 0.491–0.971 for the stromal references (Figure 5b). OX4 correlations within tumour or stroma were also low. Patient summaries of the median absolute change in compartment effect between complete and half-count data ranged from 0.020 to 0.053. Thus, a modest group-level distributional difference persisted despite unstable spot ranks under count splitting. Both half-scores used the same section-specific standardisation reference.

After depth conditioning, the mean OX4 compartment effect decreased to 0.114 and remained positive in five of six patients (Figure 5c). Restricting stroma to the IC-low/medium classes retained positive patient-average effects in all five eligible patients (0.097–0.360); A121573 lacked those classes in both sections and was retained as unavailable for this sensitivity. The absolute OX4 compartment effect occupied the 28th–96th percentiles of section-specific matched backgrounds, with only one section at or above the 95th percentile. Its median split correlation occupied the 27th–93rd percentiles (Figure 5d). These comparisons did not show a consistently exceptional compartment effect or a universal additional stability deficit for OX4 under the specified matching rules.

Figure 5. Morphological compartment effects and local count-split stability. (a) Each point is a patient's mean tumour-versus-stroma Cliff's delta across two sections; horizontal bars show the six-patient mean. (b) Each point is the average of two section-wise median correlations from 30 count splits; bars show the six-patient mean. (c) Paired patient OX4 effects before and after linear conditioning on log1p(total counts). (d) Each point is one of 12 sections, showing the percentile of the OX4 absolute compartment effect and median split correlation among 100 expression-and-detection-matched four-gene sets. Dashed lines mark the 50th percentiles. Sections are not independent patients and matched sets differ between sections. Source H&E annotations define the morphological comparison; no new human pathology review was performed. Count splits are conditional technical repetitions.

4. Discussion

In each of the six patients with histological annotations, OX4 was higher on average in tumour than in stroma. Yet the same spots changed rank substantially when we split their counts in half. A shift between the distributions of two tissue compartments can therefore remain visible even when a high-scoring spot is difficult to identify reproducibly within a section.

The size of the compartment difference helps put that result in perspective. One patient had a near-zero average difference across sections. Conditioning on total counts reduced the six-patient mean, and only one section reached the 95th percentile of its matched background. Tissue composition and RNA yield are both plausible contributors to the contrast.

Holding the scoring algorithm constant clarified how member choice changes the measurement. OX4 averages standardised transcript abundances with positive weights. Protein interactions, phosphorylation and the timing of oxidative stress are outside that definition. Substituting members changes the construct even when each candidate is described as oxeiptosis-associated. Likewise, the published lncRNA models and enrichment procedures represent distinct uses of seed lists. Direct comparisons require both the member list and the scoring algorithm.

OX4 shared variation with mitochondrial, proliferation and NRF2-related transcript references. The attenuation after GSE17538 subcohort adjustment illustrates how cohort composition can contribute to a strong transcript-level association. Institutional labels combine technical and biological differences, so the adjustment cannot distinguish their origins. Removing shared genes addresses direct overlap but leaves broader expression covariance and cell-composition effects. Accordingly, partial correlations and M2 coefficients describe conditional associations rather than causal sources of OX4 variation.

Reference-set overlap and score correlation captured different properties of the comparator panel. Some zero-overlap pairs remained strongly correlated, even though overlap and correlation were positively associated within each cohort. Both quantities should therefore be reported with their measurement universe and cohort. Because pairs share sets and patients, the 108 comparisons provide descriptive context rather than 108 independent tests.

Matched backgrounds placed the adjusted OX4 effect near the middle of distributions from four-gene sets with similar expression and variability. Adding a candidate-to-OX4 correlation constraint produced the same broad pattern. This extends the comparison logic of random-signature studies [9] to an oxeiptosis-associated candidate in colorectal cancer. The inference remains conditional on the matching rules: it is neither an equivalence result nor a test of pathway function. A calibrated score or a stress-dependent association could behave differently.

The spatial references help interpret low split correlations. Stromal markers retained spot ranks more consistently under the same processing scheme. Many four-gene sets matched to OX4 for expression and detection also had low correlations. The one GSE226997 sample with complete epithelial coverage showed that a marker set can rank spots poorly even when its tissue role is familiar. Both the absolute correlation and the choice of comparator affect the interpretation.

Original-file tracing resolved a separate uncertainty: the GSE226997 count representation did not introduce the missing KRT19 feature or alter the measured counts. Feature universes and archived processing metadata nevertheless differed among patients. Such differences should be reported alongside score coverage rather than interpreted as biological absence. The GSE226997 expression-assisted compartments could not supply an expression-independent morphology test [16]. The separate H&E-annotated dataset provided that reference [18], although its stromal classes contain varying immune mixtures and do not identify a pure cellular source of OX4.

Readers can assess a small transcript score more directly when its members, weights, standardisation population and platform coverage are specified. Reporting the outcome definition, reference-score correlations and position among matched sets then places a clinical association in context. For spatial data, a compartment contrast and a split-count correlation describe different scales of measurement. To calibrate OX4 against oxeiptosis itself would require paired RNA and stress-dependent death measurements capable of distinguishing alternative death processes.

4.1. Limitations

The study used retrospective processed bulk-expression data with incomplete treatment, molecular and cell-composition information. The primary cohorts differ in clinical composition, assay and follow-up, and residual confounding is plausible. We did not estimate tumour purity or test a cell-type-specific model. Missingness and zero follow-up records require assumptions, particularly for recurrence endpoints. M2 and some sensitivity models have limited events relative to their parameter counts, and model-based intervals may understate uncertainty from analytic choices.

Candidate and reference lists were investigator specified, and some sensitivity analyses were added after initial results without external preregistration. Background comparisons cover finite gene universes and matching choices. Within-cohort standardisation does not provide a locked assay for individual patients, and survival association does not establish functional calibration. Confidence intervals remained compatible with potentially relevant effects. The targeted literature audit was access limited and AI assisted, without completed independent human ratings. The small fourth-platform analysis was exploratory, and the two-member simulation did not reproduce an empirical cross-platform comparison (Tables S25 and S26).

The spatial extension was exploratory and followed the bulk results. It included few patients, different sampling geometries and targeted feature universes. Random count splitting reduces the counts available to each half and measures a conditional technical property, not repeat-section or biological reproducibility. The split halves share the full-data standardisation reference and arise from one observed matrix. Depth can reflect tissue composition as well as technical yield, so conditional correlations do not identify causal artifacts. Spatial matched backgrounds controlled marginal expression and detection but not all gene–gene or spatial covariance. GSE226997 did not include a separate matched-set analysis, independent pathology labels or a functional death readout. Original filtered matrices were verified, but complete sequencing archives and original high-resolution image transformations were not reconstructed.

The morphology analysis used six patients with two paired sections each. This sample limits the precision of the patient bootstrap intervals and of the discrete sign-flip tests. Annotations came from the original investigators, with no additional blinded rating in this analysis. A second independently annotated cohort would test the reproducibility of the tumour–stroma contrast. Cell composition, RNA yield and other expression programmes could each contribute to the observed difference.

5. Conclusions

Across six patients, OX4 showed a modest tumour–stroma difference despite low spot-level rank stability after count splitting. In three bulk cohorts, adjusted survival estimates varied in direction and lay near the centres of their matched backgrounds. For this four-transcript score, tissue distribution, robustness to reduced counts and clinical association each provide a distinct piece of evidence.

Supplementary Materials

Texts S1–S4 describe the methods, inputs and verification. Tables S1–S37 contain cohort summaries, model estimates, source evidence and spatial results. Code S1 contains the spatial scripts and derived tables. The accompanying supplement index identifies each file.

Author Contributions

Conceptualization, J.W. and J.Y. Methodology, J.W. Investigation and data curation, J.W., P.S., W.Z., Y.Z. and Y.T. Validation and visualization, P.S., W.Z., Y.Z. and Y.T. Funding acquisition, Y.Z. and J.Y. Project administration, J.W. and J.Y. Supervision, J.Y. Writing—original draft preparation, J.W. Writing—review and editing, J.W., P.S., W.Z., Y.Z., Y.T. and J.Y.

Funding

This work was supported by the National Natural Science Foundation of China (No. 82360467) and the Science and Technology Fund Project of Guizhou Provincial Science and Technology Program (Nos. QKH JC-ZK [2022]-358 and [2022]-348). The funders had no role in study design, data collection and analysis, the decision to publish or preparation of the manuscript.

Institutional Review Board Statement

Each contributing cohort describes its ethics approval and consent procedures in the source publication. We analysed publicly released, de-identified records. We did not contact participants, collect additional specimens or submit a new ethics application for this secondary analysis.

Clinical trial number: not applicable.

Informed Consent Statement

Consent procedures are described in the source cohort publications. The present study used public de-identified data and involved no contact with participants.

Consent for Publication

Not applicable. This manuscript contains no identifiable individual participant information.

Data Availability Statement

Bulk expression and clinical data are available through UCSC Xena (TCGA-COADREAD) and GEO (GSE39582, GSE17538 and GSE71187). Spatial inputs are available through GEO (GSE280315 and GSE226997) and the source morphology archive (https://doi.org/10.5281/zenodo.7760264). The bulk-model inputs, R scripts and expected outputs are archived at https://doi.org/10.5281/zenodo.22882797. The combined reproducibility package contains the additional spatial scripts, derived results and a portable check of the morphology summaries (https://github.com/wangjinglin0905-gif/oxeiptosis-crc-reproducibility/tree/v0.17.0). The version-specific archive is https://doi.org/10.5281/zenodo.24330126.

Acknowledgments

The authors thank the investigators and participants who generated and shared the public datasets analysed in this study. The authors take responsibility for the analyses, interpretation and text, including the use of WorkBuddy, OpenAI Codex and Kimi described in Materials and Methods. Exact tool version identifiers were not retained for all stages.

Conflicts of Interest

The authors declare no competing financial or non-financial interests.

References

1. Holze, C.; Michaudel, C.; Mackowiak, C.; Haas, D.A.; Benda, C.; Hubel, P.; Pennemann, F.L.; Schnepf, D.; Wettmarshausen, J.; Braun, M.; et al. Oxeiptosis, a ROS-induced caspase-independent apoptosis-like cell-death pathway. Nat Immunol. 2018, 19, 130–140. https://doi.org/10.1038/s41590-017-0013-y

2. Oikawa, D.; Gi, M.; Kosako, H.; Shimizu, K.; Takahashi, H.; Shiota, M.; Hosomi, S.; Komakura, K.; Wanibuchi, H.; Tsuruta, D.; et al. OTUD1 deubiquitinase regulates NF-κB- and KEAP1-mediated inflammatory responses and reactive oxygen species-associated cell death pathways. Cell Death Dis. 2022, 13, 694. https://doi.org/10.1038/s41419-022-05145-5

3. Pallichankandy, S.; Thayyullathil, F.; Cheratta, A.R.; Subburayan, K.; Alakkal, A.; Sultana, M.; Drou, N.; Arshad, M.; Tariq, S.; Galadari, S. Targeting oxeiptosis-mediated tumor suppression: a novel approach to treat colorectal cancers by sanguinarine. Cell Death Discov. 2023, 9, 94. https://doi.org/10.1038/s41420-023-01376-3

4. Niu, L.; Wu, Z. Identification and validation of oxeiptosis-associated lncRNAs and prognosis-related signature genes to predict the immune status in uterine corpus endometrial carcinoma. Aging (Albany NY). 2023, 15, 4236–4252. https://doi.org/10.18632/aging.204726

5. Wen, L.; Xu, K.; Huang, M.; Pan, Q. Identification of oxeiptosis-associated lncRNAs and prognosis-related signature to predict the immune status in gastric cancer. Medicine (Baltimore). 2024, 103, e37189. https://doi.org/10.1097/md.0000000000037189

6. Shen, G.; Liu, Z.; Wang, M.; Zhao, Y.; Liu, X.; Hou, Y.; Ma, W.; Han, J.; Zhou, X.; Ren, D.; et al. Neoadjuvant apatinib addition to sintilimab and carboplatin-taxane based chemotherapy in patients with early triple-negative breast cancer: the phase 2 NeoSAC trial. Signal Transduct Target Ther. 2025, 10, 41. https://doi.org/10.1038/s41392-025-02137-7

7. Du, J.O.; Huang, Q.K.; Sun, X.C.; Hu, S.K.; Lin, T.S. Machine learning-driven transcriptomic and single-cell profiling of programed cell death patterns in colon cancer. Sci Prog. 2026, 109, 00368504261432995. https://doi.org/10.1177/00368504261432995

8. Ma, J.Y.; Wang, Y.X.; Zhao, Z.Y.; Xiong, Z.Y.; Zhang, Z.L.; Cai, J.; Guo, J.W. Identification of key programmed cell death genes for predicting prognosis and treatment sensitivity in colorectal cancer. Front Oncol. 2024, 14, 1483987. https://doi.org/10.3389/fonc.2024.1483987

9. Venet, D.; Dumont, J.E.; Detours, V. Most random gene expression signatures are significantly associated with breast cancer outcome. PLoS Comput Biol. 2011, 7, e1002240. https://doi.org/10.1371/journal.pcbi.1002240

10. Marisa, L.; de Reyniès, A.; Duval, A.; Selves, J.; Gaub, M.P.; Vescovo, L.; Etienne-Grimaldi, M.C.; Schiappa, R.; Guenot, D.; Ayadi, M.; et al. Gene expression classification of colon cancer into molecular subtypes: characterization, validation, and prognostic value. PLoS Med. 2013, 10, e1001453. https://doi.org/10.1371/journal.pmed.1001453

11. Smith, J.J.; Deane, N.G.; Wu, F.; Merchant, N.B.; Zhang, B.; Jiang, A.; Lu, P.; Johnson, J.C.; Schmidt, C.; Bailey, C.E.; et al. Experimentally derived metastasis gene expression profile predicts recurrence and death in patients with colon cancer. Gastroenterology. 2010, 138, 958–968. https://doi.org/10.1053/j.gastro.2009.11.005

12. Liu, J.; Lichtenberg, T.; Hoadley, K.A.; Poisson, L.M.; Lazar, A.J.; Cherniack, A.D.; Kovatich, A.J.; Benz, C.C.; Levine, D.A.; Lee, A.V.; et al. An Integrated TCGA Pan-Cancer Clinical Data Resource to Drive High-Quality Survival Outcome Analytics. Cell. 2018, 173, 400–416.e11. https://doi.org/10.1016/j.cell.2018.02.052

13. Singh, A.; Boldin-Adamsky, S.; Thimmulappa, R.K.; Rath, S.K.; Ashush, H.; Coulter, J.; Blackford, A.; Goodman, S.N.; Bunz, F.; Watson, W.H.; et al. RNAi-mediated silencing of nuclear factor erythroid-2-related factor 2 gene expression in non-small cell lung cancer inhibits tumor growth and increases efficacy of chemotherapy. Cancer Res. 2008, 68, 7975–7984. https://doi.org/10.1158/0008-5472.can-08-1401

14. Rath, S.; Sharma, R.; Gupta, R.; Ast, T.; Chan, C.; Durham, T.J.; Goodman, R.P.; Grabarek, Z.; Haas, M.E.; Hung, W.H.W.; et al. MitoCarta3.0: an updated mitochondrial proteome now with sub-organelle localization and pathway annotations. Nucleic Acids Res. 2021, 49, D1541–D1547. https://doi.org/10.1093/nar/gkaa1011

15. Oliveira, M.F.; Romero, J.P.; Chung, M.; Williams, S.R.; Gottscho, A.D.; Gupta, A.; Pilipauskas, S.E.; Mohabbat, S.; Raman, N.; Sukovich, D.J.; et al. High-definition spatial transcriptomic profiling of immune cell populations in colorectal cancer. Nat Genet. 2025, 57, 1512–1523. https://doi.org/10.1038/s41588-025-02193-3

16. Park, S.S.; Lee, Y.K.; Choi, Y.W.; Lim, S.B.; Park, S.H.; Kim, H.K.; Shin, J.S.; Kim, Y.H.; Lee, D.H.; Kim, J.H.; et al. Cellular senescence is associated with the spatial evolution toward a higher metastatic phenotype in colorectal cancer. Cell Rep. 2024, 43, 113912. https://doi.org/10.1016/j.celrep.2024.113912

17. Chen, J.; Zhou, M.; Wu, W.; Zhang, J.; Li, Y.; Li, D. STimage-1K4M: A histopathology image-gene expression dataset for spatial transcriptomics. Adv Neural Inf Process Syst. 2024, 37, 35796–35823. https://doi.org/10.52202/079017-1129

18. Valdeolivas, A.; Amberg, B.; Giroud, N.; Richardson, M.; Gálvez, E.J.C.; Badillo, S.; Julien-Laferrière, A.; Túrós, D.; Voith von Voithenberg, L.; Wells, I.; et al. Profiling the heterogeneity of colorectal cancer consensus molecular subtypes using spatial transcriptomics. NPJ Precis Oncol. 2024, 8, 10. https://doi.org/10.1038/s41698-023-00488-4

19. Zou, Y.; Xie, J.; Zheng, S.; Liu, W.; Tang, Y.; Tian, W.; Deng, X.; Wu, L.; Zhang, Y.; Wong, C.W.; et al. Leveraging diverse cell-death patterns to predict the prognosis and drug sensitivity of triple-negative breast cancer patients after surgery. Int J Surg. 2022, 107, 106936. https://doi.org/10.1016/j.ijsu.2022.106936

20. Kang, P.; Chen, J.; Zhang, W.; Guo, N.; Yi, X.; Cui, T.; Chen, J.; Yang, Y.; Wang, Y.; Du, P.; et al. Oxeiptosis: a novel pathway of melanocytes death in response to oxidative stress in vitiligo. Cell Death Discov. 2022, 8, 70. https://doi.org/10.1038/s41420-022-00863-3

21. Li, Y.D.; Chen, X.; Dong, X.D.; Ding, Z.; Chen, Z.S.; Fang, S. Methyl-cantharidimide suppresses cyclin-dependent kinase 1 (CDK1) and induces oxeiptosis in liver cancer. J Exp Clin Cancer Res. 2026, 45, 90. https://doi.org/10.1186/s13046-026-03669-8
