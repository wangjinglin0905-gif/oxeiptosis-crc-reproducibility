# Text S3. Exploratory spatial measurement extension

This text describes GSE280315 and GSE226997; the subsequent morphology dataset is described separately in Text S4. This extension was added after the original bulk study. Each computational stage had a local analysis contract before its corresponding analysis; the stages were sequential and informed by earlier findings. No external preregistration, independent pathology rating or new wet experiment was performed.

## 1. Analysis units and provenance

GSE280315: five filtered 8-μm HDF5 matrices, paired position tables and overview H&E images. Section labels P1CRC/P2CRC/P5CRC/P3NAT/P5NAT indicate four donors and one tumour/normal-adjacent pair. The five sections are not five independent patients. Bins were summed into origin-aligned 32/64/128/256-μm squares. A second grid shifted by half the region width in both axes was used for a descriptive boundary sensitivity. Correlations between projected grids/scales used shared underlying bins and were not independent validation.

GSE226997: GSM7089855–58, four patients. STimage-1K4M files were fixed at Hugging Face commit b2073c6e8bb7e6d8dfd15c0db6134ba4f2a8a84f. Original count provenance was verified using the complete filtered_feature_bc_matrix.h5 members extracted from each original GEO tar.gz prefix. The four 64-MiB prefixes were retrieved as 2-MiB HTTP ranges; status 206, exact Content-Range and byte count were required. SHA256, URL and acquisition times were recorded. Full compressed archives were not obtained; their terminal gzip CRC and FASTQ processing were not verified.

HDF5 matrices contained 17,944 features in GSM7089855/56/58 and 19,144 in GSM7089857. All 294,437,824 counts, feature positions and barcode orders agreed with the derivatives. Duplicate symbols received unique suffixes without aggregation. The mapping tables identify each feature by its original Ensembl ID. OX4 and the predefined controls were unaffected by the duplicate names. Depth denominators matched original filtered-HDF5 column totals exactly. Thus, earlier technical estimates were retained without recomputation.

The HDF5 files contained a Visium Human Transcriptome Probe Set target group. Software attributes were spaceranger-1.3.1 for three samples and spaceranger-1.3.0 for GSM7089857. GEO free text mentions a different software version and inconsistent annotation-version labels. The audit records these differences and does not assign an unverified exact annotation release. All four HDF5 genome fields read GRCh38. Counts and the denominator refer to represented features in the deposited matrices, not an unrestricted transcriptome.

## 2. Eligibility and construction

GSE280315 regions required represented bins occupying at least 50% of the grid and total counts >=500. GSE226997 spots required total counts >=500 and exact expression/coordinate key correspondence. No four-positive selection, imputation, denoising or spatial smoothing was used. Absent feature columns were not converted into zero counts.

For unit i and gene g, x_ig = log(1 + 10,000 c_ig/t_i). The section-specific mean and sample SD of x_ig among eligible units were used for standardization. Spatial OX4 was the mean of the four standardized values. The bulk construct used its original cohort calibration and composite scaling; no numeric cutoff or patient score was transferred from bulk to space.

Constructs: OX4 = KEAP1/PGAM5/AIFM1/OTUD1; stromal = COL1A1/COL1A2/COL3A1/DCN; epithelial = EPCAM/KRT8/KRT18/KRT19. Controls are investigator-defined marker averages. All members were required. The five supplied GSE280315 feature/name arrays lacked exact KRT18; three original GSE226997 matrices lacked KRT19. These constructs were marked unavailable, not shortened. Gene alias/probe remediation was not performed for GSE280315.

## 3. Technical splitting

For each unit, each construct-member count was split by Binomial(c_ig, 0.5). Remaining represented counts were split by Binomial(t_i − sum_g c_ig, 0.5), preserving member/total consistency. The second half was the complement. Each half used its own total for CP10K normalization; both used the full-data gene means and SDs. This is conditional count splitting with a shared reference, not two independently measured biological replicates.

The primary statistic was Spearman correlation between half-scores across eligible units in a section. We report medians and P10/P90 of technical repetitions. These are not patient confidence intervals. Depth-conditioned sensitivity uses rank scores regressed on an intercept and both half-depth ranks; the correlation of residuals is reported. It does not identify a causal technical contribution, because depth can also track tissue structure.

The GSE280315 pilot used 30 splits for each of 5 sections × 4 scales =600 repetitions, seed 20261117 with section accession, scale and repetition included in SeedSequence. Its composite feasibility gate required >=50 eligible regions, >=80% joint member detection, median rho>=0.7, P10>=0.6, and absolute full-score/log-depth rho<=0.5. None of the 20 combinations passed all conditions. The thresholds were local decision rules, not an assay standard.

The matched GSE280315 stage used seed 20261118 and 10 repetitions per construct. The GSE226997 stage used SeedSequence([20261119, numeric_GSM, construct_index, repetition]), with construct order OX4, epithelial, stromal, and 30 repetitions per complete construct. The numeric seed is not a date of execution. GSE226997 used the median/P10 stability thresholds and separately reported detection and depth correlations; it did not apply an identical physical-region gate or run a second matched-set experiment.

## 4. Spatial matched backgrounds

Each OX4 member was matched within section to genes with mean CP10K between 0.5 and 2 times its mean and detection-rate absolute difference <=0.05. OX4 members, predefined control members, duplicate-symbol genes and MT-/RPL/RPS-prefixed genes were excluded. Candidates were ranked by squared log2 mean ratio plus squared detection-rate difference divided by 0.05 squared. Up to the nearest 30 were retained; a pool of fewer than 10 would stop that section rather than relax constraints. All sections met the pool criterion.

Each section used 100 unique four-gene sets, without duplicate members within a set. Genes could recur across sets. The reported percentile is the fraction of the 100 matched medians strictly below OX4. The stage criterion for an additional deficit was at least 3 of 5 sections below the matched P10. No section met that condition, including the depth-conditioned comparison. These ranks are descriptive, and overlapping matched sets are not independent experiments.

This stage completed 5 × (100 matched + OX4 + stromal) ×10 =5,100 splits. Matching controlled marginal expression and detection, not variance, pairwise co-expression, spatial autocorrelation or cellular source. Its interpretation is separate from the clinical survival-background comparison in the original manuscript.

## 5. Pathology scope

H&E/coordinate overlays were inspected for gross flips or shifts. All four GSE226997 coordinate centres fell inside their images. No independently annotated histological endpoint was assessed in these two datasets. Park et al. describe compartments defined jointly by expression clustering and H&E morphology; those labels cannot be treated as expression-independent morphology truth.

The GSE226997 assessment was limited to coordinate alignment. No independently rated tumour–stroma regions were available for this dataset.

## 6. Verification and deliverables

Prior pilot QA compared selected region counts with original sparse matrices and recreated a fixed split with absolute correlation difference 5.55e-17. The matched stage recreated a control's aggregation through an independent bincount implementation, with zero correlation discrepancy. GSE226997 QA recreated a fixed split using explicit formulas and rank-correlation computation, again differing by 5.55e-17. These are numerical checks, not independent human review or independent experiments.

Table S28 has all 20 scale summaries; S29 has all 510 matched-stage construct summaries; S30 has all nine calculable GSE226997 construct summaries; S31 has the four original-file audits. Full repetitions, candidate pools, memberships and region/spot tables remain in the project stage directories listed in spatial_artifact_index.csv. Existing bulk Figures 1–4 and Tables 1–3 were retained. Main Tables 4–5 were generated directly from the frozen spatial summary CSVs. The subsequent morphology analysis is displayed in main Figure 5 and detailed in Text S4.

Table S27 maps the current bibliography to its uses in the main text.

Code S1 provides the scripts, source accessions and derived tables. The portable verifier reads frozen morphology tables. Reanalysis from the original expression matrices requires configuring the historical scripts for local paths.
