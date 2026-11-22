# Reproducing the updated figures

## Figure 1

Source image, spatial counts and histological annotations are from Alberto Valdeolivas, *Profiling the Heterogeneity of Colorectal Cancer Consensus Molecular Subtypes using Spatial Transcriptomics: datasets*, [source record](https://doi.org/10.5281/zenodo.7760264), [CC BY 4.0](https://creativecommons.org/licenses/by/4.0/). The displayed section is SN048_A121573_Rep1. The image has been cropped, resized and lightened beneath the overlays. Data-derived spots and calibrated 500 µm bars have been added. Panels a–c share the same field and physical calibration. Panel b shows the two groups entering the compartment contrast. Panel c shows section-standardized OX4 at all 2,049 eligible spots. Colours saturate at ±2 SD. The table retains the full values. Panels d–f reproduce the archived quantitative panels with revised lettering.

Run from the root of the unpacked 0.20.1 archive. Supply the original named section ZIP from the source data record. The quantitative source image is bundled in figure_sources/. Use Python with NumPy, Pillow and Arial regular/bold. Supply the directory containing arial.ttf and arialbd.ttf on systems outside Windows.

```sh
python plotting/figure1.py SN048_A121573_Rep1.zip . rendered/figure1 /path/to/Arial
```

The script writes a TIFF, the matching per-spot table and a provenance record. The record includes the SHA-256 of the source ZIP and image, common crop, Space Ranger scale factors, lengths and positions of all three scale bars, and the within-section Cliff’s delta. The score uses KEAP1, PGAM5, AIFM1 and OTUD1 with equal weights. The section was selected by ID before consulting score values.

## Figure 6

The corrected plot spans −0.40 to 1.05 on the correlation axis. Four negative correlations from TCGA now fall within the displayed range. The input is the same S11 table of 108 pairs included in this release. Use Python with Matplotlib, Pillow and Arial installed.

```sh
python plotting/plot_figure6.py tables rendered/Figure_6.tif
```
