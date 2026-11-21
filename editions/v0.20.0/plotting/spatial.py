"""Render the two-cohort spatial panels from the verified section and patient tables."""

from pathlib import Path

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
from PIL import Image

ROOT = Path(__file__).resolve().parents[1]
OUT = ROOT / 'rendered_spatial'
FIG = OUT / 'figures'
REVIEW = OUT / 'sources'
FIRST, SECOND = ('Valdeolivas', 'GSE294385')
TABLES = ROOT / 'tables'
COHORTS = ["Valdeolivas", "GSE294385"]
COLORS = ["#216889", "#9A4C18"]
PLOTS = []
plt.rcParams.update({"font.family": "Arial", "font.size": 9,
                     "axes.spines.right": False, "axes.spines.top": False,
                     "axes.grid": True, "axes.axisbelow": True, "grid.alpha": .23,
                     "pdf.fonttype": 42, "svg.fonttype": "none"})


def letter(ax, label):
    ax.text(-.19, 1.06, label, transform=ax.transAxes, fontweight="bold", fontsize=12)


def save(fig, stem):
    fig.savefig(REVIEW / f"{stem}.svg", bbox_inches="tight")
    png = REVIEW / f"{stem}_proof.png"
    fig.savefig(png, dpi=600, bbox_inches="tight", facecolor="white")
    with Image.open(png) as opened:
        opened.convert("RGB").save(FIG / f"{stem}.tif", compression="tiff_lzw", dpi=(600, 600))
    plt.close(fig)


def main():
    FIG.mkdir(parents=True, exist_ok=True)
    REVIEW.mkdir(parents=True, exist_ok=True)
    fig, ax = plt.subplots(1, 3, figsize=(9.65, 3.05), layout="constrained")
    bgfig, bgax = plt.subplots(1, 2, figsize=(6.7, 2.87), layout="constrained", sharex=True, sharey=True)
    for idx, (root, label, color) in enumerate(zip([FIRST, SECOND], COHORTS, COLORS)):
        people = pd.read_csv(TABLES / ("S32_patient_effects.csv" if idx==0 else "S38_patient_effects.csv"), dtype={"patient": str})
        people = people[people.construct.eq("OX4")].sort_values("patient")
        infer = pd.read_csv(TABLES / ("S32_patient_inference.csv" if idx==0 else "S38_patient_inference.csv")).set_index("construct").loc["OX4"]
        n = len(people)
        assert n == infer.n_patients
        assert np.isclose(people.cliff_delta.mean(), infer.mean_delta)
        x = idx + np.linspace(-.13, .13, n)
        ax[0].scatter(x, people.cliff_delta, ec="white", linewidth=.45, color=color, s=28, zorder=3)
        ax[0].errorbar(idx+.24, infer.mean_delta,
                       yerr=[[infer.mean_delta-infer.bootstrap_lo], [infer.bootstrap_hi-infer.mean_delta]],
                       fmt="D", ms=4, capsize=3, color="#222222", zorder=4)
        for row in people.itertuples():
            ax[1].plot([0, 1], [row.cliff_delta, row.depth_adjusted_delta], marker="o",
                       ms=3, lw=.85, color=color, alpha=.7, label=label if row.Index == people.index[0] else None)
            PLOTS.append(dict(panel="a/b", cohort=label, unit=row.patient,
                              cliff_delta=row.cliff_delta, depth_adjusted_delta=row.depth_adjusted_delta))
        source = "target_section_results.csv" if idx == 0 else "section_effects.csv"
        section = pd.read_csv(TABLES / ("S33_section_results.csv" if idx==0 else "S39_section_effects.csv"))
        paired = section.pivot(index="sample", columns="construct", values="median_split_rho")
        assert len(paired) == [12, 8][idx] and (paired.OX4 < paired.stromal).all()
        ax[2].scatter(paired.stromal, paired.OX4, color=color, label=label, alpha=.85, s=26,
                      marker=["o", "^"][idx], edgecolors="white", linewidths=.45, zorder=3)
        for key, row in paired.iterrows():
            PLOTS.append(dict(panel="c", cohort=label, unit=key, ox4_rho=row.OX4, stromal_rho=row.stromal))
        matches = pd.read_csv(TABLES / ("S34_matched_positions.csv" if idx==0 else "S41_background_positions.csv"))
        col = "abs_delta_percentile" if idx == 0 else "absolute_delta_percentile"
        ba = bgax[idx]
        ba.scatter(matches[col], matches.split_rho_percentile, color=color, s=26, edgecolors="white",
                   linewidths=.5, marker=["o", "^"][idx], zorder=3)
        ba.axhline(50, c="#777777", ls="--", lw=.75)
        ba.axvline(50, c="#777777", ls="--", lw=.75)
        ba.set(title=f"{label} ({len(matches)} sections)", xlabel="Absolute compartment effect percentile",
               xlim=(0, 102), ylim=(0, 102), xticks=[0, 25, 50, 75, 100], yticks=[0, 25, 50, 75, 100])
        letter(ba, "ab"[idx])
        for row in matches.to_dict("records"):
            PLOTS.append(dict(panel="S1", cohort=label, unit=row["sample"],
                              delta_percentile=row[col], rho_percentile=row["split_rho_percentile"]))
    ax[0].set(xlim=(-.34,1.47), ylim=(-.04,.74), xticks=[0,1],
              xticklabels=["Valdeolivas\n6 patients", "GSE294385\n5 patients"], ylabel="Cliff’s delta")
    ax[0].axhline(0, color="#777777", linewidth=.7)
    ax[1].set(xlim=(-.15,1.15), ylim=(-.12,.74), xticks=[0,1],
              xticklabels=["Full score", "Depth residual"], ylabel="Cliff’s delta")
    ax[1].axhline(0, color="#777777", linewidth=.7)
    ax[2].set(xlim=(-.025,1.025), ylim=(-.025,1.025), xlabel="Stromal reference split ρ",
              ylabel="OX4 split ρ", xticks=[0,.25,.5,.75,1], yticks=[0,.25,.5,.75,1])
    ax[2].plot([0,1], [0,1], c="#777777", ls="--", lw=.85)
    ax[2].legend(loc="upper left", frameon=False, fontsize=7)
    for panel, axes in zip("abc", ax):
        letter(axes, panel)
    bgax[0].set_ylabel("Count-split correlation percentile")
    pd.DataFrame(PLOTS).to_csv(REVIEW / "spatial_plot_values.csv", index=False)
    save(fig, "Figure1")
    save(bgfig, "FigureS1")


if __name__ == "__main__":
    main()
