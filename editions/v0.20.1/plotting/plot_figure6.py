"""Plot all reference-set pairs from the tables in Zenodo 10.5281/zenodo.24421516.

Usage: python plot_figure6.py <unpacked-package>/tables Figure_6.tif
Dependencies: Matplotlib, Pillow and Arial.
"""
import csv
import sys
from pathlib import Path

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from PIL import Image


def draw(tables, outfile):
    with (tables / "S11_reference_pair_all108.csv").open(encoding="utf-8-sig") as f:
        rows = list(csv.DictReader(f))
    plt.rcParams.update({"font.family":"Arial", "font.size":9, "axes.spines.top":False,
                         "axes.spines.right":False, "axes.linewidth":0.7})
    fig, axes = plt.subplots(1, 3, figsize=(6.5, 3), sharey=True)
    fig.subplots_adjust(left=.12,right=.982,bottom=.19,top=.80,wspace=.23)
    for i, (ax, cohort) in enumerate(zip(axes, ("TCGA-COADREAD","GSE39582","GSE17538"))):
        group = [row for row in rows if row["cohort"] == cohort]
        x = [float(row["jaccard_measured_clean_audit"]) for row in group]
        y = [float(row["rho_clean_audit"]) for row in group]
        assert len(group) == 36 and len({(r["set_a"],r["set_b"]) for r in group}) == 36
        ax.set(xlim=(-.01,.2),ylim=(-.40,1.05),xticks=[0,.1,.2],
               yticks=[-.4,-.2,0,.2,.4,.6,.8,1],xlabel="Jaccard overlap")
        assert all(ax.get_xlim()[0]<v<ax.get_xlim()[1] for v in x)
        assert all(ax.get_ylim()[0]<v<ax.get_ylim()[1] for v in y)
        ax.scatter(x,y,s=18,color="#2166A5",alpha=.70,edgecolors="white",linewidth=.3)
        ax.axvline(.05,c="#707070",ls="--",lw=.7)
        ax.axhline(.50,c="#707070",ls="--",lw=.7)
        ax.set_title(f"{cohort}\nn = {group[0]['n']}",fontsize=9)
        ax.text(-.12,1.20,chr(97+i),transform=ax.transAxes,fontweight="bold",fontsize=11)
    axes[0].set_ylabel("Reference-score Spearman ρ")
    outfile.parent.mkdir(parents=True,exist_ok=True)
    import io
    buf=io.BytesIO()
    fig.savefig(buf,format="png",dpi=600,facecolor="white")
    plt.close(fig)
    buf.seek(0)
    with Image.open(buf) as im:
        im.convert("RGB").save(outfile,dpi=(600,600),compression="tiff_lzw")


if __name__ == "__main__":
    draw(Path(sys.argv[1]),Path(sys.argv[2]))
