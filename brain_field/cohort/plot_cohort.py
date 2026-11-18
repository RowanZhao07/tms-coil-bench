"""Plot the two drive comparisons within each subject; keep the pilot separate."""
import argparse
from pathlib import Path
import sys

import matplotlib as mpl
mpl.use("Agg")
import matplotlib.pyplot as plt
from matplotlib.lines import Line2D
from matplotlib.ticker import MaxNLocator
import numpy as np
import pandas as pd

HERE = Path(__file__).resolve().parent
IDS = ["sub01", "sub09", "sub10", "sub12", "sub15", "FieldTrip_Subject01"]
Y = np.array([5.5, 4.5, 3.5, 2.5, 1.5, 0.0])
COILS = {
    "h1": ("H1", "#147665", "o"), "double_cone": ("Double-cone", "#aa5730", "D"),
    "figure8": ("Figure-8", "#446987", "s"),
}


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--audit-tools", type=Path)
    parser.add_argument("--audit-json", type=Path)
    args = parser.parse_args()
    out = HERE / "results"
    data = pd.read_csv(out / "per_subject.csv")
    if data.duplicated(["subject", "coil"]).any() or len(data) != len(IDS) * len(COILS):
        raise ValueError("Expected exactly one measurement for each subject and coil")
    mpl.rcParams.update({"pdf.fonttype": 42, "svg.fonttype": "none", "font.family": "sans-serif",
        "font.sans-serif": ["Arial", "DejaVu Sans"], "font.size": 8, "axes.labelsize": 8,
        "axes.spines.top": False, "axes.spines.right": False, "axes.linewidth": .65,
        "xtick.labelsize": 7.5, "ytick.labelsize": 8, "legend.frameon": False})
    fig, axes = plt.subplots(1, 2, figsize=(180/25.4, 81/25.4))
    fig.subplots_adjust(left=.16, right=.98, bottom=.26, top=.88, wspace=.24)
    panels = [("fixed_target_V_m", "At 100 A/µs"), ("matched_target_V_m", "At cortical P99 = 30 V/m")]
    for index, (ax, (column, title)) in enumerate(zip(axes, panels)):
        scores = data.pivot(index="subject", columns="coil", values=column).loc[IDS]
        offsets = {"h1": 0, "double_cone": .13, "figure8": -.13}
        for y, (_, row) in zip(Y, scores.iterrows()):
            ordered = row.sort_values()
            ax.plot(ordered, [y + offsets[k] for k in ordered.index],
                    color="#bac3c7", lw=1.0, zorder=1)
        for coil, (_, color, marker) in COILS.items():
            ax.scatter(scores[coil], Y + offsets[coil], c=color, marker=marker, s=29,
                       linewidths=.35, edgecolors="white", zorder=2)
        ax.axhline(.73, color="#a3abb1", ls=(0, (3, 3)), lw=.65)
        ax.set_title(title, loc="left", x=.075, pad=9, fontsize=9)
        ax.text(0, 1, "ab"[index], va="baseline", fontweight="bold", fontsize=10,
                transform=ax.transAxes + mpl.transforms.ScaledTranslation(0, 9/72, fig.dpi_scale_trans))
        ax.set(ylim=(-.55,6), xlim=(0, scores.to_numpy().max()*1.10), yticks=Y,
               yticklabels=IDS[:-1]+["FieldTrip pilot"] if index==0 else [""]*len(IDS),
               xlabel="Left thalamus mean |E| (V/m)")
        ax.tick_params(axis="y", length=0)
        ax.xaxis.set_major_locator(MaxNLocator(nbins=5))
    handles = [Line2D([], [], linestyle="none", marker=m, color=c, markersize=5, label=label)
               for label,c,m in COILS.values()]
    fig.legend(handles=handles, loc="lower center", bbox_to_anchor=(.56,.015), ncol=3, columnspacing=2)
    fig.canvas.draw()
    if args.audit_tools:
        sys.path.insert(0,str(args.audit_tools))
        from audit_panel_alignment import require_matplotlib_panel_alignment
        require_matplotlib_panel_alignment(fig, axes=list(axes), panel_ids=["a", "b"],
                                           strict=True, json_out=args.audit_json)
    fig.savefig(out/"paired_fields.pdf", dpi=600)
    fig.savefig(out/"paired_fields.svg")
    fig.savefig(out/"paired_fields.png", dpi=600)
    fig.savefig(out/"paired_fields.tiff", dpi=600, pil_kwargs={"compression": "tiff_lzw"})
    plt.close(fig)


if __name__ == "__main__":
    main()
