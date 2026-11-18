"""Plot paired cortical exposure and threshold sensitivity for five native MRI heads."""
import argparse
import string
import sys
from pathlib import Path

import matplotlib as mpl
mpl.use("Agg")
import matplotlib.pyplot as plt
from matplotlib.lines import Line2D
from matplotlib.ticker import MaxNLocator
import numpy as np
import pandas as pd


RESULTS = Path(__file__).resolve().parent / "results"
SUBJECTS = ["sub01", "sub09", "sub10", "sub12", "sub15"]
COILS = {
    "h1": ("H1", "#147665", "o", "-"),
    "double_cone": ("Double-cone", "#aa5730", "D", "--"),
    "figure8": ("Figure-8", "#446987", "s", "-."),
}
TITLES = {"cortex_P99_30": "Cortical P99 = 30 V/m", "target_mean_1": "Thalamic mean = 1 V/m"}


def matrix(data, condition, measurement):
    by_subject = data.loc[data.condition == condition].pivot(
        index="subject", columns="coil", values=measurement
    )
    if set(by_subject.index) != set(SUBJECTS) or set(by_subject.columns) != set(COILS):
        raise ValueError("Expected all five subjects and all three coils")
    values = by_subject.loc[SUBJECTS, list(COILS)]
    if not np.isfinite(values.to_numpy()).all() or (values < 0).any().any():
        raise ValueError(f"Invalid {measurement} in {condition}")
    return values


def heading(ax, index, condition):
    ax.set_title(TITLES[condition], loc="left", x=.095, pad=9, fontsize=8.5)
    ax.annotate(string.ascii_lowercase[index], (0, 1), xycoords="axes fraction",
                xytext=(0, 9), textcoords="offset points", va="baseline",
                fontsize=10, fontweight="bold")


def legend(fig, lines=False):
    handles = [Line2D([], [], color=color, marker=marker, markersize=4.5,
                      lw=1.35, ls=style if lines else "none", label=label)
               for label, color, marker, style in COILS.values()]
    fig.legend(handles=handles, loc="lower center", ncol=3, bbox_to_anchor=(.53,.012),
               columnspacing=2.2, handletextpad=.45)


def paired_plot(conditions, thresholds):
    volumes = thresholds.loc[thresholds.threshold_V_m == 15]
    matched_cortex = matrix(volumes, "cortex_P99_30", "cortex_cm3")
    matched_target = matrix(volumes, "target_mean_1", "cortex_cm3")
    cortical_p99 = matrix(conditions, "target_mean_1", "cortex_P99_V_m")
    maximum_volume = max(matched_cortex.to_numpy().max(), matched_target.to_numpy().max())
    volume_label = "Cortex ≥15 V/m (cm³)"
    panels = [
        ("cortex_P99_30", matched_cortex, volume_label, maximum_volume),
        ("target_mean_1", cortical_p99, "Cortical P99 (V/m)", cortical_p99.to_numpy().max()),
        ("target_mean_1", matched_target, volume_label, maximum_volume),
    ]
    fig, axes = plt.subplots(1, 3, figsize=(180/25.4, 72/25.4))
    fig.subplots_adjust(left=.085, right=.986, top=.83, bottom=.27, wspace=.19)
    lanes = np.arange(len(SUBJECTS), 0, -1)
    offsets = {"h1": .17, "double_cone": 0, "figure8": -.17}
    for index, (ax, (condition, values, label, maximum)) in enumerate(zip(axes, panels)):
        for y, (_, row) in zip(lanes, values.iterrows()):
            ax.hlines(y, row.min(), row.max(), lw=.85, color="#c6cdd0", zorder=1)
        for coil, (_, color, marker, _) in COILS.items():
            ax.scatter(values[coil], lanes + offsets[coil], s=21, c=color, marker=marker,
                       edgecolor="white", linewidth=.38, zorder=3, clip_on=False)
        if condition == "target_mean_1" and label == volume_label:
            # Several H1 volumes resolve to only a few tetrahedra. Label zero and small positives.
            for y, volume in zip(lanes + offsets["h1"], values["h1"]):
                ax.annotate(format(volume, ".2g"), (volume, y), xytext=(4.5, 0),
                            textcoords="offset points", va="center", color=COILS["h1"][1],
                            fontsize=7)
        # The left margin leaves the full marker visible for zero or near-zero volumes.
        ax.set(xlim=(-.025*maximum, maximum*1.06), ylim=(.53,5.47),
               yticks=lanes, yticklabels=SUBJECTS if index == 0 else [""]*len(SUBJECTS),
               xlabel=label)
        ax.tick_params(axis="y", length=0, pad=3)
        ax.xaxis.set_major_locator(MaxNLocator(nbins=4, steps=[1, 2, 4, 5, 10], integer=True))
        heading(ax, index, condition)
    legend(fig)
    return fig, axes


def threshold_plot(curves):
    fig, axes = plt.subplots(1, 2, figsize=(180/25.4, 76/25.4))
    fig.subplots_adjust(left=.085, right=.985, top=.845, bottom=.26, wspace=.15)
    displayed = curves.loc[curves.threshold_V_m.between(1, 60)]
    if not np.isfinite(displayed.cortex_percent).all() or not displayed.cortex_percent.between(0,100).all():
        raise ValueError("Cortical percentages must be finite and between 0 and 100")
    ymax = displayed.cortex_percent.max()
    for index, (ax, condition) in enumerate(zip(axes, TITLES)):
        series = displayed.loc[displayed.condition == condition]
        if set(series.coil) != set(COILS):
            raise ValueError("Missing coil in exposure curves")
        for coil, (_, color, marker, style) in COILS.items():
            values = series.loc[series.coil == coil].pivot(
                index="threshold_V_m", columns="subject", values="cortex_percent"
            ).sort_index()
            if set(values.columns) != set(SUBJECTS) or values.isna().any().any():
                raise ValueError("Incomplete curve for a subject")
            if (values.diff().iloc[1:] > 1e-8).any().any():
                raise ValueError("Cumulative exposure must decrease with threshold")
            median = values.median(axis=1)
            ax.fill_between(values.index, values.min(axis=1), values.max(axis=1),
                            facecolor=color, alpha=.12, edgecolor="none", zorder=1)
            ax.plot(values.index, median, color=color, ls=style, lw=1.35, zorder=3)
            ax.scatter([15], [median.loc[15]], marker=marker, c=color, s=16,
                       edgecolor="white", linewidth=.3, zorder=4)
        ax.set_yscale("function", functions=(lambda v: np.sqrt(np.maximum(v, 0)), np.square))
        ticks = [v for v in [0, 1, 5, 15, 30, 60, 100] if v <= ymax*1.07]
        ax.set(ylim=(0, min(100,ymax*1.07)), xlim=(1,60), xticks=[1,7.5,15,30,45,60],
               xticklabels=["1", "7.5", "15", "30", "45", "60"], yticks=ticks,
               xlabel="Field threshold (V/m)")
        if index == 0:
            ax.set_ylabel("Cortical volume (%)\nSquare-root scale")
        else:
            ax.tick_params(axis="y", labelleft=False)
        ax.grid(axis="y", linewidth=.4, color="#dce0e1")
        ax.set_axisbelow(True)
        heading(ax, index, condition)
    legend(fig, lines=True)
    return fig, axes


def export(fig, axes, name, args):
    fig.canvas.draw()
    if args.audit_tools:
        sys.path.insert(0, str(args.audit_tools))
        from audit_panel_alignment import require_matplotlib_panel_alignment
        args.audit_dir.mkdir(parents=True, exist_ok=True)
        require_matplotlib_panel_alignment(fig, axes=list(axes),
            panel_ids=list(string.ascii_lowercase[:len(axes)]), strict=True,
            json_out=args.audit_dir / f"{name}.alignment.json")
    output = RESULTS / name
    fig.savefig(output.with_suffix(".svg"))
    fig.savefig(output.with_suffix(".pdf"), dpi=600)
    fig.savefig(output.with_suffix(".png"), dpi=600)
    fig.savefig(output.with_suffix(".tiff"), dpi=600, pil_kwargs={"compression": "tiff_lzw"})
    plt.close(fig)


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--audit-tools", type=Path)
    parser.add_argument("--audit-dir", type=Path, help="Directory for the optional alignment reports")
    args = parser.parse_args()
    if args.audit_tools and not args.audit_dir:
        parser.error("--audit-dir is required when checking panel alignment")
    mpl.rcParams.update({"pdf.fonttype": 42, "svg.fonttype": "none", "font.family": "sans-serif",
                         "font.sans-serif": ["Arial", "DejaVu Sans"], "font.size": 8,
                         "axes.spines.top": False, "axes.spines.right": False,
                         "axes.linewidth": .65, "axes.labelsize": 8, "xtick.labelsize": 7.5,
                         "ytick.labelsize": 7.5, "legend.frameon": False})
    conditions = pd.read_csv(RESULTS / "conditions.csv")
    thresholds = pd.read_csv(RESULTS / "thresholds.csv")
    curves = pd.read_csv(RESULTS / "curves.csv")
    if set(conditions.condition) != set(TITLES) or len(conditions) != 30:
        raise ValueError("Expected 30 subject-by-coil-by-condition records")
    export(*paired_plot(conditions, thresholds), "focality_paired", args)
    export(*threshold_plot(curves), "focality_thresholds", args)


if __name__ == "__main__":
    main()
