"""Redraw the MRI, sensor selection and calibration figures from saved experiments."""
from argparse import ArgumentParser
from pathlib import Path
import string
import sys

import matplotlib as mpl
mpl.use("Agg")
import matplotlib.pyplot as plt
from matplotlib.collections import LineCollection
from matplotlib.colors import LinearSegmentedColormap
from matplotlib.lines import Line2D
from matplotlib.ticker import FixedLocator, ScalarFormatter
import numpy as np
import pandas as pd
from scipy.io import loadmat

HERE = Path(__file__).resolve().parent
ROOT = HERE.parents[1]
MRI = ROOT / "placement_pilot/public_benchmark/matched_figure8"
LAB = ROOT / "external_measurements"
INK, TEAL, RUST, BLUE, GREY = "#25343d", "#087863", "#ae5529", "#37679b", "#829098"
CMAP = LinearSegmentedColormap.from_list("signed_field", ["#174972", "#619dbf", "#f8f7f3", "#cf8770", "#8e372b"])
CMAP.set_bad("#e6e9e8")
STYLE = ("Fixed6", "Dopt6", "Weighted6")
SITE_NAMES = ("Symmetric", "D-optimal", "Residual-weighted")
SITE_COLORS = (INK, RUST, TEAL)
SITE_MARKS = ("s", "^", "D")

mpl.rcParams.update({
    "font.family": "sans-serif", "font.sans-serif": ["Arial", "DejaVu Sans"],
    "font.size": 7.5, "axes.titlesize": 8.5, "axes.labelsize": 7.5,
    "axes.titlelocation": "left", "axes.labelcolor": INK, "text.color": INK,
    "xtick.color": INK, "ytick.color": INK, "axes.edgecolor": GREY,
    "axes.spines.top": False, "axes.spines.right": False, "axes.linewidth": .65,
    "xtick.labelsize": 7, "ytick.labelsize": 7, "xtick.major.size": 2.5, "ytick.major.size": 2.5,
    "xtick.major.width": .6, "ytick.major.width": .6, "legend.frameon": False,
    "legend.fontsize": 7, "legend.handlelength": 1.7, "legend.columnspacing": 1.2,
    "svg.fonttype": "none", "pdf.fonttype": 42, "ps.fonttype": 42, "savefig.facecolor": "white", "lines.linewidth": 1.3,
})


def read(folder, name):
    return pd.read_csv(folder / name)


def frame(n, height_mm, *, left=.10, right=.97, bottom=.26, top=.85, gap=.49):
    fig, axes = plt.subplots(1, n, figsize=(180 / 25.4, height_mm / 25.4), squeeze=False)
    fig.subplots_adjust(left=left, right=right, bottom=bottom, top=top, wspace=gap)
    return fig, list(axes[0])


def title(ax, letter, text):
    ax.set_title(text, loc="left", pad=8, fontweight="normal")
    ax.text(-.105, 1.055, letter, transform=ax.transAxes, fontsize=9.5, fontweight="bold",
            va="bottom", ha="right")


def key(fig, handles, y=.01, ncol=None):
    fig.legend(handles=handles, loc="lower center", bbox_to_anchor=(.5, y), ncol=ncol or len(handles))


def rule(ax, axis="y"):
    ax.set_axisbelow(True)
    ax.grid(axis=axis, color="#dce1df", linewidth=.5)


def mark(color, marker, label, *, filled=True):
    return Line2D([], [], linestyle="none", marker=marker, markersize=4.7, markeredgewidth=.9,
                  markeredgecolor=color, markerfacecolor=color if filled else "white", label=label)


def field(ax, values, mask, limit=180):
    img = ax.imshow(np.ma.array(values, mask=~mask), cmap=CMAP, vmin=-limit, vmax=limit,
                    extent=(.5, 160.5, 160.5, .5), interpolation="nearest", aspect="equal")
    ax.set(xlim=(14.5, 150.5), ylim=(151.5, 15.5), xticks=[30, 80, 130], yticks=[30, 80, 130],
           xlabel="Pixel column")
    return img


def maps(d):
    valid = d["mask"].astype(bool)
    assert valid.sum() == 10622
    pred = d["prediction_T"] * 1e6
    measured = d["measured_T"] * 1e6
    fig, axes = frame(3, 73, left=.080, right=.965, bottom=.275, top=.86, gap=.26)
    for i, (ax, data, header) in enumerate(zip(axes, [pred[..., 0], pred[..., 1], measured],
                                               ["Original pose", "Corrected pose", "Phantom MRI"])):
        im = field(ax, data, valid)
        title(ax, string.ascii_lowercase[i], header)
        if i == 0:
            ax.set_ylabel("Pixel row")
        else:
            ax.tick_params(labelleft=False)
    bar_ax = fig.add_axes([.355, .114, .34, .030])
    colorbar = fig.colorbar(im, cax=bar_ax, orientation="horizontal", extend="both",
                            ticks=[-180, -90, 0, 90, 180])
    colorbar.set_label("Signed field along scanner z (µT)", labelpad=1)
    scores = read(MRI, "placement_version_comparison.csv")
    for i, version in enumerate(("v1.0", "v1.1")):
        target = scores.query("model == 'BSM-819' and reference == 'phantom_MRI' and condition == @version").iloc[0]
        calc = 100 * np.linalg.norm((pred[..., i]-measured)[valid])/np.linalg.norm(measured[valid])
        assert np.isclose(calc, target.relative_L2_percent)
    return fig, axes


def projection(d):
    fig, axes = frame(2, 87, left=.095, right=.97, bottom=.25, top=.88, gap=.39)
    ax, profile = axes
    x = 85.033 - (np.arange(160) + .5)
    z = -85.401 + 160.5 - (np.arange(160) + 1)
    ax.contourf(x, z, d["mask"], levels=[.5, 1.5], colors=["#dde4e1"])
    ax.contour(x, z, d["mask"], levels=[.5], colors=[GREY], linewidths=.55)
    for i, color in enumerate([RUST, TEAL]):
        wires = d["registered_wires_mm"][0, i]
        segments = np.stack([wires[:, (0, 2)], wires[:, (3, 5)]], axis=1)
        ax.add_collection(LineCollection(segments, colors=color, linewidths=.43))
        ax.plot(*d["poses"][[0, 2], 3, i], marker="o", mfc="white", mec=color, ms=3.1)
    ax.set(xlim=(-105, 105), ylim=(-90, 110), xlabel="Scanner x (mm)", ylabel="Scanner z (mm)",
           xticks=[-100, -50, 0, 50, 100], yticks=[-75, -25, 25, 75, 100])
    ax.set_aspect("equal", adjustable="datalim")
    title(ax, "a", "Coil and slice in scanner space")
    row = read(MRI, "profile_row50.csv")
    series = [("v1_0_uT", RUST, "--", "Original"), ("v1_1_uT", TEAL, "-", "Corrected")]
    for column, color, ls, label in series:
        profile.plot(row.column_one_based, row[column], color=color, ls=ls, lw=1.3, label=label)
    profile.plot(row.column_one_based, row.MRI_uT, marker=".", ls="none", color=INK, ms=2.3, label="Phantom MRI")
    profile.axhline(0, color=GREY, linewidth=.5, zorder=0)
    profile.set(xlabel="Pixel column in row 50", ylabel="Signed field (µT)", xticks=[40, 60, 80, 100, 120])
    title(profile, "b", "One fixed line through the field")
    key(fig, profile.get_legend_handles_labels()[0], .018)
    return fig, axes


def replay(d):
    fig, axes = frame(2, 84, left=.105, right=.97, bottom=.24, top=.87, gap=.42)
    ax, compare = axes
    valid = d["mask"].astype(bool)
    yy, xx = np.indices(valid.shape) + 1
    ax.contourf(xx, yy, valid, levels=[.5, 1.5], colors=["#eceeeb"])
    ax.contourf(xx, yy, valid & (yy >= 90), levels=[.5, 1.5], colors=["#c4dbd4"])
    ax.contour(xx, yy, valid, levels=[.5], colors=[GREY], linewidths=.6)
    fixed = read(LAB / "pose_check", "six_points.csv")
    ax.scatter(fixed.image_column, fixed.image_row, s=16, c=INK, marker="s")
    selected = read(LAB / "pose_refinement", "calibration_points.csv").query("design == 'physics6'")
    ax.scatter(selected.column, selected.row, s=28, facecolors="none", edgecolors=RUST, marker="D", linewidth=.9)
    ax.axhline(89.5, ls="--", lw=.8, color=TEAL)
    ax.set(xlim=(16, 145), ylim=(151, 22), xticks=[30, 60, 90, 120], yticks=[30, 60, 90, 120, 150],
           xlabel="Pixel column", ylabel="Pixel row")
    ax.set_aspect("equal", adjustable="datalim")
    title(ax, "a", "Six-site layouts and the reserved region")
    scores = read(LAB / "pose_check", "pose_holdout_scores.csv").query("region == 'held_out_rows_90_end'")
    for idx, rec in enumerate(scores.itertuples()):
        values = [rec.fixed_L2_percent, rec.gain_corrected_L2_percent]
        color = (RUST, TEAL)[idx]
        compare.plot([0, 1], values, marker=("^", "o")[idx], color=color, label=rec.pose, ms=5)
    rule(compare)
    compare.set(xlim=(-.18, 1.18), ylim=(0, 220), xticks=[0, 1], xticklabels=["Fixed current", "Gain fitted"],
                ylabel="Held-out field error (%)", yticks=[0, 50, 100, 150, 200])
    title(compare, "b", "5,325 pixels held out")
    key(fig, [mark(INK, "s", "Fixed sites"), mark(RUST, "D", "Sensitivity sites", filled=False),
              Line2D([], [], color=RUST, marker="^", label="Original pose"),
              Line2D([], [], color=TEAL, marker="o", label="Corrected pose")], .022)
    return fig, axes


def residual():
    row = read(MRI, "profile_row50.csv")
    fig, ax = plt.subplots(figsize=(89/25.4, 68/25.4))
    fig.subplots_adjust(left=.20, right=.96, bottom=.33, top=.88)
    ax.axhline(0, color=GREY, lw=.65)
    for column, color, style, label in [
        ("v1_1_uT", TEAL, "-", "Corrected model"), ("author_uT", INK, "--", "Authors’ simulation")
    ]:
        ax.plot(row.column_one_based, row[column]-row.MRI_uT, color=color, ls=style, lw=1.1, label=label)
    ax.set(xlabel="Pixel column in row 50", ylabel="Prediction − MRI (µT)", xticks=[40, 60, 80, 100, 120])
    ax.set_title("Residual along the same profile", pad=7)
    rule(ax)
    key(fig, ax.get_legend_handles_labels()[0], .025)
    return fig, [ax]


def axial():
    folder = LAB / "pose_refinement"
    table = read(folder, "pose_ablation.csv")
    assert (table.test_L2_percent > 0).all()
    control = read(folder, "corrected_start_control.csv")
    ref = read(folder, "reference_only.csv")
    fig, axes = frame(2, 88, left=.10, right=.975, bottom=.29, top=.87, gap=.42)
    order = ["unchanged", "gain", "translation", "translation_gain"]
    for ax, tab, baseline in zip(axes, [table, control], [False, True]):
        for idx, (group, color, marker, label) in enumerate([
            ("fixed6", BLUE, "o", "Fixed six sites"), ("physics6", TEAL, "D", "Sensitivity-selected sites")
        ]):
            rows = tab[tab.design.eq(group)].set_index("model").loc[order]
            ax.scatter(np.arange(4) + (idx - .5)*.16, rows.test_L2_percent, s=25, marker=marker,
                       facecolor="white" if idx==0 else color, edgecolor=color, linewidth=1.1, label=label, zorder=3)
        ax.set(xticks=range(4), xticklabels=["Frozen", "Gain", "Shift", "Shift\n+ gain"],
               xlim=(-.40, 3.4), ylabel="Held-out field error (%)")
        rule(ax)
        if baseline:
            ax.set_ylim(23.4, 26.3)
            ax.yaxis.set_major_locator(FixedLocator([23.5, 24, 24.5, 25, 25.5, 26]))
        else:
            ax.set_yscale("log")
            ax.set_ylim(18, 260)
            ax.yaxis.set_major_locator(FixedLocator([20, 40, 80, 160]))
            ax.yaxis.set_major_formatter(ScalarFormatter())
            ax.yaxis.set_minor_locator(FixedLocator([]))
    reference = ref.loc[ref.reference.str.contains("v1.1") & ~ref.reference.str.contains("only"), "test_L2_percent"]
    assert len(reference) == 1
    for ax in axes:
        ax.axhline(reference.iloc[0], color=GREY, lw=.85, ls="--", zorder=1)
    title(axes[0], "a", "Starting from the old pose")
    title(axes[1], "b", "Starting from the corrected pose")
    key(fig, [mark(BLUE, "o", "Fixed six sites", filled=False), mark(TEAL, "D", "Sensitivity-selected sites"),
              Line2D([], [], color=GREY, ls="--", label="Published corrected pose")], .015)
    return fig, axes


def sensors():
    folder = LAB / "sensor_design"
    sites = read(folder, "sensor_locations.csv")
    pool = read(folder, "candidate_weights.csv")
    scores = read(folder, "layout_summary.csv").query("fit == 'Joint'")
    fig, axes = frame(2, 89, left=.09, right=.97, bottom=.29, top=.865, gap=.43)
    ax, metric = axes
    ax.scatter(pool.x_mm, pool.y_mm, color="#c6ceca", s=2)
    for method, name, color, shape in zip(STYLE, SITE_NAMES, SITE_COLORS, SITE_MARKS):
        chosen = sites[sites.layout.eq(method)]
        assert len(chosen) == 6
        ax.scatter(chosen.x_mm, chosen.y_mm, s=23, facecolors="none", edgecolors=color, marker=shape, linewidth=.95)
    ax.set(xlim=(-63, 63), ylim=(-63, 63), xlabel="Scanner x (mm)", ylabel="Scanner y (mm)",
           xticks=[-56, -28, 0, 28, 56], yticks=[-56, -28, 0, 28, 56])
    ax.set_aspect("equal", adjustable="datalim")
    title(ax, "a", "Six-site layouts")
    clouds = [("Random", "#aeb9b6", "Random"), ("Shuffled", "#a690bc", "Shuffled weights")]
    for prefix, color, label in clouds:
        part = scores[scores.layout.str.startswith(prefix)]
        assert len(part) == 24
        metric.scatter(part.mean_position_error_mm, part.mean_test_L2_pct, s=7, color=color, linewidth=0, label=label)
    for method, name, color, shape in zip(STYLE, SITE_NAMES, SITE_COLORS, SITE_MARKS):
        rec = scores[scores.layout.eq(method)].iloc[0]
        metric.scatter(rec.mean_position_error_mm, rec.mean_test_L2_pct, s=32, c=color, marker=shape, zorder=3,
                       edgecolor="white", linewidth=.45)
    metric.set(xlabel="Mean position error (mm)", ylabel="Mean held-out field error (%)",
               xlim=(-.035, 1.1), ylim=(2.39, 3.45), yticks=[2.4, 2.6, 2.8, 3., 3.2, 3.4])
    rule(metric)
    title(metric, "b", "Position recovery and field prediction")
    handles = [mark(c, s, n) for c, s, n in zip(SITE_COLORS, SITE_MARKS, SITE_NAMES)]
    handles += [mark(c, ".", n) for _, c, n in clouds]
    key(fig, handles, .015)
    return fig, axes


def decay():
    scores = read(LAB / "decay_calibration", "test_scores.csv")
    assert len(scores) == 30
    fig, axes = frame(2, 86, left=.15, right=.985, bottom=.255, top=.87, gap=.11)
    styles = [("Cool-B65", BLUE, "o"), ("Cool-D-B80", RUST, "^"), ("C-B60", TEAL, "D")]
    order = ["one_anchor", "gain_only", "free_offset", "shrunken_offset", "guarded_offset"]
    labels = ["One anchor", "Four-point gain", "Gain + offset", "Shrunken offset", "Gated offset"]
    for panel, ax in enumerate(axes):
        for idx, (coil, color, symbol) in enumerate(styles):
            chosen = scores.coil.eq(coil) & (scores.test_region.eq("30_60_common") if panel == 0 else ~scores.test_region.eq("30_60_common"))
            rows = scores[chosen].set_index("method").loc[order]
            assert len(rows) == 5
            ax.plot(rows.L2_percent, np.arange(5) + (idx-1)*.14, marker=symbol, color=color,
                    markersize=4, linewidth=.75, alpha=.95)
        ax.set(yticks=range(5), yticklabels=labels if panel == 0 else [], ylim=(4.5, -.5), xlim=(0, 28),
               xticks=[0, 5, 10, 15, 20, 25], xlabel="Relative field error (%)")
        if panel:
            ax.tick_params(axis="y", length=0)
        rule(ax, "x")
    title(axes[0], "a", "Shared test range • 30–60 mm")
    title(axes[1], "b", "Extended tests by laboratory")
    key(fig, [mark(color, symbol, coil) for coil, color, symbol in styles], .009)
    return fig, axes


def main():
    parser = ArgumentParser(description=__doc__)
    parser.add_argument("--audit-tools", type=Path, help="Directory with audit_panel_alignment.py (optional)")
    parser.add_argument("--qa-dir", type=Path, help="Output directory for alignment audit reports")
    parser.add_argument("--svg-dir", type=Path, help="Also export figures as SVG into this directory")
    args = parser.parse_args()
    if args.audit_tools:
        if args.qa_dir is None:
            parser.error("--audit-tools requires --qa-dir")
        sys.path.insert(0, str(args.audit_tools))
        from audit_panel_alignment import require_matplotlib_panel_alignment
        args.qa_dir.mkdir(parents=True, exist_ok=True)
    if args.svg_dir:
        args.svg_dir.mkdir(parents=True, exist_ok=True)
    d = loadmat(MRI / "figure8_slice_results.mat")
    jobs = [("mri_field_maps", lambda: maps(d)), ("mri_pose_profile", lambda: projection(d)),
            ("mri_profile_residual", residual),
            ("six_point_replay", lambda: replay(d)), ("axial_refinement", axial),
            ("sensor_design", sensors), ("depth_calibration", decay)]
    for name, build in jobs:
        fig, panels = build()
        fig.canvas.draw()
        if args.audit_tools and len(panels) > 1:
            require_matplotlib_panel_alignment(fig, axes=panels, panel_ids=list(string.ascii_lowercase[:len(panels)]),
                                               json_out=args.qa_dir / f"{name}.alignment.json", strict=True)
        fig.savefig(HERE / f"{name}.pdf", dpi=600, metadata={"Title": name.replace("_", " ")})
        if args.svg_dir:
            fig.savefig(args.svg_dir / f"{name}.svg")
        fig.savefig(HERE / f"{name}.png", dpi=600)
        plt.close(fig)


if __name__ == "__main__":
    main()

