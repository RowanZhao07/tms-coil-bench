"""Redraw the matched-coil comparisons from compact, tracked numerical results."""
from argparse import ArgumentParser
from pathlib import Path
import sys

import matplotlib as mpl
mpl.use("Agg")
import matplotlib.pyplot as plt
from matplotlib.collections import LineCollection, PolyCollection
from matplotlib.lines import Line2D
from matplotlib.ticker import MaxNLocator
import numpy as np
import pandas as pd


HERE = Path(__file__).resolve().parent
RESULTS = HERE / "results"
OUT = HERE.parent / "figures/brain"
COILS = ["h1", "double_cone", "figure8"]
NAMES = ["H1", "Double-cone", "Figure-8"]
COLORS = ["#087863", "#a75833", "#456687"]
MARKS = ["o", "D", "s"]
INK = "#25343d"

mpl.rcParams.update({
    "font.family": "sans-serif", "font.sans-serif": ["Arial", "DejaVu Sans"],
    "font.size": 7.5, "axes.labelsize": 7.2, "axes.titlesize": 8.5, "axes.titlelocation": "left",
    "axes.labelcolor": INK, "text.color": INK, "xtick.color": INK, "ytick.color": INK,
    "xtick.labelsize": 7, "ytick.labelsize": 7, "xtick.major.size": 2.5, "ytick.major.size": 2.5,
    "axes.spines.top": False, "axes.spines.right": False, "axes.linewidth": .65,
    "axes.edgecolor": "#849097", "pdf.fonttype": 42, "ps.fonttype": 42, "svg.fonttype": "none",
    "legend.frameon": False, "savefig.facecolor": "white",
})


def titles(fig, axes, labels):
    for idx, (ax, label) in enumerate(zip(axes, labels)):
        ax.set_title(label, pad=8, loc="left", x=.09)
        ax.text(0, 1, chr(97 + idx), va="baseline", fontweight="bold", fontsize=9.2,
                transform=ax.transAxes + mpl.transforms.ScaledTranslation(0, 8/72, fig.dpi_scale_trans))


def comparison(inputs):
    air = pd.read_csv(inputs / "air_field.csv").set_index("coil")
    brain = pd.read_csv(inputs / "summary.csv").set_index("coil")
    if not set(air.index) == set(brain.index) == set(COILS) or not brain.index.is_unique or not air.index.is_unique:
        raise ValueError("Each result table must contain exactly one record per coil.")
    if not np.allclose(air.loc[COILS, "B_roi_mean_mT"], brain.loc[COILS, "B_roi_mean_mT"]):
        raise ValueError("The magnetic-field results are inconsistent between the two tables.")
    fig, axes = plt.subplots(1, 3, figsize=(180 / 25.4, 64 / 25.4))
    fig.subplots_adjust(left=.115, right=.98, bottom=.240, top=.82, wspace=.24)
    series = [(air["B_sphere30_mT"], "Point |B| (mT)"),
              (air["B_roi_mean_mT"], "ROI mean |B| (mT)"),
              (brain["E_roi_mean_V_m"], "ROI mean |E| (V/m)")]
    for idx, (ax, (values, unit)) in enumerate(zip(axes, series)):
        values = values.loc[COILS].to_numpy(float)
        if not np.isfinite(values).all() or np.any(values < 0):
            raise ValueError(f"Invalid field in {unit}")
        for row, (value, color, mark) in enumerate(zip(values, COLORS, MARKS)):
            ax.plot([0, value], [row, row], color="#cdd4d2", lw=1.05, zorder=1, solid_capstyle="butt")
            ax.plot(value, row, color=color, marker=mark, linestyle="none", markersize=5.8)
        ax.set(yticks=[0, 1, 2], yticklabels=NAMES if idx == 0 else ["", "", ""],
               ylim=(2.52, -.52), xlabel=unit, xlim=(0, values.max() * 1.11))
        ax.tick_params(axis="y", length=0, pad=5)
        ax.xaxis.set_major_locator(MaxNLocator(nbins=4, steps=[1, 2, 2.5, 5, 10]))
    titles(fig, axes, ["Air · sphere pose", "Air · MRI pose", "Brain · MRI pose"])
    return fig, axes


def boundaries(triangles, labels):
    """Find segment endpoints on every boundary between distinct regions."""
    edges = np.sort(triangles[:, [[0, 1], [1, 2], [2, 0]]].reshape(-1, 2), axis=1)
    labels = np.repeat(labels, 3)
    order = np.lexsort((edges[:, 1], edges[:, 0]))
    edges = edges[order]
    labels = labels[order]
    same_next = np.r_[np.all(edges[:-1] == edges[1:], axis=1), False]
    inside = same_next & np.r_[labels[:-1] != labels[1:], False]
    border = ~same_next & ~np.r_[False, same_next[:-1]] & (labels != 0)
    return edges[inside | border]


def section(inputs):
    with np.load(inputs / "coronal_section.npz", allow_pickle=False) as data:
        verts, tri, tissue, roi = (data[k] for k in ("xz_mm", "triangles", "tissue", "roi"))
        fields = np.stack([data[coil] for coil in COILS])
        x, z, t1, thalamus = (data[k] for k in ("image_x_mm", "image_z_mm", "t1", "left_thalamus"))
    if not roi.any() or not np.isin(tissue[roi], [1, 2]).all():
        raise ValueError("The section must cross the brain target.")
    polys = verts[tri]
    brain = np.isin(tissue, [1, 2])
    extent = (x[0]-.5, x[-1]+.5, z[0]-.5, z[-1]+.5)
    norm = mpl.colors.Normalize(0, np.max(fields[:, brain]))
    cmap = mpl.colormaps["magma"]
    # Preserve all brain tissue and keep the physical coordinates undistorted in all panels.
    fig, panels = plt.subplots(2, 2, figsize=(180 / 25.4, 147 / 25.4))
    axes = panels.ravel()
    fig.subplots_adjust(left=.065, right=.985, bottom=.14, top=.952, wspace=.09, hspace=.30)
    tissue_edges = verts[boundaries(tri, np.where(brain, tissue, 0))]
    roi_edges = verts[boundaries(tri, roi.astype(int))]
    ax = axes[0]
    ax.imshow(t1, cmap="gray", extent=extent, origin="lower", interpolation="nearest",
              vmin=0, vmax=np.percentile(t1[t1>0], 99))
    ax.contour(x, z, thalamus, levels=[.5], colors=["#57daca"], linewidths=.9)
    for i, ax in enumerate(axes[1:]):
        ax.add_collection(PolyCollection(polys, facecolors="#dbe1df", edgecolors="none",
                                         antialiaseds=False, rasterized=True))
        ax.add_collection(PolyCollection(polys[brain], facecolors=cmap(norm(fields[i, brain])),
                                         edgecolors="none", antialiaseds=False, rasterized=True))
        ax.add_collection(LineCollection(tissue_edges, colors="#aab3ba", linewidths=.19, alpha=.5,
                                         rasterized=True))
    for i, ax in enumerate(axes):
        ax.add_collection(LineCollection(roi_edges, colors=INK, linewidths=1.4))
        ax.add_collection(LineCollection(roi_edges, colors="white", linewidths=.8))
        ax.set(xlim=extent[:2], ylim=extent[2:])
        ax.set_aspect("equal", adjustable="box")
        ax.xaxis.set_major_locator(MaxNLocator(nbins=4, steps=[1, 2, 5, 10]))
        ax.yaxis.set_major_locator(MaxNLocator(nbins=4, steps=[1, 2, 5, 10]))
        if i >= 2:
            ax.set_xlabel("MRI x (mm)")
        else:
            ax.tick_params(labelbottom=False)
        if i % 2 == 0:
            ax.set_ylabel("MRI z (mm)")
        else:
            ax.tick_params(labelleft=False)
    titles(fig, axes, ["T1 anatomy"] + NAMES)
    colorbar_ax = fig.add_axes([.115, .056, .34, .022])
    cbar = fig.colorbar(mpl.cm.ScalarMappable(norm=norm, cmap=cmap), cax=colorbar_ax, orientation="horizontal")
    cbar.locator = MaxNLocator(nbins=5)
    cbar.update_ticks()
    cbar.set_label("Electric-field magnitude (V/m)", labelpad=2)
    fig.legend(handles=[Line2D([], [], lw=1.1, color="#138e80", label="Left thalamus on T1 panel"),
                        Line2D([], [], marker="s", mfc="white", mec=INK, ms=5,
                               ls="none", label="Outlined target elements")],
               loc="lower left", bbox_to_anchor=(.55, .018), fontsize=7)
    return fig, axes


def main():
    parser = ArgumentParser(description=__doc__)
    parser.add_argument("--results", type=Path, default=RESULTS)
    parser.add_argument("--out", type=Path, default=OUT)
    parser.add_argument("--audit-tools", type=Path, help="Optional directory containing the panel-alignment audit")
    parser.add_argument("--qa-dir", type=Path)
    parser.add_argument("--svg-dir", type=Path, help="Optional destination for editable SVG exports")
    args = parser.parse_args()
    if args.audit_tools:
        if not args.qa_dir:
            parser.error("--audit-tools requires --qa-dir")
        sys.path.insert(0, str(args.audit_tools))
        from audit_panel_alignment import require_matplotlib_panel_alignment
        args.qa_dir.mkdir(parents=True, exist_ok=True)
    args.out.mkdir(parents=True, exist_ok=True)
    if args.svg_dir:
        args.svg_dir.mkdir(parents=True, exist_ok=True)
    for name, draw in [("target_comparison", comparison), ("coronal_fields", section)]:
        fig, axes = draw(args.results)
        fig.canvas.draw()
        if args.audit_tools:
            require_matplotlib_panel_alignment(fig, axes=list(axes), panel_ids=list("abcd"[:len(axes)]), strict=True,
                                               json_out=args.qa_dir / f"{name}.alignment.json")
        fig.savefig(args.out / f"{name}.pdf", dpi=600, metadata={"Title": name.replace("_", " ")})
        if args.svg_dir:
            fig.savefig(args.svg_dir / f"{name}.svg", dpi=600)
        fig.savefig(args.out / f"{name}.png", dpi=600)
        plt.close(fig)
        print(f"Wrote {name}.pdf and .png")


if __name__ == "__main__":
    main()
