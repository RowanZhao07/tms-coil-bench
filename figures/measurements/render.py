"""Redraw the published-measurement figures from the project's measured data tables.

Run from any directory. Supply --audit-script only for the optional external
layout checker; ordinary reproduction needs only NumPy, SciPy, Matplotlib and pandas.
"""
from __future__ import annotations

import argparse
import hashlib
import importlib.util
from pathlib import Path
import string

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib.lines import Line2D
from matplotlib.ticker import ScalarFormatter
import numpy as np
import pandas as pd
from scipy.io import loadmat

HERE = Path(__file__).resolve().parent
ROOT = HERE.parents[1]
DATA = ROOT / "external_measurements"
INK, BLUE, ORANGE, GREY = "#26313B", "#006B8F", "#AE5300", "#697782"
PALE = "#E7F1F4"
WIDTH_MM = 183

plt.rcParams.update({
    "font.family": "sans-serif", "font.sans-serif": ["DejaVu Sans", "Arial"],
    "font.size": 7.8, "axes.labelsize": 8,
    "axes.titlesize": 8.4, "axes.titleweight": "normal", "axes.titlepad": 8,
    "axes.spines.top": False, "axes.spines.right": False, "axes.linewidth": .65,
    "axes.edgecolor": INK, "axes.labelcolor": INK, "text.color": INK,
    "xtick.color": INK, "ytick.color": INK, "xtick.labelsize": 7.2, "ytick.labelsize": 7.2,
    "xtick.direction": "out", "ytick.direction": "out",
    "xtick.major.width": .65, "ytick.major.width": .65,
    "xtick.major.size": 2.5, "ytick.major.size": 2.5,
    "legend.frameon": False, "legend.fontsize": 7.2, "legend.handlelength": 2,
    "pdf.fonttype": 42, "pdf.use14corefonts": False, "svg.fonttype": "none",
    "savefig.facecolor": "white", "lines.linewidth": 1.25,
})


def read(folder, name):
    return pd.read_csv(DATA / folder / (name + ".csv"))


def row(titles, height=77, bottom=.25):
    fig, axes = plt.subplots(1, len(titles), figsize=(WIDTH_MM / 25.4, height / 25.4), squeeze=False)
    axes = axes[0]
    fig.subplots_adjust(left=.092, right=.975, bottom=bottom, top=.84,
                        wspace=.34 if len(titles) == 2 else .43)
    for ax, title, letter in zip(axes, titles, string.ascii_lowercase):
        ax.set_title(title, loc="left")
        ax.text(-.16, 1.115, letter, transform=ax.transAxes, va="bottom",
                fontsize=10, weight="bold", clip_on=False, gid="panel-label")
    return fig, axes


def bottom_key(fig, handles, ncol=None):
    fig.legend(handles=handles, loc="lower center", bbox_to_anchor=(.52, .012),
               ncol=ncol or len(handles), columnspacing=1.7, borderaxespad=0)


def linekey(label, color=BLUE, style="-"):
    return Line2D([], [], color=color, ls=style, label=label)


def dotkey(label, color=INK, marker="o", filled=False):
    return Line2D([], [], ls="none", marker=marker, ms=4,
                  markerfacecolor=color if filled else "white", markeredgecolor=color, label=label)


def dots(ax, x, y, color=INK, marker="o", size=15, filled=False, **kwargs):
    return ax.scatter(x, y, s=size, marker=marker, linewidths=.8,
                      facecolors=color if filled else "white", edgecolors=color, zorder=3, **kwargs)


def scan():
    t = read("measured_coil_library", "mc_b70_heldout_points")
    t = t[np.abs(t.z_mm-23.8) < .2]
    assert len(t) == 2642
    # A uniform 7 mm raster gives one equally sized cell per observation.
    # The supplementary 3 mm central raster remains in the all-point scores.
    t = t[np.isclose(t.z_mm, t.z_mm.mode().iloc[0])].sort_values(["y_mm", "x_mm"])
    assert len(t)==41**2 and len(t[["x_mm","y_mm"]].drop_duplicates())==len(t)
    truth = t[[f"measured_B{v}_uT_A" for v in "xyz"]].to_numpy()
    model = t[[f"predicted_B{v}_uT_A" for v in "xyz"]].to_numpy()
    residual = np.linalg.norm(model-truth, axis=1)
    scale = np.ceil(np.abs(np.r_[truth[:, 2], model[:, 2]]).max()/10)*10
    resid_max = float(np.ceil(residual.max()))
    fig, axes = row(["Measured $B_z$", "Predicted $B_z$", "Vector residual"], height=84, bottom=.34)
    fig.subplots_adjust(left=.079, right=.979, wspace=.38)
    for i, (ax, values) in enumerate(zip(axes, (truth[:, 2], model[:, 2], residual))):
        norm = plt.Normalize(-scale, scale) if i < 2 else plt.Normalize(0, resid_max)
        cloud = ax.pcolormesh(np.sort(t.x_mm.unique()),np.sort(t.y_mm.unique()),
                              values.reshape(41,41),shading="nearest",
                              cmap="RdBu_r" if i < 2 else "viridis",norm=norm,
                              edgecolors="face",linewidth=.15)
        ax.set(xlim=(-145,145), ylim=(-145,145), xticks=(-100,0,100), yticks=(-100,0,100),
               xlabel="x (mm)", ylabel="y (mm)" if i == 0 else "")
        ax.set_aspect("equal")
        fig.canvas.draw()
        box = ax.get_position()
        barax = fig.add_axes([box.x0, .153, box.width, .026], label=f"<colorbar-{i}>")
        cb = fig.colorbar(cloud, cax=barax, orientation="horizontal",
                          ticks=[-scale, 0, scale] if i < 2 else [0, resid_max/2, resid_max])
        cb.set_label("Flux density (µT/A)" if i < 2 else "Residual norm (µT/A)", labelpad=2)
        cb.solids.set_rasterized(False)
        cb.solids.set_edgecolor("face")
        cb.outline.set_linewidth(.45)
    return fig, axes


def waveform():
    folder = "measured_coil_library"
    scores = read(folder, "waveform_scores").sort_values("relative_L2_percent")
    manifest = read(folder, "waveform_manifest").set_index("filename")
    fig, axes = row(["Magstim DCC", "Magstim D70", "All 25 configurations"], height=78, bottom=.29)
    for ax, name in zip(axes, ("Magstim_DCC", "Magstim_D70")):
        fname = name + ".tsv"
        path = ROOT / "tmp" / folder / fname
        if not path.exists():
            raise SystemExit("Run external_measurements/measured_coil_library/prepare_inputs.py first.")
        data = path.read_bytes()
        sha = hashlib.sha1(b"blob " + str(len(data)).encode() + b"\0" + data).hexdigest()
        assert sha == manifest.loc[fname, "git_blob"]
        time, obs, fit = np.loadtxt(path, comments="#").T
        score = 100*np.linalg.norm(obs-fit)/np.linalg.norm(obs)
        assert np.isclose(score, scores.loc[scores.coil == name, "relative_L2_percent"].iloc[0])
        ax.plot(time*1e6, obs, color=INK, lw=1.05)
        ax.plot(time*1e6, fit, "--", color=ORANGE, lw=1.2)
        ax.set(xlabel="Time (µs)", xlim=(-20,500), xticks=[0,200,400], ylim=(-1.32,1.22),
               yticks=[-1,0,1], ylabel="Pickup signal (a.u.)" if ax is axes[0] else "")
    ax = axes[2]
    err = scores.relative_L2_percent.to_numpy()
    assert len(err) == 25
    dots(ax, np.arange(1,26), err, BLUE, size=11, filled=True)
    ax.axhline(np.median(err), color=GREY, ls=":", lw=1)
    ax.set(xlabel="Rank by error", ylabel="Relative L2 error (%)", xlim=(0,26),
           xticks=[1,13,25], ylim=(0,20), yticks=[0,5,10,15,20])
    bottom_key(fig, [linekey("Measured waveform", INK), linekey("Published fit", ORANGE, "--"),
                     linekey("Median error", GREY, ":")])
    return fig, axes


def france():
    t = read("lyon", "air_profiles")
    fig, axes = row(["Cool-B65", "Cool-D-B80"], height=81, bottom=.27)
    for ax, (coil, d) in zip(axes, t.groupby("coil",sort=True)):
        d = d.sort_values("distance_mm")
        near = d.distance_mm <= 60
        ax.axvspan(15,60,color=PALE, lw=0)
        ax.plot(d.distance_mm[near], d.model_scaled_T[near]*1000,color=BLUE)
        ax.scatter(d.distance_mm, d.model_scaled_T*1000, s=13, marker="D", color=BLUE, zorder=2)
        dots(ax, d.distance_mm, d.B_T*1000, size=13)
        anchor = d[d.distance_mm==10]
        dots(ax, anchor.distance_mm, anchor.B_T*1000, ORANGE, "s", 24, filled=True)
        assert np.isclose(anchor.B_T.iloc[0],anchor.model_scaled_T.iloc[0])
        ax.set(yscale="log", xlabel="Distance from contact (mm)", xlim=(-5,147), xticks=[0,20,40,60,100,140],
               ylim=(3,1400), yticks=[5,10,50,100,500,1000],
               ylabel="Peak field (mT)" if ax is axes[0] else "")
        ax.yaxis.set_major_formatter(ScalarFormatter())
        ax.tick_params(which="minor",length=1.5)
    bottom_key(fig,[dotkey("Measured"),dotkey("Model at sampled depths",BLUE,"D",True),
                    dotkey("10 mm anchor",ORANGE,"s",True)])
    return fig, axes


def germany():
    data = read("regensburg", "distance_profiles")
    fig, axes = row(["C-B60", "D-B80"], height=80, bottom=.28)
    for ax, (name, d) in zip(axes,data.groupby("coil",sort=True)):
        ax.axvspan(41,70,color=PALE,lw=0)
        ax.plot(d.distance_mm,d.distance_exp_T*1000,color=BLUE)
        if name=="C-B60":
            ax.plot(d.distance_mm,d.CCD_normal_scaled_T*1000,ls="--",color=ORANGE)
        dots(ax,d.distance_mm,d.digitized_peak_T*1000,size=5)
        ax.set(xlabel="Axial distance (mm)",ylabel="Peak Hall component (mT)" if ax is axes[0] else "",
               xlim=(-1,102),xticks=[0,25,50,75,100],ylim=(-20,1050),yticks=[0,250,500,750,1000])
    bottom_key(fig,[dotkey("Digitised Hall mean"),linekey("Near-range exponential"),
                    linekey("Same-model dipoles",ORANGE,"--")])
    return fig, axes


def usa():
    obs=read("abbasi","hall_depths")
    tx=read("abbasi","current_transfer")
    assert len(obs)==14 and len(tx)==6
    fig,axes=row(["Absolute field at measured depths","Transfer to higher current"],height=82,bottom=.29)
    for ax,limit in zip(axes,(270,210)):
        ax.plot([0,limit],[0,limit],":",color=GREY,lw=.85)
        ax.set(xlim=(0,limit),ylim=(0,limit),xlabel="Measured peak (mT)",ylabel="Predicted peak (mT)")
        ax.set_aspect("equal",adjustable="box")
    dots(axes[0],obs.Hall_mT,obs.thin_loop_mT,ORANGE,"s")
    dots(axes[0],obs.Hall_mT,obs.author_COMSOL_mT,BLUE,"o")
    for name, marker, color in [("A","o",BLUE),("B","s",ORANGE)]:
        d=tx[tx.coil==name]
        dots(axes[1],d.measured_high_drive_mT,d.predicted_high_drive_mT,color,marker,20,True)
    bottom_key(fig,[dotkey("Thin-loop model",ORANGE,"s"),dotkey("Published COMSOL",BLUE),
                    dotkey("Coil A",BLUE,"o",True),dotkey("Coil B",ORANGE,"s",True),
                    linekey("Agreement",GREY,":")])
    return fig,axes


def italy():
    t=read("independent_lab","italy_results")
    fig,axes=row(["Absolute peak field","Error by measurement position"],height=79,bottom=.28)
    ax=axes[0]
    for dx,col,color,marker in [(-.19,"measured_peak_uT",INK,"o"),(0,"predicted_peak_uT",BLUE,"s"),
                                (.19,"published_sim4life_uT",ORANGE,"^")]:
        dots(ax,t.point+dx,t[col],color,marker,18,filled=(dx==0))
    ax.set(yscale="log",ylabel="Peak field (µT)",ylim=(35,1.2e6))
    ax=axes[1]
    ax.axhspan(-5,5,color=PALE,lw=0)
    ax.axhline(0,color=GREY,lw=.7)
    dots(ax,t.point-.07,t.relative_error_pct,BLUE,"s",18,True)
    dots(ax,t.point+.07,t.published_error_pct,ORANGE,"^",18)
    ax.set(ylabel="Signed relative error (%)",ylim=(-22,26),yticks=[-20,-10,0,10,20])
    for ax in axes: ax.set(xlabel="Position in instrument report",xlim=(.5,7.5),xticks=range(1,8))
    bottom_key(fig,[dotkey("Measured"),dotkey("MATLAB rings",BLUE,"s",True),dotkey("Published Sim4Life",ORANGE,"^")])
    return fig,axes


def finland():
    d=loadmat(DATA/"independent_lab"/"finland_profiles.mat")
    profiles,coords=d["probePeaks"],d["values"]
    assert profiles.shape==(100,30) and coords.shape==(31,3)
    fig,axes=row(["Thirty targeting records","Field centroid in 31 workbook rows"],height=88,bottom=.31)
    heat=axes[0].pcolormesh(np.arange(.5,31),np.arange(-.5,100),profiles,shading="flat",
                            cmap="cividis",edgecolors="face",linewidth=.15)
    axes[0].set(xlabel="Target record",ylabel="Probe position ID",xticks=[1,10,20,30],yticks=[0,25,50,75,99])
    ax=axes[1]
    ax.axhline(0,lw=.6,c=GREY,ls=":"); ax.axvline(0,lw=.6,c=GREY,ls=":")
    points=ax.scatter(coords[:,0],coords[:,1],c=coords[:,2],cmap="viridis",s=16,
                      edgecolors=INK,lw=.4, vmin=25, vmax=30)
    ax.set(xlabel="Centroid x (mm)",ylabel="Centroid y (mm)",xlim=(-1.8,1.15),ylim=(-.52,.52),
           xticks=[-1.5,-1,-.5,0,.5,1],yticks=[-.4,-.2,0,.2,.4])
    for i,(ax,obj,label) in enumerate(zip(axes,[heat,points],["Peak pickup signal (a.u.)","Peak electric field (V/m)"])):
        box=ax.get_position()
        cb=fig.colorbar(obj,cax=fig.add_axes([box.x0+.04,.115,box.width-.08,.025],label=f"<colorbar-{i}>"),
                        orientation="horizontal",ticks=([0.05,.15,.25] if i==0 else [25,27.5,30]))
        cb.set_label(label,labelpad=2);cb.outline.set_linewidth(.5)
        cb.solids.set_rasterized(False)
        cb.solids.set_edgecolor("face")
    return fig,axes


def israel():
    t=read("measured_coil_library/roth","digitized_e_fields")
    fig,axes=row(["Shared lateral probe path","Central normal of figure-of-eight coil"],height=78,bottom=.27)
    for i,coil,color,marker in [(0,"H1_left",ORANGE,"s"),(0,"figure8_left",BLUE,"o"),
                                (1,"figure8_normal",BLUE,"o")]:
        d=t[t.coil==coil].sort_values("distance_cm")
        axes[i].plot(d.distance_cm,d.field_V_per_m,color=color,lw=.95)
        dots(axes[i],d.distance_cm,d.field_V_per_m,color,marker,size=19)
    axes[0].axvline(4.5,ls=":",color=GREY,lw=1)
    for i,ax in enumerate(axes):
        ax.set(xlabel="Distance along probe path (cm)",ylabel="Induced field (V/m)" if i==0 else "",
               xlim=(0,10),ylim=(0,170),xticks=[0,2,4,6,8,10],yticks=[0,40,80,120,160])
    bottom_key(fig,[dotkey("H1",ORANGE,"s"),dotkey("Figure-of-eight",BLUE),
                    linekey("4.5 cm comparison",GREY,":")])
    return fig,axes


FIGURES = {"m01_vector_scan": scan, "m02_pickup_waveforms": waveform,
           "m03_french_decay": france, "m04_german_decay": germany,
           "m05_current_transfer": usa, "m06_circular_field": italy,
           "m07_target_repositioning": finland, "m08_saline_electric_field": israel}


def main():
    parser=argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--audit-script", type=Path, help="Optional path to audit_panel_alignment.py")
    parser.add_argument("--qa-dir",type=Path,help="Write temporary alignment audit records here")
    parser.add_argument("--print-dir",type=Path,help="Also export SVG and LZW-compressed TIFF to this directory")
    parser.add_argument("--only",choices=FIGURES,help="Render a single figure")
    args=parser.parse_args()
    auditor=None
    if args.audit_script:
        import sys
        sys.path.insert(0,str(args.audit_script.resolve().parent))
        spec=importlib.util.spec_from_file_location("audit_panel_alignment",args.audit_script)
        auditor=importlib.util.module_from_spec(spec)
        sys.modules[spec.name]=auditor
        spec.loader.exec_module(auditor)
    if args.qa_dir: args.qa_dir.mkdir(parents=True,exist_ok=True)
    if args.print_dir: args.print_dir.mkdir(parents=True,exist_ok=True)
    for name,build in FIGURES.items():
        if args.only and args.only!=name: continue
        fig,axes=build()
        if auditor:
            auditor.require_matplotlib_panel_alignment(fig,axes=axes,panel_ids=list(string.ascii_lowercase[:len(axes)]),
                json_out=args.qa_dir/(name+".alignment.json") if args.qa_dir else None)
        fig.savefig(HERE/(name+".pdf"),metadata={"Title":name.replace("_"," "), "Creator":"Matplotlib"})
        fig.savefig(HERE/(name+".png"),dpi=600)
        if args.print_dir:
            fig.savefig(args.print_dir/(name+".svg"))
            fig.savefig(args.print_dir/(name+".tiff"),dpi=600,pil_kwargs={"compression":"tiff_lzw"})
        plt.close(fig)
        print("Rendered", name)

if __name__=="__main__":
    main()
