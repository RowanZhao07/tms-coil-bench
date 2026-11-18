"""Redraw the coil atlas and numerical validation from the saved numerical outputs."""
from argparse import ArgumentParser
import csv
from pathlib import Path
import sys

import matplotlib as mpl
mpl.use("Agg")
import matplotlib.pyplot as plt
from matplotlib.collections import LineCollection
from matplotlib.lines import Line2D
from matplotlib.ticker import MaxNLocator
import numpy as np
from scipy.io import loadmat


HERE = Path(__file__).resolve().parent
ROOT = HERE.parents[1]
RES = ROOT / "results"
PILOT = ROOT / "placement_pilot"
FEM = ROOT / "fem" / "results"
WIDTH = 180
INK, TEAL, ORANGE, LIGHT = "#273843", "#00786e", "#b85419", "#9ba9ac"
BLUE = "#426f9b"
LABELS = {
    "circular": "Circular", "h1": "H1", "h1_simplified": "Simplified H1",
    "h1_alternative": "Alternative H1", "pediatric_figure8": "Pediatric figure-eight",
    "elliptical": "Elliptical", "double_cone": "Double-cone", "pancake": "Pancake",
    "pediatric_circular": "Pediatric circular", "t_coil": "T-shaped",
    "butterfly": "Butterfly", "figure8_array": "Figure-eight array",
    "folded_figure8": "Folded figure-eight", "folded_circular": "Folded circular",
}
mpl.rcParams.update({
    "font.family": "sans-serif", "font.sans-serif": ["Arial", "DejaVu Sans"],
    "font.size": 7.4, "axes.labelsize": 7.4,
    "axes.titlesize": 8.3, "axes.titlelocation": "left", "axes.titlepad": 6,
    "axes.spines.top": False, "axes.spines.right": False, "axes.linewidth": .55,
    "xtick.labelsize": 7, "ytick.labelsize": 7, "xtick.major.size": 2.7,
    "ytick.major.size": 2.7, "xtick.major.width": .5, "ytick.major.width": .5,
    "xtick.major.pad": 2.4, "ytick.major.pad": 2.4,
    "axes.edgecolor": INK, "axes.labelcolor": INK, "text.color": INK,
    "xtick.color": INK, "ytick.color": INK, "legend.fontsize": 6.8,
    "legend.frameon": False, "pdf.fonttype": 42, "ps.fonttype": 42,
    "svg.fonttype": "none", "axes.unicode_minus": True, "image.interpolation": "nearest",
})


def rows(path):
    with path.open(encoding="utf-8-sig", newline="") as stream:
        return list(csv.DictReader(stream))


def numeric(path):
    return np.loadtxt(path, delimiter=",", skiprows=1, ndmin=2)


def values(records, key):
    return np.array([float(r[key]) for r in records])


def figure(height):
    return plt.figure(figsize=(WIDTH / 25.4, height / 25.4))


def panel(fig, mm, panel_id, title=None):
    """Position a panel in physical millimetres, measured from bottom left."""
    width, height = fig.get_size_inches() * 25.4
    x, y, w, h = mm
    ax = fig.add_axes([x / width, y / height, w / width, h / height], label=panel_id)
    if title:
        ax.set_title(title, loc="left")
    return ax


def paper_text(fig, x, y, text, **kwargs):
    width, height = fig.get_size_inches() * 25.4
    return fig.text(x / width, y / height, text, **kwargs)


def panel_heading(fig, x, y, code, text):
    paper_text(fig, x, y, code, fontsize=9, weight="bold", va="baseline")
    paper_text(fig, x + 5, y, text, fontsize=8.3, va="baseline")


def grid(ax):
    ax.set_axisbelow(True)
    ax.grid(axis="y", color="#e5e9e9", lw=.45)


def save(fig, path, args, row_groups=(), column_groups=()):
    path.parent.mkdir(parents=True, exist_ok=True)
    if args.audit_scripts:
        from audit_panel_alignment import require_matplotlib_panel_alignment
        args.qa_dir.mkdir(parents=True, exist_ok=True)
        study_axes = [a for a in fig.axes if not a.get_label().startswith("_")]
        require_matplotlib_panel_alignment(
            fig, axes=study_axes, panel_ids=[a.get_label() for a in study_axes],
            row_groups=list(row_groups), column_groups=list(column_groups),
            json_out=args.qa_dir / (path.name + "_alignment.json"),
        )
    fig.savefig(path.with_suffix(".pdf"), dpi=600, metadata={"Title": path.stem.replace("_", " ")})
    fig.savefig(path.with_suffix(".png"), dpi=600, facecolor="white")
    if args.svg_dir:
        args.svg_dir.mkdir(parents=True,exist_ok=True)
        fig.savefig(args.svg_dir / (path.name + ".svg"), dpi=600)
    plt.close(fig)


def project(xyz):
    az, el = np.deg2rad([30, 26])
    horizontal = np.array([np.cos(az), -np.sin(az), 0])
    vertical = np.array([np.sin(az)*np.sin(el), np.cos(az)*np.sin(el), np.cos(el)])
    return np.stack([xyz @ horizontal, xyz @ vertical], axis=-1)


def windings(data, scale, basis):
    start, end = data[:, :3] * scale, data[:, 3:6] * scale
    cut = np.flatnonzero(np.linalg.norm(start[1:] - end[:-1], axis=1) > 1e-5) + 1
    blocks = np.split(np.arange(len(data)), cut)
    return [basis(np.vstack((start[block], end[block[-1]]))) for block in blocks]


def add_wires(ax, lines, color=TEAL, lw=.42):
    ax.add_collection(LineCollection(lines, colors=color, linewidths=lw))


def center_wires(lines):
    xy = np.concatenate(lines)
    midpoint = .5 * (xy.min(axis=0) + xy.max(axis=0))
    return [p - midpoint for p in lines], np.ptp(xy, axis=0)


def scale_bar(ax, length, label, pos=(.05, .055)):
    xlim, ylim = ax.get_xlim(), ax.get_ylim()
    x = xlim[0] + pos[0] * np.ptp(xlim)
    y = ylim[0] + pos[1] * np.ptp(ylim)
    ax.plot([x, x + length], [y, y], lw=1, color=INK, solid_capstyle="butt")
    ax.annotate(label, xy=(x+length/2, y), xytext=(0, 2.3), textcoords="offset points",
                ha="center", va="bottom", fontsize=6.6)


def lattice(data, x_col, y_col, v_col, a=1000, b=1000):
    x, j = np.unique(data[:, x_col] * a, return_inverse=True)
    y, i = np.unique(data[:, y_col] * b, return_inverse=True)
    if x.size * y.size != data.shape[0]:
        raise ValueError("Expected a complete rectangular sampling grid.")
    field = np.empty((len(y), len(x)))
    field[i, j] = data[:, v_col]
    return x, y, field


def draw_image(ax, x, y, b, **kwargs):
    half_x, half_y = (x[1] - x[0])/2, (y[1] - y[0])/2
    image = ax.imshow(b, extent=[x[0]-half_x, x[-1]+half_x, y[0]-half_y, y[-1]+half_y],
                      origin="lower", aspect="auto", interpolation="nearest", **kwargs)
    ax.set(xlim=(x[0], x[-1]), ylim=(y[0], y[-1]))
    return image


def observations(source, points, field):
    vector=source[field]
    return np.column_stack((source[points],vector,np.linalg.norm(vector,axis=1)))


def atlas_figure(coil, source, wires, m, args):
    xy = observations(source,"xyPts","Bxy")
    xz = observations(source,"xzPts","Bxz")
    axis = observations(source,"axisPts","Baxis")
    assert np.allclose(np.linalg.norm(xy[:, 3:6], axis=1), xy[:, 6])
    assert np.isclose(1000*xy[:, 6].max(), float(m["xy20_peak_mT"]))
    axisval = 1000*axis[:, 6]
    assert np.isclose(axisval[np.isclose(axis[:, 2], .030)].item(), float(m["B30_mT"]))
    fig = figure(139)
    paper_text(fig, 5, 133, LABELS[coil], fontsize=10, weight="bold")
    panel_heading(fig, 8, 122, "a", "Winding geometry")
    panel_heading(fig, 95, 122, "b", "Field at 20 mm")
    panel_heading(fig, 8, 47, "c", "Field on the central axis")
    panel_heading(fig, 95, 47, "d", "Field along the centre plane")
    geom = panel(fig, (11, 56, 65, 60), "a")
    lines, span = center_wires(wires)
    s = max(80, span[0]/.90, span[1]/.72)
    geom.set(xlim=(-s/2,s/2), ylim=(-s*.48,s*.44), aspect="equal")
    add_wires(geom, lines)
    scale_bar(geom, 25, "25 mm")
    geom.axis("off")
    mapxy = panel(fig, (102, 59, 57, 57), "b")
    xx, yy, plane = lattice(xy, 0, 1, 6)
    i, j, section = lattice(xz, 0, 2, 6)
    vmax = 1000 * max(plane.max(), section.max())
    norm = mpl.colors.Normalize(0, vmax)
    img = draw_image(mapxy, xx, yy, 1000*plane, cmap="viridis", norm=norm)
    mapxy.contour(xx, yy, plane, levels=[.5*plane.max()], colors="white", linewidths=.7)
    mapxy.set(xlabel="x (mm)", ylabel="y (mm)", xticks=[-150,0,150], yticks=[-150,0,150])
    cax = panel(fig, (163, 63, 2.7, 44), "_scale")
    cb = fig.colorbar(img, cax=cax)
    cb.locator = MaxNLocator(nbins=4)
    cb.update_ticks()
    cb.ax.tick_params(labelsize=6.7)
    cb.ax.set_title("mT", fontsize=7.2, loc="left", pad=3)
    curve = panel(fig, (15, 12, 61, 29), "c")
    curve.plot(1000*axis[:, 2], axisval, color=TEAL, lw=1.5)
    curve.set(xlim=(10,100), xlabel="Depth (mm)", ylabel="|B| (mT)",
              xticks=[10, 30, 50, 70, 100], ylim=(0, axisval.max()*1.07))
    curve.yaxis.set_major_locator(MaxNLocator(4))
    grid(curve)
    mapxz = panel(fig, (102, 20, 57, 14.25), "d")
    draw_image(mapxz, i, j, 1000*section, cmap="viridis", norm=norm)
    mapxz.set(xlabel="x (mm)", ylabel="Depth (mm)", xticks=[-150, 0, 150], yticks=[10, 55, 100])
    mapxz.invert_yaxis()
    save(fig, HERE / "atlas" / coil, args, column_groups=[["b","d"]])


def geometry_overview(geometries, args):
    fig = figure(151)
    all_shapes = {k:center_wires(v) for k,v in geometries.items()}
    width = max(max(span[0],1.5*span[1]) for _, span in all_shapes.values()) * 1.14
    row_groups = []
    for n, (coil,(lines,_)) in enumerate(all_shapes.items()):
        row, col = divmod(n, 4)
        x, y = 4 + col*44, 115 - row*36
        ax = panel(fig, (x, y, 39, 26), coil)
        ax.set(xlim=(-width/2,width/2), ylim=(-width/3,width/3))
        add_wires(ax, lines, lw=.37)
        ax.axis("off")
        paper_text(fig, x, y+29, LABELS[coil], fontsize=7.6)
        if col == 0:
            row_groups.append([])
        row_groups[-1].append(coil)
    key = panel(fig, (92,7,39,26), "_scale")
    key.set(xlim=(-width/2,width/2), ylim=(-width/3,width/3))
    scale_bar(key, 50, "50 mm", pos=(.12,.39))
    key.axis("off")
    save(fig, HERE / "geometry_overview", args, row_groups=row_groups)


def comparison(metrics, args):
    fig=figure(126)
    ordered=sorted(metrics, key=lambda r:float(r["B30_mT"]), reverse=True)
    names=[LABELS[r["id"]] for r in ordered]
    yy=np.arange(len(names))
    ax=panel(fig,(37,15,70,94),"a")
    panel_heading(fig, 6, 119, "a", "Central-axis field at 30 mm")
    columns=["B30_mT","B30_equal_10m_mT"]
    x, y=(values(ordered, key) for key in columns)
    ax.hlines(yy,x,y,color="#b7c3c2",lw=1)
    ax.scatter(x,yy,s=15,color=TEAL,zorder=3,label="Fixed 1 kA")
    ax.scatter(y,yy,s=18,facecolor="white",edgecolor=ORANGE,linewidths=.8,
               marker="D",zorder=3,label="Scaled to 10 m")
    ax.set(yticks=yy,yticklabels=names,ylim=(13.7,-.8),xlabel="|B| (mT)",xlim=(0,310))
    ax.legend(ncol=2,loc="lower left",bbox_to_anchor=(0,1.005),borderaxespad=0,
              handletextpad=.35,columnspacing=.9)
    ax.tick_params(axis="y",length=0)
    ax2=panel(fig,(119,15,55,94),"b")
    panel_heading(fig, 117, 119, "b","Depth retention")
    xx=100*values(ordered,"B50_mT")/values(ordered,"B10_mT")
    ax2.hlines(yy,0,xx,color="#b7c3c2",lw=.8)
    ax2.scatter(xx,yy,color=BLUE,s=16)
    ax2.set(xlim=(0,42),ylim=(13.7,-.8),yticks=yy,yticklabels=[],xlabel="|B| at 50 / 10 mm (%)")
    ax2.tick_params(axis="y",length=0)
    save(fig,HERE/"coil_comparison",args,row_groups=[["a","b"]])


def placements(args):
    geom_dir=PILOT / "placement" / "target_frame_segments"
    fig=figure(97)
    cases=["h1","figure8","double_cone"]
    names=["H1","Figure-eight","Double-cone"]
    panel_ids=[]
    for n,coil in enumerate(cases):
        x=13+n*57
        ax=panel(fig,(x,25.7,45,46.3),coil)
        panel_ids.append(coil)
        # Equal physical scales and common limits for the x-z projection.
        ax.set(xlim=(-135,110),ylim=(-175,77),xticks=[-100,0,100],yticks=[-150,-100,-50,0,50])
        segments=windings(numeric(geom_dir/(coil+".csv")),1,lambda p:p[:,[0,2]])
        th=np.linspace(0,2*np.pi,301)
        ax.plot(85*np.cos(th),85*np.sin(th)-85,color=LIGHT,lw=1)
        add_wires(ax,segments,lw=.4)
        ax.scatter([0]*3,[-10,-30,-50],marker="o",color=ORANGE,s=10,zorder=4)
        ax.set_xlabel("x (mm)")
        if n==0:
            ax.set_ylabel("z (mm)")
        else:
            ax.tick_params(labelleft=False)
        panel_heading(fig,x-6,79,chr(97+n),names[n])
    handles=[Line2D([],[],color=TEAL,label="Projected winding"),
             Line2D([],[],color=LIGHT,label="85 mm sphere"),
             Line2D([],[],lw=0,marker="o",markersize=3.5,color=ORANGE,label="Common targets")]
    fig.legend(handles=handles,ncol=3,loc="center",bbox_to_anchor=(.5,.07))
    save(fig,HERE/"placement",args,row_groups=[panel_ids])


def reference_surfaces(args):
    fig=figure(76); data=rows(PILOT/"reference_surface_comparison.csv")
    for n,(coil,title) in enumerate([("h1","H1"),("figure8","Figure-eight"),("double_cone","Double-cone")]):
        rr=[r for r in data if r["coil"]==coil]
        xx=values(rr,"depth_mm")
        x=14+n*57
        ax=panel(fig,(x,21,44,39),coil)
        panel_heading(fig,x-4,68,chr(97+n),title)
        for col,colour,fmt,label in [("sphere_norm_mT",TEAL,"o-","Sphere"),("flat_norm_mT",ORANGE,"s--","Plane")]:
            ax.plot(xx,values(rr,col),fmt,c=colour,ms=3.8,lw=1.2,mfc="white" if col=="flat_norm_mT" else colour,
                    label=label)
        ax.set(xlim=(7,53),ylim=(0,80),xticks=[10,30,50],yticks=[0,20,40,60,80],xlabel="Depth (mm)")
        if n==0:ax.set_ylabel("Target |B| (mT)")
        else:ax.tick_params(labelleft=False)
        grid(ax)
        if n==2: fig.legend(*ax.get_legend_handles_labels(),loc="center",bbox_to_anchor=(.5,.065),ncol=2)
    save(fig,HERE/"reference_surfaces",args,row_groups=[["h1","figure8","double_cone"]])


def validation(args):
    fig=figure(81)
    axes=[panel(fig,(15+i*57,21,43,42),k) for i,k in enumerate(["a","b","c"])]
    names=["Single-turn analytic check","18-turn circular field","18-turn residual"]
    for i,ax in enumerate(axes):
        panel_heading(fig,8+i*57,72,chr(97+i),names[i]);grid(ax)
    c=rows(RES/"validation"/"circle_analytic_convergence.csv")
    axes[0].loglog(values(c,"segments"),100*values(c,"max_axis_relative_error"),"o-",color=TEAL,lw=1.3,ms=3.5)
    axes[0].set(xlabel="Segments per turn",ylabel="Max. axial error (%)",xticks=[48,192,768],xticklabels=["48","192","768"])
    axes[0].set_yticks([.001,.01,.1],["0.001","0.01","0.1"])
    r=rows(FEM/"multiturn_main_circular"/"main_model_comparison.csv")
    z=1000*values(r,"z_m")
    axes[1].plot(z,1000*values(r,"Bz_MATLAB_T"),color=TEAL,lw=1.3,label="Filament model")
    axes[1].scatter(z,1000*values(r,"Bz_FEM_T"),s=16,facecolors="white",edgecolors=ORANGE,label="Finite elements",zorder=3)
    axes[1].set(xlabel="Depth (mm)",ylabel="Axial field (mT)",xticks=[10,30,50],ylim=(0,290))
    fig.legend(*axes[1].get_legend_handles_labels(),loc="center",bbox_to_anchor=(.5,.055),ncol=2,handlelength=1.3)
    axes[2].axhline(0,color=LIGHT,lw=.7)
    axes[2].plot(z,values(r,"FEM_vs_MATLAB_pct"),"o-",color=ORANGE,ms=3.5,lw=1.2)
    axes[2].set(xlabel="Depth (mm)",ylabel="FEM − filament (%)",xticks=[10,30,50],ylim=(-.075,.075),yticks=[-.06,0,.06])
    save(fig,HERE/"solver_validation",args,row_groups=[["a","b","c"]])


def cone_validation(args):
    fig=figure(91)
    ax=panel(fig,(15,24,64,49),"a")
    panel_heading(fig,8,81,"a","Mesh and domain convergence")
    data=rows(PILOT/"fem3d"/"results"/"mesh_domain_convergence.csv")
    for colour,shape,h in [(BLUE,"s",8),(LIGHT,"^",4),(ORANGE,"D",2),(TEAL,"o",1.5)]:
        rr=[r for r in data if float(r["h_mm"])==h]
        ax.plot(values(rr,"half_width_mm"),values(rr,"relative_vector_l2_pct"),marker=shape,
                ms=4,lw=1.2,color=colour,label=f"{h:g} mm")
    ax.set(xlabel="Domain half-width (mm)",ylabel="Vector relative error (%)",xticks=[192,256,320,480],ylim=(0,8))
    fig.legend(*ax.get_legend_handles_labels(),loc="center",bbox_to_anchor=(.28,.065),ncol=4,handlelength=1.1,
               columnspacing=.7,handletextpad=.3,title="Mesh step",title_fontsize=7)
    grid(ax)
    ax2=panel(fig,(109,24,49,49),"b")
    panel_heading(fig,99,81,"b","Local error at 20 mm")
    d=numeric(PILOT/"fem3d"/"results"/"final_roi.csv")
    sub=d[np.isclose(d[:,2],.02)]
    xx,yy,b=lattice(sub,0,1,10)
    img=draw_image(ax2,xx,yy,b,cmap="magma",vmin=0,vmax=b.max())
    ax2.set(xlabel="x (mm)",ylabel="y (mm)",xticks=[-90,0,90],yticks=[-90,0,90])
    cb=fig.colorbar(img,cax=panel(fig,(162,28,2.7,38),"_scale"))
    cb.set_label("Pointwise error (%)")
    save(fig,HERE/"double_cone_fem",args,row_groups=[["a","b"]])


def normalized_depth(sources,args):
    families=[
        ("Circular and elliptical",["circular","pancake","pediatric_circular","elliptical"]),
        ("H-shaped windings",["h1","h1_simplified","h1_alternative"]),
        ("Two-wing and array",["pediatric_figure8","t_coil","figure8_array","folded_figure8"]),
        ("Bent windings",["double_cone","butterfly","folded_circular"]),
    ]
    fig=figure(119)
    for n,(title,coils) in enumerate(families):
        row,col=divmod(n,2); x=15+col*89; y=72-row*55; pid=chr(97+n)
        ax=panel(fig,(x,y,72,32),pid)
        panel_heading(fig,x-7,y+38,pid,title)
        for coil,colour,style in zip(coils,[TEAL,ORANGE,BLUE,INK],["-","--","-.",":"]):
            data=observations(sources[coil],"axisPts","Baxis")
            baseline=data[np.isclose(data[:,2],.01),6].item()
            ax.plot(data[:,2]*1000,100*data[:,6]/baseline,style,color=colour,lw=1.2,label=LABELS[coil])
        ax.set(xlim=(10,100),ylim=(0,104),xticks=[10,30,50,70,100],yticks=[0,50,100])
        if not col:ax.set_ylabel("Field / 10 mm field (%)")
        if row:ax.set_xlabel("Depth (mm)")
        ax.axhline(50,color="#e5e9e9",lw=.45,zorder=0)
        ax.legend(loc="upper right",handlelength=1.9,borderpad=0,labelspacing=.15)
    save(fig,HERE/"normalized_depth",args,row_groups=[["a","b"],["c","d"]],column_groups=[["a","c"],["b","d"]])


def sensitivity(args):
    fig=figure(87)
    axes=[panel(fig,(15+i*57,24,42,43),k) for i,k in enumerate(["a","b","c"])]
    titles=["Fold angle","Pediatric geometry","Midpoint integration"]
    for i,ax in enumerate(axes):
        panel_heading(fig,8+i*57,77,chr(97+i),titles[i]);grid(ax)
    data=rows(RES/"sensitivity"/"fold_angle_scan.csv")
    for key,fmt,colour,label in [("axis_B30_mT","o-",TEAL,"Axial field"),("xy20_half_area_cm2","s--",ORANGE,"Half-peak area")]:
        y=values(data,key)
        axes[0].plot(values(data,"fold_angle_deg"),100*y/y[0],fmt,c=colour,lw=1.1,ms=3,label=label)
    axes[0].set(xlabel="Fold angle (°)",ylabel="Relative to 0° (%)",xticks=[0,30,60])
    fig.legend(*axes[0].get_legend_handles_labels(),loc="center",bbox_to_anchor=(.20,.065),ncol=1,
               handlelength=1.3,labelspacing=.15)
    dat=rows(RES/"sensitivity"/"pediatric_size_equal_NI.csv")
    axes[1].axhline(1,c=LIGHT,lw=.7)
    axes[1].plot(values(dat,"depth_mm"),values(dat,"larger_over_compact"),color=TEAL,lw=1.5)
    axes[1].set(xlabel="Depth (mm)",ylabel="Larger / compact field")
    q=rows(RES/"sensitivity"/"midpoint_subdivision_convergence.csv")
    for name,style,colour,label in [("elliptical","o-",TEAL,"Elliptical"),("h1_simplified","s--",BLUE,"Simplified H1"),
                                    ("h1_alternative","^:",ORANGE,"Alternative H1")]:
        rr=[r for r in q if r["id"]==name]
        axes[2].semilogy(values(rr,"subdivisions_per_original_segment"),100*values(rr,"midpoint_relative_L2"),style,
                         c=colour,ms=2.5,lw=1.1,label=label)
    axes[2].set(xlabel="Subdivisions",ylabel="Relative error (%)",xticks=[1,16,32])
    axes[2].set_yticks([.001,.1,10],["0.001","0.1","10"])
    fig.legend(*axes[2].get_legend_handles_labels(),loc="center",bbox_to_anchor=(.84,.072),ncol=1,
               handlelength=1.3,labelspacing=.1)
    save(fig,HERE/"model_sensitivity",args,row_groups=[["a","b","c"]])


def main():
    parser=ArgumentParser(description=__doc__)
    parser.add_argument("--audit-scripts",type=Path,help="Optional directory containing audit_panel_alignment.py")
    parser.add_argument("--qa-dir",type=Path,default=ROOT/"tmp"/"core-fig-qa")
    parser.add_argument("--svg-dir",type=Path,help="Optional destination for separately editable SVG figures")
    parser.add_argument("--only",help="One coil ID or 'summary' (all summary figures)")
    args=parser.parse_args()
    if args.audit_scripts:sys.path.insert(0,str(args.audit_scripts))
    metrics=rows(RES/"metrics.csv")
    assert {r["id"] for r in metrics} == set(LABELS)
    sources={coil:loadmat(RES/"fields"/(coil+".mat"),simplify_cells=True) for coil in LABELS}
    geom={coil:windings(sources[coil]["c"]["segments"],1000,project) for coil in LABELS}
    for m in metrics:
        if args.only in (None,m["id"]):
            atlas_figure(m["id"],sources[m["id"]],geom[m["id"]],m,args)
    if args.only in (None,"summary"):
        geometry_overview(geom,args);comparison(metrics,args);placements(args);reference_surfaces(args)
        validation(args);cone_validation(args);sensitivity(args);normalized_depth(sources,args)


if __name__=="__main__":
    main()
