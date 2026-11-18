"""Register the three original windings to one MRI scalp and evaluate free-space B.

All geometry uses native right/anterior/superior coordinates in millimetres.
Run with the brain-field experiment environment. The FEM reads the
same centreline geometry and source-to-head matrices produced here.
"""
import argparse
import csv
import json
from pathlib import Path
import sys

import numpy as np
from scipy.optimize import brentq
from scipy.spatial import cKDTree
from check_primary import wire_dadt


ROOT = Path(__file__).resolve().parents[1]
PILOT = ROOT / "placement_pilot" / "placement"
sys.path.insert(0, str(PILOT))
from run_placement import CURRENT, decode, evaluate_b, transform_segments

HERE = Path(__file__).resolve().parent
CACHE = ROOT / "tmp" / "brain_field"
COILS = ("h1", "figure8", "double_cone")
GAP = 12.0
DIDT = 1e8


def dot(a, b):
    return np.einsum("ij,ij->i", a, b)


def segment_point_sq(a, b, p):
    direction = b - a
    fraction = np.clip(dot(p - a, direction) / np.maximum(dot(direction, direction), 1e-30), 0, 1)
    v = p - a - direction * fraction[:, None]
    return dot(v, v)


def inside_face(point, a, b, c, normal):
    return ((dot(np.cross(b-a, point-a), normal) >= -1e-11)
            & (dot(np.cross(c-b, point-b), normal) >= -1e-11)
            & (dot(np.cross(a-c, point-c), normal) >= -1e-11))


def segment_triangle_sq(start, end, triangles):
    """Minimum squared distance for each paired closed segment and triangle."""
    a, b, c = triangles.transpose(1, 0, 2)
    direction = end - start
    uu = dot(direction, direction)
    normal = np.cross(b-a, c-a)
    normal /= np.sqrt(dot(normal, normal))[:, None]
    h0, h1 = dot(start-a, normal), dot(end-a, normal)
    best = np.full(len(start), np.inf)

    # Endpoints projected on the triangle interior, including face crossings.
    for p, h in ((start, h0), (end, h1)):
        proj = p - h[:, None]*normal
        best = np.minimum(best, np.where(inside_face(proj, a, b, c, normal), h*h, np.inf))
    crossing = h0*h1 < 0
    alpha = np.divide(h0, h0-h1, out=np.zeros_like(h0), where=crossing)
    hit = inside_face(start+alpha[:, None]*direction, a, b, c, normal)
    best[crossing & hit] = 0

    for p, q in ((a, b), (b, c), (c, a)):
        # Bounding features for the segment-segment minimum.
        best = np.minimum(best, segment_point_sq(start, end, p))
        best = np.minimum(best, segment_point_sq(p, q, start))
        best = np.minimum(best, segment_point_sq(p, q, end))
        v, w = q-p, start-p
        vv, uv = dot(v,v), dot(direction,v)
        uw, vw = dot(direction,w), dot(v,w)
        det = uu*vv - uv*uv
        skew = det > 1e-13*uu*vv
        s = np.divide(uv*vw-vv*uw, det, out=np.zeros_like(det), where=skew)
        t = np.divide(uu*vw-uv*uw, det, out=np.zeros_like(det), where=skew)
        inside = skew & (s >= 0) & (s <= 1) & (t >= 0) & (t <= 1)
        d = w + s[:, None]*direction - t[:, None]*v
        best = np.minimum(best, np.where(inside, dot(d,d), np.inf))
    return np.maximum(best, 0)


class Surface:
    def __init__(self, triangles):
        tri = np.asarray(triangles, dtype=float)
        area2 = np.linalg.norm(np.cross(tri[:,1]-tri[:,0], tri[:,2]-tri[:,0]), axis=1)
        self.triangles = tri[area2 > 1e-10]
        self.centres = self.triangles.mean(axis=1)
        self.tree = cKDTree(self.centres)
        self.radius = np.linalg.norm(self.triangles-self.centres[:,None,:], axis=2).max()

    def clearance(self, a, b, lower_bound_only=False):
        """Exact segment-surface minimum; centroid spheres bound all excluded pairs."""
        mid = (a+b)/2
        half = np.linalg.norm(a-b, axis=1)/2
        dist, idx = self.tree.query(mid)
        lower = dist - self.radius - half
        if lower_bound_only:
            return float(lower.min())
        ix = int(dist.argmin())
        best = float(np.sqrt(segment_triangle_sq(a[ix:ix+1], b[ix:ix+1],
                                                 self.triangles[idx[ix:ix+1]])[0]))
        candidates = np.flatnonzero(lower <= best+1e-8)
        for subset in np.array_split(candidates, max(1, (len(candidates)+127)//128)):
            near = self.tree.query_ball_point(mid[subset], best+self.radius+half[subset]+1e-8)
            counts = np.fromiter(map(len, near), int)
            if not counts.any():
                continue
            ii = np.repeat(subset, counts)
            jj = np.concatenate([c for c in near if len(c)]).astype(int)
            dsq = segment_triangle_sq(a[ii], b[ii], self.triangles[jj])
            best = min(best, float(np.sqrt(dsq.min())))
        return best


def scalp_frame(skin, target):
    """Take the upper scalp intersection, smooth its normal over a 12 mm disk."""
    a,b,c = skin.triangles.transpose(1,0,2)
    v,w = b-a,c-a
    rx, ry = (target[:2] - a[:,:2]).T
    det = v[:,0]*w[:,1] - v[:,1]*w[:,0]
    valid = np.abs(det) > 1e-10
    s = np.divide(rx*w[:,1]-ry*w[:,0], det, out=np.zeros_like(det), where=valid)
    t = np.divide(v[:,0]*ry-v[:,1]*rx, det, out=np.zeros_like(det), where=valid)
    ok = valid & (s >= -1e-8) & (t >= -1e-8) & (s+t <= 1+1e-8)
    z = a[:,2] + s*v[:,2] + t*w[:,2]
    if not np.any(ok & (z > target[2])):
        raise ValueError("No superior scalp intersection above the selected target")
    origin = np.r_[target[:2], z[ok].max()]
    ids = skin.tree.query_ball_point(origin, 12)
    normals = np.cross(v[ids], w[ids])
    normals *= np.where(normals[:,2] < 0, -1., 1.)[:,None]
    nz = normals.sum(axis=0)
    nz /= np.linalg.norm(nz)
    nx = np.array([1.,0.,0.])-nz*nz[0]
    nx /= np.linalg.norm(nx)
    basis = np.column_stack((nx, np.cross(nz,nx), nz))
    assert np.isclose(np.linalg.det(basis), 1)
    return origin, basis, len(ids)


def place(skin, segments, old_matrix, origin, basis, requested):
    rot = basis @ old_matrix[:3,:3]
    shift0 = origin + basis @ old_matrix[:3,3]
    world = transform_segments(segments, rot, shift0)
    a,b = world[:,:3],world[:,3:6]
    axis = basis[:,2]
    def gap(h, lower=False):
        return skin.clearance(a+h*axis, b+h*axis, lower)
    # Begin entirely above the scalp. The distance to a fixed set is
    # 1-Lipschitz under translation, so large steps in empty space are safe.
    high = (skin.triangles.reshape(-1,3)@axis).max() - min((a@axis).min(), (b@axis).min()) + requested + 3
    while gap(high, True) > requested+1:
        high -= gap(high, True)-requested
    ghigh = gap(high)-requested
    if ghigh < -1e-7:
        raise RuntimeError("Lost the outside clearance bracket")
    for _ in range(200):
        low = high-max(1., ghigh)
        glow = gap(low)-requested
        if glow <= 0:
            break
        high, ghigh = low, glow
    else:
        raise RuntimeError("No external winding-clearance root")
    answer = brentq(lambda h: gap(h)-requested, low, high, xtol=1e-6)
    matrix = np.eye(4)
    matrix[:3,:3] = rot
    matrix[:3,3] = shift0 + answer*axis
    actual = gap(answer)
    assert abs(actual-requested)<1e-5
    return matrix, actual, answer


def write_tcd(key, source, output):
    mid = .5*(source[:,:3]+source[:,3:6])
    dl = (source[:,3:6]-source[:,:3])*source[:,6:7]
    metadata = dict(name=key, version="1.0", stimulatorList=[dict(name="Normalized current")],
                    coilElementList=[dict(type=2, name="Rigid winding", stimulator=0,
                                          points=mid.tolist(), values=dl.tolist())])
    output.mkdir(parents=True, exist_ok=True)
    (output/f"{key}.tcd").write_text(json.dumps(metadata, separators=(",",":")), encoding="utf-8")


def verify_tcd_field(key, matrix, observations, exact_mT):
    """Check source serialization, orientation, units and transform at every ROI cell."""
    coil = json.loads((CACHE/"coils"/f"{key}.tcd").read_text(encoding="utf-8"))
    elements = coil["coilElementList"][0]
    mid_m = (np.asarray(elements["points"]) @ matrix[:3,:3].T + matrix[:3,3])*1e-3
    dl_m = np.asarray(elements["values"]) @ matrix[:3,:3].T*1e-3
    mids = []
    for observation in observations:
        r_m = observation*1e-3 - mid_m
        summands = np.cross(dl_m, r_m) / np.linalg.norm(r_m, axis=1)[:,None]**3
        mids.append(1e-7 * CURRENT * 1000 * summands.sum(axis=0))
    error = np.linalg.norm(np.asarray(mids)-exact_mT) / np.linalg.norm(exact_mT)
    if error > 1e-4:
        raise ValueError(f"Wrong TCD transform, units or winding for {key}: relative B difference {error}")
    return float(error)


def casing_distance(skin, matrix):
    native = json.loads((ROOT/"reference_models"/"Brainsway_H1.tcd").read_text(encoding="utf-8"))
    shells = []
    for model in native["coilModels"]:
        vertices = decode(model["points"]) @ matrix[:3,:3].T + matrix[:3,3]
        shells.append(vertices[decode(model["faces"],np.int64)])
    shell = Surface(np.concatenate(shells))
    def edges(surface):
        v = surface.triangles
        return np.concatenate([v[:,[0,1]],v[:,[1,2]],v[:,[2,0]]])
    rim = edges(shell)
    other = edges(skin)
    return min(skin.clearance(rim[:,0],rim[:,1]), shell.clearance(other[:,0],other[:,1]))


def weighted_percentile(values, weights, percentile):
    order = np.argsort(values)
    limit = np.cumsum(weights[order])
    index = np.searchsorted(limit, percentile/100*limit[-1], side="left")
    return float(values[order[index]])


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--tcd-only", action="store_true", help="Export coil sources before the head mesh is prepared")
    args = parser.parse_args()
    sources = {key:np.loadtxt(PILOT/"source_geometry"/f"{key}.csv", delimiter=",", skiprows=1) for key in COILS}
    for key, s in sources.items():
        write_tcd(key, s, CACHE/"coils")
    if args.tcd_only:
        return
    data = np.load(CACHE/"scalp.npz")
    skin = Surface(data["vertices_mm"][data["faces"]])
    with np.load(CACHE/"target_roi.npz") as roi:
        target = roi["center_mm"]
        points, weights = roi["centroids_mm"], roi["volumes_mm3"]
        numbers = roi["element_numbers"]
    assert points.shape==(len(weights),3) and np.all(weights>0) and len(numbers)==len(weights)
    origin, basis, count = scalp_frame(skin, target)
    hist = json.loads((PILOT/"pose.json").read_text(encoding="utf-8"))
    with (PILOT/"target_fields.csv").open(encoding="utf-8") as f:
        sphere = {row["coil"]:float(row["Bnorm_mT"]) for row in csv.DictReader(f)
                  if row["case"]=="nominal" and [float(row[k]) for k in ("x_mm","y_mm","z_mm")]==[0,0,-30]}
    params = dict(system="native MRI RAS mm", target_center_mm=target.tolist(), scalp_anchor_mm=origin.tolist(),
                  superior_depth_mm=float(origin[2]-target[2]), tangent_basis_columns_RAS=basis.tolist(),
                  area_weighted_normal_patch_radius_mm=12, normal_patch_triangles=count,
                  prescribed_wire_gap_mm=GAP, coils={})
    out = HERE/"results"
    out.mkdir(exist_ok=True)
    rows = []
    for name in COILS:
        s = sources[name]
        prior = np.array(hist["coils"][name]["poses"]["nominal"]["source_to_world"])
        matrix, actual, dh = place(skin, s, prior, origin, basis, GAP)
        shell_gap = casing_distance(skin, matrix) if name=="h1" else None
        if shell_gap is not None and shell_gap < 1.e-5:
            raise ValueError("H1 casing contacts scalp at this common winding clearance")
        params["coils"][name] = dict(source_to_head_mm=matrix.tolist(), wire_gap_mm=actual,
                                     adjustment_from_sphere_mm=dh, casing_gap_mm=shell_gap,
                                     source_segments=len(s), source=f"placement_pilot/placement/source_geometry/{name}.csv")
        print(f"{name}: winding {actual:.6f} mm; housing {shell_gap}; axial shift {dh:.5f} mm", flush=True)
    # Save the common placement before the air-field integration, so the
    # independent finite-element solver can start immediately.
    (out/"placement.json").write_text(json.dumps(params, indent=2)+"\n", encoding="utf-8")
    for name in COILS:
        matrix = np.asarray(params["coils"][name]["source_to_head_mm"])
        winding = transform_segments(sources[name], matrix[:3,:3], matrix[:3,3])
        fields = 1000*evaluate_b(winding, np.vstack((target,points)))
        values = np.linalg.norm(fields,axis=1)
        # Hold the original spherical placement fixed and change only depth.
        prior = np.asarray(hist["coils"][name]["poses"]["nominal"]["source_to_world"])
        sphere_winding = transform_segments(sources[name], prior[:3,:3], prior[:3,3])
        sphere_deep = evaluate_b(sphere_winding, [0,0,-params["superior_depth_mm"]])
        tcd_error = verify_tcd_field(name, matrix, points, fields[1:])
        # The primary contribution is -dA/dt. Its norm is reported at each
        # centroid directly, before the tissue potential is introduced.
        primary = wire_dadt(points, sources[name], matrix, DIDT)
        rows.append(dict(coil=name,I_A=CURRENT,B_center_mT=values[0],B_roi_mean_mT=np.average(values[1:],weights=weights),
                         B_roi_p95_mT=weighted_percentile(values[1:],weights,95),B_sphere30_mT=sphere[name],
                         B_sphere_same_depth_mT=1000*np.linalg.norm(sphere_deep),
                         sphere_matched_depth_mm=params["superior_depth_mm"],
                         tcd_midpoint_vs_exact_relative_L2=tcd_error,
                         didt_A_per_s=DIDT,
                         E_primary_centroid_exact_mean_V_m=np.average(np.linalg.norm(primary,axis=1),weights=weights),
                         roi_volume_mm3=weights.sum(),roi_elements=len(weights)))
        np.savetxt(out/f"air_roi_{name}.csv", np.column_stack([numbers,fields[1:],values[1:]]), delimiter=",",
                   fmt=["%d"]+["%.12g"]*4, header="element_number,Bx_mT,By_mT,Bz_mT,Bnorm_mT", comments="")
    with (out/"air_field.csv").open("w", newline="", encoding="utf-8") as f:
        writer = csv.DictWriter(f, fieldnames=list(rows[0]))
        writer.writeheader()
        writer.writerows(rows)


if __name__ == "__main__":
    main()
