"""Place three fixed windings and define one left thalamus ROI in a native RAS head.

Run with simnibs_python. A common winding clearance is specified in millimetres.
An intersecting H1 casing prevents publication of a pose for that subject.
"""
import argparse
import gc
import hashlib
import json
from pathlib import Path
import sys
from tempfile import TemporaryDirectory

import nibabel as nib
import numpy as np
import simnibs

HERE = Path(__file__).resolve().parent
sys.path.insert(0, str(HERE.parent))
import define_target
import geometry

WORK = HERE.parents[1] / "tmp" / "brain_cohort"


def sha256(path):
    with path.open("rb") as stream:
        return hashlib.file_digest(stream, "sha256").hexdigest()


def label_registration(mesh, labels):
    """Report anatomical agreement independently of any electrical field."""
    data = np.asanyarray(labels.dataobj)
    inverse = np.linalg.inv(labels.affine)
    record = {}
    for side, key in [(10, "left_thalamus"), (49, "right_thalamus")]:
        indices = np.column_stack(np.where(data == side))
        if not len(indices):
            raise ValueError(f"No label {side} in anatomical volume")
        center = nib.affines.apply_affine(labels.affine, indices.mean(axis=0))
        record[f"{key}_center_ras_mm"] = center.tolist()
    if record["left_thalamus_center_ras_mm"][0] >= record["right_thalamus_center_ras_mm"][0]:
        raise ValueError("Left and right thalamus labels are not in expected RAS order")
    for tag, expected, name in [(1, [2, 41], "white_matter"), (2, [3, 42], "gray_matter")]:
        cells = np.flatnonzero((mesh.elm.elm_type == 4) & (mesh.elm.tag1 == tag))
        if not len(cells):
            raise ValueError(f"Missing tetrahedra with tissue tag {tag}")
        sample = cells[np.linspace(0, len(cells)-1, min(20000, len(cells)), dtype=int)]
        coords = mesh.nodes.node_coord[mesh.elm.node_number_list[sample, :4]-1].mean(axis=1)
        selected = define_target.anatomy_at(coords, data, inverse)
        record[f"{name}_checked_elements"] = len(sample)
        record[f"{name}_fraction_with_cerebral_label"] = float(np.isin(selected, expected).mean())
    return record


def place_coils(scalp, center, gap, subject):
    origin, frame, count = geometry.scalp_frame(scalp, center)
    original = json.loads((geometry.PILOT / "pose.json").read_text(encoding="utf-8"))
    record = {
        "subject": subject, "system": "native MRI RAS mm", "target_center_mm": center.tolist(),
        "scalp_anchor_mm": origin.tolist(), "superior_depth_mm": float(origin[2]-center[2]),
        "tangent_basis_columns_RAS": frame.tolist(), "area_weighted_normal_patch_radius_mm": 12,
        "normal_patch_triangles": count, "prescribed_wire_gap_mm": gap, "coils": {},
    }
    for name in geometry.COILS:
        csv_path = geometry.PILOT / "source_geometry" / f"{name}.csv"
        source = np.loadtxt(csv_path, delimiter=",", skiprows=1)
        coils = WORK / "coils"
        coils.mkdir(parents=True, exist_ok=True)
        with TemporaryDirectory(dir=coils) as interim:
            geometry.write_tcd(name, source, Path(interim))
            supplied = Path(interim) / f"{name}.tcd"
            shared = coils / supplied.name
            if shared.is_file():
                if sha256(shared) != sha256(supplied):
                    raise ValueError(f"Shared coil source differs from current input geometry: {name}")
            else:
                supplied.replace(shared)
        old = np.asarray(original["coils"][name]["poses"]["nominal"]["source_to_world"])
        matrix, actual, shift = geometry.place(scalp, source, old, origin, frame, gap)
        shell = geometry.casing_distance(scalp, matrix) if name == "h1" else None
        if shell is not None and shell < 1e-5:
            raise ValueError(f"{subject}: H1 casing intersects scalp at {gap:.2f} mm winding clearance")
        record["coils"][name] = {
            "source_to_head_mm": matrix.tolist(), "wire_gap_mm": actual, "casing_gap_mm": shell,
            "adjustment_from_sphere_mm": shift, "source_segments": len(source),
            "source": csv_path.relative_to(geometry.ROOT).as_posix(), "source_sha256": sha256(csv_path),
        }
    return record


def prepare(subject, mesh_path, label_path, gap, results_dir):
    cache = WORK / subject
    dst = results_dir / subject
    dst.mkdir(parents=True, exist_ok=True)
    mesh = simnibs.read_msh(str(mesh_path))
    labels = nib.load(label_path)
    registration = label_registration(mesh, labels)
    target = define_target.define(mesh, labels, subject=subject, cache=cache, results=dst)
    target.update(mesh_sha256=sha256(mesh_path), labels_sha256=sha256(label_path),
                  label_registration=registration)
    (dst / "target.json").write_text(json.dumps(target, indent=2)+"\n", encoding="utf-8")
    del mesh, labels
    gc.collect()
    with np.load(cache / "scalp.npz") as patch:
        skin = geometry.Surface(patch["vertices_mm"][patch["faces"]])
    placement = place_coils(skin, np.asarray(target["center_ras_mm"]), gap, subject)
    placement["mesh_sha256"] = target["mesh_sha256"]
    temp = dst / "placement.json.tmp"
    temp.write_text(json.dumps(placement, indent=2)+"\n", encoding="utf-8")
    temp.replace(dst / "placement.json")
    print(f"{subject}: {target['roi_element_count']} target tetrahedra, "
          f"{target['roi_volume_mm3']:.2f} mm3, scalp depth {placement['superior_depth_mm']:.2f} mm, "
          f"H1 casing clearance {placement['coils']['h1']['casing_gap_mm']:.3f} mm", flush=True)


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("subject", help="Unique identifier of the MRI subject")
    parser.add_argument("--mesh", type=Path, required=True)
    parser.add_argument("--labels", type=Path, required=True,
                        help="Native anatomical labels with FreeSurfer 10, 49, 3 and 42 (subject-RAS mm)")
    parser.add_argument("--gap-mm", type=float, default=geometry.GAP)
    parser.add_argument("--results-dir", type=Path, default=HERE/"results")
    args = parser.parse_args()
    if Path(args.subject).name != args.subject or args.subject in (".", "..") or args.gap_mm <= 0:
        parser.error("Use a single subject ID and a positive clearance")
    prepare(args.subject, args.mesh.resolve(), args.labels.resolve(), args.gap_mm, args.results_dir)
