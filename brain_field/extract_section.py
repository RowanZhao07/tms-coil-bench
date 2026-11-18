"""Intersect the solved tetrahedra with the coronal plane through the target."""
from argparse import ArgumentParser
import gc
import hashlib
import json
from pathlib import Path

import numpy as np


HERE = Path(__file__).resolve().parent
COILS = ("h1", "figure8", "double_cone")
EDGES = np.array([[0, 1], [0, 2], [0, 3], [1, 2], [1, 3], [2, 3]])


def mesh_digest(mesh):
    """Ensure that a saved element index identifies the same tetrahedron in each solve."""
    checksum = hashlib.sha256()
    for array in (mesh.nodes.node_coord, mesh.elm.node_number_list, mesh.elm.elm_type, mesh.elm.tag1):
        checksum.update(np.ascontiguousarray(array).view(np.uint8))
    return checksum.hexdigest()


def intersect_tetrahedra(nodes, connectivity, types, y_mm):
    """Return vertices (x,z), triangle connectivity and original row indices."""
    faces, owners = [], []
    end_a, end_b = EDGES.T
    for start in range(0, len(types), 100_000):
        rows = start + np.flatnonzero(types[start:start + 100_000] == 4)
        if not len(rows):
            continue
        identifiers = connectivity[rows] - 1  # SimNIBS node numbers are one-based.
        heights = nodes[identifiers, 1]
        cut = (heights.min(axis=1) < y_mm) & (heights.max(axis=1) >= y_mm)
        if not cut.any():
            continue
        rows = rows[cut]
        xyz = nodes[identifiers[cut]]
        below = xyz[:, :, 1] < y_mm
        edges_crossed = below[:, end_a] != below[:, end_b]
        for count in (3, 4):
            selected = edges_crossed.sum(axis=1) == count
            if not selected.any():
                continue
            verts = xyz[selected]
            batch, edges = np.nonzero(edges_crossed[selected])
            a, b = verts[batch, end_a[edges]], verts[batch, end_b[edges]]
            t = (y_mm - a[:, 1]) / (b[:, 1] - a[:, 1])
            polygon = (a[:, (0, 2)] + t[:, None] * (b[:, (0, 2)] - a[:, (0, 2)]))
            polygon = polygon.reshape(-1, count, 2)
            centred = polygon - polygon.mean(axis=1, keepdims=True)
            angle = np.arctan2(centred[:, :, 1], centred[:, :, 0])
            polygon = np.take_along_axis(polygon, np.argsort(angle, axis=1)[:, :, None], axis=1)
            for i in range(1, count - 1):
                triangle = polygon[:, [0, i, i + 1]]
                ab = triangle[:, 1] - triangle[:, 0]
                ac = triangle[:, 2] - triangle[:, 0]
                area_twice = ab[:, 0] * ac[:, 1] - ab[:, 1] * ac[:, 0]
                # A plane passing through just one vertex or one edge has zero area.
                keep = area_twice > 1e-10
                faces.append(triangle[keep])
                owners.append(rows[selected][keep])
    if not faces:
        raise ValueError(f"No tetrahedra intersect y={y_mm} mm.")
    faces = np.concatenate(faces)
    owners = np.concatenate(owners)
    vertices, inverse = np.unique(faces.round(8).reshape(-1, 2), axis=0, return_inverse=True)
    return vertices, inverse.reshape(-1, 3).astype(np.int32), owners


def anatomy_section(anatomy_dir, vertices, faces, tissue, y_mm):
    """Sample T1 intensities linearly and anatomical labels by nearest voxel in native RAS."""
    import nibabel as nib
    from scipy.ndimage import map_coordinates

    brain_vertices = vertices[faces[np.isin(tissue, [1, 2])]].reshape(-1, 2)
    lower = np.floor(brain_vertices.min(axis=0) - [11, 7])
    upper = np.ceil(brain_vertices.max(axis=0) + [11, 15])
    x = np.arange(lower[0], upper[0] + 1)
    z = np.arange(lower[1], upper[1] + 1)
    xx, zz = np.meshgrid(x, z)
    ras = np.stack((xx, np.full_like(xx, y_mm), zz), axis=-1)

    def sample(path, order):
        image = nib.load(path)
        ijk = nib.affines.apply_affine(np.linalg.inv(image.affine), ras)
        volume = np.asanyarray(image.dataobj)
        return map_coordinates(volume, np.moveaxis(ijk, -1, 0), order=order,
                               prefilter=False, mode="constant", cval=0)

    t1 = sample(anatomy_dir / "T1.nii.gz", 1)
    labels = sample(anatomy_dir / "segmentation/labeling.nii.gz", 0)
    return dict(image_x_mm=x, image_z_mm=z, t1=t1.astype(np.float32), left_thalamus=(labels == 10))


def main():
    parser = ArgumentParser(description=__doc__)
    parser.add_argument("--field-dir", type=Path, default=HERE.parent / "tmp/brain_field/fields")
    parser.add_argument("--roi", type=Path, default=HERE.parent / "tmp/brain_field/target_roi.npz")
    parser.add_argument("--anatomy-dir", type=Path, default=HERE.parent / "tmp/brain_field/m2m_Subject01")
    parser.add_argument("--results", type=Path, default=HERE / "results")
    args = parser.parse_args()
    from simnibs.mesh_tools.mesh_io import read_msh

    target = json.loads((args.results / "target.json").read_text(encoding="utf-8-sig"))
    center = np.asarray(target["center_ras_mm"], dtype=float)
    radius = float(target["radius_mm"])
    with np.load(args.roi, allow_pickle=False) as roi_data:
        roi_element_numbers = np.asarray(roi_data["element_numbers"])
    fields = {}
    fingerprint = None
    for coil in COILS:
        mesh = read_msh(str(args.field_dir / f"{coil}.msh"))
        digest = mesh_digest(mesh)
        if fingerprint is None:
            fingerprint = digest
            vertices, faces, owners = intersect_tetrahedra(
                mesh.nodes.node_coord, mesh.elm.node_number_list, mesh.elm.elm_type, center[1]
            )
            tissue = np.asarray(mesh.elm.tag1)[owners].astype(np.int32)
        elif digest != fingerprint:
            raise ValueError(f"Node coordinates, element ordering or tissue tags differ for {coil}.")
        field_vector = mesh.field["E"].value[owners]
        magnitude = np.linalg.norm(field_vector, axis=1)
        if not np.isfinite(magnitude).all():
            raise ValueError(f"Invalid field on the section for {coil}.")
        fields[coil] = magnitude.astype(np.float32)
        del mesh
        gc.collect()
    anatomy = anatomy_section(args.anatomy_dir, vertices, faces, tissue, center[1])
    np.savez_compressed(args.results / "coronal_section.npz", xz_mm=vertices, triangles=faces,
                        tissue=tissue, element_number=owners.astype(np.int32) + 1,
                        roi=np.isin(owners + 1, roi_element_numbers),
                        center_ras_mm=center, radius_mm=radius, **fields, **anatomy)
    print(f"Saved {len(faces):,} triangles from {np.unique(owners).size:,} intersected tetrahedra.")


if __name__ == "__main__":
    main()
