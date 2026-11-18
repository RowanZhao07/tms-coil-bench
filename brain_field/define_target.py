"""Define one MRI-space target and cortical comparator before calculating any field."""
import hashlib
import json
from pathlib import Path

import nibabel as nib
import numpy as np


ROOT = Path(__file__).resolve().parents[1]
CACHE = ROOT / "tmp" / "brain_field"
HEAD = CACHE / "m2m_Subject01"
RESULTS = ROOT / "brain_field" / "results"
RADII_MM = (3.0, 5.0, 7.5)
PRIMARY_RADIUS_MM = 5.0


def topology_hash(mesh):
    digest = hashlib.sha256()
    for values in (mesh.elm.elm_type, mesh.elm.tag1, mesh.elm.node_number_list):
        digest.update(np.ascontiguousarray(values, dtype="<i4").view(np.uint8))
    return digest.hexdigest()


def anatomy_at(centers, volume, inverse):
    voxels = np.floor(nib.affines.apply_affine(inverse, centers) + 0.5).astype(int)
    inside = ((voxels >= 0) & (voxels < volume.shape)).all(axis=1)
    labels = np.zeros(len(centers), dtype=volume.dtype)
    labels[inside] = volume[tuple(voxels[inside].T)]
    return labels


def define(mesh, labeling, *, subject="FieldTrip Subject01", cache=CACHE, results=RESULTS):
    volume = np.asanyarray(labeling.dataobj)
    inverse = np.linalg.inv(labeling.affine)
    indices = np.column_stack(np.where(volume == 10))
    if not len(indices):
        raise ValueError("Left thalamus label 10 missing")
    roi_voxels = nib.affines.apply_affine(labeling.affine, indices)
    anatomical_centroid = roi_voxels.mean(axis=0)
    center = roi_voxels[np.argmin(np.linalg.norm(roi_voxels - anatomical_centroid, axis=1))]
    gray = np.flatnonzero((mesh.elm.elm_type == 4) & (mesh.elm.tag1 == 2))
    regional, cortical = [], []
    native_thalamus_volume = 0.0
    for start in range(0, len(gray), 60000):
        ids = gray[start:start + 60000]
        corners = mesh.nodes.node_coord[mesh.elm.node_number_list[ids, :4] - 1]
        centers = corners.mean(axis=1)
        volumes = np.abs(np.linalg.det(corners[:, 1:] - corners[:, :1])) / 6.0
        if not np.all(volumes > 0):
            raise ValueError("Non-positive gray matter tetrahedron volume")
        labels = anatomy_at(centers, volume, inverse)
        thalamus = labels == 10
        native_thalamus_volume += volumes[thalamus].sum()
        radial = np.linalg.norm(centers - center, axis=1)
        near = thalamus & (radial <= max(RADII_MM))
        regional.append((ids[near]+1, centers[near], volumes[near], radial[near]))
        cortex = np.isin(labels, [3, 42])
        cortical.append((ids[cortex]+1, volumes[cortex]))
    ids, xyz, weights, radii = (np.concatenate(xs) for xs in zip(*regional))
    ref_ids, ref_weights = (np.concatenate(xs) for xs in zip(*cortical))
    chosen = radii <= PRIMARY_RADIUS_MM
    if not np.any(radii <= min(RADII_MM)) or not len(ref_ids):
        raise ValueError("Empty target or cortical reference")
    cache.mkdir(parents=True, exist_ok=True)
    results.mkdir(parents=True, exist_ok=True)
    skin = (mesh.elm.elm_type == 2) & (mesh.elm.tag1 == 1005)
    skin_nodes, skin_faces = np.unique(mesh.elm.node_number_list[skin, :3] - 1, return_inverse=True)
    if not len(skin_nodes):
        raise ValueError("Missing skin surface tag 1005")
    np.savez_compressed(cache / "scalp.npz", vertices_mm=mesh.nodes.node_coord[skin_nodes],
                        faces=skin_faces.reshape(-1, 3))
    np.savez_compressed(
        cache / "target_roi.npz", center_mm=center, centroids_mm=xyz[chosen],
        volumes_mm3=weights[chosen], element_numbers=ids[chosen],
        reference_element_numbers=ref_ids, reference_volumes_mm3=ref_weights,
        regional_element_numbers=ids, regional_volumes_mm3=weights, distances_mm=radii,
    )
    record = {
        "subject": subject, "space": "subject MRI world RAS, mm",
        "target": "left thalamus", "target_label": 10, "gm_tissue_tag": 2,
        "anatomical_centroid_ras_mm": anatomical_centroid.tolist(),
        "center_ras_mm": center.tolist(), "radius_mm": PRIMARY_RADIUS_MM,
        "center_rule": "label-10 voxel centre nearest to the centroid of all label-10 voxel centres",
        "element_rule": "tetrahedron barycentre inside sphere, mesh GM, nearest anatomical voxel has label 10",
        "cortical_reference": "mesh GM, nearest anatomical voxel has label 3 or 42 (bilateral cerebral cortex)",
        "roi_element_count": int(chosen.sum()), "roi_volume_mm3": float(weights[chosen].sum()),
        "cortical_element_count": len(ref_ids), "cortical_volume_mm3": float(ref_weights.sum()),
        "thalamus_label_voxels": len(indices), "thalamus_mesh_gm_volume_mm3": native_thalamus_volume,
        "topology_sha256": topology_hash(mesh),
        "sensitivity_radii_mm": list(RADII_MM),
    }
    (results / "target.json").write_text(json.dumps(record, indent=2) + "\n", encoding="utf-8")
    return record


if __name__ == "__main__":
    import simnibs
    mesh = simnibs.read_msh(str(HEAD / "Subject01.msh"))
    labels = nib.load(HEAD / "segmentation" / "labeling.nii.gz")
    print(json.dumps(define(mesh, labels), indent=2))
