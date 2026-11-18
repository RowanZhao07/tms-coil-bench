"""Pair the FEM and free-space measurements at the same tetrahedra."""
import csv
import json

import numpy as np
import simnibs

from define_target import CACHE, RESULTS, RADII_MM, topology_hash


def weighted_quantile(values, weights, fraction):
    order = np.argsort(values)
    cumulative = np.cumsum(weights[order], dtype=float)
    return values[order[np.searchsorted(cumulative, fraction * cumulative[-1])]]


def write_rows(filename, rows):
    with (RESULTS / filename).open("w", newline="", encoding="utf-8") as stream:
        writer = csv.DictWriter(stream, fieldnames=rows[0].keys())
        writer.writeheader()
        writer.writerows(rows)


def main():
    meta = json.loads((RESULTS / "target.json").read_text(encoding="utf-8"))
    with (RESULTS / "air_field.csv").open(encoding="utf-8") as stream:
        air = list(csv.DictReader(stream))
    rows, sensitivity, targets = [], [], []
    for coil in air:
        name = coil["coil"]
        mesh = simnibs.read_msh(str(CACHE / "fields" / (name + ".msh")))
        if topology_hash(mesh) != meta["topology_sha256"]:
            raise ValueError(f"{name}: FEM element numbers do not match reference")
        efield = mesh.field["E"].value
        with np.load(CACHE / "target_roi.npz") as roi:
            ids = roi["element_numbers"] - 1
            centers = mesh.nodes.node_coord[mesh.elm.node_number_list[ids] - 1].mean(axis=1)
            np.testing.assert_allclose(centers, roi["centroids_mm"], atol=2e-5, rtol=0)
            if int(coil["roi_elements"]) != len(ids):
                raise ValueError("Air and FEM target selections have different sizes")
            w = roi["volumes_mm3"]
            ref = roi["reference_element_numbers"] - 1
            ref_w = roi["reference_volumes_mm3"]
            e = np.linalg.norm(efield[ids], axis=1)
            cortical = np.linalg.norm(efield[ref], axis=1)
            if not (np.all(np.isfinite(cortical)) and np.all(np.isfinite(e))):
                raise ValueError("Non-finite field in brain")
            ave = float(np.average(e, weights=w))
            p99 = float(weighted_quantile(cortical, ref_w, .99))
            b = np.genfromtxt(RESULTS / f"air_roi_{name}.csv", names=True, delimiter=",")
            np.testing.assert_array_equal(b["element_number"], ids + 1)
            np.testing.assert_allclose(np.average(b["Bnorm_mT"], weights=w),
                                       float(coil["B_roi_mean_mT"]), rtol=2e-10)
            row = dict(coil=name, B_sphere30_mT=float(coil["B_sphere30_mT"]),
                       B_center_mT=float(coil["B_center_mT"]),
                       B_roi_mean_mT=float(coil["B_roi_mean_mT"]), E_roi_mean_V_m=ave,
                       E_roi_p95_V_m=weighted_quantile(e, w, .95),
                       E_reference_p99_V_m=p99, target_to_reference_pct=100 * ave / p99)
            if "D" in mesh.field:
                field = mesh.field["D"]
                if any(field is item for item in mesh.nodedata):
                    d = field.value[mesh.elm.node_number_list[ids] - 1].mean(axis=1)
                else:
                    d = field.value[ids]
                row["E_primary_roi_mean_V_m"] = np.average(np.linalg.norm(d, axis=1), weights=w)
            rows.append(row)
            extra = roi["regional_element_numbers"] - 1
            extra_w = roi["regional_volumes_mm3"]
            distance = roi["distances_mm"]
            extra_e = np.linalg.norm(efield[extra], axis=1)
            for radius in RADII_MM:
                sel = distance <= radius
                sensitivity.append(dict(coil=name, radius_mm=radius, elements=int(sel.sum()),
                                        volume_mm3=extra_w[sel].sum(),
                                        E_mean_V_m=np.average(extra_e[sel], weights=extra_w[sel])))
            for i, number in enumerate(ids + 1):
                targets.append(dict(coil=name, element_number=int(number), x_mm=centers[i, 0],
                                    y_mm=centers[i, 1], z_mm=centers[i, 2], volume_mm3=w[i],
                                    Bnorm_mT=b["Bnorm_mT"][i], Ex_V_m=efield[ids[i], 0],
                                    Ey_V_m=efield[ids[i], 1], Ez_V_m=efield[ids[i], 2],
                                    Enorm_V_m=e[i]))
        del mesh, efield
    write_rows("summary.csv", rows)
    write_rows("roi_radius.csv", sensitivity)
    write_rows("target_elements.csv", targets)
    print(json.dumps(rows, indent=2))


if __name__ == "__main__":
    main()
