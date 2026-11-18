"""Check the FEM excitation at the actual ROI nodes against analytical wire integrals."""
import csv
import json

import numpy as np

from define_target import CACHE, RESULTS, ROOT


def wire_dadt(nodes_mm, segments, transform, didt):
    start = segments[:, :3] @ transform[:3, :3].T + transform[:3, 3]
    stop = segments[:, 3:6] @ transform[:3, :3].T + transform[:3, 3]
    wire = stop - start
    lengths = np.linalg.norm(wire, axis=1)
    tangent = (wire / lengths[:, None]) * segments[:, 6:7]
    values = []
    for point in nodes_mm:
        offset = point - start
        axial = np.einsum("ij,ij->i", offset, wire) / lengths
        perpendicular = np.linalg.norm(offset - axial[:, None] * wire / lengths[:, None], axis=1)
        # Integral of dl/r along a straight segment. The common distance unit cancels.
        kernel = np.arcsinh(axial / perpendicular) - np.arcsinh((axial - lengths) / perpendicular)
        values.append(1e-7 * didt * (kernel @ tangent))
    return np.asarray(values)


if __name__ == "__main__":
    import simnibs
    placement = json.loads((RESULTS / "placement.json").read_text(encoding="utf-8"))
    rows = []
    with np.load(CACHE / "target_roi.npz") as roi:
        element_ids = roi["element_numbers"] - 1
    for name, coil in placement["coils"].items():
        mesh = simnibs.read_msh(str(CACHE / "fields" / f"{name}.msh"))
        node_ids, inverse = np.unique(mesh.elm.node_number_list[element_ids] - 1, return_inverse=True)
        segments = np.loadtxt(ROOT / coil["source"], skiprows=1, delimiter=",")
        reference = wire_dadt(mesh.nodes.node_coord[node_ids], segments,
                              np.asarray(coil["source_to_head_mm"]), 1e8)
        reference = reference[inverse.reshape(-1, 4)].mean(axis=1)
        field = mesh.field["D"]
        if any(field is item for item in mesh.nodedata):
            predicted = field.value[mesh.elm.node_number_list[element_ids] - 1].mean(axis=1)
        else:
            predicted = field.value[element_ids]
        delta = predicted - reference
        rel_l2 = 100 * np.linalg.norm(delta) / np.linalg.norm(reference)
        if not np.isfinite(rel_l2) or rel_l2 > 1.0:
            raise ValueError(f"Excitation at {name} ROI elements differs by {rel_l2}%")
        rows.append(dict(coil=name, target_nodes=len(node_ids), target_elements=len(element_ids),
                         dAdt_relative_vector_L2_pct=rel_l2,
                         dAdt_max_vector_difference_V_m=np.linalg.norm(delta, axis=1).max()))
        del mesh, field
    with (RESULTS / "source_check.csv").open("w", encoding="utf-8", newline="") as f:
        writer = csv.DictWriter(f, fieldnames=rows[0].keys())
        writer.writeheader()
        writer.writerows(rows)
    print(rows)
