"""Compute both common-drive and cortical-P99-matched fields for a native MRI head.

Run with simnibs_python. The subject argument selects the previously prepared
target, placement and input mesh. Outputs for the common drive are full meshes;
the separate scalar-potential solves at matched drive are recorded numerically.
"""
import argparse
import csv
import json
import os
import sys
import time
from pathlib import Path

HERE = Path(__file__).resolve().parent
sys.path.insert(0, str(HERE.parent))
from solve import DIDT, fem, hash_file, mesh_io, np, solve_current, solve_subject
from analyze import weighted_quantile
from define_target import topology_hash


WORK = HERE.parents[1] / "tmp" / "brain_cohort"
CORTEX_SETPOINT_V_M = 30.0


def run(subject, head, expected_hash=None):
    cached = WORK / subject
    results = HERE / "results" / subject
    target = json.loads((results / "target.json").read_text(encoding="utf-8"))
    pose_file = results / "placement.json"
    pose = json.loads(pose_file.read_text(encoding="utf-8"))
    if pose["subject"] != subject or target["subject"] != subject or pose["mesh_sha256"] != target["mesh_sha256"]:
        raise ValueError("Placement and region masks belong to different subjects or meshes")
    if expected_hash is not None and expected_hash != target["mesh_sha256"]:
        raise ValueError("Target was not defined on the specified head mesh")
    with np.load(cached / "target_roi.npz", allow_pickle=False) as roi:
        target_ids = roi["element_numbers"] - 1
        cortex_ids = roi["reference_element_numbers"] - 1
        target_volumes = roi["volumes_mm3"]
        cortex_volumes = roi["reference_volumes_mm3"]
        expected_centers = roi["centroids_mm"]
    rows = []

    def metrics(vectors):
        magnitudes = np.linalg.norm(vectors[target_ids], axis=1)
        cortical = np.linalg.norm(vectors[cortex_ids], axis=1)
        return float(np.average(magnitudes, weights=target_volumes)), float(
            weighted_quantile(cortical, cortex_volumes, .99))

    def matched_drive(coil, mesh, cond, system, free_dof, dadt, base):
        if not rows:
            if topology_hash(mesh) != target["topology_sha256"]:
                raise ValueError("Target and solver head topology differ")
            actual_centers = mesh.nodes.node_coord[mesh.elm.node_number_list[target_ids] - 1].mean(axis=1)
            np.testing.assert_allclose(actual_centers, expected_centers, rtol=0, atol=2e-5)
            if not np.all(mesh.elm.tag1[cortex_ids] == 2):
                raise ValueError("Cortical reference points must belong to gray matter")
        base_roi, base_p99 = metrics(base.field["E"].value)
        if not np.isfinite(base_p99) or base_p99 <= 0:
            raise ValueError("Cortical P99 must be positive and finite")
        gain = CORTEX_SETPOINT_V_M / base_p99
        start = time.perf_counter()
        # Build a new right-hand side from the scaled physical source and solve again.
        current = mesh_io.ElementData(gain * dadt.value, name="D", mesh=mesh)
        potential, convergence = solve_current(system, free_dof, current)
        v = mesh_io.NodeData(potential, name="v", mesh=mesh)
        rescaled = fem.calc_fields(v, "E", cond=cond, dadt=current)
        result_roi, result_p99 = metrics(rescaled.field["E"].value)
        # Verify directions and amplitudes on the target, not just the pooled norms.
        expected = base.field["E"].value[target_ids] * gain
        error = np.linalg.norm(rescaled.field["E"].value[target_ids] - expected) / np.linalg.norm(expected)
        if not np.isfinite(error) or error > 1e-6:
            raise RuntimeError(f"Matched-drive ROI vectors differ from linear scaling: {error}")
        np.testing.assert_allclose(result_p99, CORTEX_SETPOINT_V_M, rtol=1e-7)
        np.testing.assert_allclose(result_roi, gain * base_roi, rtol=1e-6)
        rows.append(dict(coil=coil, unit_didt_A_per_s=DIDT, fixed_roi_mean_V_m=base_roi,
                         fixed_cortex_p99_V_m=base_p99, gain=gain, matched_didt_A_per_s=DIDT*gain,
                         requested_cortex_p99_V_m=CORTEX_SETPOINT_V_M,
                         matched_p99_V_m=result_p99, matched_roi_mean_V_m=result_roi,
                         roi_vector_scaling_rel_L2=float(error), **convergence,
                         matched_solve_seconds=time.perf_counter()-start))
        print(f"{subject} {coil}: matched P99 {result_p99:.10f}, ROI {result_roi:.7g} V/m", flush=True)

    record = solve_subject(head, pose_file, WORK / "coils",
                           cached / "fields", target["mesh_sha256"], second_drive=matched_drive)
    record.update(subject=subject, paired_drive_CortexP99_V_m=CORTEX_SETPOINT_V_M,
                  roi_file_sha256=hash_file(cached / "target_roi.npz"), peak_includes_both_drives=True)
    table = results / "matched_check.csv"
    temporary = table.with_suffix(".tmp")
    with temporary.open("w", newline="", encoding="utf-8") as stream:
        writer = csv.DictWriter(stream, fieldnames=rows[0].keys(), lineterminator="\n")
        writer.writeheader()
        writer.writerows(rows)
    os.replace(temporary, table)
    audit = results / "solver_run.json"
    temporary = audit.with_suffix(".tmp")
    temporary.write_text(json.dumps(record, indent=2) + "\n", encoding="utf-8")
    os.replace(temporary, audit)


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("subject")
    parser.add_argument("--mesh", type=Path, help="Override the native mesh location")
    parser.add_argument("--sha256", help="Optional independently obtained native mesh checksum")
    args = parser.parse_args()
    mesh = args.mesh if args.mesh else WORK / "heads" / args.subject / (args.subject + ".msh")
    run(args.subject, mesh, args.sha256)
