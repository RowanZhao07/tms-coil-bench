"""Re-solve the stored physical source at a 1 V/m thalamus mean and audit the volume readout."""
import argparse
import csv
import gc
import json
from pathlib import Path
import sys

HERE = Path(__file__).resolve().parent
sys.path.insert(0, str(HERE.parents[1]))
from solve import cond_utils, fem, hash_file, mesh_io, np, solve_current
from analyze import weighted_quantile
from define_target import topology_hash

WORK = HERE.parents[2] / "tmp" / "brain_cohort"
SUBJECTS = ("sub01", "sub09", "sub10", "sub12", "sub15")


def solve_one(subject):
    prior = HERE.parent / "results" / subject
    audit = json.loads((prior / "solver_run.json").read_text(encoding="utf-8"))
    target = json.loads((prior / "target.json").read_text(encoding="utf-8"))
    head_file = WORK / "heads" / subject / f"{subject}.msh"
    roi_file = WORK / subject / "target_roi.npz"
    if hash_file(head_file) != audit["mesh_sha256"] or hash_file(roi_file) != audit["roi_file_sha256"]:
        raise ValueError(f"{subject}: head or region differs from the fixed-drive run")
    with np.load(roi_file) as roi:
        cortex_ids = roi["reference_element_numbers"] - 1
        cortex_volumes = roi["reference_volumes_mm3"]
        target_ids = roi["element_numbers"] - 1
        target_volumes = roi["volumes_mm3"]
    mesh = mesh_io.read_msh(str(head_file))
    cond = cond_utils.cond2elmdata(mesh, [c.value for c in cond_utils.standard_cond()])
    system = fem.TMSFEM(mesh, cond, solver_options="hypre", store_G=True)
    system.prepare_solver()
    free_dof = np.ones(system.dof_map.nr, dtype=bool)
    free_dof[system.dof_map[system.dirichlet.nodes]] = False
    rows = []
    for coil in ("h1", "double_cone", "figure8"):
        base = mesh_io.read_msh(str(WORK / subject / "fields" / f"{coil}.msh"))
        if topology_hash(base) != target["topology_sha256"] or not np.array_equal(base.nodes.node_coord, mesh.nodes.node_coord):
            raise ValueError(f"{subject}/{coil}: base field and head differ")
        fixed_E = base.field["E"].value
        target_mean = np.average(np.linalg.norm(fixed_E[target_ids], axis=1), weights=target_volumes)
        gain = 1.0 / target_mean
        reference = gain * fixed_E[cortex_ids]
        source = mesh_io.ElementData(gain * base.field["D"].value, name="D", mesh=mesh)
        del base, fixed_E
        gc.collect()
        potential, convergence = solve_current(system, free_dof, source)
        v = mesh_io.NodeData(potential, name="v", mesh=mesh)
        result = fem.calc_fields(v, "E", cond=cond, dadt=source)
        fields = result.field["E"].value
        target_read = np.average(np.linalg.norm(fields[target_ids], axis=1), weights=target_volumes)
        cortex = fields[cortex_ids]
        cortex_norm = np.linalg.norm(cortex, axis=1)
        p99 = weighted_quantile(cortex_norm, cortex_volumes, .99)
        relative = np.linalg.norm(cortex-reference) / np.linalg.norm(reference)
        reference_norm = np.linalg.norm(reference, axis=1)
        if abs(target_read-1.0) > 1e-7 or relative > 1e-6:
            raise ValueError(f"{subject}/{coil}: target or cortical vector readback failed")
        row = dict(subject=subject, coil=coil, target_mean_V_m=target_read, cortex_P99_V_m=p99,
                   didt_A_per_us=audit["didt_A_per_s"] * gain / 1e6,
                   cortex_vector_relative_L2=relative, **convergence)
        for threshold in (1, 15):
            actual = cortex_norm >= threshold
            row[f"cortex_ge_{threshold}_cm3"] = cortex_volumes[actual].sum() / 1000
            row[f"changed_elements_at_{threshold}"] = np.count_nonzero(actual != (reference_norm >= threshold))
        rows.append(row)
        print(f"{subject}/{coil}: target={target_read:.9f}, P99={p99:.5f}, V15={row['cortex_ge_15_cm3']:.6f}", flush=True)
        del result, fields, cortex, source, potential, v
        gc.collect()
    return rows


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("subject", choices=SUBJECTS)
    subject = parser.parse_args().subject
    values = solve_one(subject)
    folder = HERE / "results" / "verification"
    folder.mkdir(exist_ok=True, parents=True)
    out = folder / f"{subject}.csv"
    temporary = out.with_suffix(".tmp")
    with temporary.open("w", newline="", encoding="utf-8") as stream:
        writer = csv.DictWriter(stream, fieldnames=values[0], lineterminator="\n")
        writer.writeheader()
        writer.writerows(values)
    temporary.replace(out)
