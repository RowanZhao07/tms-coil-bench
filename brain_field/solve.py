"""Solve the three fixed-pose TMS fields on the unchanged Subject01 head mesh.

Run with the Python interpreter supplied by SimNIBS 4.6.0. Large .msh files
stay in tmp/brain_field; compact audit values are written beside the metrics.
"""
import gc
import hashlib
import importlib.metadata
import json
import os
import platform
import time
from pathlib import Path

# Set these before importing numerical libraries so the benchmark is reproducible
# on a laptop CPU. Run one right-hand side at a time.
for key in ("OMP_NUM_THREADS", "MKL_NUM_THREADS", "OPENBLAS_NUM_THREADS"):
    os.environ[key] = "8"

import fmm3dpy
import numpy as np
import simnibs
from simnibs.mesh_tools import mesh_io
from simnibs.simulation import fem
from simnibs.simulation.tms_coil.tms_coil import TmsCoil
from simnibs.simulation.tms_coil.tms_coil_element import LineSegmentElements
from simnibs.utils import cond_utils

HERE = Path(__file__).resolve().parent
WORK = HERE.parent / "tmp" / "brain_field"
RESULTS = HERE / "results"
DIDT = 1e8  # A/s
EPS = 1e-8


def hash_file(path):
    with path.open("rb") as stream:
        return hashlib.file_digest(stream, "sha256").hexdigest()


def working_set_peak_mib():
    if os.name == "nt":
        import ctypes
        from ctypes import wintypes
        class ProcessMemoryCounters(ctypes.Structure):
            _fields_ = [("cb", wintypes.DWORD), ("page_faults", wintypes.DWORD)] + [
                (name, ctypes.c_size_t) for name in ("peak_rss", "rss", "peak_paged",
                  "paged", "peak_nonpaged", "nonpaged", "pagefile", "peak_pagefile")]
        counter = ProcessMemoryCounters()
        counter.cb = ctypes.sizeof(counter)
        get_process = ctypes.windll.kernel32.GetCurrentProcess
        get_process.restype = wintypes.HANDLE
        get_info = ctypes.windll.psapi.GetProcessMemoryInfo
        get_info.argtypes = [wintypes.HANDLE, ctypes.POINTER(ProcessMemoryCounters), wintypes.DWORD]
        if not get_info(get_process(), ctypes.byref(counter), counter.cb):
            raise OSError("GetProcessMemoryInfo failed")
        return counter.peak_rss / 2**20
    import resource
    peak = resource.getrusage(resource.RUSAGE_SELF).ru_maxrss
    return peak / (2**20 if platform.system() == "Darwin" else 1024)


def primary(mesh, coil, affine):
    """SimNIBS midpoint line-current kernel, evaluated with its FMM3D library."""
    if not all(isinstance(el, LineSegmentElements) for el in coil.elements):
        raise TypeError("Expected only directed line elements")
    points = np.concatenate([el.get_points(affine) for el in coil.elements]) * 0.001
    directions = np.concatenate([el.get_values(affine) for el in coil.elements]) * 0.001
    nodes = mesh.nodes.node_coord
    fmm = fmm3dpy.lfmm3d(eps=EPS, sources=points.T, charges=directions.T,
                          targets=nodes.T * 0.001, nd=3, pgt=1)
    if fmm.ier != 0:
        raise RuntimeError(f"FMM3D status {fmm.ier}")
    values = (4e-7 * np.pi * DIDT * fmm.pottarg.T).copy()
    # Test against the native SimNIBS coil API and an explicit SI kernel at
    # reproducible head nodes, including nodes close to the deep target.
    rng = np.random.default_rng(7834)
    target = np.asarray(json.loads((RESULTS / "placement.json").read_text("utf-8"))["target_center_mm"])
    close = np.argsort(np.linalg.norm(nodes - target, axis=1))[:16]
    idx = np.unique(np.r_[close, rng.choice(len(nodes), size=64, replace=False)])
    reference = DIDT * sum(el.get_a_field(nodes[idx], affine) for el in coil.elements)
    direct = np.array([1e-7 * DIDT * (directions / np.linalg.norm(points - r * 0.001, axis=1)[:, None]).sum(0)
                       for r in nodes[idx]])
    def rel(p, q):
        return float(np.linalg.norm(p - q) / np.linalg.norm(q))
    audit = {"segments": len(points), "validation_nodes": len(idx), "fmm_requested_eps": EPS,
             "fmm_vs_native_rel_L2": rel(values[idx], reference),
             "fmm_vs_explicit_SI_rel_L2": rel(values[idx], direct),
             "fmm_vs_native_max_difference_V_per_m": float(np.linalg.norm(values[idx] - reference, axis=1).max())}
    if max(audit["fmm_vs_native_rel_L2"], audit["fmm_vs_explicit_SI_rel_L2"]) > 1e-6:
        raise ValueError(f"Unverified primary field: {audit}")
    return mesh_io.NodeData(values, name="D", mesh=mesh).node_data2elm_data(), audit


def main():
    start = time.perf_counter()
    head_path = WORK / "m2m_Subject01" / "Subject01.msh"
    spec = json.loads((HERE / "head_inputs.json").read_text(encoding="utf-8"))
    if hash_file(head_path) != spec["Subject01.msh"]["sha256"]:
        raise ValueError("Head checksum does not match head_inputs.json")
    mesh = mesh_io.read_msh(str(head_path))
    pose = json.loads((RESULTS / "placement.json").read_text(encoding="utf-8"))
    cond_list = cond_utils.standard_cond()
    cond = cond_utils.cond2elmdata(mesh, [c.value for c in cond_list])
    volume_tags = np.unique(mesh.elm.tag1[mesh.elm.elm_type == 4])
    print("Assembling the common TMS system", flush=True)
    system = fem.TMSFEM(mesh, cond, solver_options="hypre", store_G=True)
    system.prepare_solver()
    free_dof = np.ones(system.dof_map.nr, dtype=bool)
    free_dof[system.dof_map[system.dirichlet.nodes]] = False
    record = {"simnibs_version": importlib.metadata.version("simnibs"), "python_version": platform.python_version(),
              "solver": "SimNIBS TMSFEM, PETSc CG + HYPRE BoomerAMG", "platform": platform.system(),
              "threads_requested": 8, "parallel_coils": 1, "gpu_used": False,
              "runtime_packages": {key: importlib.metadata.version(key) for key in ("numpy", "scipy", "fmm3dpy", "petsc4py")},
              "ksp_relative_tolerance": float(system._solver.ksp.getTolerances()[0]),
              "D_definition": "tetrahedron mean of four nodal dA/dt vectors in V/m; primary E is -D",
              "mesh_sha256": spec["Subject01.msh"]["sha256"], "mesh_nodes": mesh.nodes.nr,
              "tetrahedra": int(np.count_nonzero(mesh.elm.elm_type == 4)), "didt_A_per_s": DIDT,
              "conductivity_S_per_m": {int(t): {"tissue": cond_list[t-1].name, "value": cond_list[t-1].value} for t in volume_tags},
              "ground_nodes_one_based": np.atleast_1d(system.dirichlet.nodes).tolist(),
              "common_setup_s": time.perf_counter()-start, "coils": {}}
    results_dir = WORK / "fields"
    results_dir.mkdir(parents=True, exist_ok=True)
    for key in ("h1", "figure8", "double_cone"):
        begin = time.perf_counter()
        affine = np.asarray(pose["coils"][key]["source_to_head_mm"])
        coil_file = WORK / "coils" / f"{key}.tcd"
        coil = TmsCoil.from_file(str(coil_file))
        print(f"Evaluating {key} coil and checking primary field", flush=True)
        dadt, audit = primary(mesh, coil, affine)
        primary_s = time.perf_counter() - begin
        rhs = system.assemble_rhs(dadt)
        potential = system.solve(rhs)
        ksp = system._solver.ksp
        status = int(ksp.getConvergedReason())
        residual = (system.A @ potential - rhs)[free_dof]
        normalized = float(np.linalg.norm(residual) / np.linalg.norm(rhs[free_dof]))
        if status <= 0 or normalized > 1e-6:
            raise RuntimeError(f"FEM failed: {key}, {status}, {normalized}")
        v = mesh_io.NodeData(potential, name="v", mesh=mesh)
        out = fem.calc_fields(v, "EeD", cond=cond, dadt=dadt)
        temp = results_dir / f"{key}.part.msh"
        mesh_io.write_msh(out, str(temp))
        os.replace(temp, results_dir / f"{key}.msh")
        record["coils"][key] = dict(audit, tcd_sha256=hash_file(coil_file), primary_seconds=primary_s,
            elapsed_seconds=time.perf_counter()-begin, free_dof_relative_residual=normalized,
            ksp_iterations=int(ksp.getIterationNumber()), ksp_converged_reason=status,
            process_peak_rss_MiB_so_far=working_set_peak_mib())
        print(key, record["coils"][key], flush=True)
        del out, dadt, v, rhs, potential
        gc.collect()
    record["process_peak_rss_MiB"] = working_set_peak_mib()
    record["elapsed_seconds"] = time.perf_counter()-start
    (RESULTS / "solver_run.json").write_text(json.dumps(record, indent=2)+"\n", encoding="utf-8")


if __name__ == "__main__":
    main()
