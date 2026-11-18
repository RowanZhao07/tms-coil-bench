"""Three-dimensional Galerkin verification of the study's double-cone coil."""
from __future__ import annotations

import argparse
import csv
import hashlib
import itertools
import json
import math
import platform
import time
from pathlib import Path

import numpy as np
import scipy
from scipy.fft import dstn, idstn
from scipy.sparse import diags, kron

HERE = Path(__file__).resolve().parent
ROOT = HERE.parents[1]
GEOMETRY = ROOT / "results" / "geometry" / "double_cone.csv"
MU0 = 4e-7 * np.pi
CURRENT = 1000.0
CORNERS = tuple(itertools.product((0, 1), repeat=3))
CASES = [(0.192, .004), (.256, .008), (.256, .004), (.256, .002),
         (.320, .004), (.320, .002), (.480, .004), (.480, .002), (.480, .0015)]


def write_csv(path, records):
    with path.open("w", encoding="utf-8", newline="") as stream:
        writer = csv.DictWriter(stream, fieldnames=list(records[0]))
        writer.writeheader()
        writer.writerows(records)


def eigenvalues(n, h):
    k = np.arange(1, n)
    v = np.cos(np.pi * k / n)
    stiffness = (2 - 2 * v) / h
    mass = (4 + 2 * v) * h / 6
    return (stiffness[:, None, None] * mass[None, :, None] * mass[None, None, :]
            + mass[:, None, None] * stiffness[None, :, None] * mass[None, None, :]
            + mass[:, None, None] * mass[None, :, None] * stiffness[None, None, :])


def solve(rhs, denominator, workers):
    """Exact algebraic solution of the tensor-product Q1 stiffness system."""
    coeffs = dstn(rhs, type=1, norm="ortho", workers=workers, overwrite_x=True)
    coeffs /= denominator
    return idstn(coeffs, type=1, norm="ortho", workers=workers, overwrite_x=True)


def check_galerkin_solver():
    # Independent explicit sparse FE matrix on a small grid checks spectral assembly.
    n, h = 10, .01
    stiff = diags([-np.ones(n-2), 2*np.ones(n-1), -np.ones(n-2)], [-1, 0, 1]) / h
    mass = diags([np.ones(n-2), 4*np.ones(n-1), np.ones(n-2)], [-1, 0, 1]) * (h/6)
    mat = kron(kron(stiff, mass), mass) + kron(kron(mass, stiff), mass) + kron(kron(mass, mass), stiff)
    rhs = np.random.default_rng(41).normal(size=(n-1, n-1, n-1))
    sol = solve(rhs.copy(), eigenvalues(n, h), 1)
    residual = np.linalg.norm(mat @ sol.ravel() - rhs.ravel()) / np.linalg.norm(rhs)
    assert residual < 1e-12
    return float(residual)


def integrate_line_load(segments, extent, h):
    """Distribute each straight segment to Q1 basis functions by exact cellwise quadrature."""
    locations, vector_weights = [], []
    for segment in segments:
        start, dl = segment[:3], segment[3:6] - segment[:3]
        a = (start + extent) / h
        b = (start + dl + extent) / h
        breaks = [0., 1.]
        for j in range(3):
            if abs(b[j] - a[j]) < 1e-13:
                continue
            lo, hi = sorted((a[j], b[j]))
            for plane in range(math.floor(lo) + 1, math.ceil(hi)):
                t = (plane - a[j]) / (b[j] - a[j])
                if 1e-12 < t < 1 - 1e-12:
                    breaks.append(t)
        breaks = sorted(set(breaks))
        for left, right in zip(breaks[:-1], breaks[1:]):
            mid = (left + right) / 2
            half = (right - left) / 2
            for sgn in (-1, 1):
                t = mid + sgn * half / math.sqrt(3)
                locations.append(a + t * (b - a))
                vector_weights.append(half * dl * segment[6])
    q = np.asarray(locations)
    weights = np.asarray(vector_weights)
    base = np.floor(q).astype(int)
    local = q - base
    idx, vals = [], []
    for corner in CORNERS:
        off = np.array(corner)
        phi = np.prod(np.where(off, local, 1-local), axis=1)
        idx.append(base+off-1)  # eliminate homogeneous boundary; source is strictly interior
        vals.append((MU0 * CURRENT * phi[:, None]) * weights)
    idx = np.vstack(idx)
    vals = np.vstack(vals)
    assert np.max(np.abs(vals.sum(axis=0))) < 1e-13  # weighted dl of the closed coil
    n = round(2 * extent / h)
    assert np.min(idx) >= 0 and np.max(idx) < n-1
    return tuple(idx.T), vals


def recovered_gradient(a, p, extent, h):
    """Mean adjacent Q1 gradients at grid vertices, then trilinear interpolation."""
    grid = (p + extent) / h
    ijk = np.floor(grid).astype(int)
    frac = grid - ijk
    grad = np.zeros_like(p)
    for corner in CORNERS:
        off = np.asarray(corner)
        node = ijk + off - 1
        w = np.prod(np.where(off, frac, 1-frac), axis=1)
        for dim in range(3):
            step = np.eye(3, dtype=int)[dim]
            grad[:, dim] += w * (a[tuple((node+step).T)] - a[tuple((node-step).T)]) / (2*h)
    return grad


def roi():
    xy = np.arange(-20, 21) * .0045
    points = [(x, y, .02) for x in xy for y in xy]
    points += [(x, 0, z/1000) for x in xy for z in range(10, 101, 2)]
    points += [(0, 0, z/1000) for z in range(10, 101)]
    return np.unique(np.round(np.array(points), 8), axis=0)


def compute(segments, p, extent, h, workers):
    n = round(2*extent/h)
    idx, v = integrate_line_load(segments, extent, h)
    lam = eigenvalues(n, h)
    grad_a = np.empty((len(p), 3, 3))
    for j in range(3):
        rhs = np.zeros((n-1, n-1, n-1))
        np.add.at(rhs, idx, v[:, j])
        a = solve(rhs, lam, workers)
        grad_a[:, j, :] = recovered_gradient(a, p, extent, h)
        del a, rhs
    return np.stack((grad_a[:, 2, 1] - grad_a[:, 1, 2],
                     grad_a[:, 0, 2] - grad_a[:, 2, 0],
                     grad_a[:, 1, 0] - grad_a[:, 0, 1]), axis=1)


def reference(p):
    # Used ONLY AFTER all PDE solutions have been obtained. No reference on boundary or RHS.
    data = {}
    for suffix in ("axis", "xy20", "xz"):
        table = np.loadtxt(ROOT / "results" / "fields" / f"double_cone_{suffix}.csv",
                           delimiter=",", skiprows=1)
        for row in table:
            data[tuple(np.round(row[:3], 8))] = row[3:6]
    return np.array([data[tuple(q)] for q in p])


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("--workers", type=int, default=4)
    args = parser.parse_args()
    result = HERE / "results"
    result.mkdir(exist_ok=True)
    sparse_resid = check_galerkin_solver()
    seg = np.loadtxt(GEOMETRY, delimiter=",", skiprows=1)
    p = roi()
    estimates = []
    print("Small assembled-system residual", sparse_resid, flush=True)
    for extent, h in CASES:
        timer = time.perf_counter()
        field = compute(seg, p, extent, h, args.workers)
        elapsed = time.perf_counter() - timer
        estimates.append((extent, h, elapsed, field))
        print(f"3D Q1 half-width={extent*1000:.0f}mm h={h*1000:.0f}mm done {elapsed:.2f}s", flush=True)

    ref = reference(p)
    assert np.isfinite(ref).all()
    rows = []
    for extent, h, elapsed, b in estimates:
        delta = b - ref
        rel = np.linalg.norm(delta, axis=1) / np.linalg.norm(ref, axis=1)
        n = round(2*extent/h)
        rows.append(dict(half_width_mm=extent*1000, h_mm=h*1000, cells=n**3,
                         interior_dof_per_component=(n-1)**3, seconds=elapsed,
                         relative_vector_l2_pct=100*np.linalg.norm(delta)/np.linalg.norm(ref),
                         median_pointwise_vector_error_pct=100*np.median(rel),
                         p95_pointwise_vector_error_pct=100*np.percentile(rel, 95),
                         max_pointwise_vector_error_pct=100*np.max(rel),
                         max_vector_difference_mT=1000*np.max(np.linalg.norm(delta, axis=1))))
    write_csv(result / "mesh_domain_convergence.csv", rows)
    b = estimates[-1][-1]
    rel = np.linalg.norm(b-ref, axis=1) / np.linalg.norm(ref, axis=1)
    names = "x_m,y_m,z_m,FEM_Bx_T,FEM_By_T,FEM_Bz_T,MATLAB_Bx_T,MATLAB_By_T,MATLAB_Bz_T,vector_difference_T,relative_vector_error_pct"
    table = np.column_stack((p, b, ref, np.linalg.norm(b-ref, axis=1), 100*rel))
    np.savetxt(result / "final_roi.csv", table, delimiter=",", header=names, comments="")
    # Intentionally include bilateral and noncentral probes on both planes.
    probes = np.array([(0,0,z) for z in (10,20,30,40,50)]
                      + [(-45,-45,20),(45,45,20),(27,-18,20),(-27,18,20),
                         (45,0,30),(-45,0,50),(0,45,20),(18,0,70)], dtype=float)
    take = [np.flatnonzero(np.all(np.isclose(p,q/1000, atol=1e-12, rtol=0), axis=1))[0] for q in probes]
    np.savetxt(result / "selected_probes.csv", table[take], delimiter=",", header=names, comments="")
    metadata = dict(source_segments=len(seg), current_A=CURRENT, mu0_H_per_m=MU0,
                    source_sha256=hashlib.sha256(GEOMETRY.read_bytes()).hexdigest(),
                    field_reference="Existing MATLAB exact polygonal-wire results; read only after FEM solves",
                    boundary="All components of A=0 on six cube faces",
                    source="Closed filament line measure, cellwise exact two-point Gauss on Q1 test functions",
                    equation="-Laplacian(A) = mu0 J, homogeneous permeability, divergence-free closed source",
                    finite_element="Tensor-product trilinear conforming hexahedra, consistent Q1 matrices",
                    curl_recovery="Average nodal gradient and trilinear interpolation from nodal values",
                    residual_vs_assembled_Q1=sparse_resid, n_roi=len(p),
                    peak_process_working_set_bytes=None, python=platform.python_version(), numpy=np.__version__,
                    scipy=scipy.__version__, fft_workers=args.workers)
    try:
        import psutil
        metadata["peak_process_working_set_bytes"] = getattr(psutil.Process().memory_info(), "peak_wset", None)
    except ImportError:
        pass
    (result / "metadata.json").write_text(json.dumps(metadata, indent=2), encoding="utf-8")
    print(json.dumps(rows, indent=2), flush=True)

if __name__ == "__main__":
    main()
