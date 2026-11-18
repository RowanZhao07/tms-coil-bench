"""Independent Magpylib 6.3.0 evaluation of MATLAB-exported directed wires.

All spatial points in the MATLAB axis, xy and xz data are checked by default.
This is independent software using the same filament model, NOT FEM.
Run with .venv/Scripts/python.exe simulation/compare_magpylib.py
"""
from __future__ import annotations
import argparse
import csv
import json
from pathlib import Path
import platform
import time
import numpy as np
import magpylib as magpy
from scipy.constants import mu_0

ROOT = Path(__file__).resolve().parents[1]


def magpy_field(s: np.ndarray, p: np.ndarray, current: float = 1000.) -> np.ndarray:
    """Only source aggregation and batching are local; kernel is Magpylib."""
    out = np.zeros_like(p)
    for p0 in range(0, len(p), 128):
        obs = p[p0:p0 + 128]
        accum = np.zeros((len(obs), 3))
        for s0 in range(0, len(s), 512):
            seg = s[s0:s0 + 512]
            n, m = len(obs), len(seg)
            h = magpy.core.current_polyline_Hfield(
                observers=np.repeat(obs, m, axis=0),
                segments_start=np.tile(seg[:, :3], (n, 1)),
                segments_end=np.tile(seg[:, 3:6], (n, 1)),
                currents=np.tile(seg[:, 6] * current, n),
            )
            accum += h.reshape(n, m, 3).sum(axis=1)
        # Use the same explicitly declared engineering mu0 as MATLAB, so the
        # comparison isolates the two implementations rather than CODATA drift.
        out[p0:p0 + len(obs)] = accum * (4 * np.pi * 1e-7)
    return out


def run(stride: int = 1):
    out = ROOT / 'results' / 'magpylib'
    out.mkdir(parents=True, exist_ok=True)
    rows = []
    start = time.perf_counter()
    sources = sorted((ROOT / 'results' / 'geometry').glob('*.csv'))
    assert len(sources) == 14, f'Expected 14 geometries, found {len(sources)}'
    for source in sources:
        s = np.loadtxt(source, delimiter=',', skiprows=1)
        for section in ('axis', 'xy20', 'xz'):
            source_field = ROOT / 'results' / 'fields' / f'{source.stem}_{section}.csv'
            v = np.loadtxt(source_field, delimiter=',', skiprows=1)[::stride]
            p, matlab = v[:, :3], v[:, 3:6]
            t = time.perf_counter()
            b = magpy_field(s, p)
            err = b - matlab
            absolute = np.linalg.norm(err, axis=1)
            ref = np.linalg.norm(matlab, axis=1)
            # Avoid falsely huge relative errors at vector cancellation nodes.
            valid = ref > max(float(ref.max()) * 1e-8, 1e-15)
            row = dict(id=source.stem, section=section, points=len(p),
                       relative_L2=float(np.linalg.norm(err) / np.linalg.norm(matlab)),
                       max_absolute_T=float(absolute.max()),
                       max_relative_nonzero=float((absolute[valid] / ref[valid]).max()),
                       seconds=time.perf_counter() - t)
            rows.append(row)
            np.savetxt(out / f'{source.stem}_{section}.csv',
                       np.column_stack([p, b, np.linalg.norm(b, axis=1), absolute]),
                       delimiter=',', header='x_m,y_m,z_m,Bx_T,By_T,Bz_T,Bnorm_T,abs_vector_error_T', comments='')
            print(f'{source.stem:22} {section:5} n={len(p):5} L2={row["relative_L2"]:.3e} time={row["seconds"]:.1f}s', flush=True)
            if row['relative_L2'] >= 1e-8:
                raise RuntimeError(f'Independent field mismatch: {row}')
    with (out / 'comparison.csv').open('w', newline='', encoding='utf-8') as f:
        writer = csv.DictWriter(f, fieldnames=list(rows[0]))
        writer.writeheader(); writer.writerows(rows)
    metadata = dict(package='magpylib', version=magpy.__version__, python=platform.python_version(),
                    numpy=np.__version__, mu0_used=4*np.pi*1e-7, scipy_codata_mu0=mu_0,
                    stride=stride, geometry_count=14, total_points=sum(x['points'] for x in rows),
                    elapsed_seconds=time.perf_counter()-start,
                    comparison='Independent software implementation of the same directed-filament model; not FEM',
                    same_geometry_limitation='Agreement does not validate uncertain dimensions, source gaps, omitted leads, or clinical effects')
    (out / 'metadata.json').write_text(json.dumps(metadata, indent=2), encoding='utf-8')
    print(json.dumps(metadata, indent=2), flush=True)


if __name__ == '__main__':
    ap = argparse.ArgumentParser()
    ap.add_argument('--stride', type=int, default=1, help='1=check every saved observation point')
    args = ap.parse_args()
    if args.stride < 1:
        raise ValueError('stride must be positive')
    run(args.stride)
