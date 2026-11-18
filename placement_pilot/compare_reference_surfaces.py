"""Compare two clearance conventions for exactly the same oriented coil.

Input is the *placed* segment export from placement/. The only allowed change
is a translation along world +z. All inputs are mm and magnetic output is mT.
The plane is a flat geometric phantom, not an estimate of anatomical scalp.
"""
from pathlib import Path
import csv
import json
import numpy as np
import magpylib

ROOT = Path(__file__).resolve().parent
GAP_MM = 12.0
RADIUS_MM = 85.0
DEPTHS_MM = (10, 30, 50)


def field_mT(segments, point_mm):
    starts, ends = segments[:, :3] / 1000, segments[:, 3:6] / 1000
    observers = np.broadcast_to(np.array(point_mm) / 1000, starts.shape).copy()
    h = magpylib.core.current_polyline_Hfield(observers, starts, ends, 1000 * segments[:, 6])
    return h.sum(axis=0) * (1000 * 4 * np.pi * 1e-7)


def sphere_gap(segments):
    center = np.array([0.0, 0.0, -RADIUS_MM])
    a, b = segments[:, :3], segments[:, 3:6]
    direction = b - a
    along = np.einsum('ij,ij->i', center - a, direction) / np.einsum('ij,ij->i', direction, direction)
    nearest = a + np.clip(along, 0, 1)[:, None] * direction
    return np.linalg.norm(nearest - center, axis=1).min() - RADIUS_MM


def main():
    rows = []
    files = sorted((ROOT / 'placement' / 'target_frame_segments').glob('*.csv'))
    if len(files) != 3:
        raise ValueError('Run the three-coil spherical placement before this comparison.')
    with (ROOT / 'placement' / 'target_fields.csv').open(encoding='utf-8-sig') as fp:
        nominal = {
            (row['coil'], -float(row['z_mm'])): np.array([float(row[k]) for k in ('Bx_mT', 'By_mT', 'Bz_mT')])
            for row in csv.DictReader(fp)
            if row['case'] == 'nominal' and float(row['x_mm']) == float(row['y_mm']) == 0.0
        }
    for file in files:
        s = np.loadtxt(file, delimiter=',', skiprows=1, ndmin=2)
        if s.shape[1] != 7 or not np.isfinite(s).all():
            raise ValueError(f'Expected six mm endpoint coordinates and one current weight in {file}')
        sphere_gap_mm = sphere_gap(s)
        if abs(sphere_gap_mm-GAP_MM) > 0.01:
            raise ValueError(f'Incompatible baseline gap for {file.stem}: {sphere_gap_mm}')
        flat = s.copy()
        dz = GAP_MM - s[:, [2, 5]].min()
        flat[:, [2, 5]] += dz
        for depth in DEPTHS_MM:
            target = [0, 0, -depth]
            bs = field_mT(s, target)
            np.testing.assert_allclose(bs, nominal[file.stem, depth], rtol=1e-8, atol=1e-8)
            bf = field_mT(flat, target)
            ns, nf = np.linalg.norm(bs), np.linalg.norm(bf)
            rows.append(dict(coil=file.stem, depth_mm=depth, sphere_radius_mm=RADIUS_MM,
                             centerline_gap_mm=GAP_MM, shift_to_flat_mm=dz,
                             sphere_Bx_mT=bs[0],sphere_By_mT=bs[1],sphere_Bz_mT=bs[2],sphere_norm_mT=ns,
                             flat_Bx_mT=bf[0],flat_By_mT=bf[1],flat_Bz_mT=bf[2],flat_norm_mT=nf,
                             flat_minus_sphere_percent=100*(nf/ns-1)))
    with (ROOT/'reference_surface_comparison.csv').open('w', newline='', encoding='utf-8-sig') as fp:
        writer=csv.DictWriter(fp, fieldnames=list(rows[0])); writer.writeheader(); writer.writerows(rows)
    print(json.dumps([r for r in rows if r['depth_mm']==30],indent=2))


if __name__ == '__main__':
    main()
