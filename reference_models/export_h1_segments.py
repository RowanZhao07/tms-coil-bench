"""Decode the unmodified SimNIBS H1 .tcd file into SI line elements.

The upstream LineSegmentElements implementation uses points as element
MIDPOINTS, despite the brief web format text referring to starting positions.
No deformation, rotation, rescaling, closure repair, or decimation is applied.
Run with Python + NumPy + SciPy. Output files are derived data under GPL-3.0.
"""
from pathlib import Path
import base64
import hashlib
import json
import zlib
import numpy as np
from scipy.io import savemat

HERE = Path(__file__).resolve().parent
SOURCE = HERE / 'Brainsway_H1.tcd'


def decode(value):
    if isinstance(value, str):
        return np.frombuffer(zlib.decompress(base64.b64decode(value)), dtype=np.float64).reshape(-1, 3)
    return np.asarray(value, dtype=float)


def main():
    obj = json.loads(SOURCE.read_text(encoding='utf-8'))
    pp, vv, groups = [], [], []
    group_info = []
    for i, element in enumerate(obj['coilElementList']):
        assert element['type'] == 2
        p, v = decode(element['points']), decode(element['values'])
        keep = np.linalg.norm(v, axis=1) > 1e-12
        pp.append(p[keep] * 1e-3)
        vv.append(v[keep] * 1e-3)
        groups.append(np.full(int(keep.sum()), i, dtype=int))
        group_info.append({'group_id': i, 'name': element['name'], 'segments': int(keep.sum()), 'deformations': element.get('deformations', [])})
    mid = np.vstack(pp)
    dl = np.vstack(vv)
    group = np.concatenate(groups)
    start, end = mid - dl / 2, mid + dl / 2
    np.savetxt(HERE / 'h1_segments_si.csv', np.column_stack([group, start, end, mid, dl]),
               delimiter=',', header='group_id,x1_m,y1_m,z1_m,x2_m,y2_m,z2_m,xmid_m,ymid_m,zmid_m,dlx_m,dly_m,dlz_m', comments='')
    savemat(HERE / 'h1_segments_si.mat', {'mid_m': mid, 'dl_m': dl, 'start_m': start, 'end_m': end, 'group_id': group})
    metadata = {
        'source_url': 'https://raw.githubusercontent.com/simnibs/simnibs/fde509f6f8de38a4b03947944cf7ee93b8bebabe/simnibs/resources/coil_models/flexible_coils/Brainsway_H1.tcd',
        'source_sha256': hashlib.sha256(SOURCE.read_bytes()).hexdigest(),
        'retrieved_date': '2026-11-18', 'license': 'GPL-3.0; see SIMNIBS_LICENSE.txt',
        'model_version': obj.get('version'), 'coordinate_unit': 'm',
        'interpretation': 'midpoint and dl from official LineSegmentElements implementation; no deformations applied; original coordinate system retained',
        'segments': len(mid), 'groups': group_info,
        'bbox_mid_m': [mid.min(0).tolist(), mid.max(0).tolist()],
        'conductor_length_m': float(np.linalg.norm(dl, axis=1).sum()),
        'sum_dl_m': dl.sum(0).tolist(),
        'closure_note': 'Source line elements are not an exactly closed circuit: summed dl is nonzero. Retain this measured/source limitation; do not silently append a fictitious return conductor.',
        'references': ['https://simnibs.github.io/simnibs/build/html/documentation/coils.html', 'https://doi.org/10.1016/j.brs.2025.05.136'],
    }
    (HERE / 'h1_metadata.json').write_text(json.dumps(metadata, indent=2), encoding='utf-8')
    print(json.dumps({k: metadata[k] for k in ['segments', 'conductor_length_m', 'sum_dl_m', 'source_sha256']}, indent=2))


if __name__ == '__main__':
    main()
