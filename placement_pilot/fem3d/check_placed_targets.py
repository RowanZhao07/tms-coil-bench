"""Verify the nominal placed-coil field at the exact spherical-phantom targets."""
import csv
import json
import numpy as np
from run_double_cone import GEOMETRY, HERE, compute, write_csv

placement = HERE.parent / 'placement'
pose = json.loads((placement/'pose.json').read_text())['coils']['double_cone']['poses']['nominal']
matrix = np.array(pose['source_to_world'])
R, translation_mm = matrix[:3, :3], matrix[:3, 3]
assert np.allclose(R.T @ R, np.eye(3)) and np.linalg.det(R) > .99999
world_mm = np.asarray(pose['targets_world_mm'])
native_m = ((world_mm - translation_mm) @ R) * .001
segments = np.loadtxt(GEOMETRY, delimiter=',', skiprows=1)
# Verify the exact geometry being tested matches the independently exported posed winding.
posed = np.loadtxt(placement/'target_frame_segments/double_cone.csv', delimiter=',', skiprows=1)
for cols in [slice(0,3), slice(3,6)]:
    expected = (1000*segments[:,cols])@R.T + translation_mm
    assert np.allclose(expected, posed[:,cols], atol=1e-7, rtol=1e-12)
assert np.array_equal(segments[:, 6], posed[:, 6])
estimates = []
for h in [.002, .0015]:
    native_b = compute(segments, native_m, .480, h, 4)
    estimates.append((h, native_b @ R.T * 1e3))
with (placement/'target_fields.csv').open(encoding='utf-8-sig') as stream:
    refs = {(float(r['x_mm']), float(r['y_mm']),float(r['z_mm'])):
            np.array([float(r[c]) for c in ('Bx_mT','By_mT','Bz_mT')])
            for r in csv.DictReader(stream) if r['coil']=='double_cone' and r['case']=='nominal'}
result=[]
for h, world_b in estimates:
    for point, b in zip(world_mm,world_b):
        ref = refs[tuple(point)]
        err = np.linalg.norm(b-ref) / np.linalg.norm(ref)
        result.append(dict(h_mm=h*1000,x_mm=point[0],y_mm=point[1],z_mm=point[2],
                           FEM_Bx_mT=b[0],FEM_By_mT=b[1],FEM_Bz_mT=b[2],
                           FEM_norm_mT=np.linalg.norm(b),ref_Bx_mT=ref[0],ref_By_mT=ref[1],
                           ref_Bz_mT=ref[2],ref_norm_mT=np.linalg.norm(ref),vector_error_pct=err*100))
write_csv(HERE/'results'/'posed_double_cone.csv', result)
print(len(world_mm),'physical targets per mesh;',max(r['vector_error_pct'] for r in result if r['h_mm']==1.5),'% maximum discrepancy on final mesh.')
