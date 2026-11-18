"""Rigid registration of three coil sources to a spherical clearance reference.

World origin is the scalp surrogate apex and +z points outwards. All geometric
quantities are in millimetres. Field code converts source/target data to SI.
Run with the coil_study project Python environment (numpy / scipy).
"""
from pathlib import Path
import base64
import csv
import hashlib
import json
import zlib

import numpy as np
from scipy.io import loadmat
from scipy.optimize import brentq, least_squares, minimize

HERE = Path(__file__).resolve().parent
PROJECT = HERE.parents[1]
CURRENT = 1000.0
RADIUS = 85.0
GAP = 12.0
AXIS_TARGETS = np.array([[0., 0., -d] for d in (10, 30, 50)])
TARGETS = np.vstack((AXIS_TARGETS, [[10, 0, -30], [-10, 0, -30],
                                    [0, 10, -30], [0, -10, -30]]))
HEADER = 'ax_mm,ay_mm,az_mm,bx_mm,by_mm,bz_mm,weight'


def decode(field, dtype=np.float64):
    if isinstance(field, str):
        return np.frombuffer(zlib.decompress(base64.b64decode(field)), dtype=dtype).reshape(-1, 3)
    return np.asarray(field, dtype=dtype)


def save_table(name, rows):
    with (HERE/name).open('w', newline='', encoding='utf-8') as fp:
        w = csv.DictWriter(fp, fieldnames=list(rows[0]))
        w.writeheader()
        w.writerows(rows)


def transform_segments(s, r, t):
    return np.column_stack((s[:, :3] @ r.T + t, s[:, 3:6] @ r.T + t, s[:, 6]))


def transform_matrix(r, t):
    matrix = np.eye(4)
    matrix[:3, :3] = r
    matrix[:3, 3] = t
    return matrix


def exact_closest_segment(a, b, centre):
    """Closest point to centre on each finite segment, including its endpoints."""
    ab = b - a
    t = np.einsum('ij,ij->i', centre-a, ab) / np.einsum('ij,ij->i', ab, ab)
    return a + np.clip(t, 0, 1)[:, None] * ab


def closest_shell_distance(triangles, centre):
    """Minimum distance from centre to triangles, including faces, edges and vertices."""
    a,b,c = triangles.transpose(1,0,2)
    candidates = [np.linalg.norm(exact_closest_segment(x,y,centre)-centre,axis=1).min()
                  for x,y in [(a,b),(b,c),(c,a)]]
    ab,ac = b-a,c-a
    normal=np.cross(ab,ac)
    normal2=np.einsum('ij,ij->i',normal,normal)
    good=normal2>1e-18
    a,b,c,normal,normal2=a[good],b[good],c[good],normal[good],normal2[good]
    signed=np.einsum('ij,ij->i',centre-a,normal)
    proj=centre-(signed/normal2)[:,None]*normal
    inside=np.ones(len(a),dtype=bool)
    for x,y in [(a,b),(b,c),(c,a)]:
        inside &= np.einsum('ij,ij->i',np.cross(y-x,proj-x),normal)>=-1e-9
    if inside.any():
        candidates.append(np.min(abs(signed[inside])/np.sqrt(normal2[inside])))
    return float(min(candidates))


def h1_source():
    folder=PROJECT/'reference_models'
    source=folder/'Brainsway_H1.tcd'
    tcd=json.loads(source.read_text(encoding='utf-8'))
    groups=tcd['coilModels']
    used=[e['elementCasing'] for e in tcd['coilElementList']]
    assert len(set(used))==len(groups)==18
    designated=np.vstack([decode(groups[i]['minDistancePoints']) for i in used])
    def residual(v):
        return np.linalg.norm(designated-v[:3],axis=1)-v[3]
    fit=least_squares(residual,[0,-60,50,90],loss='soft_l1',f_scale=5)
    fixed=next(e for e in tcd['coilElementList'] if e['name']=='H1_no_def')
    mark=decode(groups[fixed['elementCasing']]['minDistancePoints'])
    local_seed=mark[np.argmin(np.linalg.norm(mark,axis=1))]
    patch=mark[np.linalg.norm(mark-local_seed,axis=1)<=12]
    anchor=patch.mean(0)
    z_axis=anchor-fit.x[:3]
    z_axis/=np.linalg.norm(z_axis)
    # Lock azimuth using the original +x axis projected into the tangent plane.
    x_axis=np.array([1.,0,0])-z_axis*z_axis[0]
    x_axis/=np.linalg.norm(x_axis)
    r=np.vstack((x_axis,np.cross(z_axis,x_axis),z_axis))
    h=loadmat(folder/'h1_segments_si.mat')
    a,b=h['start_m']*1000,h['end_m']*1000
    s=np.column_stack((a,b,np.ones(len(a))))
    fit_rotation=r.copy()
    a0=(a-anchor)@r.T
    b0=(b-anchor)@r.T
    # A second, explicitly geometry-only registration sets the tilt of this
    # asymmetric rigid helmet. The figure-eight and double-cone have specified
    # symmetry/hinge axes. Constrain possible head centres to a 30-degree cone
    # about the marked working face. The target field never enters the fit.
    def constraint(c):
        near=exact_closest_segment(a0,b0,c)
        lengths=np.linalg.norm(c-near,axis=1)
        k=lengths.argmin()
        return lengths[k]-RADIUS-GAP,(c-near[k])/lengths[k]
    constraints=[dict(type='ineq',fun=lambda c:constraint(c)[0],jac=lambda c:constraint(c)[1]),
                 dict(type='ineq',fun=lambda c:-c[2]*np.tan(np.deg2rad(30))-np.linalg.norm(c[:2]))]
    solutions=[]
    for initial in [[0,0,-190],[-15,0,-110],[20,0,-108],[0,-15,-110],[30,-30,-100]]:
        estimate=minimize(np.linalg.norm,initial,method='SLSQP',constraints=constraints,
                          options={'ftol':1e-9,'maxiter':120})
        if estimate.success and constraint(estimate.x)[0]>-1e-6:
            solutions.append(estimate)
    assert solutions
    c_local=min(solutions,key=lambda x:x.fun).x
    z_axis=-(c_local@r)/np.linalg.norm(c_local)
    x_axis=np.array([1.,0,0])-z_axis*z_axis[0]
    x_axis/=np.linalg.norm(x_axis)
    r=np.array([x_axis,np.cross(z_axis,x_axis),z_axis])
    shells=[]
    for i in used:
        m=groups[i]
        pts=decode(m['points']); faces=decode(m['faces'],np.int64)
        assert faces.min()>=0 and faces.max()<len(pts)
        shells.append(pts[faces])
    metadata=dict(anchor_mm=anchor.tolist(), designated_points=len(designated),
                  anchor_patch_size=len(patch), fitted_sphere_centre_source_mm=fit.x[:3].tolist(),
                  fitted_sphere_radius_source_mm=fit.x[3], shell_triangles=sum(map(len,shells)),
                  fit_median_abs_residual_mm=float(np.median(np.abs(fit.fun))),
                  fit_p90_abs_residual_mm=float(np.percentile(np.abs(fit.fun),90)),
                  tilt_from_initial_marked_fit_deg=float(np.rad2deg(np.arccos(np.clip(r[2]@fit_rotation[2],-1.,1.)))),
                  geometric_multistart_converged=len(solutions),
                  initial_marked_fit_rotation=fit_rotation.tolist(),
                  source_sha256=hashlib.sha256(source.read_bytes()).hexdigest(),
                  reference_rule='Fixed H1_no_def component: average marked minDistancePoints within 12 mm of marked point closest to native origin. Initialize normal from robust fit of all 18 sets of marks. Refine its direction within 30 degrees to minimize geometric fiducial-to-apex distance under full 12-mm wire clearance at radius 85 mm. Fix that direction for all sensitivity cases. Five geometric starts; local result, not proof of global optimum; no field or clinical target used in optimization.')
    return s,r,anchor,np.vstack(shells),metadata


def figure8():
    turns=[]
    theta=np.linspace(0,2*np.pi,769)
    for side in [-1,1]:
        for radius in [4,12,20,28,36,44]:
            p=np.column_stack((side*44+radius*np.cos(theta),radius*np.sin(theta),0*theta))
            turns.append(np.column_stack((p[:-1],p[1:],np.full(768,side))))
    return np.vstack(turns)


def evaluate_b(s, observers):
    """Exact straight-filament angular factor; 1000 A in each unit-weight turn."""
    a,b,w=s[:,:3]*1e-3,s[:,3:6]*1e-3,s[:,6]
    dl=b-a; ell=np.linalg.norm(dl,axis=1); u=dl/ell[:,None]
    out=[]
    for p in np.atleast_2d(observers)*1e-3:
        d=p-a
        z=np.einsum('ij,ij->i',d,u)
        perp=d-z[:,None]*u
        h2=np.einsum('ij,ij->i',perp,perp)
        axial=(h2<1e-23)
        assert not np.any(axial & (z>=0)&(z<=ell))
        factor=z/np.sqrt(h2+z*z)-(z-ell)/np.sqrt(h2+(z-ell)**2)
        factor=w*factor/np.where(axial,1.,h2)
        factor[axial]=0
        out.append(1e-7*CURRENT*np.sum(np.cross(u,perp)*factor[:,None],axis=0))
    return np.array(out)


def register(s, r, anchor, radius, gap, shift_x=0., pitch_deg=0.):
    centre=np.array([0.,0.,-radius])
    ang=np.deg2rad(pitch_deg)
    tilt=np.array([[np.cos(ang),0,np.sin(ang)],[0,1,0],[-np.sin(ang),0,np.cos(ang)]])
    r=tilt@r
    t0=-(r@anchor)+np.array([shift_x,0,0])
    before=transform_segments(s,r,t0)
    a,b=before[:,:3],before[:,3:6]
    def clearance(shift):
        c=centre-np.array([0.,0.,shift])
        nearest=exact_closest_segment(a,b,c)
        return float(np.min(np.linalg.norm(nearest-c,axis=1))-radius-gap)
    # Start above the sphere and approach along -z to the outside contact
    # branch. Find its first crossing rather than solving an interior root.
    maxshift=max(60.,gap+40.-min(a[:,2].min(),b[:,2].min()))
    upper=maxshift
    for lower in np.arange(upper-2,-maxshift,-2):
        if clearance(lower)<=0:
            break
        upper=lower
    else:
        raise ValueError('Unable to bracket outside contact')
    shift=brentq(clearance,lower,upper,xtol=1e-11)
    t=t0+np.array([0.,0.,shift]); world=transform_segments(s,r,t)
    a,b=world[:,:3],world[:,3:6]
    near=exact_closest_segment(a,b,centre)
    clearance_by_segment=np.linalg.norm(near-centre,axis=1)-radius
    ix=clearance_by_segment.argmin(); cpoint=near[ix]
    assert np.isclose(clearance_by_segment[ix],gap,atol=1e-7)
    assert np.isclose(np.linalg.det(r),1.) and np.allclose(r@r.T,np.eye(3))
    # Independent dense subdivisions bound the exact min from above.
    sample_min=min(np.linalg.norm(a+(b-a)*k/8-centre,axis=1).min()
                   for k in range(9))-radius
    assert sample_min>=gap-1e-7 and sample_min-gap<0.025
    return world,r,t,dict(min_wire_gap_mm=clearance_by_segment[ix],
                       sampled_wire_gap_mm=sample_min,
                       effective_anchor_z_mm=float((r@anchor+t)[2]),
                       contact_x_mm=float(cpoint[0]),contact_y_mm=float(cpoint[1]),
                       contact_z_mm=float(cpoint[2]),contact_xy_offset_mm=float(np.linalg.norm(cpoint[:2])),
                       closest_wire_to_apex_mm=float(np.linalg.norm(exact_closest_segment(a,b,np.zeros(3)),axis=1).min()))


def main():
    h=h1_source()
    cone=np.loadtxt(PROJECT/'results/geometry/double_cone.csv',delimiter=',',skiprows=1)
    cone[:,:6]*=1000
    sources={'h1':h,'figure8':(figure8(),np.eye(3),np.zeros(3),None,
                    {'reference_rule':'Origin is the tangent junction of two outer loops; native y follows the shared tangent.'}),
             'double_cone':(cone,np.eye(3),np.zeros(3),None,
                    {'reference_rule':'Origin is central symmetry/hinge axis at project z=0 upper-endpoint plane. Wings descend around sphere, project y follows hinge.'})}
    for sub in ['source_geometry','target_frame_segments']:
        (HERE/sub).mkdir(exist_ok=True)
    pose_rows=[];field_rows=[];placement={}
    source_checks=[]
    cases=[('nominal',85,12,0,0),('radius80',80,12,0,0),('radius90',90,12,0,0),
           ('gap10',85,10,0,0),('gap15',85,15,0,0),('lateral_minus5',85,12,-5,0),
           ('lateral_plus5',85,12,5,0),('tilt_minus5',85,12,0,-5),('tilt_plus5',85,12,0,5)]
    for key,(s,rot,anchor,shells,metadata) in sources.items():
        np.savetxt(HERE/'source_geometry'/f'{key}.csv',s,delimiter=',',header=HEADER,comments='')
        if key in ('figure8','double_cone'):
            path=(PROJECT/'results/validation/document_foe_case.csv' if key=='figure8' else
                  PROJECT/'results/fields/double_cone_axis.csv')
            v=np.loadtxt(path,delimiter=',',skiprows=1)
            b=evaluate_b(s,1000*v[:,:3])
            err=float(np.linalg.norm(b-v[:,3:6])/np.linalg.norm(v[:,3:6]))
            assert err<1e-10
            source_checks.append(dict(coil=key,source_reference=str(path.relative_to(PROJECT)),
                                      points=len(v),relative_vector_L2=err))
        placement[key]={**metadata,'poses':{}}
        coil_cases=list(cases)
        if key=='h1':
            coil_cases.append(('initial_marked_normal',85,12,0,0))
        for case,radius,gap,dx,tilt in coil_cases:
            selected_rotation=(np.asarray(metadata['initial_marked_fit_rotation'])
                               if case=='initial_marked_normal' else rot)
            world,r,t,dist=register(s,selected_rotation,anchor,radius,gap,dx,tilt)
            target_b=evaluate_b(world,TARGETS)
            # Rigid-body covariance via inverse-mapping observers into source frame.
            source_b=evaluate_b(s,(TARGETS-t)@r)
            cov=np.linalg.norm(target_b-source_b@r.T)/np.linalg.norm(target_b)
            assert cov<1e-11
            shell_gap=None
            if shells is not None:
                c=(np.array([0,0,-radius])-t)@r
                shell_gap=closest_shell_distance(shells,c)-radius
                assert shell_gap>0, (key,case,shell_gap)
            pose_rows.append(dict(coil=key,case=case,sphere_radius_mm=radius,
                                  wire_gap_requested_mm=gap,lateral_shift_mm=dx,tilt_deg=tilt,
                                  **dist,shell_gap_mm=shell_gap,field_covariance_relative_L2=cov))
            placement[key]['poses'][case]=dict(source_to_world=transform_matrix(r,t).tolist(),
                                      sphere_centre_world_mm=[0,0,-radius],targets_world_mm=TARGETS.tolist())
            for p,v in zip(TARGETS,target_b):
                field_rows.append(dict(coil=key,case=case,x_mm=p[0],y_mm=p[1],z_mm=p[2],
                       Bx_mT=v[0]*1e3,By_mT=v[1]*1e3,Bz_mT=v[2]*1e3,Bnorm_mT=np.linalg.norm(v)*1e3))
            if case=='nominal':
                np.savetxt(HERE/'target_frame_segments'/f'{key}.csv',world,delimiter=',',header=HEADER,comments='')
    save_table('pose.csv',pose_rows);save_table('target_fields.csv',field_rows)
    orientation_cases=[]
    for case in ('initial_marked_normal','nominal'):
        pr=next(row for row in pose_rows if row['coil']=='h1' and row['case']==case)
        fr=[row for row in field_rows if row['coil']=='h1' and row['case']==case
            and row['x_mm']==0 and row['y_mm']==0]
        orientation_cases.append(dict(rule=case,
                tilt_from_initial_deg=(0. if case=='initial_marked_normal' else h[4]['tilt_from_initial_marked_fit_deg']),
                wire_gap_mm=pr['min_wire_gap_mm'],shell_gap_mm=pr['shell_gap_mm'],
                anchor_z_mm=pr['effective_anchor_z_mm'],closest_wire_to_apex_mm=pr['closest_wire_to_apex_mm'],
                B10_mT=fr[0]['Bnorm_mT'],B30_mT=fr[1]['Bnorm_mT'],B50_mT=fr[2]['Bnorm_mT']))
    save_table('h1_orientation_comparison.csv',orientation_cases)
    save_table('source_field_checks.csv',source_checks)
    (HERE/'pose.json').write_text(json.dumps({'units':'mm','matrix_convention':'world_xyz = R @ source_xyz + t for column vectors',
                                  'nominal_radius_mm':85,'nominal_wire_gap_mm':12,'current_A':CURRENT,'coils':placement},indent=2))
    for p in pose_rows:
        if p['case']=='nominal': print(p)


if __name__=='__main__':
    main()
