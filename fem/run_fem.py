"""Run genuine Gmsh/GetDP axisymmetric FEM and independent analytic benchmarks.

Run: ..\\.venv\\Scripts\\python.exe run_fem.py
Portable solver executables live outside deliverable source in ../tools/.
"""
from __future__ import annotations
import argparse, csv, json, math, pathlib, re, shutil, subprocess, time
import numpy as np
from numpy.polynomial.legendre import leggauss

ROOT = pathlib.Path(__file__).resolve().parent
TOOLS = ROOT.parent / 'tools'
MU0 = 4 * math.pi * 1e-7
R, WIRE_A, CURRENT = .05, .001, 1000.

def executable(name):
    found = list(TOOLS.rglob(name + '.exe'))
    if found: return found[0]
    found = shutil.which(name)
    if found: return pathlib.Path(found)
    raise FileNotFoundError(f'{name} not found under {TOOLS} or PATH')

def run(cmd, cwd, logfile):
    start=time.perf_counter()
    p=subprocess.run([str(x) for x in cmd], cwd=cwd, capture_output=True,
                     text=True, encoding='utf-8', errors='replace')
    logfile.write_text('COMMAND: '+subprocess.list2cmdline([str(x) for x in cmd])+'\n'+p.stdout+p.stderr,
                       encoding='utf-8')
    if p.returncode: raise RuntimeError(f'Failure {p.returncode}; inspect {logfile}')
    return time.perf_counter()-start, p.stdout+p.stderr

def geometry(L, scale):
    # Positive-radius half-disk; radius L = 0.3 / 0.5 / 0.8 m.
    # Local source and axis refinements independent of outer radius.
    return f'''R={R}; a={WIRE_A}; L={L}; s={scale};
lc=a/4*s; la=0.001*s; lo=0.025*s;
Point(1)={{0,0,0,la}};
Point(2)={{0,-L,0,lo}}; Point(3)={{L,0,0,lo}}; Point(4)={{0,L,0,lo}};
Circle(1)={{2,1,3}}; Circle(2)={{3,1,4}};
Point(5)={{0,-0.08,0,la}}; Point(6)={{0,0.08,0,la}};
Line(3)={{4,6}}; Line(4)={{6,1}}; Line(5)={{1,5}}; Line(6)={{5,2}};
Point(10)={{R,0,0,lc}};
Point(11)={{R+a,0,0,lc}}; Point(12)={{R,a,0,lc}};
Point(13)={{R-a,0,0,lc}}; Point(14)={{R,-a,0,lc}};
Circle(11)={{11,10,12}}; Circle(12)={{12,10,13}};
Circle(13)={{13,10,14}}; Circle(14)={{14,10,11}};
Curve Loop(1)={{1,2,3,4,5,6}}; Curve Loop(2)={{11,12,13,14}};
Plane Surface(1)={{1,2}}; Plane Surface(2)={{2}};
Physical Surface("air",1)={{1}}; Physical Surface("coil",2)={{2}};
Physical Curve("A0_outer",3)={{1,2}};
Physical Curve("axis_natural",4)={{3,4,5,6}};
Mesh.Algorithm=6; Mesh.MshFileVersion=2.2;
Transfinite Curve {{11,12,13,14}} = {max(9,round(17/scale))};
'''

def mesh_info(path):
    lines=path.read_text().splitlines()
    ni=lines.index('$Nodes'); n=int(lines[ni+1]); nodes={}
    for l in lines[ni+2:ni+2+n]:
        t=l.split(); nodes[int(t[0])]=np.array([float(t[1]),float(t[2])])
    ei=lines.index('$Elements'); ne=int(lines[ei+1]); area=0.; nt=0
    for l in lines[ei+2:ei+2+ne]:
        t=[int(x) for x in l.split()]
        if t[1]==2:
            nt+=1
            if t[3]==2:
                a,b,c=(nodes[k] for k in t[3+t[2]:])
                area+=abs((b[0]-a[0])*(c[1]-a[1])-(b[1]-a[1])*(c[0]-a[0]))/2
    return n,nt,area

def finite_cross_section(z):
    # Independent quadrature over a circular conductor cross section.
    q,w=leggauss(64); rho=(q+1)*WIRE_A/2; wr=w*WIRE_A/2
    theta=2*np.pi*np.arange(256)/256
    rr=R+rho[:,None]*np.cos(theta)[None,:]
    zz=rho[:,None]*np.sin(theta)[None,:]
    kern=MU0*rr**2/(2*(rr**2+(z-zz)**2)**1.5)
    return CURRENT/(np.pi*WIRE_A**2)*np.sum(kern*rho[:,None]*wr[:,None])*(2*np.pi/256)

def write_csv(path, rows):
    with path.open('w',newline='',encoding='utf-8-sig') as f:
        writer=csv.DictWriter(f,fieldnames=list(rows[0])); writer.writeheader(); writer.writerows(rows)

def main():
    ap=argparse.ArgumentParser(); ap.add_argument('--quick',action='store_true')
    args=ap.parse_args()
    gmsh,getdp=executable('gmsh'),executable('getdp')
    out=ROOT/'results'; out.mkdir(exist_ok=True)
    versions={}
    for n,e in [('gmsh',gmsh),('getdp',getdp)]:
        _,txt=run([e,'-version'],ROOT,out/(n+'_version.log')); versions[n]=txt.strip()
    cases=[(.5,2.,'mesh_coarse'),(.5,1.,'mesh_medium'),(.5,.5,'mesh_fine'),
           (.3,.5,'boundary_030'),(.8,.5,'boundary_080')]
    if args.quick: cases=cases[:1]
    records=[]; summary=[]
    for L,scale,label in cases:
        case=out/label; case.mkdir(exist_ok=True)
        (case/'circle.geo').write_text(geometry(L,scale),encoding='utf-8')
        mt,_=run([gmsh,'circle.geo','-2','-format','msh2','-o','circle.msh','-nt','1'],case,case/'gmsh.log')
        nodes,tris,area=mesh_info(case/'circle.msh')
        (case/'case_parameters.pro').write_text(f'CurrentDensity = {CURRENT/area:.17g};\n',encoding='ascii')
        shutil.copy2(ROOT/'circle_axisymmetric.pro',case/'circle.pro')
        st,txt=run([getdp,'circle.pro','-msh','circle.msh','-solve','SolveA','-pos','Export','-v','4'],case,case/'getdp.log')
        # GetDP Table columns: entity type, element id, xyz, local xyz, Bxyz.
        data=np.loadtxt(case/'axis_b.txt',ndmin=2)
        np.savetxt(case/'axis_raw_columns.csv',data,delimiter=',')
        if data.shape[0]!=5: raise RuntimeError(f'Expected five samples; got {data.shape}')
        print(label,'raw shape',data.shape,'first',data[0],flush=True)
        bz=data[:,-2]
        z=np.linspace(.01,.05,5)
        assert data.shape[1] == 11 and np.allclose(data[:,2],0) and np.allclose(data[:,3],z)
        assert np.all(np.isfinite(data)) and np.all(bz>0)
        assert np.max(abs(data[:,-3])) < 1e-10, 'Axis radial B must vanish by regularity'
        analytic=MU0*CURRENT*R**2/(2*(R**2+z*z)**1.5)
        finite=np.array([finite_cross_section(v) for v in z])
        for i in range(5):
            records.append(dict(case=label,domain_radius_m=L,mesh_scale=scale,nodes=nodes,triangles=tris,
                z_m=z[i],Bz_FEM_T=bz[i],Bz_filament_T=analytic[i],Bz_finite_cross_T=finite[i],
                relative_error_filament_pct=100*(bz[i]-analytic[i])/analytic[i],
                relative_error_finite_pct=100*(bz[i]-finite[i])/finite[i]))
        mem=re.findall(r'Mem = ([\d.]+)Mb',txt)
        s=dict(case=label,domain_radius_m=L,mesh_scale=scale,nodes=nodes,triangles=tris,
          mesh_source_area_m2=area,source_current_A=CURRENT,source_current_density_A_m2=CURRENT/area,
          gmsh_seconds=mt,getdp_seconds=st,getdp_reported_memory_MB=float(mem[-1]) if mem else None,
          max_abs_finite_error_pct=float(np.max(abs((bz-finite)/finite))*100),
          max_abs_filament_error_pct=float(np.max(abs((bz-analytic)/analytic))*100))
        summary.append(s); print(json.dumps(s),flush=True)
    write_csv(out/'axis_comparison.csv',records); write_csv(out/'convergence_summary.csv',summary)
    (out/'metadata.json').write_text(json.dumps(dict(versions=versions,R_m=R,wire_radius_m=WIRE_A,
          current_A=CURRENT,axis_evaluation_radius_m=0,formulation='axisymmetric transformed A_phi/r FEM; first-order triangles; seven-point quadrature',
          boundary='A_phi/r=0 on finite half-disk outer boundary; natural regularity on symmetry axis',
          notes='Coil J normalized using actual triangular cross-sectional mesh area; finite circular cross-section differs from filament'),indent=2),encoding='utf-8')

if __name__=='__main__': main()
