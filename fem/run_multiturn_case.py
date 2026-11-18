"""Actual 18-turn main circular model FEM, separate from the single-turn benchmark.
Each circular conductor cross section has its own physical group and exactly 1000 A.
"""
from __future__ import annotations
import json, math, pathlib, re
import numpy as np
from numpy.polynomial.legendre import leggauss
from run_fem import executable, run, write_csv

BASE=pathlib.Path(__file__).resolve().parent
ROOT=BASE.parent
OUT=BASE/'results'/'multiturn_main_circular'
MU0=4*math.pi*1e-7
LOOPS=[(float(r),float(z)) for z in [0.,-.0045] for r in np.linspace(.017,.055,9)]
A=.002; CURRENT=1000.; L=.8

def geometry():
    lines=[f'''L={L}; la=0.001; lo=0.025; lc=0.0005;
Point(1)={{0,0,0,la}};
Point(2)={{0,-L,0,lo}}; Point(3)={{L,0,0,lo}}; Point(4)={{0,L,0,lo}};
Circle(1)={{2,1,3}}; Circle(2)={{3,1,4}};
Point(5)={{0,-0.08,0,la}}; Point(6)={{0,0.08,0,la}};
Line(3)={{4,6}}; Line(4)={{6,1}}; Line(5)={{1,5}}; Line(6)={{5,2}};
Curve Loop(1)={{1,2,3,4,5,6}};''']
    for i,(r,z) in enumerate(LOOPS):
        p=100+i*5; c=100+i*4; loop=10+i
        lines.append(f'''Point({p})={{{r},{z},0,lc}};
Point({p+1})={{{r+A},{z},0,lc}}; Point({p+2})={{{r},{z+A},0,lc}};
Point({p+3})={{{r-A},{z},0,lc}}; Point({p+4})={{{r},{z-A},0,lc}};
Circle({c})={{{p+1},{p},{p+2}}}; Circle({c+1})={{{p+2},{p},{p+3}}};
Circle({c+2})={{{p+3},{p},{p+4}}}; Circle({c+3})={{{p+4},{p},{p+1}}};
Curve Loop({loop})={{{c},{c+1},{c+2},{c+3}}};
Plane Surface({loop})={{{loop}}};
Physical Surface("turn_{i+1}",{100+i})={{{loop}}};
Transfinite Curve {{{c},{c+1},{c+2},{c+3}}}=17;''')
    holes=','.join(str(10+i) for i in range(18))
    lines.append(f'''Plane Surface(1)={{1,{holes}}}; Physical Surface("air",1)={{1}};
Physical Curve("outer_u0",3)={{1,2}}; Physical Curve("axis_natural",4)={{3,4,5,6}};
Mesh.Algorithm=6; Mesh.MshFileVersion=2.2;''')
    return '\n'.join(lines)

def mesh_areas(path):
    lines=path.read_text().splitlines(); ni=lines.index('$Nodes'); n=int(lines[ni+1]); nodes={}
    for line in lines[ni+2:ni+2+n]:
        t=line.split(); nodes[int(t[0])]=np.array([float(t[1]),float(t[2])])
    ei=lines.index('$Elements'); ne=int(lines[ei+1]); areas={i:0. for i in range(100,118)}; nt=0
    for line in lines[ei+2:ei+2+ne]:
        t=[int(x) for x in line.split()]
        if t[1]==2:
            nt+=1
            if t[3] in areas:
                a,b,c=(nodes[k] for k in t[3+t[2]:])
                areas[t[3]]+=abs((b[0]-a[0])*(c[1]-a[1])-(b[1]-a[1])*(c[0]-a[0]))/2
    assert all(v>0 for v in areas.values())
    return n,nt,areas

def references(z,nradial=64,ntheta=256):
    q,w=leggauss(nradial); rho=(q+1)*A/2; wr=w*A/2
    theta=2*np.pi*np.arange(ntheta)/ntheta
    finite=0.; filament=0.
    for radius,center in LOOPS:
        rr=radius+rho[:,None]*np.cos(theta)[None,:]
        zz=center+rho[:,None]*np.sin(theta)[None,:]
        kernel=MU0*rr**2/(2*(rr**2+(z-zz)**2)**1.5)
        finite+=CURRENT/(np.pi*A*A)*np.sum(kernel*rho[:,None]*wr[:,None])*(2*np.pi/ntheta)
        filament+=MU0*CURRENT*radius**2/(2*(radius**2+(z-center)**2)**1.5)
    return finite,filament

def main():
    OUT.mkdir(parents=True,exist_ok=True)
    (OUT/'circle.geo').write_text(geometry(),encoding='utf-8')
    mt,_=run([executable('gmsh'),'circle.geo','-2','-format','msh2','-o','circle.msh','-nt','1'],OUT,OUT/'gmsh.log')
    n,nt,areas=mesh_areas(OUT/'circle.msh')
    turn_rows=[dict(turn=i+1,physical_group=100+i,radius_m=r,z_center_m=z,wire_radius_m=A,
                   mesh_area_m2=areas[100+i],current_density_A_m2=CURRENT/areas[100+i],current_A=CURRENT)
               for i,(r,z) in enumerate(LOOPS)]
    write_csv(OUT/'turn_currents.csv',turn_rows)
    tags=','.join(str(i) for i in range(100,118))
    pro=(BASE/'circle_axisymmetric.pro').read_text()
    pro=pro.replace('Coil = Region[2];',f'Coil = Region[{{{tags}}}];')
    pro=pro.replace('Domain = Region[{1,2}]',f'Domain = Region[{{1,{tags}}}]')
    pro=pro.replace('source[Coil] = (4.e-7*Pi)*X[]^2*CurrentDensity;',
        '\n'.join(f'source[Region[{i}]] = (4.e-7*Pi)*X[]^2*{CURRENT/areas[i]:.17g};' for i in range(100,118)))
    (OUT/'case_parameters.pro').write_text('// Individual turn current densities are in circle.pro.\n')
    (OUT/'circle.pro').write_text(pro,encoding='utf-8')
    st,log=run([executable('getdp'),'circle.pro','-msh','circle.msh','-solve','SolveA','-pos','Export','-v','4'],OUT,OUT/'getdp.log')
    raw=np.loadtxt(OUT/'axis_b.txt'); zz=np.linspace(.01,.05,5)
    assert np.allclose(raw[:,3],zz) and np.all(np.isfinite(raw)) and np.max(abs(raw[:,-3]))<1e-10
    matlab=np.genfromtxt(ROOT/'results'/'fields'/'circular_axis.csv',delimiter=',',names=True,encoding='utf-8-sig')
    rows=[]
    for i,z in enumerate(zz):
        finite,filament=references(z); high,_=references(z,96,384)
        assert abs(finite-high)/high<1e-9, 'Finite cross-section quadrature not converged'
        match=matlab[np.isclose(matlab['z_m'],z,rtol=0,atol=1e-10)]
        assert len(match)==1
        m=float(match['Bz_T'][0]); f=float(raw[i,-2])
        rows.append(dict(z_m=z,Bz_FEM_T=f,Bz_MATLAB_T=m,Bz_filament_exact_T=filament,
                         Bz_finite_cross_T=finite,FEM_vs_MATLAB_pct=100*(f-m)/m,
                         FEM_vs_finite_cross_pct=100*(f-finite)/finite,
                         finite_cross_vs_filament_pct=100*(finite-filament)/filament,
                         MATLAB_vs_filament_pct=100*(m-filament)/filament))
    write_csv(OUT/'main_model_comparison.csv',rows)
    mem=re.findall(r'Mem = ([\d.]+)Mb',log)
    summary=dict(coil='circular',turns=18,current_per_turn_A=CURRENT,wire_radius_m=A,domain_radius_m=L,
                 nodes=n,triangles=nt,gmsh_seconds=mt,getdp_seconds=st,getdp_reported_memory_MB=float(mem[-1]),
                 max_abs_FEM_vs_MATLAB_pct=max(abs(r['FEM_vs_MATLAB_pct']) for r in rows),
                 max_abs_FEM_vs_finite_cross_pct=max(abs(r['FEM_vs_finite_cross_pct']) for r in rows),
                 max_abs_MATLAB_vs_filament_pct=max(abs(r['MATLAB_vs_filament_pct']) for r in rows),
                 note='Independent FEM of 18 finite circular wire cross-sections; MATLAB model is 192-segment polygonal filaments. Each turn separately normalized to 1000 A. This additional case has not undergone its own mesh convergence sweep.')
    (OUT/'summary.json').write_text(json.dumps(summary,indent=2),encoding='utf-8')
    print(json.dumps(summary,indent=2),flush=True)

if __name__=='__main__': main()
