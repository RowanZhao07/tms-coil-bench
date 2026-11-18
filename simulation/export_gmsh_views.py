"""Export computed CSV magnetic fields into Gmsh's native post-processing format.

This script creates visualization files; it does not solve a finite-element PDE.
Prefer matching results/magpylib CSVs, falling back to results/fields (MATLAB).
Run using ../.venv/Scripts/python.exe export_gmsh_views.py [--validate]
"""
from __future__ import annotations
import argparse
import csv
import json
import pathlib
import shutil
import subprocess
import sys
import time
import numpy as np

ROOT = pathlib.Path(__file__).resolve().parents[1]
OUT = ROOT / 'results' / 'gmsh_views'

def locate_source(coil, plane):
    filename = f'{coil}_{plane}.csv'
    for folder, software in [('magpylib', 'Magpylib'), ('fields', 'MATLAB')]:
        path = ROOT / 'results' / folder / filename
        if path.exists():
            return path, software
    raise FileNotFoundError(filename)

def scalar_view(path, name, plane):
    table = np.genfromtxt(path, delimiter=',', names=True, encoding='utf-8-sig')
    axes = ('x_m', 'y_m') if plane == 'xy20' else ('x_m', 'z_m')
    a, b = np.unique(table[axes[0]]), np.unique(table[axes[1]])
    order = np.lexsort((table[axes[1]], table[axes[0]]))
    table = table[order]
    if len(table) != len(a)*len(b):
        raise ValueError(f'{path}: expected complete rectangular grid')
    xyz = np.column_stack([table[k] for k in ('x_m', 'y_m', 'z_m')]) * 1000
    mag = np.sqrt(sum(table[k]**2 for k in ('Bx_T', 'By_T', 'Bz_T'))) * 1000
    if not np.all(np.isfinite(xyz)) or not np.all(np.isfinite(mag)):
        raise ValueError(f'{path}: non-finite values cannot be visualized')
    xyz, mag = xyz.reshape(len(a),len(b),3), mag.reshape(len(a),len(b))
    lines = [f'View "{name} |B| [mT]; coordinates [mm]" {{']
    for i in range(len(a)-1):
        for j in range(len(b)-1):
            indices = ((i,j),(i+1,j),(i+1,j+1),(i,j+1))
            coords = ','.join(f'{v:.10g}' for ij in indices for v in xyz[ij])
            values = ','.join(f'{mag[ij]:.10g}' for ij in indices)
            lines.append(f'SQ({coords}){{{values}}};')
    lines.append('};')
    return '\n'.join(lines), dict(points=int(mag.size), quadrangles=(len(a)-1)*(len(b)-1),
                                  maximum_mT=float(mag.max()), minimum_mT=float(mag.min()))

def wire_view(coil):
    path = ROOT / 'results' / 'geometry' / f'{coil}.csv'
    t = np.genfromtxt(path, delimiter=',', names=True, encoding='utf-8-sig')
    lines = [f'View "{coil}: wire geometry [mm]" {{']
    for row in t:
        coords = ','.join(f'{1000*row[k]:.10g}' for k in ('ax_m','ay_m','az_m','bx_m','by_m','bz_m'))
        lines.append(f'SL({coords}){{1,1}};')
    lines.append('};')
    return '\n'.join(lines),len(t)

def viewer_geometry(coil):
    return f'''// Data viewer only. Imported B fields are not a Gmsh FEM solution.
Merge "{coil}.pos";
General.Axes = 3;
General.SmallAxes = 1;
General.Orthographic = 1;
General.RotationX = 55;
General.RotationY = 0;
General.RotationZ = 25;
View[0].ShowScale = 1;
View[1].ShowScale = 0;
View[2].ShowScale = 0;
View[2].LineWidth = 2;
View[2].ColorTable = {{ {{20,20,20}} }};
'''

def launcher(filename):
    return rf'''@echo off
set "GM=%~dp0..\..\tools\gmsh-4.15.2\gmsh-4.15.2-Windows64\gmsh.exe"
if not exist "%GM%" set "GM=%~dp0..\..\tools\gmsh-4.15.2-Windows64\gmsh.exe"
if exist "%GM%" (
  start "" "%GM%" "%~dp0{filename}"
) else (
  start "" gmsh "%~dp0{filename}"
)
'''

README = r'''# Gmsh 实验查看器

双击对应的 `open_线圈名.cmd`，或用 Gmsh 打开 `线圈名.geo`。同名 `.pos` 也可直接打开。三个独立图层分别是距线圈基准面20 mm的xy场、y=0的xz场、导线；在Gmsh左侧“Post-processing”中可分别显隐。

14类线圈的 **B由MATLAB或Magpylib计算，Gmsh在此仅负责导入与显示**。每个图层名称及 `manifest.json` 明确标出源软件、原始CSV路径和单位。导出脚本优先采用存在的 `results/magpylib/<coil>_<plane>.csv`，否则采用MATLAB的 `results/fields`；两个平面分别选源。未来新增Magpylib结果后重新运行脚本即可刷新。

14类查看器中的坐标单位为mm，色标为|B|的mT。绘图采用原始采样网格的四边形插值，细节分辨率受CSV采样间距限制。导线图层仅表示几何中心线，不表示导线内部有限截面场。单独打开`.pos`时，Gmsh默认也会给导线显示常数色标；打开配套`.geo`则隐藏此色标。

## 真正有限元结果

`open_fem_circle.cmd` 或 `fem_circle.geo` 打开 `fem/results/boundary_080` 中实际GetDP有限元解的B和A场。此查看器的坐标为SI米，B为T、A为T·m，和上述14类CSV查看器单位不同。圆环FEM的方程、收敛与运行日志位于 `fem/README.md`。只有该fem目录的结果来自本项目真实PDE求解，不能将14类场导入Gmsh的行为称作14类均已完成FEM。

`open_fem_main_circular.cmd` 或 `fem_main_circular.geo` 则打开主模型的18环圆形线圈实际FEM场，与MATLAB主结果的轴场对照保存在 `fem/results/multiturn_main_circular/main_model_comparison.csv`。

## 复现

在项目根目录运行：

```powershell
& .\.venv\Scripts\python.exe .\simulation\export_gmsh_views.py --validate
```

`--validate` 调用本机Gmsh解析圆形线圈查看器，再使用Gmsh Python API创建图形上下文并导出一张xy平面PNG；结果与日志写在本目录。软件二进制使用项目tools中的便携Gmsh，找不到时查询PATH；API优先使用tools/gmsh_python，也可通过pip安装gmsh==4.15.2。软件二进制不属于结果压缩包，普通查看不需要Python API。
'''

def main():
    ap=argparse.ArgumentParser(); ap.add_argument('--validate', action='store_true'); args=ap.parse_args()
    OUT.mkdir(parents=True,exist_ok=True)
    names=sorted(p.stem[:-5] for p in (ROOT/'results'/'fields').glob('*_xy20.csv'))
    manifest=[]
    for coil in names:
        parts=[]; record=dict(coil=coil,kind='imported_field_visualization_not_FEM',layers=[])
        for plane in ('xy20','xz'):
            source, software=locate_source(coil,plane)
            block, stats=scalar_view(source, f'{coil}: {plane}; source {software}', plane)
            parts.append(block); record['layers'].append(dict(plane=plane,software=software,
                source=str(source.relative_to(ROOT)).replace('\\','/'),**stats))
        block, n=wire_view(coil); parts.append(block); record['wire_segments']=n
        (OUT/f'{coil}.pos').write_text('\n'.join(parts)+'\n',encoding='utf-8')
        (OUT/f'{coil}.geo').write_text(viewer_geometry(coil),encoding='utf-8')
        (OUT/f'open_{coil}.cmd').write_text(launcher(coil+'.geo'),encoding='ascii')
        manifest.append(record)
    (OUT/'manifest.json').write_text(json.dumps(manifest,indent=2,ensure_ascii=False),encoding='utf-8')
    (OUT/'README.md').write_text(README,encoding='utf-8')
    (OUT/'fem_circle.geo').write_text('''// Actual GetDP axisymmetric FEM result. Coordinates m, B T, A T*m.
Merge "../../fem/results/boundary_080/field_b.pos";
Merge "../../fem/results/boundary_080/field_a.pos";
General.Axes = 3;
General.Orthographic = 1;
View[0].Name = "GetDP FEM: B [T], r-z coordinates [m]";
View[1].Name = "GetDP FEM: A_phi [T*m], r-z coordinates [m]";
View[1].Visible = 0;
''',encoding='utf-8')
    (OUT/'open_fem_circle.cmd').write_text(launcher('fem_circle.geo'),encoding='ascii')
    (OUT/'fem_main_circular.geo').write_text('''// Actual 18-turn main circular coil GetDP FEM. Coordinates m, B T, A T*m.
Merge "../../fem/results/multiturn_main_circular/field_b.pos";
Merge "../../fem/results/multiturn_main_circular/field_a.pos";
General.Axes = 3;
General.Orthographic = 1;
View[0].Name = "GetDP FEM: main 18-turn circular B [T]; coordinates [m]";
View[1].Name = "GetDP FEM: main 18-turn circular A_phi [T*m]";
View[1].Visible = 0;
''',encoding='utf-8')
    (OUT/'open_fem_main_circular.cmd').write_text(launcher('fem_main_circular.geo'),encoding='ascii')
    print(f'Exported {len(manifest)} coil viewers to {OUT}',flush=True)
    if args.validate:
        exes=list((ROOT/'tools').rglob('gmsh.exe'))
        exe=str(exes[0]) if exes else shutil.which('gmsh')
        if not exe: raise FileNotFoundError('Gmsh executable required for --validate')
        preview='''Include "circular.geo";
General.RotationX = 0; General.RotationY = 0; General.RotationZ = 0;
General.GraphicsWidth = 1200; General.GraphicsHeight = 900;
View[1].Visible=0; View[2].Visible=0;
'''
        (OUT/'preview_circular.geo').write_text(preview,encoding='utf-8')
        start=time.perf_counter()
        cmd=[exe,'preview_circular.geo','-0','-o','parsed_empty.msh']
        p=subprocess.run(cmd,cwd=OUT,capture_output=True,text=True,encoding='utf-8',errors='replace',timeout=120)
        log='COMMAND: '+subprocess.list2cmdline(cmd)+'\n'+p.stdout+p.stderr
        (OUT/'gmsh_validation.log').write_text(log,encoding='utf-8')
        if p.returncode or 'Error' in log: raise RuntimeError('Gmsh validation failed; inspect gmsh_validation.log')
        sys.path.insert(0,str(ROOT/'tools'/'gmsh_python'))
        import gmsh
        png=OUT/'circular_xy20_preview.png'
        gmsh.initialize()
        try:
            gmsh.logger.start()
            gmsh.option.setNumber('General.Terminal',1)
            gmsh.fltk.initialize()
            gmsh.open(str(OUT/'preview_circular.geo'))
            gmsh.fltk.update()
            gmsh.graphics.draw()
            gmsh.write(str(png))
            with (OUT/'gmsh_validation.log').open('a',encoding='utf-8') as f:
                f.write('\nGmsh Python API rendering:\n'+'\n'.join(gmsh.logger.get())+'\n')
            gmsh.logger.stop()
        finally:
            gmsh.finalize()
        if not png.is_file() or png.stat().st_size<1000: raise RuntimeError('Gmsh did not create a valid PNG')
        from PIL import Image
        pixels=np.asarray(Image.open(png).convert('RGB'))
        if float(pixels.std()) < 5:
            raise RuntimeError('Gmsh image is blank; graphics context did not render the field')
        (OUT/'validation.json').write_text(json.dumps(dict(success=True,seconds=time.perf_counter()-start,
                    preview_file=png.name,preview_bytes=png.stat().st_size,coil_viewers=len(manifest)),indent=2),encoding='utf-8')
        print(f'Gmsh parsed native views and exported {png}',flush=True)

if __name__=='__main__': main()
