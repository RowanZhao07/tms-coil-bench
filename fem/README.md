# 圆环线圈独立有限元验证（已实际执行）

## 模型与软件

- Gmsh **4.15.2** 生成三角网格；GetDP **3.5.0** 执行磁静态有限元。版本已由本机二进制 `-version` 核实，日志位于 `results/`。
- 单匝空气芯圆环：主半径 R=50 mm、导线圆截面半径 a=1 mm、电流 I=1000 A，所有区域 μ=μ0。
- 计算域是二维 r-z 半圆盘，旋转后为半径 L 的三维球域。轴对称是该圆环的精确降维，不能用于直接验证非轴对称 H1 或折叠线圈。
- 导线内电流密度均匀。为消除线性边界多边形面积导致的总电流误差，使用实际源区域三角形面积 S_mesh，令 J=I/S_mesh。导线形状仍随网格加密逼近圆截面。
- 外边界磁矢势设零；轴线上采用解的正则性。没有铁磁体、人体组织、涡流、驱动电路或热效应。

## 方程

采用轴对称磁矢势 Aφ，并引入 u=Aφ/r。由磁静态 Maxwell 方程可得：

```text
-div(r^3 grad u) = μ0 r^2 Jφ
Br = -r ∂u/∂z
Bz = 2u + r ∂u/∂r
```

使用一阶三角形节点基函数，七点三角形 Gauss 积分，GetDP 的稀疏求解器。该变量变换避免在轴上直接除以 r；轴上精确地有 Br=0、Bz=2u。Gmsh/GetDP 坐标 X 对应物理 r，Y 对应物理 z，输出表中的 By 因此是物理 Bz。

这是实际求解 PDE 的独立有限元实验；不是将 MATLAB 或 Magpylib 的 B 数据导入网格后显示。

## 实测结果与收敛

在轴上 z=10、20、30、40、50 mm 采样。有限截面参考值通过对圆截面内圆环解析轴场作独立高阶积分获得；另给出半径50 mm、1000 A细线圆环闭式解析值。

| 工况 | 球域半径/m | 节点 | 三角形 | 对有限截面参考最大误差 |
|---|---:|---:|---:|---:|
| 粗网格 | 0.5 | 3,593 | 7,014 | 0.6620% |
| 中网格 | 0.5 | 13,863 | 27,386 | 0.2968% |
| 细网格 | 0.5 | 54,527 | 108,380 | 0.2857% |
| 细网格、小外域 | 0.3 | 36,045 | 71,572 | 1.3389% |
| 细网格、大外域 | 0.8 | 85,846 | 170,780 | **0.0649%** |

0.5 m 外域下，加密到细网格后误差改善变小，说明外边界截断已成为主要误差来源。保持局部网格尺度、将外域从0.3扩大到0.5、0.8 m时，误差持续减小。最终工况对细线解析的最大误差为0.0624%；有限截面与细线不是严格同一物理模型，不能把两者差值全归为求解器误差。

最终0.8 m工况，本次实测 Gmsh约5.1 s、GetDP约9.6 s，GetDP日志记录内存约357 MB。数值会随机器负载变化；这是CPU计算，未使用RTX4060。上述内存为GetDP自身报告值，不是整机或Python/Gmsh合计峰值。

## 文件与复现

- `run_fem.py`：生成全部五组几何、网格、求解输入，实际启动软件，解析原始表格，导出比较CSV；需要Python和NumPy。
- `circle_axisymmetric.pro`：可检查的独立 PDE 弱式、边界与后处理定义。
- `results/axis_comparison.csv`：25个观测值、两类参考值及有符号相对误差。
- `results/convergence_summary.csv`：网格、域、运行时间、最大误差。
- `results/<case>/gmsh.log` 与 `getdp.log`：实际求解日志。
- `results/<case>/circle.geo/.msh/.pro/.pre/.res`：几何、网格、求解器输入与解。
- `results/<case>/field_b.pos`、`field_a.pos`：可由Gmsh打开的场图。坐标和单位为SI；完整场图包括线圈附近，峰值不可与距线圈1 cm的轴场混用。

在 PowerShell 运行：

```powershell
& .\run_fem.ps1
# 快速单工况检查：
& .\run_fem.ps1 -Quick
```

脚本从 `../tools/` 递归寻找 `gmsh.exe`、`getdp.exe`，找不到时使用PATH。本机便携软件已置于tools，最终分发压缩包应排除软件二进制；换机时从官方获取并放入tools即可。


本实验只为代表性圆环提供独立FEM证据。其他13类线圈的磁场计算需各自验证，不因本圆环测试通过而自动成为已完成FEM验证；此结果也不是物理测量或人体安全验证。

## 追加：与主MATLAB模型对应的18环圆形线圈

`run_multiturn_case.py` 已实际求解主圆线圈：两层z=0/-4.5 mm，每层9环，半径17～55 mm、间隔4.75 mm；每环1000 A；导线圆截面半径2 mm。18个截面各有独立物理组，分别按实际网格面积设置J_i=1000/S_i，详见 `results/multiturn_main_circular/turn_currents.csv`。

外域半径0.8 m，共41,493节点和82,528三角形；GetDP约4.71 s，自报内存178.766 MB。z=10～50 mm五个轴点上，FEM与主MATLAB结果最大差为 **0.05720%**，与18个有限截面圆环的独立高阶积分最大差为 **0.02961%**。高阶积分以64×256和96×384两组节点检查，二者相对差小于1e-9。

复算生成的结果、模型与求解日志位于 `results/multiturn_main_circular/`，直接比较表为 `main_model_comparison.csv`。本新增工况没有独立进行网格/边界扫参，不能将此前单匝模型的收敛结果等同为18环工况已完成同样扫参；其有限圆截面也与MATLAB的192段多边形细线模型有所区别。

复现命令：`..\.venv\Scripts\python.exe run_multiturn_case.py`。运行完整批次生成 Gmsh 视图后，可从仓库根目录打开 `results/gmsh_views/open_fem_main_circular.cmd` 查看有限元场。

## 软件来源

- Gmsh 4.15.2 Windows包：https://gmsh.info/bin/Windows/gmsh-4.15.2-Windows64.zip
- GetDP 3.5.0 Windows包：https://getdp.info/bin/Windows/getdp-3.5.0-Windows64c.zip
- Gmsh主页：https://gmsh.info/
- GetDP主页与手册：https://getdp.info/
