# 三类 TMS 线圈的摆位对照与非轴对称场验证

这个补充实验检验同一个问题：规定相同电流和靶点时，选取何种几何表面控制线圈间隙，会怎样影响三类绕组的磁场比较。计算对象是空气中的磁感应强度 B。球面定义距离和坐标，单位权重绕组电流为 1000 A，三种线圈分别保留自身的匝数和尺寸。相同三套绕组在个体 MRI 上的摆位及左丘脑电场见[同头模深部比较](../brain_field/README.md)。

## 同一线圈，两种摆位基准

以半径 85 mm 的球面模拟头皮轮廓，球顶设为原点。理想八字以双翼相切处定位，95° 双锥以对称中心和铰链定位，H1 用公开模型的工作面标记定位并固定形状。沿球顶法线平移，直至整条绕组中心线距球面最近处为 12 mm。H1 的完整外壳距球面最近处为 2.25 mm。各线圈球顶到最近绕组的直线距离分别为 H1 23.74、八字 12.00、双锥 37.34 mm。

平面对照保留各自球面方案的旋转、横向位置、电流、共同靶点。仅沿球顶法线平移，使整条绕组中心线到球顶切平面的最小净距也是 12 mm。这个特意设置的几何对照会将包覆式绕组整体抬高，可以显示把平面假设用于曲面绕组时产生的影响。旧的十四类实验还使用了其他摆位约定，数值请单独使用。

| 线圈 | 换用平面需抬高 mm | 球面下 10 mm 的 B 模 mT | 球面下 30 mm | 球面下 50 mm | 平面对照 30 mm | 30 mm 相对球面变化 |
|---|---:|---:|---:|---:|---:|---:|
| 固定形状 H1 | 111.41 | 55.89 | 42.78 | 32.31 | 8.34 | −80.5% |
| 理想八字 | 0.00 | 36.76 | 15.40 | 7.06 | 15.40 | 0.0% |
| 理想双锥 | 31.56 | 74.00 | 35.02 | 15.05 | 8.99 | −74.3% |

在这组名义球面摆位下，最强轴向靶点场随深度从 10 mm 的双锥转为 30、50 mm 的 H1。改用切平面后，30 mm 处三者的顺序翻转。H1 是固定零形变的公开几何，八字与双锥是理想绕组。这组结果适用于具体几何、激励和摆位；缺少相同能耗或同一仪器输出的约束。半径、间隙、横移与倾角的逐项扰动见[摆位说明](placement/README.md)和[全部目标值](placement/target_fields.csv)。[相同朝向的平面和球面对照](reference_surface_comparison.csv) · [摆位示意图](../figures/core/placement.png) · [原始 MATLAB 程序交叉核验](matlab_reference_check.csv)

同一个球面上也能看到排名对摆位的依赖。只用 H1 外壳标记拟合工作面法线时，H1 在 30 mm 的 B 模为 20.57 mT。以相同锚点为准，将 H1 朝向再旋转 25.54°、使工作面标记更接近球顶，并重新求 12 mm 净距后，B 模升至 42.78 mT。两种放置均未与球面相交；最近绕组到球顶的距离从 85.11 mm 缩至 23.74 mm。两种朝向都由几何确定，均未用靶场作为优化目标。双锥保持 35.02 mT，因此只改变 H1 工作面如何对准靶点，已足以交换二者在 30 mm 处的先后。[同一 H1 的朝向对照](placement/h1_orientation_comparison.csv)

## 对双锥场解的独立数值验证

对同一份 3840 段双锥绕组单独装配三维六面体有限元方程，解三个磁矢势分量并取旋度，用 9 组网格/外域条件检查误差。最终网格边长 1.5 mm、空气立方体半宽 480 mm。原坐标系 3571 个观测点的整体相对向量 L2 差为 0.401%，最大单点差为 2.479%。统一摆位下另核验 7 个靶点（含 4 个离轴点），最大相对向量差 0.535%。30 mm 轴上靶点的有限元与线元积分分别给出 34.952 和 35.024 mT。此检验覆盖同一理想绕组的场计算；商用双锥的实物绕组尺寸、外壳和电流需要另有测量支持。[方法与复现](fem3d/README.md) · [7 个靶点的逐点结果](fem3d/results/posed_double_cone.csv)

## 公开体模测量与既有研究

Petrov 等公开了商用八字线圈 MRI 体模的有符号 Bz 切片。读取发布的测量与模型数组，在同一掩膜的 10,622 个像素上用 MATLAB 重新评估，作者的单环、平面多匝、多层模型相对实测的 L2 差分别为 9.00%、7.36%、7.71%。这些数值衡量作者模型对其体模数据的符合程度。本项目的理想八字使用另一份几何。[本地复核说明与代码](public_benchmark/README.md) · [公开数据目录与适用范围](public_data.md) · 原始论文〔1〕

进一步按作者发布的两个版本重建参考件后，已查明九匝单层模型的 4.50 倍幅值差来自 v1.0 过时的神经导航配准矩阵。v1.1 更正使线圈原点沿扫描仪 z 移动 32.27 mm。保持等效电流 3.576 A 和 10,622 个像素不变，预测与实测 Bz 的均方根比由 4.5018 降至 1.0162，相对 L2 差由 352.56% 降至 7.63%。该实验检验公开八字参考件的这张切片；理想八字和双锥仍沿用各自的几何与验证条件。[图解、实验和复现说明](public_benchmark/matched_figure8/README.md) · 作者对旧摆位的说明〔2〕

已有研究系统比较了 50 种线圈的脑内电场深度和聚焦性，另有针对大型柔性线圈及双锥的真实头模摆位研究。本实验定量展示比较结论对摆位假设的依赖，并为非轴对称场给出可复算的三维误差。要论证脑内深部靶区刺激，应进一步研究有组织电导率的头模感应电场 E，明确波形或 dI/dt、靶区与非靶区指标，并核对物理模型与测量。单用本实验的 B 模得不出疗效或通用优胜线圈。50 种线圈比较〔3〕 · 柔性 H1 研究〔4〕 · 双锥摆位研究〔5〕

## 复现

源文件位于项目的 `coil_study` 目录。可直接使用现有 MATLAB 和 Python 环境。在新环境中，从 `coil_study` 运行 `py -3.13 -m venv .venv`，再运行 `& .\.venv\Scripts\python.exe -m pip install -r requirements.txt` 安装计算库。以下命令在 `coil_study` 下运行。验证使用 MATLAB R2025a，细网格有限元使用 CPU，整组有限元进程内存峰值约 5.91 GiB。

```powershell
& .\.venv\Scripts\python.exe .\placement_pilot\placement\run_placement.py
& .\.venv\Scripts\python.exe .\placement_pilot\compare_reference_surfaces.py
matlab -batch "addpath('placement_pilot'); verify_reference_surfaces"
& .\.venv\Scripts\python.exe .\placement_pilot\fem3d\run_double_cone.py
& .\.venv\Scripts\python.exe .\placement_pilot\fem3d\check_placed_targets.py
matlab -batch "addpath('placement_pilot/public_benchmark'); rescore_petrov"
matlab -batch "addpath('placement_pilot/public_benchmark/matched_figure8'); validate_petrov_figure8; plot_petrov_comparison"
```

## 论文与数据来源

1. [原始论文](https://doi.org/10.1371/journal.pone.0178952)
2. [作者对旧摆位的说明](https://github.com/pip010/plosone-data-pub01/releases/tag/v1.1)
3. [50 种线圈比较](https://pubmed.ncbi.nlm.nih.gov/22483681/)
4. [柔性 H1 研究](https://doi.org/10.1016/j.brs.2025.05.136)
5. [双锥摆位研究](https://pubmed.ncbi.nlm.nih.gov/40472866/)
