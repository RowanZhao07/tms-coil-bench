# 线圈计算与仿真平台调研

核查日期：2026-11-18。只将直接打开的官方页面作为版本依据；检索摘要可能滞后。

## 建议选择

1. **主计算：MATLAB R2025a + 分段 Biot–Savart。** 针对空气芯闭合线圈，采用解析有限直线段积分、曲线几何离散和线性叠加，既能处理平面圈，也能处理折叠、H 型与阵列结构。物理模型应随材料/频率/组织导电性变化，不能仅因名字不同便随意更换基本定律。使用均匀空气磁导率时，无须训练深度学习模型。
2. **全体线圈的独立软件交叉检查：Magpylib 6.3.0。** 官方 PyPI 页面记录发布于 2026-11-07，支持 Python 3.11+、NumPy/SciPy，基于向量化解析磁静态表达式。Circle 适合圆环独立解析校验，Polyline 可表达折线电流；输出 B 为 T，几何输入为 m，电流为 A。必须固定版本和单位。此比较是独立实现验证，并非独立 FEM 验证。
3. **有限元层：Gmsh + GetDP。** Gmsh 为网格与后处理软件，GetDP 才是 PDE/FEM 求解器。主页核实 Gmsh 稳定版 4.15.2（2026-03-24）；GetDP 主页仍将 3.5.0（2022-05-13）列为稳定版，而当前在线手册标题已出现 4.0.0，故不能把手册标题当成已安装/已正式发布的版本。正式实验应记录实际二进制 `--version`。GetDP 支持 1D/2D/轴对称/3D，以及静态、时域、频域问题，Windows 二进制包含 PETSc/MUMPS。对单机，先做圆形线圈轴对称 FEM 基准，再增加代表性 3D 模型，可降低内存和网格调试成本。
4. **需要人体感应电场时：SimNIBS 4.6.0。** 官方在线文档为 4.6.0，提供 TMS coil + 头部网格 + 电导率 + dI/dt 的有限元流程。安装约占 3 GB 磁盘。SimNIBS 的 H1/H4/H7 柔性线圈基于 X-ray 重建，可作为 H1 结构的可靠资料来源，但不能将未核实尺寸的自建近似线圈冒称为真实商业 H1。

## 单机 RTX 4060 的适配边界

Magpylib 的普通接口基于 NumPy/SciPy，选用它不要求 CUDA，也不会自动利用 RTX 4060。GetDP/SimNIBS 的上述官方基本流程以 CPU 和主内存为主要资源。能够运行的网格规模取决于 RAM、离散阶次、求解器和域大小，不能只凭显卡型号保证任意 3D 模型可运行。建议分块计算观测点、顺序处理线圈、限制并发 FEM 作业，并以实际运行时间与峰值内存作为本机证据。若 MATLAB 已装 Parallel Computing Toolbox，可另行评估 gpuArray；无需为了“最新模型”引入未经验证的神经网络近似。

## B 与 E 的区别

磁静态计算给出某个电流 I 下的 B(I)，在线性空气模型中 B 与 I 成正比。TMS 的组织电场包含由线圈矢势时间变化产生的初级场与导电组织边界产生的电势项，即 E = -∂A/∂t - ∇φ；它与 dI/dt 线性相关。没有脉冲波形/dI/dt、导电率和头模时，不能从 B 峰值或 B 半高区直接宣称刺激深度、临床疗效或儿童适用性。“儿童线圈”仅应作为几何类别名称。

## 主要官方来源

- Magpylib PyPI（6.3.0、发布日期、Python 要求）：https://pypi.org/project/magpylib/
- Magpylib current API（Circle、Polyline、SI 单位与 B 输出）：https://magpylib.readthedocs.io/en/stable/_autogen/magpylib.current.html
- Magpylib 官方仓库（计算模型）：https://github.com/magpylib/magpylib
- Gmsh 下载和版本主页：https://gmsh.info/?lang=en
- GetDP 能力、稳定下载与 Windows 求解库：https://getdp.info/
- GetDP 在线手册与磁静态教程：https://www.getdp.info/doc/texinfo/getdp.html
- SimNIBS 安装（约 3 GB 磁盘）：https://simnibs.github.io/simnibs/build/html/installation/simnibs_installer.html
- SimNIBS TMS FEM API：https://simnibs.github.io/simnibs/build/html/documentation/python_lib/auto/simnibs.simulation.fem.tms_coil.html
- SimNIBS FAQ（dI/dt、单位、E 缩放）：https://simnibs.github.io/simnibs/build/html/faq.html
- SimNIBS 官方 coil 文档源（H1/H4/H7 与测量重建模型）：https://github.com/simnibs/simnibs/blob/master/docs/documentation/coils.rst

## 尚未完成、不能声称的项目

本调研未安装任何平台、未执行 FEM、未测运行时间或峰值内存；以上仅为平台选型依据。报告中的实验结果必须来自后续实际日志。没有取得 Elmer 当前版本和硬件要求的充分官方证据，故不将其列为本次优选执行平台。
