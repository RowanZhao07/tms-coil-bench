# 从头皮上方的绕组到左丘脑的感应电场

这里比较三个固定几何模型。H1 沿用公开重建绕组的零形变状态，八字沿用台架的理想八字，双锥沿用折起两翼的参数模型。三者与球面台架实验使用同一套线段坐标、方向与安匝权重。[绕组与台架摆位](../placement_pilot/placement/README.md) · [脑内结果图册](../figures/brain/README.md)

## 先确定比较位置

头模取自 FieldTrip 的 Subject01。作者将一名受试者的 T1 加权 MRI 转为右、前、上（RAS）坐标，用 CHARM 分割组织并生成四面体。我们直接使用作者发布的同一张网格及其原生解剖分区。解剖标签 10 对应左丘脑，3 和 42 对应左右大脑皮层。原始 T1 体素约 1.25 mm，各组织电导率按 SimNIBS 的各向同性设置。

1. 求整块左丘脑标签的体素中心均值，选距离该均值最近的左丘脑体素中心作为靶心。所得个体 RAS 坐标为 $(-7.40,-20.08,15.24)$ mm。
2. 在以靶心为球心、半径 5 mm 的范围内，保留重心既属于左丘脑标签又属于网格灰质的四面体。共 62 个单元，总体积 508.75 mm³。半径 3 mm 和 7.5 mm 用作敏感性检验。各线圈使用完全相同的单元编号和体积。
3. 从靶心沿 MRI 上方投射射线至头皮。锚点比靶心高 80.00 mm。以周围 12 mm 邻域的面积加权皮面法线和解剖右方向建立切平面坐标系，继承球面实验的绕组朝向，再沿法线平移。三种绕组的中心线到真实皮面的最短距离均为 12.00 mm；H1 的公开外壳离皮面最近 2.47 mm。几何定位在读取电场前完成。[靶区与统计口径](results/target.json) · [三组刚体变换与实测几何间距](results/placement.json)

## 脑组织怎样改变电场

以相同的电流变化率 $\dot I=10^8\ \mathrm{A/s}$（100 A/µs）驱动三种线圈。在磁准静态近似下，绕组提供初级场 $-\partial\mathbf A/\partial t$，导电组织界面和头部边界还会产生由电势梯度描述的场。SimNIBS 求解

$$
\mathbf E=-\frac{\partial\mathbf A}{\partial t}-\nabla\phi,\qquad
\nabla\!\cdot(\sigma\mathbf E)=0,\qquad
\mathbf n\!\cdot\sigma\mathbf E=0\quad\text{（外边界）}.
$$

$\sigma$ 表示各组织电导率，$\mathbf A$ 是矢势，$\phi$ 是电势。有限元在每个四面体给出总场的三个分量。先求每个单元的模长再按体积平均，能防止朝向不同的电场相互抵消。另取整个双侧大脑皮层灰质按体积加权的第 99 百分位作皮层强场参照：

$$
\bar E_R=\frac{\sum_{i\in R}V_i\|\mathbf E_i\|}{\sum_{i\in R}V_i},
\qquad
\eta=100\,\frac{\bar E_R}{Q^{(V)}_{0.99}(\|\mathbf E\|_{\rm cortex})}\%.
$$

百分位按累计体积达到 99% 时的场值确定。$\bar E_R$ 越大，同一驱动下靶区平均场越强；$\eta$ 越高，则相对于该线圈自身的皮层强场，靶区保留的场越多。两者衡量不同问题。空气对照统一用 1,000 A，在相同靶区的单元重心计算 $|\mathbf B|$。磁场与电场分别排序，单位分别是 mT 和 V/m。

## 结果与数值检验

按 H1、双锥、八字排列，靶区平均初级场为 **50.18、30.86、10.20 V/m**，加入组织电势后，总场分别为 **10.72、2.23、0.65 V/m**。各单元先对矢量求和，再取模长，才能计算总场均值。三种线圈的靶区均值分别相当于各自双侧大脑皮层 P99 的 **8.69%、3.36%、2.17%**。把靶区半径改为 3 mm 或 7.5 mm，三种线圈的排序相同。[主表](results/summary.csv) · [半径与纳入体积](results/roi_radius.csv)

- **绕组是否按原样放进求解器。** 直接读取导出的三只线圈文件，与原绕组的有限直线段积分对照。空气磁场的三分量相对 L2 差最高为 0.0011%。另在靶区单元涉及的 38 个节点直接积分磁矢势，再按四个顶点平均；与有限元输入的初级场向量相对 L2 差最高为 0.00174%。[几何与磁场核验](results/air_field.csv) · [矢势核验](results/source_check.csv)
- **线性系统是否收敛。** 4,734,021 个四面体、839,118 个节点共用同一刚度矩阵。PETSc 共轭梯度与 HYPRE 多重网格对每只线圈都迭代 37 次，自由节点相对残差为 $2.20\times10^{-10}$ 至 $3.61\times10^{-10}$。这个残差度量方程求解精度；该头模的分割误差、组织电导率与绕组形状仍会影响预测。[运行记录与各组织电导率](results/solver_run.json)
- **单机开销。** 在 i7-14700HX、32 GB 内存、RTX 4060 Laptop 的 Windows 笔记本上，三组 CPU 求解从读网格至写完场文件共 118.45 s，峰值工作集为 3.39 GiB。完整场文件每只约 561 MiB。首次下载、摆位、统计与绘图耗时另计。

## 复跑与核查

从仓库根目录执行。下列命令使用 SimNIBS 4.6.0 安装附带的 Python 3.11 环境；若已激活相同环境，可将 `simnibs_python` 替换为该环境的 `python`。磁场对照、有限元和统计均可在 CPU 运行。头模、三份完整三维场和中间缓存写到 `tmp/brain_field/`，公开输入按文件大小和 SHA-256 核验。仓库保存逐单元靶区结果、运行记录及可独立重绘的截面。

~~~powershell
simnibs_python brain_field/prepare_head.py
simnibs_python brain_field/define_target.py
simnibs_python brain_field/geometry.py
simnibs_python brain_field/solve.py
simnibs_python brain_field/analyze.py
simnibs_python brain_field/check_primary.py
simnibs_python brain_field/extract_section.py
~~~

绘图脚本和排版文件见[图册与重绘命令](../figures/brain/README.md)。[同条件主结果](results/summary.csv) · [3、5、7.5 mm 半径对照](results/roi_radius.csv) · [靶区各单元的电场向量及空气磁场](results/target_elements.csv) · [源项与解析线积分核查](results/source_check.csv) · [线性求解、组织电导率与内存记录](results/solver_run.json)

## 参考文献、头模和软件来源

- FieldTrip 团队。Subject01 原始 T1、头模构建坐标、分割与有限元网格。[数据说明与教程](https://www.fieldtriptoolbox.org/getting_started/othersoftware/simnibs/) · [可下载的头模](https://download.fieldtriptoolbox.org/getting_started/simnibs/m2m_Subject01/) · [个体解剖标签定义](https://download.fieldtriptoolbox.org/getting_started/simnibs/m2m_Subject01/segmentation/labeling_LUT.txt)
- Puonti O 等，2020。*Accurate and robust whole-head segmentation from magnetic resonance images for individualized head modeling*。[论文](https://doi.org/10.1016/j.neuroimage.2020.117044)
- SimNIBS 团队。有限元、磁场源项、电场单位及 ROI 的体积加权分析。[SimNIBS 4.6.0 固定版](https://github.com/simnibs/simnibs/releases/tag/v4.6.0) · [原生求解及场后处理](https://github.com/simnibs/simnibs/blob/v4.6.0/simnibs/simulation/fem.py) · [ROI 分析](https://simnibs.github.io/simnibs/build/html/tutorial/analysis.html) · [驱动与单位](https://simnibs.github.io/simnibs/build/html/faq.html)
- Worbs T 等，2025，提供 H1 重建几何并研究头皮贴合与形变对电场的影响。Deng Z-D 等，2013，比较线圈深度与聚焦性。[H1 来源](https://doi.org/10.1016/j.brs.2025.05.136) · [深度与聚焦研究](https://doi.org/10.1016/j.brs.2012.02.005)
