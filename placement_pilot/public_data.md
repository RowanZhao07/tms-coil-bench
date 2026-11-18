# 用于统一摆位与非轴对称验证的公开数据

扩展实验已经复跑德国和法国的距离曲线、丹麦 MC-B70 的三轴扫描与 25 款线圈的拾取波形，以及以色列 H1 盐水体模的电场剖面。美国微型八字、意大利圆形与芬兰/巴西阵列的测量也已纳入。下表保留与本节统一摆位最直接相关的来源。[跨实验室数据、资料形态和适用范围](../external_measurements/README.md)

截至 2026-11-18，可以从公开材料复核一只 8 字线圈的体模实测磁场。H1、8 字与双锥的统一摆位也有公开头模及线圈模型可用。三个型号在同一头模、同一摆位协议下的三分量实测磁场，本次检索未找到可直接下载的成套数据。

| 来源 | 已核实的对象和数据 | 可支持的检验 |
|---|---|---|
| Petrov 等，2017 | Magstim Rapid² 配 MR 兼容的八字线圈、琼脂体模和 3 T MRI。公开 MAT 文件含五个 160 × 160 的矩阵。由 MRI 相位推算的 `data_MRI1` 是一个有符号的磁场分量 Bz 切片，另四个矩阵来自模型计算。| 可以逐像素检验模拟的 Bz 场形和绝对幅值。作者在 v1.1 修复了用于摆放线圈的神经导航配准矩阵，测量与预测矩阵的文件内容沿用 v1.0。[论文与等效电流](https://doi.org/10.1371/journal.pone.0178952) · [测量和对照矩阵](https://github.com/pip010/plosone-data-pub01/blob/v1.1/phantom_experiments/error_metrics/data.mat) · [作者的摆位修复说明](https://github.com/pip010/plosone-data-pub01/releases/tag/v1.1) |
| Drakaki 等，2022 | 25 款商业线圈的模型来自直接场测量。SimNIBS 仓库可见 `MagStim_D70.ccd`、`MagStim_DCC.ccd`、`MagVenture_C-B60.ccd`、`MagVenture_Cool-D-B80.ccd` 等。以 DCC 为例，CCD 文件存有 2,000 个等效磁偶极子的中心和加权方向，并给出刺激器信息。| 可用匹配型号的独立等效源预测空间场，检验自行建立的绕组几何。该对照保留了测量约束，但文件中的数据行是模型参数，不能作为 2,000 个探头观测。D-B80、Cool-D-B80 和自建理想双锥应分别标明。[论文](https://doi.org/10.1016/j.brs.2022.04.017) · [官方模型目录](https://github.com/simnibs/simnibs/tree/master/simnibs/resources/coil_models/Drakaki_BrainStim_2022) · [DCC 文件](https://raw.githubusercontent.com/simnibs/simnibs/master/simnibs/resources/coil_models/Drakaki_BrainStim_2022/MagStim_DCC.ccd) |
| Schecklmann 等，2020 | Figure 2 给出 MagVenture C-B60 8 字与 D-B80 双锥的峰值磁场距离曲线。横向霍尔探头沿中心法线以 1 mm 步距移动，刺激器使用 100% 单相脉冲。探头只读其敏感方向的分量。| 可在对准传感面后复核衰减趋势。图上的点属于发表的曲线。若数字化并使用这些点，需要记录读图误差。作者的数据声明为依申请获取，未找到公开原始读数下载。[论文与图 2、3](https://doi.org/10.1038/s41598-020-58034-2) |
| Roth 等，2007 | 在盐水头形体模内用双电极探头测得 H1、H2 和 8 字线圈的感应 **E**，论文给出空间分布。| 可检查 H1 同类设计的深度与扩散趋势。比较时复现盐水电导边界、各线圈姿态和论文归一化，再与相应的 E 分量或幅值比较。本文未提供可直接逐点复算的 H1 三轴 **B** 数据。[原论文](https://www.brainsway.com/wp-content/uploads/2018/12/Roth_Journal_of_Clinical_Neurophysiology_2007_Fig-8_H-coil_Electrical_Field_Distribution.pdf) |
| Worbs 等，2025；SimNIBS | 柔性 H1 的 `.tcd` 提供 X 射线重建的绕组、外壳与可动部件，便于保存刚体变换与形变参数；官方示例数据提供同一套 Ernie 真实头模、MNI152 和球头模。| 可用同一网格、头皮锚点和外壳净距重复摆位。H1 外壳模型和头部组织网格支持建模与复现，尚不能验证探头实测场强。作者已经研究过大型柔性线圈的自动位置与形变优化，因此“统一摆位”宜作为比较控制条件。[论文](https://doi.org/10.1016/j.brs.2025.05.136) · [项目使用的 H1 文件版本](https://github.com/simnibs/simnibs/blob/fde509f6f8de38a4b03947944cf7ee93b8bebabe/simnibs/resources/coil_models/flexible_coils/Brainsway_H1.tcd) · [官方头模说明与下载](https://simnibs.github.io/simnibs/build/html/dataset.html) |

## 将实测切片用于独立验证

作者仓库的两个版本均可下载施加 TMS 和关闭 TMS 时的刺激实验 NIfTI 图像，以及供 SCIRun 读取的 NRRD 文件。上游目录说明介绍了 Philips PAR/REC 转换流程，两个版本的文件清单中提供的是转换后的图像。本实验直接采用作者已处理成特斯拉的 `data_MRI1`。从成像相位到物理量的过程及脉冲时序由原研究给出。[公开成像文件清单](https://github.com/pip010/plosone-data-pub01/tree/v1.1/phantom_experiments/mri_data/PHILIPS_3T_Achieva_RAW_FILES) · [原始测量方法](https://doi.org/10.1371/journal.pone.0178952)

本项目用作者给定的绕组、等效电流和线圈姿态独立计算磁场，将三维磁场投影到扫描仪轴向，再在相同切片上取样。比较遵循作者 MATLAB 脚本的数组朝向和体模掩膜，共用 10,622 个有效像素。对于单层九匝模型，修复配准后，场强均方根与实测的比值由 4.5018 回到 1.0162，相对 L2 误差从 352.56% 降到 7.63%。两次计算使用相同的电流、网格和评价像素。[项目图解与重跑指南](../README.md) · [八字线圈逐模型对照](public_benchmark/matched_figure8/placement_version_comparison.csv) · [上游修正提交](https://github.com/pip010/plosone-data-pub01/commit/4cb2218)

公开验证这个 8 字参考件可以增加对同一求解流程处理非轴对称几何的信心。本项目的理想 8 字或儿童 8 字若几何、尺寸及回流导线不同，需与该参考件区分。对 H1、双锥，型号或形状未充分对齐时，适合报告间距、角度与回流导线差异对场形的影响；人体内总感应电场仍需头模型及其导电边界。
