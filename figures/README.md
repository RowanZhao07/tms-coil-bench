# 研究图件与重绘入口

四个图册共 **39 幅**，每幅均提供浏览用 PNG 和排版用 PDF。打开下表的图册，即可在每幅图下看到中文解读、单位与数据入口。不同线圈的名称、左右方向、磁场分量和距离起点均沿用相应实验定义。

| 图册 | 幅数 | 内容 |
|---|---:|---|
| [绕组、摆位与求解器](core/README.md) | 22 | 8 幅汇总图和 14 幅逐线圈图。几何投影、轴上与截面磁场、深度归一化、摆位规则、有限元收敛和参数敏感性。 |
| [公开台架测量](measurements/README.md) | 8 | 三轴扫描、拾取脉冲、霍尔距离曲线、跨电流预测、圆形线圈绝对场、阵列定位和盐水体模电场。 |
| [MRI 配准、反演与消融](methods/README.md) | 7 | 同色标场图、绕组与切片位置、剖线残差、六点检验、轴向修正、三轴选点和跨距离校准。 |
| [左丘脑靶区与头内电场](brain/README.md) | 2 | 三种线圈的空气磁场和靶区电场排序，个体 T1 解剖与三个电场的同位置冠状切面。 |

## 取图与读图

投稿排版可直接选 PDF。文字、坐标轴、标记和曲线保留为矢量对象。方法图册保留体模实测的像素网格，头内电场图保留原四面体单元的场值。宽约 180–183 mm 的横版图适合双栏，方法图册另有 89 mm 单栏剖线。600 dpi PNG 用于文档插图或需要位图的工作流。原始采样网格决定科学分辨率，导出更多像素不会产生新的测量。

图中的“相对 L2 误差”是 $\varepsilon=100\|p-y\|_2/\|y\|_2$，其中 $p$ 为预测，$y$ 为对应的实测或数值参照。若曲线在弱场区看似接近，仍应检查误差旁的物理单位；同样的绝对偏差除以较小的实测值，会得到较大的百分比。

## 从仓库数据重新绘制

在仓库根目录运行以下 PowerShell 命令。Windows 和 Python 3.13 下已验证。绘图只需 Python；MRI 切片矩阵、逐点预测、几何与汇总表随库提供。丹麦拾取波形从作者固定版本获取，准备脚本会核验文件后缓存。

~~~powershell
python -m venv .venv
& .\.venv\Scripts\python.exe -m pip install -r figures/requirements.txt
& .\.venv\Scripts\python.exe external_measurements/measured_coil_library/prepare_inputs.py
& .\.venv\Scripts\python.exe figures/core/render_core.py
& .\.venv\Scripts\python.exe figures/measurements/render.py
& .\.venv\Scripts\python.exe figures/methods/draw_methods.py
& .\.venv\Scripts\python.exe brain_field/plot_results.py
~~~

各脚本按自身位置解析输入和输出目录。需要重新解场、拟合或生成表格时，使用[主 README 的实验复跑命令](../README.md#reproduce)。

## 参考与数据来源

绘图的原始输入、测量协议和作者项目见[主文数据与文献目录](../README.md#sources)。排版尺寸参考 [Nature 图件制作指南](https://www.nature.com/documents/Final_guide_to_authors.pdf)。
