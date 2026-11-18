# SimNIBS H1 参考几何

本目录保留了 H1 官方模型、解释其线元的上游参考实现，以及本项目导出的 SI 单位线元表。参考实现中一个线元由中点及有向增量给出。若中点是 $p$、增量是 $v$，两个端点分别为 $p-v/2$ 和 $p+v/2$。导出脚本不改变模型的形变设置；输入文件的 SHA-256、分组、长度和闭合残差见[来源与几何记录](h1_metadata.json)。上游文件及其适用许可见 [SimNIBS 许可证](SIMNIBS_LICENSE.txt)。

## 上游文件和版本

以下文件取自 SimNIBS 提交 `fde509f6f8de38a4b03947944cf7ee93b8bebabe`。本地参考实现另存为 `tms_coil_element_reference.py`，仅用于核对线元约定。

- [Brainsway_H1.tcd](https://github.com/simnibs/simnibs/blob/fde509f6f8de38a4b03947944cf7ee93b8bebabe/simnibs/resources/coil_models/flexible_coils/Brainsway_H1.tcd)
- [线圈文档 coils.rst](https://github.com/simnibs/simnibs/blob/fde509f6f8de38a4b03947944cf7ee93b8bebabe/docs/documentation/coils.rst)
- [格式定义 coil_model.schema.json](https://github.com/simnibs/simnibs/blob/fde509f6f8de38a4b03947944cf7ee93b8bebabe/simnibs/resources/coil_models/coil_model.schema.json)
- [线元实现 tms_coil_element.py](https://github.com/simnibs/simnibs/blob/fde509f6f8de38a4b03947944cf7ee93b8bebabe/simnibs/simulation/tms_coil/tms_coil_element.py)
