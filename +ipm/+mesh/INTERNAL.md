# mesh：网格与运行算子

把冻结配置和可选网格覆盖装配为完整 `ops`，统一几何、导数、求积与椭圆离散数据。
导航：[包索引](../README.md) · [结构总览](../../STRUCTURE.md) ·
[第一象限架构记录](../../ARCHITECTURE_QUADRANT_LEVELSET_20260914.md)。

## 接口与数据

- `build(config,gridOverride)`：构造完整运行算子；覆盖仅替换网格分组，不改配置。
- `fdMatrix(axis,derivativeOrder,stencilWidth)`：局部多项式差分矩阵。
- `quadrature`、`quadrature6`：配套高阶权重；`quality`：网格质量诊断。
- `sixthOrderPolicy()`：固定的六阶网格准入阈值及推荐 CFL。

`ops.x` 是行向量，`ops.y` 是列向量，场为 `ny × nx`。
`ops` 保存数值运行数据，不包含整份配置；`integrationWeights` 统一用于积分和守恒校正。
基线导数和控制体积权重保持原实现；四、六阶分别使用七点、九点导数与配套求积。
高阶单边闭合不对称，不能复用基线的 SPD 假设。

象限选项仅保存非负 `x`：`quadrantDerivative` 在局部虚拟模板上折叠奇偶系数，
`Dx` 作用于偶密度，`oddDx` 作用于奇流函数/水平速度；高阶 `Tx` 对奇流函数装配。
`quadrantQuadrature` 是和旧全域正半域一致的单侧权重；
非均匀 `X(ξ)` 的逻辑度量由奇模板求导，避免对称轴单边闭合改变 WENO 物理间距。
Poisson 内部未知数仍为 `(nx-2)×(ny-2)`。

## 六阶准入

二维稀疏矩阵装配及分解前，每轴必须满足：相邻单元比 `<=1.15`、
九点归一化模板 `rcond>=1e-9`、最小权重/平均权重 `>=1e-4`、度量误差 `<=1e-6`。
自动拉伸受限；显式、定制或重网格后的不合格坐标直接拒绝，不静默修补。
这些门槛不构成任意非均匀网格的稳定性证明。
网格质量与特征单元计数只是网格准入、诊断、remesh 和 hard-stop 数据；
mesh 不实现 width/scaling feedback，也不修改精确规范返回的比例率。

本模块仅跨模块调用 [field](../+field/INTERNAL.md) 的 `weno5Geometry`；椭圆矩阵及缓存由 `build` 直接装配。
运行时的度量算子由 `field.poissonOperator` 处理；候选网格由
[remesh](../+remesh/INTERNAL.md) 提交给 `build`，本模块不调用推进、输出或研究代码。
相关验证：[基线网格](../../tests/+ipmtests/+baseline/grid.m)、[四阶核心](../../tests/+ipmtests/+fourth/core.m)、
[六阶核心](../../tests/+ipmtests/+sixth/core.m)。
