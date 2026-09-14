# field：密度、速度与输运

在给定 `ops` 和尺度上计算物理场及空间离散量；不拥有时间循环、配置默认值或文件输出。
导航：[包索引](../README.md) · [结构总览](../../STRUCTURE.md) ·
[第一象限架构记录](../../ARCHITECTURE_QUADRANT_LEVELSET_20260914.md)。

## 入口

- `initialDensity(ops,physics)`：内置初值或用户函数采样，返回 `ny × nx` 密度。
- `velocity`、`biotSavart`：由密度导数经椭圆求解得到速度、流函数及诊断信息。
- `poisson`、`poissonOperator`、`greenBoundary`：统一椭圆度量、边界值与求解路径。
- `transport(rho,u1,u2,ops,...)`：按已选 MUSCL、非均匀 WENO5、FD-WENO5 或 FD-WENO7 输运。
- `reconstructPhysical(rho,flow,scale,ops)`：统一重建物理坐标和场。

## 数值约定

物理壁面不透流；`open/closed` 控制人工边界。保守壁面模式与 `advective_upwind`
不可混同。`weno5_nonuniform` 的两层边界退回 MUSCL，不能据内部重构阶数声明整体四阶。
FD-WENO 使用映射、常量残差校正及适用的闭边界质量投影。

第一象限显式模式只积分正源并对 `X=0` 使用零流函数迹；WENO5-FD 水平 ghost
按偶密度/奇速度反射，以概念全域节点数与双分裂通量共同尺度复现旧正半轴 RHS。
Green 源压缩只保留正源，但沿用全域概念分箱编号。
`closedMassRateTarget` 默认零，非零目标仅供 WENO7 重标度装配使用。

各向异性椭圆度量是 `kappa=Cy/Cx`；每个 RK 阶段传入当前比例尺。
基线 `kappa=1` 的原有代数组合保持不变；高阶单边闭合使用 LU/direct。
`rho_physical=rho/Comega`；场始终与调用时的当前网格配套。
本模块只消费 `evolve` 传入的已解比例率并装配输运/源项；它不使用特征网格计数，
也不实现幅值、宽度或位移 restoring feedback。

模块通过传入的 `ops` 获取数值数据，不调用其他模块解析配置或推进状态。
相关验证：[椭圆](../../tests/+ipmtests/+baseline/elliptic.m)、[输运](../../tests/+ipmtests/+baseline/transport.m)、
[四阶](../../tests/+ipmtests/+fourth/core.m)、[六阶](../../tests/+ipmtests/+sixth/core.m)。
