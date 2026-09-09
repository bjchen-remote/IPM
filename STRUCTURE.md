# 结构与接口

本副本按职责维护原求解器的命名空间、数值路径和版本化契约。
`ipm.solve(opts)` 是唯一新算例入口；`ipm.solve(overrides,checkpoint)` 从同网格可信
接受步恢复；`ipm.verify(...)` 只负责验证调度。

## 目录地图

```text
ipm_structured/
├── README.md                 使用入口
├── STRUCTURE.md              本文件：边界与契约
├── CHANGELOG.md              本副本的变更及验证记录
├── +ipm/
│   ├── solve.m               求解流程
│   ├── verify.m              验证入口
│   ├── +config/              选项、默认值、校验
│   ├── +mesh/                网格、导数、求积、算子装配
│   ├── +field/               密度、椭圆求解、速度、输运
│   ├── +evolve/              状态、重标度、时间推进
│   ├── +remesh/              候选网格、迁移、验收
│   ├── +diagnostics/         特征、质量、停止条件
│   └── +output/              记录、v2 结果、保存、绘图
├── examples/                 短小、非写入算例
├── tests/                    分层验证、新旧一致性与来源清单
├── research/                 实验、分析与历史资料
└── result/                   本副本的生成输出
```

每个运行模块各有一份 `INTERNAL.md`。`+ipm` 的父目录是 MATLAB 路径入口；
不复制第二套扁平 `main_ipm` 接口，不依靠递归 `addpath` 拼装运行环境。

## 模块职责与依赖

| 模块 | 主要入口 | 直接使用的其他模块 |
|---|---|---|
| [config](+ipm/+config/INTERNAL.md) | `resolve`, `schema`, `sixthOrder` | `mesh.sixthOrderPolicy` 的常量 |
| [mesh](+ipm/+mesh/INTERNAL.md) | `build`, `fdMatrix`, `quality` | `field.weno5Geometry` |
| [field](+ipm/+field/INTERNAL.md) | `initialDensity`, `velocity`, `transport` | 无；使用传入的 `ops` |
| [diagnostics](+ipm/+diagnostics/INTERNAL.md) | `measure`, `trackFeatures`, `stopPolicy` | 无；计算输入数据的诊断量 |
| [remesh](+ipm/+remesh/INTERNAL.md) | `adapt` | `mesh`, `diagnostics` |
| [output](+ipm/+output/INTERNAL.md) | `record`, `finalize`, `validate`, `makeCheckpoint`, `restoreCheckpoint` | `config`, `mesh`, `field`, `evolve`, `diagnostics` |
| [evolve](+ipm/+evolve/INTERNAL.md) | `initialize`, `advance`, `flow` | 以上六个模块；统一协调状态 |

数值运行模块不直接引用 `tests/`、`research/` 或原版目录；研究和测试可以调用运行模块。
自定义初值通过 `initialCondition` 函数句柄注入，不把研究路径写进核心。
只有显式一致性验证读取原版。模块内部函数虽可由 MATLAB 包名访问，仍属于内部接口；
用户通常只需要 `ipm.solve`、配置辅助函数以及结果读取/绘图函数。

## 一次计算的路径

```text
平铺 opts（或 config.activeCase）
  -> evolve.initialize
       -> config.resolve       只解析一次
       -> mesh.build -> field.initialDensity -> 初始化重标度与流场
  -> solve 循环
       -> evolve.advance -> 选步长 -> 配套 RK -> remesh.adapt
       -> output.record -> diagnostics.stopPolicy
       -> output.maybeCheckpoint（仅连续可信前缀中的已记录接受步）
  -> output.finalize -> 物理量重建 -> v2 校验
  -> output.write / report / plotResult

可信 checkpoint -> output.readCheckpoint -> output.restoreCheckpoint
  -> 重建同一网格算子与流场，恢复完整 state/log/output cursor -> solve 循环
```

`solve` 只组织流程。物理模式、各向同性和各向异性动态模式共用数值路径，
各向同性是 `Cx=Cy` 的特例。步进及重网格结果经有限性检查后才保留；
非有限状态回退至上一个已接受状态，返回 `non_finite_solution`；数值模块异常继续抛出。

## 四种数据对象

| 对象 | 所有者 | 契约 |
|---|---|---|
| `config` | `config.resolve` | 完整、分组、`frozen=true`，新配置 `schemaVersion=4` |
| `ops` | `mesh.build` | 当前网格及完整数值算子，不存整份配置 |
| `state` | `evolve` | 已接受的 `rho/flow/scale/ops`、时钟、步数及冻结配置 |
| `result` | `output.finalize` | 外部 v2 数据结构，不保存运行中的求解对象 |
| `checkpoint` | `output` | schema v4 接受步事务；保存可恢复 state、log 和调度 cursor |

场数组均为 `ny × nx`，行对应 `y`，列对应 `x`。当前网格坐标与场必须成对。
重网格通过 `mesh.build(config,gridOverride)` 构建候选算子，不重新解析配置。
`ops.integrationWeights` 是积分、守恒检查与迁移校正的共同权重。

外部 `opts` 为平铺标量结构体；内部配置分为 `grid/time/physics/elliptic/transport/`
`scaling/remesh/diagnostics/output` 九组。选项全集只在
[config.schema](+ipm/+config/schema.m) 定义；未知字段和不兼容组合报错。
`config.activeCase` 是无参数生产算例，不能与配置表默认值混同。

config schema 4 固定
`scalingContract='exact_gauge_no_feedback_v1'`。`lengthScaleGain` 只能为 1，
`omegaGaugeGain`、`widthGaugeGain` 和 `travelingWaveGain` 只能为 0；
`maxDynamicRate`、`adaptiveGain` 组件及五个 width-controller 选项
(`adaptiveLengthScaling`、`widthExpansionStrength`、`widthContractionOnset`、
`widthContractionStrength`、`maxWidthRateCorrection`) 不再属于当前契约。
可选`remesh.autonomousMesh`仅显式启用，不给旧配置注入默认；省略版本仍是固定节点数/箱的version1。
显式version2允许在零时刻冻结的方向增点参考族中自动迁移，并预先登记节点资源上限；
已实际从原始t=0跨过两个方向的自然增长，但有限箱、有限族的运行尚不构成长时误差验收。
配置Nx/Ny保留初始预算，当前实际节点数、参考族成员与各时刻场/网格配对单独记录。
可选 `initialMeshObservationFallback` 为原始k8初值提供一次自动解析重测，
仅原初选因候选族耗尽或核心/前沿不足失败时使用；观察网格与实际演化网格分开记录。
该选项不注入旧配置，恢复不重测或重采样初值；实际网格仍经过原质量与迁移门。
动态比例率由所选规范的瞬时代数恒等式唯一给出；
参考量只用于精确条件、定向、正性和病态拒绝，不生成 restoring 源项。
特征的网格单元计数仍可用于 telemetry、remesh 提案/验收和 hard stop，但不能修改
`c_l`、`c_omega`、`c_r` 或时间速度。

四个数值选择分别为 `spatialDiscretization`、`transportScheme`、`timeIntegrator`、
`remeshTransferScheme`。基线是 `legacy_second_order/muscl_minmod/ssprk3/pchip`；
光滑四阶路径显式选择 `high_order/weno5_fd/ssprk54`，整体阶数证据不包含触发式迁移。
六阶严格要求 `sixth_order/weno7_fd/rk6/sixth_order` 四者配套。
两种高阶空间后端采用非对称闭合；各向异性求解只能选择 direct，不能使用 PCG。
网格准入规则与证据范围见 [mesh/INTERNAL.md](+ipm/+mesh/INTERNAL.md) 和
[研究索引](research/README.md)，合法配置不自动获得整体阶数保证。

## 对外结果

`result.schemaVersion=2`，包含 `metadata/state/grid/scale/physical/history/quality/fit/`
`config/elliptic/snapshots`。`state.rho` 是当前演化场，`physical.rho` 是物理密度。
`state.normalizedTime`、`state.canonicalTime`、`state.physicalTime` 不可混用；
对应历史字段是 `history.common.t/canonicalTau/physicalTime`。
`history` 分为 `common/gauge/mesh/anisotropic`，同组记录对齐；
启用快照时，每个 `snapshots.rho{k}` 与自己的 `x{k}/y{k}` 和时刻配套。

物理坐标重建为 `physical.x=(ops.x-Xshift)/Cx`、`physical.y=ops.y/Cy`，
物理密度为演化密度除以 `Comega`。正比例尺在状态中用对数保存；
`Lx=1/Cx`、`Ly=1/Cy`、`aspect=Cy/Cx`；
守恒源系数为 `cx+cy+comega`。重标度量与物理量的重建统一由 `field` 完成。

`output.validate` 接受 v2 结构体或含变量 `result` 的 MAT 文件；不猜测旧字段或转换 v1。
较早 v2 结果内的配置 schema 1--3 仍可原样读取分析，其缺少的选择只按旧语义校验、不补写。
请求的输出路径保留在 `config.output`，实际路径在 `metadata`；默认相对当前工作目录。
保存使用新文件名并追加清单，不覆盖旧产物。可信数据窗口和拟合结果只是数值诊断，
不能替代网格加密、边界截断检查或奇性证明。

checkpoint 不是普通 result：它额外保存 `baseX/baseY`、完整 `ops.rescaling`、内部对数
尺度、`mass0/rhoRange0`、当前 `rho/x/y`、标量历史及 `nextOutput/nextCheckpoint`。
临时 MAT 在目标目录写完后才安装为带步数的不可变文件；末条记录不在连续可信前缀时不写。
恢复白名单不包含网格、物理、数值离散、诊断阈值或重标度选项。已有 terminal v2 result
只能通过 `checkpointFromResult` 重建同网格 checkpoint，并须通过初始记录、尺度、质量、
值域和末态流场 overlap 检查；它不能补回旧 result 未保存的 bitwise 对数尺度。
当前 checkpoint/config 版本均为 4，且续算时必须重放同一
`exact_gauge_no_feedback_v1` 运行期契约。schema 3 及更早产物只读分析，不可续算或隐式升级。
autonomousMesh version2的原生检查点保存并验证实际参考族、级别、方向增点账本、触发证据和
逐记录节点数；支持从保存的实际网格恢复。该版本暂不支持result→checkpoint桥接，明确拒绝，
不能从普通结果猜补变量网格运行状态。

## 文档与历史

当前行为以本文件、模块文档和实现为准。[tests/source_map.json](tests/source_map.json)
记录完整迁移对应；[research/notes](research/notes) 原样保存原版技术文档。
历史文件中的路径、函数名和研究计划按历史上下文阅读，不作为当前执行指令。
