# output：记录、结果与产物

拥有日志、v2 外部契约及可选可视化/保存，不推进状态，也不决定数值方法。
导航：[包索引](../README.md) · [结构总览](../../STRUCTURE.md)。
接受 [evolve](../+evolve/INTERNAL.md) 的已接受状态与当前网格；
checkpoint 恢复后交回同一演化入口。

## 接口

- `initializeLog`、`record`、`recordIfNew`、`shouldRecord`：统一记录时机与历史字段。
- `finalize(state,log,stopReason)`：物理量重建、质量摘要、可信窗拟合与结果校验。
- `validate(input)`：接受 v2 结构体或含变量 `result` 的 MAT 文件；`validateV2` 校验契约。
- `snapshotAt`、`trustedMask`、`compare`：读取匹配网格的快照、筛选历史与比较结果。
- `write`、`uniquePath`：保存 MAT 与追加清单；`report`：终端摘要。
- `makeCheckpoint`、`writeCheckpoint`、`readCheckpoint`：打包、原子安装并严格校验
  schema-v4 接受步事务。
- `maybeCheckpoint`：只在已记录末态仍属于连续可信前缀时执行 cadence/exit 写入。
- `restoreCheckpoint`：按白名单覆盖时间/步数终止界、`rhoXStop` 和输出控制，重建
  same-grid 完整状态。
- `checkpointFromResult`：从 terminal-trusted v2 结果重建并 overlap 复核一次性 checkpoint；
  不换网格，也不把指数化尺度冒充 bitwise 原值。
- `plotResult`、`visualize`、`closeLive`：结果图、在线绘图/视频与资源清理。

## 数据与写入约定

v2 包含 `metadata/state/grid/scale/physical/history/quality/fit/config/elliptic/snapshots`。
第一象限结果与原生 checkpoint 的轴/场只含 `X>=0`，并保留冻结的可选
`grid.quadrantOnly` 标记；不可将旧全域检查点裁半或修改 `xlim` 后续算。
历史分组及各时钟必须对齐；启用快照时，每帧保留自己的坐标，不能套用末态网格。
所有末态场为有限实数组且尺寸匹配。加载时不猜测变量名、不升级旧结果或注入旧字段。
历史 v2 结果的配置 schema 1/2/3 可原样校验；新运行配置为
schema 4，且固定 `scalingContract=exact_gauge_no_feedback_v1`。schema 3 保留五个
已退役的宽度反馈字段，但只能用于读取/分析，不会被隐式升级。
持久化校验同时要求 `lengthScaleGain=1` 以及
`omegaGaugeGain=widthGaugeGain=travelingWaveGain=0`；特征网格计数只能解释为
诊断/remesh/hard-stop 历史，不是可重放的变化率反馈。

三个固定锚点规范的 `history.gauge` 字段按实际运行名持久化：

- template：`omegaGaugeTemplateCenter`、`omegaGaugeTemplateRadius`、
  `omegaGaugeTemplateSupportPoints`、`omegaGaugeTemplateValue`、
  `omegaGaugeTemplateCondition`、`omegaGaugeTemplateRateResidual`、
  `omegaGaugeTemplateRateScale` 和 `omegaGaugeReferenceTemplateValue`；
- bulk：`omegaGaugeBulkCenterX`、`omegaGaugeBulkCenterY`、`omegaGaugeBulkRadius`、
  `omegaGaugeBulkSupportPointsX`、`omegaGaugeBulkSupportPointsY`、
  `omegaGaugeBulkValue`、`omegaGaugeBulkCondition`、`omegaGaugeBulkEnergy`、
  `omegaGaugeBulkRateResidual`、`omegaGaugeBulkRateScale` 和
  `omegaGaugeReferenceBulkValue`；
- wall L4：`omegaGaugeWindowCenter`、`omegaGaugeWindowRadius`、
  `omegaGaugeWindowValue`、`omegaGaugeWindowCondition`、`omegaGaugeWindowEnergy`、
  `omegaGaugeWindowRateResidual`、`omegaGaugeReferenceWindowValue`、
  `omegaGaugeWindowMomentOrder`、`omegaGaugeWindowMoment`、
  `omegaGaugeWindowSupportPoints` 和 `omegaGaugeWindowRateScale`。

三组均记录通用 `omegaGaugeEnergy` 和归一化 `omegaGaugeResidual`。
`gradient_energy` 还记录 `omegaGaugeEnergyX/Y`、
`omegaGaugeForcingX/Y`、总 `omegaGaugeForcing`、绝对 forcing 及其相消比。
所有动态各向同性规范同时记录
`cOmegaOverCL=c_omega/c_l` 和
`cLMinusCOmegaOverCL=(c_l-c_omega)/c_l`；分母在舍入尺度上不可分辨于零时写为
`NaN`，不以替代值参与演化或停机。

选中 `wall_omega_quadratic_peak` 时，`history.common` 还保存五个纯代数、
只读量：`physicalQuadraticPeak=(C_l/C_omega)P`、其倒数、
`quadraticInverseSlopeAlgebraic=-(c_l-c_omega)/P`、
`quadraticLocalTerminalTime=t+P/((c_l-c_omega)physicalQuadraticPeak)` 及
`rescaledRhoXInfToQuadraticPeak=G/P`。只有所选二次峰的 `P`、尺度和
有关分母有限且为正时才写入数值；未选该规范（包括 `none`
中的反事实峰）或条件不足时写 `NaN`。`finalize` 将可信前缀原样放入
`fit.quadraticPeakTelemetry`，并明示 `empiricalClassificationPerformed=false`；这些
恒等式量不执行回归，不产生终止时间、Type-I 或 blow-up 的经验声明。
续算时若旧 schema-4 checkpoint 前缀尚无这五字段，`finalize` 只在该前缀
的缺失/`NaN` 位置上，用已存 `C_l,C_omega,c_l,c_omega,physicalTime,G,P`
重建；新记录值优先且与重建值交叉核对。`source`、`storedRecords`、
`reconstructedRecords` 和 `storedReconstructionConsistent` 保留这一来源证据，
重建不回写历史。
checkpoint 运行期 reference 白名单对应包含
`referenceAnchorWallTemplate`、`referenceAnchorWallTemplateProjection`、
`referenceAnchorBulkGradientL2`、`referenceAnchorWallWindowValue` 和
`referenceOuterWallDensityWindowValue`；不得从结果诊断反推补写。

`saveResults=false` 时不写结果文件；绘图、在线图、视频有各自独立开关。
请求路径在 `config.output`，实际路径在 `metadata`；默认相对 MATLAB 当前工作目录。
MAT 先写临时文件再安装，失败仅清理由本次尝试生成的临时文件；
已有同名产物改用带 case ID 的路径，成功保存后追加 `manifest.jsonl`。
保存成功但清单追加失败会警告，不假装整次保存未发生。

checkpoint 与 result 分离。前者保存当前 `rho/x/y`、参考轴、重标度运行期 reference、内部尺度、
初值质量/值域、标量日志、可选历史快照以及输出/checkpoint 游标；因此
`history.wallCore` 在每个标量输出时刻保存以跟踪峰为中心的最多 129 个壁面节点：
重标度/物理/随动坐标、物理 `rho`、物理 `rho_x`、核心宽度和全局 `max|rho_x|`。
这是必备的一维轻量诊断，会进入 checkpoint 签名；它不开启或代替二维 `snapshots`。
旧 schema-4 checkpoint 缺少该可选组时仍可恢复，记录从恢复后的第一个输出时刻开始。

`storeSnapshots=false` 仍可恢复。写入先在目标目录生成临时 MAT，再以带 case ID 和步数的
新路径安装，并追加 `checkpoint_manifest.jsonl`，绝不覆盖请求的基名文件。

原生 checkpoint 的 same-grid split run 必须逐位复现。checkpoint schema 4 只接受
config schema 4；schema 1--3 checkpoint 属于旧演化契约，明确拒绝而不
隐式升级。terminal result 转换同样只接受 config schema 4；历史 config schema 1/2/3
result 仍可分析，但不能桥接到当前续算。可转换的 terminal result 先从冻结
配置重建初始状态，核对 `ops.rescaling/mass0/rhoRange0` 与 time-zero 记录，再重算末态
`psi/omega/u1/u2` 并核对全部可重建标量记录；terminal record 不可信或 overlap 超差即拒绝。
旧 result 只有 `Cx/Cy/Comega`，所以其内部 `logC` 只能由 `log(C)` 恢复，续算属于受控的
机器精度桥接。恢复时先从 schema-4 冻结配置重建 fresh ops，再只按显式白名单
恢复初始 gauge/网格分辨率 reference；索引由当前轴重算，配置派生的 `rescaling`
字段不从 payload 覆盖。恢复选项白名单只包括 `finalTime/physicalFinalTime/maxSteps`、产物/可视化控制、
`caseMetadata` 和 `checkpoint`；不允许数值、物理、网格、诊断或重标度覆盖。

`result.fit` 顶层继续保留旧 `blowupFit` 字段；`maximumGrowthRate` 和
`canonicalMaximumGrowth` 子结构分别在连续可信前缀的物理时间与 canonical `tau`
上拟合 `history.common.physicalRhoXInf`。子结构显式记录 observable、输入可用性和
可信记录数；这些分类只作证据，不会改变步长、规范或停止条件。
`gaugeRateConvergence` 子结构在同一连续可信前缀上记录
`c_l/c_omega` 及两个归一化比值的固定嵌套尾窗统计。其 1% 门只产生
`provisionalStable` 诊断，并以
`diagnostic_only_requires_box_time_grid_refinement` 明示要求独立的箱体、时间步和
网格加密证据；它不反馈演化、CFL 或停机决策。`finalize` 会从冻结配置
显式传入所选 `c_omega` 是否为 prescribed：物理模式、各向异性固定速率和
各向同性 `cOmegaGauge='none'` 都不具有经验收敛资格；其他各向同性代数规范
可进入新的零极限数值门，但仍需独立加密证据。

依赖 `config` 校验选项、`field` 重建物理量、`diagnostics` 生成诊断；checkpoint 恢复还
调用 `mesh.build` 与 `evolve.flow/initialize`，但不推进时间。
相关验证：[输出契约](../../tests/+ipmtests/+baseline/output.m)、
[checkpoint](../../tests/+ipmtests/+baseline/checkpoint.m)、
[集成](../../tests/+ipmtests/+baseline/integration.m)。

## 自动网格的持久化边界

可选autonomousMesh仅属于config schema4；历史schema1–3不能携带。无该字段的旧schema4不补默认。
validateV2用共享正规化器及递归类型/形状检查校验policy，再投影真实末态和完整history调用validateController。
readCheckpoint在签名路径校验账本/节点/时钟，restoreCheckpoint再与重建flow交叉验证。
启用version1时所有初选/迁移算子均明确来自custom_axes，恢复按保存轴直接构造一次；旧未启用分支保留。
metadata中的controller由既有完整payload签名覆盖；RHS cache仍不持久化。

version2配置N保持原始值，实际N/base必须匹配referenceFamily/currentLevelId；reader在LU前逐历史与快照核对级别。
仅v2的history.common新增nodeCountX/nodeCountY/meshLevelId。刚提交事务的lastDecision保留原source级别，不能改成target。
result→CP以ipm:CheckpointResultVariableNodeFamily拒绝未实现的族谱重建；原生CP支持按实际轴一次恢复。

显式初始观察回退策略及其used/unused证据随冻结配置和controller持久化。
读回仅验证保存的规则、观察轴、真实初选和完整账本；不得为校验重新求解析初值或重新规划。
旧配置省略该选项时保持原字段集合。实际H64四步与2+2步恢复的数学状态、完整history、
snapshots和cursor严格一致；结构体字段排列顺序不属于数值差异，但字段集合、类型、形状、
每个值及NaN位置必须相同。验证证据见研究目录中的初始回退合同。
