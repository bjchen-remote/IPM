# config：配置解析

将用户的平铺 `opts` 变成唯一一份冻结分组配置，不构造场、不推进时间。
总体契约见 [STRUCTURE.md](../../STRUCTURE.md)。

## 接口

- `schema()`：选项名、分组、原始默认值、类型及合法取值的唯一清单。
- `resolve(opts)`：覆盖、规范化、推导、校验，返回 `schemaVersion=4`、`frozen=true` 的配置。
- `activeCase()`：无参数求解时使用的生产覆盖项；不等于空结构体的配置默认值。
- `longTimeProfile(settings)`：把服务器长跑的少量公开设置映射为显式 version-4
  自动网格选项；只解析配置，不创建目录或启动求解。
- `sixthOrder(overrides)`：构造完整六阶覆盖项；仍由 `resolve` 作最终校验。
- `tupleViolation(...)`：统一检查耦合数值选择，供解析与持久化结果校验共用。
- `outputPaths(opts)`：推导默认结果和视频路径，不创建文件。

## 约定

内部配置分为 `grid/time/physics/elliptic/transport/scaling/remesh/diagnostics/output` 九组。
未知字段报错；`customX/customY` 必须成对。求解初始化只调用一次 `resolve`，
重网格直接传入局部网格覆盖，不重新解析。不要向内部代码传第二份平铺选项。

`checkpoint` 是无默认值的可选 output 策略；出现时规范化为
`enabled/file/every/atExit`。`every` 是正的 canonical 时间间隔或 `Inf`；省略文件名时
从 `resultFile` 派生。未显式提供该选项的旧配置和 v2 结果不会被注入新字段。

`autonomousMesh` 是无默认值的可选 remesh 策略；未显式提供时旧配置逐值不变。
`autonomousMeshPolicy(input,choices)` 只做纯参数正规化，保存 version=1、enabled
（默认 true）、双向 targetCoreCells 和门槛、固定 qualityLimits/search 以及 native
canonical 预测窗。第二参数可选；提供时统一核对允许的数值元组，供解析和持久化校验共用。
目标为 [32,256] 内的两个整数，事务最低单元数为 target-1，触发/预测缓冲/终点
分别为 target×26/32、22/32、20/32；完整显式派生字段重读必须与目标逐值一致。
历史 core14、安全 .70、峰跳/守恒/范围、网格质量和搜索范围固定为已注册值；
目标允许不等两轴，但不保证给定节点预算可行。enabled 只允许 dynamic/isotropic/
double_odd_omega、transport_anchor 与 wall_omega_quadratic_peak、高阶/WENO5-FD/
SSPRK54/high_order 迁移，并要求 adaptiveRemesh=true、initialAnalyticRemesh=true、
adaptiveLevels=[.1 .5 .9]；disabled 不增加这些元组限制。它不更改 C、CFL 或 maxDt，
配置正规化本身不证明运行期自动网格已可用。

version5额外冻结 `nodeFamily.initialNodeCount` 和启动时的 `maximumTotalNodes`，由两者
确定方向因子/节点级，记录动态成员数的有界搜索身份；旧version1--4正规化不变。

`initialMeshObservationFallback` 是另一项无默认值的可选 remesh 策略，
由纯函数 `initialMeshObservationPolicy` 正规化。显式启用的 version1 固定为
`primitive_k8_fixed_probe_v1`：仅支持原始 `degenerate_primitive`、power8、
anchor1、double-odd 及已启用的 autonomousMesh v2/v3/v4/v5，要求对称箱 H>4、Ymax>4。
内区半宽4、间距.025、尾比1.05、20格平滑坡、最多600000观察节点及一次回退均由版本固定。
省略该选项不补字段、不改变旧配置；它不更改实际求解节点预算、C或时间步规则。

配置 schema 4 固定
`scalingContract='exact_gauge_no_feedback_v1'`。`maxDynamicRate` 已删除；
`adaptiveGain` 运行组件与五个 width-controller 选项
`adaptiveLengthScaling`、`widthExpansionStrength`、`widthContractionOnset`、
`widthContractionStrength`、`maxWidthRateCorrection` 也已退役。
`lengthScaleGain` 必须等于 1，`omegaGaugeGain`、
`widthGaugeGain`、`travelingWaveGain` 必须等于 0。这些不是可调增益，
而是持久化的 no-feedback 约束。特征网格单元计数仍可配置为诊断、remesh 和
hard-stop 门，不是比例率控制器输入。

旧 `anchor_wall_strain` 只是一条近似历史关系，schema 4 明确拒绝运行；
历史 schema 1--3 结果仍可按原字段只读。离散壁面模板规范当前只允许固定网格，
在未注册模板重采样规则前不与 adaptive remesh 混用。
幅值规范枚举包含三个固定锚点窗比较模式：
`anchor_wall_template_projection`、`anchor_bulk_gradient_l2` 和
`anchor_wall_window_l4`。它们与 L2 窗共用
`transportAnchorX=1` 和 `omegaGaugeWindowRadius`；解析器不从网格计数推导或移动窗。

外部 result-v2 读取器仍可原样校验历史 config schema 1--3；旧结构不注入
新字段、不升级，也不能用于 schema-4 checkpoint 续算。

基线保留 `legacy_second_order/muscl_minmod/ssprk3/pchip`。
六阶必须同时选择 `sixth_order/weno7_fd/rk6/sixth_order`；四阶 FD-WENO 与高阶迁移
分别要求 `high_order` 空间后端。两种高阶后端的各向异性求解均禁止 PCG。

唯一跨模块依赖是 `mesh.sixthOrderPolicy()` 的准入常量与推荐 CFL；
不得在本模块运行算子装配或研究脚本。六阶网格阈值不是用户可调的公共选项。
相关验证：[基线配置](../../tests/+ipmtests/+baseline/config.m)、[六阶核心](../../tests/+ipmtests/+sixth/core.m)。
