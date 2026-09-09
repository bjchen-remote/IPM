# 分层验证

添加本副本根目录后，统一通过 `ipm.verify(suite,mode,originalRoot)` 调度。
`suite` 默认为 `baseline`，`mode` 默认为 `quick`；无需手动添加测试包目录。

| 调用 | 覆盖 |
|---|---|
| `ipm.verify()` | 配置、网格、椭圆、输运、对称、重标度、重网格、输出、checkpoint 与集成基线 |
| `ipm.verify('fourth')` | 四阶空间/时间/物理收敛、Green 边界、迁移与 quick 稳定性扫描 |
| `ipm.verify('fourth','heavy')` | 四阶套件，扩展稳定性扫描 |
| `ipm.verify('sixth')` | 六阶空间、WENO7、RK6、网格准入、迁移与物理收敛 |
| `ipm.verify('all')` | 独立目录结构检查，以及基线、四阶、六阶；不需要原版 |
| `ipm.verify('all','heavy')` | 同上，其中四阶采用 heavy 稳定性扫描 |
| `ipm.verify('equivalence')` | 显式比较本副本与保留的原求解器 |

`quick/heavy` 仅改变四阶稳定性扫描范围，不降低其他套件的通过标准。
`all` 不隐含新旧一致性检查，保证本副本脱离原目录后仍可独立验证。

一致性检查默认读取相邻 `sixth_order_integration/`，也可以指定其他原版位置：

```matlab
ipm.verify('equivalence','quick','/absolute/path/to/sixth_order_integration');
```

该检查比较 6 个纯物理模式短程求解：二、四、六阶数值元组分别搭配
`double_odd_omega` 与 `half_plane` 对称性；另比较 PCHIP、四阶、六阶三组非平凡迁移。
配置的网格、时间、物理、椭圆、输运、重网格、诊断和输出域仍逐域精确比较；
求解结果剔除运行元数据、版本化配置和独立演进的拟合诊断。此外，
schema 4 直接以 `tau+dt` 更新规范时钟，而保留旧版把常数导数 `1`
经 RK 表重组；纯物理模式的轨迹未变，但舍入可在三个重复的规范时钟字段
`state.canonicalTime`、`history.common.canonicalTau` 和
`snapshots.canonicalTime` 中不同；物理内核投影仅额外剔除这三项。
`state.rho/psi/omega/velocity`、`state.physicalTime`、整个 `physical`
重建、`history.common.physicalTime`、其余全部历史物理观测量、所有尺度、
normalized clock、质量、快照和椭圆信息仍全部用 `isequaln` 比较。

2026-09-09 的首次生命周期后新旧检查，保留在
`result/verification/autonomous_runtime_20260909/lifecycle_and_suites_v1/equivalence_failure.json`，
因旧版缺少此前新增的 23 个 `history.common` 字段而失败；这次失败证据不改写。
新版测试先要求新增字段集合**恰好**等于以下白名单，反向缺失集合为空，
再仅从新版投影剔除这些旧版不存在的字段。未知新增项、缺失白名单项或缺失旧字段都拒绝，
不以字段交集隐式扩大范围：

- 步号：`acceptedStep`。
- 上一接受步的选择记录：`acceptedCanonicalDt`、`transportRate`、`cflDtLimit`、
  `maximumDtLimit`、`canonicalEndpointDtLimit`、`physicalEndpointDtLimit`、
  `physicalClockSpeed`、`realizedCfl`、`timestepActiveLimiter`、
  `amplitudeSourceStep`、`conservativeSourceStep`。
- 物理网格覆盖：`physicalDomainXMinimum`、`physicalDomainXMaximum`、
  `physicalDomainXRadius`、`physicalDomainYMaximum`、`physicalMinimumDx`、`physicalMinimumDy`。
- 所选二次峰的代数量：`physicalQuadraticPeak`、`inversePhysicalQuadraticPeak`、
  `quadraticInverseSlopeAlgebraic`、`quadraticLocalTerminalTime`、
  `rescaledRhoXInfToQuadraticPeak`。

这些字段由 `output.record` 追加：步长组复制 `selectTimestep/emptyTimestep` 的观测，
域/间距组来自轴与尺度，二次峰组来自局部只读代数函数。新增记录不修改本检查中
未启用自动网格的物理轨迹；`acceptedStep` 等日志也可供新控制器审计使用，
不能据此将启用控制器的运行解释为旧版算法。

每个算例还要求快照及历史覆盖每一个接受步，以旧版快照、旧版 flow/ops
独立重算下一行的 transport rate、四个步长上限、实际选择及 limiter，
并与旧版 `ipm_select_timestep` 的返回 dt 逐值比较。
物理域/间距和源项步幅从旧版轴、共有尺度与共有规范率重算；初态 NaN、
limiter 字符串、Inf 均保持 `isequaln` 语义。
独立按 dt 直接累加还校验新版三个规范时钟副本。
本套件六个物理算例未选二次峰规范，因此五个新增峰值字段须精确为 NaN；
动态有效代数值由 `ipmtests.baseline.output` 单独验证，这里不扩大其覆盖声明。
测试包含五个故意损坏的投影/派生量拒绝检查，所有共有字段仍使用原始 `isequaln`，
没有新增数值容差。

迁移比较场、指标、非规范化数值算子及分解后的求解作用。新版在
`operator.gridQuality.x/y` 中额外记录下列四个只读求积质量诊断，
因保留旧版没有对应字段，算子投影仅剔除它们：

- `quadratureWeightsStrictlyPositive`
- `minimumQuadratureWeightToControlWidthRatio`
- `maximumQuadratureWeightToControlWidthRatio`
- `maximumQuadratureWeightToControlWidthAbsoluteDeviationFromUnity`

其余 `gridQuality` 字段、实际 quadrature weights、所有差分/输运/椭圆
算子、迁移场与指标、以及 Poisson 分解后的求解作用仍精确比较。

schema 4 的动态规范化已改为无 limiter、无 restoring feedback 的精确规范合同，
因此不再与保留旧版的动态轨迹做 bitwise 等价；该合同由基线重标度专项测试。
无参数算例的未变选项与六阶辅助配置仍做精确比较。所有数值比较不使用误差容限，
测试结束恢复原路径和工作目录。

来源提交和原版全部 140 个已跟踪文件的 SHA-256 记录在 [source_map.json](source_map.json)，
其中 `functions` 是旧入口到新包名的对应表；`files` 记录迁移位置，`destination=null`
表示未原样复制。来源哈希用于追溯，不等同于新文件应保留相同字节。

```text
tests/
├── +ipmtests/
│   ├── +baseline/       原基线分层门槛
│   ├── +fourth/         四阶验证
│   ├── +sixth/          六阶验证
│   ├── +support/        制造解等测试辅助函数
│   └── equivalence.m    新旧数值一致性
└── source_map.json      来源与迁移清单
```

各套件通过断言报告失败；检查实际报告和异常，不把脚本启动或单个算例结束视为整套通过。
历史 `ipmtests.release` 仅组合基线与四阶，不等于包含六阶的 `ipm.verify('all')`。
测试可能使用临时文件验证输出行为；不要将测试目录加入生产运行依赖。
当前执行结果和环境见 [VALIDATION.md](VALIDATION.md)，变更见 [CHANGELOG.md](../CHANGELOG.md)；旧证据的适用范围见
[研究索引](../research/README.md)，尤其不可将局部算子测试扩大解释为任意问题的整体收敛保证。
