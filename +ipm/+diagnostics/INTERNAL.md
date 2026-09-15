# diagnostics：特征、质量与停止条件

由当前状态或记录序列计算诊断量，不修改密度、算子或配置，也不执行文件写入。
导航：[包索引](../README.md) · [结构总览](../../STRUCTURE.md) ·
[第一象限架构记录](../../ARCHITECTURE_QUADRANT_LEVELSET_20260914.md)。

## 接口

- `measure(rho,flow,ops,mass0)`：当前场的质量、梯度、值域、边界及分辨率指标。
- `trackFeatures`：供精确规范定位、诊断、重网格与停机使用的当前峰值/宽度特征。
- `measureFeatures`、`peakLocation`、`peakResolution`：一维特征及有效点数测量。
- `resolutionFailed(flow,ops)`：分辨率不足的离线判定；只作为质量证据，不终止步进。
- `stopPolicy(...)`：保留历史阈值分类供离线审计；求解主循环不再调用它。
- `finalQuality(history)`：最终质量摘要。
- `blowupFit(t,gradInf)`：输入时间窗内的梯度拟合，不自行声明奇性。
- `maximumGrowthRateFit(t,M)`：物理时间下比较指数、有限时幂律与
  双指数，使用 e-fold 尾窗、留出验证及采样稳定性门。
- `canonicalMaximumGrowthFit(tau,M)`：递增 canonical `tau_c` 下在同一
  `log(M)` 误差域比较线性最大值与指数最大值，并以二次
  `log(M)` 模型只作非渐近曲率守门。
- `gaugeRateConvergence(tau,cL,cOmega,options)`：在 canonical `tau_c` 末
  25%/40%/60% 三个固定嵌套窗口上，以 canonical-time 梯形权重
  记录 `c_l`、`c_omega`、`c_omega/c_l` 及
  `(c_l-c_omega)/c_l` 的加权均值、斜率、相对趋势和去趋势 RMS；
  另外预注册无量纲零极限候选 `z=c_omega/c_l -> 0`。

## 约定

特征定义要在规范、迁移与记录之间一致；不能仅为美化曲线修改阈值或掩去失败点。
当前流场及网格为输入事实；诊断函数不重新求解速度或默认配置。
没有可用正壁面特征的物理算例不应被当作动态规范初始化失败。

可信历史窗口由 `output.trustedMask` 选择，再交给时间序列拟合；
拟合与质量指标不替代空间加密、远场边界检查或爆破证明。
`output.finalize` 将连续可信前缀上的 `physicalRhoXInf` 同时交给物理时间
`maximumGrowthRateFit` 和 canonical 时间 `canonicalMaximumGrowthFit`，并保留旧
`blowupFit(physicalTime,physicalGradInf)` 结果以维持外部契约。
同一连续可信前缀还会交给 `gaugeRateConvergence`；该诊断至少
需要 16 个点和严格递增的 `tau_c`。`provisionalStable` 固定要求
`c_l` 与 `c_omega` 在所有尾窗内的绝对相对趋势均不超过 1%，
三窗加权均值的相对极差均不超过 1%，且所有尾窗的
去趋势相对 RMS 均不超过 1%。比值另有独立的
`ratioDiagnosticsStable`，不替代两个规范速率的判据。所有输出明示标记
`diagnostic_only_requires_box_time_grid_refinement`；通过该门也不是渐近常数的证明。
若任一核心速率的加权均值太接近零，使相对趋势或相对 RMS
无定义，诊断记为 `relative_statistics_undefined`，不得进入候选排名。
零极限分支不改变上述 `provisionalStable` 的版本-2 数值语义。它要求
`c_l` 非零且通过旧相对稳定门，并要求每个嵌套尾窗的
`max(abs(z))`、`abs(slope(z))*window_span`、去趋势 RMS，以及三窗
加权均值的绝对极差，全部不超过同一个 1% 阈值。这只产生数值字段
`cOmegaZeroLimitStable`。新的 `zeroAwareProvisionalStable` 还要求调用者通过
`options.selectedRatePrescribed=false` 明确声明所选速率不是预设常数。省略该选项时
provenance 不可用；传入 `true` 时则速率为 prescribed，两种情况的
`empiricalRateConvergenceEligible` 都为 false。因此精确零不会自动成为“经验收敛”。
当前判据记为 `criteriaVersion=3`；旧结果中没有该版本的
`provisionalStable` 不能直接当作新判据，必须从历史序列重算。
物理 `t` 与 canonical `tau_c` 的拟合字段不可互换解释；
`tau_c` 也不是剩余时间 `T-t`。
若且仅若重标度场的最大值也趋于常数，而
`c_l -> c_l*`、`c_omega -> c_omega*` 且
`kappa=c_l*-c_omega*>0`，则尺度恒等式预测
`max|rho_x| ~ exp(kappa*tau_c)` 且
`1/max|rho_x| ~ const*(T-t)`。因而“线性增长”不应在
canonical 时间上预设给原始最大值；必须分别比较 canonical
线性/指数模型和物理时间上的倒数线性/有限时幂律。
停止原因属于结果契约，即使提前停止也要保留最后可用状态及对应诊断。
在 schema-4 `exact_gauge_no_feedback_v1` 中，特征的位置/场值可以进入规范的精确
代数约束；但峰值、层级、面积或连通宽度的网格单元计数只能被记录、触发
remesh 或硬停机，不得作为 `c_l/c_omega/c_r` 的反馈信号。

本模块只调用 [mesh](../+mesh/INTERNAL.md) 的象限判定/奇偶导数选择，
向 [remesh](../+remesh/INTERNAL.md) 提供只读特征，不调用求解或重网格事务，
便于测试和研究复用。
相关验证：[数值核心](../../tests/+ipmtests/+baseline/numericCore.m)、[缩放](../../tests/+ipmtests/+baseline/scaling.m)、
[输出](../../tests/+ipmtests/+baseline/output.m)。

meshFeatureIntervals读取真实rho/轴/Dx及逐值相等的rho*Dx，输出正壁面核心、前沿和所选节点列的纵向宽度。
第一象限模式直接读取非负轴，不构造镜像观测；远边界诊断只包含 `X=H` 和 `Y=Ymax`，
`X=0` 是对称轴而非人工远边界。
其actualCoreCells保留原生90%计数；该纵向定义仍来自nodal peak column，不冒充连续峰位置的宽度。
这些观测仅用于网格规划与验收，不能反馈规范速率或修改演化场。
