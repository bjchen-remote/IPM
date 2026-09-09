# 当前 active datum 的四阶 256 x 256 / 384 x 384 爆破测试

日期：2026-09-03

## 结论摘要

本轮针对当前维护的 `degenerate_primitive`, `k=8` 初值

\[
\rho_0(x,y)=\frac18\log\left(1+\frac{|x|^8}{1+y^8}\right)
\]

运行了 `high_order / weno5_fd / ssprk54 / high_order` 数值组合。由于
`double_odd` 对称实现存储整个水平区间，513 x 257 个节点对应第一象限
256 x 256 个单元；这是报告中“256 x 256 网格”的含义。

最有效的数值改进不是放宽停机阈值，而是用双峰集中监视函数取代分段
`lattice` 重映射。513 x 257 lattice 网格只可信到物理时间
`0.0219116`；相同节点数下，全局解析双高斯 `A=16,w=0.55` 网格的
连续局部 trusted prefix 到达 `t=1.659748`，可信窗口延长
75.747438 倍。该算例的 raw stop 为 `t=1.797499`，原因仍是
`grid_resolution_failure`，故 raw 终点不能当作证据。大小盒的局部场
对照到 `t=0.5` 低于 0.34%；时间步对照只完成到 `t=0.12`，
其误差为 `1e-9` 量级。

在相同解析双高斯 `A=16,w=0.55` 上进一步加密到第一象限
384 x 384 单元（存储 769 x 385 节点）后，连续可信前缀到
`t=1.75158228257`，相对 256 x 256 的 `t=1.65974756404` 延长
5.533053%。计算请求到 `t=2.05`，但 raw 轨迹在
`t=1.85594004481` 就触发 `grid_resolution_failure`；这个 raw 终点同样
不是证据终点，也不能表述为已计算到 `t=2.05`。

静态正则性复核随后发现，此前归档的平滑集中网格都使用了
`legacy_abs_gaussian` 监视函数；其 `abs(x)` 因子在 `x=0` 非 C1，
不满足四阶映射应有的全局光滑性。这些 legacy 结果因而只作为
历史对照。新的 `analytic_double_gaussian` 已成为 helper 默认，
其 `A=16,w=0.7` 独立复跑与 legacy 标量轨迹高度一致。进一步把
解析 monitor 宽度改为 `w=0.55` 后，连续可信时间比 analytic
`w=0.7` 提高 0.927975%，比 legacy `w=0.7` 提高 0.814928%。

旧的点数窗口 `blowupFit` 门槛在所有已完成 legacy 与 analytic 网格上均为
`candidate=false`。针对最大值增长率的新复核改用
`gamma=d(log M)/dt`、物理时间加权、固定等距评分点和按 e-fold 定义的窗口。
在 256 x 256 基线上，`A=16,w=0.55` 的最大梯度和壁峰在三个较长窗口内均更偏向
有限时幂律；最短窗口会随输入重采样在非线性模型与不可区分之间翻转，
已按保守门槛降级。梯度和壁峰总共只覆盖
`1.686/1.738` 个 e-fold。因而现在可以报告“稳定的加速增长和有限时幂律
候选外推”，不能报告可信地区分了幂律与双指数，更不能报告找到了
可信爆破解。
新的 384 x 384 轨迹将梯度/壁峰动态范围提高到
`1.998490/2.042679` 个 e-fold，但仍未达到 3 e-fold 可信门槛；
它的最终模型序列也是最短窗口不可区分、后三个窗口偏有限时幂律。
按 canonical 重标度时间 `canonicalTau` 重新分析后，物理最大梯度和
壁峰的最后 `0.25--0.5` e-fold 局部窗支持 `tau_c`-指数增长，
明显优于 `tau_c`-线性增长；较长窗仍检出速率曲率和漂移。
这与物理时间下的近 Type-I 幂律并不矛盾，
而是双尺度动态重标度轨道的两种时钟表示。

## 关于 MATLAB 智能分解

MATLAB 对应接口是 `decomposition`。当前高阶 Poisson 路径已经使用
`decomposition(A,'lu')`，并复用该分解求解每个时间步。对实际 33 x 33
高阶矩阵测试 `decomposition(A)` 时，MATLAB 自动选择的类型仍是 `lu`；
该矩阵的相对非对称量约为 `0.941`，不能使用 Cholesky。故把显式 LU
改写成自动模式不会改变算法：对本问题，
`decomposition(A)` 与现有 `decomposition(A,'lu')` 等价，且现有分解已在
各时间步复用。分解器不是当前瓶颈，也不能修复早期失效。
129 x 129 初始失败
发生在分解之前：过度伸缩网格令四阶求积权变为负值。

## 数值范围与可信判据

本轮没有把内部 `blowupFit.candidate` 当成爆破证明。可信结论额外要求：

- 终点位于 trusted mask 内，不能在壁面振荡或分辨率停机之后外推；
- 密度范围、壁面轮廓、峰值符号、网格安全性和相邻单元比通过检查；
- 自定义映射和监视函数具有四阶离散所需的全局正则性与正求积权；
- 梯度、壁峰和宽度诊断给出一致、跨窗口稳定的奇性时间与指数；
- 时间步、空间网格和计算盒变化下结果收敛；
- 明确报告远场和开放边界警告。

内部候选门槛仅要求幂律拟合 `R^2 > 0.98` 且比指数拟合高 `0.01`；它本身
不足以支持“可信爆破”。

## 最大值对数增长率复核

新的主诊断令

\[
M(t)=\|\nabla\rho(t)\|_\infty,\qquad
\gamma(t)=\frac{d}{dt}\log M(t).
\]

简单指数、有限时幂律和双指数分别对应

\[
\gamma=\lambda,\qquad
\gamma=\frac{p}{T-t},\qquad
\gamma=\gamma_0 e^{\kappa(t-t_0)}.
\]

实现使用局部二次拟合估计 `gamma` 和 `gammaDot`，三种带宽用于敏感性
检查；模型排序仍直接拟合 `log(M)`，避免把微分噪声当作独立观测。
所有回归和误差均采用物理时间梯形权重，每个增长窗口最多用 64 个等距
物理时间评分点，blocked validation 按前 70% 物理时间训练、后 30% 验证。
非线性 profile 命中搜索边界、退化到简单指数，或 AICc/验证/增长率带宽
不能一致时，结论为 `indistinguishable`。此外把完整输入先重采样到最多
512 个等距物理时刻再重复判别；若结论改变，同样降级。窗口由最后
`0.25/0.5/1/1.5` 个 e-fold 定义，不再由最后若干 solver 记录定义。
等距抽样降低了输出密度偏置，但不会把同一条数值轨迹变成 64 个独立
实验；AICc 差只用于区间内模型排序，不解释为统计置信度或爆破概率。

256 x 256 analytic `A=16,w=0.55` 基线的连续可信前缀结果为：

| 梯度增长窗口（e-fold） | 胜出模型 | `T` | `p` | 简单指数 `lambda` |
|---:|---|---:|---:|---:|
| 0.25 | indistinguishable（原序列偏 power） | 2.002703 | 1.033165 | 2.658160 |
| 0.50 | finite-time power | 2.035373 | 1.120764 | 2.352610 |
| 1.00 | finite-time power | 2.045224 | 1.143174 | 1.810676 |
| 1.50 | finite-time power | 2.063708 | 1.178507 | 1.374928 |

`gamma` 在 `t=0.5/0.8/1.0/1.2/1.4/1.6` 约为
`0.7445/0.9420/1.1107/1.3579/1.7703/2.5653`，到可信终点三种带宽给出
`2.96936/2.96883/3.01407`；相应 `gammaDot` 为
`5.3579/5.3128/7.8730`。这排除了“当前区间内常 `gamma` 平台”作为最佳
描述，但二阶量仍有约 48% 的带宽跨度，不能据此严格排除更一般的非定常
指数增长。

壁峰总增长为 `1.738296` 个 e-fold，终点 `gamma` 为
`2.88196/2.88051/2.87833`。其原始 `0.25` e-fold 窗口偏双指数，但在
完整输入等距重采样后降为不可区分；后面三个窗口稳定选择有限时幂律，
因此整体仍是 `window_dependent`。动态 `wall_omega_peak` 规范还给出
不依赖数值微分的物理缩放率贡献

\[
\gamma_{\rm scale}
=(\mathrm{canonicalCL}-\mathrm{canonicalCOmega})\,C_l/C_\omega .
\]

其可信终点为 `2.883854`。它不包含受控重标度壁峰自身的微小变化；与
完整物理壁峰三带宽数值导数的整段相对 RMS 差为 `0.78%--1.02%`，
确认了壁峰加速率本身不是差分假象。

宽一些的 analytic `A=16,w=0.70` 网格给出相同的最终模式：最短梯度
窗口不可区分，三个较长窗口选择有限时幂律；原序列拟合参数为
`T=1.99978/2.03507/2.04724/2.06671` 和
`p=1.02945/1.12041/1.14706/1.18304`。这支持拟合对这两种解析节点分配
的内部稳定性，但两者节点数相同，不是空间加密检验。

对这个 256 x 256 基线，最强而严谨的解释是：有限时幂律在现有数据中受到更多支持，简单指数在
可信区间内被一致排在后面；但候选 `T` 比可信终点晚约
`0.34--0.40`，总动态范围不到两个 e-fold，两个最大值的短窗口又都对
输入采样敏感。
在达到至少 3 个 e-fold，并完成晚期时间步、计算盒和多层更高节点数的重合
验证前，power 与 double exponential 均不能升级为可信渐近分类。

## 384 x 384 空间加密长跑与增长曲线

本节数值直接取自归档的 q384
[`result`](../../../result/verification/fourth_order_active_analytic_A16_w055_q384_t205_dt002.mat)
和 [`assessment`](../../../result/verification/fourth_order_active_analytic_A16_w055_q384_t205_dt002_assessment.mat)。
该算例存储 769 x 385 节点，即 `double_odd_omega` 对称下第一象限
384 x 384 单元。计算盒仍为 `[-281,281] x [0,132]`，解析双高斯
monitor 仍为 `A=16,w=0.55`，但中心目标间距按网格加密缩小为
`[0.0193174285333333,0.000833333333333333]`，`fineCount=196609`。
初始 x/y 核心点数为 `108.093590/172.672499`，最大相邻单元比为
`1.191629/1.023452`，因而这是对物理中心间距的真正加密，而不是仅增加外区节点。

数值组合为
`high_order / weno5_fd / ssprk54 / high_order`，`CFL=0.25`，
`maxDt=0.002`，请求物理终点为 `t=2.05`。计算在 2982 步后于
raw `t=1.85594004481` 触发 `grid_resolution_failure`，所以它并未到达
`t=2.05`。更重要的是，1007 个记录中只有前 802 个构成连续可信前缀，
证据终点为 `t=1.75158228257`。**raw 停机时间不是证据终点**，
不允许用 `1.85594004481` 上的 under-resolved 数据推动拟合结论。
与 q256 相同 monitor 的可信终点 `1.65974756404` 相比，q384 达到
`1.05533053370` 倍，即延长 `5.533053%`。可信终点的 safety 为
`0.998638640`，x/y 核心点为 `8.010906/14.598992`。

在这个连续可信前缀上，最大梯度增长 `7.377904` 倍，即
`1.998489603` 个 e-fold；壁峰增长 `7.711241` 倍，即
`2.042679102` 个 e-fold。梯度终点的三带宽 `gamma` 为
`3.893169/3.890059/3.920758`，壁峰为
`3.772760/3.768482/3.757316`。四个增长窗口的最终判决如下；
`T,p` 只在最终选中有限时幂律时报告。

| 窗口（e-fold） | 梯度最终判决 | 梯度 `T / p` | 壁峰最终判决 | 壁峰 `T / p` |
|---:|---|---:|---|---:|
| 0.25 | indistinguishable | -- | indistinguishable | -- |
| 0.50 | finite-time power | `2.034970 / 1.117645` | finite-time power | `2.036468 / 1.093508` |
| 1.00 | finite-time power | `2.038213 / 1.127549` | finite-time power | `2.033039 / 1.083758` |
| 1.50 | finite-time power | `2.045426 / 1.145480` | finite-time power | `2.022373 / 1.058170` |

梯度的 0.25 e-fold 原序列 AICc 偏双指数，但等距重采样后不可区分，
因而被降级；壁峰的最短窗口在采样敏感性门槛前已不可区分。
两个序列的总体共识均是 `window_dependent`。完整六面板图见
[`q384 growth curve`](../../../result/verification/fourth_order_active_analytic_A16_w055_q384_t205_dt002_growth.png)，
256/384 公共时间比较见
[`resolution comparison`](../../../result/verification/fourth_order_active_A16_w055_q256_vs_q384_growth.png)；
两张图都只使用各自的 trusted prefix。

q384 的 assessment 明确给出 `rateEvidenceSufficient=false`、
`canonicalRateEvidenceSufficient=false`、
`singleRunCandidate=false` 和 `credibleCandidate=false`。原因是梯度仍不足
3 e-fold，最短窗口不可区分，拟合的 `T=2.02--2.05` 全部位于
`t=1.751582` 可信证据区间之外，而且计算先于请求终点触发了分辨率门。
此外，可信终点远场速度比仍为 `0.930036`，物理质量漂移为
`-0.946684`。因此这次加密支持“加速增长在更高节点数下继续，且较长窗口偏有限时幂律”，
但**仍没有找到可信爆破解**。一次 256 -> 384 加密也不足以建立空间收敛阶。

## Canonical `tau_c` 下的双尺度增长

本节中的 `tau_c` 专指程序保存的递增重标度时钟
`common.canonicalTau`，不是文献中常写成 `T-t` 的剩余时间。
本次 q256/q384 可信前缀的 `timeSpeed` 都为 1，所以它正是实际
canonical 规范时钟，而不是限速后的替代坐标。程序缩放约定给出

\[
q(\tau_c)=\frac{C_l}{C_\omega},\qquad
\frac{dt}{d\tau_c}=\frac1q,\qquad
M_{\rm phys}=q\,M_R.
\]

对壁面峰值，`wall_omega_peak` 规范将重标度峰
`M_R` 锁定在近乎常数；q384 上它仅从 `0.686068588`
变到 `0.686067018`。因而

\[
M_{\rm wall}=q M_{R,\rm wall},\qquad
\partial_{\tau_c}\log M_{\rm wall}
=c_l-c_\omega+\partial_{\tau_c}\log M_{R,\rm wall}
\simeq c_l-c_\omega .
\]

q384 连续可信段为 `tau_c=0--4.006902434`。在共同 `log M`
误差域中比较线性最大值 `M=A+B(tau_c-tau_0)` 与指数最大值
`M=A exp(kappa(tau_c-tau_0))`，得到：

| q384 尾窗 | 梯度 `kappa` | 壁峰 `kappa` | 梯度指数/线性验证 NRMSE | 壁峰指数/线性验证 NRMSE |
|---:|---:|---:|---:|---:|
| `Delta tau_c=0.5` | 0.51635 | 0.50317 | 0.049% / 0.606% | 0.148% / 0.484% |
| `Delta tau_c=1.0` | 0.52174 | 0.50831 | 0.211% / 2.218% | 0.259% / 2.060% |
| `Delta tau_c=1.5` | 0.52760 | 0.51380 | -- | -- |
| `Delta tau_c=2.0` | 0.53379 | 0.51896 | -- | -- |

专用 canonical-time 拟合器在同一 `log M` 误差域中额外加入
`quadratic-log` 曲率守门，避免把任意缓慢变率曲线强制二分。q384 按
`0.25/0.5/1/1.5` e-fold 的最终判决为：

| e-fold 窗 | 梯度 | 壁峰 |
|---:|---|---|
| 0.25 | exponential maximum | exponential maximum |
| 0.50 | exponential maximum | exponential maximum |
| 1.00 | curved/nonasymptotic | indistinguishable |
| 1.50 | curved/nonasymptotic | curved/nonasymptotic |

所以数据排斥“当前尾段只是 `tau_c` 线性最大值”，但还不能把
整段写成已达渐近常指数率。在 q256/q384 共同终点
`tau_c=3.400853591` 上，后 `Delta tau_c=1` 的梯度斜率为
`0.536353/0.536174`，壁峰为 `0.521907/0.521742`，显示共同区间
的空间一致性。但在 q384 更长的末段，斜率仍由长窗的
`0.52--0.54` 向短窗的 `0.50--0.52` 漂移，所以尚不能声称已经
达到严格常数平台。

尺度分解进一步给出，q384 最后 `Delta tau_c=0.5` 上

\[
C_l\sim e^{0.2302\tau_c},\qquad
C_\omega\sim e^{-0.2729\tau_c},\qquad
q\sim e^{0.5032\tau_c}.
\]

终点即时组合率 `c_l-c_omega=0.490022`。如果这个正极限
`kappa` 持续，则

\[
T-t=\int_{\tau_c}^{\infty}q(s)^{-1}\,ds
\sim\frac{1}{\kappa q(\tau_c)},\qquad
M_{\rm wall}\sim\frac{M_{R,\rm wall}}{\kappa(T-t)}.
\]

由末段和终点率估计 `T=2.007--2.016`，与物理时间拟合的
`T=2.02--2.04` 同量级。因而此处更正后的主解释是：
**物理 `t` 下的近 Type-I 有限时幂律，与 canonical `tau_c` 下的
指数轨道，是同一双尺度机制的等价描述。**

如果所说的“`tau_c` 下线性/指数增长”是指重标度最大值
`M_R`，则还有更精细的区分：

\[
M_R\sim \tau_c^m
\Rightarrow M_{\rm phys}\sim
\frac{[\log(1/(T-t))]^m}{T-t},
\]

\[
M_R\sim e^{\nu\tau_c}
\Rightarrow M_{\rm phys}\sim
(T-t)^{-(1+\nu/\kappa)}.
\]

当前 q384 重标度全局梯度在最后 `Delta tau_c=1` 只有
`nu=0.01341` 的弱增长，而重标度壁峰被规范固定。与此同时，
重标度壁面 FWHM/核宽在同一窗中约按
`exp(-0.66481 tau_c)` / `exp(-0.65796 tau_c)` 收缩，表明固定
`X` 坐标中还有未冻结的内尺度；这也正是最终触发分辨率失效的结构。

但壁峰的指数律主要是固定峰值规范与尺度比的恒等式，
不能重复计为独立证据。当前算例的 `dynamicScaleGeometry`
仍为 `isotropic`，`C_x=C_y` 也不能自动证明两个空间尺度。真正独立的
双尺度判据应来自重标度全局梯度、水平/竖直特征宽度、形状收敛及
理论预测率比的同时收敛。因此当前最准确的结论是
“**强烈且跨网格一致的 `tau_c`-指数候选，与双尺度机制相容；
独立双尺度率比尚未收敛**”，而不是已经证明爆破。

`tau_c` 下的两网格增长、尺度分解、规范率与尾部 `T` 估计见
[`canonical-tau growth figure`](../../../result/verification/fourth_order_active_A16_w055_q256_vs_q384_tau_growth.png)。

## 网格优化

下表中未显式标注 analytic 的双峰网格均是
`legacy_abs_gaussian` 归档数值；最后三行是全局解析双高斯复跑与加密。
“连续可信”表示通过 solver trusted mask，但仍需结合映射正则性、
远场和收敛对照审计。

| 存储节点 | 网格 | 初始 x 核心点 | 最大相邻比 x | `max |Delta^2 log(dx)|` | 连续可信至 |
|---|---|---:|---:|---:|---:|
| 257 x 257 | 原始 sinh | 3.24 | 约 1.12 | 小 | 仅时间步探针 |
| 257 x 257 | lattice，两次 | 9.25 | 1.50 | 较硬 | `0.0751061` |
| 513 x 257 | lattice，两次 | 27.01 | 1.50 | 0.352 | `0.0219116` |
| 513 x 257 | 大盒平滑，A=2,w=0.85 | 16.1061 | 1.153270 | 0.0188 | 无早期 lattice 振荡 |
| 513 x 257 | 小盒平滑，A=0.5,w=0.7 | 15.5816 | 1.041273 | 0.002540 | `0.5` |
| 513 x 257 | 小盒聚焦，A=1,w=0.7 | 19.5553 | 1.054873 | 0.005760 | `0.678278` |
| 513 x 257 | x/y 均衡聚焦，A=2,w=0.7 | 26.3628 | 1.075481 | 0.012430 | `0.960992` |
| 513 x 257 | 集中，A=3,w=0.7 | 32.0122 | 1.093088 | 0.019292 | `1.334758` |
| 513 x 257 | 近终端，A=12,w=0.7 | 58.2049 | 1.239133 | 0.078061 | `1.621842` |
| 513 x 257 | 稳定极限，A=16,w=0.7 | 63.6904 | 1.302353 | 0.099551 | `1.646331` |
| 513 x 257 | analytic 稳定基线，A=16,w=0.7 | 63.4079 | 1.309616 | 0.076436 | `1.644487` |
| 513 x 257 | analytic 当前最优，A=16,w=0.55 | 72.0562 | 1.301284 | 0.074148 | `1.659748` |
| 769 x 385 | analytic q384 加密，A=16,w=0.55 | 108.0936 | 1.191629 | 0.036693 | `1.751582` |

小盒为 `[-281,281] x [0,132]`。上表初始核心点采用实际 solver
诊断（近终端和稳定极限网格由同一解析初值与 factory 网格复核）。
所有接受网格均保持正四阶求积权，也没有改变任何 trusted-mask
或硬停机阈值。

更激进的 `A=20/24/28` 分支在壁面轮廓审计中因振荡被拒绝；
`A=20,w=0.7` 的 `t=0.2` 短筛通过不等于长时接受。
对 `A=16` 缩窄监视宽度也不能单调改善结果：`w=0.35` 在
raw `t=0.275560` 因壁面振荡停机（可信至 `0.269956`），
`w=0.45` 在 raw `t=0.453370` 停机（可信至 `0.449072`）。
legacy `w=0.55` 只完整确认到 `t=0.5`；在监视函数正则性问题
确认后，其冗余 legacy 长跑已中止，不用于刷新最佳可信时间。

## 解析双高斯复核

`fourth_order_smooth_peak_grid` 现默认使用全局解析的
`analytic_double_gaussian`；旧命名算例则显式锁定
`legacy_abs_gaussian` 以保留结果复现性。新的
`analytic_stable_limit_box_quadrant_256` 使用 `A=16,w=0.7`，是对最优
legacy 网格的必要重跑。

解析网格的 profile 元数据为 version 2，monitor 积分使用
`fineCount=131073`。`A=16,w=0.7` 初始 x/y 核心点为
`63.40794/115.08640`，最大相邻单元比为
`1.30961649/1.03538355`，`log(dx)`/`log(dy)` 最大二阶差分为
`0.07643618/0.00120556`。7 点 stencil 最小 rcond 为
`8.804e-6/9.565e-6`，归一化最小求积权为
`0.0063475/0.00077649`；因此初始映射通过正求积权、相邻比与
stencil 条件数门槛。

更强的解析幅度在壁面振荡门下被拒绝：

| analytic 筛选 | raw stop | 连续可信至 | 停机原因 | 四类候选 |
|---|---:|---:|---|---|
| A=24,w=0.7 | 0.132121853 | 0.128389363 | `wall_profile_oscillation` | 0/0/0/0 |
| A=20,w=0.7 | 0.279627095 | 0.274291876 | `wall_profile_oscillation` | 0/0/0/0 |
| A=18,w=0.7 | 0.369873884 | 0.366480746 | `wall_profile_oscillation` | 0/0/0/0 |

因此保留 `A=16,w=0.7` 作为解析稳定基线。它在 raw
`t=1.79154446692` 因 `grid_resolution_failure` 停机，869 个记录中前 664 个
构成连续局部 trusted prefix，可信至 `t=1.64448713611`。该端点的
梯度/壁峰增长为 `5.16185266/5.44512851`，x/y 核心点为
`8.02825/14.69442`，四类候选全部为 false。其可信时间是 513 x 257
lattice 结果的 75.051 倍，比 legacy `A=16` 仅低 0.112%。

旧的最后 40 点 `R^2` 门槛仍为 false；该点数窗口会随晚期输出密度改变，
且旧版 profile 搜索曾受有限上界影响，现不再作为主判据。按上节的
物理时间/e-fold 增长率诊断，`w=0.70` 的梯度三个较长窗口偏有限时
幂律，最短窗口因输入采样敏感而降为不可区分；总动态范围只有
`1.641296` 个 e-fold。这仍是增长外推，不是爆破证明。

在 `t=0.96/1.3/1.644` 的共同区间内，analytic/legacy 主要标量差约为
`1e-4--6e-4`，说明去掉原点非 C1 点后主轨迹未发生实质改变。
`t=0.5` 的 FWHM 曾单独显示约 9.4% 差异，但同时核心宽度差只有
`8.8e-5`。这是 FWHM 的离散交点分支/采样别名跳变，而不是主轨迹
的 9.4% 偏离；因而单一 FWHM 数值不能用来判定物理不一致。

解析 `A=16,w=0.55` 的 256 x 256 基线算例现已完成，并成为该分辨率下最优的局部数值质量
前缀。它与上述 `w=0.7` 使用同一 profile version 2 和
`fineCount=131073`；初始 x/y 核心点为 `72.056164/115.151066`，最大相邻比
为 `1.301284189/1.035383549`，曲率为 `0.074148338/0.001205563`。7 点
stencil 最小 rcond 为 `9.0403e-6/9.5650e-6`，归一化最小求积权为
`0.0054952/0.00077649`，未归一化最小求积权为
`0.006020124/0.000398820`。

该算例在 raw `t=1.79749926110` 因 `grid_resolution_failure` 停机；880 个记录中
前 681 个构成连续局部 trusted prefix，到 `t=1.65974756404`。相对
513 x 257 lattice，这是 75.747438 倍；它分别比 analytic/legacy
`A=16,w=0.7` 延长 0.927975%/0.814928%。可信端点 safety 为
`0.99808118`，x/y 核心点为 `8.01538/13.76353`，梯度/壁峰增长为
`5.397535/5.687645`，四类候选仍全部为 false。

旧的点数窗口 `R^2` 候选门槛仍为 false，但不再承担增长模型分类。
新的增长率诊断在梯度四个 e-fold 窗口给出
`T=2.00270--2.06371,p=1.03316--1.17851`；采样敏感性门槛后的判决为
`indistinguishable/power/power/power`，壁峰也是同一最终序列（其原始
最短窗口偏双指数）。所以当前是有限时幂律候选，而不是已经由多个
最大值诊断一致确认的渐近律。

`w=0.55` 相对 `w=0.7` 在 `t=0.96/1.3/1.644` 的梯度、壁峰和核心宽度
差异多为 `1e-4--8e-4`。个别早期采样中梯度相对差可达 1.1%--1.3%，
FWHM 为 6%--7%；但在匹配物理时间 `t=0.5` 时，梯度/核心宽度差仅为
`4e-6/3e-6`。这与前述 FWHM 离散交点分支一致，不应把孤立宽度
采样跳变解释为物理轨迹分叉。

可信端点的物理/重标度质量漂移为 `-0.929168/-0.346311`，远场
速度/源项比为 `0.923439/0.0144411`，密度范围违反为 `7.434e-8`，
散度和 Poisson 残差为 `4.625e-13/1.569e-11`。它是 256 x 256 基线的最优局部
数值质量前缀，但这些晚期警告和当时缺失的多网格收敛检验仍阻止将它
升级为可信爆破证据。

## 已完成长时结果

下表的核心点数和增长率均取自最后连续可信记录，不是后面的
under-resolved raw 终点。候选顺序为梯度/壁峰/逆 FWHM/逆核心宽度。
标为 smooth/focused/balanced/concentrated/near-terminal/stable-limit 的旧行
均使用 legacy monitor；显式标为 `analytic stable/optimal/q384` 的三行
使用解析双高斯。

| 算例 | 请求终点 | raw stop | 连续可信至 | 停机原因 | 可信端点 x/y 核心点 | 梯度/壁峰增长 | 候选 |
|---|---:|---:|---:|---|---:|---:|---|
| 257 square lattice | 0.10 | 0.077056 | 0.075106 | `wall_profile_oscillation` | 9.557 / 44.474 | 1.02750 / 1.02533 | 1/0/1/0，不可信 |
| 513 x 257 lattice | 0.12 | 0.023895 | 0.021912 | `wall_profile_oscillation` | 27.608 / 43.775 | 1.02420 / 1.00809 | 1/0/0/0，不可信 |
| 513 x 257 smooth large | 0.03 | 0.03 | 0.03 | `physical_final_time` | 15.983 / 43.617 | 1.00065 / 1.01146 | 0/0/0/0 |
| 513 x 257 smooth small | 0.12 | 0.12 | 0.12 | `physical_final_time` | 14.987 / 41.104 | 1.00454 / 1.04757 | 0/0/0/0 |
| 513 x 257 smooth small | 0.50 | 0.50 | 0.50 | `physical_final_time` | 12.066 / 11.987 | 1.08862 / 1.27241 | 0/0/0/0 |
| 513 x 257 smooth large | 0.50 | 0.50 | 0.50 | `physical_final_time` | 13.057 / 20.446 | 1.09074 / 1.27429 | 0/0/0/0 |
| 513 x 257 focused, A=1 | 1.00 | 0.968606 | 0.678278 | `grid_resolution_failure` | 13.128 / 8.006 | 1.25606 / 1.43544 | 0/0/0/0 |
| 513 x 257 balanced, A=2 | 1.20 | 1.20 | 0.960992 | `physical_final_time`，终态非 trusted | 13.235 / 9.294 | 1.64743 / 1.82401 | 0/0/0/0 |
| 513 x 257 concentrated, A=3 | 1.50 | 1.50 | 1.334758 | `physical_final_time`，终态非 trusted | 9.274 / 9.383 | 2.67910 / 2.88000 | 0/0/0/0 |
| 513 x 257 near-terminal, A=12 | 1.85 | 1.779629 | 1.621842 | `grid_resolution_failure` | 8.029 / 13.802 | 4.839197 / 5.113862 | 0/0/0/0 |
| 513 x 257 stable-limit, A=16 | 1.95 | 1.790679 | 1.646331 | `grid_resolution_failure` | 8.014 / 14.210 | 5.189615 / 5.474471 | 0/0/0/0 |
| 513 x 257 analytic stable, A=16 | 1.95 | 1.791544 | 1.644487 | `grid_resolution_failure` | 8.028 / 14.694 | 5.161853 / 5.445129 | 0/0/0/0 |
| 513 x 257 analytic optimal, A=16,w=0.55 | 1.95 | 1.797499 | 1.659748 | `grid_resolution_failure` | 8.015 / 13.764 | 5.397535 / 5.687645 | 0/0/0/0 |
| 769 x 385 analytic q384, A=16,w=0.55 | 2.05 | 1.855940 | 1.751582 | `grid_resolution_failure` | 8.011 / 14.599 | 7.377904 / 7.711241 | 0/0/0/0 |

`t=0.5` 小盒终点 safety 为 `0.667372`，大盒为 `0.612715`；密度范围
违反为零，散度约 `1.5e-13`。小盒远场源项比为 `0.01018`，大盒降到
`7.79e-5`。远场速度比仍分别约为 `0.700` 和 `0.919`，反映此慢衰减
datum 的长程速度；这项警告没有被隐藏。

更晚的可信端点警告明显更强。`A=12` 的物理/重标度质量漂移为
`-0.921048/-0.311130`，远场速度/源项比为 `0.916196/0.014634`；
legacy `A=16` 对应为 `-0.926362/-0.333854` 和
`0.949963/0.018510`。analytic `A=16` 的同类量为
`-0.925973/-0.331997` 和 `0.950681/0.0187067`，密度范围违反为
`2.301e-6`，散度和 Poisson 残差分别为 `5.03e-13` 和 `1.60e-11`。
analytic `A=16,w=0.55` 最佳前缀的同类量为
`-0.929168/-0.346311` 和 `0.923439/0.0144411`，密度范围违反为
`7.434e-8`，散度和 Poisson 残差分别为 `4.625e-13` 和
`1.569e-11`。
当前 dynamic/open 问题与慢衰减初值不提供
封闭有限域质量守恒，因此这些数值不能解释为爆破证据；它们恰是晚期
仍需更大计算盒和更高空间分辨率复核的原因。

## 时间步和计算盒对照

本节对照均基于同一 `legacy_abs_gaussian` 映射族。它们能检验旧网格族
内的时间步、计算盒和重叠轨迹一致性，但不能消除原点处非 C1 的
映射正则性缺口。

| 对照 | `rho` 相对 L2 | `omega` 相对 L2 | 壁峰相对差 | 壁宽相对差 |
|---|---:|---:|---:|---:|
| `t=0.005`, `dt=5e-4` 对 `1e-4` | 7.51e-16 | 1.16e-14 | 9.79e-16 | 1.09e-15 |
| `t=0.12`, `dt=2e-3` 对 `5e-4` | 3.92e-10 | 6.43e-10 | 4.53e-11 | 4.49e-10 |
| `t=0.005`, 大盒对小盒 | 2.64e-4 | 7.56e-4 | 2.38e-4 | 1.46e-3 |
| `t=0.5`, 大盒对小盒 | 1.76e-3 | 3.33e-3 | 1.24e-3 | 2.91e-3 |

`t=0.5` 的跟踪峰位置差为 `0.00114`。因此核心区对计算盒的局部敏感性
很小，但这不等同于远场速度已全局收敛。

在 x-only 聚焦与 x/y 均衡聚焦网格的共同可信窗内，`t=0.678`
的梯度、壁峰和壁宽相对差分别为 `3.06e-4`、`-2.37e-4` 和
`6.67e-4`。这证明延长主要来自更晚地耗尽核心分辨率，而非改变早期物理轨迹。
均衡与 concentrated 网格在 `t=0.96` 的同三项差异也仅为
`2.85e-4`、`3.23e-4` 和 `-7.06e-4`。

`A=16,w=0.7` legacy 网格的旧点数窗口拟合低于 `R^2` 候选门槛；这类
窗口又受非均匀输出密度和 profile 搜索范围影响，已由物理时间/e-fold
诊断替代。该轨迹仍来自非 C1 legacy monitor，只能作为历史对照，不能
单独用于可信爆破判断。

证据覆盖时间也必须分开报告：时间步对照只到 `t=0.12`，大小盒
局部场对照到 `t=0.5`，同 513 x 257 节点数的 legacy 映射重叠对照到
`t=0.96`。现已有一条 769 x 385 节点的 q384 独立空间加密轨迹，
它将连续可信终点从 q256 的 `t=1.659748` 提高到
`t=1.751582`。但在这个晚期区间仍无时间步减半或更大计算盒复核，
也没有第三个空间分辨率用于建立收敛序列。已完成的 analytic/legacy
同节点对照仍只是映射正则性复核，应与 q256 -> q384 节点加密分开解读。

## 历史结果不可混用

当前 `k=8` primitive 是 2026-08-29 才切换的 active datum，仓库历史没有
给出它的可信奇性时间。旧的 `t` 约 1.4 结果属于另一个
`initialCondition='degenerate', k=4` physical 算例；`T` 约 0.77 的强候选
又属于 tuned smooth signed-slope datum。CCF heavy-tail 是机制对照，但其
最佳历史长跑也在分辨率失效前给出 `candidate=false`。这些结果不能用来
给当前测试补一个“预期爆破时刻”。

## 复现

```matlab
cd ipm_structured
addpath(pwd)
addpath('research/experiments','research/analysis')

% 以下旧命名算例显式复现 legacy_abs_gaussian 归档结果；
% 它们不是 analytic monitor 的四阶正则性复核。
opts = fourth_order_blowup_case('smooth_box_quadrant_256',struct( ...
    'physicalFinalTime',0.5,'maxDt',0.002,'outputEvery',0.005, ...
    'resultFile',fullfile('result','verification', ...
    'fourth_order_active_smooth_box_quadrant_256_t050_dt002.mat')));
result = ipm.solve(opts);
assessment = ipm_assess_fourth_order_blowup(result);

focusedOptions = fourth_order_blowup_case('focused_box_quadrant_256');
focusedResult = ipm.solve(focusedOptions);
focusedAssessment = ipm_assess_fourth_order_blowup(focusedResult);

balancedOptions = fourth_order_blowup_case('balanced_box_quadrant_256');
balancedResult = ipm.solve(balancedOptions);
balancedAssessment = ipm_assess_fourth_order_blowup(balancedResult);

concentratedOptions = fourth_order_blowup_case( ...
    'concentrated_box_quadrant_256');
concentratedResult = ipm.solve(concentratedOptions);
concentratedAssessment = ipm_assess_fourth_order_blowup( ...
    concentratedResult);

nearTerminalOptions = fourth_order_blowup_case( ...
    'near_terminal_box_quadrant_256');
nearTerminalResult = ipm.solve(nearTerminalOptions);
nearTerminalAssessment = ipm_assess_fourth_order_blowup( ...
    nearTerminalResult);

stableLegacyOptions = fourth_order_blowup_case( ...
    'stable_limit_box_quadrant_256');
stableLegacyResult = ipm.solve(stableLegacyOptions);
stableLegacyAssessment = ipm_assess_fourth_order_blowup( ...
    stableLegacyResult);

% 新的全局解析双高斯复核路径。
analyticOptions = fourth_order_blowup_case( ...
    'analytic_stable_limit_box_quadrant_256');
analyticResult = ipm.solve(analyticOptions);
analyticAssessment = ipm_assess_fourth_order_blowup(analyticResult);

analyticOptimalOptions = fourth_order_blowup_case( ...
    'analytic_moderate_narrow_limit_box_quadrant_256');
analyticOptimalResult = ipm.solve(analyticOptimalOptions);
analyticOptimalAssessment = ipm_assess_fourth_order_blowup( ...
    analyticOptimalResult);

% 384 x 384 空间加密长跑：保存 result/assessment 并绘制增长图。
denseProducts = run_fourth_order_dense_growth();

% 主增长率输出：gamma、gammaDot、瞬时 T/p 及三模型窗口比较。
gradientRateFit = analyticOptimalAssessment.gradientGrowthRateFit;
wallRateFit = analyticOptimalAssessment.wallPeakGrowthRateFit;

% Canonical tau_c 下的线性/指数/曲率守门与双尺度图。
gradientTauFit = analyticOptimalAssessment.gradientCanonicalGrowthFit;
wallTauFit = analyticOptimalAssessment.wallPeakCanonicalGrowthFit;
[tauFigure,tauAnalysis] = ipm_plot_tau_growth(); %#ok<NASGU,ASGLU>

verify_fourth_order_blowup_setup
```

上述默认输出名用于新复跑。已归档证据使用带
`_t050_dt002`/`_t185_dt002`/`_t195_dt002` 等后缀的显式
`resultFile`；assessment 也是单独保存的，不应将默认文件名与归档运行混同。

主要实现与审计入口：

- `research/experiments/fourth_order_blowup_case.m`
- `research/experiments/fourth_order_smooth_peak_grid.m`
- `+ipm/+diagnostics/maximumGrowthRateFit.m`
- `+ipm/+diagnostics/canonicalMaximumGrowthFit.m`
- `research/analysis/ipm_assess_fourth_order_blowup.m`
- `research/analysis/ipm_plot_tau_growth.m`
- `research/analysis/verify_fourth_order_blowup_setup.m`

其中 `verify_fourth_order_blowup_setup` 已成功审计 5 个 analytic 256-cell、
1 个 dense 384-cell 与 1 个 legacy factory，包括完整四阶 tuple、profile/version、网格质量、
reserved metadata 与 continuous-prefix 辅助器契约；已保存结果的实际可信前缀
由 `ipm_assess_fourth_order_blowup` 另行判定。

主要结果位于 `result/verification/fourth_order_active_*.mat` 和
`result/verification/fourth_order_analytic_screen_*.mat`，并由
`result/verification/manifest.jsonl` 记录。当前最佳长跑与解析系列汇总为：

- [`fourth_order_active_analytic_A16_w055_t195_dt002.mat`](../../../result/verification/fourth_order_active_analytic_A16_w055_t195_dt002.mat)
- [`fourth_order_active_analytic_A16_w055_q384_t205_dt002.mat`](../../../result/verification/fourth_order_active_analytic_A16_w055_q384_t205_dt002.mat)
- [`fourth_order_active_analytic_A16_w055_q384_t205_dt002_assessment.mat`](../../../result/verification/fourth_order_active_analytic_A16_w055_q384_t205_dt002_assessment.mat)
- [`fourth_order_active_analytic_A16_w055_q384_t205_dt002_growth.png`](../../../result/verification/fourth_order_active_analytic_A16_w055_q384_t205_dt002_growth.png)
- [`fourth_order_active_A16_w055_q256_vs_q384_growth.png`](../../../result/verification/fourth_order_active_A16_w055_q256_vs_q384_growth.png)
- [`fourth_order_active_A16_w055_q256_vs_q384_tau_growth.png`](../../../result/verification/fourth_order_active_A16_w055_q256_vs_q384_tau_growth.png)
- [`fourth_order_active_A16_w055_q256_vs_q384_tau_growth_analysis.mat`](../../../result/verification/fourth_order_active_A16_w055_q256_vs_q384_tau_growth_analysis.mat)
- [`fourth_order_active_analytic_campaign_summary.mat`](../../../result/verification/fourth_order_active_analytic_campaign_summary.mat)
- [`fourth_order_active_growth_rate_comparison.mat`](../../../result/verification/fourth_order_active_growth_rate_comparison.mat)

## 当前可支持的表述

本轮已经构造出第一象限 256 x 256 和 384 x 384 单元的四阶组件测试。
全局解析双高斯 `A=16,w=0.55` 网格将连续局部 trusted 物理时间从
513 x 257 lattice 的约 `0.0219` 延长到 `1.659748`，为
75.747438 倍。它比 analytic `w=0.7` 基线延长 0.927975%，
但共同前缀上的主要标量轨迹只有小幅变化。进一步的 q384
节点加密将证据终点延长到 `1.751582`，相对 q256 增加
5.533053%。它请求到 `t=2.05`，但在 raw `t=1.855940` 触发
`grid_resolution_failure`；该 raw 终点不是证据终点。时间步、
大小盒和 legacy 族内的网格重分布对照分别只覆盖到
`t=0.12/0.5/0.96`。analytic/legacy 对照是监视函数正则性复核，
不是节点加密。q384 是首个更高节点数独立轨迹，但单次加密仍不能
建立空间收敛阶；晚期的远场与质量警告也仍然显著。

旧的四种点数窗口 `R^2` 门槛均为 `candidate=false`；新的最大值增长率
审计则显示 q256 与 q384 的三个较长梯度窗口均偏向有限时幂律。
q384 的最大梯度/壁峰动态范围是 `1.998490/2.042679` e-fold，
仍不足 3 e-fold。在 canonical `tau_c` 中，q384 的最后
`0.25/0.5` e-fold 窗选择指数最大值，较长窗却触发曲率守门；
末段组合指数率约为 `0.50`。这把原来的物理时间幂律更准确地
重解释为“近 Type-I 物理增长 / canonical-time 指数轨道”的
同一双尺度候选。壁峰关系部分由规范恒等式保证，较长窗率仍漂移，
且晚期时间步、计算盒、形状收敛与第三层空间分辨率复核不足。
因此当前证据支持“与双尺度机制相容的 canonical-`tau_c` 局部指数候选”，
不支持宣称已经找到可信爆破解。
