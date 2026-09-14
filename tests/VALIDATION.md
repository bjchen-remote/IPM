# 结构化副本验证记录

## 2026-09-14：架构整理后的最终回归

冻结 9.14 包与源码构造的 1025/1024 实际节点象限轴、锚点修正轴及其构造信息
均 `isequaln`。变更运行文件和象限测试的 Code Analyzer 为 0 条。
`ipm.verify('quadrant')` exit0（4.814 秒）；`ipm.verify('equivalence','quick')`
的 6 个旧全域物理算例与 3 次迁移逐值一致（2.627 秒）；最终
`ipm.verify('all','quick')` exit0（91.022 秒）。未重新发行、未运行服务器复制主文件预飞，
也未测百万节点 LU 与长时收敛。

## 2026-09-14：复制主文件的绝对路径契约

`server/main_quadrant_release.m` 只用自身路径确定作业目录，固定源码
`releaseRoot='/data/user/hd58131/ipm/ipm_long_time_server_20260914'`。
固定路径字面量、源文件中的真实选项构造、MATLAB R2026a `config.resolve`
与 Code Analyzer 均通过：`1025×513`、`xlim=[0,8]`、象限开启、无旧自动网格；
打包脚本 `bash -n` 通过。
修正后发行包根 `main.m` 与源码入口逐字一致，`SHA256SUMS`、ZIP 完整性及
ZIP SHA-256 校验通过；从发行包运行 `ipm.verify('quadrant')` exit0
（9.525 秒）。
旧自定位候选的数值配置预飞不能替代复制件路径测试，已从正式发行目录撤出。
真正的复制件预飞需要服务器存在上述发布绝对路径；本机尚未执行该测试，
也未运行 525825 节点的 PDE/LU。

## 2026-09-14：象限直接分配与 level-set 独立触发

横轴整数分配只评估 1 次，65×33 原生迁移与 1025 节点一维轴回归通过；
默认 90% 核心触发不受旧综合 `safetyFactor` 支配，非 90% 用户观察层回退
到直接 90% 测量；无正峰的物理态跳过反演。最终 `ipm.verify('quadrant')` exit0（6.541 秒），
`ipm.verify('equivalence','quick')` 的 6 个旧全域物理算例及 3 个迁移逐值一致
（7.380 秒）。`ipm.verify('all','quick')` 在最后的边界调整前 exit0
（141.792 秒）；完整最终版本未重复运行该全套。变更的运行文件 Code Analyzer
仅 `evolve.initialize` 留有原有 MSNU 抑制提示，其余 0 条。
1024² 二维求解/分解性能、长期再网格容量与时空收敛仍未验证。

## 2026-09-14：第一象限直接重网格

显式 `quadrantOnly` 的 `ipm.verify('quadrant')` 在 MATLAB R2026a 通过：
129×65 全域与 65×65 象限动态 RHS 最大绝对差 3.9774e-14；
65×33 单提案原生迁移、非均匀聚焦轴对照（最坏差 1.2434e-14）、
4 步短算、远边界轴分类及 schema-v4 checkpoint 验签恢复通过。
旧全域 equivalence 的 6 条物理轨迹和 3 次迁移仍逐值一致。
`ipm.verify('all','quick')` 在最后一项仅限象限自定义逻辑度量修正前 exit0，
耗时 112.5367 秒；修正后单独的象限专项再 exit0，旧全域路径未触及。
最后修改的 `mesh.build`、象限测试及短例 Code Analyzer 均 0 条；
全运行包 125 个 `.m` 文件仅有未改动 `evolve.initialize.m` 的 1 条旧 `MSNU` 抑制提示。
尚未运行百万节点 LU/长时收敛，不能把短算认证外推到服务器规模。

## 2026-09-05：schema 4 exact-no-feedback 验证

本机 MATLAB R2026a 对当前副本的 129 个 `.m` 文件执行 Code Analyzer，
结果为 0 条告警；全套 baseline 通过。
在补严 checkpoint 的退役字段拒绝表后，又对两个受影响文件执行 Code Analyzer
（0 条告警）并重跑完整 checkpoint/restart baseline；动态与四阶 split-run、
不可信终态抑制及伪 schema-4 `state.rescaling.maxDynamicRate` 负测试全部通过。
统一所有 schema-4 退役字段清单后，最终再次检查全部 129 个维护范围 `.m`
文件（0 条告警），并重跑 config/grid/elliptic/numeric-core/transport/symmetry/
scaling/remesh/smoke/output/checkpoint/integration 共 12 层 baseline，全部通过。

物理内核 equivalence 覆盖空间阶 2/4/6 分别搭配
`double_odd_omega` 和 `half_plane` 的 6 个纯物理模式短程解，
以及阶 2/4/6 的 3 个非平凡 remesh。schema 4 直接更新规范时钟后，
`state.canonicalTime`、`history.common.canonicalTau` 和
`snapshots.canonicalTime` 与旧版只有浮点重组舍入差别。新版还在
`operator.gridQuality.x/y` 中增加了四个只读求积质量诊断：
`quadratureWeightsStrictlyPositive`、
`minimumQuadratureWeightToControlWidthRatio`、
`maximumQuadratureWeightToControlWidthRatio` 和
`maximumQuadratureWeightToControlWidthAbsoluteDeviationFromUnity`。
投影仅排除这三个规范时钟位置与四个新诊断字段；物理场、
`physicalTime`、`rho/omega/velocity`、其余历史物理观测量、尺度、
求积权重、迁移场/指标、数值算子和 Poisson 求解作用均以
`isequaln` 精确一致，没有使用容差。

q256 五候选共用初值/RHS preflight 给出
`c_l=0.486855155934`，`X=1` anchor defect 为 `3.262e-17`，
五个 `omegaGaugeResidual` 均为 0：

| `cOmegaGauge` | `c_omega` | 局部 condition |
|---|---:|---:|
| `gradient_energy` | 0.242828367953 | N/A |
| `anchor_wall_window_l2` | 0.0574064170493 | 0.707558 |
| `anchor_wall_template_projection` | 0.0574064170493 | 1 |
| `anchor_bulk_gradient_l2` | 0.239324557527 | 0.711453 |
| `anchor_wall_window_l4` | 0.0794076810301 | 0.802283 |

两候选 `t=0.1` smoke 中，`gradient_energy` 和
`anchor_wall_window_l2` 均以 72 步结束，分类为 `target_reached`，
且 `hardPassed=true`、`eligible=true`；两者的对称双向 pairwise physical
covariance 门槛也通过。该 smoke 产物位于
[`q256_comega_gauge_tournaments_exact_v4`](../result/verification/q256_comega_gauge_tournaments_exact_v4/tournament_20260905T063951664Z)。

### 正式 q256 五候选规范竞赛

预注册的五候选 group-stage 在 `2337.84 s` 内完成，未触及 `7200 s`
预算边界。五例均以 432 步到达物理时间 `0.55`，均满足
`hardPassed=true`、`productionEligible=true`，没有 designed safety stop、
数值/网格门失败、规范退化或意外崩溃。10/10 对称双向物理协变比较全部完成并
通过；全局最坏差异为 physical rho relative L2 `1.510e-12`、physical omega
relative L2 `3.224e-12`、physical `max|rho_x1|` relative `6.765e-14`。

排名使用冻结的六维字典序向量，而非事后加权总分；首项是尾段
`max(abs(c_l drift),abs(c_omega drift))`：

| 排名 | `cOmegaGauge` | 首项最坏 drift | 最大 detrended RMS | 终点 `(c_l,c_omega)` |
|---:|---|---:|---:|---:|
| 1 | `gradient_energy` | 0.0390410 | 1.56868e-4 | (0.532713, 0.265508) |
| 2 | `anchor_wall_template_projection` | 0.117985 | 3.75635e-3 | (0.471819, 0.0386551) |
| 3 | `anchor_wall_window_l2` | 0.399637 | 5.42665e-3 | (0.469538, 0.0231403) |
| 4 | `anchor_bulk_gradient_l2` | 0.606956 | 9.52427e-3 | (0.499906, 0.0673326) |
| 5 | `anchor_wall_window_l4` | 0.679437 | 5.50459e-3 | (0.472483, 0.0215006) |

可信历史上的最坏 `X=1` transport-anchor 规范化缺损为 `2.533e-16`，
最坏通用 omega-gauge residual 为 `1.139e-16`；每条历史均逐点满足
`lengthScaleGain=1`、`widthControlMode=0`。五例共同终态为
physical `max|rho_x1|=0.901020`、`||grad rho||_inf=0.945037`、严格 core
计数 `75.96/72.37`、safety `0.120316`。共同远边界 velocity/omega 比
`0.9322/0.03339` 仍超过告警线，所以本结果只证明无反馈规范的 q256
相对稳定性排序，不构成爆破或 q512 生产可信性证据。

不可变证据见 [group summary](../result/verification/q256_comega_gauge_group_stage_v4/group_20260905T064403198Z/group_stage_summary.mat)、
[tournament summary](../result/verification/q256_comega_gauge_tournaments_exact_v4/tournament_20260905T064403596Z/tournament_summary.mat)
与 [rate plot](../result/verification/q256_comega_gauge_tournaments_exact_v4/tournament_20260905T064403596Z/q256_gauge_tournament_rates.png)。

### q512 大箱 exact-gauge/CFL 冻结

长度规范继续用 `transport_anchor`：固定 `X=1`，
`c_l=-U_1(1,0)`，因而 `V_1(1,0)=0`。新候选
`wall_omega_quadratic_peak` 在 wall `R_X` 当前极大节点及两个邻点上
构造非均匀三点二次插值，在其严格内部顶点用同一组权计算
`P=Q[R_X]`、`F=Q[B_X]`，并直接取 `c_omega=-F/P`。所以
`P_tau=F+c_omega P=0` 是瞬时代数约束，不是 peak feedback。全部
本轮运行都直接积分 `tau`、`timeSpeed=1`，没有 limiter 或平滑增益。
`none` 仅是预设 `c_omega=0, C_omega=1` 的 control，其恒零速率不计作
经验收敛。

初始 RHS box screen 在 q512、`H=281,1e3,1e4,1e5,1e6` 上不推进
PDE。`H=1e5 -> 1e6` 时，`gradient_energy` 的全局能量由
`1.79572e5` 增至 `1.79579e6`，forcing 由 `-4.38572e4` 增至
`-4.38591e5`；对非衰减 `degenerate_primitive` datum 明显是 extensive，因而
即使比值给出 `c_omega=0.244232359244` 也被否决。局部二次顶点在
`H=1e6` 给出 `c_l=0.488462559170`、`c_omega=0.125042865249`、
`P=0.686073738019`，并通过所有几何、条件数和残差门。正式证据为
[box-screen summary](../result/verification/q512_comega_initial_box_screens_v4/screen_20260905T110530211Z_tpbe425904_b6b4_4b3f_9aa4_c02b888cd290/comega_box_screen_summary.mat)。
该作业只筛规范，不证截断收敛；`H=1e6` 的 far-velocity ratio 仍为
`0.986819`。

随后的四组 fresh-root 短算在同一 `H=1e6` 网格上交叉
`wall_omega_quadratic_peak/none` 与 `CFL=0.25/0.50`，积分至
`t=0.03`，`maxDt=0.004`。4/4 运行和全部成对比较都通过硬门。
二次顶点两个 CFL 分别用 `40/21` 步、`184.672/106.776 s`；
`none` 用时 `185.350/106.231 s`，总墙钟加速约 `1.74x`。同规范
CFL 比较的最坏 physical rho/omega relative L2 为
`8.003e-12/2.963e-9`，physical peak 相对差小于 `4e-14`；同 CFL
跨规范最坏 rho/omega 为 `1.156e-11/2.051e-11`。因此后续长算冻结
`CFL=0.50`。证据见 [gauge-by-CFL summary](../result/verification/q512_fresh_gauge_cfl_short_screens_v4/screen_20260905T114636013Z_tp3e066220_eec8_4c94_a0c7_b02a376fc67f/gauge_cfl_screen_summary.mat)。

q256 二次顶点轨迹延至 `t=1.7` 时，终点
`(c_l,c_omega,kappa)=(0.238287,-0.266436,0.504723)`。最后 `0.125` 物理时间
窗上，物理二次顶点的逆量 affine 拟合斜率为 `-0.748374`，代数
平均斜率为 `-0.747521`，相差 `0.114%`，且 RMS/逆量下降量为
`0.113%`。这只是 inverse-peak 的 Type-I candidate：同窗口
`c_l/c_omega/kappa` 仍漂移 `19.4%/14.2%/3.08%`，且只覆盖 `1.861`
e-folds，未达 1% 稳定性与 3 e-folds 门。因此
[q256 result](../result/verification/q256_quadratic_peak_continuations_v4/continuation_20260905T120651140Z/result.mat)
不支持 blowup 或渐近常速率声明。

q512 fresh root/continuation 现在允许根作业显式选择
`gradient_energy` 或 `wall_omega_quadratic_peak`，但从 checkpoint 续算时必须
冻结并通过同一 gauge 的全链条审计。真实 q512 二步根作业
[hard pass](../result/verification/q512_fresh_exact_gauge_roots_v4/root_20260905T122642961Z_large_box_smoke_tp5f74d6c0_96d5_4c8a_b18f_88b19ea00f92/root_summary.mat)；在此基础上的最终
[continuation summary](../result/verification/q512_schema4_continuations/continuation_20260905T124215721Z_tp1b1f9942_ef22_4494_b5f7_ee4ca7b4a77c/continuation_summary.mat)
为 `completed_endpoint`，segment/total steps 为 `2/4`，墙钟 `27.9553 s`，
`postAudit.hardPassed=1`、`parentHistoryExactPrefix=1`。两次早期 postrun
rejection 定位到空 `anisotropic` 组的 `fieldnames/isfield` 空形状边界，
修复为逐字段 `cellfun` 存在性检查并将空标量组视作 exact prefix。
这是 artifact audit bug，不是数值不匹配。

发射前对当前 202 个 MATLAB 文件逐一执行 Code Analyzer，结果为
`warnedFiles=0, diagnostics=0`；`ipm.verify('baseline')` 的 12 层验证
全部通过。长算 preflight 冻结为 q512 `1025 x 513`、初始
`H=1e6`、`wall_omega_quadratic_peak`、`CFL=0.50`、`maxDt=0.004`、
`t_final=1.7`、`tau_final=4.5`、`maxSteps=4000`，其 [summary](../result/verification/q512_fresh_exact_gauge_roots_v4/root_20260905T124547741Z_large_box_campaign_tp2922f311_a512_4591_8ba1_898027e4bd92/root_summary.mat)
的 mesh/contract 门均通过。

2026-09-05 20:47 CST 已用带并发锁的后台会话
`ipm_q512_quad_h1e6_t170` 启动长算。[active root summary](../result/verification/q512_fresh_exact_gauge_roots_v4/root_20260905T124741359Z_large_box_campaign_tp359ea53d_96ea_4ca8_b160_17e0fe61f50f/root_summary.mat)
及同目录 `root_manifest.jsonl` 是权威状态记录；启动时 manifest 已写入
`root_solve_started`。`H=1e6` 仅是初始物理半宽，历史还逐点记录
随 `C_l` 变化的实时物理覆盖范围。运行中状态不作为 blowup 或
渐近结论。

## 2026-09-01：初始结构化副本验证

验证时间：2026-09-01。环境：本机 MATLAB R2026a Update 5
（`26.1.0.3346908`），从 `/private/tmp` 启动干净会话，仅添加本副本根目录。
保留版本：`sixth_order_integration@65f999b7e205dd4decef06c937d1aa003cf8b744`。

| 检查 | 本次结果 |
|---|---|
| 原版本保留 | 清单内 140 个来源文件的 SHA-256 全部未变；其余旧工作目录未编辑 |
| 核心迁移审查 | 78 个核心文件仅变更函数名、声明和调用名，数值表达式及运算顺序保留 |
| 独立目录 | 103 个迁移函数均解析至本副本；求解器的 75 个静态依赖不指向旧版、研究或测试目录 |
| `ipm.verify()` | 基础 11 层全部通过，包括输出保护、重网格回退与快照配对 |
| `ipm.verify('fourth')` | 完整四阶套件通过，包含 quick 稳定性扫描和 320 步物理演化 |
| `ipm.verify('sixth')` | 完整六阶套件通过，包括 WENO7、RK6、网格准入、迁移及物理收敛 |
| `ipm.verify('equivalence')` | 11 组完整演化结果与 3 组非平凡重网格迁移均精确等值 |
| MATLAB Code Analyzer | 新副本全部 125 个 `.m` 文件零消息；78 个对应原核心文件也为零 |
| 三个短算例 | 均以 2 步到达 `t=2e-4`，停止原因为 `physical_final_time` |
| 绘图入口 | 终态图 4 个坐标轴、单步实时图均正常；隐藏窗口测试后全部关闭 |

一致性使用 `isequaln`，不放宽数值误差：仅剔除每次生成的
`metadata.caseId` 与 `metadata.createdAt`。其余配置、场、物理重建、网格、尺度、
历史、质量、拟合、椭圆信息与快照一并比较。覆盖五种基线模式，以及四/六阶的
物理、各向同性动态和各向异性固定规范模式。
三组迁移分别检查二/四/六阶的场、指标、算子和 Poisson 分解求解结果。
原版与新版在动态模式下保留相同的末步物理时钟行为。

`ipm.verify` 完成后恢复原 MATLAB 路径；一致性检查也恢复工作目录。
独立会话中原 `main_ipm` 不在路径上，显式一致性检查结束后仍不在路径上。
本次原始日志与报告保存在
[`result/verification/structured_copy_20260901`](../result/verification/structured_copy_20260901)。

### 范围

本次验证证明重组后的副本维持已检查的数值行为，没有扩大原方法的适用范围。
未运行四阶 `heavy` 扫描、长时间生产算例、视频编码或任意非均匀网格认证。
三个小网格示例仅用于运行演示，其远场截断诊断提示仍保留，不代表生产精度。
四/六阶整体收敛证据仍限于原套件覆盖的光滑固定网格问题；
不据此声称任意动态规范、触发式重网格或非光滑解的整体高阶收敛。
参考 `iclm` 的目录组织来自可读目录树；OneDrive 占位文档正文未能读取。
