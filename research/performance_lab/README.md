# 性能实验室

这里的代码只做有界微基准和候选内核比较。生产时间推进仍只有 `ipm.solve`；
研究目录不会被求解器调用。不改变 schema-4 `exact_gauge_no_feedback_v1` 规范、
可信门、RK 系数、CFL 或当前生产 checkpoint。

## 使用

从 `ipm_structured` 的父目录启动 MATLAB 后：

```matlab
addpath('ipm_structured');
addpath('ipm_structured/research/performance_lab');

% 默认小型 WENO 微基准，不读文件、不推进 PDE。
report = ipm_perflab_benchmark_weno5();

% 必须显式传入短小 opts 或原生 schema-4 checkpoint。
% 已有 checkpoint 保持创建它时的 MATLAB 计算线程设置。
report = ipm_perflab_benchmark(opts, ...
    struct('repetitions',2,'solveSteps',2,'profileSolve',true));

% 唯一允许跨线程读取既有 checkpoint 的实验：先固定 10 线程恢复，
% 随后仅在内存中比较 flow/advance，结束或异常时恢复 10 线程。
report = ipm_perflab_thread_probe(checkpointFile, ...
    struct('threads',[10,1,2,4],'repetitions',1, ...
    'finalTime',4,'physicalFinalTime',10));
```

`ipm_perflab_benchmark` 的默认 `solveSteps=0` 只测冻结状态。
启用 `profileSolve` 会替换 MATLAB profiler buffer，计时含 restore/finalize；
`actualSolveSteps` 记录实际完成步数，不能把安全门提前停止当作完成请求步数。
`makeCheckpoint` 的测量包括事务打包、签名和校验，不含磁盘压缩、写入或 manifest。
微基准函数不保存结果；调用者可把报告写入 `result/verification/` 的新目录。

## 已执行证据（2026-09-08）

`result/verification/performance_lab_20260908/` 保存本次 MAT 证据。
所有计时都是共享机器上的局部样本，不是独占资源的稳定墙钟指标。

- `weno5_block_singlethread.mat`：在 `maxNumCompThreads(1)` 下，
  WENO5 独立 face 分块候选对 double/single、extrapolate/reflect、
  平滑/常量/零/跨幅值/间断数据的每个方向 100 项比较全部逐位一致；
  q512 两个方向的所有定时结果也逐位一致。
  `513×1025` 上 block 32 为 `15.317 ms`，原核约 `17.886 ms`，`1.168x`；
  `1025×513` 上 block 16 为 `15.024 ms`，原核约 `15.838 ms`，`1.054x`。
  每档 5 次取中位数。分块收益依赖方向、线程和数组尺寸，不据此更换生产核。
- `frozen_smoke.mat`：`17×33`、四阶 physical、固定网格、2 步
  `ipm.solve` 完成。已有 RHS cache 路径与无 cache 路径的 rho/flow/scale/cache
  全部逐位一致，`rho*Dx'` 与 `flow.source` 逐位一致。该小网格的
  cached/uncached step 中位为 `14.337/20.190 ms`，仅是 smoke 数据。
  同时验证了 profiler 和 checkpoint 打包测量入口。
- `q512_thread_probe.mat`：在 10 线程严格恢复真实 step-1818 checkpoint 后，
  每档只测一个 flow 和 advance。10/1/2/4 线程 advance 分别为
  `6.468/5.369/4.389/4.348 s`；RHS、步后 rho、scale 的比较差异均为零。
  首档 flow 计时异常偏大，不能据此宣称线程加速。
  两次/档逆序复测的已完成部分反而给出 10 线程 `4.714 s`、4 线程
  `4.925 s`。因机器内存压力，第二轮在母任务要求下已取消；完整复测没有
  产物，已有终端输出保存在 `q512_thread_probe_repeat_partial.log`。
  **线程加速证据不足，生产维持 10 线程。**

`ipm_perflab_signature_probe` 已写出，用于定位 1/2/4 线程下现有签名的
确切差异字段；`ipm_perflab_bridge_threads` 已写出，用于源环境严格恢复后
在目标线程环境用现有接口进行 0 时间重序列化、严格回读和 payload 比对。
因数值计算现已要求统一串行排队，这两个 helper 尚未执行验证，不能交给生产使用。
它们不修改签名算法或放宽兼容/可信门。

命令行 `-singleCompThread` 在这台 MATLAB R2026a Mac 启动时曾报 Qt/neon
错误；成功的单线程实验使用进程内 `maxNumCompThreads(1)`。
已有 checkpoint 签名包含浮点归约，改变线程数后直接读取可能校验失败；
本实验没有修改签名实现，也没有给它放宽容差。

## 审计结论与优先级

1. **先测线程数，而后改内核。** 同一冻结 q512 状态分别测
   `flow` 和 `advance`，才能区分 WENO、稀疏 LU 回代、多线程开销及输出成本。
   线程变化属于浮点路径变化，必须记录 rate/field 差异，不能称为逐位等价。
2. **保留已实现的三项收益。** `stepRk` 已复用接受态 RHS，SSPRK54 每连续步
   从 6 次 flow 降为 5 次；`rhs` 已按规范惰性跳过无用的物理 transport；
   high-order `biotSavart` 已复用 `poisson` 的 residual。这些不是本轮新收益。
3. **低风险消重。** `diagnostics.measure` 第一个 `rho*Dx'` 可复用
   `flow.source`；`velocity` 就以完全相同表达式生成它，接受态/流场/网格配对是
   已有契约。需要用输出和 checkpoint 测试检查全部诊断逐位一致。它只减少
   记录开销，不会减少每个 RK stage 的主成本。
4. **分块 WENO 是有限收益候选。** 高阶 transport 每个 RHS 有四次 WENO5
   调用，原核按 face 各做两个小重构。候选只把独立 face 放入第三维，保持
   每个元素的算术式及分组顺序，避免新近似、阈值或限幅。单核大网格微测
   仅支持约 5%–17% 的核级收益，完整步收益更低；应在选定线程数后复测。
5. **控制重复初始化和输出成本。** q512 固定网格 LU 已在 `mesh.build`
   缓存，isotropic 每步不会重复分解；频繁 restart/regrid 会重建它。
   `record` 动态增加很多标量列，checkpoint 又签名/校验整段历史，极长历史
   下可能增长。先用实测分离 `restore`、`record`、`makeCheckpoint` 与真实
   磁盘 I/O，不削减当前可信诊断或放宽 checkpoint cadence 来制造加速。
6. **anisotropic 属于单独研究。** 每阶段变化的 kappa 使其重新组装并
   factorize 椭圆算子；固定 Tx/Ty 的 Kronecker 部件可缓存，但复用不同
   kappa 的 LU 或把各向异性 PCG 引入 high-order 都不是等价改动。

本目录没有修改共享求解器或 CHANGELOG。后续默认10线程的真实895×386缓存
ABBA已补齐四kernel/RHS/整步bitwise，但整体收益未达5%门，见
[WENO_ACTUAL_ABBA.md](WENO_ACTUAL_ABBA.md)。没有进入新旧连续轨迹/output/checkpoint验收。

另有已通过 8 例小网格逐位验证的换网格内存候选：在创建新分解前释放旧
`ops.poisson`。详见 [REGRID_MEMORY.md](REGRID_MEMORY.md)；尚无 q512 RSS/墙钟证据。

Poisson 行均衡候选已完成 27 组小网格 MMS 和 6 组受控行单位压力测试，
结果为 `row_equilibration_small.mat`。二次幂缩放的所有 MMS 流函数/速度
最大数值差为零；极端行尺度下 RCOND 改善明显，但原方程误差仍在相近
舍入水平，尚无真实网格精度或速度提升证据。最小集成方案、完整 RHS 和
原始残差要求、适用分支和计时限制见 [ROW_EQUILIBRATION.md](ROW_EQUILIBRATION.md)。
三个相关 MATLAB 文件均已通过 `checkcode`，未修改共享求解器。

Kronecker-sum 的缓存实 Schur/快速对角化初筛已完成，结果为
`sylvester_small.mat`：72 项维护网格 MMS/Green 比较通过，但 H=1e6 原始
差分轴的拉伸阶梯揭示显著小行误差，标准 dense 候选当前拒绝进入 q512。
详见 [SYLVESTER_SCREENING.md](SYLVESTER_SCREENING.md)，含 28 项阶梯比较、
分量后向误差、已知离散解、实际微计时与缓存复杂度；三个文件 checkcode 为 0。

原生 checkpoint 时间步加密入口已通过 tiny physical/dynamic 正例和过短
canonical 上界负例，详见 [TIME_REFINEMENT.md](TIME_REFINEMENT.md)。两支
冻结同一物理终点，保持同网格、逐值历史前缀和完整可信。随后授权的真实
q512 step2226 窗口也已通过：粗/细6/11步、物理时间差1.199e-14，lab三场
L2差至多1.658e-12；原始产物与告警全部保存，进程正常退出。

只读三个 q512 v2 末态的 `alpha(Y)X` 背景输运拆分，使全域余量 rate 降低
6.20x/9.37x/10.19x，原 rate 与维护 selector 逐值一致。它没有执行密度重映射
或修改 CFL，不是整步加速结论；完整 alpha 剖面和下一步误差问题见
[ANCHOR_TRANSPORT.md](ANCHOR_TRANSPORT.md)。

冻结背景回溯插值进一步完成 96 组合成、16 组 frozen stage2 和 28 组
恒等/常数/边界检查。完整场 poly6/poly8 约17–23ms；误差、半群、逆向、
L∞与文献限定见 [BACKGROUND_REMAP.md](BACKGROUND_REMAP.md)。没有运行新的非线性PDE。

后续非线性重映射中点研究已完成两档全可信小网格的光滑段二阶验证，并保留了
65×33 二次峰模板切换导致的跨事件失阶反例。候选通过完整 `F-B` 和 remap
二阶缺损修正保持原半离散方程；原全速 WENO 余项仍有显著谱尺度，不能直接用
余速 CFL 放宽稳定步。详见 [NONLINEAR_REMAP_MIDPOINT.md](NONLINEAR_REMAP_MIDPOINT.md)。
没有修改生产推进或 C 规则，当前不批准集成或 q512 测试。

stage006 tau6.0932 的冻结平台预算筛查完成 15 项、9 项 ready：fine96/112/128
的 .85/.90/.95 锚点比例均通过，fine144/160 被原 ratio≤1.08 门拒绝。
已有 fine112/.90、fine128/.90 两档后备，详见
[PLATFORM_BUDGET_STAGE006.md](PLATFORM_BUDGET_STAGE006.md)；没有进行 LU、PDE 或事务。

当前轴单跳筛查沿用维护 equalizer/转移评分，保留 immutable root 和真实深度审计，
显式固定精确 X=1，见 [CURRENT_AXIS_SINGLE_HOP.md](CURRENT_AXIS_SINGLE_HOP.md)。
在同一最新末态的有界非单调幅度搜索找到两项完整冻结门通过的候选，并确认
sigma=.25 的维护幅度搜索漏解，见 [NONMONOTONE_CURRENT_AXIS.md](NONMONOTONE_CURRENT_AXIS.md)。
最新v4 step3326的独立恢复筛查已找到合格x，但18条y family及幅度路径均未达到原21门，
见 [V4_STAGE002_GRID_RECOVERY.md](V4_STAGE002_GRID_RECOVERY.md)。
同一latest epoch的三档静态盒阶梯与五支matched-physical协议已通过原门，
见 [LATE_BOX_LATEST_EPOCH.md](LATE_BOX_LATEST_EPOCH.md)；仍只构成短区间独立缩箱证据。
较晚实际native参考终点的三支natural-epsilon续算入口与tiny正/负验收见
[LONG_BOX_ACTUAL_ENDPOINT.md](LONG_BOX_ACTUAL_ENDPOINT.md)，保持原门、精确历史前缀与新IVP lineage。
最新fresh895的无LU成本分解、Green复用负结果及公开LU三角缓存的小网格精度拒绝见
[FRESH_RHS_COST_AND_LU_SCREEN.md](FRESH_RHS_COST_AND_LU_SCREEN.md)。未找到新的合格等价加速，不放宽bitwise或误差门。

已有缓存上的两次原flow计时入口已通过tiny逐值验证，接入网格实验已有restore窗口，
不另建大LU。互斥函数/行级成本的定义、检查及profiler开销限制见
[FROZEN_FLOW_PROFILE.md](FROZEN_FLOW_PROFILE.md)。

冻结率probe的复用条件和墙面单行负例见
[RATE_PROBE_REUSE_AUDIT.md](RATE_PROBE_REUSE_AUDIT.md)。最新step328的原组装两档扰动
把异常90%宽度响应定位到全线LF最大值折点，墙y项影响很小，见
[LATE_RATE_RESPONSE.md](LATE_RATE_RESPONSE.md)；只诊断，不改变原C、推进或失败门。
