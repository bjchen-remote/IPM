# 无 LU 启动预飞与自动续接审计

2026-09-09；只新增 performance_lab 研究入口，未修改生产 `+ipm`、当前 driver、admit 或已冻结主线。

`ipm_perflab_preflight_fresh_campaign(plan, newReportDirectory)` 已实现并完成 8 项回归。它检查启动前提，**不证明从原问题 t=0 连续稳定到长时，也不保证未来仍有可接受网格**。

## 接口和明确边界

```matlab
plan = struct('sourceRoot', absoluteFrozenSource, ...
    'outputDirectory', absoluteNewCampaignDirectory, ...
    'startKind', 'time_zero', 'entryPoint', 'ipm.solve', ...
    'options', actualFlatOptions, ...
    'additionalEntries', {{'run_fresh_profile_campaign'}});
report = ipm_perflab_preflight_fresh_campaign(plan, absoluteNewAuditDirectory);
assert(report.passed);
```

- 要求默认 10 线程、绝对路径和新的审计/计算目录。计算目录保持未创建，只在已存在父目录做临时写入/清理检查。结果落为注册 plan MAT、报告 MAT/JSON、源码 SHA256 闭包；调用结束恢复原 cwd/path。
- 暂时切换到冻结 sourceRoot，显式加入对应 longtime/grid/experiments 路径。对注册入口的维护命名空间调用和句柄进行保守传递扫描，再用 MATLAB `requiredFilesAndProducts` 补充依赖。每项 `which` 必须位于 sourceRoot；外部同名函数不能填补缺文件。实际受测闭包 109 个源码文件，没有动态调用点或外部运行依赖。
- 对任意运行时构造函数名/未注册 callback 不声称静态依赖完整。`callbackDependencies` 可以增加注册依赖；原 t=0 的初值采样只执行原生预定义初值，任意函数句柄初值会明确拒绝并要求另行无 LU 审查。
- `time_zero` 要求唯一入口 `ipm.solve`，用维护 `config.resolve` 一次校验实际 options，再用实际显式 `customX/customY` 和维护 `initialDensity` 检查初值尺寸、有限性及范围。没有调用 `initialize` 或 `mesh.build`；解析初始重网格和实际初始流场仍须独立验收。
- `fresh_checkpoint` 需要 `bootstrapResultFile`、`initialDataFile`、`qualificationFile`，调用原 `fresh_profile_campaign_review/admit` 严格读取原生 checkpoint，逐值核对状态、全部历史、时钟、尺度、快照、配置、参考量及 fresh lineage。该模式明确 `originalTimeZeroInputRecognized=false`。链式 admission 的证据目录只替换为预飞自己的新目录；原 qualification、checkpoint 和 campaign 不写入。
- 当前作用范围是 schema4/high_order/exact_gauge_no_feedback_v1 的对称固定节点数 campaign。实际轴必须保留端点、严格单调、完整偶对称和精确 ±transportAnchorX。原轴门未放宽：ratio≤1.08、log-curvature≤.01、stencil-rcond≥1e-9、min quadrature/mean≥1e-8、权重严格正、weight/control∈[.35,1.65]。
- 可选 `axisProbes` 为结构数组，每项包含 `label/functionName/arguments/geometry`。只允许维护纯轴工厂 `mesh_core_patch_axis`、`mesh_general_rounded_axis`、`ipm_gridlab_equalize_axis`。执行实际构造，以另一条原轴配对作完整轴质量检查。每个注册方向至少一项通过。**这不是 frozen-field 的 core/front、峰跳、守恒、范围或实际 transaction-ready 验收**；这些原门继续由 designer/transaction 执行。

纯轴回归示例：

```matlab
plan.axisProbes = struct('label','uniform rounded constructor', ...
    'functionName','mesh_general_rounded_axis', ...
    'arguments',{{1,2,49,1/12,8,.5,2}},'geometry','x');
```

## 实际验证

受测源码完整冻结于 `result/verification/performance_lab_20260908/campaign_preflight_v2_source`，184 文件 SHA 清单在运行后逐项验证无变化。报告 `campaign_preflight_v2/regression.json`，原始执行 profile 在同目录 `regression.mat`。

| 场景 | 实测结果 |
|---|---|
| 完整冻结源码，49×25 实际初值/原门轴/rounded 工厂 | 通过 |
| 历史 `fresh_adaptive_source_v1` 缺 `mesh_general_rounded_axis` | 启动前 `PreflightDependencies` 拒绝，准确列出缺文件 |
| 两个真实构造调用均没有可行 cell budget | 保存各自失败，并以 `PreflightAxisFamiliesExhausted` 拒绝 |
| y 轴质量失败 | `PreflightAxisQuality` 拒绝 |
| 已存在输出目录 | `PreflightOutputExists` 拒绝 |
| 相对输出路径 | `PreflightOutputPath` 拒绝 |
| 用 late-fresh driver 冒作原问题 t=0 | `PreflightTimeZeroEntry` 拒绝 |
| 实际 fresh step3668，parent-equivalent τ=10.0432007518592 | 原 14-link 链全部重读通过；完整 history/pair/reference/axis 门通过，core=25.0091876/30.6492393 |

上述全部在一个人为 cwd/path 冲突目录中调用：该目录提供同名 `mesh_general_rounded_axis`，执行即报错。完整源正确使用冻结实现；历史 v1 虽然能从 ambient 找到同名函数，仍因其不在 sourceRoot 被拒绝。各次调用后 cwd/path 精确恢复。

8/8 符合预注册预期，MATLAB exit0。执行 profile 确认 `solve/mesh.build/poisson/velocity/initialize/flow/restoreCheckpoint/writeCheckpoint` 均未调用，零 LU、零 PDE、零 native checkpoint 写入。v1 7 项旧报告保留。v2 实测后仅调整一处分析器 AGROW 注释和 if/elseif 排版；最终三个 live helper 的 checkcode 单独保存为 `campaign_preflight_final_checkcode.mat/.log`。

## 当前批次运行与恢复边界

1. `run_fresh_profile_campaign` 是晚态缩箱后的独立 fresh IVP；其 fresh time0 并非原问题物理 t=0。`mesh_design_fresh_core_patch` 同样要求 latePhysicalBoxBranch 的父 crop 来源。通用 t=0 网格应复用纯几何工厂，而不能伪造该 lineage。
2. 现 driver 完成 maxStages 后返回 `batch_review_required`；异常/候选耗尽返回 `review_required`。保存 accepted native checkpoint 很完善，但目前没有自动创建下一批、重新入场和持续运行到目标的外层状态机。无人运行需要在不改变数值门的前提下接通这一层。
3. 每段都读完整 trusted history 并核对源前缀、参考量；每次新批 admission 会重读整条已注册链。它们不需要 LU，但 I/O/哈希/审计工作随历史增长。可缓存不可变证据摘要以减少重复读取的设计，必须另有文件 SHA/身份验证，不能用上次 passed 标志代替当前原生检查。
4. 如果 transaction 已提交但随后 solve 未产出配对 result，driver 把 latestPairedResultFile 留空，同时保留有效新 CP。当前 fresh 入场接口要求 CP 与 result 配对；未来恢复器须显式处理 `transaction_committed → solve_pending`，依据 CP、transaction audit 和 manifest 恢复，不能重新做一次同时间迁移或跳过该阶段。
5. 数值候选失败与程序依赖失败应分别分类。当前 try_regrid 捕获异常并尝试后续 family，历史缺依赖因此重复四次。预飞已提前拒绝此类错误；本次没有改变 driver 的异常处理或任何候选门。
6. checkpoint 的 payload 保存完整状态、尺度、冻结参考、日志和调度 cursor；写入使用新文件及临时安装，manifest 追加失败不会删除已有效写入的 CP。无人恢复应先严格验签和原生身份/前缀，再从不可变文件续算，不覆盖旧成果。

## LU 生命周期与下一项等价优化

**研究外部 transaction。** `ipm_gridlab_regrid_checkpoint` 在旧 CP 严格 restore 后，先验证候选轴，随后 `rmfield(state.ops,'poisson')`，再保留 oldState 并 transfer/build 新缓存。旧 transfer/audit 所需的是导数、求积、几何、参考与旧流场；RHS cache 不含 ops/decomposition。返回只有新 CP/audit，局部新 LU 在返回/异常展开后释放；后续独立 solve 再 restore 会重建它。若多个 candidate 进入真实 transaction 后失败，每次会重新构建源和候选缓存。

**真正在线自适应。** `advance` 的 `acceptedState=state` 仍引用旧 ops，用于 RK/迁移后非有限回退；`remeshIfNeeded → adapt → transfer` 同时构建新 LU。只删除某个局部 ops.poisson 不会释放 acceptedState 中的引用。若上一候选 transfer 返回但 validate 拒绝，`candidateOps` 也会保留到下次赋值/函数退出，下一次 build 可能与该失败缓存重叠。后者可先独立验证及时清理失败局部变量；消除旧 accepted LU 则需要完整的 rollback 数据与按需重建策略，不能只删字段宣称解决。

**restore 的潜在双 build。** `restoreCheckpoint` 先 `mesh.build(config)`，再比较坐标；不同则用 saved.x/y 的 gridOverride 第二次 build。原生在线 remesh 保持初始冻结 config，因此恢复时常会先构建废弃初始 LU；现研究 gridlab CP 已把 custom 轴写入 config，通常只有一次。

可检验候选是：在不装配 Poisson 的情况下，复用维护轴生成算术确定旧实现是否会进入第二分支；只有确定不同，才直接使用 saved axes 构建一次。**不可对所有 CP 无条件 custom override**：未迁移的 analytic stretched/uniform 网格使用解析 metric，custom 分支使用 FD metric，legacy 的 Dx 构造也不同。坐标相同不保证这些算子及 RHS/CFL 逐值相同。

下一步等价测试须覆盖原始 uniform/stretched、自定义初态、1/多次在线 remesh、研究 gridlab CP。检查全部 ops 数组/导数/metric/求积、flow/RHS/cache、尺度、log/cursor、reference 和 metadata，随后原 same-grid split-run 验证；不能只检查坐标、全局残差或最终 profile。本次仅定位和设计，未修改/运行该生产优化。

## 调度研究 helper 的独立复核

`mesh_autonomous_decision` 不改场、CFL、时钟或历史；显式 timeUnit 匹配、只在当前 remesh epoch 估计趋势、保留每接受步检查要求均正确。已向 root 指出 v1 先按固定 .2 预测、后缩短步数预算 horizon 会造成保守提前 regrid。root 已在新的独立版本改为实际计划 horizon，并添加了该负例及 endpoint-clipping 的步长预算回归；原 v1 证据不覆盖。这里不把该调度研究或静态预飞当作自动网格的时间/空间/箱体收敛认证。
