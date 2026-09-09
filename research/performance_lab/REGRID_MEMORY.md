# 换网格时旧 Poisson 分解的存活期

建议的最小改动是在 `ipm_gridlab_regrid_checkpoint` 校验候选轴后、
`oldState = state` **之前**执行：

```matlab
state.ops = rmfield(state.ops,'poisson');
oldState = state;
```

`ops.poisson` 是 `ipm.mesh.build` 创建的 `decomposition` 对象。
`high_order` 与 `sixth_order` 使用 LU；`legacy_second_order` 使用加权矩阵的
Cholesky。`ops.A`、`ops.Tx`、`ops.Ty` 则是稀疏数值矩阵，不是分解对象。
本候选只移除 `poisson`，其余字段保留。

原函数先复制完整 `state` 给 `oldState`，随后 `ipm.remesh.transfer` 中的
`ipm.mesh.build` 构造候选分解。因此旧分解被两份状态持有，并一直存活到
事务函数返回。只在 `transfer` 内删除局部 `ops.poisson`，或复制 `oldState`
后只从 `state.ops` 删除它，都不能消除外层状态对旧分解的引用。

schema-4 checkpoint 不保存 `ops` 或 RHS cache。`restoreCheckpoint` 返回后，
其局部状态/算子变量退出作用域；调用者此时只有 `state.ops.poisson` 持有这个
恢复时新建的旧网格分解。`state.flow.poissonSolveInfo` 只有求解诊断，
`state.rhsCache` 只有坐标、场、尺度及 RHS/scale-rate，都不含分解对象。
在上述位置移除字段即可在新 `mesh.build` 前消除旧分解的存活引用。
MATLAB 内存池何时把页交还操作系统属于另一问题，不能仅凭字段删除宣称
RSS 立即下降；但内存分配器有机会重用已释放的分解存储。

## 旧字段依赖审计

| 消费者 | 从旧 `ops` 读取的字段 |
|---|---|
| `validate_axes` | `x,y,nx,ny` |
| `transfer` 所有分支 | `x,y,symmetryMode,rescaling,baseX,baseY,remeshCount` |
| `transfer` PCHIP 分支另需 | `rescalingMode`；动态模式还需 `hx` |
| `transfer` 六阶分支另需 | `integrationHx,integrationHy` |
| `audit_transaction` | `integrationWeights,Dx,nx,ny,remeshCount` |

以上均不读取旧 `poisson`。旧 `flow` 也不参与该事务的 audit。
候选的新流场仍通过原 `ipm.evolve.flow(rhoNew,opsNew,state.scale)` 计算，
它消费的是新分解；记录和 checkpoint 仍走原来的维护路径。

这一步只适用于外部 checkpoint 换网格事务。它不直接适用于
`ipm.remesh.adapt` 的运行期重网格：后者失败时要原样回退并继续演化，
因此不能提前销毁回退态的求解缓存。外部事务失败则抛出异常，源 checkpoint
始终不变，重试可重新严格恢复它。

## 独立候选与验证

`ipm_perflab_regrid_checkpoint.m` 是现维护 regrid 函数的研究副本。
除入口改名、释放旧 `poisson` 及额外的第四个输出 `memoryInfo` 外，
transfer、flow、audit、record 与 checkpoint 原文保留。
`memoryInfo` 只记录对象 class、旧字段是否已移除和原 transfer metrics。
它不进入 checkpoint、audit 或配置。

`ipm_perflab_verify_regrid_memory` 使用 `65×33` 的小型 fixture，比较
PCHIP/四阶/六阶物理模式，以及四阶动态 `wall_omega_quadratic_peak` 规范，
每种各测原网格和同端点微变形网格。逐位比较：

- 整个 checkpoint payload、signature，以及除创建时间外的整个 checkpoint；
- 完整 audit、直接调用维护 `transfer` 得到的 rho 与全部 transfer metrics；
- 源 checkpoint 不变，以及新网格分解仍存在。

测试 fixture 的分辨率计数下限与事务阈值显式用于小网格等价检查，不是生产
网格/可信门建议。物理 PCHIP 本就没有守恒校正，因此仅它的测试质量门为
`1e-3`；初始 `33×17` 的变形 fixture 也触发了原版的可信前缀拒绝，最终
改为充分分辨这些测试特征的 `65×33`，没有更改生产函数或生产验收逻辑。

真实 q512 的峰值 RSS 和墙钟收益尚未测量；应在计算队列空出后串行比较，
避免多个 MATLAB 进程的交换/压缩压力污染数据。

2026-09-08 已执行：上述 8 例全部通过，包括所有 payload/signature/audit/
transfer 的逐位比较；旧、新 `poisson` 的运行期 class 均确认为 `decomposition`。
报告位于 `result/verification/performance_lab_20260908/regrid_memory_small.mat`。
候选 regrid 文件 Code Analyzer 为 0 条；验证脚本最初有 1 条 ASGLU（只为随后
`clear` 而接收 `transferOps`），已改用 `~` 接收，不涉及任何数值表达式。
