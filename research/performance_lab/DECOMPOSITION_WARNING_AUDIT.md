# Poisson 重复条件警告：已缓存检查，未量化日志输出

2026-09-09，只读审计。按照本次任务“若原实现已缓存或已研究过则明确报告并止损”的条件，本轮没有启动 MATLAB、小矩阵、LU、flow 或 PDE，也没有修改生产和运行源。

## 直接结论

本机 MATLAB R2026a 的 `decomposition` **已经缓存条件估计，不在每次回代重估 rcond**。`research/performance_lab/FRESH_RHS_COST_AND_LU_SCREEN.md` 第25行也已记录同一结论。当前可疑成本只有重复条件分支、消息格式化、调用栈和日志输出，不能将它描述为重复昂贵条件数计算。

生产高阶 `mesh.build` 第175行使用 `decomposition(A,'lu')`，在同一 ops 生命周期重用 `ops.poisson`；`field.poisson` 的 isotropic kappa=1 路径直接使用该对象和原 Green 修正 RHS。各向异性 kappa 改变时会另建 factor，不能跨 A 当作同一审计对象。

本机可读主文档/实现路径：

- `/Applications/MATLAB_R2026a.app/toolbox/matlab/matfun/decomposition.m`：第126–132行公开 `CheckCondition` 可读写属性；第301/317–325行在第一次启用检查且缓存为空时生成 `rcondSaved`；第381行执行原 solve，第395–396行在解已算完后调用警告逻辑；第595–610行 `rcond(object)` 返回已有缓存；第697–706行警告逻辑只使用缓存。
- `/Applications/MATLAB_R2026a.app/toolbox/matlab/matfun/+matlab/+internal/+decomposition/SparseLU.m`：内部 LU 的 `solve` 调用原 UMFPACKWrapper。此内部类明确标为未公开接口；不提出调用/修改该内部对象、refinement 或因子格式。

公开属性从 true 改为 false 仅关闭后置警告检查，不改变源文件中已经执行的 solve 调用。源码层面可以推导同对象的数值求解路径不变；**本轮没有新增实测，不能把这个静态结论写成已经通过 bitwise 或速度验收**。

## 可复用的既有实际成本证据

`result/longtime/20260908_campaign_v1/fresh_tau744_baseline_rhs_cache_v1/profile/profile_report.json` 保存了默认10线程、真实895×386 step328原 flow 的排他 profiler：

| 项目 | 一次调用耗时 |
|---|---:|
| 完整 flow（含 profiler） | 639.054 ms |
| decomposition 原 solve 行381 | 465.561 ms |
| warnIfIllConditioned 总调用 | 0.07625 ms |

检查占该次完整 flow 约 **0.0119%**。该次没有执行第704行的 nearlySingular warning 发出分支，因此不含当前 v4 大量警告的实际 I/O。不能据此将 v4 日志输出成本断言为零，也不能据数十万行直接推导显著求解加速。

全生产 `+ipm` 搜索没有用 `lastwarn` 控制数值路径；仅 output 模块另有清单写入警告。不过全局 warning/lastwarn、调用栈、日志内容和“warning as error”的会话行为仍属于外部可观察差异，不能声称完全不改变运行行为。

## 若以后单独登记日志优化：最小可检验设计

这是未执行的研究测试设计，不是生产补丁或可用验收：

1. 每个真正新 factor，包括新网格、恢复时重建、失败事务回滚后的重建，保持默认 `CheckCondition=true`。以第一次实际 Poisson RHS 执行原回代，正常发出首次警告，并保存 A 的身份、节点/epoch、缓存 `rcond(factor)`、警告类别与消息、有限性状态。不得在首次检查前将构造参数设 false，也不得用全局 `warning off` 吞掉新 factor 的证据。
2. 复制同一分解对象，再只将副本的公开 `CheckCondition=false`。不重算 LU，不重设 permutation、scale 或 refinement。将副本属性复原为 true 后，以 `isequaln` 核对原对象，并核对 rcond、所有公开不变量和同一 A/RHS 逐位一致；必须实际测试对象复制/属性语义，不能只凭推测。
3. 首先登记小矩阵条件阶梯及小网格。保存每一对原/候选 backsolve、完整 `field.poisson` 返回的 psi/RHS/Green边界/operator/solveInfo，以及原 `flow` 的全部 RHS/尺度率/场；全部必须 bitwise。小残差不能代替前向误差或病态矩阵的物理解精度，也不能因 warning 变少宣称条件数改善。
4. 随机 ABBA 顺序分别量化默认警告输出与仅关闭副本重复检查，记录真实日志字节/行数和 wall time；首次警告及每个新 A 的警告另存。已有无警告的 profiler 数据不能替代这一量测。全局静音只可作为额外隔离实验且必须恢复，不作为生产候选。
5. 只有小例严格通过且重复日志成本显示有意义的收益，才考虑在已有 restore 自然窗口中借用一次真实 cache；不可为此额外建大 LU。本任务明确在识别缓存与既有记录后止损，因此这些步骤均尚未执行。

重复警告治理可能改善日志可读性和体积；当前证据不支持把它列为已经验证的主计算效率优化。当前主线仍保持所有原条件警告。
