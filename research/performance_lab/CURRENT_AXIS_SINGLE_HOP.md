# 当前轴单跳、精确 X=1 的冻结网格研究

本研究只读取真实已接受场，沿用维护 equalizer 和 frozen-pair score。未调用 Poisson/LU、时间推进或 checkpoint restore，没有更改生产 C、solver、driver、root 或 runtimeRefs。所有候选仍须独立换网格事务与短程 PDE 验证；这里的 `transactionReady` 仅表示冻结门通过。

## Stage006：tau=6.0932

入口 `ipm_perflab_screen_current_axis.m` 的五个路径参数依次是：当前 `result.mat`、与之配对的 `mesh_design.mat`、原 root 的 `mesh_root_reference.mat`、最后一次已提交的 globally rooted `moving_platform.mat`、新的输出目录。可选第六参数为 `struct('sigmas',[.35,.5],'xConstructionPaddings',1.10)`；默认筛九档 sigma、两档 padding。

来源为 `adaptive_campaign_v3/stage_006/result.mat`（step2800、remeshCount7），外部参考仍是 unregridded H1e6 root。代码逐值验证最后已提交 platform 的 candidateX/Y 等于实际 current x/y、该 platform 的 reference axes 等于原 root；并检查它由已注册 root factory 精确构造、冻结门通过。因此 **当前增量引用深度是0，候选是1**，与历史 remeshCount7 不混淆。每次 sigma/padding 都从同一份 entry current x 构造，绝不把本轮前一候选当 reference。

现有 gridlab fallback 的九档均满足原 mesh targets，但锚点是 core/front 严格交集中心 `0.9700575505527312`，不是运行规范的 `X=1`。这些原候选到1的最近节点距离为 `3.98e-5` 至 `2.25e-4`，均无精确1节点，故只保留原设计/provenance，不直接作为规范候选。

新增独立构造使用同一维护 `ipm_gridlab_equalize_axis`：reference=current entry x、focusCenters=1、anchorPhysicalPosition=1。将原来精确±1节点对应索引 snap 回±1后，用原完整 pair score 重新检查整个 x/y 轴和转移。snap 位移最大约 `1.05e-7`，源于现 equalizer 在 H1e6 域尺度上的锚点容差；所有结果均在 snap **之后** 评分。y 仍从原 root y 构造，padding=1.15，不走当前 y 的增量引用。

硬门保持：core21/21、front20，ratio≤1.08、log-curvature≤.01、stencil rcond≥1e-9、全局 quadrature ratio≥1e-8、局部 quadrature/control-width 在[.35,1.65]，峰跳≤2e-3、质量缺陷≤5e-12、range≤2e-4。分辨率舍入容差1e-6与维护 `mesh_platform_candidate` 一致；本次全部计数也满足不减容差的21/20。

| sigma | x padding1：front / ratio | x padding1.10：front / ratio | padding1.10 峰跳 | 全部门 |
|---:|---:|---:|---:|:---:|
| 0.18 | 24.65687 / 1.06787973 | 25.41293 / 1.06969089 | 6.445e-05 | 通过 |
| 0.25 | 24.72598 / 1.06764772 | 25.97480 / 1.06853851 | 1.170e-04 | 通过 |
| 0.35 | 24.75309 / 1.06750677 | 26.15848 / 1.06817105 | 5.421e-05 | 通过 |
| 0.5 | 24.75985 / 1.06740691 | 26.19637 / 1.06825983 | 6.712e-05 | 通过 |
| 0.75 | 24.75956 / 1.06734229 | 26.18745 / 1.06856247 | 7.713e-06 | 通过 |
| 1 | 24.75862 / 1.06744248 | 26.17817 / 1.06898205 | 3.430e-05 | 通过 |
| 1.5 | 24.75766 / 1.06752808 | 26.16940 / 1.06934138 | 7.423e-05 | 通过 |
| 2 | 24.75727 / 1.06756113 | 26.16588 / 1.06948034 | 6.421e-05 | 通过 |
| 3 | 24.75696 / 1.06758579 | 26.16315 / 1.06958410 | 5.178e-05 | 通过 |

x padding1 的 coreX=21；padding1.10 的 coreX=23.1；所有 coreY=24.15。18/18 候选通过。可先保留 sigma=.35、padding1.10（ratio1.06817105，峰跳5.42e-5）及 sigma=.5、padding1.10 两项。没有根据此处冻结指标推断未来 PDE、时间步收益或无限次增量重网格可行性。

同一实际场的深度负例显式设 currentDepth=maximumDepth=1，维护回退返回 `maximum_incremental_reference_depth_reached`、attempted=false、usedAsReference=false、candidateIndices=[]，未绕过引用深度门。

结果位于 `result/verification/performance_lab_20260908/current_axis_stage006_v1/`，含两个完整维护 design、18个候选 MAT、`current_axis_report.mat/json`、原生 fallback 审计和 profiler no-PDE/no-Poisson 审计。初次包装器曾误读字段名并退出；原设计已保存，修正为 `reference.currentAxisXFallback` 后严格核对 frozen options 再复用。恢复运行 exit0、checkcode0，未覆盖原失败日志。

## 集成边界

当前长期 driver 尚未接入本研究。原 regrid metadata 的 `gridLabRegrids` 没有 incremental-depth 字段；不能在下一 stage 无证据地把深度重新设0。将来使用单跳候选时，须把 depth1、entry-axis signature、globally rooted parent、原 root reference 和实际事务结果持久绑定。下一次若仍基于这个 incremental child，则现有 maximumDepth1 应拒绝；只有有证据证明当前轴再次来自 immutable root 的直接全局构造，才能重新认定 depth0。研究入口要求 globally rooted parent 逐值匹配实际轴，不能把 current 冒充 root。

## Stage009：tau=6.693201945463523 的拒绝与机制核对

最新可信来源为 `adaptive_campaign_v3/stage_009/result.mat`（step3096，remeshCount8），实际轴逐值匹配 `stage_008/platform_design/moving_platform.mat`，仍是直接 root factory 生成，故增量 depth0→1 有证据。此时 positive X=1 索引为190；root 最近1的正轴节点索引为104，其位置距1为2.220e-16（原root并非逐位精确1，保持原轴不修改）。未覆盖 stage006 结果。

按实际末态只复测 sigma=.35/.5、x padding1.10。原 focus=anchor=1 两项均触及 ratio1.08、`mesh_gate_before_target`，core分别6.69112/7.16772、front5.78909/6.10003，0/2 ready。其 mesh、局部quadrature及transfer门通过，但分辨率门明确拒绝。`current_axis_stage009_v1/` 保存两份候选、原维护 fallback、depth负例及完整无LU审计。

维护 equalizer 的 Gaussian kernel 分别使用 `focusCenters` 和 `anchorPhysicalPosition`，允许把 focus 固定在真实 core/front 交集中心0.9717792027776753，同时 anchor仍设1。此设置只更改独立构造参数，没有改 equalizer 实现。两个分离设置的实际 frozen scores如下：

| sigma | 选中A | coreX / coreY | front | ratio | 结果 |
|---:|---:|---:|---:|---:|:---:|
| 0.35 | 2.00640267 | 11.113883 / 24.150000 | 8.399540 | 1.08000000 | 分辨率拒绝 |
| 0.5 | 3.38829310 | 8.584853 / 24.150000 | 6.955785 | 1.08000000 | 分辨率拒绝 |

额外仅做网格轴的 A=0/.01/.05/.1 机制检查（无场演化），原focus1和交集focus都测试两档sigma。A=0轴相对原轴的节点∞差仅1.421e-14，精确anchor仍1，core17.8823064、front21.3041846，说明6.7并非初始化损失。低幅度先增加core；随后计数随A非单调。维护 `resolve_bracket` 在无target成功时选择 **最大幅度的admitted trial**，其约定不是最大分辨率。故失败输出可能比entry轴更差，不能当成“最佳可行分辨率”。这也意味着失败返回的 candidate 绝不能越过当前resolution门进入事务。

例如focus=交集、sigma=.35的保存路径：A0 → core17.8823/front21.3042；A.25 →19.392/21.460；A.5 →20.6229/20.256；A1 →19.966/15.518；A2 →11.158/8.425；最终A2.0064才触及ratio门。已保存的所有trial均未达到padding目标23.1/22；没有穷尽连续幅度参数，更没有证明所有当前轴family数学上不可行。分离focus的两项仍不能作为主线后备。

`current_axis_focus_stage009_v1/` 保存16条0/低幅轴记录、两个完整candidate、原维护fallback与深度负例、无Poisson/PDE profiler审计、执行recipe。`current_axis_saved_audit.mat/json` 对三组已保存结果重新逐值核对源场/轴及原root signatures，保存全幅度trace和质量/转移统计。最终checkcode0，各进程exit0。唯一启动失败是stage006初次包装器字段读取错误，原失败日志保留；后续数值候选拒绝均按原门正常完成。
