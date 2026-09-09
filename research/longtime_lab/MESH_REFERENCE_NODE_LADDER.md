# 从零时间容量停止：只在 x 加点即可取得合格研究候选

2026-09-09；这是独立无 LU 几何与冻结转移实验，没有改变生产代码、配置、旧结果或运行中的冻结源，没有创建变量 N checkpoint 或推进 PDE。

来源是原 t=0 压力运行 [from_zero_tau4_v1](../../result/verification/autonomous_runtime_20260909/from_zero_tau4_v1/report.json) 的严格末态：step1134，canonical tau=2.6886192715883674，physical t=1.613029050260455，321×161，remeshCount=13，原生核心 26.00091/32.34267，累计峰跳 0.0007106094。该有限箱压力实验安全停于 `autonomous_mesh_axis_capacity`，并未达到其 tau4 目标。

[四层实际报告](../../result/verification/autonomous_runtime_20260909/node_ladder_from_zero_1134_v2/run/report.json)、[完整性审计](../../result/verification/autonomous_runtime_20260909/node_ladder_from_zero_1134_v2/final_audit.json) 已完成；MATLAB 91130 exit0。原始双轴扩展 v1 亦 exit0 并保留。两个实际输入及全部冻结源 SHA 未变。原生 CP/result 的 rho、轴、步数、三个时钟、配置、history、caseId、remeshCount、指数尺度与位移精确配对，整个原始 history trusted。

## 容量原因与最小梯度

注册次序按总节点数递增，最多评估原搜索给出的前三对轴；一对通过明确迁移门后该层停止。未为了通过而改变原 70 个 x / 4 个 y 搜索参数或门槛。

| N | 总节点 | 合格 x / y 轴 | 结果 |
|---|---:|---:|---|
| 321×161 | 51,681 | 0/70；2/4 | x 容量耗尽 |
| 321×321 | 103,041 | 0/70；4/4 | 只加 y 不能解决 x 瓶颈 |
| **641×161** | **103,201** | **70/70；2/4** | **最小合格研究候选** |
| 641×321 | 205,761 | 70/70；4/4 | 同样通过，但此时无需双轴加点 |

原 x 搜索有 30 个组合因 rounded branch 的节点预算而不能构造；其余 40 个全部违反 adjacent ratio≤1.08。其中 38/34 个还分别违反局部 quadrature 下/上限，22 个出现 global q 与非正 quadrature，8 个违反 log-spacing curvature。统计可以重叠。原 y 的两个拒绝仅为 core target；另外两个本就可用。因此这里的停止不是箱体裁剪或放宽网格门的理由，而是固定 x 节点预算下注册族无法同时承载核心和外侧过渡。

停机过程中失败的下一尝试状态未当成有效源保存。本实验从实际已接受 step1134 重新读取并计算，旧族在该可信源本身也确实没有合格 x；没有插值/回推未接受的下一步场。

## 最小候选的量化结果

候选文件：[level_3_candidate_1.mat](../../result/verification/autonomous_runtime_20260909/node_ladder_from_zero_1134_v2/run/level_3_candidate_1.mat)，含 `candidate/score/observation/rhoNew/snapshot`。它的 `nativeTransactionSupported=false`、`nativeTransactionPerformed=false`。

- 原 hard grid 门保持：adjacent ratio≤1.08、curvature≤.01、rcond≥1e-9、global q≥1e-8、local q/control width 在 [.35,1.65] 且权重为正。
- 实际 x/y adjacent ratio=1.0246510/1.0571510；global q=.1033325/.006057865。
- 冻结原特征区间上的 core/front=35.19231/36.80000/30.76132，均超过32/32/20。
- 用真实旧 321×161 rho 经两轴高阶保守插值后，重新计算派生 source 与原单峰诊断：core=35.16929/36.02969，front=30.57193，超过实际迁移 floor31/31/front20。
- 全域 rhoX 峰相对变化1.02760485e-4≤.002；质量相对缺陷1.44767e-15≤5e-12；值域相对违规2.26230e-7≤2e-4。加到原累计预算后仍小于.02。

独立的[保存后标量审计](../../result/verification/autonomous_runtime_20260909/node_ladder_from_zero_1134_v2/run/independent_saved_scalar_audit.json) 还核对了冻结32/32/20目标。该 Python HDF5读取本身不是原生签名验证；原生严格校验已在前面的 MATLAB 实验中执行。

往返场误差原样保留：field L2=4.82566e-7，rhoX L2=3.65422e-4、LInf=5.14710e-3（约0.515%），rhoY L2=3.73031e-5、LInf=7.68376e-4，x 对称缺陷5.47e-15。`score.admissible` 只表示轴准入，不能把它解释成所有误差目标或动态比较已通过。本实验另显式检查了登记的核心/前沿/峰/质量/值域/累计门；未把场误差观察自动升级为动态资格。

## reference family 与真实场的分离

[mesh_plan_reference_node_count](mesh_plan_reference_node_count.m) 的研究接口为：

```matlab
[candidates, report] = mesh_plan_reference_node_count( ...
    featureFromActualOldView, oldAxes, targetReference, actualAnchor, policy);
```

`featureFromActualOldView` 只由严格读取的原 rho/x/y 与同一原生高阶 Dx 得到。wrapper 复用原 `plannedAxisPairs` 的 core-patch、y equalizer、全轴质量、排序和最多三对候选；仅将构造的节点预算与源场大小分开。**没有创建更高 N 的 source rho，也没有把插值数据用于提取旧源特征。** 同 N 时，候选结构和全部 x/y 评分与维护版 `plannedAxisPairs` 逐值相同，已实际验证。

reference 的根是 `checkpoint.state.baseX/baseY`，并与原 initialization.selectedBaseX/Y 逐值核对。需要加点的一轴在其 normalized index 上用 PCHIP 求中间点，旧 knot 位置显式回填原值；x 只计算非负半轴再镜像。原节点、端点、0、±transportAnchor 均精确保留。只 x 加点时 y 的 reference 数值保持原样。

这里的嵌套性只属于 **reference family**。最终候选使用 core-patch/equalizer 重新分布，所以除端点/0/±anchor 外，不保证保留所有旧实际节点。“只 x 加点”也不意味着 y 不移动：y 仍可以在原161个节点上重分配，必须评分完整的 x/y 转移。

已测试入口：

```matlab
report = mesh_probe_reference_node_ladder( ...
    actualCheckpointFile, pairedResultFile, newOutputDirectory);
```

完整冻结运行脚本在 `node_ladder_from_zero_1134_v2/run_frozen.m`。输出目录必须为新目录；参考 family 四层及搜索门在运行前登记。

## 生产接入仍缺的合同

现维护版没有变量 N 事务合同。本候选不能交给原 same-N checkpoint/restore 路径，不能改 config.nx/ny 后伪装为原生续算，也不能让旧 history 的轴/尺度发生暗中替换。

[policy-v2 合同设计](../performance_lab/AUTONOMOUS_NODE_FAMILY_V2_DESIGN.md) 负责区分不可变初始 grid、reference family 层号和当前实际 N；后续至少还需真正变量 N 原生事务、完整 reference/gauge/clock/lineage 审计、同物理时间比较、并且保留原场与率误差门。本次离线 rhoNew 不是新 native 状态，未求出新格 Poisson/flow/rates，也没有稳定性结论。

641×161约为当前895×386主线的30%节点，但 sparse LU峰值不能只按节点数线性估计。当前机器16 GiB RAM；本实验无 LU，没有据节点数启动任何大型算例。生产候选仍须内存准入和串行因子生命周期验证。

这证明增加 x 的节点预算可以解除此特定可信末态的注册几何瓶颈，并提供通过明确离线迁移门的最小候选；不证明任意后期或任意盒边界都能靠同一层级成功，也不证明已经完成从零到目标长时间的连续运行。
