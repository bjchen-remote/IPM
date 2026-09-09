# 分层搜索与方向节点族的显式组合合同（提案）

2026-09-09。本文仅设计下一候选版本；不修改维护或冻结源码，不把两份独立资格合成生产资格。配置/version/controller由root与efficiency负责，纯几何模块已单独锁定。

## 已审对象与确定缺口

[方向节点族v3候选](../../result/verification/autonomous_runtime_20260909/node_family_v3_directional_candidate_v1/README.md)有8成员：`[1,1;2,1;1,2;2,2;3,1;1,3;3,2;2,3]`。实际初始321×161对应最大308481点，全部低于310000；按乘积/注册索引排序为`[1,3,2,6,5,4,8,7]`。16条最大分轴单调路径最长恰3次增长，第三次已到[3,2]或[2,3]。未见新增方向/资源放行漏洞；完整family仍从真正初选base重建，factor3只保证root节点精确保留，不保证factor2中间节点嵌套。原转移通路允许不嵌套轴。

[独立只读审计](../../result/verification/autonomous_runtime_20260909/directional_v3_independent_review_v1_af3e6647462b/README.md)确认一个旧v2共享的持久化边界：validator只限制总`candidateIndex <= eligible成员数*3`，attemptSummary没有逐次目标成员和成员内排名。真实noLU合成反例将sourceLevel1→targetLevel2的接受编号改为24（该成员排第3，实际最大9），仍被接受；25的总上界控制被拒。[完整反例](../../result/verification/autonomous_runtime_20260909/directional_v3_independent_review_v1_af3e6647462b/probe.json)保留，不是PDE或原生checkpoint证据。运行中的原planner确实每成员最多3并按顺序枚举；此缺口是持久化不能证明这条性质，并非观察到实际违规。efficiency已另建修复候选的任务，旧v1冻结不改；修复结果应独立验证。

后补：新directional_candidate_v2的[3正13反](../../result/verification/autonomous_runtime_20260909/node_family_v3_directional_candidate_v2/run_v1/attempt_order_report.json)已实际通过，包含直接载入上述反例及补齐字段后仍超预算的反例；[三文件独立只读复核](../../result/verification/autonomous_runtime_20260909/directional_v3_independent_review_v1_af3e6647462b/v2_remedy_readonly_review.json)未见接受账本目标/级内连续编号/成员顺序的新问题，v1/v2字段不受影响。修复仍不等于原生资格，完整计划与失败记录的证据边界见下文。

当前[纯分层模块](../../result/verification/autonomous_runtime_20260909/hierarchical_planner_candidate_v1_41e91d9ebf20/README.md)只新增两个`+ipm/+remesh`函数，无research数值依赖。8057原70报告/3候选逐值不变；9205在330试轴处取得与研究实现相同的3对；500预算/Y全拒/单位协变/合法不同N均有实际noLU资格。其首候选原生迁移目前仍未执行。

## 最小显式接口

建议组合能力使用**新策略版本**（下称拟议version4，最终编号由root决定），不把旧v3悄悄改成500项搜索。v1/v2的配置规范值、planner分支、报告形状、metadata和checkpoint读取合同保持逐值原样；已启动的v3也不能按新版本重标。

version4沿方向v3的八成员/generator/cap310000/max3growth和全部质量、物理、CFL、峰/质量/范围/累计门。保留`policy.search`的原70×4/每成员3配置，另外冻结一个显式`axisSearchPolicy`描述符：

```matlab
struct('version',1, ...
  'algorithm','primary_then_lower_octave_dyadic_v1', ...
  'activation','evolved_remesh_only', ...
  'lowerExtensionOctaves',1,'refinementDepth',2, ...
  'maximumAxisTrialsPerMember',500,'maximumPairsPerMember',3, ...
  'stageCounts',[70,260,170], ...
  'stopRule','first_completed_stage_with_qualified_pair', ...
  'rankingRule','quality_margin_then_schedule_index_then_y_index')
```

初版建议只在**演化重网格**调用分层模块。t=0仍按既有原70初选、实际解析重采样门和显式一次固定观察fallback；两条规则均在t=0配置时冻结，整个运行无需人工增加列表。这减少初始观察双阶段账本的额外变化。若希望分层也用于初始原格/固定观察格，应另行显式登记作用域，重新验证实际解析core、phase索引及一次fallback资格，不能顺带启用。

纯接口不变：

```matlab
[pairs, detail, searchAudit] = ipm.remesh.hierarchicalAxisPairs( ...
    acceptedView, immutableMemberAxes, actualAnchor, policy, ...
    policy.nodeFamily.maximumTotalNodes);
```

`plannedAxisPairs`仍负责原primary；非空时前两输出直接返回其原对象，searchAudit只放第三输出。family planner用原eligible次序逐成员调用；每成员单独拥有500轴/最多3对额度，最多8成员、4000X与24个返回对，不是全局只留前三。不同N无keep；演化过滤qualified keep后保留原成员内排名，首个真实尝试可能从localPairIndex2开始。Y全拒不做额外X搜索；某阶段成功后不运行下一阶段。原生三项失败不触发第二轮隐性搜索。

候选核心预测始终来自同一真实source的配对场，reference只提供真实单位/箱/节点预算；不改箱、不把候选转成新的source、不中途采样原始初值。cap不是LU内存模型；额外4000项纯几何最坏规划成本也须记录，不能只报告首次早停开销。

## 最少持久化内容

方向v3修复中的每次尝试`targetLevelId/localPairIndex`应成为基础。组合版在此之上明确以下层次：

| 位置 | 新增或保留内容 | 验证要求 |
| --- | --- | --- |
| policy及初始化 | 完整规范化axisSearchPolicy、原参考family/generator/资源mask | 新字段只属于新版本；旧配置/CP不注入默认、不重标升级；恢复绝不采样初值 |
| controllerDecision | 原source step/clocks/remesh/level/nodeCount；有序eligible成员；每成员搜索摘要 | eligible按真实family重算，顺序精确；每成员70/330/500停止量与stage规则一致；保留primary是否非空、Y是否可行、fallbackUsed与未计算阶段 |
| 计划输出的紧凑描述符 | 全部返回对（最多24）各自的targetLevelId、localPairIndex、unchanged、searchPhase、refinementLevel、xScheduleIndex、yTrialIndex | 保留keep过滤前排名；原primary成功不得出现refined候选；X索引精确映射固定schedule，Y对应原四sigma；每成员最多3 |
| attemptSummary及接受audit | 与计划描述符相同的身份字段；原passed/错误/拒绝原因；原实际field/mesh/peak/mass/range/累计审计 | 尝试序列是过滤keep后的计划前缀，不能重复、逆序、跳过已登记提案或把最终接受项归给别的成员；接受项目标与实际新轴/family必须一致 |
| 失败记录 | 完整当前规划报告、真实尝试前缀、确切stopReason | 几何为空、field/native拒绝、预算耗尽分别记录；安全回到真正已接受CP，不扩大搜索或修改原判定 |

不要把巨大trial报告放进每个时间步history。每次规划保留完整报告作为独立artifact；checkpoint控制器保存上述有界描述符/摘要及其来源绑定。只存`passed`、总候选数或外部文件路径不足以验证规则。实际历史source未保存时，validator能验证账本自洽及身份/边界，**不能重新证明过去的primary为空、质量计算或场迁移正确**；需要源paired snapshot/CP与同冻源重播的事件资格来补足。文件哈希也不等于重新计算数值证据。

## 具体接线边界与最小验证

1. 配置正规化只在新版本接受固定axisSearchPolicy；family/growth/cap值来自明确组合注册。将现有变量N检查扩至新版本，但不改变旧分支。`initialMeshObservationPolicy`仅扩新版本兼容检查，初始数值算法维持其原作用域。
2. `planAutonomousMesh`在新版本的演化family循环调用纯模块，保存第三输出摘要和有序候选描述符。`applyAutonomousMesh`将身份带入attemptSummary及接受audit；原transfer、flow、审计公式不需要改。controller写入/读取验证器按上述身份和前缀规则核对。
3. 原v1/v2真实CP/result全文、默认四步和split保持逐值；方向v3修复版的独立资格先完成。新组合须重跑noLU：8057 primary不变、9205分层成功、Y全拒、500耗尽、成员独立预算、不同N无keep、cap拒绝、单位协变及参数索引映射；新增每成员预算/逆序/重复/身份错配/伪fallback/旧CP重标反例。
4. 原生测试另排队：至少一次真正分层选出的同N事务与一次方向factor3事务；随后相同完整终点的whole/prefix/resume，完整数学状态/history/ledger/refs/cursor配对。最后以新版本从原t0自然运行跨初始观察、普通重网格、增长与分层后备，无人工改列表。两份独立noLU报告不能替代这些组合运行。

物理箱随尺度缩小、远场覆盖消失与边界闭合误差仍独立于网格搜索；更多节点或500项候选不能消除这些科学限制。有限参数族耗尽不是“没有任何可行网格”，有限八成员也不是无限自适应。
