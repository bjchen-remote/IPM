# 在原t0冻结算法的分层参数搜索

2026-09-09。研究入口为`mesh_hierarchical_axis_search(view,referenceAxes,anchor,policy,maximumTotalNodes)`，配套纯枚举器`mesh_hierarchical_axis_schedule(search)`。两者不接受步号、时间、case标签或初始密度采样器；不修改维护`+ipm`或任何现有运行冻结源。它解决有限参数网格漏解，不解决物理域闭合或任意窄核心的无限容量。

## 从构造器约束生成范围，停止人工补整数列表

维护`roundedAxis`要求fine为整数且至少8，rounding为整数且至少2；左右外支各需至少`rounding+2`个cell。因此

`fine + 2*(rounding+2) <= (Nx-1)/2`

是可直接拒绝的必要节点条件。左式也衡量细平台与两个平滑过渡占用的最小cell预算；小值先试，可以给远尾留下更多节点，但不保证质量或core更好。实际anchor最近节点会离散切换，warp也可能推走核心细区，不能用粗点失败或分数插值证明整片参数区域不可行。

旧注册为fine32:16:128、rounding16:8:48与fraction[.5,.65]，共70项。仅在旧范围内取中点永远到不了rounding8–14，因此必须在**运行前**登记下界探索，而不是遇到晚期失败后临时指定36/40/44等数值。这里固定以下算法：

| 层 | 从原注册推导的整数集合 | 去重后额度 |
| --- | --- | --- |
| L0 | 原fine、rounding、fraction的枚举及排名逐值保留 | 70 |
| L1 | 两个下界各向下扩一个二倍尺度，受构造器8/2下界保护；上界不变；原参数间隔各二分。当前得到fine16:8:128、rounding8:4:48 | 新260，累计330 |
| L2 | 相同范围，间隔再各二分。当前fine16:4:128、rounding8:2:48 | 只取余下170，累计500 |

新增元组先去重，再按`[fine+2*(rounding+2), fine, rounding, fraction原索引]`稳定排序。L0排序完全不变。最多500个不同X元组/目标reference成员，包含构造前被必要预算拒绝的记录；不是每层500，更不是无限循环。Y保持原四个sigma、原equalizer和immutable reference，warp保持.25、x/y padding保持1.10/1.15。节点资源cap由caller显式提供，不能随搜索放大。本次两个895×386研究态的cap为345470，N和箱子始终不变。

lower-octave、二分深度2和500是明确的有限搜索设计选择，不是可行性定理。该设计受历史失败启发；历史覆盖检查不能伪称独立验证。合法原参数数组必须是规则整数间距，且能两次二分；不支持的注册形状明确拒绝，不能静默舍入、改变其原节点预算。

## 可检验的停止与资格规则

1. 先完整运行原planner。只要它返回任何qualified轴对（包括合法keep），直接返回其原候选与原报告；不触发后备或重新排名。
2. 若原Y族无合格轴，停止本X细分并保留Y失败原因。X参数搜索不能解决已知Y瓶颈；更高注册node-family是另一层调度，不能擅开新N。
3. 仅原轴对为空且Y可行，依次完整计算L1、然后必要时L2已登记预算前缀。某层有合格X后，与原合格Y配对，按原最弱无量纲质量裕量排名，稳定tie使用schedule索引和原Y索引；返回最多3对并停止，后续元组为**未计算**。
4. 全500仍无轴对时，返回`registered_search_budget_exhausted_not_global_infeasibility`。整个L2整数空间尚有未测点，更不涉及其它工厂，不能写“无可行网格”。原范围上界、warp或Y族若需要变化，须作为未来新版本预先登记。
5. 轴对仍须通过完整冻结field与真实native迁移。目标core32、前沿20、精确±anchor，以及相邻比1.08、曲率.01、rcond1e-9、全局q比1e-8、局部q/controlWidth [.35,1.65]全部保留。实际事务31、峰/质量/范围/累计门不变；本纯几何模块不作CP提交或flow。已返回轴对后field/native失败，不自动借此扩大搜索或绕开本次3对预算。

当前实现保留每个已试元组的原因、质量和覆盖报告，未试元组只在schedule中，不产生假评分。每次规划都从同一个真实paired source与immutable reference出发；不将本层候选串成下一候选的“新source”。主线程可在不同reference family之间应用相同有限算法；每级预算需独立，较小成员耗尽不能饿死较大成员，同时完整运行的总体尝试数与资源cap应预登记。

## 历史参数覆盖与真实态验证

原step9205的40项研究有10个X通过。固定新schedule在500以内包含全部10个元组：fine40/round16两项在105/106；fine44/round8在381/382；fine36/round14在387/388；fine44/round12在405/406；fine48/round14在437/438。**这些是已有报告的键覆盖，未声称本算法本次全部重算它们。** 如果L1已取得合格轴，L2应依约不运行。

规则已在真实几何运行前保存于
`result/longtime/20260908_campaign_v1/hierarchical_axis_search_v2_60e9cb6651e4/registration.json`。
同目录`historical_coverage.json`保存上述覆盖，`input_manifest.json`记录两份真实CP/result及原axis-plan。
实际运行注册两个案例，均已完成，MATLAB exit0：

- 9205：L0原70项及其空候选与原保存报告逐值相同。L1完整新增260项，累计330项，19.284秒找到合格轴并停止；L2未计算。首项为此前人工40项列表没有的fine16/round16/.5，说明实际运行没有仅重放已知成功键。
- 8057：4.090秒；原70项报告与原3个候选全部逐值相同，后备未调用、额外试轴记录为空。

[两个真实态报告](../../result/longtime/20260908_campaign_v1/hierarchical_axis_search_v2_60e9cb6651e4/report.json)保存原生签名、result/history配对、来源未改和profiler无LU/flow/PDE/transfer证据。113个冻源及全部输入SHA复验不变。

随后在独立目录固定保存的前三，不重新排名或搜索，用原完整冻结场scorer真实执行三次双向高阶转移，全部原门通过，合计1.162秒：

| 固定排名 | fine/round/fraction；原Y sigma | source特征在候选轴上的core X/Y；front | rhoX峰相对跳变 | 守恒相对缺陷 | 相对值域越界 |
| --- | --- | --- | --- | --- | --- |
| 1 | 16/16/.5；.5 | 33.9890/36.8000；28.9425 | 8.60836e-5 | 4.39852e-15 | 2.13787e-7 |
| 2 | 16/16/.65；.5 | 33.9946/36.8000；28.9729 | 8.97481e-5 | 4.24143e-15 | 2.13787e-7 |
| 3 | 32/8/.5；.5 | 35.1929/36.8000；33.6703 | 2.95189e-5 | 4.55561e-15 | 2.87639e-7 |

这些core数是原scorer的几何指标，**不是生产迁移后重测core**。三次均保存完整field/rhoX/rhoY/roundtrip观测；不能将表中的峰跳或质量门当作动态收敛证明。第一项与此前累计预算相加约.00058329665，仍须原生重算而非直接写ledger。

证据与待执行入口：

- [完整三对报告](../../result/longtime/20260908_campaign_v1/hierarchical_pairs_native_prepared_v1_6a44ec12d610/pairs/report.json)；[首候选](../../result/longtime/20260908_campaign_v1/hierarchical_pairs_native_prepared_v1_6a44ec12d610/pairs/candidate_01.mat)。
- [原生预飞](../../result/longtime/20260908_campaign_v1/hierarchical_pairs_native_prepared_v1_6a44ec12d610/native_preflight/report.json)已实际通过；其profiler无LU/flow/PDE/transfer，MATLAB session87694 exit0。
- [零时间原生迁移runner](../../result/longtime/20260908_campaign_v1/hierarchical_pairs_native_prepared_v1_6a44ec12d610/run_native_trial_01.m)仅准备，尚未执行。新kind与searchOrigin显式标记hierarchical来源/schedule77，保留原31/31实际core、front20、峰.002/累计.02及全部原生history/ref/clock/cursor/roundtrip门。不是原40族candidate。
- [字节审计](../../result/longtime/20260908_campaign_v1/hierarchical_pairs_native_prepared_v1_6a44ec12d610/final_byte_audit.json)：122个源及4个输入均未变；111个`+ipm`与既有原生试验冻源逐字相同，完整scorer也逐字相同。维护`+ipm`、主线冻源和旧判定未改。

首试`hierarchical_axis_search_v1_ad690f1961b2`仅在CodeAnalyzer前置检查失败，没有计算真实几何。第二目录只改`sum(logical)`的静态写法和有界cell追加注释；范围、顺序、预算、门与停止规则均未改变。最初进度消息把最后一个通过组合手工抄作round10，完整机器读取更正为round14，最晚索引438；枚举算法和500上限没有因此调整。

## 下一运行期的最小集成资格

已将上述数值算法移入[独立纯模块候选副本](../../result/verification/autonomous_runtime_20260909/hierarchical_planner_candidate_v1_41e91d9ebf20/README.md)，仅新增`+ipm/+remesh/hierarchicalAxisPairs.m`与`hierarchicalAxisSchedule.m`，119个原生产文件逐字保留。接口为`[candidates,report,searchAudit]=hierarchicalAxisPairs(view,referenceAxes,anchor,policy,maximumTotalNodes)`：原primary成功时前两输出完全不变，搜索审计单独放第三输出。后备报告将原70行作为逐值前缀附加试轴，候选全局索引可直接定位报告行。

该副本的[真实noLU回归](../../result/verification/autonomous_runtime_20260909/hierarchical_planner_candidate_v1_41e91d9ebf20/run/report.json)已exit0，53.075秒：8057原报告/3候选exact；9205与研究实现候选exact、330项停止；5坏输入拒绝；Y全拒不fallback；不同N的必要分支预算负例精确500项耗尽且无伪keep；二进制长度单位×2轴/选择exact；真实1134的321×161源到合法641×161参考成员仍保持原primary候选/report exact。两新文件CodeAnalyzer=0、依赖闭包全部`+ipm`、profiler无LU/flow/PDE/transfer。[SHA复验](../../result/verification/autonomous_runtime_20260909/hierarchical_planner_candidate_v1_41e91d9ebf20/final_byte_audit.json)123源+11输入未变。

纯模块没有修改维护`+ipm`，没有配置启用或controller/ledger接线；与节点族v3候选是两份独立协议，不能拼成组合版本已通过。原生迁移仍只准备未执行，完整包布局/原生持久化/split与从t0长时集成仍需新的显式资格。

上述纯回归已覆盖Y全拒、500上限、必要分支预算、资源cap与配对拒绝；后续仍需完整控制器的case/时间标签独立性、保存/恢复同源选择及失败后的安全交回资格。旧版本保持原70逐值行为，新能力只能由显式新策略版本启用。

若后续实际field/native三候选均失败，本模块仍可能安全交回；不能承诺所有未来网格问题都已自动解决。原t0到目标物理时间的无需人工改族实验及独立空间/时间/箱子误差资格仍是最终验收。固定计算箱随Cx增长缩成小物理箱、远场覆盖缺失与外尾边界模型是另一未解决能力，增加搜索次数或降低q门都不能替代它。
