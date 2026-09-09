# 与 fresh lineage 无关的配对网格设计器

入口：`[candidate,report]=mesh_design_paired_snapshot(snapshot,referenceAxes,actualAnchor,controls)`。
这是研究级、无 LU 的同箱同节点数设计；没有原生迁移、PDE、历史生成或 box promotion。

`snapshot` 包含 `rho/x/y/scale(Cx,Cy,Comega)/canonicalTime/physicalTime/trusted/provenance`；
`referenceAxes` 包含不可变真实 `x/y/provenance`。二者必须使用相同坐标单位、相同节点数和端点。
调用者负责真实原生资格；输入 `trusted` 是调用者的资格声明，函数不声称已验证签名。
原初始轴可以只有接近 anchor 的节点，但所有输出候选和可接受 keep 均严格含 `±actualAnchor`；
不会为了 snap 修改不可变参考轴。

控制范围在进入设计前固定：目标核心32/32、前沿20；fine=[32 48 64 80 96 112 128]、
rounding=[16 24 32 40 48]、core fraction=[.5 .65]，共70个x轴；y sigma=[.18 .25 .35 .5]。
x/y padding=1.10/1.15，warp radius=.25 anchor units。先完整枚举轴质量/区间覆盖，
再按最弱无量纲质量裕量稳定排序，至多3个完整冻结 pair transfer，首个全门通过即返回。
质量、峰跳、质量守恒、值域及局部q门全部固定，controls不能覆盖这些门。
改变搜索范围须在调用前明确 controls，并满足有限 axis/pair 预算；报告保存全部已尝试拒绝。

数值不准入返回有原因的 status；无效输入、缺依赖等程序错误抛出。`transactionReady`
仅表示具备进行实际事务试验的冻结资格。研究 scorer 和生产迁移不逐值等同，实际核心和
累计预算仍需在原生迁移后重新验收。

## 实际三态回归

冻结输出：[v5完整报告](../../result/longtime/20260909_paired_designer_v5/run/report.json)，
同目录各态 MAT 保存 candidate/design/snapshot/reference；137个冻结文件、精确runner与SHA均保留。
MATLAB exit0，两个新研究文件 Code Analyzer 无发现，全程无 LU/PDE。

| 输入 | x轴通过 | y轴通过 | 最终结果 |
|---|---:|---:|---|
| 原始t=0解析初值，H1e6，1025×513 | 50/70 | 4/4 | core34.70193/130.65508、front25.86644；fine32/round16/fraction.65，y sigma.18；峰跳1.06e-9 |
| 原大箱step3326，tau7.0932 | 64/70 | 0/4 | 明确`registered_y_family_infeasible`；y在global q≈1e-8处仅到core20.135～20.716，不能达到32；没有尝试PDE |
| 当前fresh step3668，等价tau10.0432，895×386 | 29/70 | 3/4 | core35.11770/36.80000、front34.60023；fine32/round24/fraction.5，y sigma.5；峰跳1.56723e-6 |

最新fresh候选 x/y 相邻比1.06667138/1.06518277，global q比1.16760e-6/4.96576e-8。
相邻比比旧64/.65族有更多余量，但y全局q仍有限；这不是任意长时间容量保证。
候选在 `run/fresh_step3668.mat` 的 `candidate` 变量。实际事务未执行，不能称动态通过。

原t=0重采样使用原native保存的physics和base轴，调用原`ipm.field.initialDensity`解析公式；
原root没有存t=0 snapshot，故不伪造t=0 checkpoint。
[额外初值审计](../../result/longtime/20260909_paired_designer_v5/initial_invariant_audit.json)
确认 physics、base轴、rhoRange0、mass0 均逐值一致。原x最近anchor为1.0000000000000002；
输入如实保留，输出候选精确±1。

五个负例（错配rho、不可信、变箱、覆盖硬门、超出注册搜索预算）全部拒绝。
最新fresh所选候选的坐标单位×.37和×4完整配对检查均通过；最大轴误差1.324e-13、
核心/前沿误差7.231e-11、全部无量纲迁移指标差1.733e-10，预设容差1e-8。
这些unit_factor文件是明确的合成单位变换fixture，不能被当作实际native候选来源。

保留的早期失败：v1只在Code Analyzer停止；v2的whichguard在设计前拒绝缺失两几何helper；
v3拒绝原初始轴非精确anchor的过强输入条件；v4停止于空struct数组报告追加。
v5修复的是输入/闭包/记录机制，所有输出数学硬门不变。

## 依赖与下一层边界

研究接口依赖 `mesh_core_patch_axis → mesh_general_rounded_axis`、
`ipm_gridlab_equalize_axis`、`ipm_gridlab_feature_profiles`、`ipm_gridlab_score_frozen_pair`
及其既有 `+ipm` 纯网格/诊断依赖。实际回归额外调用原native读取和fresh严格review，
这些读取不属于通用设计器的输入合同。

## 纯生产层迁移：已完成无LU回归

新增 `ipm.diagnostics.meshFeatureIntervals(acceptedView)` 与
`[candidates,report]=ipm.remesh.plannedAxisPairs(acceptedView,referenceAxes,anchor,policy)`。
前者要求真实 `rho/x/y/Dx/source/trusted`，强制 `source==rho*Dx'` 逐值；返回
`coreInterval/frontInterval/yCoreWidth/coreCenter/actualCoreCells/leftFrontCells`。
actualCoreCells 沿原 `peakResolution` 的 .9 X/Y 定义，竖直列仍为原 nodal 峰列。

后者接受已正规化 policy，不重新解析 config，只筛全轴质量与预测区间覆盖。
候选字段 `x/y/unchanged/xIndex/yIndex/predictedCells/quality`；按原稳定排序返回最多3个。
`predictedCells=[xCore,yCore,leftFront]`，不能作为真实迁移后的核心验收。
policy 接口为 `targetCoreCells`、`minimumFrontCells`、`qualityLimits` 和 `search`；
`search` 使用 `fineCells/roundingCells/coreFineCellFractions/ySigma/maximumAxisCandidates/`
`maximumPairCandidates/xPadding/yPadding/maximumWarpFraction`，含义与上述研究controls一致。

纯工厂迁入 `ipm.remesh.roundedAxis/corePatchAxis/equalizedAxis`。
没有迁入完整冻结pair scorer，也不调用 `flow`、`transfer`、checkpoint 或 research 文件。
[生产层实际回归](../../result/longtime/20260909_production_axis_planner_v2/run/report.json)
确认三实际态的几何区间、所有70x/4y轴评分逐值一致，t0/fresh3668的研究选中pair均为
生产候选1，full3326同样明确零候选。tiny解析例、两单位缩放、三负例全通过。
五个新增生产模块与测试 Code Analyzer 无发现，MATLAB exit0。

`requiredFilesAndProducts` 的闭包只有8个 `+ipm` 文件：上述5个，加
`diagnostics.peakResolution`、`mesh.quality`、`mesh.quadrature`，没有research依赖。
冻结副本、runner和全部SHA保存在 `20260909_production_axis_planner_v2`；v1只在测试脚本
格式检查停止，未执行设计。此次迁移保持数学运算，报告类型和数值拒绝的命名空间按生产模块命名。

在线实际transfer后的统一验收和调度由根任务集成进原求解器。此处仅完成独立模块迁移及
no-LU回归，尚不能替代从t=0无人干预的完整长时验收。
