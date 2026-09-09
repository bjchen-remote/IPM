# 自动长时网格：故障图谱与从零时间验收设计

审阅日期：2026-09-09。本文依据实际代码、保留的失败产物和已完成对照整理；本次只读审阅，没有新建椭圆算子或推进 PDE，没有修改当前或后续主线的冻结源码。

最终工程验收应是：给定原始解析初值、预先登记的箱子、误差预算和长时终点，仅调用 `ipm.solve(opts)`，程序从物理 `t=0` 自行完成初始网格、预测、候选搜索、事务验收及恢复，达到预定终点时没有人工修改网格参数、换来源或重写验收门。现有晚态自适应运行是这项工作的证据和组件，尚未达到这个验收。

有限节点数、固定箱子及正的全局求积下限不能保证任意长时间的任意窄核心都可表示。自动处理候选失败与诚实保存可信终态是必要能力；仅能安全退出仍不算“已完成从零到长时无人干预”。工程完成与奇异性、极限 PROFILE 的科学认定也必须分别记录。

状态更新：本文早期“当前代码/从 t=0 接口/CI”段保留的是首次审阅时的阶段记录。
现已实现并实测显式 v2 自动节点族；本页末尾的“第二方向自然增长”补充了真实从零运行证据，
不能再把早期“尚未实现”理解为最新版代码状态。有限注册族的运行通过仍不替代域、空间、时间误差资格。

## 当前自动网格设计的验收清单（2026-09-09更新）

| 环节 | 必须自动完成的行为 | 现有证据与缺口 |
|---|---|---|
| 原始零时刻 | 从解析初值识别真实核心；初选后重新解析测量，粗测失准时按预登记规则重测 | H32/H64/H128预选轴四步通过；自动fallback H64/H32四步、默认回归与split均通过，H16原初选成功且未使用fallback。显式可选规则已由root合入维护源码，H64从原t0长时运行已启动；未宣称长时科学资格 |
| 动态分布与增点 | 每个接受态评估两轴；先同级调整，再试资源允许的方向增点；始终保留准确anchor | v2从0自然1→2→4增长已发生；最高级后只能同级移动，有限族不能保证任意长时 |
| 候选失败 | 先无LU淘汰，再实际迁移/flow验收；失败保持来源及累计预算，不临时降门 | 原生拒绝、回退与跨增长恢复有实测；历史未保存的瞬时场不能补造 |
| 资源与调度 | 一个候选一个LU；节点/内存/步数预算明确，检查点安全落盘 | 主线+一个addon受控；节点cap不是RAM模型，4096分段预算仅研究调度 |
| 持久化 | 每帧rho与各自轴配对，级别/时钟/来源/触发/误差账本可严格核验 | 两次自然增长的实际原生split均逐位通过；空快照与真实跨尺寸快照分列；合法空迁移账本也须覆盖 |
| 科学可信度 | 同物理时刻核验场、梯度、率、箱子、空间与时间敏感性，保留旧判定 | 原B cap110k容量停止未到终点；新cap210k A/B均从零到共同物理终点，但两窗科学门仍false，CDE当次比较仍pending。固定H8在τ8的登记物理远场覆盖为0；不能以平滑Profile或小残差代替域误差控制 |
| 搜索分辨率 | 预算内既找细化程度，也找满足光滑/尾长约束的过渡长度，避免粗参数间隔造成假容量失败 | 裁箱step9205原70个X全拒；另注册40项中10个X通过，首项原生事务已过。新预登记分层算法在累计330项时自动成功，前三完整冻结场门通过；8057原70逐值保留，后备未调用。分层首项原生迁移尚未执行 |

每条失败记录都应说明触发、实际原因、修改规则、对应回归及尚未验证的范围。
安全停止是必要保护；只有目标长时区间无需人工改网格/源码/参数而完成，才满足用户的最终工程要求。

## 已登记的共同门槛

以下为本轮长时研究实际采用的门槛，不应被核心模块较宽的默认值替代：两轴相邻单元比 `<=1.08`、log spacing 曲率 `<=.01`、七点局部 stencil rcond `>=1e-9`、求积权重全正、`min(w)/mean(w)>=1e-8`、局部 `w/controlWidth` 在 `[.35,1.65]`。实际迁移峰跳 `<=.002`，累计绝对峰跳 `<=.02`，质量相对缺陷 `<=5e-12`，新增相对值域越界 `<=2e-4`，前沿至少 20 个单元。

当前 fresh 895×386 主线以核心 32/32 为设计目标，实际事务至少 31/31，每段末至少 20/20；完整历史最低 14/14、安全指标 `<.70`。这些层次不能混用：14 是可信历史底线，不是新设计目标；设计预测 32 也不等于实际迁移后达到 31。单次空间对照有各自更高的已冻结门（42 支事务 40、末态 34；56 支事务 54、末态 44）。

## 故障图谱

### F01：固定细区位置、宽度与晚时核心不匹配

- **触发**：当前核心离开原细平台；原参考轴 equalizer 继续缩窄 sigma 仍找不到合格轴；或初值型候选把 `H` 的固定比例当最小宽度。
- **原因**：固定形状/节点分配并未随核心漂移及收缩改变。`minimumWidthFraction=1e-6` 在 `H=1e6` 给出宽度 1，初始审计的严格核心却仅约 `.0472/.0139`。局部加密被外侧平台和参考节点编号限制。
- **证据**：[最初审计](mesh_findings_20260908.md)；[tau4.683 的六项更窄族](../../result/longtime/20260908_campaign_v1/mesh_narrow_stage003_v1/mesh_search_report.json)（sigma .35/.25/.18 × target21/24 全为 `infeasible_x`）；[平台节点预算与真实核心漂移](../../result/verification/performance_lab_20260908/platform_budget_stage006_v1/platform_screen_report.json)。
- **自动规则**：每个可信接受态重新测量核心、前沿和平台边缘的距离；以实际核心宽度/局部单元尺度给最小宽度；候选中心跟随核心。保留原始参考轴，但允许 anchor 所占节点编号改变。
- **回归**：解析峰在固定背景上跨越平台边缘并逐渐收缩，包含同一峰的多个节点相位；强制旧族失败后新族自动接管，检查所有轴门及对称性。
- **限制**：六项失败仅是该有限族失败，不是所有网格无解；当前新族仍有有限节点预算及正核心中心假设。

### F02：把 transport anchor 与核心强绑在同一常宽平台

- **触发**：`|peak-anchor|/coreWidth` 变大；为覆盖精确 anchor 而将常宽细区延长，挤占连接两侧粗网格所需节点。
- **原因**：规范采样点与最大梯度核心承担不同作用，二者可以分离。把 anchor 固定在细平台内部，或者把其原节点编号也固定，增加了不必要的几何约束。
- **证据**：[最新 epoch 的常宽族](../../result/longtime/20260908_campaign_v1/fresh_epoch3326_anchor_family_v1/scan_report.json) target24/32 × fine96/128/160/192/224/256 × fraction .90/.95 共 24 项全部失败；[直接核心补丁与推荐候选逐值等价](../../result/longtime/20260908_campaign_v1/direct_core_patch_equivalence_v1/report.json)；[实际核心补丁事务及对照](../../result/longtime/20260908_campaign_v1/fresh_epoch3326_core_patch_smoke_v1/smoke_report.json)。
- **自动规则**：细区围绕真实核心；通过紧支撑 C2 位移把某个过渡区节点放到精确 `±anchor`，节点编号自由。位移后必须回到真实单位重新检查整轴、成对迁移和前沿，不能只检查 bump 的解析光滑性。
- **回归**：anchor 在细区内、边缘、过渡区三类；anchor=.152142… 与 1；单位缩放 .25/.37/4；895 与其他奇数 x 节点数；过大位移导致非单调的候选必须拒绝。
- **限制**：`mesh_core_patch_axis` 假设正 anchor、正核心中心、对称 x、单一正峰以及可容纳的紧支撑 warp；不适用于核心穿原点、多峰竞争或任意规范，不能无条件泛化。

### F03：全局求积下限与远尾空间跨度耗尽节点预算

- **触发**：增加核心目标或减小核心尺度后，y 轴先撞 `min(w)/mean(w)`，或外尾连接处撞相邻比/曲率门；继续调整 x 也不能恢复完整配对。
- **原因**：`ipm.mesh.quality` 的全局比使用整轴平均权重。固定域长 L、节点 N 时，约有 `min(w)>=1e-8*L/N`，它不是局部稳定性的同义词。大箱的远尾使平均权重很大；坐标归一化同时缩放分子分母，不会消除该限制。
- **证据**：[step3326 的原根/当前 y 有限恢复搜索](../../result/verification/performance_lab_20260908/vertical_recovery_v4_stage002_v1/vertical_report.json) 没有 ready pair；[target42/56 冻结族](../../result/longtime/20260908_campaign_v1/fresh_tau724_finer_family_v3/screen/report.json)；[56 实际事务](../../result/longtime/20260908_campaign_v1/fresh_tau724_target56_long_v1/run/transaction_audit.json) 的 y ratio `1.078944679`、global q ratio `4.57611396e-8`，已接近两条边界。
- **自动规则**：预飞同时登记全局、局部 q 裕量和外尾单元预算，提前预测目标可行性；先改变节点分配及补丁形状。箱子缩小或增节点只能作为另行验证的契约能力，不能用降低原门来伪装恢复。
- **回归**：同一局部核心、多个 H 的冻结轴梯度；同单位缩放应不改变 q 比；让 y 失败而 x 通过的配对必须整体拒绝。记录具体质量门与失败轴，不只输出 `infeasible`。
- **限制**：早期研究事务只支持同N/同箱；后续显式v2已实现并实际验证注册节点族增长，但仍保持计算箱固定，不能任意增加第五成员或超资源cap。有限箱晚态成功不证明从 t=0 的小箱误差足够小。

### F04：预测核心通过，原生迁移核心却失败

- **触发**：冻结场评分达到设计目标，生产重映射后某轴低于实际事务门。
- **原因**：冻结评分器的迁移与生产高阶插值、光滑守恒校正并非逐值相同；核心和竖直取样列也随重构改变。只用预测核心验收，会提交实际欠分辨状态。
- **证据**：[tau4.883 原始失败日志](../../result/longtime/20260908_campaign_v1/mesh_platform_tau488_smoke_v1.log)：预测约 23.1/21，原生实际 23.0883/19.453，低于 20/20，未提交 CP；[之后独立注册并通过的 v2](../../result/longtime/20260908_campaign_v1/mesh_platform_tau488_smoke_v2/mesh_smoke_report.json)。
- **自动规则**：冻结设计用于排序与无 LU 淘汰；实际迁移再计算核心/前沿/峰/质量/值域以及全部局部 q 门，全部通过才写原生 CP。可预先登记设计余量（当前 x1.10/y1.15），不能在失败后修改实际门。
- **回归**：保存本次冻结-pass/native-fail 场景；tiny 原生事务覆盖生产校正显著改变窄核心的场；失败必须保持源 rho/axes/scales/history 逐值不变，无 remeshCount 增加。
- **限制**：固定 15% y 余量不是通用误差界；后续须用实际迁移统计和稳健几何决定裕量，但硬门独立。

### F05：累计迁移误差与重网格过密

- **触发**：每段无条件重网格；每次峰跳均小，却在多次迁移后累积；使用已变形当前轴无限递归当新根，参考来源逐步丢失。
- **原因**：单次守恒只约束一个积分量；小峰跳也不控制整个场、尾部或导数。多次插值可引入耗散、相位漂移和高阶导数误差。
- **证据**：[当前 fresh driver](run_fresh_profile_campaign.m) 的 `.002/.02` 实际与累计峰跳门；[跨重网格真实链](../../result/longtime/20260908_campaign_v1/fresh_qualification_chain_nolu_v1/positive_v2_stage1_final_v2/chain_report.json)；[有限 current-reference 搜索的来源深度约束](../../result/verification/performance_lab_20260908/current_axis_nonmonotone_v4_stage002_v1/search/search_report.json)。
- **自动规则**：有足够核心和预测余量时保持同网格；重网格后设滞回/最短收益区间。累计预算必须从实际逐事务峰跳重算，保持原始轴/初始参考及每次来源，不能信一个 JSON 总计。候选所有参数从同一个可信源出发，禁止本轮候选互相递归派生。
- **回归**：重复 A→B→A 的低资源事务循环；逐次检查质量/范围/峰及局部场、累计预算耗尽拒绝、参考深度和事务链缺失拒绝。
- **限制**：累计绝对峰跳 .02 不是全场累计误差证明；还需周期性同物理时刻的空间/时间/箱子敏感性对照。

### F06：nodal 竖直宽度的假跳与观察器阶数误认

- **触发**：x 网格变化使最近峰列切换；原离散 y 宽超过 .005，而同连续峰处的 y 宽差远小于门槛。
- **原因**：旧宽度在不同 x 节点列采样，混入横向偏移效应。连续峰位置与横纵耦合也会使加密误差比不均匀，不能根据局部插值阶数宣称几何通用四阶。
- **证据**：[step19→核心补丁辅助判定](../../result/longtime/20260908_campaign_v1/fresh_epoch3326_core_patch_smoke_v1/auxiliary_verdict_v1.json)；[42 终点宽度审计](../../result/longtime/20260908_campaign_v1/fresh_tau724_target42_long_v1/run/continuous_width_audit.json)：旧 y 宽差约 2.1216%，同连续峰仅 `1.2771e-5`；[4 相位/3 层 MMS](../../result/longtime/20260908_campaign_v1/continuous_width_observer_validation_v2/report.json)：最细误差 `4.25e-6`、协变 `8.88e-15`，但每次误差降至少 8 的预设总判仍 false。
- **自动规则**：保持旧离散宽和原门原样；另行预注册固定物理窗口、同连续峰定义的宽度与误差界。报告 nodal column offset、连续宽度与共同物理列的分解。核心实际事务门和原线性场门不能被新的观察器绕过。
- **回归**：峰相位 0/.23/.49/.77，尺度与幅值各跨数量级；多峰/扁峰/窗口边界退化明确拒绝；旧 false/辅助 true 并存的报告序列化不得改写旧判定。
- **限制**：C1 峰和 90% 水平集仍有退化/根选择问题；稳健 smooth 宽泛函是额外观察定义，不能倒改先前的 90% 门。

### F07：PROFILE 形状稳定，真实形状导数仍对网格敏感且非单调

- **触发**：归一化形状差很小，`U_tau`/残差场差较大；改变宽度观察器后率更稳定，残差场仍不随三档加细单调靠近。
- **原因**：高阶空间导数、峰/宽时间导数与输运项抵消对误差更敏感。标量残差范数稳定不能代替残差场稳定；观察器误差与 PDE 空间误差必须分离。
- **证据**：[56 原门判定](../../result/longtime/20260908_campaign_v1/fresh_tau724_target56_long_v1/run/final_verdicts.json) 保留 legacy=false/aux=false（原 linear rhoX L2 `.00204027>.002`，Hermite/连续宽/峰率分别通过）；[三档 smooth 率与范数](../../result/longtime/20260908_campaign_v1/acceleration_lab/smooth_cache_observer_tpc387fcb7_e2ef_4b66_8356_0bf5971dcb56/report.json)；[同物理时刻完整残差场](../../result/longtime/20260908_campaign_v1/acceleration_lab/smooth_cache_field_comparison_tpa0360d92_e1de_45f6_a48d_81a4bd434e40/report.json)。
- **实际结果**：在绝对 `t=1.964631495951036`，smooth primary 的 beta 粗/42/56 为 `-.621658299/-.621636625/-.621632706`；核心残差范数约 `.00743844/.00747899/.00750790`。但 core/holdout 场差粗→42 为 `7.9748%/7.5074%`，42→56 为 `9.7139%/14.0199%`，未单调收敛。旧 90% 观察定义的结果也保留，不与 smooth 混合拟合。
- **自动规则**：把“可以继续可信推进”和“残差已空间收敛”分成状态。保留真实 `Omega/FX`、配对尺度/轴/时钟，用固定观察协议作至少三档敏感性；不让拟合或更平滑观察器自动提升物理收敛资格。
- **回归**：[smooth tiny 入口](../acceleration_lab/ipm_accellab_test_smooth_inner.m) 的 MMS、单位/平移协变、中央时间差分二阶；实际缓存的 `Omega±h*FX` 方向检查；注入范数相近但方向相反的场，验证不能被“范数稳定”放行。
- **限制**：三档均从同一晚态迁移后演化，是晚态空间敏感性，尚非独立 t=0 网格收敛；局部细化比也不等于统一全域 h 比。当前证据不能认证极限残差或加速收敛。

### F08：小箱、fresh 初值与物理时钟混用

- **触发**：拿 fresh local physical time 与原生 absolute time 直接比较；将截箱后的数据套入父 CP 并继承不存在的历史；在不同箱子更换内层网格或自然 epsilon 后把差异全部归为边界误差。
- **原因**：晚态截箱是新的有限箱初值问题。其 `rho0=exp(-parentLogComega)*parentRho`、轴除以 parentCx、anchor 除以 parentCx，新本地时钟从 0 开始，父历史没有被求解。
- **证据**：[初始多箱内层轴审计](../../result/longtime/20260908_initial_box_axis_audit_v1.json) 显示旧多箱仅部分平台完全相同，全部 inner nodes 不一致；[tiny 五支协议](../../result/longtime/20260908_fresh_box_protocol_tiny_v2/protocol_report.json)；[最新 epoch 五支自然/公共 epsilon 对照](../../result/verification/performance_lab_20260908/late_box_v4_stage002_epoch_v1/dynamic/protocol_report.json)；[更长三箱对照](../../result/verification/performance_lab_20260908/long_box_v4_stage002_epoch_v1/dynamic/three_box_report.json)；[修正 epoch 的真实比较](../../result/longtime/20260908_campaign_v1/fresh_core_patch_long_tau724_epoch_audit_v2/epoch_observer_report.json)。
- **自动规则**：fresh 分支独立 caseId、初始采样、lineage 与本地 0 历史；`absoluteT=parentEpoch+localPhysicalTime`；`parentEquivalentTau=parentTau+lambda0*localCanonicalTime`，`lambda0=parentCx/parentComega`。比较用匹配的实际 absoluteT，率明确 fresh/parent/physical 单位。截箱先冻结场筛查，再同 inner nodes 的动态 matched-time 分支，并分离自然与共同 epsilon。
- **回归**：tiny 全箱 fresh 重启对原 CP 续算的物理协变；敏感小箱必须拒绝；故意遗漏 epoch、错误 parentCx/anchor、错初值/假继承历史都应拒绝；所有分支严格到达真实终点。
- **限制**：静态冻结 crop 不是动态边界验证；`farVelocityRatio` 不是边界误差。晚态短/长区间通过仍不能授权原 t=0 小箱。当前同网格原生续算不能直接变箱。

### F09：冻结依赖缺失被误归为候选族耗尽

- **触发**：前几段同网格正常，首次需要设计时才调用未冻结的 helper；所有候选捕获同一程序错误后最终报 family exhausted。
- **原因**：只对顶层函数做 `which` 或 Code Analyzer 不能证明运行依赖闭包。当前 driver 的候选 catch 较宽，可把 `MATLAB:UndefinedFunction` 记成普通候选失败。
- **证据**：[真实 v1 第一次失败](../../result/longtime/20260908_campaign_v1/fresh_adaptive_campaign_v1/stage_003/design_01/campaign_failure.json) 缺 `mesh_general_rounded_axis`；[保留的 v1 总状态](../../result/longtime/20260908_campaign_v1/fresh_adaptive_campaign_v1/campaign_status.json)。前两段真实完成，失败不能解释成数学上无可行网格。
- **自动规则**：启动前在冻结源内实际执行一次会走完所有几何依赖的 no-LU 设计；分别登记数值不准入、资源耗尽和程序错误。后两类不能盲目试完所有参数。冻结清单、哈希、`which` 检查及隔离 cwd 都需要保留。
- **回归**：删去 `mesh_general_rounded_axis.m` 的独立复制 fixture 应在预飞阶段报缺依赖；有完整依赖但所有轴不合格的 fixture 应报不同的数值拒绝；不得启动第一步 PDE。[performance_lab 预飞实际 v2](../../result/verification/performance_lab_20260908/campaign_preflight_v2/regression.json) 已通过 8/8，包括隔离 ambient 同名函数、历史漏依赖、全部工厂失败、轴门、输出目录和错标原 t=0；实际 step3668 的 14-link 原生链重读也通过，无 LU/PDE。
- **限制**：当前完整 v2/后续冻结解决已知漏文件，但还不是一个自动依赖生成与验证系统。预飞也不能代替真实迁移/长期运行。

### F10：低核心预测不足、有限候选与分段驱动仍需要人工接续

- **触发**：固定 `.2` 段内核心下降比拟合快；无足够历史时零衰减 fallback；全部预注册族失败；达到 maxStages/maxSteps 后仅保存状态。
- **原因**：当前 `fresh_profile_campaign_review` 只用同 remesh 最近 `.35` 等价 tau 的 log-core 斜率，无预测误差界或平台边缘预测；`run_fresh_profile_campaign` 固定四族和 895×386，并以 `batch_review_required` 结束批次。
- **证据**：[实际预测及资格 no-LU 测试](../../result/longtime/20260908_campaign_v1/fresh_profile_driver_nolu_v3/report.json)；[当前 driver](run_fresh_profile_campaign.m) 和 [review](fresh_profile_campaign_review.m)；前述真实窄族与平台族耗尽。
- **自动规则**：在原核心硬停机之前，按核心预测、峰移动、平台剩余距离及 q 裕量选择下一检查间隔；无趋势时缩短探测段，不能把零衰减解释为安全。预注册顺序包含移动 patch、不同细区节点分配、轻量有限搜索与减小段长。控制器必须在同一任务内自动滚动批次，错误和真实不可行则保存可信 CP 并报告实际限制。
- **回归**：人工可解析的平稳、指数收缩、突然加速和移动峰序列；触发前后的滞回；未到终点却达到 maxSteps 时必须失败；候选全失败仍有精确可恢复源；完全可信链可以无人干预跨批次恢复。
- **限制**：现有局部预测是规划工具，不是 PDE 保证；安全退出和短桥接不能被计为从零到预定长时的成功。新的段长算法须预注册，并与原 CFL/maxDt 保持分工。

### F11：原生签名、配对与历史来源的运行环境依赖

- **触发**：10 线程生成的旧 CP 在 4/1 线程下读签名失败；只比较 caseId 或 JSON passed；反序列化的匿名初值句柄身份不等；重网格后沿用仅支持同网格的 bootstrap 资格。
- **原因**：旧 numeric signature 包含浮点归约；句柄身份不是捕获数据相等；原生事务只有终端历史记录允许替换，普通续段要求完整前缀相等。
- **证据**：[原检查点审计约定](mesh_findings_20260908.md)；[完整资格子门负例](../../result/longtime/20260908_campaign_v1/fresh_profile_driver_nolu_v3/final_qualification_checks.json)；[缺段/错误 sourceCP 负例](../../result/longtime/20260908_campaign_v1/fresh_qualification_chain_nolu_v1/test_report.json)；[新链准入入口](fresh_profile_campaign_admit.m)。
- **自动规则**：本批旧 CP 原生 read/restore 保持 10 线程；不绕过签名。逐值核对 rho、配对 x/y、全部尺度与时钟、完整历史、数值配置、初值 `functions(handle)` 的完整捕获数据；事务 step/时钟/尺度不变、remeshCount 精确 +1。references 的数值不变，仅 origin/pin 最近节点索引可按原规则更新。标量监视器明确不具原生验证资格。
- **回归**：错 rho/axes/scale/time/初值捕获/lineage 的逐项负例；合法匿名句柄重载；非节点 pin 的最近节点规则；同网格前缀和重网格仅终端替换；链缺段不能跳过。
- **限制**：跨线程/跨平台签名可移植性仍需未来单独版本化解决，不通过调松比较误差解决；当前链帮助程序仍带 fresh 分支契约。

### F12：多 LU 内存、时间步成本与资源失败

- **触发**：旧/新网格 LU 同时常驻，多条 q512 辅助任务同时 restore；重网格加密极少数单元后 CFL 大幅收紧。
- **原因**：稀疏椭圆分解内存远大于冻结轴/场数据；最小输运局部尺度限制整体步长。把“核心点更多”作为唯一优化目标会忽略总推进成本。
- **证据**：[最初被资源协调中断的记录](mesh_findings_20260908.md)；[现事务 helper 提前释放旧 poisson](../grid_lab/ipm_gridlab_regrid_checkpoint.m)；[42 实际 267 步/697.7s](../../result/longtime/20260908_campaign_v1/fresh_tau724_target42_long_v1/run/fine_branch/branch_report.json) 与 [56 实际 355 步/925.4s](../../result/longtime/20260908_campaign_v1/fresh_tau724_target56_long_v1/run/fine_branch/branch_report.json)，两者同物理区间。
- **自动规则**：先 no-LU 筛轴/配对，再串行做实际迁移；释放旧 factor 后构建新 factor。设计评分同时记录预期时间步、q 裕量和可用资源，但成本不能覆盖可信度硬门。动态演化和瞬时缓存审计各有明确独占额度。
- **回归**：小数据监视旧 poisson 在新 build 前已释放；空资源预算必须在 LU 前拒绝；缓存 allowlist 禁止保存 factor；中断与正常退出状态分别记录。
- **限制**：上述耗时是本机同源实际观察，不是性能模型；矩阵近奇异警告需要配对真实 Poisson residual 审计，不能仅靠警告文本放行或否决。

### F13：原始粗网格观测失准，或者连可迭代的候选都没有

- **触发**：初始321×161均匀轴的单cell padding比真实核心宽；H32原top3实际解析核心不足31；H64/H128原X有限族为0，不能继续“从失败候选迭代”。
- **原因**：观察分辨率与求解节点预算被混为一件事。增加箱长却保持粗观察N，会改变局部导数、五节点前沿平滑的物理跨度及padding；这不等于真实解析初值没有合格求解网格。
- **应对规则**：只在原t=0、原初选有限族耗尽或全部候选仅实际core/front失败时，允许显式一次固定解析观察。原source轴、plan、失败原因保留；观察轴独立采样原解析初值，候选再以其自身轴采样，并通过实际native flow初始化和全family资格。不得把观测场写作PDE接受态，或在演化途中回采初值。
- **实际证据**：独立候选源码的[H64报告](../../result/verification/autonomous_runtime_20260909/initial_observation_fallback_candidate_v1/native_H64_v1/report.json)与[独立只读审计](../../result/verification/autonomous_runtime_20260909/initial_fallback_H64_independent_v2_9d44a7cbb7/report.json)：原source为321×161、core1.613515/5.849160、X0/70/Y3/4；一次559×264解析观察得到X3/70/Y3/4及3个pure候选。第1个native候选成功，初态core34.379507/36.786509、front24.963657；5帧真实场均161×321，4步至τ4e-5，remesh/growth均0。原source、observer、selected三套特征由独立解析重算逐值核对，18项门全过；另有157个源/输入文件前后SHA不变及完整冻源manifest复验。
- **另两条真实分支**：[H32](../../result/verification/autonomous_runtime_20260909/initial_observation_fallback_candidate_v1/native_H32_v1/report.json)因`actual_analytic_core_or_front`使用一次fallback后四步通过；[H16](../../result/verification/autonomous_runtime_20260909/initial_observation_fallback_candidate_v1/native_H16_v1/report.json)原初选成功，`observerUsed=false`，没有多余观测。
- **回归与限制**：12配置/上下文、14损坏ledger、8附加eligible反例及2纯合成ledger正例保留；合成不等于native。H256仍失败；k8/anchor1/双奇/t0与观察预算是显式规则范围，不是通用初值算法。任意合法solver N/预算可配置，并非永久硬编码321或离散H。完成F17后root已将显式可选规则合入维护源码；旧省略字段语义不变。[H64原t0长时注册](../../result/verification/autonomous_runtime_20260909/from_zero_H64_tau12_auto_initial_v1/registration.json)已真实启动，目标τ12不是已完成结论。四步不是大箱长时误差证明。H64独立审计首试脚本与变量同名在读取前失败，旧目录保留，未重复PDE。

### F14：计算箱固定，物理箱却随尺度缩小，登记远场覆盖变成零

- **触发**：各向同性、零平移下`Hphysical=Hcanonical/Cx`持续下降；仍以计算坐标“大H”或当前箱的固定百分比声称覆盖原物理远场。
- **实际证据**：[原生物理远場覆盖审计](../../result/longtime/20260908_campaign_v1/acceleration_lab/physical_far_coverage_eedb04e0cb5d4258be7b2a227c5e74ba/run/report.json)严格配对H8从零τ8/step4545：Cx4.626468、物理半宽1.729181、物理高度.864590、最大半径1.933283。原登记物理环带[4,6]与[6,8]均0个节点、query coverage=0。按当前箱比例取出的环带虽有节点，半径小于4，不能替代远场条件。所有native/trust门仍通过，因此两种资格必须独立。
- **应对规则**：预先登记固定物理观察窗口及覆盖/方向/原生间距要求，逐接受态或有界检查点重算实际域。覆盖不足时记录`unknown/uncovered`，不能把无数据记为零误差，也不能仅靠提高core目标或增节点修复物理域缺失。改箱是单独PDE边界契约，不能混入同箱regrid。
- **回归与限制**：相同计算箱而Cx不同的纯单位例；有核心、有近箱环带但固定远场覆盖为0的真实负例；fresh空间轴按scale/shift、时间按epoch分别变换。有限环带即使完整覆盖，也不能证明盒外加权界沿PDE传播。`A/H`条件率尺度不得直接加入C；完整外尾Dirichlet迹与固定锚点率有不同误差阶。

### F15：有限搜索的参数间隔造成假容量结论，anchor位移会推走核心细区

- **触发与原因**：fresh step9205、等价τ12.44320075原70个X全拒，其中56项触相邻比、52项core不足、32项front不足（原因可重叠）。全质量合格14项中最好core仅26.9176。某些原anchor warp在当前核心处移了1.38–2.70个core宽；仅“精确anchor已保留”不保证细区仍覆盖核心。
- **证据**：[原族逐门总结](../../result/longtime/20260908_campaign_v1/performance_lab_fresh_stage5_geometry_v1/original_axis_summary.json)；另行预登记fine36/40/44/48 × rounding8/10/12/14/16 × fraction.5/.65共40项，复用2项、仅新算38项，原warp及全部门不变，[10个X通过](../../result/longtime/20260908_campaign_v1/performance_lab_fresh_stage5_geometry_v1/passing_candidates_summary.json)，未触发另9项后备。不能把此新族成功回写为旧70族成功。
- **实际资格层次**：[前三完整冻结配对](../../result/longtime/20260908_campaign_v1/performance_lab_fresh_stage5_pairs_v1/run/report.json)均过；[首项真实native事务](../../result/longtime/20260908_campaign_v1/performance_lab_fresh_stage5_native_trial_v1/report.json)保持895×386/同箱/全部时钟，core35.124212/35.963732、front约39.974，峰跳8.3466589e-6、累计.0005055596954、质量缺陷0、相对范围越界4.4293e-7，严格CP回读、history/ref/lineage均通过；该事务没有推进PDE。
- **应对规则**：冻结有限候选族时给出系统的fine/rounding范围及稳定顺序；原族全无轴对后才进入独立注册后备，最大尝试数仍固定。先全轴和区间覆盖，后完整field，最后native事务，三种通过不可互代。新增族须按版本从未来配置启用，不在运行中偷偷扩展。
- **回归与限制**：保留原70全拒而40中可行的真实边界，以及axis pass/field fail/native fail分层负例；已成功原族不能被后备抢走默认选择。新首项是晚态恢复候选资格，不是从原t0无需人工改族的工程验收，也不证明所有轴设计均被有限搜索覆盖。
- **自动规则的独立实测**：[预登记分层搜索](AUTONOMOUS_PARAMETER_SEARCH.md)从构造器最小fine/rounding及两侧分支cell约束生成下界一倍尺度拓展和两级二分，最多500个不同X元组。9205原70报告逐值复现、L1新增260后停止，前三包含原40族之外的fine16/round16；8057保持原70/3候选逐值相同、后备不运行。其前三完整冻结场门均通过、首项noLU原生预飞通过，实际原生事务尚未运行。500耗尽只能表示登记搜索预算耗尽，不能证明无可行网格；此研究实现尚未接入从t0运行期，物理箱闭合也仍是另一能力。

### F16：增点解除容量停止，但同物理终点仍有显著箱/网格联合敏感性

- **原失败**：H16/target32/cap110000旧B保存step1520，未达登记物理终点；未提交1521的Y候选全部core不足，level4受资源mask禁止。旧失败/CP不升级。
- **新实际运行**：[全新t0、cap210000 B2](../../result/verification/autonomous_runtime_20260909/expanded_zero_controls_cap210000_v1/plan/factorial/runs/B_H16_T32/endpoint_report.json)已到t1.9140859724939359，step2046/τ4.66701456，641×321/level4、16次regrid/2次增长，core29.702984/32.670913、安全.1543043、累计峰跳.00092929665；原生全门通过。cap变化只打开预登记成员，不是任意新N，也不是RAM模型。
- **科学失败保留**：[新A/B实际比较](../../result/verification/autonomous_runtime_20260909/expanded_zero_controls_cap210000_v1/AB_actual_comparator/report.json)两者均`native_qualified`并覆盖两物理比较窗口，科学门却均false。物理rhoX最大值A4.0590846、B6.2994455（约+55.19%）；连续横宽A.01216547、B.00607397（约−50.07%）。physicalGradInf另为4.8904366/18.9873486，不能与rhoX混称“G”。
- **应对规则/回归**：最少H8/H16×target32/42及独立半步时间控制在同物理终点比较；双向linear/legacy与Hermite/continuous各自保存，D−C−B+A在固定query有符号求差。缺失C/D/E明确pending，不插值预测终点代替实际结果。把“原生达终点但科学门false”列为必要负例。
- **限制**：A/B初选轴、增长次数、实际core与物理域均不同，是箱/空间/自适应事件联合敏感性；既不能称纯箱隔离，也不能仅归因Green边界。t0壳积分揭示外尾机制，但不认证演化场尾、PDE空间阶或长时PROFILE。

### F17：研究审计/测试协议失败不能误记为数值路径失败

- **合法空账本**：旧H32预选四步已完成，研究auditor访问空`transactions.sourceLevelId`出错；独立只读修正后复用已存CP/result通过，没有重复PDE。
- **合法NaN**：[default-v1逐项诊断](../../result/verification/autonomous_runtime_20260909/initial_fallback_default_review_4dd528ad81/report.json)证明原`isequal(history)`是唯一失败项。六个初始步诊断包含合法NaN；完整class/shape/值及NaN位置相同，其余9组原isequal全true。新harness按严格class/shape与等NaN位置比较，不忽略字段；既有四步v1数据复用。
- **终点协议真的不同**：[H64 split逐组诊断](../../result/verification/autonomous_runtime_20260909/initial_fallback_split_review_e35ce47562/report.json)只见`history.common.canonicalEndpointDtLimit`前三项不同：原prefix把finalTime4e-5改为2e-5。12组native及state/physical/scale/grid/snapshots/controller/gauge均相同，但完整history=false保留。新独立suite必须保持同finalTime4e-5，仅maxSteps2截断；不能豁免这个有数值语义的诊断。
- **schema测试漏覆盖**：候选新增显式可选初选规则后，旧baseline测试仍只允许原kind及105/101/4计数而失败；真实为106选项/101默认/5可选。独立candidate_v2仅更新测试白名单和精确集合，全部`+ipm`与已锁定v1字节相同；原full-suite失败保留。更新期望不是宣布整套回归已通过。
- **后续实际通过**：[省略新选项的v1/v2默认回归](../../result/verification/autonomous_runtime_20260909/initial_fallback_native_harness_v2_fec7afeab9/default_candidate_v1_v2/report.json)原生配对与所选数学字段逐位通过；[all/equivalence运行日志](../../result/verification/autonomous_runtime_20260909/initial_fallback_full_suites_v2/run.log)及同目录MAT报告保留实际套件结果。新H64 split保持相同终点，只以maxSteps截断；其第二次审计误把rescaling字段顺序视为数值差异，[完整只读补验](../../result/verification/autonomous_runtime_20260909/initial_fallback_split_v2_readonly_review/complete_report.json)按字段集合/class/shape/值/NaN位置严格核对，完整history/snapshots/cursor全部通过并含5反例。恢复仅1次mesh.build、没有initialDensity/planner重采样；只读补验未重跑PDE或LU，两个旧split失败均保留。
- **规则与限制**：记录失败发生阶段、已保存可信产物、独立复核结果和确切未通过项；测试修复放新目录，原源码/报告不覆盖。字段顺序可作为表示差异，不能因此放宽字段集合、类、形状、值、NaN位置或有物理语义的历史记录。默认/split/套件资格不能替代长时空间/时间/域误差验证。

## 当前代码：可复用部分与不能直接泛化的部分

| 部件 | 可以复用 | 通用化前必须处理 |
|---|---|---|
| `+ipm/remesh/adapt, transfer, interpolate` | 提案→迁移→验收→提交；失败回退；高阶一致迁移与守恒校正 | 当前 `validate` 的 dynamic `peakSafe` 直接成立，核心主要只守旧值约 98%，缺当前 campaign 的绝对目标、front、累计峰跳、局部 q 门。不能只打开 `adaptiveRemesh` 就声称具有研究协议的全部门。 |
| `+ipm/evolve/initialize` | t=0 解析初值在接受的初始网格上重新采样，建立本案真实 references | 当前初始网格选择沿旧 `adapt`，尚非新 core-patch 控制器；t=0 不应转移旧样本或伪造 parent crop。 |
| `+ipm/evolve/remeshIfNeeded` | 每次接受步后触发、重建真实 flow/RHS | 只用当前 safety 和最大 remesh 次数，无提前预测/候选资源预算。 |
| `mesh_general_rounded_axis`、`mesh_core_patch_axis` | 无 LU 纯几何；任意正单位、奇数 x 节点；core/anchor 可分离 | 新候选构造器，不是历史截箱轴的逐值生成器；单峰、对称、正 anchor/中心及最小节点预算仍有约束；整轴质量必须外部复核。 |
| `ipm_gridlab_equalize_axis`、`feature_profiles`、`score_frozen_pair` | 轴均衡、实际场的核心/前沿特征、完整配对筛选 | 固定单调搜索和一套 sigma 非普适；需要有界多族搜索及错误分类；冻结 transfer 不能代替原生事务。 |
| `mesh_platform_candidate`、`mesh_design_stage` | H=1e6 历史晚态的诊断/回放证据 | 前者硬编码 anchor=1、512 正半轴预算和旧 root factory；后者要求 root dataset。不得直接用作任意 t=0 初值的默认设计器。 |
| `mesh_design_fresh_core_patch` | 实际 lineage、immutable 初始轴、精确 anchor、成对全门设计 | 必须有 `latePhysicalBoxBranch` 与父 CP/截箱样本；目标至少21、前沿20，并依赖单一 snapshot。不能拿它给原始 t=0 算例捏造 fresh lineage。 |
| `ipm_gridlab_regrid_checkpoint` | 真原生同箱同节点事务；时钟、references、历史终端替换与新签名 | 是研究 CP 工具，会合法重建网格配置；不适合每步恢复/重建 LU。通用在线控制应复用其验收逻辑并整合进原 solver 事务，而非另写 PDE 路径；增节点/变箱明确禁止。 |
| `run_fresh_profile_campaign`、review/admit/chain | 精确配对、同 remesh 预测、累计预算、明确失败保存、批次链 | 当前固定895×386、fresh时钟、caller资格、四族、分段和人工批次结束；只能作为控制器原型。当前 admit 已支持重查实际链，旧仅同网格限制不可再当成已修复后的现状。 |
| 连续与 smooth 观察器、真实 RHS cache | 分离几何采样与实际 PDE 残差；no-LU 回归 | 没有授权替换原门或控制 gauge rates；观察器通过不等于 PDE 残差收敛。 |

## 从 t=0 使用的设计接口

研究级通用配对接口现已实现并完成三实际态 no-LU 回归，见
[MESH_PAIRED_DESIGNER.md](MESH_PAIRED_DESIGNER.md)：原t=0和当前fresh态找到全门冻结候选，
旧大箱step3326在固定节点数下明确报告y族q-floor容量失败。它不要求fresh lineage，
但尚未集成下面提议的在线控制器或执行从零到长时的PDE验收。

保持唯一数值入口 `ipm.solve(opts)`。新网格控制器属于 `+ipm/remesh` 内部，数值模块不依赖 research 文件。先在独立研究原型完成以下契约，验证后再迁入核心；不能修改运行中的冻结副本。

纯规划接口可写作：

```matlab
[decision, nextMemory] = meshPlan(stateView, initialReference, policy, memory)
% decision.kind = keep | propose | shorten_review_interval | infeasible
% decision.candidates = ordered x/y pairs + construction/quality reasons
% 此处没有 LU、迁移、规范率修改或 history 写入。
```

`stateView` 提供真实当前 rho/配对轴、尺度/全部时钟、可信特征、精确 transport anchor 和可用资源；`initialReference` 保存本案 t=0 解析初值描述、实际接受初始轴与 references 来源；`policy` 冻结硬门、目标核心、前沿、候选预算、终点和允许节点/箱策略；`memory` 保存同 remesh 的最近趋势、上次接受/拒绝、实际累计峰跳及检查调度。memory 进入版本化 checkpoint 并可逐值恢复，不能在批次重启时归零。

在线执行使用现有的 `propose→transfer→validate→commit` 骨架：先筛有序候选，实际迁移后按统一硬门重算，成功才增加 remeshCount、重建真实 flow/RHS 并写审计。普通接受步保留当前状态；无解时不得反馈修改 `c_l/c_omega/c_r`。将纯几何 helper 的实现迁入适当数值模块之前，须有旧 solver 逐值一致性和迁移回归。

初始阶段先用解析初值选择适用的单峰/多特征规划策略；在最终接受的初始轴上重新采样解析初值，生成新原始 caseId、t=0/history0，不需要 parent lineage。保留全部初始网格尝试，但运行 references 以最终真实初态按原规则建立。晚态 fresh 分支继续使用其独立、明确的 lineage，二者是不同合法入口。

候选顺序应预先登记为：保持网格 → 随核心移动的 patch/anchor 分离 → 不同细区/过渡节点预算 → 有界均衡或单跳 current-axis 备选。具体 family 不应在运行中手改。检查间隔由保守预测及边缘距离决定；缩短的是网格复核调度，不擅改 CFL、maxDt 或物理时间单位。若要允许增节点，先单独实现并验证真实原生增节点事务与存储/内存预算；现阶段不能假设已有该能力。

从零时间的盒子必须在启动前用完整 inner-node 一致的盒子梯度选定并冻结。晚期自动截箱不应作为同一原生 IVP 的无痕 fallback。若最后采用显式多箱分支算法，验收描述必须写明它解决的是带受检验边界近似的分支问题。

## 最小连续集成与最终验收

| 层次 | 最少场景 | 通过条件/产物 | 当前状态 |
|---|---|---|---|
| 每次改动，无 LU | helper 闭包、策略类型、精确 anchor、对称性、单位缩放、节点预算、y 单独失败、移动/收缩趋势、不可行原因分类 | Code Analyzer；全部质量门；同一输入确定性输出；缺依赖在 PDE 前失败 | 既有几何/协变证据可复用；performance_lab fresh 预飞已 8/8，通用 t=0 设计仍需独立验证。 |
| 每次改动，小网格 | 原生事务接受/拒绝、反复往返、累计预算、解析 t=0 重采样、完整 checkpoint 恢复、同网格前缀/事务终端替换 | 真实原生签名及配对；拒绝精确不改变源；至少跨两次自动重网格与一次批次恢复 | 部分 tiny helper 已有；新的从零控制器尚未实现，不能写成已通过。 |
| 观察器与缓存，无新 LU | 相位/单位/MMS、smooth 方向差分、残差范数相同而场不同、legacy false 保留 | 预先登记误差界；正确错误分类；缓存无 factor；不能凭 norm 稳定升级残差资格 | 既有 continuous/smooth 与实际三档证据可回放。 |
| 合并前，小/中网格独立 t=0 运行 | 同原始解析初值，两档空间×两档时间；移动峰/强收缩；固定初始协议自行跨多个 regrid | 到真实共同物理终点、完整 trusted、各层门全部满足；0人工参数修改；保存全部候选拒绝与时间开销 | 待新接口完成后注册，当前晚态 smoke 不替代此层。 |
| 发布/科学验收，排队大规模 | 已选盒子至少两个独立初始网格，从 t=0 到预先登记长时终点；必要共同 epsilon/盒子梯度 | 全段自动、原生可恢复；实际物理场/率/峰/宽及真实残差场的空间时间敏感性满足独立注册阈值；资源在预算内 | 未完成，主线现有长时结果与此项分开报告。 |

最小完成日志应包含 t=0 配置/源码哈希、每次规划触发原因和预测误差、候选拒绝轴与质量门、实际事务前后场/网格/时钟、累计误差、每段真实终点及 stopReason、链审计、无人工修改记录。最终验收不能只统计重网格次数、最大梯度或 `exit0`；也不能把候选全失败后的安全退出标成目标达成。

本文所有历史链接应随文档检查实际存在；未运行的接口和 CI 项目明确为提议。已有失败文件、旧判定与冻结源码保持不变。

## 2026-09-09 新增实测：长尾节点预算与验收边界

v3 step4768后四旧候选全部在x相邻比1.08门失败，范围1.08199–1.09218。
64/96个fine与40个rounding占用尾部节点，warp之前log-spacing斜率已超过log(1.08)；
因此不能只调整锚点warp。旧y sigma .25预测core31.72也不足32，只是被更早x拒绝遮住。
通用搜索自动选择32 fine/16 rounding/.5 fraction/y sigma .5，x/y相邻比约1.06674/1.06682。
实际原生迁移最终core34.7329/36.1246、front31.0700，跳变6.77e-5、质量缺损0，全部原门通过。
证据：../../result/longtime/20260908_campaign_v1/fresh_adaptive_campaign_v4/stage_001/regrid_audit.json。

实际候选负例另暴露NaN权重造成NaN比较误通过；现验证积分权重形状、实数、有限、正性及派生指标有限。
控制器审查发现forecast可能选择unchanged网格，重复分解并消耗事务预算；已禁止演化同轴提交，
并让高于目标的短窗预测只继续观测。相关独立回归仍记录真实通过范围，详见运行期集成记录。

## 后续实证：控制器窗口与节点预算

- [长期窗口报告](AUTONOMOUS_CONTROLLER_WINDOW_BOUND.md) 明确0.35时间范围不是独立于dt的内存界，给出真实native计数、无LU微基准和固定预算采样/CP重播的待集成规则。from-zero最终验收必须包含高接受步数的同epoch窗口，不能用重网格后只剩1行的末态证明有界。
- [原t=0容量停止与单轴加点](MESH_REFERENCE_NODE_LADDER.md)：step1134/tau2.688619的321×161原x族0/70、y2/4；只加y仍失败，641×161首次通过原几何与明确离线迁移门。reference family嵌套、actual网格移动和变量N原生事务必须分别审计；未生成新CP或完成动态资格。
# 补充：原始 t=0 粗测几何使初选预测过于乐观

`initial_box_geometry_v2_v3` 的同 k8/321×161 三箱纯解析筛查中，H8/16 的原 top3
均过门；H32 三个 top3 的 Xcore 预测34.76185、重新解析实际27.63098，全部被原31门
拒绝，轴质量及所有 reference family 成员仍全部通过。H32 的源 dx=.2，单侧
`max(rawHalfWidth,local_spacing)` 把核心设计左端从原始交点1.126836移到1.0；
重新解析后物理设计宽缩小20.56%。固定五节点前沿平滑也随网格改变物理跨度。
这属于原始粗源的几何估计失配，不能直接标成固定 N 容量耗尽或箱误差。

应对规则是保留实际解析重采样门；将有界 t=0 解析再设计作为独立未来协议，不能降门
或把未验证的其它候选当成功。回归需要 H32 三候选实际拒绝、H8/16 正例、完整
family 质量以及“完成筛查但并非全部箱通过”的状态分离。详细证据与限制见
[INITIAL_BOX_ANALYTIC_SELECTION.md](INITIAL_BOX_ANALYTIC_SELECTION.md)。

## 2026-09-09 实证：第二方向自然增长与最高级后的容量

- **触发与真实事件**：`from_zero_tau8_v2_cap210000_v1` 的第二增长实际发生于 **step2477**，
  canonical τ=5.594870008425091、physical t=2.1092457725037628，remesh20→21，
  level2 的641×161→level4 的641×321。step2481/τ5.600240554933915是首个已保存增长后CP，
  不是事务发生步。该案第一次增长是step1135/τ2.6913754578032374，321×161→641×161。
- **因果规则**：保存的 `current_core_trigger` 源核心为25.99223592/31.87618618，
  X低于原26门；.2τ预测22.72173306/29.41753709。X最大衰减割线来自step2468→2469、
  Y来自2421→2422，均在事件之前及原.35τ窗口内。保存的fit/最大割线/指数预测代数均复核。
  触发的是X核心门，不应据此断言旧Y轴哪一质量门失败。
- **真实事务值**：实际attemptSummary仅candidate1且通过；新核心35.16633484/35.26912843，
  front34.24239898，峰跳8.092688637e-5，累计由.00151371138460到.00159463827097。
  原生读取器验证完整ledger、原门、时钟、源/目标level、固定参考family和初始化来源。
  没有改箱、重置时钟、放宽质量门或反馈修改规范率。
- **证据**：[完整无LU审计](../../result/verification/autonomous_runtime_20260909/second_growth_capacity_audit_df297f1061/run/report.json)、
  [保存的第二增长原事务](../../result/verification/autonomous_runtime_20260909/second_growth_capacity_audit_df297f1061/run/second_growth_raw.json)、
  [输入与冻源SHA复验](../../result/verification/autonomous_runtime_20260909/second_growth_capacity_audit_df297f1061/final_byte_audit.json)。
  session28745实际exit0，CodeAnalyzer=0，profile无mesh.build/flow/transfer/advance/solve；
  只重建小型配对Dx与轴诊断，没有椭圆LU。
- **当前容量快照**：注册时最新完整CP为step2906/τ6.200020140905038/t2.1551574783515375，
  core32.71735806/34.42634536、安全.1327739892。
  `common.physicalGradInf=max hypot(rhoX,rhoY)=12.75979941`；横向导数的
  `common.physicalRhoXInf=max abs(rhoX)`约10.2579，二者不能共用含混的“G”名称。
  [精确标量与名称说明](../../result/verification/autonomous_runtime_20260909/second_growth_capacity_audit_df297f1061/observable_names.json)
  保留二次峰的第三个独立值；旧`fresh_profile_campaign_review.physicalGradientMaximum`实际指rhoXInf，
  不能按字段英文名称误解为全梯度模长。该新增标量提取不是又一次原生签名验证；签名已由上述严格审计完成。
  最高level4有54/70个合格X轴、3/4个Y轴；top3含1个qualified keep及2个迁移提案。
  两提案预测core34.66335/36.8与35.16261/36.8，front37.49846与34.07395；**尚未实际transfer**。
  当前X/Y相邻比1.045374/1.035991，global q比.0121113/.00171484，均有原门余量。
  16个X试轴触及相邻比门，另有局部求积门拒绝；一个Y试轴不足目标core，全部失败证据保留。
- **应对规则**：继续使用冻结搜索与同level原生事务，失败保存最新可信CP，不擅加节点或换箱。
  641×321共205761点；cap210000剩4239点不能转化为未注册的新成员。
  两次增长额度已用完，且分轴单调条件下只有level4仍合资格；剩余278次remesh与
  .01834624989峰预算只允许进一步同级调整，不保证还能推进多少时间。
- **回归场景**：从真实t=0自动跨1→2→4、两次增长及其后的同级regrid，检查完整原生恢复；
  在最高级强制候选族无解时必须诚实容量停止；禁止回退节点级、伪造第五成员、超cap或重置累计预算。
  将“观察到增长后的CP步号”误当成事件步号、及稀疏history缺少事件行，列入只读审计负例。
- **未解决的证据范围**：此案`storeSnapshots=false`、`outputEvery=.01`，history事件两侧是
  step2475与2481，没有step2477本身；未保存同一步迁移前后完整rho和完整失败axisReport。
  因此本次严格核对签名、三个CP的history/ledger精确前缀和累计求和，
  但**没有独立重放历史transfer或重算该历史峰跳**；这些值来自已签名的原生实际事务审计。
  最新CP的场可独立重算当前核心/前沿与剩余几何，不能冒充过去事务瞬间的场。
  这证明有限族的两个方向增长已实际可用，不能证明无限自适应或任意长时间可解。

本次研究审计还保留三次未完成尝试：`second_growth_capacity_audit_dfb58e6785`的静态缩进告警、
`..._b6eeed4947`中错误假定增长在2481、`..._bf2515c250`中错误要求稀疏history必须含事件行。
它们是研究审计假设失败，不是生产事务失败；最终脚本明确记录可核验范围，不补写任何旧数据。

## 2026-09-09 补充：更大初始箱的观测失准与调度上限

固定321×161的初始粗观测在H64/H128都找不到X候选；预登记的有限解析观测格恢复了原门合格候选，求解节点预算不变。H256仍无候选，不能靠扩大粗网格测量误差来声称容量；详见[原始零时刻解析观测研究](INITIAL_ANALYTIC_OBSERVATION_PROBE.md)。后续预选H32/H64/H128及自动fallback H32/H64已有原生四步证据，显式可选能力在F17验证后已合入维护源码；省略字段仍保持旧默认。

另将研究调度失败与数值网格失败分开：fresh v4 stage6在step7510达到原512步上限，状态core28.79/34.08与原数值门通过，但未到原目标且stageAccepted=false保持。显式完整原生链审计后，v5用4096调度预算先补原目标，在step7545通过；没有改CFL/dt/网格门。当前v5 step8057/equivalent tau12.043200751859164通过，下一事务核心33.4496/36.9361、累计峰跳.0004972130365。调度恢复本身不算用户要求的单次从零自主网格资格。

### Latest actual capacity case: original B H16

See EXPANDED_ZERO_CONTROLS.md for the original cap110000 failure at attempted step1521 (saved step1520), the Y-core rejection, and the separately registered cap210000 five-case experiment. Old outputs are immutable. New runs all start at t0; no stopped checkpoint is upgraded.

A separate research auditor failed on the valid empty initial transaction ledger. growthCount now handles this as zero; H32 existing four-step outputs were audited without repeating PDE. H32/H64/H128 preselected native four-step qualification and the automatic H32/H64 fallback probes passed. The explicit optional fallback subsequently passed default/split/full-suite qualification and was integrated by root; the immutable candidate sources and original failures remain available. F17 distinguishes these checks from unresolved scientific accuracy.

最新版接口与必须分离的资格见[自动网格设计约束](AUTONOMOUS_MESH_ACCEPTANCE_CONSTRAINTS.md)。后补的实际证据优先于本文保留的早期阶段描述；原失败判定始终不覆盖。


## 2026-09-09 实际H8 τ11.673容量停止：有限搜索遗漏可行轴

原v2 τ8→12续算最后接受step11049（τ11.673078880623883、物理t2.30798220707281），下一未保存step11050的forecast预测X core21.9975988低于22，唯一可用level4的原X 70项全拒、Y 2/4通过，停止并保留原失败。11049峰rhoX76.86049、二维全梯度103.51727，原预定τ12未到。

记录中的X拒绝分类有重叠：66 adjacent ratio、32 core、9 front、46 local q min、30 local q max、2 global q/nonpositive。旧失败没有保存11050完整rho或每trial数值quality，不能重算该瞬间。

独立只读审查V/h8_tau11673_capacity_nolu_v1以真正保存的11049场作反事实：同641×321原70仍X0/Y2；预登记hierarchy到330项有20个X合格，首对fine16/round16/fraction.5/Ysigma.5，实际离线core33.61963/36.69159、front27.67842、峰跳2.72829e-5、质量守恒通过。fine16只是plateau构造参数，不等于实际core格数。X3成员961×321原70即有61个X；Y3成员仍需X后备。

当前源轴global q X5.54515e-4/Y6.30023e-5，远高于1e-8下界；本例不能归咎当前global q耗尽。总结为细区/过渡参数覆盖不足，需要普适有限分层搜索与分方向增长。没有修改旧CP/规则、没有把接受11049反事实称作失败11050迁移成功；新的原始t0冻结组合运行才可验证自动跨过这一障碍。
