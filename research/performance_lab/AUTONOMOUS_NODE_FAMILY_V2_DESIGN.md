# 自动加点策略 v2：只读接口审计与合同草案

状态：2026-09-09 独立设计，未实现、未启用、没有新 LU 或 PDE 计算。现有 v1 产物和续跑规则不变。几何筛查由 mesh_audit 独立进行，本文件不将其候选轴解释为已通过原生迁移。

## 已有证据与任务边界

`result/verification/autonomous_runtime_20260909/from_zero_tau4_v1/report.json` 记录原始零时刻起跑、321×161、13 次自然迁移、全历史 trusted，最终在 step 1134、tau=2.6886192715883674、物理时间 1.613029050260455 以 `autonomous_mesh_axis_capacity` 安全停止。末态 core=[26.000909525348369,32.342665683010424]，累计实际峰跳 0.00071060941279805351。该证据说明有限注册族在这个实际状态耗尽，不证明所有同 N 网格都不可行，更不证明奇异性。

目标是新运行从初始条件起，在冻结的资源和网格合同内自动增加必要方向的节点，继续同一 PDE、规范、物理时钟、case、历史和原生检查点。这个有限族仍可能耗尽；预登记有限加点能力不等于已保证任意长时运行。

## 固定 N 假设与可复用能力

以下位置以本次审计源码为准，路径相对工程根；审计伴随结果保存源码 SHA。

| 位置 | 当前行为 | v2 需要的边界 |
|---|---|---|
| `+ipm/+config/resolve.m:335` | 初始 `customX/Y` 长度须等于配置 `nx/ny`；双奇对称要求奇数 nx | 保留。配置描述初始预算，不重写为末态预算。family 是另一项经过正规化的策略/原生证据。 |
| `+ipm/+config/autonomousMeshPolicy.m:35` | 仅 version=1，固定搜索族和质量常数 | 按显式 version 分派；省略 version 仍为 1。v2 才接受加点合同；完整存储策略重读严格同类型、同值。 |
| `+ipm/+mesh/build.m:20,118,140,246` | custom 轴可输入，但计算网格、Kronecker 维数、数组分配及 ops.nx/ny 使用 override.nx/ny | v2 构造前要求 override.nx=numel(customX)、ny=numel(customY)。它是内部瞬时构造输入，不修改 config。当前代码不能只换 custom 轴长度。 |
| `+ipm/+remesh/transfer.m:16,58` | 两处 override 继承 config 中的初始 nx/ny | v2 验证 proposal 的注册级别后设置实际 override 计数；v1 原分支不动。 |
| `+ipm/+remesh/transfer.m:67` | 无条件保留旧 ops.baseX/Y | 同级维持原值；升级只能安装初始化已注册的 target member，不能把候选轴或当前轴重新当根。 |
| `+ipm/+remesh/interpolate.m:38,48,71,224` | oldAxis 约束 sampleDimension，queries 可有不同数量；两方向守恒修正使用各自旧/新权重 | 有变 N 的数组能力，仍需正式不同 N 的守恒、壁面、对称和真实峰跳回归。保留唯一 transfer 路径。 |
| `+ipm/+remesh/plannedAxisPairs.m:6,36,124` | reference 长度必须等于 current 轴；x 工厂/预算用 numel(current.x)，y 使用 reference.y | 分离“真实 source/view”与“目标计数/reference”。source.rho/source/Dx 始终按旧实际轴配对；targetNodeCount 从 family 取得，不能给 view 填虚构新轴。 |
| `+ipm/+evolve/planAutonomousMesh.m:44` | 唯一 reference=ops.baseX/Y，单层候选耗尽即停 | v2 枚举符合单调关系的注册 level；同级先试，再按实际总节点数递增。触发逻辑、预测窗、core/safety 门不变。 |
| `+ipm/+evolve/applyAutonomousMesh.m:23,53` | 初始化按初始预算选轴；成功后 selectedBase=实际 trial 基轴；之后迁移并递增 remeshCount | 初选 level 1 成功后、第一步 PDE 前生成并冻结 family。迁移成功才原子提交 level/base/ledger/cache，候选失败不推进 level。 |
| `+ipm/+remesh/auditCandidate.m:51` | 拒绝节点数或 base 任一变化；要求 sameBoxAndNodeCount=true | v2 用分别可审计的 sameBox、registeredNodeCounts、matchingReferenceMember 和 admittedLevelTransition；不能在变 N 时继续伪称 sameNodeCount。 |
| `+ipm/+remesh/controllerTelemetry.m:7` | memory.version=1，窗口只有 remesh epoch、core、safety | v1 保留；v2 增加 levelId/actualNodeCount，same-time commit 仍替换末观测，remesh epoch 更新后重置趋势窗。 |
| `+ipm/+remesh/validateController.m:48,67,78,110` | currentN=configN；initial selected/source 轴也与 current 同 N；base 必须等于初选；事务要求 sameBoxAndNodeCount | v2 单独校验 initialN、当前注册 N、初始根、当前 family member、分量单调 level 轨迹和按 epoch 匹配的快照。绝不能仅删除这些检查。 |
| `+ipm/+output/readCheckpoint.m:113,161` | base 长度=current；configN=current | base 长度=current 的规则可保留，但 base 改为有证据的 current member。configN=current 只对旧/v1保留，v2改为 configN=initialN 且 currentN=registered memberN。先无 LU 全校验。 |
| `+ipm/+output/makeCheckpoint.m:12`、`checkpointSignature.m:7` | 原样保存 config、metadata、x/y、baseX/Y、尺寸化 rho 签名和每帧签名 | 存储天然能容纳不同形状；metadata 中的 family/level 会进入现有签名。不要为 v2改变 v1签名算法或既有 payload。 |
| `+ipm/+output/restoreCheckpoint.m:35` | enabled 分支已按 saved 实际轴和计数一次 custom build，后恢复 base/runtime references 并重算 flow | 保留一次 LU 结构，在 build 前完成 v2 family/level审计。base 从已匹配 payload恢复；post-build flow/末窗口 exact 检查仍保留。 |
| `+ipm/+output/validateV2.m:31,51,614` | 场尺寸按实际 grid；快照逐帧按自身 x/y；enabled controller 当前仍固定 N | 通用 shape 检查不必弱化。v2 helper 要额外逐帧核对当时的 level/nodeCount，不能只核末帧。 |
| `+ipm/+output/record.m:25,142,317`、`finalize.m:46` | snapshots 是 cell，每帧独立 rho/x/y；末态也按 ops 实际轴；历史没有显式节点数 | v2 新历史增加 nodeCountX/Y 与 meshLevelId，始终记录实际 ops；v1 不注入字段。配置初始预算和末态实际预算分别解释。 |
| `+ipm/+output/checkpointFromResult.m:29,120` | 从配置重建 seed；末态可 build 实际N，但 base恢复为 seed初始base；seed LU同时仍活着 | 新版若支持此可选桥接，必须从已验 family 选末级 base，并修 LU 生命周期。未完成前明确拒绝 v2 结果转CP，不能伪造 lineage；原生 CP 本身足以无人续接。 |
| `+ipm/+remesh/propose.m:49` | 旧 adaptive 路径可能选 current/base，隐含同 N参考 | v2 仍使用自主 planner，禁止错误落入旧 propose/adapt 路径。旧路径不扩资格。 |

数值内核 `field.poisson`、速度重建、SSPRK、transport/RHS cache 等主要使用 ops.nx/ny、ops.x/y 和字段实际大小，不从 config 推断当前数组大小。`makeRhsCache` 存储实际轴、rho 和 remeshCount；升级后必须全部重建，不能继承旧尺寸缓存。`snapshotAt`/结果比较和可视化也以逐帧/实际轴工作，但需要一个跨尺寸运行的端到端回归，不能仅凭静态审计称已合格。

`initializeMetadata` 的 caseId 含初始 N；这是生成时身份，不应加点后重命名 case。另以 level/nodeCount 历史描述实际网格。

## 建议的预登记 family

这里只给字段草案，不是已获准的默认参数。策略显式 `version=2`，其余当前已验证数值元组与所有门不变。建议新增 `nodeFamily`：

```matlab
nodeFamily = struct( ...
  'generator','selected_base_index_pchip_v1', ...
  'cellFactors',[1 1;2 1;1 2;2 2], ...
  'ordering','node_product_then_registration_index', ...
  'componentwiseNondecreasing',true, ...
  'maximumAcceptedGrowthTransitions',2, ...
  'maximumTotalNodes',PRE_REGISTERED_VALUE, ...
  'maximumEstimatedPeakBytes',PRE_REGISTERED_VALUE, ...
  'resourceModel','FROZEN_MEASURED_MODEL_ID');
```

增长次数上界 2 是该四成员单调图的最长严格增长路径，不是新增的重网格次数额度。`maxRemeshes`、累计峰跳、步数和时间上界不因加点重置。资源值需根据 mesh 几何和实际硬件/LU证据确定；本文件不指定未经验证的内存常数。应明确失败时保持可信 checkpoint 并报资源/容量原因，不临时放宽门。

factor 乘的是区间数：Nx=(Nx0−1)fx+1，Ny=(Ny0−1)fy+1。321×161 初始预算对应：

| factor | 实际 N | 网格节点乘积 | Poisson 内点数 |
|---|---:|---:|---:|
| [1,1] | 321×161 | 51,681 | 50,721 |
| [2,1] | 641×161 | 103,201 | 101,601 |
| [1,2] | 321×321 | 103,041 | 101,761 |
| [2,2] | 641×321 | 205,761 | 203,841 |

按总节点数排名时，[1,2] 略先于 [2,1]；不要把注册行号当大小顺序。若明确只有 x 族失败而 y 原候选合格，先排除不能解决 x 容量的 y-only 成员；这个排除必须来自本轮实际轴筛查原因，不能仅由一个核心计数猜测。最终质量排序和稳定同分规则也应预登记。

## 零时刻生成及不可变性

1. 仍以 config 初始 Nx0/Ny0 完成现有 level 1 真实解析初值选轴及原门审计。`initialization.sourceBaseX/Y` 表示配置生成轴，`selectedBaseX/Y` 永远保存这次实际成功的初选轴。若 level1失败，本版本不自动改初始预算偷渡一次升级，应如实初始失败；更宽的初始注册需要另立合同。
2. 只在 level1 成功后，用 selectedBase 的归一化索引坐标生成全部成员。x 在正半轴做 PCHIP，再镜像得到精确双奇几何；y 在 [0,1] 索引上做 PCHIP。factor1必须直接复制原数组。factor2新数组的奇数位置直接逐值赋回旧节点，插入值只写偶数位置；端点、0、±anchor因此精确保留，anchor编号允许改变。
3. 每个成员直接来自同一 selectedBase，绝不从前一个候选/末态/另一个已插值成员递归生成。保存每个 member 的 factors、nodeCount、完整基轴、根节点嵌入索引、generation version、嵌入逐值检查和原 1D 质量审计。family1与selectedBase exact。跨机器若生成算术不能逐值重现，应拒绝依赖重生成的恢复；保存轴本身是权威，不能用重建近似值覆盖它。
4. 在第一接受步前冻结完整 family，记录构造发生在 original step/tau/t/remesh=0。初始化后不得追加成员或改变生成参数。几何不合格成员不能成为可进入 level；是使整个注册失败还是持久化固定 disabled mask，必须在版本中选一种，不得到末态临时解释。保守首版建议注册成员均需在 t0通过基轴质量门，否则拒绝这份 v2运行资格。
5. 保留的是 reference family 中旧节点，而非承诺每次重网格后的实际节点均为上一 actual轴的超集。现有 corePatch/equalizer会移动实际网格。若未来要求 actual网格也节点嵌套，那是另一种几何构造，本方案未验证。

## 当前级别、候选及提交

建议 memory.version=2；新增 `referenceFamily` 和 `currentLevelId`。`referenceFamily.rootX/Y` 可复用 initialization.selectedBase，`members(k).baseX/Y` 是所有级别的不可变数组。runtime `ops.baseX/Y` 仅作当前 member 的别名，其值必须与该 member 严格匹配；根不能随着别名变化被重定义。

某一已接受状态需换网格时，先查其 level，然后枚举所有 factors≥currentFactors 的成员；同级允许搬移，严格增长可只改变一个方向。`[2,1]→[1,2]` 必须拒绝，即使两者节点乘积接近或行号增加。每个 level 的几何都读取同一个实际 source.rho、source axes、source Dx/flow；target member 只提供预算和reference，不取代真实源。

为了防止前三个同级候选把增长路径永久挤掉，需要显式的二层预算：每个 eligible level 最多原 `maximumPairCandidates=3`，各轴搜索仍≤70，levels最多4；全部真实尝试最多12，按级别小节点数顺序，某一 level 无候选或全部原门拒绝后才进入下一 level。未来若需全局最多3，需要另立公平分配规则，不能假定现有截断自然给增长候选机会。每次真实 transfer 都从同一个已接受 source直接产生；不能将未提交的候选当下一次源。

以下内容在 trial 审计前必须与 source exact：config、case/run lineage、scale、三时钟、step、mass0、rhoRange0、所有原始 runtime references。只重算 originIndex/pinIndex 等当前轴索引。原 referencePeakPoints、adaptiveTarget/safety points 等不能因为新 N 增加而重置。所验证 scope仍是 transport_anchor/quadratic_peak；不要扩展到存有旧网格模板数组的其他规范。

v2实际事务的证据应明确保存：source/target levelId、source/target factors 与 N、sourceStep/三时钟/sourceRemeshCount、sameBox、registeredNodeCounts、matchingReferenceMembers、admittedLevelTransition，以及全部原 peak/mass/range/core/front/quality 审计。所有检查通过并重算完整 finite flow 后，才一次提交新 base、level、remeshCount+1、账本与新 RHS cache；新增 levelId 不是新物理 epoch。

`relativePeakJump` 仍是各自原生 Dx 的全域导数峰比较；mass仍用各自原积分权重、分母 max(abs(oldMass),eps)；累计为真实通过的所有迁移，包括同级与增长。禁止为新节点插值幅值投影、修平 C、放宽峰跳或换质量分母。PCHIP reference构造通过不代表 field transfer通过。

## 持久化和恢复

- config.grid.nx/ny 永远等于 initialN；config中的初始custom轴也保持原值。actualN属于原生状态和level，不能把配置改成actualN来骗旧检查。
- 策略 v2、memory v2、family生成版本和不可变根应严格类型/shape/value匹配。frame实际N来自其已接受epoch的 level，不接受“只要更大就可以”。
- history.mesh 添加两个 nodeCount和 levelId；window同样记录。每帧 rho/x/y 分别按自身N校验，并用 acceptedStep/remeshCount→transaction前缀确定 level；历史首帧必须 level1/configN，末帧和window与state exact。不删除或补写旧帧来对齐形状。
- native CP仍存完整原始logs、scales、cursor、rho/axes/base和metadata。现有签名已包含这些值；policy版本与checkpoint schema不是同一概念。如果新增额外顶层状态字段，签名字段也必须同步明确纳入；优先放入已经签名的完整metadata，避免漏签。是否保留CP schema4由最终合同审核决定，不能以外壳数字相同为理由绕过v2校验。
- `readCheckpoint` 先无LU校验config、family、嵌入根、各level质量/资源准入、ledger单调图、currentN/base及所有历史对应。建议调用相同版本的纯生成器从持久化 selectedBase 重算全部 family，要求 class/shape/value exact，而不是只查旧节点嵌入；若环境导致重算不同则拒绝恢复，不覆盖已存轴。v1继续完整旧规则，不注入family、不改变storedpolicy或签名。
- `restoreCheckpoint` 对已通过的v2直接按saved实际x/y/nx/ny作一次custom build，恢复匹配member base、原runtimeRefs，重算flow/cache，并严格核对terminal telemetry。不能先build初始level再发现N不同重build。
- v1 checkpoint不能被overrides升级为v2；现有restore白名单本就不允许策略覆盖。新v2必须从初始条件单独开始并全程使用冻结版本。新代码需要继续读取并按原算法续跑v1，而不是“读时迁移”。
- 原生CP是必须路径；result→CP桥接不是无人续接必需功能。若它暂不具备family-aware seed/base恢复与严格校验，应尽早显式拒绝v2，保留result可读性和nativeCP完整资格。不能把该桥接当前对尺度使用log重建的容差当native恢复exact证据。

## 内存、失败和无人运行

当前 `solve.m:45` 在 `advance` 返回后去除唯一旧 factor，再调用apply；拒绝时从旧精确A重建。这一所有权结构保留。升级比同 N 更要避免 trial、acceptedState、seed 或外部变量持有旧 `decomposition`；删除一个局部字段不证明LU已释放。

不能把节点乘积的约4倍误作 LU峰内存的约4倍：sparse fill、ordering和阵列副本可能非线性增长。只有数组密集存储随NxNy线性，family基轴自身仅O(Nx+Ny)。高阶 build 目前在 metricX/Y最终正性检查前就可能factor；所以新reference和实际候选的无LU预检应先检查所有现有1D质量、映射metric、数值可分辨间距，再做资源准入。资源准入不取代原数值门。

预登记上限需要针对目标硬件测量原矩阵组装、LU峰、flow/cache与source/trial并存、snapshot/checkpoint写入的峰和可用磁盘。动态空闲内存可以作为更保守的拒绝理由，不可以触发未登记的其他N。对真实OOM不要承诺catch后一定能恢复；分配前必须保有安全余量，定期原子native CP使意外进程退出后可继续。

候选无解、全部实际迁移拒绝、资源不足、程序依赖错误分别记原因；都不前移时钟或更新level。尚有更大注册成员且本轮预算未耗尽才自动继续枚举；整个有限族失败则保留可信状态并退出。用户不应被要求手工改网格来完成已登记动作，但最终资源有限的停止必须如实报告。

## 实现资格验证清单

1. 纯配置/无LU：v1正规化和旧CP/result字节/签名完全不变；v1→v2静默升级拒绝；未知生成器/非法factor/资源缺项拒绝。
2. 纯几何：factor1 exact；factor2旧节点嵌入exact；端点/0/±anchor和对称exact；每成员从同根、无递归；所有原质量与metric门，真实 step1134仅作冻结几何评分。
3. 小网格真实迁移：x-only/y-only/both，从同source独立，完整原质量/核心/峰/质量守恒/范围/finite flow门；v1同N原结果exact。拒绝后rho/scales/旧A作用/ledger/level保持exact。
4. lifecycle：预/后增长CP签名、每个成员单次restore、跨尺寸snapshots、所有history/window/原始references exact；连续与checkpoint分段唯一solve逐值一致；故意错member/base/N、降级、绕过一级、伪造旧节点、错epoch帧都拒绝。若允许从[1,1]直达[2,2]则它是合法边，不能把所有“跳级”误拒。
5. 自动触发：从t0到超过原step1134容量点，真实自然增长无需人工注入，保留完整可信前缀与失败；不能以手工调用事务代替自动策略验收。
6. 生产预算与数值资格：记录真实LU峰/时间和磁盘；再做从t0独立目标32/42/56或同物理终点空间加密控制。可运行更长不自动证明更准确，也不单凭梯度上升称强奇异性。

本设计交付只满足接口审计与合同明确化。具体factor族、所有成员初始质量可行性、资源上限和最终代码分派必须等待独立几何与硬件证据。
