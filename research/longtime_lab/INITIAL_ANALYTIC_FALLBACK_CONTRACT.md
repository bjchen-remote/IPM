# 一次性原始初值观测 fallback：提案与已完成资格

本文件记录显式能力的契约。当前维护版生产没有自动解析观测 fallback；三个已运行的预选原生probe
只把已选网格通过customX/customY交给现有`ipm.solve`。没有改生产默认、原v2 policy或旧checkpoint。

## 已有证据及边界

研究入口`mesh_screen_initial_box_geometry`的`fixed_analytic_probe_v1`分支已经冻结并运行。
`mesh_prepare_preselected_initial_probes`又独立重算四个源和全部八个candidate的解析rho、配对七点Dx、
rho乘Dx、轴质量、实际核心/前沿及全部四级referenceAxisFamily，结果与保存数据逐值相等。

| 半宽H / Ymax | 观测轴与节点积 | 求解参考/候选 | 第一candidate实际core X/Y | 前沿 | 纯几何结果 |
|---|---|---|---|---:|---|
|32 / 16|527×247 = 130169|321×161|34.472881 / 36.784269|25.159955|3个pair通过|
|64 / 32|559×264 = 147576|321×161|34.379507 / 36.786509|24.963657|3个pair通过|
|128 / 64|589×280 = 164920|321×161|34.306021 / 36.780977|24.509062|2个pair通过|
|256 / 128|617×295 = 182015|321×161|无候选|—|X族0候选，保留容量拒绝|

观测节点上限600000与实际求解节点资源上限110000分别登记。观测网格不构造椭圆算子，不成为求解网格，
不写任何虚构的PDE历史。候选所属初始family的资源mask均为`[true,true,true,false]`；
被资源拒绝的双轴增长成员仍通过几何门。H128候选X相邻比1.07984356已接近原1.08上限，
因此此结果既不证明更大箱存在合格网格，也不证明这些箱能稳定推进到指定长时间。

冻结证据根为`result/verification/autonomous_runtime_20260909`（以下路径相对此根）：

- 原均匀粗观测及三轮边界：`large_initial_box_bounds_v1`。
- 固定解析观测：`large_initial_analytic_probe_v1`，前三个`screen/box_k_candidate_1.mat`被选中。
- 独立再审计与原生注册：`preselected_initial_native_preparation_v1/plan/audit.json`。
  noLU预飞实际exit0，CodeAnalyzer=0、profile无build/flow/transfer/solve/advance；冻源SHA未变。

## 建议的显式接口（尚未实现）

建议增加独立可选平铺选项`initialMeshObservationFallback`，正规化后归入`config.remesh`；
不要增加到当前必须逐值一致的`autonomousMesh.version=2`结构中。省略此选项时配置中不注入字段，
继续原v2初始化、原失败及原checkpoint解释。显式启用才接受如下固定注册：

```matlab
opts.initialMeshObservationFallback = struct( ...
    'version',1, 'enabled',true, ...
    'rule','primitive_k8_fixed_probe_v1');

% Proposed pure internal planner, not an existing callable production API:
[candidates,evidence] = ipm.remesh.planInitialAnalyticFallback( ...
    config, originalAxes, originalInitialAttempt);
```

版本1的规范化值必须固定：内区半宽4、spacing=.025、尾相邻比上限1.05、log斜率坡20单元、
观测上限600000、最多一次fallback、原每轴试验上限70/4和最终pair上限3。
本次验收求解policy仍为已注册v2/target32/32/cap110000；此次试验不扩大候选族、节点预算或质量阈值。
永久helper应读取已正规化的求解policy预算；观测规则自身的固定数字不作为可任意修改的调参接口。
扩大实际验收声明需要另行注册，不能仅从参数能被解析推断已有数值证据。

**首轮验收登记范围应明确限制**：原字符串初值`degenerate_primitive`、k=8、
`double_odd_omega`、transport anchor精确1、原配置均匀321×161轴、
箱集合`[-H,H]×[0,H/2]`且H属于`[32,64,128,256]`。
H256是已测可安全拒绝的输入，不是成功范围；目前没有区间内所有H的证据。
任意函数句柄初值、不同k/anchor、fresh晚态重启不属于这个命名规则。
321×161、target32、cap110000和离散H集合是当前试验登记范围，不应成为永久公共接口的硬编码。
其他节点预算或箱子不能复用此处的“已验证”标签；是否可合法试算应由明确的规则输入域、
原轴/资源/解析核心门决定，详见下述最小接线审查。数学上能生成轴也不等于具备动态误差资格。

`originalAxes`保存原配置生成的求解参考轴，fallback期间不变。
helper不接受调用者提供的rho、演化state或checkpoint；它只能从已校验的原始解析初值自行采样。
`originalInitialAttempt`保存原初选的完整有限族、每项实际解析重采样审核和拒绝原因，
调用方须证明这是唯一新运行的step=0、canonical/physical/normalized clocks=0、remesh=0，
尚未提交任何接受态或t0历史。零时刻上下文不能由一个孤立布尔`isInitial=true`替代。

## 调用顺序与接受事务

1. 完整运行原初选。已有候选真正通过原初始化门时，原路径及浮点顺序保持，不调用fallback。
2. 仅当有限族无提案，或候选在解析重采样后全因core/front不达标时，允许一次固定解析观测。
   配置不兼容、非有限初值、参考族/资源无效、签名错误、LU/flow失败等不属于观测误差补救范围。
   对原轴族质量导致无提案，可记录为有限族耗尽；观测轴自身必须独立先过全部原质量门。
3. 从同一原解析初值生成固定观测rho和配对七点Dx；精确解析导数只作旁路误差诊断，
   不替代实际planner输入。固定reference仍是原配置的321×161轴，source观测节点可以不同。
4. 原planner排序产生最多3个pair。逐项回到实际321×161候选重新解析rho/Dx，并检查
   全部原质量、精确±1、actual core≥31/31、front≥20，以及全部参考族成员几何与冻结资源mask。
   第一项真正通过原native初始化流场/安全门后才作为selected base提交；没有反复改变观测spacing的循环。
5. 所有尝试失败则保存原拒绝与fallback拒绝，初始化报错；不伪造可恢复checkpoint、时间历史或成功t0。
   持有旧LU的实现必须在下一候选build之前清除所有alias，不让观测辅助对象携带LU。

这里的候选解析重采样是原始t0初选，不是重网格迁移；不能额外套用旧粗采样到新解析场的
守恒校正来改变原始datum。一旦已提交t0并进入演化，所有后续网格变化必须走原生transfer/audit，
绝不回调初始解析公式取代演化场。

## 持久化与恢复

新可选config字段只接受严格规范化值；缺字段旧schema4和所有旧CP保持exact，不静默补写。
须同步新增output/config/controller验证，再考虑启用；本提案不修改现有验证器。
新运行的initialization应附独立`observationFallback`证据：规则与完整固定controls、
original axes、原拒绝列表、是否实际触发、固定观测轴/质量/节点积、选中pair原序号、
解析candidate实际门和selected base family。现有sourceBase与selectedBase的原义仍保留。
保存证据不能把辅助观测nodes写入history.nodeCountX/Y、meshLevelId或solver资源账本。

恢复只验证已保存config/初始证据及真正family，再按当前saved axes恢复；
不得重新初选、重新执行fallback或从当前rho推测原始初值。故CP白名单不允许改变此选项，
也不允许把旧CP静默升级为fallback算例。

## 原生预选probe已注册，尚非集成验证

稳定研究入口：

```matlab
report = mesh_run_preselected_initial_probe(registrationFile, label);
```

冻结目录`preselected_initial_native_preparation_v1/source`，registration为同目录根下
`plan/registration.mat`。独立launch_H32.m、launch_H64.m、launch_H128.m
分别只调用一次原`ipm.solve`：默认10线程、canonical终点4e-5、maxDt1e-5、maxSteps4，
严格检查到达终点、原生result/CP/config/history、源/选中/终态轴exact、t0解析rho和所有零时钟、
qualified keep序号1、终态core≥31/31、front≥20、remesh0。
最大步数退出不等于成功。t0证据来自终态原生CP保存的初始化账本与配对快照，
不宣称另存了不存在的零步CP。

三项已由root串行完成四步，均达到canonical 4e-5并通过原生资格；H32在
`preselected_native_review_and_remaining_v2/H32_actual_review/report.json`只读重审，
H64/H128在原准备目录的`plan/runs/H*_preselected_initial/report.json`。
终态core分别为34.472526/36.784317、34.379177/36.786571、34.305711/36.781038。
H32求解早已成功，原research auditor错误访问空transactions字段导致审计退出；
修复仅在研究审计器，重审已有结果，没有重跑H32 PDE，原失败保留。
这些是外部预选轴的原生资格，不能代替下面新算法的集成运行。
将来真正集成还至少需要：无选项旧v2结果/CP bitwise回归、无选项H8/H16原成功初选不触发、
H32/64/128失败后一次补救、H256安全拒绝、错误初值/箱/clock/重复调用负例、
新CP逐值split-replay，以及单LU生命周期。四步原生预选资格不替代这些集成回归。

## 现有代码的最小接线审查（只读，2026-09-09）

当前`initialize`先构造原初态与流场，调用`planAutonomousMesh(state,true)`，
然后释放原poisson再调用`applyAutonomousMesh`，成功返回后`solve`才创建log并记录t0。
因此可以在`initialize`的两处现有报错前接一次fallback：纯计划无候选，或实际初选全失败且原因确属
可补救的解析分辨率。首版保留原始flow及成功路径顺序；“在原始LU之前做纯初选”是另一个效率重构，
不是本次能力必要条件。原初始flow本身异常时仍如实失败，不以fallback吞掉它。

### 必须修改及新增的函数

| 文件/函数 | 最小职责 |
|---|---|
|`+ipm/+config/schema.m`|登记无默认可选`initialMeshObservationFallback`，独立kind，归remesh。|
|新`+ipm/+config/initialMeshObservationPolicy.m`|纯正规化固定规则/一次预算，校验enabled自主v2及原k8/anchor1兼容性；不构造场，不改原autonomousMeshPolicy。|
|`+ipm/+config/resolve.m`|在autonomousMesh正规化之后，仅字段存在时调用新正规化器。|
|`+ipm/+output/validateV2.m`|新增kind验证，调用同正规化器并严格class/shape/value比较；schema≤3不接受新字段；无字段不补值。|
|新`+ipm/+remesh/planInitialAnalyticFallback.m`|固定观测轴生成、原initialDensity/pairedDx、现plannedAxisPairs、候选解析几何及全family审计；无build/flow/transfer。可将固定轴构造提为该函数local helper，避免额外公共API。|
|`+ipm/+evolve/initialize.m`|唯一一次fallback分派、原失败保存、真实零时刻证明、原LU别名释放；保持config和真实来源state不变，将fallback候选交回现有apply。|
|`+ipm/+evolve/applyAutonomousMesh.m`|仅在带已核对fallback证据的initial plan成功时把证据附入initialization；原initial_audit/native flow/family提交不另建数值路径。|
|`+ipm/+remesh/validateController.m`|新option与reserved initialization证据相互一致、类型和值严格验证；缺/禁用option却出现fallback账本拒绝，旧无字段分支不变。|

上述是六个现有文件加两个pure helper的最小建议，不是已提交的补丁。
`planAutonomousMesh`本身可以不改：初始化分派保留原plan的真实来源投影，只换候选/轴报告并附证据。
`controllerTelemetry`也不必认识观测网格，它始终只接真实state。
`solve`、`advance`、`transfer`、`auditCandidate`、`record`没有必要变化；原始t0选择仍不增加remeshCount。

### 三个不可绕过的实际合同

1. **来源投影不能换成观测网格。** `applyAutonomousMesh`的v2入口逐值核对
   plan.sourceStep/三个时钟/remesh/N/level和plan.coreCells与真实state。
   fallback的观测source虽然使用同一解析datum，但节点数和核心计数不同，不能冒充这里的state。
   保留原plan的来源字段，观测features仅进axisReport/initialization.observationFallback。
2. **两个阶段不能拼成六个主attempts。** `validateController`要求selectedCandidateIndex≤3，
   且initialization.attempts长度恰等于选中index，只有最后一项通过。
   主字段继续描述最终选中的局部阶段1..3；原初选的失败及完整报告放入独立originalAttempt。
   额外记录selectedPhase=`original`或`fixed_observation`，并严格核对触发条件和每阶段预算。
3. **不能只按当前failure字符串猜原因。** `applyAutonomousMesh/initial_audit`将失败聚合成
   `initial_axis_or_analytic_resolution_gate`，这一名字不能证明仅core/front失败。
   最小实现可读取现有audit的qualityPassed、exactZeroTime、coreCells、leftFrontCells、轴质量及errorIdentifier，
   同时独立核对±anchor；所有项只因core/front不足才进入这一分支。
   有异常、post_transfer_flow_resolution_or_finiteness或无有效audit时不要自动“修复”。
   如需更细reason字段，只在显式新能力分支中增加，避免改变旧CP的audit逐值结构。

所有原候选在apply前的纯解析筛查失败时，可避免这些候选的新LU，但这是新选项启用分支内的优化。
若要把未尝试native的pure拒绝记入账本，应明确audit层级，不能填进原native attempts并伪称运行了flow。
已有实现会在原family构造不通过时返回异常；这仍应是独立family拒绝，不因观测核心改善而放行。

### 原生保存/读取原则及无需修改的路径

`makeCheckpoint`已经完整打包config和runMetadata；`checkpointSignature`递归包含两者。
新增规范化字段和初始化证据因此自然受现有签名覆盖，无需改变签名算法或checkpoint schema。
不要将观测rho/Dx/整套临时数组塞入每份CP：持久化固定规则、观测轴、质量/特征和原始尝试的有界证据即可；
完整观测场可保存为独立研究artifact，不能让它变成恢复所需外部依赖。

`readCheckpoint`本来就会正规化并逐值核对config，随后调用validateController；
`restoreCheckpoint`也调用同validator。把新证据的验证封在validator中即可，两个reader通常无需改代码。
恢复的已有白名单未包含新选项，故其覆盖会自然被拒绝，不能加到白名单。
恢复继续按saved.x/y一次build，复制saved.rho/base/scale/mass0/rhoRange0/runtimeRefs/log/cursor，
重算原flow/RHS；它没有调用initialize、initialDensity或planner的理由。
重新生成**一维固定观测轴**来检查metadata是否与规则匹配，可以作为pure验证；
不要在reader里再次生成解析rho、重新评估候选或覆盖保存的任何数据。

`checkpointFromResult`已经对自主v2明确拒绝；维持原拒绝，不为了fallback增添第二种族谱重建。
也无需让`restoreRuntimeReferences`加入观测参考值：观测量不属于规范reference，
真实初始references应仍由最终选中的native轴/初值初始化并原样持久化。

### 永久规则与一次实验限制的分层

- **规则必要条件**：只用于原始t0；显式新选项且自主v2；命名规则确实对应k8 primitive、anchor1、原对称性；
  固定一次观测和有界候选；全部原质量/核心/前沿/资源门；不能接收演化场或CP。
- **当前轴公式的输入域**：对称X端点H>4、Ymax>4，尾fzero有合法括号、内区及±1精确保留，
  实际观测节点积≤600000，观测轴自身质量通过。Ymax=H/2和H是2的幂不是公式要求。
  小箱可继续原初选；若未来需要对端点≤4执行同类观测，应注册清楚裁切公式，不能悄悄改本规则。
- **可复用求解预算输入**：从config/immutable originalAxes读取Nx/Ny和policy target/cap，
  让原planner和全family门决定容量；321×161、target32、cap110000不应写死在永久helper。
  预算为其他已允许值时可返回“合法试算/无资格证据”或正常容量拒绝，不能自动贴本研究的已验证标签。
- **当前实际验收矩阵**：只有H32/64/128的321×161/target32/cap110000解析几何通过，H256失败；
  原生probe待运行，其他组合更没有动态误差资格。报告应单列`matchesValidationRegistration`和
  `runtimeGatesPassed`，不要用一个supported布尔混合数学输入合法与实验已通过。

### 必需的恢复验证点

1. 旧无option的配置正规化和实际旧v2 CP严格读取、恢复及split结果逐值不变；不能注入空fallback结构。
2. 新选项enabled但原初选成功时，明确used=false/attemptCount=0，原native场、轴、尺度与无选项相同；
   metadata显式多出新证据是预期，不把整体结构称bitwise相等。
3. 新fallback成功后，测试storeSnapshots=true和false的原生保存/重读；结果与CP严格配对，
   source base仍原配置轴，selected base/family为真实选中候选，观测N从未写入solver历史。
4. 新CP→restore→继续的split与不中断的同配置运行逐值相同：rho、ops轴/base、scale、flow、
   reference values、history、controller、调度cursor。profiler必须显示恢复仅一次build/flow，
   没有initialize、initialDensity、planInitialAnalyticFallback或plannedAxisPairs调用。
5. 对新policy缺字段/错误class、used与计数不符、错误规则/原始拒绝/观测节点/选中phase、
   伪造obsN为levelN、修改option的resume override及v2 result桥接分别拒绝；
   各拒绝发生在LU前，旧CP不受这些新字段规则影响。

本节审查后，root授权将八处变更仅实现于
`result/verification/autonomous_runtime_20260909/initial_observation_fallback_candidate_v1/source`。
维护`+ipm`及全部旧冻源未改。独立candidate已通过noLU配置、四箱几何逐值、34个枚举反例、
used/unused合成账本及旧实际v1结果/v2 CP读取；新算法的native/default/split脚本已登记，由root排队，
尚不能把预选probe通过写成该算法的native集成通过。具体清单见候选目录`IMPLEMENTATION_REVIEW.md`。
