# remesh：网格迁移事务

将网格提案、场迁移与验收分开；候选失败不能改变已接受的场或网格。
导航：[包索引](../README.md) · [结构总览](../../STRUCTURE.md) ·
[第一象限架构记录](../../ARCHITECTURE_QUADRANT_LEVELSET_20260914.md)。

## 接口与流程

`adapt(rho,ops,flow,config,amplitudeScale)` 返回 `rhoNew/opsNew/info`；
`info.applied` 表示真正提交，只有接受后 `remeshCount` 才增加。

```text
第一象限：level-set 测量 -> 单组 x/y 轴 -> 质量门 -> 原生迁移 -> 接受/原样回退
旧全域：旧候选计划 -> 原生迁移 -> 验收 -> 减弱重试/原样回退
```

提案负责候选坐标、局部分辨率与单元比；迁移只构造候选场和算子；
验收负责改善幅度、值域、峰值变化及配套守恒指标。旧全域通用路径最多尝试
初始提案和八次减弱重试，可恢复的高阶网格/迁移拒绝会重试。
第一象限直接方法不减弱重试；任一门失败即返回原场/原轴。
程序错误或不相关异常向调用者传播。

## 方法与约束

`pchip/high_order/sixth_order` 使用各自配套插值与积分权重。
六阶迁移使用八点、七次局部插值；高阶无限制插值可能在不光滑场上过冲，
所以必须保留值域验收、拒绝和精确回退，不能将局部插值阶数等同于整体演化阶数。
重网格只覆盖 `mesh.build` 的网格参数，不重新解析配置或改变数值选择。
解析初始重网格的重新采样属于 `evolve.initialize`，不是迁移插值。
特征网格单元计数可以触发提案、参与候选设计/验收并支持分辨率硬停机；
在 schema-4 `exact_gauge_no_feedback_v1` 中，它们不得反馈修改任何动态比例率。

依赖 [mesh](../+mesh/INTERNAL.md) 构建候选算子，依赖
[diagnostics](../+diagnostics/INTERNAL.md) 测量特征；由
[evolve](../+evolve/INTERNAL.md) 调度，不负责求解循环或输出。
相关验证：[基线迁移](../../tests/+ipmtests/+baseline/remesh.m)、[四阶迁移](../../tests/+ipmtests/+fourth/remap.m)、
[六阶核心](../../tests/+ipmtests/+sixth/core.m)。

## 第一象限：一次 level-set 事务

显式第一象限模式不复用下述旧对称全域候选族。`directLevelSetAdapt` 从已接受场的
90% 壁面核心/前沿边界与纵向宽度一次反演一个正 `X` 轴和一个 `Y` 轴，
最多一对几何提案、一次原生迁移；质量、守恒、峰值和值域门失败即回退。
`evolve.remeshIfNeeded` 在配置的最高观察层为 90% 时，复用当前 `flow` 的
横向/纵向核心单元数作廉价预筛，与内部几何测量共用 `directLevelSetTargets`
的目标和 80% 阈值；其他观察层直接由本方法重测 90% level set。
解析初值的重网格也直接由本方法判断；象限不再使用旧综合 `safetyFactor`
触发。没有有限核心观测的无峰物理态不尝试反演。最终以重测的 level-set 几何为准。
现阶段保持同节点数，需更大分辨率时应在初始配置给出更大的实际正半轴 `nx`。
象限 `roundedAxis` 不再枚举左右单元分配：对细区外的两段长度 `L/R` 和细格宽 `h`，
先按 `log(1+L/h):log(1+R/h)` 分配固定的剩余单元预算，再投影到单元数与
最小格宽均可行的区间；只为这一组分配求解两条单调标量斜率。
若目标细格宽大于固定节点预算可容纳的宽度，先由左右跨度和连续平衡点
相邻的两个整数分配解析求出最大可行宽度，再一次限幅；不枚举轴、不重复尝试。
`roundedGeometricAxis` 对纵向宽度另求解一条单调标量方程。
这些标量二分不是二维场候选搜索；仍须通过整轴质量与原生迁移验收。
旧全域 `roundedAxis` 的遍历分配逻辑保持不变。
`roundedAxis(...,positiveOnly=true)` 的 `nodeCount` 是实际存储的 `[0,H]` 节点数；
默认旧全域调用的 `nodeCount` 仍是 `[-H,H]` 节点数。构造信息另记概念全域节点数，
调用者无需为象限先做 `2*nodeCount-1` 换算。

## 旧对称全域：自动网格计划与审计

plannedAxisPairs/roundedAxis/corePatchAxis/equalizedAxis只构造一维几何，数值依赖不进入research。
固定70个横轴和4个纵轴候选按最差无量纲质量余量排序，最多交回3对；保持实际单位、箱和节点数。
auditCandidate在真实transfer之后、计数提交之前验证配置/参考/时钟/初始不变量、场与轴配对、
完整梯度峰、质量/值域、原生90%核心、前沿及全轴质量。积分权重和所有派生门数必须实、有限且有效。
controllerTelemetry保存同epoch接受步观察和实际跳变账本；validateController是纯持久化一致性检查。
其校验不能重建未保存的历史迁移场，因此不能用metadata中的passed替代新的原生迁移计算。
version1不允许演化同轴事务，不支持节点数增长；其合同保持。

显式version2通过referenceAxisFamily冻结方向节点族与资源mask。planner区分真实来源尺寸和目标参考尺寸；
transfer用完整族/sourceLevelId/targetLevelId纯核对后仅override实际N和base，仍走唯一原插值通路。
audit分别核对同箱、注册尺寸/参考和合法单调跃迁，保留所有原数值门；differentN不伪称sameBoxAndNodeCount。

version3把冻结方向族扩展到三倍单轴与3×2/2×3组合；version4在每个方向成员内部按固定
层级扩大参数搜索。`hierarchicalAxisSchedule`和`hierarchicalAxisPairs`只补充候选生成，
每层仍应用同一质量门，整个尝试前缀与来源身份写入controller账本并由
`validateController`重读。默认和旧version1/2不注入这些字段，也不升级旧checkpoint。
version5从初始网格尺寸与显式节点预算构造更多方向级别；每一级仅在既有质量门内按
最大相邻比、质量裕量、稳定索引排序。其参考族和搜索上限都随冻结政策持久化，
不允许续算时改预算。同一来源、政策、特征区间下完全相同的参考轴试探在一次规划中缓存，
不跨接受步或续算使用；缓存结果的候选、报告与审计需逐值等价。几何候选不是实际迁移通过或长时容量保证。
显式 `densityVersion=1` 只属于 version5 新算例：原 70 个 X 候选把最宽的
48-cell 过渡替为 8-cell 过渡；第四个 Y 候选改用单侧 rounded-geometric
密度，使外区对数格宽斜率均摊。前三个 Y 候选和缺省旧政策完全保持。
`densityVersion=2` 复用同一 Y 函数和 70 个 X 候选，冻结
`search.xPadding=1.50`（v1 是 1.10），使后期同节点网格在质量门内保有更多 X
核心单元；版本值属于检查点冻结配置，不能在恢复时改变。
均摊密度的 X 初选圆滑参数 `[8,16,24,32,40]` 在补充阶段把最大值扩至 44，
使第一补充阶段仍有注册的 260 个唯一试探，完整阶段计数保持 `70/260/170`。
R2.5 未补齐此阶段，在首次触发补充搜索时产生账本异常；修复与回归见
[搜索阶段修复](../../research/longtime_lab/BALANCED_SEARCH_STAGE_REPAIR_20260914.md)。
新密度仍通过同一整轴质量门、核心/前沿预测门与实际原生迁移审计；
它不是放宽节点预算或突破质量门的快捷路径。

`planInitialAnalyticFallback` 是显式初始观察策略的纯计划入口，仅在原始t0使用。
省略失败账本时只返回启动前的解析预飞候选，不得冒充控制器事务；传入失败账本时
仍严格要求原始初选族耗尽或真实解析核心/前沿单独拒绝。
它保留真实粗格的source/clock/step/node字段，将新候选和观察证据另列；每阶段最多3个候选，
第二阶段索引仍局部编号，原失败账本不改写。`validateController` 同时核对观察规则、
使用次数、阶段索引和实际初选轴；保存的观察轴不属于求解referenceFamily。

`continuousVerticalCore` 是不参与原生决策的只读纵向观察器。二次峰幅值规范
提供 `omegaGaugeQuadraticPeakX` 时保持原评估点；其他幅值规范（包括 R2.2 的
外壁面密度窗）没有该专有点，改用 `trackFeatures` 已计算的连续
`trackedPeakX`。旧 R2.2 实验中观察器对新 C 每步产生
`ipm:ContinuousVerticalCoreInput`，异常被隔离、原生重网格仍运行；这项后备
修复只恢复观察数据，不改变网格请求、场迁移或 C 率。
