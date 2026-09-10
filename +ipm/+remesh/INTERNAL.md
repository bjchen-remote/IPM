# remesh：网格迁移事务

将网格提案、场迁移与验收分开；候选失败不能改变已接受的场或网格。
总体契约见 [STRUCTURE.md](../../STRUCTURE.md)。

## 接口与流程

`adapt(rho,ops,flow,config,amplitudeScale)` 返回 `rhoNew/opsNew/info`；
`info.applied` 表示真正提交，只有接受后 `remeshCount` 才增加。

```text
propose（axis） -> transfer（interpolate + mesh.build）
               -> validate -> 接受 / 减弱变形并重试 / 原样回退
```

提案负责候选坐标、局部分辨率与单元比；迁移只构造候选场和算子；
验收负责改善幅度、值域、峰值变化及配套守恒指标。最多尝试初始提案和八次减弱重试。
可恢复的高阶网格/迁移拒绝会重试；程序错误或不相关异常向调用者传播。

## 方法与约束

`pchip/high_order/sixth_order` 使用各自配套插值与积分权重。
六阶迁移使用八点、七次局部插值；高阶无限制插值可能在不光滑场上过冲，
所以必须保留值域验收、拒绝和精确回退，不能将局部插值阶数等同于整体演化阶数。
重网格只覆盖 `mesh.build` 的网格参数，不重新解析配置或改变数值选择。
解析初始重网格的重新采样属于 `evolve.initialize`，不是迁移插值。
特征网格单元计数可以触发提案、参与候选设计/验收并支持分辨率硬停机；
在 schema-4 `exact_gauge_no_feedback_v1` 中，它们不得反馈修改任何动态比例率。

依赖 `mesh` 构建候选算子，依赖 `diagnostics` 测量特征；不负责求解循环或输出。
相关验证：[基线迁移](../../tests/+ipmtests/+baseline/remesh.m)、[四阶迁移](../../tests/+ipmtests/+fourth/remap.m)、
[六阶核心](../../tests/+ipmtests/+sixth/core.m)。

## 自动网格的纯计划和实际审计

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

`planInitialAnalyticFallback` 是显式初始观察策略的纯计划入口，仅在原始t0使用。
它保留真实粗格的source/clock/step/node字段，将新候选和观察证据另列；每阶段最多3个候选，
第二阶段索引仍局部编号，原失败账本不改写。`validateController` 同时核对观察规则、
使用次数、阶段索引和实际初选轴；保存的观察轴不属于求解referenceFamily。
