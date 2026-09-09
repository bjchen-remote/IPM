# 自动方向增点 v2：运行合同与验证

最终目标仍是从原始物理 t=0 自行推进到长时并验证强奇异性；有限族和短程测试不是最终资格。

显式 policy.version=2 才启用新合同，省略版本仍为1，旧检查点不能静默升级。
所有C计算、PDE、时间步规则、质量/core/front/峰跳/质量守恒/范围门保持。

## 预先登记的网格族

零时刻在解析初选轴上冻结 referenceFamily。区间倍数为 [1 1;2 1;1 2;2 2]，
各成员都从同一初选base通过索引PCHIP生成并精确保留旧节点、端点、0和锚点。
全部级别先过轴质量；显式 maximumTotalNodes 决定冻结的资源准入mask。
节点上限不是已测LU峰内存模型，也不保证任意长时可计算。

配置Nx/Ny永远指初始预算；实际轴和base必须匹配当前族成员。实际重分布轴可以移动，
不能将参考族嵌套解释为所有实际轴嵌套。当前级别先试，再按节点乘积/注册号排序，
各维只增不减，每级独立最多3提案（含qualified keep）；不同N无keep，演化不提交同轴事务。
所有候选都直接使用同一个真实来源场。

## 原子提交和恢复

solve在advance返回后释放旧LU；apply逐候选构建、原插值/守恒、实际审计和flow。
全部通过后才同步提交base、currentLevelId、remeshCount、逐事务controllerDecision/attemptSummary与新RHS cache。
审计分别说明sameBox、注册尺寸、参考匹配和合法跃迁；不同N时sameBoxAndNodeCount为false。
source的尺度、时钟、步数、初始不变量和runtime references保持。

v2 history.common另存nodeCountX/nodeCountY/meshLevelId，window保存nodeCount/levelId。
原生CP在LU前核对完整参考族、每个历史/快照所属级别、真实触发因果和累计账本；恢复按保存实际轴只构造一次。
刚提交事务后的lastDecision仍是该事务的source决策；初始化前的core与初选后core分别有记录。
这些是一致性校验，不重放未保存的历史迁移场。result→CP尚不支持v2，明确报
ipm:CheckpointResultVariableNodeFamily；原生检查点是已实现的续跑入口。

## 实际验证记录

证据根为 result/verification/autonomous_runtime_20260909。

- node_family_geometry_v2_v3：配置10正/26反、几何6正/7反、8个实际状态比较通过；v1与冻结旧版逐值一致。
- controller_v2_persistence_v3：47反例、9正边界、44个v1反例与旧实际签名通过。
- campaign_v1/acceleration_lab/variable_n_audit_ddc207b7ee144aabbce8c58b5283f3f9：50项纯迁移合同通过；明确为无LU的state-view测试。
- v2_smoke_and_suites_v1：真实零时刻、4步CP分段场/完整日志逐值一致，前4步物理字段与冻结v1逐值相同。all及equivalence通过；6个物理运行、3项重网格核逐值保持。
- from_zero_tau4_v2_cap110000_v1：已从原始t=0连续完成tau4，step1653、physical t1.9140859724939365、17次实际迁移，step1135/tau2.6913754578032374自然由321×161增长至641×161；cap110000、箱[-8,8]×[0,4]及全部门全程冻结。末core32.85919162/34.66892314、G4.890436581、累计跳.001148393161，wall686.82秒、maxRSS3.002GB。
- natural_growth_restart_v1：20份真实CP严格可读，13份共同v1前缀全部既有字段/history逐位相同；真实CP1102复跑到1177跨自然增长，场/尺度/参考/完整日志/跨尺寸快照/controller及新native CP逐位一致；6个重签语义损坏均拒绝。
- from_zero_tau8_v2_cap210000_v1：另一原始t=0、cap210000冻结运行已连续完成tau8，step4545、25次remesh、两次自然增长，末641×321、core27.38530409/29.72754570、累计jump.00179034457402；physical t2.2412067821954538，full physicalGradInf26.50250791。wall3763.11秒、maxRSS3.653GB、0swaps；仍需域/空间/时间误差资格。

初始三箱研究还发现H32粗网格的几何padding导致解析重采样后x核心27.63<31，正确拒绝。
仅零时刻有界解析重测研究在H32第二轮使核心34.6654/36.7910、前沿25.1172通过；尚未加入v2，见INITIAL_ANALYTIC_ITERATION.md。
控制器时间窗仍没有固定样本数上限，固定容量聚合窗尚未进入生产。
