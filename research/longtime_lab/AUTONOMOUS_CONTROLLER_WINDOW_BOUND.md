# 自适应控制器窗口的长期复杂度与有界采样建议

2026-09-09；只读审查及已执行的无 LU 微基准。没有修改 `+ipm`、配置或任何运行中的冻结源。当前 0.35 canonical 时间窗**没有独立于步长的记录数上限**；有限时间跨度不能代替内存上界。建议在下一版策略中显式注册有界观察池，同时保持每个接受步的即时核心数与安全门检查。

实际输出：[MATLAB 报告](../../result/longtime/20260909_controller_complexity_v1/run/report.json)、[输入及源码完整性审计](../../result/longtime/20260909_controller_complexity_v1/final_audit.json)、[研究采样模型](../../result/longtime/20260909_controller_complexity_v1/sampling_model/report.json)。MATLAB 83804 和 Python 97000 均退出 0。三份原生 CP 在 MATLAB 10 线程下严格读取；未 restore、建 LU、调用 flow 或推进 PDE。

## 当前成本来自哪里

[advance](../../+ipm/+evolve/advance.m) 在每个接受步调用 [planAutonomousMesh](../../+ipm/+evolve/planAutonomousMesh.m)。当前窗口有 m 行，提交过 R 次重网格：

| 操作 | 工作量及空间 |
|---|---|
| controllerTelemetry 追加/同一步替换 | 五个数组共六列 double；窗口数值本体为 48m 字节 |
| 逐步完整性检查 | 对六列拼接后检查 finite，另扫两列正核心；O(m) 扫描和临时数组 |
| 按 epoch 与 0.35 时间过滤 | 扫 m 行并复制五字段，O(m)；原状态/回滚状态的写时复制还会保留额外副本 |
| 完整历史峰跳累加核对 | 每步 sum 全部 R 项，O(R)；当前 maxRemeshes 是独立的有限预算 |
| 两轴预测 | 两次 degree-1 polyfit，两次 log/diff/secant 全量扫描，O(m) |
| checkpointSignature | runMetadata 的数值数组原样进入 signature；并未摘要或哈希窗口。构造主要是结构遍历/写时共享，不能误报成一次完整数值扫描 |
| 原生读写校验与存储 | isequaln、controller/window/history 校验及 MAT 序列化随数据量增长；payload 与 signature 逻辑上均持有完整窗口。makeCheckpoint 和 writeCheckpoint 又各调用一次 readCheckpoint |

每个接受步 dt≈常数时，m≈min(当前 epoch 的步数, 0.35/dt)。若一个同网格 epoch 的时间跨度始终小于 0.35，n 步累计窗口工作可达 O(n²)。固定正 minDt 和 maxSteps 会使**单个已配置运行**形式上有限，但该上界不是实用采样预算，且跨批 maxSteps 可合法增加。selectTimestep 在 dt<minDt 时停止而不推进；不能靠这道极晚的步长停止门控制内存。

以 minDt=1e-10 计，0.35/minDt 已约 35 亿行，窗口裸数值约 168 GB，尚未计 signature、临时数组、完整 history、PDE 场或 MATLAB 结构开销。实际 lifecycle 配置 minDt=1e-12，其形式上界还大两个数量级。这里是容量算式，未实际分配此规模。

## 实际窗口必须与稀疏输出日志区别

最新读取的长时分支是 fresh case 的 step4999：native tau=0.11339773068976232，物理局部时刻 0.0594080399189198，epoch=7。其完整日志只有 917 条，当前 epoch 只有 26 条。另一份实际原生重网格 CP 明确在 step4768/tau=0.11022997002750386 开始这个 epoch。因此若逐接受步记录，当前应有 **4999−4768+1=232 行**，不是 26 行。原生读取和这两个 CP 的路径均在报告中。

该分支本身尚无新生产 controller ledger；232 是由已验证原生接受步与 epoch 起点得到的每步采样计数，不能冒充实际已保存的 controller window。末步 native dt=1.3953666906311881e-5；若同网格、同 dt 维持满 0.35 native tau，估计约 25,084 行。该估计不是未来步长预测。

此 fresh case 的 native tau 与 parentEquivalentTau 不同；不能把“0.35 native tau”直接解释成“0.35 parent-equivalent tau”。这里没有改变时钟或推广 fresh case 为原 t=0 运行。

真正 from-zero lifecycle CP 的 step2 已有 remeshCount=2，实际 controller window 只有末行 1 条。该测试覆盖事务与恢复，却不会暴露长时间同 epoch 的窗口增长。最终从 t=0 的验收必须另含长期无重网格/高步数窗口场景。

## 小资源实测

同机 10 线程，固定合成平滑核心趋势；冻结调用原 controllerTelemetry、planAutonomousMesh 和 checkpointSignature。状态不请求网格设计，不包含 PDE ops。各点取 3–40 次中位数；小规模有 JIT/计时噪声，不能宣称整表严格单调。

| 行数 | 窗口裸数组 | telemetry | plan 含 telemetry | 合成 MAT 保存 |
|---:|---:|---:|---:|---:|
| 512 | 24 KiB | 0.034 ms | 0.136 ms | 2.75 ms |
| 1024 | 48 KiB | 0.092 ms | 0.141 ms | 4.55 ms |
| 4096 | 192 KiB | 0.151 ms | 0.370 ms | 9.73 ms |
| 16384 | 768 KiB | 0.496 ms | 1.624 ms | 19.34 ms |
| 65536 | 3 MiB | 2.537 ms | 4.982 ms | 78.76 ms |

65536 行的 signature 逻辑体积约 3.16 MB；签名构造约 0.082 ms，而值比较约 0.320 ms。合成 payload+signature 的压缩 MAT 为 2.47 MB；压缩率依赖数据，不能作为最坏内存上界。这些是控制器自身微基准，不能替代整步 PDE 性能测量；当前数百行还不是主计算瓶颈。

## 拟合项存在可消除的数学冗余

设 y_i=log(core_i)，t_i 严格递增。OLS 的斜率可写成全部两点割线斜率的凸组合，权重与 (t_j−t_i)² 成正比。每条两点割线又是其间相邻割线的凸组合。因此拟合斜率位于相邻割线斜率的最小值与最大值之间。

当前 decay=max(0,−OLS slope,所有相邻衰减割线)，在精确算术下等于 max(0,所有相邻衰减割线)。已做 400 个非均匀时间/随机序列拟合检查，越出凸包量为 0。未来可以省去这两次 polyfit；但浮点次序可能改变最低位，不能直接宣称与当前冻结算法 bitwise 相同，也不能静默修改旧 CP 的恢复语义。

## 建议下一版：32 个最近真样本 + 480 个时间桶

建议注册 512 行的观察池：最近 32 个真实接受步，加最多 480 个非空时间桶中的最后一个真实样本；按 acceptedStep 去重并排序。桶宽 h=0.35/(480−2)≈7.3222e-4，给两端部分桶及浮点边界留余量。桶原点为当前 remesh epoch 的实际 canonicalTime，原始时钟本身始终原样保存。空桶不填点、不插值、不构造虚假接受步。

只保留最近 512 步虽能限制内存，但在小 dt 下会丢失 0.35 时间跨度；时间桶解决这个问题。有效样本不足、时间覆盖不足或时间分辨率退化必须显式报告，不能假称完整趋势窗。

每个接受步仍执行现有 core<26、原硬核心门、安全门与事务门；这些不得采样。每次新 remesh 清空趋势池，用事务后同一真实时钟/step 的核心计数重新起始，禁止把换网格跳变作为收缩速率。

预测采用真实相邻接受步的 log-core 衰减割线，在每个桶另存两轴最大值及对应真实端点证据；每步的相邻速率计算 O(1)。只在桶的末样本也早于 t−0.35 时删除该桶。这样桶内部分过期的旧最大值可能多保留一小段时间，预测可更保守，不会漏掉完整窗口中仍有效的最大割线。不能先丢弃桶内点再从稀疏代表点重算割线，否则会漏掉短促收缩。

**512 仅指完整观察行；最大割线统计槽和端点证据也必须单独设固定预算。** 当前研究模型最多 480 个此类槽，每槽保存两轴常量大小证据，整体仍为固定空间；裸 double 数据约百 KiB量级而非“只有24 KiB”。生产实现应给 whos/序列化的实际硬预算，避免把额外摘要藏在行数之外。全量窗口/ledger/history 校验放在 checkpoint 入出边界与事务提交边界；逐步只验证新增真样本及固定池不变量。

研究参考模型 [mesh_controller_sampling_reference.py](mesh_controller_sampling_reference.py) 已执行 25,001 个合成接受步，包含步长从 5e-4 降至 1e-7、epoch 切换和 6 个序列化拆分点：最大完整观察行 510、最大桶数 479；所有保留点来自原始流，末样本始终精确，182 次衰减上包络对照全部保守，序列化恢复后每步完整内存与直跑完全相同，重复同一观察幂等。该模型不是生产 MATLAB 实现，也不是原生 CP/PDE 分支；不能据此宣称实际解或最终网格决策与旧算法一致。

## 进入从零时间最终验收前的最小补充

1. 新策略/控制器版本明确登记预算、桶边界、同 step 规则、桶内最大值失效规则。旧版本 CP 保持原算法语义，不在 restore 时隐式截断窗口。
2. MATLAB 原型覆盖逐渐变小 dt、超过预算十倍的同 epoch 记录、恰落桶/窗边界、收缩尖峰刚过期、同 step 重复、真正 remesh、初态、时间分辨率退化。固定池及额外证据槽均有硬上界。
3. 对同一接受观测流任意拆分保存/恢复，窗口、桶、最大割线证据、lastDecision 与重网格触发须逐值一致。checkpoint 保存全部必要状态，不靠外部 history 重建被丢弃的样本。
4. 保留每步即时核心/安全检查，在同一 from-zero 初值做生产 ipm.solve 直跑/split-run 生命周期验证。新采样可能更早触发网格，必须重新通过实际 transfer、累计峰跳及全部质量门；不允许用有界采样作为放宽数值门的理由。
5. 在注册长时端点记录整个运行的最大观察行数、统计槽数、控制器 bytes 和 CP 体积；不能只看某个重网格后为1的末窗口。全 history 仍随输出次数增长，transaction ledger 受 maxRemeshes 单独约束；限制窗口并不自动限制整个 solver 的所有历史存储。
