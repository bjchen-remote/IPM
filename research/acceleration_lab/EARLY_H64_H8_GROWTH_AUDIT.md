# H64 与 H8 早期自然增点：真实事务只读审计

结论：H64 较早用完原两次增点，是 **Y 预测触发与已登记有限候选容量共同作用** 的记录事实；现有证据不足分配二者的因果比例。H64 并没有比 H8 更频繁的早期迁移。不能把大箱的早增点归因于独有的预测噪声，也不能把有限候选未产出当作所有可能网格的容量不够。

各自冻结 runtime 默认10线程严格读取真实 CP：H64 step1269 / τ2.80012407015，H8 step4545 / τ8。两者所有 history trusted，均无 snapshots，history 分别281/801行，非每接受步保存。实际 LU/PDE 次数均为0，原输入 SHA 未改。这里使用的 H64 case 为901b551b，不是早期短探针。

来源、逐事件完整表、原证据：[growth_comparison.json](/Users/bjchen/Library/CloudStorage/OneDrive-北京大学/Research/IPM/ipm_structured/result/longtime/20260908_campaign_v1/acceleration_lab/early_H64_H8_growth_872fe5db5b5c44789d9168710801175a/growth_comparison.json)、[all_transactions.csv](/Users/bjchen/Library/CloudStorage/OneDrive-北京大学/Research/IPM/ipm_structured/result/longtime/20260908_campaign_v1/acceleration_lab/early_H64_H8_growth_872fe5db5b5c44789d9168710801175a/all_transactions.csv)、[H64 evidence](/Users/bjchen/Library/CloudStorage/OneDrive-北京大学/Research/IPM/ipm_structured/result/longtime/20260908_campaign_v1/acceleration_lab/early_H64_H8_growth_872fe5db5b5c44789d9168710801175a/H64/evidence.json)、[H8 evidence](/Users/bjchen/Library/CloudStorage/OneDrive-北京大学/Research/IPM/ipm_structured/result/longtime/20260908_campaign_v1/acceleration_lab/early_H64_H8_growth_872fe5db5b5c44789d9168710801175a/H8/evidence.json)。

## 四次实际增点

| 箱 | step / τ / physical t | level | 触发 | source core X / Y | predicted X / Y | 新 core X / Y |
|---|---|---|---|---|---|---|
| H64 | 319 / 0.804260918 / 0.682763659 | 1→2 | Y forecast | 31.47017/31.16996 | 27.91025/1.41010 | 34.74558/37.77500 |
| H64 | 555 / 1.463991424 / 1.073028706 | 2→4 | Y forecast | 32.62200/31.53945 | 28.58955/2.87275 | 35.40499/37.33760 |
| H8 | 1135 / 2.691375458 / 1.613866760 | 1→2 | X current<26 | 25.95699/32.33849 | 22.49762/2.40749 | 35.16712/36.03221 |
| H8 | 2477 / 5.594870008 / 2.109245773 | 2→4 | X current<26 | 25.99224/31.87619 | 22.72173/29.41754 | 35.16633/35.26913 |

两箱的 autonomous policy 逐字段完全相同：目标32、当前触发26、预测buffer22、review长度0.2、trend窗口0.35、每族最多3对、相同几何/质量/迁移门、节点cap210000。level1/2/3/4分别321×161、641×161、321×321、641×321。两个 growth 的实际完整峰值跳变：H64为9.35091e-5、8.81514e-5；H8为1.17685e-4、8.09269e-5；原质量门均通过。

## 预测敏感性可以直接测到，selector 原因仍不完整

| H64增长 | Y最大负log割线 | Y窗口OLS衰减 | 倍数 | 割线区间Δτ | 实际Y计数单步变化 | 同状态仅OLS预测Y |
|---|---:|---:|---:|---:|---|---:|
| step319 | 15.4789664 | 0.8712610 | 17.766 | 0.002859132 | 32.58041420→31.16996294 (-4.329%) | 26.185418 |
| step555 | 11.9798416 | 0.3640383 | 32.908 | 0.002726352 | 32.58657875→31.53945310 (-3.213%) | 29.324742 |

以上两组计数是保存的 step318/554 原 history 与 step319/555 **迁移前** controllerDecision 的真实计数组合。重算负log割线与记录值相符，均在同一迁移epoch，因而不是把迁移跳变当作原时间导数。只使用OLS的同状态预测分别26.1854、29.3247，高于22；这是明确标记的单状态公式对照，不是重跑另一控制器或证明可以安全跳过迁移。

记录没有保存319/555迁移前峰列和0.9交点模板；后续 history 已经过新网格，不能拿新峰列与旧列比较后宣布该单步下降由selector切换造成。H8第一次增长所用 Y 最大割线发生1087→1088；同epoch的保存区间1087→1091中峰列1→0.9964929462且Y宽下降，支持采样区间的关联，但不能精确归因1088。H8第二次增长Y最大/OLS仅1.257倍，触发仍是X当前core<26。

纵向原观测取离散墙峰列的Y截面，其列切换可以放大宽度计数的非平滑性；最大负割线会保留该尖峰。这里没有源rho时，不能以这个机制的一般可能性替代逐事件证据。

## 早期频率与有限候选容量

预登记共同早期窗口为τ∈[0,1.6]：H64有9次迁移（2次增点），H8有11次（0次增点）；所有这些早期迁移均为Y forecast，间隔中位数分别0.14140和0.09197。H8也有比H64更大的Y负割线（最高40.9937）。因此“大箱预测噪声更大”不能解释为何H64更早增点。

四次增长都接受 combined candidateIndex=1。按冻结 planner 的节点乘积/登记序号排序，第一次增长前的较小 eligible levels 为[1,3]，第二次为[2]；这些层没有返回可实际尝试的候选，才轮到增长层的第一个候选。此项是保存的 first-candidate 证据结合冻结确定排序得到的推论，非补造 axisReport。

已提交 audit 仅保存 controllerDecision、accepted candidate quality 和 attemptSummary，没有保存各层 axisReport / xTrials / yTrials。因此不能确定“某个x邻格比”或“某个y宽门”是较小族失败的唯一原因，也不能重建所有未返回候选。第一增长连保持Nx321而仅增Y的level3也没有返回候选，与横向可用候选限制相符；第二增长发生在level2、改增Y，与纵向限制相符，但都不是完整逐轴失败证明。

初选本身也不同：H64使用一次已登记固定观察fallback；H8没有启用该项。初始core H64为34.3795/36.7865，H8为34.8595/39.3226。仅对真实保存的初始selected axes统计固定区域节点：

| 固定初始区域 | H64节点 | H8节点 |
|---|---:|---:|
| positiveX0to2 | 94 | 104 |
| positiveX2to8 | 34 | 57 |
| positiveXBeyond8 | 33 | 0 |
| Y0to1 | 38 | 41 |
| Y1to4 | 29 | 120 |
| YBeyond4 | 94 | 0 |

这些是包含原端点的固定区域节点计数，不是core cell门。H64为远场分配了33个正X外区节点和94个Y>4节点；相同161正半轴节点预算与平滑过渡约束下，后续集中能力更有限是合理解释，但仍包含初选方法、非局部箱影响与随后轨迹差，不能由此作纯边界误差归因。初始c_l H64=.4817725984，H8=.4349185431，已不同；相同canonical τ不代表相同physical时刻。

## 下一项最小可证伪诊断

先增加研究边界输出：每次requested source的原rho/axes/Dx、峰列坐标/索引及Y取样列、两0.9 root模板；保留完整各层axisReport、实际top3索引、fit/max割线的端点计数与峰列、单步Δτ。继续使用原controller，不修改门或接受逻辑。这样的缓存可无LU重算连续Hermite同峰宽，与原离散计数在同一真实场上比较。

若必须复现现有事件：真实H64 step318 CP恰在第一次增长前，冻结native advance一步到319，尚未apply前缓存rho和完整meshPlan；第二次从真实step532 CP按原native走23步至555，在前置554和555缓存。这是一旦另行授权LU/PDE窗口后才执行的最小范围；本审计未执行。不能只拿step318或532替代缺失的触发源场并声称重现失败轴门。

本报告不改变正在运行的H64或H8，不放宽旧网格门，也不认定强奇异性成立。
