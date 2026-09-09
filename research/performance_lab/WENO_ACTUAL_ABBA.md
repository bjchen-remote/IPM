# 最新895×386冻结缓存的WENO实际ABBA

结论：block16、block32在默认10线程下均逐位一致，但完整RHS/一步没有达到预注册的稳定5%收益，因此不进入匹配物理时间的下一轨迹对照，不改生产WENO。

输入是最新epoch长盒对照的crop_natural原生step81，case `20260908T134443474_dynamic_isotropic_895x386_tp6fdc3a52_fed0_4876_aa77_e1c92ff19ee4`。其absolute physical time约1.9590400984028897，fresh local elapsed .004512961717177134。输入由上一三盒完整trusted/native-pair/history/lineage审计通过，不把它伪作旧full-box case。

运行冻结目录为 `result/verification/performance_lab_20260908/weno_actual_step81_v1_source`，保存117份MATLAB文件及哈希。98份+ipm沿用此前已实测代码；两个研究namespace分别仅复制9个数值调用链函数以重定向WENO（advance、stepSsprk54、stepRk、flow、rhs、isotropicGauge、assembleRhs、transport、weno5FluxDerivative）。除WENO独立face分块及调用名路由外，表达式/系数/规范/诊断不变。研究namespace没有solve或restore，LU及其它功能仍调用原维护实现，不通过切换path或修改已运行源码测速度。

`ipm_perflab_weno_actual_abba` 默认10线程，严格restore一次12.18079秒，此后所有样本共享这一份ops/decomposition与接受态RHS cache。每个完整advance都回到相同冻结输入，输出不提交，不记录新history，不保存CP。因配置禁止adaptive remesh，样本不构建其它LU。生成器 `ipm_perflab_build_weno_probe.py` 保存实际源码来源和文件哈希。

预注册block16、32；四个真实WENO调用是rho_x、rho_y、unit_x、unit_y，保留维护全局LF、epsilon、边界外推、metric和计算步距。每个操作先A/B各warm一次，然后4个ABBA或BAAB块，每块A/B各两次。固定独立RandStream seed20260908预先决定次序。两个kernel、RHS与step的warm输出及全部计时输出都与baseline reference比较；RHS包括整个flow，advance包括除共享ops对象以外所有state字段（rho、scale、flow、rhsCache、config等）。没有用容差放过roundoff。

代码先在实际tiny fresh checkpoint上以两块/操作通过全链bitwise，所有19份研究MATLAB函数Code Analyzer为0，生成器Python AST通过。实际大网格每候选每操作16个样本，共192个计时调用；另24个warm调用。计时包含两侧相同的浅层返回结构打包，精确比较在计时外。对应完整advance共36次冻结重复，不是36步连续轨迹。

预注册准入：所有操作逐位一致，且whole RHS或whole advance的4个成对块速度比均≥1.05。不得只选择最好一次或某个内核作为整体收益。全部原始order/seconds/bitwise矩阵保存在每操作JSON与benchmark_report.mat/json。

| block | 操作 | 原核中位ms | 候选中位ms | 总中位速度比 | 4块速度比范围 |
|---:|:---|---:|---:|---:|---:|
| 16 | rho_x | 13.7408 | 24.4093 | 0.562931 | 0.518659–0.614572 |
| 16 | rho_y | 12.4886 | 10.8012 | 1.15622 | 0.60471–1.27607 |
| 16 | unit_x | 14.1055 | 10.5513 | 1.33685 | 0.844925–1.64122 |
| 16 | unit_y | 12.1738 | 10.8462 | 1.1224 | 1.00259–1.37073 |
| 16 | entire_rhs | 601.119 | 651.411 | 0.922796 | 0.72411–1.26458 |
| 16 | entire_advance | 2593.36 | 2584.53 | 1.00342 | 0.980618–1.04334 |
| 32 | rho_x | 10.9512 | 21.9612 | 0.498661 | 0.496509–0.573668 |
| 32 | rho_y | 10.1357 | 18.375 | 0.551604 | 0.478862–0.556991 |
| 32 | unit_x | 11.5404 | 20.9138 | 0.551808 | 0.501333–0.724781 |
| 32 | unit_y | 10.2139 | 9.93535 | 1.02804 | 0.989034–1.0599 |
| 32 | entire_rhs | 513.595 | 512.58 | 1.00198 | 0.985707–1.0162 |
| 32 | entire_advance | 2577.74 | 2554.38 | 1.00914 | 1.00141–1.01896 |

两候选全部输出bitwise通过。block16整步只有1.003416x，四块横跨0.980618–1.043343；block32整步1.009144x，四块1.001410–1.018960。block16 RHS存在明显共享机器波动（单次 .541–1.155秒）；block32的RHS样本集中在约.492–.525秒，但整体差仅.2%。这不能支持预设的稳定5%收益。kernel形状、线程、JIT与运行时状态影响明显，旧单线程5%–17% kernel收益没有转化为当前整步收益。

session94800正常exit0，唯一我的大LU窗口已释放。117份冻结源码、source manifest和step81 CP字节哈希均未变；completion明确anyStableWholeGain=false、separateMatchedPhysicalTestLaunched=false、productionModified=false。研究候选保留为逐位等价但当前性能门拒绝的记录；没有隐藏失败样本或扩大线程/网格搜索。
