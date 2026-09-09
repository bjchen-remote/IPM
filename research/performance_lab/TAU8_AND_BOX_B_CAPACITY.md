# τ8实际末态与B容量停止：无LU审查

本轮研究只读真实原生CP/result，重建配对七点Dx、实际rho*Dx及有限候选轴。没有restore、椭圆LU、flow、PDE推进或transfer。新的轴分辨率均为当前真实特征在候选轴上的几何计数，不能当迁移后实际core或动态资格。

## τ8仍有同级容量

step4545/τ8/t=2.2412067821954538，level4=641×321，25次迁移，两次增长额度已用完。当前实际core=27.38530409/29.72754570、front31.68220、安全.1461506671。原70×4有限搜索有35个X、3个Y，105个配对；原排序top3皆为迁移提案，无qualified keep。首项几何core=34.48462320/36.80000017、front39.21142516，尚未transfer。

相较既有τ6.20002014审查的54个X，当前X可行族缩至35个。现在35个X触及相邻比门，叠加16个局部求积下门与8个上门；Y仅sigma.18不足目标。不能由两次实际观测推算未来耗尽时刻，也不能称固定N全局不可能。

|质量指标|原硬门|当前X / Y|首提案X / Y|
|---|---|---|---|
|最大相邻比|≤1.08|1.05375845 / 1.04191707|1.05601840 / 1.04356353|
|log spacing curvature|≤.01|.00326111 / .000598228|.00340206 / .000621231|
|局部七点stencil rcond|≥1e-9|9.19937e-6 / 1.00779e-5|9.14811e-6 / 1.00761e-5|
|全箱平均归一化最小q|≥1e-8|.00488682 / .000653640|.00380063 / .000499544|
|q/control width 最小|≥.35|.503467 / .636682|.495648 / .636650|
|q/control width 最大|≤1.65|1.447973 / 1.381080|1.451279 / 1.381261|

所有q严格为正。此stencil rcond是局部FD诊断，不能当椭圆A条件数或前向误差界。
当前触发26的余量为1.3853/3.7275，保存的原.2τ预测24.00187/23.86437仍高于22，lastDecision=resolved。剩余275次迁移、累计峰跳额度.01820965543不是未来时长保证。完整.35τ窗口当前374条；v2是有限时间窗，尚非固定记录数上界。

## B的实际失败与反事实严格分开

B_H16_T32保存CP为step1520/τ3.819993033851268，物理t1.8011493615315282，未到登记1.9140859724939365，原stopReason=autonomous_mesh_axis_capacity保持。实际失败发生在未提交并回退的step1521/τ3.822215542624246；保存lastFailure不是CP1520同场。

保存失败证据：forecast_core_trigger，源core27.94013845/31.61670472，预测24.38299995/2.65060876。level2有67/70个X而0/4个Y，四Y都拒于y_core_target；cap110000使level4 resourceAdmitted=false。最近Y最大割线12.39448来自1520→1521，不能在没有1521完整rho时把该诊断变化直接归因为某种峰模板。

当前CP1520独立重筛有66X/1Y、3提案，其唯一合格Y为sigma.5，core32.49975188，所有四sigma [.18,.25,.35,.5]均在相邻比约1.08处达到mesh_gate_before_target；各core为2.733692/4.529482/10.925967/32.499752。现轴Y相邻比1.079999999979997，仅约2e-11余量。它在保存态尚可行不推翻下一步真实容量失败。

同一CP1520 rho、原初选base、原全部搜索/质量门，只把另注册的cap设210000后，原已有level4可进入：66X/4Y、3提案，首项几何core35.194644/36.80000006、front32.102456；首Y相邻比1.033972。family的唯一数值变化是resourceAdmitted mask，轴/端点/全部quality逐值相同。本项不是原B恢复，也不是1521失败场的独立重演；需全新t0登记运行和实际事务证明。

## 下一版最小有限方向族提案（未实现）

先按原v2保留τ8的同level续跑。若之后需要新增方向，最小研究提案是原四成员加[3,2]，从原始t0初选selectedBase直接做归一化索引PCHIP，cap310000，最多3次增长。当前额外两方向都已注册筛查：

|相对原始节点的cell factors|N / 节点积|当前实际场的纯几何可行轴|首提案core X/Y|
|---|---|---|---|
|[3,2]|961×321 / 308481|70X / 3Y|35.181388 / 36.800000|
|[2,3]|641×481 / 308321|35X / 4Y|34.484623 / 36.800000|

factor3的所有新reference轴通过全部原门，精确保留原selected knots、端点、0、±1；factor2独立重算与原family逐值相同。factor3不声称保留factor2的全部中间reference节点，更不声称真实PDE网格嵌套。当前[3,2]把X邻比降至1.032090并使70个X都可行；[2,3]不能消除当前X瓶颈，因此不是本轮最小首选。此判断不外推未来Y场。

未来版本必须重新登记规则名、完整不可变family、cap、增长预算、排序与总pair预算。五成员最多15提案（每member原上限3）；current/target level、N/base、componentwise不减、root来源、source/target审计与完整history/CP均须严格支持新版本。当前v1/v2及其旧checkpoint不可升级；不能把这次研究的临时geometry cap写回原CP。所有core31/31、front20、峰.002/累计.02、质量5e-12、range2e-4、安全.70及轴质量门保持。

## 16 GiB的资源边界

原τ8整次运行实际maximum RSS=3653304320 bytes、peak footprint=3821654360 bytes（3.5592 GiB），0 swaps；这是205761点真实运行证据。原cap110000此前peak=3221622696 bytes。308481点没有新的LU实测，本报告不能认证其与大LU主线并行可承受；节点数增加约49.92%不等于LU内存线性增长。一次ps资源快照在sandbox内不可用，未据它作容量判断。

下一版必须先做独立资源资格，再从t0冻结启用：最少取当前最优新增几何与最差质量余量候选，串行创建一次实际原矩阵factor/flow，确认旧factor和所有alias已经释放，真实峰内存与可靠恢复/实际迁移门通过。建议预登记附加进程footprint上限5 GiB，只有“主线已测峰 + 附加5 GiB + 至少3 GiB系统/安全余量 ≤16 GiB”且实测无新增swap才能允许并行；若主线峰未测或总预算不够，只能串行资格运行，不可声称已准入。该5 GiB是待检验的资源拒绝界，不是本轮测量结果；超界应拒绝更大族。后续长跑仍遵循主线+至多一个附加LU，完整CP后释放，失败恢复旧A时不与新factor共存。

更长τ仍受同级有限几何族、真实transfer峰/质量/range/safety、累计峰预算、计数/步数、CFL成本及有限物理箱误差约束。当前计算箱[-8,8]×[0,4]映射到末态物理箱约±1.72918、Ymax .864590；几何审查没有给出域误差或奇异性资格。

## 实际验证与保存

主审查session19722 exit0：B纯审查7.63秒，τ8含两额外几何11.17秒；原同级规划2.44秒，两额外方向5.10/2.60秒。另session58489仅把既存MAT候选质量导出JSON，未重新计算几何。两个helper CodeAnalyzer均0，profile无restore/build/flow/Poisson/transfer/advance/solve调用。维护+ipm未编辑，也未重跑已有suite。

CP/result逐值核对rho、轴、全部history、config、时钟与controller/metadata；原路径latestCheckpointFile被单独准确记录，因为maybeCheckpoint先序列化immutableCP，再更新运行metadata指向新文件。除该输出路径外所有metadata逐值相同。这不是CP与result整个文件或未保存flow缓存的逐位等价证明。

v1只因过期CodeAnalyzer抑制标记停止，未读CP；v2严格读取通过后，因错误要求latestCheckpointFile也相同而停止。两次研究审计失败与旧产物保持；v3仅修正这两个审计假设，所有共有数学量继续isequaln，无容差修改。

主要文件：B/report.json、B/saved_last_failure.json、tau8/report.json、两处report.mat（完整candidate轴/每项axisReport）、axis_detail_export.json、manifest_before.json、final_audit.json。原8个输入及119冻结文件SHA最终逐个复验。

证据根：`result/verification/autonomous_runtime_20260909/tau8_and_B_capacity_nolu_v3`。
