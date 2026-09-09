# 非单调幅度搜索：Stage009 实际冻结场

最新可信源为 `adaptive_campaign_v3/stage_009/result.mat`，tau6.693201945463523、step3096、remeshCount8。实际轴与已提交stage008 globally rooted平台逐值相同，因此本研究从真实depth0构造depth1；原immutable root和runtimeRefs没有替换。只用维护equalizer、quadrature/quality和frozen-pair评分，没有LU/PDE、checkpoint restore或事务。

## 冻结方案

研究入口 `ipm_perflab_search_current_axis_amplitude(previousReportFile, outputDirectory)` 接受此前已审计的stage009 `current_axis_report.mat`。先重新核验原场、当前轴、root、globally rooted parent及深度负例。`frozen_protocol.json` 在数值筛选前落盘。

- sigma预注册为.18/.25/.35/.5/.75/1；Gaussian focus取实际core/front交集中心0.9717792027776753，anchorPhysicalPosition=1。
- 每档上界由维护equalizer首次mesh barrier取得（另外预设A≤16）；这不是全参数可行性的数学上界。
- 在[0,A上界]取33点，最多选择两处局部峰邻域，各用5点做4轮细化；重复A缓存。本轮每档恰好41个唯一A，合计246次轴构造。
- 目标为 `min(coreCells/21,frontCells/20)`，不加x padding。y始终由原root y构造，21×1.15=24.15。
- 每个A都从同一entry current x出发，通过维护equalizer的 `initialAmplitude=maximumAmplitude=A` 和不可达到的内部计数目标请求原轴构造。A=0用其baseline返回。这里的大内部目标只防止构造器提前退出；实际评分/接收门始终21/20。
- 每条轴先将entry精确±1的原索引snap回±1，再做完整axis quality和分辨率评分。mesh ratio≤1.08、curvature≤.01、stencil rcond≥1e-9、quadrature ratio≥1e-8、局部quadrature比[.35,1.65]，均保持。
- 六档各自最优中排名前三才执行全场frozen-pair转移。继续使用core21/21、front20、峰跳≤2e-3、质量缺陷≤5e-12、range≤2e-4；真实事务core floor20不降低。axis-only计数与维护全pair计数均校验相差<1e-10。

## 结果

| sigma | A上界 | 最优采样A | core / front | 最小比例 | 完整pair |
|---:|---:|---:|---:|---:|:---:|
| 0.18 | 0.9355535 | 0.328905525 | 21.804113 / 20.754160 | 1.03770798 | 通过 |
| 0.25 | 1.2840091 | 0.418807642 | 21.310408 / 20.288361 | 1.01441803 | 通过 |
| 0.35 | 2.0064027 | 0.556463241 | 20.797207 / 19.820526 | 0.99034321 | 分辨率拒绝 |
| 0.5 | 3.3882931 | 0.774277916 | 20.342312 / 19.370344 | 0.96851719 | 未评分；轴目标不足 |
| 0.75 | 5.8469200 | 1.244754446 | 20.033140 / 19.111821 | 0.95395907 | 未评分；轴目标不足 |
| 1 | 6.4914621 | 1.914474174 | 19.923167 / 19.004350 | 0.94872223 | 未评分；轴目标不足 |

两项完整pair通过：sigma=.18的峰跳1.48755e-5、ratio1.07301570；sigma=.25的峰跳2.26203e-4、ratio1.07274344。二者coreY=24.150000003。候选固定精确±1，均仍未做原生checkpoint事务或短程PDE，不能仅凭此处冻结门直接接入生产。

## 维护搜索是否漏过可行幅度

在完全相同的focus、anchor、sigma和原targets21/20下，额外直接调用维护equalizer进行轴级对照：

| sigma | 维护返回 | 维护A | core / front | 此次搜索找到了维护漏解 |
|---:|:---|---:|---:|:---:|
| 0.18 | target_met | 0.248853802 | 21.000000 / 21.552500 | 否 |
| 0.25 | mesh_gate_before_target | 1.284009059 | 13.986465 / 9.967513 | 是 |
| 0.35 | mesh_gate_before_target | 2.006402673 | 11.113883 / 8.399540 | 否 |
| 0.5 | mesh_gate_before_target | 3.388293102 | 8.584853 / 6.955785 | 否 |
| 0.75 | mesh_gate_before_target | 5.846919965 | 7.963314 / 6.600802 | 否 |
| 1 | mesh_gate_before_target | 6.491462100 | 11.487895 / 8.937482 | 否 |

因此 **sigma=.25 是实际漏解的明确反例**。原幅度递增采样未命中中间可行区，最后返回最大admitted A1.284009；新搜索A.418807642满足21/20，并通过同一完整转移门。sigma=.18的可行性本身不属于漏解，原搜索A.248853802已经满足core21/front21.5525；新搜索提高的是最小分辨率余量。

其余四档是本次有界采样的负结果，不能推广为连续参数全局不可行。最大admitted幅度并不等于最佳核心/前沿分辨率，分辨率门必须独立保留。

## 产物与限制

全部结果位于 `result/verification/performance_lab_20260908/current_axis_nonmonotone_stage009_v1/`：冻结protocol、六档轴搜索MAT与全trace、前三候选、`search_report.mat/json`、`maintained_selector_comparison.mat/json`、无Poisson/PDE profiler审计及执行recipe。

最终helper checkcode0，计算和只读对照进程exit0。最初两条isscalar静态提示在数值启动前被检查门拦下，修复后运行，旧日志保留。初版输出继承了来源审计中的旧实验参数标签；终版把本次参数与 `sourceAuditProvenance` 分开保存，原始序列化文件保存在 `raw_initial_provenance/`，并逐值核验所有数值数组、score和trace完全不变，详见 `final_provenance_audit.json`。

研究不修改维护equalizer或driver。若以后使用该候选，须将depth1与实际事务轴和原root持久绑定；禁止下一stage无证据地复位depth0。真正长时间收益尚需原生事务、短程精度/稳定性及持续运行验证。
