# 有界 t=0 解析网格迭代：实际结果

H32 的一轮初选不足，可以在固定原参考轴与原门下通过第二轮解析观察解决。
这是独立研究入口 `mesh_iterate_initial_box_geometry(previousScreenRun,outDir)` 的结果；
生产 v1/v2 没有接入此循环，没有改候选族、cap、目标、spacing 倍率或验收门。

## 预登记协议

三箱和原输入来自已完成的 `initial_box_geometry_v2_v3/run`；原始 primitive k8、
321×161、version2 policy、target32/32、cap110000 均逐值保留。最多三轮，每轮调用
现有 planner 的原 top3。一旦某候选通过实际解析 core31/31、front20、全部原轴质量、
精确±1 和完整 reference family 质量，立即停止该箱的筛查。

候选密度始终调用同一 `initialDensity` 重新解析采样，使用配对七点 `Dx` 重算 source；
不插值旧场，不求流场、不装配 LU，不创建 PDE 历史或新的物理 epoch。每个观察的
physical/canonical time 和 accepted step 都显式记录为0。所有轮的设计 reference
仍是该箱**原配置 uniform 轴**；它不会被失败候选替换。只有最终纯资格候选附带
可供未来原生初始化审查的 `referenceAxisFamily`，不称它已经在 solver 中提交。

只有本轮全部候选的拒绝原因均为实际 core/front 不足，且其轴质量、精确锚点和
完整 family 均通过时，才允许下一轮。选择分数为
`min([actualCore./[31 31], actualFront/20])`，最大者优先，原候选顺序稳定破同分。
其它原因立即停止；三轮仍失败也如实停止。

## 实测

| 箱半宽 / ymax | 接受轮次 / 候选 | 实际 core X/Y | 前沿 | X/Y 相邻比 |
|---|---|---|---:|---|
| 8 / 4 | 1 / 1 | 34.859451 / 39.322607 | 25.281313 | 1.039226 / 1.000000 |
| 16 / 8 | 1 / 1 | 34.833744 / 36.361138 | 24.803930 | 1.051377 / 1.038025 |
| 32 / 16 | 2 / 1 | 34.665402 / 36.790981 | 25.117156 | 1.061070 / 1.048065 |

所有箱第一轮的全部 planner 输出轴与前次筛查逐值一致；实际评估的第一轮解析
测量也逐值一致。H8/H16 第一候选直接通过，后续候选不再评估。

H32 第一轮的三个原失败保留：Xcore 均27.630978，Ycore 为37.887987、38.211785、
38.483430，front均22.222307。分数均为.8913218735448346，因此稳定选择候选1
作为下一轮 t=0 观察，而不是临时增大 spacing 参数或挑选较大 Ycore。

粗源原始90%宽为 .47920317，受单元下限放大后的设计宽为 .60603965，放大因子
1.26468205。第二轮 source 来自候选1的同一解析初值，原始宽和设计宽均为
.48143532，放大因子恢复为1。重新设计后第一个候选的实际 core34.665402/36.790981
和front25.117156全部通过，四个 reference family 成员也全部通过原几何门；
资源 mask 仍为 `[true true true false]`。

共执行4次纯 planner 调用和6个解析候选测量，均在预登记上限内。
最终 H32 产物为
`initial_box_iteration_v2_v1/run/box_3_round_2_candidate_1.mat`，包含解析 view、
实际 feature/measurement、candidate、原 reference、config、audit 和完整 family。

## 执行证据与限制

权威目录：`result/verification/autonomous_runtime_20260909/initial_box_iteration_v2_v1/`。
MATLAB session76205 exit0，Code Analyzer 为0；所有输入、冻结源 SHA 未变，维护研究
入口与执行冻结文件逐字节一致。见 `run/report.json`、`final_audit.json` 与
`readable_summary.json`。原一轮结果的 H32 `false` 完全保留，不覆盖或重新解释。

报告中 `allBoxesPureGeometryQualified=true` 与 `nativeInitializationQualification=false`
同时成立。新循环仅证明这三项解析初选的有界解决能力；未计算流场 safety、规范条件、
RHS 或动态误差，也未检验一般初值、多峰或更大箱。将来若迁入生产，应注册新的初始
规划版本并保留该轮次与失败证据；不能因此改称现有 v2 已具备这项能力。
