# 唯一从零时间步控制对照

本项补充 H8/H16 × target32/42 矩阵中尚未覆盖的时间离散敏感性。
只登记一个新算例 E，复用已经完成的 A，不重复跑第二个粗基准。

| 项目 | A（真实复用） | E（本次仅登记） |
|---|---|---|
| 原解析初值、箱 | primitive k8，[-8,8]×[0,4] | 完全相同 |
| 初始配置节点数 | 321×161 | 完全相同 |
| 冻结自动网格 | v2，target32/32，cap110000 | 完全相同 |
| CFL | .35 | .175 |
| maxDt | .005 | .0025 |
| 最大接受步数 | 30000 | 60000 |
| 终止条件 | 实际 canonical τ=4 | finalTime=Inf |
| 共同物理终点 | 1.9140859724939365 | 精确匹配 A 实际值 |

原空间/输运/RK/椭圆选择、WENO epsilon=1e-12、网格搜索族、质量门、事务与累计误差门不变。
新caseId、原物理/规范时间零、解析初值重新采样；不继承 A 的演化历史、不调用 checkpointFromResult。
CFL 与 maxDt 同时减半，避免仅减 CFL 在 maxDt 主导时没有实际减小步长。
maxSteps60000只是上限，不能作为成功条件；容量、时间步或可信度停止均保留原生CP并报告未完成。

准备入口 `mesh_prepare_time_control(factorialRegistration,frozenSourceDirectory,newPlanDirectory)`
严格重读 A 的真实 result/CP，核对全场、轴、尺度、时钟、history、config及原始初始化/控制器来源。
它逐值核对除时间控制和明确输出路径/描述外的全部数值配置域；
初始几何资格复用四格预飞中 A 已与真实 selectedBaseAxes 逐值相同的结果，
因为原解析物理、网格、reference、policy均逐值不变，而纯planner不读取时间步参数。
这不是 E 已通过原生初始化或时间推进的证据。

唯一 dormant runner 是登记目录 `launch.m`，调用已冻结的
`mesh_run_box_space_case(registration.mat,'E_H8_T32_HALF_TIME')`。
同一个已测单例运行器检查source SHA/依赖闭包、纯初选资格及新输出目录，
再调用唯一 `ipm.solve`，结束后严格读取真正原生末态，并与 result 和期望配置配对。
它保持10线程签名验证，拒绝未到共同物理终点、非完整trusted、原历史核心/安全门、
endpoint核心20/20、registered family或累计峰跳失败。
此次准备没有启动它，没有新增LU/PDE；root在资源窗口单独排队。

## 共同物理观察与解释

E−A在同物理终点比较，沿用 [四格预注册观察协议](FROM_ZERO_BOX_SPACE_FACTORIAL.md)：
core [0,.8]×[0,.25] 与 holdout [0,2]×[0,1]；每支必须实际覆盖，禁止外推。
双向linear与Hermite场差并列；原legacy离散宽判定原样保留，连续同峰90%宽.005辅助门独立报告。
rho L2/Inf≤1e-3，rhoX/rhoY≤2e-3，rhoX最大值/二次峰≤2e-3；
物理时间log尺度率 abs(diff)≤max(1e-6,.002|reference|)，两个方向均须满足，canonical率也并列。
固定物理peak-window [.1,.6]映回每支真实规范坐标后，才使用连续几何helper。

明确分列 `history.common.physicalGradInf=max hypot(rhoX,rhoY)` 与
`history.common.physicalRhoXInf=max abs(rhoX)`，二次拟合正峰又是第三个量。
旧 `fresh_profile_campaign_review.physicalGradientMaximum` 以及旧 comparison 接口的同名字段
实际指rhoXInf，禁止用英文名猜作full gradient。E登记中保存这些映射。

逐步CFL/实际dt、每次remesh与growth的step/真实时钟、source/target级别和累计峰跳一起保存。
E的自动事件时刻和网格轨迹可能变化；两档结果可检验时间步与自动事件的联合敏感性，
不能单靠一次减半声称固定网格RK阶、自动迁移整体阶或时间误差已可忽略。
是否达到比较门必须由实际终态报告决定，不能由初始配置一致或减半参数推断。

## 实际准备状态

`result/verification/autonomous_runtime_20260909/time_control_half_v2` 已完成 noLU 预飞，
session5726 exit0，新prepare helper CodeAnalyzer=0；baseline真实原生配对、非时间数值域逐值相同、
checkpoint策略相同、完整源SHA/求解依赖闭包均通过。profile没有solve/build/flow/transfer/advance。
唯一待排队文件是该目录 `launch.m`；`plan/registration.mat` 是其完整冻结登记。
`plan/preflight_report.json` 和 `final_preparation_audit.json` 保存实际准备结果，**E的solve尚未执行**。
v1临时预飞文件括号错误导致解析期退出，失败日志保留，没有启动任何数值运行。

补充避免不同研究接口的旧名称混淆：本次复用的
`mesh_audit_box_space_case.physicalGradientMaximum`指 **physicalGradInf**，
其`physicalRhoXMaximum`才是physicalRhoXInf；这与旧fresh review/comparator同名字段的含义不同。
`final_preparation_audit.json`逐项登记了该映射，后续比较应按明确的history原字段读取，
不能直接把不同helper的`physicalGradientMaximum`互相比较。
