# 从原解析 t=0 的箱子 × 核心分辨率对照

本注册只检验现有自动流程在共同物理终点的敏感性，不证明无穷域、空间收敛阶或奇性。
新运行入口仍只有 `ipm.solve`。准备阶段不 build、flow、LU、transfer 或推进时间。

## 最少的四项矩阵

| 标签 | 初始计算箱 | targetCoreCells | 处理 |
|---|---|---|---|
| A_H8_T32 | [-8,8] × [0,4] | 32/32 | 复用真实 baseline step1653 |
| B_H16_T32 | [-16,16] × [0,8] | 32/32 | 注册新增 |
| C_H8_T42 | [-8,8] × [0,4] | 42/42 | 注册新增 |
| D_H16_T42 | [-16,16] × [0,8] | 42/42 | 注册新增 |

全部由同一 `degenerate_primitive` k=8 解析初值、321×161 配置初始轴、零物理/规范时钟开始。
使用 A 的完整冻结 v2 生产源码。唯一变化是箱子、目标核心数、终止条件、输出和明确的新 caseMetadata；
公共离散、椭圆选择、CFL、maxDt=.005、WENO epsilon、搜索族、原质量门与 cap110000 保持。
不继承 A 的历史/参考量，不把新算例标成 checkpoint continuation。

共同终点固定为 A 实际 `physicalTime=1.9140859724939365`；A 的 canonical tau=4 只是其自身时钟。
新增三项设置 `finalTime=Inf`、`physicalFinalTime` 为上述值、maxSteps=30000。
最大步数、节点容量或任意停止原因本身不算抵达终点。

policy 的既有 target 参数同时缩放触发和端点要求，因此 T42 并非只换轴后静态加密。
T32 的事务/端点核心门是31/31、20/20；T42 是41/41、26.25/26.25。
两者全部历史核心≥14、安全系数<.70、前沿≥20、单次峰跳≤.002、累计绝对峰跳≤.02，
原质量/质量守恒/范围门不变。节点族预登记321×161、641×161、321×321、641×321；
最后一个205761节点成员不在110000资源预算内。即使无资源准入，其参考轴质量仍须通过。

## 初选与运行资格

`mesh_prepare_box_space_factorial(baselineDirectory,frozenSourceDirectory,newPlanDirectory)`
严格用10线程读取 A 的原生 CP/result，核对场、轴、尺度、时钟、完整history、config、case及controller。
然后每项仅调用生产 `plannedAxisPairs` 一轮/top3；用候选轴重新采样原解析初值和配对7点Dx，
检查实际核心/前沿、全部轴门、精确±anchor及完整referenceAxisFamily。
A 初选轴还须与真实 t=0 selectedBaseAxes 逐值相同。

这里的 pureGeometryLaunchable 只表示候选几何可行；候选流场安全、实际 native 初始化与整个动态区间尚未验证。
不使用 H32 研究中的解析迭代，不扩大搜索族，不对失败候选临时乘 spacing 系数。
每项失败审计、全部 top3 和源粗测→候选真实测量差保存。

`mesh_run_box_space_case(registrationFile,label)` 只接 B/C/D 中一个标签；校验冻结源码 SHA 和依赖闭包，
拒绝已有输出目录或未通过纯初选的标签，之后才调用一次 `ipm.solve`。
每项单独 MATLAB 进程、默认10线程，结束后退出释放全部 LU；root 只排主线 + 一个 addon。
不得把准备文件、登记或最大步数退出写作已完成。

`mesh_audit_box_space_case(resultFile,checkpointFile,targetPhysicalTime,expectedConfig)`
无 LU 重读原生终态并精确核对 result；要求 t=0 初值来源、全历史可信、原箱和注册family、
累计误差与目标相关端点门，以及实际共同物理终点。它也用于复用 A。
已有 V2 checkpointFromResult 不支持变量节点族；本协议只接真正原生 CP，绝不重建伪 CP。

## 预先固定的共同观察协议

所有比较只在相同物理时刻；同一坐标映射 `xphys=(x-Xshift)/Cx, yphys=y/Cx`，
同一物理场 `rho/Comega, (rho*Dx')Cx/Comega, (Dy*rho)Cx/Comega`。
固定两个物理窗口：core [0,.8]×[0,.25]，holdout [0,2]×[0,1]。
A 实际物理箱约±2.44475×1.22238，覆盖两窗；其余每项都须另验覆盖，缺少则保留失败，禁止外推填值。

四条边分别 B−A、D−C（箱子方向），C−A、D−B（目标方向）。
每条边每窗口并列保存：

- 双向线性插值 rho 的相对 L2/Inf≤1e-3，rhoX/rhoY≤2e-3。
- 双向 continuous-Hermite 观察及其与线性观察的差；它不替换原 linear 判定。
- rhoX 最大值、二次峰值相对差≤2e-3；full physical gradient G 另列，不能把 rhoXInf 称为 full G。
- 原 legacy nodal-column 核宽相对差≤.005，单独原样保存。
- 同一连续峰的90%横纵宽相对差≤.005，作为预注册辅助判定。
- 规范率和物理时间率均保存；主要同物理时刻率门在 d(log C)/dt 上应用
  `abs(diff)<=max(1e-6,.002*abs(reference))`，两个方向均须通过；原 canonical 率门也并列报告。

使用既有 `mesh_compare_fresh_box_data` 的 linear/legacy 门和 `compare_space_observers` 的双向linear/Hermite观察。
旧 comparator 的字段 `physicalGradientMaximum` 指 rhoX 峰，需明确映射，不能误填 full G。
率必须在调用 comparator 前变换；同一 native canonical tau 不能代替共同物理时刻。
连续几何用 `continuous_inner_geometry`，每支将同一物理 peak-window [.1,.6] 映回自己的规范坐标。
窗口必须包含唯一非退化连续正峰及两侧90%交点；不满足就拒绝此观察，不另选对自己有利的峰。
该 observer 的绝对误差及单位协变已实测，而逐相位每层误差下降≥8未全部通过；不宣称通用四阶几何。

交互项 `D−C−B+A` 必须先在相同固定物理 query nodes 上作有符号相减，再取范数；
每窗321×161均匀query点，只作为观察算子，不冒充求解网格或提高有效分辨率。
linear、Hermite分别报告，以 A 同窗范数归一；标量同时给绝对交互和按A归一的相对交互。
四条边必须在每一方向先比场，禁止只比较标量范数后再称场差。

每项一起记录实际节点数、网格轴、重网格/增长时刻、CFL/实际步长、累计峰跳、有效WENO物理幅度尺度。
初始 selectedBaseAxes、参考量、规范尺度、节点族轨迹随箱子与目标不同，故矩阵是**箱子与自适应网格联合敏感性**。
2×2差分能显示两因素及交互是否重要，不能完全拆出纯边界误差；两档分辨率也不能估计可靠收敛阶。
farVelocityRatio 不能当边界误差。

## 时间误差的下一最小项（未派发）

若 B/C/D 都通过完整区间资格，且空间/箱子差已可解释，建议最多再做 D 的一档 CFL/2，
其他全部不变，仍从原 t=0 到同一物理终点。先检查实际 dt 中 CFL 限制比例：
如果 maxDt 主导，CFL/2未必把有效时间步减半，此时应先明确新的时间步控制注册。
自动事件时刻可能随 CFL 改变，所以这一项衡量时间步与触发序列联合敏感性，不是固定网格时间阶验证。
本轮只注册三项新增运行，不预先启动或扩成五项额外计算。

## 2026-09-09 已执行的纯预飞

最终目录：`result/verification/autonomous_runtime_20260909/box_space_factorial_v3`。
MATLAB session65453 exit0；四个研究文件 CodeAnalyzer=0、冻结生产源码与 baseline 逐文件相同，SHA复验通过。
`plan/preflight_report.json` 保存原生 A 配对和四项初选；没有启动任何新增 solve。
原 CFL=.35、maxDt=.005、configured WENO epsilon=1e-12。
附加无 LU 边界检查 session50835 exit0：错误共同物理终点、错误 expectedConfig、重跑已有 A 均被拒绝，
记录在 `plan/boundary_checks.json`。Dormant runner 的真正 solve 分支尚未执行。

| 标签 | 初选实际核心 X/Y | 前沿单元 | X/Y 最大相邻比 | 纯几何 |
|---|---|---|---|---|
| A | 34.85945 / 39.32261 | 25.28131 | 1.03923 / 1.00000 | 通过，轴与真实初选逐值相同 |
| B | 34.83374 / 36.36114 | 24.80393 | 1.05138 / 1.03802 | 通过 |
| C | 45.56645 / 47.74219 | 25.50348 | 1.04966 / 1.01584 | 通过 |
| D | 45.61563 / 47.96637 | 24.14394 | 1.06002 / 1.04137 | 通过 |

全部候选的完整参考族质量通过；budget mask 均为 true/true/true/false。
该表不是动态结果，更细目标可能更早用尽方向增长预算或容量。
B/C/D 的 dormant 单例 batch 文件分别在上述目录的 `launch/`；root 可逐个排队，不能并发自动派发三项。
执行 `matlab -batch "run('.../launch/B_H16_T32.m')"`（完整路径需由 root 指定），随后退出该 MATLAB 再派下一项。

失败证据保留：v1 在无数值执行前被研究 SHA helper 的 CodeAnalyzer 风格提示拒绝；
v2 在只读 baseline 审计中因把 physicalQuadraticPeak 写到 gauge 而非 common 组报错；
v3 修正这两个研究脚本问题，生产源码、阈值和旧产物均未更改。
