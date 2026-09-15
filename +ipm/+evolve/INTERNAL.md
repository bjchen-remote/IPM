# evolve：状态与时间推进

协调完整数值生命周期，不另建求解器。外层循环由 `ipm.solve` 持有。
导航：[包索引](../README.md) · [结构总览](../../STRUCTURE.md)。
网格请求交给 [remesh](../+remesh/INTERNAL.md)，结果记录交给
[output](../+output/INTERNAL.md)。

## 路径

- `initialize(opts)`：一次配置解析、网格、初值、尺度、流场及必要的解析初始重网格。
- `advance(state)`：选步长、配套 RK、有限性检查及必要的重网格；成功 REMESH 后发布非终止质量告警。
- `flow`、`rhs`：空间场、特征、规范、缩放输运与源项的共同组装路径。
- `stepSsprk3`、`stepSsprk54`、`stepRk6`：通过 `stepRk` 共享阶段状态契约。
- `isActive`、`selectTimestep`：处理归一化、规范及物理时钟和步数限制。

`state` 保存冻结配置、当前 `ops/rho/flow/scale`、初始质量/值域、步数与时刻。
接受状态只在完整步进与迁移通过有限性检查后保留；失败时回退并报告停止原因，
不吞掉数值函数抛出的异常。最后一个有限状态供外层记录和输出。
同网格恢复由 `output.restoreCheckpoint` 重建 `ops/flow` 后交回同一个外层循环；它不另建
时间推进入口。checkpoint 必须保留参考轴和派生重标度参考量，不能只保存末态 `rho`。

固定网格上，每个显式 RK 步末为接受态计算 `flow` 时已同时生成该状态的
RHS。`state.rhsCache` 仅暂存这一份 RHS 及对应 scale-rate，下一步只在
`rho`、完整 scale、两条坐标轴、`remeshCount` 和 geometry 全部逐值精确相同时
复用为 k1。未命中时立即执行原 `flow` 路径，不改变任何后续 stage 或 RK 算术顺序。
SSPRK3、SSPRK(5,4) 与 RK6 因此分别将连续步的 RHS 求值数从 4/6/9 降为
3/5/8。接受的 remesh 使用新网格上本就需要的 `flow` 重建缓存；未接受步回退
完整 `acceptedState`，包括其旧缓存。缓存是派生瞬态数据，`makeCheckpoint` 不持久化；
`restoreCheckpoint` 在已校验的恢复态上重算 `flow/RHS` 后重建。

## 重标度

物理、各向同性动态、各向异性动态模式共用路径；`Cx=Cy` 是各向同性特例。
正比例尺存对数；物理模式所有尺度变化率为零。各向异性只支持 `fixed`、
`peak_translation` 规范，后者要求 half-plane 对称模式。
每个 RK 阶段使用自己的比例尺、椭圆度量和规范；守恒源系数为 `cx+cy+comega`。
RK6 使用六阶推进行，不用嵌入五阶行自适应控步，也不声称 SSP/TVD/保正。
canonical `tau` 是直接积分时钟；不再有 `maxDynamicRate` 限速或
状态依赖的公共时间重参数化。`canonicalRateMagnitude` 仅是原始率的诊断量，
不参与步进；CFL、终止时钟和有限性安全门仍独立生效。

冻结壁面模板投影以初始离散网格上的模板为定义，目前只准用于固定网格；配置解析会拒绝
其与 adaptive remesh 的组合，避免在尚未注册模板重采样规则时悄悄改变规范泛函。

schema 4 的运行期边界是 `exact_gauge_no_feedback_v1`：`c_l`、`c_omega`和
`c_r` 只由所选规范的瞬时代数恒等式给出。`lengthScaleGain=1`，
`omegaGaugeGain=widthGaugeGain=travelingWaveGain=0`；不得添加幅值、宽度或位移 restoring
项。旧 `adaptiveGain` width controller 不在 RHS/规范路径上。峰值、层级或连通宽度的
网格单元计数可由诊断层产生，但只能驱动 telemetry 和 remesh；网格质量不再终止演化。

### 固定锚点幅值规范

记不含幅值源的 RHS 为 `B`，`R_tau=B+c_omega R`，并用归一化紧支集
光滑 bump 权定义 `⟨·⟩_w`。壁面窗固定于 `X_a=1`、半径
`omegaGaugeWindowRadius` 和开支撑 `|X-X_a|<r`；bulk 窗为上半平面中以
`(X_a,0)` 为圆心的开半圆。窗必须严格留在正 `X` 与计算域内，每个相关坐标
至少有三个正权节点。

R2.2 的 `outer_wall_density_window_l2` 另取 `X=2,Y=0` 附近的平滑余弦平方
窗，默认支撑 `1.5<X<2.5`，以 `R` 本身而非 `R_X` 定幅值：
`M=⟨R^2⟩_w`，`c_omega=-⟨R B⟩_w/M`。在固定网格的半离散方程上
`d_tau M=0`；重网格重新计算求积权和窗值，但保留原始 reference。
窗与 `X≈1` 的局部尖峰分离，避免用可能在幂律下不可积的全域梯度量。
该规范只固定一个积分量，不推出所有外区斜率逐点不变。

- `anchor_wall_template_projection`：初始壁面梯度在窗内归一化为冻结模板
  `T=R_X(0)/||R_X(0)||_w`。令 `P=⟨R_X,T⟩_w`、`F=⟨B_X,T⟩_w`，则
  `c_omega=-F/P` 且原始残差为 `F+c_omega P`。条件数诊断为
  `|P|/(||R_X||_w ||T||_w)`；`P` 必须与冻结初始投影同号且远离零。
  因为投影使用梯度和固定紧支集权，远场密度常数不能主导该规范。
- `anchor_bulk_gradient_l2`：令
  `E=⟨R_X^2+R_Y^2⟩_w`、`F=⟨R_X B_X+R_Y B_Y⟩_w`，则
  `c_omega=-F/E` 且原始残差为 `F+c_omega E`。记录值是 `sqrt(E)`，
  condition 是它与支撑内最大点梯度的比值；冻结参考值和当前 `E` 必须为正。
- `anchor_wall_window_l4`：令
  `M4=⟨R_X^4⟩_w`、`F4=⟨R_X^3 B_X⟩_w`，则
  `c_omega=-F4/M4` 且原始残差为 `F4+c_omega M4`。因子 4 在
  `d_tau M4=4(F4+c_omega M4)` 中消去。记录值是 `M4^(1/4)`，
  energy 另记录 `⟨R_X^2⟩_w`，moment order 必须为 4。

三者的 reference 在初始化时冻结，为条件/符号/正性安全门和 checkpoint
重放证据；它们不产生追踪初值偏差的 restoring 项。

`gradient_energy` 的原有精确条件写为
`E=E_X+E_Y`、`F=F_X+F_Y`、`c_omega=-F/E`，其中
`E_X=<R_X,R_X>`、`E_Y=<R_Y,R_Y>`、
`F_X=<R_X,B_X>`、`F_Y=<R_Y,B_Y>`。运行时另行记录
`omegaGaugeEnergyX/Y`、`omegaGaugeForcingX/Y` 和总
`omegaGaugeForcing`，并用 `omegaGaugeAbsoluteForcing=<|R_X B_X+R_Y B_Y|>`
及 `omegaGaugeForcingCancellationRatio=|F|/omegaGaugeAbsoluteForcing`
辨别大计算域下两个方向、远场贡献与空间正负相消。总量、分量和绝对 forcing
必须全为有限值，否则规范显式拒绝该状态。
决定 `c_omega` 的总能量和总 forcing 仍使用原表达式与求和顺序；分量值不做平滑、
限幅或反馈。

`wall_omega_quadratic_peak` 是局部峰值候选。活动壁面 `R_X` 峰节点及左右
相邻节点在原非均匀 `X` 坐标上唯一确定二次多项式 `Q[R_X]`，其严格内点凹顶点
记为 `X_v`。同一个顶点和同一组三点 Lagrange 权重作用于 `B_X`，定义
`P=Q[R_X](X_v)`、`F=Q[B_X](X_v)`，并令 `c_omega=-F/P`。在活动三点模板
不切换的区间内，包络定理给出 `d_tau P=F+c_omega P=0`；顶点运动项因
`d_X Q[R_X](X_v)=0` 精确消失。运行记录顶点、峰值、曲率条件、插值权条件、
子网格抬升及原始率残差。非凹、非内点、跟踪窗边缘、近零或非有限模板均显式
失败；该规则不平滑、不限幅，也不追踪冻结参考值。活动模板切换仍是需要在
时间/CFL/空间加密中检查的分段光滑事件，不能由代数零残差掩盖。

`gaugeContract` 元数据 schema 2 为长度规范同时记录精确的离散
evaluation 和参数。`wall_omega_width` 与 `wall_density_width`
明示使用分段线性交点；后者对无括号、反向或近零斜率显式失败，
并记录交点区间、斜率与条件数。`omega_peak_location` 的峰点来自
`R_X` 的三点二次插值，而 `R_XX/R_XXX` 使用微分矩阵后的线性插值；
因此它只是代数 phase 候选，不宣称精确离散峰位条件。运行另记录
`R_XX(X_peak)` 归一化缺损、全域 `X R_XXX` 响应条件及包含
长度/平移/幅值项的完整 phase 残差；商式自身的零残差不是峰位证据。

本模块协调其余六个运行模块；`output` 只负责初始化元数据，数值核心不读研究或测试代码。
相关验证：[缩放](../../tests/+ipmtests/+baseline/scaling.m)、[集成](../../tests/+ipmtests/+baseline/integration.m)、
[四阶时间](../../tests/+ipmtests/+fourth/time.m)、[六阶时间](../../tests/+ipmtests/+sixth/time.m)。

## 可选自动网格政策

显式autonomousMesh启用的version1路径在零时刻重新采样解析初值，并冻结接受的base轴/初始不变量。
advance第三输出仅返回纯网格计划；solve在该调用返回后释放旧poisson引用，再执行applyAutonomousMesh。
这样旧接受步的局部引用不与新LU共存。实际迁移按原transfer执行，全部审计和新flow通过后才提交计数。
失败计划回退前一接受步；候选实际迁移全部失败则保留已通过端点门的来源步并按同一A恢复LU。
两类退出都记录原因；错误与数值拒绝分开。memory、真实迁移审计与有界时间窗放入runMetadata。
version1当前仅固定节点数/箱；时间窗有限不等于记录数有固定上限，长期加点与有界采样仍在研究。
完整运行期范围与证据见../../research/longtime_lab/AUTONOMOUS_RUNTIME_INTEGRATION.md。

显式version2在零时刻冻结referenceFamily，按当前实际场枚举分量单调、资源准入的成员。
各级独立最多3提案，不同N禁止keep；原始config N和物理epoch不变。
apply只接受来源时钟/step/level/N匹配的计划，audit和flow通过后同步提交新base/currentLevel、计数和真实决策账本。
原生跨尺寸恢复及每级历史合同见../../research/longtime_lab/AUTONOMOUS_NODE_GROWTH_RUNTIME_V2.md。
version5仍由同一 `advance → planAutonomousMesh → applyAutonomousMesh` 路径运行，
方向族由初始节点数与预设资源上限确定。规划按目标总节点数递增，成员内按最大相邻比
排序；实际迁移失败继续尝试已注册后备，不触碰C率或原接受场。

显式 `initialMeshObservationFallback` 只在原始零时刻启用。先执行原初选；仅有限族耗尽，
或原解析候选全部仅因核心/前沿分辨率失败时，才用固定解析观察格重测并重新规划一次。
原失败、原粗格测量与新观察证据分别保存，不能用观察格冒充实际演化格。
实际候选仍从原解析初值采样，经原数值门及flow验收后提交。恢复使用保存的实际轴和场，
不重新采样初值或运行初选。H32/H64实际回退、H16不回退、默认路径和原生split已验证；
长时从零试验及独立误差资格另行记录。
