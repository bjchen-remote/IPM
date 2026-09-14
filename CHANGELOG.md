# 本副本变更记录

## 2026-09-14：极端冻结态同节点核心余量 / R2.5 候选

- 新增仅对自动网格 version 5 显式生效的 `meshDensityVersion=2`：沿用 v1 的单侧几何 Y 密度与 70 个 X 候选，把冻结的 `search.xPadding` 从 1.10 改成 1.50。省略字段、`0` 与既有 `1` 的配置和旧检查点保持原值；未增加节点预算或放宽相邻比 1.08、曲率、核心、前沿及原生场迁移门。
- 验签旧 C `τ=12.400290` 的原生 `641×321→961×321` v1 候选通过，X/Y 最大比 `1.044411/1.038809`，实际核心 `35.152/36.981`。旧 C `τ=14.200201` 的同节点 v1 候选也通过，但 X 核心仅 `32.889`，高于设计目标不足 1 格；这是后期寿命余量问题，而非仅仅网格比集中。
- 更强 X 密度在旧 C `τ=14.100094/14.200201` 的 961×321 原生反事实迁移均通过，X 核心 `45.288/47.767`、X 最大比 `1.050411/1.059386`；新 C `τ=3.500136` 的 641×321 正式 v2 候选原生迁移通过，X/Y 核心 `46.456/36.158`、最大比 `1.050595/1.026118`。极端旧态 X 的格宽变比负担参与率约 79.71%，最重 5% 接口承担约 6.42%；Y 参与率约 99.17%。详细逐态 JSON 与方法边界见 `research/longtime_lab/EXTREME_FROZEN_MESH_RESERVE_20260914.md`。
- v2 从原始物理 `t=0` 自主推进到 canonical `τ=.5`，304 步、6 次自然重布、正常到达目标，末 X/Y 核心 `43.442/36.230`、最大比 `1.049928/1.043198`。另从 `τ=.25` 验签 checkpoint 续算到 `.30`，自然重布 `1→2` 且冻结版本 2 不变。公开服务器接口与完整 `ipm.verify('all','quick')`（94.6515 秒）均 exit0；新模式连续长时、空间/箱长收敛仍未认证。
- 旧数值等价快速套件 `ipm.verify('equivalence','quick')` 也通过（6.2773 秒）。

## 2026-09-14：显式同节点均摊密度函数试验

- 新增 `longTimeProfile` 的 `meshDensityVersion=1` 显式选项（只配自动网格 version 5）；缺省 0 的冻结 policy 不增加字段，既有检查点保持原密度规则。新模式将第四个 Y 候选换成单侧圆滑几何密度，让外区对数格宽增长几乎均摊；X 的 70 个既定候选预算中加入 8-cell 过渡，不增加节点或放宽质量门。
- 验签冻结态 τ3.500136、641×321：正式新模式首个改变网格候选 X/Y 最大相邻比 1.047919/1.026118，旧实际网格 1.049257/1.039181。原生迁移审计和后续 flow 均通过，实际核心 X/Y=34.946/35.840、前沿34.373、相对峰跳7.59e-5、质量缺损0。Y 接口对数负担参与率由77.19%增至99.59%。同场手选参数还达到 X=1.047314，但非正式候选；不冒充已上线结果。
- 旧 C τ12.400290 冻结态上，两策略在**同一目标 961×321** 的首个纯几何候选总最大比分别为旧1.059659、新1.044411；新候选预测核心35.144/36.800、前沿39.852，尚未做该大网格的原生迁移。旧实际同级641×321的X已接近均摊，不能把新增节点说成密度改进。
- 新模式从物理t=0到τ=.5，277步、6次自然重布、正常到达设定时间；原生检查点τ=.25→.30续跑通过，密度版本冻结。前两次实际Y选择仍为旧Gaussian候选1；新Y的正式选择和原生迁移证据来自τ3.5冻结态，尚无新模式连续长跑。定向服务器接口、`ipm.verify('all','quick')`（100.9264秒）、六物理算例/三迁移 `ipm.verify('equivalence','quick')`（6.7301秒）均exit0。图、脚本、旧浮点检查点签名对线程数敏感的首次失败与详细边界见 `research/longtime_lab/BALANCED_MESH_DENSITY_20260914.md`。

## 2026-09-14：R2.3 候选版外区截断范数与新 C 观察器修复

- 新增 `server/analyze_outer_profile` 用户入口：从验签的接受态断点只读恢复原生 RHS，在重标定固定坐标测量避开 `(1,0)` 半圆核的 bulk/壁面 `R_{Xτ}` 截断 L1/L2/L4/L∞，以及两态公共观察格的 Profile 增量；默认输出新的 JSON，不修改轨道。半径 0.05/0.1/0.2/0.4、外窗 `[0,4]×[0,2]` 与 401×201/801×401 观察格均明确记录。候选截断梯度 C 率只作为冻结态比较，不替换 R2.2 的 `(2,0)` C。
- R2.2 新 C 的连续纵向核心影子观测器每个接受步因旧二次峰规范专有位置为 NaN 而失败，原生网格决策未受影响。新版本在该专有点缺失时改用已计算的 `trackedPeakX`；旧二次峰分支保持不变。真实 `τ=3.00041744` 新 C 断点的只读恢复观测得到连续 Y 核心 34.6993 格；独立原始 `t=0` 10 步短跑到 `τ=.01`，影子观测 10 次、零失败、零请求分歧。R2.2 正在运行的轨道仍使用旧发布源码。
- 旧峰值 C 的四个验签冻结态独立归一化后，`δ=.2` bulk L2 外区相对增量为 0.00523/0.00676/0.00157；新 C 原始 `t=0` 轨道到 `τ=3.00041744` 已自主重布 21 次、网格 641×161，同一范数的 `R_{Xτ}` 为 0.12154/0.10165/0.09278/0.09822（`τ≈1.5/2/2.5/3`）。新轨道还没有外区收敛证据；H8 远边界速度/源比 0.9872/0.1599 超过 0.01 告警阈值，最终 Profile 需大箱与时空误差比较。
- 合成稳态/均匀增长及双观察格截断范数测试、合成与真实断点的观察器测试、用户分析入口均实际通过；修改后完整 `ipm.verify('all','quick')`（105.9568 秒）及 `ipm.verify('equivalence','quick')`（8.0027 秒）均 exit0。新 C 长时强奇异性及外区无穷时间极限仍未认证。

## 2026-09-14：R2.2 外壁面密度幅值规范

- 新增 `outer_wall_density_window_l2`：固定 `(2,0)` 周围半宽 0.5 的平滑壁面密度加权 L2 值，按 `c_omega=-<R,B>_w/<R,R>_w` 精确选率。它避开 `X≈1` 的局部尖峰，因此不依赖可能在最终幂律下发散的全域梯度积分。初始参考值、checkpoint 恢复、历史诊断和实际网格迁移均已接入。固定网格半离散窗值率残差为零；自动迁移的离散窗值仍可能有跳变。
- 服务器 R2.2 的新算例默认新 C，保留旧 `wall_omega_quadratic_peak` 作为明确可选对照；自动网格 version 4/5 的准入同步允许新规范。旧 R2/R2.1 检查点保持原包续算，不能用新默认规范更换冻结配置。
- 四个验签旧轨道冻结态 `τ=12.3003–14.2002` 的 `(2,0)` 密度降低 44.664%，同时外壁样点 `R_X/R(2,0)` 变化均小于 0.39%。`τ=14.2002` 原生 RHS 的新窗率为 `+0.00149598426`；三步内存短分支的窗值相对漂移 `2.22e−16`、点值变化 `2.05e−10`，没有写回旧断点。从原始 `t=0` 的 version-5 服务器配置跑 10 步到 `τ=0.01`，窗值漂移 `2.22e−16`，尚未触发后续自动重网格。独立小网格短步、记录、参考恢复及迁移测试通过。
- 修改后 `ipmtests.layout`、完整 baseline、fourth quick、sixth 及六算例/三迁移 equivalence 分别直接运行并 exit0。`ipm.verify('all','quick')` 两次在 MATLAB 进程启动阶段出现 Qt/NEON fatal error（exit137），未进入测试；同环境独立布局和所有组成套件随后通过。R2.2 尚无新规范从零到长时间、多网格/大箱收敛或最终幂律认证。

## 2026-09-13：第五版自动初始网格、原生事务与服务器包

- 服务器第五版默认 `initialNodeCount='auto'`。H8/Y4 保持原 321×161 起点；对更大的计算箱，在建立二维 Poisson 算子和 LU 之前，从原始解析 `t=0` 数据按二维节点成本筛选质量合格的 X/Y 初始数量，并将选中值冻结进配置与 checkpoint。H64/Y32、H128/Y64 仍选 321×161；H256/Y128 的固定 321×161 初始候选耗尽，自动选 401×161，原生初始化、两步推进及 checkpoint 续算到 τ=.008 均通过。固定 H256 失败日志改为可复现的简短摘要。
- H8 第五版从原始 `t=0` 自主推进至 τ=.5，共 273 步、5 次自然自动重布，无人工中途改网格。另用原生 `mesh.build`、`remesh.transfer`、审计及 flow 验证早态 321×161→641×161 方向增点；这是直接迁移试验，并非自然演化中已经触发的增点。
- 晚期可信 τ≈9.300547 冻结态复核：节点成本 205761、308481、462241、616161、821121 的合格最差相邻比分别为 1.062837、1.048762、1.036248、1.030233、1.025002；3×2 与 3×3 的独立无 LU 前向/回转评分通过。前沿与 τ≈7.400881 相同，但晚态数值更差，说明需持续动态补足瓶颈方向；不构成动态网格误差界。
- 修改后完整 `ipm.verify('all','quick')` 和六物理算例/三原生迁移的 `ipm.verify('equivalence','quick')` 均实际 exit0。H64/H128/H256 短程算例仍有未消除的计算箱告警，第五版大网格长程和独立空间/时间/箱长收敛尚未验收。

## 2026-09-13：隔离试验中的 version-5 自适应方向节点族

- 新策略从原始 `Nx×Ny` 与启动时 `maximumTotalNodes` 冻结方向增长级别，保持同箱、精确锚点、分量单调增长、完整候选质量和实际迁移验收。每级已注册候选按最大相邻比最小、质量裕量和稳定索引排序；旧 version-4 默认配置与排序不变。
- 原始 H8 `t=0` 轨道的可信 τ≈7.400881 冻结态完成无 LU 测试：预算 1000000 登记 38 级；`641×321`、`961×321`、`641×481`、`1281×321`、`961×481` 的最佳最大相邻比分别为 `1.052817709`、`1.041135950`、`1.052817709`、`1.041135950`、`1.029972698`。首个返回的改变网格候选与逐对枚举的最小最大比一致；尚未做 version-5 的原生大网格迁移或从零长跑。
- 独立冻结场前向/回转迁移评分覆盖 `641×321 → 961×321`：核心 X/Y `35.1867/36.8000`、前沿 `36.3926`，相对峰跳 `7.32e-5`、守恒缺损 `3.18e-15`、值域越界 `1.42e-9`。`961×481` 的相对峰跳 `7.638e-5`、守恒缺损 `3.535e-15`、值域越界 `1.921e-9`。这是无 LU 的研究评分，不替代生产 `mesh.build` 后的原生事务。
- 单次规划复用相同参考 X/Y 轴的纯几何试探；四组真实冻结态的候选、质量报告和审计逐值一致。1000000 节点预算的 19 个可扩展级别，完整规划计时从 121.365 秒降至 62.345 秒；310000 节点预算从 7.820 秒降至 5.659 秒。逻辑试探数仍逐级登记，缓存仅在同一来源和政策的一次规划内有效。
- 服务器接口无 LU 回归通过，覆盖 version-4 默认、version-5 1000000 节点预算与大箱初始观察；旧配置测试 `autonomousMeshConfig` 和旧可变节点账本 `autonomousControllerV2` 通过。version-5 原生 H8 `t=0` 三步到 τ=.004、原生 checkpoint 验证以及再续三步到 τ=.008 均通过。完整数值套件、version-5 演化重布和长跑仍待验证，不能据此发布为已完成的强奇异性 Profile。

## 2026-09-10：服务器接口与 version-4 自动网格

- 加入面向服务器的 `ipm.config.longTimeProfile(settings)`、
  `server/profile_settings.m` 和批处理启动脚本。公开设置
  `maximumAdjacentGridRatio=2` 精确映射到原 `remeshMaximumCellRatio`，在相邻网格比
  不超过带舍入容差的 2 时不触发 `grid_smoothness_failure`，其他数值安全停机保持。
- 将已完成全套数值回归的自动网格 version 4 分层/方向增点实现合入维护路径，保留
  已验证的 rounded-axis 浮点停滞提前返回。新增独立服务器接口配置与停止边界测试；
  新接口/停止边界及大箱观察配置实测通过；合并后的 `ipm.verify('all','quick')`
  与六算例/三迁移 `equivalence` 均实际 exit0。首次接口测试因测试自身把
  `timeIntegrator` 误读到 transport 域而失败，修正测试字段后通过；数值实现未因该失败修改。

## 2026-09-09：自动网格构造的逐值提速

- 仅在 roundedAxis 二分中点已等于端点时提前返回，原数值返回值、轴和候选选择保持。400实际构造×3长度单位的1200组完整比较通过；6轮ABBA构造器速度比1.31–1.45、中位1.415，不声称整段PDE同幅提速。
- 完整原70/分层330候选、500耗尽、单位协变和真实不同N等noLU回归62625；核比较83422；all quick与6算例/3迁移新旧一致性12177，均实际exit0。维护112数值文件与已测套件副本逐字一致，仅roundedAxis一文件合入，运行冻源未改。
- 证据：result/verification/autonomous_runtime_20260909/rounded_growth_stagnation_candidate_v1。新增research/longtime_lab/AUTONOMOUS_MESH_DESIGN_SUMMARY.md归纳从零全自动到长时的要求、实际失败与未验收边界。

## 2026-09-08 午间续算与独立误差检查

- 后续v3 stage6已到tau6.093201945464、物理t1.915696770059、G14.298746917185，
  达3.03694 e-folds但目标继续活动。真实target42重网格/短程空间对照已通过，
  old/new末态core约20/20与46/48、梯度峰相对差4.8437e-5；原线性lab差和
  新的双向Hermite辅助差同时保留，未修改空间门或声称独立初始网格收敛。
- 五支q512短箱协议实际完成，full/crop的自然与公共epsilon均通过；.001物理
  区间内梯度最大相对Inf差约4.23e-6，c_l差8.82e-5。全箱原生/新初始化协变
  已实测，仍只支持短区间晚态敏感性，不代表先前历史或无穷箱极限已验证。
- 固定候选的分块精确积分解释了旧双线性采样的虚假残差恶化；Hermite训练
  积分范数仅改善约.16%，小于native限制敏感性，尾部速度仍超门，拒绝保持。
  真实主线1427行采样将峰值漂移定位于已记录中心变化/重网格区间；中心未变
  区间的净漂移约7e-14，瞬时PPrime约3e-17。连续峰值定义继续独立研究。
- 新增保留原幅度的物理密度/梯度与重标度导数尖化图，真实原生8帧验证并查看
  图片无布局问题。tau3.002→5.693核心物理宽缩约10倍、峰曲率长约31倍，
  图中未拟合奇异时刻或极限指数。修正Lawson算法在两档可信小网格光滑段
  二阶，跨模板切换失败完整保留，未改主线CFL/数值公式。
- 冻结v3主线已到tau5.493201945464、物理t1.882633244567、G10.895228899245，
  增长2.76510 e-folds，全部历史可信。下一次移动平台事务通过，继承累计绝对峰跳
  7.45208591e-4；继续推进，尚未达到强奇异性/无穷时间极限验收。
- 从tau5.093配对场构造同箱同N的局部分辨率候选，32/42目标预测core分别
  35.2/36.8与46.2/47.6651，所有冻结场门通过；动态敏感性待测，不称全局网格收敛。
- 真实三帧Poisson保存解的逐行后向误差检查已通过：原全局残差逐值复现，
  最大分量误差约2.79--3.24e-16，没有超过1e-12的行。保持RCOND告警，
  该结果不是正向误差界，不替代空间/有限箱验证。
- 固定q512、同remeshCount3的商空间secant候选被拒，训练L2残差乘1.33452，
  尾部速度变化.00104319也超门。真实失败rho另存，不把回退baseline当候选。
  线性C1 Hermite观察器已通过非均匀三次多项式、真实幅值导数时间差分、跨cell
  连续性与加密测试；真实候选复核仍排队，原拒绝完整保留。
- 新物理fresh箱长分支从配对物理采样重新初始化参考与局部时钟，并显式记录
  绝对epoch，不伪造旧历史。5支小网格协议含原生/full/crop及公共WENO epsilon
  对照实际通过；每支终态原生CP严格配对。真实q512箱长协议另行运行，未晋升。

## 2026-09-08：长时间实验、冻结场网格与换网格内存优化

- 从原 schema4 q512 step1818 以不可修改源码副本续算两个 Delta tau=0.2 段。
  已到 tau=4.073201938625、物理 t=1.755056773798，step1936，physical max|rho_x|=
  5.519022127615。两段原历史逐值前缀保持，原生 checkpoint 与末态 rho 精确一致，
  全历史可信；末态双轴 strict core=16.5468/16.3783，safety=0.48845，进入主动换网格。
- 独立网格实验修正 H=1e6 下监视器默认最小宽度与实际窄核不匹配的问题。
  在 tau=3.6732 的冻结场上，保持 ratio<=1.08、curvature<=0.01、正求积及局部权重门，
  得到 strict core24/24 且精确保留 x=+-1 的候选；rho_x 峰跳 1.21913e-4、
  field relative L2=4.53458e-6。该初态候选的动态对照因并发内存压力中断，未记为通过。
  tau=4.0732 的同参数候选失败于 x 可行性；target21、sigma1的新候选通过全部硬门。
  真实 transaction 得core21.0294/21.8734、峰跳1.92799e-4；随后旧/新网格
  同物理时刻3/6步对照通过，物理时间差5.33e-15，最大梯度相对差2.29345e-4，
  lab窗口rho/omega/rhoY相对L2差为.000369314/.00127881/.00135956。
  自适应冻结源码主线首段已至tau4.28320194454、t1.77985698442，最大梯度6.1224393；
  第二次重网格target21、sigma.75通过，累计绝对峰跳2.91738e-4，继续向tau4.4832推进。
- `ipm_gridlab_regrid_checkpoint` 在轴校验后、复制 oldState 前释放旧
  `ops.poisson` decomposition。维护中的 transfer、flow、审计、记录和签名表达式不变，
  避免旧/新 LU 同时存活。二阶 PCHIP、四阶、六阶物理模式及四阶精确二次峰动态模式，
  各 identity/微变形共8例：完整 payload/signature、除创建时间外的 checkpoint、
  完整 audit、直接 transfer 场及 metrics 全部逐值一致，源 checkpoint 不变。
  这是小网格等价与对象生命周期证据；尚未测量 q512 内存节省或声称完整生产提速。
- 独立 WENO 分块测试在两组各100项比较及 q512 尺寸中逐位一致，单线程内核收益
  约5%--17%，未集成生产。线程数试验第一轮收益未复现；生产保留10线程。
  历史 checkpoint 的数值签名含线程相关浮点归约，改变线程数会严格校验失败；
  原10线程验证通过，未改签名算法、未降低校验标准。
- 原 schema4 q512 最后5帧的独立 secant 外推 memory1/2/4 均未被采用：某些平均
  RHS 残差降低，但局部梯度 RHS 无穷范数反升（代表例1.812倍）。规范峰值历史
  漂移8.8471e-5的原始拒绝保留，历史正幅值投影明确标为研究候选、不是物理续算。
  双轴窄核 log-width 斜率约-0.652/-0.651，后续检验平移/收缩与形状漂移的区别。
- tau4.0732单态局部平移/双轴收缩拟合，在独立外圈将梯度RHS L2残差降至
  原值约10.2%、无穷残差约9.49%（严格二次峰平移条件）；全域无穷残差反升
  88.9倍，因此仅为内层坐标诊断，不作为全域C规则或已收敛Profile。
  精确内层坐标导数已通过非均匀多项式/Gaussian时间差分、非零幅值导数、模板
  切换拒绝与tiny原RHS测试；实际q512精确坐标对照另行记录。
- paired-grid profile与冻结截箱诊断的tiny v3实际通过（Green边界，65x33/67x35）。
  原生10线程校验后恢复调用线程；全箱核心速度差为0，截箱能检测非零变化，
  过小截箱被拒绝，两个入口Code Analyzer无消息。尚无晚态箱长/分辨率验收。
- 晚态箱长分类wrapper额外tiny实际通过：full-box零差、敏感crop拒绝、无动态晋升、
  调用线程恢复；早先三个MATLAB启动失败已独立保留。轻量HDF5监视器与已原生
  校验tiny checkpoint的7个MATLAB标量逐值一致，但它本身不校验signature/可信历史。
- endpoint CFL只读定位与维护selectTimestep的率逐值一致：tau4.4832全域rate251.225，
  双峰各±3个横核宽/3个竖核宽内rate29.022、壁面9.257；最紧点在X1.14564、
  Y9.70213，source仅为峰值8.56e-7。支持研究跨y延伸的窄x网格所造成的输运限制，
  不授权放宽原全域CFL。该诊断不建立LU、不推进PDE，Code Analyzer无消息。
- tau4.2832固定内层节点的真实晚态静态箱长屏已完成：H约1e4（901x386）
  核心速度relative L2/inf变化1.37981e-4/1.20419e-4，c_l相对变化8.79963e-5，
  全部静态门通过；H约3k也过，H约1k速度门失败。没有动态箱收敛结论或主线箱修改。
- 自适应主线stage3已到tau4.68320194454、t1.82058739639，G7.4360598；stage4
  注册候选族耗尽后安全退出（原生step2226保留）。新增更强网格族仍在研究，
  同网格Delta tau=.1的安全bridge段另行启动。原生接续入口新增CP/配对result
  识别与逐值配对门，静态检查通过，尚未将该新版作为实际主线运行。
- 两个同网格bridge短段实际完成至tau4.88320194454、G8.1802859。移动细平台族
  释放精确x=1的固定节点编号限制，原位置/参考与硬门不变；无竖向余量的首个
  实际事务yCore19.453<20遭拒，未提交CP。加15%竖向候选余量后，真实
  core23.0883/22.4604、峰跳1.679586e-4；同物理时刻old/new各5步通过，
  时间差1.60e-14，峰差1.96096e-4，三场lab L2最多.00124213。
  已用186文件冻结源码adaptive_source_v3实际重启adaptive_campaign_v3，
  step2317、tau4.893201945464起步向5.0932推进；v2快照未启用。
- step2226晚态真实时间步加密粗/细6/11步通过，时间差1.20e-14，三场lab L2
  至多1.66e-12，1174行源前缀/完整可信/原生终点CP均通过。背景alpha(Y)X拆分
  的剩余率约为原率的1/10；高阶背景remap已做合成真解/半群/逆向/缺失入流测试，
  未把它当成已完成的非线性大步算法，未改变原CFL与+ipm。
- 独立row-equilibration小测未证明真实精度提升，不集成。dense Schur/Sylvester
  在admitted小网格比较通过，但raw高拉伸矩阵的全局残差5.50e-14仍可伴随
  已知解相对误差1.90e-4、分量后向误差3.99e-6，因此拒绝外推至q512生产。
- 产物位于 `result/longtime/20260908_campaign_v1/` 与
  `result/verification/performance_lab_20260908/`；新研究入口在
  `research/longtime_lab/`、`research/acceleration_lab/`、`research/performance_lab/`。
  尚未达到强奇异性或无穷时间极限证据，目标继续活动。

## 2026-09-05：schema 4 精确规范与无反馈重标度契约

- 新配置升为 schema 4，并固定
  `scalingContract='exact_gauge_no_feedback_v1'`。`lengthScaleGain` 固定为 1，
  `omegaGaugeGain`、`widthGaugeGain`、`travelingWaveGain` 固定为 0；规范只解
  瞬时代数约束，不再叠加幅值、宽度或位移 restoring 项。
- 删除 `maxDynamicRate`，退役 `adaptiveGain` 运行组件和
  `adaptiveLengthScaling`、`widthExpansionStrength`、`widthContractionOnset`、
  `widthContractionStrength`、`maxWidthRateCorrection` 五个 width-controller 选项。
  峰值、层级、面积和连通宽度
  的网格单元计数仅作 telemetry、remesh 和 hard-stop 输入，不得修改
  `c_l/c_omega/c_r` 或时间速度。
- checkpoint 读取器也将 `state.rescaling.maxDynamicRate` 作为 schema-4
  禁止字段显式拒绝，并由负测试覆盖；旧 limiter 不能借受污染的续算产物重新进入。
- 精确幅值规范接口包括 `anchor_wall_template_projection`、
  `anchor_bulk_gradient_l2` 和 `anchor_wall_window_l4`；它们的 reference、
  condition、energy 和 rate residual 属于可审计的规范证据，不是反馈增益。
- q256 规范竞赛的每个结果还必须在完整连续可信历史上证明
  `c_l=c_lNominal=canonicalCLNominal`、`lengthScaleGain=1`、
  `widthControlMode=0`、两种 width correction 恒为零，以及固定 `X=1` 处
  `V_1(1,0)=0` 的规范化残差不超过 `1e-10`。候选条件数、残差与恒等式门槛
  来自结果内冻结的 candidate specification，不能由后续注册表改动重解释。
- q256 group-stage 在首个求解前冻结完整候选、网格、物理协变性、checkpoint
  cadence 与时间预算；正常安全停机、数值/网格门失败、规范退化及意外崩溃分开
  归类，异常后从不可覆盖的带步数 checkpoint 中定位最新存档。
- 预注册 q256 五候选正式赛程在 `2337.84 s` 内完成，五例均到达 `t=0.55`
  并通过全部硬门和 10/10 对称物理协变比较。字典序排名为
  `gradient_energy`、`anchor_wall_template_projection`、
  `anchor_wall_window_l2`、`anchor_bulk_gradient_l2`、
  `anchor_wall_window_l4`；对应尾段最坏 `c_l/c_omega` 相对漂移为
  `3.904%/11.799%/39.964%/60.696%/67.944%`。物理解最坏 rho/omega
  relative L2 差异仅 `1.510e-12/3.224e-12`，但共同远边界比仍告警，
  因而该结果只作无反馈规范筛选，不作爆破或 q512 可信性声明。
- 新增局部 `wall_omega_quadratic_peak` 精确幅值规范：在 wall `R_X`
  活动极大节点的非均匀三点二次顶点上，用同一插值求值算子
  得 `P=Q[R_X]`、`F=Q[B_X]`，直接取 `c_omega=-F/P`。
  `transport_anchor` 仍严格取 `c_l=-U_1(1,0)`；`tau` 直接积分、
  `timeSpeed=1`，无 limiter/平滑/宽度反馈。`none` 只是预设
  `c_omega=0,C_omega=1` 的负对照，不把恒零速率当成经验收敛。
- q512 初始 RHS box screen 在 `H=281,1e3,1e4,1e5,1e6` 上确认
  `gradient_energy` 的全局能量/forcing 对非衰减 datum 随箱子线性增长，
  因 extensivity 被否决。局部二次规范在 `H=1e6` 给出
  `(c_l,c_omega,P)=(0.488462559170,0.125042865249,0.686073738019)`
  并通过所有硬门；但 far-velocity ratio `0.986819` 使其仍只是初始
  normalization screen。不可变证据位于
  `result/verification/q512_comega_initial_box_screens_v4/screen_20260905T110530211Z_tpbe425904_b6b4_4b3f_9aa4_c02b888cd290/`。
- `wall_omega_quadratic_peak/none x CFL=0.25/0.50` 的 q512 `H=1e6`
  四组 fresh-root 短算全部 hard pass。`CFL=0.50` 将每个算例从
  40 步降至 21 步，总墙钟加速约 `1.74x`；同规范 CFL 比较最坏
  rho/omega relative L2 为 `8.003e-12/2.963e-9`，因而冻结
  `CFL=0.50` 用于长算。证据位于
  `result/verification/q512_fresh_gauge_cfl_short_screens_v4/screen_20260905T114636013Z_tp3e066220_eec8_4c94_a0c7_b02a376fc67f/`。
- q256 二次顶点续算至 `t=1.7` 时，最后 `0.125` 物理时间窗内的
  inverse-peak affine 斜率/代数平均斜率为
  `-0.748374/-0.747521`，但 `c_l/c_omega/kappa` 仍漂移
  `19.4%/14.2%/3.08%`，且只有 `1.861` e-folds。因此只记为
  Type-I candidate，不声称 blowup 或渐近常速率。
- q512 fresh-root/continuation 现在在根作业中支持
  `gradient_energy` 与 `wall_omega_quadratic_peak`，但同一链上 gauge 必须
  逐 checkpoint 冻结。真实二步根作业与再两步续算已 hard pass；
  续算终态为 `completed_endpoint`、segment/total steps `2/4`、
  `parentHistoryExactPrefix=1`，证据位于
  `result/verification/q512_schema4_continuations/continuation_20260905T124215721Z_tp1b1f9942_ef22_4494_b5f7_ee4ca7b4a77c/`。
  早期两次 postrun rejection 是空 `anisotropic` 组触发的
  `fieldnames/isfield` 空形状 artifact-audit 边界缺陷；现已改为逐字段
  exact-prefix 检查，不是 PDE 数值失配。
- 新增独立 q256 `c_l` 赛程：冻结 `cOmegaGauge='gradient_energy'`，预注册现有六条
  精确 `lengthGauge`，保存完整配置、公式、运行成本、连续可信前缀和 canonical
  末尾 25%/40%/60% 速率统计，并对成功候选作全部双向物理场比较。裸调用只跑
  `t=1e-4` smoke；`t=0.55` 六候选正式赛程要求显式 `runMode='formal'` 与
  `confirmLongRun=true`。共同 CFL/`maxDt` 可整体覆盖但不能因候选而异。
  六候选 smoke 在 `21.33 s` 内全部通过精确合同及 15/15 pair，最坏 rho/omega
  relative L2 为 `3.884e-6/7.958e-6`；两步历史不足 16 点，故未产生速率排名，
  正式赛程尚未运行。三个新增研究函数 Code Analyzer 均为 0 条消息，formal
  preflight/长算确认门/非法 CFL 负测均通过。
- 生产 checkpoint 升为 schema 4，只接受 config schema 4 及完整运行期 reference
  白名单。schema 3 及更早结果只读分析；旧 checkpoint 不可续算、转换或隐式升级。
- 同步根 README、STRUCTURE、模块 INTERNAL 和研究笔记；历史公式保留作证据，
  但明确标为 retired historical contract。`tests/source_map.json` 保留
  `ipm_adaptive_gain.m` 的来源哈希，将其维护目的设为 `null`。

## 2026-09-05：双轴冻结场网格实验室与可续算换网格事务

- 冻结数据集升级为 schema v2：每个样本同时保留重标度 `rho/x/y` 与物理
  `physicalRho/physicalX/physicalY`、尺度时钟、可信标记和场签名；每个来源保留 metadata、
  完整配置、网格签名、节点数和停机原因。schema-v2 特征同时记录
  水平/竖直 source 包络、物理 `rho_x/rho_y` 最大值及 0.1/0.5/0.9 连通宽度；
  0.9 分量作为续算的严格 core。
- 水平解析监视器使用非对称代数尾 core、窄左前沿和宽 bridge；竖直监视器
  以 `y=0` 为边界 core，显式关闭虚假内部前沿并以宽 bridge 覆盖 source/
  envelope 尾部。新增单边 y 轴构造及与对称 x 轴一致的比率、log-spacing
  曲率、stencil 条件数和正求积权门。
- 新增全局 log-cell-width equalizer，把特征压缩分摊到大量相邻单元；它在固定
  endpoint 下搜索满足 core 目标的最小可准入幅度。网格比率上限是硬门，
  equalizer 和 stage designer 都不会为凑够点数自动放宽；不可行时返回
  在该上限下的最佳可达点数与限制门。
- 新增双轴 pair scorer：仅对 x/y 均通过网格准入的候选执行维护中的
  old-new-old 高阶迁移，对每个冻结时间审计场、壁迹、`rho_x/rho_y`、前向最大值
  变化、质量、值域、对称性以及水平/竖直严格 core 和前沿点数。
- 新增 `ipm_gridlab_design_continuation_stage`：从当前可信末态构建 x/y 候选，并总是
  将它们应用到外部 schema-v2 全局参考轴，避免逐段在当前网格上叠加局部
  密度比。默认目标为 x/y 的 0.9 core 各 18 cells、x 前沿 8 cells，
  adjacent-ratio cap 为 1.08；这些是网格准入目标，不是可信爆破证据。
- 明确区分两种 checkpoint：`ipm_gridlab_make_checkpoint`/
  `ipm_gridlab_restore_checkpoint` 是旧的实验室 schema-v1 离线状态包，不是生产动态
  restart；`ipm.output.*Checkpoint` 是完整保留
  history/cursor 的 schema-v2 可信接受步事务，由唯一入口
  `ipm.solve(overrides,checkpoint)` 同网格恢复。
- 新增 `ipm_gridlab_regrid_checkpoint`：在不推进 PDE 时间的前提下，用维护中的
  remesh transfer 将可信生产 checkpoint 换到等节点、同 endpoint 的新网格，重建
  operators/flow，保留尺度、时钟、初值不变量、history/cursor，并在相同接受
  时刻替换末行。默认事务门要求 x/y 严格 core 各至少 16 cells、
  `max|rho_x|` 相对改变不超过 `5e-3`、质量缺损不超过 `5e-12`、相对值域越界
  不超过 `2e-4`；这些门与网格质量审计全部通过时才产生新的生产
  checkpoint。
- 新增独立研究驱动 `ipm_gridlab_run_q512_campaign`，用全局参考轴自动循环
  design -> 事务式 regrid -> `ipm.solve`。默认每段 `Delta tau=0.35`（硬上限
  0.42）、可信记录态 checkpoint cadence 0.20、ratio cap 1.08、x/y 严格
  core 设计目标 21/20 和 x-front 20。硬网格门包括双轴严格正权、
  `0.35<=q_i/h_CV_i<=1.65`、相邻比不超过 1.08、log-spacing 曲率不超过
  0.01，以及全局 `min(q)/mean(q)>=1e-8` 兼容底线。原 `1e-4` 仍作为独立
  telemetry warning 阈值，低于它会在 design/regrid/terminal 中显式标记，
  并非被忽略或删除。
  Stage 3 -> 4 的全局 y 比由 `1.526e-4` 下降到 `1.084e-4`，而局部
  `q_i/h_CV_i` 比率保持稳定且在硬区间内；因此前者诊断的是双尺度全局
  span，不是局部网格突变。
  campaign 直接使用 `ipm.mesh.quality` 的同定义局部求积字段，
  记录实际选中的全局 reference provenance，并显式拥有可选的
  uniform flat-tail fallback 参数，不隐式继承 stage designer 内部默认值。
  每段末强制连续可信前缀、可信终态、双 core
  至少 14 且 safety 严格小于 0.70。换网格 `max|rho_x|` 单次跳变不超过
  `2e-3`，并按整条 provenance 以 `2e-2` 上限累计守门；result bridge 的直接
  使用与继承血统均写入不可覆盖 MAT summary 和 JSONL manifest。新阈值策略
  已实现并完成文本审阅；Code Analyzer 与动态验证仍待父任务运行，本记录不预先宣称它们已通过。

## 2026-09-04：可信接受步 checkpoint 与同网格恢复

- `ipm.solve(overrides,checkpoint)` 仍是唯一时间推进入口；新增可选 canonical cadence/
  exit checkpoint 策略，且只在最新记录属于连续可信前缀时落盘。
- schema-v2 checkpoint 保存完整 `rho/x/y`、参考轴、内部尺度、重标度参考状态、初值
  质量/值域、标量历史与 `nextOutput/nextCheckpoint` 游标；`storeSnapshots=false` 不影响
  可恢复性。临时 MAT 原子安装为带 case ID/step 的不可覆盖文件，并追加独立清单。
- resume 只允许终止界、结果/图像/视频/元数据/checkpoint 控制覆盖；同网格动态二阶和
  129 x 65 四阶短程 split-run 均逐位等同于不中断运行。不可信终态不产生 checkpoint。
- 新增 `checkpointFromResult`，把 terminal-trusted v2 结果受控转换为同网格 checkpoint；
  转换会重建并 overlap 核验 `ops.rescaling/scale/mass0/rhoRange0/flow`。旧 result 未保存
  内部 `logC`，故该桥接只承诺机器精度一致，不宣称 bitwise exact。
- 已将本次 1025 x 513 q512 可信终态转换为 step 5300 的生产 checkpoint；完整初始和
  末态 overlap 检查通过，未推进额外 PDE 时间步。

## 2026-08-31：建立结构化副本

- 从保留的 `../sixth_order_integration` 建立独立副本；来源提交为
  `65f999b7e205dd4decef06c937d1aa003cf8b744`。
- 将原来的扁平函数按职责放入 `ipm.config/mesh/field/evolve/remesh/diagnostics/output`，
  求解入口改为 `ipm.solve(opts)`；不保留第二个求解包装入口。
- 以 `ipm.verify` 调度基线、四阶、六阶及新旧一致性验证；短算例与研究实验分开。
- 新建结构总览和各模块 `INTERNAL.md`；原有技术说明与研究历史保存在 `research/`。
- 不修改数值公式、选项语义或 v2 结果契约；不复制原目录 `.git` 与生成的批量产物。
- [source_map.json](tests/source_map.json) 记录原版全部 140 个已跟踪文件的 SHA-256 与目的路径，
  以及 103 个函数的新旧名称；原版 `AGENTS.md` 的副本以本目录维护约定替代。

## 2026-09-01：验证完成

- 已确认：基线全部 11 层、四阶完整套件（quick 稳定性扫描）、六阶完整套件通过。
- 已确认：11 个短程求解和 3 个非平凡网格迁移的新旧数值逐项完全一致；
  求解结果仅排除每次运行生成的 `metadata.caseId/createdAt`。
- 已确认：原版 140 个已跟踪文件的 SHA-256 保持不变。
- MATLAB R2026a Update 5 下，103 个迁移函数解析及 75 个求解依赖检查通过。
- 全部 125 个 MATLAB 文件 Code Analyzer 零消息；三个短算例及终态/实时绘图通过。
- 干净会话只添加新根目录即可运行；验证结束后恢复原路径和工作目录。

完整记录由 [tests/VALIDATION.md](tests/VALIDATION.md) 汇总；四阶 heavy 扩展扫描未计入上述通过项。
`research/notes` 中的历史通过记录不视作本副本已通过。

## 2026-09-03：四阶 active datum 最大值增长率审计

- 新增 `ipm.diagnostics.maximumGrowthRateFit`：在物理时间上估计
  `d(log M)/dt` 与其导数，并以时间加权、固定等距评分点和 e-fold
  窗口比较简单指数、有限时幂律及双指数。
- 修正兼容 `blowupFit` 的有限时刻搜索上界，并显式记录 profile 是否
  命中边界，避免把指数极限误报为有限奇性时刻。
- `ipm_assess_fourth_order_blowup` 现同时分析最大梯度和跟踪壁峰，并用
  动态壁峰规范给出的独立物理增长率核验数值微分。
- 解析指数/幂律/双指数在均匀与强非均匀采样上的回归测试通过；完整
  baseline 和四阶算例设置审计通过。当前解析 256 x 256 单元结果支持
  加速增长及有限时幂律候选外推，但不足以支持可信爆破分类。

## 2026-09-04：384 x 384 空间加密与增长曲线

- 新增 `analytic_moderate_narrow_quadrant_384`：769 x 385 存储节点
  对应第一象限 384 x 384 单元，两个中心目标间距同比缩为 256
  算例的 2/3，并保持原四阶数值组合与全部质量门。
- 新增 `run_fourth_order_dense_growth`，自动保存原始结果、连续可信前缀
  assessment 及最大值增长图；重复运行时后续产物从求解器返回的
  实际唯一结果路径派生，避免文件错配。
- 新增单分辨率六面板增长图和 256/384 分辨率对比图，均严格截取
  `continuousTrustedPrefix`。
- 384 算例在 2982 步后于 raw 物理时间 `1.855940045` 因
  `grid_resolution_failure` 停机，可信前缀到 `1.751582283`，比 256
  算例延长 5.533%。公共可信时段的最大梯度差在粗网格终点仅
  0.079%；但动态范围仍只有 1.998 e-fold，因此结论仍是空间一致的
  有限时幂律候选，不是可信爆破解。
- 新配置静态检查、网格准入审计与完整四阶 quick 证书通过。
- 按 canonical `tau_c` 增加线性最大值、指数最大值与曲率守门
  模型，并新增 256/384 的 `tau_c` 六面板双尺度图。q384 最后
  `Delta tau_c=0.5` 上梯度/壁峰指数率为 `0.51635/0.50317`；
  线性模型的留出误差约高一个数量级。
- `ipm_assess_fourth_order_blowup` 现保存物理最大值、重标度梯度、
  `rho_x/rho_y`、尺度比及水平/竖直内宽的 canonical-time 诊断；
  384 长跑驱动器也会在 q256 基线存在时自动输出 `tau_c` 对比图。
- 明确 `dt/dtau=C_omega/C_l` 与
  `M_phys=(C_l/C_omega) M_R`：物理时间的近 Type-I 幂律与
  canonical `tau_c` 的指数轨道是同一候选机制的两种表示。
  壁峰关系部分由规范恒等式保证，所以不重复当作独立爆破证据。

## 2026-09-04：Lorentzian LAT 平滑化与 512 x 512 计算

- 为减弱固定非均匀网格的局部极端变形，在 256 x 256 第一象限粗网格上先筛选
  LAT 监视器，并选定中心 `1.14`、宽度 `0.20`、强度 `A=28` 的双
  Lorentzian 设计。水平最大相邻网格比/最大 `log(dx)` 曲率由旧 q256
  Gaussian 映射的 `1.301284189/0.0741483378` 降为
  `1.10097818/0.00240067378`。
- 优化后 q256 在 2732 步后因网格分辨率门停机，raw/trusted 物理时间为
  `1.83683959412/1.72352016681`，trusted canonical `tau=3.80105633`。
  q512 存储 `1025 x 513` 节点，对应第一象限 `512 x 512` 单元；它在
  4800 步的 `maxSteps` 预算上限停止，956/956 个输出记录全部可信，不是
  数值失效。raw=trusted 物理时间为 `1.83433595901`，canonical
  `tau=4.77268406`，终点 x/y core 为 `8.28986/12.1711`，safety 为
  `0.965034`；峰值 RSS 为 `5.423628288 GB`，swap 为 0。
- q256/q512 公共物理时段终点的梯度/壁峰相对误差为
  `8.719e-5/4.207e-5`，log-RMS 误差为 `2.327e-4/1.660e-4`。
  q512 物理时间 1.5 e-fold 幂律拟合给出梯度
  `T=2.03899451,p=1.12989826` 和壁峰
  `T=2.03623379,p=1.09028074`。最后四分之一 canonical-`tau` 区间的
  指数率为 `0.5059/0.4913`，指数拟合误差
  `0.0012/0.0016` 低于线性拟合的 `0.0129/0.0117`；但完整
  1.5 e-fold 的 canonical 判定仍是 `curved_nonasymptotic`。
- 梯度/壁峰总动态范围只有 `2.382/2.415` e-fold，低于 3 e-fold
  证据门；远场 source 最高占比 `2.716%`，远场 velocity 占比约
  `94.9%`，且缺少晚时半步长、更大计算盒与第三个匹配网格复核。因此
  细化审计保持 `credible=false`：结果支持空间一致的双时钟候选机制，
  不支持夸大为可信爆破解。
- 新增 `run_fourth_order_lattice_screen`、
  `run_fourth_order_lattice_validation` 和 `run_fourth_order_lattice_512`，分别负责
  q256 短程 LAT 筛选、选中映射的 q256 长程验证以及 q512 最终计算与
  产物归档。新增 `ipm_assess_lattice_refinement` 用于公共时段细化误差和保守可信性
  审计，`ipm_plot_canonical_maximum_models` 用于在共同 e-fold 窗口比较 canonical
  线性/指数最大值模型及增长率。

## 2026-09-04：相邻密度比均摊的协调 LAT 与严格嵌套 q512

- 新增 `coordinated_log_spacing` v5：固定峰路径细平台后，在 q256 正半轴全部
  121 个可行整数拆分中最小化左右支的最坏 `abs(Delta log(dx))`，同值时再最小化
  两侧斜率差。选中的 F92/R20 拆分为 `L/F/R=50/92/114`，左右渐近斜率
  `0.085282394/0.084378533`，平衡度 `0.989402`，最大相邻单元比
  `1.089024557`。162 个非零 active-transition 接口的均摊效率为 `0.882370`；
  另有 93 个平台内 `ratio=1` 的设计性零跳变，不混入 active 均摊结论。
- q256 的显式 92-cell 核心平台为 `[0.95,1.30]`，连同两侧等宽桥接单元后，
  实际连续平台为 94 cells、范围 `[0.9461956522,1.3038043478]`。q512 由每个
  q256 单元保宽二分构造，粗节点嵌套误差为浮点舍入量级，实际连续平台为
  188 cells；最大相邻比降至 `1.043847442`，最大 log-spacing 曲率为
  `0.001349049`。
- F92/R20 q256 长算在 raw/trusted 物理时间
  `1.88073510666/1.76292414513` 停止，对应 canonical
  `tau=5.368011277/4.096011277`。相对旧 A28 q256，可信物理时间延长
  `2.286%`、可信 tau 延长 `7.760%`；共同物理时刻的梯度只差 `0.039%`，
  横向 core 点数增加约 `42.2%`，失效瓶颈由 x 转移到 y。
- 严格匹配的 q512 用时约 7 小时 13 分，在 5300 步预算上限停止；1009/1009
  个记录全部可信，raw=trusted 物理时间 `1.85662160032`、canonical
  `tau=5.03835845624`。终点梯度/壁峰为 `10.2985731/8.70832525`，增长
  `2.51216/2.54105` e-fold；x/y core 为 `10.1425/8.23707`，safety 为
  `0.971219`。相对旧 q512，物理时间和 canonical tau 分别延长
  `1.215%/5.567%`。
- q256/q512 的 domain、数值组合、映射来源及网格嵌套全部匹配。在共同可信
  区间 `[0,1.76292414513]`，梯度 log-RMS/末端相对误差为
  `2.048e-4/7.614e-4`，壁峰为 `1.585e-4/3.429e-4`。q512 的物理时间
  0.5--1.5 e-fold 窗口继续选择有限时幂律，梯度约
  `T=2.039,p=1.13`；canonical 最近 0.25--0.5 e-fold 选择指数增长，
  但 1.0--1.5 e-fold 仍为 `curved_nonasymptotic`。
- 可信爆破判定保持 `false`：梯度动态范围仍不足 3 e-fold，终点远场速度比
  达 `0.9850`，且缺少半时间步、更大计算盒和第三个匹配网格。现有结果支持
  “物理时间有限时幂律 + canonical tau 短尾指数”的双时钟候选机制，不能当作
  已收敛的标准 IPM 爆破解。
- `ipm_plot_lattice_spacing_balance` 现明确区分 active-transition 与平台零跳变；
  设置审计直接回归 q256/q512 节点嵌套和 188-cell 连续等宽平台。


## 2026-09-08 05:29 UTC — 主线tau7.0932与连续几何/较长空间检查

- v4从tau6.6932经112/.90平台事务继续两段到step3326、tau7.09320194546、physical t1.95452713669、G22.09532515。全trusted、源前缀、末态原生配对通过；stage003全部注册mesh family耗尽，安全退出。原CP保存，目标继续active。
- 当前轴非单调搜索在v3stage009 sigma.25发现维护搜索漏过的可行幅度，独立完整pair通过；最新v4stage002 x已可行但root y sigma1不足21，未提交不ready候选。
- target42长空间支152步/709.60秒完成并匹配主线t1.86943578867，细端core40.495/44.769、safety.19756，G差1.1406e-4、双向Hermite梯度L2约7.47e-4；全部原生/前缀与补充CP-axis逐值检查通过。共享早期粗历史，不称独立网格收敛。
- 连续C1同峰横纵宽研究观测通过合成/单位协变，消除旧nearest-column纵宽约3%观测扰动。11帧观察/滚动预测已实测，双观察密度一致；固定早期拟合到最新仍胜线性，但滚动衰减率不稳定且有下界命中，未认定极限/未推进预测场。所有旧负结果保留。


## 2026-09-08 05:47 UTC — 原算法到tau7.2432，最新缩箱及宽度观测审计

- bridge003固定原网格78步/378.07秒到step3404、tau7.24320194546、physical t1.9590400984、G23.5397256；末态core16.3781/17.9260/safety.48846，原生精确配对/全可信/源前缀与既定core门通过，session退出。
- v4stage002独立[1,.1,.01]冻结缩箱原门全部通过，.01为895×386；5支.001物理native/fresh/epsilon对照在推进，前两支全箱协变过，不宣称全协议完成。
- 901×386 fresh target32真实事务及old/new5/6步可信通过，原smoke仅离散y宽差.00859686超.005，总false保留。4原生CP noLU连续峰宽度审计证明差异几乎全来自nodal峰列：同连续峰宽差2.401e-5，固定共同物理列差2.267e-6。正做独立解析观测验证；无生产门或规范修改。


## 2026-09-08 06:27 UTC — fresh细核心长推进及三盒对照、固定未来预测

- 新core patch独立实际资格通过（旧离散宽smoke=false保留，事先注册连续同峰宽aux=true）；895×386分支109步到fresh139/parent等价tau7.24320075186，G23.53790686、core32.7604/36.3166、safety.24420，全部native/历史前缀/lineage/端点门通过。与原full-box同物理t1.9590400984，G相对差7.73e-5、Hermite omegaL2约1.35e-4。旧末尾局部时钟观测报错单独修复/审计，未重算PDE或改原CP。
- 最新三natural盒长对照完整pass：full/crop各累计81步，middle新IVP79步，同绝对端点误差<=1.11e-15。full/small rhoY Inf1.679e-4，full/middle1.531e-5，middle/small1.526e-4；全箱与原native实际协变通过。long common-epsilon未扩展，dynamicBoxConvergenceEstablished仍false。
- 固定increment rank1在3过去留出帧预测误差.0141%/.0328%/.0554%，胜同窗口线性4.5--7.4倍；未来tau[7.4,9.1]模型/源码已冻结，不用新帧重拟合，不推进预测场。
- 一致C1内层几何及真实U_tau完成解析/tiny/尺度/退化测试，latest真实缓存core/外圈RMS .00778688/.01464538，只剩归一Eulerian演化3.11%/6.29%；核心绝对残差尚未下降，不能宣称突破慢收敛。C1几何最细MMS绝对误差和协变过、部分加密比例门失败，allPassed=false和所有旧失败保留，不声称普遍四阶。


## 2026-09-08 06:52 UTC — 主线tau7.6432、真实换网格及首个未来预测

- 新895×386主线完成两段365步到fresh504、parent等价tau7.64320075186、physical t1.96977603989、G27.80419522，全部可信/原生/前缀/源资格过。v1首次换网格因冻结遗漏general轴helper安全退出；保留原结果后新v2完整195文件源码及实际依赖预检，已真实regrid至core35.1964/36.7277、峰跳1.43e-5、累计.000113765，再启动12段，原PDE未改。
- 固定预测首次真正未来帧tau7.4432：rank1 L2 .00134884，linear .00398521、singleDecay .00432142；仍胜线性约2.95倍但优势下降，固定模型未重拟合。
- 默认10线程实际895×386 ABBA中WENO blocked两候选完全bitwise，整步仅.34%/.91%收益，不达5%门、不采用。独立受尾部约束的tiny JF候选约.3%残差改善并通过小实验门，非原IVP轨迹、未证明长时间提速。
- 从真实RHS补上内层形状积分漂移，瞬时时间幅宽指标随四窗口为.395至.382，不宣称空间Hölder指数。积分与密度重建差<=5.99e-7，真实方向差分导数误差<=8.43e-7；原JSON向量方向导致的测试false保留并独立修正。


## 2026-09-09 从零时间自动网格要求与长时接续（16:10 UTC）

用户将最终网格验收明确为从物理 t=0 连续稳定跑到长时间，中途无需人工暂停改网格。
已写入 OBJECTIVE_20260908_ZH.md；整理历史失败与自动化回归场景，现有晚期研究驱动不冒称已完成。
fresh v2 的12个实际段全部通过终点/原生/完整历史门，step3668、等价tau10.043200751859162、
absolute t2.0079305293351908、physical G70.729822076606、core25.00919/30.64924。
v3新freeze197文件与先前通过的admission审计SHA完全一致；重新读取v1+v2完整原生链后，
48段研究接续已启动（session60123）；首个真实事务core35.15559/37.00691、peakjump3.279234e-5，
累计绝对jump.0003249532。没有修改运行中的冻结源或主线C/输运/质量门。


## 2026-09-09 自动网格组件与因果回放（16:31 UTC）

新增 AUTONOMOUS_MESH_FAILURE_ATLAS.md，12类真实故障与从t=0验收流程，所有历史证据链接核对存在。
新增研究 mesh_autonomous_decision / mesh_test_autonomous_decision，不改生产数值路径；
current/forecast触发、同remesh趋势、冷启动一步观察、动态review时间预算、明确时钟和失败分类。
v2通过12项fixture、真实fresh742记录/6epoch及原始从零历史920记录；缩短review不修改CFL。
原v1过早forecast与端点极短dt风险由独立审查发现并修正，两个反例另行测试，冻结v1保留。
这仍是历史回放，未证明新控制器从零跑到长时间。核心集成、跨多次真实迁移及独立t0误差梯度待做。
主线新v3第1段392步/1003.10s到tau10.243200751859163、G76.10824270879445，全部实际端点门过，继续第2段。
冻结未来预测器8帧收齐且未重拟合：末帧tau8.8432 rank1/linear L2 .00756266/.01117191；
固定候选极限距离先降后增，不能作为已知极限。科学图扩展至实际tau10.0432。

## 2026-09-09 17:15 UTC — 自动网格运行期开发与4768长时接续

- 实际v3至step4768/tau10.6432/G87.9733后旧候选相邻比失效，保存最后可信原生断点。
- v4研究主线采用纯70x/4y规划，98个原数值内核文件逐字保持；17-link链重读后实际迁移通过，93443继续运行。
- 维护版新增显式可选autonomousMesh，零时刻解析初选、每接受步纯规划、外层释放LU后的实际事务及元数据账本。
  旧未启用分支不改公式或默认值。仍属开发中，未宣称从0自动跑到长时。
- 真实候选审计42项回归通过（修复NaN权重误通过）；两种321x161初值各完成0→4e-5/4步，初始样本exact。
  所有初次失败和冻结源码保留。完整合同、多迁移、CP split及all/equivalence等待本轮执行，不能记为通过。

## 2026-09-09 — 自动网格第一轮完整验证与自然容量失败

- lifecycle_v2：实际2次注入迁移、候选失败精确回退、同A分解作用及真实RHS/flow/全部ops数据重建exact，原生CP分段轨迹全部exact。
- all实际通过；原版equivalence初因已有23个只读新增日志失败，记录保留。严格对齐并独立逐值重算后6run/3remesh全过，5个损坏负例拒绝，无容差增加。
- policy 9正/30反、output/controller hook 5正/14反、controller 44反+2真实from0、plan noLU14项均通过，原始输入SHA不变。
- 从0自然压力测试无中途修改，13次重网格后在tau2.68861927因固定X节点容量停于可信step1134；未完成注册tau4。
- 主线v4已有2新终点到tau11.0432/G101.44254；版本1到长时/强奇异性目标尚未完成，开始只读/几何层验证自动加点v2。

### 2026-09-09：显式方向增点v2与真实短程恢复

预登记参考族、方向单调增点、原生迁移审计及完整跨尺寸持久化已接入；v1数值路径与读/续跑规则保持。
实际v2_smoke_and_suites_v1三项smoke/all/equivalence均通过，真实4步CP分段完整字段逐值一致，第一段物理字段与冻结v1逐值相同。
独立纯配置/几何、50项迁移合同、47持久化负例/9正边界通过。尚未据此声明自然增点或长时资格。
from_zero_tau4_v2_cap110000_v1已从原始t=0启动，全程冻结、固定小箱、显式节点cap，不手工触发迁移。
最新六帧实际Profile图至equivalent tau11.243；主线随后通过tau11.643/G125.0466。图与峰值增长均不是极限或Holder证明。
详细合同与限制见research/longtime_lab/AUTONOMOUS_NODE_GROWTH_RUNTIME_V2.md。

### 2026-09-09：真实从零自动增点与调度恢复

- v2原始t=0连续tau4成功，step1653、17 remesh、1自然X增长；20 native CP与跨增长1102→1177全状态/历史/快照重放逐位通过。全部短套件通过，仍不代表大箱/无限时误差资格。
- 第二个从零tau8运行cap210000仍运行中，已在tau5.60024实际到641×321，未改途中输入。
- 主线v4 stage6达到512步上限，step7510/equivalent tau11.83017，原stageAccepted=false保留。实际23 links/9 remesh重新审计通过；17调度反例/3正例通过。研究driver v3登记4096步预算和明确首段补齐旧目标，v5已启动；104个+ipm核与v4逐字相同，未修改任何数值门或运行中源码。

### 2026-09-09：扩大初始箱的纯解析观测研究

研究screen可显式登记有限箱列表，并另选固定解析观测格；默认原均匀观测保留。求解参考仍321×161。真实几何筛查使原先无X候选的H64/H128得到原门合格pair，H256仍失败；未构造LU/PDE或修改生产。原H32数值逐位复核通过，仅异常调用栈的路径/行号作为非数值元数据单独说明。证据及未完成native资格见INITIAL_ANALYTIC_OBSERVATION_PROBE.md。

### 2026-09-09：从零自动网格完成tau8

from_zero_tau8_v2_cap210000_v1已exit0，一次冻结运行4545步、25次remesh、两次自然方向增点到641×321，末核心27.3853/29.7275、全history可信。wall3763.11s/maxRSS3.653GB/0swaps；full physicalGradInf26.5025，尚无大箱/空间/时间独立误差资格。第二增长原生split已开始，未提前计通过。旧跨尺寸快照回归改为直接检查真实rho/x/y帧尺寸，1177/1653实际281/401帧补验noLU通过，旧报告保留。

- 2026-09-09：第二次自然增点的原生 CP2370→2481 实际 split/replay 通过；完整数学状态/log/cursor/controller及串行重建flow/physical逐位，空快照合同单独通过。τ8真实轨迹图与源/输入SHA复核完成；主线实际推进至等效τ12.2432，B大域从零对照在运行。仍未获得域/空间/时间误差资格或强奇异性结论。

- 2026-09-09: H32/H64/H128 native four-step preselected-axis probes passed. Fixed the research empty-transaction auditor and extracted a read-only preselected audit; H32 original failure retained, no PDE rerun. Registered five fresh t0 cap210000 controls after original B capacity failure. Preflight passed; maintained +ipm unchanged. See research/longtime_lab/EXPANDED_ZERO_CONTROLS.md.

## 2026-09-09：τ8与H16实际容量无LU审查

- 独立performance_lab入口仅严格读取真实CP/result并做原70×4纯几何搜索。τ8/step4545同级35X/3Y、3迁移提案；B实际失败在被回退的1521，保存1520，失败Y四族均不足core目标。另登记cap210k在保存B场上给出level4候选；未transfer或PDE，原失败不改。
- 围绕τ8原始selected base的factor3两方向纯几何通过，[3,2]可行70X/3Y，[2,3]35X/4Y。只作为未来版本/资源资格提案，未编辑+ipm、旧policy/CP不升级，不声称16GiB并行LU已认证。session19722 exit0，两个helper CodeAnalyzer0、119冻源/8输入SHA不变；早先静态告警和metadata写后路径审计假设失败保留。证据与完整质量余量见 research/performance_lab/TAU8_AND_BOX_B_CAPACITY.md。

## 2026-09-09：fresh9205 有限几何及冻结配对审计

新增 performance_lab 的 fresh_stage5_geometry/pairs/native_trial 独立研究入口，不改维护数值模块或运行冻源。原 stage005 的 70 X 全拒保留；登记 40 个 X 参数（38 新构造、2 旧失败复用）获得 10 个全质量/core/front 通过，未执行条件保护 warp。按原最弱质量余量固定排序仅评分 3 对 X/Y，3 对均过原冻结迁移门，源 step9205 的 clocks/base/runtimeRefs 不变；profile 确认无 LU/flow/PDE。几何 112 源/5 输入与配对 117 源/3 输入 SHA 未变。原生同保存时刻迁移脚本只准备，必须另排 LU 后测实际 core31/31、front20、全部原门及原生读回；此处不声称已原生迁移或完成后续 PDE。

### 后续实际原生迁移与从零对照

- candidate1 实际原生迁移 session9927 exit0：step9205 不推进时钟，核心35.1242/35.9637、front39.974、峰相对跳8.34666e-6，质量/迁移/历史/不可变参考/原生读回全部通过。118冻源和5输入SHA复核未变，实测peak footprint5.577GB。输出 `performance_lab_fresh_stage5_native_trial_v1`；尚未沿此新候选推进PDE，不构成原t0自动长时资格。
- 新cap210000 A/B分别1654/2046步到共同物理终点，B自动两次增长。实际无LU比较 session6662 exit0，两窗科学门均失败：rhoX峰差55.19%、连续x核心宽度差50.07%；C/D/E仍待完成。初始空间截断的A/H解释不能替代这些联合误差对照。
- 独立初始观测fallback候选H64/H32/H16原生四步均exit0，前两者各自动使用一次观察回退，H16原网格合格不使用。维护+ipm未合入。默认旧路径和H64 split的首轮严格审计失败保留，正在逐字段检查；不得写作兼容/恢复已通过。

### 2026-09-09：自动初始观察规则通过资格并合入

合入八个runtime文件（新增两个纯policy/planner）与schema覆盖测试，显式可选
`initialMeshObservationFallback` 支持原始k8大箱初选失败后的自动解析重测；省略时不补字段。
全部运行代码与已验证candidate_v2字节相同，合入前维护原文件逐字核对并备份。
H64/H32实际回退与H16 unused通过；all/equivalence session30711通过六种物理元组与三种迁移，
default v1/v2 session43690数学组逐位通过。H64实际2+2步恢复session15200完成，
只读完整补验53584验证全部原生数学状态/history/snapshots/cursor严格相同、只重建一次且不重跑初选。
首次合法NaN比较、prefix终点改动与结构字段排列导致的审计失败均原样保留；不放宽数值门或剔除历史字段。
维护目录postcopy49092无LU验证路径、schema、旧v2与新observer native读取通过。
证据集中在V/initial_observation_promotion_v1；新H64原始t0至tau12单次运行已启动，尚未完成长时目标。

C2从原始t0自动到共同物理终点，step1782、target42、一次增长；ABC实际观察已完成。
AC峰差.0212%但核心窗梯度场最大差约.8%–1.0%，两窗原科学门仍失败；D/E待运行。

## 2026-09-09：未来 factor3 家族与资源脚本准备

新增 performance_lab 三个 MATLAB 研究入口及独立资源 watchdog；原四 reference 成员逐字段 bitwise，仅形成未实现的 policy3/family2 草案，追加[3,2]/cap310000/growth3。复用原 tau8 既存候选、只重构一对轴并作冻结插值，全部原 field 门通过；没有重跑70×4、没有LU/PDE或旧CP升级。原生 normalizer 仍拒绝新草案。961×321 单factor/flow资源脚本仅dormant，默认false和root窗口门已预检；3文件CodeAnalyzer0、9纯watchdog测试与真实旧内存日志解析通过，123源/4原输入SHA不变。资源峰值和native3运行/恢复资格尚未获得。

- 2026-09-09：独立研究副本 `result/verification/autonomous_runtime_20260909/node_family_v3_directional_candidate_v1` 实现显式八成员 node-family v3（原4+[3,1]/[1,3]/[3,2]/[2,3]，cap310000、最多3增长、每族3候选），兼容既有可选t0解析观察；维护+ipm及旧CP不改，旧五成员草稿单独保留。实际noLU回归11869/57799均exit0：5旧原生CP/4旧result完整逐值、原4参考成员bitwise、8合成方向账本、41配置/持久化负例，补充4纯审计/telemetry正例和13审计/迁移前拒负例；20文件CodeAnalyzer0。两次补充启动因/tmp同名脚本cwd遮蔽而在helper前失败，日志保留，最后显式cwd+which修复。未执行v3 LU、原生迁移/CPsplit或PDE，尚无资源/长时资格；详见 research/performance_lab/DIRECTIONAL_NODE_FAMILY_V3.md。

- 2026-09-09：独立八成员v3候选修订2补充每级候选持久化证据，仅新增v3演化candidate/成功audit.attemptSummary的targetLevelId/localPairIndex及严格顺序/≤3/末目标核对；v1旧冻源保留，维护+ipm不改。mesh独立合成24→第3eligible成员反例保留并在新边界拒绝。无LU16633 exit0，旧全回归及4正13反审计重过，新3正13反顺序测试全过，21文件CodeAnalyzer0；125源/测试和10输入SHA一致。仍未执行v3 LU/PDE/原生迁移或restart资格；详见 result/verification/autonomous_runtime_20260909/node_family_v3_directional_candidate_v2/README.md。

- 2026-09-09：独立 `result/verification/autonomous_runtime_20260909/native_v3_H64_qualification_preparation_v1` 封存 H64 原始 t0 的 v3 八方向自然 factor3/原生 split dormant 协议（τ8、20000步、原数值门、无快照），维护 +ipm 与 v3 冻源未改。默认 false 实测 session73227 exit0，2正14反、3.5014735秒、5文件 CodeAnalyzer0；12项 Python watchdog 纯测试通过。119 runtime、输入 CP/bundles 与 SHA 全部 exact。原生 true、资源 LU/flow 和 factor3/PDE/replay 均未执行；资源门要求两方向实际证书及系统 swap 原始基线/活动增量，子进程0swaps不代表系统未换页。较早13反测试和执行helper原样保留。

- 2026-09-09：H8 原 v2 τ8→12 continuation 的真实 capacity failure 已独立 noLU 审查，原 CP/result/失败不变。`V/h8_tau11673_capacity_nolu_v1` session18959 exit0，29.0663s；接受11049场同N原X0/70，既定hierarchy330项有20X通过且首完整offline pair通过；X3原70有61X，Y3仍需hierarchy，三首pair全部offline门过。严格区分未保存11050 rho 与接受11049反事实；未迁移、未建LU、未推进PDE、未升级CP。CodeAnalyzer0，源及输入SHA exact；参数只读补充54932 exit0。

- 2026-09-09：独立 v4 持久化审查 `V/v4_independent_dispatch_review_v1` session36243 exit0，26.5882s，无LU/PDE。四个明确synthetic failure ledger的错误N/类型/epoch/step被首版完整validator接受，共有committed对照正确拒绝；原反例与源永久保留并交独立修订。实际H8接受11049配对rho/Dx搭配明确synthetic触发控制，检验4/8/7九候选真实dispatch与四项transfer前拒绝，不能称旧CP升级或native v4谱系。129源/测试、3输入SHA exact，helper CodeAnalyzer0。随后修订需区分真实accepted source与未接受next-step failure、历史保留及endpoint变更，当前不称已获得native资格。

- 2026-09-09：独立高阶原A+低阶五点Cholesky预条件GMRES小格筛查 `V/preconditioned_poisson_small_v1` session60099 exit0，49/65/129，数值7.49809s/进程18.3613s，child聚合RSS峰1.159GB<1.5GiB；无大LU/PDE/CP读写。27项通过原轴质量的问题，GMRES原方程/场门8过，加固定两缺陷修正27过，但17个修正内层flag3如实保留；65/H1e6/s16未准入压力轴修正后已知误差1.073e-8仍超门，尽管逐分量后向误差2.17e-16。预条件factor存储代理较小，tiny回代明显更慢，无生产/大场资格。123登记文件SHA exact、CodeAnalyzer0，原始失败与时序全保留；详见 research/performance_lab/LOW_ORDER_POISSON_PRECONDITIONER_SCREEN.md。

- 2026-09-09：后续Poisson有界测试仅准备于 `V/preconditioned_poisson_early_stop_preparation_v1`，未运行MATLAB/LU/PDE。冻结前轮119运行文件逐字相同；新登记129²方向拉伸六对轴、局部/远尾已知解及原MMS/Green，每阶段完整记录、原分量/行残差同时达门即停且最多两修正，与固定路径共同前缀须bitwise；旧s16失败仅只读引用。helper及两个runner CodeAnalyzer0，watchdog4纯测试及默认不启动检查通过，aggregateRSS1.5GiB/300秒并需root显式串行小窗口；新轴质量、修正精度和真实成本尚未实测。
