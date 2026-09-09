# 当前研究状态（2026-09-09）

目标ACTIVE，尚未获得经独立误差验证的强奇异性/近无穷时间Profile。硬要求：最终算法从原始物理t=0自行到长时，途中无需人工改网格、源码或参数。设计与实际失败的简明对应见[AUTONOMOUS_MESH_DESIGN_SUMMARY.md](AUTONOMOUS_MESH_DESIGN_SUMMARY.md)。
D=result/longtime/20260908_campaign_v1；V=result/verification/autonomous_runtime_20260909。旧状态按时间归档在D。

## 当前补充进展（优先于下方旧阶段描述）

- 裁剪研究v6第二段完整通过：step10439/equivalentτ12.843200751859165/绝对物理t2.025122944423388，physicalRhoXInf187.04933334、core25.14578/32.62213。下一段预测X21.13781触发重网格。00:23:51UTC root在第二段完整保存/验收后SIGTERM本人研究进程以释放大LU槽位，session29101 exit143；配对CP/result SHA保留于D/fresh_adaptive_launch_v6_candidate/resource_handoff_after_stage2_v1/result.json。这是外部研究资源调度，不是数值失败，也不是从零自动性完成。
- H64原始t0一次solve的session85571继续，最新已读manifest为step6866/τ8.600266968099579/物理t2.0210839391054036、641×321；原τ12目标及全部冻结规则不变。
- 组合v4修订2已完成第四阶/六阶/新旧一致性三项小套件（session5335 exit0，114.26秒），actualSuitesPassed=true，baseline/all仍未完成。独立正确性窗口峰组RSS1.919GB、footprint2.174GB，系统swap使用量未超过启动基线、无swapout，但发生4次系统swapin；无法归因于具体进程，原严格资源门不通过，不作性能证据。此前两次严格启动中断原样保留。
- 完整all/quick+equivalence已在独立growth_free_compute_windows_v1窗口完成session62138 exit0，106.67秒，全部数值套件通过，峰组RSS3.100GB/time-l峰RSS2.931GB、footprint2.305GB，源码与输入SHA未变。系统swap未超过启动基线、swapout零增量、swapin增8；新计算窗口通过，原strict资格false，不作性能比较。H8原始t0 τ16待当前一次late9205 baseline-only短窗释放；H8授权文件尚为false。
- H64冻结齐次速度模型的真正未来检验完成：预先选定首份τ≥8原生CP，实际step5874/τ8.000103920706282/物理t2.0089828325708456。固定墙扇区r[4,6]/[6,8]相对f/r模型误差0.778342%/0.872358%，原leading误差26.87543%/19.60957%，同一预登记2%观察门通过，源/输入SHA exact。实际session57519 exit0、无LU/PDE。仅有限扇区源场未来一致性，不是全外域边界、PDE精度或奇性资格。产物V/H64_frozen_homogeneous_future_tau8_v2；v1静态发现字段错误后未执行，原文件保留。
- 9842几何补全扫描session67675 exit0：正常分层搜索在330X返回top3，完整500X中56个通过原门，原4Y仅1个通过。当前Y全局求积比1.2699e-8接近原1e-8经验门，局部求积和stencil仍正常；不据此修改门。产物D/global_quadrature_9842_audit_v2，第一版仅CodeAnalyzer摘要告警停止，未读CP。

## 既有阶段记录与下一批

- H8 长时session24877已exit1：原τ8 CP续算至最后接受step11049/τ11.673078880623883/物理t2.30798220707281，641×321；下一未保存11050触发X forecast后原70项全拒，autonomous_mesh_axis_capacity。预定τ12未到。末rhoXInf76.86049347、二维physicalGradInf103.51727485，35remesh、累计峰跳.002413479978、完整可信链仍保留。冻结11049同N反事实分层330项找到20个X，首对完整离线门通过；未保存11050场不重构，尚无该状态原生迁移资格。
- 新A2/B2/C2均已exit0到共同物理终点1.9140859724939365，step1654/2046/1782；B2自动两次增长，末641×321，C2目标42且一次增长到641×161。ABC_actual_comparator实际配对通过，AB、AC两窗科学门均失败；D/E待排队。
- 原裁箱长时session51945已exit1：stage4 step9205/equivalent τ12.443200751859164端点通过，stage5原X族70全拒。新40族candidate1真实原生迁移已过，core35.1242/35.9637、front39.974、峰跳8.34666e-6。独立v6 driver/chain的27链、24反例及7trace边界均过noLU，source锁定；正在运行session29101执行D/fresh_adaptive_launch_v6_candidate/dormant_run.m，48个预登记阶段；从原9205按登记自行选择，不接独立试验CP。启动前210源/157输入SHA复核不变。首段自动迁移并实际PDE到step9842/equivalentτ12.643200751859165/绝对物理t2.0243654513057767，physicalRhoXInf175.15501916、core30.67025/34.41205，stageAccepted=true，1694.32秒；现stage2目标12.8432。这条研究续算不算原始t0自动性资格。
- 当前第二大LU session85571：V/from_zero_H64_tau12_auto_initial_v1/run.m，从原始t0一次ipm.solve，H64/Y32、uniform321×161、cap210000、tau12/maxSteps60000、全部规则事先冻结。真实自动初选后已自然增X与Y，step5586/tau7.800078720/物理t2.004261585，641×321；尚在运行，不保证有限族一定达到终点。
- 原cap110000 B session28196已exit1：保存step1520/tau3.819993/tphys1.801149，未达目标。下一步1521预测触发时Y候选全因core不足拒绝，回退到可信1520；失败与CP保留。旧cap110k C/D/E未运行，现由新五组cap210k队列替代，不能标失败/通过。
- 第二自然增长原生split session88314已exit0：真实2370→2481跨step2477/level2→4，全部数学状态、完整log/cursor/controller逐位，串行restore后psi/omega/velocity/physical及scales也逐位。没有虚构空快照的跨尺寸资格。
- τ8完整轨迹图session18761已exit0：V/tau8_complete_figure_v1/rendered/native_progress.png，真实801行；源105文件及输入SHA复核不变。
- H32/H64/H128原生四步全部通过。H32旧research auditor错误访问空transactions；修复后只读既有结果通过V/preselected_native_review_and_remaining_v2/H32_actual_review，不重跑PDE。H64/H128报告仍在原preselected_initial_native_preparation_v1/plan/runs。
- 新五组cap210k预检session19405 exit0：原候选轴和全family几何逐值相同，仅第四级资源mask变为准入，全部从原t0启动。119冻源/104原runtimeSHA复核通过；原门、C和CFL不变。
- 比较器最新冻源在expanded_zero_controls_cap210000_v1/source：真实新A/B/C native_qualified，D/E pending，旧B严格rejected。AC峰差.0212%，但核心窗rhoX/rhoY场Linf差约.82%/1.05%，未过原科学门；不单看峰值认定空间收敛。
- 初始观察策略已合入维护+ipm八文件（其中二个新增）与schema测试；原生H64/H32回退和H16 unused、all/equivalence30711、default43690、实际split15200后只读完整补验53584均通过。postcopy49092 exit0，实际维护路径/config/新旧native读回无LU通过；所有运行冻源未改。

机器16GiB/10核，只允许主线+一个附加大LU。不要修改任何运行中冻结源，不要伪造缺失结果或放宽数值门。

## 长时主线实际进展

v4 stage6因512步调度上限在step7510/equivalent tau11.83016857安全停止，原stageAccepted=false永久保留。
v5只改研究调度/原生链准入，104个+ipm文件与v4逐字相同；v3登记预算4096，显式allowlist7510，首段补原targetLocal=.14749081749446574。23-link原链/17调度反例/3正例通过。
v5首段35步到step7545/tau11.843200751859165，通过原门；随后24-link完整实际链重审也通过（D/v5_first_stage_chain_review_v1）。
最新严格终点step9205/equivalent tau12.443200751859164、绝对物理t2.023556209610179，max abs(rho_x)=163.90802554162542、core21.92523212/32.99159751、安全.3648764107，累计jump.0004972130365。
stage4通过原门，下一原X候选族耗尽。40项新有限参数族仅新算38项，10项全几何门通过，原warp不变；按预登记条件不执行后备9项。前三完整冻结场pair通过，第一项真实transfer和原生读回已过，但未沿新网格推进PDE。
进一步的分层参数搜索在9205实际L1即找到pair（330X/260新），19.284秒；8057原primary报告/三候选逐值不变且无fallback。它按constructor预算细化参数，未使用step/时间特判；前三完整冻结场pair也已过，首项实际原生迁移session98240已exit0：schedule77/core33.97135/36.47823/front29.3687/峰跳8.6084e-5，累计.000583297052，签名CP roundtrip通过，无PDE。实际峰内存6.248GB，122源/4输入运行后SHA不变。纯模块候选V/hierarchical_planner_candidate_v1_41e91d9ebf20仅新增两函数、119旧文件exact；session83203已通过原报告/选择、500预算耗尽、坏输入、单位协变和不同N真实reference等noLU检查，未接入policy，不能把几何成功计作实际演化。
最新图D/continuous_profile_tau1224_figure_v1/presentation_v1/profile_evidence.png来自六个真实CP；当前物理0.9核心宽度3.461686896e-5/8.014657820e-6，图源及输入SHA复核通过。

主线是已登记的独立裁箱物理IVP，equivalent tau只是比较时钟，不声称继承父问题原生history，不重置epoch或再裁箱。
C保持exact_gauge_no_feedback_v1/[1,0,0,0]/transport_anchor/wall_omega_quadratic_peak，原CFL/dt/数值门不变。RCOND警告保留，不以小残差代替前向精度。
旧fresh_profile_campaign_review.physicalGradientMaximum实际指physicalRhoXInf；从零新报告同名字段可能指physicalGradInf，比较时必须直接取明确history原字段。

## 原始t=0自动网格：tau8已完整达到

V/from_zero_tau8_v2_cap210000_v1 session74751已exit0，原始t=0一次冻结ipm.solve到tau8，无中途改动/人工触发。
step4545、物理t2.2412067821954538、final_time；25次真实remesh、2次自然增长，末641×321/level4，core27.38530409/29.72754570，累计峰跳.00179034457402。
二维全梯度physicalGradInf=26.502507910077174；全native history trusted，原始t0/配对/预定终点检查通过。
wall3763.106秒、maxRSS3.653304320GB、peakfootprint3.821654360GB、0swaps。
第一次增长step1135/tau2.6913754578032374（321×161→641×161）；第二次step2477/tau5.594870008425091（641×161→641×321）。有限族与有限箱仍不是无限时或域/空间/时间误差资格。

原v1在1134/tau2.688619固定N容量安全停止，失败保留。此前独立v2 cap110000连续到tau4：step1653、17remesh、1自然增长；1102→1177完整原生split逐位通过。
其20份CP严格可读、13份共有v1前缀逐位。原快照测试曾以history节点数代替实际帧尺寸，已修复：真实1177 CP281帧/1653 CP及result401帧，各帧rho/x/y配对且确含321×161与641×161；补验noLU通过，旧报告未改。
第二增长2370→2481预检及实际完整原生split均通过（V/second_growth_native_restart_v1），原快照为空并单独验证。最大附加LU同时1；真实runtime93.53秒/maxRSS3.4508GB。

全部v2配置/纯几何/实际迁移/持久化短测试与all/equivalence已通过。政策default仍v1，v2显式启用，旧CP不升级；result→CP对v2明确不支持，只用native CP。

## 新初始网格规则研究

V/large_initial_box_bounds_v1：相同321×161，H32初选解析core失败可通过最多3轮中的第2轮修复；H64/H128/H256粗观测原X族皆0，无候选可迭代。
V/large_initial_analytic_probe_v1：固定解析观测格（内spacing.025、半宽4、尾比≤1.05/20格平滑坡、总观测≤600000），仍设计原321×161、原门/原70×4搜索。H32/H64/H128实际候选全部core/front/family几何通过，H256仍0。
观测格是解析采样而非求解格或演化场插值；其内区FD导数误差3.116e-7，原H64/H128粗观测误差9.28%/20.56%。三项预选轴与真正自动fallback的H64/H32原生四步均过，H16无需回退；维护模块现已接入显式可选规则，兼容/恢复/完整套件均完成，H64从零长时在运行。
首次default审计误用isequal比较合法NaN；首次split人为改变prefix finalTime，造成3个终点dt诊断差异；第二次split的strict helper将rescaling字段顺序当差异。原失败均保留。完整补验不豁免history/cursor字段，严格核对字段集合/类型/形状/值/NaN位置，实际恢复只build一次且不重跑初选。
原H32全部数值字段逐位重审通过，仅异常stack的冻结路径/调用者行号差异单列；失败原样保留。详见INITIAL_ANALYTIC_OBSERVATION_PROBE.md及故障图谱。

## 下一版自动规则的独立候选

八成员方向族v3的独立v1源已过noLU：保留原四成员，追加[3,1]/[1,3]/[3,2]/[2,3]，cap310000、最多3次增长。旧默认/v1/v2、5真实CP/4真实result与旧四族保持，41反例及附加边界检查通过；尚无资源、native迁移/恢复或PDE资格。独立审查实际复现持久化attemptSummary仅总上界、缺逐级来源的证据缺口（并非runtime实际超额），故首次v3 PDE前另建v2修订，逐attempt保存targetLevelId/localPairIndex并严格核对次序/每级预算；旧失败与冻源保留。

v3 修订v2已有逐级次序/预算核验；真实H64冻结场961×321与641×481两项operator-only已顺序完成（session97156/71839），1build/1flow、有限流场、峰内存3.376/3.361GB，各窗口系统swap使用量无增长且swapins/swapouts零变化。系统已有约2.43GB swap，不能写成系统无swap。尚无v3自然factor3迁移/PDE资格。

组合v4已锁定V/hierarchical_directional_v4_candidate_v1/source（122文件）；旧v1/v2/v3 5CP/4result/plan逐值，8正36反、实际8057/9205、单位协变/500耗尽/differentN及H64原初选无LU通过。20文件CodeAnalyzer为0。仍待组合全数值回归、原始t0一次长跑、自然后备/增长及严格split；不能把独立模块资格拼作组合PDE资格。正在准备H8从零τ16/max60000与H64从零τ8两份独立协议，源仍保留旧roundedAxis。

H64两早期增长都是Y forecast触发，前步到迁移前core计数下降4.329%/3.213%；只用OLS的同状态反事实预测26.185/29.325并不触22，而max割线触发。但H8早期更多Y forecast迁移仍不增N；当前transaction未存axisReport，不能凭末态断定较小族何种质量门失败。该缺口与初始远场节点分配混合，详见research/acceleration_lab/EARLY_H64_H8_GROWTH_AUDIT.md。

roundedAxis浮点停滞优化已作为唯一runtime文件合入：1200构造axes/info逐值、8057/9205完整纯规划与全部hierarchy反例通过，构造器6轮ABBA速度比1.3107–1.4510、中位1.4149；all quick/6旧算例/3迁移12177已exit0，维护112个.m与已测源exact。只认构造器收益，无整段PDE提速声明。所有运行与v3/v4冻源未动。证据V/rounded_growth_stagnation_candidate_v1。

## 尚未突破的加速与误差问题

晚期32/42/56网格残差场差非单调，不能靠标量形状稳定认定极限。单侧空间有效指数随尺度/时间漂移，时间宽度指数不等于Holder指数。
新Poisson memory4子空间shadow在49/65各80步保持全部原数值组逐位；78尝试77通过1违反原残差门，原拒绝保留。完整oracle审计比原flow慢，未省回代、不改C/轨迹、没有生产加速，已止损。
RCOND已被decomposition缓存，不存在每RHS重估热点；警告输出成本未单独量化，未关闭警告。先前JFNK/Lawson/LF/Sylvester/行缩放等结论保持。
固定容量趋势桶窗仍仅研究通过，非原polyfit逐位等价，未接入生产。

ABC的实际共同物理时刻内形状观察16855已exit0：各自连续峰/90%宽对齐后，Hermite B−A wall L2/Inf为1.9863%/2.8353%，vertical .27457%/.47916%；C−A分别.15806%/.25160%和.04157%/.08906%。原linear结果、全部原固定窗failed gates保留；对齐不能补足域收敛，也不能忽略A/B canonical时间不同。图D/acceleration_lab/ABC_aligned_inner_shape_9f16c481a07a41e3b0207a405494dfcf/run/ABC_aligned_inner_shapes.png已目视检查。

## 原始初选的域/离散联合敏感性

A原CP1653与B已完成CP501的真实初始行皆step=t=τ=0、Cx=Cy=Cω=1，基础物理/规范/输运/椭圆配置exact。初始c_l .4349185430864652→.46173112115203213，相差6.16497%；rhoX峰相对差7.73e-7。双方选轴不同，不能将该率差单独归因边界，也不替代终点控制。证据见research/acceleration_lab/A_B_ZERO_TIME_OBSERVATION.md。

原始非紧支撑源的t0精确Green外壳积分显示H收敛主项A/H，A≈.4247699062；这是空间截断收敛，不是时间慢模态的证据。真实tau8物理半箱1.729180878，预登记r=4–6及6–8远场窗口均无覆盖；不能把A/H当此末态实际遗漏率。完整边界trace模板的纯积分MMS在研究，尚未注入PDE或修改C。详见research/acceleration_lab/PHYSICAL_FAR_COVERAGE_ACTUAL.md。

## 新增实际远场与时间外推结果

原始t0完整外域Dirichlet trace固定算子实验已实际完成H8/H32（51638/68377）。每例1build、2Poisson、0PDE；原rho/omega/网格/参考不变，原C公式重算。c_l分别 .4349185430864652→.48821245677200925、.4750938595229709→.4883710476626841；两箱差 .0401753164→.0001585909。116注册源/输入SHA保持；这只是t0外域效应实测，不能充当演化远场闭合，未改生产C。

晚期固定rank1内形状外推已做真正未来保留检验：训练equivalentτ9.243–10.643，冻结后读10.843–12.243。终点wall L2 rank1 .6671% vs linear .3115%，vertical .1576% vs .08056%；原科学门拒绝rank1，不认定突破线性慢收敛。几何量预测与已知未来宽度对齐的形状误差分开保存，详见research/acceleration_lab/LATE_INNER_TEMPORAL_HOLDOUT.md。

H8/H64最新只读轨迹图：V/H8_H64_progress_observation_v2/native_progress.png（图中H64固定到τ6.8004、非当前实时终点），源/输入实际计数见final_audit。H8 τ11.673全梯度103.51727，物理半箱1.59350；H64图端全梯度22.70992、半箱10.50867。两终点时钟不同，只是进度图，不能作域收敛比较。

组合v4修订2已锁定：V/hierarchical_directional_v4_candidate_v2/source，仅validator修failed source边界；原旧兼容/8正36反与6正13反通过，独立审查未发现剩余确定问题。保留合法旧未接受attempt跨缩短终界/同step更短接受及后续增长的语义；未保存场仍不可重演。原v1漏洞证据V/v4_independent_dispatch_review_v1保留。small_suites(fourth/sixth/equiv)准备中；baseline含真实513×257 build/solve，等待大slot，不能声称全套已过。

低阶预条件Poisson tiny60099：27合法问题加两次原残差缺陷修正后场/方程门通过，但更慢，且17个内层停滞flag保留。未过原轴门的s16压力例依然已知解1.073e-8失败。没有生产替换或大场资格；下一候选只准备合法强方向轴与原残差达门即止。详见research/performance_lab/LOW_ORDER_POISSON_PRECONDITIONER_SCREEN.md。

H64真实物理远场[4,6]/[6,8]两开发CP均有原支持，但Ω相对初始1/r主项差18–27%；条件次阶tB/r²缩小约2.3–3.5倍。进一步无拟合固定齐次速度的二维特征模型差为.40–.81%，ODE阶梯与Jacobian独立差分通过。源码与τ8未来CP规则事先锁定，prospective检验仍未执行。它不是自洽IPM、没有改C/外边界，没有证明均匀全外场传播；详见research/acceleration_lab/H64_EVOLVED_FAR_SOURCE_OBSERVATIONS.md及V/H64_frozen_homogeneous_characteristics_v1。

late9205真正JFNK研究baseline-only准备完成：D/acceleration_lab/late_jfnk_baseline_preparation_v2_8c7f1e019a8147f29dcfc3070e850a25。未来仅一次实际saved-grid build+一次flow、原record逐值、释放factor后保存无LU cache；当前0LU/PDE。齐次幅度伪解与DP[Jv]+D²P[v,F]=0已明确，不能把裸缩幅降残差当shape加速；bordered square Newton/GMRES与矩形observer最小二乘分开，候选未实现/执行。
