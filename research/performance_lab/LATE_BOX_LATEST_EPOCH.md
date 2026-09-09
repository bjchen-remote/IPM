# 最新epoch的静态盒阶梯与独立fresh协议

本实验固定来源为adaptive_campaign_v4/stage_002的原生step3326 checkpoint，canonical tau7.093201945463523、physical epoch1.9545271366857122。仅使用root分配的额外LU窗口，内部各LU/solve串行，默认10线程严格验签。原checkpoint、main bridge及其网格/历史均不修改。

在独立输出 `result/verification/performance_lab_20260908/late_box_v4_stage002_epoch_v1/` 内冻结102份MATLAB源码：+ipm精确复制parent adaptive_source_v4（复制前逐文件确认与维护+ipm一致）；fresh协议及branch两个helper取mesh此前实际通过的冻结版本，static两个helper原样复制维护版本。源文件和parent checkpoint SHA256均保存于source_manifest.json，runner从该独立source目录运行并断言解析路径，不调用共享路径的演化实现。

静态调用mesh_screen_late_box(source,newOut,[1,.1,.01])，由它内部调用screen_frozen_box，避免重复LU。固定common inner nodes、rho和原runtime reference值，只裁掉outer nodes，Green边界使用裁后的source重算。原静态全门保持：full-box control、速度L2/Inf 1e-3、cL/cOmega相对1e-3/绝对1e-6、quadratic peak/位置1e-10、Poisson residual1e-8及原mesh质量门。相邻.1与.01都passed/eligible才进入动态。

动态严格使用原mesh_run_fresh_box_protocol(source,.01,.001,newOut, controls allowLarge=true, maximumSteps=64, staticReportFile=实际JSON)。五支顺序为native continuation、full natural、crop natural、full common epsilon、crop common epsilon。native继承原同网格case/history；四个fresh均在物理坐标与物理rho上生成新case、time0和新runtime references，absolute epoch单独记录。它们不是changed-box native continuation。

五支必须到同一absolute physical target1.9555271366857122，误差≤1e-11；完整result及native checkpoint历史trusted，terminal rho/axes/step/clocks/case逐值匹配，实际terminal core≥14、safety<.70。先通过full native/fresh covariance门，再允许crop分支；natural/common-epsilon两套crop/full误差门原样。任何依赖门失败即停止后续，保留全部原始文件，不修改helper/阈值。

初始restored全盒decomposition报RCOND3.303172e-18，原警告保留，未抑制；后续结论按预注册残差、trusted、covariance及空间/物理比较门给出，不从小全局残差推出所有局部行精度。

本次静态与五支动态比较现已全部通过。该结果只支持这个epoch与这个短区间的独立缩箱比较，不能回溯验证早期full-box历史、证明无穷物理盒收敛或直接授权生产缩箱。

## 实际静态结果

| 比例 | 节点x×y | 实际H_x / H_y（parent rescaled） | 速度相对L2 | Poisson residual | 全门 |
|---:|---:|---:|---:|---:|:---:|
| 1 | 1025×513 | 1e+06 / 1e+06 | 0 | 2.0616e-09 | 通过 |
| 0.1 | 959×449 | 94880.2 / 98127.1 | 1.39271e-05 | 2.18851e-09 | 通过 |
| 0.01 | 895×386 | 9669.01 / 9963.1 | 0.000152671 | 2.04747e-09 | 通过 |

full-box对照速度与rates逐值零差；.01核心速度Inf1.33126e-4，cL相对差9.39143e-5，cOmega相对差2.97673e-11。相邻.1/.01均满足原门，因此依原流程进入动态，无阈值改动。

## 五支matched-physical结果

所有分支各完成19步，stopReason=physical_final_time，absolute终点为1.955527136685712，误差均为0。每支完整result和native checkpoint历史trusted，终点数组/网格/step/clock/case逐值配对通过；四个fresh均新case/history。lambda=32.205386661267056，公共物理比较窗为[-.304284470,.304284470]×[0,.152142235]。

| 分支 | 核心x / y | safety | solve墙钟秒（共享机器） |
|:---|---:|---:|---:|
| native_continuation | 17.6443677 / 18.5503857 | 0.45340248 | 105.338 |
| full_natural | 17.6443677 / 18.5503857 | 0.45340248 | 101.469 |
| crop_natural | 17.6443677 / 18.5503787 | 0.45340248 | 58.101 |
| full_common_epsilon | 17.6443677 / 18.5503857 | 0.45340248 | 100.853 |
| crop_common_epsilon | 17.6443677 / 18.5503787 | 0.45340248 | 60.159 |

时间只记录本次含初始化/求解的共享机器观测，不作为独占资源的速度基准。full native/fresh covariance：坐标相对Inf1.87913e-14、rho1.78198e-14、rhoX7.99965e-12、rhoY6.80867e-10、cL6.56313e-12、cOmega9.13697e-14，全部通过原门，crop在该门通过后才启动。

| 指标 | 自然epsilon full/crop | 公共epsilon full/crop |
|:---|---:|---:|
| rho相对L2 | 6.896546e-07 | 6.896356e-07 |
| rho相对Inf | 1.777121e-06 | 1.777068e-06 |
| rhoX相对L2 | 2.387542e-06 | 2.387478e-06 |
| rhoX相对Inf | 2.728113e-05 | 2.728032e-05 |
| rhoY相对L2 | 2.476031e-06 | 2.475975e-06 |
| rhoY相对Inf | 3.427341e-05 | 3.42724e-05 |
| 物理梯度峰相对差 | 2.248909e-07 | 2.248843e-07 |
| 物理二次峰相对差 | 2.248562e-07 | 2.248496e-07 |
| 物理core宽x相对差 | 2.268004e-07 | 2.267938e-07 |
| 物理core宽y相对差 | 2.183557e-07 | 2.183492e-07 |
| parent单位cL相对差 | 9.420722e-05 | 9.420149e-05 |
| parent单位cOmega相对差 | 3.827476e-08 | 3.827242e-08 |

两套比较全部passed。仅改变full-box epsilon的对照给rho L2=7.79622e-11、rhoX Inf=1.21485e-9、rhoY Inf=1.53963e-9；所见缩箱差异不是这项epsilon控制造成。global LF线最大值和外边界closure仍可随crop改变，本实验没有把它们冻结为同一算子。

crop_natural生成新case `20260908T134443474_dynamic_isotropic_895x386_tp6fdc3a52_fed0_4876_aa77_e1c92ff19ee4`，终点checkpoint step19；其初态物理盒为±1471.0649871964852、ymax1515.8089040149603，物理anchor=.15214223499880195。原始initial_physical_data、result、native CP、完整lineage和branch audit均已交给mesh作独立noLU网格准备，未把它写回原native continuation。

## 完整性与补充记录审计

session58103正常exit0，额外LU窗口已释放，没有启动额外PDE。源checkpoint bytes及全部102份冻结MATLAB文件SHA256与执行前完全相同。四个原helper Code Analyzer为0，无公式或门限编辑。补充脚本 `ipm_perflab_audit_box_records.py` 只用已维护ctypes reader读取已完整通过native配对的checkpoint标量；不替代native签名验证，不建ops、不推进PDE。

该段保存记录的Poisson residual最大2.30664e-9，包含native继承完整历史时最大2.41569e-9。保存记录均有限，段内core最低17.64437/18.55038、safety最高.45340248。Poisson residual是补充观察，trustedMask本身不实施该残差门，保存历史也不覆盖每个未记录RK子阶段。

实际主报告为 `dynamic/protocol_report.mat/json`，`completed=true`、`shortIntervalComparisonPassed=true`、`eligibleForLongerIndependentBoxTest=true`；保留 `dynamicBoxConvergenceEstablished=false`、`productionDomainModified=false`。另有静态report、源码/CP哈希manifest、completion、supplemental_saved_record_audit及完整run.log。未宣称长期盒误差已验证或直接promotion。
