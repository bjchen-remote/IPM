# V4 step3326：冻结网格恢复研究

本轮没有找到满足全部原门的xy候选。最佳x已达标；保持纵向全局quadrature门后，已扫描y family的最好core为20.7850318，低于21。没有修改生产配置、原root/runtimeRefs、门限或checkpoint，也没有做LU/PDE。

## 实际来源

从v4 `campaign_manifest.jsonl` 的最后 `review_required` 和 `endpoint` 确定 `stage_002` 的实际step3326、tau7.093201945463523、physicalTime1.9545271366857122。默认10线程调用原生 `ipm.output.readCheckpoint` 严格校验schema4 integrity，再与配对result的rho、x/y及canonicalTime逐值比较，完整结果历史trusted。

最后真正接受的网格是 `stage_001/platform_design/moving_platform_02.mat`（fine112/.90）；同目录 `moving_platform.mat` 是被拒的fine96，不作来源。代码逐值验证current轴等于该已接受platform、platform由immutable root factory精确产生、其root reference仍与原root轴相同。当前entry x/y各为增量depth0；基于current的候选各增到1。原root未被current冒充，checkpoint只读完整性校验，没有restore或修改runtimeRefs。

## 横向非单调搜索

沿用已经验证的六档sigma和33点加局部细化，共246个唯一A，仅前三项做完整frozen-pair。focus固定实际core/front交集，anchor精确1。最佳sigma=.18、A=.301671318，coreX21.577345、front20.555727，ratioX1.07756168，峰跳1.77615e-5。

另外五档x没有达到21/20。同参数维护selector的sigma=.18原本就能达到core21/front21.37252，本轮未出现新的维护搜索漏解。与原root-y sigma=.25组合时，coreY20.446874，故三个full-pair候选都被分辨率门拒绝。固定x MAT仍保存，用作独立y恢复实验的输入。

## 纵向原门扫描

先从immutable root y扫描sigma=.18/.25/.35/.5/.75/1/1.5/2/3，目标21×1.15=24.15。全部不足后，再从真实entry y逐项独立构造单跳候选。没有一个candidate被用作下一个reference。

18条原构造均未达到最低21，因此每条在0至其原构造上界之间进一步取33个A（合计594个路径点），始终调用维护equalizer的同一轴构造。没有发现比原构造更高的可接收计数。所有轴保留ratio≤1.08、curvature≤.01、stencil rcond≥1e-9、全局quadrature比≥1e-8和局部quadrature/control-width比[.35,1.65]。

| sigma | immutable-root y最好core | current-y单跳最好core |
|---:|---:|---:|
| 0.18 | 20.13532472 | 20.43612104 |
| 0.25 | 20.44687448 | 20.44687450 |
| 0.35 | 20.62238375 | 20.45261219 |
| 0.5 | 20.71645611 | 20.45559932 |
| 0.75 | 20.76117901 | 20.45699605 |
| 1 | 20.77424768 | 20.45740075 |
| 1.5 | 20.78242114 | 20.45765308 |
| 2 | 20.78503183 | 20.45773358 |
| 3 | 10.08481222 | 20.45778852 |

最好原root来源为sigma2、A11.465507、core20.7850318；最好current来源为sigma3、A.815946879、core20.4577885。root sigma3在预注册A≤16上限只到10.0848122，明确记为`maximum_amplitude_reached`；它不是质量门饱和或全参数不可行证据。其他已完成路径的终点主要触及全局minimumQuadratureWeightRatio=1e-8。

在原y sigma=.25处，ratio仅1.05312、curvature.001595、stencil rcond9.54e-6；实际限制是最小weight1.94932e-5相对全箱平均的ratio=1e-8。保留H1e6与ny513时，这个门会限制可用的近壁间距。本轮仅检查两个来源及上述Gaussian family/有限幅度；不能由这些负结果声称所有更一般网格都不可能。

没有纵向轴达到21，故未再进行无意义的全场转移或事务（0个y full-pair）。既有x full-pair的质量/转移结果仍完整保留。真实事务floor20不能被拿来绕过本轮21/21冻结设计门。

## 产物

- `current_axis_nonmonotone_v4_stage002_v1/`：严格checkpoint/result配对、单跳正/负来源审计、六档x搜索、前三full-pair、维护selector对照、纵向限制详情和无LU/PDE profiler审计。
- `vertical_recovery_v4_stage002_v1/`：冻结protocol、18条原构造、18份33点幅度路径、所有最优轴、来源轴/quality、`vertical_report.mat/json`及无LU/PDE profiler审计。
- 研究入口 `ipm_perflab_search_vertical_recovery.m`。静态检查0，运行exit0。初版reference metadata的深度字段更名为明确的entryReferenceDepth=0和candidateIncrementalDepth=0/1，原metadata保留，并逐值核验source axes/quality未变。

以上均为恢复准备证据。主线的后续短桥接或缩箱restart及matched-physical验证，由独立生产/控制窗口安排，本研究未执行。
