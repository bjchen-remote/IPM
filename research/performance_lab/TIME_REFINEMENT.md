# 同物理终点的时间步加密入口

`ipm_perflab_time_refinement(sourceFile,newOutputDirectory,settings)` 已用
真实 `ipm.solve` 和原生 checkpoint 完成 tiny 正/负验证。它比较源配置
`CFL=.5,maxDt=.004` 与独立实验分支 `.25,.002`，不修改共享求解器。

原生 resume 白名单不允许改 CFL/maxDt。本入口先严格恢复源 checkpoint，
仅改变 `state.config.time` 的这两个字段，再通过现有
`makeCheckpoint/writeCheckpoint/readCheckpoint` 写入独立 fine-seed。
除恢复 metadata 外，接受态数学字段、完整 log（包括 snapshots）和 cursor
都必须逐值相同。两个数值设置变化在独立 provenance 中明确记录，这不是
原方法逐位 restart 的声明。恢复产生的旧 LU 在两个 `ipm.solve` 前释放。

`settings.physicalFinalTime` 和 `settings.canonicalFinalTime` 必须显式提供，
在任何 PDE 推进前连同阈值、源签名和原配置保存为 `frozen_study.mat/json`。
两支共用同一物理终点、canonical 上界和步数预算。默认只允许 `65²` 以下
网格；真实计算须在母任务给出的串行窗口中显式设置 `allowLargeGrid=true`。
不自动改变线程，读取源文件前要求 10 线程。源网格必须关闭自动 remesh，
保持 source 的 remesh 配置；重网格后的原生 checkpoint 可作为新的 source。

```matlab
addpath('ipm_structured');
addpath('ipm_structured/research/performance_lab');
settings = struct('physicalFinalTime',targetPhysicalTime, ...
    'canonicalFinalTime',canonicalCeiling,'maxSegmentSteps',128, ...
    'window',[-2,2,2]);
report = ipm_perflab_time_refinement(sourceFile,newDirectory,settings);
```

报告检查每支全部历史的维护 `trustedMask`、源历史全部字段前缀的逐值一致、
数值轴/remeshCount 不变、终点 checkpoint 的严格回读及末态/历史一致。
两支必须真正到达预冻结的物理目标，时间差也必须达标；仅两支时间相同
并不足以通过。fine 还必须实际用了更多步，避免端点裁切使两个设置等效。
所有现有诊断和门保持原样，失败报告仍保存。

物理 lab 窗口的 rho/omega 比较复用 `ipm.output.compare` 的现有线性插值、
控制体权重和范数。rhoY 用维护四/六阶 `fdMatrix` 及 `Cy/Comega` 物理变换，
交给相同的比较核；不另写测试版权重/误差公式。`lab_deviations.mat` 保留
三场的 reference/candidate/difference 数组与物理坐标。另存原接口的峰对齐
比较作为辅助指标。P、曲率 C、二次峰位置/forcing、`c_l/c_omega`、
`C_l/C_omega`、rho-x/rho-y/全梯度峰均保存原值、差值和相对差。

默认研究阈值为：物理时间 `1e-10`；rho/omega/rhoY 的 lab L2 为
`5e-4/5e-3/5e-3`；梯度峰和 P 为 `5e-3`；规范率和曲率为 `1e-2`；
尺度为 `1e-3`。可以在开跑前通过 `settings.thresholds` 指定更严格的研究
阈值，运行中不根据结果调门。这些不是原生 solver 的新接受门，也不构成
无限时间或收敛阶证明。

`time_refinement_tiny_v2/tiny_time_refinement_report.mat` 的实际验证：

| 65×33 native source | coarse/fine 步数 | 两支物理时间差 | 目标时间误差 | 最大 lab L2 差 |
|---|---:|---:|---:|---:|
| physical，源已走 2 步 | 3/6 | 0 | -1.04e-17 | 1.66e-12 |
| transport_anchor + quadratic，源已走 2 步 | 4/7 | 0 | -1.94e-11 | 7.06e-13 |

动态终点的共同微小欠量来自维护 `minDt` 终点规则，未修正时钟或降 minDt；
在预冻结的 `1e-10` 内。两支前缀、网格、全历史可信、原生读写均通过。
负例故意给过短 canonical 上界：两支只走 1/1 步、时间差与场差均为零，
但距离目标仍 `0.0110009`，被正确判为失败。两个 MATLAB 文件 checkcode=0。

随后在母任务明确安排的串行窗口，用真实重网格后 step2226 checkpoint
完成了 `time_refinement_q512_step2226_v1`。目标 `t0+.001`、canonical
上界 `tau0+.03`；粗/细 6/11 步、墙钟 `48.46/82.40 s`（含原生 restore/
finalize/写入），最终物理时间分别为 `1.821587396389057` 和
`1.821587396389069`，差 `1.199e-14`。1174 行历史前缀逐值保留，完整可信、
固定网格/规范、分支与两支终点的严格原生 checkpoint 检查全部通过。

138798 个物理 lab 节点上 rho/omega/rhoY relative L2 为
`2.409e-14/1.658e-12/4.854e-13`，relative inf 为
`2.886e-14/2.868e-11/2.209e-11`。P、曲率 C、physical rho-x 峰的相对差为
`2.912e-15/-7.924e-14/3.660e-14`，`c_l/c_omega` 相对差为
`1.766e-12/7.141e-13`。这支持该短晚态窗口的时间误差很小，不能外推任意
长期误差或省略空间/远场检查。

初次 MATLAB 启动因 Qt neon 环境错误退出 137，尚未运行脚本；失败日志
保留。普通管道方式的一次重试成功，全部原 RCOND=`8.911277e-17` 告警
保留在 diary，数值计算正常退出 0。没有关闭告警或修改数值门。

全部证据位于 `result/verification/performance_lab_20260908/`，实际 q512
报告为 `time_refinement_q512_step2226_v1/time_refinement_report.mat/json`。
