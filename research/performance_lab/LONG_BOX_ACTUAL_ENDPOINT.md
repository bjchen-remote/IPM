# 最新epoch三盒在实际较晚原生终点的对照

本研究固定原生step3326、physical epoch 1.9545271366857122，以原主线bridge_segment_003的step3404为真正更晚参考：tau 7.2432019454635235、absolute physical time 1.9590400984028897。对应fresh IVP的总local physical elapsed为0.004512961717177522。目标和参考在运行前固定；不把未到目标的分支用于盒误差结论。

入口是 `ipm_perflab_extend_box_to_reference(shortReportFile,referenceResultFile,referenceCheckpointFile,newOut,controls)`。controls固定absolutePhysicalTarget、middleFraction=.1、maximumAdditionalSteps=128、allowLarge=true。已有完整短5支报告来自 `late_box_v4_stage002_epoch_v1/dynamic/protocol_report.json`，本轮独立输出为 `long_box_v4_stage002_epoch_v1/`。

仅三支natural public epsilon，按full_natural、crop_natural、middle_natural串行执行。full和.01 crop分别从先前fresh step19合法原生checkpoint延续，保留其case、时间和完整history；.1 middle由原生step3326真实物理场生成新IVP、新case与time0。三个物理盒仍保留各自原节点预算与原自然epsilon。先前common-epsilon短控制不延长，本轮不得声称长区间epsilon控制完成。

参考先只读验证result和严格native checkpoint，检查同source case/网格、source history精确前缀、共同物理终点、完整trusted、result/native完整history、rho、axes、时钟和终点尺度逐值配对。无reference restore/LU。每支只经唯一ipm.solve路径推进；full/crop检查全部原数值配置、捕获初始物理数组、runtime refs和fresh lineage保留。默认10线程，未修改任何+ipm公式或checkpoint签名。

较晚终点重新计算native/full covariance，原物理rho、rhoX、rhoY、axes与rates门不变。fresh rates除初态lambda=Cx0/Comega0恢复parent canonical单位。full covariance失败立即停止crop；full/small双向原linear比较失败立即停止middle。完整三盒做full/small、full/middle、middle/small，保留原field/梯度峰/二次峰/物理宽/rate门；actual core最低14、safety<.70、全历史trusted、physical_final_time与endpoint误差≤1e-11仍为原门。

`three_box_extension_tiny_v2` 用实际65×65局域场验证source→短5支→更晚native参考→三支延续，三个分支全部达到共同终点、full covariance与三盒比较通过，prefix/refs/lineage/数值配置通过。一步cap负例以maximum_steps退出但未达物理终点，立即写failure并且未创建crop/middle目录。tiny显式沿用维护fixture的resolution豁免，不用于生产core和safety结论。

首轮tiny_v1完整保留。其full支数学数据/历史/终点均正确，审计错误来自MATLAB反序列化匿名函数对象的身份不等；捕获physicalX/Y/rho及完整functions结构逐值相同。v2按完整functions记录比较匿名函数（并非仅func2str），同时加固result/native完整history与scale精确配对。Code Analyzer对两个research入口0警告。已作独立只读时钟、率单位和prefix审查，无额外PDE实现。

实际运行前冻结102份MATLAB源码及8个输入的SHA256，保留launch_recipe和完整日志。+ipm与先前最新epoch五支所用冻结代码相同。补充保存记录审计使用独立ctypes标量reader，无LU、无PDE，不替代native验证；保存Poisson residual不等价于未记录RK子阶段误差或componentwise精度。

实际三支现已全部通过，session63135正常exit0并释放LU。`dynamicBoxConvergenceEstablished=false` 与 `productionDomainModified=false` 保持，不作无限盒或早期原历史验证主张。

## 实际较晚终点结果

三支均以physical_final_time退出，full/crop的物理目标误差为4.440892098500626e-16，middle为1.1102230246251565e-15；早期进展消息中full误差写0在此明确更正。完整history、终点尺度与native payload精确配对，full/crop原history前缀、runtime refs、匿名初始数据捕获和数值规范全部通过。

| 分支 | 本段步数 / fresh总步数 | 末态core x / y | safety | 墙钟秒 |
|:---|---:|---:|---:|---:|
| full_natural | 62 / 81 | 16.3781211 / 17.9259703 | 0.488456519 | 288.053 |
| crop_natural | 62 / 81 | 16.3781222 / 17.9259402 | 0.488456486 | 164.46 |
| middle_natural | 79 / 79 | 16.3781212 / 17.9259676 | 0.488456516 | 274.326 |

full与crop同为追加62步，单次共享机器墙钟比1.7515；保存时间包含初始化/验证等开销，不作为独占资源benchmark。middle从epoch新起，共79步，不能直接以总墙钟和另外两支作速度比。

真正较晚native/full covariance为坐标Inf8.87520e-14、rhoInf2.45217e-13、rhoXInf3.86742e-11、rhoYInf9.95413e-10、parent单位cL相对2.87566e-8、cOmega相对2.87410e-8，原门全部通过。

| 指标 | full / .01 | full / .1 | .1 / .01 |
|:---|---:|---:|---:|
| rho相对L2 | 3.143885e-06 | 2.875872e-07 | 2.856583e-06 |
| rho相对Inf | 8.484599e-06 | 7.737905e-07 | 7.710809e-06 |
| rhoX相对L2 | 1.125827e-05 | 1.068323e-06 | 1.022626e-05 |
| rhoX相对Inf | 0.0001323013 | 1.206581e-05 | 0.0001202355 |
| rhoY相对L2 | 1.175303e-05 | 1.137061e-06 | 1.06743e-05 |
| rhoY相对Inf | 0.0001679016 | 1.531252e-05 | 0.0001525891 |
| 物理梯度峰相对差 | 1.053692e-06 | 9.609656e-08 | 9.575954e-07 |
| 物理二次峰相对差 | 1.053818e-06 | 9.610871e-08 | 9.577098e-07 |
| 物理core宽x相对差 | 1.089014e-06 | 9.931899e-08 | 9.896948e-07 |
| 物理core宽y相对差 | 9.227875e-07 | 8.414662e-08 | 8.386409e-07 |
| parent单位cL相对差 | 9.517315e-05 | 8.654077e-06 | 8.651982e-05 |
| parent单位cOmega相对差 | 1.927285e-07 | 1.111693e-08 | 2.038454e-07 |

所有原linear双向field/峰/宽/率门通过。虽然全域盒尺寸有约十倍阶梯，以上结论仍限定于固定late epoch和本次elapsed .00451296，不回溯证明早期边界误差，不改dynamicBoxConvergenceEstablished=false。common epsilon控制只在先前elapsed .001完成。

补充保存记录审计：三支全部保存标量有限，Poisson residual最大2.4917833e-9，段内core x最低16.378121、y最低17.848034，safety最高.488456519。它们不覆盖未记录RK子阶段；Poisson residual不是trustedMask门。102份运行源码和8输入SHA256全部未变。

额外原始数组lineage审计先tiny通过，再在实际三盒只读严格加载4份CP（parent+3终点），逐值核对parent case/clock/scales、box节点索引、固定fresh轴、physical anchor、捕获初始physicalX/Y/rho及report/native完整lineage，全部通过。无restore/LU/PDE。原始报告分别为dynamic/three_box_report.mat/json、original_data_lineage_audit.json、supplemental_saved_record_audit.json；tiny审计初次Code Analyzer的isscalar建议在启动数值前修复并保留失败日志。
