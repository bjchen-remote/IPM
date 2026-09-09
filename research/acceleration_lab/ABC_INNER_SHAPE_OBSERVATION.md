# ABC 同物理时刻：连续峰与双宽对齐的内形状

实际观察结果：在各自连续峰幅、峰位、90%横/纵宽对齐后，A/B墙面形状的相对L2 / 采样Inf差为 **1.9863% / 2.8353%**，纵向为 **0.27457% / 0.47916%**。差异大幅小于原峰幅和固定物理场的差，但仍显著可测；这不是原科学门通过，也不是极限Profile或强奇异性的证据。

完整输入与预登记：[registration.json](/Users/bjchen/Library/CloudStorage/OneDrive-北京大学/Research/IPM/ipm_structured/result/longtime/20260908_campaign_v1/acceleration_lab/ABC_aligned_inner_shape_9f16c481a07a41e3b0207a405494dfcf/registration.json)。输出：[report.json](/Users/bjchen/Library/CloudStorage/OneDrive-北京大学/Research/IPM/ipm_structured/result/longtime/20260908_campaign_v1/acceleration_lab/ABC_aligned_inner_shape_9f16c481a07a41e3b0207a405494dfcf/run/report.json)、[profiles.mat](/Users/bjchen/Library/CloudStorage/OneDrive-北京大学/Research/IPM/ipm_structured/result/longtime/20260908_campaign_v1/acceleration_lab/ABC_aligned_inner_shape_9f16c481a07a41e3b0207a405494dfcf/run/profiles.mat)、逐曲线CSV，以及 [PNG](/Users/bjchen/Library/CloudStorage/OneDrive-北京大学/Research/IPM/ipm_structured/result/longtime/20260908_campaign_v1/acceleration_lab/ABC_aligned_inner_shape_9f16c481a07a41e3b0207a405494dfcf/run/ABC_aligned_inner_shapes.png) / [PDF](/Users/bjchen/Library/CloudStorage/OneDrive-北京大学/Research/IPM/ipm_structured/result/longtime/20260908_campaign_v1/acceleration_lab/ABC_aligned_inner_shape_9f16c481a07a41e3b0207a405494dfcf/run/ABC_aligned_inner_shapes.pdf)。

## 真实端点与表示

分别重读ABC真实原生CP和配对result，使用原冻结runtime、默认10线程；完整原native端点/家族/expected config门及与旧comparator的audit逐值一致。相同物理时刻t=1.9140859724939365（实际三值相差≤3.4e-15），canonical time不同。连续geometry按旧comparator同一physical peak窗口[.1,.6]重算，结果与已保存geometry逐值相同，没有新挑峰窗口或修改规范。

| Case | canonical τ | physical peak x | physical C1 peak P | physical wall width | physical vertical width | native N |
|---|---:|---:|---:|---:|---:|---|
| A_H8_T32 | 4.0000000000 | 0.2986938121 | 4.0593595414 | 0.0121654717 | 0.0032340898 | 641×161 |
| B_H16_T32 | 4.6670145627 | 0.2326148204 | 6.2996479054 | 0.0060739689 | 0.0016532954 | 641×321 |
| C_H8_T42 | 3.9994186529 | 0.2988404635 | 4.0600072196 | 0.0121531618 | 0.0032268242 | 641×161 |

定义 `Uwall(xi)=H2[Omega](a+wx*xi,0)/P`，`Uvertical(eta)=H2[Omega](a,wy*eta)/P`。wx是围绕峰的connected .9P两个横向交点距离，wy是同一连续峰位置a上的首次纵向.9P下穿。P、a、wx、wy完全采用各case自己的连续Hermite几何。物理场与canonical场的幅值/长度转换在这个表示中相消；此变换只用于读数据，没有投影或改变任何原rho。

wall固定xi∈[-2,2]，vertical固定eta∈[0,3]。线性与Hermite两观察器用 **同一已登记几何和P**，只更换场的插值，因此线性曲线在零点不被强制成1。这里没有把linear观察器另找峰/交点后再择优。

## 全部三组差异

相对L2分母为第一个case的该完整曲线L2；下表Inf分母为同一曲线采样最大绝对值。反方向分母也全保存在JSON。

| second − first | 曲线 | Hermite L2 | Hermite Inf | linear L2 | linear Inf |
|---|---|---:|---:|---:|---:|
| B−A | wall | 1.986300% | 2.835325% | 1.986754% | 2.841656% |
| B−A | vertical | 0.274573% | 0.479164% | 0.274122% | 0.476623% |
| C−A | wall | 0.158065% | 0.251596% | 0.157914% | 0.252167% |
| C−A | vertical | 0.041571% | 0.089061% | 0.044333% | 0.094785% |
| C−B | wall | 1.871319% | 2.695893% | 1.871383% | 2.713217% |
| C−B | vertical | 0.242473% | 0.390103% | 0.240067% | 0.383496% |

L2在两条曲线原生cell端点变换后的并集上分段做Gauss4：固定几何下Hermite沿每段为三次多项式，差平方次数≤6，故理想算术下对**已表示的离散多项式**精确；不把这个代数事实当作原PDE解的精度界。Inf使用预登记801/601点及固定加倍到1601/1201点，全部pair最大变化为3.26867e-6；未宣称已求多项式的严格全局极值。

## 观察误差与原失败分别保留

| Case | wall linear−Hermite Inf | vertical linear−Hermite Inf |
|---|---:|---:|
| A | 0.022296% | 0.011523% |
| B | 0.023178% | 0.027558% |
| C | 0.019551% | 0.006474% |

最大自身插值差约0.02756%。A/B内形状差超过这个观察器差，也高于A/C的空间/网格策略敏感性（wall Inf0.2516%、vertical0.08906%）。这些比值提示不能把A/B剩余差完全当作本次线性插值噪声；A/C不是严格误差上界，尤其连续峰/根几何本身的空间误差尚未消失。原Hermite测试仅支持已登记绝对误差/协变，不能笼统声称所有phase四阶。

原“55%”应准确称physicalRhoXInf峰幅差：B相对A为55.1937%，quadratic峰约55.1881%。原固定物理窗中B→A的rhoX场Hermite Inf实际是153.775%，小窗[0,.8]×[0,.25]的L2为60.9606%，外窗[0,2]×[0,1]为37.0742%。它们与对齐的一维曲线范数定义不同，不能把两者相减解释为可严格归因的误差比例。旧AB和AC原linear/legacy及Hermite/continuous科学门均保持false；原报告SHA完全未变。

B在相同物理时刻已到τ4.6670，A/C约τ4.0000，且B峰位更靠近原点、physical两宽近乎减半、峰幅更高；这些几何变化造成固定物理位置比较与内坐标比较截然不同。剩余wall约2.8%、vertical约0.48%可以来自不同canonical进展、箱影响、动态网格事件与离散误差的组合。本次没有同箱多个时间的同观察序列，不能从ABC三条曲线证明剩余形状正在收敛或量化慢模态衰减率。D/E仍pending，不能补齐不存在的控制。

## 对下一步加速研究的含义

更有针对性的研究对象是这个几何对齐后的U及其真实链式Uτ，而不是把仍显著平移/收缩的原宏观场直接当作定常解。下一小候选首先需要同一条原PDE轨迹多个实际时刻的rho、真实FX、真实P′及几何导数缓存，分离形状残差和几何率；应与相同RHS预算的自然前进对照，并保留外圈/尾部velocity和共同physical-time原场硬门。

本轮没有真实RHS缓存，因此没有计算新Uτ、没有生成加速候选或把任何局部外推注入PDE；旧secant/三规范失败全部保持。两次绘图Python环境缺少matplotlib仅影响工具选择，最终用MATLAB直接画已保存曲线；一次初始启动路径笔误在研究函数前停止，记录在failed_launch.json，实际数据运行exit0。
