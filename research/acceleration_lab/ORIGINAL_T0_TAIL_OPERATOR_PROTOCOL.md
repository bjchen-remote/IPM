# 原始 t0、固定算子的完整外尾边界对照

状态：两套实际初态的无 LU 准备已通过。真实 Poisson / velocity / 原规范对照尚未运行；执行脚本默认 `executeNative=false`。不改生产代码、初值、时钟、C 公式或任何旧实验。

## 固定输入与所问问题

使用已完成 A/H8 的最早 step116 checkpoint，以及 H32 预选初态的 step1 checkpoint。两者保存了真实 t0 snapshot；严格签名、全部 trusted history、remeshCount=0、当前 axes=base axes=t0 axes 均已核对。这里观测的是首帧 t=0，不是 step116 或 step1 的演化场。两者各 321×161（51,681 个节点），分别比较自己同一网格上的两次边界选择；A 与 H32 的网格本来不同，不能把它们的差直接解释为纯箱边界误差。

[准备登记](/Users/bjchen/Library/CloudStorage/OneDrive-北京大学/Research/IPM/ipm_structured/result/longtime/20260908_campaign_v1/acceleration_lab/original_t0_tail_operator_preparation_cf06cc9724df40bf8eec8ffabda3f2dd/registration.json) 固定两份 CP、冻结 runtime、节点上限 60,000、求积阶数、窗口和预算。[准备结果](/Users/bjchen/Library/CloudStorage/OneDrive-北京大学/Research/IPM/ipm_structured/result/longtime/20260908_campaign_v1/acceleration_lab/original_t0_tail_operator_preparation_cf06cc9724df40bf8eec8ffabda3f2dd/run/report.json) 保存严格来源与全部求积门。

对每一网格，未来比较：

1. 原 `ipm.evolve.rhs` 的完整结果，作为实际数值路径 baseline。
2. 相同 rho、成对 `Omega=rho*Dx'`、算子、原 Crefs，在原 `greenBoundary` 上仅加原始解析 k8 初值在盒外的完整 Dirichlet 迹；用原 `poisson(..., boundaryOverride)` 解出 psi / u 后，调用原 `isotropicGauge` 和完整 WENO 组装。

这不是直接给 c_l 加 A/H。盒内原 FD 源、Green 求积/截断/压缩及椭圆离散全部保留，所以即使加入准确的盒外迹，也不产生全空间速度的数值真值认证。

## 完整外尾，非内点展开

原始数据为

`rho0(x,y)=log(1+x^8/(1+y^8))/8`,
`Omega0(x,y)=x^7/(1+x^8+y^8)`。

正象限外部区域为 `([0,∞)×[0,∞)) \ ([0,H]×[0,H/2])`，用原 double-odd Green 镜像核积分；负 x 边界由奇对称得到，底墙/对称轴严格为零。964 个实际边界查询包括左右边、顶边、底边及重复角点。全迹使用 Green 核，不用 `xy*A/H` 代替人工边界上的函数。

为了稳定独立核验，分别计算：

- 齐次 leading 源 `x^7/(x^8+y^8)`，径向解析原函数 + 角向 Gauss256。
- 精确余项 `-x^7/((x^8+y^8)*(1+x^8+y^8))`，以观察点为极点的完整二维积分，Gauss128/256。
- 独立用原 `Omega0` 在观察点极坐标直接积分，Gauss512。

前两项之和与第三项比较；避免把仅 leading tail 的 MMS 冒充 exact initial datum。原全迹 MMS 的首轮 256 阶变阶门失败、随后固定 512 阶追加通过均原样保留，见 [FULL_GREEN_TAIL_TRACE_MMS.md](/Users/bjchen/Library/CloudStorage/OneDrive-北京大学/Research/IPM/ipm_structured/research/acceleration_lab/FULL_GREEN_TAIL_TRACE_MMS.md)，本次没有重跑或改写它们。

| 无 LU 实测 | A / H8 | H32 |
|---|---:|---:|
| 两方法 trace 差 / max(H,abs(trace)) | 6.84314e-14 | 6.84314e-14 |
| 余项128→256差 / max(H,abs(trace)) | 2.42725e-13 | 1.56882e-18 |
| 近轴 `H*psi/(qx*qy)` 的 scaled 差 | 5.54334e-13 | 5.54445e-13 |
| 精确初值余项最大绝对迹 | 3.31344e-7 | 2.02254e-11 |
| 完整外尾最大绝对迹 | 1.80336688 | 7.21346814 |
| 964点全部求积时间（含独立高阶校验） | 16.889 s | 13.912 s |

固定门分别为 trace scaled 1e-10、近轴商 scaled 1e-9；轴上零值精确。初态 rho 与原解析采样逐值相同，纯 1D 维护 Dx/Dy 重建的 Omega 最大值与 native history 首行逐值相同，mass/range 及 axes/base 同样精确。profiler 检查无 build/restore/velocity/poisson/initializeScaling/advance/solve/transfer 调用。

这些积分差是离散求积一致性证据，不是严格区间误差界；还没有经椭圆逆算子传播成速度/规范率误差界。

## 恢复参考量而不重新初始化

`initializeScaling` 会计算原速度并冻结 strainTarget 等参考量。未来试验禁止调用它，尤其不能在 tail 分支重新初始化并改变问题。

一次 `mesh.build(config, saved custom axes)` 后，调用维护的 `restoreRuntimeReferences` 恢复源 CP 的原冻结 references，要求完整 rescaling 逐值相同。再使用纯 `initialScale` 的 t0 identity scale，绝不借用 CP 当前非零时刻的 scale。

原 baseline 必须先满足：

- Dx/Dy、quadrature、axes、原 refs 与准备 context 逐值相同。
- 原 `rhs` 算出的 c_l/c_omega/c_r 与实际首行 history 逐值相同。
- 从 baseline 实际 physical flow 用原 `isotropicGauge` 重组的完整 RHS 与原 rhs 逐值相同。

任一失败即保存失败并停止，不转成容差通过。tail branch 保持 rho/Omega/refs/峰值不变；由完整 RHS 的成对 Dx 计算真实 quadratic P′，另列原代数 residual，不能强置 P′=0。

## 无 LU 入口回归与独立审查

[回归结果](/Users/bjchen/Library/CloudStorage/OneDrive-北京大学/Research/IPM/ipm_structured/result/longtime/20260908_campaign_v1/acceleration_lab/original_t0_tail_operator_harness_a34bf878d09c40efae8d895959d3033c/run/report.json) 的两正例（省略第三参数、显式 false）与三项 test-only 损坏 context（Omega、mass0、初始 c_l）均按预期通过/拒绝。profiler 确认无 LU/PDE 路径调用。准备与回归的原 CP、context、冻结 runtime 共 114 项输入 SHA 均未改变；原准备的两个 helper SHA 亦未改变。失败反例不是修改过的原生 checkpoint。

mesh_audit 独立只读复核未发现这两个已注册 case 的算子、边界、速度或原规范路径错误；目前尚未实跑 true 分支。入口按固定 SHA 的这两份 context 使用，不能把它当作支持任意初值/任意动态几何的通用外尾求解接口。

## 固定测量与执行资源

核心窗口固定 `[0,2]×[0,1]`，holdout 固定 `[2,4]×[0,2]`，另列完整盒。输出原/新 u1、u2、维护 Dx/Dy 的四项速度梯度、c_l/c_omega/c_r、完整 RHS 差、真实 P′、全边界 native/tail/combined 数组和代数 Poisson residual。t0 的 Cx=Comega=1，因此 canonical 与 physical rates 此时一致。

每个 case 单独 MATLAB 进程，默认10线程；一次 build、同 factor 两次 solve（baseline 与 tail）、不推进时间。保存数组后清除 ops 释放 LU。记录 build/solve/gauge 时间、build前后/释放后 RSS。不能以这两个回代的计时宣称时间演化加速；三类误差——外尾积分、盒内边界/FD源、椭圆离散——仍需分开。

[冻结 harness 目录](/Users/bjchen/Library/CloudStorage/OneDrive-北京大学/Research/IPM/ipm_structured/result/longtime/20260908_campaign_v1/acceleration_lab/original_t0_tail_operator_harness_a34bf878d09c40efae8d895959d3033c) 含 `source/`、runtime/原 CP/context SHA、无 LU 回归及两份 dormant 脚本。当前维护目录的新改动不进入这套锁定 runtime。

```matlab
% 仅在 root 明确给出 LU 窗口后，单 case 调用；普通准备默认 false。
report = ipm_accellab_run_t0_tail_operator_pair(contextFile, uniqueOut, true);
```

[H8 context](/Users/bjchen/Library/CloudStorage/OneDrive-北京大学/Research/IPM/ipm_structured/result/longtime/20260908_campaign_v1/acceleration_lab/original_t0_tail_operator_preparation_cf06cc9724df40bf8eec8ffabda3f2dd/run/A_H8_initial_context.mat)；[H32 context](/Users/bjchen/Library/CloudStorage/OneDrive-北京大学/Research/IPM/ipm_structured/result/longtime/20260908_campaign_v1/acceleration_lab/original_t0_tail_operator_preparation_cf06cc9724df40bf8eec8ffabda3f2dd/run/H32_initial_context.mat)。它们是明确标记 `isNativeCheckpoint=false` 的研究数组上下文，不是伪造原生 CP。

这个试验只回答原始 t0 的固定算子对边界外尾贡献是否敏感。不验证未知演化场的外尾闭合，不把任何校正回填到 t>0，不替代 B/C/D/E 的原预登记终点比较，不宣称强奇异性或慢收敛已被突破。

## root 启动层补验

旧dormant脚本未显式cd；为避免从维护目录运行时的cwd包优先级风险，实际true执行改用D/acceleration_lab/original_t0_tail_root_launcher_v1/run_native_A_H8.m与run_native_H32.m。原数值helper及116份源/输入未改，显式冻结cwd/which/依赖/SHA的两支false补验48643已exit0，CodeAnalyzer=0，无LU/PDE。该补验不是原脚本已经跑错数值源的证据；真正算子比较仍待资源槽。
