# Kronecker-sum 求解候选初筛，2026-09-08

**当前拒绝将普通 dense 实 Schur 或特征分解候选送入 q512 冻结测试。**
小型温和网格上的速度和缓存尺寸有吸引力，但同一维护差分矩阵的拉伸阶梯
已经显示小行相对精度损失：全局原方程残差仍约 `1e-14` 时，Schur 解可有
`1.9e-4` 的已知离散解误差和 `4.0e-6` 的分量后向误差。该失败发生在
direct LU 仍明显更准确的阶梯，不能仅以全局残差合格接受候选。

这是研究筛选结果，没有修改 `+ipm`、grid_lab、PDE、Green 边界、可信门或
checkpoint；没有进行 q512 或 PDE 时间推进。全部数值实验为 `17²–65²`，
MATLAB R2026a，10 计算线程，共享机器短计时。

## 相同离散问题与实现

维护矩阵是 `A=kron(Ix,Ty)+kron(Tx,Iy)`。对内点矩阵 `X`，有

```matlab
A*X(:) = reshape(Ty*X + X*Tx.',[],1);
```

因此边界装配后完整 `B` 的问题等价于 `Ty*X+X*Tx.'=B`。候选只替换
内点线性代数，不改 `Tx/Ty`、未知量的排列、源项或边界值。

- `ipm_perflab_sylvester_factor`：缓存 direct LU、实 Schur、一般非对称
  特征分解，或供 MATLAB 原生 `sylvester` 使用的 1D dense 矩阵。
- `ipm_perflab_sylvester_solve`：接收完整内点 RHS；不生成 Green 边界。
- `ipm_perflab_verify_sylvester`：用维护 `field.poisson` 装配 RHS 和提供
  direct 解；插回相同边界后比较 psi、速度、MMS、原方程残差及分量后向误差。

实 Schur 路径分解 `Ty=Qy*Sy*Qy.'`、`Tx.'=Qx*Sx*Qx.'`，求
`Sy*Z+Z*Sx=Qy.'*B*Qx`，再还原 `X=Qy*Z*Qx.'`。三角方程调用公开
`sylvester`；本机可读实现会识别实 Schur 输入并跳过再次 Schur 化。
候选没有调用 MATLAB internal API。直接对原始 `Ty/Tx.'` 调用公开
`sylvester` 得到相同的数值差和失败行为，排除了本候选缓存变换造成问题
这一解释。

快速对角化保留一般非对称特征向量，使用 LU 解特征向量线性系统；没有
假设正交性，没有丢弃复数虚部，也没有通过显式取实部掩盖误差。
它在极端拉伸下没有解决稳定性问题，部分条目更差。

## 已维护网格与 Green/MMS 比较

9 组实际通过 `mesh.build` 的网格：四阶 `n=17,33,65` 均匀，四阶 `n=65`
拉伸 2、4；六阶 `n=33,65` 均匀，六阶 `n=65` 拉伸 2、4。每个网格包含：

1. 显式非零侧边界的光滑 MMS，底边和顶边为零。
2. `rho=-exp(-X.^2-Y.^2)`、`source=rho*Dx.'` 和维护 image-Green 边界。

4 种线性代数方法共 72 项，全部通过研究筛选比较；这些不是生产接受门。
所有边界数值完全相同。例：四阶 `65²`、拉伸 4 的 MMS，direct psi 相对
误差 `3.4933e-8`、速度相对误差 `1.2688e-6`；缓存 Schur 对 direct 的
psi/速度差 `9.55e-13/1.88e-12`。同网格 Green 问题的 psi/速度差为
`3.97e-13/3.06e-12`。微小网格通过不足以证明大尺度分离下可靠。

## H=1e6 拉伸阶梯：已知离散解

另用 `65²` sinh 轴 `x=H*sinh(s*sx)/sinh(s)`、
`y=H*sinh(s*sy)/sinh(s)`，`H=1e6`，`s=0,4,8,12,16,20,24`。
每个轴仍使用维护四阶 7 点 `fdMatrix` 形成完全相同的 Kronecker-sum。
给定一个在计算坐标上平滑、含较高频分量的内点矩阵 `Xexact`，由
`B=reshape(A*Xexact(:),...)` 构造已知离散解问题。

这是原始轴/差分矩阵的代数压力试验，**没有声称这些高拉伸微小轴全部通过
生产网格门，也没有把它们当作 H=1e6 PDE 的离散精度试验**。它比人为行
单位缩放更直接地测试实际非均匀 stencil 的 Kronecker 分解。

| 拉伸 s / 行最大系数跨度 | direct 已知解相对误差 | 实 Schur 已知解相对误差 | 实 Schur 原方程全局残差 | 实 Schur 分量后向误差 |
|---|---:|---:|---:|---:|
| 4 / 2.96 阶 | 3.39e-15 | 7.91e-13 | 3.10e-13 | 2.76e-14 |
| 8 / 6.24 阶 | 1.21e-13 | 7.67e-11 | 9.23e-14 | 1.97e-11 |
| 12 / 9.39 阶 | 2.91e-12 | 4.78e-8 | 8.48e-14 | 2.38e-9 |
| 16 / 12.51 阶 | 5.91e-9 | 1.90e-4 | 5.50e-14 | 3.99e-6 |

`s=12` 的普通特征分解误差为 `9.47e-8`，分量后向误差 `1.22e-8`；
`s=16` 已知解误差 `4.46e-4`。它没有提供有用的替代。

`s=20/24` 下 direct 对已知解也失效，因此该两档不能把 direct 当作
准确真值；MAT 文件保留完整失败数据，不据这两档判定某候选更准确。
拒绝普通 Schur 的证据已在此前的 `s=12/16` 足够明确。

## 精度损失位置与残差解释

`s=16` 的 `Tx.'` Schur 重构 `Qx*Sx*Qx.'`，全局相对矩阵误差仅
`6.47e-15`，但按各原始行的最大系数归一化后，最大重构误差为 `7.87e-4`。
因此全局范数下小扰动未保留远场小行的相对精度。这与 Schur 变换在巨大
系数跨度下丢失小尺度信息的解释一致；本实验未用高精度特征值参考确定
具体哪一个小特征值最先受损。

报告的 `left/rightComponentwiseReconstruction` 字段实际定义为按原始
行最大系数归一化的矩阵重构缺陷，包含原结构零位上产生的填充；它不是
逐元素相对误差。方程的 `componentwiseBackwardError` 则严格使用

```matlab
max(abs(A*x-b)./(abs(A)*abs(x)+abs(b)))
```

并以 `realmin` 处理零分母。原方程全局残差保持维护表达式中的 `eps`
分母下界。另有公共行归一化残差；`s=16` Schur 为 `2.61e-4`。
这些补充证据揭示了全局最大行主导的残差可能掩盖的小行误差，没有用于
放宽或替代现有生产门。

## 缓存、初始化与时间成本

令内点尺寸为 `m=ny-2`、`n=nx-2`。缓存实 Schur 的四个实 dense 矩阵
总计 `16*(m²+n²)` 字节，初始化复杂度 `O(m³+n³)`。每个 RHS 的四次
稠密基变换约 `4*m*n*(m+n)` 浮点运算，另有三角 Sylvester 求解。

| 内点尺寸 | 四矩阵缓存 | 一个场 | 仅四次变换 FLOPs |
|---|---:|---:|---:|
| 63×63 | 0.121 MiB | 0.030 MiB | 2.00e6 |
| 511×1023（q512） | 19.95 MiB | 3.99 MiB | 3.21e9 |
| 1023×2047 | 79.91 MiB | 15.98 MiB | 2.57e10 |

这不含原始 ops、构建/BLAS 工作区、RHS 及变换中间数组；不是实测 RSS。
q512 仅按维度估算，没有构造或分解该尺寸矩阵。普通 sparse LU 的实际
内存依赖排序和 fill，decomposition 内部存储在 MATLAB `whos` 中不透明，
本报告没有虚构 LU 内存值。

四阶 `65²`、拉伸 4 的一次设置与 3 次回代样本：

| 方法 | 设置 | MMS 回代中位 | Green 回代中位 |
|---|---:|---:|---:|
| 缓存 sparse LU | 13.51 ms | 1.067 ms | 1.054 ms |
| 缓存实 Schur | 1.434 ms | 0.085 ms | 0.083 ms |
| 每次原生 sylvester | 0.026 ms（仅存原矩阵） | 1.020 ms | 0.985 ms |

原生 sylvester 的分解成本位于每次调用中。小网格缓存 Schur 的收益不能
外推 q512：四次变换从约 2e6 增为 3.2e9 FLOPs，且已经有精度失败。
快速对角化计时受多线程/复数 dense 解影响明显，未作为性能推荐。

若未来重新考虑这种结构，应先得到能够保留远场小行相对精度的独立算法
及其验证；仅把原矩阵交给现有 dense Schur、普通 eig，或仅做全局缩放，
不满足本次进入真实 q512 比较的依据。当前继续保留维护 sparse LU 路径。

## 证据与资料

`result/verification/performance_lab_20260908/sylvester_small.mat` 保存
72 项维护网格比较和 28 项拉伸阶梯条目；失败条目全部保留，未生成
“全通过”的总标志。三个新增 MATLAB 文件 `checkcode` 均为 0。

MathWorks 公开 `sylvester` 的矩阵方程和 Kronecker 表达，以及只接受
dense 矩阵的接口要求。[MathWorks sylvester](https://www.mathworks.com/help/matlab/ref/sylvester.html)

MathWorks 公开实 Schur 的正交相似分解和 1×1/2×2 准三角块形式；候选
直接使用该形式及公开求解接口。[MathWorks schur](https://www.mathworks.com/help/matlab/ref/schur.html)
