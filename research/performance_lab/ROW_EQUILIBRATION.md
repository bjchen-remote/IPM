# Poisson 行均衡候选，2026-09-08

候选已通过 27 组小网格 MMS 比较和 6 组受控行单位压力测试。二次幂行缩放
在这些 MMS 中使流函数和两个速度分量的最大数值差均为零。人为引入 16 个
数量级的行单位差异后，缩放使 decomposition 的 RCOND 从 `3.74e-17` 变为
`1.37e-4`；解误差和原方程残差仍在相近的舍入误差水平。当前证据支持继续
隔离验证这一表示变换，不能宣称它已改善真实 H=1e6 网格的精度或运行速度。

本实验没有修改共享 `+ipm`、生产规范、边界条件、可信门或 checkpoint；
没有运行 q512，也没有推进 PDE。计算线程保持 10。

## 当前矩阵与缓存

`+ipm/+mesh/build.m` 对四阶和六阶分别用 7 点和 9 点物理坐标二阶差分矩阵，
抽取内点块 `Tx`、`Ty`，按 MATLAB 列优先的 `ny × nx` 排列形成

```matlab
A = kron(speye(nxi),Ty) + kron(Tx,speye(nyi));
ops.A = A;
ops.poisson = decomposition(A,'lu');
```

`ops.A`、`ops.Tx`、`ops.Ty` 是稀疏数值矩阵；`ops.poisson` 是缓存的
`decomposition` 对象。各向同性 RK stage 仅回代，不重新分解。
非均匀网格的二阶导数系数随局部长度平方的倒数变化，还受局部 stencil
比例和单侧闭合影响。因此远近网格尺度悬殊会产生很大的行系数差异；单靠
坐标范围 H 不能准确推断矩阵条件数或物理解误差。

`+ipm/+field/poisson.m` 先把两侧和上下边界贡献减进 `rhs`，随后才回代。
目前四/六阶残差使用原始 `poissonOperator` 与完整 `rhs`：

```matlab
norm(poissonOperator*interior-rhs(:),inf) / max(norm(rhs(:),inf),eps)
```

该检查必须保留。人工行缩放会改变全局范数权重，仅检查缩放后的残差可能
掩盖原方程某些行的误差。实验额外记录分量后向误差
`max(abs(A*x-b)./(abs(A)*abs(x)+abs(b)))` 和公共行归一化残差，用于比较；
这些指标没有替换生产门。

## 等价候选与最少集成位置

设 `m_i = max_j abs(A_ij)`，选择正对角矩阵 `D`。普通 row-max 取
`D_ii=1/m_i`。优先研究的二次幂版本用 `[~,e]=log2(m)`，取
`D_ii=2^(-e_i)`，使缩放后每行最大模落在 `[0.5,1)`。在不发生溢出或
下溢损失时，二次幂乘法只改变二进制指数。候选对稀疏非零数和系数往返
一致性作了显式检查；这不保证后续 LU、pivot、迭代改进或 PDE 轨道逐位不变。

`A*x=b` 与 `(D*A)*x=D*b` 数学等价，没有改未知量、网格、PDE 或源项。
实际缓存只有一个缩放后的 factor 和长度为内点数的 scale 向量，不保存
额外的缩放稀疏矩阵。研究接口：

```matlab
[cache,info] = ipm_perflab_row_factor(A,'power2_row_max');
[interior,metrics] = ipm_perflab_row_solve(cache,completeRhs,A);
```

若后续独立源码版本验证成功，最少生产改动仍在现有唯一路径：

1. `mesh.build` 的四/六阶 LU 位置保存 `ops.poissonRowScale`，直接构建
   `decomposition(spdiags(scale,...)*A,'lu')`，避免先构建未缩放 LU。
   `ops.A`、`Tx`、`Ty`、边界系数和导数矩阵全部保留原值。
2. `field.poisson` 的四/六阶 cached 回代只将已完成边界装配的 `rhs(:)`
   乘相同 scale。不得只缩放 `source`，也不得对边界贡献漏乘或重复乘。
3. 继续以原矩阵和原完整 RHS 计算残差、返回 operator 和 interiorRhs；
   全部 flow、诊断、接受门、时间推进继续走维护实现和 `ipm.solve`。

各向异性 `kappa~=1` 使用 `(1/kappa)*kron(Tx,I)+kappa*kron(I,Ty)`。
若覆盖该分支，必须从当前 kappa 的算子计算当前 scale/factor；不能复用
各向同性的 factor。当前实验仅验证 `kappa=1`，未验证各向异性。
二阶路径已有积分权重和对称化后的 Cholesky，任意左行缩放会破坏其对称性；
该路径不属于本候选范围。

构建时需要一个临时缩放稀疏矩阵和 scale 向量，但不应同时保有原、新两个
LU。现有换网格事务仍应先释放旧 `ops.poisson`。在研究 MMS 中原 oracle
factor 和小型候选 factor 有意共存，不能用该测试的内存推算生产峰值。

checkpoint 不保存 ops，恢复时会重建分解。改变分解路径可能改变恢复流场
和后续浮点轨道，因此不能把“等价方程”当作旧 checkpoint exact 契约的
自动证明。真实网格比较、恢复/打包/严格回读、短程分支和原门验证仍需在
隔离版本中完成；目前没有生产默认改动。

## 实测与限制

脚本 `ipm_perflab_verify_row_equilibration` 使用四阶均匀、四阶拉伸
`[0.9,0.7]` 和六阶均匀三类网格，各有 `17²/33²/49²` 三档。
MMS 为 `exp(0.3*x)*cos(1.2*x)*y*(1-y)`，显式非零人工侧边界验证完整
RHS 装配，物理底边和顶边为零。所有方法保持边界数值完全相同。
三种算法分别是原 LU、二次幂 row-max、普通 row-max。

`49²` 的相对 MMS 最大误差如下；二次幂候选与维护解的数值差为零，普通
row-max 只产生舍入级差异。该特定光滑函数不足以单独建立一般收敛阶。

| 网格 | 流函数相对误差 | 速度相对误差 |
|---|---:|---:|
| 四阶均匀 | 5.07e-11 | 4.54e-9 |
| 四阶拉伸 | 2.04e-10 | 1.46e-8 |
| 六阶均匀 | 1.94e-13 | 1.73e-11 |

这些小网格的行最大系数跨度只有约 `0.13–0.60` 个数量级，不能代表生产
远近尺度。另将最后一个小网格的矩阵和 RHS 左乘人为行单位矩阵，分别增加
8 和 16 个数量级跨度，保持已知离散解。该压力试验隔离了行单位因素，
没有模拟生产网格的截断误差、远场边界误差或 stencil 病态。

| 16 阶数量级行单位试验 | decomposition RCOND | 离散解相对误差 | 原方程相对残差 | 分量后向误差 |
|---|---:|---:|---:|---:|
| 原 LU | 3.74e-17 | 5.60e-15 | 5.42e-15 | 1.68e-16 |
| 二次幂 row-max | 1.37e-4 | 4.09e-15 | 5.56e-15 | 1.70e-16 |
| 普通 row-max | 1.37e-4 | 1.42e-14 | 7.70e-15 | 1.83e-16 |

全部 warning 保持开启。原 LU 在 16 阶行单位试验触发近奇异警告，但仍给出
舍入级准确的已知离散解。RCOND 的明显变化不能据此解释为物理解精度的
同等改善，也不能据此断言生产 RCOND 告警无害。

只做了每档 5 次缓存回代、一个初始化的共享机器微计时。例如 `49²` 四阶
拉伸的原/二次幂 factor 为 `6.191/6.220 ms`、回代中位为 `0.415/0.394 ms`，
二次幂预处理为 `0.366 ms`；六阶均匀回代为 `0.559/0.555 ms`。数据不足以
证明稳定加速。16 阶行单位原/二次幂回代 `0.932/0.574 ms` 中原路径含警告
输出成本，不能当作纯数值内核加速；初始化也仅有一个样本。

全部结果保存于 `result/verification/performance_lab_20260908/row_equilibration_small.mat`：
`report.allChecksPassed=true`，27 个 MMS 条目和 6 个行单位压力条目。
三个新增 MATLAB 文件 `checkcode` 均为 0。证据范围是小网格椭圆求解，
没有真实 q512 精度、内存、整步速度或长期轨道证据。

## 一手资料与解释

MathWorks 明确给出 sparse LU decomposition 包含行缩放：
`P*(R\A)*Q=L*U`。因此当前实现并非未经任何均衡的稀疏 LU，外部再做 D
可能主要改变输入表示、条件估计或警告行为，而非获得一个全新稳定机制。
[MathWorks decomposition](https://www.mathworks.com/help/matlab/ref/decomposition.html)

MathWorks 的通用 `equilibrate` 还会使用排列和列缩放；文档提醒，缩放系统
上的小残差并不自动保证原系统上的小残差。本候选只使用行缩放，保留原始
残差验证。[MathWorks equilibrate](https://www.mathworks.com/help/matlab/ref/equilibrate.html)

本机 R2026a 可读 MATLAB 层显示 sparse LU 交给 UMFPACK wrapper，rcond
也委托该 wrapper。SuiteSparse 的 UMFPACK 官方说明描述其内部行缩放，
并将 UMFPACK 自身的 RCOND 信息定义为 U 对角幅值比的粗略指标。MATLAB
闭源 wrapper 的确切映射未在本实验验证，故报告只标注公开
`rcond(decomposition)` 的返回值，不把它直接等同于高精度范数条件数。
[SuiteSparse UMFPACK 官方接口说明](https://raw.githubusercontent.com/DrTimothyAldenDavis/SuiteSparse/dev/UMFPACK/Include/umfpack.h)
