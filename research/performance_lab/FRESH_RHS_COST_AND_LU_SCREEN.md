# 最新fresh895 RHS成本与等价优化筛查

本轮未找到通过全部条件、值得进入大网格的等价加速候选。没有新增大LU或PDE：Green/导数成本直接读取最新crop step81完整result；LU候选只用≤65×33小矩阵。已释放所有研究进程。

## 可证的成本尺度

默认10线程WENO实际ABBA后段的维护whole RHS中位为513.594687毫秒。以下为同一真实场分别计时所得尺度比，不是同一个profiler采样的排他性占比，不能把差值直接称作LU耗时。

| 已测部分 | 毫秒 | 相对whole RHS尺度 |
|:---|---:|---:|
| 4个维护WENO kernel | 42.8413 | 8.341% |
| 维护Green边界 | 11.3966 | 2.219% |
| rhoX和两速度导数 | 1.92973 | 0.3757% |
| 原A乘psi残差作用 | 1.66383 | 0.324% |
| 未归属差值（含LU、规范和其它开销） | 455.763 | 88.74% |

WENO稳定段四kernel合计42.8413ms，约8.34%；此前进展中粗估41.4ms/8.1%在此更正。源码显示Lu回代是剩余部分中最合理的重点，但目前没有真实缓存的排他性solve计时，不能把88.7%全部归到LU。

`ipm_perflab_field_cost_no_lu`按维护7点fdMatrix重建一维Dx/Dy，rho*Dx与保存omega、两个psi导数与保存velocity全部逐位一致；原Kronecker A共4,450,194非零元，原残差1.8816828106664292e-9逐位一致。只组装A，不分解、不调用poisson，不演化。八轮交错正/反操作次序，全部原时样保存于field_cost_step81_v1.json。

## 已确认不能使用的伪消重

isotropic kappa=1的poissonOperator仅返回ops.A，没有每stage重建矩阵。ops.poisson是一次decomposition(A,'lu')缓存；quadratic/anchor每flow只有一个Poisson RHS。SSPRK54阶段顺序依赖新rho，不能不改方法就合批未来stage RHS。接受态RHS缓存与高阶残差复用均已存在，不重复计算第二个稀疏残差。Green源筛选、signed centroid压缩依赖当前source，不能跨stage固定核或把左右边界强制互为负号。

读取本机MATLAB R2026a公开decomposition实现：set.CheckCondition已保存rcondSaved，warnIfIllConditioned只读取该缓存；没有每次重估rcond可省。CheckCondition=false不构成量级候选，未修改它。原内部SparseLU使用UMFPACKWrapper；其原refinement行为没有被关闭或借内部API改写。

## Green精确复用：正确但更慢

`ipm_perflab_green_candidate`只在调用内复用四个平方距离；精确对称double-odd左右边界共享距离，但分别计算原ratio/log/GEMV，不以left=-right替代原浮点运算。源选择、signed centroid、权重、chunk128均保持。

48个小尺寸/均匀与拉伸/范围/符号/边界模式测试全部bitwise。真实输入用维护quadrature重建最小Green ops，先逐位匹配保存psi的全部外边界（遵循原角点覆盖顺序），再比较候选，24个6组随机ABBA计时输出全部bitwise。原Green中位11.3966ms，候选14.1114ms，速度比.807617x，六块范围.694792–.998870。即使这部分大幅优化，占whole RHS尺度也只约2.2%，因此拒绝继续。

## 公开LU三角缓存：潜在速度但不满足精确替换

可独立检验的候选使用公开 `lu(A,'vector')` 返回L/U、行列排列及D，缓存公开triangular decomposition，按照原矩阵的 `P*(D\A)*Q=L*U` 关系求解，A和b不变。此契约由[MathWorks lu文档](https://www.mathworks.com/help/matlab/ref/lu.html)给出。没有使用内部UMFPACK接口，也没有设置未公开AllowIterativeRefinement。L/U中间因子的condition warning不代替A的条件，原A rcond逐项记录。

独立 `ipm_perflab_public_lu_tiny` 测33×17、65×33和stretch0/4/12，使用维护7点矩阵与已知离散MMS，比较原native decomposition、公开三角解的bitwise/前向误差/原方程componentwise backward error，并保留4组ABBA样本。stretch12是故意病态的代数压力网格，不声称通过生产mesh质量门。

| 网格nx | stretch | bitwise | 速度比 | 相对解差 | native分量后向误差 | candidate分量后向误差 |
|---:|---:|:---:|---:|---:|---:|---:|
| 33 | 0 | 否 | 0.696926 | 1.44343e-15 | 1.39797e-16 | 5.31138e-16 |
| 33 | 4 | 否 | 0.685643 | 1.14794e-13 | 1.72959e-16 | 1.0525e-14 |
| 33 | 12 | 否 | 0.675467 | 314.554 | 1.30773e-13 | 4.08364e-11 |
| 65 | 0 | 否 | 1.58628 | 5.21932e-15 | 1.65845e-16 | 1.05119e-15 |
| 65 | 4 | 否 | 1.90723 | 8.15103e-14 | 1.83637e-16 | 2.05197e-15 |
| 65 | 12 | 否 | 2.1524 | 1.79426e-10 | 1.76148e-16 | 2.5078e-13 |

全部六例不bitwise；65×33虽有1.59–2.15x三角解速度比，不能当作合格整体增益。病态33×17例公开解相对native差达到315，candidate分量后向误差4.08e-11；65×33强拉伸例后向误差也从1.76e-16升到2.51e-13。公开因子/三角路径不保证原UMFPACK solve/refinement的浮点行为，因此按精确替换要求拒绝，不启动大LU或PDE验证。不能拿小的global residual掩盖这些差异。

下一次若主线已有可用cache窗口，最少额外工作是直接分离原ops.poisson\rhs、原Green、gauge的同一flow计时，以核实剩余热点；当前证据尚不支持改solver或停止原生refinement。dense Sylvester早先拒绝结论保持，不重新列为合格方案。
