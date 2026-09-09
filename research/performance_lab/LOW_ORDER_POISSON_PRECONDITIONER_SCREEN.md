# 原高阶 Poisson 的低阶预条件小网格筛查

本候选目前只显示了较小的因子存储潜力，没有获得速度或大场精度资格。原高阶矩阵 A、原 Dirichlet/Green RHS 不变；GMRES 解的是原 A，不把低阶矩阵当成替代 PDE。维护 +ipm、运行源、checkpoint 与全部生产门未改。

## 已有证据与本轮区别

- `FROZEN_FLOW_PROFILE.md` 的 fresh895 实际排他计时，缓存回代约占整 RHS 的72.9%。Poisson矩阵及分解已复用，不能再次把重复build当热点；每个stage有新rho，不能无改变算法地合批未来 RHS。
- `SYLVESTER_SCREENING.md` 已拒绝普通dense Schur/eig：原方程全局残差约1e-14仍伴随小行误差和1.9e-4已知解误差。本轮不使用它们。
- `FRESH_RHS_COST_AND_LU_SCREEN.md` 已拒绝公开LU三角缓存作bitwise替换。它不保证原native refinement的结果；本轮只把公开L/U的whos存储用作明示的fill代理，不用它求解。
- `ROW_EQUILIBRATION.md` 表明行单位均衡不能自动改善物理解或速度；`DECOMPOSITION_WARNING_AUDIT.md` 已证实条件估计缓存，只有警告输出的有限开销。
- 本轮检索 research 中没有找到同一高阶 A、低阶五点Poisson预条件GMRES或线方向多重网格的已执行试验。已有各向异性低阶PCG使用原isotropic factor，属于另一问题。

## 预登记算法与检查

M为相同非均匀轴上的三点二阶负Laplacian Kronecker和。以低阶单元体积 W 加权，Q=(WM+(WM)')/2，使用公开稀疏Cholesky及其置换缓存；预条件作用是近似M^{-1}r=Q^{-1}Wr。高阶A仍是维护七点非对称闭合，未被对称化。每个未知量、源和边界处理不变。

固定 GMRES restart30、tol1e-12、最多10个outer；第三方法在初次GMRES后最多做两次完整原残差 b-Ax 的GMRES缺陷修正。没有针对输出调参或扩大阈值。MATLAB的relres/resvec属于预条件系统；本报告单独计算原方程逐分量后向误差与行归一残差，不用该内置relres作精度结论。[MathWorks GMRES](https://www.mathworks.com/help/matlab/ref/gmres.html)、[MathWorks Cholesky](https://www.mathworks.com/help/matlab/ref/chol.html)

三个尺寸49²/65²/129²，每个有H8均匀、H8拉伸和H1e6拉伸三组，s=(N-1)log(1.07)/2。九对轴均实际通过原六项轴质量门（含1.08/.01/rcond1e-9/全局q1e-8/局部q[.35,1.65]）及mesh.build。这是几何与线性代数测试，不是经过核心/transfer验收的PDE轨迹。每组有已知离散解、非零边界连续MMS和原image-Green RHS。连续MMS使用显式边界oracle，速度误差是同维护Dx/Dy的原始导数，不声称该强制非零壁边界属于原IPM初值。

另有65²/H1e6/s12、s16压力轴，明确未过原轴质量门，不能当合法生产网格。它们使用相同维护七点A与已知离散解，仅审查大量级行跨度风险。已知离散右端是浮点A*x；没有高精度真解，故同时保留native误差和candidate-native差，不从单个误差值推出稳定性定理。

阈值在运行前写入registration：逐分量后向误差1e-12、行归一残差1e-10、psi相对native差1e-9、速度差1e-8、离散已知误差1e-8。native已知误差超过1e-8时不能用作合格oracle。所有失败原样保存。

## 实际结果

session60099 exit0；数值部分7.49809秒，含启动的受监控进程18.3613秒，所有LU已释放。没有PDE、大矩阵、checkpoint读写。最大系统维数16129未知量。聚合child RSS峰1,159,397,376 B，低于1.5GiB硬限；time-l峰memory footprint为1,407,191,608 B。child0swaps不代表系统未换页，本轮不据此给并行大场资源资格。

| 方法 | 27个轴质量通过的问题满足场/方程门 | 最大逐分量后向误差 |
|---|---:|---:|
| 原native LU | 27 | 2.50e-16 |
| 低阶预条件GMRES | 8 | 5.72e-11 |
| 上述GMRES加最多两次缺陷修正 | 27 | 2.69e-16 |

普通GMRES的内置flag全部为0，仍有19项原逐行误差门失败，说明不能使用预条件残差替代原问题的验收。修正后的27项场/方程门通过不等于每个内层求解都收敛：17次修正solve返回原flag3（停滞），警告和resvec全部保留，未改成0。修正后最大psi/速度对native差分别5.28e-14/4.57e-13。

| 65²强拉伸压力轴 | 方法 | 已知解相对误差 | psi对native差 | 原逐分量后向误差 |
|---|---|---:|---:|---:|
| s12 | 原LU | 8.49e-13 | 0 | 1.86e-16 |
| s12 | GMRES+两修正 | 2.09e-12 | 2.25e-12 | 1.76e-16 |
| s16 | 原LU | 4.54e-9 | 0 | 1.75e-16 |
| s16 | GMRES | 4.48e-6 | 4.48e-6 | 4.28e-12 |
| s16 | GMRES+两修正 | 1.073e-8 | 6.19e-9 | 2.17e-16 |

s16两修正仍违反预登记已知解/psi差门，而且修正内层flag包含1和3。即使分量后向误差已接近舍入极限，也不能忽略病态系统的前向差异或宣布大跨度精度合格。

## 存储与成本的实际边界

129²的低阶Cholesky结构缓存whos约9.13MB，31向量GMRES基底的简单下界约4.00MB；公开高阶L/U及排列/缩放的whos代理约83.8MB。低阶L的nnz350112，高阶公开L+U约5,181,278。该差距是稀疏fill的实际量级证据；GMRES还有其它工作区，公开LU代理也不等于native decomposition内部缓存或整个应用RSS，不能由此推定961×321/641×481实际峰值。

下面只是每网格三个RHS单次样本的中位成本，不是随机ABBA性能认证。首例有JIT启动影响，没有因不利样本选择性重测。

| H8拉伸网格 | 原回代 | GMRES | GMRES+两修正 |
|---|---:|---:|---:|
| 49² | .661ms | 3.787ms | 16.995ms |
| 65² | 1.245ms | 5.579ms | 22.712ms |
| 129² | 10.368ms | 18.545ms | 86.993ms |

因此本轮无速度通过结论；固定两次修正还在已达到舍入极限时花费额外工作。也不能用小格回代成本直接外推大场LU，但当前证据不支持启动大型operator。

## 下一项可判定的有限方案

若继续此方向，应先独立登记更强、仍实际通过六质量门的129²方向拉伸轴（例如y拉伸与x分开，以实际门结果为准），加局部与远尾离散已知模式，并保留本轮s16失败。把“达到原分量/行残差门即停止、最多两修正”的精确停止条件单独注册，不能只提高GMRES容差或让最后flag掩盖原方程误差；仍须逐项比较psi/导数/已知误差。

线方向多重网格可作为以后替换M^{-1}的独立候选，有进一步减少五点factor存储的潜力，但本轮精确低阶求解本身已不能让所有强跨度问题通过，故没有直接上MG或宣称其能修复前向误差。它应先证明对原非对称A的收敛与残差/速度门，再讨论大场工作区和实测速度。当前原native LU继续作为运行路径。

## 原始文件

`result/verification/autonomous_runtime_20260909/preconditioned_poisson_small_v1/` 保存冻结source、登记、helper、watchdog、全部原解/残差历史、逐网格矩阵与代理因子、原日志、资源时序、compact_summary和final_audit。123个登记文件SHA不变，helper CodeAnalyzer0。无旧失败被覆盖。研究入口为 `ipm_perflab_low_order_preconditioned_poisson`；使用固定目录时必须另选新的输出目录。
