# 全方向 \(\rho_x\) 衰减条件下的 IPM 爆破：现状、约束与新候选

本文固定物理约定

\[
\rho_t+u\cdot\nabla\rho=0,
\qquad
-\Delta\Psi=\rho_x,
\qquad
u=\nabla^\perp\Psi=(-\Psi_y,\Psi_x).
\tag{D0}
\]

本文所说的远场硬条件是二维一致意义下

\[
\rho_x(t,\cdot)\in C_0(\overline{\mathbb H}),
\qquad
\lim_{|(x,y)|\to\infty,\ y\geq0}\rho_x(t,x,y)=0
\quad (t<T).
\tag{D1}
\]

除非特别标为“退化边界例”，还要求真正爆破的量就是水平导数：

\[
\rho_x\not\equiv0,
\qquad
\|\rho_x(t)\|_{L^\infty}\longrightarrow\infty
\quad(t\uparrow T).
\tag{D2}
\]

这一区分是必要的；若只要求某个梯度分量爆破而允许
\(\rho_x\equiv0\)，问题会被一个与 self field 无关的被动例子平凡化。

## 1. 当前结论一览

| 类别 | 光滑初值 | \(\rho_x\in C_0\) | \(\rho_x\) 爆破 | 标准衰减 Biot--Savart | 判断 |
|---|---:|---:|---:|---:|---|
| 旋转/平移直条带精确轨道 | 是 | 否 | 是 | 否 | 沿条带切向不衰减 |
| 水平 compact front | 是 | 否 | 是 | 否 | 独立于 \(y\) |
| affine 锐扇核 | 局部是 | 否 | 是 | 否或锐扇域 | \(\rho_x\) 是空间常数 |
| singular \(K_0/K_1\) Bessel 稳态 | 否 | 是 | 否 | 可在穿孔域成立 | 每个时刻 tip 已奇异且静态 |
| Dembski 锐角域定理 | 仅 Lipschitz 于角点 | 是 | 是 | 是（锐角域） | 严格爆破，但不是全局光滑初值 |
| 边界跳跃锐扇 Riccati 解 | 上侧 (C^{1,\alpha})、角外光滑；奇延拓跳跃 | 是 | 是 | 是（锐角域） | **有限质量幂律尾的严格边界驱动爆破** |
| 全平面 analytic affine contract | 密度初态是 | 是 | 是，Type I | 否 | **精确轨道；预设奇异远场合同** |
| 半平面 odd affine contract | 是 | 是 | 是，Type I | 否 | 峰逃向 \(y=+\infty\)，不是固定点奇点 |
| 被动 transverse front | 是 | 是，且恒为零 | 否；只有 \(\rho_y\) 爆破 | 否 | 仅作逻辑边界例 |
| compact-right-fan leading tail | 渐近数据，非初值 | 是 | 未证明 | 仅满足 leading Poisson | 精确相容，含 signed return |
| 本文 \(\beta=0\) fixed point | 若存在则是 | 若存在则是 | 若存在则 Type I | 必须满足 self moment | finite-level-area 静态支已排除 |
| signed-head/edge/Casimir-cap | 目标是 | 目标是 | 目标 Type I | 目标是 | **首选非定常未闭合 gluing** |

因此，若只按字面要求 (D1)--(D2)，并允许预先指定非衰减且在 $T$ 奇异的
harmonic far-field contract，则本文已经得到密度初态为解析 Schwartz 的精确
Type-I 轨道，见第 7 节。它不是由密度 Cauchy 数据自动选出的正例。
若再要求“标准衰减 Biot--Savart、半平面、固定有限点、全局光滑初值”，目前仍
没有完成存在性证明。若只放宽角点光滑，则 Dembski 已给出紧支撑严格爆破定理。
新的进展是：任意右锐扇的衰减 leading tail、局部自生 signed head 和远端
Casimir cap 分别都能严格构造；静态 compact fixed point 与静态光滑 wall edge
则各有一个严格障碍。因此目前最短路线已经从“找单一光滑稳态”改成非定常
head/edge/cap gluing。这些合同和证据等级必须分开，不能互相偷换。

若保留本项目指定的 \(x_2=0\) 非零上侧迹，并允许第二条扇边是
不可穿透锐角壁，则有一个更强的完全正例：角点 Poisson Hessian 恒等式
严格闭合为 \(a_t=(\tan L/2)a^2\)。
`analysis/EXACT_BOUNDARY_JUMP_POWER_FAN_BLOWUP_ZH.md` 给出光滑有限质量
\(r^{-5}\) 幂律尾初值、底边跳跃、右扇正核/负 return 和有限时梯度爆破的
完整证明。它解决边界驱动锐扇合同，不解决无第二扇边的标准半平面问题。

## 2. 已知严格爆破与它没有解决的部分

Dembski 的锐角域结果
[arXiv:2511.01827](https://arxiv.org/abs/2511.01827)
构造了边界密度为零的有限时间奇点，并可把密度在空间无穷远截成紧支撑。
其主导结构为

\[
\rho(t,r,\theta)=rP(t,\theta),
\qquad
\Psi(t,r,\theta)=r^2G(t,\theta),
\tag{D3}
\]

其中角向系统出现 \((1-t)^{-1}\) 爆破。紧支撑立即给出 (D1)。但是
\(rP(\theta)\) 在顶点一般只到 Lipschitz，且方程定义在严格小于半平面的锐角域；
它不是从全局 \(C^\infty\) 半平面初值产生的结论。这个定理证明了“远端衰减和
无 boundary mass 并不排除 IPM 爆破”，同时把真正剩余问题压缩为：如何用一个
动态内尺度把顶点的方向依赖 Lipschitz 核正则化，而不破坏角向双曲机制。

作为逻辑边界，令 \(\tau=T-t\)、\(f\in C_c^\infty(\mathbb R)\)，则

\[
\rho=\rho_c+f(y/\tau),
\qquad
\Psi=-\frac{(x-x_*)y}{\tau}
\tag{D4}
\]

逐点满足 (D0)，并有

\[
\rho_x\equiv0,
\qquad
\|\rho_y\|_\infty=\frac{\|f'\|_\infty}{T-t}.
\tag{D5}
\]

它说明只写“\(\rho_x\to0\)”还不够；本文不把 (D4) 计作所求水平梯度爆破。
其压缩完全来自预先指定的奇异 harmonic strain，标准衰减 Biot--Savart 会把该
调和部分选为零。

## 3. \(\beta=0\) bounded-contrast：fixed point 审计与非定常自屏蔽路线

令

\[
L(t)=c(T-t),
\qquad
a'(t)=b,
\qquad
z=(X,Y)=\left(\frac{x-a(t)}{L(t)},\frac{y}{L(t)}\right),
\tag{D6}
\]

并取

\[
\rho(t,x,y)=\rho_c+R(z),
\qquad
\Psi(t,x,y)=L(t)P(z).
\tag{D7}
\]

则完整物理 IPM 严格等价于

\[
\boxed{
(U+cz-be_1)\cdot\nabla R=0,
\qquad
-\Delta P=R_X,
\qquad
U=\nabla^\perp P,
\qquad
P(X,0)=0 .}
\tag{B0}
\]

若 (B0) 存在一个有界光滑解，满足

\[
R_X\in C_0(\overline{\mathbb H}),
\qquad
R_X(z_*)\neq0,
\tag{D8}
\]

则 (D7) 立即给出

\[
\rho_x(t,x,y)=\frac1{c(T-t)}R_X(z),
\qquad
\|\rho_x(t)\|_\infty
=\frac{\|R_X\|_\infty}{c(T-t)}.
\tag{D9}
\]

这会同时满足：有界密度、光滑初值、全方向 \(\rho_x\) 衰减和精确 Type-I
水平梯度爆破。它利用的是 \(R\) 在无穷远可趋于依赖方向的常数；因此不要求
\(R\in L^p\)，旧的两个 Casimir 排除不适用。

### 3.1 壁面内核不是自由参数

写壁面展开

\[
R=r_0(X)+r_1(X)Y+\cdots,
\qquad
P=p_1(X)Y+\frac12p_2(X)Y^2+\cdots.
\tag{D10}
\]

在 \(r_0'(X)\neq0\) 的 wall transition 上，(B0) 逐阶强制

\[
p_1=cX-b,
\qquad
p_2=-r_0',
\qquad
r_1=-\frac{(r_0')^2}{2c}.
\tag{D11}
\]

若把爆破点和平移 gauge 置于零，并记 \(r_0'(0)=\alpha\)，则

\[
R=\alpha X-\frac{\alpha^2}{2c}Y+O(|z|^2),
\qquad
P=cXY-\frac\alpha2Y^2+O(|z|^3).
\tag{D12}
\]

这正是此前 affine 锐扇的局部 Taylor 法形。它现在只充当内核，不允许把
\(XY\) 应变延伸到无穷远。另一个直接后果是：若 wall trace 的导数在任意大
\(|X|\) 仍非零，则 (D11) 强迫 \(P_Y=cX-b\) 一直增长，与下面的屏蔽外尾冲突。
所以合格 wall front 应在远端真正变平；平坦紧支撑边缘而不是 tanh 尾最自然。

### 3.2 一个显式相容的方向常数外尾

设 \(b=0\)，并在 \(r\to\infty\) 取

\[
R(r,\theta)=F_0(\theta)+r^{-1}F_1(\theta)+\cdots,
\qquad
P(r,\theta)=rG_0(\theta)+G_1(\theta)+\cdots.
\tag{D13}
\]

Poisson 方程首阶为

\[
G_0''+G_0=\sin\theta\,F_0',
\qquad
G_0(0)=G_0(\pi)=0.
\tag{D14}
\]

避免 \(r\log r\sin\theta\) 共振的 Fredholm 条件是

\[
\int_0^\pi \sin^2\theta\,F_0'(\theta)\,d\theta=0.
\tag{D15}
\]

一个完全显式的非平凡首项是

\[
\boxed{
F_0=A\left(\cos\theta+\frac53\cos3\theta\right),
\qquad
G_0=-\frac{4A}{3}\sin^4\theta .}
\tag{D16}
\]

它逐项满足 (D14)--(D15)，并且

\[
F_0(0)=\frac{8A}{3},
\qquad
F_0(\pi)=-\frac{8A}{3},
\qquad
R_X=O(r^{-1}).
\tag{D17}
\]

更一般地，leading far end 可以任意集中到给定的右侧锐角扇。任取

\[
0<\theta_-<\theta_+<\frac\pi2,
\qquad
0\not\equiv G_0\in C_c^\infty((\theta_-,\theta_+)),
\]

并定义

\[
\boxed{F_0'(\theta)=\frac{G_0''(\theta)+G_0(\theta)}{\sin\theta}.}
\tag{D17a}
\]

则 (D14) 按定义严格成立，而且

\[
\int_0^\pi\sin^2\theta F_0'\,d\theta
=\int_0^\pi\sin\theta(G_0''+G_0)\,d\theta=0
\tag{D17b}
\]

由两次分部积分自动成立。因此

\[
R_X=-\frac{\sin\theta}{r}F_0'(\theta)+O(r^{-2})
\tag{D17c}
\]

只在 \((\theta_-,\theta_+)\) 内 active，并沿所有方向趋零。若还要求扇形两侧
的 \(F_0\) 取同一个常数，只需再加一个线性矩条件
\(\int_{\theta_-}^{\theta_+}F_0'\,d\theta=0\)。事实上两次分部积分还给出

\[
\int F_0'\,d\theta
=\int\frac{G_0''+G_0}{\sin\theta}\,d\theta
=2\int G_0(\theta)\csc^3\theta\,d\theta.
\tag{D17d}
\]

因此取两个不交的标准 compact bumps \(\phi_1,\phi_2\)，令
\[
G_0=\phi_1-
\frac{\int\phi_1\csc^3\theta\,d\theta}
{\int\phi_2\csc^3\theta\,d\theta}\phi_2,
\]
即可显式满足该条件。由于 (D17b) 的权严格为正，任何非零 \(F_0'\) 必须变号，所以
right-fan localization 与 signed angular return 是同一个结构，而不是缺陷。
这已经证明“右侧锐角集中 + 全方向 \(R_X\to0\)”在 leading inner/outer data
层面完全相容。

`analysis/ipm_plot_compact_acute_tail.m` 对一个
\(8^\circ<\theta<38^\circ\) 的 two-bump 例同时作图和审计：leading
Poisson、Fredholm、same-constant return 残差分别为
\(2.22\times10^{-16}\)、\(1.50\times10^{-18}\)、
\(5.80\times10^{-18}\)，扇外 activity 为机器意义的严格零。该脚本只验证
asymptotic data，不声称已得到 nonlinear fixed profile。

第一次输运修正为

\[
F_1=\frac1cG_0F_0',
\tag{D18}
\]

而 \(G_1\) 由 Dirichlet ODE

\[
G_1''=\cos\theta F_1+\sin\theta F_1'
=(\sin\theta F_1)'
\tag{D19}
\]

唯一确定。更高阶形成递归的无限角模级联。当前尚未证明该级数收敛，也尚未
证明它能与 (D12) 连成全半平面异宿；但 (D16) 已经证明新远端条件与椭圆方程、
壁面条件及第一输运修正之间没有 leading-order 矛盾。

若要求 \(F_0\) 严格单调，则 (D15) 不可能成立；这时必须接受
\(r\log r\sin\theta\) 的流函数共振，或引入 signed angular return。因而
“\(\rho_x\) 全局单符号、远场速度不过度增长、左右常数不同”三者一般不能同时
保留。

同样，若 \(P\) 在无穷远留下 \(B_\infty XY\)，leading transport 给出

\[
B_\infty\sin(2\theta)F_0'(\theta)=0.
\tag{D20}
\]

对非恒定 \(F_0\) 只能有 \(B_\infty=0\)。所以 affine strain 必须被外层屏蔽；
它只能是内端 Hessian，不能再作为全局 harmonic contract。

### 3.3 壁面必须是 compact-active，特征必须弯曲

壁面条件给出一个比“渐近衰减”更强的逐点恒等式：

\[
\boxed{\bigl(cX-b-P_Y(X,0)\bigr)R_X(X,0)=0.}
\tag{D20a}
\]

若外端为 \(P=rG_0(\theta)+o_{C^1}(r)\)，则壁面上 \(P_Y=O(1)\)，而
\(cX-b\) 线性增长。因此当 \(|X|\) 充分大时，必须有

\[
R_X(X,0)\equiv0,
\tag{D20b}
\]

不是仅仅趋于零。非平凡桥接若存在，只能有一个端点无限阶平坦的
\(C^\infty\) compact active wall segment；实解析 wall trace 会被解析延拓
直接排除。

令 \(a(X)=R_X(X,0)\)，在 active segment 上写

\[
P=cXY+\sum_{n\ge2}p_n(X)Y^n,
\qquad R=\sum_{m\ge0}r_m(X)Y^m.
\]

Poisson 与输运逐阶强制

\[
p_2=-\frac a2,\quad r_1=-\frac{a^2}{2c},\quad
p_3=\frac{aa'}{6c},\quad
r_2=\frac{5a^2a'}{16c^2}.
\tag{D20c}
\]

所以只要 \(a\) 不是常数，横向层级立即变成无限维；affine 解只是常斜率时的
偶然截断。用壁面足点 \(q\) 标记特征，令 \(X=X(q,Y)\)、\(R=F(q)\)，输运
还可以精确积分成

\[
\boxed{P(q,Y)=c\left(2\int_0^YX(q,s)\,ds-YX(q,Y)\right)-bY,}
\tag{D20d}
\]

记
\(J=X_q>0\)、\(K=X_Y\)、\(I=\int_0^YX(q,s)\,ds\)、
\(\mathcal Q=2I_q-YJ\)、\(H=X-YK\)，则剩余 Poisson 方程精确成为

\[
-\partial_q\!\left[c\frac{1+K^2}{J}\mathcal Q-K(cH-b)\right]
-\partial_Y\!\left[-cK\mathcal Q+J(cH-b)\right]=F'(q).
\tag{D20d1}
\]

其有效特征场为

\[
W=U+cz-be_1=\frac{2cI_q}{J}(K,1).
\tag{D20d2}
\]

所以 \(J>0\) 时 active 特征向上开放。壁面展开给出

\[
X_Y(q,0)=\frac{F'(q)}{2c},
\qquad
X_{YY}(q,0)=\frac{F'(q)F''(q)}{8c^2}.
\tag{D20e}
\]

事实上 \(F'(q)=2cX_Y(q,0)\) 是 hodograph PDE 的自然输出，不能再作为独立
Neumann 数据。若以 \(X\) 和顶端射线斜率 \(\xi(q)\) 为未知，interior PDE、
bottom value、两侧 exterior DtN 以及 top 的 value/derivative polyhomogeneous
matching 的离散方程总数恰等于未知量数；正确的自由边界 Newton 系统并不过定。

因此 wall slope 一变化，特征离开壁面时就必须弯曲；直 ridge、平移 ridge 和
共同宽度 \(q=X/s(Y)\) 都不能闭合一个有界非恒定剖面。

### 3.4 对数修正不是可选项，而是第三阶共振强制项

纯幂尝试

\[
R\sim\sum_{n\ge0}r^{-n}F_n(\theta),\qquad
P\sim\sum_{n\ge0}r^{1-n}G_n(\theta)
\]

满足三角递推

\[
G_n''+(n-1)^2G_n
=n\cos\theta F_n+\sin\theta F_n',
\tag{D20f}
\]

而输运中的 \(-cnF_n\) 对每个 \(n\ge1\) 可逆。以 (D16) 为首项继续计算，
在 \(n=3\) 得到 Poisson 右端对 Dirichlet kernel \(\sin2\theta\) 的投影

\[
\boxed{[S_3]_{\sin2\theta}
=-\frac{25}{48}\frac{A^4}{c^3}\ne0.}
\tag{D20g}
\]

该系数与前一级自由 harmonic mode 无关，所以纯 Laurent 级数严格失败。利用

\[
-\Delta\bigl(r^{-2}\log r\,\sin2\theta\bigr)
=4r^{-4}\sin2\theta,
\]

唯一的最小修正是

\[
\boxed{P_{3,\log}=
\frac{25}{192}\frac{A^4}{c^3}
r^{-2}\log r\,\sin2\theta.}
\tag{D20h}
\]

它随后强制 \(R\) 出现 \(r^{-4}\log r\) 项。允许
\(r^{-n}(\log r)^j\) 后可以形式递推到任意有限阶，但尚未证明级数收敛或与
内核全局 gluing。指数小尾 \(e^{-\mu r^\alpha}\) 位于所有代数阶之外，不能
消除这个 \(r^{-4}\sin2\theta\) Fredholm 投影。此前询问的对数修正在这里
确实有用，而且是方程强制的。

这个共振不是显式三角例 (D16) 的偶然。对任意首项令 \(h=F_0'\)、\(b=0\)，
纯幂输运递推的前三阶为

\[
F_1=\frac{G_0h}{c},\qquad
F_2=\frac{(G_0F_1)'}{2c},
\]

\[
G_1''=(\sin\theta F_1)',\qquad
G_2''+G_2=2\cos\theta F_2+\sin\theta F_2',
\]

\[
F_3=\frac{2G_0'F_2+G_0F_2'+G_1'F_1-G_2h}{3c}.
\tag{D20h1}
\]

若 \(S_n=n\cos\theta F_n+\sin\theta F_n'\)，则严格分部积分给出

\[
\boxed{
\langle S_n,\sin((n-1)\theta)\rangle
=(n-1)\int_0^\pi\sin((n-2)\theta)F_n(\theta)\,d\theta.}
\tag{D20h2}
\]

所以 \(n=2\) 的 kernel 条件自动成立，而第一个可能障碍是
\(2\int\sin\theta F_3\)。若它非零，必须加入

\[
\boxed{P_{3,\log}=-\frac1\pi
\left(\int_0^\pi\sin\theta F_3\,d\theta\right)
r^{-2}\log r\,\sin2\theta.}
\tag{D20h3}
\]

这个泛函不随 \(G_2\mapsto G_2+\lambda\sin\theta\) 改变，因为改变量正比于
leading Fredholm 矩 \(\int\sin^2\theta h=0\)。对 (D17a)--(D17d) 的一个
两 bump compact-right-fan 例，独立 \(2401\) 点 sine audit 给出上述 log 系数
约为 \(4.64\)，且随加密稳定。这是 compact acute 子族中共振非零的数值证据，
不是严格的全族非消失定理。

形式递推还有一个直接对应目标条件的结果：每个 \(F_n\) 及 density 的 log
系数仍由 \(h\) 或较低阶 compact \(F_m\) 相乘产生，因此角支撑仍留在原右锐扇；
但 \(G_n\) 由全局 Dirichlet 角 ODE 求出，共振
\(\sin((n-1)\theta)\log r\) 也会传播到所有角度。故应要求
\(\rho_x\) 锐扇集中，而不能同时要求 \(P,U\) 角紧支撑。

### 3.5 显式 Poisson seed 与第一轮数值否证

令 \(S=1+X^2+Y^2\)、\(K=-4A/3\)。下式给出一个有界光滑的显式 seed：

\[
P^\sharp=\frac{cXY+K(Y^2+Y^4)}{S^{3/2}},
\]

\[
R^\sharp=\frac{c(-X^2Y-Y^3-4Y)
+K(3XY^4+X^3Y^2-2X^5+4XY^2-4X^3-2X)}{S^{5/2}}.
\tag{D20i}
\]

直接微分可验 \(-\Delta P^\sharp=(R^\sharp)_X\)、\(P^\sharp|_{Y=0}=0\)、
\((R^\sharp)_X\to0\)，且其远端正是 (D16)。但它不满足 (D20a)：
\(P_Y^\sharp(X,0)=cX/(1+X^2)^{3/2}\)，所以它只是满足椭圆约束和内外渐近的
Newton seed，绝不是爆破解。

脚本 `analysis/ipm_search_beta0_profile.m` 对开放出流、显式 harmonic 外应变和
精确离散 Poisson 消元做了第一轮 Gauss--Newton 搜索。最高测试使用
\(129\times49\) 网格和 \(32\times14=448\) 个模态，得到

\[
\mathrm{Poisson}=1.85\times10^{-11},\quad
\mathrm{bulk/wall}=0.2104/0.1044,\quad
\max_{\rm collar}|R_X|=0.2424.
\tag{D20j}
\]

残差没有趋零，系数范数约为 \(1116\)，wall trace 出现明显高频振荡；外应变
扫描还把最优点压向退化极限。故这轮计算是清楚的 noncandidate，而不是数值
存在性证据。它支持下一轮改用 compact wall slope、hodograph 弯曲特征和显式
polyhomogeneous 对数外端，而不是继续增加有限直模。

### 3.6 标准 Green 自屏蔽：椭圆层面可构造，输运层面仍未闭合

在半平面 Dirichlet Green 归一化下，不允许另加自由的 \(BXY\) 时，局部
双曲应变不是参数，而必须满足精确逆矩条件

\[
\boxed{
c=P_{XY}(0,0)=\frac2\pi\operatorname{p.v.}
\int_{\mathbb H}\frac{XY}{(X^2+Y^2)^2}R_X\,dX\,dY.}
\tag{D20k}
\]

可分部积分时，右端也等于

\[
\frac2\pi\operatorname{p.v.}\int_{\mathbb H}
R\frac{Y(3X^2-Y^2)}{(X^2+Y^2)^3}\,dX\,dY.
\tag{D20l}
\]

若 bounded return 全部位于重整化距离 \(r\ge D\)，则

\[
|P_{XY}^{\rm far}(0,0)|\le\frac{4\|R\|_\infty}{\pi D};
\qquad
R\in L^1\Longrightarrow
|P_{XY}^{\rm far}(0,0)|\le\frac{2\|R\|_1}{\pi D^3}.
\tag{D20m}
\]

所以外逃 reservoir 不能单独产生 order-one 的 inner \(c\)；必须有 signed
density head 留在重整化 \(O(1)\)、物理 \(O(T-t)\) 的区域。

不过这个 signed head 在椭圆层面有完全显式的紧支撑构造。取

\[
Q=R_X=X\phi(r),\qquad \phi\in C_c^\infty((r_1,r_2)),
\qquad R(X,Y)=\int_{-\infty}^XQ(S,Y)\,dS.
\tag{D20n}
\]

则每条水平线上 \(\int Q\,dX=0\)，\(R,Q\) 光滑紧支撑，右侧 \(Q>0\)、
左侧为负 return，而且

\[
\boxed{P_{XY}(0,0)=\frac4{3\pi}\int_0^\infty\phi(r)\,dr.}
\tag{D20o}
\]

令 \(\int\phi=3\pi c/4\) 就精确选出 \(P=cXY+O(r^4)\)。因此 Poisson、
source-neutral 多极矩和所需正负号都不是障碍；正梯度还可用角权集中到右锐扇，
负回流放在左侧。但 (D20n) 只是一时刻的椭圆 gluing，不满足完整输运。

确实，固定 \(\beta=0\) 场
\(W=U+cz-be_1\) 有 \(\nabla\cdot W=2c\)。若一个非平凡 stationary profile
的每个非背景 level set 都有有限面积，取支撑于非背景值域的非负 \(h(R)\)，
则

\[
0=\int_{\mathbb H}\nabla\cdot(Wh(R))
=2c\int_{\mathbb H}h(R),
\tag{D20p}
\]

矛盾。故不存在非平凡 compact/finite-level-area 的 stationary
\(\beta=0\) profile；方向常数无限面积尾或非定常 return 至少要保留一个。

非定常路线的面积匹配反而是精确相容的。令 \(s=-\log(T-t)\)，则

\[
cR_s+(U+cz-be_1)\cdot\nabla R=0,
\qquad -\Delta P=R_X.
\tag{D20q}
\]

若物理初值的 \(\lambda\)-超水平集面积为 \(A_\lambda\)，而局部角尾
\(F_0(\theta)\) 在 \(r_z\sim e^s a(\theta)\) 处接 cap，则所有 level-area
Casimir 要求

\[
\boxed{A_\lambda=\frac{c^2}{2}
\int_{\{F_0(\theta)>\lambda\}}a(\theta)^2\,d\theta.}
\tag{D20r}
\]

这个 cap 的物理半径是 \(O(1)\)，对 inner \(P_{XY}\) 的影响只有
\(O(e^{-s})=O(T-t)\)；order-one 的重整化应变仍由局部 signed head 产生。
因此“内头产生爆破 + 外 cap 保存全部 Casimir”在尺度与椭圆矩上相容，剩下的
核心问题精确缩成非定常输运连接层。

## 4. 只要求 \(\rho_x\to0\) 后，哪些旧排除仍有效

仍然有效且不依赖有限能的障碍包括：

1. \(D=\rho+\mu\Psi=0\) 的纯 Bessel 比例闭包只能给静态解；参数调制必须由
   order-one \(D\) 驱动。
2. small-argument \(K_0/K_1\) 的局部单极/偶极矩不能由小 corrector 切掉；
   修正相对量仍至少为 \(1/[q\log(1/q)]\) 或 \(1/q\)。
3. 单 transverse mode、有限 Fourier/有限 ridge 和任意共同 wall envelope 的
   最高模障碍仍成立；真正二维闭合必然产生无限 cascade。
4. 纯齐次、单符号、无 sheet 的严格锐扇以及有限 \(Y\)-多项式 bridge 仍被局部
   代数排除。
5. 若 affine 共动密度被冻结且 self velocity 处处切于 density level sets，任何
   有界 regular superlevel component 都导致通量符号矛盾；二维 compact bump
   不能靠这种 tangent closure 被仿射压缩。
6. 固定点有效速度的散度为 \(2c>0\)，所以 return 不能封在闭 characteristic
   cell 中。

不再能直接使用的结论包括：

1. 依赖 \(R\in L^p\) 的 Casimir、有限势能和 source-neutrality 全局积分排除；
2. 用两个有限 \(L^p\) 范数强制面积保持 reservoir 的结论；
3. 在衰减、有限能算子域中得到的纯 Bessel 实点谱排除；
4. 依赖 \(L_x^2\) 的 Fourier multiplier 和有限能 Fredholm 选择。

换言之，外层现在可以由方向依赖常数或更慢的非 \(L^p\) 尾承担 reservoir；
但 local \(D\)-层、wall jet、Bessel 多极矩和无限模要求一个也没有消失。

## 5. 为什么一条“逐渐变宽的光滑条带”仍不自动闭合

一个自然尝试是

\[
q=\frac{X}{s(Y)},
\qquad
R=F(q),
\tag{D21}
\]

其中 \(s(Y)\to\infty\)，于是 \(R_X=F'(q)/s(Y)\) 看起来满足 (D1)。
输运方程可以精确积分出所有允许的流函数：

\[
P=cq\,[2S(Y)-Ys(Y)]+H(q),
\qquad
S'=s.
\tag{D22}
\]

把它代回 Poisson 方程后，\(P_{XX}\) 产生 \(H''/s^2\)，而其他曲率项含
\((s'/s)^2H''\) 和 \(((s'/s)^2-(s'/s)')H'\)。要与
\(F'(q)/s\) 对所有 \(Y\) 匹配，最自然的幂律展开中唯一非平凡尺度是
\(s(Y)\sim(Y+d)^2\)，但此时 \(s^{-2}\) 项强制 \(H''=0\)，继而

\[
F(q)=F_c-Cq^2.
\tag{D23}
\]

这就是一个抛物多项式核；它在 \(|q|\to\infty\) 无界，若硬截断又产生 edge
不匹配。因此共同单相位的弯曲/变宽条带仍不足以给出全局光滑 front；需要
额外角模或非切向 self transport。

## 6. \(0<\beta<1\) 分支仍未被排除，但不能是最终 hard fan

若不要求密度有界，只要求 (D1)，则

\[
R(r,\theta)\sim r^\beta F(\theta),
\qquad
0<\beta<1
\tag{D24}
\]

自动给出 \(R_X=O(r^{\beta-1})\to0\)。所以旧的 \(r^\beta\) 路线仍是合法
的慢尾 reservoir。不过纯齐次 hard fan 仍被无-sheet 与边缘匹配排除，必须带
低齐次 signed return、弯曲 edge 或无限角模修正。

对于 Type-I 内尺度

\[
\rho(t,x)=R\!\left(\frac{x}{T-t}\right),
\tag{D25}
\]

令 \(r=e^z\)、\(R=rF(z,\theta)\)、\(P=r^2G(z,\theta)\)，可得精确圆柱系统

\[
-\bigl(G_{zz}+4G_z+G_{\theta\theta}+4G\bigr)
=\cos\theta(F+F_z)-\sin\theta F_\theta,
\tag{D26}
\]

\[
(1-G_\theta)(F+F_z)+(2G+G_z)F_\theta=0.
\tag{D27}
\]

若 \(z\to+\infty\) 真趋于非平凡 1-homogeneous fan，则 \(R_X\) 趋于一个
非零角函数，而不是零。因此该 fan 最多只能是中间渐近；最外面仍需
\(\beta=0\) 方向常数层或随重整化时间外逃的 return。

## 7. 开放 harmonic contract 下的精确二维闭合

令 \(A(t)\in SL(2)\)、\(\zeta=A(t)(x-a(t))\)，并让 harmonic affine flow
满足 \(A_t+AM=0\)。写

\[
\rho=\rho_c+F(t,\zeta),
\qquad
\Psi=\Psi_H+H(t,\zeta).
\tag{D28}
\]

则完整 IPM 精确化为

\[
F_t+\{H,F\}_\zeta=0,
\qquad
-\nabla_\zeta\cdot(AA^T\nabla_\zeta H)
=(Ae_1)\cdot\nabla_\zeta F.
\tag{D29}
\]

这允许真正二维的 \(F\)，并自动包含端帽和无限 transverse cascade。对
\(A=\operatorname{diag}(\tau^{-1},\tau)\)、\(\tau=T-t\)、
\(\sigma=\tau^2/2\)、\(H=\tau\Phi\)，它精确变成

\[
F_\sigma=\{\Phi,F\},
\qquad
-(\Phi_{XX}+4\sigma^2\Phi_{YY})=F_X,
\qquad
\rho_x=\tau^{-1}F_X.
\tag{D30}
\]

对 \(Y\) 作 Fourier 变换，椭圆逆算子有显式核

\[
\widehat\Phi(X,\eta)
=-\frac12\int_{\mathbb R}\operatorname{sgn}(X-Z)
e^{-2\sigma|\eta||X-Z|}\widehat F(Z,\eta)\,dZ.
\tag{D30a}
\]

这不是形式近似，而是 sign-antiderivative 与 \(Y\)-Poisson 半群的精确组合。
因此对所有 \(\sigma\ge0\) 有一致核估计

\[
\|\partial_X^p\partial_Y^q\Phi\|_\infty
\le \frac12
\|\partial_X^p\partial_Y^qF\|_{L_X^1L_Y^\infty}.
\tag{D30b}
\]

在同时给 \(X,Y\) 正解析半径的 Banach scale 中，把衰减权重放在被输运的
\(U=\nabla F\) 上，并用 \(\mathcal K_\sigma:L_X^1\to L_X^\infty\) 控制无权
速度；不能把权重错误地放到具有 Poisson 远尾的 \(\mathcal K_\sigma U\) 上。
此时

\[
U_\sigma=\nabla\bigl[(J\mathcal K_\sigma U)\cdot U\bigr]
\]

只损失一个解析半径。具体可取

\[
\|U\|_r=\sum_\alpha\frac{r^{|\alpha|}}{\alpha!}
\left(\|\partial^\alpha U\|_\infty
+\|\partial^\alpha U\|_{L_X^1L_Y^\infty}\right),
\tag{D30b'}
\]

则对 \(0<r'<r\)，一致于 \(\sigma\ge0\)，

\[
\|\mathcal N_\sigma(U)-\mathcal N_\sigma(V)\|_{r'}
\le\frac{C}{r-r'}(\|U\|_r+\|V\|_r)\|U-V\|_r.
\tag{D30b''}
\]

标准 Ovsyannikov/Cauchy--Kowalevski 迭代因而从 \(\sigma=0\) 给出一个局部
唯一解析衰减解。若各阶导数同时属于 \(C_0\)，非线性每项含一个衰减的
\(U\) 因子，故该性质保持。Gaussian 情形的 profile 速度及所有导数有界，
特征流与恒等映射相差有界且各阶导数有界，所以 Schwartz 扰动由与逆流复合而
保持 Schwartz。

### 7.1 满足 \(\rho_x\in C_0\) 的光滑 Type-I 精确轨道

在全平面取终端剖面

\[
F_*(X,Y)=\varepsilon e^{-X^2-Y^2}.
\tag{D30c}
\]

取 \(\rho_c=0\)。上述解析局部存在性给出
\(0\le\sigma\le\sigma_0\) 的精确解。把它按 (D28)--
(D30) 拉回物理坐标，并在 \(t=0\) 取 \(T=\sqrt{2\sigma_0}\)，得到解析 Schwartz
初值和无源 IPM 的精确轨道。对每个 \(t<T\)，

\[
\rho_x(t,\cdot)\in C_0(\mathbb R^2),
\qquad
\tau\|\rho_x(t)\|_\infty
=\|F_X(\sigma)\|_\infty
\longrightarrow |\varepsilon|\sqrt{\frac2e}>0.
\tag{D30d}
\]

最大点满足
\((x,y)=(\pm\tau/\sqrt2+o(\tau),O(\tau))\)，所以峰不向空间无穷远逃逸。
但这不是孤立点集中：对任意固定有限 \(y_0\)，取
\(x=\pm\tau/\sqrt2+o(\tau)\) 时更精确地有
\[
\tau\rho_x(t,\pm\tau/\sqrt2,y_0)
\longrightarrow\mp\varepsilon\sqrt{2/e}.
\tag{D30d'}
\]
物理解在
\(x\) 方向压成宽 \(O(\tau)\) 的层，同时在 \(y\) 方向展宽到
\(O(\tau^{-1})\)；极限奇异集是竖直 sheet，而不是右侧锐角点。这仍然是按
(D1)--(D2) 字面意义的光滑 Type-I 轨道，但不能代替最终锐角目标。

半平面 Dirichlet 可取 odd terminal profile

\[
F_*(X,Y)=\varepsilon e^{-X^2}Y e^{-Y^2},
\]

奇延拓保持 \(\Phi|_{Y=0}=0\)，并有
\(\tau\|\rho_x\|_\infty\to|\varepsilon|/e\)。但是其峰位
\(y_*\sim(\sqrt2\tau)^{-1}\to\infty\)，所以它只是半平面全局范数爆破，不是
固定壁点奇点。若 wall trace 随 \(X\) 非常数，Dirichlet Green 核会产生
\(\Phi_Y|_{Y=0}\sim\tau^{-2}\mathcal H_XF\) 的新快层，不能套用 odd 定理。

### 7.2 这个正例的合同边界

(D30) 在 Fourier \(L^2\) 中仍有坏频带

\[
k_X\sim\tau^2k_Y,
\qquad
\sup\left|\tau\frac{k_Xk_Y}{k_X^2+\tau^4k_Y^2}\right|
=\frac1{2\tau}.
\tag{D30e}
\]

所以普通各向同性 Sobolev LWP 不能证明任意 \(C_c^\infty\) blob 延续到
\(\tau=0\)；上面的正例依赖解析局域化核估计。更重要的是
\(\Psi_H=xy/\tau\) 二次增长、其速度线性增长、无限能，且
\(\|\nabla u_H\|\sim\tau^{-1}\)。
它逐点没有体源，却要求预先给定
$M(t)=\operatorname{diag}(-1/\tau,1/\tau)$ 的时间奇异远场压力/速度合同；这个
未来奇异合同不由 $t=0$ 的光滑密度 Cauchy 数据决定。因此它只是开放 harmonic
far-field contract 下的精确 trajectory，不是“从光滑初值自动爆破”的标准
Cauchy 正例。标准问题仍需用外部 density return 自洽地产生同一个局部 affine
Hessian。

还有一个严格的非微扰结论。若去掉 \(\Psi_H\)，仍使用同一面积保持坐标和
\(\Psi_s=\tau\Phi\)，则方程变成

\[
2\sigma F_\sigma
=XF_X-YF_Y+2\sigma\{\Phi,F\}.
\tag{D30f}
\]

若 \(F\to F_*\) 于局部 \(C^1\)，且 self bracket 在 \(\sigma\downarrow0\) 保持
有界，则终端必须满足

\[
XF_{*X}-YF_{*Y}=0,
\qquad F_*=H(XY).
\tag{D30g}
\]

但 \(F_{*X}=YH'(XY)\in C_0\) 强制 \(H'\equiv0\)：固定任意
\(q=XY\)，令 \(Y\to\infty\)、\(X=q/Y\)，即可得到
\(YH'(q)\to0\)。所以非平凡的自洽替代绝不可能是上面 harmonic trajectory
的小扰动；它必须使 \(\{\Phi,F\}=O(\sigma^{-1})\)，即由 wall fast layer、
锐角内核或进入 leading order 的 density return 激活退化椭圆坏频带。

### 7.3 壁面 Hilbert 闭合强制动态 edge

同一面积保持坐标还能把“density 能否自生取代 \(\Psi_H\)”化成一个严格的
壁面必要条件。令 \(s=-\log\tau\)，去掉外加 harmonic strain，则完整方程是

\[
F_s+XF_X-YF_Y+\tau^2\{\Phi,F\}=0,
\qquad
-(\Phi_{XX}+\tau^4\Phi_{YY})=F_X.
\tag{D30h}
\]

半平面 Dirichlet--Neumann 数据的 Fourier 公式为

\[
\widehat{\Phi_Y}(k,0)=\frac{ik}{\tau^4}
\int_0^\infty e^{-|k|Y/\tau^2}\widehat F(k,Y)\,dY.
\tag{D30i}
\]

若 wall fast layer 中 \(F\to f(X)\) 足够光滑，则

\[
\tau^2\Phi_Y(X,0)\longrightarrow-\mathcal Hf(X),
\]

所以 stationary wall trace 在 active set 上必须满足

\[
\boxed{(X+\mathcal Hf)f'=0.}
\tag{D30j}
\]

若只有一个 active interval \((-a,a)\) 且两端接同一常数，有限 Hilbert 反演的
基本解是

\[
f(X)=C-\sqrt{a^2-X^2},\qquad |X|<a,
\tag{D30k}
\]

因为 \(\mathcal H[\sqrt{a^2-\cdot^2}_+]=X\)。但
\(f'=X/\sqrt{a^2-X^2}\) 在两端有平方根奇性。因此单区间、完全 stationary、
自屏蔽的 wall closure 不可能同时具有 \(C^\infty\) compact-active edge。
从光滑初值出发时，edge/return 必须随 \(s\) 运动并形成第二尺度。

局部符号本身是对的：

\[
\lim_{\tau\to0}\tau\Psi_{xy}(0,0)=-\Lambda f(0).
\tag{D30l}
\]

例如 \(f=-\varepsilon e^{-X^2}\) 给
\(-\Lambda f(0)=2\varepsilon/\sqrt\pi>0\)，能在中心产生正确的正双曲应变；
它只满足中心 Taylor matching，不满足整个 active interval 上
\(\mathcal Hf=-X\)。这把剩余结构定得很具体：不是再找一个光滑静态 wall
fixed point，而是找 signed head、动态平方根 edge regularization 和
\(r_z\sim e^s\) Casimir cap 的三层连接。

## 8. 当前最小的可判定任务

优先级已经改变：

1. **先解自洽的 \(\beta=0\) BVP。** 内端用 compact wall slope、
   (D20c)--(D20e) 和本征值条件 (D20k)，外端必须用含 (D20h3) 的
   polyhomogeneous 对数层级；寻找 genuinely 2D、curved-characteristic 的
   affine-to-angular-constant 异宿。
2. **改造数值 unknown，而不是盲目加模。** 以 hodograph \(X(q,Y)\) 或
   \(Q=R_X\) 为主变量，用 (D20i) 作 Poisson seed；把当前 448-mode
   noncandidate 当作回归基线，并检查 wall flatness、log coefficient 和最大残差。
3. **并行求非定常 return。** (D20p) 已排除 finite-level-area compact
   stationary profile，(D30k) 又排除单区间光滑 stationary wall edge；直接以
   \(Q=R_X\in C_0\) 为未知，保留方向常数背景
   \(B(Y)=\lim_{X\to-\infty}R\)，用 (D20r) 锁定 cap，演化完整的
   输运--Poisson 系统；检验 signed head 局部收敛、edge 动态平滑，而 return
   向 \(|z|\to\infty\) 外逃。
4. **Bessel 降为中间 Green corrector。** 它仍负责局部椭圆响应和 flat-edge
   匹配，但不再承担最终远尾，也不再假设 small-argument Bessel 主导内核。

截至目前最精确的判断是

\[
\boxed{
\begin{gathered}
\rho_x\to0\text{ 不排除爆破：锐角域有紧支撑 Lipschitz 顶点定理，}\\
\text{预设奇异 harmonic contract 下还有解析密度初态的 exact Type-I 轨道。}\\
\text{标准衰减 Biot--Savart 的光滑固定点问题仍未闭合；最短路线是}\\
\beta=0\text{ signed head 到强制对数锐角尾，再接动态 edge/Casimir cap。}
\end{gathered}}
\tag{D31}
\]
