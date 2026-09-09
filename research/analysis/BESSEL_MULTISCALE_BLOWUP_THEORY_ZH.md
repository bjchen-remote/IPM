# 利用精确 Bessel 稳态寻找多尺度 IPM 梯度爆破

日期：2026-08-28

中心点取在壁上 \((1,0)\)，并记

\[
X=x_1-1,\qquad Y=x_2,\qquad r=(X^2+Y^2)^{1/2}.
\]

本文专门回答一个问题：已经找到的穿孔域精确 Bessel 公式（等价地，带 tip
分布源的奇异稳态核），能否在理论上
导出一个真正的多尺度有限时爆破？

为避免把 formal matching 误当成存在性证明，下面始终区分三种状态：

- **精确/已证明**：由 IPM 方程、Green 恒等式或守恒律直接得到；
- **条件性**：若某个明确的重整化边值问题存在且其可解性系数为正，则缩放律严格跟随；
- **未证明**：目前还没有构造出完整的全局物理解。

## 0. 结论先行

结论不是“Bessel 稳态稍微收缩就会爆破”。目前已经能严格分成两部分：

> **精确动态内核已经找到：**在以移动参数吸收到固定物理中心后的坐标
> \((X,Y)\) 中，对任意 \(k>0\)，
> \[
> \Psi_{\rm af}=\frac{Y(X-kY)}{T-t},
> \qquad
> \rho_{\rm af}-\rho_c=\frac{2k(X-kY)}{T-t}
> \]
> 是上半平面中的完整 IPM 精确解。它满足不可穿壁、
> \(\partial_X\rho=2k/(T-t)>0\)，正相是右向锐角扇，顶点是双曲临界点；
> 奇延拓后还在 \(Y=0\) 有精确跳跃。
>
> **但这还不是有界数据爆破：**该解在空间无界、能量无穷，而且正相扇形
> 不是紧支撑扇形。它只能作为局部 inner normal form；必须在更外尺度改形并
> 连接保持全部 Casimir 的 return/reservoir。

> **严格排除：**固定 \(\mu\) 的奇异 Bessel 解不能作为
> “有界内核 + Bessel 主中间层”的小扰动基态。局域 Poisson 多极矩、
> 壁面条件以及所有 Casimir 守恒量共同排除了这一机制。
>
> **仍可利用：**Bessel 公式给出了动态内层方程的精确漂移 Helmholtz
> Green 核、右向 Gaussian 尾迹和可以投影的外层模态。把它固定为背景后，
> 总场可由一个 closure perturbation \(d\) 精确参数化：
> \[
> \Psi=\Psi_B+T_\mu d,\qquad
> \rho=\rho_B+(I-\mu T_\mu)d.
> \]
> 这给出真正的 inner--outer 椭圆桥，但 \(d\) 在内区必为 order one，不能是
> 全域小扰动。

> **后续精确突破：**把 order-one 余留真正保留后，可构造出
> 不在上述 small-residual no-go 假设内的完整轨道。其中
>
> \[
> \rho=\tanh\!\frac{x-a}{T-t},
> \qquad
> \Psi=\frac{(x-a)y}{T-t}
> -(T-t)\log\cosh\!\frac{x-a}{T-t}-\dot a\,y
> \]
>
> 给出有界密度的精确 \((T-t)^{-1}\) 梯度爆破；水平壁上还有
>
> \[
> P_*=cXY-bYe^{qX},
> \qquad
> R_*=(2-\beta)cYe^{qX},
> \qquad
> b=\frac{(2-\beta)c}{q}.
> \]
>
> 更强的是，后一 profile 族在右向单调的 \(Y\)-线性不变类与
> \(X\ge0\) characteristic sector 中有完全显式的非线性特征公式，
> 并严格局部收敛到由初始
> 壁面 jet 选出的 \(q<0\) directional Bessel 模。详细证明见
> [BESSEL_BRIDGE_EXACT_BLOWUP_ZH.md](BESSEL_BRIDGE_EXACT_BLOWUP_ZH.md)。
>
> 这些轨道需要不衰减、无穷能的 harmonic far-field contract；
> 水平壁 Bessel 族在某个远场方向还无界。因此它们已经证明
> “Bessel bridge + dynamic scaling”机制及其类内收敛，但仍不是标准
> source-neutral、衰减、有限能半平面 Cauchy 爆破定理。

还可把“光滑初值”说得更精确。对 \(0<A<1\)、\(\kappa>0\)，

\[
r_0(X)=\frac{(2-\beta)c}
{1-A+Ae^{(\kappa/A)X}}
\]

是全 \(X\) 上的 \(C_b^\infty\) 降维初值，它生成的精确特征解满足

\[
r(s,\cdot)\to(2-\beta)c e^{-\kappa X}
\quad\hbox{in }C^1_{\rm loc}(\mathbb R).
\]

但该轨道在有限重整化时间后从左无穷远失去一致有界性，而且
\(R=Yr\) 使物理初始密度随 \(Y\) 增长。因此这证明的是
“光滑但非衰减初值 \(\to\) 局部重整化 Bessel 收敛态”，而非标准
\(C_b\) 物理相空间中的吸引子。

更强的精确构造已经消除了“密度本身无界”这个缺陷。取 \(c>g>0\)，
并选 \(f\in C_c^\infty(\mathbb R)\)，使它在右支撑端点满足

\[
f(\xi)=e^{-1/(1-\xi)}\quad(\xi<1\hbox{ 且充分接近 }1),
\qquad f(\xi)=0\quad(\xi\ge1).
\]

令

\[
d=(gs)^{-1},\quad \varepsilon=1-d,\quad w=d^2,
\quad \delta=e^{-gs},\quad L=e^{-cs},\quad A=c-2/s,
\]

\[
R=\delta^{-1}f(\varepsilon+wX),
\qquad
p=\frac1{\delta w}\int_{\varepsilon+wX}^{1}f(\eta)\,d\eta,
\qquad
P=AXY+p.
\]

则 \(p_X=-R\)，而且逐项严格满足完整重整化 IPM

\[
R_s+(U+cz-ge_1)\cdot\nabla R=gR,
\qquad -\Delta P=R_X.
\]

在每个固定紧集上，

\[
(P,R)\longrightarrow(cXY+e^{-X},e^{-X})
\quad\hbox{in }C^\infty_{\rm loc}.
\]

物理密度为

\[
\rho-\rho_c=f\!\left(\varepsilon+\frac{x-a}{L_0}\right),
\qquad L_0=\frac Lw=g^2s^2e^{-cs},
\]

所以任意有限初始时刻的数据有界、光滑且在 \(x\) 方向紧支撑。这给出
“普通光滑物理密度 \(\to\) 局部 Bessel 收敛态”的严格肯定答案，并把
Bessel 核尺度 \(L\)、核到端点的 \(Lgs\) 尺度、以及 Casimir reservoir
尺度 \(L(gs)^2\) 明确分开。但密度仍不随 \(y\) 衰减，\(AXY\) 是随时间
增强的无穷能调和远场，且不满足水平不可穿壁；所以它仍不是由标准衰减
Biot--Savart 规范唯一决定的有限能 Cauchy 解。完整证明见 companion 报告
第 13.3 节。

若只问“能否由光滑初值产生 Bessel 含时极限”，还可以固定一个完全
时间无关的远场合同。取常应变
\(\Psi_\infty=\sigma(x-x_*)y\)，并让 \(F\in C_c^\infty\) 的右端具有
标准平坦边缘 \(F(\xi)=e^{1/\xi}\)，记
\(H(\xi)=\int_\xi^0F(\eta)\,d\eta\)。精确解

\[
\rho-\rho_c=F\!\left(\frac{x-x_*}{e^{-q_0-\sigma t}}\right),
\qquad
\Psi=\sigma(x-x_*)y
+e^{-q_0-\sigma t}
H\!\left(\frac{x-x_*}{e^{-q_0-\sigma t}}\right)
\]

由初始密度与该固定 contract 在 \(y\)-平移不变 classical 类中唯一决定。
在支撑尺度 \(\ell=e^{-q}\) 内再取
\(a=x_*-\ell/q\)、\(L=\ell/q^2\)、\(\delta=\ell\)，则

\[
(P,R)\longrightarrow(e^{-X},e^{-X})
\quad\hbox{in }C^\infty_{\rm loc}.
\]

这是没有 order-one harmonic residual 的纯 directional Bessel tangent。
它发生在 \(t\to\infty\)，不是有限时爆破，也尚未证明对邻近初值稳定。

横向局域化也有一个精确中间结果：把第 13.3 节的 compact front 写成
\(\xi=X-k(s)Y\)、\(k=k_0e^{-2cs}\)，并配上显式旋转调和二次场，可得到
有界、每条固定 \(x\) 或固定 \(y\) 切片都紧支撑的有限时轨道，仍严格局部
收敛到 \((cXY+e^{-X},e^{-X})\)。但总支撑是无限斜条带，二维质量和能量
仍无穷。相反，固定 transverse 单模、共同 \(y\)-cutoff 及任意有限
Fourier 闭合都可严格排除；有限能版本必须使用非可分、无限模、
order-one 二维 return。证明见 companion 报告第 14--15 节。

本文不再指定任何数值 \(\beta\)。目前没有可靠数据或 Fredholm 条件选择它；
\(0<\beta<1\) 必须作为全局匹配本征量求出。精确 affine 解本身完全不选择
\(\beta\)。

若要让 Bessel 仍参与组织外层，明确未被排除的路线是

\[
\text{主导 dynamic inner/}D\text{ 层}
\longleftrightarrow
\text{固定或调制的 Bessel--Gaussian 背景响应}
\longleftrightarrow
\text{全局 return/reservoir 层}.
\]

若还要在 \(Y=0\) 保留密度跳跃与非零 \(\partial_X\rho\)，必须再加一个
leading wall layer。

另外还有第 9 节的固定 \(\mu\)、\(D\)-主导局域 Type-II 缩放；它不是
Bessel-主导机制，并且只有在另加非局域等测度储层后才未被排除。

因此，Bessel 公式可以“用来找爆破”。对标准 localized \(K_0/K_1\)
内外层，它的正确用法仍是 Green 核和外层参数响应；但对与双曲
应变特征方向对齐的 directional mode，Bessel 成分本身可在 order-one
harmonic residual 帮助下闭合，并产生上述局部吸引奇化。两种陈述的
适用假设不同，互不矛盾。

## 1. 基本方程与精确 Bessel 族

物理 IPM 方程为

\[
\rho_t+U\cdot\nabla\rho=0,
\qquad
-\Delta\Psi=\partial_X\rho,
\qquad
U=\nabla^\perp\Psi=(-\Psi_Y,\Psi_X).
\tag{1.1}
\]

若令 \(\rho=-\mu\Psi\)，\(\mu>0\)，则非线性输运项恒等地消失，而 Poisson 方程化为

\[
L_\mu\Psi:=\bigl(\Delta-\mu\partial_X\bigr)\Psi=0.
\tag{1.2}
\]

写成 \(\Psi=e^{\mu X/2}\Phi\)，则

\[
L_\mu\Psi=e^{\mu X/2}
\left(\Delta-\frac{\mu^2}{4}\right)\Phi.
\]

因此有右向衰减的精确成员

\[
\Psi_0=-Ae^{\mu X/2}K_0\!\left(\frac{\mu r}{2}\right),
\qquad
\rho_0=\mu Ae^{\mu X/2}K_0\!\left(\frac{\mu r}{2}\right),
\tag{1.3}
\]

以及满足不可穿壁流函数条件的法向偶极支

\[
\Psi_1=-Ae^{\mu X/2}K_1\!\left(\frac{\mu r}{2}\right)\sin\theta,
\qquad
\rho_1=\mu Ae^{\mu X/2}K_1\!\left(\frac{\mu r}{2}\right)\sin\theta.
\tag{1.4}
\]

两者在去掉顶点的域内逐点满足完整物理稳态方程，不是拟合式。

但分布意义下它们携带不可忽略的 tip 源：

\[
L_\mu\Psi_0=2\pi A\delta_0,
\tag{1.5}
\]

\[
\Psi_1=-\frac2\mu\partial_Y\Psi_0,
\qquad
L_\mu\Psi_1=-\frac{4\pi A}{\mu}\partial_Y\delta_0.
\tag{1.6}
\]

\(K_0\) 是单极/点源；\(K_1\) 是法向偶极/壁上偶极源。

## 2. 闭包缺陷 \(D\)：从稳态到动态的精确变量

对任意状态定义

\[
\boxed{D:=\rho+\mu\Psi.}
\tag{2.1}
\]

当 \(\mu\) 对空间为常数时，(1.1) 精确等价于

\[
\boxed{
L_\mu\Psi=-\partial_XD,
\qquad
\rho_t+U\cdot\nabla D=0.
}
\tag{2.2}
\]

第二式只用了 \(U\cdot\nabla\Psi=0\)。它给出三个精确结论：

1. Bessel 闭包正是 \(D=0\)。
2. \(D=0\) 时 \(\rho_t=0\)，所以不能把 \(A,\mu\) 或中心改成时间函数后仍称为物理解。
3. 所有真实尺度运动都由 \(D\neq0\) 驱动。

分布恒等式 (1.5)--(1.6) 还使这一点更强。在选取
\(D(-\infty,Y)=0\) 的规范后，无源正则化至少需要

\[
D_0\simeq-2\pi A H(X)\delta(Y),
\tag{2.3}
\]

或

\[
D_1\simeq\frac{4\pi A}{\mu}H(X)\partial_Y\delta(Y).
\tag{2.4}
\]

这是沿右射线延伸的 closure layer，而不是只在顶点附近的小残差。

## 3. 第一个严格 no-go：不可压缩输运保持整个密度分布

只要解保持经典且边界无通量，对所有合适的 \(F\)，

\[
\frac{d}{dt}\int F(\rho(t,x))\,dx=0.
\tag{3.1}
\]

因此不仅 \(L^\infty\) 守恒，所有允许的 \(L^p\) 以及超水平集面积

\[
\bigl|\{x:\rho(t,x)>s\}\bigr|
\]

都守恒。这会排除一个孤立、密度差为 \(O(1)\)、在两个方向上都收缩的
localized core。若它真的出现，必须同时有：

- 连到远场的 expanding filament/reservoir，用来保持各个 level set 的面积；或
- vanishing-contrast 机制

\[
\delta(t)\to0,
\qquad
\frac{\delta(t)}{L_{\min}(t)}\to\infty,
\tag{3.2}
\]

使密度本身越来越平，但梯度越来越大。

还可以对任意候选的全局局域化 self-similar profile 做一行排除。若

\[
(U+c z)\cdot\nabla R=\beta cR,
\qquad c>0,\quad\beta>0,
\tag{3.3}
\]

且 \(R\in L^p\)且无无穷远通量，乘以 \(|R|^{p-2}R\) 后积分得

\[
-\frac{2c}{p}\int|R|^p
=\beta c\int|R|^p,
\]

所以 \(R=0\)。这说明非平凡收缩固定点必须向无穷远输出 Casimir 通量，
或者与另一个外尺度储层匹配。

## 4. 第二个严格 no-go：局域 Poisson 多极矩

### 4.1 \(K_0\) 单极矩

这一整圆论证用于内点或全平面延拓。对当前位于物理壁上的顶点，\(K_0\)
本来已被不可穿壁条件排除；下面的矩结论是对其全平面/穿壁基准的审计。
对包含顶点的圆 \(B_r\)，积分 Poisson 方程：

\[
-\int_{\partial B_r}\partial_n\Psi\,ds
=\int_{\partial B_r}\rho n_X\,ds.
\tag{4.1}
\]

常数密度背景在右端自动消去。记

\[
\delta_r=\operatorname{osc}_{\partial B_r}\rho.
\]

由 \(\Psi_0\sim A\log r\) 可得

\[
2\pi|A|+o(|A|)\lesssim r\delta_r.
\tag{4.2}
\]

### 4.2 \(K_1\) 法向偶极矩

对上半圆 \(B_r^+\)，以 \(h=Y\) 作 Green 恒等式：

\[
\int_{\partial B_r^+}
\bigl(Y\partial_n\Psi-\Psi n_Y\bigr)\,ds
=-\int_{\partial B_r^+}Y\rho n_X\,ds.
\tag{4.3}
\]

\(K_1\) 的小参数形式为

\[
\Psi_1\sim-\frac{2A}{\mu}\frac{Y}{r^2}.
\]

因此 (4.3) 的左端为 \(2\pi A/\mu+o(A/\mu)\)，并强制

\[
\frac{|A|}{\mu}\lesssim r^2\delta_r.
\tag{4.4}
\]

### 4.3 对原三层设想的结论

下面的“Bessel 主导”有一个明确的匹配范数含义：总流函数在匹配圆上保留
Bessel 的 leading normal-flux/法向偶极泛函，而总密度与修正则用该圆上的角向
oscillation 测量。若连这个 leading flux 都被修正消掉，那么 Bessel 在流函数中也已经
不再主导。

令 \(q=\mu r\ll1\)。在小 Bessel 参数区，

\[
\delta_{B,0}
\sim\mu|A|\log\frac1q
=\frac{|A|}{r}\,q\log\frac1q,
\tag{4.5}
\]

\[
\delta_{B,1}\sim\frac{|A|}{r}.
\tag{4.6}
\]

但 (4.2) 要求 \(K_0\) 修正至少为 \(|A|/r\)，即比 Bessel 分量大

\[
\frac1{q\log(1/q)}.
\]

(4.4) 要求 \(K_1\) 非对称修正至少为

\[
\frac{|A|}{\mu r^2}
=\frac{\delta_{B,1}}q.
\]

所以已经可以严格说：

\[
\boxed{
\ell\ll r\ll\mu^{-1}
\quad\text{中不存在“小修正 + 主导 small-argument Bessel”的匹配。}
}
\tag{4.7}
\]

因此先前只依据密度最大值推出的

\[
A\sim\frac1{\mu\log(1/\ell)}
\quad(K_0),
\qquad
A\sim\ell
\quad(K_1)
\]

不是自洽的 source-free gluing law。

## 5. 第三个严格 no-go：壁面跳跃不在 Bessel 闭包中

不可穿壁条件为

\[
U_2=\Psi_X=0\quad\text{on }Y=0,
\]

所以 \(\Psi\) 在壁上对 \(X\) 为常数。若再令 \(D=0\)，则

\[
\rho|_{Y=0}=-\mu\Psi|_{Y=0}
\]

也是 \(X\)-常数，因而

\[
\partial_X\rho|_{Y=0}=0.
\tag{5.1}
\]

具体地：

- \(K_0\) 有非零壁面密度迹，但 \(\Psi_X\neq0\)，它穿过壁面；
- \(K_1\) 满足 \(\Psi=0\) 与不可穿壁，但 \(\rho=\partial_X\rho=0\) 于壁上。

所以用户观察到的 \(Y=0\) 单侧非零 trace/奇延拓跳跃，以及壁上有统一符号的
\(\partial_X\rho\)，必须由 \(D\) 的 leading wall layer 承担。这不能是高阶小修正。

## 6. 全局源、Gaussian 尾迹和补偿层

\(K_0\) 在下游的上半平面横截质量满足

\[
\int_0^\infty\rho_0(X,Y)\,dY\longrightarrow\pi A.
\tag{6.1}
\]

因此它的总质量随下游长度线性发散。\(K_1\) 虽然横截质量衰减，但

\[
\int_0^\infty\rho_1(X,Y)\,dY
\sim2A\sqrt{\frac{\pi}{\mu X}},
\tag{6.2}
\]

所以总质量仍按 \(\sqrt X\) 发散。它的法向一阶矩则保持

\[
\int_0^\infty Y\rho_1(X,Y)\,dY
\longrightarrow\frac{2\pi A}{\mu},
\tag{6.3}
\]

这正是 tip dipole 的远场印记。

对有限势能、足够衰减、不可穿壁的光滑解，还有精确恒等式

\[
\frac{d}{dt}\int_{Y>0}Y\rho\,dX\,dY
=-\int_{Y>0}|U|^2\,dX\,dY.
\tag{6.4}
\]

因此任何同时满足这些有限性假设的物理稳态都必须 \(U=0\)，从而是分层的；
若再要求水平衰减，就只剩平凡解。

Bessel 解之所以不与 (6.4) 矛盾，正是因为它们有 tip 奇性、非可积尾，并且
\(K_0\) 还有壁面通量。

一个最小的 source-neutral \(K_1\) 外层模板是放置一对反号偶极：

\[
\Psi_{\rm pair}(z)
=\Psi_{1,A}(z)-\Psi_{1,A}(z-Se_1),
\tag{6.5}
\]

\[
D_{\rm strip}
=\frac{4\pi A}{\mu}
\bigl[H(X)-H(X-S)\bigr]\partial_Y\delta(Y).
\tag{6.6}
\]

它们逐分布满足

\[
L_\mu\Psi_{\rm pair}=-\partial_XD_{\rm strip}.
\tag{6.7}
\]

(6.5)--(6.7) 还不是光滑物理解，但它是正确的双端 source-neutral 外层骨架。
更精确地，它只闭合了椭圆/Poisson 方程：应同时定义

\[
\rho_{\rm pair}=-\mu\Psi_{\rm pair}+D_{\rm strip}.
\]

还没有 \(U_{\rm pair}\cdot\nabla D_{\rm strip}=0\) 的输运恒等式，因而绝不能把它称为稳态 IPM 解。
上述常数使用全平面奇延拓的 \(\delta\) 归一化；若直接把偶极分布放在半平面边界上，
系数会随单侧分布的规约而变。

## 7. Bessel 的正确角色：重整化内层的精确 Green 核

令动态 Bessel 长度为

\[
L(t)=\mu(t)^{-1}.
\]

在重整化变量 \(z=(x-a)/L\) 中，选取

\[
\rho-\rho_c=\delta(t)R(z),
\qquad
\Psi=\delta(t)L(t)P(z).
\tag{7.1}
\]

其中 \(\rho_c\) 取为被钉住物质点上的常数密度值。常数密度不进入 Poisson 方程，
也不进入输运梯度，因而令中心化闭包缺陷 \(\widetilde D=D-\rho_c\)。则

\[
\widetilde D
=\delta(t)Q(z),
\qquad
Q:=R+P.
\tag{7.2}
\]

重整化 Poisson 方程可写成两个完全等价的形式：

\[
-\Delta_zP=\partial_{z_1}R,
\tag{7.3}
\]

\[
(\Delta_z-\partial_{z_1})P=-\partial_{z_1}Q.
\tag{7.4}
\]

第二式中的基本解恰好就是同步得到的 Bessel 公式：

\[
G(z)=-\frac1{2\pi}e^{z_1/2}K_0\!\left(\frac{|z|}{2}\right),
\qquad
(\Delta-\partial_{z_1})G=\delta_0.
\tag{7.5}
\]

在上半平面取 Dirichlet 像法核

\[
G_H(z,z')=G(z-z')-G(z-z'^*),
\qquad z'^*=(z_1',-z_2'),
\tag{7.6}
\]

则

\[
P(z)=-\int_{\mathbb R_+^2}
\partial_{z_1}G_H(z,z')Q(z')\,dz'
+P_{\rm hom}(z).
\tag{7.7}
\]

公式 (7.7) 直接作为收敛卷积使用时，要求 \(Q\) 足够衰减/紧支撑，或者已在有限匹配边界截断并
保留边界项。第 8 节匹配 \(r^\beta\) 终端 cusp 时，应在有限匹配半径使用 (7.7)，
或将非衰减部分放入 \(P_{\rm hom}\) 与重整化卷积；不应把朴素全域积分无条件延伸到增长 profile。

\(K_0\) 是 (7.5) 的单极核，\(K_1\) 是它的法向导数/偶极核。这就是精确
Bessel 稳态在爆破问题中最坚固的理论用途：它精确表示了 closure defect \(Q\)
如何产生外层流函数。

## 8. 条件路线 A：收缩的内在 Bessel 尺度

### 8.1 精确重整化固定点方程

设终端密度在顶点附近有 \(0<\beta<1\) 阶 cusp，因而选

\[
\delta(t)=L(t)^\beta.
\tag{8.1}
\]

若尺度和中心满足

\[
\dot L=-c\delta,
\qquad
\dot a=\delta b e_1,
\qquad c>0,
\tag{8.2}
\]

则代入物理 IPM 后，静止重整化 profile 必须满足

\[
\boxed{
\begin{aligned}
&(U+c z-b e_1)\cdot\nabla R=\beta cR,\\
&-\Delta P=\partial_{z_1}R,\\
&U=\nabla^\perp P.
\end{aligned}}
\tag{8.3}
\]

在左右对称或已钉住中心的设定中 \(b=0\)；去掉水平对称时，\(b\) 是需由
平移零模决定的调制参数。若 \(a(t)\) 是壁上物质点，还必须有 pinning 兼容条件

\[
b=U_1(0)=-\partial_{z_2}P(0),
\tag{8.3b}
\]

此时 \(\rho_c\) 才是沿该物质轨道的常数。

要在 \(t\uparrow T\) 时匹配一个固定物理坐标中的 \(r^\beta\) 终端 cusp，内层 profile 必须在
进入下一外层之前满足

\[
R(z)\sim|z|^\beta F(\arg z),
\qquad |z|\to\infty,
\tag{8.3a}
\]

因为 \(L^\beta R((x-a)/L)\sim r^\beta F(\theta)\)。所以重整化内层本身不能全局衰减；
纯衰减 Bessel 项也不能单独提供 (8.3a)。必须由 \(Q\)、\(P_{\rm hom}\) 或更外层完成
cusp/return 匹配，再在更大的物理尺度上截断或连接等测度 reservoir。

使用 \(Q=R+P\)，(8.3) 还可写成

\[
(\Delta-\partial_{z_1})P=-\partial_{z_1}Q,
\tag{8.4}
\]

\[
-\beta cR+(c z-b e_1)\cdot\nabla R
+\nabla^\perp P\cdot\nabla Q=0.
\tag{8.5}
\]

这是一个明确、可检验、无经验拟合参数的非线性边值问题；
\(c,b\) 是由可解性条件决定的本征/调制参数。其半平面条件为

\[
P(z_1,0)=0,
\]

而不需要令 \(R(z_1,0)=0\)；实际上壁面跳跃正由

\[
Q(z_1,0)=R(z_1,0)
\]

承担。

### 8.2 一个完全显式的 affine 动态内核

不需要拟合或数值求解，(8.3) 实际上已有一个完全显式的局域固定点。令

\[
\xi=z_1-\frac bc,
\qquad
\eta=z_2,
\qquad
B=(1-\beta)c,
\qquad
k>0.
\tag{AF1}
\]

定义

\[
\boxed{
P_{\rm af}=B\eta(\xi-k\eta),
\qquad
R_{\rm af}=2Bk(\xi-k\eta).
}
\tag{AF2}
\]

则

\[
-\Delta P_{\rm af}=2Bk=\partial_\xi R_{\rm af},
\tag{AF3}
\]

而速度为

\[
U_{\rm af}
=\nabla^\perp P_{\rm af}
=B(-\xi+2k\eta,\eta).
\tag{AF4}
\]

直接乘开可得

\[
\bigl(U_{\rm af}+c(\xi,\eta)\bigr)\cdot\nabla R_{\rm af}
=\beta cR_{\rm af}.
\tag{AF5}
\]

因此 (AF2) 对每个 \(0<\beta<1\)、\(c>0\)、\(k>0\) 都是严格的重整化固定点。
它还不只是偶然猜中的多项式：它是任何足够光滑候选固定点在壁上顶点的强制
一阶 Taylor jet。把真实停滞点放在原点并归一化 \(R(0)=0\)，若

\[
R=\alpha\xi+d\eta+O(r^2),
\qquad
P=B\xi\eta-\frac\alpha2\eta^2+O(r^3),
\tag{AF5a}
\]

其中 \(P=0\) 于壁面，而 \(-\Delta P=R_\xi\) 已强制 \(-\alpha\eta^2/2\)
这一项，则把线性项代入输运方程得到

\[
B=(1-\beta)c,
\qquad
d=-\frac{\alpha^2}{2(1-\beta)c}.
\tag{AF5b}
\]

因此 \(\alpha=2Bk\) 时恰好恢复 (AF2)。若某个有界的 affine-to-cusp
heteroclinic 固定点存在，它在 \(z\to0\) 必须趋于这个双曲 jet；不是任意选择。
它有用户要求的三个几何特征：

\[
P_{\rm af}(\xi,0)=0,
\qquad
\partial_\xi R_{\rm af}=2Bk>0,
\tag{AF6}
\]

\[
R_{\rm af}>0
\quad\Longleftrightarrow\quad
\eta>0, \xi>k\eta.
\tag{AF7}
\]

正相因而是以顶点为中心、向右张开的锐角扇，其开角为

\[
\theta_*=\arctan\frac1k.
\tag{AF8}
\]

Taylor 斜率与开角还满足可检验的局部关系

\[
\tan\theta_{\rm core}
=\frac1k
=\frac{2(1-\beta)c}{\alpha}.
\tag{AF8a}
\]

所以在 \(\alpha>0\) 的右扇中，\(c>0\) 已是局部必要条件。若再把这个角与
第 8.8 节的真空外扇角 \(\theta_*=\beta\pi/(1+\beta)\) 对齐，则

\[
\frac\alpha c
=2(1-\beta)\cot\!\frac{\beta\pi}{1+\beta}.
\tag{AF8b}
\]

这只是自然的 gluing 归一化；affine 解在零线外仍有负相，而不是真空。

在顶点，\(U_{\rm af}(0)=0\)，而

\[
\nabla U_{\rm af}
=B\begin{pmatrix}-1&2k\\0&1\end{pmatrix}
\]

有特征值 \(\pm B\)。因此这个内核的顶点是一个真正的双曲型临界点，
与 \(K_0\) 尖点的旋转主部本质不同。

它还可以精确地奇延拓到全平面：

\[
\boxed{
P_{\rm odd}=B\bigl(\xi\eta-k\eta|\eta|\bigr),
\qquad
R_{\rm odd}=2Bk\bigl(\xi\operatorname{sgn}\eta-k\eta\bigr).
}
\tag{AF10}
\]

\(P_{\rm odd}\in C^1\)，且分布意义下

\[
-\Delta P_{\rm odd}
=2Bk\operatorname{sgn}\eta
=\partial_\xi R_{\rm odd},
\]

没有额外的壁上 Dirac 质量。同时

\[
[R_{\rm odd}]_{\eta=0}=4Bk\xi,
\qquad
(U_{\rm odd}+cz-be_1)_2|_{\eta=0}=0.
\tag{AF11}
\]

所以跳跃面是 characteristic 的，输运方程也不产生额外的界面 Dirac。
这是一个真正精确的带壁面跳跃弱 fixed point。

但“正相是扇形”不等于“profile 支撑在扇形内”。若将 (AF2) 零延拓为

\[
P_+=B\eta(\xi-k\eta)_+,
\qquad
R_+=2Bk(\xi-k\eta)_+,
\]

则

\[
-\Delta P_+
=2Bk\mathbf1_{\{\xi>k\eta\}}
-B(1+k^2)\eta\,\delta_{\xi=k\eta},
\tag{AF12}
\]

而 \(\partial_\xi R_+\) 没有第二项。故该零延拓不是全半平面 IPM 解。
另一方面，若把 \(0<\eta<\xi/k\) 本身当作一个独立扇形域，则两条边上
\(P=0\)，且 (AF2) 在域内是完全精确的。

把 (AF2) 反变换到物理变量。这里先把 \(b/c\) 的位移吸收到物理中心：若原来的
中心参数为 \(a(t)\)，则改用

\[
a_*(t)=a(t)+\frac bcL(t).
\]

由 \(\dot a=bL^\beta\)、\(\dot L=-cL^\beta\) 可见 \(\dot a_*=0\)，所以
\(a_*\) 正是固定的物理顶点；以下仍把相对它的坐标记作 \((X,Y)\)。令

\[
L(t)^{1-\beta}=c(1-\beta)(T-t)=B(T-t),
\]

则在上半平面得到特别简单的精确物理公式

\[
\boxed{
\Psi_{\rm af}(t,X,Y)
=\frac{Y(X-kY)}{T-t},
\qquad
\rho_{\rm af}(t,X,Y)-\rho_c
=\frac{2k(X-kY)}{T-t}.
}
\tag{AF13}
\]

令 \(\tau=T-t\)，则相应速度为

\[
U_{\rm af}=\frac1\tau(-X+2kY,Y).
\]

于是可在物理变量中直接核验

\[
-\Delta\Psi_{\rm af}=\frac{2k}{\tau}=\partial_X\rho_{\rm af},
\qquad
\partial_t\rho_{\rm af}
=\frac{2k(X-kY)}{\tau^2}
=-U_{\rm af}\cdot\nabla\rho_{\rm af}.
\tag{AF13a}
\]

它的水平导数和速度应变为

\[
\partial_X\rho_{\rm af}=\frac{2k}{T-t}>0,
\qquad
\operatorname{spec}(\nabla U_{\rm af})
=\left\{-\frac1{T-t},\frac1{T-t}\right\}.
\tag{AF14}
\]

若坚持把原参数 \(a(t)\) 本身定义为物质顶点，则 (8.3b) 对该 affine family 给出

\[
b=U_1(0)=(1-\beta)b,
\]

故 \(\beta>0\) 时必须 \(b=0\)。因此非零 \(b\) 在这里不是一个新的行波速度，
而只是把坐标原点从真实停滞点偏移了 \(-(b/c)L\)。

公式 (AF13) 是一个精确有限时解，但它在空间上为 affine，在每个 \(t<T\) 就有
无界密度/无限能远场。因此它严格证明了“所需符号、右扇、壁跳和双曲聚焦可在
IPM 中同时精确出现”，但它还不是从有界密度初值产生的梯度爆破。
它的正确用途是作为下面 Bessel/return gluing 的精确动态内核。

还必须注意一个缩放规范事实：反变换后的 (AF13) 完全不含 \(\beta\) 或 \(c\)。
同一个全局 affine 物理解可以用任意 \(0<\beta<1\) 的 \(L(t)\) 重参数化；它没有
一个由自身选出的收缩 core 尺度。因而任意指定 \(\beta\) 的写法都不等于已经
构造出该指数的有界多尺度爆破。只有外层局域化、终端 cusp 与可解性条件才可能
真正选出 \(\beta\) 和 \(L(t)\)。

### 8.3 势能投影给出的缩放系数符号

假设暂时可以对 (8.3) 在整个上半平面分部积分，且无穷远的输运与椭圆边界项均为零。
记

\[
I_1=\int_{z_2>0}z_2R(z)\,dz,
\qquad
E=\int_{z_2>0}|\nabla P(z)|^2\,dz.
\]

用 \(z_2\) 乘 (8.3) 的输运方程后积分，并使用

\[
\int R U_2\,dz=-E,
\qquad
\nabla\cdot(z_2z)=3z_2,
\]

可得到精确投影恒等式

\[
\boxed{E=(\beta+3)cI_1.}
\tag{8.5a}
\]

所以若 \(I_1>0\)，无外通量情形下必然 \(c>0\)。但第 3 节的 Casimir no-go 又表明，
一个非平凡局域化 profile 不可能真的同时满足这些零通量假设。对截断在匹配边界的
真正多尺度解，(8.5a) 应改成

\[
E+\mathcal B_{\rm out}=(\beta+3)cI_1,
\tag{8.5b}
\]

为固定符号约定，设 \(V=U+cz-be_1\)，并在截断域 \(\Omega\) 上定义

\[
\begin{aligned}
B_{\rm tr}
&=\oint_{\partial\Omega}z_2R\,V\cdot n\,ds,\\
B_{\rm ell}
&=\oint_{\partial\Omega}
\bigl(PRn_1+P\partial_nP\bigr)\,ds,\\
\mathcal B_{\rm out}&=B_{\rm tr}-B_{\rm ell}.
\end{aligned}
\tag{8.5d}
\]

壁面上因 \(z_2=P=0\) 没有贡献；(8.5d) 只记录外匹配边界上的 return/reservoir
通量。于是
“\(c>0\) 吗？”已经化成一个定量符号问题：

\[
c=\frac{E+\mathcal B_{\rm out}}{(\beta+3)I_1}.
\tag{8.5c}
\]

只要 \(\mathcal B_{\rm out}>-E\) 且 \(I_1>0\)，就强制聚焦符号 \(c>0\)；
\(\mathcal B_{\rm out}=-E\) 时只能得到 \(c=0\)。

### 8.4 为什么 \(D\) 必须与 profile 同阶

按 (8.2)，

\[
\rho_t\sim\frac{\delta^2}{L},
\qquad
U\sim\delta.
\]

而 (2.2) 中动力项的大小为

\[
U\cdot\nabla D
\sim\frac{\delta}{L}\,D_{\rm amp}.
\]

在当前的 self-induced、非退化平衡中（\(U\sim\delta\)、\(\dot L\sim-\delta\)，没有
\(O(1)\) 外加背景流，也没有 leading 正交消去），平衡时间变化自然要求

\[
\widetilde D_{\rm amp}\sim\delta,
\qquad Q=O(R).
\tag{8.6}
\]

因此在这一缩放分支上，动态解不应是
\(\widetilde D=o(\rho-\rho_c)\) 的 Bessel 小扰动。若存在外部背景流或结构性消去，
则必须重新做 dominant balance；(8.6) 不是仅由 (2.2) 无条件推出的 no-go。

纯 Bessel 闭包在重整化变量中是

\[
Q=0,
\qquad R=-P.
\]

若 \(c>0\)，(8.5) 就会要求 \(R\) 严格 \(\beta\)-齐次；Bessel 函数并不齐次，故
稳态本身不能是收缩固定点。

### 8.5 由固定点一旦存在就强制得到的缩放律

本小节沿用第 7 节的单尺度识别 \(L=\mu^{-1}\)，即 Bessel turnover 与内自相似
尺度同阶（第 8.7 节的 \(q\asymp1\) 端点）。\(L,\delta\) 与梯度率由 fixed-point
缩放决定；\(\mu,A\) 的下述幂律只属于这个最小转换分支。

由 (8.1)--(8.2)，

\[
L(t)^{1-\beta}
=c(1-\beta)(T-t).
\tag{8.7}
\]

所以

\[
L(t)\asymp(T-t)^{1/(1-\beta)},
\tag{8.8}
\]

\[
\delta(t)\asymp(T-t)^{\beta/(1-\beta)},
\tag{8.9}
\]

\[
\mu(t)=L(t)^{-1}
\asymp(T-t)^{-1/(1-\beta)},
\tag{8.10}
\]

\[
\|\nabla\rho(t)\|_\infty
\asymp\frac{\delta(t)}{L(t)}
\asymp\frac1{T-t}.
\tag{8.11}
\]

在 Bessel 转换尺度 \(r\sim L=\mu^{-1}\) 处，需要

\[
A(t)\sim\delta(t)L(t)=L(t)^{\beta+1}.
\tag{8.12}
\]

在上述单尺度分支中，这些量由尺度不变性强制。第 8.2 节已经给出满足壁面与符号条件的局部
affine 固定点；未知的不再是“局部 fixed point 是否存在”，而是能否把它延拓为
满足终端 cusp、外层通量、有限初值和全部 Casimir 条件的全局复合解，以及该
复合边值问题选出的 \(c\) 是否严格为正。

### 8.6 \(\beta\) 必须由全局匹配选择

目前没有可靠数据或解析可解性条件选择某个具体 \(\beta\)，所以不再把任何数值
指数写入主动结论。一般地，若 terminal cusp 为

\[
\rho(T,x)-\rho(T,(1,0))\sim r^\beta F(\theta),
\qquad 0<\beta<1,
\tag{8.13}
\]

则非退化扇区内

\[
|\nabla\rho(T,x)|\sim r^{\beta-1},
\]

而时间缩放仍由 (8.7)--(8.12) 给出。第 8.8 节的 leading 扇外真空匹配还给出

\[
\theta_*=\frac{\pi\beta}{1+\beta},
\qquad
\beta=\frac{\theta_*}{\pi-\theta_*}.
\tag{8.13a}
\]

但这只是特定 free-boundary 渐近类中的条件关系。只有稳定测得 \(\theta_*\) 与径向
斜率，或由完整 affine-to-cusp Fredholm 问题得到离散可解性条件，才能确定
\(\beta\)。在此之前它是未知本征量。

### 8.7 有界振幅允许的 Bessel 尺度窗口

本节额外假设外层匹配已经选出 \(\beta\)，并把相应候选内尺度记为 \(\ell\)：

\[
\tau=T-t,
\qquad
\ell\asymp\tau^{1/(1-\beta)},
\qquad
\delta=\ell^\beta\asymp\frac{\ell}{\tau}.
\tag{8.14}
\]

这个关系正是 affine-to-cusp 匹配而非 affine 解单独选出的：在 \(r\sim\ell\)
处，affine 密度差为 \(\ell/\tau\)，而 terminal cusp 要求它为 \(\ell^\beta\)；
两者相等才给出 \(\ell^{1-\beta}\asymp\tau\)。

精确 affine 法形在物理距离 \(r\) 上的密度差为

\[
|\rho_{\rm af}-\rho_c|\asymp\frac r\tau.
\tag{8.15}
\]

所以若它在半径 \(R_{\rm af}(t)\) 内近似成立，而物理密度保持一致有界，则只得到
上界

\[
R_{\rm af}(t)\lesssim\tau=\frac{\ell}{\delta}.
\tag{8.16}
\]

它可以更早过渡到第 8.8 节的 cusp；(8.16) 不是方程强制出的第二尺度。只有在
affine 密度差增长到 \(O(1)\) 才发生过渡的**最大延伸/饱和情形**，才有
\(R_{\rm af}\asymp\tau\)。

现在令 Bessel 长度 \(L_B=\mu^{-1}\)，并记 \(q=\mu\ell=\ell/L_B\le1\)。
若用 \(K_1\) 近场在 \(r=\ell\) 匹配 \(\delta\)，则

\[
A\asymp\delta\ell=ell^{\beta+1}.
\tag{8.17}
\]

而 (4.4) 要求匹配圆上的非恒定补偿 oscillation 至少为

\[
M_{\rm corr}\gtrsim\frac{A}{\mu\ell^2}
\asymp\frac{\delta}{q}.
\tag{8.18}
\]

这给出一个比单个“distinguished law”更准确的尺度窗口：

1. 若只有与内核同阶的 \(M_{\rm corr}\asymp\delta\)，则必须
   \[
   q\asymp1,
   \qquad L_B\asymp\ell,
   \qquad \mu\asymp\ell^{-1}.
   \tag{8.19}
   \]
   此时 Bessel turnover 可在内核边缘出现，但不是一个分离尺度。
2. 若允许补偿层仍保持有界 \(M_{\rm corr}=O(1)\)，则至多有
   \[
   q\gtrsim\delta,
   \qquad
   \ell\lesssim L_B\lesssim\frac{\ell}{\delta}.
   \tag{8.20}
   \]
   最大分离的条件饱和端点是
   \[
   q\asymp\delta,
   \qquad
   L_B\asymp\frac{\ell}{\delta}=\ell^{1-\beta},
   \qquad
   \mu\asymp\ell^{\beta-1}.
   \tag{8.21}
   \]

用未知的一般 \(\beta\) 写，两个端点以及共同的形式匹配振幅为

\[
\boxed{
\begin{gathered}
\ell\asymp\tau^{1/(1-\beta)},\qquad
\delta\asymp\tau^{\beta/(1-\beta)},\qquad
A\asymp\tau^{(1+\beta)/(1-\beta)},\\
L_B\asymp\tau^{1/(1-\beta)},\quad
\mu\asymp\tau^{-1/(1-\beta)}
\qquad\text{或}\qquad
L_B\asymp\tau,\quad\mu\asymp\tau^{-1}.
\end{gathered}}
\tag{8.22}
\]

最大分离端点不是已经闭合的构造。它要求 \(O(1)\) 的**非恒定**壁面/过渡层矩；
常数密度背景或常数壁跳不贡献 (4.3) 所测的水平偶极矩。更严重的是，一个位于
\(\ell\) 尺度的泛型 \(O(1)\) 补偿层会产生 \(O(1)\) 速度和
\(O(1/\ell)\) 应变，超过 affine 分支的 \(O(\delta)\) 速度与
\(O(\delta/\ell)\) 应变。除非能证明明确的矩投影/应变消去，它会把动力学改成
另一条 Type-II 分支。

等价地，若在最大分离端点仍只让 affine 内核自己供矩，则 (4.4) 只允许

\[
A\lesssim\mu\ell^2\delta\asymp\ell^{2\beta+1},
\tag{8.22a}
\]

对应 \(r=\ell\) 处仅 \(O(\delta^2)\) 的低阶 Bessel 响应。因而分离的 Bessel
层必然是强 closure/\(D\) 耦合问题，不能被当作 affine 内核的小线性尾；衰减的
\(K_1\) 项也不能独自连接向外增长的 affine 场，主匹配必须包含下面的
homogeneous cusp/\(Q\) 分量。

### 8.8 角向硬支撑锐扇的终端 cusp 与第一修正

精确 affine 解的正相虽是锐扇，却不能直接零延拓成锐扇支撑。另一方面，作为
\(|z|\to\infty\) 的 leading terminal cusp，椭圆方程确实允许一个扇外真空的
角向硬支撑锐扇（它沿半径无界，并非空间紧支撑），
而且去掉与加入 \(X\)-反射对称会给出不同结论。

先不施加 \(X\)-反射对称。在整个上半平面 \(0<\theta<\pi\) 中令

\[
m=\beta+1,
\qquad
\theta_*=\pi-\frac\pi m=\frac{\beta}{\beta+1}\pi,
\qquad
q=r\sin(\theta_* -\theta),
\tag{8.23}
\]

并取

\[
R_0=A(q_+)^\beta.
\tag{8.24}
\]

角度 (8.23) 不是预先拟合的。内部一个特解是

\[
P_{\rm part}=-\frac{A\sin\theta_*}{m}(q_+)^m.
\]

为同时令 \(P=0\) 于右壁 \(\theta=0\)，必须加入一个在扇边为零的内部调和项
\(C_*r^m\sin(m(\theta-\theta_*))\)。扇外真空区的解也必须调和，并在左壁
\(\theta=\pi\) 为零。由于 \(P,\partial_\theta P\) 要跨扇边连续且
\(C_*\neq0\)，必有

\[
\sin\bigl(m(\pi-\theta_*)\bigr)=0.
\]

对 \(1<m<2\) 与锐角 \(0<\theta_*<\pi/2\)，唯一可能是第一模
\(m(\pi-\theta_*)=\pi\)，这正是 (8.23) 的角度公式。

定义

\[
C_*=-\frac{A(\sin\theta_*)^{m+1}}
{m\sin(m\theta_*)},
\tag{8.25}
\]

以及

\[
P_0=r^m
\begin{cases}
-\dfrac{A\sin\theta_*}{m}\sin^m(\theta_* -\theta)
+C_*\sin\bigl(m(\theta-\theta_*)\bigr),
&0<\theta<\theta_*,\\[2mm]
C_*\sin\bigl(m(\pi-\theta)\bigr),
&\theta_*<\theta<\pi.
\end{cases}
\tag{8.26}
\]

直接计算得到

\[
-\Delta P_0=\partial_XR_0.
\tag{8.27}
\]

而且 \(P_0=0\) 于 \(\theta=0,\pi\)，\(P_0\) 与 \(\partial_\theta P_0\)
在 \(\theta=\theta_*\) 连续，所以扇边没有椭圆 sheet，且扇边是 streamline。
在支撑内，若 \(A>0\)，

\[
\partial_XR_0
=A\beta\sin\theta_*\,q^{\beta-1}>0.
\tag{8.28}
\]

因此由两侧 wall 与真空 harmonic gap 强制出的角度及其反解为

\[
\boxed{
\theta_*=\frac{\pi\beta}{1+\beta},
\qquad
\beta=\frac{\theta_*}{\pi-\theta_*}.
}
\tag{8.29}
\]

对任意 \(0<\beta<1\)，它都给出“以 \((1,0)\) 为顶点、向右、锐角支撑、第一象限
\(\partial_X\rho\) 统一为正”的 leading terminal profile。若把 \(A\) 变号，
导数统一变为负。

这里必须保留一个层级区别：(8.24)--(8.27) 精确解了椭圆方程，并且其
\(r^\beta\) dilation 项精确相消，但它还不是完整非线性 fixed point，因为
\(U_0\cdot\nabla R_0\) 一般不为零。远离扇边与角点的形式渐近展开中，下一项的
齐次次数为 \(2\beta-1\)，并被强制为

\[
R_1=\frac{U_0\cdot\nabla R_0}{c(1-\beta)},
\qquad
-\Delta P_1=\partial_XR_1.
\tag{8.30}
\]

对一般 \(\beta\)，反缩放后的物理贡献为

\[
\ell^\beta R_1(r/\ell)
\sim \ell^{1-\beta}r^{2\beta-1}
\asymp (T-t)r^{2\beta-1},
\tag{8.31}
\]

所以它在固定外点趋于零，并且在重整化远场中比终端 \(r^\beta\) cusp 低阶。
当 \(\beta<1/2\) 时这一修正随 \(r\) 衰减；\(\beta\ge1/2\) 时虽不衰减，
相对主项的比值仍为 \(r^{\beta-1}\to0\)。

但硬支撑一般不能逐阶保持。若还强迫 \(P_1\) 在同一真空 gap 中调和并在两端为零，
它的次数 \(2\beta\) 会要求新的共振

\[
2\beta(\pi-\theta_*)=n\pi.
\tag{8.31a}
\]

代入 (8.29) 后，左端为 \(2\beta\pi/(1+\beta)\in(0,\pi)\)，对所有
\(0<\beta<1\) 都不等于正整数倍 \(\pi\)，除非内部 edge 导数发生额外的偶然消去。
因此 (8.29) 只是 leading Poisson/free-boundary 支撑角；完整 fixed point 很可能
需要次阶 signed return layer 或轻微弯曲的 edge，不能把“全阶精确硬支撑”当作
已经证明。

事实上还可以排除把 leading 项本身误当成完整精确解。假设在
\(0<\theta<\theta_*<\pi/2\) 内有单符号的纯齐次 fixed point

\[
R=r^\beta F(\theta),
\qquad
P=r^{\beta+1}G(\theta),
\tag{8.31b}
\]

并要求 \(P=0\) 于壁面、零延拓跨扇边无 elliptic sheet。齐次 dilation 项此时
精确相消，因此输运式强制 \(J(P,R)=0\)。在连通单相区中，由齐次性必有

\[
G=C_0|F|^{(\beta+1)/\beta}
\tag{8.31c}
\]

（符号固定）。若 wall trace 非零，\(P|_{\theta=0}=0\) 已直接与这个函数关系
矛盾；故只能令 \(F(0)=0\)。记
\(m=(\beta+1)/\beta\)，Poisson 角方程为

\[
-\bigl(G''+(\beta+1)^2G\bigr)
=\beta\cos\theta\,F-\sin\theta\,F'.
\tag{8.31d}
\]

非退化端点平衡分别强制

\[
F(\theta)\sim A\theta^{2\beta}
\quad(\theta\downarrow0),
\qquad
F(\theta)\sim B(\theta_*-\theta)^\beta
\quad(\theta\uparrow\theta_*).
\tag{8.31e}
\]

对单相 \(A,B\) 同号，(8.31d) 的右端在 wall 端符号为 \(-A\)，在 edge 端
符号为 \(+B\)；而由 (8.31c) 算出的左端在两端具有同一个由 \(C_0\) 决定的
符号，矛盾。因此

\[
\boxed{
\text{单符号、纯齐次、无 sheet 的精确硬支撑锐扇不存在。}
}
\tag{8.31f}
\]

这不否定 (8.23)--(8.29) 作为 leading Poisson/free-boundary 外几何；它严格说明
完整 nonlinear tail 必须同时含低齐次修正、edge bending、外 harmonic 场或
signed return，而不能只保留一个 \(r^\beta F(\theta)\)。

这个 homogeneous tail 也精确解释了 Casimir 从哪里流出。对
\(V=U+cz-be_1\) 和任意光滑 \(F\)，固定点方程给出

\[
\oint_{\partial\Omega}F(R)V\cdot n\,ds
=c\int_\Omega\bigl(2F(R)+\beta R F'(R)\bigr)\,dz.
\tag{8.32}
\]

取 \(F(R)=|R|^p\) 时，\(r^\beta\) 尾在半径 \(R_{\rm out}\) 上的通量与
\(R_{\rm out}^{\beta p+2}\) 同阶，正好等于体积分所需的
\((\beta p+2)c\) 倍；它不是可删去的小尾。

更一般地，若某个窄角尾满足

\[
R\sim r^\gamma\mathcal F(r^\sigma\theta),
\qquad U=o(r),
\tag{8.32a}
\]

则在大圆上比较 (8.32) 的幂次会给

\[
p(\gamma-\beta)=\sigma.
\]

让两个不同的 \(p\) 同时成立，必然强制

\[
\boxed{\sigma=0,
\qquad \gamma=\beta.}
\tag{8.32b}
\]

而 Bessel--Gaussian 尾分别对应

\[
K_0:\ (\gamma,\sigma)=(-1/2,1/2),
\qquad
K_1:\ (\gamma,\sigma)=(-1,1/2).
\tag{8.32c}
\]

所以 Bessel 窄扇严格不能承担最终 fixed-point 外通量；它至多是从 growing cusp
中减去后的中间 Green 响应。

若额外施加局部 \(X\)-反射对称，把真空 gap 截在 \(\theta=\pi/2\)，则 harmonic
匹配要求

\[
m\left(\frac\pi2-\theta_*\right)=n\pi.
\tag{8.33}
\]

对 \(1<m<2\) 与锐角 \(0<\theta_*<\pi/2\)，没有可行的正整数 \(n\)。因此这种
“单个紧支撑右扇 + 真空 gap”的简单终端结构与 \(X\)-反射对称不相容；去掉该
对称后，(8.23)--(8.29) 则给出唯一的第一模角度。最自然的全局目标由此变成

\[
\boxed{
z=O(1)\text{ 的精确 affine 双曲内核}
\longrightarrow
\text{closure/Bessel 过渡}
\longrightarrow
r^\beta\text{ 角向硬支撑终端 cusp}
\longrightarrow
\text{全局等测度 reservoir}.
}
\tag{8.34}
\]

若令 affine 内核的零线角与 (8.29) 对齐，(AF8b) 给出唯一的局部
斜率--聚焦率比。但 affine 解在扇外是负相；沿固定点特征线
\(dR/ds=\beta cR\) 时符号不能自行消失。因此这个负相必须进入真正的
edge/return layer，不能直接截掉。角度对齐只是 gluing 的必要归一化，不是已经
得到的全局解。

壁面一阶方程本身没有阻止这个异宿过渡。钉住物质顶点后，记
\(R_w(X)=R(X,0^+)\)，则

\[
\bigl(U_1(X,0)-U_1(0,0)+cX\bigr)R_w'(X)
=\beta cR_w(X).
\tag{8.35}
\]

affine jet 的 \(R_w\sim\alpha X\) 要求括号中的总壁面速度为
\(\beta cX\)，而 terminal cusp 的 \(R_w\sim CX^\beta\) 要求它趋于 \(cX\)。
两端在一阶壁面动力学上相容；真正困难在非局域 Poisson/edge/return 匹配。

### 8.9 最简单 affine-to-power 桥的精确 no-go

可以把“真正困难”进一步写成一个严格的降维障碍。令

\[
s=X-kY,
\qquad Q_k=1+k^2,
\qquad P=Yh(s).
\]

Poisson 方程的一般积分为

\[
R=2kh-Q_kYh'+F(Y).
\tag{8.36}
\]

在 \((s,Y)\) 坐标中，重整化特征速度恰好是

\[
\dot s=cs-h(s),
\qquad
\dot Y=Y(h'+c).
\tag{8.37}
\]

若只看壁面且取 \(F(0)=0\)，则存在一个很诱人的 formal heteroclinic：

\[
(cs-h)h'=\beta ch.
\tag{8.38}
\]

令 \(q_h=h/(cs)\)，其正分支由

\[
\boxed{
\frac{(1-\beta-q_h)^\beta}{q_h}
=C s^{1-\beta},
\qquad 0<q_h<1-\beta
}
\tag{8.39}
\]

隐式给出，并满足

\[
h(s)=(1-\beta)cs-O(s^{1/\beta})\qquad(s\downarrow0),
\qquad
h(s)\sim C_\infty s^\beta\qquad(s\to\infty).
\tag{8.40}
\]

而且 \(0<h'<(1-\beta)c\)、\(h''<0\)，所以形式 Poisson profile 保持

\[
R_X=2kh'-Q_kYh''>0
\tag{8.41}
\]

于整个右扇。这说明 wall trace 本身确实允许 affine-to-power 过渡。

但完整二维输运把它严格排除。把 (8.38) 代回 bulk 方程，剩余残差正好是

\[
\boxed{\mathcal E=-2Q_kY(h')^2.}
\tag{8.42}
\]

任何只依赖 \(Y\) 的 \(F\) 都不能消掉它；完整变量分离进一步强制连续的 \(h'\)
为常数。因此 \(P=Yh(X-kY)\) 类中不存在精确 affine-to-power/Bessel
heteroclinic。

这个结论也不是少加一个法向项造成的。更一般地，取 \(N\ge2\) 并令

\[
P=\sum_{j=1}^N p_j(X)Y^j,
\qquad
R=\sum_{j=0}^N r_j(X)Y^j.
\]

Poisson 的最高系数给 \(r_N=-p_N'+A_N\)，而输运的 \(Y^{2N-1}\) 系数只能来自
Jacobian，并强制

\[
N(p_N'r_N-p_Nr_N')=0.
\]

所以

\[
p_N'=A_N+\lambda p_N
\]

对某个常数 \(\lambda\) 成立，最高横向模只能是 affine/指数型。可是
\(r^{\beta+1}\) 的有限角多项式尾会要求
\(p_N(X)\sim X^{\beta+1-N}\)，其指数在 \(N\ge2\) 时为负，与上述 ODE 不相容。
逐次降阶仍回到上面的 \(N=1\) no-go。所以真正的 bridge 必须含无限角模或非多项式
二维 wall layer。Bessel corrector 的具体任务，就是抵消 (8.42) 并把 inner 所需的
横向系数从 \(Q_k(1-\beta)c\) 调制到 outer 的零。

分类还留下一个精确但不局域化的对照族。令

\[
B=(1-\beta)c,
\qquad
\nu=\frac{\beta}{2-\beta},
\]

则对任意常数 \(K\)，

\[
\boxed{
P=BY(X-kY),
\qquad
R=2Bk(X-kY)+K Y^\nu
}
\tag{8.43}
\]

逐点满足完整 fixed-point 与 Poisson 方程，并保持 \(R_X=2Bk>0\)。附加项是唯一
允许的 \(X\)-独立被动 return eigenmode，可以弯曲零线；但它在 generic rays 上仍
被 affine 线性增长支配，且 \(\nu<1\) 使壁法向导数奇异，所以不能完成
bounded affine-to-cusp gluing。

### 8.10 一个精确多齐次异宿，但几何不符合目标

还有一个有用的严格对照：IPM fixed point 的确可以精确地从 affine 内渐近过渡到
\(r^\beta\) 外渐近，只是最简单的这一个没有锐角、wall jump 与所需水平 jet。
令

\[
P=cYH(X),
\qquad
R=-cYH'(X),
\qquad
a_\beta=1-\beta.
\tag{8.44}
\]

置

\[
F=X-H,
\qquad p=H',
\]

则 Poisson 方程自动成立，而 fixed-point 输运恰好化成

\[
Fp'+a_\beta p+p^2=0,
\qquad
F'=1-p.
\tag{8.45}
\]

在 \(-a_\beta<p<0\) 的分支上有第一积分

\[
\boxed{
\frac{F^{a_\beta}(-p)}{(p+a_\beta)^{a_\beta+1}}=K>0.
}
\tag{8.46}
\]

它唯一确定一个全局异宿，并满足

\[
\begin{aligned}
X\downarrow0:&\quad
H'(X)=-a_\beta+O\!\left(X^{a_\beta/(a_\beta+1)}\right),\\
X\to\infty:&\quad
-H'(X)\sim K a_\beta^{a_\beta+1}X^{-a_\beta}.
\end{aligned}
\tag{8.47}
\]

所以

\[
P\sim-a_\beta cXY,
\quad R\sim a_\beta cY
\qquad(X\downarrow0),
\tag{8.48}
\]

而

\[
P\sim-\frac{cK a_\beta^{a_\beta+1}}\beta YX^\beta,
\quad
R\sim cK a_\beta^{a_\beta+1}YX^{\beta-1}
\qquad(X\to\infty).
\tag{8.49}
\]

由于 \(p'>0\)，它有统一 \(R_X=-cYp'<0\)。把 \(P,R\) 在 \(X<0\) 置零时，
\(X=0\) 上 Poisson 方程两侧的 delta 系数完全相等，而且总法向速度为零，所以这是
一个右向 \(90^\circ\) 扇的精确弱 fixed point。

令

\[
L^{1-\beta}=c(1-\beta)(T-t),
\]

其物理解在内区满足

\[
\Psi\sim-\frac{XY}{T-t},
\qquad
\rho\sim\frac{Y}{T-t},
\tag{8.50}
\]

而外区趋于时间无关的 \(YX^{\beta-1}\) cusp。这是一个真正精确的多齐次自相似
bridge，并再次说明方程本身没有选择单个数值 \(\beta\)。但它只有直角支撑、
\(R|_{Y=0}=0\)、尖端正则性不足，且在 \(Y\) 方向无界/无限能；因此它只是严格的
near-miss，不是所需的有界数据锐扇爆破。

## 9. 条件路线 B：\(K_1\) 诱导的局域双尺度 closure layer

还有一个更偏 Type-II 的局域机制。设水平微核尺度为 \(\ell\)，法向 closure
层宽度为 \(w\)，且

\[
X=\frac{x_1-a}{\ell},
\qquad
Y=\frac{x_2}{w},
\qquad
\rho=MR,
\qquad
\Psi=M\ell P.
\tag{9.1}
\]

在 \(\ell\ll w\) 时，精确重整化 Poisson 方程为

\[
-\left(P_{XX}+\frac{\ell^2}{w^2}P_{YY}\right)=R_X.
\tag{9.2}
\]

对应的局域定常输运方程为

\[
(-P_Y+c_xX+c_a)R_X
+(P_X+c_yY)R_Y=0,
\tag{9.3}
\]

其中

\[
c_x=-\frac{w\dot\ell}{M\ell},
\qquad
c_y=-\frac{\dot w}{M},
\qquad
c_a=-\frac{w\dot a}{M\ell}.
\tag{9.4}
\]

要从偶极矩得到宽度下界，还需明确一个匹配假设：在 \(r\sim\ell\) 处，
\(K_1\) 密度 \(A/\ell\) 与振幅为 \(M\) 的微核同阶，即

\[
A\asymp M\ell,
\tag{9.4a}
\]

并且 closure 层中的相关角向 oscillation 不超过 \(O(M)\)。在这个额外假设下，
(4.4) 才给出

\[
w^2\gtrsim\frac{\ell}{\mu}.
\tag{9.5}
\]

若 (9.4a) 不成立、偶极系数更小或 leading pole 被另一层抵消，则 (9.5) 不是独立
结论。下面的显式缩放律另外假设该下界被饱和，即
\(w^2\asymp\ell/\mu\)。若 \(\mu\) 固定，则 \(c_y/c_x=1/2\)。
若带下游 Casimir 通量的完整内层边值问题
能够选出 \(c_x=c_*>0\)，那么

\[
\dot\ell=-c_*M\frac{\ell}{w}
\asymp-c_*M\sqrt{\mu\ell}.
\tag{9.6}
\]

从而形式上得到

\[
\boxed{
\ell(t)\asymp(T-t)^2,
\qquad
w(t)\asymp T-t,
}
\tag{9.7}
\]

以及

\[
|\partial_{x_1}\rho|\asymp(T-t)^{-2},
\qquad
|\partial_{x_2}\rho|\asymp(T-t)^{-1}.
\tag{9.8}
\]

速度梯度不应被统一记成 \(w^{-1}\)。从 (9.1) 逐分量得

\[
\partial_{x_1}U_2=\frac{M}{\ell}P_{XX},
\qquad
\partial_{x_1}U_1, \partial_{x_2}U_2=O\!\left(\frac{M}{w}\right),
\qquad
\partial_{x_2}U_1=O\!\left(\frac{M\ell}{w^2}\right).
\tag{9.8a}
\]

因此对泛型 \(P_{XX}=O(1)\)，\(|\nabla U|\) 的最大分量与
\(\ell^{-1}\asymp(T-t)^{-2}\) 同阶；只有横向混合分量是 \((T-t)^{-1}\) 量级。

但必须非常明确地说：(9.7)--(9.8) 只是

\[
\boxed{\text{formal local closure-layer law},}
\]

不是全局爆破解。原因有两个：

1. 在微核 \(r=\ell\) 处令 \(q=\mu\ell\ll1\)，多极矩已经证明 \(D\) 层比
   \(K_1\) 分量大 \(1/q\)。在饱和宽度 \(r=w\) 处该比值仍为
   \(1/(\mu w)\asymp1/\sqrt q\)。所以这是 \(D\)-主导爆破，不是 Bessel-主导爆破。
2. 若 (9.3) 的 profile 在无穷远衰减且没有通量，由
   \(\nabla\cdot(U+c_xXe_1+c_yYe_2)=c_x+c_y>0\)，乘以 \(|R|^{p-2}R\)
   积分会直接推出 \(R=0\)。

因此 Type-II 局域律要成为物理解，必须与一个非局域等测度 reservoir/return layer
同时构造。

当前恒等式并没有单独推出 \(K_0\) 的尺度 ODE；尤其 \(K_0\) 尖点主速度是旋转而非
双曲聚焦。因而不能仅凭其对数近场严格断言“只有无限时收缩”，也不能据此得到
有限时爆破。任何 \(K_0\) 调制律都必须从完整 closure/return 边值问题的投影重新导出。

## 10. 右向有效集中区、统一符号区与 Bessel 抛物尾

在下游 \(X\gg L=\mu^{-1}\) 且 \(Y=O(\sqrt{LX})\) 时，

\[
\rho_0
\sim A\sqrt{\frac{\pi}{LX}}
\exp\!\left(-\frac{Y^2}{4LX}\right),
\tag{10.1}
\]

\[
\rho_1
\sim A\sqrt{\frac{\pi}{L}}
\frac{Y}{X^{3/2}}
\exp\!\left(-\frac{Y^2}{4LX}\right).
\tag{10.2}
\]

所以这个外层不是紧支撑圆锥，而只有右向张开的抛物有效集中区：

\[
|Y|\sim\sqrt{LX},
\qquad
|\theta|\sim\sqrt{\frac{L}{X}}.
\tag{10.3}
\]

对 \(A>0\)，远场中

\[
\partial_X\rho_0<0
\quad\text{when}\quad
Y^2<2LX,
\tag{10.4}
\]

\[
\partial_X\rho_1<0
\quad\text{when}\quad
Y^2<6LX.
\tag{10.5}
\]

把 \(A\) 变号就会把导数符号整体翻转。因此 Bessel 外层确实能给出一个
“右向张开、其中 \(\partial_X\rho\) 统一符号”的有效集中区，但不是支撑集。

但在该区外导数必然换号。只有在整个上半平面 \(X\in\mathbb R\) 中同时要求
\(X\to\pm\infty\) 水平衰减，并要求每个固定 \(Y\) 上的 \(\partial_X\rho\) 全局
同号，才会由单调性强制 profile 平凡；这不能误套到只有 \(X>0\) 的第一象限。
对 Bessel 外响应，正确的符号条件应只施加在 (10.4) 或 (10.5) 所给出的动态
抛物区内。真正扇外真空的锐扇候选是第 8.8 节的 terminal cusp，而不是这里的
Gaussian 尾。

## 11. 把精确 Bessel 稳态作为背景：扰动桥接方程

### 11.1 一般稳态背景上的精确扰动方程

先不使用 Bessel 的线性闭包。对任意固定稳态背景
\((\rho_B,\Psi_B,U_B)\)，写

\[
\rho=\rho_B+\eta,
\qquad
\Psi=\Psi_B+\phi,
\qquad
u=\nabla^\perp\phi.
\]

则扰动严格满足

\[
\boxed{
\begin{aligned}
&\eta_t+U_B\cdot\nabla\eta
+u\cdot\nabla\rho_B
+u\cdot\nabla\eta=0,\\
&-\Delta\phi=\partial_X\eta.
\end{aligned}}
\tag{11.0}
\]

这已经是“固定外稳态 + 动态内扰动”的精确起点。它没有假设扰动在内区小；
只要外区所用的背景满足相同的物理边界/源 contract，(11.0) 就没有近似。

### 11.2 \(D_B=0\) Bessel 背景下的一变量化

先对固定 \(\mu\)、固定 tip 分布源（或穿孔域固定内边界）的 \(D_B=0\) Bessel
基态写

\[
\rho=\rho_B+\eta,
\qquad
\Psi=\Psi_B+\phi,
\qquad
d=\eta+\mu\phi.
\]

令 \(T_\mu d=\phi\) 由

\[
(\Delta-\mu\partial_X)\phi=-\partial_Xd
\tag{11.1}
\]

定义。利用 Jacobian 的反对称性，完整非线性扰动方程精确化为

\[
\boxed{
(I-\mu T_\mu)d_t
+(U_B+U_d)\cdot\nabla d=0.
}
\tag{11.2}
\]

同时总场由单个变量 \(d\) 精确重建：

\[
\boxed{
\Psi=\Psi_B+T_\mu d,
\qquad
\rho=\rho_B+(I-\mu T_\mu)d,
\qquad
U=U_B+\nabla^\perp T_\mu d.
}
\tag{11.2a}
\]

因此不应先用 cutoff 手工拼 \(\rho\) 和 \(\Psi\)；只需构造 \(d\)，椭圆桥接就由
\(T_\mu\) 自动完成。在线性化中，

\[
\boxed{
\sigma(I-\mu T_\mu)d
+U_B\cdot\nabla d=0.
}
\tag{11.3}
\]

固定的纯 \(D_B=0\) 背景还有一个谱限制。记
\(M_\mu=I-\mu T_\mu\)。在全平面，或在沿 \(Y=0\) 取 Dirichlet 条件的上半平面，
对实、衰减且 source-neutral 的非零 \(d\)，Fourier/正弦变换给

\[
\langle d,M_\mu d\rangle
=\int
\frac{|k|^4}{|k|^4+\mu^2k_X^2}
|\widehat d(k)|^2\,dk>0.
\tag{11.3a}
\]

若 (11.3) 有实特征值 \(\sigma\)，由于算子系数为实，可取实特征函数；在所有
边界都无通量时

\[
0=\sigma\langle d,M_\mu d\rangle
+\int d\,U_B\cdot\nabla d
=\sigma\langle d,M_\mu d\rangle.
\tag{11.3b}
\]

故纯固定 Bessel 背景没有任何非零的**实**点谱特征值。若它有线性增长，只能来自
真正非自伴的复特征值 \(\Re\sigma>0\)；而在正则复合背景中，(11.5) 的
\(u_d\cdot\nabla D_B\) 会破坏 (11.3b)，实增长模也可能重新出现。这个结论没有
证明纯背景完全稳定，但排除了把某个实“压缩本征值”直接从 Bessel 参数族读出。

在全平面 Fourier 变量中，

\[
\widehat{T_\mu}(k)
=\frac{ik_X}{|k|^2+i\mu k_X},
\qquad
\widehat{I-\mu T_\mu}(k)
=\frac{|k|^2}{|k|^2+i\mu k_X}.
\tag{11.4}
\]

(11.4) 对 \(k\neq0\) 使用。零频、tip source 和外部 homogeneous Bessel 模必须由
source/moment 条件及边界规范另行固定；而且 \(M_\mu^{-1}\) 在低频并非一致有界。
因此 (11.3) 的自然算子域必须是带 source-neutrality 与加权远场条件的空间，不能
把它当作朴素无约束的 \(L^2\) 本征问题。

设内扰动的长度为 \(\ell\)，\(q=\mu\ell\)。对其典型频率
\(|k|\sim\ell^{-1}\)，

\[
\widehat{I-\mu T_\mu}=1+O(q),
\qquad
\widehat{\mu T_\mu}=O(q)
\quad(q\ll1).
\tag{11.4a}
\]

所以这个算子本身给出一个清楚的桥：

- 在 \(r\sim\ell\ll\mu^{-1}\) 的高频内区，\(d\) 近似就是物理密度扰动，
  演化的 principal part 是总速度对 \(d\) 的输运；
- 在 \(r\sim\mu^{-1}\) 的外区，\(I-\mu T_\mu\) 的完整非局域结构不可忽略，
  并自动产生 Bessel 转换。

在半平面 Dirichlet contract 中，外流函数扰动就是

\[
\phi(z)=-\int\partial_XG_{\mu,H}(z,z')d(z')\,dz',
\tag{11.4b}
\]

其中 \(G_{\mu,H}\) 是 (7.5)--(7.6) 的有量纲版本。这是目前最直接的
inner-to-outer 精确桥接公式。

### 11.3 为什么这里的“扰动”在内区不能小

第 4 节的多极矩立即作用到 (11.2a)。若 \(q\ll1\)，为了把奇异 Bessel tip
改造成有界内核，修正相对 \(K_0\) 或 \(K_1\) 至少分别大

\[
\frac1{q\log(1/q)}
\qquad\text{或}\qquad
\frac1q.
\tag{11.4c}
\]

所以这是一个**外区小、内区非微扰**的 singular perturbation，而不是在全域范数
中围绕 Bessel 做小扰动。对不可穿壁的 \(K_1\) 背景尤其清楚：
\(\Psi_B=\rho_B=0\) 于壁面，而所需的非零 wall trace 全部由
\(d|_{Y=0}\) 承担；它必与主 profile 同阶。

另外，(11.1)--(11.2a) 在分布意义下默认保持 Bessel 的固定 tip source，或把
tip 当作穿孔内边界。若目标是全域无源物理解，记

\[
S_B=L_\mu\Psi_B,
\]

则总场必须改用

\[
L_\mu\phi=-\partial_Xd-S_B,
\tag{11.4d}
\]

或者等价地加入 source-canceling wall/return particular solution。这项是 leading
约束，不能藏进小误差。最稳妥的使用方式，是只在不含 tip 的外匹配环带把 Bessel
当背景，并由内层与 return 层共同给它边界数据。

### 11.4 壁面扰动的精确压缩判据

对 Dirichlet \(K_1\) 背景，若总场也满足 \(\Psi|_{Y=0}=0\)，则壁上
\(\rho=d\)，且完整非线性方程退化为

\[
d_t+U_1d_X=0
\qquad (Y=0).
\tag{11.4e}
\]

沿壁面特征 \(\dot X=U_1(t,X,0)\)，

\[
d_X(t,X(t))
=d_X(0,X_0)
\exp\!\left[-\int_0^t\partial_XU_1(s,X(s),0)\,ds\right].
\tag{11.4f}
\]

所以真正需要证明的是总壁面应变的负部产生发散积分。原始 \(K_1\) 背景在
\(X>0\) 的壁面速度为

\[
U_{B,1}(X,0)
=A e^{\mu X/2}\frac{K_1(\mu X/2)}{X}
\sim\frac{2A}{\mu X^2}
\qquad(X\downarrow0).
\tag{11.4g}
\]

若 \(A<0\)，流向 tip，但 \(\partial_XU_{B,1}>0\)，是扩张而非压缩；若
\(A>0\)，应变为负但流向右外侧。故 raw \(K_1\) 近场不能同时给出所需的
“向 tip 运动 + 负水平应变”。\(K_0\) 的 tip 主部又是旋转。于是 desired affine
压缩必须由 \(u_d\) 在内区改变 leading strain；这再次说明不能只研究一个全域小的
线性扰动。

这个方向冲突还能直接积分。只保留 (11.4g) 的 tip 主项并记
\(C=2A/\mu\)，壁面特征与梯度满足

\[
X(t)^3=X_0^3+3Ct,
\qquad
d_X(t)=d_X(0)\left(\frac{X(t)}{X_0}\right)^2.
\tag{11.4ga}
\]

当 \(A<0\) 时特征有限时撞向已有奇异 tip，但扰动梯度反而趋于零；当 \(A>0\)
时梯度只随向外运动作无限时代数增长。故 raw \(K_1\) 壁流本身没有产生所需的
有限时压缩梯度爆破。

### 11.5 固定背景谱、参数调制与 exact-outer no-go

这是应该真正计算的非对称谱问题。局域上，\(K_0\) 顶点流场的主部是旋转：

\[
\Psi_0\sim A\log r,
\qquad
U_0\sim A\frac{(-Y,X)}{r^2},
\]

它不是直接的双曲聚焦流。\(K_1\) 有开放的偶极扇区，所以比 \(K_0\) 更值得继续，
但必须把壁层、return layer 和非局域 Casimir 通量一起放入算子域。

决定 (8.2) 中 \(c>0\) 或 (9.6) 中 \(c_*>0\) 的，不是 Bessel 公式本身，而是这个
完整内外层算子的伴随 Fredholm 投影。

若允许背景参数 \(\lambda=(A,\mu,a,\ldots)\) 随时间调制，但每个冻结的
\((\rho_B^\lambda,\Psi_B^\lambda)\) 都是稳态，则不用猜测可以直接从 (11.0) 得到

\[
\boxed{
\eta_t+U_B^\lambda\cdot\nabla\eta
+u\cdot\nabla\rho_B^\lambda
+u\cdot\nabla\eta
=-\dot\lambda\cdot\partial_\lambda\rho_B^\lambda.
}
\tag{11.4h}
\]

把扰动对伴随的平移、振幅、尺度零模取正交，便得到有限维调制系统

\[
\mathcal M(\lambda,d)\dot\lambda=\mathcal F(\lambda,d).
\tag{11.4i}
\]

这才是决定 \(\dot\mu\)、中心漂移和可能聚焦率的方程；直接令稳态参数随时间变化
而不保留右端强迫是不成立的。

还有一个局部化 no-go。若在收缩 fixed-point 坐标的某个完整外开集令
\(Q=R+P=0\)，则该区同时满足

\[
(\Delta-\partial_X)P=0,
\qquad
(cz-be_1)\cdot\nabla P=\beta cP.
\tag{11.4j}
\]

平移中心后第二式令 \(P\) 为 \(\beta\)-齐次，而第一式含不同齐次次数的
\(\Delta P\) 与 \(P_X\)。对 \(0<\beta<1\)，二者只能分别为零，进而
\(P=0\)。所以不能存在“有限过渡层之外恰好等于非零 Bessel”的收缩 fixed point；
Bessel 只能是带非零但渐小 \(Q\) 的外渐近响应，或固定物理尺度的时间依赖背景。

把这个结论写成有源方程更适合实际求解。归一化 \(\mu=1\)，令
\(L=\Delta-\partial_X\)，并写

\[
P=P_B+\phi,
\qquad Q=q,
\qquad LP_B=0,
\qquad L\phi=-q_X.
\]

固定点方程精确等价于

\[
\boxed{
\begin{aligned}
&U_B\cdot\nabla q+U_\phi\cdot\nabla q
+(cz-be_1)\cdot\nabla(q-\phi)-\beta c(q-\phi)\\
&\hspace{25mm}
=(cz-be_1)\cdot\nabla P_B-\beta cP_B.
\end{aligned}}
\tag{11.4k}
\]

右端是 Bessel 背景不可消掉的 **scaling defect**。钉住 \(b=0\) 时，\(K_0\)
近极点的 defect 含 \(\beta c\log r\) 主项；一般 \(K_\nu\) pole
\(P_B\sim r^{-\nu}\) 的 defect 为

\[
-(\nu+\beta)cP_B
\tag{11.4l}
\]

同阶。故在自相似 bridge 中，\(q\) 必须与 Bessel pole 同阶地响应；不能把它当
小 remainder。

甚至“affine 内核 + 一个精确衰减 Bessel 项”的直接叠加也失败。令

\[
P_A=B Y(X-kY),
\quad R_A=2Bk(X-kY),
\quad B=(1-\beta)c,
\]

并设另一个分量满足 \(R_B=-mP_B\)、
\((\Delta-m\partial_X)P_B=0\)。若总和仍是 fixed point，交叉输运强制

\[
\beta c q\,\partial_qP_B
+\left((2-\beta)cY+\frac{2Bk}{m}\right)\partial_YP_B
=\beta cP_B,
\qquad q=X-kY.
\tag{11.4m}
\]

沿正向特征，\(q\) 与平移后的 \(Y\) 都指数趋向无穷，而 \(P_B\) 必须按
\(e^{\beta cs}\) 增长；与 Bessel 衰减矛盾。因此唯一衰减解是 \(P_B=0\)。
这排除了最简单的线性叠加，bridge 必须通过非线性、同阶的 \(q=d\) 完成。

### 11.6 光滑背景、奇异背景与内层尺度

固定物理 Bessel 场若在候选 core 点光滑，只能作为低阶外背景。取
\(d\) 的内幅度 \(\delta=\ell^\beta\)，中心随背景速度移动，并令
\(\dot\ell=-c\delta\)。在 \(z=O(1)\) 的内变量中，冻结背景留下的两个典型系数为

\[
\frac{\bar U(a+\ell z)-\bar U(a)}{\delta}
=O(\ell^{1-\beta}),
\qquad
\frac{\ell}{\delta}\nabla\bar\rho(a+\ell z)
=O(\ell^{1-\beta}).
\tag{11.4n}
\]

它们对所有 \(0<\beta<1\) 都趋于零，所以光滑固定背景不能在 leading order
关闭 affine tail 或选择 \(\beta\)；它只提供外边界条件和低阶应变。

若把 Bessel pole 正好放在 tip，固定振幅的 \(K_0/K_1\) 项反而在内变量中发散，
而且初态已经奇异。要让它在 inner equation 中保持 \(O(1)\)，必须令
\(A(t)\to0\) 并通常同时调制 \(\mu(t)\)；此时 (11.4h) 的参数残差与主非线性
同阶，已经不是“固定稳态 + 小扰动”。因此稳态背景法的严谨定位是：

\[
\boxed{
\text{Bessel 负责外层 Green/谱结构；order-one 的 }d
\text{ 负责内层压缩与尺度选择。}
}
\tag{11.4o}
\]

### 11.7 中性模与 growing mode 的 Casimir 正交条件

固定 \(\mu\) 时，Bessel 振幅和平移切向仍在 \(D=0\) homogeneous kernel 中，
只是中性参数漂移。\(\mu\)-切向满足

\[
d_\mu=-\Psi_B,
\qquad
T_\mu d_\mu=\partial_\mu\Psi_B,
\qquad
\sigma=0.
\tag{11.4p}
\]

所以显然的 scale tangent 也是中性模，不是已经找到的增长模。

此外，对任何可积扰动/有限截断域且无扰动边界通量的设置，全部 Casimir 的一阶
变分守恒。若
\(\eta=(I-\mu T_\mu)d\) 是 \(\sigma\neq0\) 的线性特征模，则对所有合适的
\(F\)，

\[
\int F'(\rho_B)\eta\,dx=0.
\tag{11.4q}
\]

由 coarea 公式，这等价于几乎每条背景 level set 上

\[
\int_{\{\rho_B=s\}}
\frac{\eta}{|\nabla\rho_B|}\,d\ell=0.
\tag{11.4r}
\]

因此任何真正 growing eigenmode 都必须沿背景等值线强烈变号；单符号 fan 或纯
振幅扰动不可能是增长模。这给谱计算提供了无限组必须显式实施的伴随正交条件。

### 11.8 正则复合背景的完整线性化

还要强调，(11.2)--(11.3) 不是正则化复合基态的完整线性化。若 gluing 后的基态
有 \(D_B\neq0\)，并且 \(U_B\cdot\nabla D_B=0\)，则精确扰动方程为

\[
(I-\mu T_\mu)d_t
+U_B\cdot\nabla d
+u_d\cdot\nabla D_B
+u_d\cdot\nabla d=0,
\tag{11.5}
\]

线性算子必须保留 \(u_d\cdot\nabla D_B\)。若 tip/source 参数本身也随时间调制，
还会出现相应的参数导数强迫。因此简单 Bessel 谱 (11.3) 只能审计固定奇异核，
不能单独决定完整 gluing 的 Fredholm 符号。

## 12. 最小未排除的多尺度几何

与所有上述约束同时相容的最小结构是：

1. **精确双曲 Taylor 内核**：在候选 \(r\sim\ell\) 区由 (AF2) 组织，
   \(Q=R+P\) 与 \(R\) 同阶，并满足 (AF5b) 的强制 jet 关系。
2. **壁面/扇边 closure 层**：承担非零 wall trace、跨壁跳跃、统一符号的
   \(R_X\)，同时消去直接截断 affine 正相产生的斜边 sheet，并把负相转入 return。
3. **terminal cusp 匹配**：无 \(X\)-反射对称时，以 (8.23)--(8.31) 的
   \(r^\beta\) leading 角向硬支撑锐扇作为重整化远场；它提供不可删除的 Casimir
   外通量，次阶则允许 signed return/edge 修正。
4. **Bessel 背景/closure 扰动桥**：在 (8.20) 允许的 \(L_B\) 窗口内，用
   (11.2a) 的单变量 \(d\) 重建总场，并由 (7.5)--(7.7) 精确表示外响应；
   \(d\) 在外区可小、在内区必须 order one。Bessel 可以有 Gaussian 抛物有效
   集中区，但不是最终支撑边界，也不能在 \(q\ll1\) 时独自主导。
5. **全局 return/reservoir/filament**：消掉 tip 源和远场矩，并保持每一个 density
   level set 的测度。

最后一层还有一个明确的速率必要条件。对一般各向异性固定点

\[
(U+c_xXe_1+c_yYe_2)\cdot\nabla R=c_\rho R,
\tag{12.1}
\]

若没有外通量，则每个 \(p\) 都满足

\[
(pc_\rho+c_x+c_y)\int|R|^p=0.
\tag{12.2}
\]

用两个不同的 \(p\) 就强制 \(c_\rho=0\) 与 \(c_x+c_y=0\)。所以任何改变振幅的
vanishing-contrast 层都不能自行局域化；一个 \(O(1)\) return filament 若要在
无额外通量下保存全部 Casimir，只能采用面积保持的双曲速率（一向收缩、另一向
等速扩张）。壁面本身因 \(V\cdot n=0\) 不能提供这部分通量。

return 也不能被关进一个有限的 material cell。若某条闭曲线（或与不可穿壁边界
合围的曲线）处处 characteristic，即 \(V\cdot n=0\)，则散度定理给

\[
0=\oint_{\partial\Omega}V\cdot n,ds
=2c|\Omega|,
\tag{12.3}
\]

与 \(c>0\)、\(|\Omega|>0\) 矛盾。因此 return 必须开向无穷远并携带通量，
或者是连续的 \(D\)-主导层；不能是封闭的 piecewise-affine 补偿泡。

若施加水平左右反射对称，简单的扇外真空单扇终端尾已被 (8.33) 排除；去掉该对称
后，需增加平移正交条件，但对 affine jet 而言非零 \(b\) 只是坐标规范。无论是否
对称，多极矩与 Casimir 障碍都不会消失。

## 13. 下一个可判定的理论任务

下一步按“先背景谱、再非线性 gluing”的顺序进行，而不预设 \(\beta\)：

1. 先选择清楚的背景 contract：固定 tip/穿孔域用于纯谱审计；真正物理解则先构造
   source-canceling 的正则复合背景 \((\rho_B,\Psi_B,D_B)\)。
2. 在该背景上实现 (11.3) 或 (11.5) 的谱问题，算子域同时包含不可穿壁、远端
   source neutrality 以及 (11.4q)--(11.4r) 的 Casimir 正交条件。
3. 检查是否存在 \(\Re\sigma>0\) 且能产生负 wall strain 的变号模。对纯
   \(D_B=0\) 背景，(11.3a)--(11.3b) 已排除非零实特征值，故只能寻找复增长模；
   对正则复合背景则必须保留 (11.5)。振幅、平移和 \(\mu\)-切向都是中性模，
   不能冒充不稳定性。
4. 若有不稳定模，沿其非线性延拓求解 (11.2)，直到 \(d\) 在内区达到 order one；
   这时才切换到 (AF5a)--(AF5b) 的 affine Taylor core 描述。
5. 在固定点层面解 (11.4k) 的强迫 heteroclinic：内端为 affine jet，外端只要求
   Bessel-weighted 渐近而不能令 \(q\) 紧支撑；同时用 (8.42) 作为首个必须被二维
   corrector 消掉的显式 forcing。
6. 更外端匹配一般 \(r^\beta\) cusp/return flux，并把 \(\beta,c,\mu\)、平移、法向
   宽度和 return length 都作为未知调制参数。
7. 对伴随核施加 (8.32)、(12.2) 与源矩条件。只有某个
   \(0<\beta<1\)、\(c>0\) 被 Fredholm 系统选出且误差可逆，(8.7)--(8.12) 才升级
   为有界数据爆破律。

因此当前真正的“一锤定音”步骤，是先判定 Bessel/复合背景是否存在合格的
compressive growing mode；若存在，再证明其非线性轨道进入 affine-to-cusp
heteroclinic。若谱中没有这种模，固定稳态背景路线就被关闭。

## 14. 与已知理论的关系

目前检索到的严格 IPM 爆破结果并没有使用本文这一 Bessel 内外层机制；
这只能说明“尚未找到相同构造”，不能当作新颖性证明。相关的方法路线包括：

- [Chae (2006)](https://arxiv.org/abs/math/0601060) 对 divergence-free transport 自相似爆破的非存在机制，
  与本文的 Casimir 排除一致。
- [Elgindi (2014)](https://arxiv.org/abs/1411.6958) 对 IPM 平滑衰减稳态的分层性与稳定性分析，
  说明 Bessel 必须依靠奇点和非标准远场逃离经典分类。
- [Elgindi--Pasqualotto (2023)](https://arxiv.org/abs/2310.19780) 的“从不稳定模到奇点”路线，
  适合在 (11.3) 真的出现聚焦不稳定模时使用。
- [Collot--Prange--Tan (2025)](https://arxiv.org/abs/2507.17381) 建立了特殊无限能 IPM
  自相似爆破的稳定性，表明“重整化稳态 + 调制 + 加权估计”是可行的，
  但其 profile 不是本文 Bessel 稳态。
- [Dembski (2025)](https://arxiv.org/abs/2511.01827) 给出边界密度为零的 IPM 有限时奇点，
  其锐扇、加权空间和不稳定模控制方法对当前问题有直接参考价值。

## 15. 最终判断

现在可以把结论说得很精确：

1. **已经得到的完整精确解：**(AF13) 是带右向正相锐扇、统一
   \(\partial_X\rho\) 符号、双曲临界点和奇延拓壁跳的 affine 有限时解；但它
   无界、无限能，而且其 \(\beta\) 只是重整化规范。
2. **已经得到的强制局部结构：**任何光滑候选 fixed point 的非退化壁上 Taylor
   jet 都必须是 (AF5a)--(AF5b)。所以 affine core 是正确的局部法形，不是拟合。
3. **已经严格排除的机制：**固定 \(\mu\)、有界平滑内核、small-argument Bessel
   主中间层与局域小 \(D\) 的简单三层爆破。
4. **Bessel 的精确用途：**对标准 localized \(K_0/K_1\)，它给出 closure
   方程的 Green 核、source/dipole 外层基与右向 Gaussian 抛物响应；
   对 directional exponential/Jordan mode，它还可与 order-one harmonic residual
   精确闭合。两者都不是最终的 source-neutral 紧支撑锐扇。
5. **背景扰动方程已经闭合：**固定 source/穿孔域中，(11.2)--(11.2a) 把整个
   inner--outer 问题化为单个 \(d\) 的精确非线性演化。可是多极矩、wall trace 与
   scaling defect 都证明 \(d\) 在内区必须 order one；简单的小扰动或 affine+Bessel
   线性叠加严格不行。纯 \(D_B=0\) 背景还没有非零实点谱；可能的增长只能是复模，
   或来自含 \(D_B\neq0\) 的正则复合背景。
6. **存在一个精确多尺度 near-miss：**(8.44)--(8.50) 精确连接 affine 内渐近与
   \(r^\beta\) 外渐近，并有统一 \(R_X\) 符号；但它是直角扇、wall trace 为零、
   无限能且尖端不够正则。它证明多齐次异宿在方程中可能存在，也证明方程没有
   自动选择某个数值 \(\beta\)。
7. **新的精确 leading 外几何：**不施加 \(X\)-反射对称时，椭圆方程与无 sheet
   匹配强制一般 \(\beta\) 的 leading 扇外真空终端扇角为
   \(\theta_*=\pi\beta/(1+\beta)\)；(8.30) 给出第一非线性修正，但完整 fixed
   point 尚未构造，也尚未选择 \(\beta\)。而 (8.31b)--(8.31f) 进一步证明：
   单符号、纯齐次、无 sheet 的硬支撑锐扇本身不可能是完整精确 fixed point；
   低齐次修正、弯曲 edge、外 harmonic 场或 signed return 至少要出现一种。
8. **仍未被排除但尚未证明的：**affine Taylor core、同阶 \(Q/D\) 壁/扇边层、
   Bessel Green 响应、\(r^\beta\) terminal cusp 与面积保持 return filament 组成的
   有界多尺度爆破。
9. **已找到的 bounded-density 精确爆破：**在允许不衰减 harmonic
   far-field contract 的全平面，\(\rho=\tanh((x-a)/(T-t))\) 是完整 IPM
   轨道，\(\|\rho\|_\infty=1\) 而 \(\|\rho_x\|_\infty=(T-t)^{-1}\)。
10. **已找到的 Bessel 收敛机制：**水平壁 directional profile 族
    \(R_*=(2-\beta)cYe^{qX}\) 在右向、单调、\(Y\)-线性不变类的
    \(X\ge0\) characteristic sector 中
    具有完全显式的非线性特征公式，并严格 \(C^1_{\rm loc}\) 收敛到
    \(q=r_0'(0)/((2-\beta)c)<0\) 的 Bessel 指数模；统一 \(R_X<0\) 符号保持。
11. **更强的光滑初值 Bessel 收敛轨道：**存在有界、水平紧支撑的
    \(C^\infty\) 物理密度，其完整精确动态重标度解局部
    \(C^\infty\) 收敛到 \((P_*,R_*)=(cXY+e^{-X},e^{-X})\)。它显式具有
    \(L\)、\(Lgs\)、\(L(gs)^2\) 三种长度；固定内区梯度是 Type-I，
    最大梯度位于外部 return 尺度并更快增长。
12. **固定合同下的纯 Bessel 长时极限：**在一个时间无关的常双曲应变
    far-field contract 下，水平紧支撑 \(C^\infty\) 初值的唯一一维演化
    在二次平坦边缘尺度上严格 \(C^\infty_{\rm loc}\) 收敛到
    \((e^{-X},e^{-X})\)，且 \(D_*=0\)。它是 \(t=\infty\) tangent，
    不是有限时爆破或稳定性定理。
13. **横向局域化的精确边界：**wall-compatible 正弦轨道可光滑收敛到
    wall-linear Bessel，但左远场无界；旋转条带轨道有界且每条横竖切片
    紧支撑，但二维总支撑无限。单 transverse mode、有限 Fourier 与小
    finite-energy envelope 均被严格排除。
14. **仍开放的物理目标：**第 9--13 项依赖无穷能调和应变，而水平壁
    Bessel 族的某个远场方向无界。标准衰减、有限能、source-neutral、
    带右向锐扇和 wall jump 的 Cauchy 爆破定理仍需二维 return gluing。
15. **含时收敛态的精确分类：**势能恒等式排除固定物理区域中向
    非零有限能 Bessel 的长时强收敛；两个不同 \(L^p\) Casimir 排除
    shrinking profile 的全局重整化强收敛。未被排除且已有显式例子的，
    只是固定重整化紧集上的局部收敛；若要求标准有限能数据并保留
    Bessel core，还必须同时携带 order-one \(D\)/return 层与外部 reservoir。

因此答案是：

\[
\boxed{
\begin{gathered}
\text{Bessel bridge 已经产生精确奇化轨道，并在一个自然不变类中}\\
\text{具有严格非线性局部收敛；尚未证明的是用有限能 return 外层}\\
\text{代替奇化的调和远场，并保留水平壁、锐扇、跳跃与全部 Casimir。}
\end{gathered}
}
\]

新的精确公式与稳定性证明集中在
[BESSEL_BRIDGE_EXACT_BLOWUP_ZH.md](BESSEL_BRIDGE_EXACT_BLOWUP_ZH.md)。严格有限能
存在性的剩余关键，是把其 directional core 放入满足全部
Casimir/source 条件的复合 Bessel/return 背景，再构造 wall/edge 与
等测度 return 层并证明二维谱稳定性。

## 16. 只要求全方向 \(R_X\to0\) 后的更新

这里必须把旧目标中的“有限能/全部 Casimir 有限”与新的硬条件分开。若现在
只要求

\[
R_X\in C_0(\overline{\mathbb H}),
\qquad R_X\not\equiv0,
\tag{16.1}
\]

则第 3、6、12 节中依赖 \(R\in L^p\)、有限势能或无穷远通量消失的全局排除
不再自动适用。方向依赖常数或 \(r^\beta\;(0<\beta<1)\) 的慢尾本身可以承担
非有限质量 reservoir。局部 Poisson 多极矩、wall \(D\)-层、纯 Bessel 静态性、
纯齐次 hard-fan 和有限 transverse-mode 障碍则原封不动地保留。

在 bounded-contrast 端点，新的精确 fixed-point 目标是

\[
(U+cz-be_1)\cdot\nabla R=0,
\qquad -\Delta P=R_X,
\qquad U=\nabla^\perp P,
\qquad P|_{Y=0}=0.
\tag{16.2}
\]

若存在有界光滑解且 \(R_X\in C_0\)、\(R_X(0,0)\ne0\)，则

\[
\rho(t,x,y)=\rho_c+R\!\left(
\frac{x-a(t)}{c(T-t)},\frac{y}{c(T-t)}\right)
\tag{16.3}
\]

是从光滑初值出发的精确 Type-I 水平梯度爆破。其 wall Taylor jet 仍被强制为
affine 核，而一个通过 leading Fredholm 审计的方向常数外端是

\[
F_0=A\left(\cos\theta+\frac53\cos3\theta\right),
\qquad
G_0=-\frac{4A}{3}\sin^4\theta,
\tag{16.4}
\]

满足 \(R_X=O(r^{-1})\)。第一次输运修正是
\(F_1=G_0F_0'/c\)。完整存在性现在等价于构造

\[
\text{affine wall jet}
\longleftrightarrow
\text{order-one infinite-mode }D\text{ layer}
\longleftrightarrow
\text{angular-constant tail}.
\tag{16.5}
\]

这条路线不再需要有限能 return，但仍不是 affine+Bessel 的线性叠加；全局
quadratic harmonic strain 还必须被外尾屏蔽。完整的条件重审、显式外端、
逐渐变宽单相位条带的排除，以及开放 harmonic contract 下的精确二维约化见
[RHOX_DECAY_BLOWUP_STATUS_ZH.md](RHOX_DECAY_BLOWUP_STATUS_ZH.md)。
