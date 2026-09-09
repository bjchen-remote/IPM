# 以 Bessel 漂移--Helmholtz 算子为桥的精确 IPM 梯度爆破

日期：2026-08-28

本文把 Bessel 从“待收缩的稳态 profile”改成两种更可靠的对象：

1. 精确 closure/Green 算子；
2. 被 order-one 余留场驱动的 directional Bessel mode。

这样确实能得到新的完整精确轨道。不过必须一开始就说明范围：本文找到的
有界密度例子位于全平面并带无穷能调和应变；满足水平不可穿壁的精确 Bessel
例子仍在远场无界。因此它们证明了机制真实存在，但尚不是标准衰减有限能半平面
初值的爆破定理。更精确地说，它们是在允许指定不衰减 harmonic
far-field contract 下的精确 IPM 轨道；该调和部分不由标准衰减
Biot--Savart 规范单独确定。

为了简化公式，显式解的终端点先写在 \((0,0)\)。将所有 \(x\)
换成 \(x_1-1\)，就把终端点平移到原问题的 \((1,0)\)；IPM 方程和
下面的所有验证在这个平移下不变。

## 0. 最短结论

最干净的有界密度例子是

\[
\boxed{
\rho(t,x,y)=\tanh\!\frac{x}{T-t},
\qquad
\Psi(t,x,y)=\frac{xy}{T-t}
-(T-t)\log\cosh\!\frac{x}{T-t}.
}
\tag{E0}
\]

它在允许上述 harmonic far-field contract 的全平面逐点满足完整 IPM：

\[
\rho_t+u\cdot\nabla\rho=0,
\qquad
-\Delta\Psi=\rho_x,
\qquad
u=\nabla^\perp\Psi,
\]

并且

\[
\|\rho(t)\|_{L^\infty}=1,
\qquad
\rho_x(t,0,y)=\frac1{T-t}.
\tag{E1}
\]

更一般地，任意 \(f\in C_b^1\) 的一维 profile 都能放进同一精确不变类；
重整化后 profile 完全冻结。当 \(f'\not\equiv0\) 时，这给出一个
无限维、类内中性稳定的 Type-I 梯度爆破族，而不是一次偶然抵消。

本文还回答了一个更强的问题：Bessel 能否从普通光滑初值中作为含时极限
出现？第 13.3 节构造了 \(x\)-紧支撑的
\(f\in C_c^\infty\)，使完整物理密度始终有界光滑，而动态重标度场满足

\[
(P,R)\longrightarrow(cXY+e^{-X},e^{-X})
\quad\hbox{in }C^\infty_{\rm loc}.
\]

该解的固定内区梯度为 Type-I，外部 return 层的最大梯度更快，并显式展示
\(L\)、\(L|\log(T-t)|\) 与 \(L|\log(T-t)|^2\) 三种长度。代价仍是密度
不随 \(y\) 衰减以及随时间增强的无穷能调和应变；所以这是开放边界
far-field contract 下的肯定例，而不是标准有限能 Cauchy 吸引子。

第 14 节进一步固定一个**时间无关**的双曲远场合同。对水平紧支撑的
\(C^\infty\) 初始密度，唯一的一维 classical 演化在 \(t\to\infty\) 的
二次平坦边缘尺度上满足

\[
(P,R)\longrightarrow(e^{-X},e^{-X})
\quad\hbox{in }C^\infty_{\rm loc},
\qquad D_*=0.
\]

所以 Bessel 确实可以是“光滑初值 + 固定 far-field contract”的含时纯
收敛态；但这只是无限时 tangent limit，尚未证明邻域稳定性，也不解决标准
有限能有限时爆破。

第 15 节还给出两个进一步的精确桥接：一个满足水平不可穿壁并从正弦法向模
光滑收敛到 wall-linear Bessel profile；另一个是有界、每条横竖切片都紧支撑
的旋转条带有限时爆破，并局部收敛到 \(e^{-X}\)。前者在左远场无界，后者
仍有无限长斜条带和无穷能。单 transverse mode、共同 \(y\)-cutoff、有限
Fourier 闭合和小有限能 tangent residual 随后都被严格排除，因而真正有限能
版本必须包含非可分、无限模、order-one 的二维 return。

它与 Bessel 的精确联系如下。固定任意无量纲 drift 参数 \(m\)，令

\[
D=R+mP.
\]

则椭圆方程严格变成

\[
(\Delta-m\partial_X)P=-\partial_XD.
\tag{E2}
\]

二维算子 \(\Delta-m\partial_X\) 的 Green 核正是指数加权的 \(K_0\)；在上述
共线不变类中，它降成一维指数 Green 核。爆破 profile 是 order-one residual
通过这个 Bessel 算子产生的精确响应。

另外还存在两个满足水平壁条件的精确族：一个给任意
\(0\le\beta<1\) 的单尺度 directional Bessel 内区奇化；另一个给出
水平展开、法向收缩、Bessel 长度更小的精确三尺度内区奇化。
它们在第 7、8 节给出。

## 1. 不预设 \(\beta\) 的精确 dynamic scaling

取

\[
(x_1,x_2)=a(t)e_1+L(t)z,
\qquad
\frac{ds}{dt}=\frac{\delta(t)}{L(t)},
\]

\[
\rho=\rho_c+\delta R,
\qquad
\Psi=\delta L P,
\qquad
u=\delta U,
\qquad
U=\nabla^\perp P.
\]

定义

\[
c=-\frac{L_s}{L},
\qquad
\gamma=-\frac{\delta_s}{\delta},
\qquad
b=\frac{a_s}{L}=\frac{\dot a}{\delta}.
\tag{E3}
\]

则物理 IPM 与下式严格等价：

\[
\boxed{
R_s+(U+cz-be_1)\cdot\nabla R=\gamma R,
\qquad
-\Delta P=R_X.
}
\tag{E4}
\]

只有在极限 \(c\to c_*>0\)、\(\gamma\to\gamma_*\) 已经由完整方程选出后，才定义

\[
\beta_*:=\frac{\gamma_*}{c_*}.
\tag{E5}
\]

本文不预设它的数值。

## 2. Bessel 背景加余留的完整非线性方程

令物理 Bessel 参数为 \(\mu_{\rm ph}(t)\)，并设

\[
m(s)=\mu_{\rm ph}(t)L(t),
\qquad
Q=R+mP.
\]

则

\[
L_mP:=(\Delta-m\partial_X)P=-Q_X.
\tag{E6}
\]

以下抽象分解先固定一个重整化区域 \(\Omega\)，以及不随 \(m\)
变化的边界、辐射和零模规范，并假设 \(L_m\) 在所用 graph space
上可逆。定义

\[
T_m:=L_m^{-1}(-\partial_X).
\]

选取任意 homogeneous Bessel 背景 \(L_mP_B=0\)，写

\[
P=P_B+T_md,
\qquad
Q=d,
\qquad
L_mT_md=-d_X,
\]

并记

\[
A_m:=I-mT_m.
\]

于是

\[
R=A_md-mP_B,
\qquad
U=U_B+U_d.
\tag{E7}
\]

没有丢弃任何项的余留方程是

\[
\boxed{
\begin{aligned}
A_md_s
&+(U_B+U_d)\cdot\nabla d
+(cz-be_1)\cdot\nabla(A_md)-\gamma A_md\\
=\;&\partial_s(mP_B)
+m(cz-be_1)\cdot\nabla P_B-\gamma mP_B
-(\partial_sA_m)d.
\end{aligned}}
\tag{E8}
\]

参数导数也可完全显式化：

\[
\partial_mT_m=-T_m^2,
\qquad
\partial_mA_m=-T_mA_m.
\tag{E9}
\]

(E9) 还要求这个 inverse contract 对 \(m\) 可微。若 \(P_B\) 取标准
\(K_0/K_1\) tip 核，则它在穿孔域是 homogeneous，但在全域分布
意义下带源 \(S_B=L_mP_B\)。全域无源方程必须改写为

\[
L_m\phi=-d_X-S_B,
\]

或先加入恰好取消 \(S_B\) 的 particular/return field。

所以“Bessel + 余留”不是近似分解；(E8) 是完整 IPM 的另一套精确变量。

若锁定重整化 Bessel 参数 \(m=m_*>0\) 及 normalized \(P_B\)，(E8) 成为自治
forced PDE。其右端含不可删除的 scaling defect

\[
m\big[(cz-be_1)\cdot\nabla P_B-\gamma P_B\big].
\tag{E10}
\]

因此 scaled 平衡点通常必须有 order-one 的 \(d_*\)；纯 Bessel 的 \(d=0\)
不是平衡点。

## 3. 重整化收敛严格推出有限时爆破

这一部分不依赖显式解。令

\[
H(s)=\frac{L(s)}{\delta(s)}.
\]

由 (E3)，

\[
H_s=-(c-\gamma)H.
\tag{E11}
\]

假设 \(L,δ>0\)，某个 classical exact orbit 在固定的重整化区域上
对全部 \(s\in[s_0,\infty)\) 存在，并且

\[
c(s)\to c_*>0,
\qquad
\gamma(s)\to\gamma_*,
\qquad
\lambda_*:=c_*-\gamma_*>0,
\]

并且某个精确轨道满足

\[
R(s)\to R_*
\quad\hbox{in }C^1_{\rm loc},
\qquad
\nabla R_*(z_0)\ne0.
\]

因为 \(dt/ds=H\)，

\[
T-t(s)=\int_s^\infty H(\sigma)\,d\sigma
\sim\frac{H(s)}{\lambda_*}.
\]

而 \(\nabla_x\rho=H^{-1}\nabla_zR\)，所以

\[
\boxed{
(T-t)\nabla_x\rho(t,a(t)e_1+L(t)z_0)
\longrightarrow
\frac{\nabla R_*(z_0)}{\lambda_*}.
}
\tag{E12}
\]

这给出了“Bessel-remainder 重整化收敛 \(\Rightarrow\) 物理梯度爆破”的精确
充分条件。(E12) 本身证明的是选定内点上的 \((T-t)^{-1}\) 渐近；
若还有全局的一致 \(C^1\) 上界，才能把整个解的梯度范数称为 Type-I。
若 \(0<\gamma_*<c_*\)，则 \(\delta\to0\)，是 density contrast 消失而梯度
发散的局部 Type-I 速率；要把对比消失提升为全局结论，还需
\(\sup_s\|R(s)\|_\infty<\infty\)。若还要爆破中心收敛到有限点，一个充分条件是

\[
\int_{s_0}^\infty L(s)|b(s)|\,ds<\infty,
\]

例如当 \(b\) 有界而 \(L\) 指数衰减时该条件自动成立。

## 4. 一个无限维精确 fixed-point 类

先取最简单的 \(\gamma=0,c=1,b=0\)。令

\[
\xi=X-kY,
\qquad
Q_k=1+k^2,
\]

并对任意 \(h\in C^2\) 定义

\[
\boxed{
P=Y\xi+h(\xi),
\qquad
R=2k\xi-Q_kh'(\xi).
}
\tag{E13}
\]

首先，

\[
-\Delta P=2k-Q_kh''=R_X.
\tag{E14}
\]

其次，若 \(V=\nabla^\perp P+(X,Y)\)，则

\[
V\cdot\nabla\xi=V_X-kV_Y=0.
\tag{E15}
\]

而 \(R\) 只依赖 \(\xi\)，故

\[
V\cdot\nabla R=0.
\tag{E16}
\]

所以 (E13) 对每个 \(h\) 都是完整非线性 fixed point。

反缩放时令 \(\tau=T-t\)，得到逐点精确物理解

\[
\boxed{
\begin{aligned}
\Psi(t,x,y)
&=\frac{y(x-ky)}{\tau}
+\tau h\!\left(\frac{x-ky}{\tau}\right),\\
\rho(t,x,y)-\rho_c
&=\frac{2k(x-ky)}{\tau}
-Q_kh'\!\left(\frac{x-ky}{\tau}\right).
\end{aligned}}
\tag{E17}
\]

允许中心 \(a(t)\) 时，只需把 \(x\) 换成 \(x-a(t)\)，并在 \(\Psi\) 中加入
\(-\dot a(t)y\) 的匀速 harmonic gauge。

## 5. 有界密度的完整精确爆破

在 (E17) 取 \(k=0\)。对任意 \(f\in C_b^1(\mathbb R)\)（若要任意阶
classical 规则性，可取 \(f\in C_b^\infty\)），令

\[
h'=-f.
\]

则

\[
\boxed{
\rho(t,x,y)=\rho_c+f\!\left(\frac{x}{T-t}\right)
}
\tag{E18}
\]

是精确解的密度分量，并且

\[
\rho_x(t,x,y)
=\frac1{T-t}f'\!\left(\frac{x}{T-t}\right).
\tag{E19}
\]

因此只要 \(f'\not\equiv0\)，就在有限时产生梯度爆破，同时

\[
\|\rho(t)-\rho_c\|_{L^\infty}=\|f\|_{L^\infty}
\]

严格保持。可以取 \(f\in C_c^\infty(\mathbb R)\)，得到水平方向紧支撑的平滑
密度；也可以取 \(f=\tanh\)，得到 (E0)。后者还有

\[
\rho(t)\to\operatorname{sgn}x,
\qquad
\rho_x(t)\rightharpoonup2\delta_{x=0}
\]

（逐个水平截面）。

这个解不是由密度值变大造成的；其机制是 density-free harmonic strain

\[
\Psi_H=\frac{xy}{T-t}
\]

把一个有界 density front 压到宽度 \(T-t\)，而 self-induced 速度始终沿密度
等值线，因而不破坏该压缩。

对 (E0) 直接计算

\[
U(t,0,0)=0,
\qquad
\nabla U(t,0,0)
=\frac1{T-t}
\begin{pmatrix}
-1&0\\
-1&1
\end{pmatrix}.
\tag{E19a}
\]

因此固定中心的顶点确实是一个双曲临界点，两个特征值为
\(\pm(T-t)^{-1}\)。它是有界被输送量的真正梯度爆破，但驱动它的
调和应变在无穷远不衰减。

## 6. Bessel Green 桥在精确族中的显式形式

回到 (E13)，固定 \(m\ne0\)，并定义

\[
D=R+mP.
\]

则

\[
D=2k\xi+mY\xi+\widetilde d(\xi),
\qquad
\widetilde d=mh-Q_kh'.
\tag{E20}
\]

同时

\[
(\Delta-m\partial_X)P=-D_X.
\tag{E21}
\]

给定 order-one 一维 residual \(\widetilde d\)，函数 \(h\) 由一阶 Green
方程精确决定：

\[
h'-\lambda h=-\frac {\widetilde d}{Q_k},
\qquad
\lambda=\frac m{Q_k},
\tag{E22}
\]

\[
\boxed{
h(\xi)=Ae^{\lambda\xi}
+\frac1{Q_k}\int_\xi^\infty
e^{\lambda(\xi-\sigma)}\widetilde d(\sigma)\,d\sigma,
}
\tag{E23}
\]

其中辐射方向应随 \(\operatorname{sgn}m\) 选择。

齐次项

\[
P_B=Ae^{\lambda\xi},
\qquad
R_B=-mP_B
\]

严格满足

\[
(\Delta-m\partial_X)P_B=0,
\qquad
D_B=0.
\tag{E24}
\]

它是 drift--Helmholtz/Bessel 算子的 directional mode。该二维算子的
基本解为 \(e^{mX/2}K_0(|m|r/2)\)；在当前共线不变子空间中，它的
解析精确降成指数 Green 公式 (E23)。这不要求把单个实指数模
本身误认为全局衰减的 \(K_0\)。

若 \(\widetilde d\) 紧支撑，适当选择 (E23) 的辐射端后，\(h\) 可以在一侧恰为 Bessel
齐次模、在另一侧恰为零。这是一个 order-one residual 完成的半显式精确 bridge，
不是“小 Bessel + 小修正”。

## 7. 任意 \(0\le\beta<1\) 的水平壁 Bessel-mode 爆破

上面的 bounded-density 共线族不满足水平 Dirichlet wall；但还有一个不同的
精确族确实满足 \(P|_{Y=0}=0\)。取

\[
0\le\beta<1,
\qquad
c>0,
\qquad
q\ne0,
\]

\[
b=\frac{(2-\beta)c}{q},
\]

这里的 \(b\) 同时就是 (E4) 中平移漂移 \(-be_1\) 的参数。

并定义

\[
\boxed{
P_*=cXY-bY e^{qX},
\qquad
R_*=(2-\beta)cY e^{qX}.
}
\tag{E25}
\]

它满足

\[
-\Delta P_*=(R_*)_X,
\qquad
(U_*+cz-be_1)\cdot\nabla R_*=\beta cR_*.
\tag{E26}
\]

验证中的关键抵消是

\[
(U_*+cz-be_1)\cdot\nabla R_*
=(2c-bq)R_*=\beta cR_*.
\]

而且

\[
P_*=0,
\qquad
R_*=0
\quad\hbox{于 }Y=0.
\]

把

\[
P_B=-bY e^{qX},
\qquad
P_H=cXY
\]

分开，则

\[
(\Delta-q\partial_X)P_B=0,
\qquad
R_*=-qP_B,
\]

而完整 closure residual 恰好是

\[
\boxed{
D_*=R_*+qP_*=qcXY.
}
\tag{E27}
\]

所以这里的机制非常明确：Bessel directional mode 携带密度；order-one harmonic
residual 产生双曲压缩；平移和 dilation defect 精确平衡。

令

\[
L^{1-\beta}=c(1-\beta)(T-t),
\qquad
\delta=L^\beta,
\qquad
a(t)=a_*-\frac bcL(t).
\]

取 \(a_*=1\)，中心就在 \(t\uparrow T\) 时收敛到用户指定的
\((1,0)\)。

令

\[
X=\frac{x-a(t)}{L(t)},
\qquad
Y=\frac y{L(t)}.
\]

则完整物理解显式为

\[
\boxed{
\begin{aligned}
\rho(t,x,y)-\rho_c
&=(2-\beta)cL^{\beta-1}y
  e^{q(x-a)/L},\\
\Psi(t,x,y)
&=L^\beta\!\left[
  \frac{c(x-a)y}{L}-by e^{q(x-a)/L}
  \right].
\end{aligned}}
\tag{E27a}
\]

若取 \(q<0\)，它在重整化右半区 \(X>0\) 向右指数衰减。这里
\(q\) 是带方向的 drift 参数，\(|q|/L\) 才是物理逆长度；
\(q<0\) 是标准方向分支的 \(X\)-反射。此时

\[
R_*>0,
\qquad
\partial_XR_*=q(2-\beta)cY e^{qX}<0
\quad(Y>0).
\tag{E28}
\]

其内区 density contrast 为 \(O(L^\beta)\)，但

\[
\partial_y\rho
=\frac{2-\beta}{(1-\beta)(T-t)}e^{qX},
\tag{E29}
\]

且在固定内坐标 \(X=\xi\)、\(Y=\eta>0\)，也就是物理点
\(x=a+L\xi\)、\(y=L\eta\)处，

\[
\partial_x\rho
=\frac{q(2-\beta)}{(1-\beta)(T-t)}
\eta e^{q\xi}.
\tag{E30}
\]

因此它是带水平不可穿壁、右向衰减分支和统一 \(\rho_x\) 符号的
精确内尺度 \((T-t)^{-1}\) 梯度增长解。其缺点是另一水平端指数增长、
法向线性增长和无限能；它的全局 \(L^\infty\) 范数从一开始就不有限，
所以不是“有界初始范数首次发散”的爆破。

## 8. 一个水平壁上的精确三尺度 Bessel 爆破

还可以让 Bessel mode 真正呈现三种不同尺度。记 \(\tau=T-t\)，取

\[
F=C e^{\alpha\tau(x-a)}
\sinh\!\left(\frac{\Gamma y}{\tau}\right),
\]

\[
\mu(t)=\alpha\tau+
\frac{\Gamma^2}{\alpha\tau^3},
\qquad
\alpha,\Gamma,C>0,
\tag{E31}
\]

并定义

\[
\boxed{
\rho=F,
\qquad
\Psi=-\frac{(x-a)y}{\tau}-\dot a\,y-\frac{F}{\mu(t)}.
}
\tag{E32}
\]

调和部分的两个物质不变量是

\[
I=\tau(x-a),
\qquad
J=\frac y\tau.
\]

所以它精确输运 \(F=Ce^{\alpha I}\sinh(\Gamma J)\)。另一方面

\[
\Delta F=\mu F_x,
\]

故 \(\Psi_B=-F/\mu\) 满足 Bessel closure，其速度与 \(\nabla F\) 正交。
这逐项证明了 (E32) 是完整 IPM 解。

它满足

\[
\Psi|_{y=0}=0,
\qquad
\rho|_{y=0}=0,
\]

\[
\rho_x=\alpha\tau\rho>0\quad (y>0),
\qquad
\rho_y(x,0,t)
=\frac{C\Gamma}{\tau}e^{\alpha\tau(x-a)}.
\tag{E33}
\]

三个尺度为

\[
L_x\sim\tau^{-1},
\qquad
w_y\sim\tau,
\qquad
L_B=\mu^{-1}
=\frac{\alpha\tau^3}{\Gamma^2+\alpha^2\tau^4}
\sim\tau^3,
\tag{E34}
\]

若取特征尺度 \(L_x=(\alpha\tau)^{-1}\)、\(w_y=\tau/\Gamma\)，则

\[
\frac{L_x}{\mu}
=\frac{\tau^2}{\Gamma^2+\alpha^2\tau^4}
\sim w_y^2
\qquad(\tau\downarrow0).
\tag{E35}
\]

因此 Bessel 算子确实能够参与一个精确的多尺度内区梯度增长。
这里的障碍是 \(y\to\infty\) 的指数增长、无锐扇和零 wall density trace。
特别地，对任意固定 \(y>0\)，\(\sinh(\Gamma y/\tau)\) 还会使密度本身
增长；因此这不是 bounded-density pure-gradient blow-up。

## 9. 精确稳定性：directional Bessel profile 族的扇区内局部吸引

### 9.1 有界密度共线类：中性但不吸引

对第 4--6 节的共线类，允许 \(h=h(s,\xi)\)。由于 (E15)，完整
重整化方程只剩

\[
\partial_s\partial_\xi h=0.
\]

模掉流函数的空间常数 gauge 后，

\[
\boxed{h_s=0.}
\tag{E36}
\]

所以该精确不变流形内的 profile distance 保持不变：这是非线性
中性/orbital stability，但没有向单一 \(h_*\) 的吸引。这个结论不涉及
离开共线流形的二维扰动。

### 9.2 水平壁 directional 类的完整非线性降维

现在固定

\[
0\le\beta<1,
\qquad c>0,
\qquad \Lambda=(2-\beta)c,
\]

并在 \(X\ge0\) 考虑自然的 \(Y\)-线性不变类

\[
P=cXY+Yh(s,X),
\qquad
R=Yr(s,X).
\tag{E36a}
\]

这里 \(X=0\) 是物质 separatrix，不是原 IPM 问题的水平物理壁。
以下定理是 \(X\ge0\) 右侧 characteristic sector 内的结论。要将它
视为整个 \(X\)-直线上的经典解，还需指定一个与右迹及所需高阶 jet
匹配的光滑 \(X<0\) 延拓。不能把右支简单零延拓。

选取 Poisson 分支 \(r=-h_X\)，并令 \(v=h+b\)。代入完整的
(E4)，不做线性化，精确得到

\[
\boxed{
r_s-vr_X=r(r-\Lambda),
\qquad
v_X=-r.
}
\tag{E36b}
\]

取右向 separatrix/gauge 条件

\[
v(s,0)=0,
\qquad
r(s,0)=\Lambda.
\tag{E36c}
\]

这是当前相对 phase space 中的边界归一化。具体地，可固定
\(b\) 并在每个 \(s\) 施加 \(h(s,0)=-b\)，或调制 \(b=b(s)\) 来持续强制
\(v(s,0)=0\)；它不是由 (E36b) 自动保持的流函数常数 gauge。
在该归一化下，边界方程为
\(r_s(s,0)=r(s,0)(r(s,0)-\Lambda)\)，所以 \(r(0)=\Lambda\) 的
codimension-one 条件才在演化中保持。
令特征流为

\[
\frac{dX}{ds}=-v(s,X),
\qquad
X(0,\alpha)=\alpha,
\]

并记

\[
C(\alpha)=\frac{\Lambda}{r_0(\alpha)}-1.
\]

沿特征线，\(u=1/r\) 满足线性方程 \(u_s=\Lambda u-1\)；
流映射 Jacobian \(J=X_\alpha\) 满足 \(J_s=rJ\)。因此有完全显式的
非线性公式

\[
\boxed{
r(s,X(s,\alpha))
=\frac{\Lambda}{1+C(\alpha)e^{\Lambda s}},
}
\tag{E36d}
\]

\[
\boxed{
X_\alpha(s,\alpha)
=e^{\Lambda s}
\frac{1+C(\alpha)}{1+C(\alpha)e^{\Lambda s}}.
}
\tag{E36e}
\]

此外，由链式法则还有精确导数公式

\[
r_X(s,X(s,\alpha))
=-\frac{\Lambda C'(\alpha)}
{(1+C(\alpha))(1+e^{\Lambda s}C(\alpha))}.
\tag{E36e'}
\]

### 9.3 严格局部收敛到 Bessel 指数模

假设 \(r_0\in C^1([0,\infty))\)，

\[
r_0(0)=\Lambda,
\qquad
0<r_0(\alpha)<\Lambda\quad(\alpha>0),
\tag{E36f}
\]

且对某个 \(\kappa>0\)、\(\theta>0\)，初值在 \(\alpha=0\) 附近满足
可逐项微分的展开

\[
C(\alpha)=\kappa\alpha+O(\alpha^{1+\theta}),
\qquad
C'(\alpha)=\kappa+O(\alpha^\theta),
\qquad
\kappa=-\frac{r_0'(0)}{\Lambda}.
\tag{E36g}
\]

在 (E36d)--(E36e) 中置 \(\alpha=e^{-\Lambda s}\zeta\)，得

\[
X(s,e^{-\Lambda s}\zeta)
\longrightarrow
\frac1\kappa\log(1+\kappa\zeta),
\]

\[
r(s,X(s,e^{-\Lambda s}\zeta))
\longrightarrow
\frac{\Lambda}{1+\kappa\zeta}.
\]

更明确地，令 \(E=e^{\Lambda s}\)，则

\[
X_E(\zeta):=X(s,\zeta/E)
=\int_0^\zeta
\frac{1+C(u/E)}{1+EC(u/E)}\,du
\longrightarrow
\frac1\kappa\log(1+\kappa\zeta)
\]

在每个有界 \(\zeta\) 区间上以 \(C^1\) 收敛。极限导数
\((1+\kappa\zeta)^{-1}\) 在这些区间上一致离开零，所以逆映射也
\(C^1\) 收敛。又因 (E36f) 给出 \(C\ge0\)，(E36e) 给出
\(X_\alpha\ge1\)，故流映射覆盖整个 \([0,\infty)\)。

消去 \(\zeta\) 后便得严格的局部收敛

\[
\boxed{
r(s,X)\longrightarrow r_*(X):=\Lambda e^{-\kappa X}
\quad\hbox{in }C^1_{\rm loc}([0,\infty)).
}
\tag{E36h}
\]

更高阶带余项控制的初值给出相应的 \(C^k_{\rm loc}\) 收敛。
极限的 directional Bessel exponent 不是预先拟合的常数，而是由初始
边界 jet 精确选出：

\[
\boxed{
q_\infty=-\kappa
=\frac{r_0'(0)}{\Lambda}<0.
}
\tag{E36i}
\]

由 \(v_X=-r\)、\(v(0)=0\)，极限流为

\[
v_*(X)=\frac{\Lambda}{q_\infty}
\bigl(1-e^{q_\infty X}\bigr).
\]

若取 \(b=\Lambda/q_\infty\)，重建的 \(h_*=v_*-b\) 正是
\(-b e^{q_\infty X}\)，因而极限就是 (E25)。一般 \(b\) 只使极限
多一个常数 \(Y\)-harmonic 项，它改变空间均匀的水平速度，并非
无物理效果的流函数常数 gauge。只有同时调制 \(b\)/中心运动时，
才能把这个均匀速度吸收。
因为 \(q_\infty\) 随初始 jet 改变，一般结论是对整个指数族的
orbital asymptotic stability；只有在固定 \(r_0'(0)=\Lambda q\) 的子空间中，
才是向某个预指定 (E25) profile 收敛。

若还有 \(r_0'<0\)，则 (E36d)--(E36e) 直接给出

\[
r_X(s,X)<0,
\qquad
R_X(s,X,Y)=Yr_X(s,X)<0
\quad(Y>0),
\tag{E36j}
\]

即用户要求的统一 \(\partial_XR\) 符号在整个稳定类中保持。
由 \(c-\gamma=(1-\beta)c\) 和 (E12)，还得到

\[
\boxed{
(T-t)\partial_y\rho
\longrightarrow
\frac{2-\beta}{1-\beta}e^{-\kappa X}
}
\tag{E36k}
\]

这里采用第 7 节的 \(L,\delta,a\)，并在物理点
\(x=a(t)+L(t)X\)、\(y=L(t)Y\) 上取极限；当 \(b\) 为常数时
\(a(t)\to a_*\)。(E36k) 在每个固定内坐标紧集上成立，但只是
内区局部梯度奇化；该 \(Y\)-线性族的全局范数从初始时就无界。
这样，(E25) 就不再是孤立的
显式抵消：它在这个右向、单调、\(Y\)-线性不变类中有严格的
非线性局部轨道吸引盆。对每个固定物质标号 \(\alpha>0\)，(E36e) 还给出
\(X(s,\alpha)\to\infty\)。因而收敛是辐射型/局部的，它从未声称
全局 Casimir 范数向单一 profile 收敛，所以不与无扩散守恒矛盾。

该吸引盆还有一个清楚的边界：若某处 \(r_0>\Lambda\)，则 \(C<0\)，
(E36d) 的分母在有限重整化时间消失，规则轨道退出上述稳定盆。
若 \(r_0'(0)=0\)，则 \(\kappa=0\)，指数极限退化，必须重新选择更高阶
缩放。

### 9.4 三尺度族的最大共同-\(\mu\) 模态类：严格中性

固定 (E31) 的 \(\mu(\tau)\)，并只在不显含 \(\tau\) 的
\(C^2\) 函数 \(F(I,J)\) 类中写

\[
\rho=F(I,J),
\qquad
\Psi=\Psi_H-\frac{F}{\mu(\tau)}.
\]

要求 closure 对所有 \(\tau\) 成立，等价于

\[
F_{II}=\alpha F_I,
\qquad
F_{JJ}=\frac{\Gamma^2}{\alpha}F_I.
\tag{E36l}
\]

全部解为

\[
F=e^{\alpha I}
\bigl(A_+e^{\Gamma J}+A_-e^{-\Gamma J}\bigr)
+B_0+B_1J.
\]

再施加 wall trace \(F(I,0)=0\)，只剩

\[
\boxed{
F=Ae^{\alpha I}\sinh(\Gamma J)+B_1J.
}
\tag{E36m}
\]

这个分类在“同一 \(\mu(\tau)\)、\(F\) 不显含 \(\tau\)”的类中是精确的。
\(A,B_1\) 在完整非线性演化中保持，因为
harmonic 流保持 \(I,J\)，而 Bessel self-velocity 与 \(\nabla F\) 正交。
因此 (E32) 在这个最大的共同-\(\mu\) exact modal 类中中性稳定，
但精确的 \(A\) 与 \(B_1\) 零模禁止它渐近吸引到单一 profile。

### 9.5 任意二维扰动仍是开放问题

对一般二维复合 fixed point，令 \(d=d_*+g\)。固定 \(m\) 和全部
调制参数时，线性化的正确广义谱 pencil 是

\[
\big[\sigma A_m+\mathcal L_*\big]g=0.
\tag{E37}
\]

若以真实密度扰动 \(\eta=A_mg\) 为状态，生成元是
\(-\mathcal L_*A_m^{-1}\)；在无界全/半平面的自然空间中，
\(A_m^{-1}\) 在低频无界。所以一般稳定性还必须依次建立：

1. source-neutral graph domain 上的闭生成元与局部适定性；
2. 投影掉平移、缩放模并处理无限组 Casimir 伴随约束后的
   generalized resolvent gap；
3. 有限维不稳定子空间、可逆调制矩阵和非线性 tame estimate。

只在这些前提下，一个有界、coercive 的 \(S=S^*>0\) 满足

\[
\mathcal G_*^*S+S\mathcal G_*\le-2\kappa S,
\qquad
\kappa>0,
\tag{E38}
\]

才可与 Lyapunov--Perron 构造一起给出 center-stable manifold。本文已经
证明 (E25) 在上述 \(Y\)-线性类中的非线性收敛，但没有证明 (E25)
或 (E32) 对任意二维扰动稳定。Bessel mass operator 的低频缺口只是
无界域自然空间的陈述；在有界穿孔域或显式 IR cutoff 下可恢复 coercivity。

## 10. 为什么标准 \(K_0/K_1\) 不能直接替换 directional mode

这里需要分清“可调制的运动”和真正障碍。均匀平移速度可以用
pole 中心的特征运动吸收；conformal dilation 也可部分进入振幅和
\(\mu\) 调制。不可由有限参数吸收的是 trace-free 双曲应变产生的
新同阶角模。

例如在 pole 附近，令 \(S=\operatorname{diag}(1,-1)\)。对

\[
P_{K_0}\sim-A\log r,
\]

有

\[
(Sx)\cdot\nabla P_{K_0}
\sim-A\cos(2\theta).
\tag{E38a}
\]

但 \(K_0\) 的振幅、\(\mu\) 与中心切向在这一阶只产生
\(\log r\)、径向/常数项或 \(r^{-1}\) 的 \(n=1\) 模，无法产生
order-one 的 \(n=2\) 模。同样，对与水平壁兼容的法向 \(K_1\) 主部

\[
P_{K_1}\sim\frac{\sin\theta}{r},
\]

直接计算得

\[
(Sx)\cdot\nabla
\left(\frac{\sin\theta}{r}\right)
=-\frac{\sin(3\theta)}r.
\tag{E38b}
\]

这个同阶 \(n=3\) 模也不在振幅、方向、中心和 \(\mu\) 的有限
切空间中。因此，在没有同阶角向 residual 时，标准局域
\(K_0/K_1\) 族不能在非退化双曲应变下由有限参数保持闭合。
directional exponential 能闭合，是因为它只含一个波矢方向，并可与应变
特征方向对齐。

因此 actual \(K\)-profile 的下一步不是有限参数调制，而是无限角模、order-one
\(D\)-层。它必须在 pole/core 附近以同阶抵消 strain，并在外部恢复 source
neutrality。标准 localized \(K\) 族的完整阶次/no-go 审计见
[BESSEL_MULTISCALE_BLOWUP_THEORY_ZH.md](BESSEL_MULTISCALE_BLOWUP_THEORY_ZH.md)。

## 11. 从精确机制到有限能水平壁爆破还缺什么

第 5 节证明了“有界 density front + \((T-t)^{-1}\) harmonic strain”能精确爆破；
第 7、8 节证明了 Bessel closure mode 与 order-one residual 能在水平壁精确共存。
剩下的物理 gluing 目标由此变得具体：

1. 用一个 source-neutral 的 Bessel/return 外层，在内区产生
   \(cXY+o(r^2)\) 的 harmonic 双曲应变；
2. 用同阶二维 wall corrector 把第 5 节非恒定的 wall streamfunction 改成
   Dirichlet 常数，同时保留压缩率；
3. 让 Casimir 通过开放的 return filament 输往重整化无穷远；
4. 在 inner 与 outer 两个空间上证明 block spectral gap，并控制中间 no-neck 区。

若 inner/outer 线性块在同一投影 energy 下分别 coercively
dissipative，gap 为 \(\kappa_{\rm in},\kappa_{\rm out}>0\)，outer resolvent
在 source-neutral 子空间可控，且耦合算子
\(\mathcal C_{\rm io},\mathcal C_{\rm oi}\) 有界，那么一个直接的 small-gain
充分条件是

\[
\|\mathcal C_{\rm io}\|\,
\|\mathcal C_{\rm oi}\|
<\kappa_{\rm in}\kappa_{\rm out}.
\tag{E39}
\]

等价地，需要 Schur complement

\[
\mathcal L_{\rm in}
-\mathcal C_{\rm io}\mathcal L_{\rm out}^{-1}
 \mathcal C_{\rm oi}
\]

在 inner 中性/压缩模上给出正确符号。

## 12. 最终判断

现在已经可以严格回答“Bessel 能否作为桥梁帮助发现新爆破机制”：

\[
\boxed{
\begin{gathered}
\text{能。order-one closure residual 通过 Bessel 漂移--Helmholtz Green 算子}\\
\text{产生可被 harmonic/composite strain 精确压缩的 density profile；}\\
\text{directional wall profile 族还在自然不变类的右侧扇区中}\\
\text{具有显式的非线性局部轨道收敛机制。}
\end{gathered}}
\tag{E40}
\]

已证明的三层结论是：

1. (E0)/(E18) 是允许不衰减 harmonic far-field contract 时的
   bounded-density pure-gradient blow-up；
2. (E25) 是水平壁上的精确 directional Bessel 内区奇化，并且
   (E36d)--(E36k) 证明该一参数族在右向单调 \(Y\)-线性类的
   \(X\ge0\) characteristic sector 中局部轨道吸引；
3. (E32) 是水平壁上的精确三尺度 Bessel 内区奇化，但其最大
   共同-\(\mu\) 模态类只有中性稳定，没有向单一 profile 的吸引。
4. (E56)--(E67) 从有界、水平紧支撑的 \(C^\infty\) 物理密度出发，
   精确产生局部 \(C^\infty\) 收敛到 \(e^{-X}\) directional Bessel 核的
   三长度爆破轨道；它依赖开放边界的无穷能 harmonic strain。
5. (E80)--(E90) 在固定、时间无关的双曲远场合同下，从水平紧支撑
   \(C^\infty\) 初值唯一演化并长时局部收敛到纯
   \((P_*,R_*)=(e^{-X},e^{-X})\)；它是无限时边缘 tangent，而非有限时爆破。
6. (E93)--(E110) 分别给出 wall-compatible 正弦收敛轨道、横竖切片紧支撑
   的旋转条带有限时 Bessel 爆破，并证明所有固定 transverse 单模/有限
   Fourier 修补都不足以产生全局有限能解。

这些解不与先前的 localized \(K_0/K_1\) no-go 矛盾。(E25) 在壁上
\(R=R_X=0\)，且另一水平端增长；(E32) 的 closure/harmonic strain
是主阶并在远场无界。它们都避开了“localized \(K\) 背景 + 全域小余留”
的前提。连续参数族也不选择某个数值 \(\beta\)；真正有限能全局问题的
\(\beta\) 仍必须由 source-neutral return/Fredholm 条件选出。

尚未完成的是：用 source-neutral、有限能的外层替换无穷能
harmonic strain，并同时满足水平壁、右向锐扇、wall jump 与全部
Casimir。这已经被缩减为一个明确的二维 Bessel/return gluing 与谱稳定问题。

## 13. Bessel 是否可以是光滑初值的含时收敛态？

答案取决于“收敛”的拓扑。现在可以严格分成三层：

1. 在固定物理坐标中，标准有限能无源解不能长时强收敛到
   非零 Bessel 稳态；
2. 缩放后的全局 \(L^p\cap L^q\) 强收敛也被 Casimir 排除；
3. 但在固定重整化紧集上的局部收敛确实可以发生；事实上，下面给出
   从有界、水平紧支撑的 \(C^\infty\) 密度出发并精确收敛到
   directional Bessel 核的完整轨道。

### 13.1 从全 \(X\) 逐点光滑约化初值出发的精确收敛轨道

沿用第 9 节的

\[
\Lambda=(2-\beta)c,
\qquad
P=cXY+Yh,
\qquad
R=Yr.
\]

取任意

\[
\kappa>0,
\qquad
0<A<1,
\qquad
k=\frac\kappa A,
\]

并定义全直线上的光滑有界约化初值

\[
\boxed{
C(\alpha)=A(e^{k\alpha}-1),
\qquad
r_0(\alpha)
=\frac{\Lambda}{1-A+Ae^{k\alpha}}.
}
\tag{E41}
\]

它满足

\[
r_0\in C_b^\infty(\mathbb R),
\qquad
r_0(0)=\Lambda,
\qquad
r_0'(0)=-\Lambda\kappa.
\tag{E42}
\]

令

\[
E=e^{\Lambda s},
\qquad
D(s,\alpha)=1+EC(\alpha),
\qquad
J(s,\alpha)=E\frac{1+C(\alpha)}{D(s,\alpha)}.
\tag{E43}
\]

当 \(EA\le1\) 时取 \(I_s=\mathbb R\)。当 \(EA>1\) 时取

\[
I_s=(\alpha_c(s),\infty),
\qquad
\alpha_c(s)=\frac1k
\log\!\left(1-\frac1{EA}\right).
\tag{E44}
\]

然后以特征参数形式定义

\[
\boxed{
X(s,\alpha)=\int_0^\alpha J(s,a)\,da,
\qquad
r(s,X(s,\alpha))=\frac{\Lambda}{D(s,\alpha)},
\qquad
v(s,X(s,\alpha))=-X_s(s,\alpha).
}
\tag{E45}
\]

对每个 \(s\ge0\)，\(\alpha\mapsto X(s,\alpha)\) 都是从 \(I_s\) 到
\(\mathbb R\) 的光滑严格递增微分同胚。事实上，

\[
(\log J)_s=\frac{\Lambda}{D}=r,
\]

因此

\[
X_s=-v,
\qquad
v_X=-\frac{X_{s\alpha}}{X_\alpha}=-r,
\qquad
\frac{dr}{ds}\bigg|_{\alpha}=r(r-\Lambda).
\tag{E46}
\]

这恰好就是 (E36b)。令 \(h=v-b\)，则 (E45) 给出完整
重整化 IPM 的精确解。当 \(s=0\) 时 \(J=1\)、\(X=\alpha\)，
所以它确实从 (E41) 的逐点 \(C^\infty\) 约化初值出发，而不是把极限
profile 直接当作初值。这里完整密度 \(R=Yr\) 仍随 \(Y\) 线性增长，
因此不能把“约化初值有界”误读成二维初值属于 \(C_b\)。

令 \(K\Subset(-1/\kappa,\infty)\) 为任意紧的重标号区间，并取
\(\alpha=\zeta/E\)、\(\zeta\in K\)。则

\[
X(s,\zeta/E)
\longrightarrow
\frac1\kappa\log(1+\kappa\zeta),
\qquad
r(s,X(s,\zeta/E))
\longrightarrow
\frac{\Lambda}{1+\kappa\zeta}.
\]

逆变换后得到全直线上的局部收敛

\[
\boxed{
r(s,\cdot)
\longrightarrow
\Lambda e^{-\kappa X}
\quad\hbox{in }C^1_{\rm loc}(\mathbb R).
}
\tag{E47}
\]

取

\[
q_\infty=-\kappa,
\qquad
b=\frac{\Lambda}{q_\infty}=-\frac{\Lambda}{\kappa},
\]

则极限的 \((P,R)\) 逐字等于 (E25)。映射
\(X=\kappa^{-1}\log(1+\kappa\zeta)\) 把
\((-1/\kappa,\infty)\) 覆盖到全直线。因此，在允许下面所述
远场合同的逐点光滑类中，directional Bessel profile 确实是一个
含时重整化收敛态。需要注意：\(s\ge s_c\) 后这个 Eulerian 分支只在每个
有限 \((s,X)\) 上光滑，左远场非有界且特征流不完备；它不是标准的
全局 classical \(C_b\) Cauchy 轨道，详见第 13.4 节。

### 13.2 更强的正例：有界、水平衰减密度向 Bessel 模精确收敛

上一例的物理密度因 \(Y\)-线性因子而无界。其实还可以构造
一个密度全局有界、水平衰减、初始梯度有限的精确收敛轨道。

取常数

\[
0<g<c,
\qquad
f(\xi)=\operatorname{sech}\xi,
\]

并定义

\[
L(s)=e^{-cs},
\qquad
\varepsilon(s)=gs,
\qquad
\delta(s)=\operatorname{sech}(gs),
\]

\[
b=g,
\qquad
\Gamma(s):=-\frac{\delta_s}{\delta}
=g\tanh(gs).
\tag{E48}
\]

令

\[
\boxed{
R(s,X,Y)=\frac{f(X+\varepsilon(s))}{\delta(s)},
\qquad
p(s,X)=\frac1{\delta(s)}
\int_{X+\varepsilon(s)}^\infty f(\eta)\,d\eta,
\qquad
P=cXY+p.
}
\tag{E49}
\]

由定义 \(p_X=-R\)，因此

\[
-\Delta P=-p_{XX}=R_X.
\]

另一方面，

\[
U=\nabla^\perp P=(-cX,cY-R),
\]

所以完整重整化输运速度为

\[
U+cz-be_1=(-g,2cY-R).
\]

由于 \(R\) 不依赖 \(Y\)，而且

\[
R_s=gR_X+\Gamma(s)R,
\]

完整输运方程逐项变成

\[
\boxed{
R_s+(U+cz-be_1)\cdot\nabla R
=R_s-gR_X
=\Gamma(s)R.
}
\tag{E50}
\]

因此 (E49) 是非自治 dynamic scaling 下的完整精确 IPM 轨道。
由 \(\operatorname{sech}\xi\sim2e^{-\xi}\) 与
\(\int_\xi^\infty\operatorname{sech}\eta\,d\eta\sim2e^{-\xi}\)，有

\[
\boxed{
(R(s),p(s),\Gamma(s))
\longrightarrow
(e^{-X},e^{-X},g)
\quad\hbox{in }C^k_{\rm loc}
\quad\hbox{for every }k.
}
\tag{E51}
\]

这里使用了辐射规范 \(p(+\infty)=0\)。收敛仅为固定 \(X\)-紧集上的
局部收敛，并非全局 \(L^\infty\) 收敛；事实上

\[
\sup_X R(s,X)=\cosh(gs)\longrightarrow\infty.
\]

有限 \(s\) 时 \(\Gamma(s)\) 不是常数，而极限 scaling exponent 为

\[
\beta_*:=\lim_{s\to\infty}\frac{\Gamma(s)}c=\frac gc\in(0,1).
\]

极限 Bessel 部分 \(P_B=e^{-X}\)、\(R_B=e^{-X}\) 满足

\[
(\Delta+\partial_X)P_B=0,
\qquad
R_B=P_B.
\]

即带符号 drift 参数 \(m=-1\) 的 directional Bessel mode。完整极限的
closure residual 为

\[
D_*=R_*+mP_*=-cXY,
\tag{E52}
\]

所以它再次是“Bessel mode + order-one harmonic residual”的复合态。

物理时间由

\[
\frac{dt}{ds}=\frac L\delta=e^{-cs}\cosh(gs)
\]

给出，因此收缩在有限时刻 \(T\) 结束，并且

\[
T-t(s)=\frac12\left[
\frac{e^{-(c-g)s}}{c-g}
+\frac{e^{-(c+g)s}}{c+g}
\right].
\tag{E53}
\]

取 \(a_s=bL=gL\)，即

\[
a(s)=a_*-\frac gcL(s),
\]

物理密度恰好是

\[
\boxed{
\rho(t,x,y)-\rho_c
=\operatorname{sech}\!\left(
\frac{x-a(t)}{L(t)}+gs(t)
\right).
}
\tag{E54}
\]

它对每个 \(t<T\) 都属于 \(C_b^\infty\)，且在 \(x\) 方向衰减，

\[
\|\rho(t)-\rho_c\|_\infty=1,
\qquad
\|\rho_x(t)\|_\infty=\frac1{2L(t)}.
\]

在每个固定内坐标 \(X\) 上，

\[
\boxed{
(T-t)\rho_x(t,a(t)+L(t)X,y)
\longrightarrow
-\frac{e^{-X}}{c-g}.
}
\tag{E55}
\]

但全局梯度最大值位于 \(X=-gs+O(1)\)，而且更精确地

\[
\|\rho_x(t)\|_\infty
\sim\frac12\,[2(c-g)(T-t)]^{-c/(c-g)}.
\tag{E55a}
\]

所以它比 \((T-t)^{-1}\) 更快发散。这是一个同时
含有对数漂移内核和 Bessel 局部极限的多尺度轨道。

这个例子已经肯定回答：Bessel 可以是有界光滑密度初值的
含时重整化收敛态。但它的密度不依赖 \(y\)，调和项 \(cXY\)
在无穷远不衰减，因而横向质量与速度能量无穷；另外
\(P(X,0)=p(X)\) 不是常数，所以它也不满足水平不可穿壁。
它是全平面/开放边界 far-field contract 下的精确解，不是标准
有限能半平面 Cauchy 轨道。

### 13.3 更强的正例：水平紧支撑光滑密度向 Bessel 核精确收敛

上一例已经给出水平衰减；现在可以把它加强为每个 \(t<T\) 都在
\(x\) 方向紧支撑。取常数

\[
c>g>0,
\]

并选 \(f\in C_c^\infty(\mathbb R)\)、\(f\ge0\)。要求存在 \(r_0>0\)，使

\[
f(\xi)=\exp\!\left[-\frac1{1-\xi}\right]
\quad(1-r_0<\xi<1),
\qquad
f(\xi)=0\quad(\xi\ge1),
\tag{E56}
\]

而在更左侧作任意光滑紧支撑截断。固定

\[
s_0>\max\!\left\{\frac2c,\frac1{gr_0}\right\},
\]

并只考虑 \(s\ge s_0\)。定义

\[
d=\frac1{gs},
\qquad
\varepsilon=1-d,
\qquad
w=d^2,
\qquad
\delta=e^{-gs},
\qquad
L=e^{-cs},
\tag{E57}
\]

\[
A(s)=c+\frac{w_s}{w}=c-\frac2s>0,
\qquad
b=\frac{\varepsilon_s}{w}=g,
\qquad
\Gamma=-\frac{\delta_s}{\delta}=g.
\tag{E58}
\]

令

\[
\boxed{
R(s,X,Y)=\delta^{-1}f(\varepsilon+wX),
\qquad
p(s,X)=\frac1{\delta w}
\int_{\varepsilon+wX}^{1}f(\eta)\,d\eta,
\qquad
P=A(s)XY+p(s,X).
}
\tag{E59}
\]

由 \(p_X=-R\) 立刻得到 \(-\Delta P=R_X\)。同时

\[
U=\nabla^\perp P=(-AX,AY-R),
\qquad
(U+cz-be_1)_X=\frac{2X}{s}-g,
\]

而

\[
R_s=gR+\left(g-\frac{2X}{s}\right)R_X.
\]

故没有舍去任何项地得到

\[
\boxed{
R_s+(U+cz-be_1)\cdot\nabla R=gR,
\qquad
-\Delta P=R_X.
}
\tag{E60}
\]

特别地，\(A_s\) 不会产生额外项，因为 IPM 的瞬时椭圆速度恢复式和
输运方程中本来就没有 \(P_s\)。

对每个固定的 \(X\)-紧集，充分大的 \(s\) 使
\(\varepsilon+wX\) 落在 (E56) 的精确右边缘区。此时

\[
R(s,X)=\exp\!\left[-\frac{X}{1-dX}\right]
\longrightarrow e^{-X}
\quad\hbox{in }C^\infty_{\rm loc},
\tag{E61}
\]

并且由所选规范

\[
p(s,X)=\int_X^{1/d}R(s,Z)\,dZ
\longrightarrow e^{-X}
\quad\hbox{in }C^\infty_{\rm loc}.
\tag{E62}
\]

所以这是从普通光滑、水平紧支撑密度出发的完整精确收敛：

\[
\boxed{
(P,R)\longrightarrow(P_*,R_*)
=(cXY+e^{-X},e^{-X})
\quad\hbox{in }C^\infty_{\rm loc}.
}
\tag{E63}
\]

其 Bessel 部分对应 \(m=-1\)，即

\[
(\Delta+\partial_X)e^{-X}=0,
\qquad
D_*=R_*-P_*=-cXY.
\]

而且 \(R_{*,X}=-e^{-X}<0\)；因此在每个固定内区紧集上，充分大时间后
\(R_X\) 都有统一负号。由于非零紧支撑函数的导数必须在外部 return
区域换号，这个符号结论不能扩成全空间结论。

物理重构同样完全显式。取

\[
t(s)=T-\frac{e^{-(c-g)s}}{c-g},
\qquad
a(s)=a_*-\frac gc e^{-cs},
\qquad
L_0(s)=\frac Lw=g^2s^2e^{-cs},
\tag{E64}
\]

并记

\[
\xi=\varepsilon+\frac{x-a}{L_0},
\qquad
H(\xi)=\int_\xi^1f(\eta)\,d\eta,
\qquad
\varkappa(t)=A(s)\frac\delta L.
\]

则完整物理解为

\[
\boxed{
\rho(t,x,y)=\rho_c+f(\xi),
\qquad
\Psi(t,x,y)=\varkappa(t)(x-a)y+L_0(t)H(\xi).
}
\tag{E65}
\]

直接计算给出

\[
u=(-\varkappa(x-a),\varkappa y-f(\xi)),
\qquad
-\Delta\Psi=\frac{f'(\xi)}{L_0}=\rho_x,
\]

以及

\[
\frac{\dot L_0}{L_0}=-\varkappa,
\qquad
\dot a=L_0\dot\varepsilon=g\delta.
\]

故 \(\xi_t=\varkappa(x-a)/L_0\)、\(u_1\xi_x=-\varkappa(x-a)/L_0\)，
从而 \(\rho_t+u\cdot\nabla\rho=0\) 逐点成立。固定 Bessel 内坐标 \(X\) 时

\[
\boxed{
(T-t)\rho_x(t,a(t)+L(t)X,y)
\longrightarrow-\frac{e^{-X}}{c-g}.
}
\tag{E66}
\]

另一方面，全局最大梯度位于外部 return 尺度而非固定 \(X\) 内核：

\[
\|\rho_x(t)\|_\infty
=\frac{\|f'\|_\infty}{L_0}
\asymp\frac{e^{cs}}{s^2}.
\tag{E67}
\]

这里同时出现三种可辨认长度：Bessel 核尺度 \(L\)、核到右边缘的距离
\(Lgs\)，以及装载 Casimir 的外部支撑尺度 \(L/w=L_0\)。紧支撑 return
层在重整化坐标中逃向无穷远，正因如此局部 Bessel 收敛不违反全局输运不变量。

这个构造给出了目前最强的肯定答案：\(\rho(t)-\rho_c\) 对每个
\(t<T\) 都有界、\(C^\infty\) 且在 \(x\) 方向紧支撑，包括任意固定
\(s_0\) 时刻的初值。但它仍不依赖 \(y\)，所以二维 \(L^p\) 质量无穷；
\(P(X,0)=p(X)\) 非常数，故不满足水平不可穿壁；\(A(s)XY\) 还给出
线性、随时间增强的无穷能远场应变。它因此是开放边界/外加 harmonic
far-field contract 下的精确梯度爆破，而不是标准衰减 Biot--Savart
规范下的有限能 Cauchy 爆破。

### 13.4 第一个光滑收敛例子为什么仍非标准有界 Cauchy 轨道

定义

\[
s_c=\frac1\Lambda\log\frac1A.
\tag{E68}
\]

在 \(s<s_c\) 时，\(r(s,\cdot)\) 全局有界。当 \(s\uparrow s_c\) 时，
左端 \(X\to-\infty\) 的范数失去一致有界性。对 \(s\ge s_c\)，
\(r(s,X)\) 在每个有限 \((s,X)\) 上仍是 \(C^\infty\)，但不再属于
\(C_b\)；一部分左侧特征在有限重整化时间逃到 \(-\infty\)。

这不是 (E41) 的偶然缺陷。指数极限要求

\[
r_0'(0)=-\Lambda\kappa<0.
\]

因此任意两侧 \(C^1\) 延拓都在某个 \(X<0\) 处有
\(r_0>\Lambda\)。沿完整特征，

\[
r(s)=
\frac{\Lambda r_0}
{e^{\Lambda s}\Lambda-(e^{\Lambda s}-1)r_0}
\]

会在有限 \(s\) 使分母消失。若假设 \(r\) 全时间属于 \(L^\infty\)，
则 \(v_X=-r\) 使速度全局 Lipschitz、特征完备，这与上述有限时
logistic 发散矛盾。所以在精确 \(Y\)-线性类中，不存在既两侧
光滑、又全时间保持 \(C_b\) 的指数吸引轨道。

进一步，对任意初始尺度 \(L_0,δ_0>0\)，物理初值为

\[
\rho_0-\rho_c
=\delta_0\frac y{L_0}
r_0\!\left(\frac{x-a_0}{L_0}\right).
\tag{E69}
\]

它逐点 \(C^\infty\)，但因 \(r_0(0)=\Lambda\ne0\)，在无界半平面上
密度随 \(y\) 线性增长；又因 \(r_0'\not\equiv0\)，\(\rho_{0x}\) 也随
\(y\) 增长。所以 (E41) 回答了“光滑初值”问题，但没有回答
“全局有界、衰减、有限能初值”问题。

### 13.5 固定物理坐标中的长时 Bessel 强吸引被势能排除

令 \(\vartheta=\rho-\rho_\infty\)，其中常数背景已被减去。假设解对
全部 \(t\ge0\) 存在，扰动势能有限且有统一下界，\(u\in L^2\)，并且
边界与无穷远没有通量。对这样的标准无源、无穿壁半平面解，有

\[
\frac d{dt}\int_\Omega y\vartheta(t,x,y)\,dx\,dy
=-\int_\Omega|u(t,x,y)|^2\,dx\,dy.
\tag{E70}
\]

若左边势能有有限下界，则

\[
\int_0^\infty\|u(t)\|_2^2\,dt<\infty.
\]

因此，若在某个固定物理紧集 \(K\) 上，当 \(t\to\infty\) 时
\(u(t)\to u_B\) strongly in \(L^2(K)\)，则必须 \(u_B|_K=0\)。
对全局强极限，并假设密度/流函数的收敛足以通过椭圆关系取极限，这给出

\[
u_B=0,
\qquad
\partial_x\rho_B=0.
\tag{E71}
\]

再加水平衰减，只剩零或纯分层态，不可能是非零 Bessel 稳态。这个论证
只排除 \(t\to\infty\) 的固定物理坐标强吸引，不排除有限时动态重标度局部极限。
同一事实也可由静态椭圆方程直接看出。若

\[
(\Delta-m\partial_X)P_B=0,
\qquad
P_B\in H_0^1(\Omega)
\]

为实函数，并且具有足以令水平漂移边界项消失的衰减/截断，则

\[
0=\int_\Omega P_B(\Delta-m\partial_X)P_B
=-\int_\Omega|\nabla P_B|^2,
\tag{E72}
\]

所以 \(P_B=0\)。非零 \(K_0/K_1\) 或 directional mode 必须通过
pole source、边界/无穷远通量或无穷能避开 (E72)。

### 13.6 全局重整化强收敛被两个 Casimir 同时排除

令 \(\vartheta=\rho-\rho_c\)，并假设

\[
\vartheta(t,x)
=\delta(t)R_t\!\left(\frac{x-a(t)e_1}{L(t)}\right).
\]

假设 \(\delta,L>0\)。若对两个不同的
\(1\le p\ne q<\infty\)，初始的两个 Casimir 都有限且非零，并且

\[
R_t\longrightarrow R_B\ne0
\quad\hbox{globally strongly in }L^p\cap L^q,
\tag{E73}
\]

则不可压输运的 Casimir 守恒给出

\[
\delta^pL^2\|R_t\|_p^p=C_p,
\qquad
\delta^qL^2\|R_t\|_q^q=C_q.
\tag{E74}
\]

其中 \(R_B\) 的两个范数都非零且有限，则两式强制

\[
\delta(t)\to\delta_*>0,
\qquad
L(t)\to L_*>0.
\tag{E75}
\]

所以任何 \(L\to0\) 或 \(\delta\to0\) 的爆破都不可能由单一 Bessel
profile 的全局 \(L^p\cap L^q\) 强收敛描述。所有 level-set 面积必须
留在外层 reservoir/return filament 中。

### 13.7 \(K_0/K_1\) 还有额外的 source-moment 障碍

对每个 source-neutral 光滑解，若完整圆盘 \(B_r\) 严格包含于定义域
（或已作全平面延拓），则包围 tip 的轮廓满足

\[
-\oint_{\partial B_r}\partial_nP
=\oint_{\partial B_r}R\,n_X,
\tag{E76}
\]

半圆上还有

\[
\oint_{\partial B_r^+}
(Y\partial_nP-Pn_Y)
=-\oint_{\partial B_r^+}YR\,n_X.
\tag{E77}
\]

纯 \(K_0\) trace 在完整圆的 (E76) 中留下非零 monopole。对满足
Dirichlet 壁条件的法向 \(K_1\) 模，(E77) 在半圆中留下非零 normal
dipole；纯 \(K_0\) 本身已经不满足该壁条件。一般壁上半圆若没有这些
条件会多出壁边界项，不能直接套用 (E76)。因此 source-neutral 解不能
在相应圆/半圆上，以能控制 \(P,\partial_nP,R\) trace 的拓扑收敛到纯
\(K_0/K_1\)。严格能推出的是：每条包围轮廓上至少一个 trace 必须有
order-one 修正；它可以是弥散修正或内部补偿场，恒等式本身并没有证明
一条几何上局域的 layer 必须穿过每条轮廓。

### 13.8 若保留 Bessel core，标准有限能候选所需的复合形式

上述 no-go 不排除下面的复合局部收敛；反过来，若坚持保留一个局部
Bessel core，则这些成分是由现有恒等式所要求的必要结构：

\[
\boxed{
\text{local rescaled Bessel response}
+\text{order-one }D\text{/wall layer}
+\text{source-neutral return}
+\text{outer Casimir reservoir}.
}
\tag{E78}
\]

它只能在固定重整化紧集上收敛。对 \(K_0/K_1\)，紧集还必须
避开 tip 与 closure ray；对 directional mode，指数增长端必须在
进入物理远场前由二维 return 层截断。外层同时承担全部 Casimir、
source moments 和势能通量。

目前尚未证明这个有限能复合态存在，也没有排除完全不同的非 Bessel
机制。因此对用户的新问题，最精确的答案是：

\[
\boxed{
\begin{gathered}
\text{Bessel 确实可以是光滑初值的局部重整化含时收敛态；}\\
\text{但它不可能是标准有限能无源动力学的全局强吸引子。}\\
\text{若标准有限能候选保留 Bessel core，则必须加入 order-one return/outer reservoir。}
\end{gathered}
}
\tag{E79}
\]

## 14. 固定远场合同下，Bessel 确实是光滑初值的长时收敛态

第 13.2--13.3 节的有限时收缩依赖随时间增强的 affine harmonic strain。
如果问题只问“Bessel 能否由某个光滑初值在含时演化中出现”，则可以进一步
把远场合同固定为时间无关，并得到一个更接近通常 Cauchy 含义的肯定答案。

### 14.1 时间无关双曲应变下的精确紧支撑轨道

取 \(\sigma>0\)、固定点 \(x_*\)，以及 \(F\in C_c^\infty(\mathbb R)\)、
\(F\ge0\)。假设存在 \(r_0>0\)，使它在右支撑端点满足

\[
F(\xi)=e^{1/\xi}\quad(-r_0<\xi<0),
\qquad
F(\xi)=0\quad(\xi\ge0),
\qquad
H(\xi)=\int_\xi^0F(\eta)\,d\eta.
\tag{E80}
\]

再固定充分大的 \(q_0>0\)，并令

\[
q(t)=q_0+\sigma t,
\qquad
\ell(t)=e^{-q(t)},
\qquad
\xi=\frac{x-x_*}{\ell(t)},
\]

并定义

\[
\boxed{
\rho(t,x,y)=\rho_c+F(\xi),
\qquad
\Psi(t,x,y)=\sigma(x-x_*)y+\ell(t)H(\xi).
}
\tag{E81}
\]

这里的远场合同

\[
\Psi-\sigma(x-x_*)y\in L^\infty(\mathbb R^2),
\qquad
\lim_{x\to+\infty}
[\Psi-\sigma(x-x_*)y]=0
\tag{E82}
\]

完全不随时间改变。两个满足同一 Poisson 方程与 (E82) 的流函数之差是
有界整调和函数；由 Liouville 定理及右端规范只能为零。因此一旦固定
\(\sigma\) 和这个远场合同，速度由当前密度唯一恢复，而不是逐时任意选择。
在 \(y\)-平移不变 classical 类中，一维线性特征还给出该轨道的唯一性；
若要声称完整二维非衰减函数类的 well-posedness，则仍需另行指定并研究
相应 Banach 空间。

直接计算

\[
u=(-\sigma(x-x_*),\sigma y-F(\xi)),
\qquad
-\Delta\Psi=\frac{F'(\xi)}\ell=\rho_x,
\]

对应 Darcy 压力甚至可以取成时间无关的

\[
p_{\rm Darcy}
=\frac\sigma2\big[(x-x_*)^2-y^2\big]-\rho_cy.
\]

以及 \(\xi_t=\sigma\xi\)、\(u_1\xi_x=-\sigma\xi\)。故

\[
\boxed{
\rho_t+u\cdot\nabla\rho=0,
\qquad
-\Delta\Psi=\rho_x
}
\tag{E83}
\]

对所有有限 \(t\ge0\) 逐点严格成立。

现在令

\[
d=q^{-1},
\qquad
\delta=\ell=e^{-q},
\qquad
a=x_*-\ell d,
\qquad
L=\ell d^2.
\tag{E84}
\]

重整化时间仍按 \(ds/dt=\delta/L\) 定义。于是

\[
\frac{ds}{dt}=q^2,
\qquad
s-s_0=\frac1{3\sigma}(q^3-q_0^3).
\]

在 \(x=a+LX\)、\(y=LY\) 上有

\[
\xi=-d+d^2X.
\]

对任意固定 \(X\)-紧集，充分大时间后落入 (E80) 的精确平坦边缘区，因而

\[
\boxed{
R(t,X)=\delta^{-1}[\rho(t,a+LX,LY)-\rho_c]
=\exp\!\left[-\frac{X}{1-dX}\right]
\longrightarrow e^{-X}
}
\tag{E85}
\]

于 \(C^\infty_{\rm loc}\) 成立。相应流函数恰好为

\[
P(t,X,Y)=\frac{\Psi(t,a+LX,LY)}{\delta L}
=\sigma(-d+d^2X)Y+p(t,X),
\]

\[
p(t,X)=\int_X^{1/d}R(t,Z)\,dZ
\longrightarrow e^{-X}
\quad\hbox{in }C^\infty_{\rm loc}.
\tag{E86}
\]

所以

\[
\boxed{
(P,R)\longrightarrow(P_*,R_*)=(e^{-X},e^{-X})
\quad\hbox{in }C^\infty_{\rm loc}.
}
\tag{E87}
\]

这一次极限没有 order-one harmonic residual：取 \(m=-1\)，则

\[
(\Delta+\partial_X)P_*=0,
\qquad
R_*=P_*,
\qquad
D_*=R_*-P_*=0.
\]

原因是所选 Bessel 子尺度比主支撑尺度小两个对数因子；固定的物理应变在
该子尺度中消失。事实上相应的重整化率为

\[
c=\sigma d^2(1+2d),
\qquad
\gamma=\sigma d^2,
\qquad
b=\sigma(d+d^2),
\tag{E88}
\]

它们都趋于零，而 \(\gamma/c\to1\)。这不是第 3 节的有限时
\(c_*-\gamma_*>0\) 情形，而是一个无限物理时间的纯 Bessel tangent state。

它仍有三种长度

\[
L=\frac\ell{q^2},
\qquad
Lq=\frac\ell q,
\qquad
Lq^2=\ell.
\tag{E89}
\]

固定 Bessel 内区与全局最大梯度分别满足

\[
\frac1{q^2}\,
\rho_x(t,a(t)+L(t)X,y)
\longrightarrow-e^{-X},
\qquad
\|\rho_x(t)\|_\infty=e^{q(t)}\|F'\|_\infty.
\tag{E90}
\]

因此所有有限时间的密度都光滑、有界并在 \(x\) 方向紧支撑，但梯度只在
\(t\to\infty\) 时无界。这个极限是精确轨道/tangent limit；它没有证明
Bessel 对邻近初值稳定，也不能被称为有限时爆破吸引子。

### 14.2 为什么有限时例子必须把远场协议算作数据的一部分

更一般地，给任意正光滑函数 \(\lambda(t)\)，令

\[
\kappa(t)=\frac{\dot\lambda(t)}{\lambda(t)},
\qquad
\xi=\lambda(t)(x-x_*),
\]

则

\[
\boxed{
\rho_\lambda=\rho_c+F(\xi),
\qquad
\Psi_\lambda=\kappa(t)(x-x_*)y
+\lambda(t)^{-1}H(\xi)
}
\tag{E91}
\]

对每个 \(\lambda\) 都是逐点精确 IPM 解。若 \(\lambda_1,\lambda_2\) 在
\(t=0\) 具有相同的值与一阶导数、但以后不同，则两条轨道具有完全相同的
\((\rho_0,\Psi_0,u_0)\)，未来却不同。因此只允许 \(O(r^2)\) 调和增长、
却不指定其 leading coefficient 的“全平面问题”甚至不是唯一的 Cauchy
问题；完整的 \(\kappa(t)\) 协议也是数据。

标准衰减/Riesz--Biot--Savart 规范排除这个调和分支。在当前
\(y\)-不变类中它选出 \(u_1=0\)，所以密度 profile 冻结，不会产生上述
收缩。第 13.2--13.3 节因此应严格称为

\[
\boxed{
\text{无 bulk source、但带预先指定的含时 affine pressure/far-field protocol 的}\\
\text{完整 classical IPM 梯度爆破解；不是 density-only 的标准无外场 Cauchy 解。}
}
\tag{E92}
\]

第 14.1 节则只需固定的时间无关 contract，因而确实回答了“某个光滑初值
能否产生 Bessel 含时收敛态”；它仍因 \(y\)-不衰减和 affine strain 而具有
无穷横向质量与无穷能。标准有限能、source-neutral 二维初值能否产生同样的
局部 Bessel tangent/attractor，仍是第 13.8 节的 gluing 开放问题。

## 15. 向壁条件与横向局域化继续推进

第 13--14 节的有界正例不满足水平壁条件，而且最简单的例子不依赖 \(y\)。
下面先给两个还能保持精确闭合的改进，再给出解释“为什么仍未达到有限能”
的严格低阶障碍。

### 15.1 满足水平壁的正弦正则化轨道

取

\[
c>0,\qquad q<0,\qquad 0\le\gamma<c,\qquad
K=2c-\gamma,\qquad b=\frac Kq,
\]

以及 \(0<k_0<|q|\)、\(k(s)=k_0e^{-2cs}\)。定义

\[
\phi_k(Y)=\frac{\sin(kY)}k,
\qquad
\boxed{
R=Ke^{qX}\phi_k(Y),
\qquad
P=cXY-\frac{q}{q^2-k^2}R.
}
\tag{E93}
\]

因为

\[
R_X=qR,\qquad
\Delta R=(q^2-k^2)R,
\]

所以 \(-\Delta P=R_X\)。Bessel 部分
\(P_B=-qR/(q^2-k^2)\) 与 \(R\) 成比例，故
\(\nabla^\perp P_B\cdot\nabla R=0\)。另一方面

\[
R_s=2cR-2cYR_Y,
\qquad
\nabla^\perp(cXY)+cz-be_1=(-b,2cY).
\]

于是

\[
\boxed{
R_s+(U+cz-be_1)\cdot\nabla R
=(2c-bq)R=\gamma R,
\qquad
-\Delta P=R_X.
}
\tag{E94}
\]

而且

\[
P(X,0)=R(X,0)=P_X(X,0)=0,
\]

所以不可穿壁严格成立。随着 \(s\to\infty\)，

\[
\boxed{
(P,R)\longrightarrow
\left(cXY-\frac KqYe^{qX},\,KYe^{qX}\right)
=(cXY-bYe^{qX},\,KYe^{qX})
}
\tag{E95}
\]

于 \(C^\infty_{\rm loc}\) 收敛；这正是 \(K=(2-\beta)c\)、
\(\gamma=\beta c\) 时的 wall directional Bessel profile。有限 \(s\)
本身也有精确 Bessel closure：

\[
m(s)=\frac{q^2-k^2}{q},
\qquad
R=-mP_B,
\qquad
(\Delta-m\partial_X)P_B=0,
\qquad
m(s)\to q.
\tag{E96}
\]

这个轨道从每个有限 \(s\) 的逐点光滑数据出发，而且在 \(Y\) 方向有界。
但 \(q<0\) 强制 \(X\to-\infty\) 时指数增长；在有限时收缩
\(L=e^{-cs}\)、\(\delta=e^{-\gamma s}\) 下，其固定 \(X\) 的全 \(Y\)
振幅还按

\[
\sup_Y|\delta R|
=\frac K{k_0}e^{(2c-\gamma)s}e^{qX}
\]

增长。因此它是右扇区/开放入流边界上的精确 wall core，而不是整个半平面
上的 \(C_b\) 密度 Cauchy 轨道。

这个缺陷在自然单模类中无法修补。若

\[
P=cXY+h(s,X)\phi_k(Y),
\qquad
R=r(s,X)\phi_k(Y),
\]

则完整方程严格蕴含

\[
r_X=k^2h-h_{XX},
\qquad
r_s-br_X+(2c-\gamma)r
+(h_Xr-hr_X)\cos(kY)=0.
\tag{E97}
\]

因此 \(h_Xr-hr_X=0\)。非零 directional 分支上
\(r=\lambda(s)h\)，从而

\[
h_{XX}+\lambda h_X-k^2h=0.
\]

两个特征根一正一负，所以全 \(X\) 有界解只能 \(h=0\)，继而
\(r_X=0\)。有限个正弦法向模也可从最高 \(2N\) 次谐波开始逐级作同一
Wronskian 归纳。故真正全局有界的 wall 轨道若存在，必须使用无限法向模
和 order-one wall/return residual。

还有一个不依赖变量分离的瞬时版本：若 homogeneous Bessel 部分满足

\[
(\Delta-m\partial_X)P_B=0
\quad\hbox{in }\mathbb H,\qquad
P_B|_{Y=0}=0,\qquad
P_B\in L^\infty(\mathbb H),
\]

则常系数漂移椭圆方程的有界 Dirichlet 唯一性给出 \(P_B\equiv0\)。
所以非零 wall Bessel 本来就只能作为局部极限，并伴随全局无界/非紧致
行为或 order-one residual。

### 15.2 有界且每条横竖切片紧支撑的旋转条带爆破

现在回到第 13.3 节的 \(c>g>0\)、平坦边缘
\(f\in C_c^\infty\) 与同一个充分大的 \(s_0\)，并沿用

\[
d=(gs)^{-1},\qquad \varepsilon=1-d,\qquad
w=d^2,\qquad
r(s,\xi)=e^{gs}f(\varepsilon+w\xi).
\]

则

\[
r_s=gr+\left(g-\frac{2\xi}{s}\right)r_\xi.
\]

取任意 \(k_0\ne0\)，并令

\[
k(s)=k_0e^{-2cs},\qquad
\xi=X-kY,\qquad
\eta=kX+Y,\qquad
Q=1+k^2,
\tag{E98}
\]

并定义

\[
\mathcal B=\frac{k_s}{Q},
\qquad
\mathcal A=c+\frac{kk_s}{Q}-\frac2s,
\]

\[
P_H=\frac{\mathcal A}{Q}\xi\eta
+\frac{\mathcal B}{2Q}(\xi^2-\eta^2),
\qquad
h(s,\xi)=\frac1Q\int_\xi^{gs}r(s,z)\,dz,
\]

\[
\boxed{
R(s,X,Y)=r(s,\xi),
\qquad
P(s,X,Y)=P_H(s,X,Y)+h(s,\xi).
}
\tag{E99}
\]

\(P_H\) 是调和二次多项式，而 \(h_\xi=-r/Q\)。由于
\(|\nabla\xi|^2=Q\)，立即得到

\[
-\Delta P=R_X.
\]

直接计算还给出

\[
\left[
\partial_s+(\nabla^\perp P+cz-ge_1)\cdot\nabla
\right]\xi
=\frac{2\xi}{s}-g.
\]

\(\nabla^\perp h\) 与 \(\nabla\xi\) 正交，所以

\[
\boxed{
R_s+(\nabla^\perp P+cz-ge_1)\cdot\nabla R=gR,
\qquad
-\Delta P=R_X.
}
\tag{E100}
\]

当 \(s\to\infty\) 时 \(k\to0\)、\(P_H\to cXY\)，第 13.3 节的平坦边缘
计算逐字给出

\[
\boxed{
(P,R)\longrightarrow(cXY+e^{-X},e^{-X})
\quad\hbox{in }C^\infty_{\rm loc}.
}
\tag{E101}
\]

物理重构为

\[
L=e^{-cs},\qquad \delta=e^{-gs},\qquad
t=T-\frac{e^{-(c-g)s}}{c-g},\qquad
a=a_*-\frac gc e^{-cs},
\]

\[
\boxed{
\rho-\rho_c
=f\!\left(
\varepsilon+w\,\frac{x-a-k(s)y}{L}
\right),
\qquad
\Psi=\delta LP.
}
\tag{E102}
\]

所以 \(\rho\) 对每个 \(t<T\) 都全局有界、\(C^\infty\)，并且固定 \(x\)
时在 \(y\) 中紧支撑、固定 \(y\) 时也在 \(x\) 中紧支撑。例如

\[
\int_{\mathbb R}|\rho-\rho_c|^p\,dy
=\frac{L}{w|k|}
\int_{\mathbb R}|f(\zeta)|^p\,d\zeta<\infty.
\tag{E103}
\]

固定 Bessel 内区与全局最大梯度分别为

\[
\boxed{
(T-t)\rho_x(t,a+LX,LY)
\longrightarrow-\frac{e^{-X}}{c-g},
\qquad
(T-t)\rho_y(t,a+LX,LY)\longrightarrow0,
}
\]

\[
\|\nabla\rho(t)\|_\infty
=\frac{w}{L}\sqrt{1+k^2}\,\|f'\|_\infty
\sim\frac{e^{cs}}{g^2s^2}\|f'\|_\infty.
\tag{E104}
\]

但支撑是沿 \(\eta\) 方向无限延伸的斜条带，所以
\(\rho-\rho_c\notin L^p(\mathbb R^2)\)，且调和应变仍使
\(\nabla\Psi\notin L^2\)。这个例子把“完全不依赖 \(y\)”加强成“每条坐标
切片都有紧支撑”，但仍没有完成二维有限质量 return gluing。

### 15.3 横向单模、有限 Fourier 与小有限能余留的严格障碍

先考虑任意固定的 stationary separated directional profile

\[
R=e^{qX}G(Y),
\qquad
P=A_0XY+e^{qX}F(Y),
\qquad q\ne0.
\tag{E105}
\]

比较完整方程中的 \(Xe^{qX}\) 与 \(e^{2qX}\) 系数，并使用椭圆方程，得到

\[
A_0=c,
\qquad
FG'-F'G=0,
\qquad
F''+q^2F+qG=0,
\qquad
2cYG'=(\gamma+bq)G.
\tag{E106}
\]

在每个非零连通分支上 \(G\) 与 \(F\) 成比例。两个 ODE 同时成立只允许

\[
\begin{cases}
F,G\ \text{为常数},&\gamma+bq=0,\\
F,G\propto Y,&\gamma+bq=2c,
\end{cases}
\qquad
G=-qF.
\tag{E107}
\]

所以不存在非零的光滑、衰减、紧支撑或周期 transverse 单包络；两个幸存者
恰好是 whole-plane constant directional mode 与 wall-linear mode。

更朴素地，若

\[
R=r(X)\eta(Y),
\qquad
p=\phi(X)\eta(Y),
\qquad
\eta\in H^2(\mathbb R)\cap L^2(\mathbb R),
\]

则 \(-\Delta p=R_X\) 的变量分离强制

\[
\eta''=\lambda\eta.
\tag{E108}
\]

全直线上没有非零 \(L^2\) 解。因此不能简单给已有 Bessel profile 乘一个
共同的光滑 \(y\)-cutoff。

在固定周期圆柱上，若解始终只有有限个 \(y\)-Fourier 模且各水平系数属于
\(L_x^2\)，最高模 \(N\) 的 \(2N\) 非线性系数强制

\[
p_N'r_N-p_Nr_N'=0.
\]

于是 \(p_N=C(t)r_N\)；再由

\[
-p_N''+(Nk)^2p_N=r_N'
\tag{E109}
\]

作 \(x\)-Fourier 变换，乘子
\(C[\zeta^2+(Nk)^2]-i\zeta\) 迫使 \(r_N=0\)。向下归纳后所有非零
transverse 模都消失。固定 strip 的共同有限 sine 展开有同一最高模障碍。
这个结论不覆盖无限 Fourier cascade、增长频谱或一般二维 return。

最后，若

\[
P=A_0XY+p,\qquad p\in H^1(\mathbb R^2),
\]

则 \(A_0\ne0\) 已经使 \(\nabla P\notin L^2\)；有限能外层必须 order one
地取消 affine strain，而不可能只是小 envelope。若还试图令自诱导流处处
切于密度等值线，即 \(p=\mathcal F(R)\)，并假设足够衰减，则

\[
\int_{\mathbb R^2}|\nabla p|^2
=\int_{\mathbb R^2}pR_X
=\int_{\mathbb R^2}\partial_X\mathcal H(R)=0.
\tag{E110}
\]

故 \(p=0\)、\(R_X=0\)。真正有限能的精确 Bessel 轨道若存在，必然是
非可分、无限模、带 order-one 二维 return/outer reservoir 的结构。这里的
论证没有排除这种结构；它正是剩余的 gluing 与稳定性问题。

## 16. 术语与文献边界

本文的 \(e^{-X}\) 不是径向 \(K_\nu(r)\) 剖面。更准确地说，它是

\[
(\Delta+\partial_X)e^{-X}=0
\]

的 directional drift--Helmholtz homogeneous mode；经指数共轭后，该算子
变成 modified Helmholtz/Bessel 算子。因此“directional Bessel mode”是
简写，最严格的名称应是
**directional modified-Helmholtz/Bessel-operator tangent**。

截至 2026-08-28 的一次非穷尽主文献检索，没有发现与
“\(C_c^\infty\) 平坦支撑边缘在双重随动尺度下精确趋于 \(e^{-X}\)”逐字
相同的二维 IPM 结果，也没有发现对这个 \(e^{-X}\) 极限的扰动吸引域证明。
因此本文只能说“目前检索未见同一机制”，不能据此作新颖性声明。

最接近但不同的结果是：

- Castro--Córdoba--Gancedo--Orive 的
  [无限能特殊类显式爆破](https://arxiv.org/abs/0806.1180)：没有平坦
  紧支撑边缘或 \(e^{-X}\) Bessel tangent；
- Collot--Prange--Tan 的
  [无限能 IPM 自相似爆破稳定性](https://arxiv.org/abs/2507.17381)：
  真正证明了特殊降维类中的 rescaled attractor，但极限是周期余弦，
  不是空间收缩的 Bessel 模；
- Dembski 的
  [锐角域无 boundary mass 奇性](https://arxiv.org/abs/2511.01827)：
  与锐扇双曲压缩更接近，但角点初值只到 Lipschitz，profile 不是本文的
  flat-edge tangent；
- Córdoba--Martínez-Zoroa 的
  [带光滑源的有限能光滑 IPM 奇性](https://arxiv.org/abs/2410.22920)：
  有真正多尺度与 affine 内外桥接，但方程含外源；
- Elgindi 的
  [分层稳态渐近稳定性](https://arxiv.org/abs/1411.6958)：
  说明标准衰减/无通量相空间与本文非衰减 affine far-field contract
  属于不同问题。

所以目前可严格使用的表述是：

\[
\boxed{
\begin{gathered}
\text{本文构造了显式精确轨道，其光滑平坦边缘在随动双尺度下}\\
\text{局部收敛到 directional modified-Helmholtz/Bessel tangent；}\\
\text{这尚不是经过扰动稳定性证明的 Bessel attractor。}
\end{gathered}
}

## 17. 新的二维 \(\rho_x\in C_0\) 审计

若把远场要求明确为

\[
\lim_{|(x,y)|\to\infty}\rho_x(t,x,y)=0
\quad\text{沿所有路径成立},
\tag{E111}
\]

则本报告已有的非退化精确爆破轨道都不合格：水平 front 独立于 \(y\)，旋转
compact strip 沿切向无限延伸，wall-sine 轨道在另一端增长，而 affine 核的
\(\rho_x\) 是空间常数。singular \(K_0/K_1\) 的导数虽在远端衰减，却是静态
punctured-domain 解，不是光滑 Cauchy 爆破。

放弃有限 \(L^p\) 后，最短的新候选是 bounded-contrast profile

\[
(U+cz-be_1)\cdot\nabla R=0,
\qquad -\Delta P=R_X,
\qquad P|_{Y=0}=0,
\tag{E112}
\]

其内端保留本报告的 affine Taylor jet，外端改成方向依赖常数。一个显式相容
首项为

\[
F_0=A\left(\cos\theta+\frac53\cos3\theta\right),
\qquad
G_0=-\frac{4A}{3}\sin^4\theta,
\qquad R_X=O(r^{-1}).
\tag{E113}
\]

它已经通过 leading Poisson、Dirichlet Fredholm 条件和第一次输运修正，但尚未
与内核连成完整二维 fixed point。因而 (E112) 是下一次数值 Newton/谱计数的
正确目标；directional Bessel 只保留为有限匹配环带内的 Green corrector。
详见 [RHOX_DECAY_BLOWUP_STATUS_ZH.md](RHOX_DECAY_BLOWUP_STATUS_ZH.md)。

## 18. 强制对数外端与开放 harmonic sheet 轨道

(E112) 的壁面条件进一步给出精确乘积恒等式

\[
(cX-b-P_Y)R_X=0\qquad (Y=0).
\tag{E114}
\]

若外端只有 $P=O(r)$，则 sufficiently far wall tails 上必须有
$R_X\equiv0$；合格候选只能使用端点无限阶平坦的 compact-active wall
segment。以 (E113) 为首项作 Laurent 递推时，第三阶 Poisson 右端对
$\sin2\theta$ 的 Fredholm 投影为

\[
-\frac{25}{48}\frac{A^4}{c^3}\ne0.
\tag{E115}
\]

所以纯幂外端不成立，最小修正被强制为

\[
P_{3,\log}
=\frac{25}{192}\frac{A^4}{c^3}
r^{-2}\log r\,\sin2\theta.
\tag{E116}
\]

指数小尾不能消除这个代数阶共振。正确数值 unknown 必须包含弯曲
hodograph 特征、无限 transverse hierarchy 和 polyhomogeneous log tail。
第一轮 $129\times49$、448-mode 搜索虽把 Poisson 残差压到
$1.85\times10^{-11}$，bulk/wall transport RMS 仍为 $0.2104/0.1044$，故明确
不是 fixed profile。

另一方面，若允许预先指定

\[
\Psi_H=\frac{xy}{\tau},\qquad
\tau=T-t,qquad X=\frac{x}{\tau},\qquad Y=\tau y,
\]

并令 $\sigma=\tau^2/2$、$\Psi_s=\tau\Phi$，则二维问题精确约化为

\[
F_\sigma=\{\Phi,F\},
\qquad
-(\Phi_{XX}+4\sigma^2\Phi_{YY})=F_X,
\qquad
\rho_x=\tau^{-1}F_X.
\tag{E117}
\]

椭圆逆的 partial-$Y$ Fourier 核是

\[
\widehat\Phi(X,\eta)
=-\frac12\int\operatorname{sgn}(X-Z)
e^{-2\sigma|\eta||X-Z|}\widehat F(Z,\eta)\,dZ.
\tag{E118}
\]

它在 $L_X^1$--analytic scale 中关于 $\sigma\ge0$ 一致有界，标准
Ovsyannikov 迭代可从 $F_*=\varepsilon e^{-X^2-Y^2}$ 构造解析衰减轨道，且

\[
\tau\|\rho_x(t)\|_\infty
\longrightarrow |\varepsilon|\sqrt{2/e}.
\tag{E119}
\]

但这是 collapsing-$x$/expanding-$y$ 的竖直 singular sheet，不是锐角点；并且
$\Psi_H$ 二次增长，其未来奇异的远场合同不由初始密度决定。它只是逐点无体源的
exact trajectory，不能当作标准自治 Cauchy 爆破。

这个差别不能用小 corrector 消除。去掉 $\Psi_H$ 后，同一坐标给出

\[
2\sigma F_\sigma=XF_X-YF_Y+2\sigma\{\Phi,F\}.
\tag{E120}
\]

若 self bracket 有界，终端只能是 $F_*=H(XY)$；再由
$F_{*X}=YH'(XY)\in C_0$ 得 $H'\equiv0$。所以任何标准 self-screened
替代都必须在 leading order 激活 $O(\sigma^{-1})$ 的 return/fast layer，而
不是对 harmonic sheet 轨道作微扰。
