# 边界跳跃、有限质量幂律尾与右锐扇 IPM 有限时爆破

## 0. 结论和合同

本文给出一个不需要 Bessel 尾、不需要外加 harmonic strain 的
边界驱动爆破机制。物理域是以 \((1,0)\) 为顶点的右开锐角楔形

\[
 \Omega_L=\{(x_1,x_2):X=x_1-1>0,\ 0<x_2<X\tan L\},
 \qquad 0<L<\frac{\pi}{2}.
\tag{1}
\]

在两条边上施加不可穿透条件 \(u\cdot n=0\)，等价地取
\(\psi=0\)。始终选取由 Dirichlet Green 核给出的有限能/Friedrichs
分支；不允许另加在无穷远增长的 harmonic strain。方程是

\[
 \rho_t+u\cdot\nabla\rho=0,
 \qquad -\Delta\psi=\partial_X\rho,
 \qquad u=\nabla^\perp\psi=(-\psi_{x_2},\psi_X).
\tag{2}
\]

我们在底边 \(x_2=0\) 保留非零上侧迹。若按当前数值合同作奇延拓，

\[
 \widetilde\rho(X,-x_2)=-\widetilde\rho(X,x_2),
\tag{3}
\]

则 \(x_2=0\) 上有指定跳跃，但本文的经典方程始终只在
\(\Omega_L\) 的上侧求解。

得到的精确结论是：存在明确的光滑、非负、有限质量幂律尾初值，
其唯一最大经典 IPM 解必在有限时间内发生梯度爆破。这里“精确解”
是由明确 Cauchy 数据唯一确定的真实 PDE 解，不是把一个近似 profile
称作解；有限质量解一般没有闭式时空公式。密度
\(L^1\) 质量和 \(L^\infty\) 振幅保持有界，每个经典时刻都有
\(\nabla\rho(x)\to0\ (|x|\to\infty)\)。正的水平梯度核位于右锐扇顶点，
幂律远尾自动提供同一扇内的负 return。

本文的“经典/强解”是指密度在单侧闭楔形中具有 \(C^{1,\alpha}\) 角点
正则性、速度至少为 \(C^{1,\delta}\)，并在顶点外光滑，其中可取
\(0<\delta<\min\{\alpha,\pi/L-2\}\)。一般 Dirichlet 解可含
\(r^{\pi/L}\sin(\pi\theta/L)\) 角点模，所以除非角度特殊，不能把速度在
顶点误称为 \(C^\infty\)。这不影响二阶迹、流映射或 Riccati 恒等式。

## 1. 完全显式的 affine 内核

记

\[
 m=\tan L,
 \qquad q=X-\frac{x_2}{m}
 =r\frac{\sin(L-\theta)}{\sin L}.
\tag{4}
\]

在 \(0<\theta<L\) 中 \(q>0\)，而在上方斜边 \(x_2=mX\) 上
\(q=0\)。令

\[
 \rho_A(t,X,x_2)=a(t)q,
\tag{5}
\]

\[
 \psi_A(t,X,x_2)
 =a(t)\left(\frac m2Xx_2-\frac12x_2^2\right).
\tag{6}
\]

则 \(\psi_A=0\) 在两条边上，并且

\[
 -\Delta\psi_A=a(t)=\partial_X\rho_A.
\tag{7}
\]

速度为

\[
 u_A=a(t)\left(x_2-\frac m2X,\frac m2x_2\right).
\tag{8}
\]

直接代入输运方程得到唯一标量 ODE

\[
 \boxed{a'(t)=\frac m2a(t)^2.}
\tag{9}
\]

因此，对 \(a_0>0\)，

\[
 a(t)=\frac{a_0}{1-\frac12ma_0t},
 \qquad T_A=\frac{2}{ma_0}.
\tag{10}
\]

且

\[
 \partial_X\rho_A=a(t)>0,
 \qquad
 \partial_{x_2}\rho_A=-\frac{a(t)}m.
\tag{11}
\]

这是一条完全显式的含时 IPM 爆破轨道。它的唯一缺点是
\(\rho_A\) 在无穷远线性增长。下一节证明，有限时爆破并不依赖这个
无限质量外尾。

## 2. 角点 Hessian 恒等式

设 \(\rho\) 是经典解，且上方斜边上 \(\rho=0\)。该条件由边界输运
保持。记

\[
 a(t)=\partial_X\rho(t,0,0).
\tag{12}
\]

斜边的切向相容性给出

\[
 \partial_{x_2}\rho(t,0,0)=-\frac{a(t)}m.
\tag{13}
\]

在顶点附近，

\[
 \partial_X\rho=a(t)+O(r^\alpha).
\tag{14}
\]

先只使用边界的二阶 jet。由 \(\psi(X,0)=0\)，

\[
 \psi_X(0)=\psi_{XX}(0)=0.
\tag{14a}
\]

再对 \(\psi(X,mX)=0\) 求两次 \(X\) 导数，并在顶点取值，得到

\[
 \psi_{XX}(0)+2m\psi_{Xx_2}(0)+m^2\psi_{x_2x_2}(0)=0.
\tag{14b}
\]

Poisson 方程在顶点给出

\[
 -\psi_{XX}(0)-\psi_{x_2x_2}(0)=a(t).
\tag{14c}
\]

联立 (14a)--(14c) 已经形式上给出
\(\psi_{Xx_2}(0)=ma(t)/2\)。下面的角点展开说明这些二阶迹确实存在，
也排除了从远场进入的低阶 harmonic 模。

定义二次多项式

\[
 \psi_2=a(t)\left(\frac m2Xx_2-\frac12x_2^2\right).
\tag{15}
\]

它严格满足两条 Dirichlet 边界和
\(-\Delta\psi_2=a(t)\)。对 \(w=\psi-\psi_2\)，

\[
 -\Delta w=O(r^\alpha),
 \qquad w|_{\partial\Omega_L}=0.
\tag{16}
\]

锐角 Dirichlet Laplacian 的第一个非零齐次调和模次数是

\[
 \lambda_1=\frac\pi L>2.
\tag{17}
\]

更准确地，对任意

\[
 0<\sigma<\min\left\{\alpha,\frac\pi L-2\right\}
\tag{17a}
\]

角点 Schauder 展开给出

\[
 w=O(r^{2+\sigma}),\qquad D^2w=O(r^\sigma).
\tag{17b}
\]

因此 \(D^2w(0,0)=0\)，从而得到不受远场影响的精确局部恒等式

\[
 \boxed{\psi_{Xx_2}(t,0,0)=\frac m2a(t).}
\tag{18}
\]

这一步是整个机制的核心。在半平面中 \(Xx_2\) 是可自由加入的
调和二次模；在严格锐角中，第二条不可穿透边把它唯一锁定为
\((m/2)a\)。

## 3. 壁面 Riccati 爆破

底边迹记为

\[
 f(t,X)=\rho(t,X,0^+).
\tag{19}
\]

由 \(\psi(X,0)=0\) 得 \(u_2(X,0)=0\)，因此

\[
 f_t+u_1(t,X,0)f_X=0.
\tag{20}
\]

顶点上 \(u_1=0\)。对 (20) 求 \(X\) 导数，并使用
\(u_{1X}=-\psi_{Xx_2}\)，便得

\[
 a'(t)=-u_{1X}(t,0,0)a(t)
 =\psi_{Xx_2}(t,0,0)a(t).
\tag{21}
\]

结合 (18)，

\[
 \boxed{a'(t)=\frac m2a(t)^2,
 \qquad
 a(t)=\frac{a_0}{1-\frac12ma_0t}.}
\tag{22}
\]

这个 ODE 对任意合格的外层都成立：外层可以是紧支撑、指数尾或
幂律尾，只要它不改变顶点的一次 jet 和边界相容性。

因为 \(\partial_X\rho\) 在单侧角点连续且 \(a(t)>0\)，对每个
\(t<T_*\) 都存在 \(\varepsilon(t)>0\) 使

\[
 \partial_X\rho(t,X,x_2)>0
 \quad\text{当 }(X,x_2)\in\Omega_L, r<\varepsilon(t).
\tag{22a}
\]

所以统一正号不是只在一条边上成立，而是在右锐扇顶点附近的二维核心中成立。

此外，底边特征流在每个经典时刻都是 \([0,\infty)\) 上保持端点的
微分同胚。于是由 (24) 型正迹出发，

\[
 f(t,X)>0\quad(X>0),
 \qquad
 [\widetilde\rho(t)]_{x_2=0}=2f(t,X)\ne0.
\tag{22b}
\]

奇延拓不只是画图约定：将 \(\psi\) 同样作奇延拓，则 \(u_1\) 为偶函数、
\(u_2\) 为奇函数，所以界面法向速度 \(u_2|_{x_2=0}=0\)。输运方程可能
产生的界面测度系数正是 \(u_2[\widetilde\rho]\)，因而为零；Poisson 方程
右端是切向导数 \(\partial_X\widetilde\rho\)，也不会微分出界面 Dirac。
故它给出双锐扇 \(\{|x_2|<X\tan L\}\) 内逐侧经典、跨中线有跳跃的
精确弱解。

## 4. 明确的有限质量幂律尾初值

取

\[
 \boxed{
 \rho_0(X,x_2)
 =a_0\frac{X-x_2/m}{(1+X^2+x_2^2)^3},
 \qquad (X,x_2)\in\Omega_L .}
\tag{23}
\]

它在闭楔形的上侧是 \(C^\infty\) 的，在斜边为零，在底边的上侧迹为

\[
 f_0(X)=a_0\frac{X}{(1+X^2)^3}>0
 \quad (X>0).
\tag{24}
\]

它在开楔形中严格为正，所以相对于物理域的支撑集恰是闭右扇
\(\overline{\Omega_L}\)；“集中于扇形”在这里是几何域局部化加扇内
代数衰减，并非扇内紧支撑，也不是把同一函数置零后冒充半平面解。

所以奇延拓的跳跃量是

\[
 [\widetilde\rho_0]_{x_2=0}
 =2a_0\frac{X}{(1+X^2)^3}.
\tag{25}
\]

在极坐标中，

\[
 \rho_0
 =a_0\frac{r\sin(L-\theta)}
 {\sin L\,(1+r^2)^3}
 \sim a_0\frac{\sin(L-\theta)}{\sin L}\,r^{-5}.
\tag{26}
\]

因此远场严格由右锐扇内的 \(r^{-5}\) 幂律主导，并且

\[
 |\nabla\rho_0|=O(r^{-6})\longrightarrow0.
\tag{27}
\]

质量还可精确计算：

\[
 \int_{\Omega_L}\rho_0\,dx
 =a_0\left(\int_0^L\frac{\sin(L-\theta)}{\sin L}\,d\theta\right)
 \left(\int_0^\infty\frac{r^2}{(1+r^2)^3}\,dr\right)
 =\boxed{\frac{\pi a_0}{16}\tan\frac L2}.
\tag{28}
\]

同时

\[
 \|\rho_0\|_\infty
 =\frac{125a_0}{216\sqrt5}<\infty.
\tag{29}
\]

指数 \(5\) 不是特殊常数。更一般地，对任意 \(p>3/2\)，

\[
 \rho_{0,p}=a_0\frac{q}{(1+r^2)^p}
 \sim a_0\frac{\sin(L-\theta)}{\sin L}\,r^{-(2p-1)}
\tag{29a}
\]

仍有相同的角点斜率和 Riccati 爆破，并且

\[
 \int_{\Omega_L}\rho_{0,p}
 =a_0\tan\frac L2\,
 \frac{\sqrt\pi}{4}\frac{\Gamma(p-\tfrac32)}{\Gamma(p)}<\infty,
 \qquad
 |\nabla\rho_{0,p}|=O(r^{-2p}).
\tag{29b}
\]

因此可实现任意可积密度幂律 \(r^{-s}\)，\(s=2p-1>2\)。本文继续用
\(p=3\) 是因为质量和峰值都化为初等常数。

这个幂律尾不是装饰。壁面上

\[
 f_0'(X)=a_0\frac{1-5X^2}{(1+X^2)^4}.
\tag{30}
\]

所以顶点附近 \(\rho_X>0\)，而 \(X>1/\sqrt5\) 后的幂律尾是负
return。这正好实现

\[
 \int_0^\infty f_0'(X)\,dX=0,
\tag{31}
\]

即有限质量密度必需的 source-neutrality。不可能在整个无界扇形中
保持 \(\rho_X\) 单符号；能严格保留的是爆破核内的统一正号。

## 5. 有限时爆破定理

先记录本节使用的延拓引理。取
\(0<\alpha<\min\{1,\pi/L-2\}\)，用
\(\mathring C^\alpha\) 表示适配角点缩放的 Hölder 范数；(23) 的远场
多项式权作为额外传播性质处理，不放进延拓范数。其定义为

\[
 \|h\|_{\mathring C^\alpha}
 :=\|h\|_\infty+
 \sup_{x\ne y}\frac{\big||x|^\alpha h(x)-|y|^\alpha h(y)\big|}
 {|x-y|^\alpha}.
\tag{31a0}
\]

严格锐角 Dirichlet 楔形中的强解若满足

\[
 \int_0^T\|\nabla\rho(t)\|_{L^\infty}\,dt<\infty,
\tag{31a}
\]

则可延拓越过 \(T\)。证明是标准 BKM 论证，但这里需要说明它为何不被角点
破坏。令 \(\kappa=\pi/L>2\)。映射 \(z\mapsto z^\kappa\) 给出楔形
Dirichlet Green 函数

\[
 G_L(z,\zeta)=\frac1{2\pi}
 \log\left|\frac{z^\kappa-\overline{\zeta}^{\,\kappa}}
 {z^\kappa-\zeta^\kappa}\right|.
\tag{31b}
\]

Elgindi--Jeong 的 sector Poisson Lemma 3.2/3.5 对任意有界右端给出唯一
Dirichlet 解及端点对数 Calderón--Zygmund 估计；其 Remark 3.12 说明在
严格锐角时不需要额外反射对称，且密度本身不需要边界条件。Dembski 的
Proposition 2.1 与 (2.5) 将同一估计用于 IPM 的局部适定性与 Hölder 延拓。
记

\[
 g(t)=\|\nabla\rho(t)\|_\infty,
 \qquad A(t)=\|\nabla\rho(t)\|_{\mathring C^\alpha}.
\tag{31c0}
\]

把上述已证明估计用于 \(f=\partial_X\rho\)，得到

\[
 \|\nabla u\|_\infty
 \le C_L g(t)\left[1+
 \log\!\left(1+c_L\frac{A(t)}{g(t)}\right)\right],
\tag{31c}
\]

其中 \(g=0\) 时按连续极限理解。对梯度输运方程作尺度不变 Hölder
差商估计（即 Dembski (2.5)），同时得到

\[
 A'(t)\le C_L\bigl(\|\nabla u\|_\infty+g(t)\bigr)A(t).
\tag{31d}
\]

利用

\[
 g\log(1+cA/g)\le C_L\{1+g\log(e+A)\}
\tag{31e}
\]

并对 \(Z=\log(e+A(t))\) 使用 Gronwall，(31a) 便排除
有限时角点 Hölder 范数失规；再代回 (31c) 还得到
\(\int_0^T\|\nabla u\|_\infty dt<\infty\)，故局部理论确实可以延拓。
这就是所需的 acute-wedge BKM 引理。非零底边迹只进入 Dirichlet Poisson
右端，不要求把有跳跃的奇延拓放进强解空间；所有估计均在上侧物理楔形内完成。

令 \(\rho(t)\) 是由 (23) 出发的唯一最大经典解，其最大存在时间为
\(T_*\)。锐角 IPM 的局部适定性给出 \(T_*>0\)。在所有
\(t<T_*\) 上，(22) 严格成立。因此

\[
 T_*\le T_A=\frac{2}{a_0\tan L}<\infty.
\tag{32}
\]

若解存在到 \(T_A\)，则顶点水平梯度给出精确 Type-I 率

\[
 \partial_X\rho(t,1,0^+)
 =\frac{2}{\tan L}\frac1{T_A-t}.
\tag{33}
\]

若在 \(T_A\) 之前出现更早的失规，IPM 的 BKM 型经典延拓准则则给出

\[
 \int_0^{T_*}\|\nabla\rho(t)\|_{L^\infty}\,dt=\infty.
\tag{34}
\]

因而无论哪个分支发生，都有

\[
 \boxed{
 T_*<\infty,
 \qquad
 \limsup_{t\uparrow T_*}\|\nabla\rho(t)\|_\infty=\infty.}
\tag{35}
\]

对所有 \(t<T_*\)，不可压流映射保持分布函数，所以

\[
 \|\rho(t)\|_{L^1}=\|\rho_0\|_{L^1},
 \qquad
 \|\rho(t)\|_{L^\infty}=\|\rho_0\|_{L^\infty}.
\tag{36}
\]

同一流映射是 \(\Omega_L\) 的微分同胚；由于 \(\rho_0>0\) 于开楔，
\(\rho(t)>0\) 也在整个开楔保持。因此正密度的几何支撑始终是右锐扇，
不会从壁面漏出。

对应双楔奇延拓的绝对质量为上式的两倍，仍然有限；其有符号总质量则因
奇对称为零。

远场衰减还有一个不依赖 Biot--Savart 渐近展开的直接证明。令

\[
 K(t)=\int_0^t\|\nabla u(s)\|_{L^\infty}\,ds<\infty,
\tag{36a}
\]

并记流映为 \(\Phi_t\)。两条不平行的 Dirichlet 边在顶点给出
\(u(t,0)=0\)，所以

\[
 e^{-K(t)}|z|\le |\Phi_t(z)|\le e^{K(t)}|z|,
 \qquad
 \|D\Phi_t^{\pm1}\|_\infty\le e^{K(t)}.
\tag{36b}
\]

输运的拉格朗日公式为

\[
 \rho(t,\Phi_t(z))=\rho_0(z),\qquad
 \nabla\rho(t,\Phi_t(z))=D\Phi_t(z)^{-T}\nabla\rho_0(z).
\tag{36c}
\]

因此 (26)--(27) 的幂阶在每个经典时刻严格保持：

\[
 \rho(t,x)=O_t(|x|^{-5}),\qquad
 \nabla\rho(t,r,\theta)=O_t(r^{-6})\to0.
\tag{37}
\]

这还是一个真正的幂律主项，而不只是上界。流映保持两条壁，并且双
Lipschitz 映射把到边界的距离保持到常数倍。因此对任意
\(0<\varepsilon<L/2\)，存在 \(c_{t,\varepsilon},C_{t,\varepsilon}>0\)，使充分大
的 \(r\) 上

\[
 c_{t,\varepsilon}r^{-5}
 \le \rho(t,r,\theta)\le
 C_{t,\varepsilon}r^{-5},
 \qquad \varepsilon\le\theta\le L-\varepsilon.
\tag{37a}
\]

这里不宣称固定 Euler 角度上的首项系数不变；流可以扭曲角变量。所需的
可积幂律阶和全方向远场梯度消失则是严格的。

## 6. 与 Bessel、半平面和“扇形支撑”的关系

1. **Bessel 不再是最终外尾。**  
   本解的远场是可积的 \(r^{-5}\) 幂律，密度和速度都没有 Bessel
   方向模带来的无限横向质量。

2. **扇形是真实物理域。**  
   这里的第二条斜边是不可穿透壁，不能将解简单置零延拓成
   整个半平面；否则流函数的法向导数跳跃会产生人工 sheet source。

3. **\(x_2=0\) 的跳跃是本机制所需的 boundary mass。**  
   如果强迫底边密度也为零，则 (24) 的壁面斥力方向不存在，问题回到
   Dembski 的 no-boundary-mass Lipschitz 角向机制。

4. **不能要求整个无界扇形中 \(\rho_X\) 都为正。**  
   有限质量和 \(\rho\to0\) 强制水平导数总积分为零。本解保证顶点
   爆破核中 \(\rho_X>0\)，并把必需的负 return 放到同一右扇的幂律远尾。

## 7. 证据等级

- (5)--(11) 是逐项可验证的显式精确解。
- (18)--(22) 是锐角 Dirichlet 正则性给出的精确局部闭合，不是渐近拟合。
- (23) 是一个具体的光滑有限质量初值；(28)--(31) 是精确计算。
- (32)--(36) 结合标准局部适定性和 IPM 延拓准则，给出真正的有限时
  梯度爆破，而不是一个数值 candidate。
- (5)--(11) 是闭式时空公式；(23) 所产生的有限质量解是唯一 Cauchy
  演化意义下的精确解，而不是另一个闭式 ansatz。这两种“精确”不能混淆。
- 初始密度 (23) 在单侧闭楔形中是 \(C^\infty\)；受非整数角点模影响，
  相应速度一般只保证上述 \(C^{1,\delta}\) 强正则性。奇延拓是双楔形中的
  分片强、跨中线跳跃弱解，不是全平面光滑解。
- 该结论不解决“整个无边界平面上 \(C_c^\infty\) 初值爆破”；它解决的正是
  本文固定的“底边跳跃 + 右锐扇壁 + 有限质量幂律尾”合同。

## 8. 相关严格文献边界

- Kevin H. Dembski,
  [*Singularity Formation in the Incompressible Porous Medium Equation without Boundary Mass*](https://arxiv.org/abs/2511.01827),
  2025/2026。该文在锐角域对无 boundary mass 的 Lipschitz 角向核证明爆破，
  并明确说明有 boundary mass 时可用 ODE 机制。本文把后一机制写成
  (18)--(22) 的显式 Riccati 闭合，并给出 (23) 的有限质量幂律尾。
- Tarek M. Elgindi, In-Jee Jeong,
  [*Finite-time Singularity Formation for Strong Solutions to the Boussinesq System*](https://arxiv.org/abs/1708.02724),
  Lemma 3.2/3.5 给出严格锐角 sector 的唯一 Dirichlet Poisson 解、
  \(\mathring C^\alpha\) 椭圆估计和 (31c) 所用的端点对数估计；
  Remark 3.12 允许在严格锐角中去除额外对称假设。
- Diego Córdoba, Francisco Gancedo, Rafael Orive,
  [*Analytical behavior of two-dimensional incompressible flow in porous media*](https://doi.org/10.1063/1.2404593),
  提供标准 IPM 局部适定性与延拓/爆破准则背景。

因此这里不声称“boundary mass 导致 ODE 爆破”是首次发现；文献已经明确
指出这一机制。本文的新工作是把当前几何中的系数 \(\tan L/2\)、持续跳跃、
有限质量 \(r^{-5}\) 尾和全方向远场梯度衰减放进同一个可逐项核验的合同中。
