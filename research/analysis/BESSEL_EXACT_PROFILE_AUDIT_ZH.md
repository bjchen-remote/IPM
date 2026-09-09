# 带 tip 源的精确 Bessel profile：独立核验与继续推导

日期：2026-08-27

本文核验同步得到的 Bessel 解，并进一步提取它与真实动态爆破之间缺少的
方程。坐标以候选点 \((1,0)\) 为中心：

\[
X=x_1-1,\qquad Y=x_2,\qquad r=(X^2+Y^2)^{1/2}.
\]

## 1. 精确性核验

物理 IPM 采用

\[
\rho_t+U\cdot\nabla\rho=0,
\qquad -\Delta\Psi=\partial_X\rho,
\qquad U=\nabla^\perp\Psi.
\]

令 \(\rho=\lambda\Psi\)。由于

\[
U\cdot\nabla\rho
=\lambda\nabla^\perp\Psi\cdot\nabla\Psi=0,
\]

稳态 Poisson 方程等价于

\[
\Delta\Psi+\lambda\partial_X\Psi=0.
\]

写成

\[
\Psi=e^{-\lambda X/2}\Phi,
\qquad k=\frac{|\lambda|}{2},
\]

则严格化为

\[
(\Delta-k^2)\Phi=0.
\]

因此

\[
\Psi_\nu=e^{-\lambda X/2}K_\nu(kr)
\bigl[a_\nu\cos(\nu\theta)+b_\nu\sin(\nu\theta)\bigr],
\qquad \rho_\nu=\lambda\Psi_\nu
\]

在不含 Bessel 极点的域内确实逐点满足完整物理稳态方程。这不是拟合或
渐近式。

最重要的右向 \(K_0\) 成员是

\[
\boxed{
\Psi_0=-Ae^{\mu X/2}K_0\!\left(\frac{\mu r}{2}\right),
\qquad
\rho_0=\mu Ae^{\mu X/2}K_0\!\left(\frac{\mu r}{2}\right),
\qquad A,\mu>0.
}
\]

它满足 \(\rho_0=-\mu\Psi_0\)。不可穿壁的最低阶成员为

\[
\boxed{
\Psi_1=-Ae^{\mu X/2}K_1\!\left(\frac{\mu r}{2}\right)\sin\theta,
\qquad
\rho_1=\mu Ae^{\mu X/2}K_1\!\left(\frac{\mu r}{2}\right)\sin\theta.
}
\]

二者还有精确关系

\[
\rho_1=-\frac{2}{\mu}\partial_Y\rho_0.
\]

## 2. “衰减”成立到什么程度

当 \(X\to+\infty\)、\(Y=O(\sqrt X)\) 时，

\[
\rho_0\sim A\sqrt{\frac{\pi\mu}{X}}
\exp\!\left(-\frac{\mu Y^2}{4X}\right),
\]

\[
\rho_1\sim A\sqrt{\pi\mu}\,
\frac{Y}{X^{3/2}}
\exp\!\left(-\frac{\mu Y^2}{4X}\right).
\]

所以两者都在 \(r\to\infty\) 时一致趋于零，并集中到角宽
\(\theta=O(r^{-1/2})\) 的右向缩角区。这严格满足“远场函数值衰减”，
但不等价于有限质量或有限能量：

\[
\rho_0\in L^p(\text{远场})\iff p>3,
\qquad
\nabla\rho_0\in L^p(\text{远场})\iff p>\frac32,
\]

\[
\rho_1\in L^p(\text{远场})\iff p>\frac32.
\]

特别地，\(K_0\) 密度不是远场 \(L^1\) 或 \(L^2\)；\(K_1\) 密度在远场
属于 \(L^2\)，但其中心行为 \(\rho_1\sim2A\sin\theta/r\) 使全局
\(L^2\) 在顶点对数发散。

中心处

\[
\rho_0\sim\mu A\log\frac1r,
\qquad |\nabla\rho_0|\sim\frac{\mu A}{r},
\]

而 \(K_1\) 更强：

\[
\rho_1\sim2A\frac{Y}{r^2},
\qquad |\nabla\rho_1|=O(r^{-2}).
\]

因此这些解是“远场衰减但中心带极点”的精确外层/中间层。\(K_0\) 在全平面
分布意义还携带顶点点源；\(K_1\) 满足齐次壁面条件，但携带边界偶极奇性。

## 3. 水平导数并非全域统一符号

对 \(K_0\)，令 \(z=\mu r/2\)，则

\[
\partial_X\rho_0
=\frac{\mu\rho_0}{2}
\left[1-\frac{K_1(z)}{K_0(z)}\cos\theta\right].
\]

它在正水平射线附近为负，在离开该射线后为正；零线满足

\[
\cos\theta=\frac{K_0(z)}{K_1(z)},
\qquad
Y^2\sim\frac{2X}{\mu}.
\]

\(K_1\) 成员也换号，其远场零线为

\[
Y^2\sim\frac{6X}{\mu}.
\]

故 Bessel 解满足右向缩角和远场衰减，但不满足“整个第一象限
\(\partial_{x_1}\rho\) 统一符号”。它与此前线性扇形解解决的是不同约束。

## 4. 新的精确闭包缺陷方程

这是继续推导的关键。对任意含时状态定义

\[
\boxed{D:=\rho+\mu\Psi,}
\]

其中 \(\mu>0\) 暂取常数。因为 \(\rho=-\mu\Psi+D\)，完整 IPM 可精确改写为

\[
\boxed{
(\Delta-\mu\partial_X)\Psi=-\partial_XD,
\qquad
\rho_t+U\cdot\nabla D=0.
}
\tag{D1}
\]

第二式利用了 \(U\cdot\nabla\Psi=0\)。这说明：

1. 精确 Bessel 流形就是 \(D=0\)。
2. 若 \(D=0\)，输运方程立即给出 \(\rho_t=0\)。在采用衰减
   Biot--Savart 规范、没有额外调和背景流时，不能通过令
   \(A=A(t)\)、\(\mu=\mu(t)\)、中心 \(a=a(t)\) 或核心尺度
   \(\ell=\ell(t)\)，把同一 Bessel 闭式族变成非平凡含时爆破解。
3. 真正的尺度选择和核心收缩必定由 \(D\neq0\) 驱动。换言之，gluing
   误差不是可以最后再消掉的小修正；它正是调制动力的来源。

若允许 \(\mu=\mu(t)\)，仍可定义 \(D=\rho+\mu(t)\Psi\)，且空间方程
保持 (D1) 的形式；时间方程中的
\(-\mu'\Psi-\mu\Psi_t+D_t\) 必须由 \(U\cdot\nabla D\) 平衡。
因此 \(\mu'\) 只能由非零闭包缺陷的可解性条件产生。

## 5. 局域多极矩障碍：有界内核不能小扰动地切掉奇点

上一版曾用

\[
s=(r^2+\ell^2)^{1/2},\qquad
\Phi_\ell=K_0(ks),\qquad k=\frac\mu2
\]

测试将 \(K_0\) 奇点正则化。虽然

\[
(\Delta-k^2)\Phi_\ell
=\ell^2\left[-\frac{k^2K_0(ks)}{s^2}
+\frac{2\partial_sK_0(ks)}{s^3}\right]
\tag{D2}
\]

对每个固定 \(r>0\) 点态趋于零，但它并不在 \(L^1\) 中变小。实际上该残差为负，且

\[
\int_{\mathbb R^2}(\Delta-k^2)\Phi_\ell\,dx
=-2\pi k\ell K_1(k\ell)\longrightarrow-2\pi.
\tag{D3}
\]

因此只看离开顶点的点态残差会误判 gluing 误差。完整 Poisson 方程会保留
单极/偶极矩，并给出更强的必要条件。

对内点圆 \(B_r\)，积分 \(-\Delta\Psi=\partial_X\rho\) 得

\[
-\int_{\partial B_r}\partial_n\Psi\,ds
=\int_{\partial B_r}\rho n_X\,ds.
\tag{M0}
\]

右端中常数背景自动消去。若
\(\delta_r=\operatorname{osc}_{\partial B_r}\rho\)，则 \(K_0\) 对数通量强制

\[
2\pi|A|+o(|A|)\lesssim r\delta_r.
\tag{M1}
\]

对上半圆中满足 \(\Psi=0\) 壁面条件的 \(K_1\) 支，以调和函数
\(h=Y\) 作 Green 恒等式：

\[
\int_{\partial B_r^+}
\bigl(Y\partial_n\Psi-\Psi n_Y\bigr)\,ds
=-\int_{\partial B_r^+}Y\rho n_X\,ds.
\tag{M2}
\]

左端抽取法向偶极系数，因而

\[
\frac{|A|}{\mu}\lesssim r^2\delta_r.
\tag{M3}
\]

令 \(q=\mu r\ll1\)。\(K_0\) 的 Bessel 密度振幅为

\[
\delta_B\sim\frac{|A|}{r}\,q\log\frac1q,
\]

而 (M1) 要求修正量至少为 \(|A|/r\)；修正与 Bessel 分量之比至少为

\[
\frac{1}{q\log(1/q)}\longrightarrow\infty.
\]

\(K_1\) 中 \(\delta_B\sim |A|/r\)，但 (M3) 要求一个至少为
\(\delta_B/q\) 的非对称修正。因此

\[
\boxed{
\ell\ll r\ll\mu^{-1}
\quad\text{中的 small-argument }K_0/K_1
\quad\text{不能是主导中间层。}
\]

这严格否定了固定 \(\mu\) 下“有界 core + 对数/偶极 Bessel 主层”的
小扰动式匹配。特别地，仅由 \(L^\infty\) 匹配得到的
\(A\sim1/[\mu\log(1/\ell)]\) 不是可行的全局缩放律。

## 6. 分布源、非局域闭包层与外层补偿

在全平面分布意义下，两个 Bessel 分支满足

\[
(\Delta-\mu\partial_X)\Psi_0=2\pi A\delta_0,
\]

\[
\Psi_1=-\frac2\mu\partial_Y\Psi_0,
\qquad
(\Delta-\mu\partial_X)\Psi_1
=-\frac{4\pi A}{\mu}\partial_Y\delta_0.
\tag{S1}
\]

结合 (D1)，并选取 \(D(-\infty,Y)=0\) 的规范，局域无源正则化必须产生

\[
D_0\simeq-2\pi A H(X)\delta(Y),
\]

或

\[
D_1\simeq\frac{4\pi A}{\mu}H(X)\partial_Y\delta(Y).
\tag{S2}
\]

所以 \(D\) 不是只支撑在 \(r=O(\ell)\) 的小残差，而是沿右射线延伸的
monopole/dipole 闭包层。若要求远场衰减，它还需在远端连接一个反号
return layer。Gaussian 尾迹继承同一个 tip 源或偶极矩，不会自动把它消掉。

简单 \(K_0(\sqrt{r^2+\ell^2})\) 正则化尤其明确。当 \(\mu\ell\ll1\) 时，
由 (D1) 恢复的 \(D_\ell\) 在下游满足

\[
D_\ell(+\infty,Y)
\sim-
\frac{\pi q_0\ell^2}
{\mu(Y^2+\ell^2)^{3/2}},
\qquad q_0=\mu A,
\]

从而

\[
\|D_\ell\|_\infty\asymp\frac{q_0}{\mu\ell}.
\]

若仍强行取 \(q_0\log(1/(\mu\ell))=O(1)\)，该闭包缺陷反而发散。

## 7. 当前最强结论

同步得到的 Bessel 族仍是重要的**穿孔域精确稳态公式/带分布源的奇异
稳态核/漂移 Helmholtz 外层核**，但它的正确角色已经改变：

- \(D=0\) 的 Bessel 族在物理时间中只能静态；
- 固定 \(\mu\) 的奇异 Bessel 分量不能主导一个平滑有界内核的小尺度匹配；
- Bessel 主导最早只能从其自然尺度 \(r_B\sim\mu^{-1}\) 开始；
- 若爆破梯度由自然尺度 \(r_B\sim\mu^{-1}\) 的 Bessel 响应贡献，则一般必要条件是

\[
\mu(t)\to\infty,
\qquad
0\leq\delta_B(t)\lesssim M,
\qquad
\mu(t)\delta_B(t)\to\infty.
\tag{W}
\]

其中 \(M=\|\rho^{\rm in}\|_\infty\) 是有界初始物理解的密度上界，而不是
奇异 \(K_0\) profile 的 \(L^\infty\) 范数。

此时 \(A\sim\delta_B/\mu\)，且
\(|\nabla\rho_B|_{r\sim\mu^{-1}}\sim\mu\delta_B\)。若不允许一个扩张的等测度
filament/reservoir，则还必须选择 \(\delta_B(t)\to0\) 的 vanishing-contrast 子窗口。
但爆破动力必须来自
与主 profile 同阶的 dynamic inner/\(D\) 层；Bessel 只能作为
\(r\gtrsim\mu^{-1}\) 的外响应。

因此最小的未排除结构不是三层，而是

\[
\boxed{
\text{非 Bessel 动态内核/壁层}
\longleftrightarrow
\text{收缩的 Bessel--Gaussian 外响应}
\longleftrightarrow
\text{全局反号补偿层/等测度储层}.
}
\]

完整的多尺度爆破方程、条件缩放律与待解的 Fredholm 符号问题记在
`analysis/BESSEL_MULTISCALE_BLOWUP_THEORY_ZH.md`。
