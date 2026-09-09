# Polar-profile equations centered at the wall point \((1,0)\)

This note is the active mathematical research contract for profile fitting.
It does not assume that the saved terminal fields have converged. In
particular, a good one-dimensional wall fit is not accepted unless the
two-dimensional polar field collapses at the same exponent.

## 1. Centering and pinning

We use the following sign convention for the physical/rescaled IPM system:

\[
\partial_\tau R
+\bigl(U+c_l\boldsymbol X+c_r\boldsymbol e_1\bigr)\cdot\nabla R
=c_\omega R,
\qquad
-\Delta\Psi=\partial_{X_1}R,
\qquad
U=\nabla^\perp\Psi
=\bigl(-\partial_{X_2}\Psi,\,\partial_{X_1}\Psi\bigr).
\]

Place the candidate singular point at \(a=(1,0)\) and write
\(z=\boldsymbol X-a\). A stationary point on the wall must satisfy

\[
c_r=-U_1(a)-c_l,
\]

because the total horizontal transport velocity must vanish at \(a\). Define
\(W=U-U_1(a)\boldsymbol e_1\) and subtract the density value at the point,
\(Q=R-R(a)\), so that a constant density background does not contaminate a
homogeneous cusp. The
centered steady equation is

\[
\bigl(W+c_l z\bigr)\cdot\nabla Q=c_\omega Q.
\tag{P1}
\]

This is the equation to fit. Merely translating a symmetric numerical field
from center zero to display center one does not introduce \(c_r\) or change the
equation.

## 2. Polar equations

Let
\[
z=(r\cos\theta,r\sin\theta),
\qquad 0<\theta<\pi,
\]
in the upper half-plane. Then

\[
\partial_{X_1}
=\cos\theta\,\partial_r-\frac{\sin\theta}{r}\,\partial_\theta,
\qquad
\Delta
=\partial_{rr}+\frac1r\partial_r+\frac1{r^2}\partial_{\theta\theta},
\]

and

\[
U_r=-\frac1r\Psi_\theta,
\qquad
U_\theta=\Psi_r,
\qquad
W_r=-\frac1r\Psi_\theta-U_1(a)\cos\theta,
\qquad
W_\theta=\Psi_r+U_1(a)\sin\theta.
\]

Consequently the steady profile system is

\[
\bigl(W_r+c_l r\bigr)Q_r+\frac{W_\theta}{r}Q_\theta
=c_\omega Q,
\tag{P2}
\]

\[
-\left(\Psi_{rr}+\frac1r\Psi_r+\frac1{r^2}\Psi_{\theta\theta}\right)
=\cos\theta\,Q_r-\frac{\sin\theta}{r}Q_\theta.
\tag{P3}
\]

The impermeable-wall condition gives
\(\Psi(r,0)=\Psi(r,\pi)=0\). It does not impose a Dirichlet condition on the
one-sided upper trace of \(Q\).

The full gradient \(g=\nabla Q\), rather than
\(\omega=\partial_{X_1}Q\) alone, obeys

\[
\bigl(W+c_l z\bigr)\cdot\nabla g
=\left[(c_\omega-c_l)I-(\nabla U)^{\mathsf T}\right]g.
\tag{P4}
\]

Thus \(\partial_{X_1}\rho\) is not a closed observable: stretching rotates and
couples its two components. This is why profile tests must retain both
\(Q_r\) and \(Q_\theta/r\) (or both Cartesian derivatives).

## 3. Symmetric and asymmetric problems

For reflection symmetry about the vertical line through \((1,0)\),

\[
Q(r,\pi-\theta)=Q(r,\theta),
\qquad
\Psi(r,\pi-\theta)=-\Psi(r,\theta).
\]

It is therefore sufficient to solve on \(0<\theta<\pi/2\), subject to

\[
\Psi(r,0)=\Psi\!\left(r,\frac\pi2\right)=0,
\qquad
Q_\theta\!\left(r,\frac\pi2\right)=0.
\]

The wall trace at \(\theta=0\) is evolved by (P2); it is not prescribed to
vanish.

The asymmetric problem must be solved on \(0<\theta<\pi\), with only
\(\Psi(r,0)=\Psi(r,\pi)=0\), together with pinning and amplitude
normalizations. The even part of \(Q\) generates the odd part of \(\Psi\),
whereas the odd part of \(Q\) generates the even part of \(\Psi\). No
condition may be imposed at \(\theta=\pi/2\). The asymmetric
problem should only be attempted after the symmetric residual and its odd
linearized mode have been measured.

## 4. Local cusp and the new angular equation

Suppose that, near the pinned point,

\[
Q(r,\theta)=r^\beta F(\theta)+\text{lower-order terms},
\qquad
\Psi(r,\theta)=r^{\beta+1}G(\theta)+\text{lower-order terms},
\qquad
0<\beta<1.
\]

Then

\[
\partial_{X_1}Q
=r^{\beta-1}\left[\beta\cos\theta\,F(\theta)
-\sin\theta\,F'(\theta)\right],
\qquad
\alpha=1-\beta.
\]

No numerical value of \(\alpha\) or \(\beta\) is supported by reliable data or
by a solvability condition. Any former numerical guess is therefore discarded;
below \(\beta\in(0,1)\) remains an unknown matching eigenvalue.
Substitution in (P3) gives the angular Poisson equation

\[
-G''-(\beta+1)^2G
=\beta\cos\theta\,F-\sin\theta\,F'.
\tag{A1}
\]

Near \(r=0\), nonlinear transport is of order \(r^{2\beta-1}\), whereas all
terms involving \(c_l\) and \(c_\omega\) are of order \(r^\beta\). Since
\(\beta<1\), the leading angular profiles must first satisfy

\[
-\beta G'F+(\beta+1)GF'=0.
\tag{A2}
\]

On every interval where \(G\) has a fixed sign, (A2) integrates explicitly:

\[
F=C\lvert G\rvert^{\beta/(\beta+1)}.
\tag{A3}
\]

Equivalently, the leading density is constant on streamlines,

\[
Q=C\lvert\Psi\rvert^{\beta/(\beta+1)}.
\tag{A4}
\]

Equations (A1) and (A3), together with the appropriate angular boundary
conditions, form the reduced nonlinear angular boundary-value problem. They
are the first place to search for a semi-analytic profile;
\(c_l,c_\omega,c_r\) enter only at
the next order after pinning.

## 5. Fixed-angle obstruction at a nonzero wall cusp

There is an important compatibility test. In a fixed angular sector, suppose
that \(F\) is continuous up to the wall, \(G\) is nonzero immediately in the
interior, and the streamfunction satisfies \(G(0)=0\). Equation (A3) then
forces \(F(0)=0\).
Therefore the coefficient of the wall singularity

\[
\partial_{X_1}Q(r,0)\sim\beta F(0)r^{\beta-1}
\]

vanishes. The same conclusion follows on an interval where \(G=0\) by using
(A1). Hence a continuous, single-scale separated cusp cannot simultaneously
have a nonzero \(r^{-\alpha}\) wall trace and satisfy the half-plane Dirichlet
condition for \(\Psi\).

This is not a nonexistence theorem for the observed blow-up. It narrows the
possible structures to at least one of:

1. an angular boundary discontinuity at \(\theta=0\);
2. a shrinking fan \(\theta_f(r)\sim r^\sigma\), so that
   \(Q\neq r^\beta F(\theta)\);
3. a second inner scale or boundary layer that changes the leading balance;
4. an intermediate, not asymptotic, wall power law.

The saved profiles visibly have radius-dependent angular cutoffs, so the next
ansatz to test is

\[
Q(r,\theta)=r^\beta H\!\left(\frac{\theta}{r^\sigma}\right).
\tag{M1}
\]

This includes \(\sigma=0\), the ordinary fixed-fan case.
`ipm_fit_polar_profile` therefore records \(\theta_f(r)\) at several
thresholds and fits \(\sigma\); a single-scale profile is accepted only when
the wall exponent, angular-\(L^2\)
exponent, and polar collapse agree.

## 6. Next-order selection

For a fixed-angle leading solution, the first correction that can balance the
modulation terms has degree one:

\[
Q=r^\beta F+rH+\cdots,
\qquad
\Psi=r^{\beta+1}G+r^2K+\cdots.
\]

It satisfies

\[
-K''-4K=\cos\theta\,H-\sin\theta\,H',
\tag{N1}
\]

\[
-G'H+(\beta+1)GH'-\beta K'F+2KF'
+(c_l\beta-c_\omega)F=0.
\tag{N2}
\]

Solvability of this linearized angular system, together with the pinning and
amplitude gauges, is the mechanism that can select the scaling-rate ratio.
For a shrinking fan, (M1) must be substituted into (P2)–(P3) before using
(N1)–(N2); the powers change, and \(\sigma\) becomes an additional selected
exponent.

## 7. An exact decaying Bessel family on punctured physical steady domains

An explicit family exhibits stretched-exponential far-field decay if center
regularity and homogeneous wall data are not required. This family solves the
unrescaled physical steady IPM equation, not the constant-rate rescaled
steady equation with \(c_l\neq0\). Introduce the centered coordinates

\[
X=x_1-1,
\qquad
Y=x_2,
\qquad
r=\sqrt{X^2+Y^2},
\qquad
\theta=\operatorname{atan2}(Y,X).
\]

Choose a real \(\lambda\neq0\) and impose the streamline relation

\[
\rho=\lambda\Psi.
\tag{B1}
\]

Then \(U\cdot\nabla\rho=0\) identically because
\(U=\nabla^\perp\Psi\). The Poisson equation reduces to

\[
\Delta\Psi+\lambda\,\partial_X\Psi=0.
\tag{B2}
\]

Writing

\[
\Psi=e^{-\lambda X/2}\Phi,
\qquad
k=\frac{\lvert\lambda\rvert}{2},
\]

transforms (B2) into the modified Helmholtz equation

\[
(\Delta-k^2)\Phi=0.
\tag{B3}
\]

Consequently, on the punctured upper half-plane, every separated mode

\[
\Theta_\nu(\theta)
=a_\nu\cos(\nu\theta)+b_\nu\sin(\nu\theta),
\]

\[
\Psi_\nu(r,\theta)
=e^{-\lambda r\cos\theta/2}K_\nu(kr)\Theta_\nu(\theta),
\qquad
\rho_\nu=\lambda\Psi_\nu.
\tag{B4}
\]

is an exact physical steady solution away from the tip; globally it is a
punctured-domain formula or a singular steady carrying the distributional tip
source audited below. Here \(K_\nu\) is the modified Bessel
function of the second kind, and \(\nu\geq0\) may be any real angular order when
no endpoint condition is imposed. Any convergent linear combination of modes
with the same \(\lambda\) remains exact because both (B1) and (B2) are retained.

Arbitrary \(a_\nu,b_\nu\) allow nonhomogeneous wall traces. If an impermeable
homogeneous-wall condition is required on \(0<\theta<\pi\), the special modes

\[
\Theta_n(\theta)=\sin(n\theta),
\qquad n=1,2,\ldots.
\tag{B5}
\]

vanish at both wall rays. This boundary specialization is optional for the
mathematical family.

For a right-facing acute wedge \(0<\theta<\theta_*<\pi/2\), the wall-fed
positive angular mode is

\[
\nu=\frac{\pi}{2\theta_*},
\qquad
\Theta(\theta)=\cos(\nu\theta).
\tag{B5a}
\]

It is nonzero on the right wall \(\theta=0\) and vanishes at the fan ray
\(\theta=\theta_*\). Choosing \(\lambda<0\) gives \(g_\lambda(0)=0\) and
\(g_\lambda(\theta)<0\) for every interior angle \(\theta>0\); the slow
algebraic tail is therefore aligned with the positive horizontal wall, and
the field decays
exponentially as it leaves that ray. The comparison mode

\[
\nu=\frac{\pi}{\theta_*},
\qquad
\Theta(\theta)=\sin(\nu\theta).
\tag{B5b}
\]

vanishes at both wedge rays. Both modes are exact in the open wedge. Extending
either one by zero outside the wedge is only a fan-domain/distributional
construction unless its normal derivative is also matched across the fan
ray; no global half-plane exactness is claimed for that zero extension.

For fixed \(0<\theta<\pi\), the large-radius expansion is

\[
\rho_\nu(r,\theta)
\sim
\lambda\sqrt{\frac{\pi}{\lvert\lambda\rvert}}\,
\Theta_\nu(\theta)r^{-1/2}e^{g_\lambda(\theta)r}
\left[
1+\frac{4\nu^2-1}{4\lvert\lambda\rvert r}
+O(r^{-2})
\right],
\]

\[
g_\lambda(\theta)
=-\frac{\lvert\lambda\rvert+\lambda\cos\theta}{2}.
\tag{B6}
\]

Thus this is precisely a stretched-exponential profile with

\[
\alpha=1,
\qquad
\beta=-\frac12,
\qquad
f(\theta)\propto\Theta_\nu(\theta),
\qquad
g(\theta)=g_\lambda(\theta)\leq0.
\tag{B7}
\]

For \(\lambda>0\), \(g_\lambda\) vanishes only on the left wall ray; for
\(\lambda<0\), it vanishes only on the right wall ray. A nonzero trace on
that slow ray still decays like \(r^{-1/2}\). A vanishing sine trace produces
a thin angular boundary layer and also tends to zero uniformly.

### 7.1 Global right-focused members

The most important global member takes \(\lambda=-\mu\), with \(\mu>0\),
has angular order zero, and uses the sign that makes the density positive.
For \(A>0\),

\[
\Psi_0(X,Y)=-Ae^{\mu X/2}K_0\!\left(\frac{\mu r}{2}\right),
\qquad
\rho_0(X,Y)=\mu Ae^{\mu X/2}K_0\!\left(\frac{\mu r}{2}\right),
\qquad
r=\sqrt{X^2+Y^2}.
\tag{B7a}
\]

It is an exact positive physical steady solution on the punctured plane:

\[
U\cdot\nabla\rho_0=0,
\qquad
-\Delta\Psi_0=\partial_X\rho_0.
\tag{B7b}
\]

Its restriction to the upper half-plane has a nonhomogeneous wall trace,
which is allowed for this research member. Its far field is

\[
\rho_0(r,\theta)
\sim A\sqrt{\frac{\pi\mu}{r}}\,
\exp\!\left[-\frac{\mu r}{2}(1-\cos\theta)\right].
\tag{B7c}
\]

For every fixed \(\theta>0\), the decay is exponential, whereas on the
positive horizontal ray it is \(r^{-1/2}\). Near that ray,

\[
\rho_0(r,\theta)
\sim A\sqrt{\frac{\pi\mu}{r}}\,e^{-\mu r\theta^2/4},
\qquad
\Delta\theta\sim(\mu r)^{-1/2}.
\tag{B7d}
\]

Thus the solution is globally defined rather than cut off at a fan boundary,
but it concentrates into a shrinking right-facing angular sector. For every
fixed acute angle \(\theta_*>0\),

\[
\frac{\rho_0(r,\theta_*)}{\rho_0(r,0)}
\sim
\exp\!\left[-\frac{\mu r}{2}(1-\cos\theta_*)\right]
\longrightarrow0.
\tag{B7e}
\]

If an impermeable homogeneous wall is wanted, the positive upper-half-plane
member is

\[
\Psi_1(r,\theta)
=-Ae^{\mu X/2}K_1\!\left(\frac{\mu r}{2}\right)\sin\theta,
\qquad
\rho_1(r,\theta)
=\mu Ae^{\mu X/2}K_1\!\left(\frac{\mu r}{2}\right)\sin\theta.
\tag{B7f}
\]

It vanishes at both wall rays and concentrates just above the positive wall,
with peak angle \(\theta_{\mathrm{peak}}\sim\sqrt{2/(\mu r)}\).

The pole can instead be placed below the physical domain. For \(h>0\), set

\[
r_h=\sqrt{X^2+(Y+h)^2},
\qquad
\Psi_h=-Ae^{\mu X/2}K_0\!\left(\frac{\mu r_h}{2}\right),
\qquad
\rho_h=\mu Ae^{\mu X/2}K_0\!\left(\frac{\mu r_h}{2}\right).
\tag{B7g}
\]

The modified-Helmholtz pole is then at \((0,-h)\). Hence
\((\rho_h,\Psi_h)\) is an exact smooth solution of the interior physical
equations throughout the
closed upper half-plane, with the same leading right-focused far field.
However, \(\Psi_h(X,0)\) varies with \(X\), so
\(U_2=\partial_X\Psi_h\) is nonzero on the wall. It is the restriction of a
whole-plane through-flow solution, not
an impermeable half-plane solution. It proves smooth interior geometry only;
it does not prove dynamical selection or time-dependent contraction toward
the singular boundary-pole member.

### 7.2 Dynamic-scaling status

Every physical solution can be represented in the dynamic coordinates, but
the global Bessel members are fixed profiles only for the zero-rate choice

\[
c_l=c_\omega=c_r=0.
\tag{B7h}
\]

Write \(\boldsymbol\xi=(X,Y)\). Because \(U\cdot\nabla\rho=0\), a
nonzero-rate fixed profile would have to satisfy

\[
c_l\boldsymbol\xi\cdot\nabla\rho+c_r\partial_X\rho=c_\omega\rho.
\tag{B7i}
\]

For (B7a),

\[
\boldsymbol\xi\cdot\nabla\log\rho_0
=-\frac12-\frac{\mu r}{2}(1-\cos\theta)+o(1),
\tag{B7j}
\]

which is not constant in \((r,\theta)\). After setting \(c_l=0\),
\(\partial_X\log\rho_0\) is also nonconstant. Thus (B7i) forces (B7h).

Under arbitrary coordinate scales it instead gives an exact nonstationary
rescaled orbit

\[
R(\boldsymbol\xi,\tau)
=C_\omega(\tau)\,
\rho_\mu\!\left(\frac{\boldsymbol\xi}{C_l(\tau)}\right),
\qquad
\mu_{\mathrm{rescaled}}(\tau)=\frac{\mu}{C_l(\tau)}.
\tag{B7k}
\]

The intrinsic Bessel length therefore drifts unless \(C_l\) is constant. This
family is a physical steady solution and a possible outer-profile model, not
yet a nontrivial dynamic-scaling attractor with \(c_l\neq0\). A smooth,
bounded IPM density also preserves its range for as long as the solution
remains classical, so convergence from smooth data to
the singular-density member can only be contemplated as an outer or
gradient-profile limit with a separate bounded inner core.

For a globally decaying rescaled steady profile there is also a necessary
mass condition. Let \(V=U+c_lz\) and suppose that
\(V\cdot\nabla Q=c_\omega Q\). Then

\[
\nabla\cdot(QV)=(2c_l+c_\omega)Q.
\tag{B7l}
\]

Integrating gives

\[
(2c_l+c_\omega)\int Q\,\mathrm d\boldsymbol X=0
\]
whenever the flux at infinity vanishes. Therefore a nonnegative, finite-mass
global profile with \(c_\omega=0\) and \(c_l\neq0\) is impossible. A
localized profile must instead use the mass-compatible rate
\(c_\omega=-2c_l\), have zero signed mass, retain a
nonintegrable outer tail, or be only an inner profile matched to a second
outer scale. This is a necessary condition, not an existence theorem.

The center is intentionally allowed to be singular:

\[
K_\nu(kr)\sim2^{\nu-1}\Gamma(\nu)(kr)^{-\nu}
\quad(\nu>0),
\qquad
K_0(kr)\sim-\log(kr).
\tag{B8}
\]

Hence the \(\nu=1\) homogeneous-wall member is a boundary-dipole-type solution.
Regularity at the center is not part of this exact-family contract. Requiring
simultaneously a regular center, homogeneous wall data, and decay at infinity
would leave only the zero solution within the linear relation (B1).

More generally, suppose an unrescaled physical steady state is sufficiently
smooth that the pressure, wall, puncture, and far-field boundary terms vanish,
and suppose that the Darcy dissipation is finite. Darcy's law and the stationary
potential-energy identity then give

\[
\int\lvert U\rvert^2\,\mathrm d\boldsymbol x
=-\int\rho U_2\,\mathrm d\boldsymbol x,
\qquad
0=\frac{\mathrm d}{\mathrm dt}
\int Y\rho\,\mathrm d\boldsymbol x
=\int\rho U_2\,\mathrm d\boldsymbol x.
\tag{B8a}
\]

Thus \(U=0\), hence \(\partial_X\rho=0\); horizontal decay forces the state to be
trivial. The admissible homogeneous-wall Bessel mode evades the identity
through its boundary pole, while the wall-fed modes additionally violate the
zero-normal-flow boundary contract. The singularity is not a removable
plotting defect.

The family records a concrete role for angularly dependent exponential
decay through \(g(\theta)\). A genuinely angle-dependent \(\alpha(\theta)\)
is not needed for (B4); it remains a separate WKB or sector-matching ansatz. Also,
the term \(c_l r\,\partial_r\rho\) rules out (B6) as the far field of a
nontrivial constant-rate rescaled steady state with \(c_l\neq0\).

`ipm_plot_exact_bessel_family` evaluates (B4) on a punctured polar wedge. It
uses the analytic derivative

\[
\rho_X
=\lambda A e^{-\lambda X/2}
\left[
-\frac{\lambda}{2}K_\nu\Theta
+k\cos\theta\,K_\nu'\Theta
-\frac{\sin\theta}{r}K_\nu\Theta'
\right].
\tag{B9}
\]

The routine reports transformed-Helmholtz Poisson, exact physical-transport, and
sector-edge normal-flow residuals. The default \(\theta_*=60^\circ\),
\(\lambda=-1\) plot compares (B5a) and (B5b), with signed logarithmic panels
exposing the nodal curve of \(\rho_X\) that is hidden by the center singularity
on a linear color scale. Setting \(\theta_*=\pi\) and orders \([0,1]\) compares
the global through-flow member (B7a) and the impermeable member (B7f) on the
actual
upper half-plane.

### 7.3 Parabolic two-scale limit and the smooth-data gap

The right-focused Bessel tail has a rigorous anisotropic blow-down that is
different from the single-scale orbit in (B7k). Write

\[
\xi=C_x(t)(x-a(t)),\qquad \eta=C_y(t)y,\qquad
R=C_\rho(t)\rho,\qquad
P=\frac{C_\rho C_y^2}{C_x}\Psi,
\]

and use the advective clock

\[
\frac{\mathrm dt}{\mathrm ds}=\frac{C_\rho C_y}{C_x^2}.
\]

With

\[
\delta=\left(\frac{C_x}{C_y}\right)^2,\qquad
c_j=\partial_s\log C_j,
\]

the exact rescaled equations are

\[
R_s+(-P_\eta+c_x\xi+c_r)R_\xi
 +(P_\xi+c_y\eta)R_\eta=c_\rho R,
\qquad
-(\delta P_{\xi\xi}+P_{\eta\eta})=R_\xi,
\qquad
\delta_s=2(c_x-c_y)\delta.
\tag{B10}
\]

For the Bessel branch, choose \(C_y^2=\mu C_x\). The exact relation
\(\rho=-\mu\Psi\) then becomes \(P=-R\). As
\(\varepsilon=\mu/C_x\to\infty\), the centered \(K_0\) and \(K_1\sin\theta\)
members converge on compact subsets of \(\xi>0\) to

\[
R_{0,*}=A\sqrt{\frac{\pi}{\xi}}
e^{-\eta^2/(4\xi)},\qquad
R_{1,*}=A\sqrt\pi\,\eta\xi^{-3/2}e^{-\eta^2/(4\xi)}.
\tag{B11}
\]

They solve \(R_\xi=R_{\eta\eta}\). Their nontrivial \(R_\xi=0\) curves are
\(\eta^2=2\xi\) and \(\eta^2=6\xi\), with finite-\(\varepsilon\) corrections

\[
\eta^2=2\xi-2\varepsilon^{-1}+O(\varepsilon^{-2}),\qquad
\eta^2=6\xi-6\varepsilon^{-1}+O(\varepsilon^{-2}).
\tag{B12}
\]

This heat equation is not the generic \(\delta=0\) IPM limit: the nonlinear
transport in (B10) remains order one unless the additional Bessel closure
\(P=-R\) is selected. For fixed physical \(\mu\), the limit requires
\(C_x,C_y\to0\), so both physical lengths expand and
\(L_y\sim L_x^{1/2}\). It is a downstream spatial far-field blow-down, not a
time-contracting cusp. A fixed external pole at depth \(h\) shifts the scaled
normal coordinate by \(C_yh=O(\varepsilon^{-1/2})\); without recentering this
is the leading convergence error.

A genuinely shrinking acute fan instead has

\[
\eta=y/x^p,\qquad p=1+\sigma>1,\qquad
\rho=x^\beta F(\eta),\qquad \Psi=x^\gamma G(\eta).
\tag{B13}
\]

The leading Poisson and transport balances give

\[
\gamma=\beta+2p-1,\qquad
-G''=\beta F-p\eta F',\qquad
\beta+\sigma=1,
\tag{B14}
\]

and an anisotropic fixed profile must satisfy

\[
-\beta G'F+\gamma GF'+(\beta c_x-c_\rho)F=0,
\qquad c_y=pc_x.
\tag{B15}
\]

Thus \(c_\rho=\beta c_x\) is only the special branch on which the leading
physical transport vanishes separately. More generally, a continuous nonzero
wall trace \(F(0)=F_0\) with \(G(\eta)=a\eta+\cdots\) is formally compatible
only if

\[
c_\rho=\beta(c_x-a),\qquad G''(0)=-\beta F_0.
\tag{B16}
\]

This supplies a concrete shooting problem, not an existence theorem. No
current theorem connects a standard unforced impermeable \(C^\infty\) IPM
initial datum to either (B11) or a solution of (B14)--(B16).

The reproducible MATLAB screens in experiments/bessel_two_scale_case.m
therefore use smooth near-Bessel data only as a falsifiable baseline. At
257-by-129 resolution through physical time 0.05, matched 3% perturbations
have K0 full/tail difference ratios 0.998510/0.998175 and K1 ratios
0.999992/1.00011. These measurements show short-time shape persistence, not
attraction. The full derivation and audit are recorded in
analysis/reports/IPM_BESSEL_TWO_SCALE_REPORT.tex.

### 7.4 Multipole correction to the former three-layer proposal

The earlier amplitude-only matching laws for a bounded core and the
small-argument Bessel singularity are not sufficient. Integrating the exact
Poisson equation over a circle gives, for the \(K_0\) coefficient,

\[
2\pi |A|+o(|A|)
\lesssim r\,\operatorname{osc}_{\partial B_r}\rho.
\tag{B17}
\]

For the impermeable \(K_1\) branch, Green's identity on a half-disk with the
harmonic test function \(h=Y\) gives

\[
\frac{|A|}{\mu}
\lesssim r^2\,\operatorname{osc}_{\partial B_r^+}\rho.
\tag{B18}
\]

When \(q=\mu r\ll1\), (B17) forces the core/closure correction to exceed the
\(K_0\) density component by at least \([q\log(1/q)]^{-1}\), while (B18)
forces the correction to exceed the \(K_1\) component by at least \(q^{-1}\).
Consequently,

\[
\boxed{
\ell\ll r\ll\mu^{-1}
\quad\Longrightarrow\quad
\text{the Bessel singularity cannot be the leading smooth-core layer.}
}
\tag{B19}
\]

This supersedes the former candidate laws
\(q_0\sim1/\log(1/\ell)\) and \(q_1\sim\ell\) as source-free gluing laws.
The Bessel family remains exact on the punctured/sourceful contract and useful, but its valid role is now the
drift-Helmholtz Green kernel and the outer response beginning at
\(r\simeq\mu^{-1}\). A dynamical construction must instead contain a
leading closure/wall layer, a shrinking \(\mu^{-1}\) Bessel--Gaussian outer
response, and a nonlocal source-neutral return/reservoir layer. The complete
proof, conditional fixed-point equations, and the general-\(\beta\) scaling laws
are in [BESSEL_MULTISCALE_BLOWUP_THEORY_ZH.md](BESSEL_MULTISCALE_BLOWUP_THEORY_ZH.md).
The subsequent exact directional-bridge families, the bounded-density
harmonic-strain example, and the nonlinear characteristic proof of local
convergence to a right-decaying Bessel exponential are in
[BESSEL_BRIDGE_EXACT_BLOWUP_ZH.md](BESSEL_BRIDGE_EXACT_BLOWUP_ZH.md). The same
report now also gives an exact orbit from bounded \(C^\infty\) density compactly
supported in \(x\), with local \(C^\infty\) convergence to \(e^{-X}\) and three
separated physical lengths; a fixed time-independent-strain long-time orbit
converging to a pure Bessel tangent; and wall-sine/rotating-strip regularizations
that identify the precise transverse-localization obstruction. These examples
use a nondecaying/infinite-energy harmonic far-field contract and do not
supersede the finite-energy/source-neutral obstruction (B19).

## 8. Numerical decision rules

- Do not infer an exponent from the \(R^2\) value of the wall-trace fit alone.
- Require a plateau of the local logarithmic slope over increasing annuli.
- Require agreement between the wall and angular-\(L^2\) radial exponents.
- Require collapse of \(r^\alpha\omega(r,\theta)\) when \(\sigma=0\), or of
  the corresponding \(\theta/r^\sigma\) curves for a shrinking fan.
- Compare symmetric and asymmetric residuals on the same compactified polar
  grid. Accept asymmetry only if it persists under reflection-symmetric grid
  refinement and lowers the full PDE residual.
- Treat the present saved endpoints as transient data until these checks are
  stable in time and resolution.
