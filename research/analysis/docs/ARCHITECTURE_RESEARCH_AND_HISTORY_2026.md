# IPM solver architecture

> Archived 2026-08-30. This document preserves the former solver and research
> narrative, including paths that no longer exist in the active tree. The
> maintained executable contract is [STRUCTURE.md](../../../STRUCTURE.md).
>
> **Retired historical contract annotation (2026-09-05).** Terms such as
> "active", "maintained", and "default" below describe the pre-schema-4
> solver at the date of the evidence. The current config/checkpoint schema is
> 4 with `scalingContract='exact_gauge_no_feedback_v1'`:
> `lengthScaleGain=1`, omega/width/travel restoring gains are zero,
> `maxDynamicRate` and the adaptive width controller are retired, and feature
> cell counts are diagnostics/remesh/hard-stop inputs only. Local annotations
> identify the old controller formulas while leaving them intact as evidence.

## Scope and active mode

The program solves the two-dimensional incompressible porous-media equation on
a rectangular truncation of the upper half-plane. The active
`symmetryMode='double_odd_omega'` retains the prescribed first-quadrant trace,
extends `omega=d_x rho` oddly across both coordinate axes, and uses the full
upper-half computational rectangle only to represent the paired `x_1` lobes.
`symmetryMode='half_plane'` retains the older one-sided formulation for
controlled comparisons. The maintained no-argument case is the symmetric
finite-degenerate strong-tail run
`rho_0=|x|^8/(1+|x|^8+y^8)`, with partial dynamic scaling, a fixed
wall-transport anchor at `X=1`, the wall-omega-peak amplitude gauge, nodal
nonuniform WENO transport, and the strict-`1.5` lattice mesh;
`ipm.config.activeCase.m` owns those overrides.
The generic completion defaults in
`ipm_defaults.m` remain the symmetric physical equation with
`c_l=c_omega=c_r=0`, so explicit option structures do not silently inherit an
experiment. Generic quantitative physical runs use the conservative far-field
choice `[-2000,2000] x [0,2000]` with `1025 x 513` points; the speed-oriented
active `k=8` case uses `513 x 257` on `[-1e5,1e5] x [0,1e5]`. A matched
`1025 x 513` rerun remains required for quantitative derivative claims.
Independent hyperbolic-sine maps in `x` and `y` choose their stretch
automatically so the first center/wall faces are approximately `0.01`; the
outer faces are deliberately coarse because Green data, rather than a false
zero condition, close the artificial boundary. `gridStretchAutomatic=true`
records that the stretch was derived from the domain and point counts; reusing
`result.opts` after changing either therefore recomputes the map. Set it false
only when `gridStretch` is intentionally prescribed numerically.
`rescalingMode='dynamic'` uses the common grid, Poisson solve,
velocity map, transport operator, and SSPRK3 step with two normalization
rates. Double-odd symmetry fixes the horizontal frame, so `c_r=0` and
`X_shift=0`; the comparison half-plane mode still supports traveling-wave
translation. `rescalingMode='physical'` is the exact no-`c` path and
uses adaptive physical coordinates only to represent contracting scales.
The opt-in `dynamicScaleGeometry='anisotropic'` is a correctness prototype
with separate horizontal and vertical coordinate factors. It is never selected
by default and therefore does not alter the maintained isotropic call graph.
Unless the caller explicitly requests adaptive remeshing, this geometry turns
it off: scale tracking, rather than repeated interpolation of the physical
field, is the primary two-scale mechanism.

Arrays are `ny`-by-`nx`: rows are `x_2=y`, columns are `x_1=x`.

## Physical mathematical contract

```text
rho_t + div(rho*u) = 0
-Delta psi = d_x rho
u = grad_perp psi = (-d_y psi, d_x psi)
```

This is equivalent to Darcy's law

```text
rho_t + u dot grad(rho) = 0
div(u) = 0
u + grad(p) = -(0,rho).
```

The rotation and sign convention are fixed as a single contract:
`-Delta psi=d_x rho` and `u=(-psi_y,psi_x)`. `ipmtests.baseline.suite` compares both signed
velocity components against a manufactured solution, so a sign change or a
90-degree rotation mistake fails the test.

In the active symmetry class,

```text
omega(-x_1,x_2)=-omega(x_1,x_2),
omega(x_1,-x_2)=-omega(x_1,x_2),
rho(-x_1,x_2)=rho(x_1,x_2),
rho(x_1,-x_2)=-rho(x_1,x_2).
```

Consequently `psi` is odd in both axes, `u_1=-psi_y` is odd/even, and
`u_2=psi_x` is even/odd in `(x_1,x_2)`. Thus `u(0,0)=0` and `u_2=0` on
`x_2=0`. The stored bottom row is the upper one-sided density trace; the user
explicitly accepts that its odd continuation has a jump at `x_2=0`, so the
upper trace is not overwritten by zero.

Production Green data integrate only first-quadrant source nodes and use the
four-image quadrant kernel

```text
G_Q = (1/(4*pi))*log((r_x^2*r_y^2)/(r_0^2*r_xy^2)),
r_0^2  =(x-xi)^2+(y-eta)^2,
r_x^2  =(x+xi)^2+(y-eta)^2,
r_y^2  =(x-xi)^2+(y+eta)^2,
r_xy^2 =(x+xi)^2+(y+eta)^2.
```

This is the odd image in each axis with the double-reflected positive image.
The older single-`x_2` half-plane kernel remains an explicit regression mode.
On the three artificial sides, production mode evaluates Green data before
the sparse interior solve.
When the retained source exceeds `greenMaxSources`, the code aggregates it in
uniform-reference grid blocks, separately for positive and negative strength,
and preserves each block's signed strength and physical centroid. It never
keeps only the largest sources and discards the rest.
`farBoundaryMode='dirichlet_zero'` is retained only for manufactured tests and
controlled comparisons. Production conclusions still require domain-
expansion checks.

Density transport uses the physical wall as a closed flux boundary. The left,
right, and top truncation sides are open by default and use the boundary nodal
state. This is necessary both for physical outflow and for a dynamically
expanding coordinate box. `transportBoundaryMode='closed'` is retained for
manufactured conservation tests; it is not used for dynamic dilation.
The default `wallTransportMode='conservative_flux'` evolves the bottom nodal
row as the half-control-volume average and therefore preserves the global
finite-volume mass contract. The opt-in
`wallTransportMode='advective_upwind'` instead applies the exact one-sided
trace equation `rho_t+u_1*rho_x=0` on that row. It removes the `O(h)` drift of
a smooth wall extremum at a stagnation point, at the cost of losing exact
finite-resolution mass conservation on one row. This option is reserved for
wall-trace convergence screens and must be compared across grids; it does not
silently alter the maintained production path.

## Initial-data contract

The no-argument maintained entry selects `initialCondition='degenerate_primitive'` with
`degeneratePower=8` through `ipm.config.activeCase.m`. The `ipm_defaults` fallback
selected by a nonempty partial option structure remains the same family with
`k=4`. The active first-quadrant upper trace is retained and extended
evenly in `x_1` so its derivative is odd:

```text
rho_0(x,y) = (1/k)*log((1+|x|^k+y^k)/(1+y^k)),
omega_0(x,y) = d_x rho_0
             = sign(x)*|x|^(k-1)/(1+|x|^k+y^k).
```

For the active `k=8` case, writing `D=1+|x|^8+y^8` gives the exact field

```text
rho_0 = (1/8)*log(D/(1+y^8)),
q_0 = omega_0 = sign(x)*|x|^7/D.
```

It has one sign on each open half-axis: `q_0>0` for `x>0` and `q_0<0` for
`x<0`. At the origin `q_0(x,0)~x^7`, so the datum has a finite seventh-order
degeneracy but no dead interval. On every generic interior ray its derivative
tail is `O(r^-1)`, and on the wall it is `q_0(x,0)~1/x`; in particular
`q_0 -> 0` in all far-field directions. The other gradient component is also
generically `O(r^-1)`, so the two-dimensional gradient `L2` energy has a
logarithmic far-field divergence. The active truncation is therefore
`[-1e5,1e5] x [0,1e5]`, and box comparisons remain part of every quantitative
claim. The PDE does not provide a maximum principle for `rho_x`; the solver
therefore records `positiveHalfPlaneOmegaMinimum` and
`positiveHalfPlaneOmegaNegativeRatio` instead of enforcing the initial sign by
projection. The active density grows like `log|x|` on the wall. This is the
necessary tradeoff for a one-signed inverse wall tail, since integrating a
positive `1/x` derivative cannot produce a bounded primitive.

The mathematical continuation of `rho` across `x_2=0` is odd, while the
stored `y=0` row is its nonzero upper trace `log(1+|x|^k)/k`. This is the
explicitly accepted jump at the symmetry axis. Point values on that axis do
not enter the Green integral; retaining the upper trace avoids replacing the
requested first-quadrant datum by a different wall datum. The wall-density
maximum therefore lies at the artificial right boundary and is not a useful
feature to pin. Initialization, `c_l` feedback, diagnostics, and video windows
track the positive-lobe wall maximum and level widths of `omega=d_x rho`;
the negative lobe is generated by symmetry. The former bounded rational datum
remains available as `initialCondition='degenerate'`. The compatibility option
`initialCondition='degenerate_density'` treats that older derivative profile
directly as an even upper density trace only for controlled comparisons.

`initialCondition='ccf_heavy_tail'` is an explicit historical/opt-in
regularized family for testing the CCF-motivated wall derivative. With
`z=(|x|-a)_+`, regularization `epsilon`, amplitude `A`, cutoff scale `R`, and
the upper-half-plane vertical envelope `V(y)=1/(1+(y/L)^q)`, define
`F=2*A*((z^2+epsilon^2)^(1/4)-epsilon^(1/2))`. The bounded family is

```text
S = 2*A*sqrt(R),
rho_0 = S*tanh(F/S)*V(y),
d_x rho_0 = sign(x)*A*z/(z^2+epsilon^2)^(3/4)*sech(F/S)^2*V(y).
```

It follows the proposed `1/sqrt(z)` tail for `z << R` and then saturates the
density monotonically instead of putting its maximum on the artificial box
boundary. Setting `R=Inf` recovers the uncut primitive, whose wall derivative
tends pointwise away from `|x|=a` to
`sign(x)*A*1_{|x|>a}/sqrt(|x|-a)` as `epsilon -> 0`. Every finite-`epsilon`
simulation has a resolved initial maximum. This separation is required:
the unregularized target is already singular initially and cannot by itself
provide evidence of finite-time blow-up. The named CCF experiment defaults are
`a=1`, `epsilon=0.08`, `A=1`, `R=16`, `L=2`, and even `q=4`.
Quantitative conclusions
still require both doubled-box and halved-`epsilon` comparisons; uncut
`R=Inf` tests additionally require explicit far-tail convergence. The
`x_1`-even stored density again gives odd
`omega` in `x_1`; its odd `x_2` continuation retains the accepted wall jump.
The active `k=8` case and the CCF experiment cases enable
`initialAnalyticRemesh`: before the first time step, the existing paired-peak
monitor proposes a feature-centered grid and the prescribed analytic density
is evaluated again on those new nodes. This
is deliberately different from interpolating a poorly resolved initial peak.
The option is false in the generic `ipm_defaults` fallback and true in
`ipm.config.activeCase`; in either case it changes only how the prescribed analytic
datum is sampled before evolution, not the PDE or solver path.

Because the accepted odd `x_2` continuation jumps, its full-plane
distributional `d_x2 rho` contains an axis-supported singular measure already
at initial time. The implemented PDE, derivatives, and growth diagnostics are
for the upper one-sided classical trace; they deliberately do not discretize
that symmetry-axis delta. Therefore full-plane gradient blow-up cannot be
inferred from this discontinuous extension. The tracked question is growth of
the regular first-quadrant/upper-trace gradients away from the prescribed jump.

The explicit CCF datum is separable, `rho_0=P(|x|)V(y)`, with `P(x)->S`. It
satisfies the intended all-direction contract `d_x rho_0 in C_0`: along a
fixed interior ray the vertical envelope decays, along a fixed-height
horizontal path `P'(x)->0`, and along mixed paths both effects apply. The full
gradient has a different far-field contract. For fixed `y>0`,
`d_y rho_0 -> S*V'(y)` as `|x|->Inf`; hence the vertical-gradient energy grows
linearly with the horizontal half-box, not logarithmically. At the historical
CCF `X=128` truncation, continuous separated quadrature gives `E_x=13.6589` and
`E_y=4591.6021`, so `99.7034%` of the initial full-gradient energy is the far
horizontal `d_y rho` reservoir. The global gradient-energy `c_omega` gauge is
therefore a valid reversible coordinate normalization but not a local
wall-`rho_x` amplitude or steady-state criterion. Claims about that opt-in
CCF family must
report the two energy components and use gauge-free physical wall
observables; a localized horizontal reservoir would be a different
mathematical initial condition and is not introduced silently.

In the former one-sided half-plane mode there is also a horizontal-frame
issue. Writing the far angular limit as
`rho_infinity=f(theta)` gives `d_x rho_infinity=g(theta)/r`. A leading
streamfunction ansatz `psi=r h(theta)` satisfies
`-(h''+h)=g`; the Dirichlet operator has the resonant `sin(theta)` mode. The
resulting truncated velocity contains a slowly box-dependent, nearly uniform
horizontal component (logarithmic in the resonant asymptotics). Thus the
unbounded velocity for this nondecaying datum requires a horizontal-frame
normalization and laboratory peak position is not a robust box-convergence
observable. Gradient norms and peak-aligned local shapes are translation
invariant. In active double-odd mode the paired `x_1` images cancel that
uniform horizontal frame component at the symmetry origin; `c_r` is therefore
identically zero. The solver does not silently cut off the datum or alter the
physical no-`c` equation.

## Dynamic-rescaling contract

> **Retired historical contract boundary.** The coordinate identities in this
> section remain useful derivations, but the pre-schema-4 rate limiter,
> adaptive-width correction, and nonzero amplitude/width/travel restoring
> formulas are not executable contracts. Current rates solve instantaneous
> exact gauge equations without those feedback terms; see
> [STRUCTURE.md](../../../STRUCTURE.md).

The implementation adapts the two-normalization dynamic-rescaling strategy
used by De Huang and collaborators to the homogeneity of IPM. Define

```text
X = C_l(tau) x + X_shift(tau) e_1
R(X,tau) = C_omega(tau) rho(x,t)
c_l = d_tau log(C_l)
c_omega = d_tau log(C_omega)
c_r = d_tau X_shift - c_l X_shift.
```

These formulas first define the canonical rescaled clock and rates. The
implemented `tau` coincides with this clock while `timeSpeed=1`; when the
common rate limiter activates it is the normalized clock described below.

Because the IPM velocity map `rho -> u` is homogeneous of degree zero in
space, the rescaled equation and physical clock are

```text
R_tau + (U + c_l X + c_r e_1) dot grad(R) = c_omega R
U = grad_perp (-Delta)^(-1) d_X1 R
dt/dtau = C_omega/C_l.
```

In code the rescaled field is still named `rho`. Do not copy the Boussinesq
factor `(c_l+c_omega)R` or its physical clock into this IPM equation.
The semi-discrete implementation sends the full combined velocity
`V=U+c_l X+c_r e_1` through one selected nonuniform conservative flux and
assembles

```text
R_tau = -div_X(R V) + (2*c_l+c_omega)R.
```

It does not discretize `c_l X dot grad R+c_r d_X1 R` with a separate centered
derivative. The former split form violated the density maximum principle by
about `1e-2` at physical time two. Open artificial fluxes plus the conservative
`2*c_l R` correction preserve a constant state to roundoff.
The methodological reference is Chen--Huang--Li, *Novel Self-similar
Finite-time Blowups with Singular Profiles of the 1D Hou--Luo Model and the
2D Boussinesq Equations* (arXiv:2604.01868). That paper supplies the gauge
strategy; the IPM equation and scaling law above are derived separately.

### Opt-in anisotropic two-scale geometry

Set `rescalingMode='dynamic'` and
`dynamicScaleGeometry='anisotropic'`. Define

```text
X = C_x(tau) x + X_shift(tau),     Y = C_y(tau) y,
R(X,Y,tau) = C_omega(tau) rho(x,y,t),
kappa = C_y/C_x = L_x/L_y,
P = C_omega*C_y*psi,
U = (-P_Y,P_X) = (C_omega*u_1, C_omega*kappa*u_2).
```

The exact elliptic and transport contracts are

```text
-div(diag(kappa^(-1),kappa) grad(P)) = R_X,
R_tau +(U_1+c_x X+c_r)R_X +(U_2+c_y Y)R_Y = c_omega R,
R_tau = -div(R V)+(c_x+c_y+c_omega)R,
V=(U_1+c_x X+c_r,U_2+c_y Y),
dt/dtau = C_omega/C_x,
kappa_tau=(c_y-c_x)kappa,
X_shift_tau=c_x X_shift+c_r.
```

Every SSPRK3 stage passes its own `kappa` to the Poisson/Green/velocity map;
using only the step-start aspect is not consistent at third order. The
discrete operator is assembled exactly as

```text
A_kappa = kappa^(-1) kron(T_x,I_y)+kappa kron(I_x,T_y).
```

Its artificial-boundary image kernel replaces every Euclidean squared
distance by `kappa*dX^2+dY^2/kappa`. The determinant is one, so the logarithmic
kernel coefficient and source quadrature weight do not acquire another
prefactor. `anisotropicPoissonSolver='direct'` is the small-grid correctness
default. The solve is isolated behind `ipm_anisotropic_poisson_solve`; the
opt-in `pcg` path uses the existing isotropic factor as a preconditioner, and
the interface can later accept an aspect-keyed factor/PCG cache without
changing callers.

> **Retired historical contract (nonzero travel restoring).** The following
> `travelingWaveGain` formula is retained to explain the old anisotropic
> comparison. Schema 4 fixes that gain to zero; no peak-position restoring
> term may be added to `c_r`.

The first manufactured-evolution gauge is explicit:
`anisotropicFixedCX`, `anisotropicFixedCY`,
`anisotropicFixedCOmega`, and `anisotropicFixedCR` specify the four canonical
rates. `anisotropicGaugeMode='fixed'` uses all four directly. For a half-plane
run, `anisotropicGaugeMode='peak_translation'` keeps `c_x,c_y,c_omega` fixed
but recomputes

```text
c_r = -U_1(X_*,0)-c_x X_*
      -travelingWaveGain*|c_x|*(X_*-X_pin).
```

Thus the existing peak-translation and `X_shift` frame are retained without
pretending that one horizontal condition determines two rates. At a fixed
state the constraint Jacobian is the one-row matrix `[X_*,1]` acting on
`[c_x,c_r]`: its joint rank is one and its joint condition number is infinite.
After `c_x` is supplied by another gauge (or prescribed), solving for `c_r`
has reduced condition number one. These rank, coefficients, raw peak velocity,
and controlled residual are saved. Double-odd symmetry still enforces
`c_r=0`; `peak_translation` is rejected there.
This opt-in `peak_translation` branch intentionally remains a transport-lock
compatibility mode; it does not claim to freeze the exact quadratic maximum.
`fixed` and `peak_translation` are the only supported anisotropic gauge
modes. Option resolution rejects every other value before numerical setup.

Physical reconstruction is axis-specific:

```text
x=(X-X_shift)/C_x,       y=Y/C_y,
rho=R/C_omega,           psi=P/(C_omega*C_y),
u_1=U_1/C_omega,         u_2=U_2/(C_omega*kappa),
rho_x=(C_x/C_omega)R_X,  rho_y=(C_y/C_omega)R_Y,
mass=integral(R)dX dY/(C_omega*C_x*C_y).
```

Consequently the physical gradient energy is
`C_omega^(-2)[(C_x/C_y)E_X+(C_y/C_x)E_Y]`; the implementation records both
axis pieces and recomputes the anisotropic effective length from the full
physical gradient. The local Casimir identity is
`d_tau[J Phi(R/C_omega)]+div[J V Phi(R/C_omega)]=0` with
`J=(C_x C_y)^(-1)`. Since physical `rho_x` is only the finite scalar
`C_x/C_omega` times `R_X` at each preterminal time, the required
all-direction `R_X -> 0` tail is preserved by the coordinate change; the
prototype does not install a radial cutoff or a physical remap at every step.

This geometry is not an extra IPM scaling symmetry. A finite positive-kappa
extended fixed point must have `c_x=c_y`; unequal limiting rates imply an
anisotropic relative orbit/nonautonomous family, or a degenerate
`kappa -> 0/infinity` limit. The fixed-rate and peak-translation paths verify
evolution and continuation infrastructure. They do not
by themselves establish convergence to a new steady profile; that requires
conditioned long-time rate/profile convergence and grid/box replication.

> **Schema-4 reading note.** In the historical list below, every term
> proportional to `k_r` is retired because `travelingWaveGain=0`. The exact
> no-feedback identities are the zero-target equations obtained without that
> restoring term.

The three maintained rate conditions are:

1. In active double-odd mode, set the fixed computational wall point
   `X_a=transportAnchorX=1` to be a stagnation point of the combined transport:

   ```text
   U_1(X_a,0)+c_lNominal*X_a=0,
   c_lNominal=-U_1(X_a,0)/X_a.
   ```

   `X_a` is independent of `pinX`, which remains the initial wall-`omega` peak
   used by tracking and remeshing. The fixed-point row has scalar response
   `X_a=1` and is therefore nonsingular. This condition supplies only a
   coordinate stagnation point; its dynamical stability must be audited from
   `d_X(U_1+c_l X)`. For the exact active `k=8` logarithmic primitive the
   zero-step canonical derivative at `X=1` is approximately `-0.182`, so the
   prescribed anchor is attracting. The next positive-side zero is near
   `X=1.333` and is repelling, initially delimiting the anchored transport
   basin. The amplitude gauge only applies a common positive time
   reparameterization and cannot change these signs. This is an instantaneous
   structural check, not yet a proof that the basin or a two-scale profile
   persists dynamically.
2. In active double-odd mode impose `c_r=0`, because horizontal translation
   would destroy the symmetry center and `u(0,0)=0`. In comparison half-plane
   mode, locate the tracked wall peak to subgrid accuracy and restore it with
   `c_r=-U_1(X_*,0)-c_l X_*-k_r |c_lNominal|(X_*-X_pin)`.
3. The maintained `cOmegaGauge='wall_omega_peak'` sets the semi-discrete rate
   of the tracked positive wall-`omega` maximum to zero. The historical
   CCF/gradient-energy gauge and other amplitude gauges remain explicit
   comparisons and do not alter the fixed transport-anchor equation.

> **Retired historical contract (position restoring).** The later sentence in
> this paragraph that locks the tracked feature with
> `-k_r|c_lNominal|(X_*-X_pin)` describes the old nonzero travel gain; schema 4
> permits no such correction.

For `lengthGauge='transport_anchor'`, the active canonical rate is

```text
c_l = c_lNominal = -U_1(1,0).
```

Validation requires `lengthScaleGain=1`, and grid-cell feedback cannot add a
width correction; otherwise the exact zero would be lost. Grid-cell counts
only request a remesh. In double-odd mode symmetry fixes `X_shift=0`. Only the
comparison half-plane mode uses translation; there `c_r` gives the tracked
feature the restoring transport velocity
`-k_r|c_lNominal|(X_*-X_pin)`, so its center is locked while its width is
controlled. The active exact residual is `U_1(1,0)+c_l=0`, together with
`U_1(0,0)=c_r=X_shift=0`. The common positive `timeSpeed` multiplies the full
canonical vector field, so the bounded field still satisfies
`V_1^(tau)(1,0)=0`. The CFL uses the full velocity `U+c_l X+c_r e_1`.

`lengthGauge='transport_anchor'` is the active double-odd default.
`lengthGauge='omega_peak_location'` remains an explicit alternative that
freezes the actual wall-derivative extremum equation; `symmetry_peak` remains
the origin-to-moving-peak passive secant comparison. Neither alternative is
silently combined with the fixed `X=1` condition.
`lengthGauge='local_strain'` remains the automatic half-plane default. The
explicit comparison `lengthGauge='wall_density_width'` instead sets
`c_lNominal=-(U_1(X_0.6)-U_1(X_0.2))/(X_0.6-X_0.2)`, where `X_q` is the first
wall crossing of `R/q_max=q`. The ratio cancels the amplitude scaling and the
velocity difference cancels a uniform frame drift. It is not the default:
stress testing showed that preserving this broader material-density width did
not preserve the much narrower 90% `R_X` core in a multiscale run.

> **Retired historical contract (adaptive width controller).** The complete
> controller and formulas below are retained as provenance only. Schema 4 has
> no `adaptiveLengthScaling`, `maxWidthRateCorrection`, or `adaptiveGain`
> runtime path. These cell counts may only support telemetry, remeshing, or a
> hard stop; they cannot correct `c_l`.

In half-plane comparison modes, `adaptiveLengthScaling=true` can add a rate
correction that opens the positive
wall-`omega` peak when it contains too few actual nonuniform cells and retracts
it before outer levels expand indefinitely. At
each `10%`, `50%`, and `90%` level, only the connected component containing
the tracked positive peak is counted. Boundary-crossing cells receive
their linearly interpolated fractional count. The normalized area count is
computed on the connected `10%` component. Thus a long tail or a separated
secondary scale cannot make the locked peak look adequately resolved. This
extracts the old
`delta0/delta1/delta2` mesh criterion into a stateless scaling feedback:

```text
adaptive targets = min(absolute targets,
                       targetPeakExpansion*initial connected counts)
areaFactor = adaptiveTargetPeakPoints/connectedAreaGridPoints
levelFactor = max(adaptiveTargetLevelPoints./levelGridPoints)
resolutionFactor = max(areaFactor,levelFactor)
expansionRatio = max(current connected counts/initial connected counts)
safetyFactor = max(minimumPeakPoints/connectedAreaGridPoints,
                   minimumLevelPoints./levelGridPoints)

if safetyFactor > 1:
    dimensionlessCorrection = widthExpansionStrength*log(safetyFactor)
else if expansionRatio > widthContractionOnset:
    dimensionlessCorrection = widthContractionStrength
                              *log(widthContractionOnset/expansionRatio)
else:
    headroom = log(widthContractionOnset/expansionRatio)
               /log(widthContractionOnset)
    dimensionlessCorrection = widthExpansionStrength
                              *max(log(resolutionFactor),0)*headroom

dimensionlessCorrection = clamp(dimensionlessCorrection,
                -maxWidthRateCorrection,maxWidthRateCorrection)
widthCorrection = |c_lNominal|*dimensionlessCorrection.
```

Outside the hard safety floor, outer expansion has priority over an
under-resolved inner core. Once the connected area or a level falls below the
absolute safety counts `12` and `[16,8,4]`, emergency opening (`mode=2`)
overrides contraction. The ordinary positive response fades continuously to
zero at the contraction onset. Defaults use targets at most `1.5` times the
initial counts, absolute target caps `32` and `[64,32,16]`, onset `1.35`,
dimensionless expansion strength `1`, contraction strength `28`, correction
cap `16`, and baseline gain one. This rate feedback is deliberately disabled
for the active double-odd extremum-location gauge; its cell counts only drive
the discrete remeshing fallback.

The reconstructed physical nodes `(X-X_shift)/C_l` always move as `C_l`
changes. In addition, `adaptiveRemesh=true` provides a bounded local fallback
when the configured safety trigger is exceeded. In physical mode one monitor
redistributes `x_1` nodes around both `+x_*` and `-x_*`, while a second
one-sided monitor redistributes `x_2` nodes around the wall. Each axis is
triggered by its own area/support/core deficit, and both may move in one
accepted remesh. Physical density is
transferred by pure range-preserving PCHIP so the accepted one-sided wall trace
and local derivative shape are not altered by a quadrature-mass correction;
remap mass drift is exposed diagnostically. The dynamic comparison path keeps
its bounded horizontal row-integral correction. Arbitrary-grid
derivatives and the Poisson factorization are rebuilt. The axis generator uses
all connected 10/50/90-percent crossing intervals. Monitor weights are
proportional to requested cells per physical width, so the narrowest scale is
resolved without spending most nodes on the outer support. For the paired
horizontal feature, a graded origin-to-peak bridge prevents the coarse central
hole produced by isolated peak monitors. Local auxiliary nodes resolve
sub-base-grid intervals before inversion of the monitor CDF. Every request
starts from the original stretched grid, rather than the previous adaptive
grid, so concentration cannot compound artificially. A logarithmic amplitude
search uses the weakest deformation that reaches buffered safety `0.4`;
minimum spacing and a production adjacent-cell-ratio cap of `6` are hard
constraints. The trigger is `0.75`, ahead of the strict safety boundary. A
redistribution is accepted only if it reduces the relevant local spacing by at
least 2% and, after rebuilding `D_x`, does not reduce either connected 90%
core count by more than 2% on a moved axis. On the unchanged cross-axis, the
floor is capped by its already configured target, so a more accurate peak
index cannot reject a candidate while that cross-axis remains above target.
Physical remaps must additionally change the wall
`omega` maximum by no more than 3%; otherwise the old field and operators are
retained. The minimum spacing is `1e-4` times the original center target and at
most 300 changes are accepted. Domain endpoints, `x_2=0`, and the exact
`x_1=0` node are retained. History records minimum spacings, the origin gap,
maximum adjacent-cell ratios, and diagnostic global max/min spacing ratios in
both axes. Only adjacent ratios are admissibility, stopping, and trusted-fit
contracts; global ratios of thousands or more are allowed.
`remeshOuterAnchorFactor=Inf` keeps this established global redistribution as
the default. CCF long-time cases set it to `3`: the candidate is built only on
the reference subgrid ending beyond the outermost connected monitor interval
plus three of its widths, and is spliced back into unchanged far nodes. Once
the existing adjacent-cell ratio exceeds 90% of its hard cap, two halo nodes
are added at each splice and only the exterior transition cell widths are
smoothed to the soft ratio `sqrt(hard cap)`. Feature cells are not smoothed.
The hard cap remains unchanged, and all nodes outside the halo anchor remain
bitwise fixed. The anchored candidate is still subject to monotonicity,
minimum-spacing, adjacent-ratio, connected-count, peak-transfer, and
paired-symmetry checks.
The default `remeshReferenceMode='original'` remains non-compounding. CCF
long-time cases opt into `current`: each new local monitor starts from the
last accepted grid, while the outer anchor fixes the far nodes. This permits
successive representation of nested core scales that are absent from the
original sinh grid. It is not a relaxation of acceptance: minimum spacing,
cell ratio, density range, paired symmetry, core preservation, and physical
peak-change checks apply after every transfer.
If the first monitor amplitude violates the cell-ratio bound before reaching
all targets, the amplitude search now bisects back to the strongest valid
incremental improvement instead of returning the unchanged grid. This is
essential for `current` mode, where several accepted local changes may be
needed to represent a nested scale; every individual change still passes the
same hard checks.
`remeshAxisProfile='lattice'` is the concentrated alternative used by the CCF
long cases. It extracts the useful contract from the legacy
`latgene/latgene_diff` family: a nearly uniform fine plateau on each connected
level interval and a narrow smooth logistic transition to the far grid.
Nested plateaus combine by a maximum, so broad support levels do not dilute
the inner core. The legacy coordinate/metric inconsistency in its additional
cubic correction is intentionally not copied. The long profile uses
transition fraction `1/8`, paired bridge floor `0.005`, bridge power `4`,
target point caps `128` and `[192,128,96]`, and initial expansion cap `8`; the
buffered remesh safety is `0.08`, so the innermost level receives its full
96-cell candidate reserve instead of merely four times the safety floor. The
outer anchor is `3`, the proactive trigger is `0.15`, and the hard adjacent
ratio is `1.5`. Up to 12 analytic pre-remesh passes are allowed; each pass
re-evaluates the prescribed formula rather than interpolating the preceding
samples. Global spacing ratios remain unconstrained.
> **Retired historical contract (rate cap).** The `0.5` common-rate cap in the
> next paragraph belongs to the removed `maxDynamicRate` normalization. It is
> retained only to identify the historical run.

The active speed-oriented `k=8` run uses the same lattice construction on
`513 x 257`, but lowers the peak/level targets to `64` and `[96,64,48]`, the
hard floors to `12` and `[16,10,8]`, and the analytic pre-remesh cap to six
passes. Its deliberately partial, robust modulation is:
the fixed condition `U_1(1,0)+c_l=0` supplies `c_l`, the tracked positive
wall-omega-peak condition supplies `c_omega`, and double-odd symmetry supplies
`c_r=0`. Its common rate
reparameterization cap is `0.5`. Both dynamic and physical modes compute the
wall-normal connected core and may remesh `x_2`; modulation never disables the
two-axis safety contract.
This keeps the useful multilevel criterion from `Data.admesh/findcri` without
importing its stateful Boussinesq call graph. Remeshing is a numerical fallback,
not a substitute for grid/domain convergence.

Write the rescaled right-hand side without its amplitude source as `B`, so
`R_tau=B+c_omega R`. The active `cOmegaGauge='wall_omega_peak'` uses the
tracked positive wall maximizer `X_*` and imposes

```text
d_tau[d_X1 R](X_*,0) = d_X1 B(X_*,0)+c_omega*d_X1 R(X_*,0) = 0
c_omega = -d_X1 B(X_*,0)/d_X1 R(X_*,0).
```

The denominator is the positive peak itself, rather than a cancellative
strain observation. The code records its location, its size relative to the
global `|d_X1 R|` maximum, and the peak-rate residual. The global wall maximum
remains a separate diagnostic, so a tracked/global mismatch is visible rather
than silently changing the normalization target.

The comparison `cOmegaGauge='gradient_energy'` preserves

```text
E_grad = integral (|d_X1 R|^2+|d_X2 R|^2) dX
c_omega = -<grad R,grad B>/E_grad.
```

This comparison uses both density derivatives and all resolved scales; its
denominator is positive unless `R` is spatially constant. The former
`cOmegaGauge='strain_point'` preserves `d_X1 U_1(0,0)` and is retained only as
an explicit comparison. Its condition number
`|d_X1 U_1(0,0)|/||d_X1 U_1||_infinity` is recorded because cancellation or a
zero crossing creates a gauge pole. `cOmegaGauge='none'` sets `c_omega=0`, the
natural amplitude gauge suggested by the IPM density maximum principle.

The code records the fixed transport-anchor coordinate, canonical and bounded
anchor velocities, comparison peak-gauge residuals and conditioning,
gradient-gauge residual, point-strain condition, nominal and corrected rates,
pin residuals, wall-peak width, and multiscale gradient
diagnostics. The latter include total gradient energy, an effective gradient
length, and disjoint energy fractions inside half the wall-peak scale, between
half and twice that scale, and outside twice that scale. The
`farBoundarySourceRatio` reports the largest outer-band `|omega|` relative to
its domain maximum. `farBoundaryVelocityRatio` retains the raw speed ratio and
`frameVelocityU1` records the tracked horizontal speed. For the nondecaying
angular datum the raw velocity ratio is not expected to vanish and is not used
as an automatic tolerance; neither quantity replaces a peak-aligned
doubled-box comparison.

The physical reconstruction from a dynamic result is

```text
x_physical = (X-X_shift e_1)/C_l
rho_physical = R/C_omega
u_physical = U/C_omega
grad_x rho_physical = (C_l/C_omega) grad_X R
mass_physical = integral(R)dX/(C_l^2*C_omega).
```

The distinction is mandatory in diagnostics. `flow.source`, `result.omega`,
`history.wallPeak`, and `history.trackedWallPeak` are computational
`d_X1 R`; the first is global and the second remains a tracked diagnostic. The
active fixed-transport gauge does not keep that derivative peak stationary;
its drift relative to `X=1` must be reported. `result.physicalOmega`,
`history.physicalOmegaInf`, `history.physicalWallOmegaPeak`, and
`history.physicalTrackedWallOmegaPeak` multiply by `C_l/C_omega` and are the
quantities used to judge physical gradient growth.
The live physical-omega field and wall trace display this reconstructed
derivative on the moving peak-aligned physical grid; the first panel displays
rescaled `R_X1` on the fixed computational `X` grid. Titles report physical
and rescaled peaks, the FWHM, and the connected 90% core width.
For the active double-odd solution, live/video and final field plots pass only
the independent physical side `x_1>=0` to the graphics objects. Peak-relative
coordinates may still be negative to the left of the positive peak, but they
never represent the mirrored physical half `x_1<0`. This is an observation
restriction only: the current transport and Poisson solve still retain the
full paired `[-L,L]` grid.
`physicalTrackedPeakX` reports the laboratory-frame peak position, whereas
`comovingPhysicalX=(X-X_*)/C_l` supplies a peak-aligned physical coordinate for
domain/refinement comparisons.

### Dynamic-time normalization

> **Retired historical contract (`maxDynamicRate`).** This entire subsection
> documents the removed pre-schema-4 time reparameterization. Schema 4 directly
> integrates canonical `tau`, has no configurable common rate cap, and does not
> apply the factor `a` below. The equations and measurements remain here only
> as historical evidence.

The peak, width, and position gauges first produce canonical rates
`hat(c_l)`, `hat(c_omega)`, and `hat(c_r)`. Simultaneous growth of the first
two contains a common time-parametrization component; clipping them
independently would change the rescaled trajectory. The solver instead uses a
single positive speed

```text
M = max(|hat(c_l)|, |hat(c_omega)|,
        |hat(c_r)|/max(|X_pin|,1))
a = min(1,maxDynamicRate/M)
c_l = a*hat(c_l),  c_omega = a*hat(c_omega),  c_r = a*hat(c_r)
R_tau = a*R_hat_tau
dt/dtau = a*C_omega/C_l.
```

The default `maxDynamicRate=2` therefore bounds the evolved `c_l` and
`c_omega` without changing the orbit in field/scale space. The same `a`
multiplies physical transport, dilation, translation, amplitude source, all
scale ODEs, and the physical clock. `history.c_l/c_omega/c_r` are the bounded
evolution rates; `history.canonicalCL/canonicalCOmega/canonicalCR` retain the
unmodified gauge rates, and `history.timeSpeed` records `a`. This is a time
reparameterization, not evidence that canonical multiscale rate growth has
disappeared.
`history.transportAnchorX`, `canonicalTransportAnchorVelocity`, and
`transportAnchorVelocity` replay the fixed-point condition without allowing a
small `timeSpeed` to conceal a canonical residual.

The solver integrates a separate canonical clock with
`d(canonicalTau)/d(normalizedTau)=a`. User options `finalTime`, `outputEvery`,
and `maxDt` are canonical-time quantities; the internal normalized step is
divided by the current `a` and remains subject to the full CFL. Thus changing
only `maxDynamicRate` changes work and normalized elapsed time, not the target
point on the canonical rescaled trajectory. `result.tau` is canonical,
`result.normalizedTau` exposes the internal clock, and both clocks are stored
in history/snapshot metadata.
`physicalFinalTime` optionally terminates at an exact reconstructed physical
time; if it is supplied without `finalTime`, the canonical limit is set to
infinity. The last step is restricted by both clocks. This is the required
path for gauge, grid, and box comparisons at a common physical instant.

## Data flow

```text
options -> operators -> initial rho -> optional gauge initialization
                                      |
                                      v
rho -> d_x rho -> Green far data -> Poisson -> signed divergence-free velocity
 |                                               |
 +---- selectable full-velocity conservative flux <----+
                         |
              optional (2*c_l+c_omega) source
                         |
                         v
       SSPRK3 for rho, log scales, physical time, and canonical time
                         |
             2-D safety audit -> optional bounded axis remesh
                         |
                         v
             diagnostics / snapshots / result
                         |
                         v
                  live figure / MP4
```

## Public interface

- `ipm.solve(userOpts)` is the only simulation entry point.
- `ipm.solve()` and `ipm.solve([])` use the stateless overrides returned by
  `ipm.config.activeCase`; explicit structures retain their stated configuration.
- `ipmtests.baseline.suite()` is the verification entry point.
- `practical_main.m` is a minimal wrapper around `ipm.solve`.
- `ipm_resolve_options(userOpts)` completes, derives, and validates one frozen
  flat option structure for a solver run. `ipm_defaults(userOpts)` is the
  temporary legacy flat-option wrapper.

## Internal responsibilities

- `ipm.config.activeCase.m`: maintained no-argument `k=8` rational-tail,
  dynamic/lattice option overrides, including
  `lengthGauge='transport_anchor'` at `X=1`,
  `cOmegaGauge='wall_omega_peak'`, and box `[-1e5,1e5] x [0,1e5]`; it contains
  no solver state or alternate time integrator. The CCF configurations live in
  named experiment factories and are opt-in.
- `ipm.mesh.build.m`: pure construction of a uniform or
  center-clustered stretched grid from already resolved options,
  reference-coordinate derivatives divided by analytic mapping Jacobians,
  matching control-volume/Jacobian weights, and reusable weighted sparse
  `-Delta` factorization; it also stores the separated `T_x,T_y` factors used
  by the opt-in anisotropic operator.
- `ipm_build_operators.m`: temporary compatibility wrapper that resolves legacy
  partial flat options before calling the pure builder. The `ipm.solve` and
  remesh paths do not use this wrapper.
- `ipm.field.greenBoundary.m`: signed-strength/centroid-preserving block compression
  and quadrature of the active four-image quadrant Green formula (or explicit
  legacy half-plane kernel) on the three artificial boundaries, with the exact
  determinant-one anisotropic metric when an aspect is supplied.
- `ipm_anisotropic_poisson_operator.m`: assemble exact `A_kappa` and its
  weighted symmetric form.
- `ipm_anisotropic_poisson_solve.m`: isolated direct correctness solve or
  opt-in isotropically preconditioned PCG bridge for larger two-scale screens.
- `ipm.field.poisson.m`: apply boundary data and perform one isotropic or
  aspect-aware interior Poisson solve.
- `ipm.field.biotSavart.m`: source-to-velocity map.
- `ipm.field.velocity.m`: complete IPM map `rho -> d_x rho -> u`.
- `ipm.field.initialDensity.m`: built-in or user-supplied density.
- `experiments/ccf_heavy_tail_case.m`: physical/dynamic, box, and onset-
  regularization comparisons for the CCF-like tail; it only returns options
  to `ipm.solve` and is not another solver entry point.
- `experiments/ccf_gaussian_trough_case.m`: smooth Schwartz wall-trace
  screens with exact initial self-strain, pointwise wall transport, and a
  refined case; it only returns options.
- `ipm.evolve.initializeScaling.m`: fixed normalization indices and target.
- `ipm.evolve.adaptiveGain.m` **[retired historical component; no maintained
  destination]**: additive outer-priority connected-width feedback for `c_l`.
- `ipm.remesh.adapt.m`: bounded local `x_1`/`x_2` redistribution,
  range-preserving physical state transfer, optional bounded dynamic-row
  conservation, and operator rebuild.
- `ipm.remesh.axis.m`: stateless nested-level monitor, graded symmetry
  bridge, target-count amplitude search, far anchoring, on-demand exterior
  transition smoothing, and grid-smoothness constraints.
- `ipm.diagnostics.peakLocation.m`: three-point quadratic subgrid maximum on a nonuniform
  line.
- `ipm.diagnostics.peakResolution.m`: connected-primary-peak resolution on the actual
  nonuniform wall grid, including interpolated level-crossing bounds.
- `ipm.field.transport.m`: selectable conservative nonuniform flux, closed on
  the physical wall and configurable on the three artificial sides; an
  explicit opt-in replaces only the bottom row by the one-sided advective
  trace update. MUSCL is the fallback and the maintained case uses nodal WENO
  on complete stencils.
- `ipm.field.weno5Geometry.m`: precompute left/right five-node Lagrange,
  candidate, divided-difference, and smoothness geometry on each current axis.
- `ipm.field.weno5Reconstruct.m`: vectorized scale-invariant nonlinear WENO blend
  for independent nodal grid lines.
- `ipm_rescaling_rates.m`: `c_l`, `c_omega`, `c_r`, the fixed wall-transport
  anchor, and gauge residuals.
- `ipm_anisotropic_transport_rhs.m`: the anisotropic rate-dependent
  semi-discrete transport/source assembly.
- `ipm.evolve.rhs.m`: physical velocity, common transport, and optional rescaling
  source assembly.
- `ipm.evolve.stepSsprk3.m`: field, scale, physical-clock, and canonical-clock
  integration; the anisotropic branch passes the stage-local aspect.
- `ipm.evolve.timeSpeed.m`: common trajectory-preserving dynamic-time speed.
- `ipm_time_speed_anisotropic.m`: the four-rate trajectory-preserving time
  speed used only by the anisotropic geometry.
- `ipm.diagnostics.measure.m`: invariants, boundary influence, and growth norms.
- `ipm.diagnostics.blowupFit.m`: optional physical-time growth fit; never a proof.
- `ipm.output.compare.m`: peak-aligned physical grid/domain comparison with
  laboratory peak shift reported separately.
- `ipm.output.visualize.m`: real-time fields, resolution/gauge traces, and MP4
  frames; no numerical-model logic.
- `ipmtests.baseline.transport.m`: independent smooth-order, conservation, jump, and
  SSPRK3 advection regressions for both transport choices.
- `analysis/ipm_fit_polar_profile.m`: read-only wall/polar profile extraction,
  radial-power fits, angular-fan scaling, and a two-dimensional single-scale
  consistency decision.
- `analysis/ipm_audit_active_ccf_wall_instantaneous.m`: read-only reconstruction
  of the physical wall observables `Q=max rho_x`, `M=-u1_x`, and the gauge-free
  Riccati coefficient `a=M/Q` for the historical opt-in CCF branch. Its authoritative
  current-WENO initial record takes zero time steps and checks the
  semi-discrete wall identity; historical records are explicitly tagged by
  their pre-WENO contracts and resolution gates.
- `analysis/ipm_analyze_active_ccf_material_fold.m`: read-only inverse-wall-map
  fit in density label `lambda=rho(x,0)`, reporting
  `A=min(x_lambda)=1/Q`, the moving label `lambda_c`, local curvature `B`,
  `deltaLambda=sqrt(A/B)`, and edge ratio `chi=lambda_c/deltaLambda`. It marks
  under-resolved saved fields untrusted rather than extrapolating through the
  grid failure.
- `analysis/ipm_audit_active_ccf_initial_width_rates.m`: current-WENO,
  zero-step, one-retained-physical-RHS tangent audit on the frozen CCF
  initial grid. After the initialization/flow solves, it differentiates
  horizontal and moving-peak vertical `rho_x` widths using virtual
  `rho +/- dt*rhs` fields only, without recomputing a virtual-field RHS. It
  also reports a stable subgrid
  `lambda_c'` and explicitly rejects curvature-derived tangents that fail its
  fit-level spread gate. The reported gamma values are instantaneous
  derivatives with respect to `A=1/Q`, never asymptotic similarity exponents.
- `analysis/ipm_plot_active_ccf_boundary_velocity.m`: zero-step plot and MAT
  audit of the historical opt-in CCF initial dynamic-transport field
  `transportU=(timeSpeed*u1+c_l*x+c_r,timeSpeed*u2+c_l*y)`. It decomposes the
  wall velocity into physical and `c_l*x` parts, locates the negative trough
  and right zero, compares physical and dynamic-frame strain, and retains
  physical elliptic metrics as a separate reference.
- `analysis/ACTIVE_CCF_DEFAULT_BLOWUP_AUDIT_ZH.md` and
  `analysis/ACTIVE_CCF_HEAVY_TAIL_THEORY_AUDIT_ZH.md`: the numerical-evidence
  boundary and the analytic DtN/material-label reduction for the former
  default CCF datum. They are retained as historical branch reports, not as
  documentation of the current no-argument `k=8` run.
- `analysis/FIXED_LABEL_GREEN_CORE_TAIL_CLOSURE_ZH.md`: independent derivation
  of the causal fixed-label Green current/production rows, core-tail cutoff
  with transition-annulus term, and the first-exit/tube inequality still
  needed for a conditional Riccati closure. Selector-dependent historical
  label diagnostics remain exploratory and are not proof inputs.
- `analysis/ipm_analyze_ccf_gaussian_trough.m`: read-only gradient, exact
  origin-strain, value-range, far-boundary, and candidate-fan diagnostics for
  the smooth Gaussian wall screen.
- `experiments/right_fan_soft_tail_case.m`: analytic Schwartz physical initial
  data made from a rotated anisotropic Gaussian trough and a broad weak return
  with exact signed-mass cancellation. Coarse, refined, wide, forced-threshold,
  and angle-scan variants are mechanism screens only.
- `analysis/ipm_analyze_right_fan_soft_tail.m`,
  `analysis/ipm_compare_right_fan_soft_tail.m`, and
  `analysis/ipm_compare_right_fan_angle_scan.m`: read-only right-fan growth,
  value-range, core-resolution, box/grid, and angular-selection audits. The
  standard 34/38.357/42 degree scan is monotone through 42 degrees and shows
  no two-scale separation; it does not select the homogeneous-cusp angle.
- `experiments/signed_slope_mixture_case.m`: smooth source-neutral physical
  initial data built from two continuous slope windows with exact signed
  coefficients satisfying `alpha=0`, so its infinite-half-plane wall velocity
  is pure compressive Burgers. Screen/refined and wide-domain variants remain
  mechanism tests, not invariant-manifold or blow-up claims.
- `analysis/ipm_analyze_signed_slope_mixture.m`: read-only exact-wall-closure,
  Burgers slope-growth, acute-fan, far-boundary, value-range, and source-budget
  diagnostics for the physical signed-slope screens. Since the signed mass is
  designed to be nearly zero, conservation is normalized by the initial
  weighted `L1` mass rather than by the signed mass. Every remeshed field is
  differentiated on its own saved `snapshotX/snapshotY` grid.
- `experiments/signed_slope_tuned_case.m`: initial-growth-selected analytic
  signed-slope data with windows `2.5/3.5` and half-width `0.15`. It includes
  raw screens, adaptive continuations, reproducible `385x193`/`513x257`
  analytic-premesh conservative candidates, a half-CFL/time-step check, the
  rejected non-bound-preserving WENO check, and a
  `[-640,640]x[0,128]` wide-box case. All variants preserve the exact initial
  wall law, source neutrality, and all-direction `rho_x -> 0`.
- `analysis/ipm_analyze_signed_slope_two_scale.m`: read-only per-snapshot
  two-scale audit. It centers on the positive wall `rho_x` peak, independently
  normalizes horizontal and vertical half-height widths, compares normalized
  fields on a common `(xi,eta)` grid, and fits reciprocal peak time, separate
  width exponents, core resolution, and profile drift. An optional common
  horizontal/vertical core-point floor excludes marginal terminal snapshots.
- `analysis/ipm_analyze_wall_riccati_feedback.m`: read-only physical
  half-plane audit of the exact wall law
  `q_t+u_1 q_x=-u_{1x}q`, `q=rho_x`. It rebuilds every selected Poisson flow
  on its saved grid, measures the nonlinear coefficient
  `a=-u_{1x}/q_peak`, checks the actual semi-discrete identity, and compares
  `t+1/(a q_peak)` with the reciprocal Type-I time.
- `analysis/ipm_analyze_wall_strain_budget.m`: Poisson-linear decomposition
  of the wall compression into sign/side, horizontal, vertical, and original-
  fan source blocks on selected saved grids. Homogeneous remote boundary data
  make every partition exactly additive; the actual Green-boundary correction
  is reported separately.
- `analysis/ipm_analyze_boundary_layer_reduction.m`: direct same-Poisson test
  of the Kiselev--Sarsam rectangular source
  `q_model(x,y)=q(x,0) 1_{0<y<h}`. It fits `h/L_x` against the complete local
  wall-strain, peak-frame velocity, and wall-evolution profiles. The present
  wide-box result rejects every single rectangular `H_h` layer despite an
  accidental peak-scalar Hilbert match.
- `analysis/ipm_analyze_multilayer_boundary_reduction.m`: alternating-node
  train/validation audit for twelve disjoint vertical slabs. It rejects a
  separable multi-`H_h` closure: the local validation errors remain
  `0.614/0.382/0.205` for strain, peak-frame velocity, and wall evolution,
  with large alternating coefficients and worse outer validation.
- `analysis/ipm_analyze_sheared_two_scale_rank.m`: low-rank audit of left and
  right `q=rho_x` traces in `zeta=(x-x_peak)/L_x-s*y/L_y`. It scans the shear,
  reports rank-2/rank-3 shape and Green-feedback errors, mode contributions,
  time drift, and cross-run subspace angles without claiming an invariant
  finite-dimensional manifold.
- `analysis/ipm_analyze_reflection_dominance.m`: exact positive-kernel split of
  the wall strain into positive and negative parts of the reflected difference
  `q(x_peak+z,y)-q(x_peak-z,y)`. It audits normalized windows and the complete
  saved physical box; finite-window positivity is not substituted for global
  reflection-cone invariance.
- `analysis/ipm_analyze_wall_interval_collapse.m`: read-only material-wall
  audit. It follows the `t=0` density labels of a connected positive wall
  branch without reconnecting across an untrusted snapshot, checks the exact
  label jump against `integral_I q`, verifies both endpoint-velocity and
  integrated-strain versions of `W'`, reports the label-average Riccati
  coefficient, moving-peak label drift, and conditional interval deadlines.
  Its nested fixed-label rows are localization probes only; they are not a
  claim that a time-dependent shrinking interval is material.
- `analysis/ipm_analyze_anisotropic_blowup_continuation.m`: read-only audit of
  a saved anisotropic continuation. It replays each stored grid/state with an
  independent aspect-aware Poisson solve, checks
  `(c_x-c_omega)-a_h*Omega_h`, reconstructs the physical clock and both scale
  factors, splits the full-box anisotropic Green feedback, and measures
  moving-peak normalized profile drift. Optional sheared rank analysis is
  diagnostic and does not add another Poisson/evolution path.
- `analysis/ipm_analyze_anisotropic_peak_rate_chain.m`: replay the stored
  canonical rates of a legacy anisotropic result, reconstruct its old
  semi-discrete RHS, and compare the exact quadratic peak/vertical chain
  against history without contaminating the audit with current gauge rates.
- `analysis/verify_ipm_anisotropic_scaling_contract.m`: analytic coordinate,
  Green-metric, conservative-form, and three-point quadratic peak Gateaux
  identities, including exact `c_omega` cancellation from vertex motion and
  a deliberately nonzero transport-minus-actual correction.
- `analysis/ANISOTROPIC_BLOWUP_CONTINUATION_AUDIT_ZH.md` and
  `analysis/anisotropic_blowup_continuation_report_template.tex`: the exact
  mapping, acceptance gates, and result-table contract for that continuation.
  Placeholders are not numerical evidence and must remain visibly unfilled
  until the long run passes the independent audit.
- `analysis/SHEARED_TWO_TRACE_RICCATI_REDUCTION_ZH.md`: derives the sheared
  two-trace reduction target, exact reflection and interval Green functionals,
  the aspect-uniform negative-tail lemma, the horizontal-zero-mode no-go for a
  universal pointwise reflection cone, and the fixed-density-label Riccati
  contract. It distinguishes a moving peak label from its limiting material
  label and records the remaining Casimir/core-tail/reservoir closure gaps.
- `analysis/WALL_RICCATI_STRAIN_CRITERION_ZH.md`: derives the exact half-plane
  left-right Green-kernel formula, a conditional finite-time wall blow-up
  proposition, its anisotropic normalized shape functional and formal
  `kappa -> 0` wall-Hilbert limit, the exact relation to CCF/1-D IPM blow-up
  theorems, the rejected single-layer reduction, and the converged three-run
  strain budget.
- `analysis/SMOOTH_SIGNED_SLOPE_BLOWUP_CANDIDATE_ZH.md`: authoritative formula,
  convergence ledger, singular-fit table, two-scale evidence, rejected WENO
  check, wide-box and strict-core comparisons, wall Riccati mechanism, and
  remaining numerical-proof gates for the strongest current smooth blow-up
  candidate.
- `analysis/ipm_analyze_slope_mixture_lift.m`: independent stationary
  signed-slope lift audit. It constructs smooth slope bumps with exact
  `alpha=0`, the beta-one-third Burgers bridge, and a first exact `A`-null
  normal-jet corrector, then reports DtN, wall-jet, bulk-residual, angular,
  neutrality, and all-direction `R_X` decay defects. Its current result is a
  wall-exact but bulk-noncandidate profile.
- `analysis/ipm_search_signed_slope_bulk_gluing.m`: constrained continuous-
  slope bulk-residual screen using analytic infinite-half-plane velocities.
  It enforces the trace, pure-Burgers DtN, and first normal-jet moments exactly
  and audits grid, signed-variation, and normalization dependence. A strict
  third-wall-jet identity proves this pure-mixture class cannot contain the
  desired smooth exact profile; low-residual/high-variation outputs are seeds
  for `A`-null correction, never candidates by themselves.
- `analysis/ipm_audit_signed_slope_bulk_continuation.m`: compares odd,
  one-sided, and smoothed one-sided pure-mixture continuations and reruns
  bump-count/variation-penalty checks with honest absolute, reaction, and
  `R_X` residual norms. The smooth odd audit gives `0.6537`, `0.6435`, and
  `0.5521`, with outside-fan `R_X^2` fraction `0.00294`; the residual remains
  order one despite the much smaller term-relative normalization.
- `analysis/ipm_search_anull_bulk_corrector.m`: finite-periodic-box
  Gauss--Newton screen of exact-DtN `Q2`--`Q4` normal correctors with
  mean-zero horizontal Gaussian centers. Its best unregularized audit is
  absolute/reaction/`R_X` residual `0.4783/0.4603/0.4051`, coefficient norm
  `349`, corrector ratio `||Q||/||R||=0.1685`, and outside-fan fraction
  `0.1118`; a `1e-2` penalty still gives absolute RMS `0.4863` and outside
  fraction `0.0775`. This is a rejected bulk screen, not evidence of global
  `R_X -> 0`.
- `analysis/ipm_search_anull_fan_corrector.m`: analysis-only
  characteristic-oriented `A`-null fan screen. It uses
  `zeta=kappa+mu*abs(k)+i*k*sigma` so each horizontal bump follows the
  straight coordinate `X-sigma*Y`; this is not a dynamically curved-
  characteristic basis. It preserves wall trace, first normal jet, and DtN,
  adds an explicit outside-fan `R_X` penalty, and reports absolute,
  reaction-normalized, `R_X`-normalized, term-relative, maximum, coefficient,
  and grid-refinement diagnostics. On the finest grids, broad gives
  `0.359438/0.249961/0.446312` with outside fraction `0.008242`, narrow gives
  `0.299396/0.208582/0.374605` with `0.005416`, and near-tip `n<=5` gives
  `0.314621/0.219182/0.390475` with `0.005631`. Their maximum residuals remain
  `5.57/5.38/5.16`, and narrow/tip coefficient norms are about `428`; all
  three are explicit noncandidates.
- `analysis/SIGNED_SLOPE_BULK_GLUING_AUDIT_ZH.md`: derivation and numerical
  evidence ledger for the signed-slope pure-mixture third-jet obstruction and
  the `Q2`--`Q4` `A`-null continuation. It separates the analytic all-direction
  decay of the base mixture/angular tail from the still-open full-space
  cap/decay problem for the corrector.
- `analysis/ANULL_FAN_CORRECTOR_AUDIT_ZH.md`: authoritative mathematical and
  numerical ledger for the straight characteristic-oriented fan screen,
  including machine-precision wall/DtN invariants, the rejected core-only
  pilot, broad/narrow/tip honest metrics, saved MAT/PNG artifacts, and the
  explicit stationary-noncandidate decision.
- `analysis/WANG_PROFILE_SEED_AUDIT_ZH.md`: official-asset and sign audit for
  the Wang IPM profiles. It records the later stable
  `lambda=1.0285722760323 +/- 1e-13`, preserves the first paper's distinct
  `1.0285722760222`, fixes `curl(U)=-H_X`, and documents that no official
  code, checkpoint, profile array, parameter tree, or SI is currently
  published. Only the per-field `Omega` ansatz is public; `H/U` envelopes and
  the exact normalization pin are not.
- `analysis/WANG_SYMMETRY_BREAKING_CONTINUATION_ZH.md`: analysis-only contract
  for reprojection, the full Poisson-coupled even/odd Jacobian, translation-
  bordered symmetry breaking, gauge phases, exact-PDE pseudo-arclength, and
  the anisotropic shape-space orbit. It reports no computed odd kernel,
  continuation branch, or two-scale profile.
- `analysis/ipm_train_compactified_ipm_profile.m`: analysis-only frozen-base
  additive multistage surrogate. Each zero-output correction is trained while
  all earlier networks remain frozen; derivatives and nonlinear residuals are
  evaluated on the composite field. It reports structured tail shells and
  honest factored/unfactored residuals and is not in the production call graph.
- `analysis/WANG_COMPACTIFIED_MULTISTAGE_CORRECTION_AUDIT_ZH.md`: schema,
  sign, frozen-parameter, shell, finite-box Green, and medium-training ledger
  for the compactified surrogate. The first correction is explicitly a
  raw-residual negative screen, not a Wang seed or an IPM profile.
- `analysis/verify_ipm_profile_fit.m`: synthetic acute-fan power-law accuracy
  regression with an arbitrary test exponent; it is outside the simulation
  call graph and does not select a physical scaling.
- `analysis/ipm_mesh_exact_fan_profile.m`: exact local polar mesh of the
  degree-one acute-fan steady state, with separate upper/lower wall traces
  and analytic physical/rescaled residual audits.
- `analysis/ipm_plot_exact_bessel_family.m`: exact physical steady Bessel
  modes on a punctured polar sector or the full upper half-plane, including
  analytic `rho_x`, signed-log structure plots, PDE residuals, and an explicit
  edge-normal-flow audit.
- `experiments/bessel_two_scale_case.m`: reproducible physical-mode K0
  external-pole and regularized K1-core screen/refinement cases; it only
  supplies options and initial-data handles to `ipm.solve`.
- `experiments/bessel_two_scale_initial_data.m` and
  `experiments/bessel_k1_two_scale_initial_data.m`: smooth near-Bessel
  initial densities, with the optional K1 wall-monitor component explicit.
- `analysis/ipm_analyze_bessel_two_scale.m` and
  `analysis/ipm_analyze_bessel_k1_two_scale.m`: read-only tail, width,
  heat-kernel, and `rho_x` zero-line diagnostics.
- `analysis/ipm_compare_bessel_two_scale_pair.m`: matched perturbed/reference
  difference ratios; a ratio near or below one is not called attraction
  without longer-time and refinement evidence.
- `analysis/verify_bessel_two_scale_initial_data.m`: pure K1-core wall and
  far-exterior regression, run with the monitor disabled.
- `analysis/reports/IPM_BESSEL_TWO_SCALE_REPORT.tex`: complete derivation,
  literature boundary, numerical screen, and next-step two-scale proposal.
- `analysis/RHOX_DECAY_BLOWUP_STATUS_ZH.md`: authoritative audit under the
  all-direction `rho_x in C_0` contract. It separates rigorous corner-domain,
  open harmonic-far-field, and standard decaying Biot--Savart conclusions;
  records the exact general-`beta` amplitude/length scaling, affine-cusp
  hodograph, CCF-to-IPM DtN gap, smooth-angular/soft-edge far tail, signed-head
  moments, dynamic cap, finite-mode no-go, and the beta-zero logarithmic
  benchmark.
- `analysis/ipm_plot_compact_acute_beta_tail.m`: reproducible smooth-angular
  general-`beta` tail supported in an arbitrary right acute sector. It solves
  the nonresonant Dirichlet angular Poisson problem spectrally and audits the
  all-direction `R_X~r^(beta-1)` decay, sector support, and projection error;
  it is leading asymptotic data, not a nonlinear fixed profile.
- `analysis/ipm_audit_farfield_resonance.m`: derives and evaluates the first
  nonlinear far-field correction, its near-Dirichlet resonance, and the
  beta-to-one-half logarithmic Fredholm limit.
- `analysis/ipm_audit_resonance_balanced_tail.m`: deterministic 22-mode
  right-fan/weak-return construction. It minimizes outside-fan `Q0^2` while
  imposing the leading cap-drift and first near-resonance moments; its
  independent-grid outside fraction is `0.00397014115`.
- `analysis/ipm_audit_balanced_tail_recursion.m`: independent angular-grid
  recurrence through `F6/G6`, including full stationary-transport residuals
  at fixed radii. The sixth-order residual is `1.0114e-2` at `r=30` and
  `2.8623e-4` at `r=100`; this is an outer asymptotic contract, not a global
  fixed profile.
- `analysis/ipm_scan_balanced_tail_beta.m`: reoptimizes the same 22-mode
  right-fan/weak-return seed while continuing `beta`. It is a far-tail
  conditioning/Pareto screen only: the stronger order separation at smaller
  `beta` is not an exponent-selection result.
- `analysis/ipm_audit_beta_third_log_tail.m`: isolates the first exact
  rational Dirichlet resonance. The generic Burgers value `beta=1/3`
  resonates at `(order,mode)=(5,2)` and the codimension-four value `beta=1/5`
  at `(4,2)`; both force a logarithmic stream correction.
- `analysis/ipm_audit_log_transseries_tail.m`: closes the full far-field
  recurrence in polynomials of `log(r)`, audits Poisson and transport
  coefficient defects, and evaluates optimally truncated residuals. Its
  dual-grid `beta=1/3` contract is trustworthy only through about `N=7--8`:
  the fine-grid residual is `0.07123` at `r=3,N=7` and `0.01737` at
  `r=4,N=8`. Coarse `N=20` improvements are explicitly rejected after angular
  refinement.
- `analysis/ipm_optimize_beta_third_unitary_overlap.m` and
  `analysis/ipm_audit_beta_third_unitary_overlap_candidate.m`: constrained
  response-space optimization and independent audit of the actual `N=8`
  log-transseries overlap. The accepted strong branch preserves cap, `M1`,
  normalization, outside-fan energy `0.0061`, and all-direction
  `rho_x=O(r^-2/3)`. On `8001x1800`, its residual at
  `r=2,3,4,5,7,10` is
  `0.014741/0.0016888/0.00036331/0.00011036/1.8313e-5/2.7290e-6`,
  about `47.5x` below the old branch for `r>=3`; the high-order coefficient
  RMS and spectral tail also decrease. It remains a formal `N<=8` far-tail
  candidate and changes the old gluing traces by `35%--37%`.
- `analysis/BETA_THIRD_UNITARY_OVERLAP_OPTIMIZATION_AUDIT_ZH.md`: authoritative
  strong/Pareto branch ledger, dual-grid residual and spectrum tables,
  resonance logs, near-skew-generator interpretation, `r=2` term-ratio audit,
  and the warning that every inner/DtN trace must be regenerated.
- `analysis/ipm_audit_beta_third_n11_continuation.m`,
  `analysis/ipm_audit_beta_third_optimal_truncation.m`, and
  `analysis/BETA_THIRD_N11_CONTINUATION_AND_TRACE_API_ZH.md`: continuation of
  the fixed strong branch through the `(n,m)=(11,6)` resonance, spectral
  cutoff audit, optimal-truncation/Gevrey screen, and the authoritative N8/N11
  ledger. The resonant log is at roundoff, while nonresonant coefficient
  growth makes N8 optimal at `r=2`; N11 remains improving for `r>=3`.
- `analysis/ipm_export_beta_third_far_tail_traces.m` and
  `analysis/ipm_evaluate_beta_third_far_tail.m`: standard regenerated gluing
  API for modal, fixed-X/fixed-Y, radial, and characteristic trace data. The
  signed density requires sign charts or `log(abs(R))`; old boundary data are
  explicitly nonreusable.
- `analysis/RESONANCE_BALANCED_FARFIELD_AUDIT_ZH.md`: authoritative
  derivation, sign caveat, logarithmic limit, deterministic seed, dual-grid
  recursion ledger, high-order small-divisor warning, and gluing boundary.
- `analysis/LOG_TRANSSERIES_FARFIELD_AUDIT_ZH.md`: authoritative exact-
  resonance recurrence, implementation audit, nonresonant regression,
  dual-grid optimal-truncation ledger, and the explicit retraction of the
  under-resolved high-order result.
- `analysis/ipm_horizontal_overlap_audit.m`: analysis-only near-horizontal
  response-chart test. The current one-mode audit preserves overlap states
  and invariants accurately but has order-one PDE residual and merely moves
  the defect across the `3--10` degree sector; it is a registered
  noncandidate.
- `analysis/SOFT_TAIL_HODOGRAPH_AUDIT_ZH.md`: derivation and independent-grid
  ledger for finite-`Y` far-tail traces, the first `F1/G1` correction, vertical
  hodograph continuation, and horizontal-chart overlap. It records that
  `beta=1/3` and `0.21` make the present vertical graph fold (`K<=0`) near the
  wall, so a fold-aware multi-chart representation is required.
- `analysis/CASIMIR_CAP_GLUING_AUDIT_ZH.md`: complete level-area/Casimir audit
  of the moving outer reservoir. It proves the rescaled cap amplitude/radius/
  area laws `L^-beta/L^-1/L^-2`, rejects stationary compact caps and radial
  cutoffs, and gives a conditional level-wise Hamiltonian correction through
  six orders. The unresolved contract is that the same generators must be
  produced by the density's IPM Poisson velocity.
- `analysis/ipm_audit_beta_third_s5_cancellation.m` and
  `analysis/BETA_THIRD_S5_CANCELLATION_AUDIT_ZH.md`: exact kinetic-energy
  identity for the first `(n,m)=(5,2)` resonance. It proves that `S5` is a
  sixth-degree functional of `F0`, not a consequence of cap plus `M1`, and
  supplies a hard-constraint counterexample. The balanced seed's near
  `[1,-4,6,-4,1]` blocks reflect a near-skew orbit, not a universal identity.
- `analysis/ipm_optimize_beta_third_resonances.m` and
  `analysis/BETA_THIRD_DOUBLE_RESONANCE_PARETO_AUDIT_ZH.md`: exact `S5/S8`
  feasibility/Pareto audit. The two constraints have no local rank no-go, but
  exact cancellation moves energy into nonresonant high modes and worsens the
  true `N=8` overlap; it is retained as a negative design lesson.
- `analysis/ipm_audit_folded_lagrangian_charts.m` and
  `analysis/FOLDED_LAGRANGIAN_MULTICHART_AUDIT_ZH.md`: canonical two-chart
  derivation and regression. A sign change of one projection derivative is a
  chart fold, not a physical singularity; only the full Lagrangian Jacobian is
  invariant. The affine regression keeps full `J>0` while one chart folds and
  satisfies canonical/Poisson/Stokes checks.
- `analysis/ipm_screen_folded_lagrangian_gluing.m` and
  `analysis/FOLDED_LAGRANGIAN_GLUING_SCREEN_RESULTS_ZH.md`: 28-parameter
  `G/H` generating-chart screen with one shared Eulerian Poisson solve,
  multiplicity, Casimir, full-`J`, seam, and all-direction-decay gates. All
  gates pass, but validation residual remains `2.52591`, only `0.29%` below
  the seed. At fixed labels the `inner_R` block is orthogonal to every pure
  geometry mode; signed-Burgers core density and Casimir-null Eulerian collar
  modes must be opened next.
- `analysis/ipm_audit_burgers_casimir_schur_reachability.m` and
  `analysis/BURGERS_CASIMIR_SCHUR_REACHABILITY_AUDIT_ZH.md`: exact quantile
  matching of the signed-Burgers core to the N7/N8 tail level distribution,
  followed by a Casimir-null Hamiltonian-collar Schur audit. Exact Casimir
  matching exists, so the old fixed-affine `inner_R` obstruction is not a
  value-distribution obstruction. However, roughly `44.1%` rank inversions and
  `60--64` tail turns give an endpoint-order no-go for every
  orientation-preserving multiplicity-one boundary-to-boundary relabeling;
  the admissible collar improves the independent residual by under `1%`.
- `analysis/ipm_plot_compact_acute_tail.m`: reproducible construction and
  plot of a two-bump, same-constant-return beta-zero leading tail whose
  `R_X` is active only in a prescribed right acute angular sector. It audits
  the leading Poisson equation, Fredholm and return moments, angular support,
  and all-direction `R_X=O(r^-1)` decay; it is asymptotic data, not a nonlinear
  fixed profile.
- `analysis/ipm_search_beta0_profile.m`: finite-dimensional beta-zero fixed-
  profile search with exact discrete Poisson elimination, affine wall-jet
  pinning, explicit harmonic outer strain, inflow-only density data, outer
  `R_X` matching, an analytic LM Jacobian, and box-tail diagnostics. Its
  current wide-box output is a documented noncandidate, not a blow-up claim.
- `analysis/ipm_search_beta_hodograph.m`: general-`beta` single-chart
  affine-cusp-to-smooth-acute-tail prototype. It eliminates the kinematic and
  transport identities through `nu=1/mu`, optimizes only the remaining
  Poisson compatibility, reconstructs `P`, and audits geometry and
  all-direction decay. The exact characteristic invariant proves that a
  regular positive chart cannot include the flat zero of a compact angular
  tail without `J=X_q` degenerating; a global candidate needs separate edge
  charts. Returned results therefore always mark `isGlobalCandidate=false`.
- `analysis/BETA_HODOGRAPH_EDGE_PROPOSITION_ZH.md`: proof of the exact
  characteristic invariant and the induced affine-to-tail matching law. A
  single regular positive chart cannot attach two `C-infinity` flat angular
  edges to nonzero affine characteristics.
- `analysis/ipm_audit_beta_edge_charts.m`: reproducible slope/angle edge-chart,
  overlap, and beta-continuation audit. The current beta-about-0.70 interior
  residual valley is explicitly nonconverged, and the direct two-edge-chart
  splice fails the overlap and Poisson tests.
- `analysis/BESSEL_EXACT_PROFILE_AUDIT_ZH.md`: independent exactness,
  decay, sign, distributional-source, and local multipole audit of the
  Bessel family.
- `analysis/BESSEL_MULTISCALE_BLOWUP_THEORY_ZH.md`: corrected no-go theorem
  for a Bessel-leading smooth-core layer, the exact closure-defect Green-kernel
  formulation, the exact affine hyperbolic blow-up/Taylor jet, the leading
  general-`beta` supported terminal-cusp angle, and conditional multiscale
  gluing laws.
- `analysis/BESSEL_BRIDGE_EXACT_BLOWUP_ZH.md`: exact dynamically rescaled
  Bessel-remainder equation; bounded-density full-plane gradient blow-up under
  a prescribed nondecaying harmonic far-field contract; exact directional and
  three-scale wall modes; and an explicit nonlinear characteristic proof that
  the right-decaying directional profile family is orbitally attracting on
  compact subsets of the right characteristic sector in the natural monotone
  `Y`-linear invariant class. It also gives an all-`X`, pointwise-smooth exact
  trajectory from a `C_b^infinity` reduced initial profile to that Bessel
  family; an even stronger trajectory whose physical density is bounded,
  `C^infinity`, and compactly supported in `x`, and whose rescaled core
  converges locally in `C^infinity` to `exp(-X)`; a time-independent
  constant-strain contract whose unique translation-invariant classical orbit
  has a pure long-time Bessel tangent; wall-compatible sine regularizations
  and bounded rotating-strip solutions with compact horizontal and vertical
  slices; together with the
  potential-energy, two-Casimir, source-moment, and global-boundedness
  obstructions, plus separated-transverse and finite-Fourier no-go results,
  that prevent these open-boundary harmonic-strain examples from becoming
  decaying finite-energy source-neutral half-plane Cauchy attractors.
- `analysis/POLAR_PROFILE_RESEARCH.md`: centered polar steady equations,
  symmetric/asymmetric boundary contracts, fixed-angle versus shrinking-fan
  alternatives, and the exact far-decaying Bessel family.

Keep these functions stateless. Add a file only for a distinct responsibility.

## Numerical and verification contracts

- First derivatives, Poisson, and the composed signed IPM velocity operator
  are second order on uniform and smoothly stretched grids.
- Multiplying the pointwise Poisson operator by the matching Jacobian/control-
  volume weights must produce a symmetric stiffness matrix to roundoff.
- Time stepping is SSPRK3 with `cfl<=0.5`.
- `transportScheme='weno5_nonuniform'` is active in `ipm.config.activeCase`.
  It keeps the solver-wide nodal-value meaning of `rho`, uses three quadratic
  candidates plus the global quartic interpolant on every complete five-node
  stencil, and rebuilds all geometry after an adaptive remesh. The two outer
  face layers use MUSCL-minmod. `transportScheme='muscl_minmod'` remains the
  generic robust fallback.
- On smooth stretched grids the WENO face reconstruction is fifth order. The
  complete conservative nodal flux remains second order because it uses the
  existing dual-volume divergence; this matches the second-order derivative,
  Poisson, and velocity contracts and must not be advertised as a fifth-order
  IPM solver. Treating nodal values as finite-volume cell averages is forbidden:
  a trial did so and created a false wall-`rho_x` secondary peak.
  The reconstruction follows [Martí, Mulet, Yáñez, and Zorío, *Efficient WENO
  Schemes for Nonuniform Grids*](https://link.springer.com/article/10.1007/s10915-024-02558-6),
  using its point-value algorithm and precomputed fixed-grid geometry.
- Physical transport is conservative. Weighted mass is conserved to roundoff
  only when the three artificial sides use `transportBoundaryMode='closed'`;
  the production open box changes mass by its resolved boundary flux.
- The discrete velocity divergence and wall-normal velocity vanish up to
  roundoff.
- Constant density gives zero source, velocity, and physical RHS.
- Constant density under pure dynamic dilation is preserved by the open
  conservative
  boundary flux plus `2*c_l R` to roundoff.
- The bounded rational upper trace, its `x_1` evenness, nonzero wall trace,
  analytic odd `omega=d_x rho`, and numerical derivative convergence are
  checked.
- The active logarithmic primitive is checked against its exact density and
  `omega=sign(x)|x|^7/(1+|x|^8+y^8)`, including one-sided sign, analytic wall
  peak, the exact transport-anchor zero, and a negative initial combined-
  velocity slope at `X=1`.
- Double-odd mode checks `rho` even, `omega` and `psi` odd, the corresponding
  odd/even velocity parities, even RHS, paired remesh symmetry, zero origin
  velocity, `c_r=0`, and `X_shift=0` through a time step.
- A constructed main peak plus separated secondary peak checks that adaptive
  resolution counts only the main peak's connected component and uses the
  true nonuniform cell widths.
- An exact quadratic on a nonuniform three-point stencil checks subgrid peak
  localization.
- Adaptive-remesh tests check strict grid monotonicity, an exact origin node,
  range and row-mass preservation, reduced peak spacing, exact quadratic
  differentiation on the resulting arbitrary grid, and rejection of repeated
  no-improvement redistribution.
- A synthetic twofold expansion checks that the additive width correction is
  negative. A conflicting-scale case, with an over-expanded outer level and
  under-resolved inner levels, checks that outer contraction has priority.
- A synthetic large-rate case checks that one common time factor enforces the
  dynamic-rate cap while retaining the canonical rate ratios.
- A one-step deliberately capped dynamic run checks that canonical and physical
  clocks advance positively and more slowly than normalized time.
- Synthetic `(1-t)^(-2)` and exponential histories check that the optional
  blow-up fitter identifies the former and rejects the latter.
- Identical peak-aligned fields with different laboratory peak positions check
  that `ipm.output.compare` reports zero shape error and the translation
  separately.
- Dynamic mode checks finite rates, the configured additive width correction,
  zero local-strain length-gauge residual, the active zero-origin symmetry residual,
  a positive wall-omega peak denominator, and zero semi-discrete peak rate.
  The full-gradient-energy gauge is tested separately for positive energy and
  zero semi-discrete energy rate. The legacy point-strain mode is also tested
  separately but is not default.
- The optional wall-density-width length gauge is checked for a finite rate
  and a zero two-crossing width residual, but is not used by production runs.
- In explicit half-plane comparison mode, adding a constant horizontal velocity
  must leave canonical `c_l,c_omega` unchanged and shift only `c_r` oppositely.
- In that comparison mode, multiplying a translated rescaled state by a positive constant must multiply
  all three canonical rates linearly and its canonical RHS quadratically.
- Peak, gradient-energy, and `c_omega=0` right-hand sides must exactly match
  the selected full-velocity conservative assembly for their canonical rates.
- Reusing automatically generated options after changing the domain or point
  count must recompute `gridStretch`; explicit numeric stretch remains fixed.
- Both transport choices use each local face distance and the matching
  nonuniform control-volume widths. Closed-boundary weighted mass, a constant
  state, a sharp jump, and short SSPRK3 translation are checked separately.
- Green data are exactly zero on the physical wall. Both the four-image
  quadrant kernel and legacy single-image half-plane kernel are checked
  against a unit weighted source on all three artificial sides.
- Smooth rational data check compressed Green boundary values and the induced
  core velocity against an uncompressed source quadrature.
- Localized data check `farBoundaryVelocityRatio` against
  `farBoundaryTolerance`. The uncut rational datum instead checks its source
  tail and requires peak-aligned box expansion because its angular far velocity
  remains `O(1)`. Enlarging only the box without retaining center resolution is
  not an acceptable convergence test.

`ipmtests.baseline.suite` enforces these items with analytic manufactured fields for
`D_x`, `D_y`, the Poisson inverse, and the complete velocity map, plus
  algebraic, divergence, mass, constant-state, initial-data, dynamic-gauge,
  Green-boundary, and doubled-domain checks.

### Current verification record (2026-08-28)

- The selected nodal WENO reconstruction has final smooth refinement order
  `4.701`; the complete conservative nodal transport has order `2.000`.
  A nonuniform-grid step translation has L1 errors `3.438e-2` for WENO and
  `5.280e-2` for MUSCL, WENO TV ratio `1.0001`, and no material range
  overshoot. Closed WENO mass rate and constant-state residual are at
  roundoff.
- A matched `513 x 257` CCF/dynamic/lattice run through physical time `0.01`
  completed 51 CFL-controlled steps and three analytic remeshes with a single
  wall-`rho_x` maximum, TV ratio exactly `1` at recorded precision, density
  range violation `5.84e-12`, x/y adjacent ratios `1.492/1.015`, and zero
  translation. This is an integration regression, not a blow-up result.
- `ipmtests.baseline.suite` passes with the transport suite included, and MATLAB Code
  Analyzer reports zero messages for the beta-zero fixed-profile search and
  compact-acute-tail construction scripts. The compact-tail audit separately
  closes its Poisson/Fredholm/return constraints at roundoff.

- The active double-odd contract now passes direct parity tests: initial and
  evolved `rho` are `x_1` even; `omega,psi,u_1` have the required odd parity;
  `u_2` and the RHS have the corresponding even parity. Reported defects are
  between zero and `3e-15`; origin speed is `4.6e-17`, and both `c_r` and
  one-step `X_shift` are exactly zero. Paired adaptive remeshing preserves
  symmetry to reported precision.
- The four-image Green kernel agrees with direct quadrant evaluation to
  `1.27e-16`; the retained half-plane kernel agrees to `1.58e-16`. Compressed
  quadrant Green data/core velocity differ from uncompressed quadrature by
  `1.23e-3/7.46e-6` in the regression.
- A short end-to-end active-mode run to physical time `0.1` retained density
  evenness `8.9e-16`, omega oddness `4.9e-15`, origin speed `2.0e-15`, zero
  translation, and zero density-range violation.
- Native MATLAB `ipmtests.baseline.suite` passes. Final refinement orders are `1.999`
  for Poisson, `1.979/2.040` for `D_x/D_y`, `2.006` for the signed
  Biot--Savart map, and `2.058` for the even-trace/odd-derivative rational
  `rho -> rho_x` test. The maximum Poisson algebraic residual is
  `5.43e-12`; divergence is `2.44e-15`; wall-normal velocity is exactly
  zero at reported precision.

### Active physical convergence record (2026-08-22)

- Every accepted physical run has `c_l=c_omega=c_r=0`, `C_l=C_omega=1`,
  and `X_shift=0` exactly in its saved history. The no-`c` equation, not a
  reconstructed dynamic trajectory, supplies all numbers in this subsection.
- A clean `257 x 129`, box-100 run became under-resolved near `t=1.376`; its
  late-time peak and fitted exponent differ materially from finer runs, so it
  is rejected as blow-up evidence.
- At `t=1.4` on box 200, refining `513 x 257 -> 1025 x 513` changes the wall
  peak by `0.63%` and FWHM by `3.36%`. The peak position changes by `5.4%`,
  and horizontal/wall-normal 90% widths still change by about `19%/36%`.
  Thus the broader peak begins to converge while the innermost multiscale
  structure does not.
- On `1025 x 513`, expanding box 200 to 400 at fixed center spacing changes
  peak/FWHM/position by `1.64%/2.27%/0.86%` at `t=1.4`; through `t=1.38`
  all three differences are about `1.1%` or less. The wall-normal 90% width
  still changes by `7.5%` at `t=1.4`.
- On the `513 x 257`, box-200 comparison, halving `maxDt` from `5e-4` to
  `2.5e-4` changes peak/FWHM by `0.27%/0.57%` at `t=1.36`,
  `1.60%/3.46%` at `t=1.38`, and `3.88%/3.86%` at `t=1.4`.
  Consequently terminal high-accuracy runs use `2.5e-4`; a visually smooth
  curve at the larger step is not treated as time converged.
- The time-refined box-400 `1025 x 513` run reaches trusted
  `t=1.41001039` with peak `37.8769`, horizontal/wall-normal core counts
  `6.07/7.99`, range violation `7.7e-7`, and exact zero modulation. It then
  stops at a deliberately excluded under-resolved state. The 40/80-point
  full-gradient windows give `T about 1.466`, `p about 1.32`, but 20/160-point
  windows drift to approximately `(1.60,3.97)` and `(1.45,1.12)`; the
  exponent is not converged. Wall `rho_x1` and `rho_x2` also favor different
  terminal scales, consistent with possible multiscale growth.

The longer grid/box/time-two records below were produced before activating
the double-odd `x_1` image and are retained as explicit `half_plane`
comparison history. They do **not** establish convergence or blow-up behavior
for the new active symmetry class; matched double-odd long runs must be
recomputed.
- Closed-box MUSCL mass rate is `2.00e-17`. Open pure dilation preserves a
  constant state with infinity error `6.66e-15`. The peak/energy/none
  canonical full-flux assembly residuals are all below `2.6e-16`.
- Multiplying a translated state by three scales the three canonical rates by
  three and the canonical RHS by nine with errors `3.09e-14` and
  `5.08e-14`. Adding horizontal velocity seven leaves `c_l,c_omega`
  unchanged and subtracts seven from `c_r`.
- Green data reduce the localized doubled-box core-velocity change from
  `6.67e-2` (zero far data) to `2.51e-3`. Signed 3000-source aggregation
  has boundary/core-velocity errors `8.16e-4/9.07e-6` in the MATLAB
  regression. On boxes 2000 and 4000, its initial core-velocity error versus
  full source quadrature is about `1e-5`; increasing to 12000 sources changes
  that only slightly.
- An independent Python mirror gives the same stretched-grid orders and Green
  errors. With the active open full-velocity flux on a `257 x 129`, box-500
  run to physical time two, density remained approximately
  `[0,0.9828]`, rescaled peak variation was `3.7e-5` relative, connected
  expansion stayed below `1.33`, and canonical rates stayed below `0.75`.
  A final gauge-covariant `129 x 65`, box-2000 mirror to physical time one
  kept the rescaled peak fixed to `3.3e-15` relative, expansion below `1.134`,
  canonical rates below `0.708`, and density in approximately
  `[0,1.000091]`.
- Real MATLAB `257 x 129`, box-500 runs at physical time two give
  `rho in [0,1.0000187]`. Halving `maxDt` from `0.005` to `0.0025`
  changes peak-aligned physical density by `1.37e-6` relative L2. Tightening
  `maxDynamicRate` from two to `0.2` changes it by `4.76e-9`; only the
  normalized elapsed time changes.
- At the same physical time, `wall_omega_peak`, `gradient_energy`, and
  `none` gauges differ in peak-aligned physical density by about
  `1.47e-6` relative L2 and in physical `rho_x` by about `7.7e-6`.
  Rescaled peaks differ as intended. This is the primary long-run evidence
  that amplitude gauge no longer contaminates the physical solution.
- Grid study on box 500: `129 x 65` is rejected (about `2.1%` density and
  `13%` derivative-field error versus `257 x 129`). `257 x 129` versus
  `513 x 257` still differs by `0.85%` in density and `4.7%` in
  `rho_x`. Raising only the vertical count from `257` to `513` at
  `nx=513` changes density by `0.16%`, `rho_x` by `0.91%`, the wall
  peak by `0.09%`, and its width by `0.33%`.
- Far-box study at `513 x 257` shows that fixed point count cannot be enlarged
  indefinitely: source-tail ratios fall from about `1.4%` (500) to
  `0.87%` (1000), `0.47%` (2000), and `0.24%` (4000), but the box-4000
  outer face grows to about 132 and corrupts local extrema. Green compression
  is not the cause.
- The matched high-resolution comparison `1025 x 513`, boxes 2000 and 4000,
  is the active quantitative far-field record. At physical time two the
  peak-aligned density difference is `6.36e-4` relative L2 and the
  `rho_x`-field difference is `6.36e-3`. The single wall maximum and
  half-height width still differ by about `9%` and `16%`; therefore the
  overall derivative field is near percent-level convergence, but a pointwise
  extremum is not yet converged enough to claim blow-up.
- Across all completed long MATLAB runs, physical density stayed within
  approximately `[0,1+2e-5]`, connected peak expansion stayed in roughly
  `1.25--1.37`, and no coupled unbounded growth of `c_l,c_omega` appeared.
  These are stability observations, not a proof of regularity or blow-up.
- The new bounded-remesh stress series isolates the remaining grid limit. On
  `257 x 129`, box 500, center target `0.02`, a half-target spacing floor
  preserves monotonic nodes, exact `X=0`, and density range to about `3e-5`,
  but the 90% derivative core still reaches the two-cell stop near physical
  time `1.94`. Increasing only the emergency `c_l` response made canonical
  `c_l,c_omega` exceed five and failed earlier; allowing quarter-target local
  spacing did not finish interactively because of the global explicit CFL.
  Therefore bounded remeshing prevents artificial cumulative contraction and
  exposes failure cleanly, but does not by itself resolve the observed late
  multiscale core. Local time stepping or a multi-block AMR transport/Poisson
  method is required before extending quantitative blow-up claims past this
  point.
- A final default-resolution attempt to extend the common physical endpoint
  from two to three remained numerically alive but did not finish within the
  four-hour test window and had no checkpoint output. It was interrupted
  cleanly; physical time two remains the furthest saved quantitative endpoint.
  The rapidly rising cost of `t=2 -> 3` is an open long-run monitoring item,
  not evidence of blow-up by itself.
- MATLAB Code Analyzer is clean after the final numerical changes. The default
  `1025 x 513`, box-2000 path is intentionally conservative and expensive;
  `513 x 257` is suitable for exploratory videos but not for pointwise
  blow-up claims.
## Blow-up observation

Physical IPM transport preserves the range of density, so density-amplitude
growth is not by itself a physical blow-up signal. The primary signal is
`history.physicalGradInf`, evaluated against `history.physicalTime`. Also
inspect `physicalGradientEnergy`, `physicalGradientEffectiveLength`, the
inner/near/outer gradient-energy fractions, time-step collapse, resolution
dependence, gauge residuals, and distance from artificial boundaries. For the
active discontinuous odd `x_2` extension, these are upper one-sided regular
gradients and exclude the prescribed axis delta already present at `t=0`.
They cannot be described as creation of a full-plane singularity. For the
uncut degenerate data, neither global gradient energy nor its resulting
`c_omega` should be compared across boxes without an explicit box-refinement
study. A
divergent point-strain `c_omega` accompanied by a vanishing `strainCondition`
is a gauge singularity, not by itself a physical blow-up. A candidate power
law must outperform exponential growth and persist under grid and domain
refinement. The fit is a diagnostic only.
Convergence of the peak-aligned density, or even the full `rho_x` field, does
not imply convergence of its single largest grid value. The earlier matched
large-box agreement belongs to half-plane comparison mode and must be
regenerated for the four-image kernel. Point extrema must not be used alone
as a blow-up claim.
The reconstructed `physicalRhoMin/physicalRhoMax` and
`physicalRangeViolation` explicitly audit the maximum principle. Evolution
stops with `maximum_principle_violation` if the violation reaches the default
`rangeStopTolerance=5e-3`. All completed real MATLAB time-two studies stayed
within about `2e-5` of the initial density range.
Wall-profile oscillation uses the connected positive trace. An adjacent-node
split of a flat main peak is not a failure by itself: multiple discrete maxima
stop or invalidate a run only when the total-variation ratio also exceeds
`oscillationTVTolerance` (default `1.02`). A negative wall lobe above
`positiveWallNegativeTolerance` remains an independent failure.

### Wall density-label collapse research contract

For every classical half-plane solution, the bottom trace obeys
`rho_bar_t+v*rho_bar_x=0` and
`q_t+partial_x(v*q)=0`, where `q=rho_bar_x` and `v=u_1|wall`.
Inside any positive wall branch, use the transported density itself as the
material coordinate:

```text
ell = rho_bar(x,t),        x = X(t,ell),        Q(ell,t)=q(X(t,ell),t),
X_t = v(X,t),              X_ell = 1/Q,
a(ell,t) = -v_x(X,t)/Q = -X_tell,
Q_t = a Q^2,               (1/Q)_t = -a.
```

These differential coordinates are used only on the open branch `Q>0`.
The complete-component zero endpoints are recovered by an improper limit of
interior labels; `X_ell=1/Q` is not treated as finite at a zero endpoint.
The two zero endpoints of a complete positive component are only one useful
choice. For any two fixed labels `ell_- < ell_+` within the branch,

```text
Delta_ell = ell_+-ell_-,
W = X(t,ell_+)-X(t,ell_-) = integral_[ell_-,ell_+] 1/Q dell,
B_I = [v(X(t,ell_-))-v(X(t,ell_+))]/Delta_ell
    = average_[ell_-,ell_+] a,
W' = -Delta_ell*B_I.
```

Thus a uniform future lower bound `B_I >= B_0 > 0` is a conditional finite-
time collapse theorem for every fixed-label subinterval, not only for zero
endpoints. It is not an unconditional result inferred from one late-time
frame. The analyzer fixes every label from snapshot 1, uses that initial
`Delta_ell` in all invariant and deadline rows, verifies both
`Delta_ell=integral_I q` and the two independent forms of `W'`, and stops
rather than reconnecting a component across a failed trust gate. A reported
deadline must remain labeled conditional on persistence of its measured
lower bound.

The pointwise reflected-difference condition is deliberately weaker evidence,
not an invariant-cone contract. Adding a horizontal density zero mode `h(y)`
leaves `q`, the Poisson velocity, the reflected difference, the current
material interval, and its current `B_I` unchanged, but changes the exact
reflection forcing and generally `B_I'` through the missing vertical trace.
Therefore unrestricted full-density IPM admits no universal `q`-only or bare
pointwise-reflection cone. A valid forward-invariance statement must also
control the core production, the signed Green tail, and the relevant
equimeasurable/Casimir leaf.

A moving wall maximum is not a fixed label. On one `C1`, nondegenerate peak
branch with no peak switching,

```text
ell_*' = v_xx(x_*,t)*q_*^2/q_xx(x_*,0,t).
```

Under the conditional scales `q_* ~ (T-t)^(-1)`,
`L_x ~ (T-t)^(a_x)`, with the corresponding first- and second-derivative
bounds, this drift is `O((T-t)^(a_x-2))`. It is time-integrable only when
`a_x>1`, yielding a limiting fixed label `ell_infinity`; the intended local
proof route is then the exact fixed-label identity
`1/Q(t)=1/Q(t_0)-integral_[t_0,t] a ds` along
`X(t,ell_infinity)`. This remains conditional until the derivative bounds and
positive time integral are proved. Converting label drift into a physical
offset additionally requires `Q >= c*q_*` on the full label segment joining
`ell_*(t)` to `ell_infinity`; it then gives only
`|X(t,ell_infinity)-x_*(t)|/L_x=O(1)`, not convergence of that normalized
offset.

The same scaling gives a core label mass
`q_* L_x ~ (T-t)^(a_x-1) -> 0`. Consequently the limits of a nested label
width `f -> 0` and `t -> T` need not commute. A time-dependent
`f(t)` does not have fixed material endpoints and cannot inherit the interval
theorem. Only under the explicit extrapolation that the observed tail lower
bound persists for the entire future classical interval does the current
full-component row give a deadline near `1.69`; versus the candidate time near
`0.772`, it is retained only as a deliberately weak sufficient screen. No
nested-label number is promoted before its analyzer is complete.

### Anisotropic continuation audit contract

The anisotropic continuation audit consumes a saved `ipm.solve` result and is
strictly read-only. It uses the exact convention
`kappa=C_y/C_x=L_x/L_y`, `dt/dtau=C_omega/C_x`, and independently recomputes
the aspect-aware Poisson/velocity map on every audited saved grid. It must
check the physical and gauge versions of the wall identity

```text
(c_x-c_omega)-a_h*Omega_h+Gamma_h/Omega_h = 0,
a_h=-U_1X/Omega_h,
```

where `Gamma_h=Omega_tau+v_* Omega_X` is the moving-peak amplitude residual;
`Gamma_h=0` is used only when the selected amplitude gauge actually enforces
it.
It also reports the full-box positive/negative anisotropic Green split,
physical-clock reconstruction, aspect degeneration, independent gauge and
Poisson residuals, resolution/maximum-principle gates, and peak-aligned
profile drift. A rate plateau alone, or a profile made stationary by a gauge,
is not convergence to a new dynamic steady state. Promotion requires all
stored fields and contemporaneous grids, trusted late-time gates, positive
Riccati feedback, mutually consistent remaining-time estimates, and grid/box
replication.

### Stationary bulk-screen research contract

Stationary signed-slope and `A`-null searches are analysis-only and never enter
`ipm.solve` or the production solver call graph. A bulk screen must report, on
the same full audit box, absolute RMS, residual/reaction, residual/`R_X`,
term-relative RMS, absolute maximum, outside-fan `R_X^2` fraction, coefficient
norm, corrector/base ratio, wall/DtN/Poisson invariants, and at least one
unoptimized grid refinement. A small term-relative value is not a candidate
criterion because it can be caused by cancellation of large transport terms.
Core-only fits that fail on the box top or periodic boundary are rejected.

The current fan-localized `A`-null basis is straight
characteristic-oriented: `sigma` tilts a bump along `X-sigma*Y`, while `mu`
and `kappa` set normal decay. It must not be described as a curved-
characteristic closure. The authoritative broad/narrow/tip finest-grid
absolute/reaction/`R_X` triples are
`0.359438/0.249961/0.446312`,
`0.299396/0.208582/0.374605`, and
`0.314621/0.219182/0.390475`; their outside fractions are
`0.008242/0.005416/0.005631`. Because the best absolute and
reaction-normalized residuals remain about `0.30` and `0.21`, maximum
residuals remain above `5`, and tip enrichment regresses, all three are
noncandidates. Any successor must change the bulk geometry, use genuinely
curved characteristics or multiple edge charts, or start from a different
bulk seed; adding more straight Taylor basis functions alone is not evidence
of convergence.

### Wang-seed continuation research contract

This path is analysis-only and does not change `ipm.solve` or the production
solver. The authoritative stable seed value is the later
`lambda=1.0285722760323 +/- 1e-13`, with
`beta=lambda/(1+lambda)=0.507042459460253`; the first paper's
`1.0285722760222` is retained as a distinct historical table value. The
paper's printed `curl(U)=+H_X` line conflicts with Darcy's law, its explicit
`Omega` ansatz, and `H_X=-Omega`; every import must use
`curl(U)=-H_X=Omega` and recompute `-Delta(P)=H_X` locally.

No official code, checkpoint, profile array, parameter tree, or SI is
currently available. The public per-field architecture specifies only
`Omega=-(X/s)*exp(N_Omega)*q`; `H`, `U1`, and `U2` envelopes and the exact
normalization pin remain unknown. Raster figures are not machine-precision
seeds. A continuation run therefore cannot begin until an auditable `H`
array or executable parameter tree with normalization metadata is supplied.

After seed reprojection, symmetry breaking must use the full Poisson-coupled
Jacobian, split into even and odd reflection blocks. Horizontal drift `b` is
an unknown and a weighted translation phase borders the odd block so the
translation tangent is not reported as a fan mode. Branch switching is
allowed only after a refinement-stable zero singular value and transversality
check; otherwise a forced homotopy is only diagnostic unless its forcing
returns to zero. Exact-PDE pseudo-arclength Newton is a proposed next stage,
not an existing result.

The anisotropic contract includes `delta_s=2*(c_x-c_y)*delta`. Hence a
stationary finite nonzero `delta` forces `c_x=c_y`; fixed-aspect-ratio states
with unequal rates are frozen diagnostics, not two-scale fixed points. A true
two-scale search must retain the `nu*partial_kappa(R)` term with
`kappa=-log(delta)` and solve a shape-space orbit from the isotropic Wang
seed toward the degenerate fan chart. No such orbit has been computed.

The local compactified multistage surrogate does not relax this seed boundary.
Its first frozen-base correction changes the expanded-tail raw transport RMS
from `0.157614` to `0.158534`, although the factored combined RMS drops by
`6.73%`. The structured tail has the requested decay and zero sampled
wrong-sign energy, but finite-box Green velocity mismatches remain
`0.7268/0.7096` on `L=20/40`. It is an infrastructure regression and negative
screen only; a successor must encode the streamfunction/Poisson nonlocal
closure rather than infer it from local div/curl penalties.

### Polar-profile research contract

Profile analysis is read-only and stays outside the time integrator. The
candidate wall point is labeled `(1,0)`, while an existing double-odd result
may retain source center zero and be translated only for display. The pinning
identity is `c_r=-U_1(1,0)-c_l`; after subtracting the center velocity and the
density value at the point, the steady equation is

```text
(W+c_l z) dot grad Q = c_omega Q,
-Delta Psi = partial_z1 Q,
W = grad_perp Psi-U_1(1,0)e_1.
```

The exact polar system and its boundary conditions are maintained in
`analysis/POLAR_PROFILE_RESEARCH.md`. A wall log-log line is never accepted
as a profile by itself. `acceptedSingleScale` additionally requires agreement
between the wall and angular-L2 radial exponents and a small full-polar
collapse error. Fan angles are measured at several relative thresholds and
fitted to `theta_f(r)~r^sigma`; threshold-dependent `sigma` is evidence
against a settled single two-scale law.

For a local cusp `Q~r^beta F(theta)`, `0<beta<1`, nonlinear transport is more
singular than the modulation terms. Its leading angular equations imply
`F=C|G|^(beta/(beta+1))`, where `Psi~r^(beta+1)G`. With a continuous angular
profile and `Psi=0` on the wall this forces `F(0)=0`; it cannot produce a
nonzero wall `partial_x1 Q~r^(beta-1)`. The active alternatives are therefore
a shrinking angular fan `F(theta/r^sigma)`, a second inner/boundary scale, an
accepted angular discontinuity, or a transient wall power law. This is a
compatibility reduction for a pure homogeneous *local* cusp, not a blow-up or
nonexistence proof. The active general-`beta` construction instead has an
affine tip and uses the `r^beta` law at rescaled infinity.

The exact general-`beta` two-scale contract is

```text
rho = rho_c + L^beta R(s,(x-a)/L,y/L)
psi = L^(1+beta) P
L_t = -c L^beta,  a_t = b L^beta,  s_t = L^(beta-1)
R_s + (U+cZ-b e_1).grad(R) - beta*c*R = 0
-Delta P = R_X,  U = grad_perp P.
```

For `0<beta<1` and `L(T)=0`,

```text
L^(1-beta) = c*(1-beta)*(T-t)
rho_x = R_X/(c*(1-beta)*(T-t)).
```

Thus every bounded nonzero stationary `R_X` gives an exact Type-I physical
gradient rate; `beta` selects the length, density amplitude, and terminal
Hölder cusp, not the time exponent. The forced affine tip is

```text
R = A*X - A^2*Y/(2*(1-beta)*c)
P = (1-beta)*c*X*Y - A*Y^2/2
P_XY(0) = (1-beta)*c.
```

Its characteristics form `Y~X^((2-beta)/beta)`. The exact graph-hodograph
unknowns `J=X_q`, `K=X_Y`, characteristic speed `m`, and `R` satisfy

```text
J_Y=K_q,  (m*J)_Y=2*J,  m*R_Y=beta*R
(1+K^2)*m_q-J*K*m_Y+m*(K*K_q-J*K_Y)=-R_q/c.
```

This affine-cusp to smooth-tail heteroclinic is the primary fixed-profile BVP.

At rescaled infinity, a smooth angular tail
`R~r^beta F0(theta)`, `P~r^(1+beta) G0(theta)` obeys

```text
G0''+(1+beta)^2 G0 = -beta*cos(theta)*F0+sin(theta)*F0'
R_X=r^(beta-1)*(beta*cos(theta)*F0-sin(theta)*F0').
```

Because `1+beta` lies strictly between the first two Dirichlet frequencies,
the angular Poisson problem is uniquely solvable for every smooth `F0`
compactly supported in any prescribed right acute sector, and `R_X -> 0` in
all directions. The older hard fan `R=A*(q_+)^beta` fails this contract:
`R_X~q^(beta-1)` is singular on the edge and does not decay along paths with
fixed edge-normal distance. A narrowing admissible edge must be soft; for
normal width `w(r)=r^gamma`, `0<gamma<1`, its gradient decays like
`r^(-gamma*(1-beta))` while its angular width tends to zero. Representing this
requires angular complexity growing like `r^(1-gamma)`; fixed finite angular
modes cannot close the desired profile.

For the wall trace `f=R(X,0)`, define

```text
Ahat[R](k)=|k|*integral_0^infinity exp(-|k|Y) Rhat(k,Y) dY.
```

The exact trace equation is

```text
f_s+(c*X-b+H(A[R]))*f_X-beta*c*f=0.
```

It reduces to CCF only under the additional DtN/columnar condition
`A[R]~=f`. Within that reduced model, `theta=-f` and
`lambda=beta/(1-beta)` map the standard CCF similarity law exactly. CCF
finite-time blow-up from smooth one-dimensional data is rigorous; its lift to
two-dimensional IPM is not. The direct numerical task is to control or measure
`A[R]-f` across the evolving normal layer.

The exact degree-one acute-fan benchmark uses centered coordinates
`X=x_1-1`, `Y=x_2`, `q=sin(theta_*)X-cos(theta_*)Y` and

```text
R=A_0 q,
Psi=A_0 sin(theta_*)/(2 cos(theta_*)) Y q,
a_0=A_0 sin(theta_*)^2/(2 cos(theta_*)).
```

It satisfies `U dot grad R=-a_0 R`, so the constant-rate steady condition is
`c_l-c_omega=a_0`. The mesh benchmark selects the spatial-only gauge
`c_l=a_0`, `c_omega=0`. Here `c_omega` scales density; the physical gradient
scale grows at rate `c_l-c_omega`. In coordinates centered at `(1,0)`,
`c_r=0`; in the original unshifted coordinate with the point fixed at
`x_1=1`, the equivalent value is `c_r=-c_l`. This benchmark is an exact
infinite-energy fan-domain/distributional solution, not a claim that the
localized degenerate-data run selects it.

The polar research contract also records an exact far-decaying family for
the unrescaled physical steady equation. In centered polar coordinates, set
`rho=lambda*Psi`, `k=|lambda|/2`, and

```text
Psi_nu = exp(-lambda*r*cos(theta)/2) K_nu(k*r)
         * (a_nu*cos(nu*theta)+b_nu*sin(nu*theta)),
rho_nu = lambda*Psi_nu.
```

Then `U dot grad rho=0` and `-Delta Psi=partial_X rho` exactly on the
punctured upper half-plane. The angular order may be any real `nu>=0` when
wall traces are unrestricted; `sin(n*theta)`, integer `n>=1`, gives the
optional homogeneous impermeable-wall specialization. Modes with the same
`lambda` may be linearly superposed. For fixed interior angle the tail is
`f(theta) r^(-1/2) exp(g_lambda(theta) r)`, with
`g_lambda=-(|lambda|+lambda*cos(theta))/2<=0`; on the single wall ray where
`g_lambda=0`, a nonzero trace still decays algebraically. Center singularity
is accepted for this research family. It is not a `c_l~=0` rescaled steady
profile and is not part of the active time integrator.

The maintained right-wedge visualization takes `lambda<0` and
`Theta=cos(pi*theta/(2*theta_*))`, so the trace is nonzero on the positive
horizontal wall, zero on the acute fan ray, algebraically slow along the
wall, and exponentially damped away from it. The comparison sine mode
vanishes on both wedge rays. These are exact solutions inside the open wedge;
their zero extensions are distributional fan-domain fields unless the normal
derivative is matched at the fan ray.

The distinguished global member instead takes `lambda=-mu<0`, `nu=0`:

```text
Psi_0 = -A exp(mu X/2) K_0(mu r/2),
rho_0 =  mu A exp(mu X/2) K_0(mu r/2).
```

It is positive and exact on the punctured plane, with far field
`rho_0~A*sqrt(pi*mu/r)*exp(-mu*r*(1-cos(theta))/2)`. Its angular width about
the positive horizontal ray is `(mu*r)^(-1/2)`, so it is a global solution
that becomes asymptotically concentrated in every fixed right-facing acute
sector rather than a hard wedge cutoff. The `K_1*sin(theta)` member gives the
optional homogeneous-wall upper-half-plane version. Replacing `r` by
`sqrt(X^2+(Y+h)^2)`, `h>0`, moves the pole below the wall and yields an exact
interior solution smooth throughout the upper half-plane, but its varying
wall streamfunction gives nonzero normal velocity; it is not an impermeable
half-plane solution.

These global Bessel members are dynamic-scaling fixed profiles only for
`c_l=c_omega=c_r=0`. With nonconstant coordinate scales they trace an exact
but nonstationary Bessel-parameter orbit `mu_rescaled=mu/C_l`; their intrinsic
length rules out a nonzero-`c_l` fixed profile. They establish a steady outer
geometry, not attraction from smooth data. Any such selection claim must
retain a bounded inner core, respect the physical density maximum principle,
and demonstrate convergence in time and resolution.

There is nevertheless a distinct anisotropic far-field limit. With
`X=C_x*x`, `Y=C_y*y`, `R=C_rho*rho`,
`P=(C_rho*C_y^2/C_x)*psi`, and the advective clock
`dt/ds=C_rho*C_y/C_x^2`, define `delta=(C_x/C_y)^2`. The exact transformed
Poisson equation is

```text
-(delta P_XX+P_YY)=R_X.
```

On the Bessel branch, `C_y^2=mu*C_x` and `P=-R`. The centered K0 and K1
members then converge as `epsilon=mu/C_x -> infinity` to, respectively,

```text
R0_*=A sqrt(pi/X) exp(-Y^2/(4X)),
R1_*=A sqrt(pi) Y X^(-3/2) exp(-Y^2/(4X)),   X>0.
```

Their nontrivial `R_X=0` curves are `Y^2=2X` and `Y^2=6X`. This is a
spatial downstream heat-kernel blow-down: for fixed `mu`, both physical
lengths expand (`L_y~sqrt(L_x)`), so it is not a time-contracting cusp.
The general `delta->0` IPM system retains order-one nonlinear transport and
does not become the heat equation unless the extra closure `P=-R` is selected.

The decaying `K_nu` branch has a center pole, while the center-regular
modified-Bessel branch grows at infinity. More generally, the Darcy
dissipation/potential-energy identity makes every sufficiently decaying,
smooth, finite-dissipation physical steady state with an impermeable wall
trivial. For a globally decaying rescaled steady state, an independent
necessary condition is

```text
(2 c_l+c_omega) integral Q = 0.
```

Hence `c_omega=0`, `c_l~=0`, nonnegative finite mass, and global decay are
not compatible. A target bounded-density gradient profile must use zero
signed mass, `c_omega=-2 c_l`, a nonintegrable/matched outer field, or a
genuinely multiscale inner/outer construction. The Bessel modes are exact
exterior benchmarks, not such a cusp profile by themselves.

There is now a stronger local correction. Circle flux for K0 and the
half-disk normal-dipole moment for K1 imply that, when `q=mu*r << 1`, any
bounded smooth-core/closure correction is larger than the corresponding
Bessel density by at least `1/(q*log(1/q))` or `1/q`, respectively. Thus a
small-argument Bessel singularity cannot be the leading layer on
`ell << r << 1/mu`. The Bessel response can begin only at its natural scale
`r~1/mu`; a viable dynamic construction needs a same-order closure/wall
layer and a nonlocal source-neutral return/reservoir layer. The former
amplitude-only K0/K1 gluing laws are superseded by this multipole audit.

The same audit now contains an exact dynamic affine normal form. After fixing
the physical wall tip and writing `X=x_1-x_*`, it is

```text
psi = Y*(X-k*Y)/(T-t),
rho-rho_c = 2*k*(X-k*Y)/(T-t),  k>0.
```

It has `rho_X=2*k/(T-t)>0`, an impermeable wall, a right-spread positive
sector, a hyperbolic tip, and an exact odd weak extension with a wall jump.
It is unbounded and infinite-energy, so it is a local normal form rather than
a bounded-data blow-up theorem. Moreover, its apparent self-similar exponent
`beta` is a rescaling gauge: the physical affine formula is independent of
`beta`. Any smooth nondegenerate wall fixed point is nevertheless forced to
have this affine first Taylor jet.

Without local `x_1` reflection symmetry, the leading homogeneous
`r^beta` Poisson/free-boundary tail can have a supported right fan only at
`theta_*=pi*beta/(1+beta)`. No reliable data or solvability condition currently
selects a numerical `beta`; it remains a matching eigenvalue. This angle is a
leading-order condition inside the hard-edge ansatz, not an admissible
all-orders profile: in addition to the first nonlinear correction of degree
`2*beta-1`, the hard edge itself violates all-direction `R_X in C_0`. A smooth
angular tail or radius-dependent soft edge must replace it and generally needs
a signed return/edge layer. Direct K1 matching has the bounded-amplitude window
`ell <= L_B <= ell/delta`, with `delta=ell^beta`; the maximally separated
endpoint requires an order-one nonconstant compensating moment and may change
the inner dynamics unless its induced strain cancels. The open theorem is an
affine-to-cusp heteroclinic with closure/wall and area-preserving Casimir
return layers, not a shrinking Bessel steady state.

The active background-field formulation is exact on a fixed-source or
punctured-domain Bessel contract. With `D_B=0`, perturbation
`d=eta+mu*phi`, and
`(Delta-mu*partial_x)(T_mu d)=-partial_x d`, reconstruct

```text
psi = psi_B + T_mu d,
rho = rho_B + (I-mu*T_mu)d,
(I-mu*T_mu)d_t + (U_B+U_d).grad(d) = 0.
```

This operator is asymptotically local on an inner scale `ell << 1/mu` and
becomes the full Bessel Green response at `r~1/mu`, so it is the preferred
elliptically exact inner/outer bridge. It is a singular perturbation: the
multipole bounds and nonzero wall trace force `d` to be order one in the
inner region. A source-free physical solution must additionally cancel the
Bessel tip distribution with a leading wall/return particular solution.
The obvious amplitude, translation, and `mu` tangent modes are neutral, while
any growing integrable mode must have zero weighted mean on almost every
Bessel density level set because all Casimirs are conserved. Positivity of
the real quadratic form of `I-mu*T_mu`, together with skew transport, also
rules out every nonzero real eigenvalue on the pure `D_B=0` background under
source-neutral no-flux conditions. Any pure-background instability must be a
genuinely complex mode; a regular composite `D_B!=0` background may recover
real growth through its extra linear coupling. The next decisive calculation
is therefore the compressive spectrum of the regular composite background,
followed by order-one nonlinear continuation.

Two further exact audits constrain that continuation. No finite polynomial in
the wall-normal coordinate can connect the forced affine wall jet to a
general `r^beta` tail; a genuine bridge needs infinitely many angular modes or
a non-polynomial wall layer. A single-sign, purely homogeneous, sheet-free
acute fan is also impossible as a complete nonlinear fixed point; the
`r^beta` fan can only be a leading outer asymptotic with lower-degree,
edge/harmonic, or signed-return corrections. There is an exact 90-degree weak
heteroclinic
`P=c*y*H(x), R=-c*y*H'(x)` connecting an affine inner asymptotic to a
homogeneous outer cusp for every `0<beta<1`, but it has zero wall trace,
insufficient tip regularity, and infinite energy. It is a useful near-miss,
not the target bounded-data fan.

## Result format

`ipm.solve` returns final rescaled/computational fields (`rho=R`, `omega=d_X1
R`, `rescaledOmega`, `psi`, `u1`, `u2`) plus reconstructed
`physicalRho`, `physicalOmega=(C_l/C_omega)d_X1 R`, and
`physicalU=U/C_omega`, grid, physical `t`, canonical rescaled `tau`, internal
`normalizedTau`, `C_l`, `C_omega`, `X_shift`, physical grids, options, history,
snapshots, stop reason, and candidate blow-up fit. In physical mode the scales
and both rescaled clocks agree. `finalQuality` gives the terminal safety,
density range, wall-profile oscillation/sign, adjacent-cell ratios, and
two-axis core counts in one compact audit structure. Sparse operators and
factorizations are not saved. Default physical and dynamic output files are
`result/production/ipm_double_odd_physical.mat` and
`result/production/ipm_double_odd_dynamic.mat`, so one mode does not silently
overwrite the other. Active-mode `X_shift` and all `c_r` histories are zero by
contract.
For anisotropic output, `C_x,C_y,L_x,L_y,aspect`, split canonical/advanced
`c_x,c_y`, `aspectRate`, the exact conservative source, aspect-tagged Poisson
solve information, axis-specific physical derivatives/velocities, and split
physical gradient energies replace scalar-scale inference. Its history retains
the common peak-frame and condition slots for result-schema compatibility; a
NaN in a nonselected gauge is intentional. For `peak_translation`,
`peakTranslationCXCoefficient`,
`peakTranslationCRCoefficient`, `peakTranslationRank`, and
`peakTranslationCondition` still describe the legacy transport row
`[X_*,1]`.
Physical grids include laboratory `physicalX`, peak-aligned
`comovingPhysicalX`, `physicalY`, and `physicalTrackedPeakX`.
Full field snapshots are disabled by default to avoid more than a gigabyte of
duplicate data in a standard long video run; set `storeSnapshots=true` to save
them. When enabled, `snapshotTimes` uses physical time, `snapshotTau` uses
canonical rescaled time, `snapshotNormalizedTau` exposes the internal clock,
and the `snapshotX`/`snapshotY` cell arrays store the contemporaneous
computational grid for every field. Post-processing must use these per-frame
grids after adaptive remeshing; a legacy multi-snapshot remeshed result without
them is rejected rather than silently differentiated on the final grid.
History includes both
clocks, physical and rescaled mass, area/support/core resolution counts, the
additive width correction, its local rate scale, and control mode,
the hard two-axis safety factor, accepted remesh count, minimum `dX`, minimum
`dY`, peak and origin `dX`, maximum adjacent-cell ratios, diagnostic global
max/min spacing ratios, horizontal and wall-normal core counts and widths,
origin velocity plus `rho`-even/`omega`-odd defects,
full-gradient energy/effective length and scale
fractions, FWHM plus `trackedWallCoreWidth`, point-strain conditioning,
wall-support local-maximum count, total-variation ratio, positive-side
negative-omega leakage, bounded and canonical scaling rates,
the common time speed, and all gauge residuals.
Options retain `gridStretchAutomatic` and the selected `transportScheme`;
`physicalFinalTime` and `finalTime`
record the physical/canonical stopping contracts. With open artificial
transport boundaries, physical mass is a moving/truncated-box diagnostic
rather than an invariant. `blowupFit.trustedThroughTime` is the last fitted
state satisfying resolution, density-range, one-peak/total-variation,
positive-wall-sign, and maximum-cell-ratio contracts;
under-resolved terminal states are excluded from every reported fit. The
primary fit uses 40 points, while `windowPoints`, `windowT`, `windowP`,
`windowR2Power`, `windowR2Exponential`, and `windowCandidate` record the
20/40/80/160-point stability audit when those windows are available.
The explicit historical CCF long-mesh audit additionally exports selected history
columns and an explicit `trusted` mask to
`result/verification/ccf_tail_long_physical_history.csv`; this CSV is a
derived review artifact, while the MAT result remains authoritative.
The optional polar analyzer returns a separate report with `wall` fit data,
`polar.angularL2Fit`, thresholded `polar.fanPower`,
`polar.separabilityError`, the exponent translation `beta=1-alpha`, and
`acceptedSingleScale`. These analysis reports do not alter or replace a
`ipm.solve` result.
The wall-interval analyzer returns a separate `diagnostics` structure with the
selected/trusted indices, snapshot-1 material labels and jump, endpoint and
width histories, three independent compression representations, material-
speed defects, moving-peak label and drift rows, limiting-label fit, nested
fixed-label audits, and explicitly conditional collapse deadlines. Its
`densityJumpRelativeDrift` and `integratedQRelativeToInitialJump` are always
normalized by snapshot 1; later snapshots never redefine the invariant.
The anisotropic continuation analyzer similarly returns an independent
`audit` with immutable input contract, history/quality gates, per-snapshot
scale, clock, Riccati, Green, Poisson, and profile-drift rows, plus an optional
sheared-rank block. Neither structure mutates or replaces the authoritative
`ipm.solve` result, and neither may silently fill a report placeholder after a
failed gate.
The characteristic-oriented `A`-null fan search returns an independent
analysis structure containing `options`, `basisAudit`, normalized basis
labels/scales, all regularization/penalty `runs`, `paretoIndices`,
`bestIndex`/`best`, and per-grid `refinement` metrics. Its broad, narrow, and
tip MAT/PNG pairs live under `result/verification` and are evidence artifacts,
not production results or accepted profiles. Every summary must preserve the
absolute/reaction/`R_X` normalizations, outside-fan fraction, maximum residual,
and explicit `noncandidate` classification; the narrow PNG is the
authoritative current visualization. Its first-panel title has been corrected
to `characteristic-oriented A-null R_X`; the underlying basis remains the
straight `X-sigma*Y` parameterization.
The exact Bessel analyzer similarly returns its polar grid, analytic fields,
per-mode PDE residuals, and per-edge normal-flow residuals without entering
the simulation call graph.
The near-Bessel analyzers return fitted `mu` (and K0 pole depth), normalized
cross-section and full-kernel errors, parabolic width exponents, `rho_x`
zero-line errors, initial/final field changes, and gradient ratios. Matched
pair comparison reports full-field and tail difference ratios. Screen MAT
files and diagnostic PNGs live under `result/verification`; they are derived
evidence and may be regenerated from `bessel_two_scale_case` through the sole
simulation entry point `ipm.solve`.
The default physical run writes
`result/production/ipm_double_odd_physical.mp4`; an explicit dynamic run writes
`result/production/ipm_double_odd_dynamic.mp4`. Each updates
a live figure at every recorded output time. Rate panels use canonical
rates against canonical time, so the common limiter cannot hide genuine raw
`c_l`, `c_omega`, or `c_r` growth; `timeSpeed`, safety factor, accepted remesh
count, and peak spacing remain visible in the live/video titles.

## Simplicity rules

1. Keep generic explicit-option completion physical and free of every
   rescaling term; keep the maintained no-argument dynamic experiment isolated
   in `ipm.config.activeCase`.
2. Keep legacy Boussinesq code out of the active call graph.
3. Share the physical numerical operators between both modes.
4. Prefer explicit structures and pure functions over stateful classes.
5. Keep plotting and fitting outside the numerical core.
6. Treat larger domains and finer grids as explicit verification runs.
7. Keep active IPM source in the root, reproducible option sets in
   `experiments/`, accepted outputs in `result/production`, convergence and
   regression data in `result/verification`, rejected/history data in
   `result/archive`, and the old solver in `legacy_boussinesq/`.

## Change log

> **Historical chronology.** Every controller/rate-limit entry below predates
> schema 4. It records what was tested at the time and does not reinstate
> `maxDynamicRate`, adaptive width feedback, or nonzero restoring gains.

- 2026-08-30: Established a single option-resolution boundary on the maintained
  solver path. `ipm.solve` now resolves once, pure operator construction and
  remesh reuse that frozen flat structure, and legacy partial-option callers
  remain supported by wrappers. Rescaling initialization now preserves every
  configured field while overlaying initial-field tracking state, so repeated
  initialization no longer needs a manual field repair. Anisotropic gauge
  configuration now exposes only the implemented `fixed` and half-plane
  `peak_translation` modes; option resolution rejects all other values.
- 2026-08-29: Enlarged the maintained initial physical box from
  `[-256,256] x [0,256]` to `[-1e5,1e5] x [0,1e5]` and raised
  `physicalFinalTime` from `0.5` to `2`, without changing the `513 x 257`
  point count or the analytic datum. The automatic lattice keeps center
  spacing `1e-2`; its maximum adjacent-cell ratio is about `1.056`, below the
  active `1.5` bound. `ipmtests.baseline.suite` and Code Analyzer pass.
- 2026-08-29: Reduced the no-argument active grid from `1025 x 513` to
  `513 x 257` for interactive turnaround, while retaining the
  `[-1e5,1e5] x [0,1e5]` tail box. Peak/level targets are now `64` and
  `[96,64,48]`, hard floors are `12` and `[16,10,8]`, and at most six analytic
  pre-remesh passes are used. Quantitative pointwise claims still require the
  explicit matched `1025 x 513` refinement.
- 2026-08-29: Raised the active canonical `maxDt` ceiling from `1e-4` to
  `1e-2`. A zero-step audit on the remeshed `513 x 257` datum gives
  `dt_CFL=1.9415e-3` in the internal clock (`1.0070e-3` physical) at
  `cfl=0.35` for the then-active bounded rational comparison, whereas the old
  cap allowed only `1e-4` physical. The maintained logarithmic-primitive case
  likewise chooses its step stagewise from CFL; a smoke run on the former
  `[-256,256] x [0,256]` box reached physical time `0.005` in three steps.
  Thus `maxDt=1e-2` is a safety ceiling,
  not a prescribed fixed step.
- 2026-08-29: Changed the no-argument `ipm.solve()` datum from the regularized
  CCF branch, first to the bounded rational comparison and then to
  `initialCondition='degenerate_primitive'`, `degeneratePower=8`, namely
  `rho_0=log((1+|x|^8+y^8)/(1+y^8))/8`, on `[-1e5,1e5] x [0,1e5]`. Its exact
  `q_0=rho_x=sign(x)|x|^7/(1+|x|^8+y^8)` has one sign on each half-axis,
  seventh-order finite origin degeneracy, and both generic-ray and wall
  inverse-distance decay. Thus `q_0 -> 0`, while the full-gradient `L2` energy
  and the wall primitive have logarithmic far-field growth. The latter is the
  unavoidable cost of a one-signed `1/x` tail. The active gauges are the exact
  `transport_anchor` at `X=1` with `lengthScaleGain=1` and
  `cOmegaGauge='wall_omega_peak'`. A nonempty partial option structure still
  receives the generic `k=4` fallback. CCF data and their gradient-energy
  normalization remain reproducible only through explicit named experiment
  options; they are no longer silently selected by `ipm.solve()`. The exact
  anchor residual is zero and its initial canonical slope is approximately
  `-0.182`; hence `X=1` is attracting for this datum. A separate repelling zero
  near `X=1.333` initially bounds the anchored transport basin. New full
  positive-half-plane minimum and negative-ratio diagnostics expose any later
  loss of the initial `rho_x>0` sign without modifying the evolved field.
  On the former `[-256,256] x [0,256]` box, a physical-time `0.005` smoke run
  takes three steps, keeps the negative-ratio
  diagnostic at zero, changes `c_l` from `0.487243929` to `0.487306433`, and
  changes the physical wall peak from `0.686006441` to `0.687247420`; this is a
  consistency check only, not blow-up or two-scale convergence evidence.
- 2026-08-29: Replaced the then-maintained `symmetry_peak` half-strength length
  condition by `lengthGauge='transport_anchor'` at the fixed computational wall
  point `X=1`. Double-odd symmetry still enforces `c_r=X_shift=0`, so the
  canonical equation is `c_l=-U_1(1,0)`; validation fixes
  `lengthScaleGain=1` and excludes adaptive width-rate corrections. The common
  positive `timeSpeed` preserves the zero exactly. New canonical/bounded anchor
  diagnostics are stored in flow/history and verified initially and after an
  SSPRK3 step. On the historical post-remesh CCF datum,
  `timeSpeed=0.579556`, `c_l=0.5`, `V_1(1,0)=3.19e-17`, and
  `d_X V_1(1,0)=-0.961733`. The wall velocity is positive at `X=0.5`, negative
  at `X=1.239698`, and returns to zero at `X=1.479396`; hence
  `[0.5,1.239698]` is an explicit inward-pointing support witness. This gauge
  no longer claims to freeze the moving wall-`rho_x` peak. `ipmtests.baseline.suite` and the
  MATLAB Code Analyzer pass.
- 2026-08-29: Refocused the then-active blow-up analysis on the exact no-argument
  `ipm.solve()` CCF heavy-tail datum and made the generic `ipm_defaults`
  fallback distinction explicit. Corrected the far-field contract: the
  finite-regularization datum is `C1`, has `rho_x` continuous/Lipschitz and in
  `C0`, but retains a fixed-`y` `rho_y` reservoir, so its gradient energy grows
  linearly with horizontal box length and is `99.7034%` normal-gradient energy
  on the maintained box. Added three bounded, read-only diagnostics: the
  current-WENO zero-step wall Riccati audit, the physical density-label inverse
  map/fold audit, and a one-RHS tangent audit of horizontal and vertical
  widths from one retained post-initialization RHS. The latter gives initial
  FWHM `gamma_x/gamma_y = 2.3474/0.3803`,
  while its aspect-rate direction disagrees with the old pre-WENO long
  transient; these are therefore recorded as local derivatives, not
  similarity exponents. The same tangent accepts `lambda_c'=-0.621263` but
  rejects `B'`, `deltaLambda'`, and `chi'` because their fit-level rates vary
  by more than 100% and `B'` changes sign. The associated theory target is the
  normalized right-minus-left Green moment: its odd direct kernel explains
  why finite edge coupling can supply leading compression while a symmetric
  interior fold cancels. The corresponding proof target now has two distinct
  gates on a causally fixed material label: a positive current core-tail
  margin and a first-exit (or tube) production inequality. The current margin
  alone is not a closure. No PDE step, parameter scan, or solver path was
  added by these audits.
- 2026-08-29: Added a zero-step boundary-velocity plot for the exact maintained
  CCF initial grid and corrected it to use the dynamic transport field rather
  than physical `u` alone. The later fixed-`X=1` gauge above supersedes the
  first negative-to-positive version: the current plot shows an attracting
  positive-to-negative anchor, a trough `-0.102030` at `X=1.168891`, and the
  outer zero `X=1.479396`.
  At the wall-`rho_x` peak, physical/dynamic strains are
  `1.432211/0.330046`; only the former enters the physical Riccati coefficient
  `a=0.653291`. This remains an initial coordinate/elliptic-contract check,
  not two-scale or blow-up evidence.
- 2026-08-29: Generalized the exact material-wall collapse theorem from the
  two zero endpoints of a complete positive component to every pair of fixed
  density labels inside that component. In label coordinates,
  `X_t=v`, `X_ell=1/Q`, `Q_t=a Q^2`, and the interval coefficient is exactly
  the label average of `a`; a moving maximum instead has the independent drift
  `ell_*'=v_xx*q_*^2/q_xx`. Recorded the conditional `a_x>1` limiting-label
  route, the collapse of core label mass
  `q_* L_x ~ (T-t)^(a_x-1)`, and the resulting noncommutation of the
  `f -> 0` and `t -> T` limits. The limiting-label interpretation is restricted
  to one non-switching `C1` peak branch; its physical-offset estimate also
  requires `Q >= c*q_*` on the connecting label segment and gives only an
  `O(1)` normalized offset, not convergence. The full-component deadline near
  `1.69` assumes the observed tail lower bound persists into the uncomputed
  future and is retained as a weak sufficient screen, not a prediction of the
  candidate time near `0.772`. Added the read-only material-interval and
  anisotropic-continuation analyzers to the file/result contracts. Also made
  explicit the horizontal-zero-mode no-go: bare pointwise reflection
  positivity is not a universal forward-invariant cone without bulk
  production, signed-tail, and Casimir-leaf control. No solver path changed.
- 2026-08-28: Replaced the rejected single-rectangle interpretation by a
  quantified reduction hierarchy. A twelve-slab separable model still fails
  validation (`0.614/0.382/0.205` locally), whereas sheared left/right traces
  at slope `s=0.30` are stably rank 2--3 across the 385, 513, and wide 513
  runs: outer-window rank-2/rank-3 shape errors are about `4.8%/1.55%`,
  rank-2 Green-feedback errors are `4.4%--5.2%`, and cross-grid subspace angles
  are `0.2--1.6` degrees. Added the exact reflection-difference audit. At
  `t=0.74` the complete-box negative/positive feedback ratios are
  `0.0312/0.0377/0.0355`, with the wide-box ratio decreasing from `0.1477` at
  `t=0.65`; the residual negative part is a far return/cap. Documented the
  aspect-uniform tail lemma
  `a_tail^- <= C/(pi*delta) R^(-delta)` and the independent dynamic Riccati
  rate defect `(c_x-c_omega)-a_h*Omega_h`. These results identify a
  reflection-positive sheared two-trace core plus a moving 2-D reservoir as
  the next closure target, but do not establish forward invariance.
- 2026-08-28: Added the opt-in anisotropic dynamic-coordinate prototype with
  exact `A_kappa`, determinant-one Green metric, stage-local SSPRK aspect,
  split physical reconstruction/diagnostics, and an unchanged isotropic
  default path. Added fixed-rate and half-plane peak-translation modes, with
  the explicit rank-one `[X_*,1]` no-go for jointly determining `c_x,c_r`.
- 2026-08-28: Expanded the tuned smooth signed-slope candidate to
  `[-640,640]x[0,128]` at `513x257`. Its last trusted time remains
  `0.750040037`, with global-`rho_x`/wall/full-gradient growth
  `15.2735/39.8590/18.5965`, zero range violation, far-source ratio
  `7.9e-6`, and `-0.93%` signed-mass drift relative to initial L1. At common
  `t=0.74`, standard/wide-box peaks differ by `0.61%`, half-widths by `1.67%`,
  peak location by `9e-6`, and normalized two-scale profiles by `0.64%`.
  Under the stricter six-points-in-both-cores gate, standard/wide Type-I times
  are `0.773187/0.772385`, peak exponents `1.0420/1.0384`, and width exponents
  `1.6945/1.5336` versus `1.6918/1.5523`. Added a direct wall-Riccati audit:
  the three main grids/boxes give positive median
  `-u_{1x}/q_peak=0.260--0.270`, instantaneous Riccati times consistent with
  the reciprocal fit, and only `2%--3%` median semi-discrete wall-identity
  defect. The exact Green-kernel pairing yields a conditional wall blow-up
  criterion and a formal `kappa=L_x/L_y -> 0` wall-Hilbert functional; at
  common `t=0.74`, actual/wall-limit coefficients differ by
  `7.95%/5.88%/0.37%` on coarse/fine/wide runs. A linear source budget is
  independently stable: right/left positive-`q` blocks contribute about
  `+3.6/-2.6` times the net strain, `|x-x_*|<L_x` supplies about `61%`, and
  `y>16 L_y` still supplies about `12%`. A stronger same-Poisson reduction to
  one rectangular Kiselev--Sarsam `H_h` layer fails: at wide-box `t=0.74`, the
  best local strain/peak-frame-velocity/wall-evolution relative errors are
  `0.868/0.900/0.300`, and fitting only the peak scalar makes profile errors
  worse. This rules out finite-box forcing and a one-layer closure, and
  identifies a core-plus-reservoir 2-D strain feedback, but does not prove its
  analytic positive lower bound.
- 2026-08-28: Continued the optimized beta-one-third tail through the next
  `(n,m)=(11,6)` resonance. The new resonant log coefficient stays near
  `1e-15`; the actual limitation is nonresonant broad-spectrum growth.
  `N=11` residuals at `r=2,3,4,5,7,10` are
  `0.02036/0.001038/0.0001257/2.445e-5/2.071e-6/1.512e-7`: N8 is optimal at
  `r=2`, while N11 still improves for `r>=3`. A four-point large-order screen
  is compatible with, but does not identify, Gevrey-1 behavior and an
  empirical `exp(-7 r^(2/3))` optimal-remainder scale. Added regenerated N8
  and N11 trace APIs; signed trace labels require sign charts.
- 2026-08-28: Opened a quantile signed-Burgers core and a Casimir-null
  Hamiltonian collar against both N7 and N8 tails. The core can match the full
  one-dimensional level-area distribution exactly, disproving a Casimir-value
  no-go. Nevertheless, `44.1%` rank inversions and `60--64` tail turns violate
  endpoint order for every orientation-preserving multiplicity-one
  boundary-to-boundary relabeling. The surviving five-dimensional collar
  Schur space has Fredholm orthogonal ratio about `0.75` and improves the
  validation residual by under `1%`; the next admissible freedom is a 2-D
  reservoir/braiding or moving interface, not another 1-D relabeling.
- 2026-08-28: Replaced exact `S5=S8=0` resonance cancellation by constrained
  optimization of the real `N=8` overlap residual. Simultaneous exact
  cancellation is locally reachable but worsens the nonresonant coefficients;
  the accepted small-nonzero-log branch instead improves the independent
  `8001x1800` residual by about `47.5x` for every tested `r>=3`, gives
  `0.014741` already at `r=2`, and reduces rather than inflates `F4--F8` RMS.
  A near-skew/Hamiltonian generator supplies efficient tangent coordinates but
  is not a hard equation for the best branch. The new tail changes full gluing
  traces by `35%--37%`; old inner/DtN results cannot be reused. Continuation
  through the next `(n,m)=(11,6)` resonance is required before any global claim.
- 2026-08-28: Replaced the invalid `K<=0` profile no-go by the invariant
  fold-aware Lagrangian criterion. Two generating charts can cross a projection
  fold while the full Jacobian, canonical area, Poisson solve, and Stokes seam
  remain regular. A 28-parameter dual-chart gluing screen passes invertibility,
  multiplicity, Casimir, seam, and all-direction-decay gates but lowers the
  independent residual only `0.29%` to `2.52591`. The fixed-label transport law
  makes the affine core density mismatch invisible to every pure geometry
  mode, so the next Schur system must open signed-Burgers core density and
  Casimir-null Eulerian collar sources rather than add chart modes.
- 2026-08-28: Fixed adaptive full-field provenance by storing `snapshotX` and
  `snapshotY` beside every saved field; remeshed legacy histories without
  per-frame grids are now rejected by field-level analysis. Added a regression
  for snapshot/grid counts, shapes, and the terminal grid. Using that contract,
  promoted the initial-growth-tuned analytic signed-slope datum to a strong
  smooth numerical blow-up candidate under the hard all-direction
  `rho_x -> 0` condition. Analytic pre-remeshing raises the initial vertical
  90%-core from `1.32` points to `8.34/10.02` on `385x193`/`513x257` grids.
  With the default conservative wall update, the last trusted fine state at
  `t=0.750041767` has global-`rho_x`/wall/full-gradient growth
  `14.8297/38.9072/18.7215`, zero range violation, `6.72/4.95` horizontal/
  vertical core points, and `0.51%` signed-mass drift relative to initial L1.
  Coarse/fine Type-I wall fits on `t>=0.55` give `T=0.772025/0.771813`;
  separate width exponents are `1.710/1.551` and `1.776/1.607`, while
  normalized two-scale profiles differ by about `1%--4%` at common times.
  Halving both maximum time step and CFL on the coarse grid gives
  `T=0.771795`, width exponents `1.717/1.564`, and only percent-level growth
  changes at common time.
  Final resolution-failed spikes are excluded. Nonuniform WENO5 fails the
  maximum-principle gate at `t=0.01`, so it is not supporting evidence.
- 2026-08-28: Closed the rational far-tail recurrence in the log-polynomial
  transseries class. At `beta=1/3` the exact resonances are
  `(n,m)=(5,2),(8,4),(11,6),...`; at `beta=1/5` the first is `(4,2)`.
  Corrected the transport coefficient back-substitution from `p` to `p+1`,
  added absolute/relative coefficient defects, and recomputed every quoted
  number. A same-seed `5001/1000` versus `8001/1800` audit accepts only about
  `N=7--8`: `N=7,r=3` has residual `0.07123` and `N=8,r=4` has `0.01737`.
  The apparent coarse-grid `N=20` core penetration is withdrawn because the
  fine-grid coefficients and residuals diverge. Also registered the exact
  Casimir-cap scaling and conditional level-wise Hamiltonian gluing, together
  with the one-mode horizontal-overlap failure and the resulting need for a
  fold-aware multi-chart core/cap fixed point.
- 2026-08-28: Upgraded the all-direction-decaying smooth-angular tail from a
  leading-order ansatz to a six-level nonlinear far-field recurrence. A
  deterministic weak return using `0.3970%` outside-fan `Q0^2` energy cancels
  the cap-drift and first near-resonance moments, suppresses the first stream
  correction by about 3800, and yields decreasing full transport residuals
  for `r>=30`. Registered the beta-one-half `r log(r)` Fredholm limit, the
  non-positive-definite resonance-moment caveat, and the high-order small
  divisor. Added analytic Schwartz right-fan initial data and honest
  grid/box/angle scans; the simple Gaussian family has early compression but
  no two-scale separation and does not select 38.357 degrees.
- 2026-08-28: Added the straight characteristic-oriented `A`-null fan audit
  and integrated its broad/narrow/tip evidence without changing the production
  solver. The exact-DtN basis uses
  `zeta=kappa+mu*abs(k)+i*k*sigma` and follows `X-sigma*Y`; it is not a
  curved-characteristic parameterization. Finest-grid honest
  absolute/reaction/`R_X` residuals are
  `0.359438/0.249961/0.446312` (broad),
  `0.299396/0.208582/0.374605` (narrow), and
  `0.314621/0.219182/0.390475` (tip), with outside-fan fractions
  `0.008242/0.005416/0.005631` and maximum residuals
  `5.5732/5.3818/5.1617`. The narrow coefficient norm is `428`, and adding
  `n=5` near-tip modes regresses. All three outputs are explicit stationary
  noncandidates; the result is only a Pareto improvement that preserves
  wall trace, first jet, and DtN while keeping fan leakage below one percent.
  Added the analysis script/audit to the file map and documented the
  independent MAT/PNG result contract; full-space capping, genuinely curved
  characteristics, and multi-chart bulk gluing remain open.
- 2026-08-28: Promoted the Wang--Leger upper-half-plane IPM profiles to the
  strongest existing unforced smooth self-similar numerical benchmark under
  the active `R_X -> 0` contract. Their isotropic ansatz
  `rho=(1-t)^lambda*H(x/(1-t)^(1+lambda))` gives the exact candidate identity
  `rho_x=-(1-t)^(-1)*Omega`, while the numerical envelope hard-codes
  `Omega=O(r^(-1/(1+lambda)))`. The later stable value is
  `lambda=1.0285722760323 +/- 1e-13`, mapping to
  `beta=lambda/(1+lambda)=0.507042459460253`; the first paper's distinct
  `1.0285722760222` is retained as the old table value. Later reported 2-D
  IPM residuals are
  about `1e-11`--`1e-13`. This remains numerical evidence rather than a
  theorem or a smooth finite-energy basin: the profiles use even density/odd
  vorticity on both horizontal sides and one isotropic scale, not the target
  one-sided right-acute-fan two-scale geometry.
  Added the two Wang audits to the file map. They record the printed curl-sign
  conflict (`curl(U)=-H_X` is authoritative), the absence of official code,
  checkpoint, profile arrays, and SI, and the fact that only the `Omega`
  per-field ansatz is public. They also define a proposed translation-bordered
  odd Newton/pseudo-arclength route and the exact anisotropic stationary no-go:
  unequal rates require a shape-space orbit, not frozen finite-`delta` steady
  states. No odd kernel, continuation branch, or two-scale orbit is claimed,
  and no production solver code was changed.
- 2026-08-28: Reworked the compactified Wang-type analysis prototype into a
  true frozen-base additive multistage correction stack. New corrections have
  a strictly zero output layer, all prior parameters remain unchanged, and
  the composite field is differentiated before evaluating nonlinear residuals.
  The first 387-parameter correction lowers factored combined RMS by `6.73%`
  but raises raw combined RMS by `0.585%`; raw transport is
  `0.157614 -> 0.158534`, while finite-box Green mismatches remain
  `0.7268/0.7096`. Registered the script and audit as an analysis-only
  infrastructure verification/raw-residual negative screen.
- 2026-08-28: Audited the signed-slope bulk continuations with non-cancelling
  residual norms. The smooth odd pure mixture has absolute/reaction/`R_X`
  residual `0.6537/0.6435/0.5521` and outside-fan fraction `0.00294`, so the
  small term-relative metric is rejected as a large-term cancellation.
  Exact-DtN `A`-null `Q2`--`Q4` corrections lower those residuals to
  `0.4783/0.4603/0.4051`, but require coefficient norm `349`, have
  `||Q||/||R||=0.1685`, and raise the outside-fan fraction to `0.1118`.
  A `1e-2` penalty still leaves absolute RMS `0.4863` and outside fraction
  `0.0775`. This crosses the third wall-jet obstruction numerically but remains
  an order-one, fan-delocalized noncandidate. The screen uses a finite periodic
  `X` box and mean-zero Gaussian centers; only the analytic signed-slope base
  and smooth angular tail currently establish all-direction `R_X -> 0`.
  Full-space corrector decay/capping and two-dimensional bulk gluing remain
  open.
- 2026-08-28: Derived the exact continuous-slope identity
  `A R_w=alpha*f+gamma*Hf`. A signed combination of a negative wide-angle
  return and positive narrow-angle head enforces `alpha=0`, `gamma>0`, hence
  the full half-plane Green operator gives the exact wall law
  `f_t-gamma*f*f_x=0`. Smooth generic Burgers breaking selects `beta=1/3`;
  smooth slope endpoints localize the density-gradient activity in a right
  acute sector and give all-direction `R_X=O(r^-2/3)`. Added the exact
  `A`-null normal-jet basis and the wall-jet hierarchy. The stationary lift
  matches the target wall/jet residuals to `2.1e-16/1.75e-10` but retains full
  bulk relative RMS `0.705`, so no two-dimensional profile is claimed. A
  source-neutral wide-box `513x257` physical screen reaches `t=0.331625` with
  wall-slope growth `1.553` versus Burgers `1.496`, acute-fan fraction `0.829`,
  and far velocity/source ratios `0.0039/5.73e-4`; the full gradient is only
  `0.982` of initial and wall-closure defect grows `0.0171 -> 0.3568`. This
  isolates the remaining problem as full bulk gluing, not wall closure.
- 2026-08-28: Derived a third-wall-normal-jet no-go for the pure continuous-
  slope manifold. With moments `M_n=int sigma^n*w`, exact first and second
  jet matching require `M1=1/(2*gamma)` and
  `M2=(1-gamma^2)/(5*gamma^2)`. The two independent cubic-bridge terms at the
  next order reduce to the impossible compatibility
  `-(1+4*gamma^2)/(5*gamma^2)=0`. A constrained analytic-velocity search can
  lower the term-normalized bulk residual to about `1e-2` only with large
  signed variation; this cannot converge to a classical pure-mixture
  solution. The exact Burgers wall branch survives, but `A`-null higher-jet
  or genuinely non-mixture bulk correctors are mandatory. The first missing
  jet is explicit:
  `q3=27*(1+4*gamma^2)*f^2*(f')^5/(70*gamma^3)`. A generalized exact
  `A`-null kernel with any even `lambda(k)>=kappa>0` preserves the selected
  jet and DtN identity while making the zero Fourier mode decay normally.
  The one-sided cubic bridge is not `C1` at its tip; a strict two-scale search
  should use the smooth odd bridge plus an outer signed cap/head, or explicitly
  introduce a third dynamic edge scale.
- 2026-08-28: Strengthened the general-beta hodograph obstruction to the
  explicit affine-to-tail matching law and audited separate edge charts. A
  beta scan found an interior residual valley near `beta=0.70`: the
  `129x161`, 169-parameter state has Poisson relative RMS `1.96e-2`, but its
  coefficient norm is 248 and an unoptimized `161x201` projection regresses
  to `3.02e-2`. The direct two-edge splice has overlap `R_X` defect `1.414`,
  while the small-angle chart drives `min(J)` to `6.38e-9` and Poisson
  residual to `2.67e3`. Both remain explicit noncandidates.
- 2026-08-28: Promoted the exact `0<beta<1` amplitude/length scaling to the
  primary smooth-data blow-up route. Derived the universal Type-I identity
  `rho_x=R_X/(c*(1-beta)*(T-t))`, the exact affine-cusp hodograph system, the
  nonresonant smooth-angular Poisson tail in an arbitrary right acute sector,
  the exact wall DtN formula and conditional CCF closure, the signed-head
  self-strain and moment-derivative constraints, and the dynamic Casimir cap.
  The hard `q_+^beta` fan was rejected because its gradient is singular on the
  edge and fails all-direction decay along fixed-normal-distance paths. A
  finite-angular-mode highest-mode argument also rules out the intended smooth
  decaying closure; the active fixed-profile search must use full hodograph or
  adaptively increasing angular complexity. The numerical CCF exponent is now
  used only as an initialization, never as a selected two-dimensional IPM
  exponent. Added the first general-`beta` hodograph residual prototype. Its
  refined single interior chart preserves transport/compatibility to about
  `2e-4`, keeps `J>0`, matches the prescribed outer density to roundoff, and
  reduces relative Poisson RMS to `6.41e-2`, but a maximum residual near the
  chart/inner-transition boundary remains order one. It is explicitly a
  noncandidate. The exact invariant
  `R*exp(-beta*t/2)/(mu*J)^(beta/2)=constant` shows why a compact-angular tail
  requires degenerate edge charts rather than a single regular chart.
- 2026-08-28: Added an opt-in one-sided advective wall-trace update after
  isolating an `O(h)` extremum drift in the conservative bottom half-cell.
  The default conservative path is unchanged. Added smooth Gaussian-trough
  CCF-wall screen/refinement options, exact initial strain
  `u1_x(0,0)=-4/(3*sqrt(pi))`, a dedicated analysis/plot routine, and a
  stagnation-point transport regression. Terminal states that fall between
  scheduled output times are now appended to history/snapshots, with a
  regression that the saved history reaches `result.physicalTime`. At
  physical time `0.6`, the accepted `257x129`/`513x257` monotone screens have
  zero density-range violation. The refined run gives gradient/wall-slope
  growth `1.2909/1.1799`, origin strain `-0.7513 -> -1.2375`, far-velocity
  ratio `0.0040`, and no remesh. This is early compressive evidence only; the
  short-window singular fit is not accepted as blow-up evidence.
- 2026-08-28: Fixed the far-field target to all-direction
  `rho_x in C_0` and separated three contracts that had previously been mixed:
  compactly supported Lipschitz corner blow-up, an analytic Type-I exact
  trajectory under a prescribed future-singular harmonic strain, and the still-open smooth
  half-plane problem with decaying Biot--Savart normalization. Derived the
  beta-zero wall product identity, compact-active wall requirement, infinite
  wall Taylor hierarchy, and exact hodograph reduction. Continued the
  direction-constant far tail through its third Fredholm resonance and proved
  that an `r^-2 log(r) sin(2 theta)` streamfunction term is forced; an
  exponentially small tail cannot remove the algebraic obstruction. Added a
  falsifiable beta-zero BVP search. Its 129-by-49, 448-mode wide-box run solves
  Poisson to `1.85e-11` but retains bulk/wall transport RMS
  `0.2104/0.1044`, outer-collar maximum `|R_X|=0.2424`, and an oscillatory
  return, so no fixed profile is claimed. Constructed an arbitrary prescribed
  right-acute compact-angular leading tail using two signed smooth bumps. The
  reproducible audit gives Poisson residual `2.22e-16`, Fredholm residual
  `1.50e-18`, same-constant return residual `5.80e-18`, and exactly zero
  activity outside the selected sector. Higher-order transport preserves the
  density-gradient angular support formally, while the elliptic velocity and
  forced log modes propagate through all angles; no global fixed profile is
  inferred from this leading-tail construction.
- 2026-08-28: Replaced the arbitrary local affine strain in the standard
  decaying half-plane contract by its exact Dirichlet Green inverse-moment
  condition `c=P_XY(0,0)`. An explicit compact signed annular density head
  realizes any prescribed positive `c`, proving elliptic self-screening and
  source-neutral sign balance are constructible. A divergence identity rules
  out nontrivial finite-level-area stationary beta-zero profiles, while the
  exact level-area law permits a nonstationary cap at rescaled radius
  `exp(s)` and fixed physical radius. The cap changes inner strain only by
  `O(exp(-s))`. An independent wall Dirichlet-to-Neumann limit gives
  `(X+Hf)f'=0`; its one-interval solution has square-root edge singularities.
  The remaining smooth-Cauchy mechanism is therefore a dynamic signed-head,
  fast-edge, and Casimir-cap gluing problem, not a smooth compact stationary
  wall fixed point.
- 2026-08-28: Classified the sense in which Bessel can be a time-dependent
  limit of smooth data. Constructed an explicit all-horizontal-line,
  pointwise-`C^infinity` characteristic orbit whose dynamically rescaled field
  converges in local `C^1` to the right-decaying directional Bessel family.
  Strengthened this to an exact orbit from bounded `C^infinity` physical
  density compactly supported in `x`, with local `C^infinity` convergence to
  `exp(-X)` and three explicit lengths `L`, `L*g*s`, and `L*(g*s)^2`.
  The compact return layer carries the transported outer reservoir while
  escaping every fixed rescaled compact set.
  Added a stronger Cauchy-contract distinction: finite-time front contraction
  requires a prescribed time-dependent affine far-field protocol, whereas a
  fixed time-independent constant strain produces a unique
  translation-invariant smooth orbit whose double edge zoom converges locally
  in `C^infinity` to the pure `exp(-X)` Bessel mode at infinite physical time.
  Added a wall-compatible sine regularization and a bounded rotating-strip
  finite-time orbit with compact slices. Proved that fixed transverse
  separation, a common vertical cutoff, and finite Fourier closures cannot
  supply the missing finite-energy localization.
  A non-exhaustive primary-literature audit found no exact precedent for the
  compact flat-edge double zoom to `exp(-X)` and no stability theorem for that
  limit. The precise term is directional modified-Helmholtz/Bessel-operator
  tangent, not radial `K_nu` profile or perturbatively proved attractor.
  Proved that its far-left loss of uniform boundedness and its `Y`-linear
  physical growth are unavoidable within the invariant class. Recorded sharp
  no-go results for fixed-coordinate finite-energy attraction, global
  two-Casimir rescaled convergence, and pure K0/K1 source-moment convergence;
  only a local composite Bessel/core plus order-one return limit remains open
  for standard finite-energy smooth data. The finite-time compact-density orbit
  still has infinite transverse mass and an externally selected, time-dependent
  harmonic strain; the fixed-contract long-time orbit removes the time
  dependence but not the infinite-energy affine far field. Neither is a
  standard decaying Biot--Savart Cauchy flow.
- 2026-08-28: Found exact Bessel-bridge singularity mechanisms after retaining
  the order-one closure residual. Added the bounded-density `tanh` front under
  nondecaying harmonic strain, the exact horizontal-wall directional Bessel
  family for general `0<=beta<1`, and a three-scale wall mode. For the
  directional family, reduced the full nonlinear rescaled IPM equation in the
  monotone `Y`-linear class to an explicit characteristic/logistic system and
  proved sector-local `C^1` orbital convergence to the right-decaying
  exponential selected by the initial boundary jet, with the sign of `rho_x`
  preserved. Recorded that
  all current exact wall mechanisms remain infinite-energy or far-field
  unbounded and therefore do not settle the finite-energy sharp-fan Cauchy
  problem.
- 2026-08-28: Discarded the formerly inserted numerical beta: neither reliable
  data nor the current analytic solvability conditions select it. All active
  blow-up scaling, cusp-angle, and matching statements now retain unknown
  `0<beta<1`; the exponent must be selected by the full perturbation/gluing
  problem.
- 2026-08-28: Recast Bessel as a fixed outer background and derived the exact
  one-variable closure-perturbation evolution and reconstruction. Recorded
  inner non-smallness, tip-source cancellation, neutral tangent modes,
  level-set Casimir orthogonality for growing modes, exclusion of nonzero real
  point spectrum for the pure background, the self-similar scaling defect,
  direct affine-plus-Bessel and finite wall-normal-mode no-go results, and the
  exact but geometrically incorrect 90-degree heteroclinic.
- 2026-08-27: Found and audited the exact affine hyperbolic IPM blow-up normal
  form and proved that it is the forced nondegenerate wall Taylor jet. Added
  the no-`x_1`-symmetry leading supported cusp angle
  `theta_*=pi*beta/(1+beta)`, its degree-`2*beta-1` first correction and
  all-orders hard-support caveat, the bounded K1 scale window, and the
  area-preserving Casimir return
  condition. Recorded that affine `beta` is a gauge and that bounded-data
  affine-to-cusp gluing remains open.
- 2026-08-27: Corrected the former Bessel-leading three-layer proposal using
  exact K0 monopole and K1 normal-dipole moment bounds. Added the closure-defect
  Green-kernel formulation, source-neutral return-layer requirement, and the
  conditional general-`beta` multiscale law; marked the older amplitude-only
  matching laws as superseded.
- 2026-08-27: Added smooth near-Bessel K0 external-pole and regularized
  K1-core physical-mode experiments plus read-only tail/`rho_x` diagnostics.
  The accepted 257-by-129, `t=0.05` screens retain their prescribed two-scale
  shapes, but matched 3% perturbations have K0 full/tail difference ratios
  `0.998510/0.998175` and K1 ratios `0.999992/1.00011`; this records short-time
  persistence, not attraction. The pure K1-core datum has zero wall trace and
  a `0.006959` far-exterior relative error; evolution adds an explicit 0.5%
  signed-Gaussian wall monitor solely for the generic wall-peak safety audit.
  Documented the exact anisotropic equations, K0/K1 heat-kernel blow-downs,
  and the strict distinction between far-field scaling and physical-time
  contraction.

- 2026-08-27: Added the decay/regularity obstruction to the exact Bessel
  family contract. A sufficiently decaying smooth physical steady state with
  an impermeable wall is trivial by the Darcy energy identity; the admissible
  `K_nu` modes survive only through a boundary pole. Also recorded the global
  rescaled mass condition `(2 c_l+c_omega) integral Q=0`, which rules out a
  nonnegative finite-mass decaying profile with `c_omega=0` and `c_l~=0`.

- 2026-08-27: Promoted the globally defined right-focused Bessel solution to
  the polar research contract. Recorded the positive punctured-plane `K_0`
  member, its `r^(-1/2) exp[-mu*r*(1-cos(theta))/2]` tail and
  `(mu*r)^(-1/2)` angular width, the homogeneous-wall `K_1*sin(theta)` member,
  and the smooth interior through-flow obtained by placing the pole below the
  wall. Proved that these profiles are dynamic-scaling fixed points only
  at zero rates; arbitrary scaling produces the parameter orbit
  `mu_rescaled=mu/C_l`. Explicitly separated steady existence from dynamical
  selection out of smooth bounded data.

- 2026-08-27: Added `ipm_plot_exact_bessel_family` and corrected the default
  visualization to the requested right-facing acute wedge. The default
  `theta_*=60 degrees`, `lambda=-1` gallery compares the nonhomogeneous
  wall-fed cosine mode `nu=pi/(2 theta_*)` with the two-edge-zero sine mode
  `nu=pi/theta_*`. It plots density and analytic `rho_x` with logarithmic
  dynamic-range compression, plus detailed signed `rho_x` and radial cuts.
  The maximum transformed-Helmholtz Poisson residual is `4.49e-15` relative
  and physical transport cancels exactly. The same function now accepts
  `theta_*=pi` for a direct upper-half-plane comparison and reports edge
  normal velocity, separating the wall-fed mode from the impermeable sine
  mode. The upper-half-plane Poisson residual is `2.11e-15`; relative edge
  flux is `1.00` for the wall-fed mode and `1.09e-16` for the sine mode.
  Documented that zero extension
  outside the wedge requires a distributional interface source unless the
  normal derivative is matched.

- 2026-08-27: Recorded the exact far-decaying physical steady Bessel family
  `rho=lambda*Psi`, obtained by conjugating
  `Delta Psi+lambda*partial_X Psi=0` to modified Helmholtz. The documented
  family permits nonhomogeneous wall traces, arbitrary real angular order,
  center singularity, and same-`lambda` modal superposition. Its fixed-angle
  far field has `alpha=1`, `beta=-1/2`, and
  `g(theta)=-(|lambda|+lambda*cos(theta))/2`; integer sine modes retain the
  optional homogeneous impermeable wall. Explicitly separated this physical
  steady family from constant-rate rescaled states with `c_l~=0`.

- 2026-08-27: Added a MATLAB local polar mesh for the exact degree-one fan
  steady state about `(1,0)`. Upper/lower density and streamfunction branches
  are meshed separately so the wall density jump remains visible; the third
  panel shows the total rescaled transport `U+c_l z`. For the default
  `theta_*=60 degrees`, `A_0=1`, the chosen gauge is
  `c_l=0.75,c_omega=0`; analytic steady, physical-transport, fan-normal,
  wall-normal, Poisson, and divergence residuals are at or below `5e-16`.
  Recorded explicitly that only `c_l-c_omega=a_0` is invariant and that
  `c_r=0` assumes centered coordinates.

- 2026-08-26: Added a read-only polar-profile research layer centered at a
  wall point displayed as `(1,0)`. Derived the centered steady IPM equations
  and separate symmetric-quadrant/asymmetric-half-plane boundary contracts.
  The local cusp balance gives a reduced nonlinear angular problem and a
  fixed-angle compatibility obstruction for a nonzero wall power law,
  motivating the measured shrinking-fan ansatz `theta_f(r)~r^sigma`. Added
  wall/polar exponent agreement, angular collapse, fan-power diagnostics, and
  a synthetic arbitrary-exponent regression. The saved physical and dynamic
  endpoints both fail the new single-scale consistency decision and remain
  transient/nonconverged data.

- 2026-08-26: Added selectable `weno5_nonuniform` transport and made it the
  maintained no-argument choice. The implementation follows the arbitrary
  nonuniform point-value algorithm of Martí--Mulet--Yáñez--Zorío (2024):
  three quadratic candidates, one global quartic, simple slope indicators,
  and grid geometry precomputed after every remesh. A rejected cell-average
  trial violated the program's nodal data contract and generated a 2.4% wall
  derivative TV increase; the retained nodal version completes the matched
  real IPM regression without secondary peaks. Added `ipmtests.baseline.transport`
  and recorded fifth-order face reconstruction, second-order full transport,
  conservation, sharp-step stability, and WENO/MUSCL error comparisons.

- 2026-08-26: Fixed the no-argument entry contract. `ipm.solve()` previously
  completed an empty structure through the generic defaults and therefore ran
  the degenerate physical case even though the maintained case was CCF
  heavy-tail with partial dynamic scaling. It now reads the stateless
  `ipm.config.activeCase` overrides; explicit option structures retain generic
  behavior. `long_dynamic` reuses the same override source to prevent drift.

- 2026-08-26: Added the legacy-inspired `lattice` adaptive-axis profile for
  CCF long runs. Flat fine-grid plateaus and narrow tanh transitions replace
  the broad overlapping monitor tails, while the paired symmetry and existing
  arbitrary-grid Poisson/difference path are unchanged. Long-run target caps
  are now `128` and `[192,128,96]` with expansion cap `8` and buffered safety
  `0.08`. Clarified and tested
  that `remeshMaximumCellRatio` bounds adjacent cells only; global spacing
  ratios are recorded separately and never reject or truncate a run.
  A strict-cap regression also fixed the case where the first
  count-satisfying trial violated the adjacent ratio: the search now retains
  the strongest legal sub-target deformation. The then-maintained CCF long cap is
  `1.5`, with repeated analytic pre-remeshing and what was then a stable
  partial-dynamic comparison (`symmetry_peak`, length gain `0.5`,
  gradient-energy amplitude gauge). The 2026-08-29 fixed-anchor entry
  supersedes that active gauge. Dynamic mode retains wall-normal diagnostics
  and remeshing.

- 2026-08-24: Extended the high-resolution CCF physical mesh audit with
  on-demand transition smoothing. The default/original-reference path is
  unchanged. In the opt-in anchored current-reference path, the two-node halo
  and exterior log-width smoother activate only after 90% of the hard
  adjacent-cell-ratio budget is spent. The accepted `1025 x 513`,
  `[-128,128] x [0,64]` run follows the previous remesh schedule through
  `t=0.3460`, then lowers the saturated ratio `8 -> 6.20` and reaches
  `t=0.392521` in 4174 steps with 26 remeshes. It stops at the horizontal
  resolution floor (`4.98` 90%-core cells); the strict all-checks trusted
  window ends at `t=0.368068` with `10.31` core cells. Terminal TV is `1`,
  far-source ratio is `1.07e-3`, range violation is `1.51e-3`, and physical
  PCHIP remap/transport mass drift is `-6.46e-3`. The three trusted-window
  fits (wall peak, inverse FWHM, inverse core width) all return
  `candidate=false`; the run is a resolved concentration/mesh-limit record,
  not a blow-up claim. The compact history is
  `result/verification/ccf_tail_long_physical_history.csv`.

- 2026-08-24: Tested and rejected post-PCHIP physical mass corrections.
  Global support corrections contaminated the far source; separate row/column
  corrections preserved quadrature mass but produced order-`10^2` artificial
  cross-direction gradients at a clean endpoint. Physical remaps therefore
  retain range-preserving PCHIP, exact wall-trace preservation, peak-change
  and gradient/range diagnostics, while exposing mass drift. The verification
  fixture now also bounds physical-remap gradient amplification. Dynamic mode
  retains its established bounded horizontal row-integral correction.

- 2026-08-24: Added opt-in far-grid anchoring and CCF long-time mesh cases.
  Their peak and 10/50/90-percent targets are `48` and `[96,48,24]`, safety
  floors are `16` and `[24,14,10]`, and the outer anchor factor is `3`.
  Existing production defaults retain their previous targets and global
  redistribution.
  The ratio-6 long run reached `t=0.2901` before its core fell to `8.08`
  cells, safety rose to `1.238`, both cell-ratio caps saturated, and TV reached
  `1.027`. Far-source and mass diagnostics remained small. The CCF-only mesh
  therefore uses a longer anchor transition and adjacent-ratio cap `8`; the
  production cap remains `6`.

- 2026-08-24: Added the regularized `ccf_heavy_tail` initial-data family whose
  wall derivative converges to the proposed
  `1_{x_1>1}(x_1-1)^(-1/2)` CCF profile before a monotone saturation scale.
  Added explicit onset, amplitude, cutoff, vertical-envelope, doubled-box,
  halved-regularization, and optional analytic pre-remesh contracts plus a
  second-order `D_x` regression away from the deliberately nonsmooth onset.
  A `1025 x 513`, doubled-box physical verification reaches `t=0.25` with
  wall-peak growth `2.065`, FWHM ratio `0.171`, core-width ratio `0.0836`,
  mass drift `3.72e-6`, density-range error `7.95e-4`, and no oscillation.
  The far-source ratio rises to `0.0134`, and power fits do not separate from
  exponential growth or share a stable terminal time; this is recorded as
  strong CCF-like concentration, not a blow-up claim.

- 2026-08-22: Replaced isolated quartic peak monitors by a target-driven
  nested 10/50/90-percent monitor with density-per-width weights and a graded
  origin-to-peak bridge. Added buffered per-axis triggering, simultaneous
  two-axis updates, minimum-spacing and adjacent-cell-ratio constraints, and
  the stateless `ipm.remesh.axis` responsibility. The accepted 1025-by-513
  physical run reaches `t=1.44` with peak `106.329`, gradient `143.169`,
  horizontal/vertical 90-percent core counts `5.02/4.47`, density-range error
  `1.04e-7`, one support maximum, total-variation ratio `1`, and cell ratios
  `5.95/5.97`. At `t=1.4` it differs from the previous 1025 baseline by
  `0.36%` in peak and `2.89%` in FWHM. The 513 time-step halving changes
  peak/FWHM by at most `0.67%/1.57%` through `t=1.428`; 513-to-1025 spatial
  differences grow to `12.0%/26.8%` there, so the strongest endpoint is
  resolved and non-oscillatory but is not claimed spatially converged.
- 2026-08-22: Added wall-profile oscillation and grid-smoothness contracts to
  history, stopping, and the trusted fit mask. Classified legacy source,
  reproducible cases, accepted outputs, verification data, and rejected
  experiments into dedicated directories without deleting history.

- 2026-08-22: Restricted live/video and final field observation to the
  independent `x_1>=0` side of the active double-odd solution. Documented that
  the numerical state still stores 512 mirrored negative nodes, one origin,
  and 512 positive nodes at `nx=1025`; converting the solve itself to a
  first-quadrant grid remains a separate Poisson/transport boundary change.

- 2026-08-22: Restored the symmetric physical equation as the active default;
  physical initialization records peak/core references for mesh tracking but
  all three modulation rates, both scale factors, and `X_shift` remain exact
  identity values. Remesh requests are checked after every SSPRK step and
  blow-up fits use only states with `safetyFactor<=1`.
- 2026-08-22: Extended non-compounding adaptive redistribution to both axes.
  The paired `x_1` monitor follows the wall peak and the one-sided `x_2`
  monitor follows its wall-normal profile; both core counts enter the stop
  criterion. Axis-wise updates avoid cross-axis overexpansion, local monitor
  quadrature resolves sub-base-grid windows, wall-normal transfer preserves
  the discontinuous upper wall trace, and physical zero-`c` regressions were
  added.
- 2026-08-22: Added physical-remap peak consistency: a candidate changing the
  wall `rho_x1` maximum by more than 3% is rejected. A proactive safety trigger
  of `0.7` was experimentally rejected after 34 remeshes failed to extend the
  trusted endpoint and degraded the terminal fit for the former isolated-peak
  monitor. That historical result does not apply to the replacement nested
  monitor, whose accepted trigger is `0.75`.
  The physical default `maxDt/outputEvery` are now `2.5e-4/0.002` based on the
  saved time-step study.
- 2026-08-22: Blow-up diagnostics now report 20/40/80/160-point window fits in
  addition to the primary 40-point result. The time-refined `1025 x 513`,
  box-400 physical run reaches trusted `t=1.41001039` with wall peak `37.8769`
  and horizontal/wall-normal core counts `6.07/7.99`, then stops before using
  an under-resolved terminal state. The 40/80-point full-gradient fits give
  `T about 1.466`, `p about 1.32`, while 20/160-point fits drift materially;
  this is recorded as a multiscale growth candidate, not a converged exponent.
- 2026-08-22: Updated live/MP4 output for physical mode. It follows the
  wall-normal width, updates both nonuniform axes after remeshing, and displays
  full-gradient growth, wall peak, two-axis safety, minimum spacings, and
  density-range error instead of three identically zero modulation curves.

- 2026-08-18: Established the physical half-plane IPM path, reusable sparse
  operators, conservative transport, verification, diagnostics, and this
  maintained architecture record.
- 2026-08-18: Fixed the signed rotation contract and added manufactured tests
  for `D_x`, `D_y`, Poisson, and the complete velocity operator.
- 2026-08-18: Made the requested one-sided `k=4` density the default; fixed
  its exact `x>0` support and nonzero wall trace.
- 2026-08-18: Added optional IPM-specific dynamic rescaling with both `c_l`
  and `c_omega`, physical-time reconstruction, semi-discrete gauge checks,
  and physical-gradient blow-up diagnostics.
- 2026-08-18: Switched the active default to dynamic scaling to prevent rapid
  extremum contraction from outrunning the fixed grid. The no-`c` physical
  mode remains available explicitly and is still tested independently.
- 2026-08-18: Added configurable length over-zoom and set
  `lengthScaleGain=2` by default because center pinning alone did not spread
  the neighborhood of the maximum sufficiently.
- 2026-08-18: Added horizontal translation rate `c_r` and integrated
  `X_shift` so over-zoom can spread the peak neighborhood while the traveling
  wave remains locked at `X_pin`.
- 2026-08-18: Added live MP4 output, adaptive peak-width feedback for `c_l`,
  a far-boundary contamination ratio, and a nonuniform second-order stretched
  grid.
- 2026-08-18: Extracted the legacy nonuniform-grid contracts without importing
  its stateful Boussinesq call graph: reference derivatives divided by the map
  Jacobian, matching weighted Poisson stiffness, image-Green artificial
  boundary data, and `10%/50%/90%` support criteria. The initial
  `[-500,500] x [0,500]`, `513 x 513` box retains the old mesh's extended
  far-field design and uses automatic `0.05` center spacing.
- 2026-08-18: Replaced discontinuous level-cell counts by linearly
  interpolated level widths and raised the outer/middle/core targets to
  `64/32/16`; the default initial core therefore activates stronger spreading
  instead of waiting for visible under-resolution.
- 2026-08-18: Replaced the fragile point-strain `c_omega` default by a
  positive full-gradient-energy gauge. Retained point strain and `c_omega=0`
  as explicit comparison modes and added multiscale gradient-energy and gauge-
  conditioning diagnostics.
- 2026-08-18: Added the bounded rational default
  `rho_0=chi(x>0)x^k/(1+x^k+y^k)` and retained the earlier logarithmic
  primitive as an explicit comparison. Reworked adaptive length feedback to
  count the main wall-`omega` peak's connected component using actual
  nonuniform cells, raised the maximum gain from eight to sixteen, and added
  direct derivative and separated-secondary-peak verification.
- 2026-08-19: Replaced linear adaptive-gain amplification by a configurable
  square-root response. This retains connected-peak resolution feedback while
  reducing excessive `c_l X` CFL cost and physical-clock slowdown.
- 2026-08-19: Corrected live and final visualization to display reconstructed
  physical `rho_x=(C_l/C_omega)R_X` instead of silently labeling the decreasing
  rescaled derivative as physical. Added explicit result fields, wall maxima,
  and a regression check showing positive initial physical peak growth.
- 2026-08-19: Added the wall-omega-peak amplitude gauge and made it default so
  the positive rescaled `R_X` maximum has zero semi-discrete rate. Retained the
  full-gradient-energy and point-strain gauges as independently verified
  multiscale and legacy comparisons.
- 2026-08-19: Replaced one-way `c_l` over-zoom by bidirectional connected-width
  feedback. Adaptive targets are now bounded multiples of the initial primary
  peak, the baseline gain is one, and gains below one retract over-expanded
  level sets while `c_omega` and `c_r` continue locking amplitude and position.
- 2026-08-19: Added a common dynamic-time normalization with default maximum
  evolved rate two. It multiplies the complete vector field and clock, so
  growing canonical `c_l` and `c_omega` are recorded but no longer make the
  evolved rates diverge or alter the gauge trajectory through separate caps.
- 2026-08-20: Replaced the unstable multiplicative width gain by an additive,
  outer-priority multiscale controller. Its positive response fades before the
  target cap, its negative response overrides an under-resolved inner core,
  and the tracked peak is restricted to the initial locked-wave neighborhood.
- 2026-08-20: Added an explicit canonical rescaled clock so the dynamic-rate
  cap is only a common time reparameterization. Added physical velocity,
  physical mass, both rescaled clocks, and their snapshot/history contracts.
- 2026-08-20: Recorded long-step, grid, rate-cap, and rational-datum domain
  stress tests. Small boxes are quantitatively rejected for the uncut
  direction-dependent far-field datum.
- 2026-08-20: Replaced destructive top-magnitude Green-source truncation by
  signed block aggregation that preserves total strength and centroid. Added
  direct boundary and core-velocity compression regressions.
- 2026-08-20: Corrected `c_r` from a zero-velocity condition at a fixed grid
  node to a subgrid tracked-peak restoration gauge. Long tests distinguish its
  actual peak offset/residual from the length-gauge residual.
- 2026-08-20: Changed physical field, wall, and final plots from misleading
  fixed `X` axes to peak-aligned physical coordinates while retaining one
  explicit rescaled-control panel.
- 2026-08-20: Audited far-velocity warnings for the angular datum. The raw
  ratio and tracked frame speed are recorded, but only source tails and
  peak-aligned box expansion are treated as actionable convergence tests.
- 2026-08-20: Added a stateless result comparator for peak-aligned density and
  physical-gradient convergence across grids or boxes.
- 2026-08-20: Replaced the box-dependent point-velocity length gauge by the
  tracked local-strain gauge. Constant horizontal velocity now changes only
  `c_r`, eliminating artificial coupled growth of `c_l` and `c_omega`.
- 2026-08-20: Added reconstructed density-range diagnostics and a configurable
  maximum-principle stop, so amplitude overshoot cannot be mistaken for IPM
  blow-up.
- 2026-08-20: Decoupled MP4 frames from full-field snapshot storage. Snapshots
  are opt-in, preventing default `513^2` long runs from retaining roughly one
  gigabyte of redundant arrays.
- 2026-08-20: Replaced separately centered dynamic-scaling derivatives by one
  full-velocity nonuniform MUSCL flux with open artificial density boundaries
  and the conservative `2*c_l R` correction. Added constant-dilation and exact
  canonical-flux regressions; long-run maximum-principle error fell from order
  `1e-2` to order `1e-5` or less.
- 2026-08-20: Made width and traveling-wave restoration feedback amplitude-
  gauge covariant by multiplying their dimensionless errors by the local
  nominal strain rate. Added translated-state covariance and real MATLAB
  peak/energy/none gauge comparisons.
- 2026-08-20: Added `physicalFinalTime`, automatic stretch provenance/reuse,
  and common-physical-time comparisons. Raised the conservative default to
  `1025 x 513` on `[-2000,2000] x [0,2000]` after matched grid/domain studies.
  The derivative field is near percent-level converged at physical time two,
  while the point maximum and half-height width remain explicitly untrusted as
  standalone blow-up evidence.
- 2026-08-20: Added a hard connected-peak safety floor and bounded local
  `x`-remeshing fallback. Redistribution is referenced to the original
  stretched grid, transfers density by PCHIP without new extrema, enforces a
  row-wise conservative bounded correction and center-spacing/CFL floor,
  rejects changes below 2%, rebuilds the arbitrary-grid Poisson/derivative
  operators, and records remesh/spacing diagnostics.
- 2026-08-20: Corrected physical/canonical terminal-time classification when
  the floating-point final remainder is smaller than `minDt`.
- 2026-08-20: Added an amplitude- and frame-invariant wall-density-width
  length gauge as an explicit comparison. Stress tests rejected it as the
  production default because a broad material-density width did not control
  the narrower derivative core; local strain remains active.
- 2026-08-20: Made double-odd `omega` symmetry active. The first-quadrant
  density trace is extended evenly in `x_1` and discontinuously/oddly in
  `x_2` as explicitly requested. Green data now use the four-image quadrant
  kernel and integrate only first-quadrant sources. The induced `rho`,
  `omega`, `psi`, and `u` parities, paired remeshing, zero origin velocity, `c_r=0`, and
  `X_shift=0` are verified; the older half-plane kernel remains an explicit
  comparison mode.
- 2026-08-20: Replaced the active double-odd local-strain length gauge by the
  symmetry-anchored origin-to-positive-peak secant gauge. A failed run showed
  that the old gauge kept the rescaled peak amplitude within about 4% but let
  its half-height width contract from `0.807` to `0.151`, exhausted six
  remeshes, and stopped with only about two 90%-core cells. The new gauge uses
  the already verified `U_1(0,0)=0` contract to remove this length drift while
  keeping `c_r=X_shift=0`. The live/video first panel now displays rescaled
  `R_X1` directly beside physical `rho_x1`.
- 2026-08-21: Replaced the passive peak-velocity/secant gauge by the active
  `R_X1` extremum-location gauge obtained from the evolution of
  `R_X1X1(X_*)=0`. On a `513 x 257` regression through physical time `1.45`,
  the peak position stayed near its initial value, the rescaled peak stayed
  near `1.06`, canonical rates remained bounded, and `timeSpeed` stayed one.
  The simultaneous FWHM/core contraction is recorded as a multiscale
  diagnostic rather than hidden by forcing one affine scale to lock both.
- 2026-08-21: Decoupled active `c_l` from grid-cell feedback. Core safety now
  triggers a paired monitor centered on the derivative peaks and sized by the
  connected 90% width. A remesh is rolled back if recomputing `R_X1` makes
  the core resolution worse. Added `trackedWallCoreWidth` to
  flow/history/video diagnostics and explicit accepted/rejected transfer
  regressions.
