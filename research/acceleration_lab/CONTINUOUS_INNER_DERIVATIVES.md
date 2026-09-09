# Consistent continuous inner geometry and its actual time derivative

These are independent, read-only diagnostics of the original PDE. They do
not change `c_l`, `c_omega`, the field, the accepted trajectory, or the
decision on any earlier secant candidate. In particular the actual
continuous peak derivative is measured even when the production quadratic
peak derivative vanishes.

## One interpolation operator throughout

Write `Omega = rho*Dx'` and `F = rhs*Dx'` on one fixed native grid. The wall
operator `H1` is the C1 cubic Hermite interpolant with maintained `Dx`
slopes. The tensor operator `H2` uses the four nodal jets

\[
\Omega,\quad \Omega D_x^T,\quad D_y\Omega,
\quad D_y(\Omega D_x^T).
\]

At fixed query coordinates these operators are linear in the field data,
so their true time derivatives are `H1[F]` and `H2[F]`. Their spatial
derivatives are derivatives of those same polynomials. Differentiating a
different interpolant or interpolating a nonlinear PCHIP slope rule would
not justify this identity. Once peak and widths are determined from the
field, the *complete* normalized coordinate map is nonlinear in the field.

Let `a` be the unique positive, nondegenerate interior maximum of the wall
interpolant on the registered window, and let `P=H1[Omega](a)`. On a fixed
peak-cell branch,

\[
P'=H_1[F](a),\qquad
a'=-\frac{\partial_xH_1[F](a)}{\partial_{xx}H_1[\Omega](a)}.
\]

For level `q=.9`, choose the nearest connected left and right roots
`ell<a<r` of `H1[Omega]=qP`. Set `wx=r-ell`. Their rates are

\[
\ell'=\frac{qP'-H_1[F](\ell)}{\partial_xH_1[\Omega](\ell)},\quad
r'=\frac{qP'-H_1[F](r)}{\partial_xH_1[\Omega](r)},\quad
w_x'=r'-\ell',\quad \beta=w_x'/w_x.
\]

The vertical width `wy` is the first positive downcrossing of
`H2[Omega](a,y)=qP`, using exactly the same `a`. Its derivative is

\[
w_y'=\frac{qP'-H_2[F](a,w_y)-a'\partial_xH_2[\Omega](a,w_y)}
 {\partial_yH_2[\Omega](a,w_y)},\qquad \gamma=w_y'/w_y.
\]

The translation term in this equation is needed even if the wall maximum
is stationary in `x`: away from the wall the trace generally has nonzero
horizontal derivative. It was explicitly exercised in the coupled
manufactured test, not assumed negligible.

For `X=a+wx*xi`, `Y=wy*eta` define

\[
U(\xi,\eta)=\frac{H_2[\Omega](X,Y)}{P},\qquad
G=U_\tau=\frac{H_2[F]+(a'+w_x'\xi)H_{2,x}[\Omega]
 +w_y'\eta H_{2,y}[\Omega]}{P}-\frac{P'}P U.
\]

All terms are retained. The rates are determined by peak and level-set
identities, with no least-squares fit. In the isotropic production
coordinates, positive physical inner compression rates per canonical time
are `c_l-beta` and `c_l-gamma`, since physical widths are `wx/C_l` and
`wy/C_l`.

## Events and conditioning

The peak helper retains the local cell scale, normalized curvature,
stationarity error, active maxima, cell fractions and local coordinate ULP
tolerances. Multiple or nearly tied maxima, boundary/flat peaks, absent or
nontransverse connected roots, large root residuals, and threshold-flat
cells intersecting the chosen component are rejected explicitly.

A C1 field has continuous first derivatives at a native cell interface,
so a transverse level crossing may retain a first derivative there. Peak
motion involves a second derivative, which may jump. If the peak lies on
or within its recorded local tolerance of an interface, the helper saves
both one-sided curvatures and phase rates and withholds a single phase rate.
It does not snap the peak to a favorable cell or infer C2 smoothness from
C1 interpolation. Genuine multi-peak exchanges require a separate
directional/event analysis. A remesh also changes the operator and cannot
be silently treated as a differentiable time step of this fixed-grid map.

These qualifications preclude a generic fourth-order claim for peak
positions or coupled widths. Even for smooth functions, a cubic-Hermite
derivative-root position typically has an `h^3` error; coupling to the
vertical trace can transmit that position error at first order. Neither a
single-phase observed ratio nor exact polynomial reproduction certifies
uniform order across all phases or events.

## Executed tiny validation

`ipm_accellab_test_continuous_inner` used a nonuniform 49 by 41 polynomial
field

\[
\Omega=A[(1-q_x^2+\delta q_x^3)(1-q_y^2)+\kappa q_xq_y^2],\quad
q_x=(x-a)/L_x,\quad q_y=y/L_y,
\]

with `a(0)=.197`, `a'=.07`, `Lx(0)=.55`, `Ly(0)=.4`,
`(log Lx)'=-.18`, `(log Ly)'=-.23`, `(log A)'=.11`,
`delta(0)=.035`, `delta'=.6`, and `kappa=.35`. The evolving asymmetry creates
a nonzero residual shape derivative after normalization. The largest
analytic geometry-rate error was `4.996e-16`; the analytic `G` error was
`2.665e-15`. Adding `.37*Omega` to the forcing shifted the measured peak
rate and left the geometric rates and normalized `G` invariant to
`1.999e-15`. Reciprocal horizontal/vertical coordinate unit changes by
`1e3` and `1e-3` agreed to `4.885e-15`.

For time increments `1e-3, 5e-4, 2.5e-4, 1.25e-4`, central analytic-path
differences of the entire recomputed pullback gave errors
`8.4259e-8, 2.1069e-8, 5.2804e-9, 1.3124e-9`, with all relevant geometry
templates fixed. This establishes the expected second-order time-difference
check of the first-derivative formula, not the time integrator's fourth
order.

The original 65 by 33 PDE and its original SSPRK54 stages supplied a
separate real-forcing check. At each increment the test took two forward
steps and evaluated `(-3U0+4Uh-U2h)/(2h)`. Its errors were
`6.15785e-6, 1.54009e-6, 3.85058e-7, 9.63398e-8`, again with fixed geometry
templates. The initial actual `P'=.00608866833856` was retained. This cost
one tiny operator build, eight tiny PDE steps and 41 tiny RHS evaluations;
there was no large-grid solve. The report is
`continuous_inner_tiny_tpcd2a4343_0e61_4a3f_8df8_b57c64465b52/` beneath the
campaign acceleration output directory.

The event tests withheld the phase rate at an exact native peak interface
while retaining the correct one-sided `.07` translation rates. A separate
smooth polynomial with horizontal and vertical roots on native interfaces
retained the correct implicit rate (`G` error `1.111e-14`). This smooth
special case is not a general C2 event result. Two additional threshold
plateau negative tests passed in
`continuous_flat_events_tp844d19c0_ba9e_4f44_8b2c_4b67e5459555/`. Those
adversarial tests use explicitly labeled fixed synthetic linear jet maps,
not PDE derivative matrices; they test event rejection only.

The earlier manufactured choice `a(0)=.173`, `delta'=.025` is preserved as
a failed time-difference experiment: its left root crossed a cell on three
increments, and errors of order `1e-11` were dominated by roundoff. Its
complete reproduced measurements remain in
`continuous_inner_tiny_tp0f73864d_df82_42b1_9a92_0b79455829fa/`; the initial
failure is also retained. The stronger interior-branch manufactured signal
was registered separately. No derivative formula or pass tolerance was
changed to turn that failed experiment into a pass.

## Actual original native cache at tau 4.78320194454

Only the saved accepted native baseline `fields{1}` at step 2270 was used
from `fixed_observers_tpf210ce8d_3aee_4a8c_93dd_c1add1567b86/`. Its actual
original `F_X`, paired axes and maintained derivatives came from the
previous strict five-checkpoint audit and fresh original RHS comparison.
The second, rejected secant field was not evaluated here. No new checkpoint
restore, Poisson solve, RHS evaluation or PDE step occurred. This is an old
native state, distinct from the later fresh finite-box case.

| Quantity | Continuous C1 value |
|---|---:|
| `P` | .686253805674884 |
| `P'` | 3.03354868260e-5 |
| `P'/P` | 4.42044715455e-5 |
| `a` | .971970246205102 |
| `a'` | -.00172747461049 |
| `wx`, `wy` | .0227212617972, .00679057269250 |
| `beta`, `gamma` | -.654182059295, -.728497126938 |
| `c_l-beta`, `c_l-gamma` | .820892929076, .895207996719 |

The vertical translation forcing was `-.0106469082926`, alongside the
fixed-X forcing `-.0419337359062`; omitting it would materially change the
vertical rate. The selected peak-cell fraction was `.66508`, normalized
curvature `-.00292722`, and its nearest-node distance `.000435115`, versus
a local interface tolerance `2.808e-14`. Root slopes were well resolved,
all root residuals at most `1.433e-15`, and both threshold-flat lists empty.
The recorded Poisson residual `4.79415e-10` and original near-singular
warning (`RCOND=8.911277e-17`) are preserved, not reassessed by this no-LU
observation.

On the fixed inner rectangle `[-2,2] x [0,3]`, the core is
`[-1,1] x [0,1.5]`; the rest is the holdout. Every native mapped cell and
rectangle boundary is split before four-point tensor Gauss quadrature.
`G` is at most bicubic per cell, so this integrates `G^2` to roundoff for
the selected interpolant. These are RMS L2 values; reported infinity
quantities remain sampled maxima rather than certified suprema.

| Region | Normalized Eulerian RMS | Complete `G` RMS | Remaining fraction |
|---|---:|---:|---:|
| Core | .164910449185 | .007075535373 | .0429053 |
| Holdout | .178121681177 | .018656742600 | .104742 |
| Full inner rectangle | .174912446789 | .016539991944 | .0945615 |
| Wall core | .162090681070 | .006200609162 | .0382540 |
| Wall holdout | .214219401683 | .022068959487 | .103020 |

Thus the exact coordinate motion accounts for about 95.7% of the core RMS
evolution and 89.5% of the holdout RMS evolution at this state. The complete
`G` divided by profile RMS is `.00862045` in the core and `.0301233` outside
it per canonical time. These ratios describe the present derivative; they
do not establish its decay, a fixed point, a slow eigenmode, a uniform
error bound, a Holder exponent or a singularity. The holdout is a local
surrounding rectangle, not the far-field matching or velocity gate needed
for a candidate trajectory.

The original 65 by 49 and 129 by 97 node-masked trapezoids are retained.
Their core weights included full boundary-node weights from the whole
rectangle, so their effective core areas exceeded 3. Their differences
from split-cell Gauss must not all be attributed to interpolation error.
`ipm_accellab_continuous_sampling_from_cache` separately clips trapezoid
weights to the same core and holdout boundaries, saving both original and
corrected measurements without recomputing any field or geometry.
For the old state, boundary-clipped core RMS values are `.007093946390`
and `.007080387532`, respectively 0.2602% and 0.06858% above the Gauss
reference. The uncorrected values were `.007518338950` and
`.007299151597`, with effective areas `3.158203125` and `3.07861328125`.
Thus most of that particular apparent core discrepancy came from boundary
weighting. The separate audit is in
`continuous_sampling_tp0a734f71_5d32_4912_b979_432a1cd77979/`. This result
does not retroactively change any old accelerator metric or rejection.

The full actual report and small observation fields are in
`continuous_inner_cached_tp340defbd_db76_46a9_a80c_3946318e1510/`. Executed
sources were frozen there before the later threshold-flat guard addition.
That guard cannot affect this frame because both flat-cell lists are empty.

## Entry points and next comparison

`ipm_accellab_continuous_inner_rates(Omega,F,x,y,Dx,Dy,window,options)`
returns the geometry, genuine derivatives, conditions and event signatures.
`ipm_accellab_continuous_residual(...)` evaluates the no-fit decomposition,
same-cell Gauss norms and both sampling densities.
`ipm_accellab_probe_continuous_saved` is specifically registered for the
audited old step-2270 cache. `ipm_accellab_probe_continuous_fresh_cache`
uses an explicitly paired fresh-case RHS cache, lineage, actual anchor
window and clock/unit conversion; it does not restore a native checkpoint.

For a fresh physical restart with parent constants `Cx0`, `Comega0`, put
`lambda0=Cx0/Comega0`. The coordinate correspondence is

\[
X_{parent}=C_{x0}X_{fresh},\quad
\Omega_{parent}=\Omega_{fresh}/\lambda_0,\quad
\Delta\tau_{parent}=\lambda_0\tau_{fresh}.
\]

Consequently widths and positions multiply by `Cx0`; `beta`, `gamma`,
`c_l`, `c_omega` and `U_tau` divide by `lambda0`; `a'` multiplies by
`Cx0/lambda0`; `P'` divides by `lambda0^2`. Absolute physical time adds the
recorded parent epoch. These unit relations do not prove two discretized
finite-box cases have identical trajectories. Actual source, window,
lineage, spatial differences and matching-time checks remain explicit.

## Executed later fresh-cache observation

The separately qualified 895 by 386 fresh finite-box case supplied a
complete `true_rhs_cache.mat` at its accepted step 139, local canonical
time corresponding to parent-equivalent `tau=7.24320075185916` and absolute
physical time `1.95904009840289`. The source continuation's trusted-history,
native pairing, exact prefix, unchanged axes, lineage/reference, matching
absolute-time and resolution gates passed. The parent process saved this
cache after its one explicitly authorized restore; this observer performed
zero further restores, operator builds, RHS evaluations or time steps.

Its `physicalAnchorX=.152142234998802` fixes the registered current
computational window `[0,.304284469997604]`. The selected peak agrees
bitwise with the global positive-half maximum. The lineage gives
`Cx0=6.57279683059654` and `lambda0=32.2053866612671`; all comparisons below
are in parent units. The source Poisson residual `6.67177e-9` is retained.

| Quantity | Later fresh observation, parent units |
|---|---:|
| `P`, `P'` | .686076040790558, -2.93723915462e-5 |
| `P'/P` | -4.28121517148e-5 |
| `a`, `a'` | .974363586787440, .00154354675378 |
| `wx`, `wy` | .00466832105026, .00127810721300 |
| `beta`, `gamma` | -.625235744978, -.662096563652 |
| Native `c_l`, `c_omega` | .0696779224956, -.350846929268 |
| `c_l-beta`, `c_l-gamma` | .694913667474, .731774486147 |
| Core, holdout `G` RMS | .007786883124, .014645378150 |
| Core, holdout `G / normalized Eulerian RHS` | .0311343, .0629433 |
| Core, holdout `G / U` | .00939542, .0232596 |

In fresh units the vertical fixed-X forcing was `-93.3884674` and its
translation correction `+44.6043084`; their sum, and actual `PPrime`,
were used. The peak-cell fraction `.57242`, normalized curvature
`-.00091752`, normalized root slopes (minimum magnitude `.00564155`),
and root residuals at most `8.798e-15` resolve the local branch. There were
no threshold-flat cells; nearest-node distance `9.272e-6` greatly exceeded
the local interface tolerance `4.325e-15`.

The later residual fractions are smaller, but the absolute core residual
is `.00779`, versus `.00708` at the old state. The holdout residual is
smaller. This pair has different grids, finite boxes and initial-value
lineage, so it does not establish temporal convergence or identify a slow
linear mode. The appropriate next measurement is the same diagnostic on
additional naturally produced true-RHS caches with their paired grids,
clock conversion and event records, accompanied by spatial and box
controls. No projected or predicted field is accepted by this analysis.

The full later report and fields are in
`continuous_inner_fresh_tpaa551c1e_6ecb_414f_89c7_4c610f2fd4b9/` beneath the
campaign acceleration output directory.
Its separately retained same-boundary sampling audit is
`continuous_sampling_tp7bcac57e_5924_41b2_b774_ba02a47a0f5f/`: the 65 by 49
and 129 by 97 core norms differ from Gauss by `+.14178%` and `+.03452%`,
and holdout norms by `-.02943%` and `-.0002024%`. Those are quadrature
differences of one fixed interpolant, not spatial PDE error estimates.
