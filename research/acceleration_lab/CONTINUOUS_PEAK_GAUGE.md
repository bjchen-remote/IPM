# Independent continuous C1 peak amplitude rule

This is a new research rule. It does not replace the maintained C rules,
change any native schema4 checkpoint, or establish a new accepted trajectory.
The initial work consists only of manufactured fields and saved 65 x 33
native frames. There are no new Poisson solves or time steps in this report.

## Definition and instantaneous algebra

On a fixed grid, write the wall gradient as
`v = (rho * Dx.')(1,:)`. Use the maintained seven-point `Dx` again for its
nodal slopes `d = v * Dx.'`. On each interval `[x_i,x_{i+1}]`, define the
cubic Hermite polynomial with endpoint values `v_i,v_{i+1}` and slopes
`d_i,d_{i+1}`. This defines one C1 function `H_h[v]` on the whole grid.
The interpolation map is linear in `v`; its slopes are not limited, clipped
or obtained from nonlinear PCHIP.

Choose one fixed closed search interval `I` and define

```
P_H(rho) = max_{a in I} H_h[(rho * Dx.')(1,:)](a).
```

For the positive wall branch in the saved tiny experiment, `I=[0,x(end)]`
is fixed before inspecting candidate rates. Each cubic contributes all
interior stationary roots and its clipped endpoints. The maximum therefore
does not depend on selecting the largest sampled node and fitting only that
node's three-point stencil.

For any two grid vectors `v,w`,

```
|max H_h[v] - max H_h[w]| <= ||H_h[v-w]||_infinity
                          <= C_h ||v-w||_infinity.
```

Here `C_h` is finite for the fixed interpolation operator and interval; it
is not asserted uniform over arbitrary increasingly singular grids. Also,
`P_H(s*rho)=s*P_H(rho)` for `s>0`. These properties show continuity of the
value functional and positive homogeneity. They do not imply that the peak
location or every higher derivative is smooth.

Suppose the maximizing point `a_H` is unique, isolated and interior, with
resolved negative one-sided curvature. For an evolving grid field,
linearity gives `partial_tau H_h[v]=H_h[v_tau]`. At the maximum the spatial
derivative is zero. Consequently the movement of `a_H` contributes no
first-order change to the peak value, and

```
P_H' = H_h[(rho_tau * Dx.')(1,:)](a_H).
```

This conclusion remains valid when an isolated maximum crosses a shared
native node: the Hermite function and its first spatial derivative agree
from both cells. Their second derivatives need not agree; neither a smooth
peak acceleration nor a universal time-integrator order follows from C1.

Let the independently assembled maintained split be
`rho_tau=F0(rho)+cOmega*rho`, with the instantaneous length and translation
rates held at their prescribed values in `F0`. For the present
`transport_anchor` length rule, `c_l=-U1(anchor,0)/anchor` depends on the
instantaneous field, not on the choice of `cOmega`. The saved symmetric
experiment has `c_r=0`. Put

```
B_H = H_h[(F0 * Dx.')(1,:)](a_H),
cOmega_H = -B_H/P_H.
```

Then `P_H'=B_H+cOmega_H*P_H=0` on this instantaneous smooth branch.
No target peak, normalization projection, relaxation gain, rate cap or
modified clock appears. Differentiating `cOmega_H` itself is unnecessary
for this first derivative: it is the scalar already present in `rho_tau`.

The same continuum coordinate-change interpretation would require
evolving `d(log C_omega)/d tau=cOmega_H` consistently with the full field
and the unchanged physical time rule. Different finite-grid gauge
discretizations are not automatically identical physical PDE trajectories.
That claim requires independent time and spatial comparisons.

## Nonsmooth and invalid cases

`ipm_accellab_hermite_peak(values,forcing,x,Dx,window,options)` records
stationarity residual, selected cell/fraction, active positions and values,
both one-sided curvatures at an active node, tolerances and validity reasons.
The default rejects nonpositive, flat/degenerate, boundary and multiple or
nearly tied maxima. Search-boundary or remesh changes are not hidden.

For finitely many exactly tied isolated maxima, the forward directional
derivative is the maximum of the forcing values at the active sites. Since
all active peak values equal `P_H`, a one-sided algebraic rate is possible:
`cOmega_H=-max_active B_H/P_H`. An optional synthetic-only test exercises
that rule. It is not the default, and an active set selected within a finite
floating-point tolerance is not a rigorous certificate of an exact tie.
Flat intervals would require a different active-set treatment and are
rejected. A tie can make the rate nondifferentiable even though the maximum
value remains continuous.

Changing the mesh changes the functional. This construction does not remove
interpolation errors or peak jumps introduced by remeshing. It also does
not guarantee a uniform lower curvature bound near a singular limit.

## Manufactured validation

`ipm_accellab_test_hermite_peak` passed on a nonuniform fixed grid. A moving
quadratic with changing amplitude crosses native cells; errors in peak
value, location and actual peak derivative were at most `5.565e-15`.
For a sampled moving Gaussian, the discrete peak derivative was checked
against central differences of the discrete maximum, without replacing it
by the continuum amplitude. Errors at time increments
`1e-2/1e-3/1e-4/1e-5` were
`1.3992e-6/1.3996e-8/1.4026e-10/9.0643e-12`.
Positive-homogeneity and forcing-linearity errors were at most `1.125e-15`.
The three default degeneracy guards and an explicit two-peak forward
directional derivative test also passed. The latter's first-order forward
difference errors decreased `2.482e-4 -> 2.484e-6`.

Report:
`result/longtime/20260908_campaign_v1/acceleration_lab/hermite_peak_mms_tpd02bec49_a100_462a_80c0_a41f02adea7a/report.json`.

## Saved original 65 x 33 switch frames

The input is **v2** of
`result/verification/performance_lab_20260908/quadratic_switch_native65_v2/switch_frames.mat`.
The source exporter separately assembled `physicalRhs` with all rates zero,
`baseRhs` with the native length/translation rates and zero amplitude rate,
and `rhs` with the original rule. It verified bitwise native replay of both
terminal rho and the complete trusted history. This analysis independently
checked, for all three frames, exact saved `Dx/Dy`,
`rhoX == rho*Dx.'`, and `rhs == baseRhs+cOmega*rho`.
The earlier v1 export is not used for this algebra.

| tau | Original quadratic P | Continuous Hermite P | Original cOmega | New algebraic cOmega_H | H peak derivative on original RHS |
| --- | --- | --- | --- | --- | --- |
| .04925 | .857239194891 | .858018906341 | .220273299661 | .215261977210 | +.004299809408 |
| .04950 | .858636300092 | .858018551770 | .212381504391 | .215172215890 | -.002394482239 |
| .04975 | .858636300092 | .858017955489 | .212313467295 | .215082366517 | -.002375765249 |

The Hermite locations were `.686393893538/.686286629767/.686179350325`,
with selected curvatures near `-3.384`. After constructing the proposed
instantaneous full field `baseRhs+cOmega_H*rho` and applying maintained
`Dx` to that full field, the Hermite peak derivative was at most `1.388e-16`.
The separate full-field differentiation versus linear derivative split
disagreed by at most `8.912e-16`, as expected from floating-point operations.
This is an instantaneous identity, not a new integration result.

To isolate continuity of the *functional*, the experiment also fixed the
offline path `v(lambda)=(1-lambda)*v_before+lambda*v_after` between the first
two saved wall fields. This prescribed path is not a PDE trajectory.
The original largest sampled node changes at `lambda=.297148672354`.
At equal distances `epsilon` to either side:

| epsilon | Original quadratic value difference | Hermite maximum value difference |
| --- | --- | --- |
| 1e-2 | .00139734451669 | -7.24936610741e-9 |
| 1e-4 | .00139735634795 | -7.24936777274e-11 |
| 1e-6 | .00139735646626 | -7.24975635080e-13 |
| 1e-8 | .00139735646744 | -7.32747196253e-15 |

Thus the old selector yields a finite value jump on this fixed field path,
while the continuous maximum approaches a common value. This does not
identify every event in the large-grid run or assert new time integration
conserves the proposed observable.

The original quadratic-gauge trajectory does **not** conserve Hermite P.
Its measured interval slopes were `-.001418280423/-.002385127207`.
Across the old selector event, the first interval's slope differs from the
mean endpoint instantaneous derivatives by `-.002370944007`; endpoint
trapezoids do not resolve that discontinuous-rate event. The next, smooth
interval differs by only `-3.463e-9`. Neither the old rowwise zero derivative
nor the new instantaneous cancellation can substitute for a finite-time
conservation check.

Report:
`result/longtime/20260908_campaign_v1/acceleration_lab/hermite_peak_switch_tp8f640587_3267_475f_bd2d_a1ad4acfef24/report.json`.

## Original tiny time-integration proposal and subsequent execution

The following proposal was subsequently authorized and executed as an
independent research experiment. See
[CONTINUOUS_PEAK_TIME_TRIAL.md](CONTINUOUS_PEAK_TIME_TRIAL.md) for the
complete finite-time results, unresolved full-field convergence order,
physical-coordinate comparison, local-tolerance audit and later q512
instantaneous evaluation. The proposal below records its original scope.

The current evidence justifies preparing a separate tiny experiment, subject
to the root's computation schedule. Use the same 65 x 33 initial state,
fixed grid, transport anchor and open-boundary discretization. Integrate
the new full RHS and all scale variables consistently; do not project rho
or reset P after any stage. Retain the original quadratic experiment as
an independent reference with its original result and rejection records.

Before running, freeze the time horizon and step-size sequence; an initial
proposal is `tauEnd=.08`, `dt=[.0005,.00025,.000125]`, which straddles the
known old selector time. This is a proposal, not a performed or accepted
trial. Save every stage's Hermite maximum, PPrime, cOmega, active cell,
one-sided curvature and all failed validity guards. Report actual peak
drift over time and step refinement, rather than merely the rowwise
algebraic residual. Record original quadratic P as a secondary observer.
No general formal time order is assumed across C1 cell changes.

The first decision is whether the new rule removes a discontinuous C-rate
event without introducing a degenerate peak, unstable RHS or unrefined
finite-step drift. Only after that would a separate physical-coordinate
comparison at common physical times and further spatial refinement be
meaningful. These checks still would not establish strong singularity or
the long-time limiting profile.
