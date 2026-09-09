# Projection-free Jacobian-free profile experiment

The new tail-constrained candidate passed all predeclared **tiny independent
profile-candidate** tests on 49 by 25 and 65 by 33 grids. Its actual core
shape residual decreased about 0.3%, and this small benefit survived
unchanged-PDE evolution to a common physical time. It is not a production
trajectory, an acceleration of the original initial-value solution, or
evidence that long-time linear convergence has been overcome.

The first, unconstrained candidate remains rejected. Although it reduced
the core residual by 2.5--2.8%, it changed the sampled tail density by about
4.4%, tail velocity by about 0.145%, and physical peak by about 0.3%.
Its complete results and executed source are preserved in
`jf_profile_tp40145ebe_18f5_4340_ad83_acfac42a0e9e/`. The second trial added
constraints to the construction; it did not relax these failed acceptance
gates or relabel that result.

## Construction and mathematical scope

The source is the original Gaussian density on `[-4,4] x [0,4]`, evolved
to canonical time `.08` using the unchanged quadratic-peak/transport-anchor
PDE. The already saved native-rule run at `dt=.000125` supplies rho and
all five scale variables. Each grid builds one tiny maintained operator;
all probes call the original Poisson/RHS with that same operator.

At the source rho, let `F(rho)` denote the complete original density RHS.
Build at most three whole-field directions spanning
`F, J_F F, J_F^2 F`, with centered Jacobian-vector products and two-pass
weighted orthogonalization. Each direction has the same weighted L2 norm
as rho. Thus the coefficient norm equals the relative field-step norm.
No local field patch, P normalization, or state projection is applied:

\[
\rho_{candidate}=\rho+Vc.
\]

The shape observer is the complete continuous C1 `G=U_tau` described in
`CONTINUOUS_INNER_DERIVATIVES.md`. Its genuine `PPrime` and geometric rates
are recomputed from every trial's original RHS. The fitting residual is
`G/P`, with the **candidate's actual P**, evaluated at a fixed baseline
split-cell Gauss observation grid. This extra division is deliberate:
for isotropic physical reconstruction,

\[
U_t=(C_l/C_\omega)G,\qquad
\omega_{physical,peak}=(C_l/C_\omega)P,
\]

so `G/P = U_t / omega_physical,peak`. It prevents a pure overall-amplitude
change from automatically looking like improved intrinsic shape dynamics.
Hard acceptance still requires actual `U_t` to decrease. Independent
`.999` and `1.001` amplitude controls showed the raw shape rate scaled by
those factors while `G/P` changed by at most `3.3e-13` on these fixtures.
This numerical homogeneity check is specific to these samples.

The Jacobian of this *complete* shape observer is differenced in each
whole-field direction at relative field increments `1e-5` and `5e-6`.
Both the continuous geometry-cell signature and the original quadratic
peak selector must remain unchanged across the probes. The relative
Frobenius difference of the two Jacobians was below `9.9e-9`. Singular
values below `1e-8` of the largest are excluded; all three directions
were retained in these runs.

For the second protocol, form a positive coefficient metric

\[
Q=I/r^2+T_\rho^TT_\rho+T_u^TT_u+p_P^Tp_P,\quad r=.002.
\]

Here the tail-density and two-component tail-velocity derivative rows use
the fixed physical observation nodes, their fixed weights, baseline tail
norms, and 80% of the unchanged `.001` acceptance tolerances. `p_P` is the
actual physical-peak derivative divided by 80% of its `.001` matching
tolerance. The candidate solves the three-dimensional constrained
Gauss--Newton model

\[
\min_c\|R+J_R c\|^2\quad\text{subject to}\quad c^TQc\le1.
\]

Cholesky `Q=L' L` and `q=Lc` reduce this to a Euclidean trust-region
least-squares problem. The scalar regularization multiplier is determined
by bisection on its norm equation. This optimizes coefficients before
constructing a field; it does not project the resulting rho, its P, or its
time derivative. Since every contribution to `Q` is nonnegative, the
ellipsoid bounds each *linearized* tail/peak matching quantity separately.
The real candidate is then checked with a fresh original Poisson/RHS;
the linearized constraint is never treated as a substitute for that check.

## Predeclared acceptance tests

Both actual physical-time core residual and intrinsic core `G/P` must
decrease at least 0.1%. Holdout L2 and sampled maxima, wall L2 and sampled
maxima, and native-FD core/holdout residuals must not increase. In addition
to candidate-native masks, the new constrained trial checks native
residuals on the fixed baseline masks, retaining both observations.
The final residual norm is independently integrated after splitting the
candidate's own native mapped cells and the same inner rectangle
boundaries; it does not reuse the fitting grid as its only metric.

Physical matching uses the explicitly finite rectangle `[-3,3] x [0,3]`,
sampled on 49 by 25 common laboratory nodes. The tail consists of those
nodes with `|x|>=2.5` or `y>=2`. This is **not the full computational box**.
Tail density relative L2, source infinity over baseline physical peak,
and velocity relative L2 each have a `.001` hard limit. Physical peak
relative change must stay below `.001`; density relative infinity on the
observation window must stay below `.003`.

Every candidate, including a rejected instantaneous one, was also advanced
with the original SSPRK54 PDE from its unchanged physical clock to
`t=.10`, alongside the source baseline. Rho and all five scale variables
evolved together. Two maximum canonical steps, `.001` and `.0005`, were
used. The actual RK states bracketing `.10` supplied C1 temporal
interpolation, using density and scale rates divided by `dt/dtau`, to
observe both branches at the same physical time. No P correction or
clock relabeling occurred. These are independently labeled research
initial-value branches with finite/CFL checks, not an assertion of a
fully trusted native continuation or checkpoint.

The same residual and matching tests must pass after both relaxations.
Time refinement must change core residual by less than `1e-4` relatively
and physical density by less than `1e-6`. Finally a cost control advances
the original source naturally for at least the candidate construction's
RHS budget, and requires the candidate's raw/intrinsic core norms to be
at least 0.1% smaller. This last comparison is at distinct physical times;
it is a cost diagnostic, not the common-time PDE validation.

## Actual constrained results

The second output is
`jf_profile_tpe83cbd24_ced9_4032_800c_a94a7a505d7b/` under the campaign
acceleration output directory. All initial and final hard tests passed on
both grids. There was no damping scan or subsequent acceptance threshold
change.

| Measurement | 49 by 25 | 65 by 33 |
|---|---:|---:|
| Relative whole-field step | .000523450 | .000523487 |
| Instantaneous physical core residual ratio | .99680539 | .99712767 |
| Instantaneous intrinsic core ratio | .99638036 | .99667922 |
| Instantaneous holdout ratio | .99783557 | .99801956 |
| Tail density relative L2 | .000544574 | .000521588 |
| Tail velocity relative L2 | .000343369 | .000349155 |
| Physical peak relative change | .000426570 | .000449951 |
| Observation-window density relative infinity | .000508860 | .000509145 |
| Core ratio after relaxation to physical `t=.10` | .99700193 | .99695889 |
| Holdout ratio after relaxation | .99833721 | .99826449 |
| Tail velocity after relaxation | .000345680 | .000351438 |
| Linear-model/actual decrease agreement | .998449 | .998973 |

The source actual continuous `PPrime` values were `.0106784758922` and
`.000151657932918`; the candidates' were `.0108051316986` and
`.000220938441008`. They were never forced to zero. The physical peak
matching constraint bounds a static candidate discrepancy; it is not a
restoring term in the evolution rule.

Each construction used 20 fresh RHS evaluations, including the independent
amplitude controls. The natural-progress control used four original RK
steps, or 21 RHS evaluations from a cold first-stage cache. Candidate
core residuals divided by those naturally advanced norms were `.9929551`
and `.9943521`. A whole-field tangent projection estimates canonical
time shifts `.0005790` and `.0005830`; these numbers are diagnostics only,
and no such shifts were applied to the candidates' clocks.

The two-grid baseline core residuals were `.31307` and `.31629` per
physical time, differing about 1%, which is larger than the 0.3% candidate
benefit. Moreover, the natural baseline core residual **increased slightly**
over the short cost-control interval. This is not an already linearly
converging long-time profile. Passing both grids demonstrates a reproducible
small constrained optimization step and matching-time persistence; it does
not demonstrate a converged spatial error, a better asymptotic rate, or a
cost advantage for long-time production. The constrained family remains a
small research result, not a reason to consume a large-grid window yet.

The `.001` versus `.0005` relaxation comparison changed physical density
by at most `1.126e-12` relatively and core residual by at most `2.031e-11`.
These small differences bound sensitivity to the two tested time steps;
they do not establish a general fourth-order temporal error law. Complete
failed and passing reports, candidate fields, original source fields,
Jacobians, coefficient metrics, and endpoint fields are retained. There
were no large-grid operator builds, large-grid PDE steps, native
checkpoints, mainline code edits, or changes to any older secant decision.

Run the fixed constrained protocol with:

```matlab
addpath('ipm_structured/research/acceleration_lab');
root = 'ipm_structured/result/longtime/20260908_campaign_v1/acceleration_lab';
ipm_accellab_test_jf_profile(fullfile(root, ...
  'hermite_peak_integration_tp995ad07e_6c1e_4745_8219_f53d760e2e6c'), ...
  root,'tail_ellipsoid');
```

This uses unique output directories and remains limited to the registered
49/65 original-equation fixtures. `unconstrained` retains the original
construction as a distinct protocol; it cannot overwrite or reaccept an
earlier run.
