# Experiments

The active numerical implementation lives in `+ipm/`, and `ipm.solve` is the
only simulation entry point. This directory contains reproducible option sets
and experiment notes. Examples below assume the current MATLAB directory is
`ipm_structured/` and its root has been added with `addpath(pwd)`.

```matlab
addpath('research/experiments');
result = ipm.solve(strong_mesh_case('production'));
```

- `production`: accepted 1025-by-513 physical run on the box `[-400,400]`
  by `[0,400]`, with `maxDt=2.5e-4`.
- `spatial_513`: 513-by-257 comparison with the same box and time step.

Rejected parameter sweeps remain historical data in the original solver's
`result/archive`; their bulk artifacts were not copied here. They are not
alternate solver entry points.

The regularized CCF-like heavy-tail family is a historical, opt-in comparison
configured by `ccf_heavy_tail_case`; it is not the no-argument `ipm.solve`
datum:

```matlab
addpath('research/experiments');
physical = ipm.solve(ccf_heavy_tail_case('screen_physical'));
dynamic = ipm.solve(ccf_heavy_tail_case('screen_dynamic'));
```

`box_large` repeats the physical test on a doubled box, while `epsilon_fine`
halves the onset regularization together with the target center spacing.
Comparing these cases is mandatory: the limiting wall derivative is singular
at `x_1=1` and decays only as `x_1^(-1/2)`.

`long_physical` and `long_dynamic` freeze the historical CCF setup explicitly
and use the doubled box through physical time
`0.5`, use the legacy-inspired flat-core `lattice` profile with 96 candidate
cells on the 90% core, and anchor the far grid outside the monitored feature.
Only adjacent cell ratios are constrained; global ratios are diagnostic and
may exceed `1e3`. The hard adjacent ratio is `1.5`. `long_dynamic` combines
this grid with the fixed `X=1` transport-anchor length gauge and a
full-gradient-energy amplitude gauge. They remain verification runs until the
long-time quality and refinement comparisons in the notes are accepted.

These cases enable analytic initial pre-remeshing around the paired off-origin
peaks. See `ccf_heavy_tail_notes.md` for the verified 2026-08-24 runs and the
explicit distinction between strong concentration and a blow-up claim.

## Smooth near-Bessel two-scale screens

`bessel_two_scale_case` supplies two physical-mode upper-half-plane tests to
the same `ipm.solve` entry point:

```matlab
addpath('research/experiments','research/analysis');

k0Reference = ipm.solve(bessel_two_scale_case('screen_unperturbed'));
k0Perturbed = ipm.solve(bessel_two_scale_case('screen_perturbed'));
k0Diagnostics = ipm_analyze_bessel_two_scale(k0Perturbed);
k0Comparison = ipm_compare_bessel_two_scale_pair(k0Reference,k0Perturbed);

k1ReferenceOptions = bessel_two_scale_case('screen_k1_perturbed',struct( ...
    'caseMetadata',struct('parameters',struct( ...
        'perturbationAmplitude',0)), ...
    'resultFile','result/verification/bessel_k1_two_scale_screen_unperturbed.mat'));
k1Reference = ipm.solve(k1ReferenceOptions);
k1Perturbed = ipm.solve(bessel_two_scale_case('screen_k1_perturbed'));
k1Diagnostics = ipm_analyze_bessel_k1_two_scale(k1Perturbed);
k1Comparison = ipm_compare_bessel_two_scale_pair(k1Reference,k1Perturbed);
```

Every Bessel-specific research parameter lives under
`caseMetadata.parameters`; solver options remain strict and typo-checked.

The K0 density has a pole below the wall and is smooth on the computational
half-plane, but its analytic streamfunction has through-wall flow. The solver
therefore recomputes an impermeable velocity; this is not an exact steady-pair
regression. The pure regularized K1-core datum has exactly zero wall density.
Physical evolution no longer depends on a scaling-gauge peak, so maintained K1
cases keep the monitor amplitude at zero and disable the inapplicable
positive-wall-omega sign gate. Run
`verify_bessel_two_scale_initial_data` to audit the homogeneous-wall datum.

The 257-by-129, `t=0.05` screens are deliberately labeled screening evidence:
the matched 3% perturbations are essentially neutral over this short interval.
Use `refined_perturbed`, `refined_k1_perturbed`, and longer user overrides
before making any attraction claim.

## Fourth-order active-datum campaign

`fourth_order_blowup_case` replaces the maintained active case's numerical
tuple by `high_order/weno5_fd/ssprk54/high_order` without changing its datum,
dynamic gauge, Green artificial boundary, or open outer transport.  The
`pilot_256` and `pilot_256_dt_half` cases use 257 nodes, hence 256 stored cells
per coordinate direction, and provide the required first time-step check:

```matlab
addpath('research/experiments','research/analysis');
pilot = ipm.solve(fourth_order_blowup_case('pilot_256'));
half = ipm.solve(fourth_order_blowup_case('pilot_256_dt_half'));
pilotAudit = ipm_assess_fourth_order_blowup(pilot);
halfAudit = ipm_assess_fourth_order_blowup(half);
verify_fourth_order_blowup_setup
```

`research/analysis/verify_fourth_order_blowup_setup.m` audits the five
analytic 256-cell factories, the dense 384-cell factory, and one retained
legacy factory: the complete fourth-order tuple, profile/version metadata,
grid quality, reserved metadata, and the continuous-prefix helper contract.
Saved-run assessments apply that helper separately to each result.

The high-order Poisson path already builds `decomposition(A,'lu')` once and
reuses it at every step.  On the actual nonsymmetric high-order matrix,
MATLAB's automatic `decomposition(A)` also selects LU.  The automatic and
explicit forms are therefore equivalent here; the factorization is not the
resolution or wall-oscillation bottleneck.

The initialization scan selected `premeshed_256`: two analytic `lattice`
passes from the current grid, outer-anchor factor 2, adjacent-cell cap 1.5,
and `maxRemeshes=2` so the accepted initial grid is frozen during evolution.
It improves the initial horizontal 90% core from about 3.24 to 9.25 points;
it does not relax the trusted-mask or hard-stop thresholds.

For the symmetry-reduced interpretation of 256-by-256 cells in the first
quadrant, use `premeshed_quadrant_256`.  It stores 513-by-257 nodes over the
full horizontal/upper-half-plane box; the same two accepted pre-remesh passes
give about 27 horizontal 90% core points while retaining the 1.5 adjacent-cell
cap.

All archived non-analytic smooth-grid names--the `smooth_*`, `focused`,
`balanced`, `concentrated`, `near_terminal`, `stable_limit`, `narrow_limit`,
`moderate_narrow_limit`, and `limit` cases--explicitly retain the
`legacy_abs_gaussian` monitor for reproducibility.  A later regularity audit
found that its use of `abs(x)` is not differentiable at the origin.  These
results remain useful numerical trajectory and mapping-stability comparisons,
but they are not evidence from a globally smooth fourth-order map.  The helper
now defaults to `analytic_double_gaussian`; analytic cases are versioned under
separate names and must be rerun independently.

`smooth_quadrant_256` replaces the piecewise lattice transition by the legacy
Gaussian monitor centered on the two analytic wall-derivative peaks.  Its
513-node horizontal grid has about 16 initial 90% core points, adjacent ratio
about 1.153, and maximum `log(dx)` curvature about 0.0188.  It is the preferred
archived mapping-stability comparison when the more concentrated lattice grid
develops secondary WENO peaks.

`smooth_box_quadrant_256` repeats that construction on
`[-281,281] x [0,132]` with a milder monitor.  Its initial horizontal core is
15.582 points, adjacent ratio 1.04127, and `log(dx)` curvature 0.002540.
The box put the far-source ratio just below 1% initially, but it rises to
0.01018 by `t=0.5`.  The far-velocity ratio remains large for this slowly
decaying datum, so this is a source-truncation comparison rather than a fully
far-field-converged claim.

`focused_box_quadrant_256` keeps the same box and center spacing but raises
the legacy monitor strength from 0.5 to 1.0.  This gives 19.555 initial
horizontal core points while keeping the adjacent ratio at 1.05487.  It is
the fixed-513-by-257 long-time continuation case; it does not relax any
quality or trusted-evidence threshold.

`balanced_box_quadrant_256` follows the observed long-time limiting
direction: it uses monitor strength 2 in x and halves the near-wall y spacing
to 0.01.  The resulting adjacent ratios remain about 1.075 in x and 1.026 in
y, its initial x/y core counts are 26.363/64.961, and all high-order
quadrature weights remain positive.  This is the archived fixed-size test used
when the x-only focused case exhausts its y core.

`concentrated_box_quadrant_256` is the next legacy fixed-size stress test.  It
uses monitor strength 3 and near-wall y spacing 0.005, while retaining adjacent
ratios near 1.093/1.029, initial x/y core counts 32.012/84.860, and positive
high-order quadrature.  Its requested `t=1.5` run is locally trusted only
through `t=1.334758`; the under-resolved terminal state is not evidence.

`near_terminal_box_quadrant_256` is an aggressive legacy diagnostic follow-up
to the formal `T` near 2 extrapolation.  It uses monitor strength 12 and near-wall y
spacing 0.00125; its initial x/y core counts are 58.205/115.086, its adjacent
ratios are 1.23913/1.03538, and its high-order quadrature remains positive.
The run stops at raw `t=1.779629` and has a continuous locally trusted prefix
through `t=1.621842`, where gradient/wall growth is 4.839197/5.113862 and all
four blow-up diagnostics are false.

`stable_limit_box_quadrant_256` is the strongest retained same-width legacy
map, with `A=16,w=0.7`.  Its initial x/y core counts are 63.690/115.174,
adjacent ratios 1.30235/1.03538, and horizontal `log(dx)` curvature 0.09955.
It stops at raw `t=1.790679`; its continuous locally trusted prefix ends at
`t=1.646331`, 75.1 times the 513-by-257 lattice evidence time.  Gradient/wall
growth is 5.189615/5.474471, while all four legacy point-window candidates
remain false.  Point-count windows over-weight densely written late records;
the growth-rate assessment below supersedes them for model classification.
This run also has the legacy-map regularity defect and is unsuitable as a
blow-up claim.

More aggressive legacy amplitude branches `A=20/24/28` were rejected by the
wall-profile oscillation audit.  The `A=20,w=0.7,t=0.2` short screen is not a
long-time acceptance.  Narrow `A=16` monitors are also non-monotone:
`narrow_limit_box_quadrant_256` (`w=0.35`) stops at raw `t=0.275560`
(trusted through 0.269956), and the `w=0.45`
stops at raw `t=0.453370` (trusted through 0.449072).  The `w=0.55` screen is
confirmed only through `t=0.5`;
`moderate_narrow_limit_box_quadrant_256` is retained for configuration
provenance, but its redundant legacy long run was terminated after the
regularity defect was identified.  `limit_box_quadrant_256` is the rejected
`A=28,w=0.35` limit branch, not a preferred production case.

`analytic_stable_limit_box_quadrant_256` replaces the legacy `abs(x)` monitor
by the globally analytic two-Gaussian profile at the same `A=16,w=0.7` and is
the required independent rerun.  Its profile-v2 grid uses `fineCount=131073`;
the initial x/y core counts are 63.408/115.086, adjacent ratios
1.30962/1.03538, and horizontal/vertical `log(dx)` curvatures
0.076436/0.001206.  Minimum seven-point stencil rconds are
`8.804e-6/9.565e-6`, and normalized quadrature margins are
`0.0063475/0.00077649`.

The completed analytic baseline stops at raw `t=1.791544` with
`grid_resolution_failure`; 664 of 869 records form the continuous locally
trusted prefix through `t=1.644487`, 75.051 times the 513-by-257 lattice
evidence time and only 0.112% below the legacy A=16 prefix.  Endpoint
gradient/wall growth is 5.161853/5.445129, with 8.028/14.694 x/y core points.
All four legacy point-window candidates are false.  The maintained assessment
now also applies `ipm.diagnostics.maximumGrowthRateFit`: it estimates
`gamma=d(log M)/dt` with three local-polynomial bandwidths, uses physical-time
weights, and compares exponential, finite-time power, and double-exponential
models on fixed physical-time scoring points and e-fold windows.  It also
downgrades decisions that change after uniform resampling of the full input.
For this `w=0.70` run, the shortest gradient window is sampling-sensitive and
the other three prefer finite-time power, with raw-series fits
`T=1.99978--2.06671` and `p=1.02945--1.18304`; the available gradient range is
only `1.64130` e-fold, so this is not a credible asymptotic classification.
The assessment also applies
`ipm.diagnostics.canonicalMaximumGrowthFit` to `common.canonicalTau`.  That
separate diagnostic compares a maximum linear in canonical tau against an
exponential maximum on one log-value error domain, with a quadratic-log
curvature guard; its fields must not be interpreted as physical singular
times or physical-time rates.

The analytic amplitude screen rejects stronger mappings through the wall
profile gate: `analytic_ratio_limit_box_quadrant_256` (`A=24`) has raw/trusted
times `0.132121853/0.128389363`,
`analytic_strong_limit_box_quadrant_256` (`A=20`) at
`0.279627095/0.274291876`, and
`analytic_intermediate_limit_box_quadrant_256` (`A=18`) at
`0.369873884/0.366480746`.  Each stops for `wall_profile_oscillation`, and
each has four false candidates.  The registered
`analytic_moderate_narrow_limit_box_quadrant_256` uses `A=16,w=0.55` and is
the completed current optimum.  Its profile-v2, `fineCount=131073` grid has
initial x/y core counts `72.056164/115.151066`, adjacent ratios
`1.301284189/1.035383549`, and `log(dx)`/`log(dy)` curvatures
`0.074148338/0.001205563`.  Minimum seven-point stencil rconds are
`9.0403e-6/9.5650e-6`; normalized and unnormalized minimum quadrature
weights are `0.0054952/0.00077649` and `0.006020124/0.000398820`.

The narrow analytic run stops for `grid_resolution_failure` at raw
`t=1.79749926110`; 681 of 880 records form its continuous locally trusted
prefix through `t=1.65974756404`.  This is 75.747438 times the 513-by-257
lattice evidence time, 0.927975% beyond the broad analytic baseline, and
0.814928% beyond the legacy `A=16,w=0.7` prefix.  Endpoint safety is
`0.99808118`, x/y core counts are `8.01538/13.76353`, and gradient/wall
growth is `5.397535/5.687645`.  All four candidates remain false.

Its old point-count `R^2` gates remain false, but those gates are no longer the
primary model classifier.  The growth-rate diagnostic gives endpoint gradient
`gamma=2.96936/2.96883/3.01407` across the three bandwidths.  The final
`0.25/0.5/1/1.5` e-fold gradient decisions are
`indistinguishable/power/power/power`, with raw-series power fits
`T=2.00270/2.03537/2.04522/2.06371` and
`p=1.03316/1.12076/1.14317/1.17851`.  The wall peak has the same final
decision sequence; its raw shortest-window preference is double exponential,
but that preference is sampling-sensitive.  The wall scale-gauge endpoint
rate contribution is `2.883854`, while the
three differentiated rates are `2.88196/2.88051/2.87833`; their full-prefix
relative RMS discrepancies are `0.78%--1.02%`.

The registered node-refinement follow-up is
`analytic_moderate_narrow_quadrant_384`.  It keeps the analytic
`A=16,w=0.55` datum and box, but stores 769-by-385 nodes: because the full
horizontal domain contains both symmetric halves, this is a
384-by-384-cell grid in the first quadrant.  Both target-center-spacing
components are exactly two thirds of the retained 256-cell case
(`0.0193174285333333` in x and `0.000833333333333333` in y), and the case
requests physical time `2.05`.  Run the complete solve, trusted-prefix
assessment, and maximum-growth plot with:

```matlab
addpath('research/experiments','research/analysis');
products = run_fourth_order_dense_growth();
```

The requested final time is a computational target, not an evidence cutoff.
The runner saves the raw result separately, evaluates only the continuous
locally trusted prefix, and exports the corresponding growth-curve PNG.  A
512-by-512 first-quadrant case is not the first long run on the current 16 GB
host: its direct-LU memory and measured short-run cost make the 384-cell case
the practical first refinement.  This is a run-order decision, not acceptance
or rejection of the 512-cell refinement.

The completed 384-cell run took 7647.862 seconds and stopped for
`grid_resolution_failure` at raw physical `t=1.85594004481`; 802 of 1007
records form its continuous trusted prefix through `t=1.75158228257`.  This
extends the 256-cell trusted endpoint by 5.533%, while the gradient and wall
peak grow by factors 7.377904 and 7.711241.  At the coarse-grid trusted
endpoint on the common interval, the two gradient maxima differ by 0.079%.
The 384-cell gradient has 1.998490 e-folds and endpoint logarithmic rate
`gamma=3.89317/3.89006/3.92076` across the three bandwidths.  Its 0.5, 1.0,
and 1.5 e-fold windows select finite-time power fits, with gradient
`T=2.03497/2.03821/2.04543` and `p=1.11764/1.12755/1.14548`; the 0.25
e-fold window remains sampling-sensitive and inconclusive.  Consequently this
is a spatially consistent finite-time-power candidate, not a credible blow-up
solution: the gradient range is still below the three-e-fold evidence gate,
there is no late half-step or larger-box replication, and the trusted-endpoint
far-velocity ratio is about 0.93.  The individual and resolution-comparison
plots are archived as
`fourth_order_active_analytic_A16_w055_q384_t205_dt002_growth.png` and
`fourth_order_active_A16_w055_q256_vs_q384_growth.png`.

The same trusted prefixes must also be inspected in the increasing canonical
rescaled clock `history.common.canonicalTau`; this is not the remaining-time
variable `T-t`.  On the 384-cell run, canonical tau reaches 4.006902.  Over
the last `Delta tau=0.5`, the physical gradient and wall peak behave as
`exp(0.51635*tau)` and `exp(0.50317*tau)`, while held-out errors for a maximum
that is linear in tau are roughly one order of magnitude larger.  At the
common q256/q384 canonical-tau endpoint, the fitted exponential rates agree
to about `2e-4` in absolute slope.  Since
`dt/dtau=C_omega/C_l` and the wall gauge keeps the rescaled wall peak nearly
constant, this canonical-tau exponential and the near-Type-I physical-time
power law are two clock descriptions of the same candidate mechanism, not
competing classifications.  The dedicated e-fold-window guard selects a
local exponential on the last 0.25 and 0.5 e-fold windows, while the longer
windows remain curved/nonasymptotic or indistinguishable.  The wall relation
is partly a gauge identity;
independent evidence must come from the global rescaled gradient, feature
widths, shape convergence, and a third grid/time-step/box study.  The
canonical-time comparison is archived as
`fourth_order_active_A16_w055_q256_vs_q384_tau_growth.png` with its sibling
`_analysis.mat` data file.

Analytic/legacy scalar diagnostics differ by about `1e-4--6e-4` at
`t=0.96/1.3/1.644`.  A lone 9.4% FWHM difference at `t=0.5` accompanies only
an `8.8e-5` core-width difference and is a discrete crossing-branch/sampling
alias, not a 9.4% trajectory discrepancy.  Width conclusions must therefore
use the connected core and neighboring diagnostics rather than one aliased
FWHM sample.

Between analytic `w=0.55` and `w=0.7`, gradient, wall, and core-width
differences at `t=0.96/1.3/1.644` are mostly `1e-4--8e-4`.  Some early
samples show 1.1%--1.3% gradient and 6%--7% FWHM branches, but at matched
physical time `t=0.5` the gradient/core-width differences are only
`4e-6/3e-6`.  This again identifies the isolated FWHM discrepancy as a
crossing/sampling alias rather than a physical trajectory split.

The legacy trajectories have time-step evidence only through `t=0.12`,
large/small-box evidence through `t=0.5`, and same-size mapping overlap through
`t=0.96`.  The analytic cases agree closely through their common prefix, but
near the optimum `t=1.659748` there is still no half-step, larger-box, or
higher-node-count corroboration.  The optimum trusted endpoint has
far-velocity/source ratios `0.923439/0.0144411`, physical/rescaled mass drift
`-0.929168/-0.346311`, range violation `7.434e-8`, divergence `4.625e-13`,
and Poisson residual `1.569e-11`; these are material late-time limitations.
The archived optimum and analytic campaign summary are
[`fourth_order_active_analytic_A16_w055_t195_dt002.mat`](../../result/verification/fourth_order_active_analytic_A16_w055_t195_dt002.mat)
and
[`fourth_order_active_analytic_campaign_summary.mat`](../../result/verification/fourth_order_active_analytic_campaign_summary.mat).
The focused fit comparison is archived as
[`fourth_order_active_growth_rate_comparison.mat`](../../result/verification/fourth_order_active_growth_rate_comparison.mat).
The
dynamic/Green/open problem remains outside the
global fixed-grid physical fourth-order certificate.  The smooth fixed-map
cases themselves disable adaptive remapping; trigger-driven high-order
remapping is an additional limitation only where enabled.  Optimize the last
continuous trusted time, not the later hard stop time, and do not relax range,
oscillation, core-resolution, or cell-ratio thresholds.

The analytic `w=0.55` case supplies a reproducible finite-time-power candidate
for the maximum gradient, not a credible blow-up solution.  Simple exponential
growth is not the best description on the present trusted interval, but power
and double exponential are not yet credibly separated: both shortest windows
are sampling-sensitive, the data span only `1.69--1.74` e-fold, and
late-time step, box, and node-refinement checks remain incomplete.  No credible
blow-up claim is currently supported.

### Lorentzian LAT optimization and 512-cell refinement

The 2026-09-04 LAT study first optimized the fixed analytic map on the
256-by-256 first-quadrant grid.  The selected profile is a symmetric double
Lorentzian centered at solver-grid `x=1.14`, with width `0.20` and strength
`A=28`.  The relevant horizontal mapping diagnostics improve from the old
q256 Gaussian map's maximum adjacent-cell ratio and maximum `log(dx)`
curvature, `1.301284189/0.0741483378`, to
`1.10097818/0.00240067378`.  These are local smoothness diagnostics for the
fixed LAT map, not global `max(dx)/min(dx)` ratios.

The optimized q256 run takes 2732 steps and stops at raw physical time
`1.83683959412` through the grid-resolution gate.  Its continuous trusted
prefix ends at physical time `1.72352016681` and canonical
`tau=3.80105633`.  The final refinement stores 1025-by-513 nodes, which is a
512-by-512-cell grid in the first quadrant.  It reaches the configured
4800-step `maxSteps` budget with raw=trusted physical time
`1.83433595901`, canonical `tau=4.77268406`, x/y core counts
`8.28986/12.1711`, and safety `0.965034`.  All 956 of 956 recorded states
are trusted: this is a budget stop, not a numerical failure.  Peak resident
memory was `5.423628288 GB`, with zero swap.

Over the common physical-time interval, q256/q512 endpoint relative errors
are `8.719e-5` for the maximum gradient and `4.207e-5` for the tracked wall
peak; their log-RMS errors are `2.327e-4/1.660e-4`.  On q512, the physical
1.5-e-fold power fits are
`T=2.03899451,p=1.12989826` for the gradient and
`T=2.03623379,p=1.09028074` for the wall peak.  Over the final quarter of
the canonical-time span, the exponential rates are `0.5059/0.4913`; errors
`0.0012/0.0016` are lower than the corresponding linear-model errors
`0.0129/0.0117`.  The full 1.5-e-fold canonical windows nevertheless remain
`curved_nonasymptotic`, so the short-window exponential behavior is not an
asymptotic classification.

The gradient/wall dynamic ranges are only `2.382/2.415` e-fold, below the
three-e-fold evidence gate.  The maximum far-source ratio is `2.716%`, the
far-velocity ratio is about `94.9%`, and no late-time half-step, larger-box,
or third matched-grid replication is available.  The refinement assessment
therefore records `credible=false`.  The result supports spatial agreement
and the candidate dual-clock interpretation--near-Type-I power growth in
physical time together with locally exponential maximum growth in canonical
tau--but it is not a credible blow-up solution.

The reproducible workflow is split among these runners and analysis files:

- `run_fourth_order_lattice_screen.m` runs the short q256 LAT candidates and
  records ranking diagnostics.
- `run_fourth_order_lattice_validation.m` performs the selected q256 long
  validation and trusted-prefix analysis.
- `run_fourth_order_lattice_512.m` runs the final q512 case and archives its
  single-grid and q256/q512 comparison products.
- `ipm_assess_lattice_refinement.m` computes common-interval refinement
  errors and applies the conservative credibility audit.
- `ipm_plot_canonical_maximum_models.m` compares canonical linear and
  exponential maximum models on common e-fold windows and plots their rates.

The principal raw results are
[`fourth_order_lat_validation_analytic_lorentzian_a28_w020_c114_quadrant_256_q256_t210.mat`](../../result/verification/fourth_order_lat_validation_analytic_lorentzian_a28_w020_c114_quadrant_256_q256_t210.mat)
and
[`fourth_order_active_analytic_lorentzian_A28_w020_c114_q512_t210_dt002.mat`](../../result/verification/fourth_order_active_analytic_lorentzian_A28_w020_c114_q512_t210_dt002.mat).
The q256/q512 growth plots are archived as
`fourth_order_active_lorentzian_A28_q256_vs_q512_growth.png`,
`fourth_order_active_lorentzian_A28_q256_vs_q512_tau_growth.png`, and
`fourth_order_active_lorentzian_A28_q256_vs_q512_canonical_models.png`.

## Scope

Only case factories and initial-data functions present in this directory are
part of the executable experiment contract. Historical designs for unshipped
case families remain research notes, not runnable examples. Add a case here
only together with its implementation and a focused verification.
