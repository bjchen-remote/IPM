# Exact-gauge laboratory

This directory isolates amplitude- and length-gauge experiments from the
maintained `ipm.solve` entry point.  In the amplitude-gauge tournament the
length normalization is frozen at

```text
V_1(1,0)=U_1(1,0)+c_l=0,
c_l=-U_1(1,0).
```

All new experiments require configuration schema 4 and
`scalingContract='exact_gauge_no_feedback_v1'`. The four frozen scaling
sentinels are

```text
lengthScaleGain / widthGaugeGain / omegaGaugeGain / travelingWaveGain
                    1 / 0 / 0 / 0.
```

The retired fields `maxDynamicRate`, `adaptiveGain`, `adaptiveLengthScaling`,
`widthExpansionStrength`, `widthContractionOnset`, `widthContractionStrength`, and
`maxWidthRateCorrection` must be absent, not merely assigned harmless-looking
values. Canonical `tau` is the integration clock and every recorded
`timeSpeed` must equal one. CFL selection and the mesh/trust/finite-state
gates remain active; they constrain numerical accuracy and are not modulation
feedback.

### Exact moving quadratic-peak amplitude gauge

The large-box candidate `wall_omega_quadratic_peak` tracks the active interior
maximum of the wall trace `R_X(X,0)`.  On the possibly nonuniform active
three-node stencil, let `Q` denote evaluation at the strictly interior vertex
of the unique interpolating quadratic.  The instantaneous algebraic rule is

```text
P = Q[R_X](X_vertex),
F = Q[B_X](X_vertex),
c_omega = -F/P,
P_tau = F + c_omega P = 0.
```

Here `B` is the base density right-hand side before the amplitude-gauge term.
The same vertex and interpolation weights are used for `P` and `F`; concavity,
vertex interiority, curvature/weight conditions, active-node margins, stencil
identity, and the normalized rate residual are retained in history.  Thus this
is a local discrete phase condition, including on nonuniform grids, rather than
a peak controller.  The length rule remains the independent transport anchor

```text
c_l = -U_1(1,0),  V_1(1,0) = 0.
```

Both rates are recomputed from the current state.  `tau` is integrated
directly with `timeSpeed=1`; there is no dynamic limiter, smoothing gain, or
width feedback.  The `none` amplitude gauge instead prescribes exactly
`c_omega=0` and `C_omega=1`.  Its quadratic rate is counterfactual telemetry,
so a constant selected zero may never be reported as empirical rate
convergence.

## Registered candidates

`ipm_gaugelab_candidate_specifications` is the single tournament whitelist.
The preregistered group consists of two incumbents and three new candidates:

- `gradient_energy` preserves the global gradient energy.
- `anchor_wall_window_l2` uses the positive wall functional
  \(Q_2=(\langle R_X^2\rangle_w)^{1/2}\).
- `anchor_wall_template_projection` freezes the initial normalized wall-slope
  template \(\phi_0\), sets \(P=\langle R_X,\phi_0\rangle_w\), and chooses
  \(c_\omega\) so that \(P_\tau=0\). Its tabulated template is defined on the
  initial mesh, so this candidate is explicitly fixed-mesh only.
- `anchor_bulk_gradient_l2` uses
  \(E=\langle |\nabla R|^2\rangle_W\), \(Q=E^{1/2}\), and
  \(c_\omega=-\langle\nabla R,\nabla B\rangle_W/E\).
- `anchor_wall_window_l4` uses
  \(M_4=\langle |R_X|^4\rangle_w\), \(Q_4=M_4^{1/4}\), and
  \(c_\omega=-\langle R_X^3B_X\rangle_w/M_4\).

The compact wall support is fixed before all runs at center `X=1`, radius
`0.75`, hence open support `(0.25,1.75)`. The bulk bump is fixed at
`(X,Y)=(1,0)`, radius `0.75`, with half-disk support inside
`(0.25,1.75) x [0,0.75)`. The center, radius, and positive-support point
counts are recorded at every trusted sample and must remain constant.

`anchor_wall_slope` and `wall_omega_peak` remain registered only as historical
comparators. `anchor_wall_strain` labels old diagnostic artifacts but is rejected
by the schema-4 runtime because its former discrete relation was not exact.
`none` is an optional negative control and is never eligible to win.

## Hard acceptance gates

Every candidate starts from the same physical initial datum on the frozen
coordinated F80/R20 q256 mesh (`513 x 257` stored nodes). That mesh must obey:

```text
maximum adjacent-cell ratio       <= 1.08
maximum log-spacing curvature     <= 0.01
minimum seven-point stencil rcond >= 1e-9
minimum global quadrature ratio   >= 1e-8
local quadrature/control ratio    in [0.35,1.65]
```

All localized gauges fail closed unless their full trusted histories contain
the fixed geometry and support diagnostics, a condition of at least `0.25`,
a constant positive reference, a positive same-sign current value, and a
normalized rate residual at most `1e-10`. The normalization is
`abs(RateResidual)/RateScale`; the raw residual is reported but is not compared
to a dimensionful universal threshold. The generic `omegaGaugeResidual` must
also be at most `1e-10`.

The wall-moment identity

```text
WindowValue^WindowMomentOrder = WindowMoment
```

must hold to relative error `1e-10`; the recorded order must equal two or four
as registered. For the L2 wall gauge, `WindowEnergy=WindowMoment` is checked
as well. The bulk identity `BulkValue^2=BulkEnergy` has the same tolerance.
Missing fields, nonpositive scales, sign changes, changing support counts, or
identity failures eliminate the candidate.

The common gates are unchanged: continuous trusted prefix through the target,
fourth-order spatial/time/transfer tuple, fixed mesh with zero remesh count,
at least 14 strict core cells per direction, and terminal safety strictly below
`0.70`.  In addition, the entire trusted history must itself certify the exact
length gauge and absence of feedback: `lengthScaleGain=1`,
`widthControlMode=0`, both width-correction histories are zero,
`c_l=c_lNominal=canonicalCLNominal`, the recorded anchor remains `X=1`, and
both the length-gauge residual and transported anchor velocity satisfy the
registered normalized `1e-10` tolerance.  Thus the artifact, rather than only
the input configuration, proves

```text
V_1(1,0)=0.
```

The physical covariance gate is
computed for every successful pair in both directions and reduced to the
absolute worst value. It checks common physical time, laboratory peak shift,
physical fields, wall peak/width, and physical `max|rho_x1|`. At least two
successful production candidates are required; no display baseline influences
eligibility or rank.

The ranking tail is selected by physical time (last 40%, capped at `0.25`),
while canonical tau remains the regression coordinate. Ranking prioritizes
the worst relative drift of `c_l,c_omega`, then detrended oscillation,
`c_l-c_omega` drift, physical covariance, and cost. A directly normalized
observable is never treated as independent shape evidence.

## Preregistered q256 group

The only mutable option is console verbosity:

```matlab
addpath('research/experiments','research/analysis','research/gauge_lab');
products = ipm_gaugelab_run_q256_group_stage_v4(struct('verbose',true));
```

The protocol is saved before any solve.  It freezes the complete candidate
specifications, mesh gates, physical-covariance gates, checkpoint cadence and
launch budget, and the returned tournament must match those structures with
`isequaln`.  A saved tournament result embeds its own candidate specification,
which is used when it is reassessed; later registry edits cannot silently
reinterpret an immutable artifact.  It runs exactly the five candidates
listed above from the original datum to physical time `0.55`, with
`maxDt=2e-3`, `CFL=0.25`, fixed radius `0.75`, and a two-hour aggregate launch
budget. A conservative forecast prevents starting a candidate that cannot
fit in the remaining budget; a solve already in progress is never killed at
an arbitrary state. Each candidate writes to a fresh timestamped directory,
uses its own result and same-gauge checkpoint, and never overwrites an earlier
artifact.

Designed gauge degeneracy, numerical/grid safety stops, wall-clock boundary,
and unexpected infrastructure errors are recorded as different failure
classes.  When a solve throws after writing a checkpoint, the runner discovers
the latest immutable, suffixed checkpoint in that candidate directory rather
than assuming that the requested base name was written literally.

Tournament results are stored under
`result/verification/q256_comega_gauge_tournaments_exact_v4/`; registration
summaries are stored under
`result/verification/q256_comega_gauge_group_stage_v4/`.

## 2026-09-05 formal q256 result

The preregistered five-candidate group completed in `2337.84 s`, below its
`7200 s` budget. All five runs reached physical time `0.55` in 432 steps,
passed their hard audits, and passed all 10 symmetric pairwise physical-
covariance comparisons. The ranking is lexicographic in the six components
described above; it is not a weighted or post-hoc scalar score:

| rank | `cOmegaGauge` | worst tail drift of `c_l,c_omega` | worst detrended RMS | terminal `(c_l,c_omega)` |
|---:|---|---:|---:|---:|
| 1 | `gradient_energy` | 0.0390410 | 1.56868e-4 | (0.532713, 0.265508) |
| 2 | `anchor_wall_template_projection` | 0.117985 | 3.75635e-3 | (0.471819, 0.0386551) |
| 3 | `anchor_wall_window_l2` | 0.399637 | 5.42665e-3 | (0.469538, 0.0231403) |
| 4 | `anchor_bulk_gradient_l2` | 0.606956 | 9.52427e-3 | (0.499906, 0.0673326) |
| 5 | `anchor_wall_window_l4` | 0.679437 | 5.50459e-3 | (0.472483, 0.0215006) |

Across the full trusted histories, the worst normalized `X=1` transport-anchor
defect was `2.533e-16` and the worst generic omega-gauge residual was
`1.139e-16`; `lengthScaleGain` was identically one and `widthControlMode` was
identically zero. The worst pairwise differences over all candidates were
`1.510e-12` in physical rho relative L2, `3.224e-12` in physical omega
relative L2, and `6.765e-14` in physical `max|rho_x1|`. Thus the candidates
changed gauge coordinates but did not measurably change the physical solution.

The common terminal physical values were `max|rho_x1|=0.901020`,
`||grad rho||_inf=0.945037`, strict core-cell counts `75.96/72.37`, and safety
factor `0.120316`. The common far-boundary velocity and omega ratios were
`0.9322` and `0.03339`, both above the registered warning level. Consequently,
`gradient_energy` is the q256 screening winner, but its roughly 3.9% tail
trend is not evidence of an asymptotic constant, and this fixed-box screen is
not a credible blow-up computation or by itself a q512 launch decision.

The immutable [group summary](../../result/verification/q256_comega_gauge_group_stage_v4/group_20260905T064403198Z/group_stage_summary.mat),
[tournament summary](../../result/verification/q256_comega_gauge_tournaments_exact_v4/tournament_20260905T064403596Z/tournament_summary.mat),
and [rate plot](../../result/verification/q256_comega_gauge_tournaments_exact_v4/tournament_20260905T064403596Z/q256_gauge_tournament_rates.png)
contain the complete registered protocol and per-candidate evidence.

The checkpoint reader also rejects `maxDynamicRate` if it is injected into a
purported schema-4 runtime-rescaling state. This closes the artifact boundary:
the removed limiter cannot be revived by configuration, result history, or
checkpoint restore.

## Preregistered q256 length-gauge comparison

`ipm_gaugelab_run_q256_cl_tournament` is a separate comparison in which
`cOmegaGauge='gradient_energy'` is frozen and the only protected-configuration
difference is `scaling.lengthGauge`.  It does not reuse the amplitude
tournament's conclusion as evidence that `c_l` is optimal.  The registry
`ipm_gaugelab_length_candidate_specifications` contains the six exact
schema-4 algebraic rate rules already implemented by
`ipm.evolve.isotropicGauge`:

| `lengthGauge` | exact algebraic rule | preregistered feature parameter |
|---|---|---|
| `transport_anchor` | `c_l=-U1(X_anchor,0)/X_anchor` | `X_anchor=1` |
| `symmetry_peak` | `c_l=-(U1(X_peak,0)-U1(0,0))/X_peak` | peak tracking factor `3` |
| `local_strain` | `c_l=-d_X U1(X_peak,0)` | peak tracking factor `3` |
| `omega_peak_location` | `c_l=physicalRhs_XX/(X_peak R_XXX)` | peak tracking factor `3` |
| `wall_omega_width` | `c_l=-physical_width_rate/current_width` | connected `0.9` level from `[0.1,0.5,0.9]` |
| `wall_density_width` | velocity difference divided by crossing separation | normalized-density levels `[0.2,0.6]` |

All six candidates share the F80/R20 q256 mesh, the fourth-order numerical
tuple, one common CFL/`maxDt`, fixed-mesh policy, and the four no-feedback
sentinels.  Their defaults are `CFL=0.25` and `maxDt=2e-3`; the caller may
override them jointly for the whole tournament with `cfl<=0.5` and positive
`maxDt`.  Both values are frozen in the protocol and cannot differ between
candidates.  The full resolved config for every candidate and the exact
formula string are saved before any solve.  The q256 axis misses `X=1` by `0.001875`,
so the transport-anchor entry explicitly retains the existing linear
interpolation error; the comparison does not conceal it or move the anchor.
Gauge-contract metadata schema 2 also records the exact evaluation functional
and its registered parameters.  In particular, both width gauges identify
piecewise-linear crossings and their selected levels; the density-width path
fails closed on an unbracketed, reversed, or unresolved crossing and records
both crossing slopes, brackets, and conditions.

`omega_peak_location` remains an algebraic screening candidate rather than an
exact discrete peak-phase condition: its `X_peak` is a quadratic peak of
`R_X`, while `R_XX` and `R_XXX` are differentiated fields sampled by linear
interpolation.  Every record therefore includes normalized `R_XX(X_peak)`, the
global `X R_XXX` response condition, and the full phase residual including
length, translation, and amplitude terms.  Its quotient residual alone is not
accepted as evidence that the phase defect vanished.

A bare call is bounded to physical time `1e-4` and at most four steps.  A
candidate subset is allowed only in smoke mode, which is useful for testing
the runner and the bidirectional comparison path:

```matlab
addpath('research/experiments','research/gauge_lab');
smoke = ipm_gaugelab_run_q256_cl_tournament();
```

The six-candidate `t=0.55` formal run requires two independent explicit
choices; `preflightOnly=true` registers and audits the same formal protocol
without solving:

```matlab
preflight = ipm_gaugelab_run_q256_cl_tournament(struct( ...
    'runMode','formal','preflightOnly',true));

formal = ipm_gaugelab_run_q256_cl_tournament(struct( ...
    'runMode','formal','confirmLongRun',true, ...
    'cfl',0.5,'maxDt',4e-3));
```

Each completed result records elapsed time, accepted steps, continuous trusted
prefix, raw and normalized exact-gauge residuals, and the maintained
`gaugeRateConvergence` diagnostic.  The latter reports canonical-time weighted
statistics on the fixed trailing 25%, 40%, and 60% windows for `c_l`,
`c_omega`, `c_omega/c_l`, and `(c_l-c_omega)/c_l`.  Its 1% flag remains
diagnostic only and requires the tail trend, nested-window mean spread, and
detrended relative RMS all to remain below 1%.  Every pair of successful
candidates is compared in both
directions at the common physical endpoint on the inner physical window; at
least two successful candidates are required.  A designed degeneracy is
listed in `comparisonExcludedCandidates` and does not erase valid pairwise
evidence among the remaining candidates.

Formal comparison eligibility is anchored to the preregistered incumbent
`transport_anchor`: an alternative must pass its bidirectional physical-field
comparison against that reference, and every pair inside the resulting
eligible subset must also pass.  A noncovariant outlier is reported and
excluded instead of invalidating all otherwise consistent candidates.  The
summary separately retains the all-pairs flag, the eligible-subset flag, and
both candidate lists; no failed comparison is erased.

The physical wall-width comparison is recomputed from every terminal physical
omega field at one common half-maximum level.  The generic result comparator's
tracked width cannot be used unchanged here because `wall_omega_width` reports
its registered 0.9 gauge level while the other runs report the diagnostic 0.5
level; mixing those definitions produces a false covariance failure.

Formal ranking is lexicographic in worst `c_l/c_omega` tail trend, nested-window
mean spread, detrended RMS, bidirectional physical covariance normalized by
each preregistered hard limit, and wall time. Undefined relative statistics
are ineligible for ranking.
A winner is therefore only a q256 fixed-box screening choice.  The local
derivative and peak-curvature rules can be resolution sensitive, while the two
width rules rely on interpolated crossings.  None of these normalizations
provides independent shape evidence, and no formal c_l tournament result has
yet been run or recorded here.  Box, time-step, and grid refinement remain
required before using a selected rule in a long q512 campaign.
Automatic q512 promotion additionally requires the winner itself to pass the
three rate-stability gates, every successful pair to pass physical covariance,
and the exclusions to be reported explicitly.  If no candidate satisfies
those conditions, the table may name a conditional screening leader but the
q512 launcher must retain its prior frozen rule.

The 2026-09-05 implementation check ran the formal preflight and the complete
default smoke, not the formal tournament.  All six smoke candidates reached
the common `t=1e-4` endpoint with exact-contract and continuous-trust checks
passing; all 15 bidirectional pairs passed.  The worst inner-window rho/omega
relative L2 differences were `3.884e-6/7.958e-6`, the worst physical
`max|rho_x|` relative difference was `8.016e-7`, and the common half-maximum
width difference was `3.850e-6`.  The smoke took `21.33 s`.  Its two-step
histories are intentionally too short for the 16-point rate-convergence
diagnostic, so it produces no ranking or stability claim.

## 2026-09-05 q512 quadratic-gauge and CFL freeze

`ipm_gaugelab_screen_q512_comega_box` performed the initial-RHS-only box
screen on the anchor-aligned q512 factory at half-widths `281`, `1e3`, `1e4`,
`1e5`, and `1e6`, without taking a PDE step.
On the final `1e5 -> 1e6` comparison, global gradient energy grew from
`1.79572e5` to `1.79579e6` and its forcing from `-4.38572e4` to
`-4.38591e5`.  Their ratio happens to give `c_omega=0.244232359244` at
`H=1e6`, but both constituents are extensive for this nondecaying
`degenerate_primitive` datum.  Consequently `gradient_energy` was rejected as
the large-box normalization, not promoted because its quotient looked stable.

The local quadratic functional was box-stable over the same comparison:
`P=0.686073738019`, while `c_l` changed from `0.488459589077` to
`0.488462559170` and `c_omega` from `0.125042864250` to
`0.125042865249`.  All registered curvature, weight, vertex, active-stencil,
runtime-agreement, and exact residual gates passed.  `none` passed only as the
prescribed control.  The formal artifact is the [q512 initial box-screen
summary](../../result/verification/q512_comega_initial_box_screens_v4/screen_20260905T110530211Z_tpbe425904_b6b4_4b3f_9aa4_c02b888cd290/comega_box_screen_summary.mat).
This is a normalization screen, not truncation evidence; in particular, the
`H=1e6` far-velocity ratio was still `0.986819`.

`ipm_gaugelab_run_q512_gauge_cfl_screen` then performed a four-run fresh-root
screen crossing
`wall_omega_quadratic_peak` and `none` with `CFL=0.25` and `0.50` on the same
q512 `H=1e6` grid through physical time `0.03` (`maxDt=0.004`).  All four hard
audits and both same-gauge CFL comparisons passed.  For the quadratic gauge,
`CFL=0.25/0.50` used `40/21` steps and `184.672/106.776 s`; the control used
`185.350/106.231 s`.  The aggregate wall-time speedup was about `1.74`.
The worst same-gauge CFL errors were `8.003e-12` for physical rho relative L2
and `2.963e-9` for physical omega relative L2, with physical-peak relative
difference below `4e-14`.  The worst same-CFL cross-gauge rho/omega errors
were `1.156e-11/2.051e-11`.  Therefore `CFL=0.50` is the frozen long-campaign
choice.  Complete evidence is in the [q512 gauge-by-CFL screen
summary](../../result/verification/q512_fresh_gauge_cfl_short_screens_v4/screen_20260905T114636013Z_tp3e066220_eec8_4c94_a0c7_b02a376fc67f/gauge_cfl_screen_summary.mat).

The q256 quadratic-gauge trajectory was separately extended to physical time
`1.7`.  Its endpoint rates were
`(c_l,c_omega,kappa)=(0.238287,-0.266436,0.504723)`, and the selected physical
quadratic peak was `4.41394`.  Over the final physical-time window of width
`0.125`, the affine fit to its inverse had slope `-0.748374`, versus mean
instantaneous algebraic slope `-0.747521`; their mismatch was `0.114%`, and
the fit RMS normalized by the inverse drop was `0.113%`.  This is a useful
Type-I candidate signature for the inverse peak, not literal linear growth of
the peak itself.  In the same window, however, `c_l/c_omega/kappa` still
changed by `19.4%/14.2%/3.08%`, and the trajectory contained only `1.861`
e-folds, below the preregistered three-e-fold minimum.  The [q256 result](../../result/verification/q256_quadratic_peak_continuations_v4/continuation_20260905T120651140Z/result.mat)
therefore supports neither a blow-up claim nor asymptotically constant rates.

A real two-step q512 quadratic root and two-step continuation also froze the
restart chain.  The [root summary](../../result/verification/q512_fresh_exact_gauge_roots_v4/root_20260905T122642961Z_large_box_smoke_tp5f74d6c0_96d5_4c8a_b18f_88b19ea00f92/root_summary.mat)
hard-passed at two total steps.  The final [continuation
summary](../../result/verification/q512_schema4_continuations/continuation_20260905T124215721Z_tp1b1f9942_ef22_4494_b5f7_ee4ca7b4a77c/continuation_summary.mat)
is `completed_endpoint`, with two segment steps/four total steps,
`27.9553 s` wall time, `postAudit.hardPassed=1`, and
`parentHistoryExactPrefix=1`.  Two earlier development replays were rejected
only because empty `anisotropic` history groups exposed an empty-shape
`fieldnames`/`isfield` prefix-audit bug.  The audit now checks child fields one
name at a time and treats two empty scalar groups as an exact prefix; no field
value mismatch or numerical failure was involved.

After all 202 MATLAB files returned zero Code Analyzer diagnostics and the
12-layer baseline suite passed, the registered long root was launched at
2026-09-05 20:47 CST in detached session `ipm_q512_quad_h1e6_t170`.  It uses
the selected quadratic gauge on `1025 x 513` nodes, initial half-width
`H=1e6`, `CFL=0.50`, `maxDt=0.004`, physical endpoint `t=1.7`, canonical
ceiling `tau=4.5`, and a 4000-step ceiling.  The [preflight
summary](../../result/verification/q512_fresh_exact_gauge_roots_v4/root_20260905T124547741Z_large_box_campaign_tp2922f311_a512_4591_8ba1_898027e4bd92/root_summary.mat)
passed the mesh and exact-contract gates.  The active [root
summary](../../result/verification/q512_fresh_exact_gauge_roots_v4/root_20260905T124741359Z_large_box_campaign_tp359ea53d_96ea_4ca8_b160_17e0fe61f50f/root_summary.mat)
and adjacent manifest are the authoritative status records.  `H=1e6` is the
initial physical half-width; every record also stores the evolving physical
coverage after division by `C_l`.  A running registration is not a numerical
result or a blow-up claim.

## Fresh schema-4 q512 roots and the box ladder

`ipm_gaugelab_run_q512_root` is the only registered fresh-root launcher in
this directory. It never accepts a result or checkpoint input and therefore
cannot silently continue the historical schema-2 q512 campaign. Every root
starts from the analytic datum at `t=tau=0`, freezes the mesh, and fixes

```text
lengthGauge / transportAnchorX / cOmegaGauge
transport_anchor / 1              / gradient_energy (default), or
                                    wall_omega_quadratic_peak
```

together with `exact_gauge_no_feedback_v1` and the
`high_order/weno5_fd/ssprk54/high_order` tuple. Checkpointing is enabled with
`atExit=true`, and the checkpoint must validate as schema 4 with no retired
field in either its configuration or runtime rescaling state.

Add the experiment and gauge directories to the MATLAB path. A read-only
registration of the matched q512 refinement is:

```matlab
addpath('research/experiments','research/gauge_lab');
root = ipm_gaugelab_run_q512_root(struct('preflightOnly',true));
```

The default executed call is a bounded smoke test: physical time `1e-4`,
canonical time at most `1e-3`, and at most two steps. Smoke mode rejects a
request beyond physical time `1e-3`, canonical time `1e-2`, or 10 steps, so a
bare call cannot accidentally launch the multi-hour calculation:

```matlab
root = ipm_gaugelab_run_q512_root();
```

The matched mode uses the strict two-times refinement of the tournament mesh,
with `1025 x 513` stored nodes on `[-281,281] x [0,132]`. Large-box mode uses
the anchor-aligned F80/R20 q512 factory and defaults to the initial physical
box `[-1e6,1e6] x [0,1e6]`:

```matlab
root = ipm_gaugelab_run_q512_root(struct( ...
    'boxMode','large_box','preflightOnly',true));
```

`boxHalfWidth` selects another large-box rung and applies the same half-width
to x and y. Each rung is a separate fresh root, not a regrid or restart. A
typical truncation study registers, then smokes, `1e4`, `1e5`, and `1e6`
before spending on a long segment. The value is the initial physical extent;
the terminal physical extent is the saved computational axis divided by
`C_l` and must be read from the result rather than inferred from the rung
label.

A long segment requires both an explicit mode and an explicit acknowledgement:

```matlab
root = ipm_gaugelab_run_q512_root(struct( ...
    'boxMode','large_box','boxHalfWidth',1e6, ...
    'runMode','campaign','confirmLongRun',true, ...
    'cOmegaGauge','wall_omega_quadratic_peak', ...
    'cfl',0.50));
```

Campaign mode defaults to a canonical discovery ceiling of `tau=12`, a
physical safety ceiling of `t=4`, and `maxSteps=30000`.  The earlier `t=2.1`
ceiling was removed because the current q256 tail extrapolation places it near
`tau=3.05`, before three e-folds of gradient amplification.  That extrapolation
selects a minimally informative window only; it is not asymptotic evidence.

`cfl`, `maxDt`, horizons, record cadence, and checkpoint cadence are recorded
but do not enter either gauge formula. A faster CFL choice must be run as a
separate registered root and compared at common physical time; it is not a
modulation parameter.

Before solving, the launcher saves `root_summary.mat` and appends
`root_manifest.jsonl` under a unique directory in
`result/verification/q512_fresh_exact_gauge_roots_v4/`. The preregistration
contains the exact formulas, full frozen configuration, requested horizons,
and mesh audit. Both axes must satisfy the same hard gates used by the q256
laboratory: adjacent ratio `1.08`, log-spacing curvature `0.01`, seven-point
stencil rcond `1e-9`, global quadrature ratio `1e-8`, and local
quadrature/control ratio `[0.35,1.65]`.

After an executed run, the summary records total wall seconds, accepted steps,
and wall seconds per step. It rejects a result unless the terminal record is in
the continuous trusted prefix, the stop reason is the requested physical or
canonical endpoint (with the solver's own `minDt` stopping tolerance),
`timeSpeed` is identically one, direct and
canonical rates agree, both exact gauge residuals are at most `1e-10`, feedback
histories remain zero, the mesh never changes, and the terminal schema-4
checkpoint matches the result. These are artifact-integrity gates, not evidence
that `c_l`, `c_omega`, or the physical maximum has reached an asymptotic law.

In particular, `gradient_energy` is a global normalization for the uncut
degenerate primitive and can depend on the truncation box. Claims that its
rates converge, that physical `max|rho_x|` grows exponentially in canonical
time, or that its inverse is affine toward a physical-time endpoint therefore
require agreement across
the fresh box ladder and a matched
resolution/time-step comparison. Raw far-boundary velocity remains telemetry
for this nondecaying datum; source-tail and matched-box field comparisons carry
the truncation decision.

### Registered schema-4 q512 continuation

`ipm_gaugelab_run_q512_continuation` is the only registered continuation
launcher in this laboratory. Its first argument must pass
`ipm.output.readCheckpoint` as a trusted schema-4 accepted-step checkpoint and
must retain the fresh q512 root lineage. It additionally requires the fixed
`1025 x 513` mesh, zero remesh count, `transport_anchor` at `X=1`,
the root's selected `gradient_energy` or `wall_omega_quadratic_peak` rule,
`exact_gauge_no_feedback_v1`, and the
`high_order/weno5_fd/ssprk54/high_order` tuple. Historical schema-1--3
checkpoints, retired feedback fields, a different gauge, or an unregistered
q512 state are rejected rather than upgraded. The large-box factory retains
its anchor-aligned `X=1` node; the matched-box factory uses the registered
fixed-coordinate linear interpolation and is not required to contain that
node exactly.

The default is a no-PDE preflight. `finalTime`, `physicalFinalTime`, and
`maxSteps` may remain unchanged or increase, but may never decrease, and all
three limits must leave the stored state open. CFL, `maxDt`, `minDt`, output
cadence, checkpoint cadence, grid, physical data, diagnostics, and every
numerical choice are copied exactly from the checkpoint. The only other public
changes are result/checkpoint/video paths and the explicit plot, live-plot,
video, and verbosity controls. In particular, this runner has no `cfl`,
`maxDt`, `lengthGauge`, `cOmegaGauge`, regrid, or checkpoint-cadence option.

```matlab
preflight = ipm_gaugelab_run_q512_continuation( ...
    root.checkpointFile,struct( ...
    'parentRootSummary',root.summaryFile, ...
    'physicalFinalTime',2e-4,'finalTime',2e-3,'maxSteps',4));
```

Set `execute=true` only after inspecting the saved
`continuation_summary.mat`. Evolution then occurs solely through
`ipm.solve(overrides,checkpoint)`. The summary stores the complete validated
parent checkpoint signature, its step and clocks, a compact persisted restore
chain, the exact frozen configuration audit, the parent-history prefix check,
the child checkpoint signature, gauge/rate/growth diagnostics, and wall time
per newly accepted step. `parentRootSummary` is required provenance.  It must
report a completed, hard-passed direct parent (fresh root or continuation),
retain the same frozen gauge and restore chain, and reference a checkpoint
with the identical signature.

A terminal physical/canonical endpoint is classified as
`completed_endpoint`. Reaching the requested cumulative `maxSteps` while both
time horizons remain genuinely open is a valid checkpointed campaign segment
and is classified as `completed_intermediate`; it is not reported as a final
endpoint or asymptotic result. Any other stop reason, a discontinuous trusted
prefix, changed protected configuration, lost parent-history prefix, or child
checkpoint mismatch fails the post-run gate.

A bounded large-box q512 split/resume exercise is:

```matlab
addpath('research/experiments','research/gauge_lab');
root = ipm_gaugelab_run_q512_root(struct( ...
    'boxMode','large_box','boxHalfWidth',1e6, ...
    'cOmegaGauge','wall_omega_quadratic_peak', ...
    'cfl',0.50,'verbose',false));
continued = ipm_gaugelab_run_q512_continuation( ...
    root.checkpointFile,struct('parentRootSummary',root.summaryFile, ...
    'execute',true,'physicalFinalTime',2e-4, ...
    'finalTime',2e-3,'maxSteps',4,'verbose',false));
```

Run the two calls serially. They add at most a few accepted q512 steps and are
only a checkpoint integrity exercise, not a CFL, truncation, late-time, or
blow-up validation. The no-PDE rejection harness
`ipm_gaugelab_test_q512_continuation` covers legacy checkpoint schemas and
attempts to change protected CFL, time-step, grid, physical, gauge, and
checkpoint-cadence fields.

### Box-ladder orchestrator

`ipm_gaugelab_run_q512_box_ladder` applies the root launcher independently to
the matched reference and, by default, an anchor-aligned large-factory control
at half-width `281` followed by `1e3`, `1e4`, `1e5`, and `1e6`. Its default is
preflight-only, so the following performs six mesh and
contract registrations but no PDE evolution:

```matlab
addpath('research/experiments','research/gauge_lab');
ladder = ipm_gaugelab_run_q512_box_ladder();
```

The large-box rungs can be changed with an increasing vector, for example
`struct('boxHalfWidths',[1e4,1e5,1e6])`. The matched `[-281,281] x [0,132]`
reference and the large-factory half-width-281 control are always included.
There is no restart, checkpoint-input, or long-
campaign option: each rung must pass through `ipm_gaugelab_run_q512_root` as a
fresh schema-4 transaction.

Executing all bounded smokes requires the explicit switch:

```matlab
ladder = ipm_gaugelab_run_q512_box_ladder(struct( ...
    'executeSmoke',true,'cfl',0.25));
```

The orchestrator inherits the root smoke caps (`physicalFinalTime<=1e-3`,
`canonicalFinalTime<=1e-2`, and `maxSteps<=10`). It writes a ladder summary and
append-only manifest before invoking any rung, then records each root's own
summary, result, and checkpoint paths. An executed rung contributes:

- initial and terminal direct/canonical `c_l`, `c_omega`, and `kappa`;
- gradient-energy X/Y parts, forcing X/Y parts, signed and absolute forcing,
  the signed/absolute forcing cancellation ratio,
  `forcing+c_omega*energy`, and its normalized algebraic residual;
- initial and actual terminal physical-domain endpoints;
- initial and terminal far-source and far-velocity values and ratios;
- physical `max|rho_x|`, accepted steps, wall seconds, and seconds per step;
- the complete tail-rate and canonical/physical maximum-growth diagnostics
  (short smoke histories are expected to report these as insufficient).

For every accepted smoke result, the terminal physical field is compared with
the matched result on the intersection of their physical domains and the
registered cap `[-2,2] x [0,2]`. Interpolation is performed in both directions;
the summary retains both directional records and their componentwise absolute
worst. Results must share a physical endpoint within the registered roundoff
tolerance before this comparison is available.

The same bidirectional comparison is also recorded for successive members of
the anchor-aligned large-box family (`281 -> 1e3 -> ... -> 1e6`). Comparing
matched 281 against large-factory 281 is a combined
factory/domain-geometry/core-phase control (their y extents differ);
successive large-family records expose the remaining box trend. These are not
pure truncation comparisons, because a fixed number of stored nodes requires
the outer-node distribution to change as the box grows.

The ladder applies a frozen screening gate to the final large-family pair that
ends at half-width `1e6` (by default `1e5 -> 1e6`). It checks the bidirectional
core-field errors, terminal `max|rho_x|`, `c_l`, `c_omega`, gradient energy and
forcing, forcing-cancellation ratio, exact `X=1` nodes on every large-box
factory rung, root post-audits,
and the terminal far-source/far-velocity ratios. The field and rate limits
match the short CFL screen in scale; both far-boundary ratios must be at most
`1e-2`. `ladder.hardPassed` is true only when all roots, comparisons, and this
campaign-box gate pass.

This remains a one-step screening gate, not a certificate of truncation
convergence: `gradient_energy` integrates over the whole box and a fixed number
of stored nodes redistributes the outer grid as the box grows. The complete
table is retained so the later campaign can distinguish initial box safety
from genuine late-time box and grid refinement evidence.

### q512 CFL-ladder orchestrator

`ipm_gaugelab_run_q512_cfl_ladder` preregisters two independent q512 roots at
`cfl=0.25` and `cfl=0.50`. Both use the same box, horizon, maximum time step,
fixed grid, exact gauge contract, and numerical tuple. The default box is the
anchor-aligned large-box factory at half-width `1e6`. The candidate list is
fixed rather than caller-selectable, and the runner accepts no checkpoint or
restart input.

The default call performs both root preflights but no PDE evolution:

```matlab
addpath('research/experiments','research/gauge_lab');
ladder = ipm_gaugelab_run_q512_cfl_ladder();
```

Execution requires the explicit switch below. The orchestrator then passes
`runMode='campaign'` and `confirmLongRun=true` to each fresh-root launcher. Its
own caps limit the experiment to physical time `0.03`, canonical time `0.05`,
`maxDt<=0.004`, and at most 1024 steps; the registered defaults use 512
steps, output cadence `1e-10`, and checkpoint cadence `0.01`. The output
cadence equals the current solver `minDt`, so every accepted step is retained
for timestep and rate statistics.

```matlab
ladder = ipm_gaugelab_run_q512_cfl_ladder(struct('execute',true));
```

Each executed root records its terminal physical/canonical clocks,
`max|rho_x|`, `c_l`, `c_omega`, their two dimensionless ratios, accepted steps,
wall time, seconds per step, and distributions of the accepted time step,
realized CFL, and active time-step bound (CFL or `maxDt`, never a gauge-rate
limiter). The summary verifies that those records cover
every accepted solver step. It also retains the complete
`gaugeRateConvergence`, physical maximum-growth, and canonical maximum-growth
diagnostics.

The terminal fields are compared on the intersection of both physical domains
and the registered cap `[-2,2] x [0,2]`. Interpolation is run in both
directions, and strict preregistered limits are applied to the componentwise
worst field errors, the common physical endpoint, terminal rates and rate
ratios, and rate-window agreement when both rate diagnostics are valid.
Performance statistics are measurements, not numerical pass gates.

This is deliberately labeled `short_window_cfl_consistency_only_not_late_time_evidence`.
A pass cannot establish that either gauge rate is asymptotically constant or
that `max|rho_x|` follows a late-time growth law; those claims still require
independent box, time-window, and spatial-refinement agreement.

## Historical artifacts

Schema-3 tournaments and the old two-candidate playoff remain useful for
reconstructing the development history, but they are now unconditionally
`read_only_history`. Stored legacy `eligible=true` flags do not confer schema-4
production eligibility, and schema-3 checkpoints cannot be resumed. Calling
`ipm_gaugelab_run_q256_playoffs` now stops before creating output and directs
new work to the schema-4 group runner.

For context only, the earlier limiter-free schema-3 screen at physical time
`0.55` ranked `gradient_energy`, `anchor_wall_slope`,
`anchor_wall_window_l2`, then `wall_omega_peak`. All reconstructed nearly the
same terminal physical field, but the common far-velocity ratio was about
`0.932`; this laboratory result was never a credible blow-up claim.
