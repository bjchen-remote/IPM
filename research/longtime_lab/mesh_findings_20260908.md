# Schema-4 large-box continuation mesh findings

This note concerns only the `wall_omega_quadratic_peak` / `transport_anchor`
root at `H=1e6`, stored nodes `1025 x 513`. It does not inherit the old
schema-3 campaign trajectory or its inferred asymptotics.

## Accepted source and first continuation interval

The native source is
`result/verification/q512_fresh_exact_gauge_roots_v4/root_20260905T124741359Z_large_box_campaign_tp359ea53d_96ea_4ca8_b160_17e0fe61f50f/checkpoint_20260905T204759339_dynamic_isotropic_1025x513_tp8c8b7442_a9d7_4682_b700_7cf8af6fe3f7_step0000001818.mat`.

It has `tau=3.6732019386250148`, physical `t=1.6999999999998241`,
`max|physical rho_x|=4.5151615670278424`, 920 continuously trusted records,
zero previous regrids, x/y strict 0.9 cores `21.55981491 / 20.38194189`,
and safety factor `0.39250431`.

Over the preceding `Delta tau=0.35`, fitted core-count slopes were
`-15.7348 / -13.2022` cells per tau. These are short-window forecasts,
not convergence laws. They justify an initial `Delta tau=0.2` segment and
fresh measurement, rather than assuming the fixed grid remains resolved
through tau 5. Maintain native checkpoint cadence at most 0.2 and a
conservative terminal gate of both cores >=14 and safety <0.70. Plan a new
grid before either core falls below about 18. A solver's much lower
`resolutionFailed` hard stop is not the conservative campaign gate.

## Measured design defects and isolated remedies

1. Both analytic seed builders default to `minimumWidthFraction=1e-6` of
   the far endpoint. At H=1e6 this imposes a width of 1, while the measured
   strict cores are only `0.04716210 / 0.01393320`. The isolated design uses
   `minimumWidthFraction=1e-10` and retains the existing minimum of two
   local cell widths. No PDE formula or production default changed.
2. The automatic y analytic candidate compressed the current 20.38-cell
   core to 44.84 cells although the target was 24. A minimum-amplitude
   original-reference y equalizer attained 24 cells with minimum dy
   `0.00051205`, versus `0.00026708` for that analytic candidate, and
   retained more global quadrature margin. This is useful mesh economy;
   its effect on actual time-step cost still requires PDE measurement.
3. Automatic x equalization anchors the intersection of core/front
   intervals. For this exact-gauge run, the independent candidate instead
   anchors at x=1, where `transport_anchor` samples the velocity. The
   equalizer's existing root tolerance scales as `1024*eps(H)` and admitted
   a `5.30e-8` coordinate drift. The independent candidate records that
   drift, sets the original reference nodes exactly to x=+1 and x=-1, and
   repeats the complete mesh and frozen-field pair scoring after that edit.
   No mesh gate is relaxed.

## Executed frozen-field checks

`mesh_design_stage` and `mesh_balance_y_candidate(...,1)` generated
`result/longtime/20260908_mesh_design_initial_v1/mesh_anchor1_exact_candidate.mat`.
The variable `candidate` contains `candidateX`, `candidateY`, the actual
reference provenance, equalization records, and the complete pair score.
The original unregridded schema-4 root result supplies the immutable global
reference axes; earlier H=281 datasets fail the required endpoint match.

The exact-anchor balanced pair passed full-q512 frozen transfer and mesh
checks: strict x/y cores 24/24, relative spatial rho-x peak jump
`1.21912929e-4`, round-trip rho relative L2 `4.53458378e-6`, and y minimum
quadrature-weight/mean-weight `8.38008293e-8`. The x=1 node is exact.
The original automatic design was feasible in 4.35 seconds. All three
research MATLAB scripts checked so far have zero Code Analyzer findings
after removal of an obsolete suppression.

All candidates retain adjacent ratio <=1.08, log-spacing curvature <=0.01,
seven-point stencil rcond >=1e-9, positive quadrature with global ratio
>=1e-8, and local q/control-width ratios within [0.35,1.65]. The old
global-ratio warning 1e-4 is already exceeded by the original H=1e6 grid;
it must not be substituted for the explicit 1e-8 compatibility floor.

Frozen transfer is not a PDE validation. `mesh_test_candidate` performs an
audited native regrid, then advances the old mesh by Delta tau=0.01 and
advances the new mesh to exactly the old run's physical end time. It reports
laboratory-frame local rho, omega and rho-y errors, a separate peak-aligned
comparison, the physical-clock match, both gauge-rate changes, spatial rho-x
change, and trusted status. The canonical endpoints may differ slightly
because the physical end time is matched. Run artifacts, not this plan,
determine whether that test passes.

The initial transaction/short-run test was started under MATLAB's original
10-thread setting, then explicitly interrupted to relieve system-wide memory
pressure while other q512 work was running. It exited with status 130 before
writing a completed transaction or either PDE result. Its directory is
`result/longtime/20260908_mesh_regrid_test_v1` and its launch log is the
adjacent `.log`. The frozen checks above passed; the transaction and PDE
comparison remain unverified and must run serially on the latest accepted
state. Do not report the interrupted test as a numerical rejection or pass.

## Continuation integration and limits

The maintained `ipm_gridlab_run_q512_campaign` does accept schema 4 and
provides useful transactional jump, core, trust and quadrature gates.
Its old default seed width fractions and original-reference dataset must
be adapted explicitly for the H=1e6 root. Its default algorithm also regrids
at every segment boundary; same-grid segments remain preferable while the
core margin is sufficient. Both its design and regrid preserve node counts
and computational endpoints. Increasing q or changing H needs a separate
audited branch rather than silently editing a checkpoint.

The core `restoreCheckpoint` accepts a horizon equal to the stored time;
therefore zero-time regrid at the physical endpoint does not require changing
the accepted physical time. A following `ipm.solve` opens the horizon using
legal overrides. The more restrictive registered fixed-grid gauge-lab
wrapper additionally refuses to lower the previously requested finalTime;
direct `ipm.solve` can legitimately use smaller intermediate segment horizons
as long as they remain later than the accepted state.

Every later candidate must be rebuilt/scored on its actual source state.
`mesh_test_candidate` checks exact equality of the source rho and axes with
the frozen design before any transaction. Long-time follow-up should retain
the per-regrid rho-x jump <=0.002 and cumulative jump <=0.02 budgets, and
continue matched-time grid, time-step and far-box comparisons. Reaching a
large gradient remains evidence to examine, not a reason to label a
singularity before those comparisons.
