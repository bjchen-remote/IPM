# Longer fixed-inner-node three-box protocol

The first q512 interval .001 passed all five registered native/full/crop
cases. Its actual results are in [mesh_fresh_box_protocol_plan.md](mesh_fresh_box_protocol_plan.md).
This next protocol keeps the same step2031 source and advances to total
local physical elapsed time **.01**, absolute time 1.7898569844190184.
It preserves every initial inner node and independently compares full,
fraction .1 and fraction .01 boxes under natural and common WENO epsilon.
No q512 run was started while preparing this longer protocol.

## Choosing the interval

`mesh_predict_box_window.py` reads only scalar histories with the ctypes HDF5
monitor. It performs no signature or trusted-prefix validation. Its saved
forecast is `20260908_campaign_v1/fresh_box_ladder_prediction_v1.json`.
It uses only the source's remeshCount=1 records and fits log core counts
against canonical time on .08/.12/.2 windows. The maximum measured decay,
safety slope and c-l-minus-c-omega are each multiplied by a declared 1.25
planning factor. This is a heuristic margin, not a rigorous future bound.

For lambda=Cx/Comega and an inflated positive kappa=c-l-c-omega, the planning
map is `Delta tau=-log(1-kappa*lambda*Delta t)/kappa`, obtained by freezing
the rate difference. Core forecasts use exponential decay in that canonical
interval. They do not treat the physical clock as canonical time.

| Total physical elapsed | Predicted canonical interval | Predicted cores x/y | Predicted safety | Planning margin |
| --- | --- | --- | --- | --- |
| .005 | .0452371 | 17.5533 / 17.9602 | .454731 | Pass |
| .01 | .0917806 | 16.8730 / 17.3165 | .471351 | Pass |
| .02 | .1891048 | 15.5345 / 16.0441 | .506104 | x below planning buffer16 |

Select .01 first. The preparation buffer is both predicted cores>=16 and
safety<.65. Actual production-case gates remain both cores>=14 and
safety<.70, a complete trusted history, and a normally validated terminal
native checkpoint exactly paired with the result. Reassess a .02 total
interval from the actual .01 histories before requesting another window.

## Six necessary branches

| Box fraction | Natural epsilon | Common epsilon |
| --- | --- | --- |
| 1, 1025x513 | Continue existing full_natural fresh native CP | Continue existing full_common_epsilon CP |
| .1, 963x449 | New physical IVP from original step2031 | New physical IVP from original step2031 |
| .01, 901x386 | Continue existing crop_natural fresh native CP | Continue existing crop_common_epsilon CP |

The four continued seeds have an existing, verified fresh history from
local physical time0 through .001. They keep their case IDs, original
runtime references, initial paired axes/data and full genuine history;
same-grid `ipm.solve` adds the remaining .009. They are not recreated with
invented history or changed domains. The two middle boxes each receive a
new case ID and new local-time0 history, and integrate .01 directly from the
same original parent late state.

The saved initial physical data can be used for inspection but is not a
trusted checkpoint substitute. Re-reading and validating the original
parent checkpoint before constructing two new IVPs costs no elliptic
factorization and provides stronger provenance than trusting an unsigned
initial-data MAT file. No operator restore is needed to multiply/crop those
initial samples.

The common epsilon remains 1/385^2=6.7465002529937596e-6 in both axes for all
three boxes. It also exceeds the middle box's natural y floor
1/448^2=4.9824617346938776e-6. All six natural/common cases are necessary:
the tiny epsilon effect observed over .001 does not justify dropping a
control at .01. Report full/middle, middle/small and full/small comparisons
within each epsilon family, plus three within-box epsilon-only comparisons.
Global Lax-Friedrichs splitting and outer closures remain box-dependent.

An extra original-parent native continuation would be a seventh case. It is
not part of the minimum ladder: the exact physical covariance and the actual
q512 .001 overlap calibrate the fresh coordinates and rates, while all box
comparisons use the same fresh-coordinate formulation. The report explicitly
sets `nativeCovarianceCheckedAtTarget=false`; it cannot present the prior
.001 test as a measured native/fresh comparison at .01. Add that seventh
sentinel only if a new inconsistency calls for it.

## Entry point and resource registration

```matlab
campaign=fullfile(pwd,'result/longtime/20260908_campaign_v1');
report=mesh_extend_fresh_box_ladder( ...
    fullfile(campaign,'fresh_box_tau428_protocol_v1/protocol_report.json'), ...
    fullfile(campaign,'fresh_box_ladder_prediction_v1.json'), ...
    fullfile(campaign,'fresh_box_tau428_ladder_v1'), ...
    struct('targetElapsedPhysicalTime',.01,'middleFraction',.1, ...
        'maximumAdditionalSteps',128,'allowLarge',true));
```

Prepared launch script: `/private/tmp/ipm_fresh_box_tau428_ladder_v1.m`.
The script must run only after the coordinator allocates a q512 window.
The wrapper checks the exact original source, prior full covariance,
trusted/passed seeds, matching static screens for both crops, sufficient
common epsilon and the registered forecast margin before any PDE.

There are six serial PDE calls, six elliptic assemblies and zero extra
reference-audit restores. All reads and solves use 10 threads. Four calls
resume native fresh checkpoints; two initialize the middle box. Each call
has a cap of128 additional accepted steps, but only a physical endpoint
within1e-11 counts as success. Approximate step counts are 40–50 total per
branch, based on the observed five accepted steps over .001; this is not
a guaranteed timestep schedule.

The short protocol took about54/27 seconds for full/small fresh cases with
five accepted steps, including initialization and saving. Scaling those
measurements for longer integration gives a planning range around20–40
minutes for six cases under similar load; middle-box timing is unmeasured.
This estimate is not a resource bound or benchmark. Peak memory stays near
one full solver plus compact observations: every solve returns without live
operators, and its result is reduced/saved before the next assembly. Native
validation reads do not build an extra LU.

All endpoint diagnostics and short-protocol comparison limits remain
unchanged. Results preserve the physical laboratory window, absolute epoch,
case lineage, full history and seed provenance. A failed branch saves its
available artifacts and stops subsequent work. A passing ladder supports a
longer independent test; it still does not validate the earlier parent
history, establish an infinite-box limit, or change the production domain.

The six-case tiny entry point actually passed with MATLAB exit0 in
`result/longtime/20260908_fresh_box_ladder_tiny_v1`, with an adjacent log.
Four existing65x33/49x25 histories continued from elapsed.02 to.03 without
changing IDs; two57x29 middle-box IVPs started at the same original epoch.
All six endpoint/trust/native-pairing tests passed and both natural/common
sensitive-small-box comparisons failed as intended. Code Analyzer returned
no findings for the new driver and comparison helper. The q512 source
bundle and hashes are frozen in
`20260908_campaign_v1/fresh_box_tau428_ladder_source_v1`.

## Reference and coordinate interpretation

A later no-LU ctypes read of the already native-validated short-run files is
saved in `fresh_box_reference_monitor_v1.json`. That read itself does not
validate signatures. All four full/small fresh branches share pinX
.2111192001667, tracking half-width1.600755013663 and identical safety
references. Transformed original-parent references instead have pinX
.2767186463902 and half-width9.048117939002. Fresh diagnostic reference counts
also differ from those at the original root. Full/small strainTarget differs
by about1.54e-4; in the selected quadratic/anchor gauge it is used only for
reported strainError, not to set the rates. The active peak windows coincide
across fresh boxes, but their difference from the original root is another
reason not to extrapolate measured short native covariance indefinitely.

Fixed inner nodes mean the common initial physical samples and then fixed
computational axes. Dynamic rescaling continues, so the reconstructed
physical axes and physical box boundaries evolve with each branch's Cx/Cy.
The box labels H~1e6/1e5/1e4 refer to original parent rescaled units. The
comparison window is fixed in laboratory physical coordinates; its endpoint
comparisons therefore use the paired current axes and interpolation. This
protocol does not impose stationary physical walls at the initial bounds.

## Actual .01 execution

After the coordinator allocated one q512 window, all six frozen cases
completed with MATLAB exit0. The window was released. Artifacts are in
`20260908_campaign_v1/fresh_box_tau428_ladder_v1`, including
`ladder_report.json`, each branch's result/native checkpoint/paired physical
observations, and an adjacent run log.

All six cases stopped for `physical_final_time`, with absolute-target error
at most2.51e-14, full trusted result and checkpoint histories, exact terminal
rho/axes/step/clocks/case-ID pairing, and preserved real fresh histories.
Minimum terminal cores were17.146948/17.933558 and maximum safety.46655534.

| Natural-epsilon comparison | rho relative L2 | rho-x relative Inf | rho-y relative Inf | Physical gradient peak relative difference | c-l relative difference |
| --- | --- | --- | --- | --- | --- |
| Full / H~1e5 | 4.07622e-7 | 3.55997e-6 | 4.12232e-6 | 1.44629e-7 | 8.09112e-6 |
| H~1e5 / H~1e4 | 4.11804e-6 | 3.59651e-5 | 4.16465e-5 | 1.46114e-6 | 8.17884e-5 |
| Full / H~1e4 | 4.52566e-6 | 3.95250e-5 | 4.57688e-5 | 1.60577e-6 | 8.98788e-5 |

The common-epsilon full/small values were rho L2=4.52546e-6,
rho-x/rho-y Inf=3.95233e-5/4.57668e-5, physical peak difference1.60569e-6,
and c-l difference8.98740e-5. Natural full/small c-omega relative difference
was1.21546e-7. The within-box epsilon-only maximum gradient Inf changes were
3.07e-11/9.09e-9/2.04e-9 for full/middle/small. These controls remain much
smaller than their corresponding box differences at .01; this is measured
at this interval, not inferred from the shorter run.

Observed solve wall times were205.4/152.3/109.8/193.2/153.2/107.4 seconds,
total921.35 seconds (15.36 minutes), excluding driver overhead outside the
per-case timer. Continued full/small cases each added36 accepted steps to
their prior5; newly initialized middle cases used40 steps. This RK
partition difference includes the earlier .001 endpoint clipping. The
physical endpoints match, but temporal refinement is still required to
separate time-discretization effects from box effects.

The supplemental no-LU `fresh_box_ladder_record_audit_v1.json` examined all
saved scalar records after native validation: minimum recorded cores match
the terminal minima, maximum recorded safety is.46655534, and maximum saved
Poisson residual is4.88e-10. The ctypes read itself performs no signature
or trusted-mask validation. The residual threshold1e-8 is a post-run
diagnostic, not a hidden runtime trustedMask gate; these extrema concern
saved records rather than every unrecorded RK substage.

Both natural/common three-box comparisons passed the registered finite-
interval thresholds. The report explicitly leaves
`nativeCovarianceCheckedAtTarget=false`,
`dynamicBoxConvergenceEstablished=false` and
`productionDomainModified=false`. Increasing the box reduces these measured
differences, but three points with moving rescaled boundaries, differing RK
partitions and no independent time/inner-grid refinement do not establish
an infinite-box error bound or validate the parent's earlier history.

Any future grid redistribution on these fresh branches must preserve their
physical anchor1/parentCx and transform the immutable root grid reference
to the same units. The existing platform helper with a hard-coded anchor1
cannot be applied directly to these physical-coordinate checkpoints.
