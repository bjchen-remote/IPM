# Arbitrary-anchor grid redistribution for a verified fresh finite box

The small-box natural-epsilon branch at local physical elapsed .01 has
901 x 386 nodes, anchor `0.21696981920359623 = 1/parentCx`, and its own
native history beginning at local time zero. Its absolute epoch remains
1.7798569844190184. Its actual computational endpoints are
x=+-2124.6819635072989 and ymax=2166.8115091110594. This is a new physical
initial-value branch from original step2031, not a changed-domain version
of the parent's native history.

## Construction and provenance

`mesh_fresh_grid_design(checkpointFile,resultFile,initialDataFile,outDir,controls)`
normally reads both current fresh and original parent checkpoints at ten
MATLAB threads. It checks the entire trusted histories, exact current
rho/axis/step/clock/case-ID pairing, and reconstructs the saved initial
physical axes and density bitwise from the parent crop. It does not restore
operators or build an LU. The immutable reference is those actual initial
fresh axes; current runtime references come from this fresh native case.
They are neither the old H1e6 root reference nor generated factory axes.

`mesh_anchor_grid_candidate(dataset,reference,anchor,candidateFile,controls)`
normalizes lengths by the positive anchor. `mesh_general_rounded_axis`
constructs a direct rounded-log candidate with an arbitrary odd node budget,
exact +-anchor nodes and the true box endpoints. Its bridge formulas and
balanced split criterion follow the existing rounded factory, but it has
no q256/power-of-two restriction or historical/nested-factory reproduction
claim. Only candidate construction is performed in normalized units.

The y equalizer uses the immutable fresh initial y axis divided by anchor.
If the candidate y count differs, a separately identified PCHIP resampling
in node-index coordinates supplies its reference budget; the original
reference remains stored unchanged. Every candidate is then multiplied
back and undergoes the complete real-unit paired-field score.

Hard gates are unchanged: adjacent ratio<=1.08, log-spacing curvature<=.01,
stencil rcond>=1e-9, global q/mean(q)>=1e-8, every q positive and local
q/control-width in [.35,1.65], exact +-anchor and box endpoints. Production
design targets cannot fall below21/21 and the left-front target is20.
Frozen transfer gates retain peak change<=.002, mass defect<=5e-12 and
range violation<=2e-4. Construction exceptions are saved as failure.mat/json;
mesh/transfer rejections retain their full scored candidate.

The frozen pair score uses its existing separable transfer. In particular
its vertical correction is not the production smooth-bubble correction.
Its predicted core counts are not native post-transfer measurements.

## Actual no-LU results

The completed source is
`fresh_box_tau428_ladder_v1/crop_natural`, native step41, local elapsed .01.
Initial data are from
`fresh_box_tau428_protocol_v1/crop_natural/initial_physical_data.mat`.
All paths below are relative to `result/longtime/20260908_campaign_v1`.
The MATLAB process exited0; log is `../20260908_fresh_anchor_design_v1.log`.

| Nodes | Target x/y | Frozen core x/y | Left-front cells | Peak relative jump | Current native count contract |
| --- | --- | --- | --- | --- | --- |
| 901 x 386 | 24/24 | 26.4 / 27.6 | 29.2945 | 3.80925e-5 | Supported |
| 901 x 386 | 32/32 | 35.2 / 36.8 | 39.0594 | 3.80984e-5 | Supported |
| 901 x 386 | 42/42 | 46.2 / 48.3 | 49.9030 | 3.23234e-5 | Supported |
| 1025 x 513 | 24/24 | 26.4 / 27.6 | 29.2945 | 3.80925e-5 | Rejected: changed counts |
| 1025 x 513 | 32/32 | 35.2 / 36.8 | 39.0594 | 3.80984e-5 | Rejected: changed counts |
| 1025 x 513 | 42/42 | 46.2 / 48.3 | 50.1431 | 3.23234e-5 | Rejected: changed counts |

All six pass every frozen gate, with mass defect<=8.74e-15 and range
violation<=7.41e-7. Files are
`fresh_anchor_design_{901,1025}_target{24,32,42}_v1/candidate.mat`.
The target24 same-count directory also stores the fully checked dataset and
immutable reference in design.mat. `fresh_anchor_design_summary_v1.json`
is a lightweight read of saved candidate scalars; it does not itself perform
native signature validation.

Registered controls are fineCells96, anchorFineCellFraction.85,
roundingCells40, horizontal padding1.10, vertical padding1.15 and
yEqualizationSigma.25. Neither padding lowers an acceptance target.
The exact engine call for the proposed first transaction candidate is:

```matlab
candidate=mesh_anchor_grid_candidate(design.dataset,design.reference, ...
    design.anchor,candidateFile,struct('nodeCountX',901,'nodeCountY',386, ...
    'targetXCoreCells',32,'targetYCoreCells',32));
```

## What the native contract currently supports

`ipm_gridlab_regrid_checkpoint` preserves both node counts and box endpoints.
Its maintained transfer has no anchor=1 assumption; origin/pin indices are
recomputed and the native reference values are preserved. The tiny test
actually passed its same-grid transaction at anchor .92188077881673891.
A changed count was explicitly rejected with `ipm:gridlab:RegridNodeCount`.
No production helper or core solver was modified for this design.

The existing transaction does not explicitly require anchor to be a node
or impose the local q/control-width bounds. A future wrapper must retain
the candidate's complete axis gates before invoking it. The old
`mesh_test_candidate` expects the original campaign design schema and a
fixed laboratory window; it is not a drop-in fresh-box smoke wrapper.
The changed-count candidates have `transactionReady=false`. Supporting
their native transaction would require a separately reviewed opt-in change,
tiny transfer/restore tests and a new production validation window.

## Tiny evidence and next dynamic experiment

`20260908_anchor_grid_tiny_v3` passed with exit0 and all four helper files
were Code Analyzer clean. Scaling an analytic paired fixture by1, .25 and4
produced exactly equal normalized candidate axes, dimensionless quality
metrics and relative full-pair metrics. The additional non-binary scale .37
had normalized-axis error<=1.43e-14, quality error<=2.50e-14 and relative
paired-metric error<=5.83e-15. The native
nonunit-anchor identity transaction preserved all rescaling references.
The native changed-count rejection was detected. This fixture uses a
declared core floor2 only for its 65 x 33 native transaction; production
mesh limits were retained. The synthetic covariance candidates were fully
mesh-scored but were not themselves promoted by the target-resolution gate.
The separate `20260908_fresh_anchor_lineage_tiny_v1` exercised the real
fresh/parent provenance wrapper and scored an expanded candidate, without
promoting its changed count to a native transaction. Failed tiny v1 is kept;
it used an incorrect score field name, corrected before v2.

The next proposed test is the same-count target32 candidate. Require the
actual native transaction to achieve cores>=31/31, plus the unchanged
peak/mass/range/axis gates, exact anchor, unchanged reference values and
unchanged clocks. Then evolve old/new grids to the same physical endpoint
(suggest first an additional .001), at ten threads, separately and serially.
Require complete trusted histories, native/result pairing, actual endpoint
arrival, both cores>=14 and safety<.70 throughout the recorded interval.
Compare rho, rho-x, rho-y, physical gradient peak and physical rates in the
same laboratory overlap window used by the finite-box protocol, with the
existing epsilon unchanged. This isolates grid transfer and density; it
does not repeat or replace the three-box boundary experiment.

The native transaction costs an old-grid restore then a new factorization;
its existing early old-LU release applies. Each of the two short solves has
one factorization, and returned states/observations must be reduced before
the next. No q512 LU, native grid transaction or PDE was run while preparing
these candidates. A successful dynamic test would authorize an independent
fresh-box continuation experiment, not changing the main H1e6 domain or
declaring an infinite-time PROFILE.


## Matched-time wrapper validation and allocated run

`mesh_test_fresh_grid` now implements that two-branch protocol. It revalidates
current/parent native lineage and the saved initial data, rescoring the
complete paired candidate before any LU. It explicitly preserves all
runtime reference values except recomputed origin/pin node indices, and
checks the actual transferred left front before writing the regrid CP.
The existing same-count native transaction still enforces its measured core,
peak, mass and range gates. Complete recorded histories, native/result
pairings, fresh time-zero lineage and actual physical endpoint arrival are
required for both short branches. The unchanged linear comparison gate and
additional bidirectional C1 Hermite observer are both saved.

The wrapper tiny v1 and v2 passed with exit0 in
`result/longtime/20260908_fresh_grid_wrapper_tiny_v{1,2}`. The latter supplies
x as a column and y as a row deliberately; canonical x-row/y-column
normalization preserves their values before native endpoint checks. Both
65 x 33 old/new branches took two accepted steps, matched physical time
exactly, preserved lineage, passed normally read native pairings, and agreed
in both observers to <=3.46e-14. Both wrapper files were Code Analyzer clean.

The coordinator then allocated the actual 901 x 386 target32 transaction
and .010-to-.011 physical elapsed comparison. The v1 wrapper attempt stopped
before any LU or PDE because vector-oriented endpoint `isequal` compared a
candidate y row against a native y column. That evidence is retained in
`fresh_grid_tau428_target32_smoke_v1/preflight_failure.json`; no state or
candidate values were changed. The fix was covered by tiny v2 and frozen
separately in `fresh_grid_tau428_target32_source_v2/source_manifest.json`.
The v2 actual run was then started in its own output directory; completion
must be established from its smoke report and process exit, not this plan.


The v2 actual run subsequently completed and MATLAB exited0. The native
transaction passed with actual core35.2387648/36.9891749, left front37.7827563,
peak jump3.80984e-5, mass defect2.072e-16 and range violation7.409e-7.
Old/new grids then took5/6 accepted steps (24.81/27.18 seconds) to the same
physical endpoint within3.39e-14. Both entire recorded histories, normal
native CP/result pairings, fresh lineage and core14/safety.70 passed. The
new terminal core was35.0252255/36.8919815; the old was17.0396878/17.8636737.

**The original combined smoke gate rejected this candidate.** The only
exceeded registered comparison was the legacy physical y-core width:
relative difference.00859685545 exceeds.005. Linear rho/rho-x/rho-y Inf
were.000298528/.000616301/.000588755, physical gradient peak difference
4.26323e-5 and quadratic peak difference7.51033e-6. Parent-canonical
c-l/c-omega relative differences were7.67629e-7/.000492887, within their
registered limits. Additional Hermite observer maxima were rho Inf2.137e-6,
omega Inf1.83548e-4 and rho-y Inf2.70886e-5. These smaller observer errors
were retained separately and did not replace the width rejection. The
saved `smoke_report.json` remains `passed=false`; process exit0 means the
experiment/report finished, not that its numerical acceptance gate passed.


The bounded follow-up `mesh_diagnose_fresh_grid_width` completed with exit0,
no LU or PDE/RHS calls, using the coordinator's validated frozen geometry
in `continuous_forecast_tau509_709_v1/source`. Its common physical peak
window is [0,.43393963840719246], twice the initial physical anchor, mapped
separately to each state's computational coordinates by its actual Cx.
All four native signatures and complete trusted histories were checked at
10 threads. The legacy width was reconstructed independently from its
actual nodal peak column and matched the recorded value. The original smoke
file and false decision were not modified.

| Pair | Legacy y-width relative difference | C1 crossing at original nodal columns | C1 width at each continuous peak | C1 width at one common physical x |
| --- | --- | --- | --- | --- |
| Zero-time source / regrid | .00853927 | .00853518 | 2.33401e-5 | 1.45361e-6 |
| Matched old / new endpoint | .00859686 | .00859201 | 2.40113e-5 | 2.26709e-6 |

At the matched endpoint, the old nodal peak column539 is +4.94370e-5 in
physical x from the continuous peak (+.13337 local cells), whereas new
column602 is +7.62566e-5 (+.42502 local cells). Thus a finer grid need not
place its maximum node closer to the continuous maximum. Merely improving
the y crossing interpolation leaves almost the entire .86% discrepancy.
The legacy width signed difference is +1.64420605e-5; changing the vertical
crossing to C1 contributes -9.12839e-9 to that difference, and moving to the
continuous peak column contributes -1.63877225e-5. The remaining continuous
width difference is +4.52096e-8. These corrections sum algebraically and
were asserted by the audit. Almost the entire observed legacy mismatch is
therefore attributable to peak-column selection in this frozen comparison.

The actual continuous physical y widths are .00188284904746 and
.00188289425702. At the shared physical x=.206971391221 they are
.00188286951890 and.00188287378753. These are independent observation results,
not revised native diagnostics or a retrospective pass. Zero-time mismatch
is already nearly as large as the matched endpoint mismatch, so this audit
also distinguishes the grid observation effect from its brief evolution.

The audit records scales and rate units explicitly. At the endpoint the
old/new fresh Cx values are1.01964732450/1.01964732559 and fresh Comega values
.970614580443/.970613259223. Dividing fresh canonical rates by the original
lambda0=8.92162799960 yields parent-canonical c-l
.190299532207/.190299386127 and c-omega
-.298420809739/-.298567897609. Multiplying fresh canonical c-l by the current
fresh Cx/Comega yields the actual physical-time log-length rates
1.78354883230/1.78354989290. Reconstructed parent Cx/Comega and physical-time
amplitude rates are also saved, so canonical-rate normalization is not
confused with a derivative in physical time.

Full evidence: `fresh_grid_tau428_target32_smoke_v2/continuous_width_audit_v1.json`.
The observer supports designing a separate, prospectively registered
continuous-geometry validation experiment. It does not lower the existing
.5% gate, upgrade the failed smoke, or authorize the main domain change.
