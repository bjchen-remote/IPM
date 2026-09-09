# Long-time schema-4 profile work

The production entry point remains `ipm.solve`. These research scripts design
and audit grids, continue accepted checkpoints, and inspect paired fields.
Checkpoint validation for this campaign uses the original 10-thread MATLAB
environment. The archived native checkpoint signature is checked normally;
changing the thread count must not be used to bypass it.

## First feasible grid at tau 4.073201938625015

Source: `result/longtime/20260908_campaign_v1/baseline_segment_002/result.mat`
and its native step-1936 checkpoint. All references come from the immutable
unregridded schema-4 H=1e6 root, not the older campaign's H=281 dataset.

The bounded Gaussian equalizer search stopped at its first fully admitted,
exact-anchor candidate:

| x/y target | x sigma | Outcome |
| --- | --- | --- |
| 24/24 | 1.5 | Original stage design: infeasible x |
| 21/21 | 1.5 | No x candidate met both the core and front targets |
| 21/21 | 2 | No x candidate met both the core and front targets |
| 21/21 | 1 | Feasible after independent pair scoring and exact x=+/-1 anchoring |

All attempts retain the x-front target of 20 and the same hard mesh gates:
adjacent ratio <=1.08, log-spacing curvature <=0.01, seven-point stencil
rcond >=1e-9, positive quadrature with global minimum/mean >=1e-8, and local
quadrature/control-width ratios inside [0.35,1.65]. The failed full designs
retain each candidate, its attained counts, mesh diagnostics and rejection
status under `result/longtime/20260908_campaign_v1/mesh_tau407_search_v1`.

The accepted axes are in
`mesh_tau407_search_v1/trial_03/mesh_anchor1_exact_candidate.mat` under that
campaign directory, variable `candidate`. Strict cores are
`21.0000000327 / 21.0000000116`, the frozen rho-x maximum jump is
`1.92799036e-4`, field round-trip relative L2 is `5.65324490e-6`, and the y
global quadrature ratio is `7.67141265e-8`. The coordinator subsequently
passed its native transaction and matched-physical-time smoke test in
`mesh_tau407_smoke_v1/mesh_smoke_report.json`: transferred cores
21.0294/21.8734, physical endpoint mismatch 5.33e-15, both histories trusted,
physical rho-x maximum difference 2.29345e-4, c-l difference 1.21486e-6,
and c-omega difference 4.54839e-4 (all relative magnitudes). Laboratory-frame
local rho/omega/rho-y L2 differences were 3.69314e-4/1.27881e-3/1.35956e-3.
This short test supports that transaction; later stages retain their gates.

Reproduce the selected family on a later state with:

```matlab
controls = struct('targetXCoreCells',21,'targetYCoreCells',21, ...
    'xEqualizationSigma',1);
mesh_design_stage(newResult,rootResult,newDirectory,controls);
mesh_balance_y_candidate(designFile,referenceFile,newCandidateFile,1);
```

`mesh_design_stage` exposes only targets and candidate-family controls; hard
gates remain fixed. `mesh_balance_y_candidate` inherits the selected x
equalizer's sigma, power and optional corridor geometry. It no longer resets
a successful family to sigma=1.5. Its exact anchor correction is recorded,
and the complete corrected axes are scored again. Controls selected for this
source are not proof that the same family remains feasible at later times.

## Adaptive continuation gates

Target 24/24 cells when feasible; the accepted 21/21 fallback uses a
transaction core floor of 20/20 (`target-1`, never below 16). A stage of at
most Delta tau=0.2 with native checkpoints every 0.1 is reasonable for the
measured source. Trigger a design when either current core is below 20.
Before launch, also estimate both endpoint cores and safety from recent
records with the same `remeshCount`. If a 0.2 segment predicts a core below
16 or safety above 0.65, first improve the grid or shorten the segment.
Do not fit a core-loss rate across an instantaneous positive regrid jump.
The actual terminal gates remain both cores >=14 and safety <0.70.

Sum the absolute per-regrid rho-x maximum jumps from all inherited
`gridLabRegrids` provenance, not only the current invocation. Retain a
per-regrid limit of 0.002 and a cumulative budget of 0.02. Check the predicted
candidate jump before transfer and the actual audited jump before continuing;
to require budget approval before installing a checkpoint, create the
transaction in memory and write it only after the cumulative check. An
exhausted budget requires a new validated branch/resolution strategy.

## Moving fine platform after the original x family is exhausted

At the stage-3 endpoint tau=4.6832, target21/24 with Gaussian sigma
.35/.25/.18 all failed. Their complete designs are in
`20260908_campaign_v1/mesh_narrow_stage003_v1`; the same stage-2 scan had
passed target21 with sigma .35 and .25. At stage 3, the maximum admitted x
core counts were only 19.4328/18.4852/17.8567, limited by ratio=1.08.
Original-root corridor shoulders16/24/32 also failed. A separate anchor-first
scan over Gaussian .18–.6 and corridor shoulders16–96 likewise failed. No
target or hard gate was reduced, and all these failures remain saved.

The root-anchor equalizer keeps x=1 at positive-half node104. The new
`mesh_platform_candidate` keeps x=1 exact while allowing its index to change.
It calls the existing rounded-log grid constructor with a smaller fine
platform centered to include the evolving core and front. The unchanged
factory is first required to reproduce the immutable original root x axis
exactly. The root dataset and original-root y equalizer remain unchanged;
no candidate/current axis becomes a recursive reference.

The smallest currently validated parameter registration is:

```matlab
% Source design retains targetXCoreCells=targetYCoreCells=21 and xFront=20.
controls=struct('positiveFineCells',96,'anchorFineCellFraction',.85, ...
    'resolutionPadding',1.10);
candidate=mesh_platform_candidate(designFile,referenceFile,outputFile,controls);
```

The requested fine spacing is
`min(coreWidth/targetX,frontWidth/frontTarget)/1.10`. Of its 96 positive-half
fine cells, 82 lie below x=1 and 14 above. Every corrected axis is scored
again after exact +/-1 snapping, with all prior mesh gates plus explicit
target21/front20, peak-jump2e-3, conservation5e-12 and range2e-4 checks.
The source Gaussian design may have `status=infeasible_x`: it supplies the
paired frozen data and immutable references, not an accepted candidate. A
driver must therefore attempt the platform independently of that old status.

The first platform candidate passed frozen tests at both tau4.7832 and
tau4.8832. The latter is
`20260908_campaign_v1/mesh_platform_bridge2_v1/candidate_01.mat` with
`source_design/mesh_design.mat` and `source_design/mesh_root_reference.mat`.
It moves the positive anchor index from104 to180 and gives core23.1/21,
front26.6656, x ratio1.062876914, x curvature0.001908415, y global quadrature
ratio4.68752e-8, peak jump1.67958612e-4, field L2 3.47272e-6, conservation
defect1.72326e-14 and range violation1.64388e-7. The native transaction later
rejected this initial candidate because the transferred y core was 19.453,
below the unchanged 20-cell gate. The coordinator added
`verticalResolutionPadding=1.15`, increasing only the candidate's y
equalization request, and tested a separate frozen version before running
adaptive v3. The original candidate remains evidence of a rejected
transaction; it must not be reused as a production acceptance. The
coordinator owns the current platform helper and frozen v3 manifests.

For the next driver, try this one platform registration first and retain old
families as separately logged fallbacks. Do not register untested platform
controls as accepted choices. Frozen target/quality success at two times is
not evidence that the same platform will remain feasible indefinitely.

## Diagnostic interpretation and verification status

`analyze_profile_series` differentiates each rho on its own saved x axis and
retains paired x/y trace coordinates. Inner-profile comparison removes each
sample's displacement, core widths and amplitude; it measures normalized
shape, not convergence in the original rescaled coordinates. The selected
vertical column and its displacement from the subgrid peak are reported.
Physical core widths and sample case IDs are also explicit.

`screen_frozen_box` keeps common rho/grid samples and restores the maintained
runtime-reference whitelist without feedback. It checks that every stored
normalization/resolution reference survives and that the crop retains the
entire comparison core and a three-node outer stencil margin. The original
LU is released before building a cropped operator. Recomputed outer
derivative closures and artificial boundary data still change the numerical
problem. The current Green boundary uses the cropped source to recompute its
traces; it does not retain the omitted source contribution. A stable frozen crop does not validate prior
dynamics or prove that the original large box is free of truncation error.

The two diagnostics passed the tiny 65x33/67x35 paired-grid fixture with
Green boundaries in `result/longtime/20260908_diagnostic_tiny_v3`: original
10-thread native checkpoint validation, caller-thread restoration to 2,
paired x/y arrays, full-box core velocity overlap exactly zero, preserved
runtime references, and rejection of an undersized crop. A 0.75 crop gave
core velocity L2 difference 0.136569, demonstrating sensitivity rather than
box convergence. Both scripts passed Code Analyzer. Earlier v1/v2 failures
and logs remain; the passing fixture and report are saved in v3. The
tau-4.0732 design search used only frozen-field and one-dimensional grid work;
its two mesh helper scripts also had zero Code Analyzer findings.

## Independent finite-box work and monitoring

`mesh_compare_initial_box_summaries` reads the old initial-RHS screen only.
The resulting `20260908_initial_box_axis_audit_v1` confirms identical peak
platforms within 1e-12 but different complete inner x/y grids across H.
Consequently its very small initial rate differences are candidate evidence,
not fixed-inner-grid or late-time box convergence. The experiment ladder,
exact proposed source checkpoint, and separately preregistered static and
dynamic criteria are in [mesh_finite_box_plan.md](mesh_finite_box_plan.md).

`mesh_screen_late_box` prepares a static report and never installs a grid or
advances the PDE. Its additional relative-infinity/peak telemetry and static
classification wrapper passed the coordinator's tiny run in
`20260908_box_wrapper_tiny_v1`: full control passed, the sensitive .75 crop
was rejected, no dynamics was promoted, and caller threads were restored.
Code Analyzer passed. Earlier Qt/neon failures occurred before MATLAB
executed and remain recorded separately. The tiny result does not establish
a q512 late-box pass; that independent run retains its own report.

`read_checkpoint_scalars.py CHECKPOINT.mat` uses the local MATLAB HDF5
library through Python ctypes and reads only clocks and final numeric
common/mesh/gauge records. It loads no rho, assembles no operator, writes no
checkpoint, and performs **no signature or trusted-prefix validation**. This
is progress monitoring only; native restoration still requires normal
10-thread MATLAB validation. System and bundled Python produced identical
step-1936 telemetry in about 0.03 seconds. Seven tiny-checkpoint clocks,
rates, core counts and peak values also matched the successful native MATLAB
diagnostic output exactly. Unsupported fields are listed explicitly, and
nonfinite numeric values are labeled.

## Fresh late-state finite-box dynamics

The actual stage-1 step-2031 static screen is saved in
`20260908_campaign_v1/late_box_tau428_v1`. Fractions through .003 passed;
.001 failed both core velocity tests. The conservative .01 candidate retains
901x386 of the original 1025x513 nodes, with core velocity relative L2/Inf
1.37981e-4/1.20419e-4. It is eligible for an independent dynamic experiment;
the static result does not establish box convergence.

`mesh_fresh_box_branch` samples the saved late state exactly on its own
physical axes, divides the density by the original amplitude scale, and
starts a new `ipm.solve` IVP with local physical time zero. Its physical
transport anchor is the parent anchor divided by the original Cx. Each
branch has a new case ID, new native history and a separate absolute epoch;
no changed-domain checkpoint is fabricated. Runtime references are freshly
initialized and explicitly labeled. `mesh_audit_fresh_references` can compare
the transformed parent references against them at a frozen child endpoint,
using one additional restore and no PDE.

The physical covariance/reference/sensitive-crop tiny test passed in
`20260908_fresh_box_tiny_v1`: full-box rho/rho-x/rho-y relative infinity
differences were 6.75e-16/1.34e-13/1.10e-13, transformed-reference RHS and
rates agreed exactly in both full/crop cases, and the intentionally small
crop produced rho L2 difference 1.37e-3. Those reference comparisons are
instantaneous tests, not separately evolved reference trajectories.

`mesh_run_fresh_box_protocol` implements five serial solves: native parent
continuation, full/crop fresh with natural epsilon, and full/crop fresh with
one common epsilon. WENO's effective regularizer is
`max(wenoEpsilon,1/(N-1)^2)`, so changing the total node count changes inner
weights even when inner nodes coincide. The shared-epsilon pair controls
this difference; it does not hold the global Lax-Friedrichs line maxima or
outer closures fixed. Full-box native/fresh covariance gates run before any
cropped solve. All five production endpoints require fully trusted histories,
both actual cores >=14 and safety <.70. The 20-cell transaction gate does not
apply because these branches perform no regrid or transfer.

The five-case wrapper itself passed `20260908_fresh_box_protocol_tiny_v1`:
five distinct IDs (one inherited native ID plus four new IVPs), covariance
rho-x difference 1.34427e-13, both sensitive-crop comparisons rejected,
Code Analyzer clean, and no dynamic box-convergence promotion. Tiny fixture
resolution exemptions are reported and cannot apply to a >20000-node case.
The registered q512 source, thresholds, units and exact call are in
[mesh_fresh_box_protocol_plan.md](mesh_fresh_box_protocol_plan.md). Preparing
that registration did not launch a q512 PDE. The coordinator subsequently
allocated one q512 window, and the frozen protocol completed all five cases
with MATLAB exit0 in `20260908_campaign_v1/fresh_box_tau428_protocol_v1`.
Every branch reached the physical target within 2.45e-14, retained complete
trusted result/checkpoint histories, paired its terminal native rho/axes/
clocks/case ID exactly with its result, and passed core>=14/safety<.70.
Full fresh/native covariance rho-x/rho-y differences were 2.64e-13/1.14e-10.

For elapsed physical time .001, natural-epsilon full/crop differences were
rho L2=4.64e-7, rho-x L2/Inf=7.62e-7/3.68e-6 and rho-y
L2/Inf=7.11e-7/4.23e-6. The physical gradient peak changed by 1.55e-7;
c-l/c-omega relative changes in parent units were 8.82e-5/1.74e-8.
The common-epsilon pair was almost identical, while changing epsilon alone
in the full box gave rho-y relative Inf difference 4.60e-12. Terminal cores
were at least 18.12965/18.54205 and safety at most .441267. These are short
interval results; they nominate a longer independent box experiment and do
not establish dynamic box convergence or change the production domain.

The proposed next interval is .01 total physical elapsed, with three boxes
full/.1/.01 and both natural/common epsilon. `mesh_predict_box_window.py`
uses only same-remesh scalar tails with a declared 1.25 planning factor;
it predicts cores16.873/17.316 and safety.47135 at .01. At .02 its predicted
x core15.534 falls below the preparation buffer16, so .01 is selected first.
This ctypes forecast does not validate signatures or guarantee future
resolution. The actual native/trusted/core14/safety.70 gates still govern.

`mesh_extend_fresh_box_ladder` reuses four existing fresh native checkpoints
and their real time-zero histories, adding only the two middle-box IVPs.
It preserves natural/common epsilon separation over the longer interval.
The minimum protocol has six serial solves; the earlier native/fresh
calibration is reused and explicitly not claimed as remeasured at the new
target. Its tiny fixture passed in `20260908_fresh_box_ladder_tiny_v1`:
four histories/IDs were continued, two new middle-box histories were added,
all six terminal native/result pairings and physical endpoints passed, and
the sensitive small box remained rejected under both epsilon families.
The forecast, resource estimate and prepared call are in
[mesh_fresh_box_ladder_plan.md](mesh_fresh_box_ladder_plan.md).

The coordinator then allocated a window and the .01 six-case ladder actually
completed with MATLAB exit0 in
`20260908_campaign_v1/fresh_box_tau428_ladder_v1`. Both natural/common
three-box comparisons passed. Natural full/small rho L2 was4.53e-6,
rho-x/rho-y Inf3.95e-5/4.58e-5 and physical peak difference1.61e-6;
the middle box's errors were about one order smaller. Within-box epsilon-only
gradient Inf differences were at most9.09e-9. All six native/result pairings,
complete trusted histories and physical endpoints passed; minimum terminal
cores17.146948/17.933558 and maximum safety.46655534 retained ample margin.
The measured solve total was15.36 minutes. This result is a fixed-interval
box test, not an infinite-box/asymptotic-profile conclusion; .01 native
covariance was not remeasured. Supplemental scalar-record audits also passed
the post-run residual check and clearly distinguish it from trustedMask.

The independent finite-box redistribution tools `mesh_fresh_grid_design`,
`mesh_anchor_grid_candidate` and `mesh_general_rounded_axis` now support any
positive anchor and arbitrary odd horizontal node budgets. They keep the
actual fresh initial axes as immutable reference, validate those samples
against the original native parent, and score every complete candidate in
the real fresh units. Tiny v3 actually passed unit covariance (including non-binary scale .37)
and the existing native transaction at nonunit anchor; changed native node counts
were correctly rejected. The six no-LU real designs at901x386/1025x513 and
targets24/32/42 all passed the unchanged frozen hard gates. At the frozen-design stage, only same-count designs were eligible for
a future native transaction; expanded-node designs remain unsupported. Exact files, evidence and the proposed
target32 matched-time smoke are in
[mesh_fresh_grid_plan.md](mesh_fresh_grid_plan.md).


The coordinator subsequently allocated the 901x386 target32 native smoke.
Its v2 run completed/exit0 with actual transferred cores35.24/36.99 and
both matched-time branches trusted, paired and resolved. The combined
numerical gate remains **false**: legacy vertical physical core width
changed.8597%, above its.5% limit. All other registered field/peak/rate
comparisons passed; extra Hermite observations remain separate and do not
replace that rejection. A pre-LU v1 endpoint check failed only on vector
orientation and is retained, with the fix covered by tiny v2. The exact
outcomes and additional no-LU width diagnosis are tracked in the dedicated
fresh-grid plan above.


A subsequent no-LU audit of source/regrid and old/new terminal native CPs
identified the width discrepancy as a nodal peak-column observation effect:
C1 crossing alone leaves the.859% difference, while measuring at continuous
wall peaks gives y-width difference2.40e-5 and using one common physical
column gives2.27e-6. Zero-time legacy mismatch was already.854%. All original
rejection evidence remains unchanged. The audit uses the root's frozen,
validated C1 geometry and a common physical window based on the actual
nonunit anchor, and reports rate/scale transformations explicitly.

For the newer step3326 fresh-box epoch, all24 requested anchor-centered
constant-patch designs failed:96 fine cells missed the core/front, while
larger patches exhausted the adjacent-ratio budget. The independent
`mesh_core_patch_axis` / `mesh_fresh_core_patch_candidate` family instead
centers the fine patch on the measured core and places the exact anchor on
an existing transition node with a compact C2 coordinate warp. Eight real
895x386 frozen candidates passed every unchanged hard gate, with scale
covariance checked and no LU/PDE. The proposed target32/fine64/.65 candidate
predicts cores35.23/36.8 and front43.62; it still requires an actual native
transaction and matched-time validation. All failures, source lineage,
candidate paths and limits are in
[mesh_fresh_core_patch_plan.md](mesh_fresh_core_patch_plan.md).


The latest core-patch candidate then completed an allocated actual native
transaction and .0003 matched-physical-time test in
`fresh_epoch3326_core_patch_smoke_v1` (MATLAB exit0). Actual transaction
core35.21/37.10 and final new-grid core35.01/36.85 passed; all native/prefix/
lineage/linear-field/rate/peak gates passed. Two decisions remain separate:
original legacy smoke=false because nodal y-width differed1.481%, while
the pre-registered continuous same-peak width auxiliary=true (matched
wall/y differences1.062e-4/1.190e-5 at the unchanged.005 threshold). The
original false reports are preserved. This is still a short independent
finite-box experiment, with no main-domain or asymptotic-profile promotion.

The independent `run_fresh_profile_campaign` driver is implemented and its
public native/qualification preflight has passed actual no-LU checks on
fresh steps19/30/139. The proposed bootstrap is the completed fine step139
at equivalent tau7.243200751859, core32.76/36.32; the next .2 stage predicts
28.96/33.84, so no immediate regrid is requested. The unchanged legacy
false and prospective auxiliary true require a caller ledger. All five
auxiliary subgates are independently enforced and negative-tested. Final
Code Analyzer found zero issues; the campaign itself has not been run.
Its immutable reference, fixed895x386 box, target32/native31 transfer gates,
clock definitions, actual reports and invocation are documented in
[mesh_fresh_profile_campaign_plan.md](mesh_fresh_profile_campaign_plan.md).

The independent step139 finer-grid screen subsequently completed12
target42/56 × fine64/96/128 × fraction.65/.5 trials on the same895x386 box.
Seven candidates passed all original frozen gates; both selected fine64/.65
families passed axis and complete paired coordinate-unit tests (maximum
paired invariant difference2.03e-12). Five candidates were rejected by the
unchanged1.08 x cell-ratio limit. No native transaction or PDE was run.
Exact candidates, margins, failure evidence and the complete frozen source
closure are recorded in
[mesh_fresh_resolution_family_plan.md](mesh_fresh_resolution_family_plan.md).

The selected target42 candidate then completed one actual267-step fine
branch to the existing coarse step328 absolute time1.964631495951036.
Native/field/rate/peak gates and the registered continuous-width plus
Hermite auxiliary gates passed; the original legacy y-width decision
remains false. Paired true RHS caches for fine step406 and coarse step328
were saved serially, and both extra MATLAB/LU processes exited0. The
dedicated finer-family plan above records the actual limits and results.

Matched-epoch true RHS analysis then exposed32.53% core derivative sensitivity
between the coarse and target42 branches. The registered target56 third
level completed355 steps to the same actual endpoint, with core53.2867 /
62.1355, safety.150131 and all native gates passed. Its original comparison
and combined auxiliary decision remain false: linear rhoX L2=.00204027
exceeds the original.002 gate. Continuous widths, Hermite fields, peaks and
rates pass separately. Its complete true RHS cache is saved, and MATLAB
exited0 with the extra LU released. This prepares the separate three-level
residual trend analysis; it does not establish spatial convergence.

Future batch restart across remeshed grids now has an independent no-LU
native qualification-chain checker. Actual v1stage2 and v2stage1 chains,
including the step504 transaction, passed; missing-stage and wrong-source
negative fixtures were rejected. No running driver was modified. See
[mesh_fresh_qualification_chain_plan.md](mesh_fresh_qualification_chain_plan.md).
