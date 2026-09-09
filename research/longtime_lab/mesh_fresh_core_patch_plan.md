# Late fresh-box core patch with an exact transition anchor

The latest audited small-box natural branch starts at absolute physical
epoch1.9545271366857122 from parent native step3326. The completed .001
fresh branch has895 x 386 nodes and anchor.15214223499880195. Its strict
frozen core center is.14820156372828294, so a fine interval forced to reach
the anchor has to span much more than the .0007640731 core width.
The actual fresh initial physical axes and density are kept as immutable
references and reconstructed bitwise from the normally read parent CP.
This lineage is unrelated to the earlier fresh case at epoch1.77985698.

## Completed constant-patch failure screen

`mesh_scan_fresh_anchor_family` tested the requested24 combinations:
target24/32, fine cells96/128/160/192/224/256, anchor fractions.90/.95.
All failed. With96 fine cells, the complete mesh was admissible but the
actual feature lay outside its useful fine interval: predicted x core was
only6.37–11.77 and left front5.83–9.45. At128 cells, the maximum adjacent
ratio already reached1.0827–1.0853, exceeding1.08. Larger fine counts reduced
the outer-cell budget further;224 reached1.1406–1.1449, and256 additionally
failed the quadrature margin. No hard gate was weakened. The y axis could
reach27.6/36.8 predicted core throughout.

Every trial, full rejected score or construction failure, is retained in
`result/longtime/20260908_campaign_v1/fresh_epoch3326_anchor_family_v1`.
MATLAB session26371 completed with exit0. No LU or PDE was run.

## Independent candidate geometry

`mesh_core_patch_axis` constructs a rounded-log grid about the measured
core center rather than the transport anchor. Its short uniform patch
uses64 or96 cells. It retains the unchanged factory's rounded transitions
and balanced split criterion, applied in units divided by the actual
transport anchor. It is a candidate factory, never a reconstruction of the
historical crop axes.

The nearest existing transition node to the transport anchor receives a
compact displacement. In normalized coordinates, the bump is
`delta*(1-u^2)^3` for `abs(u)<1` and zero outside, centered at that original
node, with radius.25 anchor units. Its value and first two derivatives
vanish at the support endpoints. The node is then set to the exact anchor;
negative x is constructed by odd reflection. Domain endpoints, zero, node
counts and exact +-anchor are retained. Monotonicity is required. This
smooth warp is not assumed to preserve quality: all axes and complete paired
fields are rescored in their original fresh units afterward.

`mesh_fresh_core_patch_candidate` reuses only the separately admitted y
candidate and identical immutable fresh reference from the constant-patch
screen. It retains all original ratio1.08, curvature.01, rcond1e-9,
q-floor1e-8, positive/local q-control [.35,1.65], target/front and transfer
gates. No core solver, existing mesh helper, original candidate or failed
smoke result was edited.

## Actual frozen results

Eight candidates (target24/32, fine64/96, core fractions.50/.65) all passed.
All place the exact anchor outside the fine interval in a transition.
Output root is
`result/longtime/20260908_campaign_v1/fresh_epoch3326_core_patch_family_v1`;
the log is `result/longtime/20260908_fresh_epoch3326_core_patch_v1.log`.
MATLAB session80395 exited0; this was entirely no-LU/no-PDE work.

| Target | Fine cells / core fraction | Predicted core x/y | Left-front cells | Peak relative jump | Max adjacent x ratio |
| --- | --- | --- | --- | --- | --- |
| 24 | 64 / .50 | 26.4049 / 27.6 | 32.6851 | 1.10977e-4 | 1.064237 |
| 24 | 64 / .65 | 26.4355 / 27.6 | 32.7293 | 7.30456e-5 | 1.064404 |
| 24 | 96 / .50 | 26.3772 / 27.6 | 32.6466 | 1.27568e-5 | 1.072571 |
| 24 | 96 / .65 | 26.4148 / 27.6 | 32.7000 | 3.65170e-5 | 1.072224 |
| 32 | 64 / .50 | 35.2406 / 36.8 | 43.3741 | 1.24056e-4 | 1.066632 |
| 32 | 64 / .65 | 35.2329 / 36.8 | 43.6196 | 9.94476e-5 | 1.066575 |
| 32 | 96 / .50 | 35.2482 / 36.8 | 43.6403 | 4.92107e-5 | 1.074732 |
| 32 | 96 / .65 | 35.1912 / 36.8 | 43.5594 | 1.25592e-4 | 1.074574 |

The proposed first actual test is
`target32_fine64_fraction65/candidate.mat`, with x curvature.0017320,
x/y global quadrature ratios6.59489e-6/2.97561e-7, conservation defect8.34e-15
and range violation1.615e-8. Its maximum coordinate displacement is
6.08635e-5 in fresh units. Choosing64 fine cells leaves more outer-grid
margin than96 while all frozen transfer errors remain below their limits.
This is a proposed candidate, not an accepted native regrid.

The new factory's one-dimensional scale test used factors1, .25, .37 and4.
The largest normalized-axis discrepancy was2.274e-13 and the largest
dimensionless quality discrepancy4.330e-14. Exact anchors and outside-patch
placement were asserted in every case. Both new files were Code Analyzer
clean. These geometry tests do not replace native transfer or PDE tests.

Exact candidate generation, using the same saved immutable source:

```matlab
candidate=mesh_fresh_core_patch_candidate( ...
    fullfile(oldFamily,'target24_fine96_fraction90/design.mat'), ...
    fullfile(oldFamily,'target32_fine96_fraction90/candidate.mat'), ...
    newCandidateFile,struct('positiveFineCells',64, ...
    'coreFineCellFraction',.65,'anchorWarpRadius',.25));
```

## Required next evidence

The new candidate is compatible with the tested same-node
`mesh_test_fresh_grid` interface. A future allocated run must still require
the actual native transaction cores>=31/31, actual left front>=20, unchanged
peak/mass/range and full axis gates, exact anchor/reference/clock/lineage
preservation, and old/new evolution to one matched physical endpoint.
Both complete native/result histories and core14/safety.70 remain required.

The earlier epoch target32 smoke remains false because its legacy nodal
y-width differed.8597%, despite continuous same-peak width differing only
2.40e-5. A future independent protocol may preregister same-continuous-peak
widths at the same.005 threshold as an additional observation while retaining
the legacy difference and all original linear field/rate/peak gates. This
must not retroactively change the earlier rejection. The coordinator's
subsequent geometry MMS found small absolute errors/covariance errors but
failed its universal per-refinement ratio target; no generic fourth-order
geometry claim is supported.


## Completed native transaction and prospectively registered observation

The coordinator explicitly allocated the target32/fine64/coreFraction.65
native window. Before any PDE, an independent auxiliary protocol was saved
as `fresh_epoch3326_core_patch_smoke_v1_auxiliary_registration.json`.
It requires continuous same-peak wall/y widths to agree within.005 at both
zero time and matched endpoint, while preserving every original native,
linear-field, rate and peak gate. It retains the original legacy verdict.
The frozen geometry is exactly from
`continuous_width_observer_validation_v2/source`; its small absolute and
covariance errors are recorded along with its failed universal per-phase
refinement-ratio criterion. No fourth-order claim is made.

The actual run completed with MATLAB exit0 (session20947) and released its
LU. The original helper and all numerical source files were frozen in
`fresh_epoch3326_core_patch_smoke_source_v1` with131 SHA256 entries.
The immutable source is the latest fresh crop native step19. The target was
an additional physical elapsed.0003, absolute time1.9558271366857122.

The measured native transaction passed cores35.2110089/37.0953495 and
left front42.6694856, exceeding the required31/31 and20. Peak jump was
9.9447634e-5, mass defect2.127e-16 and range violation1.615e-8. Its wall
cost was22.38 seconds. Old/new branches took7/11 accepted steps and
29.18/38.93 seconds. Their physical endpoints differed by2.50e-16; all
terminal native/result pairings, full trusted histories, genuine fresh
lineage and recorded core14/safety.70 gates passed. The new terminal core
was35.0090961/36.8506381 (old17.5369921/18.3595101).

The two separately retained decisions are:

| Decision | Result | Evidence |
| --- | --- | --- |
| Original legacy smoke | False | Legacy y-width difference.0148131 exceeds.005; all other original gates pass |
| Prospectively registered C1 auxiliary | True | Matched continuous wall/y differences1.06216e-4/1.19011e-5; zero-time differences6.20798e-7/4.80543e-6; all original non-width gates pass |

The original linear rho L2/Inf differences are3.62287e-4/4.00816e-4,
rho-x1.19961e-3/9.47056e-4 and rho-y9.95997e-4/6.84688e-4.
The rho-x L2 comparison has less margin than the other field comparisons,
so a longer interval or spatial ladder remains necessary before a stronger
accuracy claim. Physical gradient peak difference is1.25949e-4 and
quadratic peak difference1.63796e-5; parent-canonical c-l/c-omega relative
differences are2.57297e-6/5.98681e-4. All are below the original limits.

The independent four-CP audit again identifies the legacy y-width shift
with the nodal peak-column choice: zero-time legacy difference was already
.0147768, versus C1 same-peak4.80543e-6. At the matched endpoint a common
physical column reduces the C1 y-width difference to1.92371e-6. The signed
legacy width shift2.97672915e-6 decomposes into a y-crossing correction
1.02050e-9 and peak-column correction-2.97527194e-6, leaving a continuous
width difference2.47771e-9. The old false results and this run's original
false result remain unchanged.

All actual reports are in `fresh_epoch3326_core_patch_smoke_v1`:
`smoke_report.json`, `transaction_audit.json`,
`linear_hermite_observers.json`, `continuous_width_audit_v1.json`, and
`auxiliary_verdict_v1.json`. A normally validated new-grid checkpoint and
paired result reside in `new_grid` at native fresh step30. This is an
independent short finite-box grid experiment; it does not itself modify
the main H1e6 domain, prove box convergence at later times, or establish an
infinite-time profile.
