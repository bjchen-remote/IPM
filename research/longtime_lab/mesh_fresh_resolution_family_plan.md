# Finer frozen candidates at fresh step139

The actual source is `fresh_core_patch_long_tau724_v1/result.mat` and its
paired native step139 under `result/longtime/20260908_campaign_v1`.
Absolute physical time is1.9590400984028897; parent-equivalent tau is
7.243200751859. The source is the separately qualified fresh895x386 case.
Its immutable reference remains the actual initial axes and samples in
`performance_lab_20260908/late_box_v4_stage002_epoch_v1/dynamic/crop_natural`.
No original or running campaign source was modified.

`mesh_scan_fresh_resolution_family` tested target42/56, fine64/96/128 and
fraction.65/.5, retaining fixed box, counts, exact nonunit anchors, front20,
all original axis/local-q gates and all frozen transfer gates. Successful
axes also underwent independent construction checks with coordinate-unit
factors.25,.37,4. The first admitted fine64/.65 candidate per target received
complete paired frozen transfer tests on separately labelled synthetic
coordinate transforms. These tests do not assert PDE or gauge covariance,
and synthetic products are not native transaction inputs.

## Actual results

The complete output is
`result/longtime/20260908_campaign_v1/fresh_tau724_finer_family_v3/screen/report.json`.
The scan ran from a129-file immutable source containing the established v4
solver/grid_lab and complete new construction chain. Every checked `which`
resolved within that source; all native reads used10 threads. Code Analyzer
reported zero issues. MATLAB exited0 with no elliptic LU or PDE work.

| Target | Fine cells | Fraction | Frozen result | X/Y core | Left front |
|---|---:|---:|---|---|---:|
| 42 | 64 | .65 | passed | 46.1503 / 48.3000 | 55.3865 |
| 42 | 64 | .50 | passed | 46.1826 / 48.3000 | 51.3498 |
| 42 | 96 | .65 | passed | 46.1904 / 48.3000 | 56.2272 |
| 42 | 96 | .50 | passed | 46.1830 / 48.3000 | 55.5339 |
| 42 | 128 | .65/.50 | rejected x ratio | — | — |
| 56 | 64 | .65 | passed | 60.5169 / 64.4000 | 67.8910 |
| 56 | 64 | .50 | passed | 61.5287 / 64.4000 | 56.2603 |
| 56 | 96 | .65 | passed | 61.6043 / 64.4000 | 73.4938 |
| 56 | 96 | .50 | rejected x ratio | — | — |
| 56 | 128 | .65/.50 | rejected x ratio | — | — |

All five rejected axes exceed the unchanged1.08 adjacent-cell ratio gate:
target42/fine128 gives1.089123/1.088498; target56/fine96/.50 gives1.080154;
target56/fine128 gives1.090886/1.091309. Transfer is not attempted after
axis rejection. Full rejection records and scored candidates are retained.

Selected actual-unit candidate files, relative to the screen directory:

- `target42_fine64_fraction65/candidate.mat`: peak jump8.24951e-5,
  x/y ratio1.068891/1.067705, global q ratios4.67405e-6/1.51146e-7.
- `target56_fine64_fraction65/candidate.mat`: peak jump8.28131e-5,
  x/y ratio1.070960/1.078945, global q ratios3.50546e-6/4.57611e-8.
  The y-axis is nearer its smoothness and q-floor limits.

All21 admitted-axis construction unit tests passed. Maximum normalized
axis discrepancy was2.65e-16 and maximum dimensionless quality discrepancy
5.85e-12. All six complete selected paired unit tests passed unchanged
frozen gates; maximum normalized axis difference was7.75e-13 and paired
dimensionless metric difference2.03e-12. The tolerances were registered
before testing as5e-11 and1e-8 respectively.

These are candidates for future actual same-node transactions. Actual
post-transfer core, temporal accuracy and long-interval spatial sensitivity
have not been inferred from frozen scores. A future transaction should use
the actual step139 source, the original target-dependent native core floors
(41 or55), front20, original field/peak/rate gates and separately registered
continuous-width observation alongside the unchanged legacy width decision.

## Reproducibility and retained failures

The v3 parent directory contains source hashes, the full frozen dependency
tree and the exact runner. Its log is
`result/longtime/20260908_fresh_tau724_finer_family_v3.log`.
Earlier v1 stopped at an unused-variable Code Analyzer finding, before any
candidate. Earlier v2 produced a valid target42 candidate and then stopped
at MATLAB empty-struct array bookkeeping; both failures remain unchanged.
The fixes affected only this independent scan helper. The construction chain
is `mesh_design_fresh_core_patch` → `mesh_fresh_core_patch_candidate` →
`mesh_core_patch_axis` → `mesh_general_rounded_axis`. The last helper has no
further external research or numerical dependencies: it uses MATLAB
builtins and its local growth/geometric/bridge functions.

## Actual longer target42 comparison

The coordinator subsequently allocated one finer branch from actual step139
to the already completed coarse step328 absolute time1.964631495951036.
`mesh_run_fresh_resolution_branch` froze133 source files, checked dependency
paths, registered every gate and a separate continuous-width/Hermite
auxiliary verdict before computation. The coarse reference was not rerun.

Output: `fresh_tau724_target42_long_v1/run` under the campaign root.
The actual target42 native transaction passed core46.1381/47.8213,
front54.9836 and peak jump8.24951e-5. The fine branch advanced267 steps in
697.713s to step406 and the exact registered absolute time. Its final
core40.6851/45.3495 and safety.196632 passed the predeclared34/.30 gates;
full native pairing, history prefix, lineage, references and fixed axes passed.

The original linear-field/rate/peak gates passed. Worst protocol linear
rhoX L2/Inf were1.697e-3/9.877e-4, rhoY L2/Inf1.409e-3/7.249e-4;
G/P relative differences were4.768e-5/4.087e-6, and parent c_l/c_omega
relative differences1.018e-5/9.311e-5. The original legacy decision remains
false because nodal y-width changed2.1216%; it already differed2.5776%
at the zero-time transaction. Continuous same-peak wall/y widths differed
4.009e-4/1.277e-5 at the terminal comparison, within the unchanged.005
auxiliary bound. Both directions of Hermite field comparisons passed, with
worst rhoX Inf4.470e-4. The separately registered combined auxiliary verdict
is true. These are late-state spatial plus transfer/evolution sensitivities,
not independent-initial-grid convergence or a claim of generic fourth-order
continuous geometry.

Fine `run/true_rhs_cache.mat` was saved from one post-run native restore.
After the fine process exited0 and released its factor, the existing coarse
step328 was restored once in the separate
`fresh_tau744_baseline_rhs_cache_v1` directory. Its `true_rhs_cache.mat` also
contains actual rho/rhoRate/full flow/config/rescaling and an explicit
transport-array allowlist; no elliptic factor is serialized. Baseline
native review and cache audit are alongside it. The two actual caches share
the exact absolute physical epoch and permit later no-LU C1 derivative
sensitivity checks. The baseline factor additionally supported two original
flow evaluations for the independently tiny-tested performance profiler;
all rhoRate/scaleRate/full flow/cache/state values matched exactly. Both
MATLAB processes exited0 and released all extra LU memory.

## Registered third level for unresolved derivative sensitivity

The two actual RHS caches exposed a larger spatial sensitivity in the shape
derivative than in the shape itself. At absolute time1.964631495951036, the
coordinator measured coarse/42 differences of32.53% in the core and22.40%
in the holdout for the normalized shape time derivative, versus.0120% and
.01795% for the shape. These observations motivate the third level; they
do not establish residual convergence.

`fresh_tau724_target56_long_v1` freezes a separate133-file source and the
already admitted target56/fine64/.65 candidate. The source is the same
step139, and the physical target is the actual coarse step328 endpoint.
Before computation it registered actual transaction core54/54 and front20,
terminal core44/44 and safety<.23, at most512 added steps, unchanged CFL,
all original axis/transfer/field/rate/peak gates, and separate legacy and
continuous-width/Hermite decisions. These actual protocol floors supersede
the earlier unrun41/55 suggestion above; none was changed after seeing a
transaction. The y-axis ratio1.078944679 and global q ratio4.576114e-8 retain
the original1.08 and1e-8 limits.

The actual transaction passed core60.5010/64.0774, front67.4781 and peak
jump8.281305e-5. The branch completed355 steps in925.387s to native step494,
with absolute physical time1.964631495951036 exactly. Final core53.2867 /
62.1355 and safety.150131 passed the predeclared gates. The entire history,
native/result pairing, source prefix, lineage, runtime references and fixed
axes passed. MATLAB exited0 and the extra LU allocation was released.

The original legacy comparison and the registered combined auxiliary
comparison both remain **false**. Besides the legacy y-width difference
1.74062%, protocol linear rhoX relative L2 is.00204027136, above its unchanged
.002 gate. The other protocol field metrics pass: rho L2/Inf.00045369 /
.00068512, rhoX Inf.00114474, rhoY L2/Inf.00171024 / .00090477. Peak and rate
gates pass. G/P relative differences are4.88684e-5 / 3.61983e-6; parent
c_l/c_omega differences are9.97020e-7 / 6.64969e-5 in absolute units.

The separate continuous same-peak wall/y width differences are1.85296e-4 /
1.44315e-5, within.005, and both directions of Hermite field comparisons
pass. Their worst gradient relative L2/Inf are6.95678e-4 / 5.31344e-4.
These observations do not replace the failed original linear field gate;
`final_verdicts.json` and `continuous_auxiliary_verdict.json` preserve that
distinction. They also do not decide convergence of the shape derivative.

The true RHS cache was saved from one11.338s post-run native restore with
zero extra time steps. It contains actual rho, rhoRate, complete flow,
configuration, rescaling and an explicit transport array allowlist alongside
Omega/FX and provenance, without elliptic factors. `fine_branch/branch_report.json`
contains its strict native review. The three actual-epoch caches are now
available to the coordinator's separate no-LU residual trend analysis.
The coarse reference was not rerun. `completed_cache_inventory.json` is
only an HDF5 inventory check, not a new native signature validation.

The133-file source hash recheck passed. A hash comparison against the
target42 freeze found only the independent protocol wrapper changed;
all production solver and numerical helper files are identical. Exact
registration, source hashes, comparisons, runner and logs remain in the
new target56 output tree, leaving earlier false outcomes and frozen sources
unchanged.
