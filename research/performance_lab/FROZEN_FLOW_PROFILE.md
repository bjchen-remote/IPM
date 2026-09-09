# Existing-cache flow profile

`ipm_perflab_profile_frozen_state(state, newOutputDirectory)` accepts an already
restored native state with an accepted-state RHS cache. It performs two original
`ipm.evolve.flow` evaluations: one ordinary wall-clock measurement and one MATLAB
function/line profile. It does not restore a checkpoint, build an operator,
factor a matrix, advance the PDE, or write a checkpoint. Each flow must make
exactly one native `assembleRhs` and one cached decomposition solve.

The helper requires 10 computational threads and the maintained isotropic
`transport_anchor` / `wall_omega_quadratic_peak` path. It refuses to interrupt
an active MATLAB profiler. A prior inactive profiler buffer, if present, is
saved before collecting this experiment. The output directory must be new.

Both returned density rates, derived scale rates, complete flow structures and
complete RHS caches must match the caller's cached values with `isequaln`.
The complete input state is also checked for value equality. The test does not
claim to inspect undocumented mutable state inside MATLAB's decomposition object.

The report partitions the profiled flow time into disjoint call-graph costs:
cached decomposition backslash, Green boundary, remaining Poisson work,
remaining velocity work, transport assembly, remaining isotropic gauge work,
feature tracking, remaining RHS work, and the flow wrapper. These differences
are valid for the checked single-call path. Inclusive function and line times
must not be summed across nested calls. `functionsBySelfSeconds` also contains
profiler/harness entries outside the flow; the disjoint partition uses only the
native flow's total time. MATLAB's native sparse solve may appear as self time
of its internal `SparseLU.solve` wrapper rather than as a separate builtin row.

`profile_raw.mat` retains the original function table. `profile_report.json`
and `.mat` retain source paths, all equality checks, call counts, the partition,
and sorted function/line records. Profiling adds overhead, especially to the
many small WENO calls. These times locate costs; they do not establish a speedup
and should not be rescaled mechanically to the separate ordinary timing.

The first test used the existing 49×49 crop checkpoint at fresh step 9 from
`result/verification/performance_lab_20260908/three_box_extension_tiny_v2/positive`.
The harness performed one tiny restore; the helper performed none. Every
equality check passed, the nine required call counts were one, and `checkcode`
reported zero issues. Ordinary flow time was 0.0168325 s and profiled time was
0.0429124583 s, a 2.54938 overhead ratio. This tiny cost partition is not evidence
about the production grid. Records are in
`result/verification/performance_lab_20260908/profile_frozen_tiny_v1/`.

The mesh experiment then froze this helper and called it in its already
scheduled baseline step-328 restore, after saving the cache and before releasing
the same state and LU. That MATLAB process exited normally. No separate large
restore was performed for profiling.

## Actual 895×386 baseline cache

Records are in
`result/longtime/20260908_campaign_v1/fresh_tau744_baseline_rhs_cache_v1/profile/`.
The frozen source is the original v4 source. The fresh clocks at accepted step
328 are canonical 0.010867710115605512 and physical 0.010104359265323784;
these are local fresh clocks, not the parent absolute clocks.

Both complete flow/cache comparisons and the state equality passed. Every
required root function was called once. Ordinary flow took 0.5618838333 s;
profiled flow wall time was 0.6391122917 s, 13.7446% longer. The disjoint native
flow partition sums to 0.6390543745 s:

| Cost | Profiled seconds | Profiled flow share |
| --- | ---: | ---: |
| Cached decomposition backslash | 0.465949447 | 72.9123% |
| Transport assembly | 0.143734009 | 22.4917% |
| Green boundary | 0.015016668 | 2.3498% |
| Remaining Poisson work | 0.003607625 | 0.5645% |
| Remaining velocity work | 0.002099958 | 0.3286% |
| Remaining isotropic gauge work | 0.002741792 | 0.4290% |
| Feature tracking | 0.004087000 | 0.6395% |
| Remaining RHS work | 0.001765000 | 0.2762% |
| Flow wrapper | 0.000052875 | 0.0083% |

The leading hotspot is the one native `SparseLU.solve` call: 0.465506530 s of
self time. Its call path is `poisson.m:41`, `decomposition.m:381`, then
`SparseLU.m:36`; the final line calls MATLAB's internal sparse solver. The
decomposition wrapper itself accounts for only about 0.443 ms beyond that
native solve. There is no repeated factorization in this flow.

The four WENO calls contain 5,132 `weno_z_left` calls. Their profiled inclusive
time is 0.140933426 s; this fine-grained path is particularly exposed to
profiler overhead. The new partition replaces the earlier 88.74% *unassigned*
cost estimate. Neither that old remainder nor the new 72.91% profiled share is
an unprofiled LU percentage. The actual trace establishes the native cached
solve as the dominant cost without qualifying the previously rejected public
LU or dense Schur replacements.
