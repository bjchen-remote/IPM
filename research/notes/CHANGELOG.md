# Changelog

> Archive note (2026-09-05): this chronology stops before the structured
> solver's config/checkpoint schema-4
> `exact_gauge_no_feedback_v1` contract. Any `maxDynamicRate`, adaptive-width,
> amplitude-restoring, width-restoring, or travel-restoring behavior below is a
> retired historical contract. Current changes are recorded in
> [`../../CHANGELOG.md`](../../CHANGELOG.md).

Implementation chronology before 2026-08-30 is preserved in
[`analysis/docs/ARCHITECTURE_RESEARCH_AND_HISTORY_2026.md`](analysis/docs/ARCHITECTURE_RESEARCH_AND_HISTORY_2026.md).
This file records maintained-program changes only.

## 2026-08-30

### Integrated explicit sixth-order feature

- Ported the isolated sixth-order method from research commit `15a10d9` onto
  the grouped-configuration shared solver as an explicit feature, without
  changing the resolver baseline or the established fourth-order selection.
- Added the strict four-selector tuple `sixth_order/weno7_fd/rk6/sixth_order`.
  Each WENO7, RK6, and sixth-order-remap selector is bound bidirectionally to
  sixth-order spatial operators, including for initially fixed-grid runs.
- Added nine-point differential/Poisson closure, positive mapped degree-seven
  quadrature, mapped WENO-Z7 with degree-six ghost continuation, the
  eight-stage Dormand--Prince RK6(5)8M path, and width-eight conservative
  smooth remapping through the existing operator, transport, RK, remesh, and
  result lifecycles.
- Added the private `ipm_sixth_order_grid_policy` and executable axis-quality
  checks for adjacent-cell ratio, stencil conditioning, positive quadrature
  margin, and mapping-metric consistency. Sixth-order metric, quadrature, and
  four admission failures are recoverable remesh-candidate failures; retry and
  rollback remain transactional.
- Added independent integrated sixth-order core/time/physical certificates:
  two-tail D1/D2/nonzero-Dirichlet Poisson/velocity refinement, both WENO LF
  splits, nonzero-rate anisotropic RK6 production-stage coupling,
  source-aware explicit closed rescaled-mass assembly,
  analytic-sinh remap mass/trace preservation, inadmissible-candidate exact
  rollback, and closed-boundary physical mass drift. Documented the broader
  isolated Green/remap/stability evidence, two
  noncertifying pure-time switching diagnostics, and strict exclusions in
  `SIXTH_ORDER_RESEARCH.md`. `verify_ipm_release` retains its existing
  second/fourth-order meaning; sixth order has a separate all-runner.
- Kept the continuing 12-hour grid/CFL optimizer on
  `research/sixth-order-global6`. No incomplete checkpoint ranking is recorded
  as the maintained optimum; the integrated helper currently recommends
  `CFL=0.40` from the completed finite stability sweep.
- Passed the unchanged 11-layer `verify_ipm()` suite, the complete existing
  fourth-order quick certificate, the new integrated sixth-order certificate,
  a 29-file zero-warning MATLAB Code Analyzer run, and cross-revision exact
  result comparisons for representative second- and fourth-order nonuniform
  runs.

### Shared-core simplification and architecture contract

- Merged the duplicate SSPRK3/SSPRK(5,4) scale packing, stage RHS, and final
  flow reconstruction into `ipm_step_rk`; each tableau retains its original
  explicit arithmetic expression and its stable entry point.
- Made `integrationWeights` a mandatory runtime-operator field, removed the
  redundant Green-weight alias and unused high-order Poisson weights, and
  replaced the private remap quadrature copy with the public high-order rule.
- Centralized the three coupled high-order selector constraints in one pure
  validator shared by configuration resolution and persisted-result checking.
- Made transport report its conservative/advective form to source assembly,
  eliminating scheme-name inference and locking the rule that the spatial
  divergence source is added exactly once.
- Clarified remeshing as a candidate/validate/commit transaction: candidate
  construction preserves the accepted remesh count, and only the orchestrator
  increments it after acceptance. Removed unused candidate retry state and
  redundant grid overrides.
- Rewrote the architecture, quick-start, and fourth-order certification
  documents to distinguish resolver baseline, no-argument active case, valid
  numerical tuples, and the narrower globally certified fourth-order scope.
- Confirmed cross-revision bitwise identity for legacy physical and dynamic
  evolution, accepted PCHIP remeshing, and the following SSPRK3 step; every
  compared field had maximum numerical difference zero.
- Passed `verify_ipm_release('heavy')` after simplification, including the
  480-step physical stability run, and obtained zero Code Analyzer findings
  across all 108 MATLAB files.

### Integrated opt-in fourth-order path

- Added an opt-in, end-to-end smooth-solution path combining seven-point local
  polynomial derivatives and Poisson closure, degree-five cell quadrature,
  mapped finite-difference WENO-Z5 transport, and five-stage fourth-order
  SSPRK(5,4). The maintained second-order numerical path remains the default
  and retains its prior regression gates.
- Ported the research prototype from `cf13af6` onto the grouped-configuration
  solver after `b5837f1` without restoring deleted flat-option wrappers.
  Spatial, transport, time, and remesh order choices are independent explicit
  options validated by `ipm_resolve_config`.
- Advanced the nested frozen configuration to schema version 2 so new results
  persist those numerical choices. The version-2 result reader retains a
  narrow compatibility path for earlier nested config-version-1 results.
- Added a direct nonsymmetric high-order elliptic solve with full boundary
  coupling and high-order Green weights. PCG is rejected for this prototype
  instead of silently applying an incompatible symmetric solver.
- Added six-point degree-five local barycentric remapping with explicit
  tensor-mass and column-mass corrections that preserve the physical wall
  and top trace.
- Added independent gates for all-node spatial convergence, genuine temporal
  refinement including actual nonzero-scale isotropic and anisotropic
  SSPRK(5,4) production steps, coupled manufactured evolution, unforced
  production self-convergence, Green quadrature, repeated and integrated
  remapping, CFL limits, nonsmooth stress profiles, and long physical-IPM
  evolution. See `HIGH_ORDER_RESEARCH.md` for measured evidence and the
  deliberately uncertified scope.

### Maintained solver refactor

- Declared the logarithmic-primitive `degenerate_primitive`, power-eight datum
  as the no-argument active case and removed the bounded-rational wording.
- Kept the former CCF heavy-tail datum as an explicit named experiment so its
  long dynamic case no longer changes when the active default changes.
- Added one-pass option resolution and a pure resolved-operator builder.
- Unified isotropic and anisotropic Green, Poisson, velocity, time-speed,
  rescaled transport, and SSPRK3 primitives. Preserved explicit `kappa=1` and
  zero-rate oracle paths.
- Removed the incomplete triangular anisotropic research prototype. The
  maintained anisotropic modes are `fixed` and `peak_translation`.
- Split feature extraction and isotropic gauge selection into independent
  modules.
- Reduced `main_ipm` to lifecycle orchestration and extracted initialization,
  step selection, advance/remesh policy, state recording, stopping, physical
  reconstruction, result finalization, persistence, plotting, and reporting.
- Introduced result schema version 2 with grouped history; disabled features no
  longer allocate bulk `NaN` history fields.
- Added per-result case metadata, collision-free atomic result persistence,
  append-only JSON-lines manifests, and non-overwriting video paths.
- Decoupled physical evolution from dynamic-gauge origin and peak assumptions;
  peak-resolution history and terminal quality are optional when no feature
  exists. K1 Bessel cases therefore no longer add an artificial wall monitor
  and explicitly disable their inapplicable positive-wall sign gate.
- Split verification into configuration, grid, elliptic, numeric-core,
  transport, symmetry, scaling, remesh, smoke, output, and integration layers.
- Moved experiment-specific Bessel parameters into structured case metadata.
- Archived the former combined architecture/research/change log without
  deleting it and replaced it with the current executable contract.

### Verification status

- The complete `verify_ipm` suite passed under MATLAB R2026a, including all
  five short solver modes, two physical edge cases, and the output/manifest
  contract.
- Exact isotropic-oracle errors were zero for Green, metric operator, and
  SSPRK3; direct and PCG Poisson errors were below `4e-15`.
- WENO reconstruction order was `4.701`, the closed-boundary mass rate was at
  roundoff, and the physical/active/half-plane/anisotropic-fixed/
  peak-translation smoke runs completed in `2/3/3/3/3` steps; the asymmetric
  origin-free and peak-free physical cases completed in `2/2` steps.
- MATLAB Code Analyzer reported no warnings across maintained root, analysis,
  and experiment MATLAB files; `git diff --check` passed.

### V2-only simplification batch

- Retired the transitional `ipm_upgrade_result_v1`,
  `ipm_result_legacy_view`, and `ipm_history_legacy_view` paths after migrating
  maintained analysis and verification code to the grouped version-2 result.
  `ipm_validate_result` now rejects missing, version-1, and other schema tags;
  MAT input must contain a variable named `result`.
- Removed the flat configuration and construction wrappers
  `ipm_config_legacy_view`, `ipm_defaults`, `ipm_build_operators`, and
  `ipm_rescaling_rates`.
- Consolidated the remaining configuration layer into the declarative
  `ipm_option_schema` and single `ipm_resolve_config` boundary; removed
  `ipm_resolve_options` and `ipm_group_options`. The public entry still accepts
  flat overrides, but production and verification code use only the resulting
  nine-domain frozen configuration, with no `state.opts` copy.
- Changed remesh reconstruction to pass a local grid override to the pure
  grouped-config operator builder instead of mutating or re-resolving options.
  The adaptive-axis amplitude cap is now an explicit algorithm parameter.
- Kept requested result/video paths frozen in `result.config.output` and actual
  collision-free artifact paths in `result.metadata`. Result validation derives
  its exact configuration field sets from the declarative schema without
  invoking resolution or defaults.
- Removed the positional time-speed test signature. `ipm_time_speed` now has
  only the common rate-vector API used by production and verification callers.
- Centralized nonnegative profile masking, peaks, widths, and resolution counts
  in the pure `ipm_measure_features` primitive shared by rescaling
  initialization and runtime tracking. Flattened the internal tracker-to-gauge
  measurement record without changing the result/flow schema, crossing rules,
  or evaluation order. The maintained solver fell from 5,716 to 5,653 lines
  in this batch; the measurement core itself fell from 378 to 317 lines.
- Added direct masked-profile and wall-omega-width gauge regressions to the
  scaling verification layer.
- Deduplicated the two identical fixed-iteration searches that retain the
  strongest valid adaptive-axis fallback. The count-only bracket, final
  valid-and-count refinement, 32-evaluation order, midpoint arithmetic, and
  strict best-score comparison remain unchanged. `ipm_adaptive_axis` fell from
  433 to 407 lines.
- Removed the unused 20-file, 3,203-line `legacy_boussinesq` solver tree and
  the three-line `practical_main` convenience wrapper. `main_ipm` is now the
  sole executable entry advertised by the active tree; the retired sources
  remain recoverable from Git and `safety/pre-refactor-20260830`.
- Removed the unused `degenerate_density` compatibility datum, the Bessel
  experiment's flat metadata aliases, and the final stale
  `ipm_resolve_options` error-source labels. The terminal audit confirmed one
  configuration resolver, one rate-vector time-speed API, result schema 2, and
  no deleted adapter names in tracked executable code.
- Across the complete simplification, maintained solver code fell from 6,062
  to 5,613 lines and all tracked MATLAB from 15,144 to 11,769 lines.
  Verification grew from 2,666 to 2,992 lines to preserve v2 rejection,
  feature-measurement, width-gauge, and maintained-mode regressions.
