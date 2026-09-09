# IPM solver architecture

> Maintained-copy note (2026-09-05): this file preserves the source solver's
> architecture evidence and historical names. It is no longer the current
> executable source of truth. The structured copy uses config/checkpoint schema
> 4 and `scalingContract='exact_gauge_no_feedback_v1'`; rate limiting and width,
> amplitude, or travel restoring feedback described by older passages is
> retired. See [`../../STRUCTURE.md`](../../STRUCTURE.md).

This is the source of truth for the maintained MATLAB solver. It records
executable ownership and numerical contracts, not development history or
empirical evidence. History belongs in [`CHANGELOG.md`](CHANGELOG.md);
fourth- and sixth-order evidence and exclusions belong in
[`HIGH_ORDER_RESEARCH.md`](HIGH_ORDER_RESEARCH.md) and
[`SIXTH_ORDER_RESEARCH.md`](SIXTH_ORDER_RESEARCH.md).

## Hard invariants

- `main_ipm` is the only production solver entry point.
- Configuration is resolved exactly once per run. Internal code never rebuilds
  flat options or invokes defaults again.
- Physical, isotropic dynamic, and anisotropic dynamic modes share one
  elliptic, transport, exact-gauge, Runge--Kutta, remesh, and result path.
- Isotropic scaling is the exact specialization `C_x=C_y`.
- The resolver baseline remains the second-order regression oracle. Fourth
  order is explicit opt-in; sixth order additionally requires its complete
  four-selector tuple. Neither method is inferred from a single `order` switch.
- Adding an opt-in method cannot alter default second-order or existing
  fourth-order arithmetic. Shared-path changes require cross-revision bitwise
  regression of both established tuples.
- A remesh is a transaction: a candidate cannot mutate the accepted state
  before validation succeeds.
- Numerical claims require the verification gates in this document.

## Ownership and lifecycle

```text
flat public overrides (or ipm_active_case)
  -> ipm_initialize_run
       -> ipm_resolve_config                 resolve exactly once
       -> ipm_build_operators_resolved       build runtime operators
       -> ipm_initial_density
       -> ipm_initialize_rescaling
  -> main_ipm loop
       -> ipm_select_timestep
       -> ipm_advance_one_step
            -> ipm_step_ssprk3 / ipm_step_ssprk54 / ipm_step_rk6
                 -> ipm_step_rk              shared stage-state contract
                      -> ipm_flow_at_state
                           -> ipm_rhs
                                -> ipm_velocity
                                     -> ipm_poisson / ipm_green_boundary
                                -> feature tracking and gauge selection
                                -> ipm_assemble_rescaled_rhs
            -> ipm_remesh_if_needed
                 -> propose -> transfer candidate -> validate -> commit/retry
       -> ipm_record_state
       -> ipm_stop_policy
  -> ipm_finalize_result
       -> ipm_reconstruct_physical
       -> ipm_validate_result
  -> ipm_write_result / ipm_report_result
```

`main_ipm.m` owns control flow only. Numerical formulas stay in leaf modules,
so adding a mode cannot silently create another solver. After a step or remesh
proposal returns normally, `ipm_state_is_finite` checks the evolution fields
and logarithmic scale state. A non-finite proposal is rolled back to the last
accepted state and returns `non_finite_solution`; exceptions from numerical
modules propagate. Physical reconstruction is performed later during recording
and finalization.

Four structures have distinct ownership:

| Structure | Owner | Contract |
|---|---|---|
| `config` | `ipm_resolve_config` | Frozen, grouped, complete; persisted in the result. |
| `ops` | `ipm_build_operators_resolved` | Runtime grid/numerical data only; never stores the full config. |
| `state` | lifecycle | Current accepted `rho`, flow, scale, clocks, step, and the frozen config. |
| `result` | `ipm_finalize_result` | Versioned external data contract; contains no live solver objects. |

## Configuration and numerical selectors

`ipm_option_schema` is the sole inventory of public option names, domains,
raw defaults, types, and allowed values. `ipm_resolve_config` performs

```text
defaults -> overrides -> normalization -> derivation -> validation -> freeze
```

and returns nine domains: `grid`, `time`, `physics`, `elliptic`, `transport`,
`scaling`, `remesh`, `diagnostics`, and `output`. Unknown and removed fields
are errors. Remeshing passes a local grid override directly to the pure
operator builder; it does not resolve configuration again.

Four commonly confused configurations are deliberately different:

| Meaning | Spatial | Transport | Time | Remap |
|---|---|---|---|---|
| Resolver baseline | `legacy_second_order` | `muscl_minmod` | `ssprk3` | `pchip` |
| No-argument active case | baseline spatial | `weno5_nonuniform` | `ssprk3` | `pchip` |
| Certified fourth-order tuple | `high_order` | `weno5_fd` | `ssprk54` | fixed grid; high-order transfer tested separately |
| Sixth-order tuple | `sixth_order` | `weno7_fd` | `rk6` | `sixth_order` |

The no-argument active case also selects dynamic rescaling, a `513 x 257`
large-domain grid, analytic/adaptive remeshing, and normal output behavior. It
is a production case, not the schema baseline or a safe smoke test.

The four selectors remain independent for baseline/fourth-order controlled
experiments. Valid does not imply a global-order certificate. Sixth order is
different: its spatial, transport, time, and remap selectors form one strict
tuple even on a fixed grid, so a later remesh cannot silently reduce order.
The resolver enforces:

- `weno5_fd` requires `spatialDiscretization='high_order'`;
- `remeshTransferScheme='high_order'` requires high-order spatial operators;
- `sixth_order` spatial requires `weno7_fd`, `rk6`, and `sixth_order` remap,
  and each of those three selectors requires `sixth_order` spatial;
- both full-order spatial backends require the direct anisotropic solver, not
  PCG.

`ipm_numerical_tuple_violation` owns all coupled-selector rules for both
run-time resolution and persisted-result validation. The result validator
checks the tuple without applying defaults or resolving configuration a second
time. `ipm_sixth_order_options` is the preferred explicit constructor for the
complete sixth-order selection; it does not bypass resolution.

## Coordinates and scale state

The dynamic map is

```text
X = C_x(t) x + X_shift(t),
Y = C_y(t) y,
rho_physical = rho_rescaled / C_omega(t).
```

Positive scales are stored logarithmically. Isotropic state stores
`logC_l=log(C_x)` and has `C_y=C_x`; anisotropic state additionally stores
`logC_y`. Both store `logC_omega`, `X_shift`, physical time, and canonical
time. The anisotropic elliptic metric is

```text
kappa = C_y / C_x.
```

`ipm_time_speed` accepts the complete rate vector, translation scale, and cap.
It returns a positive time speed and the canonical rate magnitude. The reported
conservative source coefficient is always

```text
c_x + c_y + c_omega.
```

Physical mode sets every scale rate exactly to zero and does not require an
origin node or a positive wall feature. Dynamic isotropic rates are owned by
`ipm_isotropic_gauge`. Dynamic anisotropic mode supports only `fixed` and
`peak_translation`; the latter requires half-plane symmetry. Every
anisotropic RK stage uses its own current `kappa`.

## Runtime numerical contract

Arrays are `ny x nx`: rows are `y`, columns are `x`. A production `ops`
structure is complete, not a partial compatibility object. In particular it
always supplies geometry, derivative matrices, boundary modes,
`spatialDiscretization`, `dynamicScaleGeometry`, and `integrationWeights`.
The legacy integration weights equal `ops.weights` exactly. Backend-only data
such as `poissonWeights` exist only when that backend uses them.

### Spatial and elliptic backends

| Backend | Derivatives and quadrature | Elliptic solve |
|---|---|---|
| `legacy_second_order` | Existing mapped/nonuniform second-order operators and control-volume weights | Weighted SPD operator, cached Cholesky at `kappa=1`, direct or PCG anisotropic solve |
| `high_order` | Seven-point local-polynomial derivatives and degree-five cell quadrature | Full-order one-sided closure, nonsymmetric operator, LU/direct solve |
| `sixth_order` | Nine-point local-polynomial derivatives and positive mapped degree-seven quadrature | Paired full-order one-sided closure, nonsymmetric operator, LU/direct solve |

`ipm_poisson_operator(ops,kappa)` is the single metric operator. The exact
`kappa=1` legacy branch retains its original boundary updates, weighted RHS,
and arithmetic grouping. Green geometry uses

```text
kappa * dx^2 + dy^2 / kappa.
```

Both full-order Green paths share `integrationWeights`, but their global-order
claims are restricted to the smooth uncompressed source classes documented in
`HIGH_ORDER_RESEARCH.md` and `SIXTH_ORDER_RESEARCH.md`.

Before a sixth-order operator is accepted, each axis is measured by
`ipm_grid_quality` and checked against the private, fixed policy returned by
`ipm_sixth_order_grid_policy`. The current limits are adjacent-cell ratio
`1.15`, normalized nine-point stencil `rcond` floor `1e-9`, normalized positive
quadrature-weight floor `1e-4`, and mapping-metric consistency error ceiling
`1e-6`. These are internal numerical admission rules, not public option fields
and not a theorem for arbitrary nonuniform grids. Automatic sinh stretch is
capped before construction; explicit/custom/remeshed axes fail explicitly.

`ipm_velocity` computes `rho_x`, solves Poisson, and differentiates the stream
function to produce `u=(-psi_y,psi_x)`. It is common to every mode and backend.

### Transport form and source assembly

`ipm_transport_rhs` owns the discrete transport form and returns a descriptor
alongside the RHS:

- `muscl_minmod` and `weno5_nonuniform` return a finite-volume conservative
  bulk divergence. `weno5_nonuniform` uses WENO only where its full stencil is
  available and falls back to MUSCL in two boundary layers.
- `weno5_fd` returns mapped finite-difference WENO-Z5 in advective form. It
  uses linewise-global Lax--Friedrichs splitting, degree-five ghost
  extrapolation, and a discrete constant-state/free-stream correction.
- `weno7_fd` returns the paired mapped finite-difference WENO-Z7 advective
  form. It uses four cubic substencils, WENO-Z power four, linewise-global
  Lax--Friedrichs splitting, degree-six ghost continuation, and the same
  constant-state/free-stream and closed-mass contracts.
- `advective_upwind` makes the physical-wall row advective even when the bulk
  scheme is conservative. The default `conservative_flux` leaves the wall in
  the bulk form.

`ipm_assemble_rescaled_rhs` consumes the returned form descriptor; it never
guesses semantics from a scheme name. For a finite-volume bulk it adds
`(c_x+c_y)rho` exactly once to convert conservative divergence to the desired
advective spatial form. WENO-FD already performed that conversion through its
free-stream correction, so it receives no second spatial source. Both forms
then add `c_omega rho` exactly once. An advective wall row is excluded from the
finite-volume spatial correction.

For closed WENO-FD with the conservative wall, the residual quadrature mass
defect is projected onto the constant mode. The ordinary physical target is
zero. An explicit closed sixth-order rescaled assembly instead supplies the
source-aware target `(c_x+c_y)M`, so the complete RHS has weighted rate
`(c_x+c_y+c_omega)M`. Production dynamic gauges retain their established open
coordinate-transport assembly. The explicit advective wall trace is
authoritative and is not modified by any projection.

### Time integration

`ipm_step_ssprk3`, `ipm_step_ssprk54`, and `ipm_step_rk6` are thin stable entry
points. `ipm_step_rk` owns the single scale-vector encoding, stage RHS,
stage-local metric, and final flow reconstruction. It deliberately keeps the
SSPRK3 and SSPRK(5,4) update expressions separate and explicit; converting
them to a generic stage loop or low-storage formula could change the
second-order roundoff oracle. Only the eight-stage Dormand--Prince RK6(5)8M
path uses the generic explicit-tableau stage loop. Its sixth-order row advances
the state, its embedded fifth-order row is diagnostic only, and its SSP
coefficient is zero.

The lifecycle dispatches only from `config.time.timeIntegrator`. There is no
geometry-specific integrator and no duplicate evolution chain.

## Feature and gauge ownership

`ipm_measure_features` is a policy-free one-dimensional measurement primitive.
`ipm_initialize_rescaling` owns initial reference peaks, widths, and locations.
At runtime, `ipm_track_features` owns wall/vertical measurements and passes one
record to the selected gauge. Gauge modules choose rates; they do not duplicate
elliptic, transport, or time integration.

## Remesh transaction

Adaptive remeshing follows one atomic sequence:

```text
trigger
  -> propose grid
  -> validate candidate geometry, transfer field, and build candidate operators
  -> validate core counts, physical peak, and physical-density range
  -> accept: increment remeshCount once and install candidate
     reject: retain rho and ops exactly, damp proposal, and retry
```

Retries use indices `0:8`; a retry damps the original proposal amplitude by
`0.8^retry`. Only the following known numerical grid failures are recoverable:
`HighOrderQuadratureWeights`, `HighOrderGridMetric`,
`FiniteDifferenceScale`, `HighOrderRemeshBubble`, `SixthOrderGridMetric`,
`SixthOrderQuadratureWeights`, `SixthOrderQuadratureMetric`,
`SixthOrderGridCellRatio`,
`SixthOrderGridStencilCondition`, `SixthOrderGridQuadratureMargin`, and
`SixthOrderGridMetricConsistency` (all under the `ipm:` namespace). Other
exceptions propagate because they indicate programming or configuration
errors. Candidate operators retain the current `remeshCount`; only the
orchestrator commits the increment after validation.

The fourth-order transfer uses local six-point barycentric interpolation and
degree-five quadrature. The sixth-order transfer uses width-eight,
degree-seven interpolation and the exact old/new one-dimensional factors of
the solver's paired mapped tensor norm. Both use a
horizontal constant correction and wall-normal bubble correction. Range
overshoot is measured in physical-density units. A rejected candidate never
changes the accepted state.

Initial analytic remeshing uses the same grid transaction, but after accepting
a grid it re-evaluates the analytic datum there instead of installing an
interpolated field. Threshold-triggered adaptive evolution is not currently a
global fourth- or sixth-order claim.

## Result contract

The external result has top-level `schemaVersion=2` and grouped metadata,
state, grid, scale, physical fields, histories, quality, fit, configuration,
elliptic metadata, and snapshots. Disabled history groups stay empty rather
than allocating bulk `NaN` arrays. Every nonempty history group is aligned to
`history.common`.

New nested `result.config` uses configuration schema 4, records all four
numerical selectors, and fixes
`scalingContract='exact_gauge_no_feedback_v1'`. It contains neither
`maxDynamicRate` nor the retired width controller. The reader also accepts
otherwise valid result-v2 files with nested config schema 1--3 for analysis.
For schema 1, absence of `timeIntegrator`,
`spatialDiscretization`, and `remeshTransferScheme` is interpreted during
validation as SSPRK3, legacy spatial operators, and PCHIP. The structure is
returned unchanged; no missing fields are injected and no schema is upgraded.
A stored schema-1 `transportScheme` remains authoritative, but cannot claim
the later `weno5_fd` or `weno7_fd` paths. Legacy config schemas remain readable
for analysis but are not convertible to current schema-4 checkpoints.

MAT input must contain a variable named `result`. Missing top-level versions,
result schema 1, aliases, and inferred variables are rejected. Requested
output paths remain frozen in config; collision-free actual paths live in
metadata. Persistence never overwrites an existing MAT or video artifact.

## Verification gates

The runners are separate by design:

- `verify_ipm()` runs the 11 maintained layers: config, grid, elliptic,
  numeric core, transport, symmetry, scaling, remesh, smoke, output, and
  integration. It does not run the fourth-order certificate.
- `verify_ipm_high_order_all()` runs high-order spatial/coupled, temporal,
  unforced physical convergence, Green, remap, and stability gates.
- `verify_ipm_sixth_order_all()` independently runs sixth-order core/WENO,
  including two-tail D1/D2/nonzero-Dirichlet Poisson/velocity refinement and
  both WENO LF splits; RK algebra/nonlinear-ODE and nonzero-rate anisotropic
  production-stage coupling; width-eight interpolation/conservation, grid
  admission, analytic-sinh integrated remap conservation, invalid-candidate
  exact rollback, and mass-conservative unforced uniform/stretched physical
  convergence.
- `verify_ipm_release()` deliberately retains its established second/fourth
  semantics. Sixth-order release work invokes its independent all-runner in
  addition. The broader Green/repeated-remap/stability research campaigns remain
  recorded separately in `SIXTH_ORDER_RESEARCH.md`.

Every numerical or architectural change must pass MATLAB Code Analyzer,
`git diff --check`, `verify_ipm()`, and the relevant fourth-/sixth-order layers.
A release touching shared numerical code must pass
`verify_ipm_release('heavy')` and, when the sixth-order feature is affected,
`verify_ipm_sixth_order_all()`. Shared-path refactors additionally
require cross-revision bitwise comparisons of the established second- and
fourth-order accepted states and results.

Static success never substitutes for a numerical run. If MATLAB cannot
execute, numerical verification must be reported as unexecuted.

## File ownership and change rules

- Lifecycle: `main_ipm`, initialization, step selection, advance, record, and
  stop-policy modules.
- Configuration: `ipm_option_schema`, `ipm_resolve_config`,
  `ipm_numerical_tuple_violation`, the opt-in
  `ipm_sixth_order_options`, and private `ipm_sixth_order_grid_policy`.
- Operators/elliptic: resolved builder, `ipm_grid_quality`, derivative and
  fourth-/sixth-order quadrature helpers, Poisson, Green, and velocity modules.
- Dynamics: RHS, transport/source assembly, SSPRK/RK6 entry points and shared
  RK core, scale-speed logic, and WENO5/WENO7 workers.
- Features/gauges: measurement, runtime tracking, isotropic gauge, and the two
  maintained anisotropic gauges.
- Mesh: trigger, orchestrator, proposal, transfer, validation, adaptive-axis,
  and width-six/width-eight high-order interpolation modules.
- Certification: `verify_ipm`, `verify_ipm_high_order_*`, and the independent
  `verify_ipm_sixth_order_*` runners; empirical narratives live in the two
  research-certificate documents rather than this contract.
- Output: record/reconstruction/finalization, result validation, persistence,
  reporting, and visualization.

When a difference can be represented by a metric, rate vector, boundary mode,
or callback, extend the common primitive. Do not add a parallel solver, a
second configuration representation, an alias field for the same numerical
quantity, or a compatibility fallback for incomplete runtime `ops`.
