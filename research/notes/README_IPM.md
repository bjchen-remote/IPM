# IPM solver

> Maintained-copy note (2026-09-05): the structured solver now emits config
> schema 4 with `scalingContract='exact_gauge_no_feedback_v1'` and checkpoint
> schema 4. `maxDynamicRate` and the adaptive width controller are retired;
> `lengthScaleGain=1`, while omega/width/travel restoring gains are zero.
> Feature cell counts are diagnostics/remesh/hard-stop inputs only. Older text
> below is retained as source-repository evidence; where it conflicts, the
> current contract is [`../../STRUCTURE.md`](../../STRUCTURE.md).

This repository contains one maintained MATLAB solver for the two-dimensional
incompressible porous media equation on a truncated upper half-plane:

```text
rho_t + u . grad(rho) = 0,
u = grad^perp psi,
-Delta psi = d_x rho.
```

It supports physical evolution, isotropic dynamic rescaling, a small
anisotropic correctness mode, the maintained second-order method, and explicit
smooth-solution fourth- and sixth-order paths.

## Safe quick start

Use an explicit short, non-writing case when checking an installation:

```matlab
opts = struct( ...
    'nx',65,'ny',33,'xlim',[-4,4],'ymax',4, ...
    'rescalingMode','physical', ...
    'physicalFinalTime',1e-3,'maxDt',1e-4, ...
    'adaptiveRemesh',false,'initialAnalyticRemesh',false, ...
    'storeSnapshots',false,'saveResults',false, ...
    'makePlots',false,'livePlot',false,'writeVideo',false, ...
    'verbose',true);
result = main_ipm(opts);
```

`main_ipm()` with no arguments is not a smoke test. It loads
`ipm_active_case.m`: a `513 x 257` dynamic-rescaling production case on a very
large stretched domain, with analytic/adaptive remeshing and the normal output
settings. Its datum is

```text
rho(X,Y) = log(1 + |X|^8/(1+Y^8)) / 8.
```

This is the power-eight logarithmic primitive, not the former bounded rational
density.

## Numerical choices

The resolver baseline is the unchanged second-order tuple:

```matlab
opts.spatialDiscretization = 'legacy_second_order';
opts.transportScheme = 'muscl_minmod';
opts.timeIntegrator = 'ssprk3';
opts.remeshTransferScheme = 'pchip';
```

The no-argument active case overrides `transportScheme` with
`weno5_nonuniform`. That scheme has fifth-order interior reconstruction but
uses MUSCL in two boundary layers; it is not the globally fourth-order path.

For a currently certified smooth fixed-grid fourth-order run, use the complete
combination explicitly:

```matlab
smoothRho = @(X,Y) 2 + ...
    ((pi^2.*Y.*(1-Y)+2)./pi).*cos(pi.*(X+1));
opts = struct( ...
    'nx',65,'ny',65,'xlim',[-1,1],'ymax',1, ...
    'gridMode','uniform','gridStretchAutomatic',false, ...
    'gridStretch',[0,0], ...
    'spatialDiscretization','high_order', ...
    'transportScheme','weno5_fd', ...
    'timeIntegrator','ssprk54', ...
    'farBoundaryMode','dirichlet_zero', ...
    'transportBoundaryMode','closed', ...
    'wallTransportMode','conservative_flux', ...
    'adaptiveRemesh',false,'initialAnalyticRemesh',false, ...
    'initialCondition',smoothRho, ...
    'finalTime',1e-2,'physicalFinalTime',Inf, ...
    'maxDt',2.5e-4,'maxSteps',42,'cfl',0.45, ...
    'saveResults',false,'storeSnapshots',false, ...
    'makePlots',false,'livePlot',false,'writeVideo',false, ...
    'verbose',true);
result4 = main_ipm(opts);
```

The certified grid may also be smoothly mapped. Start with moderate grids:
the current production convergence campaign extends through `129 x 129`, and
the high-order nonsymmetric LU closure is not intended as the default
`1025 x 513` large-grid solver.

High-order remapping is a separate experimental selector:

```matlab
opts.remeshTransferScheme = 'high_order';
```

Its interpolation, conservation, rejection, retry, and rollback contracts are
tested, but threshold-triggered adaptive evolution is not yet included in the
global-fourth-order claim.

The sixth-order feature is a separate, strict opt-in. Use its helper so all
four numerical selectors and the conservative fixed-grid settings are selected
together:

```matlab
smoothRho = @(X,Y) 2 + ...
    ((pi^2.*Y.*(1-Y)+2)./pi).*cos(pi.*(X+1));
opts6 = ipm_sixth_order_options(struct( ...
    'nx',129,'ny',129,'xlim',[-1,1],'ymax',1, ...
    'gridMode','stretched','gridStretchAutomatic',false, ...
    'gridStretch',[0.75,0.50], ...
    'initialCondition',smoothRho, ...
    'finalTime',0.04,'physicalFinalTime',Inf, ...
    'maxDt',0.04/64,'maxSteps',66, ...
    'saveResults',false,'storeSnapshots',false, ...
    'makePlots',false,'livePlot',false,'writeVideo',false));
result6 = main_ipm(opts6);
```

The helper selects `sixth_order/weno7_fd/rk6/sixth_order`, disables adaptive
and analytic remeshing, and uses the recommended production `cfl=0.4`. Each
axis must pass a private internal grid-quality policy before a complete
operator is accepted or the expensive 2-D Poisson matrix is assembled and
factored; unsafe fixed or remeshed axes are rejected. Current evidence covers
smooth fixed uniform and admitted sinh grids only. It does not prove stability
for arbitrary nonuniform grids, SSP/TVD/positivity behavior, global order for
trigger-driven remeshing or dynamic gauges, or pure sixth-order time
convergence of the complete switching production WENO-IPM right-hand side.
See `SIXTH_ORDER_RESEARCH.md` for the exact evidence and exclusions.

## Modes and configuration

Physical evolution uses `rescalingMode='physical'`. Isotropic dynamic
rescaling uses `rescalingMode='dynamic'` with
`dynamicScaleGeometry='isotropic'`. The anisotropic correctness mode uses
`dynamicScaleGeometry='anisotropic'` and supports only `fixed` or
`peak_translation`; the latter requires half-plane symmetry.

Configuration follows one path:

```text
schema defaults -> flat overrides -> normalization -> derivation
                -> validation -> frozen groups
```

`ipm_option_schema` declares every public field. `ipm_resolve_config` returns
nine grouped domains and is called exactly once per run. Internal code never
keeps a second flat option structure. Unknown and removed options fail clearly.

The maintained structured copy resolves config schema 4 and freezes
`scalingContract='exact_gauge_no_feedback_v1'`. It has no
`maxDynamicRate`, `adaptiveLengthScaling`, `maxWidthRateCorrection`, or width
controller. `lengthScaleGain` is exactly one and the omega/width/travel
restoring gains are exactly zero. Counts of feature cells can be recorded or
used by remeshing and resolution stops, but cannot enter a scaling-rate
equation.

The four numerical selectors remain independent for baseline/fourth-order
controlled experiments, but the resolver rejects incompatible combinations.
WENO5-FD and fourth-order remapping require the paired high-order spatial
operators. Sixth-order spatial, WENO7-FD, RK6, and sixth-order remapping are
strictly bound as one complete tuple. Both full-order anisotropic Poisson paths
currently require the direct solver.

## Output contract

`main_ipm` returns result schema version 2 with grouped state, grid, scale,
physical fields, histories, quality, fit, frozen config, elliptic metadata, and
snapshots. Validate before analysis:

```matlab
result = ipm_validate_result(result);
t = result.history.common.physicalTime;
rhoPhysical = result.physical.rho;
```

New nested `result.config` uses configuration schema 4, fixes
`scalingContract='exact_gauge_no_feedback_v1'`, and omits the retired rate and
width controllers. Older result-v2 files with nested config schema 1--3 remain
readable for analysis, but cannot seed a current schema-4 checkpoint. The three
selectors absent from config schema 1 are
interpreted for validation as legacy spatial operators, SSPRK3, and PCHIP, but
the returned old structure is not modified and the fields are not injected.
Top-level result schema 1 and unversioned or aliased MAT files remain rejected.

Requested paths stay in frozen config. Collision-free paths actually used by
a run are stored in metadata, and existing MAT/video artifacts are never
overwritten. Retention rules are in `result/README.md`.

## Verification

Run the maintained second-order suite:

```matlab
verify_ipm();
```

Run the independent fourth-order certificate:

```matlab
verify_ipm_high_order_all();
verify_ipm_high_order_all('heavy');
```

Run both release gates together:

```matlab
verify_ipm_release('heavy');
```

Run the independent sixth-order certificate separately:

```matlab
verify_ipm_sixth_order_all();
```

The certified scopes, thresholds, current evidence, and explicit exclusions
are documented separately in `HIGH_ORDER_RESEARCH.md` and
`SIXTH_ORDER_RESEARCH.md`. Architectural ownership and the
transport/source/remesh contracts are in `ARCHITECTURE.md`.

## Repository map

- `main_ipm.m`: sole production solver entry.
- `ipm_*`: maintained numerical and lifecycle components.
- `verify_ipm*.m`: layered regression and accuracy certificates.
- `ARCHITECTURE.md`: executable ownership and change rules.
- `HIGH_ORDER_RESEARCH.md`: fourth-order certification scope and evidence.
- `SIXTH_ORDER_RESEARCH.md`: sixth-order certification scope, grid policy,
  evidence, and exclusions.
- `CHANGELOG.md`: maintained-program chronology.
- `experiments/` and `analysis/`: reproducible cases and post-processing.
- `result/`: ignored generated artifacts and append-only manifests.

Read `ARCHITECTURE.md` before changing formulas or adding a mode. Differences
that fit a metric, rate vector, boundary mode, or callback belong in the common
path, not a parallel solver.
