# Frozen-field analytic grid laboratory

This directory is deliberately separate from the production time integrator
and from `research/experiments/fourth_order_smooth_peak_grid.m`. Every local
entry point has the `ipm_gridlab_` prefix. Grid design and frozen-field
scoring advance no PDE time; a grid affects a continuation only after an
explicit audited checkpoint regrid transaction.

## Frozen data and feature models

`ipm_gridlab_collect_dataset` creates a schema-version-2 frozen dataset from
one or more version-2 results. Stored snapshots are used when available;
otherwise a result contributes only its terminal accepted field. Each sample
keeps the rescaled `rho/x/y`, reconstructed `physicalRho/physicalX/physicalY`,
the scale clocks and factors, physical/canonical time, trust flag, and field
signature. Each source also keeps its metadata, resolved configuration, stop
reason, node count, and grid signature. Consequently a saved dataset records
where every field and coordinate came from instead of treating an image as a
restart object.

`ipm_gridlab_feature_profiles` also emits schema version 2. It evaluates the
maintained finite-difference derivatives and records horizontal wall/source
features, vertical source/envelope features, rescaled and physical
`rho_x/rho_y` maxima, and connected-component widths at levels 0.1, 0.5, and
0.9. The 0.9 component is the strict safety core used by continuation design.

The analytic monitor families are intentionally different in the two
directions:

- `ipm_gridlab_seed_models` builds an x monitor from an asymmetric
  generalized-Student core, a narrow steep-front component, and a broad
  algebraic-tail bridge.
- `ipm_gridlab_seed_y_models` builds a boundary-centred y monitor. It disables
  the fictitious interior front and uses a broad bridge covering the union of
  the vertical source and envelope tails.

Both use smooth contrast caps. `ipm_gridlab_axis` admits odd symmetric x axes;
`ipm_gridlab_axis_positive` applies the same construction and diagnostics to a
one-sided y axis beginning at zero. A constant monitor reproduces the chosen
reference map exactly. Analytic monitors otherwise multiply an external
reference map, and inadmissible requests may be projected toward it by a
convex step in log cell width.

## Coordinated density equalization

`ipm_gridlab_equalize_axis` is the global log-cell-width fallback used when a
local analytic monitor cannot meet the requested core counts. It applies a
broad feature-centred ramp

```text
log(h_new) = log(h_reference) - A exp(-(d/sigma)^p)
```

and renormalizes the widths without moving the domain endpoints. The smallest
admissible amplitude meeting all requested cell counts is selected. This
spreads compression across many interfaces, so the adjacent density change is
shared rather than concentrated at one transition.

All axes are gated by maximum adjacent-cell ratio, maximum log-spacing
curvature, stencil conditioning, and positive quadrature margin. The adjacent
ratio cap is authoritative: neither the equalizer nor the stage designer ever
raises it automatically to make an infeasible target pass. An infeasible run
returns the strongest admitted trial and its limiting gate for inspection.

## Frozen-field scoring and stage design

`ipm_gridlab_score_frozen_pair` admits x and y together, then performs the
maintained separable high-order old-to-new-to-old transfer for every frozen
time. It reports field, wall trace, `rho_x`, `rho_y`, conservation, value-range,
x-symmetry, forward maximum-change, and strict horizontal-core,
vertical-core, and front-resolution errors. Mesh failure short-circuits the
transfer. `ipm_gridlab_evaluate_candidates` retains all rejected candidates
and ranks admitted x candidates on worst-time rather than mean-time behavior;
the pair scorer is the final two-axis check.

`ipm_gridlab_design_continuation_stage(current,reference,options)` packages the
whole offline design step. It freezes the current trusted terminal state,
extracts strict 0.9 cores, builds x/y analytic families, uses the global
equalizer when needed, applies them to axes from an external schema-v2
reference dataset, and pair-scores the selected axes. Using an external
global reference prevents successive continuation stages from ratcheting
local density ratios on top of the previous regrid.

The default design targets are at least 18 cells in each strict x/y core and 8
cells across the x front, with an adjacent-cell ratio cap of 1.08. They are
resolution/admission targets, not evidence of blowup. If either direction is
infeasible at the requested cap, the returned design records the best
attainable counts and has no usable selected pair; the caller must change the
node budget, target, model, or explicitly choose a different cap.

Minimal offline use:

```matlab
addpath('research/grid_lab')
dataset = ipm_gridlab_collect_dataset({result256a,result256b}, ...
    struct('trustedOnly',true,'outputFile','frozen_fields.mat'));
profiles = ipm_gridlab_feature_profiles(dataset);
xSeeds = ipm_gridlab_seed_models(profiles);
ySeeds = ipm_gridlab_seed_y_models(profiles);

design = ipm_gridlab_design_continuation_stage( ...
    currentResult,dataset,struct( ...
        'ratioCap',1.08, ...
        'targetXCoreCells',18, ...
        'targetYCoreCells',18, ...
        'minimumXFrontCells',8));
```

## Two distinct checkpoint contracts

There are two checkpoint-like objects in this directory tree, and they are not
interchangeable.

1. `ipm_gridlab_make_checkpoint` and `ipm_gridlab_restore_checkpoint` implement
   the original laboratory schema-v1 state package. It omits history and
   output cursors, rebuilds operators and flow for offline inspection, and is
   **not** the production dynamic-restart contract. Keep it only for isolated
   laboratory compatibility.
2. `ipm.output.makeCheckpoint/readCheckpoint/writeCheckpoint/restoreCheckpoint`
   implement the production schema-4 accepted-step checkpoint. (The compact
   solver result remains schema v2.) It contains
   the exact `rho/x/y`, base/reference axes, complete scale and rescaling state,
   initial invariants, trusted scalar history and optional snapshots, plus
   output/checkpoint cursors. Cadence and exit files are written atomically
   only from the continuous trusted prefix. `storeSnapshots=false` does not
   reduce restart completeness.

Production same-grid continuation is now supported by the single maintained
entry point:

```matlab
continued = ipm.solve(overrides,productionCheckpoint);
```

Only later terminal horizons and output controls may be overridden; numerical,
physical, grid, diagnostic, and scaling choices remain frozen. Native
same-grid split runs have been verified to reproduce uninterrupted runs
bitwise.
`ipm.output.checkpointFromResult` is a controlled bridge for an older trusted
terminal result, but reconstructed logarithmic scale variables make that
bridge machine-precision rather than bitwise exact.

## Transactional checkpoint regrid

`ipm_gridlab_regrid_checkpoint` is the only path here that turns a selected
grid into a production continuation checkpoint. It restores a trusted
production schema-4 checkpoint, advances zero PDE time, transfers `rho` with
the maintained remesh implementation, rebuilds operators and flow, preserves
scales, clocks, normalization references, initial invariants, history, and
cursors, then replaces the terminal history row at the identical accepted
time. It writes a new production checkpoint only after the transaction audit
passes.

The default transaction gates require at least 16 strict core cells in x and
y, relative `max|rho_x|` change no larger than `5e-3`, mass defect no larger
than `5e-12`, and relative range violation no larger than `2e-4`, in addition
to the mesh-quality limits (including adjacent ratio 1.08). Passing these
checks establishes a controlled numerical regrid; it does not by itself make
a subsequent singularity classification credible.

```matlab
[checkpoint,audit,fileName] = ipm_gridlab_regrid_checkpoint( ...
    sourceCheckpoint,design.selected.x,design.selected.y, ...
    struct('outputFile','regridded_checkpoint.mat'));
continued = ipm.solve(overrides,fileName);
```

The returned `audit` is part of the provenance and should be retained with the
continuation result.

## Conservative q512 campaign driver

`ipm_gridlab_run_q512_campaign` is the unattended driver for the staged q512
continuation after the manually inspected stage 3. It accepts a trusted
terminal production checkpoint (preferred) or terminal result, plus the
initial schema-v2 reference dataset. A result is converted only through
`ipm.output.checkpointFromResult`; the campaign summary and JSONL manifest
retain both direct bridge use and inherited bridge ancestry. Thus a later
native checkpoint is exact for its stored accepted state without being
mislabelled bitwise-equivalent to the trajectory before an earlier result
bridge.

By default, stage numbering starts at 4, each solve advances
`Delta tau=0.35` (with a hard maximum of 0.42), and trusted recorded
checkpoints are requested every
`Delta tau=0.20` and at exit. Direct x/y candidates are built from the
original reference axes with adjacent-ratio cap 1.08, x/y strict-core design
targets 21/20, x-front target 20, maximum log-spacing curvature 0.01,
global `min(q)/mean(q)` compatibility floor `1e-8`, and the direct y
equalizer `sigma=0.25`, `p=2`. The older `1e-4` global-ratio level remains a
separate warning threshold in telemetry; it is no longer confused with the
hard compatibility floor and is not silently ignored.
If that direct original-reference y family becomes infeasible, the explicitly
enabled final candidate is a fixed uniform global reference with
`sigma=0.35`, `p=3`; it is still subject to every unchanged gate and its
actual reference provenance is recorded in the design, result metadata,
summary, and JSONL manifest. The verified stage-4 design used the direct
original reference, not this fallback.

The campaign additionally compares every nodal quadrature weight with its
local control width
`h_CV=[dx_1/2,(dx_{i-1}+dx_i)/2,dx_N/2]`. The input, selected design,
committed regrid, and solve terminal must all have strictly positive weights
and `0.35 <= q_i/h_CV_i <= 1.65` in both axes. These local minima/maxima and
their node indices use the same definitions exposed by `ipm.mesh.quality`
and are retained in MAT summaries, JSONL records, and regrid audits. Regrid
candidate axes must also retain adjacent ratio at most 1.08 and log-spacing
curvature at most 0.01. Regrid jumps in spatial `max|rho_x|` are limited to
`2e-3` apiece and `2e-2` cumulatively; the
cumulative sum includes prior `gridLabRegrids` provenance rather than
restarting at zero.

Stage 3 to stage 4 provides the reason for separating the two quadrature
tests: the global y ratio fell from `1.526e-4` to `1.084e-4` while the local
`q_i/h_CV_i` extrema remained stable inside their hard interval. The falling
global value therefore diagnoses the growing dual-scale span between the
smallest and mean weights, not an abrupt local quadrature distortion. It is
still retained and compared with the `1e-4` warning level at design, regrid,
and terminal telemetry; only the much lower `1e-8` value is the hard
compatibility stop.

The input checkpoint and every solve terminal must retain an entirely
continuous trusted history, a trusted terminal record, at least 14 cells in
both strict 0.9 cores, and terminal safety strictly below 0.70. A trusted
native exit checkpoint is retained even if one of these stricter campaign
gates then stops automatic continuation. In particular, the known stage-2
y-core of about 13.82 does not pass this source gate; this driver is intended
to start from an accepted stage-3 terminal that has first passed it. The
wall-clock limit is a soft phase-boundary budget because design, regrid, and
`ipm.solve` currently have no mid-phase interruption callback.

Each invocation creates a unique campaign directory containing immutable
numbered MAT summaries, an append-only JSONL manifest, each full grid design,
regrid audit/checkpoint, solver result, and native cadence/exit checkpoints.
The numbered MAT summaries are authoritative; like the production output
manifest, a JSONL append failure emits a warning without invalidating an
already installed artifact.
The default target `physicalRhoXInf >= 1e3` refers to the accepted terminal
spatial maximum. Reaching it is only a numerical campaign milestone: far-box,
half-step, matched-grid, growth-fit, and other blowup credibility tests remain
separate.

The revised threshold policy and the opt-in incremental/root-corridor paths
have passed MATLAB Code Analyzer with zero findings. Each enabled fallback is
still required to pass a full-resolution design preflight and a two-step
transactional PDE smoke test before a formal continuation is launched.

```matlab
addpath('research/grid_lab')
campaign = ipm_gridlab_run_q512_campaign( ...
    stage3TerminalCheckpoint,initialReferenceDataset,struct( ...
        'wallClockBudgetSeconds',8*3600, ...
        'targetPhysicalRhoX',1e3));
```

## 2026-09-05 paused q512 continuation record

The supervised campaign was paused on request. The latest complete-history
result is stage 10 at canonical time `tau=9.17`, physical time
`t=2.00573673829`, and `max|rho_x1|=50.339074356`. Its 1836 history records
form one continuous trusted prefix. Terminal strict-core counts were
`17.7041/15.6549` and the safety factor was `0.511022`. The next stage had
advanced in memory, but interruption deliberately left no partial result;
its latest recoverable native checkpoint is the trusted stage-11 checkpoint
at `tau=9.25008277220`, `t=2.00681338692`, and
`max|rho_x1|=51.9107985671`.

Stage 10 and stage 11 used opt-in anchored incremental x equalization from
the accepted entry axis, with explicit generations one and two. Every
candidate was audited on the complete resulting axis. The stage-11 selected
design had x core/front counts `22.0000/28.9137`, maximum adjacent ratio
`1.07852052`, log-spacing curvature `0.00444213`, minimum stencil rcond
`8.67e-6`, global quadrature ratio `1.17e-4`, and local quadrature/control
ratio `[0.42266,1.48082]`. The transactional smoke test measured a per-regrid
rho-x jump `4.85e-5`, cumulative jump `0.00200786`, mass defect `5.08e-16`,
and transfer-induced range violation `9.99e-9`; all hard gates passed.

The immutable-root analytic family remained infeasible for stage 11. A new
opt-in analytic soft-corridor family was implemented and preserved for
experimentation, but the tested 16--192-cell shoulders also did not reach the
required x-core/front counts before an unchanged mesh gate or amplitude
boundary. It was therefore not used in a production regrid. This negative
result is important: merely increasing monitor strength or corridor width is
not evidence that the original-reference family can support the later dual
scale.

The trusted stage-10 growth analysis is stored beside the result as
`stage_010_growth_campaign_paused.png` and
`stage_010_growth_campaign_paused_analysis.mat`. The rho-x maximum accumulated
`4.29555` e-folds, but the canonical fit remains `window_dependent`: the two
shorter windows prefer an exponential maximum, while the two longer windows
trigger the curved-nonasymptotic guard. For the diagnostic threshold 1000,
canonical linear/exponential/curvature-guard crossings span respectively
`tau=64.33--111.28`, `16.27--16.87`, and `19.47--20.77`. These are model
extrapolations, not computed states.

A separate tail fit of `dt/dtau` gives an exponential decay rate
`0.387--0.422` and a conditional physical-time limit
`T=2.03717--2.04100`. Thus a truly linear maximum in canonical time would
diverge only as `tau -> infinity`; if the fitted time map also persisted, it
would correspond to a very weak logarithmic divergence as physical time
approached that finite `T`. The far-boundary velocity ratio remained about
one, and no half-step, larger-box, or third matched-grid continuation was
performed, so this conditional time is not a credible blow-up time.
