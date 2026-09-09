# Profile acceleration experiments

This directory tests convergence acceleration while the registered PDE
trajectory continues independently. It calls the maintained `ipm.evolve.flow`
and changes no production solver, gauge formula, time speed, or checkpoint.
An extrapolated density is a stationary-profile **candidate**, with no assigned
physical or canonical evolution time. It must never be substituted into the
accepted checkpoint chain.

## The question being tested

Let `F(R)` be the complete schema-4 isotropic rescaled spatial RHS on one fixed
grid. A smooth stationary profile in this gauge satisfies `F(R)=0`. For recent
states `R_j`, evaluate every `F_j=F(R_j)` using the current original RHS and form
`D_R=[R_1-R_0,...]` and `D_F=[F_1-F_0,...]`. With one fixed linear residual
metric `W`, the reduced secant proposal is

```text
gamma = argmin || W F_k - W D_F gamma ||_2, with SVD truncation/regularization
R_trial = R_k - damping * D_R gamma
```

The least-squares coefficients predict cancellation of resolved slow modes.
They do not establish actual nonlinear improvement. The code reevaluates the
original `F` after every proposed change and nonlinear gauge projection, and
requires a real residual decrease. This is a reduced-space secant experiment,
related to Anderson/minimal-polynomial acceleration; it does not assert exact
equivalence to any fixed-point Anderson variant. The source concepts are
[Walker and Ni (2011)](https://users.wpi.edu/~walker/Papers/Walker-Ni%2CSINUM%2CV49%2C1715-1735.pdf)
and [Scieur, d'Aspremont, and Bach (2016)](https://arxiv.org/abs/1606.04133).
The additional IPM projection and acceptance gates here require their own
numerical validation; the cited convergence theorems are not being claimed
for this nonuniform WENO problem.

`ipm_accellab_profile_probe` registers only the current exact transport anchor
at `X=1` and positive quadratic wall-gradient peak. It preserves the terminal
peak value by the algebraic retraction `R <- R * P_terminal/P(R)`, where
`P(a R)=a P(R)` for positive `a`. This fixes the amplitude coordinate of an
independent candidate. It is not an added feedback source in the evolution
equations. Active peak-stencil switches and WENO stencil changes limit local
smoothness and can invalidate the secant model; fresh residual rejection is
therefore essential.

The metric concatenates full-domain density RHS, core gradient RHS, and wall
gradient RHS with fixed baseline normalizations and positive quadrature
weights. Core windows default to `[0,4] x [0,2]`. Independent gates also require
the full-domain and core infinity norms, and each component L2 norm, to remain
nonincreasing. A large nondecaying outer density cannot conceal deterioration
of the inner gradients by dominating a single global norm. The coefficient
norm and relative jump have hard bounds. All projection/evaluation exceptions
are recorded as rejected trials; invalid input histories are errors.

## A smooth stationary limit is a hypothesis

The current quadratic gauge fixes the amplitude of `R_X`, while its narrow
core can keep contracting. `F(R)` or its differentiated residual can fail to
approach zero in this fixed isotropic coordinate even when a meaningful
singular or multiscale limit exists. A rejected smooth-profile candidate does
not reject such a limit, and an accepted candidate does not establish it.

`ipm_accellab_shape_history` therefore separately records the 90%-level core
width and 50%-level width on the wall and vertical line, the quadratic peak
curvature radius, and log core-width slopes against canonical time. It samples
`R_X/P` in inner coordinates using the respective instantaneous core widths.
The vertical width is measured at the active wall peak node; the vertical
profile itself is interpolated to the quadratic vertex. These linear
interpolations are read-only diagnostics, without a spatial-order claim.
The report retains both macro/core width ratios and inner-profile changes.
It does not convert a few width samples into an asymptotic scaling law.

## Running bounded tests

From the `ipm_structured` directory:

```matlab
addpath('research/acceleration_lab');
report = ipm_accellab_test();
```

The synthetic problem has contraction factors `0.9999` and `0.91`; six samples
are spaced by 100 ordinary iterations. It checks removal of resolved slow
modes, rank-deficient rejection, freshly evaluated nonlinear rejection,
post-projection reevaluation, and an independent veto. The small IPM check
uses `65 x 33` nodes and three ordinary accepted SSPRK(5,4) steps, then tests
the real RHS and fixed-gauge adapter. It does not start any production run.

On 2026-09-08, the synthetic residual ratio was `1.8674e-15`. The unaccelerated
slow component would need 138149 iterations for a factor of one million.
With samples only one step apart, the coefficient norm was `3.3902e4`, so the
default `1e4` safeguard correctly refused extrapolation. Sparse sampling
resolved the same slow modes without relaxing that bound.

The default tiny IPM probe was rejected by the coefficient guard. An explicitly
registered extra test with memory 1, coefficient ceiling `1e8`, and damping
`[0.01,0.001,0.0001]` reached genuine projected RHS evaluations; all trials were
rejected by the independent residual gates. Its baseline full/core/wall L2
residuals were about `0.1207/0.4531/0.2062`, with gauge residual
`-1.56e-16` and Poisson residual `1.14e-13`. This validates rejection behavior,
not an IPM speedup. The four initial `.m` files returned zero Code Analyzer
diagnostics after the final tiny run.

## Real checkpoint probes

```matlab
% files: chronological cell array of 2--6 validated same-root checkpoints
experiment = ipm_accellab_probe_checkpoints(files,struct( ...
    'allowHistoryProjection',true, ...
    'memories',[1,2,4], ...
    'outputRoot','result/longtime/20260908_campaign_v1/acceleration_lab'));
```

The wrapper validates every original checkpoint, reconstructs operators only
once, checks identical grid/runtime references/run metadata, and shares those
operators across the stored fields. Use the original 10-thread configuration
for the current q512 checkpoints: their existing signature can depend on the
floating-point reduction layout and is not automatically portable to a
single-thread read. No signature check is bypassed.

The raw quadratic-peak drift is always saved. By default drift beyond `1e-5`
rejects the history. Only explicit `allowHistoryProjection=true` permits
positive amplitude retraction of every history field to the terminal peak;
the original rejection, original drift, and all scale factors are retained.
Every projected field then gets a fresh original-RHS evaluation. It is labelled
an independent research history, not a projected physical time trajectory.

Registration, trial density candidates, full audits, options, errors, and
summary MAT/JSON files go into a unique output directory. The function never
replaces old evidence. A profile claim would still require held-out history
windows, a separately registered short relaxation using unmodified dynamics,
and matched grid/box checks. Decreased residual alone is not a singularity
claim or proof of faster production time integration.

## Actual q512 probe, 2026-09-08

The five root checkpoints at steps 1587/1660/1729/1794/1818 cover canonical
time `3.0023359 -> 3.6732019`. The original peak drift was `8.8471e-5`, so the
raw `1e-5` amplitude-history gate rejected them. Explicit positive projection
used factors from `1.00008848` down to `1`; this is saved as an independent
research history. The only excluded lineage field is the changing nonnumerical
`runMetadata.latestCheckpointFile`; all axes, complete saved rescaling, and
remaining run metadata must compare exactly. The earlier rejection caused by
checking that path too is separately preserved.

All three memory choices 1/2/4 were rejected after actual original-RHS
evaluations. For memory 4 and damping 0.5, full/core/wall L2 residuals fell to
`0.868/0.775/0.848` of baseline, but the core gradient RHS infinity norm grew
to `1.812` times baseline. Even the smaller tested steps made that local
maximum worse. Accepting only an average residual would have falsely labelled
this as a successful acceleration. The baseline full/core/wall relative L2
residuals are `0.26824/0.32832/0.31895`; they are not near zero.

The wall/vertical 90% core widths fell from `0.07306/0.02201` to
`0.04716/0.01393`, with fitted logarithmic slopes `-0.6521/-0.6513` per
canonical time. The wall core/FWHM ratio changed only `0.21164 -> 0.21314`.
After instantaneous core-width rescaling, the final two wall and vertical
profiles differ by `0.272%` and `1.272%` over `Delta tau=0.07305`. This short
window motivates an inner-coordinate motion diagnostic. It is not evidence
that the final shape or contraction rates have converged.

Complete [MAT/JSON evidence](../../result/longtime/20260908_campaign_v1/acceleration_lab/probe_tpc0982a58_aee1_4234_aad4_7ab261010094/experiment.json)
includes every trial's individual residual components and rejection. The
executed wrapper took `38.794 s`, built one operator set, and advanced zero PDE
steps. All five MATLAB files present at that point had zero Code Analyzer
diagnostics.

## Modulated inner residual: a separate, falsifiable diagnosis

Let `Omega=R_X`, `a` be its active three-node quadratic wall peak, and let the
positive measured inner lengths be `ell_x,ell_y`. In coordinates
`xi=(X-a)/ell_x`, `eta=Y/ell_y`, the canonical derivative of the unnormalized
gradient is

```text
Omega_tau|xi,eta = F_X + aPrime R_XX
                      + beta (X-a) R_XX + gamma Y R_XY,
beta = (log ell_x)' , gamma = (log ell_y)' .
```

`ipm_accellab_modulated_probe` evaluates the original `F` once, then fits those
three motion columns by weighted least squares. The training rectangle is
preregistered as `|xi|<=1, 0<=eta<=1.5`; a surrounding held-out rectangle
`|xi|<=2, 0<=eta<=3` excludes the training nodes. The fit has no knowledge of
the residual on that outer ring. Full-domain/core/wall/training/held-out L2
and infinity residuals are reported separately. Columns are normalized before
SVD, and rank/condition failures are explicit. A tiny-grid smoke can expand
both rectangles explicitly to ensure enough grid nodes; those options are
saved and are not a q512 registration.

An independent strict peak speed is derived using the same three active
nonuniform nodes and the same polynomial vertex throughout:

```text
Q[Omega]'(a)=0
aPrime = -Q[F_X]'(a) / Q[Omega]''(a).
```

This follows by differentiating the polynomial stationarity condition on a
fixed active stencil. A second fit fixes that strict `aPrime` and solves only
for `beta,gamma`, allowing comparison with the fully free fit. This is an
exact statement about the three-node polynomial; it is not identical to
sampling matrix-differentiated `R_XX` at `a`, and it does not remove the
nonsmoothness of active-stencil switches. The diagnostic records curvature,
active-node margin, vertex margin, and the phase residual.

Because physical lengths are `ell_x/C_l` and `ell_y/C_l` in isotropic scaling,
their logarithmic rates per canonical time are `beta-c_l` and `gamma-c_l`;
positive physical compression rates are therefore `c_l-beta` and
`c_l-gamma`. Conversion to derivatives per physical time additionally requires
the saved clock Jacobian. The routine reports only canonical-time rates.

These parameters currently describe a read-only coordinate motion. They do
not replace `c_l`, `c_omega`, `c_r`, modify the RHS, or supply synthetic
accepted times. A good training fit alone may only remove a local geometric
component and says nothing about held-out shape drift or singularity.
Connecting this to a new exact evolution gauge would require separately
registered peak and width functionals, differentiating all of them with the
same discrete operators, solving their coupled instantaneous conditions,
respecting symmetry/translation compatibility, and validating physical
covariance against the existing trajectory. No such C-rule modification is
implemented here.

`ipm_accellab_test_modulation` ran on 2026-09-08. An exact polynomial example
with `(aPrime,beta,gamma)=(0.031,-0.65,-0.67)` recovered the rates with error
`4.496e-16` and the strict nonuniform quadratic peak speed with error
`7.841e-16`. A deliberately contaminated outer ring remained visible despite
training residual cancellation, and duplicate geometric columns were
rejected by the rank gate. A `65 x 33` original-RHS probe also passed without
taking a PDE step. Its free-fit scaled condition was `3.166`; these smoke
rates are not the late q512 contraction rates. The saved evidence is
`result/longtime/20260908_campaign_v1/acceleration_lab/modulated_tiny_20260908.mat`.
The first actual q512 modulated probe is recorded below.

The optional `referenceLogWidthRates=[wallSlope,verticalSlope]` together with
its `referenceCanonicalWindow=[tauStart,tauEnd]` records an explicit comparison
against already measured width slopes, such as the `-0.6521/-0.6513` window
above. These values never enter the least-squares system or the PDE RHS.

### Actual single-state q512 inner-motion probe

The baseline segment-002 checkpoint at step 1936 was independently validated
and restored once, then `ipm_accellab_modulated_probe` evaluated one original
RHS. The process exited after `23.556 s` without advancing a PDE step. The
state has `tau=4.0732019386`, physical time `1.7550567738`,
`c_l=0.2103293070`, and `c_omega=-0.2861489653`.

| Diagnostic | Free translation and scales | Strict quadratic translation |
|---|---:|---:|
| `aPrime` | -0.00570124 | -0.00499849 |
| `beta` | -0.619147 | -0.602286 |
| `gamma` | -0.708103 | -0.665416 |
| Scaled matrix condition | 3.906 | 1.900 |
| Training L2 residual ratio | 0.03504 | 0.04387 |
| Outer holdout L2 residual ratio | 0.09964 | 0.10190 |
| Outer holdout infinity residual ratio | 0.10235 | 0.09493 |
| Full-domain infinity residual ratio | 91.353 | 88.927 |
| Physical compression x / canonical time | 0.82948 | 0.81262 |
| Physical compression y / canonical time | 0.91843 | 0.87575 |

The fit used 792 training nodes; the disjoint outer ring contained 1743
nodes. A similar reduction on that held-out region shows that the local
motion fit captures substantially more than a training-point cancellation.
The strict-translation holdout residual is still about `0.01655` in the
registered amplitude-normalized L2 units, so it is not a zero-residual
stationary profile. The reference width slopes `-0.6521/-0.6513` came from the
earlier window `tau=[3.0023359286,3.6732019386]`, not an independent derivative
estimate at this later time.

The severe full-domain deterioration explicitly rules out using this local
fit as a global C-rule change. The result supports a locally shrinking frame
as the next profile-analysis coordinate. It does not show global stationarity,
an accepted accelerated trajectory, asymptotic rate convergence, or strong
singularity. Other time windows, grid refinement, and preregistered window
variations remain necessary.

Complete [single-state evidence](../../result/longtime/20260908_campaign_v1/acceleration_lab/modulated_q512_tpdb39fb5f_1136_46a8_b9b1_ba459207a48e/report.json)
contains the free and fixed-phase residuals on full/core/wall/training/holdout
regions, the polynomial peak-speed identity, conditions, node counts, and
the unchanged gauge/Poisson residuals.

### Exact instantaneous inner coordinates

Enable `includeExactCoordinates=true` in `ipm_accellab_modulated_probe` to
compare fitted rates with a fully specified discrete-coordinate derivative.
This reuses the same already evaluated `F_X`; it adds no PDE step or RHS
evaluation. `ipm_accellab_exact_inner_rates` can also be called directly on
paired `Omega=R_X`, `forcing=F_X`, axes, and the active peak index.

The coordinate definitions are:

1. `P=Q[Omega](a)` and `PPrime=Q[F_X](a)` on the same active three-node
   quadratic wall template. `PPrime` is retained exactly as evaluated; it is
   not replaced by zero even if the production gauge should conserve `P`.
2. Horizontal coordinates `xLeft,xRight` are the connected, strict interior
   piecewise-linear crossings of `0.9*P` surrounding the active peak. Each
   crossing obeys

   ```text
   crossingPrime = (0.9*PPrime - interp(F_X)) / segmentSlope
   wx = xRight-xLeft
   beta = (xRightPrime-xLeftPrime)/wx.
   ```

3. The vertical trace at moving `X=a` uses the same three columns and
   quadratic Lagrange weights at every y-row, not a new interpolation rule.
   Its material forcing is

   ```text
   verticalOmega = Q[Omega](a,Y)
   verticalForcing = Q[F_X](a,Y) + aPrime*d_X Q[Omega](a,Y).
   ```

   The first downward piecewise-linear crossing of `0.9*P` above `Y=0`
   defines `wy`. Its implicit derivative gives `gamma=wyPrime/wy`.

The peak stencil, three crossing brackets and their coordinates, segment
slopes and relative drop, interpolation fractions, margins to node/template
boundaries, and raw implicit residuals are all recorded. A supplied
`exactCoordinateOptions.previousSignature` yields an explicit unchanged or
switched status; without it the status is `not_compared`. A near-tied active
peak or near-node crossing is rejected as a switching boundary, and an
unbracketed, flat, reversed, or nonfinite crossing is rejected. In the
combined probe, these expected coordinate failures remain visible as an
invalid exact-coordinate report alongside the completed LS diagnostics.

The exact rates are substituted into the same residual and exactly the same
observation masks as the previously registered free/fixed-translation fits.
The old mask widths use the maintained feature diagnostic's node-peak
threshold; they are not silently replaced by the new quadratic-P widths.
Both widths are saved, making this comparison explicit. The rate formula
itself uses the new exact coordinates. These remain local coordinate
diagnostics and never modify `C`, a clock, a density, or an accepted state.

`ipm_accellab_test_exact_inner` passed on 2026-09-08. On nonuniform grids, a
moving/scaling tilted polynomial and a tilted Gaussian gave final central
time-difference errors `3.944e-11` and `2.888e-11` for `[P,a,wx,wy]`
derivatives. The polynomial has known
`aPrime=0.031`, `beta=-0.65`, `gamma=-0.67`, and amplitude rate `0.07`.
Its quadratic peak motion and vertical derivative are exact up to roundoff;
the horizontal *piecewise-linear* width has its own discrete derivative.
Its beta errors against the continuous ideal on 81/161/321/641 x-nodes were
`0.0418/0.0026/0.0040/0.0045`: the error is not monotone with grid phase, and
no spatial order is claimed from that sequence. The time-difference check
validates the exact discrete derivative rather than asserting equality with
the continuous ideal at finite spacing.

A pure-amplitude test with nonzero `PPrime` recovers zero translation and
width rates, which detects an erroneous forced `PPrime=0`. Tests also detect
real later-time template changes and reject crossings exactly at nodes. A
tiny original-RHS comparison on the identical LS/holdout masks and the prior
modulation tests both passed. Full evidence is in
`result/longtime/20260908_campaign_v1/acceleration_lab/exact_inner_tiny_tpd9a53ab7_75c5_429c_9731_2084147db6e4/`.
The exact-coordinate q512 comparison is recorded next.

### Actual exact-coordinate q512 comparison at step 2031

The native adaptive stage-001 checkpoint at
`tau=4.2832019445`, physical time `1.7798569844`, was validated with the
original thread configuration, restored once, and evaluated once more by
the original RHS. It advanced no PDE step and exited after `33.202 s`.

| Rate / residual | Free LS | Fixed peak-speed LS | Exact crossing rates |
|---|---:|---:|---:|
| `aPrime` | -0.00411786 | -0.00353118 | -0.00353118 |
| `beta` | -0.620710 | -0.604757 | -0.667908 |
| `gamma` | -0.710509 | -0.668517 | -0.676751 |
| Outer holdout L2 ratio | 0.09205 | 0.09431 | 0.15920 |
| Outer holdout infinity ratio | 0.09437 | 0.08770 | 0.17927 |

Without any fitted coordinate rate, exact level coordinates remove about
84% of the held-out L2 evolution and 82% of its maximum. The fixed-phase LS
horizontal rate differs from the exact width rate by `0.06315`, about 9.5%
of the exact rate; it cannot be described as the same exact width gauge. The
vertical difference is `0.008234`. The lower LS residual is a best local
geometric fit, not evidence that its coefficients are these exact coordinate
derivatives.

The genuine `PPrime` was `-1.3248e-15`. Horizontal left/right and vertical
crossing template margins were `0.07489/0.31457/0.24528`; the active peak
selection margin was `0.001293`. This single sample has no preceding template
signature and correctly records `templateStatus=not_compared`. Physical
inner compression rates were `0.86455/0.87339` per canonical time.

The original Poisson solve warned twice (restore and probe) that the matrix
was ill conditioned, with reported `RCOND=1.479364e-16`; its backward
residual was `4.0279e-10`. These warnings were retained, not suppressed.
The exact-coordinate validity flag refers to its local functionals and
conditions, not a forward-error guarantee for this elliptic solve. Spatial
and truncation checks remain necessary. As before, the local coordinate
motion worsens full-domain maximum residual (here by a factor about 106),
so it is not a global C-rule proposal.

Complete [exact-rate evidence](../../result/longtime/20260908_campaign_v1/acceleration_lab/exact_inner_q512_tpc9fd0098_f80e_450c_aacf_50d097ce075b/report.json)
contains both LS alternatives on identical masks, exact implicit identities,
peak and crossing conditions, the actual `PPrime`, and all held-out metrics.

### An exact discrete observer of the inner shape

`ipm_accellab_pullback` defines a local observer of a complete full-domain
state, using the exact coordinates above and one fixed bilinear interpolant
`I_h` on the native axes:

```text
Omega = R_X,                  F = original full-domain RHS(R)
U(xi,eta) = I_h[Omega](a + wx*xi, wy*eta)/P

G(xi,eta) = [ I_h[F_X]
              + (aPrime + wxPrime*xi)*d_X I_h[Omega]
              + wyPrime*eta*d_Y I_h[Omega] ]/P
            - (PPrime/P)*U.
```

`G` is the instantaneous derivative of this particular discrete observer
along the original semidiscrete PDE, on fixed peak/crossing/observer cell
templates. The spatial derivatives here are derivatives of **the same
bilinear interpolant**. Replacing them by interpolated finite-difference
derivatives would differentiate a different observer. The previously reported
native finite-difference modulated residual remains a separate check. `PPrime`
is the actually evaluated quadratic functional, including any nonzero value;
neither the sampled template nor the forcing is corrected to make it vanish.
All observer cell indices, interpolation fractions, and moving-cell margins
are returned. Cross-cell finite differences are excluded from derivative
tests and the excluded fraction is reported.

The complete physical state is still required: `U` alone does not determine
the nonlocal Poisson velocity or the matching tail. Therefore
`Q(Flow_h(R)) = Q(R) + h*G(R) + O(h^2)` is a verified observer of the original
flow, **not yet a closed autonomous map on U**. A theorem about fixed-point
acceleration on a shape quotient would additionally require a consistent
full-field extension and control of its outer data. The current experiment
retains that full field explicitly and makes no such theorem claim.

### Guarded full-field secant for the inner shape

`ipm_accellab_quotient_secant` uses complete same-grid native `rho` histories.
It forms a reduced secant direction from differences of these full fields
and of their freshly computed `G` vectors. Every history and candidate is
explicitly multiplied by a positive scalar to have the terminal quadratic
peak amplitude before evaluating its original RHS. This is an independent
profile retraction, not an additional term in a time-stepping equation. Each
candidate then receives a fresh full-domain Poisson solve and original RHS.
No arbitrary two-dimensional inner patch is pasted into a global state.

The observer axes are fixed at 65 x-nodes on `[-2,2]` and 49 y-nodes on `[0,3]`.
The training region is `|xi| <= 1, eta <= 1.5`; its complement in that
rectangle is an outer holdout. These masks and quadrature weights are frozen
before testing candidates. All the following conditions are required:

- Fresh weighted training L2 decreases by at least
  `minimumDecrease*damping`, with `minimumDecrease=1e-3`.
- Training infinity, outer holdout L2/infinity, and wall infinity do not
  increase, allowing only `1e-12` relative plus `1e-14` absolute roundoff.
- Both core and outer holdout L2/infinity of the native finite-difference
  exact-motion residual do not increase. This independently checks that the
  interpolated observer has not hidden a local spike.
- Outside the fixed native box `|X|<2.5, Y<2`, the weighted relative L2
  density difference and velocity difference are each at most `1e-3`;
  the maximum source difference `|OmegaCandidate-OmegaBase|/P` is at most
  `1e-3`. Full-field construction plus these explicit tail checks is the
  matching protocol, not an exact outer Dirichlet condition.
- Quadratic amplitude mismatch is at most `1e-8`, the global relative field
  jump is at most `0.05`, the coefficient norm at most `1e4`, and the secant
  SVD relative floor is `1e-8`. Exact coordinate admissibility is mandatory.

Every attempt retains its numerical audit, any exception, the resolved
options, named failed gates, and required training norm. The return value
`accepted_stationary_candidate` only means that an **independent full-field
profile candidate** passed these tests. It is never saved as a native
checkpoint or inserted into a physical-time trajectory.

### Small-grid evidence and its limits

`ipm_accellab_test_quotient` passed on 2026-09-08. The synthetic observer
derivative had central time-difference error `2.217e-10`, and a pure amplitude
mode with nonzero `PPrime` was removed to roundoff. A 129 x 65 original-PDE
history supplied five full states at canonical times `0:0.016:0.064`.
For original SSPRK54 short flows of length `0.001, 0.0005, 0.00025`, the
maximum observer difference-quotient errors were
`2.764e-3, 1.380e-3, 6.899e-4`. Respectively 98.46%, 98.46%, and 100% of
observer nodes retained their cells. This validates the first-order
difference-quotient limit, not a spatial accuracy order or a stationary
shape assumption.

Memories 1 and 2 rejected every tested damping: large steps exceeded tail
tolerances, while smaller steps increased the training infinity residual.
Memory 4 with damping `1/16` passed. Its relative field jump was `0.0012562`;
tail density/source/velocity differences were
`9.593e-4 / 9.404e-4 / 5.065e-4`.

| Residual | Candidate / baseline immediately | After 4 original steps of 0.002 |
|---|---:|---:|
| Training L2 | 0.96941 | 0.96896 |
| Training infinity | 0.99972 | 0.96794 |
| Outer holdout L2 | 0.93384 | 0.91319 |
| Outer holdout infinity | 0.92430 | 0.84970 |
| Wall infinity | 0.89676 | 0.89373 |
| Native core L2 | 0.95357 | 0.95264 |
| Native core infinity | 0.99690 | 0.99726 |
| Native outer L2 | 0.93323 | 0.93341 |
| Native outer infinity | 0.92997 | 0.93121 |

The short relaxation uses two independent research initial-value problems
and the unmodified original PDE. The residual improvement persists over this
short check. However, both branches' absolute training L2 increased over
that interval (baseline `0.34627 -> 0.35594`, candidate
`0.33568 -> 0.34489`). This is evidence for a small guarded shape-residual
improvement, not established convergence acceleration, a stable profile,
or strong singularity. Early small-box data cannot establish late q512
behavior. The initial tail gates were not reapplied as a claim of exact tail
matching throughout the independent relaxation.

The initial complete report is in
`result/longtime/20260908_campaign_v1/acceleration_lab/quotient_tiny_tp216bf4bb_9da2_47ff_b4ed_e47a5b80d3bc/`.
The follow-up
`result/longtime/20260908_campaign_v1/acceleration_lab/quotient_relaxation_tp877bd519_ee6a_483f_8069_03fd853fd670/report.mat`
also saves complete source, history, and candidate density arrays as research
artifacts, plus configuration, axes, and runtime references. It has no native
checkpoint schema and cannot be resumed as a trusted production record.

### Fixed real-checkpoint experiment entry point

The next experiment fixes memory 4 and the **single damping 1/16** before
reading the q512 history. `ipm_accellab_quotient_checkpoints(files,user)` takes
four or five chronological native checkpoint paths. Four frames give three
effective differences, explicitly recorded. It does not automatically mix
grids, bypass signatures, skip rejected files, or scan parameters.

Each frame passes the maintained `ipm.output.readCheckpoint` signature check
under its original reduction/thread configuration. All native/current/base
axes, complete rescaling runtime references, remesh count, initial mass/range,
and numerical configuration must be exactly equal. Only terminal horizons
and output controls may differ in config. An explicit allowlist covers
non-numerical resume bookkeeping, output paths, and case metadata; every
difference and its actual values is saved. The remaining case/solver/gauge
lineage must be identical. Unexpected differences reject before restoration.

Only the final frame is restored. Its operator structure is shared with
the other saved full fields and scales. The maintained restore routine may
itself build the configured mesh and then rebuild for native remeshed axes;
the wrapper reports **one restore call**, not an uninstrumented claim of one
LU factorization. For five frames the fixed candidate normally needs six
fresh RHS evaluations after restoration: baseline, four earlier profiles,
and one candidate. A guard may stop earlier.

The wrapper creates a unique directory containing preregistration, successful
or failed experiment report, actual raw P drift, explicit history amplitude
factors, and `candidate.mat` as an independent research full field. Native
checkpoints are never written. Use only in a separately agreed q512 memory
window:

```matlab
experiment = ipm_accellab_quotient_checkpoints(files,struct( ...
    'outputRoot','result/longtime/20260908_campaign_v1/acceleration_lab', ...
    'expectedRemeshCount',3));
```

`ipm_accellab_test_quotient_checkpoints` passed on native 129 x 65 checkpoints.
It reproduced the fixed candidate with one restore and six fresh RHS calls.
A separately signed frame with changed WENO epsilon was rejected before
restore, and the failed report was preserved. Evidence is in
`result/longtime/20260908_campaign_v1/acceleration_lab/quotient_native_tiny_tpd9b11acd_3088_480d_91a2_9ca99e189537/`.

### Actual fixed q512 secant: rejected

The five preregistered native frames at steps
`2135/2182/2226/2248/2270`, canonical times
`4.4832019445/4.5848561299/4.6832019445/4.7336235901/4.7832019445`,
all passed their original ten-thread native signatures and the complete
same-grid/numerical-configuration checks with `remeshCount=3`. All actual
output/horizon and resume metadata differences were recorded. The maximum
raw relative P drift was only `2.8644e-14`; this history does not have the
earlier amplitude-drift problem.

One restore and six additional original RHS evaluations took `39.88 s`.
The fixed memory-4, damping-1/16 candidate was **rejected**, with no further
parameter scan. Its coefficient norm was `0.74969`, rank 4, and global
relative field jump `2.68096e-4`.

| Registered metric | Candidate / baseline | Outcome |
|---|---:|---|
| Bilinear pullback training L2 | 1.33452 | fails |
| Bilinear pullback training infinity | 1.85629 | fails |
| Bilinear pullback outer L2 | 1.00385 | fails |
| Bilinear pullback outer infinity | 0.99215 | passes |
| Native FD core L2 | 0.99970 | passes, tiny change |
| Native FD core infinity | 0.99943 | passes, tiny change |
| Native FD outer L2 | 0.99971 | passes, tiny change |
| Native FD outer infinity | 0.99943 | passes, tiny change |

Tail density/source/velocity changes were
`2.57925e-4 / 5.27083e-5 / 1.04319e-3`: the velocity change exceeds its
`1e-3` gate by 4.32%. Baseline and candidate actual PPrime were
`-1.21084e-15` and `9.27860e-16`. The unchanged original Poisson solve warned
at every evaluation with `RCOND=8.911277e-17`; baseline and candidate
backward residuals were `4.79415e-10` and `5.41604e-10`. Warnings were retained
in the full launch log and each audit. No accelerated q512 trajectory was
accepted, and small-grid success did not transfer to this fixed experiment.

The discrepancy between nearly unchanged native FD residuals and the
bilinear observer's increase motivates measuring observation-discretization
effects. It does **not** establish interpolation error as the cause and
does not permit replacing the acceptance metric after seeing the result.

Complete output is in
`result/longtime/20260908_campaign_v1/acceleration_lab/quotient_checkpoints_tp02bd67ac_8d70_431c_ad1a_9eb3d603ed96/`.
The generic API returns the baseline on rejection, so `candidate.mat` stores
that fallback as `rhoCandidate`. The actual rejected field was subsequently
materialized from the recorded full-field history, coefficients, and fixed
damping, with no LU, RHS, or time step. Its field-jump reproduction was
checked to `1e-13` and it is saved separately in
`materialized_rejection_tp71dbbaee_5fc3_4906_ab0f_629eb7677b54/rejected_trial.mat`.
It retains the original rejected decision and is not a native checkpoint.

### Conditional amplitude-width index, without a guessed singular time

At step 2031 let `C_l=exp(logC_l)`, `C_w=exp(logC_omega)`,
`ell_x=wx/C_l`, and `G_star=P*C_l/C_w`. Integrating a normalized inner wall
shape over a fixed inner interval gives

```text
Delta rho = G_star * ell_x * A(tau),    A = integral phi(xi,tau) dxi.
```

If the shape integral A is approximately constant, the measured rates give

```text
(log G_star)'    = c_l-c_w+PPrime/P =  0.4912892178
(log ell_x)'     = beta-c_l        = -0.8645496753
(log Delta rho)' = beta-c_w+PPrime/P = -0.3732604575

alpha_x = (beta-c_w+PPrime/P)/(beta-c_l) = 0.4317397463.
```

Equivalently the instantaneous `log G_star` versus `log ell_x` slope is
`alpha_x-1=-0.5682602537`. The actual `PPrime/P=-1.9305e-15` is used.
If A drifts, the density-increment index has the extra correction
`(APrime/A)/(beta-c_l)`, which this diagnostic has not measured.

The existing paired profile series gives adjacent-window proxy indices
`0.425014` over `tau=[4.07320,4.28320]` and `0.423930` over
`[4.28320,4.48320]`. A three-point log-log fit spanning both gives `0.424494`;
the preceding three-point fit gives `0.428964`. All adjacent windows and
the full eight-point fit (`0.433294`) are saved, rather than selecting only
the most favorable agreement. The series uses global grid `max|rho_x|` and
the maintained node-threshold width, whereas the instantaneous formula uses
the tracked quadratic peak and its exact width. At step 2031 the two wall
width definitions differ relatively by `5.986e-5`; the series also spans
remesh transactions. The roughly `0.00725` difference from the centered fit
is retained without claiming a statistical error bar.

This is a temporal amplitude-width scale index under a shape assumption.
It is neither a measured finite-time spatial Holder exponent nor an
established limiting Holder exponent. A limiting claim additionally needs
controlled profile/center convergence, stable scale/window behavior and
discretization/box checks. No singular time T is fitted. Even under the
hypothetical laws `G~(T-t)^(-p)` and `ell~(T-t)^q`, only
`alpha=1-p/q` follows; neither p, q, nor T is determined here.

`python3 research/acceleration_lab/ipm_accellab_scale_index.py` reruns this
standard-library-only diagnostic from saved JSON files; it opens no native
checkpoint and creates no RHS evaluation. The source hashes, exact formula,
all windows and limitations are preserved in
`result/longtime/20260908_campaign_v1/acceleration_lab/scale_index_57eb56a217734e5d828465ef428181f3/`.

### Linear C1 Hermite observer sensitivity experiment

`ipm_accellab_tensor_hermite` defines a tensor cubic Hermite interpolant.
The nodal jets are exactly `field`, `field*Dx'`, `Dy*field`, and
`Dy*(field*Dx')` for the supplied maintained native derivative operators.
Adjacent cells share these same jets. Their value and both first spatial
derivatives agree on interfaces, giving C1 continuity on a fixed grid.
The mixed jet uses the same prescribed operator composition everywhere.

The field-to-interpolant map is linear. Therefore its fixed-position time
derivative is the same Hermite operator applied to `F_X`, including its
jets from the same Dx/Dy. `ipm_accellab_pullback_hermite` adds moving-query
spatial derivatives taken from this **same Hermite polynomial** and retains
the genuine `-PPrime/P*U` term. No nonlinear PCHIP slope rule is used.
Peak and width coordinate templates retain their existing switching rules;
C1 observation interpolation does not make those separate functionals smooth
at a switching event or make different remeshed grids equivalent.

`ipm_accellab_compare_saved_candidate` is prepared for a separately granted
large-grid window. It keeps the existing baseline and rejected candidate
fixed, restores the baseline once and evaluates two original RHS values.
It records bilinear and Hermite metrics on the identical registered observer
nodes/masks, all prior native FD metrics, fresh tail differences, and actual
PPrime. It stores Omega and F_X as independent fields so further observer
diagnostics can run without another LU or RHS. It explicitly preserves
`originalDecision=rejected_unchanged`; it cannot accept the old candidate
under a newly chosen metric or rescan its damping.

`ipm_accellab_test_hermite` passed without a Poisson operator, RHS evaluation,
or PDE step. On nonuniform grids, a bicubic polynomial's value/Dx/Dy errors
were at most `2.761e-14` and the field-linearity error was `2.387e-14`.
The moving/scaling Gaussian observer, with its actual nonzero PPrime, had
central time-difference errors decreasing from `5.032e-8` at `h=1e-3` to
`4.357e-11` at `h=3e-5`, while its exact peak/width templates stayed fixed.
A pure amplitude mode left only `1.388e-16` residual.

One-sided values and both first derivatives approached matching interface
limits as the distance decreased: maximum discrepancies were
`7.684e-5/7.684e-7/7.684e-9` at distances
`1.859e-5/1.859e-7/1.859e-9`. A separate manufactured moving query actually
crossed an interface in every central difference; its derivative error
decreased `1.445e-8 -> 1.395e-11`. C1 continuity alone does not guarantee a
globally smooth second derivative, so a universal second-order central-time
claim at interfaces is not inferred from that one example.

For a smooth nonpolynomial function on 21/41/81/161 x-nodes and
17/33/65/129 y-nodes, Hermite value errors were
`3.773e-5/2.329e-6/1.198e-7/8.138e-9`, giving about fourth-order values;
Dx errors showed about third order. The corresponding bilinear value errors
were `8.593e-3/2.190e-3/4.968e-4/1.300e-4`. This quantifies interpolation
accuracy on the prescribed smooth example; it does not establish the source
of the q512 residual discrepancy.

The initial C1 smoke caught a row-query/column-axis implicit-expansion bug.
It was corrected by explicitly preserving query-array shapes and the entire
suite was rerun. That failed test is preserved in
`result/longtime/20260908_campaign_v1/acceleration_lab/hermite_failure_e2d30c9ecbd74397acdee5ad06f3464a/`;
the successful full report is in
`result/longtime/20260908_campaign_v1/acceleration_lab/hermite_tiny_tpf9d57151_9002_4672_886f_1420e9ed2f4f/`.

### Fixed-comparison source audit and no-LU follow-up

The source audit was strengthened before the actual comparison. A matching
scalar relative jump does not identify an entire rejected field.
`ipm_accellab_verify_rejected_source` now rereads all five original native
payloads with their strict signature checks, checks full native/base axes,
numerical configuration, runtime gauge references, case lineage and recorded
step/time identity, then reconstructs the recorded full-array secant formula.
The saved rejected array must agree elementwise within 64 ulps at its global
scale; exact equality and maximum discrepancy are reported. The recorded
residual audit and failed gates must also match. This happens before restore,
with only one-dimensional derivative operators. After the one restore, those
operators must equal the maintained `ops.Dx/ops.Dy` exactly.

The subsequent field-only follow-up is
`ipm_accellab_analyze_saved_observers(observerFieldsFile,comparisonReportFile,outputRoot)`.
It opens the saved `Omega`, `F_X`, exact coordinates, and derivative matrices;
it opens no native checkpoint, creates no Poisson operator and evaluates no
new RHS. It always retains the original comparison and rejection.

The fixed protocol reports all of the following, without selecting a new
best metric or parameter:

1. Observer node densities 65 x 49, 129 x 97, and 257 x 193, each at zero
   and half-step interior phase. Region and wall endpoints are always
   included. Both original masked-node norms and geometric-region trapezoid
   norms are saved, exposing boundary-weight effects separately.
2. A split-native-cell Gauss L2 reference on the exact geometric training
   and holdout rectangles. In a fixed native cell, the bilinear observer G
   is bilinear and the Hermite observer G is at most bicubic: moving-query
   rates are affine, multiplying derivatives of their same interpolants.
   Thus G squared has degree at most six in each axis. Four-point tensor
   Gauss integrates it exactly up to floating-point roundoff when native
   cell and core/holdout boundaries are split. This is an integration
   reference for the chosen interpolant, not for an unknown continuum PDE
   profile. Reported infinity norms remain finite sampled maxima, without a
   certified supremum error bound.
3. Native data restriction by strides 1, 2 and 4, retaining the *original*
   exact coordinate values and rates. Restricted Omega and F_X get
   corresponding fixed seven-node derivative matrices, while stride 1 uses
   the saved maintained matrices. Interpolant differences are integrated
   on the full native partition, which also splits the coarser cell
   boundaries. This measures observation-operator sensitivity to discarded
   samples. It does not recompute a coarse PDE, improve the original data,
   or establish PDE spatial convergence.

`ipm_accellab_test_observer_sampling` passed on a manufactured bicubic forcing
with analytically integrated squared norm: the split-cell Hermite L2 error
was `6.661e-16`. All six density/phase variants and all three restrictions
were retained, and all fixed-field scaling ratios matched their prescribed
value. No LU, RHS or PDE was used. The report is in
`result/longtime/20260908_campaign_v1/acceleration_lab/sampling_tiny_tpd3e2b5b1_1cfb_4cf4_bc45_ca5a2a47ec27/`.

### Completed fixed q512 observer comparison and decision

The authorized fixed comparison has now run and exited. Its five-source
strict native audit reconstructed the rejected full rho array **bitwise**;
the maintained derivative matrices matched exactly. On the unchanged
65 x 49 observer, bilinear metrics reproduced the original report to
`1.735e-17`. The Hermite training/holdout L2 ratios were
`.99842016/.99946990`. The fresh tail velocity change remained
`.00104318876`, exceeding the original `.001` gate. The existing 1/16
candidate remains rejected, including its complete earlier failed metrics.
Report: `fixed_observers_tpf210ce8d_3aee_4a8c_93dd_c1add1567b86/` beneath
the campaign acceleration output directory.

The subsequent field-only split-cell Gauss analysis gave bilinear
training/holdout L2 ratios `.99977942/.99931757` and Hermite ratios
`.99838177/.99946059`. Against the *same bilinear interpolant's* Gauss norm,
the original coarse-node baseline norm underestimated training L2 by
17.2468%, while the candidate norm overestimated it by 10.4603%. These
opposite sampling errors quantitatively account for the original coarse
ratio `1.33452`; all original measurements remain in the report. Hermite
ratios were much less sensitive across the six fixed sampling variants.
This diagnoses observation and integration sensitivity; it does not remove
the independent tail failure or prove PDE spatial convergence.

The Hermite Gauss norm decrease was `9.94396e-6` (0.1618%). Restricting the
same baseline data to every second native node changed the observed
training residual by `.000346355`, about 34.8 times that decrease. This is
an interpolation sensitivity estimate, not a rigorous fine-grid error
bound. Thus the tiny secant improvement does not justify consuming another
q512 computation window. A separate fixed 1/32 proposal is registered as
`prepared_deferred_not_run`; it cannot relabel the old rejection. Complete
sampling results are in
`observer_sampling_tpb762df4b_07a1_4d12_a848_197962764128/`, and the numerical
decision and unchanged gates are saved separately in
`observer_error_assessment_3ac826714a364d6bb1e49120f3428274/report.json`
and `assessment.md`.

### Independent continuous peak gauge research

`ipm_accellab_hermite_peak` defines the maximum of one fixed C1 cubic
Hermite interpolant of wall Omega, with slopes from maintained `Dx`.
`ipm_accellab_hermite_amplitude` derives the pure instantaneous amplitude
rate using independently assembled base forcing. No target-P feedback or
trajectory is introduced. The default explicitly rejects multiple, nearly
tied, flat and boundary maxima.

The manufactured test and verified native 65 x 33 switch-frame study have
passed without a Poisson solve or time integration. On a prescribed linear
field path through the old sampled-node selector tie, the old quadratic
peak difference tends to `.00139735646744`, while the Hermite maximum
difference tends to roundoff. The new instantaneous full-RHS derivative
cancels to at most `1.388e-16`, but the original trajectory does not conserve
the new observable. The later finite-step study is recorded separately below.
The derivation, actual results, invalid cases, remesh limitations and a
separate proposed tiny integration protocol are in
[CONTINUOUS_PEAK_GAUGE.md](CONTINUOUS_PEAK_GAUGE.md).

The subsequent authorized tiny integration, coarser time calibration,
common-physical-time comparison, large-window tolerance audit and latest
q512 instantaneous evaluation are complete. Their results and limits are
recorded in [CONTINUOUS_PEAK_TIME_TRIAL.md](CONTINUOUS_PEAK_TIME_TRIAL.md).
The C1 rule preserved its actual tiny peak to approximately `1e-14`, while
the old 65 x 33 quadratic rule retained a `1.63e-3` selector drift. A stable
fourth-order full-field convergence rate was not established. At fixed
physical time, the two rules' density difference was at most `1.5e-11` on
65 x 33, far below the two-grid spatial difference; rescaled-field error
and physical-field error must be distinguished.

The old window-scaled helper is frozen as
`ipm_accellab_hermite_peak_window_scaled_v0`. The current helper uses local
cell/polynomial scales and coordinate ULPs after a registered audit preserved
nine old false rejections on `1e6` windows. All 2,246 accepted primary C1
frames retained their peak values, locations and cells bitwise. The
step-3326 q512 probe found the new instantaneous amplitude rate
`-.350115814304`, versus native `-.350430869832`, with true proposed HPrime
`2.776e-16`; it did not integrate a new trajectory or change production C.

### Consistent continuous inner geometry and true shape derivative

The new read-only `ipm_accellab_continuous_inner_rates` uses the same C1
wall peak, connected `.9` wall roots, and same-peak tensor vertical trace.
It retains the real `PPrime` and derives all geometry rates from the
original forcing. The full normalized pullback derivative includes
translation in both the observation coordinates and the vertical root
equation. It fits no geometric rate and changes no evolution equation.

Nonuniform polynomial tests, coordinate/amplitude covariance, interface
and plateau guards, and a real 65 x 33 original-PDE time-difference check
have passed. The actual old native step-2270 cache and a separately
qualified later fresh finite-box cache were then analyzed with no new
Poisson or RHS. Split-cell Gauss norms separate the core and surrounding
holdout exactly for the chosen interpolant. Original node-masked weights
and separately corrected same-boundary trapezoid results are both kept.

At old canonical tau `4.78320`, the complete residual was 4.29% of the
normalized Eulerian evolution in the core and 10.47% outside. The later
fresh case at parent-equivalent tau `7.24320` retained 3.11% and 6.29%,
after explicit clock and coordinate unit conversion. Its absolute core
shape residual did not decrease relative to the old state; different
grids, boxes and IVPs also prevent estimating a temporal convergence rate
from this pair. Equations, entry points, actual results, preserved failed
tests and event limits are in
[CONTINUOUS_INNER_DERIVATIVES.md](CONTINUOUS_INNER_DERIVATIVES.md).

### Tiny Jacobian-free profile candidate

`ipm_accellab_test_jf_profile` now executes a projection-free whole-field
Krylov/Gauss--Newton candidate, with genuine continuous shape derivatives,
independent Gauss/holdout/native-FD checks, physical tail/peak hard gates,
and original-PDE relaxation to the same physical time on 49/65 grids.
The first unconstrained trial remains rejected. A separately registered
tail-constrained coefficient ellipsoid passed the small tests and retained
about 0.3% core residual improvement after relaxation. It did not modify
the original initial-value trajectory or establish long-time acceleration;
the observed two-grid baseline difference is larger than that benefit.
See [JF_PROFILE_TRIAL.md](JF_PROFILE_TRIAL.md) for the precise finite
observation window, costs, gates, original failures and limits.

### Actual continuous geometric rate constraints

`ipm_accellab_probe_four_gauge` tests independent tiny half-plane cases. A
linear continuum-generator estimate fails to enforce the actual WENO
constraints; a small damped Newton solve using the complete transport
assembly succeeds at the instantaneous RHS level on 49/65 grids. The
isotropic three-constraint alternative preserves prospective fixed-LU reuse
and reports its unconstrained vertical motion. No time integration or
production rule changed. Nonzero discrete coordinate-covariance defects and
the initial configuration/export failures remain explicit in
[CONTINUOUS_MULTI_GAUGE_RHS.md](CONTINUOUS_MULTI_GAUGE_RHS.md).

`ipm_accellab_observe_audited_cache` and
`ipm_accellab_compare_c1_cache_reports` read already-completed RHS caches and
their strict source reviews without another LU or PDE step. The comparison
uses a union of both native partitions for exact bicubic Gauss L2 integration
in each complete normalized chart, and converts actual derivatives to common
physical time before subtracting them. It preserves the source mesh decisions.

The independently registered smooth connected integral-width observer uses
theta=.5,m=4 with fixed .3/.7 holdouts. Its true derivative integrates FX and
the actual PPrime directly; the vertical trace also includes aPrime. It makes
the width-rate measurements much more stable on the three late spatial
caches, but complete residual-field differences remain 7--14% and do not
decrease monotonically. It cannot accept prior failed candidates. See
[SMOOTH_WIDTH_OBSERVER.md](SMOOTH_WIDTH_OBSERVER.md) for the derivation,
unchanged original evidence, fixed units, and actual no-LU validation.

The actual late smooth-three-rate solve remains rejected: its last rate
candidate missed the unchanged `5e-11` constraint gate, increased primary
core/holdout residual norms by factors `20.31/6.81`, and failed the full-field
and tail coordinate-covariance checks. Its one-sided directional checks
also exposed an unreliable local linear model. All trials and the older
90%-width failure are retained in
[SMOOTH_THREE_RATE_PROTOCOL.md](SMOOTH_THREE_RATE_PROTOCOL.md).

A separate face-local LF alpha prototype passed the registered mapped
kernel tests and then reproduced the actual late cached full 2D RHS
bitwise in its original global mode. At unchanged native rates its local
variant changed the primary smooth residual field by `5.35%/1.31%` in
core/holdout, with slightly increased norms. The 3.16 s cached operator
audit required no LU or time integration. It establishes discretization
sensitivity, not acceleration or production qualification; numerical
boundary changes and all earlier rejections remain explicit in
[LOCAL_LF_2D_AUDIT.md](LOCAL_LF_2D_AUDIT.md).

### Late cropped chronological holdout (fixed rank one rejected)

[LATE_INNER_TEMPORAL_HOLDOUT.md](LATE_INNER_TEMPORAL_HOLDOUT.md) records a separate training-only run on eight tau9.2432–10.6432 native endpoints and a frozen-model evaluation on eight later tau10.8432–12.2432 endpoints. All parameters were registered before loading the new future arrays; later states had prior scientific plots, so it is explicitly a chronological holdout rather than a new prospective collection. Rank1 loses to linear trend at the later wall/vertical horizons and fails its registered screen. Final wall/vertical L2 are0.6671%/0.1576%, versus linear0.3115%/0.08056%. Geometric scalars forecast better, but phase error is1.55 actual wall widths. No predicted state, time, limit or PDE acceleration is accepted. Four shape-model coefficient forecasts agree with an independent direct-matrix evaluation to9.99e-16. Figures, full scalar/shape tables, native metadata and hashes are retained.

[H64_DIRECTIONAL_OPERATOR_PREPARATION.md](H64_DIRECTIONAL_OPERATOR_PREPARATION.md) records two operator-only input bundles on the same actual H64 step2353 field, members7/8 at961×321 and641×481. The first original planned pair passes each frozen-input gate; default-false helper tests do not build factors. Native transfer, resource and PDE qualification remain unavailable until separately authorized real tests.

### Actual late JF qualification preparation only

[JF_LATE_CROPPED_QUALIFICATION_PROTOCOL.md](JF_LATE_CROPPED_QUALIFICATION_PROTOCOL.md) confirms that the old tail-ellipsoid method has only49/65 Gaussian evidence, not late cropped evidence. Strict no-LU step9205 preflight found the legacy49×25 physical matching grid has zero wall-core samples and only one vertical-core sample. Its actual C1 shape/config/lineage and legacy density/source samples are saved; velocity/F/PPrime are explicitly absent. A new one-candidate protocol preserves the old20-RHS construction and gates, adds source-resolving matching/phase and directional nonsmoothness qualifications, and leaves every large execution flag false. No old failure or main C rule changes.

### Actual original-t0 complete exterior-tail operator comparison

[T0_FULL_TAIL_OPERATOR_RESULTS.md](T0_FULL_TAIL_OPERATOR_RESULTS.md) analyzes root's completed H8/H32 fixed-operator trials (one LU, two Poisson, zero PDE per case). Actual Delta c_l matches independent exact exterior integration to1.43e-10/8.81e-10. Corrected c_l cross-case difference shrinks253× but retains essentially the original interior-method gaps; no pure-Poisson attribution or evolved-tail closure follows. Saved-field harmonic checks and independent smooth interior Green potentials/velocities passed. The continuous tail c_omega response has leading isotropic-strain cancellation, so A/H cannot be directly added to either original rule. All original experiments, inputs and first analysis reader failure remain preserved.

[LATE_JFNK_BASELINE_AND_METHOD.md](LATE_JFNK_BASELINE_AND_METHOD.md) separates a genuine square/bordered matrix-free Newton search from the rectangular inner observer fit. Conditional F[a rho]=a²F yields a bare-Newton amplitude shortcut, so actual peak/tail and G/P gates are essential. The quadratic constraint derivative contains D²P[v,F], and J does not generally preserve its tangent space. A default-false step9205 helper strictly reads the unchanged native source; its guarded future branch uses one saved-grid build and one original flow, then releases the factor before an exact array-cache readback and C1 residual observation. No candidate or PDE is enabled. A possible20-call future ledger remains a separately authorized protocol, not executed acceleration evidence.

[FULL_QUADRANT_FROZEN_TAIL_PROTOCOL.md](FULL_QUADRANT_FROZEN_TAIL_PROTOCOL.md) prepares a fixed r[4,5], full-first-quadrant extension on actual H64 steps2213/4291/5874. It uses global band denominators and positive native weights, separate supported Hermite queries and an explicit analytic x=0 trace, and preserves all original wall-sector numbers. The original32-point model must first pass a new full-angle quadrature audit; no order or stable-formula substitution is allowed after failure. This preparation has not started MATLAB or loaded CP arrays. It also derives an exact angular potential for the existing1/r² manufactured perturbation as a minimal later full-boundary Green-trace check, without boundary injection or an evolved-closure claim.
