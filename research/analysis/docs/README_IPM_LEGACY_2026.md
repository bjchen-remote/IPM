# IPM solver on the upper half-plane

> Archived 2026-08-30. This document preserves the former solver interface and
> file map. The maintained contract is
> [STRUCTURE.md](../../../STRUCTURE.md).

The underlying physical equation is

```text
rho_t + u dot grad(rho) = 0
u = grad_perp (-Delta)^(-1) d_x1 rho
```

on a nonperiodic rectangular truncation of the upper half-plane. The old
Boussinesq files are not called. The no-argument run uses the finite-degenerate
strong-tail datum

```text
rho_0 = |x|^8/(1+|x|^8+y^8)
```

on a speed-oriented `513 x 257` grid over `[-1e5,1e5] x [0,1e5]`, with the
exact `X=1` transport anchor,
the wall-omega-peak amplitude gauge, nodal nonuniform WENO transport, and a
strict-`1.5` two-axis lattice grid. An explicit option structure can still
select the exact
physical equation with `c_l=c_omega=c_r=0`. Both paths use
`symmetryMode='double_odd_omega'`: only the first-quadrant source is integrated,
and odd images across both axes generate the other three quadrants.

## Run and verify

```matlab
report = ipmtests.baseline.suite();
result = ipm.solve();
```

`ipmtests.baseline.suite` is the full layered runner. Its report groups configuration,
numerical-core, transport, four-case smoke, and integration checks. For a
quick end-to-end check use `ipmtests.baseline.smoke`; use `ipmtests.baseline.integration`
to run the former monolithic manufactured-solution and structural suite by
itself.

Read-only polar profile research is kept outside the solver:

```matlab
addpath('analysis');
fit = ipm_fit_polar_profile( ...
    'result/production/ipm_double_odd_dynamic.mat', ...
    struct('makePlot',true, ...
    'plotFile','result/verification/polar_profile_dynamic.png'));
```

The computation remains centered at the symmetry origin; the analysis only
labels that point as `(1,0)`. See `analysis/POLAR_PROFILE_RESEARCH.md` for the
polar steady equations and the symmetric/asymmetric profile conditions.

To mesh the exact local degree-one acute-fan benchmark about `(1,0)`:

```matlab
addpath('analysis');
fan = ipm_mesh_exact_fan_profile();
```

The default uses a 60-degree fan and writes
`result/verification/exact_fan_local_mesh.png`. Its `c_omega=0` is the
spatial-only gauge: the invariant condition is
`c_l-c_omega=A_0*sin(theta_*)^2/(2*cos(theta_*))`.

The exact physical steady equation also has a punctured-domain Bessel family
that decays at infinity:

```matlab
addpath('analysis');
bessel = ipm_plot_exact_bessel_family();
```

This writes `result/verification/exact_bessel_family_gallery.png` and
`result/verification/exact_bessel_rho_x_detail.png`. The center singularity
is essential: this is an analytic exterior/fan benchmark, not a smooth
localized blow-up profile or a dynamic-scaling steady state.

For the actual upper half-plane, compare the global wall-fed mode with the
impermeable mode by using

```matlab
halfPlane = ipm_plot_exact_bessel_family(struct( ...
    'thetaStar',pi,'modeOrders',[0,1], ...
    'modeCosine',[1,0],'modeSine',[0,1], ...
    'detailIndex',2, ...
    'galleryFile','result/verification/exact_bessel_halfplane.png', ...
    'detailFile','result/verification/exact_bessel_halfplane_rho_x.png'));
```

The `nu=0` member has wall through-flow. The `nu=1` sine member has zero
normal wall velocity but necessarily retains the pole at `(1,0)`.

For the smooth near-Bessel dynamical screens:

```matlab
addpath('experiments','analysis');

k0 = ipm.solve(bessel_two_scale_case('screen_perturbed'));
k0Diagnostics = ipm_analyze_bessel_two_scale(k0);

k1 = ipm.solve(bessel_two_scale_case('screen_k1_perturbed'));
k1Diagnostics = ipm_analyze_bessel_k1_two_scale(k1);

verify_bessel_two_scale_initial_data();
```

The K0 case moves its pole below the wall, making the density smooth in the
closed computational half-plane; `ipm.solve` still replaces the analytic
through-flow velocity by its impermeable Poisson solve. The regularized K1
main component has zero wall trace, while evolution adds an explicit 0.5%
signed-Gaussian wall monitor for the current generic safety audit. These are
short-time shape-retention tests, not evidence that smooth data are attracted
to a Bessel two-scale state. The complete derivation and numerical audit are
in `analysis/reports/IPM_BESSEL_TWO_SCALE_REPORT.tex` (compiled copy:
`output/pdf/ipm_bessel_two_scale_report.pdf`).

This opens a live evolution figure and writes
`result/production/ipm_double_odd_dynamic.mp4`. With no argument (or `[]`),
`ipm.solve` obtains its maintained overrides from `ipm.config.activeCase.m`:
`initialCondition='degenerate_primitive'`, `degeneratePower=8`, box
`[-1e5,1e5] x [0,1e5]`, grid `513 x 257`, physical final time `2`,
`lengthGauge='transport_anchor'` at `X=1` with
`lengthScaleGain=1`, `cOmegaGauge='wall_omega_peak'`, nodal nonuniform WENO
transport, and the strict-`1.5` lattice mesh. Disable
either output with:

```matlab
opts = ipm.config.activeCase();
opts.livePlot = false;
opts.writeVideo = false;
result = ipm.solve(opts);
```

The active `maxDt=1e-2` is only a loose canonical-time ceiling. The actual
step is selected stage by stage by `cfl=0.35` and then clipped at the requested
terminal time. A logarithmic-primitive smoke run on the former
`[-256,256] x [0,256]` box reached physical time `0.005` in three steps;
`maxDt` was not used as a fixed step size.

For comparisons of the maintained `k=8` datum at the same physical instant,
start from the complete active option set:

```matlab
opts = ipm.config.activeCase();
opts.physicalFinalTime = 2;
result = ipm.solve(opts);
```

If `finalTime` is omitted, this disables the canonical-time stop. This is the
recommended path for gauge, grid, and box comparisons.

Full field snapshots are not retained by default because the video does not
need them. Add `'storeSnapshots',true` only when the in-memory fields are
required for later analysis.

For a small generic-fallback debug run (not the maintained `k=8` active
case):

```matlab
opts = struct('nx',65,'ny',33,'finalTime',0.1, ...
    'saveResults',false,'makePlots',false);
result = ipm.solve(opts);
```

The generic fallback selected by a nonempty partial option structure is

```text
rho_0 = |x_1|^k/(1+|x_1|^k+x_2^k),   x_2>=0, k=4.
```

It is not the no-argument maintained datum. It is the even `x_1` extension of
the retained first-quadrant trace, so
`omega=d_x1 rho` is odd in `x_1`. Its continuation across `x_2=0` is odd and
therefore discontinuous, as explicitly requested; the stored bottom row is
the nonzero upper trace. The four-image kernel makes `omega` and `psi` odd in
both axes, yielding `u(0,0)=0`.

Because this accepted jump already contributes a distributional axis delta to
the full-plane `rho_x2`, reported gradients are upper one-sided regular
gradients and do not include that prescribed singular measure.

Change `k` with:

```matlab
opts = struct('initialCondition','degenerate','degeneratePower',6);
result = ipm.solve(opts);
```

The maintained no-argument datum is the logarithmic primitive family with
`k=8`:

```text
D = 1+|x|^8+y^8,
rho_0 = (1/8)*log(D/(1+y^8)),
q_0 = d_x rho_0 = sign(x)*|x|^7/D.
```

Thus `q_0>0` throughout `x>0, y>=0` and `q_0<0` throughout
`x<0, y>=0`. The origin has a finite seventh-order degeneracy,
`q_0(x,0)~x^7`, rather than a dead interval. Along a generic interior ray
`q_0=O(r^-1)`, and on the wall `q_0(x,0)~1/x`. Hence `q_0` tends to zero in
every far-field direction, but the full gradient has a generic `r^-1` tail
and its two-dimensional `L2` energy diverges logarithmically. The price of a
one-signed `1/x` wall tail is unavoidable logarithmic growth of `rho_0` along
the wall: a bounded primitive cannot have a positive nonintegrable derivative.
During evolution,
`history.positiveHalfPlaneOmegaMinimum` and
`history.positiveHalfPlaneOmegaNegativeRatio` expose any loss of the initial
positive-half-plane sign without projecting the numerical solution onto a
different PDE. An explicit
physical-coordinate replay changes only the rescaling mode:

```matlab
opts = ipm.config.activeCase();
opts.rescalingMode = 'physical';
result = ipm.solve(opts);
```

The former CCF-like datum is now a named historical/opt-in experiment rather
than the no-argument default:

```matlab
addpath('experiments');
opts = ccf_heavy_tail_case('long_dynamic');
result = ipm.solve(opts);
```

Before its monotone saturation scale, that family's wall derivative
approximates `1_{x_1>1}/sqrt(x_1-1)`. Set `ccfTailCutoffScale=Inf` only for an
explicit uncut-tail comparison. Because the unregularized limit is already
singular at the initial time, use the supplied `box_large` and `epsilon_fine`
cases before interpreting growth as an evolution-generated singularity. Its
historical `gradient_energy` amplitude gauge is a comparison normalization,
not the active default.

## Dynamic rescaling and physical tracking

`ipm.solve()` uses the maintained dynamic `k=8` rational setting by default. To
evolve the same datum in physical coordinates, change only one active option:

```matlab
opts = ipm.config.activeCase();
opts.rescalingMode = 'physical';
result = ipm.solve(opts);
```

For a smaller dynamic debug run that still starts from the maintained `k=8`
datum:

```matlab
opts = ipm.config.activeCase();
opts.nx = 129;
opts.ny = 65;
opts.xlim = [-8,8];
opts.ymax = 8;
opts.saveResults = false;
result = ipm.solve(opts);
```

The IPM-specific rescaled equation is

```text
R_tau + (U+c_l X+c_r e_1) dot grad(R) = c_omega R.
```

The code evolves its conservative form with one selected full-velocity
nonuniform conservative flux,

```text
R_tau = -div(R*(U+c_l X+c_r e_1))+(2*c_l+c_omega)R.
```

The physical wall is closed; the left, right, and top artificial density
boundaries are open. This preserves constants under coordinate dilation and
avoids the negative density produced by the former centered scaling split.
The maintained case uses `transportScheme='weno5_nonuniform'`; the generic
fallback is `muscl_minmod`. WENO gives fifth-order smooth face reconstruction
on arbitrary nonuniform nodal grids, while the complete conservative operator
remains second order under the solver's nodal/control-volume contract. Its
grid-dependent coefficients are rebuilt automatically after every remesh.

The maintained double-odd length gauge fixes the combined wall transport at
the computational point `X=1`:

```text
U_1(1,0)+c_l=0.
```

This uses `lengthScaleGain=1` and admits no adaptive width-rate correction;
connected `10%/50%/90%` counts request remeshing instead. Double-odd symmetry
fixes the horizontal frame:

```text
c_r = 0,  X_shift = 0,  u(0,0)=0.
```

The former peak-restoring `c_r` remains only in explicit
`symmetryMode='half_plane'` comparison runs. `transportAnchorX` must lie
strictly inside the positive computational domain. The principal active
options are:

```matlab
opts = struct('lengthGauge','transport_anchor', ...
    'transportAnchorX',1,'lengthScaleGain',1, ...
    'cOmegaGauge','wall_omega_peak');
result = ipm.solve(opts);
```

Changing `cOmegaGauge` changes the rescaled representation and internal clock,
but not the reconstructed physical solution. The common positive time limiter
multiplies both physical transport and every canonical rate, so it preserves
the `X=1` zero exactly.

For the exact active `k=8` logarithmic-primitive datum, the zero-step
canonical slope is approximately `-0.182`, so `X=1` is attracting. The next
positive-side zero is near `X=1.333` and is repelling; initially these two
zeros delimit the anchored right-side transport basin. Changing the amplitude
gauge cannot reverse the signs because it only supplies a common positive
time reparameterization.

For an explicit multiscale comparison, `lengthGauge='wall_density_width'`
uses the wall-density 20% and 60% crossings instead of local strain. It is not
the default: tests show that freezing this broader material transition can
still leave the narrower 90% `rho_x` core under-resolved.

Adaptive feedback is enabled by default. Its targets are capped relative to
the initial connected peak, so it does not command indefinite spreading:

```matlab
opts = ipm.config.activeCase();
opts.targetPeakPoints = 64;
opts.targetPeakExpansion = 8;
result = ipm.solve(opts);
```

The default amplitude gauge preserves the tracked positive wall-omega peak:

```text
c_omega = -d_X1 B(X_*,0)/d_X1 R(X_*,0),
R_tau = B+c_omega R.
```

The tracking window stays around the initial locked peak so a remote secondary
peak cannot take over the gauge. For comparison, the full two-dimensional
gradient gauge is available with `cOmegaGauge='gradient_energy'`. The old
point-strain choice is available with `cOmegaGauge='strain_point'`; monitor
`history.strainCondition`, since a value approaching zero means its
`c_omega` pole is a gauge failure. `cOmegaGauge='none'` sets `c_omega=0`.

The active target counts are `[96,64,48]` for the three connected levels,
with hard safety floors `[16,10,8]` and connected-area floor `12`. The partial
generic fallback retains `[64,32,16]`, `[16,8,4]`, and floor `12`. In physical mode,
adaptive monitors track both the paired `x_1` peaks and the wall-normal
`x_2` profile. The monitor combines the connected 10/50/90-percent intervals,
weights them by requested cells per physical width, and uses a graded
origin-to-peak bridge to avoid an unresolved central gap. Each axis triggers
from its own deficit and both axes may move together. Every candidate starts
from the original stretched axes; locally inserted quadrature nodes resolve
windows narrower than the base center spacing. Pure PCHIP transfer preserves
the density range and one-sided wall trace. A
physical remap is rejected if it loses more than 2% of either 90% core or
changes the wall `rho_x1` maximum by more than 3%.

```matlab
opts = struct('adaptiveRemesh',true,'remeshSafetyTrigger',0.75, ...
    'remeshTargetSafety',0.4,'remeshMaximumCellRatio',6, ...
    'remeshMinimumSpacingFactor',1e-4,'maxRemeshes',300);
result = ipm.solve(opts);
```

This exposes a hard resolution limit without accepting interpolation-induced
derivative spikes. Dynamic mode retains its bounded horizontal row-mass
correction as a separate comparison path.
`result.history` records
`peakGridPoints`, `supportGridPoints`, `coreGridPoints`,
`resolutionFactor`, `widthCorrection`, `widthRateScale`, `widthControlMode`,
`safetyFactor`, `remeshCount`, `minimumDx`, `minimumDy`, `peakDx`, `originDx`,
and `maximumCellRatioX/maximumCellRatioY` (strictly adjacent-cell ratios),
plus diagnostic-only `globalGridRatioX/globalGridRatioY`,
plus horizontal and wall-normal connected core counts and widths,
`wallPeakLocalMaxima`, `wallPeakTVRatio`, `positiveWallNegativeRatio`,
`originVelocity`, `rhoEvenDefect`, `omegaOddDefect`,
`c_lNominal`, `c_l`, `omegaGaugeResidual`, and all gauge residuals. Check
`history.trackedPeakOffset`, `travelingWaveResidual`, `wallPeakX`, and
especially the rescaled `history.wallPeakWidth` to confirm that the extremum
neighborhood is opened on the computational grid. `physicalWallPeakWidth`
records the corresponding reconstructed physical width. `history.c_r` and
`history.X_shift` remain zero in active mode and audit the symmetry-fixed
frame. `result.finalQuality` collects the terminal safety, range,
oscillation/sign, grid-ratio, and two-axis core contracts. The physical
reconstruction is

```text
x = (X-X_shift)/C_l,  rho_physical = R/C_omega,
u_physical = U/C_omega,  dt/dtau = C_omega/C_l.
```

`result.t` is physical time, `result.tau` is canonical rescaled time, and
`result.normalizedTau` is the internal rate-limited clock. `opts.finalTime`,
`outputEvery`, and `maxDt` are canonical-time quantities. In physical mode
`C_l=C_omega=1` and all clocks agree. The physical high-accuracy defaults are
`maxDt=2.5e-4` and `outputEvery=0.002`, selected by the saved time-step study.

## Accuracy and blow-up checks

The speed-oriented no-argument run uses `513 x 257`; it is suitable for
interactive exploration but is not by itself adequate evidence for a
pointwise derivative maximum. Quantitative claims should repeat the active
case at `1025 x 513` (and with a matched box refinement). Independent
hyperbolic-sine maps automatically target `dx=dy=0.01` near the origin/wall,
while the artificial boundary remains far away. Set `gridStretch` explicitly
to override the automatic map and set `gridStretchAutomatic=false`, or set
`gridMode='uniform'` for comparison. Automatically generated options retain
`gridStretchAutomatic=true`, so changing `nx`, `ny`, `xlim`, or `ymax` in a
reused `result.opts` recomputes the stretch instead of silently reusing the
old map.

The implementation follows the old grid/Poisson contracts without depending
on its Boussinesq classes: differentiation is done on a uniform reference
coordinate and divided by the analytic map Jacobian; the Poisson stiffness
uses the same Jacobian/control-volume weights; and the three artificial sides
receive four-image quadrant Green data computed from first-quadrant sources.
The legacy half-plane kernel remains available with
`symmetryMode='half_plane'`. Only the physical wall has `psi=0`.

`ipmtests.baseline.suite` tests second-order convergence of both derivatives, Poisson, the
complete signed velocity map, and conservative nodal transport. The included
transport suite separately checks fifth-order WENO face reconstruction,
closed mass conservation, constant states, jump behavior, and a nonuniform
SSPRK3 step-advection comparison against MUSCL. Arbitrary-grid/remesh
regression also checks
quadratic derivative exactness, monotonic nodes, range preservation, peak
refinement, and rejection of repeated no-improvement remeshing. It also checks
both Green boundary kernels, paired parity preservation, zero origin velocity,
conservation,
incompressibility, wall no-penetration, the even upper trace, the dynamic
gauges, and sensitivity to doubling the box.

The primary growth signal is `result.history.physicalGradInf` versus
`result.history.physicalTime`. `result.blowupFit` excludes states that fail
resolution, density-range, one-peak/total-variation, positive-wall-sign, or
grid-smoothness contracts, then compares a power law with exponential growth.
The primary result uses 40 points; `windowPoints`,
`windowT`, `windowP`, and the corresponding R-squared/candidate arrays expose
20/40/80/160-point window drift. A numerical fit is not proof: repeat with
finer grids and larger domains and check the far-boundary diagnostics and all
gradient components.
Also require `physicalRangeViolation` to remain small; the solver stops at the
default `rangeStopTolerance=5e-3` because physical IPM preserves density range.
Compare two completed runs in the common peak frame with, for example,

```matlab
comparison = ipm.output.compare(resultSmall,resultLarge,[-2,2,2]);
```

This reports local density/physical-omega errors separately from laboratory
peak translation.
The existing matched time-two box/grid records were produced with the former
half-plane kernel. They are comparison history, not convergence evidence for
the active four-image problem. Regenerate matched double-odd runs before
drawing a blow-up conclusion, and do not infer blow-up from a maximum trace
alone.
For possible multiscale concentration, also compare
`physicalGradientEnergy`, `physicalGradientEffectiveLength`, and
`gradientInnerFraction/gradientNearFraction/gradientOuterFraction`; these use
both components of `grad rho`, not `d_X1 rho` alone.
The uncut rational datum has a direction-dependent far field and converges
slowly under box expansion. Check `farBoundarySourceRatio`, the raw velocity
ratio, and `frameVelocityU1`; the raw velocity ratio need not decay for this
datum. Do not treat a single ratio as proof: repeat with a larger box while
retaining the same center spacing. The double-odd images cancel the uniform
horizontal component at the symmetry origin, which is audited through
`u(0,0)=0`, `c_r=0`, and `X_shift=0`. Use `comovingPhysicalX` for local-shape
comparisons.

Important options are defined in `ipm_defaults.m`. Physical results save by
default to `result/production/ipm_double_odd_physical.mat`; explicit dynamic
runs use `result/production/ipm_double_odd_dynamic.mat` unless `resultFile` is
supplied. Reproducible accepted-grid option sets live in `experiments/`; old
Boussinesq source is isolated in `legacy_boussinesq/`, and result files are
classified under `result/production`, `result/verification`, and
`result/archive`.
