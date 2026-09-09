# CCF-like heavy-tail experiment record

Date: 2026-08-24

## Datum and verification

This is the maintained named CCF comparison family,
`initialCondition='ccf_heavy_tail'`; it is no longer the no-argument
`ipm.solve()` datum. Its resolved
onset has `epsilon>0`, its density saturates monotonically after the scale
`R`, and its pre-cutoff wall derivative approximates
`1_{x_1>1}/sqrt(x_1-1)`. The uncut choice remains available explicitly with
`ccfTailCutoffScale=Inf` but was not used for the evolution claims below.

`ipmtests.baseline.suite` passed. On the fixed smooth portion of the tail, the last three
reported nonuniform `D_x` convergence orders were `4.279, 2.834, 2.441`.
The full Poisson, Biot-Savart, conservation, symmetry, remesh, and Code
Analyzer checks also passed.

The optional analytic pre-remesh evaluates the formula again after moving the
nodes; it does not interpolate a coarse initial peak. Initial diagnostics are:

| case | epsilon | grid | core points | safety | peak dx | far-source ratio |
|---|---:|---:|---:|---:|---:|---:|
| screen | 0.08 | 513 x 257 | 7.865 | 0.509 | 0.01764 | 0.0201 |
| box_large | 0.08 | 1025 x 513 | 8.855 | 0.452 | 0.01586 | 0.00346 |
| epsilon_fine | 0.04 | 1025 x 513 | 8.395 | 0.476 | 0.00846 | 0.0125 |

## Completed evolution comparisons

The accepted verification-scale physical run is `box_large` through
physical time `0.25`:

| diagnostic | final / initial or final value |
|---|---:|
| wall `rho_x` peak | 2.0651 |
| wall FWHM | 0.17109 |
| wall 90% core width | 0.083607 |
| full gradient maximum | 1.14775 |
| final grid safety | 0.52571 |
| remesh count | 8 |
| mass drift | 3.72e-6 |
| density-range violation | 7.95e-4 |
| far-source ratio | 0.00346 -> 0.01339 |
| trusted through time | 0.2481 |

The same small-box dynamic comparison reaches physical time `0.25`. Its
physical peak grows by `1.8701`, while the rescaled peak changes by only
`1.0113`; physical FWHM and core width fall to `0.21962` and `0.21651` of
their initial values. It uses one analytic initial remesh plus three evolution
remeshes. However, its far-source ratio grows from `0.0201` to `0.1518`, so it
is not a box-converged result. The `-0.329` physical-mass change is dominated
by flux through the changing finite physical truncation in dynamic
coordinates and is not treated as a blow-up diagnostic.

## Interpretation

The large-box physical result shows strong, nonoscillatory wall concentration
similar in shape to the CCF scenario: the derivative peak roughly doubles
while two connected widths contract by factors of about 6 and 12. A single
dynamic length scale locks the rescaled amplitude but does not lock the much
narrower core, which is consistent with multiscale concentration.

This is not yet numerical evidence of finite-time blow-up. On the trusted
large-box interval, a power fit to the wall peak gives `T=0.496`, `p=1.036`,
and `R^2=0.99836`, but the exponential fit has `R^2=0.99713`; the required
separation is absent and fitted exponents drift across windows. Inverse-FWHM
and inverse-core fits likewise do not give a common stable terminal time.
The next decisive comparisons are the resolved `epsilon_fine` physical run
and a large-box dynamic run, followed by time-step halving if their trusted
intervals agree.

## Long-mesh revision

The first anchored long run stopped at `t=0.084` because two adjacent nodes
on a flat main peak were counted as two local maxima. Their values differed by
only `0.33%`, the connected-profile TV ratio was `1.00012`, mass drift was
`3.95e-7`, and the far-source ratio had fallen to `0.00299`. The stop was
therefore classified as a diagnostic false positive. The maintained contract
now requires multiple maxima and TV ratio above `1.02` simultaneously; the
negative-lobe stop remains independent.

The corrected ratio-6 long mesh then reached `t=0.2901`. It stopped with two
significant maxima and TV ratio `1.0272`; simultaneously the 90% core had only
`8.08` cells, safety was `1.238`, and both adjacent-cell ratios had saturated
at `6`. Far-source ratio decreased to `0.00147`, mass drift was `-2.61e-5`,
and density-range violation was `0.00110`, confirming that the far anchor
worked and the failure was local resolution. The next long-mesh revision uses
a longer anchor transition and an experiment-only cell-ratio cap of `8`.

The ratio-8 run reached `t=0.3207` cleanly (`TV=1`, far-source ratio
`0.00125`, mass drift `-5.97e-6`) but stopped when the horizontal core reached
`4.99` cells. Offline reconstruction showed why no rescue grid was proposed:
restarting every monitor from the original sinh axis predicted only `0.49`
cells in the already nested 90% core. Long CCF cases therefore use an
outer-anchored `current` reference so local refinement can accumulate while
the standard production path remains non-compounding.

An unanchored endpoint probe improved the predicted core to `15.44` cells but
moved outer-band nodes by as much as `112` and raised the far-source ratio
from `0.00125` to `0.866`, so unanchored accumulation was rejected. An
intermediate monotone coordinate-rescaling anchor was tested and was later
superseded by the maintained local splice below. A horizontal-only candidate
may now change the measured vertical core
down to its configured target rather than being rolled back for a sub-target
2% cross-axis change; a moved axis still retains the strict 98% preservation
rule.

The maintained local anchor now uses a subgrid splice rather than a global
coordinate rescaling. When the current adjacent-cell ratio is below `7.2`,
the established grid path is unchanged. Above that threshold, two halo nodes
are admitted and only the exterior transition is relaxed toward the soft ratio
`sqrt(8)`; the feature cells and far nodes remain fixed. At the old clean
endpoint, two offline rescue passes increased the horizontal 90% core from
`5.00` to `8.31` cells, changed the full gradient by only `0.91%`, kept the
far-source ratio at `0.00111`, and preserved the wall trace to interpolation
roundoff. Their PCHIP mass cost was `0.30%`.

The final from-zero `long_physical` run uses `1025 x 513` nodes on
`[-128,128] x [0,64]`, CFL `0.35`, and maximum step `1e-4`. Its remesh schedule
matches the previous clean run through `t=0.3460`; the first on-demand rescue
then lowers the saturated horizontal cell ratio from `8` to `6.195`. The run
reaches `t=0.392521297` in 4174 steps with 26 accepted remeshes, compared with
the earlier limits `0.3207` and `0.386853`. It stops by
`grid_resolution_failure` with `4.976` horizontal core cells. The last state
satisfying every strict trusted-fit condition is `t=0.368068092`, where the
core has `10.307` cells and safety is `0.9703`.

At the hard terminal point the wall peak has grown by `4.624`, FWHM is
`0.0216954`, the 90% core width is `0.00223527`, full-gradient growth is
`2.674`, TV ratio is exactly `1`, far-source ratio is `0.001074`, range
violation is `0.001511`, and mass drift is `-0.006463`. Fits on the trusted
window give local power parameters `(T,p)=(0.5206,1.1438)` for the peak,
`(0.5433,2.8654)` for inverse FWHM, and `(1.1478,9.9486)` for inverse core
width, but all three return `candidate=false`. The incompatible terminal times
and unstable windows do not establish finite-time or multiscale blow-up.

Exact post-PCHIP physical mass corrections were tested and rejected. A global
support correction raised the far-source ratio to `0.147`; direction-separated
row/column corrections produced more than `200x` artificial cross-direction
gradient amplification at the clean endpoint. The maintained physical path
therefore accepts the measured sub-percent mass drift instead of hiding it
with a nonsmooth correction. A genuinely conservative two-dimensional remap
would be a separate numerical method and is the next prerequisite for pushing
substantially past the present mesh limit.

## Legacy-lattice focused-grid trial (2026-08-26)

The next long-run configuration replaces the broad overlapping nested
monitor tails by `remeshAxisProfile='lattice'`. This keeps the useful part of
the old `latgene` construction: a nearly uniform fine plateau over each
connected feature interval and a narrow tanh transition to the fixed far
grid. It does not copy the inconsistent extra cubic term that appears in
`latgene_diff` but not in `latgene`.

The first profile-only comparison requested peak/10%-50%-90% caps
`128/[192,128,96]`, initial
expansion cap `8`, and buffered safety `0.08`. The actual analytic initial
re-evaluation has `71.185` horizontal 90%-core cells. Its global horizontal
spacing ratio is `1.481e3`, while the maximum adjacent ratio is only `2.730`;
the global ratio is diagnostic and never limits or invalidates the run.

Matched `1025 x 513` physical comparisons through `t=0.1` used the same
targets and differed only in the axis profile:

| profile | x/y 90% core cells | x/y adjacent ratio | x/y global ratio |
|---|---:|---:|---:|
| lattice | 33.443 / 55.380 | 2.730 / 3.205 | 1.481e3 / 4.052e2 |
| nested | 30.345 / 34.512 | 8.000 / 8.000 | 2.827e3 / 4.598e2 |

Both runs reached the requested time in 1000 CFL-controlled steps with two
accepted remeshes and nearly identical wall peaks (`2.605` versus `2.609`).
This established the profile choice, but its adjacent ratios above `2.7` were
subsequently rejected as too rough for production.

### Strict-ratio and partial-dynamic revision

The maintained long cases now impose a hard adjacent-cell ratio of `1.5`, a
proactive safety trigger of `0.15`, and up to 12 analytic pre-remesh passes.
Fixing the strict-cap amplitude fallback was essential: when the first
count-satisfying monitor violated `1.5`, the old search discarded every
sub-target deformation. It now retains the strongest legal deformation below
that trial. Two analytic passes give an actual initial horizontal 90%-core
count of `93.562`, adjacent ratio `1.500000`, and global ratio `1.575e3`.

The results in this paragraph used the pre-anchor symmetry-peak secant with
gain `0.5`; they are retained as historical evidence and are not current
`ipm.solve()` output. The maintained `long_dynamic` now uses the fixed `X=1`
transport anchor with gain `1` and the full-gradient-energy `c_omega` gauge.
Dynamic mode computes and remeshes
the wall-normal core as well. At physical time `0.1` the matched hybrid run
has:

| diagnostic | value |
|---|---:|
| horizontal / vertical 90% core cells | 73.563 / 73.999 |
| horizontal / vertical adjacent ratio | 1.500000 / 1.243019 |
| horizontal / vertical global ratio | 2.746e3 / 4.927e2 |
| `C_l / C_omega` | 1.05305 / 1.03183 |
| instantaneous bounded `c_l / c_omega` | 0.500 / 0.304 |
| physical / rescaled wall peak | 2.5991 / 2.5472 |
| density-range violation / wall TV | 0 / 1 |

This shows partial scale absorption rather than a locked or overexpanded
level set. It is a short stability result, not yet a replacement for the old
long physical trusted window.

## Nodal nonuniform WENO integration (2026-08-26)

The maintained case now selects `weno5_nonuniform`; generic option completion
retains `muscl_minmod`. The WENO method reconstructs nodal values, consistent
with the derivative, Poisson, and PCHIP-remap contracts. A finite-volume
cell-average trial was rejected after it created a second wall-derivative peak
and TV ratio `1.024` by physical time `0.0042`.

Independent nonuniform tests give smooth face-reconstruction order `4.701`
and complete conservative-operator order `2.000`. In a step-advection test the
WENO/MUSCL L1 errors are `3.438e-2/5.280e-2`; WENO TV is `2.0002` and its range
is `[0,1]` to roundoff. A matched `513 x 257` active-case run to physical time
`0.01` completed 51 CFL-controlled steps with three analytic remeshes, one
wall maximum, TV ratio `1`, density-range violation `5.84e-12`, and adjacent
x/y ratios `1.492/1.015`. This only validates integration stability; long-time
grid/domain/time refinement remains required for any blow-up conclusion.
