# Fourth-order method and certification

> Contract note (2026-09-05): the maintained structured solver uses config and
> checkpoint schema 4 with `exact_gauge_no_feedback_v1`. Feature cell counts
> may drive diagnostics/remeshing/hard stops, not scaling rates. This certificate
> remains a fixed-grid physical/smooth-solution claim and does not certify any
> dynamic gauge, restoring controller, or trigger-driven remesh path.

This document defines the empirical fourth-order claim for the integrated IPM
solver. It is deliberately narrower than the set of configurations accepted by
the resolver. Generic architecture belongs in `ARCHITECTURE.md`; implementation
history belongs in `CHANGELOG.md`.

The separately selected WENO7/RK6 method has its own independent scope,
evidence, and exclusions in
[`SIXTH_ORDER_RESEARCH.md`](SIXTH_ORDER_RESEARCH.md); nothing in that
certificate widens or replaces the fourth-order claim below.

## Certified scope

The current end-to-end global space-time order claim covers:

- smooth solutions;
- fixed uniform or smoothly mapped grids;
- `spatialDiscretization='high_order'`;
- `transportScheme='weno5_fd'`;
- `timeIntegrator='ssprk54'`;
- homogeneous Dirichlet stream-function data;
- closed outer transport with the conservative physical wall; and
- physical evolution.

The production convergence norm includes every node, including one-sided
boundary closures. The unforced physical gate advances `main_ipm` on nested
grids to the same physical time; it does not use analytic forcing or a
finest-grid surrogate solution.

Separate fixed-grid time-order tests replay the actual isotropic and
anisotropic production steps with open outer transport and smooth,
non-switching scale configurations. They certify fourth-order RK coupling for
density and scale state, including stage-local anisotropic metrics; they do not
claim global spatial convergence for dynamic-rescaling runs.

A reproducible fixed-grid option block is in `README_IPM.md`. High-order
remapping is tested separately and is not required when adaptive remeshing is
disabled.

## Implemented method

- Seven-point local-polynomial first and second derivatives are exact through
  degree six, including boundary rows.
- Tensor quadrature integrates a local degree-five interpolant on every cell.
  The same public quadrature is used by diagnostics, Green integration, and
  conservative high-order remapping.
- The full-order one-sided Poisson closure is nonsymmetric and therefore uses
  LU/direct solves. High-order anisotropic PCG is rejected.
- Transport is mapped finite-difference WENO-Z5 with linewise-global
  Lax--Friedrichs splitting, degree-five ghost extrapolation, a grid-scaled
  smoothness regularizer, and a discrete free-stream correction.
- SSPRK(5,4) uses five explicit stages and a single shared scale-state codec.
  Every anisotropic stage evaluates its own metric aspect.
- Six-point local barycentric remapping uses a horizontal constant correction
  and a wall-normal bubble correction. The latter conserves each column while
  preserving both wall and top traces.
- A high-order remesh candidate whose physical-density range violation exceeds
  `rangeStopTolerance` is rejected. Retry and rollback occur before the state
  is installed.

The transport/source form is important: WENO-FD already returns a discrete
advective form through its free-stream correction. The assembler adds only the
amplitude source to that result; it does not add the spatial divergence source
a second time.

## Verification matrix

Every order gate requires finite errors and both final refinement segments to
exceed order `3.8`. Spatial, physical, Green, and remap gates additionally
require the complete error sequence to decrease monotonically; temporal gates
instead enforce roundoff-floor and fine-reference-adequacy checks. A single
finest-pair cancellation cannot pass the certificate.

| Claim | Executable gate | Required evidence |
|---|---|---|
| Full-node derivatives, Poisson, velocity, WENO, and coupled MMS | `verify_ipm_high_order` | Last two orders above 3.8 on uniform and stretched grids |
| Fourth-order time evolution | `verify_ipm_high_order_time` | Nonlinear ODE, forced IPM, and actual isotropic/anisotropic production steps above 3.8 |
| Unforced global physical convergence | `verify_ipm_high_order_physical_convergence` | Nested `main_ipm` all-node L2/Linf orders above 3.8 |
| Smooth Green boundary quadrature | `verify_ipm_high_order_green` | Uniform/stretched, two symmetries, `kappa=0.5,1,2`, all required tail segments above 3.8 |
| High-order transfer | `verify_ipm_high_order_remap` | Raw/conservative/repeated convergence, conservation, trace, rejection, and exact rollback |
| Empirical stability | `verify_ipm_high_order_stability` | Smooth/nonsmooth CFL sweep, mapped probe, and long physical evolution |
| Complete high-order claim | `verify_ipm_high_order_all('heavy')` | Every row above passes |
| Whole release | `verify_ipm_release('heavy')` | Maintained second-order and complete high-order suites both pass |

The manufactured coupled solution is

```text
psi = A(t) sin(pi(x+1)) y(1-y),
rho = B - A(t)/pi [pi^2 y(1-y)+2] cos(pi(x+1)),
A(t) = 1 + 0.1 sin(2 pi t).
```

It satisfies `rho_x=-Delta psi`, has zero normal velocity on all sides, and has
nontrivial tangential transport at the physical wall. Its forcing is evaluated
analytically at every RK stage rather than reconstructed from a discrete
residual.

## Current evidence

Evidence below was reproduced on 2026-08-30 with MATLAB R2026a. Values are
informative snapshots, not normative constants; the executable gates above are
the authority.

| Layer | Observed result |
|---|---|
| Spatial core | Minimum required tail order `3.912`; minimum final derivative/Poisson/velocity order `5.112` |
| Nonzero Dirichlet Poisson | Minimum final all-node order `6.246` across square/rectangular grids and three metric factors |
| Mapped WENO critical point | Final L2/Linf orders about `5.031/4.932` |
| SSPRK(5,4) conditions | Eight order-condition residuals at or below `1.11e-15` |
| Forced temporal IPM | Tail L2/Linf orders `3.970--4.009` |
| Actual production RK | Density/scale tail orders `3.977--4.134` for nonzero-scale isotropic and anisotropic fixed-rate paths |
| Unforced physical `main_ipm` | L2 orders `6.836/5.815`; Linf orders `6.719/5.565` on `17/33/65/129` grids |
| Smooth uncompressed Green | Minimum required tail `4.017`; minimum final sampled-boundary order `4.657` |
| Remap | Minimum raw/conservative/repeated tail orders `5.148/5.148/3.868`; integrated L2/Linf/wall about `6.14/5.95/5.94` |
| Heavy stability | All-profile contiguous CFL through `2.40`; recommended production CFL `0.50` |
| 480-step physical run | Range/mass/C2/C3 drift `1.266e-3 / 1.472e-15 / 5.155e-4 / 1.587e-3` |

The time tests also replay the actual production `ipm_step_ssprk54` with
nonzero scale rates. Independent tableau oracles check density, scale packing,
stage-local aspect, and final flow reconstruction at roundoff.

## Valid but not globally certified

The resolver intentionally permits controlled mixed-method experiments. These
do not inherit the global-fourth-order claim merely because one high-order
selector is present. In particular:

- `weno5_nonuniform` has a two-layer MUSCL boundary fallback;
- high-order spatial operators with SSPRK3 are not fourth order in time;
- SSPRK(5,4) with legacy spatial operators remains spatially second order;
- open transport boundaries are covered by several operator/time tests but not
  by the unforced global physical certificate;
- high-order remapping passes isolated and integrated transfer gates, while
  trigger-driven adaptive evolution has no global-order certificate.

## Deliberately excluded scope

- General compressed Green data: source cutoff, endpoint omission, logarithmic
  target/source singularity, and fixed-3000 compression need a dedicated
  high-order product-integration treatment. The current Green gate uses smooth
  sources vanishing to sixth order and disables cutoff/compression.
- State-dependent peak/crossing gauge switches: these introduce nonsmooth
  active-set changes. Current scale-order gates use smooth fixed behavior.
- Discontinuous data: a discontinuity cannot converge globally at fourth
  order. Top-hat and high-frequency cases are stability/boundedness probes.
- Very large high-order grids: the nonsymmetric LU closure is a correctness
  backend, not the final scalable elliptic solver.
- Arbitrarily long evolution: the finite stability campaign is empirical and
  is not a mathematical stability theorem.

The next architectural step for broader certification is a compatible
high-order symmetric elliptic discretization and a dedicated singular Green
quadrature, followed by a full adaptive-evolution refinement campaign.

For the distinct strict sixth-order tuple, use the independent
[`SIXTH_ORDER_RESEARCH.md`](SIXTH_ORDER_RESEARCH.md) certificate.
