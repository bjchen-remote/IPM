# Sixth-order method and research certificate

> Contract note (2026-09-05): the maintained structured solver uses config and
> checkpoint schema 4 with `exact_gauge_no_feedback_v1`. Feature cell counts
> may drive diagnostics/remeshing/hard stops, not scaling rates. The sixth-order
> certificate remains restricted to its stated physical fixed-grid scope and
> does not certify dynamic gauges or any retired restoring controller.

This document defines the empirical sixth-order claim for the integrated,
explicitly selected IPM feature on branch `feature/sixth-order-v1`. It is
independent of the maintained second- and fourth-order certificates. Generic
ownership belongs in [`ARCHITECTURE.md`](ARCHITECTURE.md); the fourth-order
certificate remains in [`HIGH_ORDER_RESEARCH.md`](HIGH_ORDER_RESEARCH.md).

The numerical evidence below was reproduced on 2026-08-30 by the isolated
research source at commit `15a10d9`. The integration branch reproduced its
integrated core/time/physical gates on the same date: minimum spatial/WENO7
orders `6.479/5.545`, rooted-tree residual `1.110e-16`, nonlinear-ODE minimum
tail `5.661`, zero production-vs-oracle RK6 anisotropic coupling error, and
source-aware closed rescaled-mass residual `1.995e-16`, with
uniform/fixed-sinh physical minima `7.231/6.821`. A
continuing 12-hour grid search on
`research/sixth-order-global6` is separate research work: no unfinished search
ranking or candidate is presented here as a final optimum.

## Explicit opt-in and certified scope

The resolver baseline and no-argument active case are unchanged. Sixth order
is enabled only by the complete tuple

```text
spatialDiscretization = sixth_order
transportScheme       = weno7_fd
timeIntegrator        = rk6
remeshTransferScheme  = sixth_order
```

`ipm_numerical_tuple_violation` requires all four selectors together; a lone
sixth-order selector is rejected. The preferred helper also installs the
fixed-grid physical-evolution settings used by the certificate:

```matlab
smoothRho = @(X,Y) 2 + ...
    ((pi^2.*Y.*(1-Y)+2)./pi).*cos(pi.*(X+1));
opts = ipm_sixth_order_options(struct( ...
    'nx',129,'ny',129,'xlim',[-1,1],'ymax',1, ...
    'gridMode','stretched','gridStretchAutomatic',false, ...
    'gridStretch',[0.75,0.50], ...
    'initialCondition',smoothRho, ...
    'finalTime',0.04,'physicalFinalTime',Inf, ...
    'maxDt',0.04/64,'maxSteps',66, ...
    'saveResults',false,'storeSnapshots',false, ...
    'makePlots',false,'livePlot',false,'writeVideo',false));
result6 = main_ipm(opts);
```

The helper defaults to a fixed uniform grid, homogeneous Dirichlet
stream-function data, closed conservative transport, physical evolution, no
analytic or adaptive remesh, and `cfl=0.4`. The caller still owns the domain,
resolution, initial condition, and terminal clocks. Changing the boundary,
scaling, grid, or remesh settings can move a valid run outside this
certificate.

The current global space-time claim is limited to smooth solutions on fixed
uniform grids or admitted fixed smooth sinh maps, with homogeneous Dirichlet
stream-function data, closed conservative transport, and physical evolution.
Every reported convergence norm includes the one-sided boundary nodes.

## Implemented method

- Nine-point local-polynomial first and second derivatives, including
  one-sided boundary rows, are exact on nodal polynomials through degree eight.
  The paired full-order Poisson closure is nonsymmetric and uses LU/direct
  solves.
- Mapped degree-seven quadrature integrates in a uniform logical coordinate,
  multiplies by the map metric, normalizes constants exactly, and supplies the
  common integration weights used by elliptic/Green work, closed-transport
  mass projection, norms, and remap corrections.
- Mapped finite-difference WENO-Z7 uses four cubic substencils, linear weights
  `[1,12,18,4]/35`, `tau7=abs(beta0-beta3)`, power `q=4`, linewise-global
  Lax--Friedrichs splitting, and a degree-six four-ghost-node continuation.
  It subtracts the mapped constant-state residual and, for a closed boundary,
  projects the remaining weighted mass defect onto the constant mode.
- The explicit eight-stage Dormand--Prince RK6(5)8M tableau propagates with its
  sixth-order row. Its embedded fifth-order row is retained for diagnostics;
  it is not used as an adaptive controller. The method has SSP coefficient
  zero and therefore makes no SSP, TVD, or positivity claim.
- The sixth-order remap uses width-eight, degree-seven local barycentric
  interpolation in both coordinates. Constant and smooth-bubble corrections
  use the paired mapped quadrature. Unlimited transfer may overshoot a jump;
  range validation, rejection, retry, and exact rollback remain authoritative.

## Private grid-admission policy

The quality limits are an internal numerical policy owned only by
`ipm_sixth_order_grid_policy`; they are not public configuration fields.
Before accepting a complete sixth-order operator or assembling/factoring its
expensive 2-D Poisson matrix, each axis must satisfy:

| Internal diagnostic | Current policy |
|---|---:|
| Maximum adjacent-cell ratio | `<= 1.15` |
| Minimum normalized nine-point stencil `rcond` | `>= 1e-9` |
| Minimum quadrature weight / mean weight | `>= 1e-4` |
| Mapping-metric consistency error | `<= 1e-6` |

Automatic sinh stretching is capped before operator construction. Explicit,
custom, and remeshed axes are rejected if they fail the policy; they are not
silently repaired. The CFL helper recommendation is `0.4`, obtained by
applying a margin to the finite tested sweep. These checks and experiments are
not a stability theorem for arbitrary nonuniform grids.

## Verification entry points

Run the unchanged maintained and fourth-order gates independently, then the
sixth-order certificate:

```matlab
verify_ipm();
verify_ipm_high_order_all('heavy');

verify_ipm_sixth_order_core();
verify_ipm_sixth_order_time();
verify_ipm_sixth_order_physical_convergence();
verify_ipm_sixth_order_all();
```

The integrated all-runner combines two-tail D1/D2/nonzero-Dirichlet
Poisson/velocity refinement, mapped WENO and both LF linear limits, RK algebra,
nonlinear-ODE refinement, nonzero-rate anisotropic production-stage coupling,
the explicit closed rescaled source/mass contract, grid admission,
width-eight interpolation/conservation, an
analytic-sinh-to-custom integrated mass/trace transfer, inadmissible-candidate
exact rollback, and mass-conservative uniform/stretched unforced physical
convergence. The broader
narrow-Green, repeated/refinement/nonsmooth-remap, and empirical-stability
campaigns in the evidence below were executed in the isolated research source
at `15a10d9`; they are evidence, not callable integration-branch gates.
Shared-path changes also require
cross-revision bitwise checks for the default second-order and complete
fourth-order tuples, MATLAB Code Analyzer, and `git diff --check`.

## Recorded research evidence

The executable thresholds require finite, monotonically decreasing errors and
both required tail segments above their declared order floors. These values are
evidence snapshots from `15a10d9`, not normative constants.

| Layer | Recorded result |
|---|---|
| Nine-point spatial core | Uniform minimum tail `6.729`; fixed mild/stress sinh minima `6.485/6.479` |
| WENO-Z7 critical points | Uniform all-node minima `5.628/5.632`; mapped `0.75` boundary/interior minima about `5.545/6.372` |
| Transport invariants | Mapped free-stream residual `1.608e-13`; relative closed-mass residual `4.033e-17` |
| RK6 algebra and smooth time tests | Rooted-tree residual `1.110e-16`; exact semidiscrete WENO MMS tails `5.955--5.981` |
| Unforced physical `main_ipm` | Minimum all-node order `7.231` uniform and `6.821` on fixed sinh `[0.75,0.50]` |
| Narrow Green data | Minimum heavy tail order `6.905750` |
| Width-eight remap | Minimum line/repeated/tensor tail orders `6.924376/6.256759/7.998956` |
| Quick empirical stability | Tested through CFL `0.50`; recommended production CFL `0.40` |

The narrow Green gate uses smooth sources vanishing to eighth order, disables
source cutoff and compression, and checks both symmetries and
`kappa=0.5,1,2`. The remap gate includes 24 smooth round trips, conservation,
boundary traces, discontinuity rejection, and rollback. Neither gate widens
the fixed smooth global certificate.

## Two noncertifying pure-time diagnostics

Two failures remain visible and must not be reinterpreted as passes:

- A nonlinear mapped-WENO fine-reference diagnostic gives first tail orders
  `5.40674` in L2 and `5.27623` in Linf, with finest reference/error ratios
  `0.20095/0.69535`; it reports `sixthOrderCertified=false`.
- The complete fixed-space production WENO-IPM diagnostic gives L2 tail orders
  `0.83541/4.88152` and Linf tail orders `1.02648/4.60958`; it also reports
  `sixthOrderCertified=false` despite finite, monotone, conservative output.

Thus the tableau and the smooth exact-semidiscrete WENO coupling are certified
at sixth order, but pure temporal sixth-order convergence of the complete
production right-hand side with state-dependent LF/WENO switching is not.

## Deliberately uncertified scope

- Arbitrary nonuniform meshes, arbitrary custom maps, arbitrarily long runs,
  and a general nonuniform-grid stability theorem.
- General compressed Green data or singular product integration outside the
  narrow smooth vanishing-source fixture.
- Threshold-triggered adaptive-remesh global order. Fixed smooth transfer and
  safe rejection do not certify solution-dependent trigger times or maps.
- Nonsmooth global convergence, maximum-principle preservation, positivity,
  TVD, or SSP behavior. WENO-Z7 and polynomial remapping are unlimited.
- State-dependent peak/crossing gauge switches and global dynamic-rescaling
  order.
- A scalable symmetric high-order elliptic solver. The current one-sided
  closure is a direct-solve correctness backend.
- Pure sixth-order time convergence of the complete switching production RHS,
  as recorded explicitly above.

## Twelve-hour search status

The resumable `12h` optimizer continues on the separate
`research/sixth-order-global6` branch. It may inform a later recommended fixed
grid/CFL design only after its full budget, summary validation, and winning
neighborhood rerun complete. Until then the integrated feature keeps the
conservative `cfl=0.4` helper default and the fixed research evidence above;
no live rank-one checkpoint is a maintained optimum.

## References

1. Y. Shen and G. Zha, "Improved Seventh-Order WENO Scheme," *48th AIAA
   Aerospace Sciences Meeting*, AIAA 2010-1451, 2010,
   doi:10.2514/6.2010-1451.
2. P. J. Prince and J. R. Dormand, "High order embedded Runge--Kutta
   formulae," *Journal of Computational and Applied Mathematics*, 7(1),
   67--75, 1981.
