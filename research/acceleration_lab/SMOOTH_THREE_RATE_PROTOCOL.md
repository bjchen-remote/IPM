# Fixed-flow smooth three-rate experiment, 2026-09-09

This independent experiment uses the original global-LF, conservative-wall
baseline step328 arrays. It changes no density, scale, velocity, wall mode,
spatial discretization, native checkpoint or production C rule. The original
90%-width failed rate solve stays rejected.

The primary constraints are the actual complete-RHS values
`[PPrime/P, aPrime/w, (log Wx)Prime]`, with theta=.5,m=4 smooth Wx and the
previously frozen baseline width calibration. The residual observer uses
the same frozen horizontal and vertical unit factors. Theta=.3/.7 are
observation holdouts only; no rate solve or parameter search uses them.

Unknown increments are q=[deltaCL,deltaCW,deltaVPeak/w]. Actual arguments to
the maintained assembleRhs are cl=cl0+q1, cw=cw0+q2,
cr=cr0+q3*w-q1*a. At q=0 this reproduces the native arguments exactly. The
linear continuum generators supply one initial estimate only. All rate
probes recompute the complete nonlinear original transport at the same
cached physical velocity. No Poisson or PDE time step is permitted.

The original strict infinity root tolerance remains 5e-11. Up to 10 Newton
iterations use forward and backward one-sided secant matrices with relative
step 1e-5, in fixed forward-then-backward order. Matrix rcond must be at least
1e-8. Such columns are not automatically a realizable generalized Jacobian.
Each proposed direction is checked against an actual combined directional
one-sided response at h and h/2; relative model discrepancy must be <=.01.
Only then may a power-of-two backtracking step be accepted by actual norm
decrease (Armijo coefficient 1e-4; at most 10 halvings). Every attempted
direction and rejection is saved. A failed iteration stops the solve and
still receives the full fixed-final-rate audit. No tolerance is relaxed.

Native and final q receive both step sizes of the coordinatewise one-sided
audit. Complete x/y global-alpha selectors and signs are recorded. A
forward/backward mismatch or selector change is reported as nonsmoothness
evidence, not hidden by a centered Jacobian. Centered values are diagnostic
only and never establish differentiability or generate a search direction.

The fixed native, saved old failed q, and new terminal q are compared under
the identical smooth charts and quadrature. Necessary instantaneous gates:

- original 5e-11 actual constraint bound;
- every theta's core/holdout/full-local and wall Gauss L2 and sampled infinity
  residuals nonincreasing relative to native (relative allowance 1e-10);
- primary complete native-FD modulated residual L2/infinity nonincreasing on
  the whole native box, fixed core and holdout masks;
- reconstructed physical density-RHS discrepancy, after removing the
  continuum generator difference, <=1e-3 of baseline complete native-RHS
  infinity norm on the whole box and on the fixed physical outer tail
  `|x|>=2.5 or y>=2.5` within `[-3,3] x [0,3]`.

Raw density-RHS and differentiated-RHS norms are also saved; their decline
alone is not meaningful after changing coordinate rates. The physical
density, physical velocity and scale are identical by construction, so
their tail differences are zero but are explicitly **not** independent
fresh-flow or changed-state tail tests. No necessary-gate outcome can
qualify production: a separate half-plane IVP, matched physical-time,
boundary, resolution, cost and freshly solved velocity-tail checks remain
absent. Differences must be compared with the measured 7--14% smooth
residual spatial sensitivity, not advertised as sub-percent acceleration.

## Actual outcome

The registered solve completed with **rejection** after 46 original transport
assemblies, zero Poisson and zero PDE steps. Zero increments reproduced the
cached native complete RHS bitwise. The linear seed had constraint vector
`[.0135432,-.218519,-.110609]`; one checked, actual descent step reached
`[3.62182e-6,-6.79289e-4,-3.15537e-5]`. The next forward/backward models
disagreed with the combined directional h/h2 response by 39%/82%, so neither
direction was followed. The 5e-11 root gate remains failed.

At the retained terminal rates, parent units are
`cl=.689917616,cw=.273185131,cr=-.610748197`, with actual smooth
`gamma=-.0447358264`. A small remaining width rate does not certify the full
shape: primary core/holdout residuals are **20.3149/6.8084** times native.
Both holdout theta values and native-FD local norms also fail. The differential
physical-RHS coordinate-chain defect is .35284 in whole-box relative
infinity norm; its fixed physical-tail counterpart is .003684, failing both
registered 1e-3 bounds.

For an identical smooth-observer comparison, the previous 90%-width failed
rates have core/holdout residual ratios 4.5533/1.7832. They remain rejected;
changing the observer has not accepted that earlier result.

The native forward/backward matrix discrepancy is .00599, reproduced at
h/2 (.00587). The retained terminal wall global alpha grows from 2.79482 to
618.301 (221.23 times), with a different maximum location. Its one-sided
h-versus-h/2 matrix differences are .00396/.00268, and the combined-direction
test rejects the model. None of these matrices is called a certified smooth
Jacobian. The test establishes failure of the actual evaluated rates; it
does not prove that every possible root of the new constraints fails.

All full fields, theta holdouts, attempts, selectors and source hashes are
in `smooth_three_rates_tpe7307950_c92f_4619_9e9a_7a9e014cdc79` under the
campaign acceleration directory. Its original registration, old failures
and source checkpoints were not overwritten.
