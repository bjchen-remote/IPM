# Experimental amplitude gauge at the wall point (2,0)

The corrected outer anchor is **(2,0)**. The implementation uses a smooth
wall window centered at `X=2`, rather than a single interpolated value, to
reduce sensitivity to an adaptive mesh crossing the anchor. Its half-width
in these experiments is `0.5`; its support is `1.5<X<2.5` at `Y=0`.

For the existing transport-anchor length gauge, write the rescaled density
equation as `R_tau = B + c_omega R`, where `B` is the maintained discrete RHS
with its amplitude term omitted. For fixed normalized nonnegative wall
quadrature weights `w_i`, the new selection is

```text
M(R) = sum_i w_i R(X_i,0)^2
c_omega = -sum_i w_i R(X_i,0) B(X_i,0) / M(R).
```

Thus `dM/dtau=0` exactly for the semidiscrete RHS on each fixed grid. The
window weights are recomputed from the current grid, and the initial window
RMS is stored as a frozen reference through remeshing and checkpoints. A
remesh transfer can cause a small jump in the observable; this gauge does not
cancel transfer error. The point-rule comparator is
`c_omega=-B(2,0)/R(2,0)`; it was not selected as the implementation because
linear interpolation at a single point can switch cells after remeshing.

If the wall gradient near `X=1` eventually behaves like
`|X-1|^(-p)`, a wall gradient `q`-moment crossing that point can diverge
when `p*q>=1`. No upper bound on `p` has been established. The selected
density window has positive distance from the putative local singularity;
its integral remains finite if the singularity stays localized and the
outer density stays finite. This is a conditional integrability argument,
not a proof of the asymptotic profile.

The screening used four signed original-`t0` H8 checkpoints, at
`tau=12.3003,13.0002,13.9001,14.2002`. The old quadratic wall-gradient
peak rule kept its peak near `0.685`, while `R(2,0)` fell from `0.01569693`
to `0.00868605` (44.664%). On the same orbit, the ratios
`R_X(X,0)/R(2,0)` at `X=1.4,1.6,2,2.4,2.8` changed by at most 0.39%.
Renormalizing those stored states by their value at `(2,0)` would therefore
make the sampled outer slopes nearly stationary and increase the inner peak
relative to that outer amplitude by about `1.807`. This is a postprocessed
comparison, **not** an integrated new-gauge orbit.

At the native `tau=13.9001` state, the old/point/window rates were
`-0.305432 / +0.00163487 / +0.00163871`. At `tau=14.2002` they were
`-0.30419842 / +0.00149257 / +0.00149598426`. The last state's
instantaneous sampled outer `R_X` log rates under the window selection were
between `-0.00102` and `-0.00128`, while the active inner peak-node log
rate was `+0.29806`. These are single-state RHS derivatives: they support
the proposed separation of an almost stationary outer wall field from a
still-evolving core, but do not prove a power-law limit or a stable
long-time new-gauge trajectory.

The selector is `outer_wall_density_window_l2`. R2.2 now makes it the
default for **new** long-time server runs; R2/R2.1 keep the old quadratic
peak gauge for their checkpoints. A flat options user can set explicitly:

```matlab
opts.cOmegaGauge = 'outer_wall_density_window_l2';
opts.omegaGaugeWindowRadius = 0.5;
```

An old peak-gauge checkpoint cannot be resumed under a different frozen
gauge configuration. Do not treat the frozen-state RHS audit as such a
resume. The candidate is not yet qualified for a scientifically certified
unattended long run.

R2.2's server configuration completed a bounded original-`t0` version-5
autonomous-mesh run: 10 steps to `tau=0.01` on the selected `321 x 161`
starting grid, with `2.22e-16` relative window drift, minimum window
condition `0.3345`, 38 support nodes, and zero normalized gauge residual.
No subsequent adaptive remesh occurred in that short interval. Separately,
three `2e-5` SSPRK54 steps from the signed `tau=14.2002` state were
evaluated **in memory** under the new gauge: relative window drift was
`2.22e-16`, point-value drift `2.05e-10`, and the inner peak increased
`1.7878e-5`. This branch is not an accepted new-`t0` checkpoint chain.
The mesh-transfer check from the bounded test preserved the frozen
reference exactly, but long-time adaptive remesh behavior under this gauge
is still untested.

Evidence:

- `evidence/outer_wall_gauge_frozen_20260914.json` (four signed checkpoints)
- `evidence/outer_wall_gauge_rhs_tau139_20260914.json` and
  `evidence/outer_wall_gauge_rhs_tau142_20260914.json` (native RHS audits)
- `evidence/outer_wall_density_gauge_test_20260914.json` (integrated
  candidate formula, tiny step, history, reference restore, mesh transfer,
  and latest signed-state RHS)

The bounded MATLAB check is `test_outer_wall_density_gauge(checkpointFile,
outputFile)`. It advances a tiny `129 x 65` test state, then only evaluates
the latest signed checkpoint's candidate RHS. It does not advance the
production checkpoint or restart a paused campaign.
`run_outer_wall_auto_smoke(outputFile)` reproduces the short original-`t0`
server run; `probe_outer_wall_late_branch(checkpointFile,outputFile)` makes
the three-step in-memory late probe. Both wrote their local outputs under
`result/verification/outer_wall_release_20260914/`.
