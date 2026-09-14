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

The optional selector is `outer_wall_density_window_l2`. The current
long-time preset and existing release remain on the quadratic peak gauge.
For a **new experimental run from time zero**, construct the usual flat
long-time options, then set:

```matlab
opts.cOmegaGauge = 'outer_wall_density_window_l2';
opts.omegaGaugeWindowRadius = 0.5;
```

An old peak-gauge checkpoint cannot be resumed under a different frozen
gauge configuration. Do not treat the frozen-state RHS audit as such a
resume. The candidate is not yet qualified for an unattended long run.

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
