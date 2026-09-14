# Outer-cutoff convergence diagnostic for the new C run

The candidate singular point is fixed at `(X,Y)=(1,0)` in the **rescaled**
coordinates. This study therefore measures `R_{X tau}` on the current
rescaled grid, rather than a physical `rho_{xt}` evaluated at a moving
physical location. On a fixed grid, the solver supplies its accepted-state
semidiscrete `R_tau`; multiplying it by the same `Dx'` used for `R_X` gives
the mixed derivative without a noisy finite difference between checkpoints.

Use a nonnegative smooth cutoff `chi_delta` that vanishes for
`r = sqrt((X-1)^2+Y^2) <= delta`, rises as
`sin^2(pi*(r-delta)/(2*delta))` on `delta<r<2*delta`, and equals one beyond
`2*delta`. Because the center lies on `Y=0`, the removed region is a
half-disk. A fixed outer window is one on `0<=X<=3, 0<=Y<=1`, then tapers
to zero at `X=4` and `Y=2`; it avoids the distant computational boundary.
The tested radii are `0.05, 0.1, 0.2, 0.4`, and both the two-dimensional
bulk and one-dimensional wall norms use `p=1,2,4,infinity`.

For finite `p`, the primary observable is

```text
E_{p,delta}(tau) = ||chi_delta R_{X tau}(tau)||_{L^p},
G_{p,delta}(tau) = ||chi_delta R_X(tau)||_{L^p},
e_{p,delta}(tau) = E_{p,delta}(tau)/G_{p,delta}(tau).
```

The same quantities are evaluated on the wall. Positive trapezoid cell
weights define all discrete norms, so an `L^p` value is never a signed
high-order quadrature cancellation. `p=infinity` is the maximum of the
cutoff field on the sampled grid. The radius sensitivity is essential:
a decay visible only after a very large exclusion zone would not establish
convergence immediately outside the singularity.

If `R_X(tau)` has a fixed-coordinate continuum interpretation, then

```text
||chi_delta [R_X(tau_2)-R_X(tau_1)]||_p
    <= integral_{tau_1}^{tau_2} E_{p,delta}(s) ds.
```

Thus `E -> 0` alone is not enough; a finite tail integral would provide
a Cauchy bound. On the actual adaptive mesh, remesh transfers can add
jumps not represented by the continuous RHS integral. A convergence
claim also needs those jumps controlled, direct paired-profile differences
on a common fixed observation grid, and spatial/time-step/box refinements.
The report also records consecutive Cauchy increments after linearly
interpolating both native `R_X` fields to a common `401×201` observation
mesh on `[0,4]×[0,2]`. The last pair is checked again on `801×401` nodes.
Those increments include intervening remesh effects; observing two meshes
agree does not establish native spatial convergence. Finite-window slopes,
time integrals and increment sums are diagnostics only.

For a possible later *alternative gauge*, write `R_tau=B+c_omega R` and
let `Q_p=||chi_delta R_X||_p^p`. For `p=2` or `4`, the exact fixed-grid
condition `d_tau Q_p=0` selects

```text
c_omega^(p,delta)
  = - integral chi_delta^p |R_X|^(p-2) R_X B_X
      / integral chi_delta^p |R_X|^p.
```

The analysis computes this rate for bulk and wall cutoffs but **does not
apply it to the trajectory**. It is compared with the released `(2,0)`
outer-density-window C. `p=1` and `p=infinity` remain diagnostics because
their derivative-based gauge conditions are nonsmooth at zeros or maxima.
No cutoff radius or `p` is promoted to a new C rule until it stays
well-conditioned and its outer derivative/candidate rate is robust over
the long new-C trajectory, grid changes, and larger boxes.

Implementation:

- `outer_cutoff_norms(state)` is a read-only per-state calculation.
- `analyze_outer_cutoff_checkpoints(runDirectory,outputFile,maxFiles)`
  restores the latest accepted checkpoints under the matching release,
  evaluates native `R_{X tau}` and common-grid Cauchy increments, and writes
  a non-overwriting JSON report.
- `outer_cutoff_pair_norms(first,second)` defines the common-grid increment.
- `test_outer_cutoff_norms()` checks exact stationary and uniform-growth
  identities and a paired-grid identity on a small analytic synthetic field.

The live R2.2 `t=0` run is separate from these files and uses the immutable
release under `../releases/ipm_long_time_server_r2_2_20260914`.

## Executed preliminary screen

The signed old quadratic-peak-C checkpoints at `tau=12.3003, 13.0002,
13.9001, 14.2002` were each divided by their own `(2,0)` wall-density
window RMS before comparing the shape of `R_X`. This is **not** a new-C
trajectory. For `delta=0.2`, bulk `L²` relative increments on the three
unequal time intervals were `0.00523, 0.00676, 0.00157`; the corresponding
`delta=0.1` numbers were `0.00552, 0.00703, 0.00167`. The last-pair bulk
`L²` observation-grid change was `1.08%` for `delta=0.2` and `2.99%` for
`delta=0.1`; higher `p`, especially infinity, was more grid sensitive.
These figures motivate `delta=0.2`, bulk `L²` as the primary shape screen,
with `delta=0.1`, wall `L²`, and bulk `L⁴/L∞` retained as sensitivity tests.
They do not establish a tail rate or convergence.

The independent original-`t=0` R2.2 new-C run reached accepted checkpoint
`tau=3.00041744` on `641×161` after 21 automatic remeshes. At `delta=0.2`,
bulk `L²` `E` was `0.12154, 0.10165, 0.09278, 0.09822` on the checkpoints
near `tau=1.5, 2, 2.5, 3`; the three relative profile increments were
`0.07249, 0.06614, 0.06703`. The final increment changed by only
`3.89e-5` relatively between the two observation meshes. Bulk `L⁴` and
`L∞` show stronger growth near the excluded core; at `delta=0.1` bulk
`L∞` `E` rose from `0.19285` to `0.378995` over this window. The wall
`L²` `E` fell from `0.08187` to `0.05451` for `delta=0.1`. Thus early wall
decay alone cannot be used to claim convergence of the two-dimensional
outer profile. At the same checkpoint the H8 far-boundary velocity/source
ratios were `0.9872/0.1599`, both above the `0.01` warning threshold, so a
larger-box comparison is necessary for quantitative Profile claims.

Evidence is under `result/verification/outer_c_long_20260914/`:
`cutoff_old_frozen_v2.json`, `cutoff_tau3_v1.json`, `status_tau3.log`, and
the untouched live run's `run/checkpoint_manifest.jsonl`. Later reports
must use new filenames; these early results are not overwritten.
