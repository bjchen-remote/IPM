# Reuse during frozen-field rate probes

The maintained `weno5_fd` path cannot precompute density-only positive and
negative derivatives and then combine arbitrary trial velocities. In
`+ipm/+field/transport.m`, all four flux derivatives take both the density
(or unit field) and the actual transport velocity. In
`+ipm/+field/weno5FluxDerivative.m`, the velocity determines the per-line
Lax–Friedrichs speed, extended physical flux, positive/negative split fluxes,
normalization scales and nonlinear WENO weights. Fixed density and grid do
not make those weights independent of trial `cX`, `cY` or `cR`. The unit-field
free-stream correction has the same dependence. Holding its weights or LF
speed fixed would define a different operator.

There are narrower reuse opportunities for the current frozen-field Newton
probes. They follow directly from the existing code; no candidate has been
implemented or numerically qualified by this audit.

* With identical `rho`, grid, physical velocity, boundary options and
  `cX/cY/cR`, a change of `cOmega` only changes the final
  `rhs = baseRhs + cOmega*rho` in `assembleRhs`. Reusing the identical
  `baseRhs` and preserving the zero branch and multiplication/addition order
  is an exact-expression candidate. The three-rate finite-difference
  Jacobian can reuse that base for both amplitude perturbations.
* A `cR`-only change leaves the y transport velocity unchanged. Its density
  and unit-field y flux derivatives can be reused. A `cX`-only or `cY`-only
  change in the four-rate experiment similarly leaves the other axis
  unchanged. Recombination must preserve the original arithmetic grouping,
  wall replacement and any boundary correction.
* For any such reuse, the final combined RHS must still be differentiated
  with the maintained `Dx` before applying the geometric functional. The
  identity `(B+c*rho)*Dx' = B*Dx'+c*(rho*Dx')` holds in real arithmetic but is
  not a bitwise replacement for the existing calculation. An analytic
  amplitude Jacobian column can be a separate solver choice; it does not
  establish exact equality with the existing finite-difference probe.

The current experiment in
`research/acceleration_lab/ipm_accellab_probe_four_gauge.m` already reuses one
physical flow for all rate probes. This audit does not add a Poisson solve,
modify that experiment, change production code, or claim a measured speedup.

## Wall-only constraint follow-up

For `weno5_fd` with `advective_upwind`, the original transport replaces the
bottom row by `-(drhoX-rho.*divergenceX)`. All LF speeds, split scales and WENO
weights in this x derivative are per x-line; the regularizer depends on Nx,
not Ny. Thus the three amplitude/peak-position/wall-width constraints have a
mathematically local wall expression. Their Hermite functionals use the wall
forcing and its maintained x derivatives. The vertical-width rate must still
be measured from a genuine full forcing after solving the three constraints.

The direct single-line candidate `ipm_perflab_wall_rhs` was compared with
original full 2D `assembleRhs` in 54 no-LU synthetic cases: 33×17, 49×25 and
65×33; stretch parameters 0, 4 and 10; three rate vectors; open and closed
boundaries, all with the advective wall. These are expression tests, not PDE
or mesh-admissibility tests. Both files passed `checkcode` with zero issues.
All 54 transport velocities matched, but no complete wall RHS matched
bitwise. Differences occurred only in the first/last three x nodes, with
maximum absolute difference 2.4726887204451486e-12. Full forcing/forcing-jet
checks consequently failed as well. Interior RHS nodes were all exact.
Records and the tested negative source are preserved in
`result/verification/performance_lab_20260908/wall_rhs_tiny_v1/`.

The observed boundary pattern is consistent with the changed matrix
multiplication shape in the degree-five ghost extension: the native call
multiplies Ny×6 by 6×3, whereas the direct wall call multiplies 1×6 by 6×3.
One cannot assume those reduction paths are bitwise identical. A candidate
retaining those original full-size boundary products while reconstructing
only the wall line has been prepared in `ipm_perflab_wall_rhs_ghosts`,
`ipm_perflab_weno5_wall_ghosts`, and `ipm_perflab_test_wall_rhs_ghosts`.
**This correction has not been executed or qualified.** It has no production
caller. The original negative result was not overwritten.

The actual late baseline step-328 cache was then confirmed by the acceleration
experiment to use `conservative_flux`, `weno5_fd`, and `open`. The advective
wall simplification therefore does not apply. The original wall RHS retains
its y flux derivative: zero nodal vertical velocity does not imply a zero
nonlinear discrete derivative. The y LF speed and split normalization depend
on each complete y-column. The actual rate Newton experiment consequently
continues to use complete original assembly.

A possible separate conservative/open trace candidate would retain all those
y-column reductions and reconstruct only the two y faces used at the wall,
plus the x wall derivative. It must retain all four density/unit terms and
the original arithmetic grouping. This is only a code-structure opportunity;
it has not been implemented or tested. Closed conservative transport would
also require its whole-domain projection and is outside this proposal.
