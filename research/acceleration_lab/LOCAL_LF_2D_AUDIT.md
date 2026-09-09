# Native-rate two-dimensional local-alpha audit

This uses only the completed baseline step328 rho/physicalflow/operators and
original scale rates. No root solve, PDE step, Poisson build, wall-mode change
or main-trajectory replacement is permitted. The sole numerical variation is
face-local alpha with conditional-global split normalization.

The actual 1/16-line cost preflight estimates 15.43 seconds for the four
complete local kernels, within a preregistered 180-second ceiling. That is an
estimate, not an optimized production implementation. The first cost wrapper
wrongly indexed one-dimensional broadcast metrics; its failure is preserved.

The full line-mode research assembly must reproduce the cached native RHS
bitwise before the local mode runs. Both modes keep the same degree-5 ghost
extrapolation, mapping, free-stream correction, original conservative wall
and its y contribution, and amplitude-source grouping. Closed boundary mass
projection is not substituted for the original open boundary.

For each density/unit flux in each axis, verify the mapped telescoping
identity between weighted derivative and the two numerical boundary faces,
with maximum relative discrepancy 5e-11 (normalized by the maximum line flux
or one). Record all numerical boundary faces and their changes, including
the lower y ghost face. That face is not the physical y=0 plane, so a nonzero
value is not silently reset to zero. Report the unchanged physical normal
wall velocity and the advective/free-stream source contributions separately;
open advective-form evolution does not assert zero total quadrature mass.

Require exact reflection pairing of x nodes and report the actual even/odd
defects of cached rho, velocities/source and both RHS values. The maintained
time integrator does not project rho onto bitwise symmetry, so that property
is measured rather than assumed. Local RHS symmetry
must stay below 5 times the native absolute symmetry defect plus 1e-8 times
the baseline RHS infinity norm. A failure is retained; it does not change the
original data or symmetry choice.

Save full F/FX differences, whole-box/core/holdout norms, the physical
locations of largest changes, axiswise transport contributions and boundary
mass balances. Recompute actual PPrime, aPrime, smooth beta/gamma and complete
U_tau under the previously frozen theta=.5 calibration, with .3/.7 holdouts.
The same rho gives identical observation queries; assert exact common Gauss
nodes/weights before subtracting residual fields. Report both norms and
field differences, never just a scalar decrease. Compare their size with the
existing 7--14% same-observer spatial sensitivity.

Passing assembly, flux and symmetry contracts qualifies only this operator
audit. It does not establish physical-time equivalence, justify a new C rule,
approve a changed spatial scheme, or accept a profile accelerator. Every old
failed root and spatial verdict remains unchanged.

## Completed original-rate result

The unique output is
`result/longtime/20260908_campaign_v1/acceleration_lab/local_lf_2d_native_rate_tp603409fb_dbb1_4611_b648_b0ad46d11c38`.
`compact_summary.json` contains the scalar audit, while `fields.mat` retains
both complete RHS arrays, their derivatives, every axis/face intermediate,
and all three observers' residual fields. `source/` freezes the calculation
and preregistration; `source_manifest.json` also hashes the maintained
kernel, transport and assembly references. The MATLAB batch exited normally.

At baseline step328, absolute physical time `1.964631495951036`, the full
line-mode RHS reproduced the cached native RHS bitwise. It took 0.311 s;
the four local-alpha kernels and assembly took 3.162 s, below the registered
180 s ceiling. These are one-call timings, not production benchmarks. No
Poisson evaluation, rate solve or time step was performed. Rates remained
`[c_l,c_omega,c_r]=[2.087177884069464,-11.323896768347396,0]` in fresh units.

The maximum mapped flux-telescoping relative discrepancies were
`1.65e-14` (global) and `1.89e-14` (local). The even-RHS absolute defect
decreased from `2.60e-10` to `4.56e-11`, within the registered bound.
Both physical and transport normal velocities at the wall remained exactly
zero. The lower numerical y-face density flux changed by at most `1.46e-12`,
the upper face by `0.0243`, and the two x boundary faces by about `215.48`.
The latter are actual changes in this different open-boundary discretization;
they are retained, not declared negligible. The open-domain quadrature mass
RHS changed from `-2.2154216e7` to `-2.2197463e7`. Only the stated mapped flux
identity was qualified, not unchanged physical boundary flux or conserved
total density mass.

| Registered observer | Core residual-field L2 difference / global | Holdout difference / global | Core norm local / global | Holdout norm local / global |
| --- | ---: | ---: | ---: | ---: |
| Primary theta=.5 | 5.35098% | 1.31468% | 1.00094717 | 1.00011899 |
| Holdout theta=.3 | 5.14777% | 1.28710% | 1.00098316 | 1.00015998 |
| Holdout theta=.7 | 5.38237% | 1.31861% | 1.00075500 | 1.00001523 |

Every residual norm rose slightly. The primary core and holdout residual
norms in parent canonical time changed from `0.00743843524/0.01409529503`
to `0.00744548074/0.01409697221`. The same exact Gauss queries, weights and
fixed width calibrations were used for both operators. The real peak
derivative was `-0.08071730691` (global) and `-0.08141036105` (local) in fresh
units; it was never set to zero. Primary beta in parent time changed from
`-.6216582987` to `-.6216508475`, and gamma from `-.6555655378` to
`-.6555345020`.

The native-node density RHS relative L2 difference is `8.47e-6` in the core,
`4.43e-6` in the holdout and `0.00594` across the whole box. For its maintained
x derivative the corresponding values are `0.001582`, `0.000812` and
`0.03033`. The largest density RHS change (`1.08086`) lies at a far corner,
computational `(X,Y)=(-1471.065,1461.679)`. The largest omega RHS change
(`1.98552`) lies near the positive core, `(X,Y)=(.1499566,.0017880)`.
This distinguishes a small raw core forcing change from a larger change in
the remaining shape derivative after geometric cancellation.

The audit passes its assembly and flux contracts but shows no reduction of
the shape residual at the original rates. Its primary 5.35% core field
difference remains comparable to, and smaller than, the earlier 7--14%
same-observer spatial sensitivity. Neither operator is thereby identified
as the continuum truth. Large proposed geometric rates, physical-time
accuracy, spatial refinement of this new operator, and tail/velocity gates
for an evolved candidate remain untested. The earlier smooth-three-rate,
90%-width and secant rejections remain in force.

The reproducible entry point is
`ipm_accellab_audit_local_lf_cache(cacheFile,smoothReportFile,costReportFile,outputRoot)`;
the exact four arguments are retained in `registration.json`. The helper
creates a unique output directory and checks native parity before local work.
