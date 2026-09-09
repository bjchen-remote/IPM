# Continuous geometric gauges: actual transport constraints

This independent research experiment changes no production rule or native
checkpoint. It initializes separate half-plane Gaussian cases on 49 x 25 and
65 x 33 nodes. The box and maintained Green boundary construction are the
same physical problem for the baseline and rate probes of each case. It does
not inherit the old double-odd trajectory: a translating gauge is incompatible
with that symmetry restriction.

Let Omega = rho Dx', and use the complete differentiable C1 peak and connected
90% widths from `ipm_accellab_continuous_inner_rates`. At one fixed rho and
fixed physical velocity, define the linear-in-forcing functional

```
L(F) = [P'/P, a'/wx, wx'/wx, wy'/wy]'.
```

Every component uses the actual F Dx'; in particular P' is measured, never
assigned zero. The maintained transport assembly evaluates
`F(cx,cy,cw,cr)` with velocities `u1+cx*X+cr` and `u2+cy*Y`, and source
`cw*rho`. At fixed data the WENO splitting is nonlinear in these rates.
Consequently the continuum density generators `[-X*rhoX,-Y*rhoY,rho,-wx*rhoX]`
and their 4 x 4 functional matrix provide only the initial rate estimate.
They do not define the accepted RHS.

The actual constraints are solved by damped Newton using centered differences
of the maintained transport assembly, relative step 1e-5, at most 12
iterations and 10 halving backtracks. The registered residual infinity bound
is 5e-11, and the Jacobian reciprocal condition bound is 1e-8. The final
Jacobian is recomputed at h/2 and must agree to relative Frobenius 1e-5.
Each rate probe reuses the same Poisson-derived physical velocity. There is
one initial tiny physical flow per grid and no PDE step.

The isotropic alternative solves the first three constraints with cx=cy,
retaining its actual vertical width rate. Its smaller matrix uses the same
generator construction and actual WENO Newton procedure. This would preserve
kappa=1 and the maintained cached Poisson LU. In the high-order direct path,
`ipm.field.poisson` uses the cached operator at kappa=1; a different aspect
builds a fresh LU decomposition on each call. Thus the four-constraint
anisotropic option currently has a major prospective cost beyond its small
rate solve.

## Executed instantaneous results

| Quantity | 49 x 25 | 65 x 33 |
|---|---:|---:|
| Linear seed's actual WENO constraint infinity norm | .1332632 | .0874879 |
| Four-constraint final infinity norm | 7.86e-14 | 2.50e-13 |
| Four-constraint final Jacobian rcond | .07810 | .07937 |
| Jacobian h versus h/2 relative difference | 4.96e-9 | 8.86e-9 |
| Newton iterations | 4 | 3 |
| Four-constraint transport assemblies, including validation | 54 | 45 |
| Three-constraint transport assemblies, including validation | 42 | 35 |
| Three-constraint actual remaining gamma | +.6496151 | +.6522324 |

For four constraints the final rates `(cx,cy,cw,cr)` are
`(.1954965,-.4536249,.0805319,.5148362)` and
`(.1879327,-.4642032,.0748909,.5248105)` respectively. The three-constraint
solutions have the same cx,cw,cr but cy=cx. This initial Gaussian has a large
remaining vertical expansion; a small late-state beta-gamma observed elsewhere
cannot be substituted for this measurement.

The four-constraint actual P' values are -8.33e-15 and +1.04e-14. Signed linear
field paths `rho +/- h F` recompute all geometry without projection. At
h=.001,.0005,.00025 the maximum derivative errors are
`2.07e-8,5.17e-9,1.29e-9` and `2.50e-8,6.25e-9,1.56e-9`; the registered
geometric signatures remain fixed. This is a directional derivative check,
not PDE time integration or a claim of finite-step gauge conservation.

The whole-box relative infinity defect in reconstructing physical density
RHS by the discrete coordinate chain rule is .00855/.00406 (four constraints),
with relative L2 .00548/.00197. The nonzero discrepancy is preserved. Exact
continuous coordinate covariance does not make a WENO discretization commute
exactly with coordinate transformations. These results pass the instantaneous
constraint protocol only. Physical-time comparison, evolving scale/aspect,
tail velocity, boundary and long-time shape tests are outstanding.

## Reproducibility and preserved failures

Run `ipm_accellab_probe_four_gauge(frameFile,outputRoot)` with the verified
`quadratic_switch_native65_v2/switch_frames.mat` source configuration. The
first attempt rejected both grids before any flow build because the old
`transport_anchor` label is illegal with half-plane symmetry. The corrected
independent anisotropic-fixed setup uses the legal unused `local_strain`
label. The original failed directory is
`four_gauge_rhs_tpbe51780d_a6a5_4360_b680_56f8936a9f0f`.

The actual computations completed in
`four_gauge_rhs_tpb03ac447_9727_480f_8ebb_e2012bbfa0d0`, saving each full MAT
report and three-constraint F3. JSON export then failed on a function handle;
the original process failure and its top-level failed report remain. Because
the old field save followed JSON export, four-constraint full RHS arrays were
not saved. No computation was repeated merely to repair the presentation.

`ipm_accellab_recover_four_gauge` read those exact MAT reports and retained both
cases in `four_gauge_recovered_tpc0c3a9aa_ffa4_4cd7_b86b_dbcc3f514030`, with
literal function-handle descriptions in JSON and original callable data in
the copied MAT files. Current source saves fields before JSON and uses
`ipm_accellab_json_view`; these are export changes after the numerical run.
Current recovery derives each pass flag from the saved result rather than
assuming success. All paths above are under the campaign acceleration output
directory; no original report was overwritten.
