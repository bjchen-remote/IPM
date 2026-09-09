# Continuous peak: completed tiny integration and q512 instantaneous probe

The independent C1 amplitude rule preserved its actual peak to about
`1e-14` in the registered tiny runs, while the original 65 x 33 quadratic
rule retained a `1.63e-3` peak drift across its selector change. This is
evidence about the normalization and scale variables. At the same physical
time, both rules produced nearly identical physical fields. A stable
fourth-order full-field convergence rate was **not** established. No
production C rule, schema, native checkpoint or mainline trajectory changed.

## Actual integration, source identity and scope

The two grids were 49 x 25 and 65 x 33 on `[-4,4] x [0,4]`, with the
original Gaussian datum and maintained WENO5/Poisson operators. The
registered primary steps were `.0005/.00025/.000125`, all ending at
canonical time `.08`. The native branch called maintained
`ipm.evolve.stepSsprk54`. The C1 branch called maintained `field.velocity`
and `assembleRhs`, retained the exact transport-anchor length rate and
symmetric zero translation rate, and used the independent C1 amplitude
rate. It did not call the old quadratic selector to define its new RHS.

All five components
`[logCl,logComega,physicalTime,Xshift,canonicalTime]` were advanced with rho
through the same SSPRK54 A/b stage arithmetic. The canonical constant-rate
update was reconstructed as in the maintained step. There was no P
projection, feedback, scale reset, altered physical clock or new native CP.
Every C1 stage's peak, actual derivative, active cell, curvatures, validity,
rates and CFL was retained. Native diagnostics were accepted-state only,
because the maintained step internals were left unchanged.

Before integration, four steps in native mode gave bitwise agreement with
the maintained method for rho and all five scale components on both grids.
The new direct operator path also reproduced the original independently
assembled base RHS bitwise. More strongly, all 321 rho frames of the
65 x 33 native `.00025` run matched the previously verified original native
replay prefix bitwise, including the known switch.

The C1 MAT records explicitly identify
`researchEquationId=independent_hermite_C1_peak_amplitude_fixed_window_v1`.
Their `operatorSourceConfig` describes the maintained operator source;
it does not pretend the new RHS is a production quadratic-gauge equation.

## Primary finite-time results

| Grid | dt | Max relative C1 H drift | Max relative native quadratic P drift |
| --- | --- | --- | --- |
| 49 x 25 | .0005 | 2.589e-15 | 9.182e-15 |
| 49 x 25 | .00025 | 1.048e-14 | 5.819e-15 |
| 49 x 25 | .000125 | 6.990e-15 | 6.854e-15 |
| 65 x 33 | .0005 | 4.142e-15 | 1.631e-3 |
| 65 x 33 | .00025 | 4.789e-15 | 1.630e-3 |
| 65 x 33 | .000125 | 1.540e-14 | 1.630e-3 |

All C1 stages had one valid isolated peak; none failed the peak or native
CFL guard. Maximum instantaneous H derivative was `3.608e-16`. The C1
trajectory still showed a discontinuity in the *secondary quadratic
observer* on the 65 grid, but that observer did not set the C1 RHS.
Neither grid's Hermite maximum crossed a Hermite cell during this short
time interval. Thus this is a time-integrated old-selector test, not a
time-integrated Hermite-cell-crossing test. The latter still requires
separate validation; manufactured crossing tests alone do not replace it.

On the 65 grid, adjacent native rho infinity differences were
`8.025e-7` and `4.804e-7`, while the C1 differences were only
`2.442e-15` and `2.554e-15`. On the 49 grid both rules' differences were
around `1e-14` to `1e-13`. The fine levels therefore cannot measure a
fourth-order C1 convergence rate.

Primary output:
`hermite_peak_integration_tp995ad07e_6c1e_4745_8219_f53d760e2e6c/` under
`result/longtime/20260908_campaign_v1/acceleration_lab/`.

## Separately registered coarser calibration: no fourth-order claim

The parent then authorized the additional steps `.004/.002/.001` to the
same `.08` endpoint. Every run completed and respected native CFL; the
largest C1 CFL was `.18143` on the 65 grid and `.13603` on the 49 grid.
The original fine protocol and all its results were retained separately.

C1 peak drift at `.004/.002/.001` was
`1.419e-13/8.414e-15/2.071e-15` on 49 x 25 and
`9.708e-14/6.213e-15/2.589e-15` on 65 x 33. Although this observable's first
reduction is consistent with fourth-order behavior before the roundoff
floor, it does not establish fourth-order accuracy of the full field.

The C1 adjacent rho infinity differences were
`1.122e-12 -> 3.698e-13` on 49 x 25 and
`2.049e-13 -> 5.018e-14` on 65 x 33. These correspond to self-ratios of
about 1.60 and 2.03 in log2, not four. Comparison against the retained
`.000125` numerical reference also did not produce stable fourth order:
C1 L2/infinity orders were about `1.451/1.164`, then `1.160/1.290` on 49,
and `3.237/2.347`, then `1.854/1.332` on 65. Reference sensitivity is
reported with each error. The errors are small, but their observed orders
must not be relabeled as a fourth-order pass. The cause is not established
by these runs; this report does not assign it to WENO nonsmoothness,
Poisson arithmetic, or roundoff without further evidence.

Coarse output:
`hermite_peak_integration_tpb5a31df0_fed3_4cd0_854e_bdb5e39b6ae9/`.
Its field-only reference assessment is
`fine_reference_tp40a2e9da_7770_490f_a84d_9101d021ab19/` within that output.

## Common physical time and the spatial error scale

A separate observer evaluated all primary runs at fixed physical time
`.075` on a fixed 97 x 49 observation grid in `[-3,3] x [0,3]`. Cubic
temporal Hermite used two genuine endpoint RHS evaluations per trajectory,
with derivatives divided by the actual `t_tau`; all five scale variables
were interpolated consistently. The spatial observer was the already
verified linear C1 tensor Hermite with maintained Dx/Dy. Reconstruction was

```
rho_physical(x,y,t) = R(Cl*x+Xshift,Cl*y,tau)/Comega.
```

The maintained-gradient observer was `Cl/Comega` times interpolated
`R*Dx.'`; the spatial derivative of the same density interpolant was
reported separately. Neither is silently substituted for the other.

At 65 x 33, the C1-minus-native physical density infinity difference was
`1.464e-11/2.819e-12/9.172e-13` across the three primary time steps. At
49 x 25 it was about `3e-15`. By comparison, the 49-minus-65 spatial
density difference was about `1.573e-4` in infinity norm and `1.045e-5`
in RMS L2; the maintained-gradient infinity difference was `6.510e-4`.
Those are two-grid differences, not a continuum error bound.

Thus a large difference in gauge-coordinate rho cannot be read directly as
the same physical PDE error. The original selector's amplitude/clock error
largely cancels in consistent physical reconstruction in this example.
The new rule's demonstrated benefit is continuous normalization and more
stable rescaled variables; it is not a demonstrated `1e-7` improvement in
physical density accuracy.

The alternative piecewise-linear temporal observer was retained and
differed from the cubic observer by about `1.1e-8` to `1.3e-8` at the
coarsest primary step and `3.8e-10` to `7.6e-10` at the finest. This is a
temporal interpolation sensitivity measurement, not an error bound for the
cubic observer. All comparisons used the predeclared cubic observer and
were not selected after seeing the answers.

Physical-observer output is
`physical_comparison_tpcfcd731a_3639_456b_801c_5db1fb7943cc/` within the
primary integration output. It used 24 tiny bracket RHS evaluations and
two tiny operator initializations, but no additional PDE steps.

## Local tolerance audit before the large-box probe

The initial helper multiplied stationarity residual by the entire search
window width and derived peak-position tolerance from its distant endpoint.
This can falsely reject an unchanged local peak merely because the box is
extended. The old helper is preserved as
`ipm_accellab_hermite_peak_window_scaled_v0.m`.

Among 36 registered narrow-Gaussian/window examples, the old policy had
nine false rejections when the window endpoint became `1e6`. For example,
at center `.70123`, width `.027`, the stationary residual was only
`3.659e-16`, yet multiplying it by `1e6` failed the old criterion. The raw
maximum value was unchanged.

The revised validity metric uses the cell coordinate and the local
polynomial coefficient scale:

```
stationarity = abs(H_x)*cellWidth/localPolynomialScale,
curvature    = H_xx*cellWidth^2/localPolynomialScale.
```

Position grouping uses the candidate's local cell/coordinate magnitude and
ULPs, without a unit-size floor or distant endpoint. Active-value tolerance
scales with the field amplitude, preserving amplitude-unit covariance.
The numerical relative thresholds were not increased. Positivity,
negative resolved curvature, isolated uniqueness, search-boundary and
degenerate/multiple-peak guards remain explicit. The maximum and amplitude
rate formulas are unchanged; no physical residual or PPrime is corrected.

All far-endpoint examples then passed with bitwise unchanged peak values,
positions and forcing values on the same grid. Coordinate unit factors
`1e-6/1/1e6`, amplitude factors `1e-12/1/1e12`, and the complete manufactured
peak/degeneracy suite passed. All 2,246 accepted frames from the six prior
C1 primary trajectories had bitwise unchanged H values, positions and
selected cells under the new policy. This rechecks accepted observations;
it is not a replay of every earlier internal stage.

Audit output:
`hermite_peak_tolerances_tp19cb277a_af5f_43a6_a9f6_766150e5508b/`.

## Latest q512 instantaneous evaluation

In the separately authorized resource window, native step 3326 at
`tau=7.093201945463523`, physical time `1.9545271366857122`, was strictly
validated using the original ten computational threads. The grid was
1025 x 513 on `[-1e6,1e6] x [0,1e6]`. One restore supplied its fresh
original Poisson/RHS cache, and one independent transport assembly supplied
F0. There was no extra Poisson solve after restore, time step or checkpoint.

| Quantity | Measured value |
| --- | --- |
| Native quadratic P | .686092598153450 |
| Continuous Hermite P | .686087978203658 |
| Hermite peak position | .974129855139977 |
| Native cOmega | -.350430869832058 |
| Research cOmega_H | -.350115814304440 |
| HPrime on original RHS | -.000216155809965 |
| HPrime after independent `F0+cOmega_H*rho` | 2.776e-16 |
| Relative full RHS L2 change | .000882414304 |
| Hermite curvature | -25817.8912941 |
| Normalized local curvature | -.00304864556 |

The full positive-half and original tracking windows selected the same
peak. The old window policy happened also to pass this particular state;
the real probe is not relabeled as an old-policy failure. The synthetic
false rejections above remain the direct evidence for the tolerance fix.

The original Poisson warning `RCOND=3.303172e-18` and relative algebraic
residual `2.061596e-9` were preserved in `runtime.log` and the report. The
new algebraic cancellation does not improve or certify that underlying
Poisson solve. The process exited and its LU window was explicitly released
to the mesh experiment.

Output:
`hermite_peak_q512_tp5ba1ee09_6ad9_432d_a19b_7eff48d389df/`.
This minimal output retained wall Omega/forcing fields and x, not the full
two-dimensional RHS. It cannot support a latest-state two-dimensional
shape residual without an existing full-field cache or another authorized
RHS opportunity. No additional restore is planned for that purpose.
