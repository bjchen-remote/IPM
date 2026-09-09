# Late frozen-field rate response

The baseline step-328 diagnostic identifies a real line-global x LF fold and
strong amplification by the 90%-root-width functional. A smooth integral-width
control reduces the anomalous width response by about 642 times. The
conservative wall y term is small in this fixed field. The original failed
three-rate Newton result remains failed. No replacement rate, PDE, wall mode,
or production algorithm is selected here.

## Inputs and contract

The input is the already audited 895×386 cache
`result/longtime/20260908_campaign_v1/fresh_tau744_baseline_rhs_cache_v1/true_rhs_cache.mat`,
at absolute physical time 1.964631495951036. Its stored transport operators
contain no Poisson matrix or factor. The density, physical velocity, original
native scales, axes, derivatives, geometry and observation window are fixed.
This experiment does not revalidate a checkpoint; it verifies the paired
cache and reproduces its complete RHS bitwise using the original assembler.

The native wall is `conservative_flux`; the transport is `weno5_fd`, high
order, open. A separate research copy changes only the wall mode to
`advective_upwind`. It is an alternative discretization, not a native restart.
Its differences from the native RHS are confirmed to occur only on the wall.

The three-rate coordinates match the acceleration experiment exactly:
`q = [c_l, c_omega, (c_r+c_l*a)/w_x]`, with the original peak location and
90%-width held fixed. Thus the first direction changes `c_r` by `-a dc_l`.
A fourth direction tests raw isotropic dilation at fixed `c_r`. The paired
continuous-form generators use the maintained discrete derivatives:

* Centered dilation: `-(X-a).*rhoX-Y.*rhoY`.
* Amplitude: `rho`.
* Normalized translation: `-w_x*rhoX`.
* Raw dilation: `-X.*rhoX-Y.*rhoY`.

These are discrete evaluations of continuum coordinate generators, not exact
spatial-invariance identities. For example, their measured centered-dilation
width response is 1.00187406 rather than exactly one.

Relative perturbations are 1e-5 and 5e-6, scaled by `max(1,abs(q_i))`. Both
signs and both wall modes are saved. The C1 interface is the same
`ipm_accellab_continuous_inner_rates` used by the original experiment:
`[PPrime/P, aPrime/w_x, beta, gamma]`. All arrays retain the original native
units. The fixed comparison core is `|X-a|<=3w_x`, `0<=Y<=3w_y`.

The source contains 112 frozen files and three input SHA256 records. Two
research files passed `checkcode` with zero issues. The observational WENO
copy first passed four tiny exact-output tests. The actual run made 36
original `assembleRhs` calls, zero Poisson calls and zero PDE steps. All
observed four-kernel reassemblies matched their native RHS exactly. Source
inputs remained value-identical; all MATLAB processes used 10 threads.

## Conservative wall effect is small here

At the original rates `[2.0871778841, -11.3238967683, 0]`, changing only the
wall mode changes the wall density RHS by at most 2.4634472e-10. The native
minus alternative geometric response is
`[-2.33455e-8, -5.42571e-8, 1.17301e-7, 1.04232e-7]`.
The original beta is -19.9976607848 in native units. The corresponding wall
mode difference at the unaccepted rate vector is similarly small:
beta changes by 1.38983e-7 and the density RHS by at most 3.15301e-10.

The y term remains essential away from the wall: its contribution to native
gamma is 7.65644. The small *wall* effect does not justify deleting y transport
or changing the native wall mode. It simply excludes this wall contribution
as an explanation of the order-one three-constraint discrepancy in this field.

## The x dilation response is amplified

For the smaller perturbation, the original centered-dilation response and its
generator comparison are:

| Constraint | Actual WENO response | Generator response |
| --- | ---: | ---: |
| PPrime/P | -0.99967234 | -1.00000016 |
| aPrime/w_x | -0.00525584 | 0.00000291 |
| beta | 2.71952083 | 1.00187406 |
| gamma | 0.98653112 | 1.00000824 |

The beta defect is 1.71764677. At the larger perturbation it is 1.71765046.
The advective-wall calculation gives the same conclusion. The y contribution
to the beta derivative is at most about 1.6e-6 in these probes; the x term
provides essentially the entire amplified response.

The centered-dilation density response differs from its generator by only
0.693197% in the fixed-core weighted L2 norm and 0.453462% on the core wall.
The derivative and 90%-width functional amplify this modest field difference.
Amplitude behaves as expected, with a core relative field defect below
1.9e-11. The normalized-translation field defect is about 6.46e-6; its measured
phase response is close to one.

The four-row, three-column central-difference Jacobians change by relative
1.7876e-5 (native) and 8.7781e-6 (alternative) between the two perturbation
sizes. Their three-constraint reciprocal condition numbers are about 0.268.
This is diagnostic information, not a relaxation or a pass of the original
Newton acceptance threshold.

## A line-global LF fold is directly observed

On all 135 sampled core-wall faces, density and unit-field positive/negative
WENO weights are exactly `[0.1,0.6,0.3]` in the baseline, all 16 perturbed
inputs, and the unaccepted vector. The 64 perturbed weight matrices were
checked independently from saved MAT arrays. There is no local nonlinear
weight change on these sampled faces.

The x-line split normalization is large: approximately 3.8848e5 for the
original wall density and 1.8664e7 at the unaccepted rates. The wall LF speed
increases from 2.79482159 to 137.54553224, a factor of about 49.2. Its maximizer
moves from x=+0.1471009717 to x=-0.1472320182. These values are in the saved
native coordinate system.

At the original symmetric state, opposite x locations compete for the global
maximum. The two registered perturbation sizes reproduce these one-sided
alpha slopes:

| Direction | Forward alpha slope | Backward alpha slope |
| --- | ---: | ---: |
| Centered dilation | 29.75118 | -0.11957555 |
| Normalized translation | 0.06311436 | -0.06311424 |
| Raw dilation, fixed c_r | 14.81581512 | 14.81581512 |

For centered dilation and translation, the active maximum changes between
x=±0.1471009717. The raw-dilation probe does not switch. A centered finite
difference therefore averages a line-global LF maximum fold; agreement of two
centered differences alone does not establish a smooth Jacobian.

Re-reading the saved density RHS fields, with no new transport evaluation,
confirms that the geometric response itself retains this fold. Centered
dilation has forward/backward beta slopes 4.45111/0.98793; normalized
translation has slopes 0.00984965/-0.00478521. Both perturbation sizes reproduce
the gaps, in either wall mode.

A cross-direction check gives quantitative local attribution. Estimate beta's
LF sensitivity from the translation one-sided gap, then multiply it by the
independently measured centered-dilation alpha slope. This predicts a beta
generator defect of 1.71773517 versus the observed 1.71764677: relative error
5.15e-5. The other perturbation size gives 4.62e-5. The two directions estimate
the same beta-per-alpha coefficient, approximately 0.115939, to within a few
parts per million. This strongly identifies the global LF response in these
frozen samples. It is not a validation of a modified alpha, local normalization,
time integrator, continuum limit or alternate C rule.

## Smooth-width control separates the amplification

An additional saved-field-only control uses the acceleration experiment's
independently tested `ipm_accellab_smooth_inner_rates` with theta=0.5 and power
4. It integrates the connected positive part of the normalized C1 profile;
the integral derivative contains the genuine PPrime and moving-peak terms,
but does not divide by a threshold-root slope. This is a new observation,
not a replacement of the original experiment or its gate. The original rate
directions, step sizes, physical flow and old translation normalization are
kept fixed. It makes zero transport, Poisson or PDE evaluations. The complete
PPrime and phase one-sided response rows reproduce the original functional
bitwise.

At the smaller perturbation, the native centered-dilation beta response is
0.9973243013, versus its smooth-width generator value 1.0000016225. The defect
is -0.0026773213, about 642 times smaller in magnitude than the original
90%-width defect of +1.7176467729. The smooth forward/backward beta slopes are
0.9946173269/1.0000312757: the LF fold is still visible, with gap 0.0054139488.
The corresponding old gap is 3.4631795037, about 640 times larger. Both step
sizes and both wall modes give this conclusion.

The baseline smooth beta is -20.0207458807 in native units, compared with
-19.9976607848 for the old width. Close baseline width rates therefore do not
imply equally robust responses to changes of the coordinate rates. The data
support both an actual LF-induced RHS fold and a large measurement-functional
amplification; the order-one 90%-width response must not be presented as an
equally large error of the entire PDE profile.

## Records

All original observations are in
`result/verification/performance_lab_20260908/late_rate_response_step328_v1/`:
`report.json/.mat`, baseline fields, eight complete positive/negative field
pairs, unaccepted-rate fields, all LF/scale/weight arrays, and the original
registration. Supplemental saved-array reports are
`saved_weight_summary.json`, `all_probe_x_weight_audit.json`,
`saved_alpha_one_sided_summary.json`, `one_sided_geometry.json`, and
`cross_direction_lf_gap_attribution.json`. The additional smooth control is
`smooth_one_sided_geometry.json/.mat`, with a separately frozen
`smooth_postprocess_source`. The original acceleration failure
and its 5e-11 acceptance gate have not been modified.
