# Smooth connected integral-width observer

This is a new auxiliary measurement of already saved true RHS fields. The
original 90%-crossing-width measurements and every failed accelerator/gauge
decision remain unchanged. Primary theta=.5, power m=4, and the theta=.3/.7
holdouts are registered together; no threshold is selected after seeing data.

For the continuous C1 wall maximum P at a, define

```
s = [(H[Omega]/P - theta)/(1-theta)]_+
Wx = integral_component s^4 dX
Wx' = integral_component 4 s^3/(1-theta)
                * (H[FX]/P - H[Omega]*P'/P^2) dX.
```

Only the above-threshold connected component containing the peak is used.
The endpoint integrand vanishes, so no crossing-speed or root-slope division
enters Wx'. The true P' is retained. The vertical width uses the connected
component starting at the wall on the same continuous peak trace. Its forcing
is H2[FX](a,Y)+a'*H2[Omega]_x(a,Y), and its y jet includes the corresponding
mixed derivative. The peak speed still uses the continuous stationarity
identity, with native-cell interface/degenerate-peak guards. Flat threshold
components and topology-changing threshold contacts remain explicit failures.

Gauss8, split at native cells and the two threshold endpoints, integrates the
degree-12 value and directional derivative polynomials to roundoff. This is
exactness for the chosen piecewise-cubic observer, not a PDE convergence claim.
`ipm_accellab_smooth_inner_rates` returns raw W,W',beta/gamma and event data.

For comparing complete inner residuals, each theta receives a single constant
horizontal/vertical unit multiplier from baseline step328: old 90% width
divided by raw W. Those constants are frozen across the 42/56 cases and any
future time; they do not affect beta. They make the baseline observation
rectangle coincide with the original one and prevent a smaller window from
creating an apparent improvement. Raw uncalibrated physical widths and their
spatial differences are saved separately.

## Actual checks and results

The bicubic translated/scaled MMS retains nonzero amplitude rate .11.
Analytic errors across all thresholds were at most 1.03e-15; coordinate units,
translation and amplitude covariance errors were at most 6.03e-15. Central
parameter-time differences at h=.02/.01/.005 decrease by factors near four.
Independent Gaussian data on 49/65 grids have measured errors; the maximum
fine-grid error among the registered quantities is 3.591e-4, below the
registered 1e-3 bound. This is not a claim of a universal spatial order.
The first test window omitted a theta=.3 component endpoint and correctly
failed; it remains saved. The corrected common [-1,1.4] analytic window
encloses every registered component.

The three actual caches at physical t=1.964631495951036 give the primary
parent-time horizontal rates:

| Case | beta | gamma | core residual L2 | holdout residual L2 |
|---|---:|---:|---:|---:|
| baseline step328 | -.621658299 | -.655565538 | .007438435 | .014095295 |
| target42 step406 | -.621636625 | -.655484686 | .007478989 | .014177303 |
| target56 step494 | -.621632706 | -.655496869 | .007507898 | .014176156 |

Every true cached-RHS directional check uses Omega +/- h*FX and recomputes
peak and raw widths. At h=1e-5/5e-6/2.5e-6 the errors decrease approximately
quadratically; the maximum final relative error is 1.80e-8, with stable peak
cells. These are functional directional checks, not extra PDE steps.

Although beta is much less sensitive than the old crossing-width beta, the
full residual fields still differ appreciably. Union-native-cell Gauss
integration on the fixed charts gives primary core/holdout relative field
differences 7.97%/7.51% for baseline-to-42, 9.71%/14.02% for 42-to-56, and
12.39%/11.45% for baseline-to-56. Both holdout thresholds give similar values
and are all retained. The nearly unchanged residual norms do not establish
convergence of the residual fields. These differences are still far larger
than the previous 0.3% candidate benefit, and are not monotone between the
successive locally refined grids.

The original observer gives 32.53%/22.40% baseline-to-42 and 16.13%/16.22%
42-to-56 differences. Both are valid recorded observations of different
functionals; this new measurement does not replace them or quantify the
unknown continuum error. The target56 legacy/auxiliary spatial verdicts
remain false. No original PDE, C rule, wall mode, source checkpoint, or
acceptance gate changed.

Evidence under the campaign acceleration directory:

- `smooth_inner_tiny_tp907fdef5_4c6e_473e_ac37_5ef41a0c1d01`
- `smooth_cache_observer_tpc387fcb7_e2ef_4b66_8356_0bf5971dcb56`
- `smooth_cache_field_comparison_tpa0360d92_e1de_45f6_a48d_81a4bd434e40`

Entry points are `ipm_accellab_test_smooth_inner`,
`ipm_accellab_smooth_cache_reports`, and
`ipm_accellab_compare_smooth_cache_fields`. All executed work used saved arrays
or manufactured data, with zero Poisson builds, zero PDE steps and zero
amplitude projections.
