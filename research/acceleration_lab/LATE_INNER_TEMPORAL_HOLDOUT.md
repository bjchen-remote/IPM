# Late inner-shape temporal holdout: fixed rank one rejected

The registered late rank-one predictor failed its finite-prediction screen. It beats persistence throughout this holdout, but the fixed linear trend becomes more accurate at later horizons. There was no PDE injection, accepted predicted time, extrapolated checkpoint, limit estimate, C-rule change, or large LU.

## What was already known and what is different here

The old `continuous_increment_dmd_tau509_709_v1` was an exploratory rank1/2/3 screen on an already-inspected series. Its subsequent `prospective_continuous_forecast_v1` genuinely froze rank1 before eight future observations at parent-equivalent tau7.4432–8.8432. Rank1 L2 error grew0.135%–0.756%, while remaining1.39–2.95 times smaller than fixed-linear error. Distance to its speculative candidate limit fell to0.104% then rose to0.456%. That successful finite forecast did not validate the candidate limit. The previous memory4/damping1/16 fresh-RHS secant remains rejected, including its tail-velocity gate failure. This new experiment does not reopen that result.

The present cohort is later and entirely on the one cropped-case native lineage. Eight endpoints at tau9.2432–10.6432 (steps2443,2704,3048,3369,3668,4060,4426,4768) are training only. The next eight registered endpoints at tau10.8432–12.2432 (steps5214,5631,6020,6525,6998,7545,8057,8650) are withheld from fitting. The premature step7510 checkpoint is excluded by the preselected equal-time endpoint rule; step7545 completes its actual intended endpoint. Every spacing is0.2 in parent-equivalent canonical units, verified against the strict native clocks and direct fresh-parent lineage.

These later CPs already existed and some had earlier scientific plots, so this is a chronological holdout, **not a newly prospective data collection**. The protocol, CP paths, source hashes, fixed rank and all thresholds were written before the new future arrays were loaded. Training ran in a separate process that loaded only training CPs. The serialized model SHA `98fdc9f4396c49ac92673b563686ce17798740a9b01476203fe063c115677aa7` was recorded before a separate evaluation process opened future CPs. Its SHA and all registered inputs remained unchanged afterward. No rank, window, clock, baseline, regularizer or coefficients were changed on the future errors.

## Fixed representation and model

Each real frame uses its own maintained paired derivative Omega=rho*Dx' and the previously tested continuous C1 wall peak and connected0.9 widths. Wall U=H2[Omega](a+wx*xi,0)/P is observed on xi[-2,2]; vertical U=H2[Omega](a,wy*eta)/P on eta[0,3]. This uses the **actual future geometry as an observation**; it is not a full physical-field prediction using predicted geometry. Actual amplitude, phase and widths are therefore reported separately.

Each trace is fit independently with rank-one increment DMD, no ridge, no eigenvalue clipping, fixed801/601 trapezoid training samples, and the predeclared relative signal floor1e-12. For X=[delta1,...,delta6], Y=[delta2,...,delta7], X=USV', B=Y*v1/sigma1, lambda=u1'*B. Future increments start at B*(u1'*delta7) and multiply their coordinate by lambda. The resulting prediction is a fixed linear combination of the eight training fields. Persistence is the last training value; linear trend is ordinary least squares using all eight training times. Four separate scalar models use log physicalP, physicala, log physicalwx and log physicalwy, with the same fixed rank and baselines.

Hermite is primary. Bilinear observations/models use exactly the same measured C1 peak and widths, with every result retained. Exact represented L2 errors use Gauss4 on the union of training and heldout native pullback cells: the squared piecewise-cubic prediction error has degree at most6. This is quadrature exactness for the represented functions, not a PDE accuracy statement. Infinity errors are uniform1601/1201 samples; the801/601 differences are recorded. The reduced1D traces are independently compared with the original full tensor-Hermite observer.

The predeclared screen requires positive stable training modes, rank1 beating both baselines in L2 and sampledInf for both primary traces at every future time, and all L2<=0.5%, Inf<=1%. It is only a finite-prediction research screen. It is false; no gates were relaxed.

## Complete primary L2 errors in percent

| Future tau | Wall persistence | Wall linear | Wall rank1 | Vertical persistence | Vertical linear | Vertical rank1 |
|---|---:|---:|---:|---:|---:|---:|
| 10.8432 | 0.16899 | 0.06288 | 0.04758 | 0.03988 | 0.03303 | 0.02616 |
| 11.0432 | 0.36772 | 0.08913 | 0.04118 | 0.08383 | 0.06107 | 0.06180 |
| 11.2432 | 0.60758 | 0.08500 | 0.08479 | 0.12478 | 0.08650 | 0.09807 |
| 11.4432 | 0.87537 | 0.12121 | 0.22627 | 0.14040 | 0.08705 | 0.11104 |
| 11.6432 | 1.04661 | 0.15178 | 0.27911 | 0.17121 | 0.10263 | 0.14017 |
| 11.8432 | 1.28636 | 0.17771 | 0.41297 | 0.19278 | 0.10923 | 0.16076 |
| 12.0432 | 1.51772 | 0.23154 | 0.55833 | 0.19756 | 0.10000 | 0.16510 |
| 12.2432 | 1.69548 | 0.31150 | 0.66708 | 0.19015 | 0.08056 | 0.15760 |

Final sampledInf errors (persistence / linear / rank1) are2.2323% /0.4695% /0.7801% on the wall and0.3178% /0.1149% /0.2365% vertically. The final rank1 L2 errors are2.14 and1.96 times the respective linear errors. The last two wall-rank1 L2 errors exceed the fixed0.5% screen.

The training rank1 modes have lambda0.880167 (wall) and0.602485 (vertical), and capture94.17%/97.61% of the leading training increment singular energy. Yet their increment-evolution fit residuals are49.69%/82.06%. Low spatial rank is not equivalent to a well-modeled decay law. The bilinear vertical mode changes to0.693746, so the small vertical increment dynamics are observer-sensitive. All coefficients and singular values are retained; no alternate rank was tested after seeing these results.

## Geometry is more predictable, but phase is not accurate at the inner scale

At the final actual tau12.2432, physicalP=153.3271203, a=0.1251960806, wx=3.461686896e-5 and wy=8.014657820e-6. Rank1 amplitude error is0.09952% compared with linear3.9273%; width errors are0.62657%/0.41318% compared with linear6.9422%/8.0989%. The peak-position absolute rank1 error is5.37355e-5, compared with linear0.00243373, but this still equals1.5523 actual wall widths. Reporting only the small relative macroscopic position error would hide an important inner-profile mismatch. Every scalar prediction/error at every horizon is saved.

The maximum same-frame Hermite–bilinear sampledInf discrepancy is0.000290376498; the maximum Inf-sampling refinement change across all observers/models/horizons is4.74894983e-05. The reduced trace vs full tensor-Hermite discrepancy is at most1.67421632e-13. These observer checks are appreciably smaller than the late wall prediction error, but they do not certify native spatial convergence. Future frames cross actual remesh epochs7–10, with native core counts as low as22.815/32.243 at step8057. The original box/mesh sensitivity and the cropped-IVP distinction remain. This is not the original-t0 A/B comparison and not an estimate of the universal infinite-time profile.

## Consequence for the next acceleration experiment

The geometry variables contain useful finite-horizon regularity, whereas this fixed single-decay representation of the remaining inner shape does not remain best on later real data. A controller can only use such geometric forecasts after a separate safety study; this report does not authorize changing its decisions. A fixed-point or secant acceleration of the raw shrinking profile is still unsupported. The next defensible research experiment should preserve the actual geometry/phase error budget and use a newly prospective later validation interval for any new shape-clock or multi-mode proposal. Retuning ranks on these eight failed holdouts and calling them new validation would not supply that evidence.

The full artifacts live in `late_inner_temporal_holdout_2eff8944330844d182e8f929bf6ebf42`: registration and hashes, train/frozen_model, complete actual native metadata, all model/observer error CSVs, saved curves, algebra check, and the static figure. No predicted field is labeled as a native state.
