# Fixed full-quadrant extension of the frozen homogeneous model

This is preparation only. No MATLAB, checkpoint-array observation, model run, Poisson or PDE was executed while preparing this protocol. The actual prospective tau8 wall-sector test has already completed and is preserved. The new full-angle experiment is an **extension selected after development and the tau8 result**, not another prospective test. It does not approve an evolved exterior-source closure or a boundary/C modification.

## Fixed source and observation contract

Use exactly original-H64 case `20260909T060628886_dynamic_isotropic_321x161_tp901b551b_f5f7_4cac_bdea_0d56d90b4559`, native steps2213/4291/5874, with the immutable runtime from `from_zero_H64_tau12_auto_initial_v1/source`. Their saved physical times are1.8142778041740544,1.973851420910533,2.0089828325708456. Read each original CP strictly with default10 threads; require original t0 history, identical case and expected step, trusted prefix, zero spatial shift, the original H64×32 box and isotropic high-order choices. Rebuild only the paired one-dimensional7-point derivatives and quadrature; never restore or build a Poisson operator. Physical axes are x=X/Cx,y=Y/Cx and Omega=exp(logCx−logComega)(rho Dx′), with the original terminal physicalRhoXInf product checked bitwise.

The fixed quarter-annulus is r in[4,5], theta in[0,pi/2]. Register17 radial and129 angular support queries, including both axes. Do not adjust the band, angle, times or model after looking at fields. The previously saved physical box heights of the first two cases are7.0821 and5.2543, so the second case can approach its artificial top boundary; actual derivative support, not nominal endpoint size, must decide eligibility. This is motivation for an explicit support check, not a prediction that it will pass or fail.

Primary error samples are original positive-x native nodes in the band under the old source support mask: rows1:Ny−3 and columns4:Nx−3. Use the original positive tensor quadrature weights divided by Cx². The restriction of native weights to a curved band is a finite sampled quadrature, not exact annular integration. Save its total weight alongside the analytic quarter-annulus area9pi/4, all point counts and eight fixed equal-angle population counts. Eligibility requires every nominal query to lie in the old source support rectangle, at least5 distinct x and3 distinct y native coordinates, at least one node in each of the eight angular bins and at most10000 native model points. Missing support remains an explicit rejection; values on a supported subset do not establish full-band coverage.

The primary norms are

    E_inf_abs = max |Omega_actual−Omega_model|,
    E_rms_abs = sqrt(sum(w delta²)/sum(w)),
    E_inf_rel = E_inf_abs / max |Omega_leading|,
    E_L2_rel  = sqrt(sum(w delta²)/sum(w Omega_leading²)).

All maxima/sums use the same stated native subset. Omega_leading=f(theta)/r. The denominators are **band-global**, finite and positive. No pointwise division by f occurs near x=0. Keep the absolute quantities and denominators, not only relative errors. The registered exploratory observation gates are E_inf_rel<=.02 and E_L2_rel<=.02 plus support, finite/source checks and model numerical qualification; these are new full-band tests and do not replace the old sector's pointwise2% test.

Exact x=0 native samples are reported separately. Odd symmetry gives the continuous model trace Omega_model(0,y)=0; this is a labeled analytic boundary identity. The locked numerical characteristic helper requires X>0 and is neither called on x=0 nor modified to support it. Strictly positive near-axis points still call the unchanged helper, and every integration or domain failure is retained. The overall sampled infinity gate includes the measured actual axis discrepancy against zero, using the positive band's global leading denominator. Absence of native axis samples cannot silently produce this combined infinity qualification. Axis zeros do not conceal errors at nearby positive x.

The fixed polar query grid also provides a separate Hermite/bilinear observation. Use the same physical Omega and paired physical derivatives Cx Dx,Cx Dy. Because Hermite X jets differentiate Omega after rho was already differentiated, a distinct composed-support guard uses x7..x(Nx−6) and y0..y(Ny−3). Unsupported queries stay NaN. Its positive weights are explicitly independent trapezoid r dr dtheta weights, **not** the original solver quadrature. Save both interpolants, their difference and all subset/global norms; only complete query support can qualify a full Hermite comparison. This observation does not upgrade spatial accuracy merely because its model difference is small.

## Original model and old sector results stay fixed

Copy the original `frozen_homogeneous_characteristics.m` byte-for-byte (SHA c37ed7fc883ce49678bfb75cdefaa8c5ac4610825d140a5c12225324508af6d0). Use its fixed32-point angular Gauss rule, original exact k8 density, backward characteristic/Jacobian formula and RK steps32/64/128. The difference checks use global leading infinity denominators on each tested set with the fixed1e-8/1e-9 bounds. Do not adapt the time-step count after failure.

Before source comparison, run a **new full-angle numerical audit of that unchanged32-point formula** on the fixed129 angles. Compare Ic,Is,p,p′,U1,U2 with separate MATLAB integral evaluations at AbsTol=RelTol=1e-13, and retain the fixed1e-11 maximum absolute discrepancy gate. The earlier passing audit covered only theta<=asin(.25) and gives no whole-angle guarantee. If the full-angle check fails, keep all values and mark the model's full-angle numerical qualification false. Source/support diagnostics may still be saved, but no full-angle observation is qualified. Do not increase Gauss order or substitute a stable identity in the locked model to pass this experiment.

For each case, load its original saved wall-sector predictions and metrics from the development or prospective MAT record. Require that their original x,y,actual Omega and leading Omega arrays agree bitwise with the corresponding mask of the new strict CP-derived field. Copy the old model arrays and metrics verbatim; explicitly report `modelRerun=false`. The model SHA is unchanged. Thus the original pointwise-denominator errors and gates remain exactly their original records, rather than being recalculated with the new global denominators. This is an immutable evidence pairing, not a claim that a fresh sector model run was performed.

The stable identities U1=−integral(theta,pi/2,cos(s)f(s)ds), U2=integral(0,theta,sin(s)f(s)ds) follow from the same continuous model and avoid head-integral subtraction or p/p′ cancellation near the axes. They are a possible **separate future numerical implementation**, not part of the present locked model. If needed after this audit, they require their own registration and must preserve the original32-rule failure and prospective sector results.

## Executable preparation

`ipm_accellab_full_quadrant_extension(registrationFile,outDir,executeObservation)` defaults to false. The prepared launchers explicitly select the frozen runtime/helper directory and cwd without a competing +ipm. The true branch also requires `IPM_FULL_QUADRANT_OBSERVATION_WINDOW=root_authorized_noLU_observation`; no permission to run it was inferred from preparation. Exact file manifests include the three native CPs, original model/sector MAT records, runtime, helper and registration. SHA and package resolution are checked before and after observations. Call profiling must exclude mesh.build, restoreCheckpoint, field.velocity/poisson and evolve.flow/advance. Candidate/PDE/native-checkpoint creation is absent. Failures and per-case arrays go only to a new directory.

The helper is newly prepared and has **not yet been run or CodeAnalyzer-checked** because the root reserved the MATLAB queue for other tests. No numerical pass, support count or full-angle discrepancy is asserted here.

## Minimal subsequent complete-boundary trace MMS

Only after the new model audit is understood, prepare a separately registered pure-integral trace check. Reuse the already verified **complete** quarter-plane Green kernel and exterior-of-rectangle geometry in `ipm_accellab_green_full_tail_trace`; do not replace the full boundary by an inner-point xy A/H expansion. Keep the original42 boundary targets (corners, near-wall/near-axis points and both sides/top) and three interior holdouts, with independent origin/target-centered quadrature and at least the previously required256/512-order comparison. There is no solver injection or LU in this stage.

A particularly small additional manufactured check is already compatible with the existing perturbation source:

    Omega_1(x,y) = d_x[x²/(x²+y²)^(3/2)]
                 = (2 cos(theta)−3 cos(theta)^3)/r².

Its full-quadrant potential with the current −Laplacian sign is

    psi_1(theta) = 1/3 − 2 theta/(3 pi) − cos(theta)/4 − cos(3 theta)/12.

It vanishes on both axes and satisfies −psi_1″=2cos(theta)−3cos(theta)^3. An independent interior Green integral plus the existing complete exterior perturbation trace must reproduce psi_1 on every registered boundary/holdout target. Treat the origin and boundary log singularities with separately split/regularized quadrature; no arbitrary finite radial cutoff may be mistaken for the infinite integral. Preserve the old full-trace gates: max(H,|psi|)-normalized discrepancy1e-10 and nonzero near-axis H psi/(qx qy) discrepancy1e-9, along with independently tightened interior quadrature. This manufactured source is not an unforced IPM solution; it tests the full1/r² trace and its units, not evolved closure.

For isotropic zero-shift scaling, a physical source maps as Omega_can(S)=Comega/Cx * Omega_phys(S/Cx). A leading f/r trace therefore scales as Comega, while a physical B(theta)/r² perturbation scales as Comega*Cx. The full boundary potential of the latter is generally order Comega*Cx, whereas its fixed-inner-point length-rate effect can scale as Hc^-2. These different quantities must not share the same error label. The exact frozen-characteristic k8 source also depends on physical time and the finite initial length scale; its complete evolving trace cannot in general be reused by multiplying one leading template solely by Comega.

Only a later independently registered trace quadrature would substitute the full frozen-characteristic source into the exterior kernel. That would remain a **conditional model trace**, even if it matches all three finite bands. These finite observations do not determine the actual source throughout the unobserved exterior, certify a self-consistent velocity, or authorize adding any correction directly to C.
