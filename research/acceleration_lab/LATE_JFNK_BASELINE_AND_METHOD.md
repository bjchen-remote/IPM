# Late9205: original-RHS baseline before matrix-free acceleration

This is a new, bounded preparation. The executable entry point only reads the actual native source by default. Its separately guarded true branch makes one saved-grid factor and one original `ipm.evolve.flow` call, releases the factor, and saves/readbacks the true arrays. It never constructs a candidate, steps the PDE or writes a native checkpoint. No true branch has been run in this preparation.

The source is the same895×386, remesh10, step9205 state qualified in `JF_LATE_CROPPED_QUALIFICATION_PROTOCOL.md`: native tau.1661212412279467, parent-equivalent tau12.443200751859164, absolute physical time2.0235562096101787. The original frozen runtime is `fresh_adaptive_source_v5`; its exact transport-anchor/quadratic-peak rule, conservative wall flux, global-LF WENO and Green boundary are unchanged. Existing tiny three-direction Gauss–Newton evidence is not a full JFNK solve. The failed late temporal rank1 holdout, secant tail gate and nonmonotone spatial residual results remain failures.

## A bare Newton step has an amplitude shortcut

Write F(rho) for the complete discrete rho derivative with the original grid, references and instantaneous rule. Hold the stored scale tuple fixed during a perturbation. For positive amplitude a, assume all relevant numerical branches stay fixed and homogeneous:

- Omega=D_x rho and the Green/Poisson velocity are linear in rho, so u[a rho]=a u[rho].
- The transport anchor gives c_l=−I[u1](x_anchor)/x_anchor, hence c_l[a rho]=a c_l[rho]. Translation is the original zero/symmetry-compatible rate; it is not an independent fitted phase correction.
- The complete geometric transport base B is degree2 if the LF split, WENO reconstruction, boundary treatment and correction retain their homogeneous branches. In particular the density flux scales by a², while the unit/free-stream correction scales by a and is multiplied by a rho.
- A fixed three-node quadratic peak has P[a rho]=a P[rho], the same vertex and evaluation weights. Its original c_omega=−DP[B]/P therefore scales by a. Consequently F[a rho]=B[a rho]+c_omega[a rho]a rho=a²F[rho].

Euler differentiation would give J_F(rho)[rho]=2F(rho). Thus delta=−rho/2 solves the bare Newton linear equation J_F delta=−F. The candidate rho/2 then has F/4 but also half the actual peak, velocity and density amplitude. Its intrinsic shape evolution is not accelerated. For the complete C1 geometry, U is amplitude invariant, G=U_tau scales by a, and P scales by a; **G/P remains invariant**. Reporting raw F or G alone would mislabel amplitude reduction as progress.

This is a conditional formula, not an asserted exact property of this actual floating-point operator. Fixed reference validity floors, degenerate or tied selectors, Green support/compression, epsilon/clamping, LF maxima, splitScale zero/floor branches and rounding may break exact homogeneity. The complete original flow must measure it. At a fixed factor, the minimal useful check needs **two additional fresh RHS calls**, at(1+e)rho and(1−e)rho with e=1e-5 fixed in advance. Save both actual fields and rates, compare F± with(1±e)²F, compare(F+−F−)/(2e) with2F, and retain the selector/validity evidence. Baseline plus these controls costs3 calls; the present baseline-only runner does not execute the two controls. A single positive probe tests finite homogeneity but does not establish a two-sided directional derivative. Failure to meet a future registered discrepancy gate is a rejection, not permission to assume the identity or adjust e.

## The actual constraints and their linearization

The transport anchor fixes an instantaneous velocity balance, I[u1](x_anchor)+c_l x_anchor=0. It does **not** pin the actual peak position. Its differential is I[delta u1](x_anchor)+delta c_l x_anchor=0, with delta c_l recomputed by the original rule at every fresh perturbed field. Keeping c_l or c_omega fixed would differentiate a different operator. The saved pin/tracking references and anchor interpolation support remain fixed; the active peak index does not become a user-chosen new phase point.

On a fixed nondegenerate quadratic stencil, let Q_rho(a) be the quadratic interpolant of D_x rho at its interior maximum. Then

    DP[v] = Q_v(a),
    a_v = −Q_v'(a)/Q_rho''(a),
    D²P[v,w] = −Q_v'(a) Q_w'(a)/Q_rho''(a).

Here Q_v uses the same three nodes and paired maintained D_x; Q_rho'' is twice the polynomial's stored quadratic coefficient. The original exact amplitude rule gives DP[F] approximately zero within that branch. Differentiating this constraint gives

    DP[J_F v] + D²P[v,F] = 0.

Therefore J_F does **not** generally map ker(DP) back into ker(DP). Ignoring the curvature term, freezing the quadratic weights across a true field change, or setting observed PPrime to zero would be mathematically wrong. At nodal selector changes, the different three-point quadratics can have a finite peak jump; the smooth formula cannot be continued through the event. The baseline records the actual full-product DP[F], the smaller stencil-product floating path and the original algebraic residual separately. They need not be bitwise identical. It also records the active selection and vertex margins. These are instantaneous measurements, not proof of P conservation across events.

A future vector-only projection Pi v=v−rho DP[v]/P can remove the radial direction from search vectors because DP[rho]=P. It must never be applied to the actual state or used to reset P. Pi J Pi is a projected search operator, not a claim that the original J preserves the tangent space. A clearer square formulation is the bordered Newton search

    [ J_F   rho ] [delta] = [−F]
    [ DP     0  ] [ mu  ]   [ 0].

Its multiplier has inverse-time units. This is a constrained linear search with a radial multiplier; it is not the original unbordered Newton equation unless mu=0. The full true candidate residual, amplitude/phase changes and every matching gate must still decide whether the direction helps. A small projected or bordered residual alone is insufficient. A finite linear tangent step can change P at second order and may change its selector; there is no nonlinear retraction/projection to hide that change. The future candidate is an independent research profile, never an accepted state of the original trajectory.

## Square GMRES and a rectangular observer are different problems

F has one component per native density node. A genuine matrix-free GMRES uses Arnoldi products from this square Jacobian, or from the explicitly labeled(N+1)-dimensional bordered search. It can be truncated at a small Krylov dimension; this is an actual Krylov iteration, not a promise of convergence of a full nonlinear solve. Record its entire Hessenberg matrix, orthogonality defect, projected residual and fresh full residual. No Sylvester/Schur or alternate flow formula is introduced.

The complete local C1 residual G=U_tau maps a global density field into values on a registered inner observation rectangle. Its Jacobian is rectangular. Ordinary square GMRES cannot simply be applied to it. Normal equations or LSQR would require an adjoint/J-transpose mechanism that this cheap forward-RHS proposal does not possess. The earlier three-vector observer Gauss–Newton fit and its tail ellipsoid remain correctly labeled reduced least-squares searches. They are not renamed JFNK.

Use fixed source physical units and fixed positive native quadrature weights. F is d rho/d native tau at fixed computational coordinates. Multiplication by Cx/Comega converts the clock to physical time but does not turn it into the physical Eulerian density derivative; the coordinate and amplitude chain terms are still needed. G/P=U_physicalTime/(actual physical C1 Omega peak) is a useful amplitude-invariant shape diagnostic. Report raw F, G and G/P separately with actual PPrime; never optimize a candidate-dependent denominator and omit the raw quantities.

The whole box and its density unknowns remain present. No Dirichlet density values are invented: the original transport ghost/boundary machinery and source-dependent Poisson boundary act on every perturbed field. Any symmetry restriction on search vectors must follow the actual source symmetry and be stated as a linear domain restriction. The source's transport boundary mode is not swapped. Tail/phase constraints do not authorize modifying those boundaries.

## A possible later20-call test, not yet enabled

The present code authorizes only the first line below. A candidate protocol must be frozen and separately authorized after actual baseline observation/cost results:

| Operation | Fresh original RHS calls |
|---|---:|
| This baseline |1|
| Same source rebaseline after a later factor rebuild, require true array parity |1|
| Positive/negative amplitude controls, fixed1e-5 |2|
| At most4 Arnoldi directions: +h,−h,+h/2 per direction |12|
| Combined proposed direction at +h,+h/2 to test Jv consistency |2|
| Exactly one fresh candidate and its identical-state RHS reproducibility check |2|
| Total across both phases |20|

For every direction the +h versus +h/2 forward difference provides a two-step discrepancy; +h versus −h supplies the one-sided gap and the centered product used by Arnoldi. This is **not** the old four-call two-centered-epsilon stencil and must not be represented as the same test. The proposed fixed gates remain2% for two-step discrepancy and one-sided gap, with actual branch changes explicitly rejected. All finite-difference directions are normalized in a source-fixed weighted density norm and use h=1e-5 relative to the source density norm. The candidate direction also has a fresh two-step combination check against the Arnoldi linear combination, so incompatible one-sided LF branches cannot be assembled into an assumed smooth Jacobian. The total budget is a cap: breakdown, invalid geometry, absent matching coverage or failed checks stop the experiment early. No rank, damping, epsilon or window scan follows failure.

The prospective search retains the previous trust radius.002, singular floor1e-8, actual peak change.001, full native densityInf.003, source/peak and tail density/source/velocity.001, baseline phase.1 source wall widths, minimum raw and intrinsic core decrease.001 and prediction agreement.1. Legacy physical matching nodes remain explicitly coarse; source-fixed native-node full-box and split-cell near-tail observations are additional hard checks, not replacements. Core/holdout/wall and sampled-maximum refinement must all be reported under both source-fixed and candidate-native geometry. A small linear multiplier or constraint root is not an acceptance gate by itself.

Even passing these instantaneous numerical gates would not establish a resolved scientific gain. Late matched-grid residuals have shown substantial nonmonotone spatial sensitivity. A single step9205 baseline cannot set a converged residual floor, and its roughly22/33 core counts do not cure that limitation. A future independent matched physical-time and spatial control is necessary before claiming acceleration. The20-call cost is already comparable to several original RK steps, so an eventual cost comparison must include natural original-PDE progress. No such comparison or candidate exists here.

## Baseline implementation and guards

Call `ipm_accellab_late_jfnk_baseline(registrationFile,outDir,executeBaseline)`. The default is false. The registration freezes the exact CP, original runtime and the seven observer/runner files by SHA256. Both before and after execution the helper checks hashes and `which` paths, and requires cwd without a competing `+ipm`. A strict10-thread native read verifies the actual signature and trusted history. Direct fresh-parent lineage and exact clocks/references are checked from their original files.

The true branch requires an additional explicit root resource-window environment marker. It starts **directly** with `mesh.build(savedConfig,savedGridOverride)`, restores original runtime references using the maintained allowlist and asserts full equality. It never calls initializeScaling, initialize or restoreCheckpoint. The latter would first build the configured grid and can build a second factor when its axes differ from the saved remesh grid. Paired7-point D_x/D_y, the original quadrature weights and saved base axes are independently checked. The only fresh evaluation is original `evolve.flow`; its Omega must match rho*D_x' exactly.

A read-only original `output.record` replay uses the saved timestep telemetry and compares every terminal common/gauge/mesh/anisotropic field bitwise. It does not append to the actual history or produce a checkpoint. The source checkpoint signature is recomputed unchanged after the work. The cache contains source state/history/signature, complete rhoRate/flow/Omega/FX and an operator array view with A and the factor removed. All cache fields are written to a temporary file and read back exactly, including function-handle definitions where present, before atomic installation. No decomposition is serialized. The factor is cleared before cache output and the no-LU continuous geometry/residual observation. Profiling must show one build, one flow/RHS/velocity/Poisson call and no lifecycle/step calls; default mode must show zero of each.

Failures preserve the unique output directory, arrays where available and a structured failure report. Geometry invalidity is retained as unavailable/failed observation, never repaired by P projection. A small Poisson residual is not a forward-velocity certificate. The helper's factor-cleared flag and call counts are not a substitute for an external peak-memory/swap watchdog. Actual large-case cost, parity, F, DP[F] and residual values remain pending until the separately granted execution window.

## Additional scientific prerequisite: is a steady target appropriate?

The original transport-anchor/quadratic-amplitude gauge need not make the desired late profile a nonzero smooth stationary density. Its anchor is a transport-velocity constraint, while the actual core can continue shrinking and the peak can move. An instantaneous small or zero DP[F] does not settle this question. Neither does a small finite-grid F after Newton iteration. This is an unresolved target question, not a claim that a steady limit is impossible.

Before any actual candidate, use the real baseline P, rho amplitude/range, full F, peak position/rate, both width log-rates, complete U_tau and G/P, each in stated native/parent/physical units. Compare geometry motion with the remaining core/holdout/wall shape drift. Retain the existing chronological and spatial evidence; a single frame cannot establish existence or nonexistence of a smooth steady profile. Persistent narrowing might correspond to a non-smooth or multi-scale limit whose representation is poorly suited to F=0 in this gauge. It cannot be erased by a fitted phase/width correction or treated as a stationary residual success.

Until target suitability and spatial resolvability are established, a JFNK result would be labeled only a constrained, finite-grid residual experiment with a changed research profile. It would not identify the original continuous infinite-time limit, prove singularity or qualify an accepted point on the original trajectory. The current frozen baseline runner remains unchanged; its preparation package has a separate `steady_target_admissibility.json` with status pending and the required actual fields. This addendum adds a scientific qualification boundary without changing C, the original data, or an earlier rejection.
