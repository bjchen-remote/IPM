# Late cropped tail-ellipsoid JF candidate: preparation, not execution

No `tail_ellipsoid` full-field candidate has been tested on the actual late cropped case. The only executed constrained records are `jf_profile_tpe83cbd24_ced9_4032_800c_a94a7a505d7b/n49` and `/n65`: Gaussian density on[-4,4]×[0,4], native quadratic/transport-anchor rule through canonical.08, followed by independent candidate/baseline relaxation to physical.10. They showed about0.3% residual improvement, smaller than the roughly1% baseline residual difference between those two grids. The original unconstrained failure and the unrelated late memory4/damping1/16 secant tail-velocity rejection remain unchanged.

The old method is a three-direction whole-field constrained Gauss–Newton construction using F,JF,J²F. It is not a full-dimensional JFNK solve of the stationary physical PDE. Its candidate changes the initial field while retaining its clocks and original scaling references; even a successful matching-time relaxation is an independent profile IVP, not a later state of the original initial-value solution.

## Actual no-LU source preflight

A new preregistered pure check used the immutable `fresh_adaptive_source_v5` and strictly read step9205 with default10 threads. Native signature, full history trust and direct fresh-parent lineage passed. No restore, velocity, Poisson, fresh RHS, candidate or PDE step was evaluated.

- Cropped case895×386, remesh10, native canonical time0.1661212412279467.
- Parent-equivalent tau12.443200751859164; absolute physical time2.0235562096101787; physical epoch1.9545271366857122.
- Native core21.9252321/32.9915975. These trusted counts are not a certificate of converged shape residual.
- Physical box[-1233.8787893,1233.8787893]×[0,1271.4084501].
- Actual C1 physical peak163.9233059 atx0.1249087373; wall/vertical0.9 widths3.122475282e-5/7.191828984e-6.
- Original dynamic isotropic, exact_gauge_no_feedback_v1, weno5_fd, conservative_flux; transport anchor0.152142235.

The old physical matching rectangle[-3,3]×[0,3] lies inside the box, but its49×25 uniform grid has **zero samples inside the actual wall0.9 core and only the wall node inside the vertical core**. Its spacing is4003 and17381 respective widths. It still measures the stated distant tail and remains useful as a legacy gate. It cannot bound the density discrepancy throughout this late observation window. Actual legacy density/Omega and actual normalized129×97 inner-shape arrays were saved; velocity, F and PPrime are explicitly absent, not inferred from history.

The source and results are frozen at `result/longtime/20260908_campaign_v1/acceleration_lab/jf_late_qualification_preparation_de2ed09db8c2476482dc1806dd111d86`. No original CP or schema is modified.

## Why the old runner cannot simply be called on9205

`ipm_accellab_test_jf_profile.m` explicitly initializes the Gaussian fixture, asserts Nx<=65/Ny<=33 and no remeshing, uses its saved tiny rho/z, and hardcodes a physical target.10 and fixed coarse matching nodes. Its addpath logic also targets the working maintained package. It must not be used as a large-case shortcut. A future harness must instead start from a strict immutable native CP, restore its existing operators/references exactly once, explicitly fix cwd/path and verify which/source hashes. It must never call initializeScaling or refreeze the strain/gauge references.

The complete original F(rho;fixed scale,refs,ops) must be used, including newly computed Poisson velocity, original instantaneous C and velocity-dependent global-LF WENO assembly at every perturbed rho. No global-LF or boundary-mode replacement is part of this experiment. Continuous C1 P,a,wx,wy and their genuine derivatives must be recomputed from this F; the original quadratic gauge PPrime and C1 PPrime are distinct. Neither is forced to zero by a projection or a fabricated observation.

No future field or future peak position enters construction. Fixed objective queries, norms, matching masks and normalization denominators are registered from the **source** alone. Candidate-local geometry is recomputed from the candidate for an additional independent residual check. Small oracle-aligned shape forecast error from the temporal study is not used as evidence for a predicted physical field.

## Proposed bounded instantaneous protocol (all execution disabled)

Keep the old construction parameters unchanged: three whole-field weighted-orthonormal directions, difference increments1e-5 and5e-6, trust radius.002, singular floor1e-8, Q tail/peak budget fraction.8, tail density/source/velocity hard tolerances.001, physical peak change.001, observation densityInf.003, minimum raw and intrinsic core decrease.001, Jacobian two-epsilon discrepancy<=.02 and prediction-agreement>=.1. These retain the old tiny gates; they do not reinstate any previous rejected candidate.

Before any expensive trial, freeze a complete new research registration including:

1. Exact CP, runtime and observer hashes; source case/clock/lineage/refs; one-factor resource reservation and measured actual baseline RHS cost. There is no existing9205 true-RHS cache in this preparation.
2. Same source-fixed inner rectangles as before: core[-1,1]×[0,1.5], full local[-2,2]×[0,3], with split native-cell Gauss4 for true U_t and G/P, wall and holdout norms, both fixed-baseline and candidate-native FD masks. Retain every original no-increase gate and the.1% core-decrease gate.
3. Preserve the legacy[-3,3]×[0,3]49×25 physical tail/peak checks **and add** source-fixed native-node full-box densityInf and physical-source matching checks. Their relative thresholds remain.003 for densityInf and.001 for sourceInf normalized by actual baseline physical peak. Do not label the old uniform-window value “full field.”
4. Add a registered source-inner matching rectangle[-4,4]×[0,6], sampled via split native cells; its near-tail region is|xi|>=2 or eta>=3. Require density L2, sourceInf/physicalPeak and velocity L2 within the same.001 hard matching tolerances, and preserve a separate baseline-fixed phase error in units of actual source wall width. The proposed phase bound is.1 source widths; it is an additional future gate, not a change to any past gate. No candidate-dependent physical sampling grid is allowed to make matching easier.
5. Record each original quadratic selector, continuous peak/root cells and margins, and directional WENO/global-LF active-point evidence. The existing twelve observer ±probes already supply forward/backward differences: report their relative gap as well as the two-epsilon centered gap. A proposed.02 directional gap gate rejects an unmodeled nonsmooth response. A centered difference across changing max selectors is not automatically a smooth Jacobian. Do not discard one side or silently reduce epsilon until a preferred result appears.
6. Construct at most one candidate from the frozen ellipsoid. Fresh original Poisson/RHS and all actual matching/residual gates decide it; no damping/rank/epsilon/window scan follows rejection. Positive/negative amplitude controls remain independent controls, not rho projections. Save candidate/rho/basis/J/Q, all failures, actual PPrime and full fields with an unmistakable non-native research label. No native checkpoint or accepted time is produced.

A compact phase0 source-cache audit should precede the full20-RHS construction: verify original native flow/RHS pairing, finite fresh C1 derivatives and root transversality, sampled-infinity refinement, masks and current spatial-error evidence. The present geometry-only preflight cannot establish any of those forcing/residual quantities. The previously measured nonmonotone multi-grid U_tau sensitivity is materially larger than the old0.3% tiny benefit; it prevents treating a similarly small late residual reduction as scientifically resolved without an independently matched spatial comparison. Core count alone is insufficient for that decision.

## Cost and later physical-time qualification

The old construction has an exact budget of20 fresh RHS evaluations: baseline1, two centered basis products4, three observer directions at two epsilons12, candidate1 and amplitude controls2. In isotropic geometry the same factor is reused. This count excludes strict restore/build overhead and any independent source spatial comparison. No9205 seconds or peak-memory bound has been measured; it would be misleading to extrapolate an H8/H64 or small-grid number. Three full-field direction arrays at895×386 cost about8.29MB of raw doubles, but the LU/Poisson data dominate resource admission. Only a root-authorized single-factor slot with measured memory and stop conditions can enable execution.

Twenty RHS calls are already about the original four-step budget (21 with a cold first cache). Any apparent gain must be compared with natural original-PDE progression using the same actual RHS budget, preserving real clocks. That distinct-time control is not a replacement for a common-physical-time comparison.

If an instantaneous discovery trial passes, freeze a separate short relaxation protocol before running it. Keep source and candidate at the same actual physical clock/scale tuple, advance rho and every scale variable with the original PDE/C rule, and use the same **absolute** physical target, converted to the fresh local clock by the saved epoch. A suitable bounded target can be registered from the source's original selectTimestep result before candidate evaluation, with two time-step capsdt0 anddt0/2 and a finite step/RHS cap. Bracketing original RK states may be observed with the existing consistent physical-time Hermite rule; no time relabeling or P correction is allowed. Preserve actual finite/CFL/trust limitations, each hard residual/tail gate, two-step sensitivity, and an independently matched spatial check. These would remain labeled independent profile IVPs. They do not provide a shortcut to the original solution's infinite-time state.

## Current decision

Do not spend a20-RHS large-case window based only on the old tiny pass. The new late source has no fresh forcing qualification here, the old physical matching samples miss its core, and a0.3% benefit would be below currently unresolved shape-residual spatial sensitivity. The useful next step is an already-authorized solver boundary cache or a separately granted one-state baseline RHS audit, followed by the registered one-candidate construction only if its observation and resource prerequisites are satisfied. This document and the no-LU arrays are preparation; all large-case JF execution flags remain false.
