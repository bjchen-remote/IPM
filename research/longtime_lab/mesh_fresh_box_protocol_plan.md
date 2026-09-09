# Registered independent late finite-box protocol

The proposed experiment starts from the accepted stage-1 native checkpoint
at tau=4.283201944539746, step2031 and physical epoch
1.7798569844190184. It investigates the effect of removing outer nodes from
this late state. Its outcome cannot validate the parent's preceding dynamics
or establish convergence as the physical box tends to infinity.

```matlab
campaign=fullfile(pwd,'result/longtime/20260908_campaign_v1');
sourceFile=fullfile(campaign,'adaptive_campaign/stage_001', ...
    'checkpoint_20260905T204759339_dynamic_isotropic_1025x513_tp8c8b7442_a9d7_4682_b700_7cf8af6fe3f7_step0000002031.mat');
controls=struct('allowLarge',true,'maximumSteps',64, ...
    'staticReportFile',fullfile(campaign,'late_box_tau428_v1/late_box_report.json'));
report=mesh_run_fresh_box_protocol(sourceFile,.01,.001, ...
    fullfile(campaign,'fresh_box_tau428_protocol_v1'),controls);
```

The corresponding prepared launch script is
`/private/tmp/ipm_fresh_box_tau428_protocol_v1.m`. It changes into the maintained
project directory explicitly. It has not been launched by the mesh agent.

## Registered cases and units

Every branch advances the same physical interval .001 to absolute time
1.7808569844190183. The five cases run sequentially in this order:

| Case | Initial mesh | Public epsilon | Native-history meaning |
| --- | --- | --- | --- |
| native_continuation | 1025x513, original coordinates | Inherited | Original same-grid checkpoint, old case/history |
| full_natural | 1025x513, physical coordinates | Inherited | New late-state IVP, new ID/history |
| crop_natural | 901x386, physical coordinates | Inherited | New late-state finite-box IVP |
| full_common_epsilon | 1025x513, physical coordinates | 1/385^2 | New late-state IVP, epsilon control |
| crop_common_epsilon | 901x386, physical coordinates | 1/385^2 | New finite-box IVP, epsilon control |

Fraction .01 means original rescaled x half-width 9792.523085957766 and y
maximum 9986.695463288403. All retained initial inner nodes and rho samples
are exactly the same subset of the parent physical data. The actual physical
bounds are these values divided by the initial Cx; H~1e4 refers to the parent
rescaled units, not the fresh physical coordinate system.

At the source, Cx=exp(1.5279970171255546),
Comega=exp(-.6604814239795477), and lambda=Cx/Comega=8.921627999603581.
Fresh coordinates, density and anchor are parent X/Cx, R/Comega and 1/Cx.
Fresh rates divided by lambda are comparable with parent canonical rates.
The initial estimate of the parent canonical interval is .00892163; this
estimate is not used to impose a different endpoint. Fresh local elapsed
time starts at zero, and absolute epoch is separate lineage metadata.

All numerical choices and diagnostics come from the parent frozen schema4
quadratic/transport-anchor contract. Fresh maxDt/minDt/outputEvery are divided
by lambda for covariance. Source CFL=.5 and maxDt=.004 remain the reference;
fresh maxDt=.004/lambda. There is no CFL halving in this first protocol.
Every branch has a hard budget of 64 accepted steps (native maxSteps is
parentStep+64); a step-limit stop before the requested endpoint fails.
Automatic remeshing is disabled, and no transfer or domain-changing native
restore occurs. Fresh runtime references are explicitly reinitialized.

Natural effective WENO epsilon values are
9.5367431640625e-7 / 3.814697265625e-6 in x/y on the full mesh, and
1.2345679012345679e-6 / 6.7465002529937596e-6 on the cropped mesh.
The common public epsilon 6.7465002529937596e-6 fixes both axes in both
control branches. The full-natural/full-common difference is reported
separately. Global linewise Lax-Friedrichs maxima, Green traces from cropped
source and outer derivative closures can still change.

## Gates and interpretation

The matching static report must have the exact source filename, fraction,
node counts and a passed/eligible entry. The archived checkpoint undergoes
normal native validation at 10 MATLAB threads. All five production branches
must reach the same endpoint within 1e-11, retain entirely trusted histories
and finite physical fields, and have both actual terminal cores >=14 and
safety <.70. Source monitoring gives cores18.2408/18.6087 and safety.43858;
that low-cost ctypes read is not a native signature validation. The 20-cell
transaction gate is unchanged and irrelevant to these no-transfer branches.

The full-native/full-fresh gate runs before any cropped PDE: coordinate
relative infinity <=1e-7, physical rho/rho-x/rho-y relative infinity <=1e-5,
and covariantly compared rate relative differences <=1e-4. These declared
research overlap tolerances do not replace checkpoint signature checks.

Natural and common-epsilon full/crop pairs are compared bidirectionally in
the fixed physical laboratory window [-2,2] x [0,1] divided by source Cx,
without shifting or inner-profile normalization. Both directions use paired
axes, linear interpolation and trapezoidal L2 integration. The short-interval
candidate thresholds are rho L2/Inf<=1e-3, both gradient L2/Inf<=2e-3,
physical gradient maximum/quadratic peak relative change<=2e-3, physical
strict core-width changes<=.005, and canonical rates expressed in parent
units differing by at most max(1e-6,.002*abs(full rate)). All raw metrics and
failures remain saved. Passing only nominates a longer independent box test;
the report never promotes production use or claims dynamic box convergence.

Further evidence should extend the physical interval, add an intermediate
larger box, halve CFL, and refine the shared inner grid. The controlled
epsilon comparison must remain separate from natural scheme behavior.
Independent fresh roots from the original initial state are still needed
to address earlier finite-box history. Freed points can be redistributed
only in a separate grid experiment after fixed-inner-node box effects are
quantified.

## Resources and actual tiny evidence

The protocol runs exactly five PDE cases and five elliptic assemblies, one
at a time, at 10 threads. No additional parent-reference endpoint restore is
included. Each ipm.solve returns a result without live operators; the wrapper
saves it and retains only compact physical fields/scalars before the next
solve. No original LU is held while constructing a fresh branch. Four
compact observations consume roughly 50 MB of field arrays; full case peak
memory should be budgeted from the existing q512 solver, not inferred from
the cropped node ratio. Cropping retains 66.14% of nodes, which is not a
measured LU-memory reduction. Every case saves an immutable native checkpoint,
complete result, lineage and compact observations. A failed gate leaves all
completed/partial case artifacts and prevents subsequent dependent cases.

The preliminary 65x33 covariance/reference fixture passed in
`result/longtime/20260908_fresh_box_tiny_v1`. Full fresh/native rho/rho-x/rho-y
differences were 6.75e-16/1.34e-13/1.10e-13; transformed-parent vs fresh
runtime references gave zero instantaneous RHS/rate changes in full/crop.
Those audits used two extra tiny restores; the q512 protocol omits them.

The five-case entry point then passed a separate tiny run in
`result/longtime/20260908_fresh_box_protocol_tiny_v1`, with its adjacent log.
It reproduced rho-x covariance 1.34427e-13, rejected the deliberately sensitive
.75 crop in both natural/common epsilon pairs, preserved five distinct IDs,
and passed Code Analyzer. The tiny fixture's resolution exception is labeled
and cannot apply to q512. No q512 fresh-box PDE was executed while preparing
this registration.

## Actual q512 execution after allocation

The coordinator subsequently allocated a single window. The frozen helpers
and hashes are in `20260908_campaign_v1/fresh_box_tau428_protocol_source_v1`;
the result is `fresh_box_tau428_protocol_v1/protocol_report.json` under that
campaign, with an adjacent `.log`. The process exited normally with code0
and released the window. Every branch reached the physical endpoint within
2.45e-14, stopped for `physical_final_time`, retained complete trusted result
and native checkpoint histories, and paired rho/axes/step/clocks/case ID
exactly between the terminal checkpoint and result. The post-registration
native-pairing check was also separately tested in tiny v2 before this run.

| Quantity | Natural-epsilon crop/full | Common-epsilon crop/full |
| --- | --- | --- |
| Physical rho relative L2 | 4.64273e-7 | 4.64270e-7 |
| Physical rho-x relative L2 / Inf | 7.61915e-7 / 3.68232e-6 | 7.61912e-7 / 3.68230e-6 |
| Physical rho-y relative L2 / Inf | 7.11025e-7 / 4.22715e-6 | 7.11022e-7 / 4.22713e-6 |
| Physical gradient maximum relative difference | 1.54980e-7 | 1.54979e-7 |
| c-l relative difference in parent units | 8.81875e-5 | 8.81861e-5 |
| c-omega relative difference in parent units | 1.74049e-8 | 1.74068e-8 |

Full fresh/native physical covariance gave rho/rho-x/rho-y relative Inf
1.45e-15/2.64e-13/1.14e-10 and rate differences 2.71e-12/1.55e-12.
Changing epsilon alone in the full box gave rho L2=2.62e-14 and a largest
reported gradient Inf change of 4.60e-12. Thus the selected regularizer
control is negligible on this short interval, without establishing that it
will remain negligible on longer integrations or different profiles.

All final cores were at least 18.12965/18.54205, and safety was at most
.441267. The observed solve wall times in registered order were
66.2/54.0/26.6/48.2/28.5 seconds. These shared-machine observations include
initialization/factorization and are not an isolated speed benchmark.
The original full-box decomposition's small-RCOND warnings remain in the
log; none was suppressed, and the registered trusted/covariance gates
passed. A subsequent no-LU scalar-history audit
`fresh_box_protocol_record_audit_v1.json` found maximum saved Poisson
residual5.87e-10. That is a supplemental post-run check: trustedMask itself
does not enforce a Poisson-residual threshold. Only the declared
short-interval comparison passed. The report
explicitly leaves `dynamicBoxConvergenceEstablished=false` and
`productionDomainModified=false`.

A reasonable next independent ladder is elapsed physical time .005, then
.02 if the same actual core/safety gates remain valid, with an intermediate
larger box and CFL/inner-grid refinement. Each new run must report its
absolute epoch and enforce the physical target; a max-step stop is not a
successful extension. Independent roots from the original initial state
remain necessary to assess prior box history.
