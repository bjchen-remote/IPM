# H64 directional operator inputs: preparation only

This experiment is a pure, non-native resource-input preparation. Both inputs come from the same strictly validated original-t0 H64 native checkpoint, step2353, canonical time4.600206263141828 and physical time1.8352642229337874. It uses the immutable `node_family_v3_directional_candidate_v2/source`; it does not upgrade the source's version2 config or ledger. The future version3 family is generated independently from the actual saved initialization selectedBase. Its first four complete members match the native family bitwise. Members7/8 have factors[3,2]/[2,3], node counts961×321/641×481, with explicit310000 cap.

The complete original70×4 bounded planner was executed once for each direction. Only its first selected pair was scored and interpolated. Both original frozen-input gate sets passed; other candidates were not used after scoring. All axis trials, whole-grid quality, fields and rejected-case diagnostics are retained. No failure or acceptance gate was changed.

| Direction | Interpolated X/Y core | Left front | Peak relative jump | Conservation defect | Range violation |
|---|---:|---:|---:|---:|---:|
| 961×321 |35.1599073 /37.0807325|34.1661566|8.23370672e-6|3.33126017e-15|9.48519240e-9|
| 641×481 |35.1353796 /36.9880135|32.3278737|8.75041376e-6|2.34422012e-15|1.24210173e-9|

Preparation took9.4194s. The profiled call graph contains no restore, mesh build, Poisson, flow, native transfer or PDE advance. The source field, config, clocks, references, ledger, native signature and original files remain unchanged.

The density inputs use the maintained six-point constant-conservative X-then-Y interpolation. They are explicitly **not** the native transfer's wall-trace-preserving Y bubble transaction. Whole-field rhoY round-trip relativeInf is0.003055/0.002513 and rhoY peak difference0.002749/0.002326; these are preserved as additional diagnostics and not claimed to certify physical-space accuracy. Passing the listed preparation gates gives no native migration, PDE, LU resource or long-time qualification.

Artifacts are under `result/verification/autonomous_runtime_20260909/h64_directional_operator_preparation_eddb13c7332b4da5b574e818ea103b7b`. Each `run/member_*/resource_bundle.mat` is a research bundle, not a checkpoint. Source and input hashes are recorded before/after preparation. `run_member_7_dormant.m` and `run_member_8_dormant.m` both pass `false` to the frozen resource helper. In this mode the helper only reads the bundle and writes an unqualified report.

A future explicitly authorized and externally watched resource call may build one target factor and execute one original frozen `ipm.evolve.flow`, with the original config and runtime references (only origin/pin indices paired to the target axes). It never restores an old factor, advances time, commits a transfer, or writes a native checkpoint. The new factor is cleared before output. Its report distinguishes actualLUExecuted, flowEvaluated, operatorFlowFinite and operatorCleared; resourceQualified remains false until a separate watchdog/time-l audit supplies measured peak footprint, swap and resource verdict. Neither this dormant implementation nor the existing H8 single-direction probe certifies either H64 direction.
