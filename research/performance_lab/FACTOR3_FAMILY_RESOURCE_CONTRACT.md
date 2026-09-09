# Future five-member node family and single-factor resource qualification

This is a research preparation, not a native policy or a checkpoint upgrade. The proposed future policy version is3 and reference-family version2, with generator `selected_base_index_pchip_integer_v2`. It must be explicitly registered before an original t=0 solve. Current version1/version2 configs, old checkpoints, their interpretation, and all production +ipm files remain unchanged. The existing production normalizer was actually tested and rejects the future draft.

## Exact finite contract

The registered members are the original four `[1,1],[2,1],[1,2],[2,2]` plus only `[3,2]`. For original321×161 they contain51681,103201,103041,205761,308481 nodes; the final member is961×321. Node cap310000 and maximum accepted growth transitions3 are explicit. Members retain registration IDs1..5; node-product then registration-index order is1,3,2,4,5. Transitions are componentwise nondecreasing: current4 may propose4/5, and current5 may propose only5. Per-level candidate cap3 gives at most15 selected pairs across the entire family. Growth count is recomputed from the accepted ledger, never trusted from a proposal.

Each reference axis comes directly from the original accepted initialization selectedBase, using normalized-index PCHIP. Original knots are written back exactly; X is built on the positive half then reflected. Endpoints,0 and±anchor remain exact. The actual original four members—including complete quality reports, indices, factors, node counts and resource flags—match the previous native family bitwise. The factor3 X reference also matches the already-saved tau8 research reference bitwise. It does not contain all factor2 intermediate knots, and actual regridded PDE axes are not nested.

Only future-policy version, reference generator/version, one extra factor row, cap210000→310000, and growth budget2→3 change. Every other numerical, geometry, quality, target, trigger, peak, mass, range and safety parameter is unchanged. A future implementation must independently support the five-level ledger, exact immutable family regeneration, history/window/snapshot node-count pairing, rollback/one-factor lifecycle, and native signed checkpoint restore. This preparation does not implement those persistence/runtime boundaries. In particular, the original config.nx/ny must remain the initial321/161; actual N can only come from a registered member. No policy is inferred or upgraded from an older checkpoint.

## Actual frozen input and noLU results

The input is the original t0 run's trusted step4545/tau8 checkpoint, at physical time2.2412067821954538, original box[-8,8]×[0,4]. This is a late operator probe, not a new t0 trajectory. The saved previous factor3 first candidate was reused; there was NO repeated70×4 search. Exactly one pair was reconstructed from its recorded parameters and matched bitwise: X fine32, rounding16, fraction.5, warp.25, coreCenter.9781786284671806, spacing9.493485383140184e-5; Y sigma.5.

A six-point constant-conservative X-then-Y interpolation, using the original maintained interpolation helper, supplies the operator's frozen density. There is no amplitude normalization or projection. This is deliberately labeled an offline interpolated field; it differs from the native transfer's wall-trace-preserving Y bubble correction and is not native migration evidence. The original frozen-pair gates all pass. Interpolated core is35.202757926/37.489428886, front36.429939084, peak change9.262629088e-05, conservation defect3.940276354e-15, range3.688012384e-09. All full diagnostics are retained.

The family, one-candidate reconstruction and field preparation took1.139276s. Three MATLAB files had zero CodeAnalyzer findings; default allowLarge=false probe passed. The profile excludes restore, LU, flow, native transfer, PDE solve and any full planner search. 123 frozen files and four original inputs retained their SHA hashes. Nine pure Python watchdog tests passed, and the parser was additionally checked against the actual old tau8 time -l log.

## Dormant resource measurement

`run_operator_dormant.m` remains unexecuted. It strictly reads the original CP once, builds only one actual961×321 operator/factor with original `ipm.mesh.build`, and evaluates exactly one original frozen `ipm.evolve.flow`. It does not restore an old factor, call native transfer/solve, advance a clock, write a checkpoint or manufacture an accepted state. Source config remains exact; the explicitly experimental gridOverride changes only nx/ny/custom axes for the low-level operator. Original rescaling references and scale are retained, with only origin/pin indices paired to the new axes. LU is cleared before field-output serialization. Original warnings are retained; the saved Poisson residual is not a forward-accuracy certificate.

Use only the frozen watchdog after root grants a later slot:

```sh
python3 '/Users/bjchen/Library/CloudStorage/OneDrive-北京大学/Research/IPM/ipm_structured/result/verification/autonomous_runtime_20260909/factor3_family_resource_preparation_v1/source/ipm_perflab_factor3_watchdog.py' --plan '/Users/bjchen/Library/CloudStorage/OneDrive-北京大学/Research/IPM/ipm_structured/result/verification/autonomous_runtime_20260909/factor3_family_resource_preparation_v1/resource_plan.json' --context /absolute/path/to/root_granted_resource_context.json
```

`root_resource_context.template.json` has rootGrantedWindow=false and no main reservation, so it refuses to launch. The caller must provide an actual completed time -l peak log and an explicitly reserved current-main peak larger than that historical observation. Historical peak is not asserted to bound future memory. The gate requires main reservation + own watched6GiB + system3GiB ≤ physical RAM. The own process group is sampled every.25s and terminated on observed RSS>6GiB, wall time>300s, or an observed increase in displayed system swap. Only that newly spawned process group can be terminated. Missing monitor data or missing final memory metrics fail closed.

Qualification additionally requires exit0, finite operator flow, actual peak footprint≤5GiB, zero child swaps, no observed displayed swap increase, and unchanged source/input SHA. The full time -l output and raw flow are retained. RSS-group monitoring and footprint are distinct metrics; discrete sampling is not a byte-perfect hard memory guarantee. The source/source-family LU peak has NOT been measured at961×321. The known205761-node run had3.821654360GB peak footprint and zero swaps; no linear-in-N extrapolation is used for admission.

Even a successful single frozen matrix probe will only qualify this bounded resource case. It will not certify all future pivot patterns, H16/H32/H64 families, native variable-N transfer, restart, PDE time accuracy, box convergence or long-time stability. If this case exceeds the fixed resource gate, its failure is retained; the proposed family is not silently admitted.
