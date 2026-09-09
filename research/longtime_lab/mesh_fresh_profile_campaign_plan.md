# Independent fresh-box profile campaign

`run_fresh_profile_campaign(bootstrapResultFile, initialDataFile, outputDirectory,
maxStages, qualificationFile)` is an audited research driver using only
`ipm.solve`. It has not been executed as a campaign. Its native preflight and
decision calculation have been tested on actual completed data without LU.

The selected bootstrap is `fresh_core_patch_long_tau724_v1/result.mat` under
`result/longtime/20260908_campaign_v1`, paired with fresh native step139.
Its actual initial samples remain
`result/verification/performance_lab_20260908/late_box_v4_stage002_epoch_v1/dynamic/crop_natural/initial_physical_data.mat`.
The fresh case is 895x386. Its local physical zero is parent absolute
epoch1.9545271366857122; no parent history is inherited. Reported labels use
`absolutePhysicalTime = epoch + localPhysicalTime` and
`parentEquivalentTau = parentCanonicalTime + lambda0 * localCanonicalTime`.
Parent-unit rates divide fresh rates by lambda0; equivalent tau is a
coordinate label, not evidence of an evolved full-box parent trajectory.

## Registered progression and acceptance

- Each stage advances equivalent tau .2, or fresh canonical time
  .006210141244493645 for this lineage. Maximum additional steps512;
  CFL, all timestep numerics, physics, transport, scaling and elliptic
  choices are preserved. Periodic CP spacing is equivalent tau .1.
- Regrid when either actual core is below26, or a log-core fit to the
  latest .35 equivalent tau on the same remesh count predicts a next-stage
  core below22. An unavailable fit is explicitly marked. Stage endpoints
  require cores20/20, safety<.70 and actual arrival at the canonical target;
  hitting maximum steps is not success. The full recorded history must
  remain trusted with cores>=14 and safety<.70.
- Direct core-patch candidates retain targets32/32 and front20. The
  registered order is fine64/fraction.65, fine64/.5, fine96/.5, fine96/.65.
  The immutable reference is the real fresh initial axes. The computational
  box and both node counts remain fixed; exact positive/negative anchor
  nodes remain required. No failed anchor-centered x family is constructed.
- Actual native transfer requires cores31/31, front>=20, peak jump<=.002,
  cumulative absolute peak jumps<=.02, mass defect<=5e-12 and relative range
  violation<=2e-4. Full axis gates remain ratio1.08, curvature.01,
  stencil rcond1e-9, positive global q ratio1e-8 and local q/control
  [.35,1.65]. Every post-transfer gate precedes CP commit. All rejected
  candidates are retained; exhaustion stops with the source CP available.
- Full native/result rho, axes, step, clocks, history, scale, physical
  samples, snapshots and configuration pair exactly. Captured anonymous
  initial-condition arrays are checked against actual saved samples and
  native parent reconstruction; complete `functions` data is compared
  across MAT reloads. Same-grid history prefixes and runtime references
  are retained; regridding may replace only the terminal history record
  and update the two grid-index reference fields.

## Caller qualification and read-only entry points

The caller must provide a JSON ledger with kind
`fresh_profile_campaign_bootstrap_qualification_v1` and fields:

```
bootstrapResultFile
bootstrapCheckpointFile
qualifiedGridCheckpointFile
legacySmokeReportFile
legacySmokePassed: false
continuousAuxiliaryVerdictFile
continuousObserverAuxiliaryPassed: true
```

The qualified grid is the actual fine step30 CP from
`fresh_epoch3326_core_patch_smoke_v1/new_grid`. Step139 must preserve its
same-grid history prefix, lineage, runtime references and numerical data.
The ledger points to the unchanged smoke report and `auxiliary_verdict_v1.json`.
Every native/linear-field/peak/rate/continuous-width auxiliary subgate must
be true; the original legacy smoke must remain false. This does not claim
generic fourth-order geometric convergence: the observer's stricter
per-phase refinement-ratio test remains false. Test fixture ledgers are
explicitly rejected by the actual campaign driver.

With the final frozen solver/research paths active and exactly10 threads,
preflight can be run without calling the driver:

```matlab
[review, cp, contract] = fresh_profile_campaign_review(bootstrap, initial);
qualification = fresh_profile_campaign_qualification(ledger, bootstrap, cp);
```

The coordinator's intended batch call after review is:

```matlab
campaign = run_fresh_profile_campaign(bootstrap, initial, newOutput, 12, ledger);
```

Twelve completed stages would reach equivalent tau9.643200751859 from this
bootstrap. A batch boundary requests review and never marks the asymptotic
scientific goal complete.

## Actual no-LU verification, 2026-09-08

`result/longtime/20260908_campaign_v1/fresh_profile_driver_nolu_v3/report.json`
passed actual step19/30/139 native/paired/lineage/physics-capture/axis checks,
history-prefix validation and prediction decisions. Step19 requests regrid;
fine step30 and step139 do not. Step139 has equivalent tau7.243200751859,
absolute physical time1.9590400984028897, core32.760389/36.316631 and
safety.244197. Next-stage predicted cores are28.961167/33.840324.
False legacy claims, altered references, physics and history were rejected.

`final_qualification_checks.json` in the same directory separately passed
the actual qualification and rejected each of five individually falsified
auxiliary subgates. Final Code Analyzer reported zero issues in all five new
files. Both MATLAB processes exited0 and built no elliptic operators; no
campaign PDE or regrid was run in this implementation task. Logs are
`result/longtime/20260908_fresh_profile_driver_nolu_v3.log` and
`result/longtime/20260908_fresh_profile_qualification_final.log`.
The directory retains scripts, exact reviewed source and SHA256 records.

Freeze the five new `run_fresh_profile_campaign` / `fresh_profile_*` files
with the coordinator's established schema4 solver and grid_lab tree, plus
`mesh_design_fresh_core_patch.m`, `mesh_fresh_core_patch_candidate.m`,
`mesh_core_patch_axis.m` and `mesh_general_rounded_axis.m`.
