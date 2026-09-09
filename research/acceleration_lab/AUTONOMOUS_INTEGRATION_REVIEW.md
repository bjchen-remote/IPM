# Maintained autonomous mesh integration: bounded independent review

The corrected controller and planner passed **14/14 no-LU regression cases** on
2026-09-09. This review does not qualify a physical trajectory, a native mesh
transaction, or the lifetime/bitwise action of a rebuilt Poisson factor.

The executable is `ipm_accellab_test_controller_plan(fixtureFile,outputRoot)` in
this directory. The completed run used the immutable source copy
`result/longtime/20260908_campaign_v1/acceleration_lab/controller_plan_source_2132bca282a5491d9dab305aa3385444`
and wrote
`result/longtime/20260908_campaign_v1/acceleration_lab/controller_plan_regression_tpd3f5ca5d_39c1_4a3f_87fc_d5c95998d506`.
The latter contains registration, JSON/MAT reports, and a source-manifest link.
MATLAB session 31666 exited successfully. No production files were changed by
this review or its regression.

## Defect found and corrected

The earlier planner could request an evolved remesh from a rapid forecast even
when the current core count exceeded the target. An unchanged qualifying pair
could then be ranked first and submitted, consuming a transaction and resetting
the trend without changing either axis. The maintained implementation now has
three independent guards:

1. A forecast requests work only if **the same axis** has predicted count below
   the buffer and current count below the target.
2. Evolved plans discard candidates marked unchanged.
3. `applyAutonomousMesh` rejects an evolved candidate unless at least one actual
   coordinate array changes, even if its `unchanged` flag falsely says otherwise.

The regression uses the saved `probe_v4/run/smooth_broad_fixture.mat` initial
samples and selected axes. It verifies that those samples exactly equal the
configured analytic initial datum and constructs only a pure seven-point `Dx`.
All core counts, steps, and clock increments used below are **synthetic unit
fixtures**, not measurements from that trajectory. The fictitious remesh ledger
is also explicitly labelled synthetic.

| Fixture, with consecutive clock increment 1e-4 | Result |
|---|---|
| Both axes 42, 41, 40, 39 | No request, although forecast is below buffer |
| Both axes 34, 33, 32, 31 | Forecast requests; one unchanged pair excluded, two genuinely changed candidates retained |
| X follows 42, 41, 40, 39; Y remains 31 | No request: forecast and deficit are on different axes |
| Same actual axes, `unchanged=true` | Rejected before build with `ipm:AutonomousMeshUnchangedTransaction` |
| Same actual axes, `unchanged=false` | Rejected before build with the same identifier |
| Identical time and accepted step | Replaces the existing final telemetry row; window length remains four |
| Same time but a different accepted step | Rejected |
| Backward time | Rejected |
| Counter change without a transaction ledger | Rejected |
| Changed policy | Rejected |
| Budget not reconstructible from the ledger | Rejected |
| Same-time synthetic new remesh epoch | Retains only the new epoch's row; old trend is cleared |
| MATLAB profiler | No mesh build, flow, transfer, advance, or checkpoint restore calls |
| Source data | Initial field, axes, and remesh count remain unchanged |

The apply-negative fixture intentionally lacks fields required by transfer.
Thus removing the early coordinate guard cannot accidentally start an expensive
build during this test. The expected error identifier is also checked, so a
later unrelated error cannot count as a pass.

## Other review conclusions

`advance` returns the candidate plan as its third output. Its local accepted
state goes out of scope before `solve` removes the source Poisson factor and
calls `applyAutonomousMesh`. The plan and RHS cache contain arrays rather than
factor objects. Initialization similarly clears its local operator alias before
releasing the source factor. These code paths support the intended lifetime
discipline; this test does not measure process memory or decomposition handles.

Each candidate starts from the original source state. Failed trial data are
cleared before the next candidate. The accepted source is assigned only after
the candidate audit, fresh flow finiteness/resolution gates, transaction-memory
update, and RHS cache construction succeed. No additional definite failed-step
commit defect was identified in the reviewed source.

The two failure boundaries intentionally differ. Planning exceptions, hard-floor
violations, and a plan with a stop reason roll back to the pre-RK accepted state.
If a valid plan's actual transfers all fail, the caller retains the already
completed safe post-RK state and stops, rebuilding that source factor from its
stored matrix. Neither path may be described merely as “all failures roll back
one step.”

The zero-time mesh selection re-evaluates the analytic initial density and
establishes initial mass, range, and gauge references on the selected grid.
Its selected base axes then become the immutable reference axes. Later transfers
preserve those references, initial mass/range, config, clocks, and metadata and
are audited before incrementing the remesh count. Zero-time resampling is a
different operation from migration and does not claim inherited mass transfer.

The initialized metadata now records `operatorGridRepresentation=custom_axes`,
and `gaugeContract` is recomputed after selection to record the actual anchor
node. This matters because equal coordinates do not imply equal mapped metrics:
an analytically generated axis and a supplied custom axis can use different
metric construction paths. Restart must replay the saved representation, not
infer it merely from coordinate equality.

## Remaining bounded validation

These are recommended tests, not newly demonstrated failures:

- Force a pre-RK planning failure and a post-RK all-candidate rejection
  separately; compare the retained field, complete clocks, accepted step,
  history, transaction ledger, and RHS cache to the appropriate source state.
- For a tiny rejected transfer, rebuild the original factor and compare its
  solve action, complete `rhoRate`/`flow`, `A`, metrics, quadrature, and cache at
  the same thread count. Static inspection only establishes the same stored
  matrix and LU kind, not bitwise equality of the rebuilt action.
- For an enabled zero-time selection, test restart equality of the custom
  metrics as well as axes, references, flow, and RHS. Keep disabled legacy
  restart parity as a separate test.
- Monitor actual repeated near-target transactions separately from this unit
  regression. The no-op defect is fixed; this test does not assert that every
  genuinely different short-lived mesh is economically useful.

The external research v4 long run uses its separately frozen numerical source
and driver. These maintained-integration review results must not be confused
with its independently observed native transaction qualification.
