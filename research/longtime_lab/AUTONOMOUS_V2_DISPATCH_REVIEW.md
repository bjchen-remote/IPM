# V2 runtime dispatch review

The reviewed production files are `evolve.planAutonomousMesh`,
`evolve.applyAutonomousMesh`, `remesh.controllerTelemetry`, `evolve.advance`,
`solve`, and the initial caller `evolve.initialize`. This was a read-only
review plus pure planning tests. No production file was edited and no
`applyAutonomousMesh`, flow, time advancement or LU was executed by this review.

No actionable defect was found in this revision's integration. The frozen
source and maintained source were still byte-identical after the checks;
the hashes are in
`result/verification/autonomous_runtime_20260909/v2_plan_dispatch_review_v1`.

The code preserves these concrete boundaries:

- Eligible levels are componentwise nondecreasing, resource admitted and
  geometrically admitted. Sorting by node product and registration index
  places the current level first. Each level calls the pure planner with its
  own three-proposal budget; exhausted small levels do not consume larger
  levels' budgets. The budget includes keep, preserving v1 semantics.
- A different-N reference cannot produce an old-axis keep. Evolution removes
  qualified same-N keeps before attempting actual transactions. Each remaining
  candidate retains its own target family index.
- Initial family registration uses the actually accepted analytic zero-time
  axes, after the initial field/flow/resolution checks. Family construction
  failure rejects initialization. It does not install a partial family or
  create a fresh physical epoch.
- `advance` only returns a pure plan. Its old accepted-state alias is gone
  before `solve` removes the source `poisson` field. `applyAutonomousMesh`
  asserts that removal and clears each rejected trial before another candidate
  is built. On complete failure `solve` rebuilds the original source factor
  from its exact stored matrix. Initial selection also clears the local `ops`
  alias before releasing the source factor. This lifetime conclusion is static;
  the review does not claim a measured memory bound.
- A v2 application checks plan/source step, clocks, remesh epoch, node family,
  dimensions and flow core observations before trying candidates. It submits
  a new level and remesh count only after transfer audit and new-flow gates.
  The audit retains the pre-transaction decision and ordered attempt summary.
  A planning stop rolls back to the previous accepted state and preserves its
  failed future decision separately; a candidate rejection preserves the
  already accepted source step and its failure/axis evidence.
- V2-only additions remain guarded. The observed v1 state and plan below are
  unchanged. The finite-time telemetry window still has no fixed sample-count
  bound; that previously documented limitation remains unresolved here.

## Executed no-LU evidence

`mesh_test_v2_plan_dispatch` completed with MATLAB exit 0. The actual step-1134
checkpoint was strictly read with ten threads. Its observed state and pure
plan were exactly equal to the frozen v1 implementation. Four separate,
explicitly synthetic controller inputs tested source levels 1–4. They are
not native solver trajectories and do not claim physical consistency of their
assigned trigger counts.

| Source level | Eligible order | Proposed counts per level |
|---|---|---|
| 1 | `[1 3 2 4]` | `[0 0 0 3]` |
| 2 | `[2 4]` | `[0 3]` |
| 3 | `[3 4]` | `[0 3]` |
| 4 | `[4]` | `[3]` |

All returned candidates matched their member dimensions, were non-keep
proposals, and preserved order. The source rho, operators and clocks were
unchanged by planning; `lastDecision` exactly matched the returned pure plan
projection. The results directly cover the case in which every smaller level
fails while a larger level still receives its full proposal budget.

`run/report.json` and `final_audit.json` in the review directory contain the
actual results and immutable input/source checks. Real transaction, rollback,
checkpoint roundtrip and natural from-zero growth qualification remain the
separate runtime integration tests; no claim here substitutes for them.
