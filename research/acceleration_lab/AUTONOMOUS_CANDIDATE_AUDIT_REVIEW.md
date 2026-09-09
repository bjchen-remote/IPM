# Read-only review of the autonomous candidate audit

Reviewed on 2026-09-09 before production integration. No solver file was
changed, no MATLAB process was started, and no transfer or LU was built.
The references are `+ipm/+remesh/auditCandidate.m`, `transfer.m`, `adapt.m`,
`research/grid_lab/ipm_gridlab_regrid_checkpoint.m`, the current fresh
campaign transaction, and the new autonomous configuration scope. Line
numbers below refer to the version read during this review.

## Concrete discrepancies

1. **Mass normalization changes the old gate.** `auditCandidate:61–64`
   uses `max(abs(oldMass),integral(abs(rho)))`, with a `realmin` floor.
   The native gridlab transaction's `audit_transaction` uses
   `max(abs(oldMass),eps)`. With signed cancellation the new denominator
   can be arbitrarily larger. For example, old mass `.01`, weighted L1
   mass `1`, and mass defect `1e-13` give `1e-11` under the original
   formula (fails `5e-12`) but `1e-13` under the new formula (passes).
   Keep the original `massRelativeDefect` gate; an L1-normalized defect can
   be reported separately. No defect is asserted for the actual
   sign-definite cached branch, where the normalizers nearly coincide.

2. **The fixed core level requires configuration admission.** The new
   audit and `meshFeatureIntervals` explicitly use levels `[.1,.5,.9]`.
   Maintained `transfer` computes its terminal core using the last entry
   of `ops.rescaling.adaptiveLevels`, as do native flow diagnostics.
   At review time the enabled autonomous scope did not restrict that
   otherwise user-configurable vector. A valid custom vector ending in
   `.8` would therefore give different meanings to audit and evolution
   core counts. The previous q512 driver explicitly rejects vectors
   other than `[.1,.5,.9]` in its `validate_q512_config` scope.
   The version-one policy should retain that restriction.

3. **Complete transaction provenance needs a caller or audit check.**
   The new function compares `scale`, step, normalized time, the entire
   rescaling struct except rebuilt indices, and the base axes. These
   checks are sound. It does not inspect `mass0`, `rhoRange0`,
   `runMetadata.caseId` / fresh lineage, or non-grid configuration.
   Changing only those fields in an otherwise passing complete state
   leaves its answer unchanged. The original checkpoint transaction
   preserves them by copying the full state, and the fresh campaign
   additionally checks lineage and clocks after the transaction. Keep
   those checks in the new caller, or include them in the complete-state
   audit. Derived candidate flow need not equal old flow after a transfer;
   it must instead be recomputed before committing history/checkpoint data.

## Checks that agree with the maintained transaction

- The jump uses the whole-field `max(abs(rho*Dx'))`, independently computed
  on each grid. It is not reduced to a wall peak or a C1 selector.
  The source must have a positive finite peak. The new denominator omits
  the native `eps` floor, which is stricter for extremely tiny amplitudes;
  it does not weaken the actual late-state gate.
- Range violation and its normalization agree with the original gridlab
  transaction. Actual x/y core counts reproduce `transfer` when
  `adaptiveLevels=[.1,.5,.9]`. The y count is the nodal vertical column at
  the wall's selected peak, not a continuous y-width claim.
- `meshFeatureIntervals` uses the same five-point smoothing, seven-point
  maintained derivative, positive-half feature selection, half-height
  front width, and fractional cell counts as the gridlab feature helper
  in the admitted high-order scope. The actual front floor is additional
  to the old generic checkpoint gate and must remain distinct from
  padded design-interval counts.
- Recomputed quality uses the admitted high-order quadrature and
  seven-point stencil. Limits `1.08/.01/1e-9/1e-8` match the recent fresh
  campaign, with extra local quadrature/control-width bounds `.35..1.65`.
  They are intentionally not the older generic gridlab defaults
  `.03` for spacing curvature and `1e-4` for the global quadrature ratio.
- Candidate `remeshCount` must equal the old count before acceptance:
  maintained `transfer` preserves it; `adapt` increments only after
  acceptance. The old gridlab checkpoint wrapper internally increments
  a copied state before its final audit, but it exposes/writes that copy
  only after all gates pass. The new pre-commit API should not inherit
  that wrapper's intermediate `+1` convention.
- Exact paired axes, endpoint identity, base axes, zero y origin, and
  both exact transport anchors are stricter than the generic transaction's
  tolerance-based axis checks and consistent with this registered policy.
  Rebuilding only origin/pin indices is allowed; all numerical reference
  values must stay fixed.

## Bounded regressions suggested before integration

Use a valid prepared candidate-state fixture without a Poisson object;
the audit only needs arrays. Retain the whole source and verify it remains
unchanged after each call. The basic passing fixture must meet the actual
registered core/front/quality gates, rather than lowering thresholds for
the test.

1. Passing zero-change state and actual transferred candidate: independently
   recompute the original mass, full-field peak and range formulas and
   require the matching audit values.
2. A signed-cancellation density mass example between the two formulas'
   acceptance boundaries: the native formula must reject it.
3. Mutation of each clock, `mass0`, initial range, base axis, numeric gauge
   reference, or lineage: reject, or document the dedicated caller check
   responsible for that mutation.
4. A changed runtime origin/pin index without the proper nearest node,
   missing exact anchor, changed endpoint, and premature remesh-count
   increment: each must reject. A commit increments once; a rejected
   attempt does not increment or consume the cumulative peak budget.
5. Off-wall derivative-peak growth with unchanged wall trace: the global
   peak gate must reject even when wall core counts remain adequate.
6. A padded proposed interval meeting the design count while the actual
   transferred `.9` core or half-height front misses its floor: reject
   using actual metrics.
7. Enabled autonomous configuration with nonstandard adaptive levels:
   reject before building operators.
8. Ensure nonfinite mass/quality/geometry measurements cannot pass through
   comparison semantics. The normal maintained builder already supplies
   finite positive weights, but a complete-state audit should fail closed
   if an invalid array fixture reaches it.

These are transaction checks only. Passing them does not establish
long-time profile convergence, validate a new C rule, or remove any
earlier spatial or acceleration failure.
