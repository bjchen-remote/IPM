# Registered variable-node transfer and candidate audit

The version-2 branches of `ipm.remesh.transfer` and
`ipm.remesh.auditCandidate` passed **50/50 bounded no-LU tests** on 2026-09-09.
MATLAB session 65856 exited 0. This is interface and array-state validation;
**no successful transfer, Poisson build, PDE step, or native transaction was
performed by this suite**. The actual transaction remains a separate test.

## Transfer interface

The existing transfer function and interpolation/conservation path are retained:

```matlab
proposal = struct('x',candidateX,'y',candidateY, ...
    'referenceFamily',memory.referenceFamily, ...
    'sourceLevelId',memory.currentLevelId, ...
    'targetLevelId',candidate.nodeFamilyIndex);
[rhoNew,opsNew,metrics] = ipm.remesh.transfer(rho,ops,config,proposal);
```

For an enabled version-2 policy, validation precedes interpolation and build.
The function reproduces the complete family with
`referenceAxisFamily(rootX,rootY,anchor,policy)` and compares every stored
member exactly. Roots must match the frozen configuration's original node
counts, box, and anchor. The actual source counts, base axes, field shape,
and endpoints must match the named source member. The target must be admitted
by the original resource and quality gates and cannot decrease either cell
factor. Its actual axes must match the target counts, preserve the box,
symmetry, and exact anchors, and remain finite and strictly increasing.

Only `gridOverride.nx/ny/customX/customY` changes for the numerical build.
`config` retains its original node counts. The result adopts the selected
member's base axes. Runtime scaling references are copied from the source;
only origin and pin indices are recomputed, as on the original path.
`remeshCount` is unchanged until the caller commits. The transfer function
does not own or fabricate a transaction ledger.

## Precommit audit

`auditCandidate(source,candidate,policy,priorAbsolutePeakJump)` reads the
version-2 registration and current level from **source** metadata. Candidate
metadata must still be exactly equal to the source metadata. The target is
identified uniquely from its base axes and node counts; a candidate does not
get to assert its own level in metadata before acceptance.

The version-2 audit adds:

```text
version = 2
sourceLevelId, targetLevelId
sourceNodeCount, targetNodeCount
sourceCellFactors, targetCellFactors
sameBox
registeredNodeCounts
matchingReferenceMembers
admittedLevelTransition
```

`sameBoxAndNodeCount` retains its literal meaning and is false for a genuine
change in node count. The four version-2 registration gates must still be
true. Registered base axes must reproduce from the original zero-time
selected axes, including the original node counts and box.

The audit reconstructs the level chain and accepted growth count from the
source transaction ledger. Each transition is monotone in both cell factors;
each source/target is admitted; at most two accepted growth transitions are
allowed. The ledger length must match the precommit remesh count. The
original `sum([transactions.relativePeakJump])` grouping must reproduce the
stored cumulative budget and the passed prior budget. No separate candidate
growth counter is trusted.

Clock, scale, config, metadata, initial mass/range, runtime-reference, quality,
whole-field peak, actual core/front, mass-defect, and range gates remain in
place. A node-count change does not relax any numerical threshold. Version 1
retains its original audit structure and same-box/same-node-count path.

## Executed evidence

The executable is `ipm_accellab_test_variable_n_audit.m`. It takes an existing
configuration fixture and an unused output directory. All tested states and
the small ledger are explicitly **synthetic API-contract fixtures**, not
checkpoints or a replayed physical trajectory. The passing polynomial is

```text
rho(X,Y) = (X^2 - X^4/32) (1 - Y^2/16),
box = [-4,4] x [0,4], initial N = 321 x 161, anchor = 1.
```

The registered cap is 110000 nodes. The same original policy admits level 1
and the individual X/Y refinements, and rejects the combined refinement.
Candidate fields were evaluated analytically, using paired pure derivative
matrices and quadrature; the successful transfer function was not called.

| Audit fixture | Actual core X / Y | Peak jump | Mass defect | Result |
|---|---:|---:|---:|---|
| 321 x 161 identity | 47.92491 / 50.59406 | 0 | 0 | Pass |
| 641 x 161 analytic samples | 95.85363 / 50.59406 | 2.21216e-5 | 1.96272e-15 | Pass |
| 321 x 321 analytic samples | 47.92491 / 101.19212 | 0 | 1.78429e-15 | Pass |
| 641 x 321 | — | — | — | Original cap rejects |

The suite also verifies same-level operation after a synthetic prior growth;
coarsening and cross-family shrinkage; altered configuration, clocks, initial
references, metadata, base axes, level identifiers, family members, budget,
ledger, weights, derivative data, box, anchor, and complex field rejection.
Sixteen invalid transfer contracts are rejected before interpolation/build.
The v1 audit is compared against a frozen pre-change function on identity and
nonzero amplitude-defect fixtures; the complete outputs are `isequaln`.
Profiler evidence confirms no build, flow, advance, restore, or remesh
interpolation call. The input fixture remains unchanged.

The first two trial fixtures used `X^4/16` and `X^4/24` instead of `X^4/32`.
They left an unwanted positive wall-gradient competitor at the negative box
boundary. The original global positive-peak/core diagnostic rejected their
small connected cores (3.28170 and 1.80615 on the root grid). Family and mass
gates had passed. Their failed runs and the explicit first-fixture diagnosis
remain saved; no production threshold was changed to make them pass.

## Artifacts and limitations

All paths below are relative to
`result/longtime/20260908_campaign_v1/acceleration_lab/`:

- Final report: `variable_n_audit_ddc207b7ee144aabbce8c58b5283f3f9/report.json`
  and `ready_source.json`.
- Frozen final source:
  `variable_n_audit_source_ddc207b7ee144aabbce8c58b5283f3f9/`.
- Preserved failures: `variable_n_audit_c6462cdddbbc4803993ed4f9f280e2fa/`
  and `variable_n_audit_5522ba3914e54f59a7c414f939f60b19/`.
- Original fixture diagnosis:
  `variable_n_initial_fixture_diagnosis_5522ba3914e54f59a7c414f939f60b19/`.
- Pre-change production copies:
  `variable_n_transfer_baseline_bbec274d9e554fdf920769a43eee34ea/`.

The two tested production files exactly match the final frozen copies:

```text
transfer.m       ba59ef0e4a59141bfeff0d9c202f7945d5245537258ca8597dbd0c6ce35b503a
auditCandidate.m 62255b1551886dcdbeab904575d197aa6fd5c082bb0f8f66845e24315426a189
```

The remaining actual test must use the original transfer, fresh flow/RHS,
native acceptance and persistence contracts, and a strict restart. The
no-LU results above do not establish transfer accuracy, covariance of the
physical trajectory, lifetime of new factors, or bitwise native continuation.
