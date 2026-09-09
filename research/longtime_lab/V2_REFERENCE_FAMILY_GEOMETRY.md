# Version 2 reference families and pure planning

The optional policy version 2 registers an immutable finite family of reference
axes. It preserves every version 1 gate, trigger, target-derived threshold,
search value, and scaling rule. This change covers pure configuration and axis
planning only. A candidate remains subject to the native transaction, paired
state/history, field, and controller admission contracts.

```matlab
policy = ipm.config.autonomousMeshPolicy(struct( ...
    'version', 2, 'nodeFamily', struct('maximumTotalNodes', 110000)));
family = ipm.remesh.referenceAxisFamily(selectedBaseX, selectedBaseY, anchor, policy);
member = family.members(2);
[candidates, report] = ipm.remesh.plannedAxisPairs(acceptedView, ...
    struct('x', member.baseX, 'y', member.baseY), anchor, policy);
```

`version=2` and `nodeFamily.maximumTotalNodes` must be explicit. Version 1 is
still the default when the version is omitted and rejects `nodeFamily`.
The cap is a positive finite integer node budget. It is not a RAM estimate or
permission to allocate a particular LU factorization. All other family fields
are fixed, normalized values:

| Field | Registered value |
|---|---|
| `generator` | `selected_base_index_pchip_v1` |
| `cellFactors` | `[1 1; 2 1; 1 2; 2 2]` |
| `ordering` | `node_product_then_registration_index` |
| `componentwiseNondecreasing` | `true` |
| `maximumAcceptedGrowthTransitions` | `2` |

The family function returns `version=1` (the family representation version),
`generator`, `rootX`, `rootY`, `anchor`, and `members`. Axes are double row-x and
column-y; the scalar anchor must be a positive double interior x node. Roots
must have exact x reflection symmetry, origin, both anchors and at least seven
nodes per axis. Its members retain registration indices; the execution order
is a separate sort by node product, then registration index.

Each member has `index`, `cellFactors`, `nodeCount=[Nx Ny]`, `baseX`, `baseY`,
logical `resourceAdmitted`, logical `qualityPassed`, and the complete original
`ipm.mesh.quality` reports `xQuality` and `yQuality`.

Factor 1 retains the input axis bitwise. Factor 2 evaluates PCHIP on normalized
node index, explicitly restores the original odd-index knots, and fixes zero.
For x it operates only on the positive half and mirrors that half exactly.
The box, old knots and both anchors are therefore exact in every reference
member. These are immutable **reference** knots: a later redistributed
candidate need not retain every old interior node.

Every axis factor is checked with seven-point conditioning and the original
quadrature, adjacent ratio, log-spacing curvature, global and local weight
limits. This includes factors used only by resource-ineligible members. Any
quality failure rejects the entire family (`ipm:AutonomousMeshReferenceQuality`);
no partially accepted family is returned. A root exceeding the cap is rejected
with `ipm:AutonomousMeshResourceCap`. Successful reconstruction is deterministic
and suitable for strict persisted family comparison by the controller reader.

For a selected root of 321 × 161 and cap 110000:

| Member | Cell factors | Nodes | Total | Resource admitted | Execution order |
|---|---|---|---:|---|---:|
| 1 | `[1 1]` | 321 × 161 | 51681 | yes | 1 |
| 2 | `[2 1]` | 641 × 161 | 103201 | yes | 3 |
| 3 | `[1 2]` | 321 × 321 | 103041 | yes | 2 |
| 4 | `[2 2]` | 641 × 321 | 205761 | no | 4 |

`plannedAxisPairs` continues to extract features from the actual paired source
`rho/x/y/Dx/source`. Version 2 uses the target reference node counts for axis
construction. It never upsamples the source to manufacture a higher-resolution
history. The target must have the same box and fit the explicit cap. Registry
membership, componentwise monotonicity, growth accounting and the selected
family index belong to the runtime transaction/controller contract; the pure
planner does not infer them from node counts.

The existing candidate fields are unchanged: `x`, `y`, `unchanged`, `xIndex`,
`yIndex`, `predictedCells`, `quality`. There are at most three candidates per
call. A qualified old-axis keep may be included only when source and target
have the same node counts. A different-N request can return no candidate;
it must never return a source-axis keep labeled as a larger target. Adding
nodes is not a general guarantee of feasible geometry.

The three-proposal limit includes a qualified keep. Evolution excludes that
keep after planning, so the current member can have only two remaining actual
transfer attempts. This preserves the version-one top-three semantics; it is
not a promise of three transfers for every member. Each eligible larger
member receives its own independent limit, and a smaller member cannot consume
that later member's proposal budget.

## Regression evidence and scope

The independent test is `mesh_test_node_family_v2`. Its frozen version-one
comparison functions differ from their originals only in entry-point names.
Exception stacks are omitted when comparing reports because test wrapper
filenames and call-site line numbers necessarily differ; all numerical arrays,
candidate coordinates, quality reports, rankings, rejection reasons and policy
fields are compared exactly.

The actual source is the strictly read and result-paired from-zero stress
checkpoint at step 1134, canonical time 2.6886192715883674. The four reference
members and all axis trials/candidates are compared to the completed
`node_ladder_from_zero_1134_v2/run` research experiment. This tests a true
321 × 161 source against larger reference counts without LU or PDE work.
The former experiment sorted cases by node product, so research levels
`[1,3,2,4]` correspond to registration members `[1,2,3,4]`.

Version-one planning is also compared on the three existing independent
snapshots: original t=0, full-box step 3326, and fresh step 3668. Additional tests
cover exact knot retention, both unit scalings `.37` and `4`, full-family
reconstruction, resource masks, normalized configuration, invalid policy and
axis rejection, and the same-N keep exclusion for different-N targets.

Frozen attempts are preserved under
`result/verification/autonomous_runtime_20260909/node_family_geometry_v2_v*`.
Attempt v1 stopped at a Code Analyzer formatting issue in the test; attempt v2
completed actual and legacy comparisons, then exposed an overly strong test
assumption that a broad analytic fixture must always admit a larger grid.
The production behavior was correct: that different-N request returned no
candidate. Attempt v3 corrects the test to require exclusion of the old keep
and leaves feasibility to the actual source positive cases.

Attempt v3 actually completed with MATLAB exit 0: 10 configuration positives,
26 configuration negatives, 6 geometry positives, 7 geometry negatives, and
8 actual-state comparisons passed. All four maintained files had zero Code
Analyzer findings. The strict step-1134 read/pairing, full immutable input hashes,
frozen source hashes, and maintained-to-tested source equality also passed;
see `node_family_geometry_v2_v3/run/report.json` and `final_audit.json`.
The pure production closure has no research dependency, and the recorded
profile contains no solver, flow, mesh build, or time advancement call.

These are geometry/configuration results. They do not certify a changed-N
native transaction, dynamically matched field convergence, an LU memory bound,
or completion of the from-zero long-time acceptance target.
