# Future fresh campaign qualification across native regrids

`fresh_profile_qualification_chain(baseQualificationFile, initialDataFile,
campaignDirectories, targetResultFile, newOutputDirectory)` is an independent
no-LU evidence checker. It does not modify or run any current frozen driver.
The initial caller ledger remains the actual step139 qualification with
legacy smoke=false and prospective continuous auxiliary=true.

The checker walks the explicitly ordered campaign directories and their
numbered stages. Every source path must equal the previous accepted native
endpoint; every result is revalidated with its exact native checkpoint,
full history, clocks, scales, physical samples, captured initial data,
lineage and unchanged numerical contract. It requires each completed
segment to reach its real equivalent-tau .2 endpoint within512 steps,
with core>=20 and safety<.70. The complete recorded history retains the
original trusted/core14/safety.70 requirements. A final failed design with
no native checkpoint advancement may be left behind when the next explicit
campaign resumes the identical preceding checkpoint. An incomplete native
advancement cannot be skipped.

Native regrids require remeshCount+1, unchanged step/scale/clocks/case/base
axes/initial invariants and runtime reference values. Only the old terminal
history row may be replaced. Both exact transport anchors remain on the
axis; originIndex and pinIndex must obey the maintained nearest-node rule
in `restoreRuntimeReferences`. The tracking pin is not assumed to equal
the exact transport anchor. Full grid quality and local-q gates are
recomputed. From actual native rho/axes, the checker independently recomputes
rhoX peak jump, mass, range, production nodal x/y core and left front.
Actual core31/front20, jump.002, mass5e-12, range2e-4 and cumulative absolute
jump.02 gates must pass; measured values must agree with the saved native
transaction audit and terminal history. Same-grid segments preserve all
native regrid provenance exactly. JSON `passed` alone grants no qualification.

## Executed no-LU evidence

Under `result/longtime/20260908_campaign_v1/fresh_qualification_chain_nolu_v1`:

- `positive_v1_stage2/chain_report.json`: real two-segment chain through
  step504, equivalent tau7.643200751859, existing remeshCount1, passed.
- `positive_v2_stage1_final_v2/chain_report.json`: final helper's real
  three-segment chain through step729, equivalent tau7.843200751859,
  remeshCount2, passed including the actual step504 native regrid.
- `test_report.json`: missing numbered stage and wrong source checkpoint
  negative fixtures were both rejected at their specific chain gates.
- Final Code Analyzer returned zero issues. MATLAB exited0 without restore,
  elliptic factorization, PDE steps or native writes. Exact final source,
  SHA256 records and audit runners are retained in this directory.

The earlier extra index audit in `positive_v2_stage1_final` failed because
it incorrectly assumed pinX must itself be a grid node. That failure is
preserved. The final check follows the actual maintained nearest-node
contract while retaining the separate exact transport-anchor requirement.

A future caller can use this helper to rebuild the evidence chain before
qualifying a new batch. The current frozen driver still accepts only its
original same-grid bootstrap schema; this report does not bypass that
contract or automatically promote a remeshed checkpoint. A subsequent
driver revision must explicitly call the chain checker rather than trusting
a saved success flag. No parent-native history, independent initial-grid
convergence or infinite-time limit is claimed.
