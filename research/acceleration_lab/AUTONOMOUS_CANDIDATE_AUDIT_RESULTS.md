# Actual checkpoint and negative regression results

`ipm_accellab_test_candidate_audit` performs strict native reads of the
step3668 checkpoint before the regrid and the separately committed native
regrid checkpoint. Both original files remain unchanged. The test constructs
a clearly marked, synthetic **precommit view** using the post-regrid rho,
axes, seven-point maintained Dx, high-order tensor quadrature and rescaling
indices. Its configuration, metadata, initial invariants, scale and clocks
come from the old state, and its remesh count is the old count, 5.
This view is not a checkpoint and is never signed or written as one.

Only `readCheckpoint`, one-dimensional derivative/quadrature/quality helpers,
the candidate audit and algebraic mutations are called. There is no
`restoreCheckpoint`, `mesh.build`, flow, Poisson evaluation, transfer or PDE
step. Both MATLAB batches used the default native ten-thread environment.

## Native metric parity

The original checkpoint transaction had count 5 before and 6 after, with
unchanged step and clocks. The new precommit audit reproduces all six
registered old audit values **bitwise**, in both regression versions:

| Quantity | Old native audit and new audit |
| --- | ---: |
| Actual x core cells | 35.155590271574454 |
| Actual y core cells | 37.00690537837285 |
| Relative whole-box omega peak jump | 3.2792341418625766e-5 |
| Native relative mass defect | 0 |
| Relative range violation | 9.414698872055325e-9 |
| Actual left-front cells | 39.38194080486816 |

The identity fixture passes the unchanged registered policy. Deliberate
off-wall peak growth with the exact old wall trace is rejected by the
whole-box peak gate. A narrow analytic field with fictitious design scores
of 999 cells is rejected by the actual core and front floors.

## Preserved first failure and corrected rerun

Version one is preserved at
`result/longtime/20260908_campaign_v1/acceleration_lab/candidate_audit_regression_tp475c0737_5832_465b_a37b_af6613e5f99e`.
It passed 25 of 27 tests. Setting a single source or candidate integration
weight to NaN produced a NaN mass defect but returned `passed=true` because
the comparison `NaN > tolerance` is false. The original MAT report retains
the NaNs; JSON encodes them as null. This result remains `allPassed=false`.
The source snapshot and all original failure cases are retained.

Root then added checks for real finite fields, paired positive finite
integration weights, finite real paired gradients and finite derived
metrics. The independent second run used a complete 117-file source freeze:
`candidate_audit_v2_source_575ebfe6dca84b5e9ee2a8382fa2094f` under the same
acceleration result root. Its manifest covers the maintained package and
the test entry point. The second output is
`candidate_audit_regression_tpe4b8a291_1083_4302_b726_5609a098575d`.

All **42/42** second-run tests pass, while native metric parity remains
bitwise. The additional tests cover incorrect weight shapes, zero, negative,
complex and infinite weights on each side; nonfinite derivative matrices
on each side; complex rho; and mass overflow from finite positive weights.
All these invalid inputs are rejected. Existing tests also cover clocks,
initial mass/range, configuration, metadata, base axes, numeric gauge
references, incorrect origin/pin indices, absent exact anchor, changed box,
premature count increment, cumulative peak-budget exhaustion and unsupported
adaptive levels.

The signed-cancellation test gives native mass defect
`1.0000057917e-10` and hypothetical L1-normalized defect
`1.0004138522e-13`; the corrected native `5e-12` gate rejects it. This verifies
that restoring the old denominator prevents the specific weakening found
in the read-only review.

The checkpoint file SHA256 values before the first run and after each run
match. Both output directories contain `checkpoint_file_immutability.json`;
the second also includes a compact scalar result and a link to its frozen
source manifest. MATLAB sessions 18909 and 14899 exited normally.

The entry point accepts
`(preFile,postFile,nativeAuditFile,outputRoot)`; exact original arguments are
saved in each `registration.json`. The frozen second-run entry point can be
called again to produce a new unique output directory. Earlier audit and
scientific failures are not overwritten. These tests qualify the stated
candidate audit cases, not the complete future autonomous solve/commit
pipeline or long-time convergence.
