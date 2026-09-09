# Same-policy paired proposal after the step4768 family stop

The fourth stage of `fresh_adaptive_campaign_v3` safely stopped after all
four registered legacy axis families failed their frozen mesh gate. The
trusted last native state is step4768, parent-equivalent tau
`10.643200751859164`, absolute physical time `2.0131543687297113`, with
physical gradient maximum `87.97334482777114`.

The independent wrapper `ipm_accellab_design_after_family_failure` read the
four original candidate files and retained their failed decisions. It then
strictly validated the native checkpoint, its entire result/history pairing,
the fresh lineage and original parent-crop samples at ten threads. No LU,
native transfer, flow evaluation, or PDE step was performed.

## Why the four fixed families failed

All four fail the x adjacent-cell-ratio limit `1.08`, with ratios
`1.082241162`, `1.081990340`, `1.092175761`, and `1.091579570`. All other
recorded x/y quality gates pass. The paired scorer rejects these axes before
performing its complete field-transfer score: the false resolution/transfer
flags are therefore untested downstream gates, not evidence of separate
measured transfer failures.

Their rounded branch log-spacing parameters reach `.07818` for 64 fine
cells or `.08721` for 96 fine cells, both above `log(1.08)=.07696104` even
before the compact anchor adjustment. Keeping 40 rounding cells with those
fine-cell counts leaves too few cells for the long coarse branches at this
new core width. The exact-anchor adjustment is not the sole cause.

The current positive core center is `.1487934777423412`, whereas the fixed
transport anchor is `.15214223499880195`, about 34.65 measured core widths
away. The old common y equalizer sigma `.25` also reaches the ratio cap and
provides only `31.7200` predicted core cells on the current feature, below
the unchanged design target 32. This secondary geometric limitation was
not reached by the original x-quality-first paired rejection. Registered
sigma `.18` provides only `21.2934`; `.35` and `.5` both reach the padded
target `36.8` without saturating the ratio limit.

## New frozen proposal, still requiring a native transaction

The verified generic designer used exactly the policy's 70 x-axis choices,
four y sigmas and at most three full paired scores. Sixteen x axes and two
y axes were admitted. Its first ranked pair passed the full frozen score
after 4.156 s of design work:

- x: 32 fine cells, 16 rounding cells, fine-cell fraction `.5`;
- y: sigma `.5`;
- x/y maximum adjacent ratios: `1.066740804/1.066824363`;
- predicted frozen x/y core cells: `34.72313649/36.80000003`;
- predicted frozen left-front cells: `30.64282824`;
- relative peak change `6.77013e-5`, conservation defect `7.93981e-15`,
  and range violation `5.89807e-6`.

The whole box, 895-by-386 node count, original base axes, exact positive and
negative anchors, and cumulative peak budget are preserved. All quality and
field gates retain their previous thresholds. The smaller fine/rounding
allocation leaves more cells for the coarse branches; its maximum branch
log-spacing parameter is about `.064576`.

The output is
`result/longtime/20260908_campaign_v1/acceleration_lab/generic_design_step4768_tp1a85f477_ce45_4c47_b765_b060caeb17fd/candidate.mat`.
It contains `candidate`, `sourceCheckpoint`, `sourceReview`, `sourceSignature`,
`contract`, `invariantAudit`, and `nativeTransactionOptions`. The exact
original transaction interface is:

```matlab
q = load(candidateFile);
[checkpoint,audit] = ipm_gridlab_regrid_checkpoint( ...
    q.sourceCheckpoint, q.candidate.candidateX, ...
    q.candidate.candidateY, q.nativeTransactionOptions);
```

The stored options retain actual core floor 31/31, peak-jump bound `.002`,
mass bound `5e-12`, range bound `2e-4` and unchanged mesh limits. The output
filename is intentionally empty: this research wrapper performed no native
transaction or checkpoint commit. The caller must still check actual
post-transfer geometry, flow, references, clocks, history/prefix and the
cumulative budget before committing or continuing.

`old_failures.json`, `source_review.json`, `design.json`, and
`compact_summary.json` preserve every old/new result. The linked 216-file
source freeze is
`generic4768_source_626d392bc64642ccb7af8c1903268514`. The MATLAB batch exited
normally. Frozen readiness is not acceptance of a new native trajectory,
and does not establish a converged profile or singularity.
