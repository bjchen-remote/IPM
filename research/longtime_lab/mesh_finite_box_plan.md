# Independent finite-box experiments

These branches investigate the H=1e6 spatial-span bottleneck without changing
the main domain or weakening its mesh gates. The global y quadrature floor is
`min(q_y) >= 1e-8 * H / ny`. At ny=513 it is about 1.95e-5 for H=1e6 and
1.95e-7 for H=1e4. A smaller admissible box therefore gives more room for
inner contraction as well as freeing tail nodes. This algebra is not an
accuracy argument for choosing the smaller box.

## Existing initial-RHS evidence

Source: `result/verification/q512_comega_initial_box_screens_v4/`;
screen ID `screen_20260905T110530211Z_tpbe425904_b6b4_4b3f_9aa4_c02b888cd290`.
Read-only stored-axis audit: `result/longtime/20260908_initial_box_axis_audit_v1.mat`
and `.json`; its MATLAB execution and Code Analyzer passed.

| H | c_l | c_omega | P | Relative c_l difference from H=1e6 |
| --- | --- | --- | --- | --- |
| 1e3 | 0.488133800036337 | 0.125042863987632 | 0.686073738019381 | 6.73049e-4 |
| 1e4 | 0.488429879504927 | 0.125042863798633 | 0.686073738019381 | 6.69031e-5 |
| 1e5 | 0.488459589077373 | 0.125042864250448 | 0.686073738019381 | 6.08049e-6 |
| 1e6 | 0.488462559169782 | 0.125042865249242 | 0.686073738019381 | 0 |

All runs used the schema-4 exact quadratic-peak/transport-anchor contract and
Green boundaries. The x platform [0.9475,1.2975] has 161 shared nodes within
absolute tolerance 1e-12, giving the same initial P and peak stencil. However,
X in [-2,2] contains 681/651/627/611 nodes at the four H values, and Y in [0,1]
contains 195/167/147/132. Only y=0 is shared within 1e-12. Minimum dy varies
from 0.000625049420 to 0.000625137032. Equal minimum spacing and the protected
wall stencil therefore do not establish equal full inner grids.

H=1e4 is the first conservative candidate, with H=1e5 as an intermediate
control and H=1e3 a further exploratory candidate. Initial c-omega differences
from H=1e6 are only about 1e-8 relative. These values justify tests, not a
smaller production domain: the old experiment took zero PDE steps, and outer
discretization changed with H. Its far-velocity ratio is a velocity amplitude
diagnostic and must never be relabeled a boundary error.

## Stage A: preserve all current inner nodes

Proposed later source is the coordinator's accepted adaptive stage-1 checkpoint:

```text
result/longtime/20260908_campaign_v1/adaptive_campaign/stage_001/
checkpoint_20260905T204759339_dynamic_isotropic_1025x513_tp8c8b7442_a9d7_4682_b700_7cf8af6fe3f7_step0000002031.mat
```

Low-cost unvalidated monitoring reads tau=4.283201944539746 and physical
t=1.7798569844190184. The actual screen must validate this native checkpoint
normally at 10 threads. It retains all saved no-feedback runtime references,
scales and common rho samples. It must run serially after allocation of a
q512 factorization window; it has not been launched.

```matlab
addpath('research/longtime_lab');
mesh_screen_late_box(sourceCheckpoint,newOutputDirectory, ...
    [1,.3,.1,.03,.01,.003,.001]);
```

The wrapper preregisters its criteria before building any cropped operator.
Actual endpoints are existing nodes and need not be equal in x/y:

| Fraction | nx / ny | Actual x half-width / y maximum | Removed x / y nodes |
| --- | --- | --- | --- |
| 1 | 1025 / 513 | 1000000 / 1000000 | 0 / 0 |
| .1 | 963 / 449 | 97346.2883 / 98151.7053 | 62 / 64 |
| .01 | 901 / 386 | 9792.52309 / 9986.69546 | 124 / 127 |
| .001 | 837 / 322 | 948.826061 / 977.943844 | 188 / 191 |

At fraction .01, the rectangle retains 66.14% of original grid samples.
This is a node-count calculation, not a measured memory or speed saving.
The complete ladder is saved in `20260908_late_box_crop_plan_stage1.json`.

Static candidate criteria are core velocity relative L2 and relative infinity
differences <=1e-3 on X=[.5,1.5], Y=[0,1]; each rate absolute difference must
be <=max(1e-6, 1e-3 times its full-box magnitude). Require P relative change
and peak-position absolute change <=1e-10, Poisson residual <=1e-8, preserved
references and all existing mesh gates (ratio1.08, curvature.01, rcond1e-9,
global q ratio1e-8, local q/control [.35,1.65]). The full-box overlap control
must reach 1e-11 for both velocity norms and approximately 1e-10 relative for
rates. A crop is nominated only when it and its next larger ladder member
pass. Every failure is retained; static nomination never approves dynamics.

Even this fixed-inner-node crop recomputes outer derivative closures and Green
traces from truncated sources. It does not retain the omitted source's Green
contribution, reconstruct prior smaller-box dynamics, or bound H=1e6 error.

## Stage B: matched physical time, then fresh roots

1. For an admitted crop, build an explicit independent branch preserving its
   inner nodes and the complete normalization contract. The native same-grid
   resume whitelist does not authorize changing domain bounds. A cropped
   checkpoint must therefore use an explicitly reviewed research transaction
   with new provenance; do not alter the original checkpoint or bypass its
   signature. The current static script deliberately does not implement this.
2. Compare that branch with the same full-box source at control Delta tau
   .01, .05, then .2, ending the crop run at the full-box run's exact physical
   time. Use fresh paired physical axes. Report laboratory-frame rho, rho-x,
   rho-y L2/Inf, physical gradient maximum, quadratic P, both rates, strict
   physical widths, Poisson residual, all trust gates, and physical coverage.
   Separately show any aligned/inner-normalized comparison so it cannot hide
   displacement, amplitude or contraction errors. Halve the CFL in both
   branches to distinguish boundary differences from temporal error.
3. Suggested preregistered short-branch candidate budgets are rho L2/Inf
   <=1e-3, both gradient L2/Inf and gradient maximum differences <=2e-3,
   rate differences <=max(1e-6,2e-3 times reference magnitude), and physical
   strict-width differences <=5e-3. These are research selection budgets,
   not a global convergence certificate; persistent drift across the three
   durations requires the larger box or further resolution.
4. Start independent fresh roots for H near 1e4, 1e5, 1e6 (include 1e3 only
   if preceding evidence supports it). Construct these by truncating one
   common root axis outside a protected inner rectangle with full stencil
   margins, rather than regenerating each H with the old box factory. Permit
   different total node counts to isolate box effects. Follow at matched
   physical times and repeat with independent inner refinement.
5. Only after fixed-inner-node box effects are controlled, redistribute the
   released point budget on the selected box. Compare original and redistributed
   grids separately through the existing audited transfer/short-time path.
   Domain truncation and extra inner resolution need separate evidence.

Current implementation status: the two basic diagnostics passed tiny v3.
The coordinator then successfully ran the static-classification wrapper and
its extra report fields in `20260908_box_wrapper_tiny_v1`: the full control
passed, the sensitive crop failed, threads were restored, Code Analyzer
passed, and no dynamic result was promoted. Earlier Qt/neon startup failures
remain recorded. This tiny evidence does not establish a q512 frozen-crop,
cropped-dynamic-branch or fresh-root box convergence result.

Subsequent execution update: the coordinator did run the q512 static screen
on stage-1 step2031, saved in `20260908_campaign_v1/late_box_tau428_v1`.
Fractions through .003 passed, while .001 failed both velocity criteria.
Fraction .01 gives velocity L2/Inf changes 1.37981e-4/1.20419e-4. This
supersedes the unexecuted status in the Stage-A planning text above, without
turning its static result into dynamic evidence. The next five-case dynamic
protocol and its actual tiny covariance tests are described in
[mesh_fresh_box_protocol_plan.md](mesh_fresh_box_protocol_plan.md).
