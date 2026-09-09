# Face-local mapped LF research protocol (kernel prototype only)

The current native `weno5FluxDerivative` splits mapped density Jq using one
alpha=max|a/J| per entire line. A large remote computational transport speed
therefore changes flux splitting near the core. This is a possible mechanism
for coordinate inconsistency, not yet its identified cause: the actual late
case also has `conservative_flux` wall transport and a retained normal
derivative contribution at the wall. The preceding failed three-rate solve
must remain rejected. The independent response decomposition now identifies
a real global-alpha max switch in the native RHS. Its large effect on the
original crossing-width functional is reduced roughly 640 times by the
smooth integral-width functional. Neither width response alone quantifies
full PDE coordinate error.

## Discrete construction to test

At face i+1/2, the existing positive reconstruction uses five extended nodes,
and the negative reconstruction uses the reversed five nodes displaced by
one. Their union has six nodes. Define alpha_face=max|a/J| on that exact union,
including required ghost nodes. Use this *same scalar for the face* in both
stencils:

```
F_plus(j;face)  = (a_j q_j + alpha_face J_j q_j)/2
F_minus(j;face) = (a_j q_j - alpha_face J_j q_j)/2.
```

Perform the existing left-biased WENO-Z reconstruction on the positive
stencil and reversed negative stencil; sum once to obtain a single face
flux. The derivative remains `diff(faceFlux)/(dXi*J_node)`. Reusing that
single face value at adjacent nodes preserves the mapped telescoping flux
identity irrespective of nonlinear weights. Independent left/right alpha
values at the same face are forbidden. Nonpositive extrapolated ghost J is
an explicit rejected mapping, not clipped or replaced by an absolute value.

This is a different spatial discretization. Smooth consistency follows by
applying the two reconstructions with a facewise constant alpha; their
leading reconstructed alpha*Jq terms cancel. It does not automatically give
a global fifth-order error estimate at nonsmooth max/absolute-value switches,
critical points, or extrapolation boundaries. Empirical order and stability
must be measured and any lower order retained.

## A separate normalization choice must not be hidden

The native normalized WENO-Z weights use the largest split flux on the whole
physical line, with the existing sqrt(eps)*commonScale floor and regularizer
max(epsilonFloor,1/(N-1)^2). There is no single array of split flux when alpha
depends on the face. Two distinguishable studies would be:

1. **Change alpha only.** For each face's alpha, compute both split scales
   over all physical nodes on the line, then use those scales for that
   face's two reconstructions. This preserves the original normalization
   definition conditional on alpha. It is expensive (quadratic work per
   line) and is only a small research isolation. It retains long-range
   dependence through normalization.
2. **Localize alpha and normalization.** Use the maximum split magnitudes on
   the same six-node face union, preserving the sqrt(eps) shared scale floor
   and the registered normalized epsilon. This is an additional numerical
   change and must have its own label and comparison.

The original native global-alpha/global-scale operator remains the baseline.
On constant velocity and constant J, the first variant should reproduce the
native face flux and derivative up to the same arithmetic paths. The second
variant need not reproduce it because its weight scales differ. Neither
variant can be called an acceleration of the original semidiscrete equation.

## Gates before any PDE integration

- Record every face alpha and its six-node speed bound; both reconstructions
  must use it, including exterior faces. Verify the weighted telescoping
  boundary identity against exported face fluxes.
- Test constant-state free-stream cancellation for positive, negative, zero,
  and small amplitudes on nonuniform smooth mapped grids. The complete
  advective form retains D(aq)-qD(a), using the identical new operator and
  boundary data for q and the unit field. This is not a license to remove
  the original conservative-wall y derivative.
- Preserve original extrapolation/reflection and wall/bulk assembly; test
  both wall modes explicitly. A reflected impermeable case and an open
  extrapolated case need separate boundary evidence. No new ghost limiter
  or normal-velocity modification is included.
- Use smooth mapped manufactured fluxes and velocity sign changes, both
  full-domain and interior errors, at 49/65/97/129 one-dimensional sizes.
  Polynomial reproduction alone is insufficient for the nonlinear WENO
  derivative. Measure convergence; do not predeclare fifth order.
- Extend a remote region while keeping the same local samples and mapping
  to expose alpha versus normalization dependence. The normalized epsilon
  itself depends on N, so repeat with the same explicit effective epsilon
  when isolating locality, and retain native-N results separately.
- At fixed tiny IPM rho/physicalflow, compare rate-response coordinate-chain
  defects for native and research operators, along the identical amplitude,
  translation and dilation directions. No different mode gets a favorable
  parameter search or changed tolerance.

Only after these pass should a separately registered tiny 49/65 original
physical-time comparison be attempted, with all density and scale variables
evolved, original physical box/Green boundary, and explicit method labels.
Full-field and inner C1 errors, outside-core velocity, fixed physical tails,
constant-state and finite-state gates, timestep refinement and spatial
refinement all remain necessary. No candidate state may be inserted in the
main checkpoint chain. A reduced one-state coordinate defect alone is not
physical-time validation, faster shape convergence, or a singularity result.

## Executed no-LU kernel prototype

`ipm_accellab_local_lf_derivative` and `ipm_accellab_test_local_lf` implement
and test the alpha-only conditional-global-scale variant. The optional local
normalization branch is present with a distinct explicit option but has not
been run or qualified. No production function was modified.

The copied line-alpha branch reproduced the maintained kernel bitwise on all
four MMS grids. The face-local conditional-global branch also reproduced the
native constant-speed, constant-metric result bitwise. Constant-state
free-stream errors were at most 3.27e-14; weighted flux-boundary identities
were within 5.56e-17, including a reflected endpoint fixture. A positive
nodal metric with nonpositive extrapolated ghost values was explicitly
rejected without clipping.

For the prescribed smooth mapped MMS with a velocity zero crossing, the
49/65/97/129 local-alpha full errors were
`9.258e-5,2.147e-5,2.574e-6,5.654e-7`, with interior errors
`3.678e-5,8.958e-6,1.190e-6,2.832e-7`. The measured orders were about five
on this example. Native full errors were
`1.020e-4,2.359e-5,2.835e-6,6.481e-7`; no general-order assertion follows
from one profile.

On one fixed grid, changing only velocity outside |x|>.55 changed the
derivative in |x|<.3 by 4.160e-6 for native global-alpha and 8.527e-14 for
the local-alpha variant. The latter still has global weight normalization;
this single result does not prove strict locality for arbitrary fields.

Evidence is saved under
`local_lf_kernel_tpc2e4230c_fec1_4890_95be_ec5874a6fa01` in the campaign
acceleration output. These are small array tests, with zero Poisson, zero
PDE steps and no physical-time validation. Two-dimensional wall assembly,
true IPM spatial/temporal convergence and geometric covariance at changed
rates are outstanding; this kernel is not a production-ready replacement.
