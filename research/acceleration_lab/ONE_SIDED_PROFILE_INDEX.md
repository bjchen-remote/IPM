# Fixed-radius one-sided Profile index

Six native finite-time profiles at parent-equivalent tau 8.4432 through 11.2432
were examined on the preregistered physical radii
`r = logspace(-5,-2,81)`. All six checkpoints passed strict native validation
at ten threads. Their continuous geometry and complete saved wall curves
reproduced bitwise from their own fields, axes, and scales. No LU, fresh RHS,
or PDE step was computed. The diagnostic process, session 65885, exited 0.

An earlier `ipm_accellab_scale_index.py` examines a conditional instantaneous
temporal amplitude-width index; `measure_continuous_amplitude_width_index.m`
also includes measured inner-shape drift. Neither computes this fixed-radius
one-sided spatial observable. Their results are not substituted here.

## Definition and fixed admission

Let `H[rho]` be the same frozen, linear tensor cubic Hermite polynomial used
for both field values and its spatial derivative. Its nodal jets come from
the maintained derivative matrices. With computational peak position `a`,
physical scales `Cx,Comega`, and `q=a+Cx*r`, the observation is

```text
Delta(r) = (H[rho](q,0) - H[rho](a,0))/Comega
Delta_r(r) = Cx/Comega * partial_X H[rho](q,0)
a_eff(r) = r * Delta_r(r) / Delta(r).
```

The center and physical 90% wall width are the original saved continuous
geometry. In particular, the numerator is **not** a separately interpolated
`rhoX` field: it differentiates the same polynomial defining the denominator.
The recorded continuous peak of interpolated `rhoX` and the density
polynomial's derivative at the center can differ slightly; their ratio is
saved rather than silently identifying them.

An index is reported only where all of the following predeclared conditions
hold: the query is inside the native domain; the density increment is positive
and finite and its derivative is finite; `r >= 4*physicalWallWidth`; there are
at least eight native knots between center and query; and the widest native
cell intersecting that interval is at most `r/4`. Every radius, ingredient,
flag, exclusion reason, node count, and cell-width ratio is saved. Excluded
indices remain NaN. No best window or fitted constant exponent was selected.
The negative side is not evaluated or combined with this observable.

An analytic cubic on a nonuniform 49-by-33 grid, tested at nine combinations
of physical coordinate/amplitude units and all 81 radii, checks these formulas.
Maximum effective-index error was `9.835e-11`; relative derivative error was
`7.786e-15`, below the predeclared `5e-9` test threshold. This verifies the
arithmetic and unit conversion, not the spatial error of the actual profiles.

## Complete finite-time observations

All 81 radii lie in-domain and have positive finite density increments in all
six cases. The four-width exclusion determines the final admitted counts;
the independent node/cell conditions also reject some smaller radii.

The two comparison radii in the last columns are the fixed grid indices 61
and 81, not selected by their resulting index values.

| Parent-equivalent tau | Admitted radii | Smallest admitted r | a_eff at r=0.00177827941 | a_eff at r=0.01 |
|---:|---:|---:|---:|---:|
| 8.443201 | 24/81 | 0.001372461 | 0.626360 | 0.515830 |
| 9.243201 | 30/81 | 0.000817523 | 0.574288 | 0.495472 |
| 10.043201 | 36/81 | 0.000486968 | 0.536225 | 0.480614 |
| 10.643201 | 40/81 | 0.000344747 | 0.514616 | 0.472061 |
| 10.843201 | 41/81 | 0.000316228 | 0.508357 | 0.469549 |
| 11.243201 | 44/81 | 0.000244062 | 0.497300 | 0.465083 |

At the latest time the admitted curve runs from 0.643465 at its smallest
admitted radius to 0.465083 at 0.01. There remains substantial radius
dependence. These curves describe the crossover outside the current resolved
core; they do not determine a unique constant spatial exponent.

Each finite-time density remains a smooth numerical profile. Where its
center derivative is nonzero, the small-radius Taylor limit of this ratio
is 1. The deliberately excluded inner radii therefore matter when discussing
limits: a finite-time outer-core observation does not establish the
infinite-time limit or justify exchanging time and radius limits. Centers
also move with the measured peak; these are not increments about one fixed
spatial point throughout the sequence.

## Separate temporal scale ratios

The adjacent-sample slope is
`Delta log(G*w) / Delta log(w)`, using the continuous physical peak as primary
`G`, and the actual physical wall 90% width as `w`. A secondary column uses
the native whole-field physical `rhoX` maximum to make the meaning of `G`
explicit. This temporal ratio is dimensionless, but is not the spatial
`a_eff(r)` above.

| Adjacent tau interval | Continuous-peak G | Native-grid G | Remesh count increase |
|---|---:|---:|---:|
| 8.443201–9.243201 | 0.379816 | 0.379978 | 1 |
| 9.243201–10.043201 | 0.375625 | 0.375263 | 1 |
| 10.043201–10.643201 | 0.374953 | 0.375007 | 1 |
| 10.643201–10.843201 | 0.371166 | 0.370989 | 1 |
| 10.843201–11.243201 | 0.372717 | 0.372726 | 0 |

Regridding changes and spatial discretization error have not been removed or
bounded in these ratios. An interval with unchanged remesh count still lacks
an independent spatial-error bound. The small difference between the two G
definitions is not such a bound. Neither table proves a limiting Hölder
exponent or establishes strong singularity.

## Reproducibility

The executable is `ipm_accellab_one_sided_profile_index.m`; its inputs are the
existing `rendered/profiles.mat` and a previously unused output directory.
It uses the six checkpoints and five frozen helpers registered by
`continuous_profile_tau1124_figure_v1/registration.json`. All upstream file
hashes were verified before execution.

Artifacts under
`result/longtime/20260908_campaign_v1/acceleration_lab/`:

- `one_sided_index_source_b52d3184d03d4189b12f69da5e75a60e/`: executable,
  frozen upstream helpers, and preregistration/hashes.
- `one_sided_index_b52d3184d03d4189b12f69da5e75a60e/`: `indices.json`,
  `indices.mat`, six complete 81-row CSV files, inclusion masks and the
  original figure.
- Its `presentation_v2/one_sided_indices.png` is a theme-only redraw from
  the saved `indices.mat`, with readable title and legend colors. It performs
  no checkpoint read or numerical recomputation.

The primary and secondary temporal slopes and all 486 radius rows remain
available, including every excluded observation. No production file,
trajectory, acceptance threshold, or original result was changed by this
diagnostic.
