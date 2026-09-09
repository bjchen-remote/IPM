"""Read saved JSON diagnostics; infer no singular time or limiting regularity.

Run with Python's standard library. This never opens a native checkpoint,
builds operators, changes a PDE state, or evaluates a new RHS.
"""
import argparse
import hashlib
import json
import math
from pathlib import Path
import uuid


def fit(records):
    x = [math.log(r["physicalWallWidth"]) for r in records]
    y = [math.log(r["physicalRhoXInf"]) for r in records]
    xm, ym = sum(x) / len(x), sum(y) / len(y)
    xx = sum((v - xm) ** 2 for v in x)
    assert xx > 0
    slope = sum((a - xm) * (b - ym) for a, b in zip(x, y)) / xx
    errors = [b - (ym + slope * (a - xm)) for a, b in zip(x, y)]
    return {
        "tau": [r["tau"] for r in records],
        "sampleCount": len(records),
        "logGradientVersusLogPhysicalWidthSlope": slope,
        "conditionalIncrementScaleIndex": 1 + slope,
        "maximumLogGradientFitResidual": max(map(abs, errors)),
        "fittingMethod": "unweighted_OLS_no_independence_or_errorbar_claim",
    }


def analyze(exact_file, series_file, output_root):
    exact = json.loads(exact_file.read_text())
    series = json.loads(series_file.read_text())
    e = exact["exactInnerCoordinates"]
    assert e["valid"] and e["level"] == 0.9
    records = series["records"]
    assert all(r["fullHistoryTrusted"] for r in records)
    assert all(a["tau"] < b["tau"] for a, b in zip(records, records[1:]))
    assert all(r["physicalWallWidth"] > 0 and r["physicalRhoXInf"] > 0 for r in records)
    beta, cl, cw = e["logScaleXRate"], exact["cl"], exact["cOmega"]
    p_rate = e["peakPrime"] / e["peak"]["value"]
    width_rate = beta - cl
    gradient_rate = cl - cw + p_rate
    increment_rate = beta - cw + p_rate
    assert width_rate < 0 and all(map(math.isfinite, [width_rate, gradient_rate, increment_rate]))
    alpha = increment_rate / width_rate
    focus = min(range(len(records)), key=lambda i: abs(records[i]["tau"] - exact["canonicalTime"]))
    assert abs(records[focus]["tau"] - exact["canonicalTime"]) < 1e-10
    windows = {"all_saved_samples": fit(records)}
    if focus >= 2:
        windows["three_points_ending_at_exact_sample"] = fit(records[focus - 2:focus + 1])
    if 0 < focus < len(records) - 1:
        windows["three_points_centered_on_exact_sample"] = fit(records[focus - 1:focus + 2])
    for window in windows.values():
        window["minusConditionalInstantaneousIndex"] = window["conditionalIncrementScaleIndex"] - alpha
    pairs = [fit([a, b]) for a, b in zip(records, records[1:])]
    sources = [{"path": str(p.resolve()), "sha256": hashlib.sha256(p.read_bytes()).hexdigest()}
               for p in [exact_file, series_file]]
    result = {
        "kind": "conditional_instantaneous_inner_increment_scale_index",
        "status": "read_only_diagnostic_completed",
        "sources": sources,
        "exactSample": {
            "tau": exact["canonicalTime"], "physicalTime": exact["physicalTime"],
            "beta": beta, "cl": cl, "cOmega": cw,
            "peak": e["peak"]["value"], "actualPeakPrime": e["peakPrime"],
            "actualPeakRelativeRate": p_rate,
            "logPhysicalQuadraticPeakGradientRate": gradient_rate,
            "logPhysicalInnerWidthRate": width_rate,
            "logDensityIncrementProxyRate": increment_rate,
            "logGradientVersusLogPhysicalWidthSlope": gradient_rate / width_rate,
            "conditionalIncrementScaleIndex": alpha,
            "pairedNodeWidthRelativeDifferenceFromExactQuadraticWidth":
                records[focus]["strictWallWidth"] / e["wallCoreWidth"] - 1,
            "nonzeroExactMotionHoldoutL2Ratio": e["metrics"]["holdout"]["l2Ratio"],
            "sourcePoissonBackwardResidual": exact["poissonResidual"],
        },
        "adjacentWindowFits": pairs,
        "registeredWindowFits": windows,
        "densityIncrementProxies": [
            {"tau": r["tau"], "physicalGradientTimesPhysicalWidth":
             r["physicalRhoXInf"] * r["physicalWallWidth"]} for r in records],
        "formula": "alpha_x=(beta-cOmega+PPrime/P)/(beta-cl)=1+d(logGstar)/d(logPhysicalWx)",
        "shapeCondition": "For a fixed normalized interval, DeltaRho=Gstar*physicalWx*A(tau); the formula assumes APrime/A is negligible. Its exact correction is (APrime/A)/(beta-cl).",
        "observableDifference": "Exact formula uses the tracked quadratic peak and its exact level width; paired series uses global grid max |rho_x| and maintained node-threshold width. The series spans remeshes.",
        "limitations": [
            "A temporal amplitude-width slope is not a measured spatial Holder exponent at any finite time.",
            "No limiting profile, convergence of the physical center, fixed limiting exponent, or singularity is established.",
            "No singular time T or powers of T-t are fitted or assumed.",
            "A hypothetical G~(T-t)^(-p), width~(T-t)^q would only give alpha=1-p/q; p, q and T remain undetermined.",
            "Shape drift, different peak/width observables, template switches, remesh transfer, grid and box errors remain relevant.",
            "The exact sample's original Poisson ill-conditioning warning remains in its source report and companion runtime_warnings.json.",
        ],
        "newPdeSteps": 0, "newRhsEvaluations": 0,
        "limitingHolderExponentEstablished": False, "strongSingularityEstablished": False,
    }
    destination = output_root / ("scale_index_" + uuid.uuid4().hex)
    destination.mkdir(parents=True, exist_ok=False)
    (destination / "report.json").write_text(json.dumps(result, indent=2, ensure_ascii=False) + "\n")
    centered = windows.get("three_points_centered_on_exact_sample", windows["all_saved_samples"])
    note = f"""# Conditional inner density-increment scale index

The exact sample at tau={exact['canonicalTime']:.10f} gives alpha_x={alpha:.8f}.
The centered three-sample fit gives {centered['conditionalIncrementScaleIndex']:.8f}.
These are finite-window diagnostics, not a limiting Holder exponent.

Let C_l=exp(logC_l), C_w=exp(logC_omega), w_x be the computational
exact 0.9 peak width, and ell_x=w_x/C_l its physical width. The physical
quadratic-peak gradient is G_star=P*C_l/C_w. Along a fixed normalized
inner interval, integrating a shape phi yields

    Delta rho = G_star * ell_x * A(tau),    A = integral phi dxi.

If A is constant or changes negligibly, the logarithmic rates per canonical
time are

    (log G_star)' = c_l-c_w+P'/P = {gradient_rate:.10f}
    (log ell_x)'  = beta-c_l    = {width_rate:.10f}
    (log Delta rho)' = beta-c_w+P'/P = {increment_rate:.10f}.

Their ratio gives alpha_x={alpha:.8f}; equivalently the gradient-width
slope is alpha_x-1={alpha-1:.8f}. The actual P'/P={p_rate:.4e} is used.
For drifting shape, add (A'/A)/(beta-c_l), which this proxy does not measure.

The existing series fits global grid max |rho_x| against its recorded
physical node-threshold width. It spans mesh transactions and does not use
identical observables to the exact instantaneous formula. At the exact
sample the recorded node width differs from the quadratic-level width by
{result['exactSample']['pairedNodeWidthRelativeDifferenceFromExactQuadraticWidth']:.5e} relatively.
The roughly 0.007 difference from the centered fit is retained as a
diagnostic discrepancy, not converted to a confidence interval.

The index concerns how amplitude and width co-vary along this evolving
solution. Concluding a limiting spatial Holder exponent would additionally
require a controlled limiting shape and spatial center, scale and window
stability, and discretization/box checks. No T is guessed. Even under a
hypothetical time-power ansatz, only alpha=1-p/q would follow; neither p nor
q is separately determined. No strong singularity or completed long-time
goal is claimed here.
"""
    (destination / "derivation.md").write_text(note)
    print(json.dumps({"outputDirectory": str(destination.resolve()),
                      "instantaneousIndex": alpha,
                      "centeredFitIndex": centered["conditionalIncrementScaleIndex"]}))
    return result


if __name__ == "__main__":
    project = Path(__file__).resolve().parents[2]
    campaign = project / "result/longtime/20260908_campaign_v1/acceleration_lab"
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--exact", type=Path, default=campaign / "exact_inner_q512_tpc9fd0098_f80e_450c_aacf_50d097ce075b/report.json")
    parser.add_argument("--series", type=Path, default=project / "result/longtime/20260908_profile_progress_v2/analysis/profile_metrics.json")
    parser.add_argument("--output-root", type=Path, default=campaign)
    args = parser.parse_args()
    analyze(args.exact, args.series, args.output_root)
