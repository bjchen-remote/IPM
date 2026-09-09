#!/usr/bin/env python3
"""Low-cost planning forecast only; no native signature/trust validation.

Read same-remesh scalar histories through the monitoring HDF5 helper. The
inflated tail slopes are a heuristic scheduling margin, not error bounds or
permission to bypass endpoint/native checkpoint/core gates in a real solve.
"""

import argparse
import json
import math
from pathlib import Path

from read_checkpoint_scalars import HDF5


def slope(x, y):
    mx, my = sum(x) / len(x), sum(y) / len(y)
    denominator = sum((v - mx) ** 2 for v in x)
    if denominator <= 0:
        raise ValueError("No positive time span in same-grid tail")
    return sum((u - mx) * (v - my) for u, v in zip(x, y)) / denominator


def predict(checkpoint, intervals, library):
    h = HDF5(library, checkpoint)
    try:
        def hist(group, name):
            return h.numeric("/checkpoint/payload/log/history/" + group + "/" + name)[0]

        tau = hist("common", "canonicalTau")
        remesh = hist("mesh", "remeshCount")
        x = hist("mesh", "coreGridPoints")
        y = hist("mesh", "verticalCoreGridPoints")
        safety = hist("mesh", "safetyFactor")
        cl, cw = hist("common", "c_l"), hist("common", "c_omega")
        epoch = hist("common", "physicalTime")[-1]
        log_cl = h.numeric("/checkpoint/payload/state/scale/logC_l")[0]
        log_cw = h.numeric("/checkpoint/payload/state/scale/logC_omega")[0]
        time_factor = math.exp(log_cl - log_cw)
    finally:
        h.close()
    if len({len(v) for v in [tau, remesh, x, y, safety, cl, cw]}) != 1:
        raise ValueError("Unpaired scalar histories")
    fits = []
    for width in [.08, .12, .2]:
        indices = [i for i in range(len(tau))
                   if remesh[i] == remesh[-1] and tau[-1] - width <= tau[i] <= tau[-1]]
        if len(indices) < 8:
            continue
        values = [v[i] for v in [tau, x, y, safety, cl, cw] for i in indices]
        if not all(math.isfinite(v) for v in values) or min(x[i] for i in indices) <= 0 or min(y[i] for i in indices) <= 0:
            raise ValueError("Nonfinite/nonpositive forecast samples")
        z = [tau[i] for i in indices]
        fits.append({"tauWindow": width, "samples": len(indices),
                     "xLogCoreDecay": max(0, -slope(z, [math.log(x[i]) for i in indices])),
                     "yLogCoreDecay": max(0, -slope(z, [math.log(y[i]) for i in indices])),
                     "safetySlope": max(0, slope(z, [safety[i] for i in indices])),
                     "maximumCLMinusCOmega": max(0, max(cl[i] - cw[i] for i in indices))})
    if not fits:
        raise ValueError("Insufficient same-grid tail for a forecast")
    inflation = 1.25
    decay = [inflation * max(f[name] for f in fits) for name in ["xLogCoreDecay", "yLogCoreDecay"]]
    safety_slope = inflation * max(f["safetySlope"] for f in fits)
    kappa = inflation * max(f["maximumCLMinusCOmega"] for f in fits)
    forecasts = []
    for elapsed in intervals:
        if elapsed <= 0 or not math.isfinite(elapsed):
            raise ValueError("Intervals must be finite and positive")
        if kappa * time_factor * elapsed >= 1:
            forecasts.append({"elapsedPhysicalTime": elapsed, "predictionAvailable": False})
            continue
        dtau = -math.log1p(-kappa * time_factor * elapsed) / kappa if kappa else time_factor * elapsed
        cores = [v * math.exp(-g * dtau) for v, g in zip([x[-1], y[-1]], decay)]
        predicted_safety = safety[-1] + safety_slope * dtau
        forecasts.append({"elapsedPhysicalTime": elapsed, "predictionAvailable": True,
                          "absolutePhysicalTarget": epoch + elapsed, "canonicalInterval": dtau,
                          "coreCells": cores, "safety": predicted_safety,
                          "planningMarginPassed": min(cores) >= 16 and predicted_safety < .65})
    admitted = [f["elapsedPhysicalTime"] for f in forecasts if f.get("planningMarginPassed", False)]
    return {"kind": "same_grid_scalar_planning_forecast", "checkpointFile": str(checkpoint),
            "nativeSignatureValidated": False, "trustedHistoryValidated": False,
            "pdeAdvanced": False, "sourceCanonicalTime": tau[-1], "sourcePhysicalEpoch": epoch,
            "sourceRemeshCount": remesh[-1], "sourceCoreCells": [x[-1], y[-1]],
            "sourceSafety": safety[-1], "canonicalTimeFactor": time_factor,
            "fits": fits, "heuristicInflation": inflation, "inflatedDecay": decay,
            "inflatedCLMinusCOmega": kappa, "forecasts": forecasts,
            "selectedElapsedPhysicalTime": max(admitted) if admitted else None,
            "actualRunGates": {"minimumCoreCells": 14, "maximumSafetyExclusive": .70,
                               "nativeCheckpointAndResultPairing": True, "entireHistoryTrusted": True,
                               "physicalEndpointRequired": True},
            "interpretation": "Monitoring-based planning only. Extrapolation can fail; native validation and actual runtime gates remain mandatory."}


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("checkpoint", type=Path)
    parser.add_argument("output", type=Path)
    parser.add_argument("--intervals", type=float, nargs="+", default=[.005, .01, .02])
    parser.add_argument("--library", default="/Applications/MATLAB_R2026a.app/bin/maca64/libhdf5.310.dylib")
    args = parser.parse_args()
    if args.output.exists():
        raise ValueError("Refusing to replace an existing forecast")
    report = predict(args.checkpoint, args.intervals, args.library)
    args.output.write_text(json.dumps(report, indent=2) + "\n")
    print(json.dumps({"selectedElapsedPhysicalTime": report["selectedElapsedPhysicalTime"],
                      "forecasts": report["forecasts"]}, indent=2))


if __name__ == "__main__":
    main()
