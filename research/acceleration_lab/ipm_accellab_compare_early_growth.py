"""Array-only postprocess of strict native growth exports; no reconstructed PDE."""
import csv
import hashlib
import json
import math
import statistics
import sys
from pathlib import Path


def array(x):
    return x if isinstance(x, list) else [x]


def main(root):
    report = {"kind": "actual_H64_H8_early_growth_comparison_v1", "earlyWindow": [0, 1.6],
              "noLU": True, "noPDE": True, "allGrowthTransactionsReported": True, "cases": {}}
    for label in ("H64", "H8"):
        d = json.loads((root/label/"evidence.json").read_text())
        m, h = d["memory"], d["history"]
        c, g, hm, p = h["common"], h["gauge"], h["mesh"], m["policy"]
        members = m["referenceFamily"]["members"]
        rr = []
        for j, step in enumerate(c["acceptedStep"]):
            rr.append({"step": step, "tau": c["canonicalTau"][j], "physicalTime": c["physicalTime"][j],
                       "epoch": hm["remeshCount"][j], "core": [hm["coreGridPoints"][j], hm["verticalCoreGridPoints"][j]],
                       "quadraticStencilCenter": g["omegaGaugeQuadraticStencilCenterX"][j],
                       "wallNodalPeakX": c["wallPeakX"][j],
                       "wallCoreWidth": hm["trackedWallCoreWidth"][j], "verticalCoreWidth": hm["trackedVerticalCoreWidth"][j]})
        events = []
        for t in m["transactions"]:
            dec = t["controllerDecision"]
            cur, pred = dec["coreCells"], dec["predictedCoreCells"]
            current = [cur[a] < p["regridCoreTrigger"][a] for a in range(2)]
            forecast = [pred[a] < p["predictedCoreBuffer"][a] and cur[a] < p["targetCoreCells"][a] for a in range(2)]
            assert dec["reason"] == ("current_core_trigger" if any(current) else "forecast_core_trigger")
            fitpred = [cur[a]*math.exp(-p["maximumReviewInterval"]*max(0, dec["decayEvidence"][a]["fitDecay"])) for a in range(2)]
            evidence = []
            for a, de in enumerate(dec["decayEvidence"]):
                span = de["secantToTime"]-de["secantFromTime"]
                before = next((r for r in reversed(rr) if r["step"] <= de["secantFromStep"]), None)
                after = next((r for r in rr if r["step"] >= de["secantToStep"]), None)
                same = before is not None and after is not None and before["epoch"] == after["epoch"] == t["sourceRemeshCount"]
                exactfrom = next((r for r in rr if r["step"] == de["secantFromStep"] and r["epoch"] == t["sourceRemeshCount"]), None)
                exactpair = None
                if exactfrom is not None and de["secantToStep"] == t["sourceStep"]:
                    observed = -math.log(cur[a]/exactfrom["core"][a])/span
                    assert abs(observed-de["maximumSecant"]) <= 1e-10
                    exactpair = {"beforeSavedHistory": exactfrom, "afterSavedDecisionCore": cur,
                                 "recomputedSecant": observed, "afterPeakSelectorUnavailable": True}
                evidence.append({**de, "canonicalSpan": span,
                    "impliedRelativeCoreChange": math.expm1(-de["maximumSecant"]*span),
                    "maximumToFitRatio": de["maximumSecant"]/max(de["fitDecay"], 1e-300),
                    "historyBracket": [before, after], "sameEpochBracket": same,
                    "sampledNodalPeakChanged": before["wallNodalPeakX"] != after["wallNodalPeakX"] if same else None,
                    "sampledQuadraticCenterChanged": before["quadraticStencilCenter"] != after["quadraticStencilCenter"] if same else None,
                    "exactSavedCountPair": exactpair, "selectorCausationProven": False})
            source = members[t["sourceLevelId"]-1]
            eligible = [x for x in members if x["resourceAdmitted"] and x["qualityPassed"] and
                        all(b >= a for a, b in zip(source["cellFactors"], x["cellFactors"]))]
            eligible.sort(key=lambda x: (math.prod(x["nodeCount"]), x["index"]))
            earlier = [x["index"] for x in eligible[:next(k for k,x in enumerate(eligible) if x["index"] == t["targetLevelId"])]]
            events.append({"ordinal": t["sourceRemeshCount"]+1, "sourceStep": t["sourceStep"],
                "sourceTau": t["sourceCanonicalTime"], "sourcePhysicalTime": t["sourcePhysicalTime"],
                "sourceLevel": t["sourceLevelId"], "targetLevel": t["targetLevelId"],
                "sourceNodeCount": t["sourceNodeCount"], "targetNodeCount": t["targetNodeCount"],
                "growth": t["sourceLevelId"] != t["targetLevelId"], "reason": dec["reason"],
                "currentTriggerAxes": current, "forecastTriggerAxes": forecast,
                "sourceCore": cur, "predictedCore": pred, "decay": dec["decayEstimate"],
                "fitOnlyPredictionAtSameState": fitpred,
                "fitOnlyWouldTriggerAtSameState": any(fitpred[a] < p["predictedCoreBuffer"][a] and cur[a] < p["targetCoreCells"][a] for a in range(2)),
                "counterfactualTrajectoryComputed": False, "newAuditCore": t["coreCells"],
                "leftFrontCells": t["leftFrontCells"], "peakJump": t["relativePeakJump"],
                "candidateIndex": t["candidateIndex"], "attemptSummary": t["attemptSummary"],
                "eligibleOrder": [x["index"] for x in eligible],
                "earlierEligibleLevelsWithoutReturnedCandidate": earlier if t["candidateIndex"] == 1 else None,
                "absenceInferredFromFrozenEnumerationAndFirstCandidate": t["candidateIndex"] == 1,
                "exactAxisRejectionReasonsSaved": "axisReport" in dec,
                "xQuality": t["xQuality"], "yQuality": t["yQuality"], "decayEvidence": evidence})
        init = m["initialization"]
        first = {group: {key: array(value)[0] for key, value in h[group].items()} for group in ("common", "mesh", "gauge")}
        axes = {"x": init["selectedBaseX"], "y": init["selectedBaseY"]}
        bands = {"positiveX0to2": sum(0 <= x <= 2 for x in axes["x"]),
                 "positiveX2to8": sum(2 < x <= 8 for x in axes["x"]),
                 "positiveXBeyond8": sum(x > 8 for x in axes["x"]),
                 "Y0to1": sum(0 <= y <= 1 for y in axes["y"]),
                 "Y1to4": sum(1 < y <= 4 for y in axes["y"]),
                 "YBeyond4": sum(y > 4 for y in axes["y"])}
        early = [e for e in events if e["sourceTau"] <= 1.6]
        report["cases"][label] = {"checkpointFile": d["checkpointFile"], "strictNativeRead": True,
            "step": d["step"], "tau": d["canonicalTime"], "snapshotCount": d["snapshotCount"],
            "historyRows": len(rr), "historyEveryAcceptedStep": len(rr) == d["step"]+1,
            "initialFirstHistory": first, "initialAudit": init["audit"],
            "initialObserverUsed": init.get("observationFallback", {}).get("used", False),
            "fixedInitialNodeBands": bands, "policy": p,
            "earlyTransactionCount": len(early), "earlyGrowthCount": sum(e["growth"] for e in early),
            "earlyMedianTransactionSpacingTau": statistics.median([b["sourceTau"]-a["sourceTau"] for a,b in zip(early,early[1:])]),
            "earlyForecastTriggerAxes": [e["forecastTriggerAxes"] for e in early],
            "growthEvents": [e for e in events if e["growth"]], "allTransactions": events}
    before = json.loads((root/"before_sha.json").read_text())
    changed = [f for f,s in before.items() if hashlib.sha256(Path(f).read_bytes()).hexdigest() != s]
    assert not changed
    report["originalInputSHAUnchanged"] = True
    report["auditedOriginalFileCount"] = len(before)
    assert report["cases"]["H64"]["policy"] == report["cases"]["H8"]["policy"]
    report["registeredControllerPoliciesExactlyEqual"] = True
    with (root/"growth_comparison.json").open("x") as f:
        json.dump(report,f,ensure_ascii=False,indent=2,allow_nan=False)
    with (root/"all_transactions.csv").open("x") as f:
        w=csv.writer(f);w.writerow(["case","ordinal","step","tau","physical_t","source_level","target_level","reason","source_core_x","source_core_y","predicted_core_x","predicted_core_y","new_core_x","new_core_y","peak_jump"])
        for label,d in report["cases"].items():
            for e in d["allTransactions"]:
                w.writerow([label,e["ordinal"],e["sourceStep"],e["sourceTau"],e["sourcePhysicalTime"],e["sourceLevel"],e["targetLevel"],e["reason"],*e["sourceCore"],*e["predictedCore"],*e["newAuditCore"],e["peakJump"]])
    print(json.dumps({label: {k:d[k] for k in ("earlyTransactionCount","earlyGrowthCount","earlyMedianTransactionSpacingTau","fixedInitialNodeBands")} for label,d in report["cases"].items()},ensure_ascii=False))


if __name__ == "__main__":
    main(Path(sys.argv[1]).resolve())
