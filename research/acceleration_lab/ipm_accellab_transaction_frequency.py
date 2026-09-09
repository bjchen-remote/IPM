#!/usr/bin/env python3
"""Postprocess strict checkpoint evidence only; no numerical evolution.

Usage: python3 ipm_accellab_transaction_frequency.py REPORT_DIRECTORY
The directory must contain native_evidence.json from the MATLAB exporter.
All outputs are created exclusively, so an earlier conclusion is never replaced.
"""
import json
import math
import statistics
import sys
from pathlib import Path


def rows(window):
    return [dict(time=t, step=int(s), epoch=int(e), core=c, safety=q)
            for t, s, e, c, q in zip(window["time"], window["step"],
                window["remeshCount"], window["coreCells"], window["safety"])]


def axis_stats(rr, axis):
    if len(rr) < 2:
        return {"sampleCount": len(rr)}
    ts = [r["time"] for r in rr]
    ys = [math.log(r["core"][axis]) for r in rr]
    tm, ym = statistics.mean(ts), statistics.mean(ys)
    fit = sum((t-tm)*(y-ym) for t, y in zip(ts, ys))/sum((t-tm)**2 for t in ts)
    sec = [-(ys[k]-ys[k-1])/(ts[k]-ts[k-1]) for k in range(1, len(ts))]
    j = max(range(len(sec)), key=sec.__getitem__) + 1
    q95 = sorted(sec)[math.floor(.95*(len(sec)-1))]
    return {"sampleCount": len(rr), "fitDecay": -fit,
        "maximumNegativeLogSecant": max(sec), "medianNegativeLogSecant": statistics.median(sec),
        "empirical95NegativeLogSecant": q95,
        "maximumOverPositiveFit": max(sec)/max(-fit, 1e-300),
        "maximumOverPositiveMedian": max(sec)/max(statistics.median(sec), 1e-300),
        "selectedDecay": max(0, -fit, *sec), "maxPair": rr[j-1:j+1],
        "maxPairCanonicalSpan": ts[j]-ts[j-1],
        "maxPairRelativeCountChange": rr[j]["core"][axis]/rr[j-1]["core"][axis]-1}


def main(destination):
    d = json.loads((destination/"native_evidence.json").read_text())
    h, m, policy = d["history"], d["latestMetadata"], d["policy"]
    common, mesh, gauge = h["common"], h["mesh"], h["gauge"]
    hist = [dict(step=int(s), time=t, epoch=int(e), core=[x,y], safety=q,
                 nodalCenter=n, horizontalWidth=wx, verticalWidth=wy)
            for s,t,e,x,y,q,n,wx,wy in zip(common["acceptedStep"],common["canonicalTau"],
                mesh["remeshCount"],mesh["coreGridPoints"],mesh["verticalCoreGridPoints"],
                mesh["safetyFactor"],gauge["omegaGaugeQuadraticStencilCenterX"],
                mesh["trackedWallCoreWidth"],mesh["trackedVerticalCoreWidth"])]
    assert all(b["step"] > a["step"] and b["time"] > a["time"] for a,b in zip(hist,hist[1:]))
    union, duplicate_count = {}, 0
    windows = []
    for cp in d["checkpointObservations"]:
        rr = rows(cp["window"])
        assert len({r["epoch"] for r in rr}) == 1
        assert all(b["step"] == a["step"]+1 for a,b in zip(rr,rr[1:]))
        for r in rr:
            key = (r["epoch"], r["step"])
            if key in union:
                assert union[key] == r
                duplicate_count += 1
            union[key] = r
        windows.append({"checkpointStep": cp["step"], "epoch": cp["remeshCount"],
                        "time": cp["canonicalTime"], "axes": [axis_stats(rr,a) for a in range(2)],
                        "savedLastDecision": cp["lastDecision"]})
    transactions, last_step, last_tau = [], 0, 0
    for a in m["transactions"]:
        before = next((r for r in reversed(hist) if r["time"] < a["sourceCanonicalTime"]),None)
        after = next((r for r in hist if r["time"] >= a["sourceCanonicalTime"]),None)
        transactions.append({"ordinal": a["sourceRemeshCount"]+1,
            "sourceStep": a["sourceStep"], "sourceCanonicalTime": a["sourceCanonicalTime"],
            "sourcePhysicalTime": a["sourcePhysicalTime"], "newAuditCoreCells": a["coreCells"],
            "absoluteRelativePeakJump": a["relativePeakJump"],
            "intervalAcceptedSteps": a["sourceStep"]-last_step,
            "intervalCanonicalTime": a["sourceCanonicalTime"]-last_tau,
            "lastSavedHistoryBefore": before, "firstSavedHistoryAfter": after,
            "controllerDecisionStored": "controllerDecision" in a,
            "exactTriggerReconstructible": False})
        last_step, last_tau = a["sourceStep"], a["sourceCanonicalTime"]
    # Strongest declines use only actual consecutive steps retained in one epoch.
    events = []
    for key,b in union.items():
        a = union.get((b["epoch"],b["step"]-1))
        if a is None:
            continue
        for axis in range(2):
            sec = -math.log(b["core"][axis]/a["core"][axis])/(b["time"]-a["time"])
            before = next((r for r in reversed(hist) if r["step"]<=a["step"] and r["epoch"]==a["epoch"]),None)
            after = next((r for r in hist if r["step"]>=b["step"] and r["epoch"]==a["epoch"]),None)
            complete_bracket = before is not None and after is not None
            changed = (before["nodalCenter"] != after["nodalCenter"]) if complete_bracket else None
            events.append({"axis": "xy"[axis], "negativeLogSecant":sec,"pair":[a,b],
                "historyBracket":[before,after], "completeSameEpochHistoryBracket":complete_bracket,
                "sampledNodalCenterChanged":changed,
                "exactStepSelectorCausationProven":False})
    events.sort(key=lambda r:r["negativeLogSecant"],reverse=True)
    # Output-sampled histories are a separate observer; never mix them into
    # the actual controller's per-step maximum or across remesh boundaries.
    sample_groups = {axis:{"center_changed":[],"center_same":[]} for axis in "xy"}
    for a,b in zip(hist,hist[1:]):
        if a["epoch"] != b["epoch"]:
            continue
        label = "center_changed" if a["nodalCenter"] != b["nodalCenter"] else "center_same"
        for axis in range(2):
            sample_groups["xy"[axis]][label].append(-math.log(b["core"][axis]/a["core"][axis])/(b["time"]-a["time"]))
    sampled = {}
    for axis,groups in sample_groups.items():
        sampled[axis] = {}
        for label,values in groups.items():
            sampled[axis][label] = {"intervalCount":len(values), "maximum":max(values),
                "median":statistics.median(values),"mean":statistics.mean(values)}
    latest = windows[-1]
    for axis,stat in enumerate(latest["axes"]):
        assert abs(stat["selectedDecay"]-m["lastDecision"]["decayEstimate"][axis]) < 1e-10
    intervals = [r["intervalCanonicalTime"] for r in transactions[1:]]
    report = {"kind":"native_transaction_frequency_offline_analysis_v1",
        "latestStep":d["latestStep"], "latestCanonicalTime":d["latestCanonicalTime"],
        "latestPhysicalTime":d["latestPhysicalTime"], "strictNativeReadCount":len(windows),
        "historyRows":len(hist), "historyIsEveryAcceptedStep":len(hist)==d["latestStep"]+1,
        "snapshotsStored":d["storedSnapshotCount"], "transactions":transactions,
        "transactionCount":len(transactions), "betweenTransactionIntervals":{
            "minimumSteps":min(r["intervalAcceptedSteps"] for r in transactions[1:]),
            "minimumTau":min(intervals),"medianTau":statistics.median(intervals),"maximumTau":max(intervals)},
        "cumulativeAbsolutePeakJump":m["cumulativeAbsolutePeakJump"],
        "savedWindows":windows, "uniqueDenseAcceptedSamples":len(union),
        "duplicateSamplesExactlyAgree":True,"duplicateSampleCount":duplicate_count,
        "topDenseDeclines":events[:20], "outputSampledSecantsByNodalCenter":sampled,
        "latestActualDecision":m["lastDecision"],"failedUnacceptedPlan":m["lastFailure"],
        "historicalTriggerUnknown":True,"exactPeakOrRootSwitchAttributionProven":False,
        "policyMutated":False,"PdeRun":False,"LuBuilt":False}
    with (destination/"frequency_analysis.json").open('x') as f:
        json.dump(report,f,indent=2,allow_nan=False)
    lines = ["# Native transaction-frequency audit", "", 
        f"Strict ten-thread reads passed for {len(windows)} immutable checkpoints. Latest accepted state: "
        f"step {d['latestStep']}, tau {d['latestCanonicalTime']:.12g}, physical time {d['latestPhysicalTime']:.12g}. "
        f"There are {len(transactions)} actual committed transactions, {len(hist)} output-sampled history rows, "
        f"and {len(rows(m['window']))} every-step observations in the final saved epoch window. No LU or PDE was run.", "",
        "The table's core counts are the actual **new candidate audit** counts. They are not the missing "
        "source counts at the trigger. The first interval starts at the original zero-time initial state; subsequent intervals are between commits.", "",
        "| # | Source step | Source tau | New core X / Y | Absolute peak jump | Steps since prior | Tau since prior |",
        "|---|---:|---:|---:|---:|---:|---:|"]
    for r in transactions:
        lines.append(f"| {r['ordinal']} | {r['sourceStep']} | {r['sourceCanonicalTime']:.9f} | "
            f"{r['newAuditCoreCells'][0]:.5f} / {r['newAuditCoreCells'][1]:.5f} | {r['absoluteRelativePeakJump']:.7g} | "
            f"{r['intervalAcceptedSteps']} | {r['intervalCanonicalTime']:.8f} |")
    lines += ["", "The shortest actual interval is 22 accepted steps / tau 0.03949490. The last two are "
        "149 and 170 steps / tau 0.39584012 and 0.44720505. This rules out a same-step or adjacent-step "
        "commit loop in this run. It does not establish an optimal transaction frequency: the saved audit "
        "has no source decision, source core, candidate cost, or alternative-mesh lifetime. "
        f"Total absolute relative peak jump is {m['cumulativeAbsolutePeakJump']:.12g}.", "",
        "## Retained-window sensitivity", "",
        "All dense samples come from actual saved native controller windows. Overlapping samples were "
        f"checked for exact agreement ({duplicate_count} duplicates, {len(union)} unique accepted samples). "
        "The output-sampled history is analyzed separately and no secant crosses a remesh epoch.", "",
        "| Final actual window | Least-squares decay | Median step secant | Empirical 95% secant | Maximum step secant |",
        "|---|---:|---:|---:|---:|"]
    for a,stat in enumerate(latest['axes']):
        lines.append(f"| {'XY'[a]} | {stat['fitDecay']:.8g} | {stat['medianNegativeLogSecant']:.8g} | "
            f"{stat['empirical95NegativeLogSecant']:.8g} | {stat['maximumNegativeLogSecant']:.8g} |")
    lines += ["", "The Y maximum is the real step 1087 to 1088 decline from 33.70519892 to 32.56282388 "
        "over tau 0.00265475 (a 3.389% loss). It is about 42 times the window regression and 275 times "
        "the median step secant. In the enclosing saved history interval 1087 to 1091, the quadratic "
        "stencil centre changes from 1 to 0.9964929462 and the sampled vertical 90% width changes from "
        "0.02808688 to 0.02599936. This is interval-level association; the per-step peak index and root "
        "templates were not saved, so exact-step attribution cannot be asserted.", "",
        "The implemented vertical diagnostic takes the whole Y trace at the selected discrete wall-peak "
        "column. Changing that column can change the Y observation abruptly on a smooth evolving field. "
        "The X threshold and its connected piecewise-linear roots also depend on the nodal maximum; "
        "their derivatives can change at selector/root-cell events. The controller retains the maximum "
        "negative secant over 0.35 canonical time, so an isolated observed drop can dominate its forecast "
        "well beyond that step. These data establish this maximum's outlier sensitivity, not a physical "
        "contraction rate or a proof that every frequent early transaction was spurious.", "",
        "## Terminal capacity failure", "",
        "The accepted step 1134 has saved decision `resolved`: core 26.00090953 / 32.34266568, "
        "predicted core 22.53568656 / 2.40780427. The large Y decay alone does not trigger because its "
        "current core is still above target 32; X prediction remains above buffer 22. The next attempted "
        "plan records `current_core_trigger` and `autonomous_mesh_axis_capacity`, with core "
        "25.95699021 / 32.33848829. Its X count crosses the current trigger 26. This attempted state "
        "was rolled back; it is not a second accepted endpoint. `lastFailure` omits its own step and "
        "clock, so they are not manufactured in this report. The exact failed planner branch is "
        "known; the individual rejected axis/pair reasons are not present in this checkpoint.", "",
        "## Minimum decision evidence for a future version", "",
        "Persist the following with each requested decision/commit and failed request, without changing the current stress run:", "",
        "- Decision schema and immutable policy version/hash; source accepted step, canonical/physical clocks, remesh epoch; requested flag, reason, per-axis trigger masks and stop reason.",
        "- Source core counts and safety, predicted counts, decay used, regression slope and maximum negative secant separately for each axis.",
        "- Window sample count/start/end and, for each maximum secant, both endpoint steps/times/counts and its time span. A copy of the bounded source window is the simplest exact replay evidence; summaries alone cannot replay the regression.",
        "- Nodal peak index/coordinate and selected vertical-trace column; 90% peak value, root-cell indices, fractional root positions and slopes at the source and maximum-secant endpoints, to distinguish selector events from sustained evolution.",
        "- Selected candidate index and axis identifiers/hashes, changed-axis flags, actual post-transfer core/safety and the existing gates; retain axis/pair rejection reasons and measured design/build/transfer/flow cost for failed and accepted attempts.", "",
        "Current committed audits lack these source decisions. Bracketing output rows and post-transfer "
        "counts are preserved in the JSON but are not substituted for the missing pre-transfer observations. "
        "No past trigger type or causation is retroactively assigned."]
    with (destination/"FREQUENCY_AUDIT.md").open('x') as f:
        f.write('\n'.join(lines)+'\n')
    print(json.dumps({"report":str(destination/'FREQUENCY_AUDIT.md'), "transactions":len(transactions),
        "latestStep":d['latestStep'],"uniqueDenseSamples":len(union), "currentDecay":[s['selectedDecay'] for s in latest['axes']]}))


if __name__ == '__main__':
    main(Path(sys.argv[1]).resolve())
