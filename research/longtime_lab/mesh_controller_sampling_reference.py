#!/usr/bin/env python3
"""Bounded deterministic sampling model, not an IPM implementation or run.

Only synthetic observations are used. It demonstrates actual-sample retention,
conservative rate envelopes and serialized split replay. No solver data change.
"""
import copy
import json
import math
from collections import deque
from pathlib import Path
import argparse

W = .35
RAW = 32
BINS = 480
H = W / (BINS - 2)  # leave room for both partial bins/rounding boundaries


def update(memory, row):
    # row = native canonical time, accepted step, remesh epoch, coreX/coreY, safety
    t, step, epoch, *_ = row
    if memory is None or memory['epoch'] != epoch:
        memory = {'kind': 'research_model_only', 'epoch': epoch, 'origin': t,
                  'raw': [], 'bins': []}
    raw, bins = memory['raw'], memory['bins']
    previous = raw[-1] if raw else None
    if previous is not None and step == previous[1]:
        assert row == previous  # idempotent repeated observation; no new sample
        return memory
    assert previous is None or (t > previous[0] and step > previous[1])
    raw.append(row)
    memory['raw'] = raw = [v for v in raw[-RAW:] if v[0] >= t-W]
    bucket = math.floor((t-memory['origin']) / H)
    assert bucket < 2**52
    if not bins or bins[-1]['index'] != bucket:
        bins.append({'index': bucket, 'last': row, 'decay': [0., 0.], 'pairs': [None, None]})
    b = bins[-1]
    b['last'] = row
    if previous is not None:
        for axis in range(2):
            decay = -(math.log(row[3+axis])-math.log(previous[3+axis]))/(t-previous[0])
            if decay > b['decay'][axis]:
                b['decay'][axis] = decay
                b['pairs'][axis] = [previous, row]
    memory['bins'] = [b for b in bins if b['last'][0] >= t-W]
    assert len(memory['bins']) <= BINS
    assert len(samples(memory)) <= RAW+BINS
    return memory


def samples(memory):
    actual = {v[1]: v for v in memory['raw']}
    actual.update({b['last'][1]: b['last'] for b in memory['bins']})
    return [actual[k] for k in sorted(actual)]


def decay_bound(memory):
    return [max((b['decay'][axis] for b in memory['bins']), default=0.) for axis in range(2)]


def stream():
    t = 0.
    yield [t, 0, 0, 40., 44., .3]
    for step in range(1, 25001):
        dt = .0005 if step <= 1000 else (1e-5 if step <= 21000 else 1e-7)
        t += dt
        epoch = int(step >= 18000)
        yield [t, step, epoch, 40*math.exp(-.5*t-.002*math.sin(30*t)),
               44*math.exp(-.35*t+.001*math.sin(100*t)), .3]


def run(out):
    out.mkdir(parents=True, exist_ok=False)
    memory = split = None
    full = deque()
    actual = {}
    maximum_rows = maximum_bins = 0
    maximum_excess = 0.
    checks = 0
    boundary_checks = 0
    splits = {501, 1000, 9876, 17999, 18000, 21987}
    for row in stream():
        actual[row[1]] = row
        if full and full[-1][2] != row[2]:
            full.clear()
        full.append(row)
        while full[0][0] < row[0]-W:
            full.popleft()
        memory = update(memory, row)
        split = update(split, row)
        if row[1] in splits:
            split = json.loads(json.dumps(split, separators=(',', ':')))
            boundary_checks += 1
        assert split == memory
        assert samples(memory)[-1] == row
        assert all(actual[v[1]] == v and v[0] >= row[0]-W for v in samples(memory))
        maximum_rows = max(maximum_rows, len(samples(memory)))
        maximum_bins = max(maximum_bins, len(memory['bins']))
        if row[1] % 137 == 0 and len(full) >= 4:
            oracle = [0., 0.]
            previous = full[0]
            for current in list(full)[1:]:
                for axis in range(2):
                    oracle[axis] = max(oracle[axis], -(math.log(current[3+axis])-math.log(previous[3+axis]))/(current[0]-previous[0]))
                previous = current
            bounded = decay_bound(memory)
            assert all(a >= b for a, b in zip(bounded, oracle))
            maximum_excess = max(maximum_excess, max(a-b for a,b in zip(bounded,oracle)))
            checks += 1
    repeated = copy.deepcopy(memory)
    assert update(repeated, memory['raw'][-1]) == memory
    report = {'kind': 'bounded_sampling_research_model_v1', 'productionImplementation': False,
              'noPDE': True, 'samplesProcessed': 25001, 'rawBudget': RAW, 'bucketBudget': BINS,
              'recordBudget': RAW+BINS, 'window': W, 'bucketWidth': H,
              'maximumRetainedRows': maximum_rows, 'maximumRetainedBuckets': maximum_bins,
              'rateEnvelopeChecks': checks, 'allRateEnvelopesConservative': True,
              'maximumExcessDecay': maximum_excess, 'serializedSplitBoundaries': boundary_checks,
              'everyStepMemoryIdenticalAfterSplit': True, 'sameObservationIdempotent': True,
              'allRetainedPointsActual': True, 'remeshEpochResetTested': True,
              'finalSerializedBytes': len(json.dumps(memory).encode()), 'allPassed': True}
    (out/'report.json').write_text(json.dumps(report, indent=2))
    (out/'final_model_memory.json').write_text(json.dumps(memory, indent=2))
    print(json.dumps(report))


if __name__ == '__main__':
    parser = argparse.ArgumentParser()
    parser.add_argument('output_directory', type=Path)
    run(parser.parse_args().output_directory)
