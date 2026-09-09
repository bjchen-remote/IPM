"""Root-gated, single-process-group MATLAB resource measurement. Dormant by default."""
import argparse
import hashlib
import json
import os
from pathlib import Path
import re
import signal
import subprocess
import sys
import time

GIB = 1024**3


def require(condition, message):
    if not condition:
        raise RuntimeError(message)


def validate_context(context, ram):
    require(context.get('rootGrantedWindow') is True, 'Root has not granted this LU window.')
    require(type(context.get('mainActive')) is bool, 'mainActive must be explicit.')
    reserve = context.get('mainReservedPeakBytes')
    require(type(reserve) is int and reserve >= 0, 'An explicit main-process reservation is required.')
    evidence = context.get('mainPeakEvidenceBytes')
    require(type(evidence) is int and 0 <= evidence <= reserve, 'Measured main peak must fit its reservation.')
    if context['mainActive']:
        require(evidence > 0 and bool(context.get('mainPeakEvidenceFile')), 'Active main needs actual peak evidence.')
        require(reserve > evidence, 'Main reservation must exceed observed historical peak; it is not an upper-bound proof.')
    else:
        require(reserve == 0 and evidence == 0, 'Inactive main must use a zero reservation.')
    require(reserve + 6*GIB + 3*GIB <= ram, 'Main reservation + watched6GiB + system3GiB exceeds physical RAM.')


def verify_manifest(plan):
    manifest = json.loads(Path(plan['executionManifest']).read_text())
    source = Path(plan['sourceRoot'])
    for item in manifest['source']:
        require(hashlib.sha256((source/item['path']).read_bytes()).hexdigest() == item['sha256'], 'Frozen source changed: '+item['path'])
    for item in manifest['inputs']:
        require(hashlib.sha256(Path(item['path']).read_bytes()).hexdigest() == item['sha256'], 'Input changed: '+item['path'])


def swap_used():
    text = subprocess.check_output(['/usr/sbin/sysctl', '-n', 'vm.swapusage'], text=True)
    match = re.search(r'used\s*=\s*([\d.]+)([KMG])', text)
    require(match is not None, 'Cannot read displayed swap usage.')
    return float(match.group(1)) * {'K':1024, 'M':1024**2, 'G':GIB}[match.group(2)]


def group_rss(pgid):
    listing = subprocess.check_output(['/bin/ps', '-axo', 'pgid=,rss='], text=True)
    return sum(int(parts[1])*1024 for line in listing.splitlines()
               if len(parts := line.split()) == 2 and int(parts[0]) == pgid)


def parse_time_output(text):
    result = {}
    for key, label in [('maximumRssBytes', 'maximum resident set size'),
                       ('peakFootprintBytes', 'peak memory footprint'), ('childSwaps', 'swaps')]:
        hits = re.findall(r'^\s*(\d+)\s+'+label+r'\s*$', text, re.MULTILINE)
        require(len(hits) == 1, 'Missing or ambiguous /usr/bin/time metric: '+label)
        result[key] = int(hits[0])
    return result


def terminate_owned_group(process):
    if process.poll() is not None:
        return
    os.killpg(process.pid, signal.SIGTERM)
    try:
        process.wait(timeout=10)
    except subprocess.TimeoutExpired:
        os.killpg(process.pid, signal.SIGKILL)
        process.wait(timeout=10)


def run(plan_file, context_file):
    require(sys.platform == 'darwin', 'This runner measures Darwin byte-valued time -l output.')
    plan = json.loads(Path(plan_file).read_text())
    context = json.loads(Path(context_file).read_text())
    ram = int(subprocess.check_output(['/usr/sbin/sysctl', '-n', 'hw.memsize'], text=True))
    validate_context(context, ram)
    if context['mainActive']:
        measured = parse_time_output(Path(context['mainPeakEvidenceFile']).read_text())
        require(measured['peakFootprintBytes'] == context['mainPeakEvidenceBytes'], 'Main historical peak does not match its actual time -l log.')
    require(plan['nodeCount'] == [961,321] and plan['maximumTotalNodes'] == 310000, 'Unexpected registered probe size.')
    verify_manifest(plan)
    out = Path(plan['operatorOutputDirectory'])
    log_file = Path(plan['processLog'])
    verdict_file = Path(plan['resourceVerdictFile'])
    require(not out.exists() and not log_file.exists() and not verdict_file.exists(), 'Use an unused immutable output.')
    before_swap = swap_used()
    registration = {'planFile':str(Path(plan_file).resolve()), 'context':context,
                    'physicalRamBytes':ram, 'mainReservationIsNotAProvenUpperBound':True,
                    'qualificationFootprintLimitBytes':5*GIB, 'watchdogGroupRssLimitBytes':6*GIB,
                    'systemReserveBytes':3*GIB, 'wallLimitSeconds':300, 'displayedSwapBeforeBytes':before_swap,
                    'nativeTransaction':False, 'PDETimeAdvance':False}
    Path(plan['resourceRegistrationFile']).write_text(json.dumps(registration, indent=2)+'\n')
    env = os.environ.copy()
    env['IPM_FACTOR3_RESOURCE_WINDOW'] = 'root_authorized_watched_single_factor'
    command = ['/usr/bin/time', '-l', '/Applications/MATLAB_R2026a.app/bin/matlab',
               '-batch', "run('"+plan['matlabScript']+"')"]
    start = time.monotonic()
    reason = None
    rss_peak = 0
    process = None
    try:
        with log_file.open('w') as log:
            process = subprocess.Popen(command, cwd=plan['sourceRoot'], stdout=log,
                                       stderr=subprocess.STDOUT, env=env, start_new_session=True)
            while process.poll() is None:
                rss = group_rss(process.pid)
                rss_peak = max(rss_peak, rss)
                if rss > 6*GIB:
                    reason = 'own_process_group_rss_above6GiB'
                elif time.monotonic()-start > 300:
                    reason = 'wall_time_above300s'
                elif swap_used() > before_swap:
                    reason = 'observed_system_swap_usage_increased'
                if reason:
                    terminate_owned_group(process)
                    break
                time.sleep(.25)
            code = process.wait()
    except BaseException as error:
        if process is not None:
            terminate_owned_group(process)
        reason = reason or ('watchdog_exception: '+str(error))
        code = -1
    result = {'exitCode':code, 'wallSeconds':time.monotonic()-start,
              'watchdogStopReason':reason, 'sampledOwnGroupRssPeakBytes':rss_peak,
              'includesTimeWrapperGroupRss':True, 'operatorCompleted':False,
              'resourceEnvelopePassed':False, 'nativePDEQualification':False,
              'oldCheckpointUpgradeSupported':False}
    try:
        verify_manifest(plan)
        result['allSourceAndInputSHAUnchanged'] = True
        result.update(parse_time_output(log_file.read_text()))
        result['displayedSwapAfterBytes'] = swap_used()
        report = json.loads((out/'report.json').read_text())
        result['operatorCompleted'] = (report['operatorBuilt'] and report['flowEvaluated']
                                      and report['factorReleasedBeforeOutput'] and report['operatorFlowFinite'])
        result['resourceEnvelopePassed'] = (code == 0 and reason is None and result['operatorCompleted']
                                           and result['peakFootprintBytes'] <= 5*GIB
                                           and result['childSwaps'] == 0
                                           and result['displayedSwapAfterBytes'] <= before_swap)
    except Exception as error:
        result['verdictFailure'] = str(error)
    verdict_file.write_text(json.dumps(result, indent=2)+'\n')
    print(json.dumps(result, indent=2))
    return 0 if result['resourceEnvelopePassed'] else 1


def self_test():
    valid = {'rootGrantedWindow':True, 'mainActive':True, 'mainReservedPeakBytes':5*GIB,
             'mainPeakEvidenceBytes':4*GIB, 'mainPeakEvidenceFile':'synthetic_fixture_only'}
    validate_context(valid, 16*GIB)
    tests = [dict(valid, rootGrantedWindow=False), dict(valid, mainReservedPeakBytes=8*GIB),
             dict(valid, mainPeakEvidenceBytes=None), dict(valid, mainReservedPeakBytes=4*GIB),
             dict(valid, mainActive=False), dict(valid, mainActive=1)]
    for case in tests:
        try:
            validate_context(case, 16*GIB)
        except RuntimeError:
            pass
        else:
            raise AssertionError('Invalid resource context passed.')
    metrics = parse_time_output('123 maximum resident set size\n456 peak memory footprint\n0 swaps\n')
    require(metrics == {'maximumRssBytes':123, 'peakFootprintBytes':456, 'childSwaps':0}, 'Metric parser mismatch.')
    try:
        parse_time_output('123 maximum resident set size\n0 swaps\n')
    except RuntimeError:
        pass
    else:
        raise AssertionError('Missing footprint silently accepted.')
    print(json.dumps({'pureSelfTestsPassed':9, 'spawnedProcesses':False, 'ranMATLAB':False}))
    return 0


if __name__ == '__main__':
    parser = argparse.ArgumentParser()
    parser.add_argument('--self-test', action='store_true')
    parser.add_argument('--plan')
    parser.add_argument('--context')
    args = parser.parse_args()
    if args.self_test:
        sys.exit(self_test())
    require(args.plan is not None and args.context is not None, 'An explicit plan and root resource context are required.')
    sys.exit(run(args.plan, args.context))
