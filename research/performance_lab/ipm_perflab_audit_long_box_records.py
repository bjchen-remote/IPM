#!/usr/bin/env python3
"""Read saved scalars after the actual-endpoint three-box native audit passes.

This is supplemental: no signature validation, operator build or PDE step.
"""
import argparse
import glob
import hashlib
import importlib.util
import json
import math
import sys
from pathlib import Path


def sha256(file):
    digest = hashlib.sha256()
    with Path(file).open('rb') as stream:
        for block in iter(lambda: stream.read(1024 * 1024), b''):
            digest.update(block)
    return digest.hexdigest()


def run(directory, reader_file):
    sys.dont_write_bytecode = True
    directory = Path(directory).resolve()
    reader_file = Path(reader_file).resolve()
    spec = importlib.util.spec_from_file_location('ipm_long_box_scalar_reader', reader_file)
    reader = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(reader)
    library = sorted(glob.glob('/Applications/MATLAB_*.app/bin/maca64/libhdf5.*.dylib'), reverse=True)[0]
    report = json.loads((directory / 'dynamic/three_box_report.json').read_text())
    assert report['completed'] and all(b['audit']['passed'] for b in report['branches'])
    rows = []
    for branch in report['branches']:
        hdf = reader.HDF5(library, branch['checkpointFile'])
        try:
            def vector(path):
                value, _ = hdf.numeric('/checkpoint/payload/log/history/' + path)
                return value if isinstance(value, list) else [value]
            steps = vector('common/acceptedStep')
            values = {name: vector(path) for name, path in {
                'coreX': 'mesh/coreGridPoints', 'coreY': 'mesh/verticalCoreGridPoints',
                'safety': 'mesh/safetyFactor', 'poissonResidual': 'common/poissonResidual'}.items()}
            assert all(len(v) == len(steps) for v in values.values())
            assert all(math.isfinite(z) for v in [steps, *values.values()] for z in v)
            initial_step = steps[-1] - branch['segmentAcceptedSteps']
            indices = [i for i, step in enumerate(steps) if step >= initial_step]
            assert indices and steps[-1] >= initial_step
            rows.append({'label': branch['label'], 'checkpointFile': branch['checkpointFile'],
                         'initialStep': initial_step, 'terminalStep': steps[-1],
                         'entireHistorySavedRecords': len(steps), 'segmentSavedRecords': len(indices),
                         'minimumSegmentRecordedCore': [min(values[n][i] for i in indices) for n in ['coreX', 'coreY']],
                         'maximumSegmentRecordedSafety': max(values['safety'][i] for i in indices),
                         'maximumSegmentRecordedPoissonResidual': max(values['poissonResidual'][i] for i in indices),
                         'maximumEntireHistoryRecordedPoissonResidual': max(values['poissonResidual']),
                         'nativeFullHistoryPreviouslyExact': branch['audit']['nativeFullHistoryExact'],
                         'allReadValuesFinite': True})
        finally:
            hdf.close()
    manifest = json.loads((directory / 'source_manifest.json').read_text())
    changed_sources = [r['path'] for r in manifest['sourceFiles']
                       if sha256(directory / r['path']) != r['sha256']]
    changed_inputs = [r['path'] for r in manifest['inputFiles'] if sha256(r['path']) != r['sha256']]
    assert not changed_sources and not changed_inputs
    audit = {'kind': 'supplemental_actual_endpoint_three_box_saved_record_audit', 'rows': rows,
             'sourceFilesUnchanged': True, 'inputFilesUnchanged': True,
             'sourceFileCount': len(manifest['sourceFiles']), 'inputFileCount': len(manifest['inputFiles']),
             'readerFile': str(reader_file), 'readerSHA256': sha256(reader_file),
             'signatureValidatedByThisRead': False, 'operatorBuilt': False, 'pdeAdvanced': False,
             'interpretation': 'Saved finite scalar observations after strict native validation and complete history pairing. '
                               'No claim about unrecorded RK stages; Poisson residual is not a trustedMask condition.'}
    destination = directory / 'supplemental_saved_record_audit.json'
    assert not destination.exists()
    destination.write_text(json.dumps(audit, indent=2))
    print(json.dumps(audit, indent=2))


if __name__ == '__main__':
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('directory')
    parser.add_argument('reader_file')
    args = parser.parse_args()
    run(args.directory, args.reader_file)
