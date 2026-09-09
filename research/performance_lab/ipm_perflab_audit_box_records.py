#!/usr/bin/env python3
"""Supplemental scalar reads only, after the native five-case protocol ends.

No MATLAB execution, signature validation, PDE or mutation. The native
protocol has already validated and paired every terminal checkpoint.
"""
import argparse
import glob
import hashlib
import importlib.util
import json
import math
import sys
from pathlib import Path


def run(directory, reader_file):
    sys.dont_write_bytecode = True
    directory = Path(directory).resolve()
    reader_file = Path(reader_file).resolve()
    spec = importlib.util.spec_from_file_location('ipm_box_scalar_reader', reader_file)
    reader = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(reader)
    library = sorted(glob.glob('/Applications/MATLAB_*.app/bin/maca64/libhdf5.*.dylib'), reverse=True)[0]
    report = json.loads((directory / 'dynamic/protocol_report.json').read_text())
    assert report['completed']
    rows = []
    for branch in report['branches']:
        assert branch['audit']['passed']
        hdf = reader.HDF5(library, branch['checkpointFile'])
        try:
            prefix = '/checkpoint/payload/log/history/'
            def vector(path):
                value, _ = hdf.numeric(prefix + path)
                return value if isinstance(value, list) else [value]
            steps = vector('common/acceptedStep')
            corex = vector('mesh/coreGridPoints')
            corey = vector('mesh/verticalCoreGridPoints')
            safety = vector('mesh/safetyFactor')
            residual = vector('common/poissonResidual')
            assert all(len(v) == len(steps) for v in [corex, corey, safety, residual])
            assert all(math.isfinite(z) for v in [steps, corex, corey, safety, residual] for z in v)
            initial_step = report['registration']['parentStep'] if branch['label'] == 'native_continuation' else 0
            indices = [i for i, step in enumerate(steps) if step >= initial_step]
            assert indices
            rows.append({'label': branch['label'], 'checkpointFile': branch['checkpointFile'],
                         'allSavedRecords': len(steps), 'branchSavedRecords': len(indices),
                         'minimumBranchRecordedCore': [min(corex[i] for i in indices), min(corey[i] for i in indices)],
                         'maximumBranchRecordedSafety': max(safety[i] for i in indices),
                         'maximumBranchRecordedPoissonResidual': max(residual[i] for i in indices),
                         'maximumEntireHistoryRecordedPoissonResidual': max(residual),
                         'nativePairingPreviouslyPassed': branch['audit']['nativeTerminalPaired'],
                         'allReadValuesFinite': True})
        finally:
            hdf.close()
    manifest = json.loads((directory / 'source_manifest.json').read_text())
    changed = [name for name, digest in manifest['files'].items()
               if hashlib.sha256((directory / 'source' / name).read_bytes()).hexdigest() != digest]
    checkpoint_unchanged = hashlib.sha256(Path(manifest['sourceCheckpoint']).read_bytes()).hexdigest() == manifest['sourceCheckpointSHA256']
    assert not changed and checkpoint_unchanged
    audit = {'kind': 'supplemental_post_run_saved_scalar_audit', 'rows': rows,
             'sourceFilesUnchanged': True, 'sourceCheckpointBytesUnchanged': True,
             'sourceFileCount': len(manifest['files']),
             'scalarReaderFile': str(reader_file),
             'scalarReaderSHA256': hashlib.sha256(reader_file.read_bytes()).hexdigest(),
             'signatureValidatedByThisRead': False, 'operatorBuilt': False, 'pdeAdvanced': False,
             'interpretation': 'Finite saved scalar records, after native checkpoint validation and pairing. '
                               'Not unrecorded RK stages, and Poisson residual is not a trustedMask condition.'}
    destination = directory / 'supplemental_saved_record_audit.json'
    assert not destination.exists()
    destination.write_text(json.dumps(audit, indent=2))
    print(json.dumps({'rows': rows, 'sourceFilesUnchanged': True, 'sourceCheckpointBytesUnchanged': True}, indent=2))


if __name__ == '__main__':
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('directory')
    parser.add_argument('reader_file')
    arguments = parser.parse_args()
    run(arguments.directory, arguments.reader_file)
