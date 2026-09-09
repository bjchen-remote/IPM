#!/usr/bin/env python3
"""Freeze minimal research call chains, changing only WENO entry routing."""
import argparse
import hashlib
import json
import re
import shutil
from pathlib import Path

CHAIN = ['evolve.advance', 'evolve.stepSsprk54', 'evolve.stepRk', 'evolve.flow',
         'evolve.rhs', 'evolve.isotropicGauge', 'evolve.assembleRhs',
         'field.transport', 'field.weno5FluxDerivative']


def build(native_source, research_directory, output):
    native_source, research_directory, output = map(Path, [native_source, research_directory, output])
    assert not output.exists()
    output.mkdir()
    shutil.copytree(native_source / '+ipm', output / '+ipm')
    originals = {}
    for entry in CHAIN:
        domain, name = entry.split('.')
        path = native_source / '+ipm' / ('+' + domain) / (name + '.m')
        originals[entry] = path.read_text()
    for block in [16, 32]:
        package = f'perflab_b{block}'
        for entry in CHAIN:
            domain, name = entry.split('.')
            destination = output / ('+' + package) / ('+' + domain) / (name + '.m')
            destination.parent.mkdir(parents=True, exist_ok=True)
            if entry == 'field.weno5FluxDerivative':
                text = (research_directory / 'ipm_perflab_weno5_candidate.m').read_text()
                text = text.replace('ipm_perflab_weno5_candidate(', 'weno5FluxDerivative(', 1)
                text = text.replace('blockSize = 16;', f'blockSize = {block};', 1)
            else:
                text = originals[entry]
                for target in CHAIN:
                    text = re.sub(r'\bipm\.' + re.escape(target) + r'\b', package + '.' + target, text)
            destination.write_text(text)
    shutil.copy2(research_directory / 'ipm_perflab_weno_actual_abba.m', output)
    manifest = {'kind': 'minimal_nine_function_weno_routing_probe',
                'nativeSource': str(native_source.resolve()), 'chain': CHAIN,
                'candidateBlockSizes': [16, 32], 'candidateSolveOrRestoreCopied': False,
                'files': {str(p.relative_to(output)): hashlib.sha256(p.read_bytes()).hexdigest()
                          for p in sorted(output.rglob('*.m'))}}
    (output / 'source_manifest.json').write_text(json.dumps(manifest, indent=2))
    print(json.dumps({'source': str(output.resolve()), 'files': len(manifest['files'])}))


if __name__ == '__main__':
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('native_source')
    parser.add_argument('research_directory')
    parser.add_argument('output')
    args = parser.parse_args()
    build(args.native_source, args.research_directory, args.output)
