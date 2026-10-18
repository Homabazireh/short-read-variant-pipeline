#!/usr/bin/env python3
"""Compile local DSL2 workflows with Nextflow inspect; execute no bioinformatics tasks."""
from __future__ import annotations

import json
from pathlib import Path
import shutil
import subprocess
import tempfile

from common import CONTAINERS


def main():
    root = Path(__file__).resolve().parents[1]
    nextflow = shutil.which('nextflow')
    if not nextflow:
        raise SystemExit('Install pinned Nextflow and put it on PATH first')
    with tempfile.TemporaryDirectory(prefix='pipeline-compile-') as name:
        temp = Path(name)
        fasta = temp / 'reference.fa'
        fasta.write_text('>chr20\nACGT\n')
        fai = temp / 'reference.fa.fai'
        fai.write_text('chr20\t4\t7\t4\t5\n')
        bed = temp / 'regions.bed'
        bed.write_text('chr20\t0\t4\n')
        cache = temp / 'cache'
        cache.mkdir()
        parameters = {
            'jobs': [], 'samples': [], 'reference_fasta': str(fasta),
            'reference_fai': str(fai), 'reference_sha256': 'compile-check-only',
            'reference_contigs': {'chr20': 4}, 'outdir': str(temp / 'output'),
            'containers': CONTAINERS, 'regions': str(bed),
            'annotation': {'vep_cache': str(cache), 'clinvar': None, 'gnomad': None},
        }
        configs = []
        for kind in ('annotation', 'benchmark'):
            configs.append((root / 'workflows' / kind / 'main.nf',
                            root / 'workflows' / kind / 'nextflow.config', parameters))
        if (root / 'main.nf').exists():
            configs.append((root / 'main.nf', root / 'nextflow.config', parameters))
        sv = root / 'configs' / 'shortread_sv.nf'
        if sv.exists():
            samples = temp / 'alignments.csv'
            samples.write_text('sample,bam,bai\n')
            configs.append((sv, None, {**parameters, 'input': str(samples),
                                      'fasta': str(fasta), 'fai': str(fai)}))
        for i, (workflow, config, params) in enumerate(configs):
            params_path = temp / f'params-{i}.json'
            params_path.write_text(json.dumps(params))
            command = [nextflow]
            if config:
                command += ['-C', str(config)]
            command += ['inspect', str(workflow), '-params-file', str(params_path), '-format', 'json']
            completed = subprocess.run(command, cwd=temp, capture_output=True, text=True)
            if completed.returncode:
                raise SystemExit(f'{workflow.name} compilation failed:\n{completed.stdout}\n{completed.stderr}')
            print(f'Compiled: {workflow.relative_to(root)}')
    return 0


if __name__ == '__main__':
    raise SystemExit(main())
