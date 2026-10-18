#!/usr/bin/env python3
"""Standard-library preflight and provenance shared by annotation/benchmarking."""
import argparse
import atexit
import csv
import gzip
import hashlib
import json
import os
import re
import shlex
import shutil
import signal
import socket
import subprocess
import time
from datetime import datetime, timezone
from pathlib import Path

CONTAINERS = {
    'bcftools': 'quay.io/biocontainers/bcftools@sha256:065115d2fbdcba9b9a9d360fe9f024ba829de2e64309a637a274ed47e0af17db',
    'vep': 'ensemblorg/ensembl-vep@sha256:f354dd8d09073e4d943acbbd02f5eb234a9d9e9d444371c1c349910f2123de11',
    'happy': 'jmcdani20/hap.py@sha256:0812d37e7210011e407914deb3da2094ec8258077b21d3211e694e5ca303b489',
    'truvari': 'quay.io/biocontainers/truvari@sha256:022be7a2dfe5decf1aff8d007a914f56ed68d0bc224888e1764820922bc11935',
}
VERSIONS = {'bcftools': '1.22', 'vep': '116.0', 'happy': '0.3.12', 'rtg_tools': '3.10.1', 'truvari': '5.3.0'}


def sha256(path):
    h = hashlib.sha256()
    with open(path, 'rb') as handle:
        for block in iter(lambda: handle.read(8 * 1024 * 1024), b''):
            h.update(block)
    return h.hexdigest()


def text_open(path):
    with open(path, 'rb') as handle:
        compressed = handle.read(2) == b'\x1f\x8b'
    return gzip.open(path, 'rt') if compressed else open(path)


def local_path(value, base=None, directory=False):
    if not value:
        raise ValueError('Missing required file/directory path')
    path = Path(value).expanduser()
    if not path.is_absolute():
        path = (Path(base) if base else Path.cwd()) / path
    path = path.resolve()
    if any(c in str(path) for c in '\n\r\t'):
        raise ValueError('Paths may not contain tabs or newlines')
    if not (path.is_dir() if directory else path.is_file()):
        raise ValueError(f'Missing {"directory" if directory else "file"}: {path}')
    return str(path)


def checked_file(value, expected, base=None):
    path = local_path(value, base)
    if not re.fullmatch(r'[0-9a-fA-F]{64}', expected or ''):
        raise ValueError(f'A real SHA256 checksum is required for {path}')
    if sha256(path) != expected.lower():
        raise ValueError(f'SHA256 mismatch: {path}')
    return path


class Reference:
    """Read an uncompressed indexed FASTA without third-party dependencies."""
    def __init__(self, path):
        self.path = path
        if str(path).endswith('.gz'):
            raise ValueError('Use an uncompressed reference FASTA with .fai index')
        self.index = {}
        with open(local_path(str(path) + '.fai')) as handle:
            for line in handle:
                c, n, offset, bases, width, *_ = line.rstrip().split('\t')
                values = tuple(map(int, (n, offset, bases, width)))
                if values[0] < 1 or values[2] < 1 or values[3] < values[2] or c in self.index:
                    raise ValueError('Invalid/duplicate FASTA index entry')
                self.index[c] = values
        if not self.index:
            raise ValueError('Empty FASTA index')
        self.handle = open(path, 'rb')

    def bases(self, chrom, start, length):
        if chrom not in self.index:
            raise ValueError(f'Contig {chrom!r} is absent from the reference')
        n, offset, line_bases, line_width = self.index[chrom]
        if start < 0 or length < 1 or start + length > n:
            raise ValueError(f'Variant outside reference: {chrom}:{start + 1}')
        result = []
        while length:
            take = min(length, line_bases - start % line_bases)
            self.handle.seek(offset + start // line_bases * line_width + start % line_bases)
            result.append(self.handle.read(take).decode('ascii'))
            start += take
            length -= take
        return ''.join(result).upper()


def load_resources(path):
    path = Path(path).resolve()
    with open(path) as handle:
        resources = json.load(handle)
    if resources.get('schema_version') != 1 or resources.get('genome_build') != 'GRCh38':
        raise ValueError('Resources require schema_version=1 and genome_build=GRCh38')
    resources['reference_fasta'] = checked_file(resources.get('reference_fasta'), resources.get('reference_sha256'), path.parent)
    resources['reference_fai'] = local_path(resources['reference_fasta'] + '.fai')
    resources['reference_fai_sha256'] = sha256(resources['reference_fai'])
    resources['manifest_sha256'] = sha256(path)
    return resources, Reference(resources['reference_fasta'])


def validate_vcf(path, sample, reference):
    """Validate sample, coordinates, declared build/lengths, order and every REF."""
    count, header, seen_header = 0, [], False
    contigs, seen_contigs, previous = {}, set(), None
    with text_open(path) as handle:
        for line in handle:
            if line.startswith('##'):
                if line.lower().startswith(('##reference=', '##assembly=')):
                    if re.search(r'grch37|hg19|chm13|t2t', line, re.I):
                        raise ValueError(f'Incompatible declared reference in {path}: {line.strip()}')
                if line.startswith('##contig='):
                    name = re.search(r'(?:<|,)ID=([^,>]+)', line)
                    size = re.search(r'(?:<|,)length=(\d+)', line, re.I)
                    if name:
                        chrom = name.group(1)
                        if chrom not in reference.index:
                            raise ValueError(f'VCF header contig {chrom} absent from reference: {path}')
                        if size and int(size.group(1)) != reference.index[chrom][0]:
                            raise ValueError(f'VCF/reference contig length mismatch: {chrom} in {path}')
                        contigs[chrom] = True
                continue
            if line.startswith('#CHROM\t'):
                header = line.rstrip().split('\t')
                if header[9:] != [sample]:
                    raise ValueError(f'{path}: expected single sample {sample!r}, found {header[9:]}')
                seen_header = True
                continue
            if line.startswith('#') or not line.strip():
                continue
            if not seen_header:
                raise ValueError(f'VCF #CHROM header missing: {path}')
            fields = line.rstrip().split('\t')
            if len(fields) != len(header):
                raise ValueError(f'VCF field count mismatch: {path}')
            chrom, pos, _, ref, alt = fields[:5]
            pos = int(pos)
            if not re.fullmatch(r'[ACGTNacgtn]+', ref):
                raise ValueError(f'Invalid VCF REF: {chrom}:{pos}')
            if reference.bases(chrom, pos - 1, len(ref)) != ref.upper():
                raise ValueError(f'REF mismatch: {chrom}:{pos} in {path}')
            if previous:
                if chrom == previous[0] and pos < previous[1] or chrom != previous[0] and chrom in seen_contigs:
                    raise ValueError(f'Unsorted VCF: {path}')
            seen_contigs.add(chrom)
            previous = (chrom, pos)
            if 'GT' not in fields[8].split(':'):
                raise ValueError(f'Genotype GT missing: {chrom}:{pos} in {path}')
            count += 1
    if not seen_header:
        raise ValueError(f'VCF #CHROM header missing: {path}')
    return {'path': str(path), 'sha256': sha256(path), 'records': count, 'sample': sample,
            'reference_check': 'all_REF_alleles_and_declared_contig_lengths'}


def read_bed(path, reference):
    intervals = {}
    with text_open(path) as handle:
        for line in handle:
            if not line.strip() or line.startswith(('#', 'track', 'browser')):
                continue
            fields = line.split()
            if len(fields) < 3:
                raise ValueError(f'Invalid BED line in {path}')
            chrom, start, end = fields[0], int(fields[1]), int(fields[2])
            if chrom not in reference.index or not 0 <= start < end <= reference.index[chrom][0]:
                raise ValueError(f'BED/reference incompatibility in {path}: {line.strip()}')
            intervals.setdefault(chrom, []).append((start, end))
    for chrom in intervals:
        merged = []
        for start, end in sorted(intervals[chrom]):
            if merged and start <= merged[-1][1]:
                merged[-1] = (merged[-1][0], max(end, merged[-1][1]))
            else:
                merged.append((start, end))
        intervals[chrom] = merged
    if not intervals:
        raise ValueError(f'Empty BED: {path}')
    return intervals


def intersect_beds(first, second):
    result = {}
    for chrom in first.keys() & second.keys():
        i = j = 0
        a, b = first[chrom], second[chrom]
        while i < len(a) and j < len(b):
            start, end = max(a[i][0], b[j][0]), min(a[i][1], b[j][1])
            if start < end:
                result.setdefault(chrom, []).append((start, end))
            if a[i][1] < b[j][1]:
                i += 1
            else:
                j += 1
    return result


def read_calls(path, benchmark=False):
    path = Path(path).resolve()
    required = {'sample', 'label', 'vcf', 'variant_type'}
    if benchmark:
        required |= {'truth_vcf', 'truth_bed'}
    with open(path, newline='') as handle:
        reader = csv.DictReader(handle)
        if not required <= set(reader.fieldnames or []):
            raise ValueError(f'Input CSV requires columns: {", ".join(sorted(required))}')
        rows = list(reader)
    if not rows:
        raise ValueError('Input CSV contains no callsets')
    seen = set()
    for row in rows:
        if None in row or any(row.get(k) is None for k in required):
            raise ValueError('Malformed CSV row')
        row.update({k: v.strip() for k, v in row.items()})
        for key in ('sample', 'label'):
            if not re.fullmatch(r'[A-Za-z0-9][A-Za-z0-9_.-]*', row[key]):
                raise ValueError(f'{key} must contain only letters, digits, dot, underscore, hyphen')
        if row['variant_type'] not in ('small', 'sv'):
            raise ValueError('variant_type must be small or sv')
        row['key'] = '__'.join(row[k] for k in ('sample', 'label', 'variant_type'))
        if row['key'] in seen:
            raise ValueError(f'Duplicate callset: {row["key"]}')
        seen.add(row['key'])
        for key in ('vcf', 'truth_vcf', 'truth_bed', 'evaluation_regions', 'stratifications', 'sv_policy'):
            if row.get(key):
                row[key] = local_path(row[key], path.parent)
    return rows


def argument_parser(description):
    p = argparse.ArgumentParser(description=description)
    p.add_argument('--input', required=True, help='Callset CSV; paths relative to its directory')
    p.add_argument('--resources', required=True, help='Checksum-pinned resources JSON')
    p.add_argument('--outdir', required=True)
    p.add_argument('--profile', choices=['apptainer', 'docker'], default='apptainer')
    p.add_argument('--executor', choices=['local', 'slurm'], default='local')
    p.add_argument('--cpus', type=int, default=8)
    p.add_argument('--memory-gb', type=int, default=32)
    p.add_argument('--resume', action='store_true')
    p.add_argument('--dry-run', action='store_true', help='Validate all inputs and print invocation without launching tools')
    return p


def validate_output_destination(outdir, resume):
    outdir = Path(outdir)
    if (outdir / '.pipeline.lock').exists():
        raise ValueError('Output directory is locked by another run; inspect .pipeline.lock before recovering a stale run')
    existing = (outdir / 'run.params.json').is_file()
    if resume and not existing:
        raise ValueError('--resume requires an existing run.params.json')
    if not resume and outdir.exists() and any(outdir.iterdir()):
        raise ValueError('Output directory is not empty; choose a new directory or resume the matching run')


def acquire_run_lock(outdir):
    lock = Path(outdir) / '.pipeline.lock'
    token = json.dumps({'host': socket.gethostname(), 'pid': os.getpid(), 'started_unix': time.time()})
    try:
        descriptor = os.open(lock, os.O_WRONLY | os.O_CREAT | os.O_EXCL, 0o600)
    except FileExistsError:
        raise ValueError('Another process acquired the output-directory lock') from None
    with os.fdopen(descriptor, 'w') as handle:
        handle.write(token + '\n')
    def release():
        try:
            if lock.read_text().strip() == token:
                lock.unlink()
        except FileNotFoundError:
            pass
    atexit.register(release)
    return release


def run_managed(command, cwd, cleanup_timeout=10):
    """Keep descendants in a managed group; finish shutdown before releasing locks."""
    process = subprocess.Popen(command, cwd=cwd, start_new_session=True)
    try:
        code = process.wait()
        if code:
            raise subprocess.CalledProcessError(code, command)
    except BaseException:
        try:
            os.killpg(process.pid, signal.SIGTERM)
        except ProcessLookupError:
            pass
        deadline = time.monotonic() + cleanup_timeout
        while time.monotonic() < deadline:
            process.poll()
            try:
                os.killpg(process.pid, 0)
            except ProcessLookupError:
                break
            time.sleep(0.05)
        try:
            os.killpg(process.pid, signal.SIGKILL)
        except ProcessLookupError:
            pass
        process.wait()
        raise


def run_workflow(kind, args, params, provenance):
    if args.cpus < 1 or args.memory_gb < 1:
        raise ValueError('cpus and memory-gb must be positive')
    root = Path(__file__).resolve().parents[1]
    outdir = Path(args.outdir).resolve()
    implementation = {str(p.relative_to(root)): sha256(p) for p in [
        root / 'scripts/common.py', root / f'scripts/{"annotate" if kind == "annotation" else "benchmark"}.py',
        root / 'workflows' / kind / 'main.nf', root / 'workflows' / kind / 'nextflow.config']}
    fingerprint = hashlib.sha256(json.dumps(provenance, sort_keys=True).encode()).hexdigest()
    params.update(outdir=str(outdir), max_cpus=args.cpus, memory_gb=args.memory_gb,
                  executor_name=args.executor, containers=CONTAINERS,
                  implementation_sha256=implementation, input_fingerprint=fingerprint)
    command = ['nextflow', '-C', str(root / 'workflows' / kind / 'nextflow.config'),
               'run', str(root / 'workflows' / kind / 'main.nf'),
               '-profile', args.profile, '-params-file', str(outdir / 'run.params.json'),
               '-work-dir', str(outdir / 'work'), '-with-trace', str(outdir / 'trace.tsv'),
               '-with-report', str(outdir / 'execution-report.html'),
               '-with-timeline', str(outdir / 'timeline.html')]
    if args.resume:
        command.append('-resume')
    provenance.update(schema_version=1, workflow=kind, command=command, containers=CONTAINERS, versions=VERSIONS,
                      profile=args.profile, executor=args.executor, cpus=args.cpus, memory_gb=args.memory_gb,
                      implementation_sha256=implementation, input_fingerprint=fingerprint,
                      hostname=socket.gethostname(), tested_nextflow_version='25.10.4')
    if args.dry_run:
        print(json.dumps({'command': shlex.join(command), 'parameters': params,
                          'preflight': provenance, 'execution': 'not_run'}, indent=2))
        return None
    if shutil.which('nextflow') is None or shutil.which(args.profile) is None:
        raise ValueError(f'Install Nextflow and {args.profile}, then place both on PATH')
    version = subprocess.run(['nextflow', '-version'], check=True, capture_output=True, text=True)
    if not re.search(r'\bversion\s+25\.10\.4(?=\s|$)', version.stdout + version.stderr):
        raise ValueError('This release requires tested Nextflow 25.10.4; set NXF_VER=25.10.4')
    validate_output_destination(outdir, args.resume)
    outdir.mkdir(parents=True, exist_ok=True)
    release_lock = acquire_run_lock(outdir)
    previous_path = outdir / 'run.params.json'
    if previous_path.exists():
        with previous_path.open() as handle:
            previous = json.load(handle)
        if not args.resume:
            raise ValueError('Output directory already contains a run; use --resume or a new directory')
        if previous != params:
            release_lock()
            raise ValueError('Resume parameters changed; choose a new output directory')
        archive = outdir / 'execution-history' / datetime.now(timezone.utc).strftime('%Y%m%dT%H%M%S.%fZ')
        archive.mkdir(parents=True)
        for name in ('trace.tsv', 'execution-report.html', 'timeline.html', 'provenance.json'):
            if (outdir / name).exists():
                shutil.copy2(outdir / name, archive / name)
        provenance['previous_execution_archive'] = str(archive)
    for name, content in params.get('generated_files', {}).items():
        target = Path(name)
        if not target.resolve().is_relative_to(outdir):
            release_lock()
            raise ValueError('Generated input path escaped output directory')
        target.parent.mkdir(parents=True, exist_ok=True)
        if not target.exists() or target.read_text() != content:
            target.write_text(content)
    previous_path.write_text(json.dumps(params, indent=2) + '\n')
    provenance['started_unix'] = time.time()
    provenance['status'] = 'running'
    manifest = outdir / 'provenance.json'
    manifest.write_text(json.dumps(provenance, indent=2) + '\n')
    started = time.monotonic()
    try:
        run_managed(command, cwd=outdir)
        provenance['status'] = 'workflow_completed'
    except BaseException as error:
        provenance['status'] = 'interrupted' if isinstance(error, KeyboardInterrupt) else 'failed'
        raise
    finally:
        provenance['elapsed_seconds'] = time.monotonic() - started
        manifest.write_text(json.dumps(provenance, indent=2) + '\n')
        if provenance['status'] != 'workflow_completed':
            release_lock()
    return outdir
