#!/usr/bin/env python3
"""Explicit, audited HG003 integration-fixture preparation; never run on study data.

Original files are retained. Existing derivatives must match their source identity
and recorded output checksums; untracked/partial outputs are never overwritten.
"""
from __future__ import annotations

import argparse
import csv
import hashlib
import io
import json
from pathlib import Path
import re
import subprocess
import tempfile
from fetch_public_validation import (acquisition_lock, digest, immutable_text,
                                     record_derived, validate_derived)


def contig_fingerprint(contigs):
    return hashlib.sha256(json.dumps(contigs, sort_keys=True).encode()).hexdigest()


def validate_alignment_references(fields, contigs):
    if fields[2] not in contigs and fields[2] != '*':
        raise ValueError('Read mapped to an excluded contig')
    if fields[6] not in contigs and fields[6] not in ('*', '='):
        raise ValueError('Read mate mapped to an excluded contig')
    # SA and OA reference names are not encoded in BAM core reference IDs, but
    # they still describe real alignments and must remain meaningful afterward.
    for tag in fields[11:]:
        if tag.startswith(('SA:Z:', 'OA:Z:')):
            for alignment in tag[5:].rstrip('\n;').split(';'):
                if alignment and alignment.split(',')[0] not in contigs:
                    raise ValueError('Supplementary/original alignment tag refers to an excluded contig')
        if tag.startswith('CC:Z:') and tag[5:].strip() not in contigs and tag[5:].strip() != '=':
            raise ValueError('Next-alignment tag refers to an excluded contig')


def prepare_hifi(base, contigs):
    source_path = base / 'HG003.hifi.subset.bam'
    output_path = base / 'HG003.hifi.compatible.bam'
    index = Path(str(output_path) + '.bai')
    marker = base / 'hifi-preparation.json'
    if not source_path.exists():
        return
    identity = {'original_sha256': digest(source_path), 'reference_contigs_sha256': contig_fingerprint(contigs)}
    if not validate_derived(marker, identity, [output_path, index]):
        with tempfile.TemporaryDirectory(prefix='.prepare-hifi-', dir=base) as directory:
            temporary = Path(directory) / output_path.name
            source = subprocess.Popen(['samtools', 'view', '-h', str(source_path)], stdout=subprocess.PIPE, text=True)
            target = subprocess.Popen(['samtools', 'view', '-b', '-o', str(temporary), '-'], stdin=subprocess.PIPE, text=True)
            removed, samples, count = [], set(), 0
            try:
                for line in source.stdout:
                    if line.startswith('@SQ\t'):
                        fields = dict(item.split(':', 1) for item in line.rstrip().split('\t')[1:])
                        if fields['SN'] not in contigs:
                            removed.append(fields['SN'])
                            continue
                        if int(fields['LN']) != contigs[fields['SN']]:
                            raise ValueError('Retained HiFi/reference contig length mismatch')
                    elif line.startswith('@RG\t'):
                        parts = line.rstrip().split('\t')
                        for i, item in enumerate(parts):
                            if item.startswith('SM:'):
                                if item[3:] not in ('HG003', 'HG3-bc20'):
                                    raise ValueError('Unexpected sample; this helper is only for the published HG003 fixture')
                                samples.add(item[3:])
                                parts[i] = 'SM:HG003'
                        line = '\t'.join(parts) + '\n'
                    elif not line.startswith('@'):
                        fields = line.rstrip().split('\t')
                        validate_alignment_references(fields, contigs)
                        count += 1
                    target.stdin.write(line)
                target.stdin.close()
                source_code, target_code = source.wait(), target.wait()
                if source_code or target_code:
                    raise RuntimeError('SAM re-encoding failed')
            except BaseException:
                source.kill()
                target.kill()
                source.wait()
                target.wait()
                raise
            if not samples or not count:
                raise ValueError('Published HiFi fixture lacks expected sample/read records')
            subprocess.run(['samtools', 'index', str(temporary)], check=True)
            subprocess.run(['samtools', 'quickcheck', '-v', str(temporary)], check=True)
            temporary.replace(output_path)
            Path(str(temporary) + '.bai').replace(index)
        record_derived(marker, identity, [output_path, index],
                       purpose='audited bounded integration fixture only', original=str(source_path),
                       derived=str(output_path), removed_unused_reference_contigs=removed,
                       original_samples=sorted(samples), documented_sample='HG003', records_preserved=count,
                       method='SAM re-encoding after core and alignment-tag contig checks; documented sample alias standardized')
    text = io.StringIO(newline='')
    writer = csv.writer(text)
    writer.writerow(['sample', 'platform', 'reads', 'chemistry', 'basecaller', 'model'])
    writer.writerow(['HG003', 'hifi', output_path, '', '', ''])
    immutable_text(base / 'hifi.samples.csv', text.getvalue())


def prepare_wes(base, contigs):
    original = base / 'idt_capture_novogene.grch38.bed'
    if not original.exists():
        return
    target = base / 'idt_capture_novogene.grch38.autosomal.bed'
    samples = base / 'wes.samples.csv'
    output_samples = base / 'wes.autosomal.samples.csv'
    marker = target.with_suffix('.provenance.json')
    identity = {'original_sha256': digest(original), 'samples_sha256': digest(samples),
                'reference_contigs_sha256': contig_fingerprint(contigs)}
    if validate_derived(marker, identity, [target, output_samples]):
        return
    excluded, retained = 0, 0
    output = io.StringIO()
    with original.open() as source:
        for line in source:
            if not line.strip() or line.startswith(('#', 'track', 'browser')):
                continue
            fields = line.split()
            chrom, start, end = fields[0], int(fields[1]), int(fields[2])
            if chrom not in contigs or not re.fullmatch(r'chr(?:[1-9]|1[0-9]|2[0-2])', chrom):
                excluded += 1
                continue
            if not 0 <= start < end <= contigs[chrom]:
                raise ValueError('Capture interval exceeds reference contig bounds')
            output.write(f'{chrom}\t{start}\t{end}\n')
            retained += 1
    if not retained:
        raise ValueError('No autosomal capture targets retained')
    with samples.open() as handle:
        rows = list(csv.DictReader(handle))
    if not rows:
        raise ValueError('WES fixture sample sheet is empty')
    for row in rows:
        row['target_bed'] = str(target)
    sample_text = io.StringIO(newline='')
    writer = csv.DictWriter(sample_text, fieldnames=list(rows[0]))
    writer.writeheader()
    writer.writerows(rows)
    immutable_text(target, output.getvalue())
    immutable_text(output_samples, sample_text.getvalue())
    record_derived(marker, identity, [target, output_samples],
                   purpose='reference-compatible autosomal public WES integration fixture',
                   original=str(original), excluded_intervals=excluded, retained_intervals=retained)


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--data', required=True, type=Path)
    args = parser.parse_args()
    if args.data.is_symlink():
        raise ValueError('Fixture directory cannot be a symbolic link')
    base = args.data.resolve()
    acquisition = json.loads((base / 'acquisition.json').read_text())
    if acquisition.get('sample') != 'HG003':
        raise ValueError('Only the documented public HG003 fixture is supported')
    with acquisition_lock(base, 'acquisition.json', acquisition):
        reference = json.loads((base / 'resources.json').read_text())
        fasta = Path(reference['reference_fasta'])
        if digest(fasta) != reference['reference_sha256']:
            raise ValueError('Reference checksum changed')
        contigs = {line.split('\t')[0]: int(line.split('\t')[1])
                   for line in Path(str(fasta) + '.fai').read_text().splitlines()}
        prepare_hifi(base, contigs)
        prepare_wes(base, contigs)


if __name__ == '__main__':
    main()
