#!/usr/bin/env python3
"""Create a tiny deterministic plumbing fixture, never a biological benchmark."""
import argparse
import csv
import hashlib
import json
import random
import sys
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parents[2] / 'scripts'))
from benchmark import TRUVARI_PARAMETERS


def generate(folder):
    if folder.exists() and any(folder.iterdir()):
        raise ValueError('Fixture destination must be empty')
    folder.mkdir(parents=True, exist_ok=True)
    rng = random.Random(115)
    sequence = ''.join(rng.choice('ACGT') for _ in range(10000))
    fasta = folder / 'SYNTHETIC_NOT_GRCH38.fa'
    fasta.write_text('>chr20\n' + sequence + '\n')
    Path(str(fasta) + '.fai').write_text('chr20\t10000\t7\t10000\t10001\n')
    def checksum(path):
        return hashlib.sha256(path.read_bytes()).hexdigest()
    (folder / 'resources.json').write_text(json.dumps({
        'schema_version': 1, 'genome_build': 'GRCh38', 'reference_fasta': str(fasta),
        'reference_sha256': checksum(fasta),
        'warning': 'Synthetic random sequence. GRCh38 label tests validation only; NOT real GRCh38.'}, indent=2) + '\n')
    header = ('##fileformat=VCFv4.2\n##reference=GRCh38\n'
              '##source=SYNTHETIC_UNIT_FIXTURE_NOT_BIOLOGICAL_GRCH38\n'
              '##contig=<ID=chr20,length=10000>\n'
              '##FORMAT=<ID=GT,Number=1,Type=String,Description="Genotype">\n'
              '#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\tSYNTHETIC\n')
    def snp(pos, gt='0/1'):
        ref = sequence[pos - 1]
        return pos, ref, {'A': 'C', 'C': 'G', 'G': 'T', 'T': 'A'}[ref], gt
    def insertion(pos, extra):
        return pos, sequence[pos - 1], sequence[pos - 1] + extra, '0/1'
    def deletion(pos, length, gt='0/1'):
        return pos, sequence[pos - 1:pos + length], sequence[pos - 1], gt
    def vcf(name, records):
        path = folder / name
        path.write_text(header + ''.join(f'chr20\t{pos}\t.\t{ref}\t{alt}\t50\tPASS\t.\tGT\t{gt}\n'
                                         for pos, ref, alt, gt in records))
        return path
    small_truth = vcf('small.truth.vcf', [snp(100), snp(300), insertion(500, 'ACG'), snp(800)])
    small_query = vcf('small.query.vcf', [snp(100), insertion(500, 'ACG'), snp(800, '1/1'), snp(1000)])
    sv_truth = vcf('sv.truth.vcf', [insertion(2000, 'A' * 99), deletion(4000, 100), insertion(6000, 'CT' * 75)])
    sv_query = vcf('sv.query.vcf', [insertion(2000, 'A' * 100), deletion(4000, 100, '1/1'), insertion(8000, 'ATG' * 20)])
    empty = vcf('empty.vcf', [])
    bed = folder / 'all.bed'
    bed.write_text('chr20\t0\t10000\n')
    policy = {'schema_version': 1, 'truvari_version': '5.3.0',
              'truth_release': 'synthetic_plumbing_fixture_v1_NOT_GIAB',
              'release_instructions_url': 'https://github.com/ACEnglish/truvari/wiki/bench',
              'rationale': 'Synthetic INS/DEL fixture checks matches, genotype errors, empty callsets and cross-size-bin stratification; not biological accuracy.',
              'reference_sha256': checksum(fasta), 'truth_vcf_sha256': checksum(sv_truth),
              'truth_bed_sha256': checksum(bed), 'matching': TRUVARI_PARAMETERS}
    policy_path = folder / 'sv-policy.json'
    policy_path.write_text(json.dumps(policy, indent=2) + '\n')
    rows = []
    for kind, truth, query in [('small', small_truth, small_query), ('sv', sv_truth, sv_query)]:
        for label, callset in [('synthetic_perturbed', query), ('synthetic_empty', empty)]:
            rows.append({'sample': 'SYNTHETIC', 'label': label, 'vcf': str(callset), 'variant_type': kind,
                         'truth_vcf': str(truth), 'truth_bed': str(bed),
                         'sv_policy': str(policy_path) if kind == 'sv' else ''})
    with (folder / 'benchmark.csv').open('w', newline='') as handle:
        writer = csv.DictWriter(handle, fieldnames=list(rows[0]))
        writer.writeheader()
        writer.writerows(rows)
    (folder / 'README.txt').write_text('Synthetic plumbing fixture. The GRCh38 label only exercises validation; the seeded random sequence is not real GRCh38 or GIAB. Never interpret these scores as biological accuracy.\nExpected perturbed results: SNP TP=1 FN=2 FP=2; INDEL TP=1 FN=0 FP=0; SV TP=2 FN=1 FP=1, one genotype error. Empty queries have recall=0 and undefined precision/F1.\n')


if __name__ == '__main__':
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--outdir', required=True)
    args = parser.parse_args()
    try:
        generate(Path(args.outdir).resolve())
    except (ValueError, OSError) as error:
        sys.exit(str(error))
