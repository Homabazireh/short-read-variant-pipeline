#!/usr/bin/env python3
"""Benchmark callsets against pinned truth using hap.py/vcfeval or Truvari."""
import csv
import html
import json
import math
import re
import sys
from pathlib import Path

from common import (argument_parser, intersect_beds, load_resources, local_path, read_bed,
                    read_calls, run_workflow, sha256, text_open, validate_vcf)

TRUVARI_PARAMETERS = {'refdist': 500, 'pctseq': 0.7, 'pctsize': 0.7, 'pctovl': 0.0,
                      'sizemin': 50, 'sizefilt': 50, 'sizemax': -1, 'pick': 'single',
                      'passonly': True, 'no_ref': 'a', 'extend': 0,
                      'scope': 'sequence_resolved_biallelic_INS_DEL'}
METRIC_COLUMNS = ['sample', 'label', 'variant_type', 'subset', 'type', 'filter',
                  'size_bin', 'stratification_method',
                  'requested_input_bed_bases', 'native_adjusted_confident_bases',
                  'truth_tp', 'query_tp', 'fn', 'fp', 'query_unassessed',
                  'precision', 'recall', 'f1', 'genotype_errors', 'native_report']


def count(value, field):
    try:
        number = float(value)
    except (TypeError, ValueError):
        raise ValueError(f'Missing/non-numeric benchmark count {field}: {value!r}') from None
    if not math.isfinite(number) or number < 0 or not number.is_integer():
        raise ValueError(f'Invalid benchmark count {field}: {value!r}')
    return int(number)


def metrics(truth_tp, query_tp, fn, fp):
    precision = query_tp / (query_tp + fp) if query_tp + fp else None
    recall = truth_tp / (truth_tp + fn) if truth_tp + fn else None
    f1 = None if precision is None or recall is None else (
        2 * precision * recall / (precision + recall) if precision + recall else 0.0)
    return {'truth_tp': truth_tp, 'query_tp': query_tp, 'fn': fn, 'fp': fp,
            'precision': precision, 'recall': recall, 'f1': f1}


def parse_happy(path):
    with open(path, newline='') as handle:
        rows = list(csv.DictReader(handle))
    extended_tp = {}
    extended = Path(str(path).replace('.summary.csv', '.extended.csv'))
    if extended.is_file():
        with extended.open(newline='') as handle:
            for row in csv.DictReader(handle):
                if all(row.get(k) == '*' for k in ('Subtype', 'Subset', 'Genotype', 'QQ')):
                    key = (row['Type'], row['Filter'])
                    if key in extended_tp:
                        raise ValueError(f'Ambiguous aggregate rows in {extended}')
                    extended_tp[key] = row
    result = []
    for row in rows:
        matching = extended_tp.get((row['Type'], row['Filter']))
        if 'QUERY.TP' not in row:
            if matching is None or any(count(row[k], k) != count(matching[k], k) for k in ('TRUTH.TP', 'TRUTH.FN', 'QUERY.FP', 'QUERY.UNK')):
                raise ValueError('hap.py summary/extended aggregate rows are missing or inconsistent')
            row['QUERY.TP'] = matching['QUERY.TP']
            row['QUERY.TP_source'] = 'native_extended_csv_aggregate'
        values = {k: count(row.get(k), k) for k in ['TRUTH.TP', 'QUERY.TP', 'TRUTH.FN', 'QUERY.FP']}
        result.append({'type': row['Type'], 'filter': row['Filter'],
                       **metrics(values['TRUTH.TP'], values['QUERY.TP'], values['TRUTH.FN'], values['QUERY.FP']),
                       'query_unassessed': count(row['QUERY.UNK'], 'QUERY.UNK'),
                       'genotype_errors': count(row['FP.gt'], 'FP.gt') if row.get('FP.gt') not in (None, '', 'NA', '.') else None,
                       'native_adjusted_confident_bases': count(matching['Subset.IS_CONF.Size'], 'Subset.IS_CONF.Size')
                           if matching and matching.get('Subset.IS_CONF.Size') not in (None, '', 'NA', '.') else None,
                       'native_report': str(path), 'native_metrics': row})
    if not result:
        raise ValueError(f'Empty hap.py summary: {path}')
    return result


def parse_truvari(path):
    with open(path) as handle:
        raw = json.load(handle)
    values = {k: count(raw.get(k), k) for k in ('TP-base', 'TP-comp', 'FN', 'FP')}
    return [{'type': 'INS_DEL_ge50', 'filter': 'PASS',
             **metrics(values['TP-base'], values['TP-comp'], values['FN'], values['FP']),
             'query_unassessed': None,
             'genotype_errors': count(raw['TP-comp_FP-gt'], 'TP-comp_FP-gt') if 'TP-comp_FP-gt' in raw else None,
             'native_report': str(path), 'native_metrics': raw}]


SV_BINS = [(50, 100, '50-99'), (100, 500, '100-499'), (500, 1000, '500-999'),
           (1000, 10000, '1000-9999'), (10000, None, '10000+')]


def truvari_size_strata(folder):
    """Stratify completed matches by each representation's own SV type/length."""
    tally = {(kind, label): {'truth_tp': 0, 'query_tp': 0, 'fn': 0, 'fp': 0, 'genotype_errors': 0}
             for kind in ('INS', 'DEL') for _, _, label in SV_BINS}
    for basename, count_key in [('tp-base', 'truth_tp'), ('tp-comp', 'query_tp'), ('fn', 'fn'), ('fp', 'fp')]:
        path = Path(folder) / (basename + '.vcf.gz')
        if not path.exists():
            path = Path(folder) / (basename + '.vcf')
        with text_open(path) as handle:
            for line in handle:
                if line.startswith('#'):
                    continue
                fields = line.rstrip().split('\t')
                ref, alt = fields[3:5]
                if not re.fullmatch(r'[ACGTNacgtn]+', alt):
                    raise ValueError('Unexpected non-sequence allele in Truvari size stratification')
                delta = len(alt) - len(ref)
                if not delta or min(len(ref), len(alt)) != 1:
                    raise ValueError('Unexpected non-INS/DEL allele in Truvari size stratification')
                label = next((name for low, high, name in SV_BINS if low <= abs(delta) and (high is None or abs(delta) < high)), None)
                if label is None:
                    raise ValueError('Unexpected sub-50 bp allele in Truvari output')
                bucket = tally[('INS' if delta > 0 else 'DEL', label)]
                bucket[count_key] += 1
                if count_key == 'query_tp':
                    info = dict(item.split('=', 1) if '=' in item else (item, True) for item in fields[7].split(';'))
                    if 'GTMatch' in info:
                        bucket['genotype_errors'] += int(info['GTMatch']) != 0
    results = []
    for (kind, label), value in tally.items():
        results.append({'type': kind, 'filter': 'PASS', 'size_bin': label,
                        'stratification_method': 'postmatch_each_representation',
                        **metrics(value['truth_tp'], value['query_tp'], value['fn'], value['fp']),
                        'query_unassessed': None, 'genotype_errors': value['genotype_errors'],
                        'native_report': str(folder)})
    return results


def load_sv_policy(path, row, record, reference_sha):
    if not path:
        raise ValueError('SV benchmarking requires --sv-policy or an sv_policy CSV column; review the truth release instructions first')
    path = local_path(path)
    with open(path) as handle:
        policy = json.load(handle)
    if policy.get('schema_version') != 1 or policy.get('truvari_version') != '5.3.0':
        raise ValueError('SV policy requires schema_version=1 and truvari_version=5.3.0')
    identity = {'reference_sha256': reference_sha, 'truth_vcf_sha256': record['truth']['sha256'],
                'truth_bed_sha256': record['truth_bed']['sha256']}
    for key, expected in identity.items():
        if policy.get(key) != expected:
            raise ValueError(f'SV policy {key} does not match this reference/truth bundle')
    if not policy.get('truth_release') or not policy.get('release_instructions_url') or not policy.get('rationale'):
        raise ValueError('SV policy requires truth_release, release_instructions_url, and rationale')
    if row.get('truth_release') and row['truth_release'] != policy['truth_release']:
        raise ValueError('SV policy truth_release differs from the CSV declaration')
    matching = policy.get('matching')
    if not isinstance(matching, dict) or set(matching) != set(TRUVARI_PARAMETERS):
        raise ValueError('SV policy matching must specify every parameter in the supplied policy template')
    for key in ('pctseq', 'pctsize', 'pctovl'):
        if not isinstance(matching[key], (int, float)) or not 0 <= matching[key] <= 1:
            raise ValueError(f'Invalid SV matching threshold: {key}')
    if not isinstance(matching['refdist'], int) or matching['refdist'] < 0:
        raise ValueError('SV policy refdist must be a nonnegative integer')
    for key in ('sizemin', 'sizefilt', 'sizemax', 'pick', 'passonly', 'no_ref', 'extend', 'scope'):
        if matching[key] != TRUVARI_PARAMETERS[key]:
            raise ValueError(f'SV policy {key} unsupported by this sequence-resolved >=50 bp workflow')
    return matching, {'path': path, 'sha256': sha256(path), 'definition': policy}


def strata(path, reference):
    if not path:
        return []
    output, seen = [], {'overall'}
    with open(path) as handle:
        for line in handle:
            if not line.strip() or line.startswith('#'):
                continue
            name, bed = line.rstrip().split('\t')
            if not re.fullmatch(r'[A-Za-z0-9][A-Za-z0-9_.-]*', name) or name in seen:
                raise ValueError(f'Invalid/duplicate stratification name: {name}')
            seen.add(name)
            bed = local_path(bed, Path(path).parent)
            output.append((name, read_bed(bed, reference), {'path': bed, 'sha256': sha256(bed)}))
    return output


def prepare_jobs(rows, reference, outdir, evaluation_regions=None, dry_run=False, sv_policy=None, reference_sha=None):
    jobs, evidence = [], []
    validation_cache = {}
    for row in rows:
        truth_sample = row.get('truth_sample') or row['sample']
        if not re.fullmatch(r'[A-Za-z0-9][A-Za-z0-9_.-]*', truth_sample):
            raise ValueError('Invalid truth_sample')
        for field, sample in [('vcf', row['sample']), ('truth_vcf', truth_sample)]:
            cache_key = (row[field], sample)
            if cache_key not in validation_cache:
                validation_cache[cache_key] = validate_vcf(row[field], sample, reference)
        truth_regions = read_bed(row['truth_bed'], reference)
        allowed = truth_regions
        eval_path = evaluation_regions or row.get('evaluation_regions')
        record = {'sample': row['sample'], 'label': row['label'],
                  'query': validation_cache[(row['vcf'], row['sample'])],
                  'truth': validation_cache[(row['truth_vcf'], truth_sample)],
                  'truth_bed': {'path': row['truth_bed'], 'sha256': sha256(row['truth_bed'])},
                  'truth_release': row.get('truth_release') or 'unspecified_user_supplied_truth',
                  'truth_release_notes': row.get('truth_release_notes') or None}
        matching = None
        if row['variant_type'] == 'sv':
            matching, record['sv_policy'] = load_sv_policy(sv_policy or row.get('sv_policy'), row, record, reference_sha)
            record['truth_release'] = record['sv_policy']['definition']['truth_release']
            record['truth_release_notes'] = record['sv_policy']['definition']['release_instructions_url']
        if eval_path:
            eval_path = local_path(eval_path)
            allowed = intersect_beds(allowed, read_bed(eval_path, reference))
            record['evaluation_regions'] = {'path': eval_path, 'sha256': sha256(eval_path)}
        if not allowed:
            raise ValueError(f'No confident evaluation regions remain for {row["key"]}')
        subsets = [('overall', allowed, None)]
        for name, regions, metadata in strata(row.get('stratifications'), reference):
            region_subset = intersect_beds(allowed, regions)
            if not region_subset:
                raise ValueError(f'Stratification {name} has no overlap with evaluation regions')
            subsets.append((name, region_subset, metadata))
        record['subsets'] = []
        for name, regions, metadata in subsets:
            key = row['key'] + '__' + name
            bed = outdir / '_inputs' / (key + '.bed')
            content = ''.join(f'{chrom}\t{start}\t{end}\n' for chrom in reference.index
                              for start, end in regions.get(chrom, []))
            region_count = sum(end - start for spans in regions.values() for start, end in spans)
            jobs.append({**row, 'key': key, 'subset': name, 'truth_sample': truth_sample,
                         'evaluation_bed': str(bed), 'evaluated_bases': region_count, 'matching': matching})
            record['subsets'].append({'name': name, 'bases': region_count, 'source': metadata,
                                      'evaluation_bed': str(bed),
                                      'evaluation_bed_contents': content})
        evidence.append(record)
    return jobs, evidence


def write_reports(outdir, jobs, provenance):
    records = []
    for job in jobs:
        folder = outdir / 'callsets' / job['key']
        native = folder / ('happy.summary.csv' if job['variant_type'] == 'small' else 'truvari/summary.json')
        rows = parse_happy(native) if job['variant_type'] == 'small' else parse_truvari(native)
        if job['variant_type'] == 'sv':
            rows.extend(truvari_size_strata(folder / 'truvari'))
        for row in rows:
            row.update({k: job[k] for k in ('sample', 'label', 'variant_type', 'subset', 'evaluated_bases')})
            row['requested_input_bed_bases'] = job['evaluated_bases']
            records.append(row)
    report = {'schema_version': 1, 'metrics': records,
              'notes': ['Truth and query true positives use their own representations and denominators.',
                        'Undefined ratios are null/blank, never perfect scores.',
                        'Strata and PASS/ALL overlap; do not sum rows.',
                        'SV size strata partition completed matches by each representation; cross-bin matches remain matched.',
                        'requested_input_bed_bases (legacy evaluated_bases) is the input BED intersection; hap.py adjusts insertion boundaries, reported separately as native_adjusted_confident_bases.',
                        'Unassessed calls are excluded from precision and are not false positives.',
                        'SV allele matching and genotype-error counts are separate; Truvari TP is not genotype concordance.'],
              'provenance': provenance}
    (outdir / 'metrics.json').write_text(json.dumps(report, indent=2, allow_nan=False) + '\n')
    with open(outdir / 'metrics.csv', 'w', newline='') as handle:
        writer = csv.DictWriter(handle, fieldnames=METRIC_COLUMNS, extrasaction='ignore')
        writer.writeheader()
        writer.writerows(records)
    columns = [c for c in METRIC_COLUMNS if c != 'native_report']
    markup = ['<!doctype html><html lang="en"><meta charset="utf-8"><title>Variant benchmark report</title>',
              '<style>body{font:15px system-ui;margin:2rem}table{border-collapse:collapse}th,td{padding:.4rem;border:1px solid #bbb}th{background:#eef}code{white-space:pre-wrap}</style>',
              '<h1>Variant benchmark report</h1>',
              '<p>Measured results from pinned tools. Consult provenance.json, trace.tsv and native reports for parameters, resource use and matching details.</p>',
              '<ul>' + ''.join('<li>' + html.escape(n) + '</li>' for n in report['notes']) + '</ul>',
              '<table><thead><tr>' + ''.join('<th>' + html.escape(c) + '</th>' for c in columns) + '</tr></thead><tbody>']
    for row in records:
        markup.append('<tr>' + ''.join('<td>' + html.escape('NA' if row.get(c) is None else
                       f'{row[c]:.6f}' if isinstance(row[c], float) else str(row.get(c, ''))) + '</td>' for c in columns) + '</tr>')
    markup.append('</tbody></table></html>')
    (outdir / 'report.html').write_text('\n'.join(markup) + '\n')
    (outdir / 'outputs.json').write_text(json.dumps({'schema_version': 1, 'outputs': [
        {'path': str(outdir / name), 'sha256': sha256(outdir / name)}
        for name in ('metrics.json', 'metrics.csv', 'report.html', 'trace.tsv')]}, indent=2) + '\n')


def main(argv=None):
    parser = argument_parser(__doc__)
    parser.add_argument('--regions', help='BED restricting evaluation, intersected with truth BED; never trims truth VCF')
    parser.add_argument('--sv-policy', help='Release-specific checksummed SV benchmark policy JSON (required for SV jobs)')
    args = parser.parse_args(argv)
    resources, reference = load_resources(args.resources)
    rows = read_calls(args.input, benchmark=True)
    outdir = Path(args.outdir).resolve()
    if (outdir / 'run.params.json').exists() and not args.resume and not args.dry_run:
        raise ValueError('Output already exists; use --resume or a new directory')
    jobs, evidence = prepare_jobs(rows, reference, outdir, args.regions, args.dry_run, args.sv_policy, resources['reference_sha256'])
    provenance = {'resources': resources, 'inputs': evidence, 'truvari_defaults': TRUVARI_PARAMETERS,
                  'happy_parameters': {'engine': 'vcfeval', 'confidence_mode': 'truth_intersection_evaluation_regions',
                                       'adjust_conf_regions': True, 'truth_trimming': False}}
    params = {'jobs': jobs, 'reference_fasta': resources['reference_fasta'],
              'reference_fai': resources['reference_fai'], 'truvari_defaults': TRUVARI_PARAMETERS,
              'generated_files': {subset['evaluation_bed']: subset['evaluation_bed_contents']
                                  for row in evidence for subset in row['subsets']}}
    completed = run_workflow('benchmark', args, params, provenance)
    if completed:
        try:
            write_reports(completed, jobs, provenance)
            provenance['status'] = 'completed'
        except BaseException:
            provenance['status'] = 'report_failed'
            raise
        finally:
            (outdir / 'provenance.json').write_text(json.dumps(provenance, indent=2) + '\n')
    return 0


if __name__ == '__main__':
    try:
        sys.exit(main())
    except (ValueError, OSError) as error:
        sys.exit(f'Benchmark preflight/report error: {error}')
