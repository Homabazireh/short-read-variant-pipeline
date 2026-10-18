#!/usr/bin/env python3
"""Launch pinned offline VEP annotation and retain all transcript consequences."""
import csv
import json
import re
import sys
from pathlib import Path

from common import (argument_parser, checked_file, load_resources, local_path,
                    read_calls, run_workflow, sha256, text_open, validate_vcf)


def annotation_resources(resources, base, reference=None):
    annotation = dict(resources.get('annotation') or {})
    if annotation.get('vep_version') != 116:
        raise ValueError('annotation.vep_version must be 116, matching the pinned VEP container')
    cache = Path(local_path(annotation.get('vep_cache'), base, directory=True))
    release_dir = cache / 'homo_sapiens' / '116_GRCh38'
    if not release_dir.is_dir():
        raise ValueError(f'Expected Ensembl cache directory: {release_dir}')
    manifest = checked_file(annotation.get('cache_manifest'), annotation.get('cache_manifest_sha256'), base)
    listed = set()
    with open(manifest) as handle:
        for line in handle:
            checksum, name = line.rstrip().split(maxsplit=1)
            name = name.lstrip('*')
            target = (cache / name).resolve()
            if not target.is_relative_to(cache.resolve()):
                raise ValueError('Cache manifest contains a path outside the cache directory')
            checked_file(str(target), checksum)
            listed.add(target)
    actual = {p.resolve() for p in release_dir.rglob('*') if p.is_file()}
    if not actual or not actual <= listed:
        raise ValueError('Cache checksum manifest must cover every file in homo_sapiens/116_GRCh38')
    annotation.update(vep_cache=str(cache), cache_manifest=manifest)
    for name in ('clinvar', 'gnomad'):
        source = annotation.get(name)
        if source is None:
            annotation[name] = None
            continue
        source = dict(source)
        if source.get('genome_build') != 'GRCh38' or not source.get('version'):
            raise ValueError(f'{name} requires genome_build=GRCh38 and a specific version')
        source['path'] = checked_file(source.get('path'), source.get('sha256'), base)
        source['index'] = local_path(source['path'] + '.tbi')
        source['index_sha256'] = sha256(source['index'])
        fields = source.get('fields')
        if not isinstance(fields, list) or not fields or not all(re.fullmatch(r'[A-Za-z0-9_.]+', f) for f in fields):
            raise ValueError(f'{name}.fields must be a nonempty list of VCF INFO field names')
        found_fields = set()
        shared_contigs = set()
        with text_open(source['path']) as handle:
            for line in handle:
                if line.startswith('##INFO=<ID='):
                    found_fields.add(line.split('ID=', 1)[1].split(',', 1)[0])
                if line.startswith('##contig=') and reference is not None:
                    chrom = re.search(r'(?:<|,)ID=([^,>]+)', line)
                    length = re.search(r'(?:<|,)length=(\d+)', line, re.I)
                    if chrom and chrom.group(1) in reference.index:
                        shared_contigs.add(chrom.group(1))
                        if length and int(length.group(1)) != reference.index[chrom.group(1)][0]:
                            raise ValueError(f'{name} annotation/reference contig length mismatch: {line.strip()}')
                if line.startswith('#CHROM'):
                    break
            else:
                raise ValueError(f'{name} annotation VCF lacks #CHROM header')
        if not set(fields) <= found_fields:
            raise ValueError(f'{name} VCF lacks requested INFO fields: {sorted(set(fields) - found_fields)}')
        if reference is not None and not shared_contigs:
            raise ValueError(f'{name} annotation has no reference-matching contigs')
        source['shared_reference_contigs'] = sorted(shared_contigs)
        annotation[name] = source
    return annotation


def annotation_to_tsv(vcf, output, meta, annotation):
    """Explode CSQ records without discarding alternative transcripts or alleles."""
    csq_fields, records, consequences = [], 0, 0
    fixed = ['sample', 'label', 'variant_type', 'CHROM', 'POS', 'ID', 'REF', 'ALT', 'QUAL',
             'FILTER', 'FORMAT', 'GENOTYPE', 'annotation_status', 'clinvar_status',
             'gnomad_status', 'interpretation_scope', 'INFO_JSON']
    with text_open(vcf) as source, open(output, 'w', newline='') as target:
        writer = None
        for line in source:
            if line.startswith('##INFO=<ID=CSQ,'):
                match = re.search(r'Format: ([^"\n]+)', line)
                if not match:
                    raise ValueError('Cannot parse VEP CSQ format declaration')
                csq_fields = match.group(1).rstrip('>').split('|')
            if line.startswith('#CHROM'):
                writer = csv.writer(target, delimiter='\t')
                writer.writerow(fixed + ['CSQ_' + f for f in csq_fields])
            if line.startswith('#'):
                continue
            if not writer:
                raise ValueError('Annotated VCF lacks header')
            fields = line.rstrip().split('\t')
            info = dict(item.split('=', 1) if '=' in item else (item, True)
                        for item in fields[7].split(';') if item != '.')
            rows = info.pop('CSQ', '').split(',')
            records += 1
            for row in rows:
                values = row.split('|') if row else [''] * len(csq_fields)
                if len(values) != len(csq_fields):
                    raise ValueError('VEP CSQ field count differs from its header')
                statuses = []
                for name in ('clinvar', 'gnomad'):
                    prefix = name.upper()
                    if meta['variant_type'] == 'sv':
                        statuses.append('not_applicable_small_variant_resource')
                    elif annotation.get(name) is None:
                        statuses.append('unavailable')
                    elif any(k == prefix or k.startswith(prefix + '_') for k in info):
                        statuses.append('matched')
                    else:
                        # VEP can expose allele-specific custom fields within CSQ.
                        matches = [v for k, v in zip(csq_fields, values) if k.startswith(prefix)]
                        statuses.append('matched' if any(matches) else 'no_match')
                writer.writerow([meta['sample'], meta['label'], meta['variant_type']] + fields[:7]
                                + fields[8:10] + ['annotated' if row else 'no_consequence_returned']
                                + statuses + ['gene_overlap_and_supported_consequences' if meta['variant_type'] == 'sv'
                                              else 'transcript_consequences']
                                + [json.dumps(info, sort_keys=True)] + values)
                consequences += bool(row)
    return {'variants': records, 'transcript_consequences': consequences}


def main(argv=None):
    args = argument_parser(__doc__).parse_args(argv)
    resources, reference = load_resources(args.resources)
    annotation = annotation_resources(resources, Path(args.resources).resolve().parent, reference)
    rows = read_calls(args.input)
    provenance = {'resources': resources, 'annotation_resources': annotation, 'inputs': []}
    for row in rows:
        result = validate_vcf(row['vcf'], row['sample'], reference)
        row['records'] = result['records']
        provenance['inputs'].append(result)
    params = {'jobs': rows, 'reference_fasta': resources['reference_fasta'],
              'reference_fai': resources['reference_fai'], 'annotation': annotation}
    outdir = run_workflow('annotation', args, params, provenance)
    if outdir:
        outputs = []
        for row in rows:
            folder = outdir / 'callsets' / row['key']
            vcf = folder / (row['key'] + '.annotated.vcf.gz')
            tsv = folder / (row['key'] + '.annotations.tsv')
            counts = annotation_to_tsv(vcf, tsv, row, annotation)
            outputs.append({**row, **counts, 'annotated_vcf': str(vcf), 'tsv': str(tsv),
                            'vcf_sha256': sha256(vcf), 'tsv_sha256': sha256(tsv)})
        (outdir / 'outputs.json').write_text(json.dumps({'schema_version': 1, 'outputs': outputs}, indent=2) + '\n')
        provenance['status'] = 'completed'
        (outdir / 'provenance.json').write_text(json.dumps(provenance, indent=2) + '\n')
    return 0


if __name__ == '__main__':
    try:
        sys.exit(main())
    except (ValueError, OSError) as error:
        sys.exit(f'Annotation preflight/report error: {error}')
