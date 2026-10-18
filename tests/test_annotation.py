"""Synthetic VEP output/parser and pinned-resource validation tests."""
import csv
import json
import sys
import tempfile
import unittest
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parents[1] / 'scripts'))
from annotate import annotation_resources, annotation_to_tsv
from common import checked_file, sha256


class AnnotationTests(unittest.TestCase):
    def setUp(self):
        self.temp = tempfile.TemporaryDirectory()
        self.addCleanup(self.temp.cleanup)
        self.root = Path(self.temp.name)
        self.vcf = self.root / 'synthetic.vep.vcf'
        self.tsv = self.root / 'annotations.tsv'
        self.meta = {'sample': 'SYNTHETIC', 'label': 'unit', 'variant_type': 'small'}
        self.header = ('##fileformat=VCFv4.2\n'
                       '##INFO=<ID=CSQ,Number=.,Type=String,Description="Consequence annotations. Format: Allele|Consequence|Feature|PICK|GNOMAD_AF">\n'
                       '#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\tSYNTHETIC\n')

    def test_all_transcripts_and_source_availability_are_preserved(self):
        self.vcf.write_text(self.header + 'chr20\t2\t.\tC\tT\t50\tPASS\t'
                            'CSQ=T|missense_variant|ENST1|1|0.01,T|intron_variant|ENST2||0.01;CLINVAR_CLNSIG=Benign\tGT\t0/1\n')
        counts = annotation_to_tsv(self.vcf, self.tsv, self.meta, {'clinvar': {}, 'gnomad': {}})
        with self.tsv.open() as handle:
            rows = list(csv.DictReader(handle, delimiter='\t'))
        self.assertEqual(counts, {'variants': 1, 'transcript_consequences': 2})
        self.assertEqual([r['CSQ_Feature'] for r in rows], ['ENST1', 'ENST2'])
        self.assertEqual(rows[1]['CSQ_PICK'], '')
        self.assertEqual(rows[0]['clinvar_status'], 'matched')
        self.assertEqual(rows[0]['gnomad_status'], 'matched')
        self.assertEqual(json.loads(rows[0]['INFO_JSON'])['CLINVAR_CLNSIG'], 'Benign')

    def test_no_match_unavailable_and_empty_are_distinct(self):
        self.vcf.write_text(self.header + 'chr20\t2\t.\tC\tT\t50\tPASS\tCSQ=T|intron_variant|ENST1||\tGT\t0/1\n')
        annotation_to_tsv(self.vcf, self.tsv, self.meta, {'clinvar': None, 'gnomad': {}})
        with self.tsv.open() as handle:
            row = next(csv.DictReader(handle, delimiter='\t'))
        self.assertEqual(row['clinvar_status'], 'unavailable')
        self.assertEqual(row['gnomad_status'], 'no_match')
        self.vcf.write_text(self.header)
        result = annotation_to_tsv(self.vcf, self.tsv, self.meta, {})
        self.assertEqual(result['variants'], 0)
        self.assertEqual(len(self.tsv.read_text().splitlines()), 1)

    def test_malformed_consequence_and_checksum_fail(self):
        self.vcf.write_text(self.header + 'chr20\t2\t.\tC\tT\t50\tPASS\tCSQ=T|truncated\tGT\t0/1\n')
        with self.assertRaisesRegex(ValueError, 'field count'):
            annotation_to_tsv(self.vcf, self.tsv, self.meta, {})
        with self.assertRaisesRegex(ValueError, 'SHA256 mismatch'):
            checked_file(str(self.vcf), '0' * 64)

    def test_cache_version_and_manifest_coverage_are_enforced(self):
        cache = self.root / 'cache'
        release = cache / 'homo_sapiens/116_GRCh38'
        release.mkdir(parents=True)
        data = release / 'info.txt'
        data.write_text('SYNTHETIC TEST FIXTURE ONLY\n')
        manifest = self.root / 'cache.sha256'
        manifest.write_text(f'{sha256(data)}  homo_sapiens/116_GRCh38/info.txt\n')
        resources = {'annotation': {'vep_version': 116, 'vep_cache': str(cache),
                                   'cache_manifest': str(manifest), 'cache_manifest_sha256': sha256(manifest)}}
        self.assertEqual(annotation_resources(resources, self.root)['vep_version'], 116)
        (release / 'untracked.txt').write_text('changed')
        with self.assertRaisesRegex(ValueError, 'cover every file'):
            annotation_resources(resources, self.root)
        resources['annotation']['vep_version'] = 115
        with self.assertRaisesRegex(ValueError, 'must be 116'):
            annotation_resources(resources, self.root)


if __name__ == '__main__':
    unittest.main()
