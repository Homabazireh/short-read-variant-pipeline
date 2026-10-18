"""Synthetic safety checks for public-data helpers; no scientific benchmark claims."""
import gzip
import json
from pathlib import Path
import sys
import tempfile
import unittest

sys.path.insert(0, str(Path(__file__).resolve().parents[1] / 'scripts'))
from fetch_public_validation import (acquisition_lock, digest, ensure_reference_fasta,
                                     fetch, immutable_text, record_derived, validate_derived)
from prepare_public_fixture import prepare_wes, validate_alignment_references


class PublicFixtureSafetyTests(unittest.TestCase):
    def setUp(self):
        self.temp = tempfile.TemporaryDirectory()
        self.addCleanup(self.temp.cleanup)
        self.root = Path(self.temp.name)

    def test_unrelated_destination_is_not_modified(self):
        out = self.root / 'existing'
        out.mkdir()
        important = out / 'important.txt'
        important.write_text('user data')
        with self.assertRaisesRegex(ValueError, 'Nonempty directory'):
            with acquisition_lock(out, 'acquisition.json', {'sample': 'HG003'}):
                self.fail('Unrelated directory accepted')
        self.assertEqual(list(out.iterdir()), [important])
        self.assertEqual(important.read_text(), 'user data')

    def test_concurrent_acquisition_and_identity_change_are_rejected(self):
        out = self.root / 'new'
        spec = {'sample': 'HG003', 'region': 'chr20:1-10'}
        with acquisition_lock(out, 'acquisition.json', spec):
            with self.assertRaisesRegex(ValueError, 'Another acquisition'):
                with acquisition_lock(out, 'acquisition.json', spec):
                    self.fail('Concurrent writer accepted')
        with acquisition_lock(out, 'acquisition.json', spec):
            pass
        with self.assertRaisesRegex(ValueError, 'Existing metadata differs'):
            with acquisition_lock(out, 'acquisition.json', {**spec, 'region': 'chr20:11-20'}):
                self.fail('Changed acquisition accepted')

    def test_reference_is_derived_from_compressed_bytes_and_resume_is_immutable(self):
        compressed = self.root / 'reference.fa.gz'
        with gzip.open(compressed, 'wb') as handle:
            handle.write(b'>chr1\nACGT\n')
        fasta = self.root / 'reference.fa'
        expected = ensure_reference_fasta(compressed, fasta)
        before = fasta.stat().st_mtime_ns
        self.assertEqual(ensure_reference_fasta(compressed, fasta), expected)
        self.assertEqual(fasta.stat().st_mtime_ns, before)
        fasta.write_text('>chr1\nTCGT\n')
        with self.assertRaisesRegex(ValueError, 'differs from the verified compressed'):
            ensure_reference_fasta(compressed, fasta)
        self.assertEqual(fasta.read_text(), '>chr1\nTCGT\n')

    def test_untracked_corrupt_reference_is_never_authenticated(self):
        compressed = self.root / 'reference.fa.gz'
        with gzip.open(compressed, 'wb') as handle:
            handle.write(b'>chr1\nACGT\n')
        fasta = self.root / 'reference.fa'
        fasta.write_text('>chr1\nAAAA\n')
        with self.assertRaises(ValueError):
            ensure_reference_fasta(compressed, fasta)
        self.assertFalse((self.root / 'reference.fa.derived.json').exists())

    def test_derived_outputs_require_matching_sources_and_hashes(self):
        output = self.root / 'derived.bam'
        marker = self.root / 'derived.json'
        source = {'source_sha256': '0' * 64}
        self.assertFalse(validate_derived(marker, source, [output]))
        output.write_bytes(b'synthetic-unit-test-content')
        with self.assertRaisesRegex(ValueError, 'Untracked/partial'):
            validate_derived(marker, source, [output])
        record_derived(marker, source, [output])
        self.assertTrue(validate_derived(marker, source, [output]))
        with self.assertRaisesRegex(ValueError, 'identity changed'):
            validate_derived(marker, {'source_sha256': '1' * 64}, [output])
        output.write_bytes(b'corrupted')
        with self.assertRaisesRegex(ValueError, 'checksum changed'):
            validate_derived(marker, source, [output])

    def test_download_resume_rejects_changed_source_and_preserves_original_metadata(self):
        source = self.root / 'source.bin'
        source.write_bytes(b'original-source')
        target = self.root / 'download.bin'
        initial = fetch(source.as_uri(), target)
        self.assertEqual(fetch(source.as_uri(), target), initial)
        target.unlink()
        source.write_bytes(b'changed-source')
        with self.assertRaisesRegex(ValueError, 'differs from its original provenance'):
            fetch(source.as_uri(), target)
        self.assertFalse(target.exists())
        self.assertEqual(json.loads((self.root / 'download.bin.source.json').read_text()), initial)

    def test_csv_crlf_resume_does_not_rewrite(self):
        path = self.root / 'samples.csv'
        text = 'sample,platform\r\nHG003,hifi\r\n'
        immutable_text(path, text)
        before = path.stat().st_mtime_ns
        immutable_text(path, text)
        self.assertEqual(path.stat().st_mtime_ns, before)
        self.assertEqual(path.read_bytes(), text.encode())

    def test_alignment_tags_cannot_point_to_removed_contigs(self):
        fields = ['read', '0', 'chr1', '1', '60', '4M', '*', '0', '0', 'ACGT', 'IIII']
        validate_alignment_references(fields + ['SA:Z:chr1,5,+,4M,60,0;'], {'chr1': 10})
        for tag in ('SA:Z:decoy,5,+,4M,60,0;', 'OA:Z:decoy,5,+,4M,60,0;', 'CC:Z:decoy'):
            with self.subTest(tag=tag), self.assertRaisesRegex(ValueError, 'excluded contig'):
                validate_alignment_references(fields + [tag], {'chr1': 10})

    def test_wes_fixture_resume_checks_derived_bed_and_csv_without_overwrite(self):
        (self.root / 'idt_capture_novogene.grch38.bed').write_text('chr1\t0\t5\nchrX\t0\t5\n')
        (self.root / 'wes.samples.csv').write_text('sample,assay,target_bed\nHG003,WES,original.bed\n')
        prepare_wes(self.root, {'chr1': 10, 'chrX': 10})
        output = self.root / 'idt_capture_novogene.grch38.autosomal.bed'
        before = output.stat().st_mtime_ns
        prepare_wes(self.root, {'chr1': 10, 'chrX': 10})
        self.assertEqual(output.stat().st_mtime_ns, before)
        self.assertEqual(output.read_text(), 'chr1\t0\t5\n')
        output.write_text('chr1\t5\t10\n')
        with self.assertRaisesRegex(ValueError, 'checksum changed'):
            prepare_wes(self.root, {'chr1': 10, 'chrX': 10})
        self.assertEqual(output.read_text(), 'chr1\t5\t10\n')


if __name__ == '__main__':
    unittest.main()
