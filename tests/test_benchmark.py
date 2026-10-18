"""Synthetic fixture tests verify plumbing; these are NOT GIAB accuracy results."""
import csv
import json
import os
import subprocess
import sys
import tempfile
import time
import unittest
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parents[1] / 'scripts'))
from benchmark import TRUVARI_PARAMETERS, load_sv_policy, metrics, parse_happy, parse_truvari, prepare_jobs, truvari_size_strata
from common import Reference, acquire_run_lock, intersect_beds, read_bed, read_calls, run_managed, validate_output_destination, validate_vcf


class SyntheticGRCh38LabelFixture(unittest.TestCase):
    """chr20 is deliberately 200 synthetic bases, not a biological GRCh38 slice."""
    def setUp(self):
        self.temp = tempfile.TemporaryDirectory()
        self.addCleanup(self.temp.cleanup)
        self.root = Path(self.temp.name)
        self.fasta = self.root / 'synthetic.fa'
        self.fasta.write_text('>chr20\n' + 'ACGT' * 50 + '\n')
        Path(str(self.fasta) + '.fai').write_text('chr20\t200\t7\t200\t201\n')
        self.reference = Reference(str(self.fasta))
        self.addCleanup(self.reference.handle.close)
        self.header = ('##fileformat=VCFv4.2\n##reference=GRCh38\n'
                       '##contig=<ID=chr20,length=200>\n'
                       '##FORMAT=<ID=GT,Number=1,Type=String,Description="Genotype">\n'
                       '#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\tSYNTHETIC\n')
        self.vcf = self.root / 'synthetic.vcf'
        self.vcf.write_text(self.header + 'chr20\t2\t.\tC\tT\t50\tPASS\t.\tGT\t0/1\n')

    def test_every_ref_is_checked_and_empty_vcf_is_valid(self):
        self.assertEqual(validate_vcf(self.vcf, 'SYNTHETIC', self.reference)['records'], 1)
        self.vcf.write_text(self.header)
        self.assertEqual(validate_vcf(self.vcf, 'SYNTHETIC', self.reference)['records'], 0)
        self.vcf.write_text(self.header + 'chr20\t2\t.\tA\tT\t50\tPASS\t.\tGT\t0/1\n')
        with self.assertRaisesRegex(ValueError, 'REF mismatch'):
            validate_vcf(self.vcf, 'SYNTHETIC', self.reference)

    def test_wrong_sample_build_and_contig_length_rejected(self):
        with self.assertRaisesRegex(ValueError, 'expected single sample'):
            validate_vcf(self.vcf, 'OTHER', self.reference)
        self.vcf.write_text(self.header.replace('GRCh38', 'GRCh37'))
        with self.assertRaisesRegex(ValueError, 'Incompatible declared reference'):
            validate_vcf(self.vcf, 'SYNTHETIC', self.reference)
        self.vcf.write_text(self.header.replace('length=200', 'length=201'))
        with self.assertRaisesRegex(ValueError, 'length mismatch'):
            validate_vcf(self.vcf, 'SYNTHETIC', self.reference)

    def test_regions_intersect_without_trimming_truth_vcf(self):
        truth_bed = self.root / 'truth.bed'
        truth_bed.write_text('chr20\t0\t100\nchr20\t80\t180\n')
        eval_bed = self.root / 'evaluation.bed'
        eval_bed.write_text('chr20\t40\t60\n')
        input_csv = self.root / 'input.csv'
        input_csv.write_text('sample,label,vcf,variant_type,truth_vcf,truth_bed,evaluation_regions\n'
                             'SYNTHETIC,unit,synthetic.vcf,small,synthetic.vcf,truth.bed,evaluation.bed\n')
        rows = read_calls(input_csv, benchmark=True)
        original = self.vcf.read_bytes()
        jobs, evidence = prepare_jobs(rows, self.reference, self.root / 'out')
        self.assertEqual(jobs[0]['evaluated_bases'], 20)
        self.assertEqual(evidence[0]['subsets'][0]['evaluation_bed_contents'], 'chr20\t40\t60\n')
        self.assertFalse(Path(jobs[0]['evaluation_bed']).exists())
        self.assertEqual(self.vcf.read_bytes(), original)
        self.assertEqual(evidence[0]['truth']['records'], 1)

    def test_bad_regions_and_duplicate_callsets_fail(self):
        bed = self.root / 'bad.bed'
        bed.write_text('20\t0\t2\n')
        with self.assertRaisesRegex(ValueError, 'incompatibility'):
            read_bed(bed, self.reference)
        csv_path = self.root / 'input.csv'
        csv_path.write_text('sample,label,vcf,variant_type\n' + 'SYNTHETIC,x,synthetic.vcf,small\n' * 2)
        with self.assertRaisesRegex(ValueError, 'Duplicate callset'):
            read_calls(csv_path)

    def test_output_lock_unrelated_content_and_invalid_resume_rejected(self):
        out = self.root / 'new-output'
        with self.assertRaisesRegex(ValueError, 'existing run'):
            validate_output_destination(out, resume=True)
        out.mkdir()
        validate_output_destination(out, resume=False)
        release = acquire_run_lock(out)
        try:
            with self.assertRaisesRegex(ValueError, 'locked'):
                validate_output_destination(out, resume=False)
            with self.assertRaisesRegex(ValueError, 'Another process'):
                acquire_run_lock(out)
        finally:
            release()
        (out / 'unrelated.txt').write_text('do not overwrite')
        with self.assertRaisesRegex(ValueError, 'not empty'):
            validate_output_destination(out, resume=False)

    def test_sv_policy_requires_exact_reference_and_truth_identity(self):
        record = {'truth': {'sha256': 'a' * 64}, 'truth_bed': {'sha256': 'b' * 64}}
        path = self.root / 'policy.json'
        policy = {'schema_version': 1, 'truvari_version': '5.3.0', 'truth_release': 'fixture-v1',
                  'release_instructions_url': 'https://example.org/fixture', 'rationale': 'Synthetic unit test',
                  'reference_sha256': 'c' * 64, 'truth_vcf_sha256': 'a' * 64,
                  'truth_bed_sha256': 'b' * 64, 'matching': TRUVARI_PARAMETERS}
        path.write_text(json.dumps(policy))
        self.assertEqual(load_sv_policy(path, {}, record, 'c' * 64)[0]['sizemin'], 50)
        with self.assertRaisesRegex(ValueError, 'reference_sha256'):
            load_sv_policy(path, {}, record, 'd' * 64)
        with self.assertRaisesRegex(ValueError, 'truth_release differs'):
            load_sv_policy(path, {'truth_release': 'other-v1'}, record, 'c' * 64)


class BenchmarkParserTests(unittest.TestCase):
    def setUp(self):
        self.temp = tempfile.TemporaryDirectory()
        self.addCleanup(self.temp.cleanup)
        self.root = Path(self.temp.name)

    def test_happy_keeps_distinct_tp_denominators_and_unassessed_calls(self):
        summary = self.root / 'happy.summary.csv'
        summary.write_text('Type,Filter,TRUTH.TP,QUERY.TP,TRUTH.FN,QUERY.FP,QUERY.UNK,FP.gt\n'
                           'SNP,PASS,8,10,2,5,30,1\nINDEL,ALL,0,0,0,0,0,NA\n')
        rows = parse_happy(summary)
        self.assertAlmostEqual(rows[0]['precision'], 10 / 15)
        self.assertAlmostEqual(rows[0]['recall'], 8 / 10)
        self.assertEqual(rows[0]['query_unassessed'], 30)
        self.assertIsNone(rows[1]['precision'])
        self.assertIsNone(rows[1]['recall'])
        self.assertIsNone(rows[1]['f1'])

    def test_truvari_allele_and_genotype_counts_are_separate(self):
        report = self.root / 'summary.json'
        report.write_text(json.dumps({'TP-base': 8, 'TP-comp': 10, 'FN': 2, 'FP': 5,
                                      'TP-comp_FP-gt': 3, 'precision': 0.66667}))
        row = parse_truvari(report)[0]
        self.assertEqual(row['truth_tp'], 8)
        self.assertEqual(row['query_tp'], 10)
        self.assertEqual(row['genotype_errors'], 3)
        self.assertAlmostEqual(row['precision'], 2 / 3)
        self.assertIsNone(row['query_unassessed'])

    def test_happy_0312_summary_uses_explicit_extended_query_tp(self):
        summary = self.root / 'happy.summary.csv'
        summary.write_text('Type,Filter,TRUTH.TP,TRUTH.FN,QUERY.FP,QUERY.UNK,FP.gt\nSNP,PASS,8,2,5,30,1\n')
        (self.root / 'happy.extended.csv').write_text('Type,Subtype,Subset,Filter,Genotype,QQ,TRUTH.TP,QUERY.TP,TRUTH.FN,QUERY.FP,QUERY.UNK,Subset.IS_CONF.Size\nSNP,*,*,PASS,*,*,8,10,2,5,30,10005.000000\n')
        row = parse_happy(summary)[0]
        self.assertEqual(row['query_tp'], 10)
        self.assertEqual(row['truth_tp'], 8)
        self.assertEqual(row['native_adjusted_confident_bases'], 10005)

    def test_nan_missing_and_negative_counts_are_errors(self):
        report = self.root / 'summary.json'
        for bad in (None, 'NA', 'nan', -1, 1.5):
            report.write_text(json.dumps({'TP-base': bad, 'TP-comp': 1, 'FN': 0, 'FP': 0}))
            with self.assertRaises(ValueError):
                parse_truvari(report)

    def test_zero_callset_does_not_become_perfect_accuracy(self):
        result = metrics(0, 0, 10, 0)
        self.assertEqual(result['recall'], 0)
        self.assertIsNone(result['precision'])
        self.assertIsNone(result['f1'])
        self.assertEqual(metrics(0, 0, 10, 10)['f1'], 0)

    @unittest.skipUnless(hasattr(os, 'killpg'), 'POSIX process-group lifecycle')
    def test_failed_parent_cleans_term_resistant_child_without_killing_unrelated_job(self):
        pidfile = self.root / 'child.pid'
        child = ('import os,signal,time,pathlib,sys; '
                 'signal.signal(signal.SIGTERM,signal.SIG_IGN); '
                 'pathlib.Path(sys.argv[1]).write_text(str(os.getpid())); time.sleep(60)')
        parent = ('import subprocess,sys,time,pathlib\n'
                  'subprocess.Popen([sys.executable,"-c",sys.argv[1],sys.argv[2]])\n'
                  'while not pathlib.Path(sys.argv[2]).exists(): time.sleep(0.01)\n'
                  'sys.exit(7)\n')
        unrelated = subprocess.Popen([sys.executable, '-c', 'import time; time.sleep(60)'], start_new_session=True)
        child_pid = None
        try:
            with self.assertRaises(subprocess.CalledProcessError):
                run_managed([sys.executable, '-c', parent, child, str(pidfile)], self.root, cleanup_timeout=0.2)
            child_pid = int(pidfile.read_text())
            deadline = time.monotonic() + 2
            while time.monotonic() < deadline:
                stat = Path(f'/proc/{child_pid}/stat')
                if not stat.exists() or stat.read_text().split()[2] == 'Z':
                    break
                time.sleep(0.02)
            else:
                self.fail('Descendant survived process-group cleanup')
            self.assertIsNone(unrelated.poll())
        finally:
            unrelated.terminate()
            unrelated.wait(timeout=5)
            if child_pid is not None:
                try:
                    os.kill(child_pid, 9)
                except ProcessLookupError:
                    pass

    def test_sv_postmatch_size_bins_preserve_cross_bin_matches(self):
        header = '##fileformat=VCFv4.2\n#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\n'
        for name, size in [('tp-base', 99), ('tp-comp', 100), ('fn', 150), ('fp', 60)]:
            (self.root / (name + '.vcf')).write_text(header + 'chr20\t100\t.\tA\tA' + 'C' * size + '\t50\tPASS\tGTMatch=1\n')
        rows = {r['size_bin']: r for r in truvari_size_strata(self.root) if r['type'] == 'INS'}
        self.assertEqual(rows['50-99']['truth_tp'], 1)
        self.assertEqual(rows['50-99']['query_tp'], 0)
        self.assertEqual(rows['100-499']['truth_tp'], 0)
        self.assertEqual(rows['100-499']['query_tp'], 1)
        self.assertEqual(rows['100-499']['genotype_errors'], 1)
        self.assertEqual(rows['50-99']['recall'], 1)
        self.assertEqual(rows['100-499']['precision'], 1)


if __name__ == '__main__':
    unittest.main()
