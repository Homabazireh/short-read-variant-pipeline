"""Fast validation tests; these do not claim sequencing pipeline accuracy."""
import argparse
import csv
import gzip
import importlib.util
import json
import io
import os
import sys
import time
from contextlib import redirect_stdout
from pathlib import Path
import tempfile
import unittest
from unittest.mock import patch, MagicMock, call

ROOT = Path(__file__).resolve().parents[1]
spec = importlib.util.spec_from_file_location('shortread', ROOT / 'scripts' / 'run_short_reads.py')
sr = importlib.util.module_from_spec(spec)
spec.loader.exec_module(sr)


class ShortReadTests(unittest.TestCase):
    def setUp(self):
        self.temp = tempfile.TemporaryDirectory()
        self.base = Path(self.temp.name)
        self.fasta = self.base / 'reference.fa'
        self.fasta.write_text('>chr1\n' + 'A' * 1000 + '\n>chrX\n' + 'A' * 100 + '\n')
        Path(str(self.fasta) + '.fai').write_text('chr1\t1000\t6\t1000\t1001\nchrX\t100\t1013\t100\t101\n')
        self.resource = self.base / 'resources.json'
        self.resource.write_text(json.dumps({'schema_version': 1, 'genome_build': 'GRCh38',
                                             'reference_fasta': str(self.fasta), 'reference_sha256': sr.sha256(self.fasta)}))
        self.reads = []
        for mate in (1, 2):
            path = self.base / f'read{mate}.fastq.gz'
            with gzip.open(path, 'wt') as handle:
                handle.write(f'@pair/{mate}\nACGT\n+\nIIII\n')
            self.reads.append(path)
        self.row = {'sample': 'HG002', 'assay': 'WGS', 'lane': 'L001', 'fastq_1': str(self.reads[0]),
                    'fastq_2': str(self.reads[1]), 'target_bed': ''}
        self.samples = self.base / 'samples.csv'
        self.write_samples([self.row])
        self.args = sr.parser().parse_args(['--input', str(self.samples), '--resources', str(self.resource), '--outdir', str(self.base / 'results')])

    def tearDown(self):
        self.temp.cleanup()

    def write_samples(self, rows):
        sr.write_csv(self.samples, sr.FIELDS, rows)

    def resources(self, caller='deepvariant'):
        return sr.load_resources(self.resource, caller)

    def test_default_is_autosomal_deepvariant_without_bqsr(self):
        resource = self.resources()
        group = sr.make_groups(sr.load_samples(self.samples), resource, 100)[0]
        params = sr.make_params(self.args, resource, group, self.base / 'results')
        self.assertEqual(group['calling_intervals'], [('chr1', 0, 1000)])
        self.assertEqual(params['aligner'], 'bwa-mem2')
        self.assertEqual(params['tools'], 'deepvariant')
        self.assertEqual(params['skip_tools'], 'baserecalibrator')
        self.assertTrue(params['igenomes_ignore'])
        self.assertIsNone(params['custom_config_base'])
        self.assertTrue(params['save_output_as_bam'])
        self.assertFalse(params['wes'])
        self.assertNotIn('use_gatk_spark', params)  # Upstream schema expects string, not bool.
        self.assertNotIn('tiddit', params['tools'])

    def test_wes_is_separate_and_padding_is_clipped_and_merged(self):
        bed = self.base / 'capture.bed'
        bed.write_text('chr1\t20\t40\nchr1\t80\t90\nchr1\t990\t1000\nchrX\t0\t50\n')
        row = {**self.row, 'assay': 'WES', 'target_bed': str(bed)}
        self.write_samples([row])
        group = sr.make_groups(sr.load_samples(self.samples), self.resources(), 100)[0]
        self.assertEqual(group['reporting_intervals'], [('chr1', 20, 40), ('chr1', 80, 90), ('chr1', 990, 1000)])
        self.assertEqual(group['calling_intervals'], [('chr1', 0, 190), ('chr1', 890, 1000)])
        self.assertTrue(sr.make_params(self.args, self.resources(), group, self.base)['wes'])

    def test_regions_scope_intersects_autosomes(self):
        bed = self.base / 'scope.bed'
        bed.write_text('chr1\t100\t200\nchrX\t0\t100\n')
        group = sr.make_groups(sr.load_samples(self.samples), self.resources(), 100, bed)[0]
        self.assertEqual(group['calling_intervals'], [('chr1', 100, 200)])
        self.assertEqual(group['reporting_intervals'], [('chr1', 100, 200)])

    def test_wes_scope_keeps_reporting_distinct_from_padding(self):
        capture = self.base / 'capture.bed'
        capture.write_text('chr1\t100\t200\n')
        scope = self.base / 'scope.bed'
        scope.write_text('chr1\t150\t250\n')
        self.write_samples([{**self.row, 'assay': 'WES', 'target_bed': str(capture)}])
        group = sr.make_groups(sr.load_samples(self.samples), self.resources(), 100, scope)[0]
        self.assertEqual(group['reporting_intervals'], [('chr1', 150, 200)])
        self.assertEqual(group['calling_intervals'], [('chr1', 150, 250)])

    def test_capture_with_incompatible_reference_contig_fails(self):
        capture = self.base / 'capture.bed'
        capture.write_text('chr1\t100\t200\nchr17_KI270909v1_alt\t1\t20\n')
        with self.assertRaisesRegex(sr.ValidationError, 'outside reference bounds'):
            sr.target_intervals(capture, self.resources()['lengths'], self.resources()['autosomes'])

    def test_reference_checksum_and_stale_index_rejected(self):
        self.fasta.write_text('>chr1\nAAAA\n')
        with self.assertRaisesRegex(sr.ValidationError, 'SHA-256 mismatch'):
            self.resources()
        data = json.loads(self.resource.read_text())
        data['reference_sha256'] = sr.sha256(self.fasta)
        self.resource.write_text(json.dumps(data))
        with self.assertRaisesRegex(sr.ValidationError, 'index contigs/lengths'):
            self.resources()

    def test_wrong_fai_offset_is_rejected_even_if_lengths_match(self):
        fai = Path(str(self.fasta) + '.fai')
        fai.write_text(fai.read_text().replace('chr1\t1000\t6\t', 'chr1\t1000\t7\t'))
        with self.assertRaisesRegex(sr.ValidationError, 'offsets/line widths'):
            self.resources()

    def test_missing_known_sites_fails_gatk(self):
        with self.assertRaisesRegex(sr.ValidationError, 'known_sites_vcf'):
            self.resources('haplotypecaller')

    def test_known_sites_mismatched_reference_rejected(self):
        vcf = self.base / 'known.vcf.gz'
        with gzip.open(vcf, 'wt') as handle:
            handle.write('##fileformat=VCFv4.2\n##contig=<ID=chr1,length=999>\n#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\n')
        Path(str(vcf) + '.tbi').write_bytes(b'index-fixture')
        data = json.loads(self.resource.read_text())
        data['known_sites_vcf'] = [{'path': str(vcf), 'sha256': sr.sha256(vcf)}]
        self.resource.write_text(json.dumps(data))
        with self.assertRaisesRegex(sr.ValidationError, 'contig mismatch'):
            self.resources('haplotypecaller')

    def test_pair_mismatch_rejected(self):
        with gzip.open(self.reads[1], 'wt') as handle:
            handle.write('@wrong/2\nACGT\n+\nIIII\n')
        with self.assertRaisesRegex(sr.ValidationError, 'paired-read identifiers'):
            sr.load_samples(self.samples)

    def test_duplicate_lanes_rejected(self):
        self.write_samples([self.row, self.row])
        with self.assertRaisesRegex(sr.ValidationError, 'Duplicate sample/lane'):
            sr.load_samples(self.samples)

    def test_wes_requires_target_and_wgs_rejects_target(self):
        self.write_samples([{**self.row, 'assay': 'WES'}])
        with self.assertRaisesRegex(sr.ValidationError, 'target_bed'):
            sr.load_samples(self.samples)
        self.write_samples([{**self.row, 'target_bed': 'capture.bed'}])
        with self.assertRaisesRegex(sr.ValidationError, 'WGS rows'):
            sr.load_samples(self.samples)

    def test_invalid_header_and_path_rejected(self):
        self.samples.write_text('sample,fastq_1,fastq_2\n')
        with self.assertRaisesRegex(sr.ValidationError, 'header must be exactly'):
            sr.load_samples(self.samples)
        with self.assertRaisesRegex(sr.ValidationError, 'unsupported characters'):
            sr.checked_path('/tmp/data $(whoami).fastq.gz', self.base, 'input')

    def test_unsafe_outdir_and_resume_mismatch_rejected(self):
        output = self.base / 'results'
        output.mkdir()
        (output / 'run_manifest.json').write_text(json.dumps({'fingerprint': 'one'}))
        with self.assertRaisesRegex(sr.ValidationError, 'nonempty'):
            sr.validate_output(output, False, 'one')
        with self.assertRaisesRegex(sr.ValidationError, 'Resume inputs'):
            sr.validate_output(output, True, 'two')
        sr.validate_output(output, True, 'one')

    def test_command_pins_revision_and_preserves_resume(self):
        group = {'id': 'wgs'}
        self.args.resume = True
        command = sr.nextflow_command(self.args, group, self.base)
        self.assertEqual(command[command.index('-r') + 1], sr.SAREK_REVISION)
        self.assertIn('-resume', command)
        self.assertEqual(command[command.index('-profile') + 1], 'apptainer')
        self.assertNotIn('--max_cpus', command)
        self.assertIn('process.resourceLimits', sr.execution_config(self.args))

    def test_collect_uses_upstream_manifest_and_checks_all_samples(self):
        base = self.base / 'group'
        (base / 'sarek' / 'csv').mkdir(parents=True)
        vcf = self.base / 'nonstandard-output-name.vcf.gz'
        vcf.write_bytes(b'output-fixture')
        Path(str(vcf) + '.tbi').write_bytes(b'index-fixture')
        manifest = base / 'sarek' / 'csv' / 'variantcalled.csv'
        sr.write_csv(manifest, ['patient', 'sample', 'variantcaller', 'vcf'],
                     [{'patient': 'HG002', 'sample': 'HG002', 'variantcaller': 'deepvariant', 'vcf': str(vcf)}])
        rows = sr.collect_calls(base, {'HG002'}, 'deepvariant')
        self.assertEqual(rows[0]['vcf'], str(vcf))
        with self.assertRaisesRegex(sr.ValidationError, 'exactly one'):
            sr.collect_calls(base, {'HG002', 'HG003'}, 'deepvariant')

    def test_dry_run_makes_no_output_and_needs_no_runtime(self):
        args = ['--input', str(self.samples), '--resources', str(self.resource), '--outdir', str(self.base / 'dry'), '--dry-run']
        with patch('sys.stdout') as stdout, patch.object(sr.shutil, 'which', return_value=None):
            self.assertEqual(sr.main(args), 0)
        self.assertFalse((self.base / 'dry').exists())

    def test_interruption_waits_for_nextflow_before_releasing_control(self):
        process = MagicMock()
        process.pid = 123456
        process.stdout.__iter__.side_effect = KeyboardInterrupt
        with patch.object(sr.subprocess, 'Popen', return_value=process) as popen, patch.object(sr.os, 'killpg') as kill:
            with self.assertRaises(KeyboardInterrupt):
                sr.run_nextflow(['nextflow'], self.base, {}, self.base / 'run.log')
        self.assertTrue(popen.call_args.kwargs['start_new_session'])
        self.assertEqual(kill.call_args_list, [call(123456, sr.signal.SIGTERM), call(123456, sr.signal.SIGKILL)])
        self.assertEqual(process.wait.call_args_list, [call(timeout=10), call()])
        process.stdout.close.assert_called_once()

    def test_log_write_failure_reaps_real_child_process(self):
        log = MagicMock()
        log.open.return_value.__enter__.return_value.write.side_effect = OSError('simulated log write failure')
        capture = io.StringIO()
        with redirect_stdout(capture), self.assertRaisesRegex(OSError, 'simulated log write failure'):
            sr.run_nextflow([sys.executable, '-u', '-c',
                             'import os,time; print(os.getpid(),flush=True); time.sleep(10)'],
                            self.base, os.environ.copy(), log)
        pid = int(capture.getvalue().strip())
        try:
            with self.assertRaises(ProcessLookupError):
                os.kill(pid, 0)
        finally:
            try:
                os.kill(pid, sr.signal.SIGKILL)
            except ProcessLookupError:
                pass

    def test_log_failure_stops_descendant_when_leader_exits_on_sigterm(self):
        log = MagicMock()
        log.open.return_value.__enter__.return_value.write.side_effect = OSError('simulated log write failure')
        child_script = "import os,signal,time; signal.signal(signal.SIGTERM,signal.SIG_IGN); open('child-ready','w').write(str(os.getpid())); time.sleep(30)"
        parent_script = (
            "import json,os,subprocess,sys,time\n"
            "from pathlib import Path\n"
            "child=subprocess.Popen([sys.executable,'-c',sys.argv[1]])\n"
            "Path('child-pids.json').write_text(json.dumps([os.getpid(),child.pid]))\n"
            "deadline=time.monotonic()+5\n"
            "while not Path('child-ready').exists() and time.monotonic()<deadline: time.sleep(0.01)\n"
            "assert Path('child-ready').exists(), 'Test child did not start'\n"
            "print('child ready',flush=True)\n"
            "time.sleep(30)\n"
        )
        pids = self.base / 'child-pids.json'
        try:
            with redirect_stdout(io.StringIO()), self.assertRaisesRegex(OSError, 'simulated log write failure'):
                sr.run_nextflow([sys.executable, '-u', '-c', parent_script, child_script],
                                self.base, os.environ.copy(), log)
            parent_pid, child_pid = json.loads(pids.read_text())
            with self.assertRaises(ProcessLookupError):
                os.kill(parent_pid, 0)
            # A killed orphan may remain a zombie until init reaps it, but must
            # have stopped executing before the launcher releases its run lock.
            deadline = time.monotonic() + 2
            while time.monotonic() < deadline:
                stat = Path(f'/proc/{child_pid}/stat')
                if not stat.exists() or stat.read_text().split()[2] == 'Z':
                    break
                time.sleep(0.01)
            else:
                self.fail('Nextflow descendant remained running after cleanup')
        finally:
            if pids.exists():
                try:
                    os.killpg(json.loads(pids.read_text())[0], sr.signal.SIGKILL)
                except ProcessLookupError:
                    pass

    def test_previous_measurements_are_archived_before_resume(self):
        trace = self.base / 'sarek.trace.tsv'
        trace.write_text('task\tpeak_rss\nindex\t70 GB\n')
        sr.archive_run_files(self.base, ['sarek.trace.tsv', 'missing.html'])
        archived = list((self.base / 'history').glob('*/sarek.trace.tsv'))
        self.assertEqual(len(archived), 1)
        trace.write_text('new execution')
        self.assertIn('70 GB', archived[0].read_text())

    def test_failure_manifest_is_written_while_lock_is_held(self):
        output = self.base / 'locked-run'
        original_write = sr.write_json
        checked = []
        def write_and_check(path, value):
            if path.name == 'run_manifest.json' and value.get('status') == 'failed':
                with (output / '.run.lock').open('a') as probe:
                    with self.assertRaises(BlockingIOError):
                        sr.fcntl.flock(probe, sr.fcntl.LOCK_EX | sr.fcntl.LOCK_NB)
                checked.append(True)
            original_write(path, value)
        version = argparse.Namespace(stdout='nextflow version 25.10.4 build 11177', stderr='')
        with patch.object(sr.shutil, 'which', return_value='/bin/tool'), \
             patch.object(sr.subprocess, 'run', return_value=version), \
             patch.object(sr, 'run_nextflow', side_effect=OSError('execution failure')), \
             patch.object(sr, 'write_json', side_effect=write_and_check), patch('sys.stderr'):
            code = sr.main(['--input', str(self.samples), '--resources', str(self.resource), '--outdir', str(output)])
        self.assertEqual(code, 2)
        self.assertEqual(checked, [True])
        self.assertEqual(json.loads((output / 'run_manifest.json').read_text())['status'], 'failed')
        with (output / '.run.lock').open('a') as probe:
            sr.fcntl.flock(probe, sr.fcntl.LOCK_EX | sr.fcntl.LOCK_NB)

    def test_sv_rejects_wes_only(self):
        bed = self.base / 'capture.bed'
        bed.write_text('chr1\t10\t100\n')
        self.write_samples([{**self.row, 'assay': 'WES', 'target_bed': str(bed)}])
        with patch('sys.stderr'):
            self.assertEqual(sr.main(['--input', str(self.samples), '--resources', str(self.resource), '--outdir', str(self.base / 'dry'), '--sv', '--dry-run']), 2)


if __name__ == '__main__':
    unittest.main()
