#!/usr/bin/env python3
"""Validated, resumable nf-core/sarek launcher. Python standard library only."""
from __future__ import annotations

import argparse
import csv
import fcntl
import gzip
import hashlib
import itertools
import json
import os
from pathlib import Path
import re
import shutil
import signal
import subprocess
import sys
from datetime import datetime, timezone

ROOT = Path(__file__).resolve().parents[1]
SAREK_VERSION = "3.10.0"
SAREK_REVISION = "8ccac7ad37b05dd792447763bf9671b719824587"
NEXTFLOW_VERSION = "25.10.4"
FIELDS = ["sample", "assay", "lane", "fastq_1", "fastq_2", "target_bed"]
SAFE_ID = re.compile(r"^[A-Za-z0-9][A-Za-z0-9_.-]*$")
SAFE_PATH = re.compile(r"^/[A-Za-z0-9_./+:-]+$")


class ValidationError(ValueError):
    pass


def sha256(path):
    digest = hashlib.sha256()
    with Path(path).open("rb") as handle:
        for chunk in iter(lambda: handle.read(4 * 1024 * 1024), b""):
            digest.update(chunk)
    return digest.hexdigest()


def checked_path(value, base, label, suffixes=()):
    if not isinstance(value, str) or not value:
        raise ValidationError(f"{label}: a nonempty file path is required")
    path = Path(value).expanduser()
    if not path.is_absolute():
        path = base / path
    path = path.resolve()
    # Upstream modules interpolate filenames into shell scripts without quoting.
    if not SAFE_PATH.fullmatch(str(path)):
        raise ValidationError(f"{label}: unsupported characters in path {path}; use simple paths without spaces or shell metacharacters")
    if not path.is_file() or path.stat().st_size == 0:
        raise ValidationError(f"{label}: missing or empty file: {path}")
    if suffixes and not any(str(path).endswith(s) for s in suffixes):
        raise ValidationError(f"{label}: expected suffix from {suffixes}: {path}")
    return path


def verify_hash(path, expected, label):
    if not isinstance(expected, str) or not re.fullmatch(r"[a-fA-F0-9]{64}", expected):
        raise ValidationError(f"{label}: a 64-character SHA-256 is required")
    actual = sha256(path)
    if actual.lower() != expected.lower():
        raise ValidationError(f"{label}: SHA-256 mismatch for {path}")
    return actual


def read_fai(path):
    lengths = {}
    with path.open() as handle:
        for line in handle:
            fields = line.rstrip("\n").split("\t")
            if len(fields) < 5:
                raise ValidationError(f"Malformed FASTA index: {path}")
            try:
                length, offset, width, byte_width = map(int, fields[1:5])
            except ValueError as exc:
                raise ValidationError(f"Noninteger FASTA index field: {path}") from exc
            if fields[0] in lengths or length <= 0 or offset < 0 or width <= 0 or byte_width < width:
                raise ValidationError(f"Invalid/duplicate contig in FASTA index: {fields[0]}")
            lengths[fields[0]] = length
    autosomes = {k: v for k, v in lengths.items() if re.fullmatch(r"(?:chr)?(?:[1-9]|1[0-9]|2[0-2])", k)}
    if not autosomes:
        raise ValidationError("Reference index has no recognized GRCh38 autosomes (1..22 or chr1..chr22)")
    if any(k.startswith("chr") for k in autosomes) and any(not k.startswith("chr") for k in autosomes):
        raise ValidationError("Reference mixes chr-prefixed and unprefixed autosomes")
    return lengths, autosomes


def validate_fai_against_fasta(fasta, fai):
    """Check names, lengths, byte offsets, and line widths against actual FASTA."""
    observed, layout = {}, {}
    current, position = None, 0
    with fasta.open("rb") as handle:
        for line in handle:
            if line.startswith(b">"):
                current = line[1:].split()[0].decode("ascii")
                if current in observed:
                    raise ValidationError(f"Duplicate FASTA contig: {current}")
                observed[current] = 0
            elif line.strip():
                if current is None:
                    raise ValidationError("Reference is not an uncompressed FASTA")
                bases = line.rstrip(b"\r\n")
                if current not in layout:
                    layout[current] = (position, len(bases), len(line))
                observed[current] += len(bases)
            position += len(line)
    if observed != fai:
        raise ValidationError("FASTA index contigs/lengths do not match the reference; rebuild with samtools faidx")
    with Path(str(fasta) + ".fai").open() as handle:
        for line in handle:
            fields = line.rstrip("\n").split("\t")
            if tuple(map(int, fields[2:5])) != layout.get(fields[0]):
                raise ValidationError("FASTA index offsets/line widths do not match the reference; rebuild with samtools faidx")


def vcf_contigs(path):
    contigs = {}
    with gzip.open(path, "rt") as handle:
        for line in handle:
            if line.startswith("##contig=<"):
                name = re.search(r"ID=([^,>]+)", line)
                length = re.search(r"length=(\d+)", line)
                if name:
                    contigs[name[1]] = int(length[1]) if length else None
            if line.startswith("#CHROM\t"):
                return contigs
    raise ValidationError(f"Missing VCF header: {path}")


def load_resources(path, caller):
    path = Path(path).resolve()
    with path.open() as handle:
        resource = json.load(handle)
    if resource.get("schema_version") != 1 or resource.get("genome_build") != "GRCh38":
        raise ValidationError("Resources must declare schema_version: 1 and genome_build: GRCh38")
    fasta = checked_path(resource.get("reference_fasta"), path.parent, "reference_fasta", (".fa", ".fna", ".fasta"))
    digest = verify_hash(fasta, resource.get("reference_sha256"), "reference_sha256")
    fai = checked_path(str(fasta) + ".fai", path.parent, "reference FASTA index")
    lengths, autosomes = read_fai(fai)
    validate_fai_against_fasta(fasta, lengths)
    known = []
    if caller == "haplotypecaller":
        records = resource.get("known_sites_vcf")
        if not isinstance(records, list) or not records:
            raise ValidationError("HaplotypeCaller with BQSR requires known_sites_vcf: [{path, sha256}, ...] for this GRCh38 reference")
        for number, entry in enumerate(records, 1):
            if not isinstance(entry, dict):
                raise ValidationError("Each known_sites_vcf entry must contain path and sha256")
            vcf = checked_path(entry.get("path"), path.parent, "known_sites_vcf", (".vcf.gz",))
            vcf_hash = verify_hash(vcf, entry.get("sha256"), "known-sites SHA-256")
            index = checked_path(str(vcf) + ".tbi", path.parent, "known-sites tabix index")
            contigs = vcf_contigs(vcf)
            if not contigs or not set(contigs).intersection(autosomes):
                raise ValidationError(f"Known-sites VCF lacks matching autosome headers: {vcf}")
            for contig, length in contigs.items():
                if contig not in lengths or (length is not None and length != lengths[contig]):
                    raise ValidationError(f"Known-sites contig mismatch: {vcf}: {contig}")
            known.append({"path": str(vcf), "sha256": vcf_hash, "index": str(index), "index_sha256": sha256(index)})
    return {"schema_version": 1, "genome_build": "GRCh38", "reference_fasta": str(fasta), "reference_sha256": digest,
            "reference_fai": str(fai), "reference_fai_sha256": sha256(fai), "lengths": lengths,
            "autosomes": autosomes, "known_sites_vcf": known}


def first_reads(path, limit=100):
    """Sample records for early pairing checks; FASTP validates the complete stream."""
    with gzip.open(path, "rt") as handle:
        for _ in range(limit):
            header = handle.readline().rstrip("\r\n")
            if not header:
                return
            seq = handle.readline().rstrip("\r\n")
            plus = handle.readline().rstrip("\r\n")
            qual = handle.readline().rstrip("\r\n")
            if not header.startswith("@") or not plus.startswith("+") or not seq or len(seq) != len(qual):
                raise ValidationError(f"Invalid FASTQ record: {path}")
            identifier = re.sub(r"/[12]$", "", header[1:].split()[0])
            yield identifier


def load_samples(path):
    path = Path(path).resolve()
    rows, seen, used_reads, per_sample = [], set(), set(), {}
    with path.open(newline="") as handle:
        reader = csv.DictReader(handle)
        if reader.fieldnames != FIELDS:
            raise ValidationError("Sample-sheet header must be exactly: " + ",".join(FIELDS))
        for number, row in enumerate(reader, 2):
            if None in row or any(value is None for value in row.values()):
                raise ValidationError(f"Malformed CSV row {number}")
            row = {k: v.strip() for k, v in row.items()}
            for key in ("sample", "lane"):
                if not SAFE_ID.fullmatch(row[key]):
                    raise ValidationError(f"Row {number}: invalid {key}; use letters, digits, _, . and -")
            key = (row["sample"], row["lane"])
            if key in seen:
                raise ValidationError(f"Duplicate sample/lane: {key}")
            seen.add(key)
            if row["assay"] not in ("WGS", "WES"):
                raise ValidationError(f"Row {number}: assay must be WGS or WES")
            reads = [checked_path(row[k], path.parent, k, (".fastq.gz", ".fq.gz")) for k in ("fastq_1", "fastq_2")]
            if reads[0] == reads[1] or any(str(p) in used_reads for p in reads):
                raise ValidationError(f"Row {number}: FASTQ paths must be distinct and cannot be reused across lanes")
            count = 0
            for left, right in itertools.zip_longest(first_reads(reads[0]), first_reads(reads[1])):
                if left != right:
                    raise ValidationError(f"Row {number}: paired-read identifiers/counts differ in first 100 reads")
                count += 1
            if not count:
                raise ValidationError(f"Row {number}: no FASTQ records")
            used_reads.update(str(p) for p in reads)
            row["fastq_1"], row["fastq_2"] = map(str, reads)
            if row["assay"] == "WES":
                row["target_bed"] = str(checked_path(row["target_bed"], path.parent, "WES target_bed", (".bed",)))
            elif row["target_bed"]:
                raise ValidationError("WGS rows must leave target_bed empty")
            group = (row["assay"], row["target_bed"])
            if row["sample"] in per_sample and per_sample[row["sample"]] != group:
                raise ValidationError(f"All lanes of {row['sample']} must use the same assay and target BED")
            per_sample[row["sample"]] = group
            rows.append(row)
    if not rows:
        raise ValidationError("Sample sheet is empty")
    return rows


def merge_intervals(intervals, lengths):
    order = {chrom: number for number, chrom in enumerate(lengths)}
    result = []
    for chrom, start, end in sorted(intervals, key=lambda x: (order[x[0]], x[1], x[2])):
        if result and result[-1][0] == chrom and start <= result[-1][2]:
            result[-1] = (chrom, result[-1][1], max(end, result[-1][2]))
        else:
            result.append((chrom, start, end))
    return result


def target_intervals(path, lengths, autosomes):
    intervals = []
    with Path(path).open() as handle:
        for number, line in enumerate(handle, 1):
            if not line.strip() or line.startswith(("#", "track ", "browser ")):
                continue
            parts = line.split()
            try:
                chrom, start, end = parts[0], int(parts[1]), int(parts[2])
            except (ValueError, IndexError) as exc:
                raise ValidationError(f"Invalid BED row {number}: {path}") from exc
            if chrom not in lengths or start < 0 or end <= start or end > lengths[chrom]:
                raise ValidationError(f"BED row {number} is outside reference bounds: {path}")
            if chrom in autosomes:
                intervals.append((chrom, start, end))
    if not intervals:
        raise ValidationError(f"No autosomal target intervals in {path}")
    return merge_intervals(intervals, lengths)


def intersect_intervals(left, right, lengths):
    by_contig = {}
    for chrom, start, end in right:
        by_contig.setdefault(chrom, []).append((start, end))
    result = []
    for chrom, start, end in left:
        for rstart, rend in by_contig.get(chrom, []):
            if rstart >= end:
                break
            if rend > start:
                result.append((chrom, max(start, rstart), min(end, rend)))
    return merge_intervals(result, lengths)


def make_groups(rows, resources, padding, regions=None):
    grouped = {}
    for row in rows:
        grouped.setdefault((row["assay"], row["target_bed"]), []).append(row)
    groups = []
    for (assay, target), samples in sorted(grouped.items()):
        name = "wgs" if assay == "WGS" else "wes-" + hashlib.sha256(target.encode()).hexdigest()[:12]
        reporting = target_intervals(target, resources["lengths"], resources["autosomes"]) if target else [(k, 0, v) for k, v in resources["autosomes"].items()]
        calling = merge_intervals([(c, max(0, start - padding), min(resources["lengths"][c], end + padding)) for c, start, end in reporting], resources["lengths"]) if target else reporting
        if regions:
            scope = target_intervals(regions, resources["lengths"], resources["autosomes"])
            reporting = intersect_intervals(reporting, scope, resources["lengths"])
            calling = intersect_intervals(calling, scope, resources["lengths"])
            if not reporting or not calling:
                raise ValidationError(f"Region scope does not overlap {name} reporting/calling intervals")
        groups.append({"id": name, "assay": assay, "target_bed": target,
                       "target_sha256": sha256(target) if target else None,
                       "samples": samples, "calling_intervals": calling, "reporting_intervals": reporting})
    return groups


def make_params(args, resources, group, outdir):
    base = outdir / "groups" / group["id"]
    params = {"input": str(base / "inputs" / "sarek.csv"), "outdir": str(base / "sarek"), "step": "mapping",
              "fasta": resources["reference_fasta"], "fasta_fai": resources["reference_fai"],
              "igenomes_ignore": True, "genome": "GRCh38", "custom_config_base": None,
              "aligner": "bwa-mem2", "tools": args.caller,
              "skip_tools": "baserecalibrator" if args.caller == "deepvariant" else "haplotypecaller_filter",
              "intervals": str(base / "inputs" / "calling.bed"), "wes": group["assay"] == "WES",
              "trim_fastq": True, "save_output_as_bam": True, "save_reference": True,
              "seq_platform": "ILLUMINA", "split_fastq": 0,
              "publish_dir_mode": "copy", "joint_germline": False, "download_cache": False}
    if args.caller == "haplotypecaller":
        params["known_indels"] = str(base / "inputs" / "known-sites" / "*.vcf.gz")
        params["known_indels_tbi"] = str(base / "inputs" / "known-sites" / "*.vcf.gz.tbi")
    return params


def nextflow_command(args, group, outdir, supplemental=False):
    base = outdir / "groups" / group["id"]
    name = "sv" if supplemental else "sarek"
    command = ["nextflow", "-log", str(base / f"{name}.nextflow.log"), "run"]
    if supplemental:
        command += [str(ROOT / "configs" / "shortread_sv.nf")]
    else:
        command += ["nf-core/sarek", "-r", SAREK_REVISION]
    command += ["-profile", args.profile, "-c", str(base / "execution.config"),
                "-params-file", str(base / f"{name}.params.json"), "-work-dir", str(base / f"work-{name}"),
                "-with-report", str(base / f"{name}.report.html"), "-with-trace", str(base / f"{name}.trace.tsv"),
                "-with-timeline", str(base / f"{name}.timeline.html")]
    if args.resume:
        command.append("-resume")
    return command


def execution_config(args, supplemental=False):
    config = [f"process.executor = '{args.executor}'", f"process.resourceLimits = [cpus: {args.cpus}, memory: {args.memory_gb}.GB]",
              "process.errorStrategy = 'terminate'", "process.maxRetries = 0",
              "trace.overwrite = true", "report.overwrite = true", "timeline.overwrite = true"]
    if args.executor == "local":
        config += [f"executor.cpus = {args.cpus}", f"executor.memory = '{args.memory_gb} GB'"]
    if args.queue:
        if not SAFE_ID.fullmatch(args.queue):
            raise ValidationError("Slurm queue must be a simple identifier")
        config.append(f"process.queue = '{args.queue}'")
    if supplemental:
        config += [f"process.cpus = {args.cpus}", f"process.memory = {args.memory_gb}.GB",
                   "profiles {", "  docker { docker.enabled = true; apptainer.enabled = false }",
                   "  apptainer { apptainer.enabled = true; apptainer.autoMounts = true; docker.enabled = false }", "}"]
    return "\n".join(config) + "\n"


def write_json(path, value):
    temporary = path.with_suffix(path.suffix + ".tmp")
    temporary.write_text(json.dumps(value, indent=2, sort_keys=True) + "\n")
    temporary.replace(path)


def write_csv(path, fields, rows):
    with path.open("w", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=fields)
        writer.writeheader()
        writer.writerows(rows)


def validate_output(outdir, resume, fingerprint):
    if outdir.is_symlink():
        raise ValidationError("Output directory cannot be a symbolic link")
    if outdir.exists():
        if not outdir.is_dir():
            raise ValidationError("Output path is not a directory")
        manifest = outdir / "run_manifest.json"
        if resume:
            if not manifest.is_file():
                raise ValidationError("--resume requires an existing run_manifest.json from this launcher")
            previous = json.loads(manifest.read_text())
            if previous.get("fingerprint") != fingerprint:
                raise ValidationError("Resume inputs, references, workflow code, or analysis settings differ; use a new output directory")
        elif any(outdir.iterdir()):
            raise ValidationError("Output directory is nonempty; use --resume for this exact run or choose a new directory")
    elif resume:
        raise ValidationError("--resume requires an existing output directory")


def inventory_files(base):
    return [str(p.resolve()) for p in sorted(base.rglob("*")) if p.is_file()]


def collect_calls(base, expected, caller):
    source = base / "sarek" / "csv" / "variantcalled.csv"
    if not source.is_file():
        raise ValidationError(f"Sarek exited without required output manifest: {source}")
    calls = []
    with source.open(newline="") as handle:
        reader = csv.DictReader(handle)
        if not {"sample", "variantcaller", "vcf"}.issubset(reader.fieldnames or []):
            raise ValidationError(f"Unexpected Sarek output manifest schema: {source}")
        for row in reader:
            if row["variantcaller"] != caller:
                continue
            sample = row["sample"]
            if sample not in expected:
                raise ValidationError(f"Unexpected sample in Sarek manifest: {sample}")
            vcf = checked_path(row["vcf"], source.parent, "Sarek output VCF", (".vcf.gz",))
            checked_path(str(vcf) + ".tbi", source.parent, "Sarek output VCF index")
            calls.append({"sample": sample, "label": caller, "vcf": str(vcf), "variant_type": "small"})
    if sorted(c["sample"] for c in calls) != sorted(expected):
        raise ValidationError("Sarek output manifest must contain exactly one final small-variant VCF per sample")
    return calls


def prepare_sv_input(base, expected):
    # TIDDIT uses duplicate-marked (pre-BQSR) BAMs, whose upstream manifest is stable.
    source = base / "sarek" / "csv" / "markduplicates_no_table.csv"
    if not source.is_file():
        raise ValidationError(f"Missing Sarek duplicate-marked alignment manifest: {source}")
    rows = []
    with source.open(newline="") as handle:
        for row in csv.DictReader(handle):
            sample = row.get("sample")
            if sample not in expected:
                raise ValidationError(f"Unexpected alignment sample: {sample}")
            bam = checked_path(row.get("bam"), source.parent, "TIDDIT input BAM", (".bam",))
            bai = checked_path(row.get("bai"), source.parent, "TIDDIT input BAI", (".bai",))
            rows.append({"sample": sample, "bam": str(bam), "bai": str(bai)})
    if sorted(r["sample"] for r in rows) != sorted(expected):
        raise ValidationError("Exactly one duplicate-marked BAM per sample is required for TIDDIT")
    write_csv(base / "inputs" / "sv.csv", ["sample", "bam", "bai"], rows)


def collect_sv_calls(base, expected):
    source = base / "sv" / "calls.csv"
    if not source.is_file():
        raise ValidationError(f"Missing TIDDIT calls manifest: {source}")
    with source.open(newline="") as handle:
        calls = list(csv.DictReader(handle))
    if sorted(r.get("sample", "") for r in calls) != sorted(expected):
        raise ValidationError("TIDDIT output manifest does not match requested samples")
    for row in calls:
        vcf = checked_path(row.get("vcf"), source.parent, "TIDDIT VCF", (".vcf.gz",))
        checked_path(str(vcf) + ".tbi", source.parent, "TIDDIT VCF index")
    return calls


def archive_run_files(directory, names):
    existing = [directory / name for name in names if (directory / name).is_file()]
    if not existing:
        return
    stamp = datetime.now(timezone.utc).strftime("%Y%m%dT%H%M%S.%fZ")
    history = directory / "history" / stamp
    history.mkdir(parents=True, exist_ok=False)
    for path in existing:
        shutil.copy2(path, history / path.name)


def run_nextflow(command, cwd, env, log):
    with log.open("a") as handle:
        process = subprocess.Popen(command, cwd=cwd, env=env, stdout=subprocess.PIPE,
                                   stderr=subprocess.STDOUT, text=True, start_new_session=True)
        try:
            for line in process.stdout:
                print(line, end="", flush=True)
                handle.write(line)
                handle.flush()
            code = process.wait()
            if code:
                raise subprocess.CalledProcessError(code, command)
        except BaseException:
            # Keep the run lock held until Nextflow and its process group have stopped.
            try:
                os.killpg(process.pid, signal.SIGTERM)
            except ProcessLookupError:
                pass
            try:
                process.wait(timeout=10)
            except subprocess.TimeoutExpired:
                pass
            # The leader can exit promptly while a descendant ignores SIGTERM.
            # Stop remaining group members even after the leader was reaped.
            try:
                os.killpg(process.pid, signal.SIGKILL)
            except ProcessLookupError:
                pass
            process.wait()
            raise
        finally:
            process.stdout.close()


def parser():
    p = argparse.ArgumentParser(description=__doc__)
    p.add_argument("--input", required=True, type=Path, help="sample,assay,lane,fastq_1,fastq_2,target_bed CSV")
    p.add_argument("--resources", required=True, type=Path)
    p.add_argument("--outdir", required=True, type=Path)
    p.add_argument("--profile", choices=["apptainer", "docker"], default="apptainer")
    p.add_argument("--executor", choices=["local", "slurm"], default="local")
    p.add_argument("--queue", help="Slurm partition (otherwise scheduler default)")
    p.add_argument("--cpus", type=int, default=8, help="CPU cap per task and total for local executor")
    p.add_argument("--memory-gb", type=int, default=96, help="Requested RAM ceiling per task and aggregate locally; local execution does not enforce RSS")
    p.add_argument("--caller", choices=["deepvariant", "haplotypecaller"], default="deepvariant")
    p.add_argument("--sv", action="store_true", help="Run TIDDIT with assembly on WGS groups")
    p.add_argument("--wes-padding", type=int, default=100, help="WES calling padding in bp; reporting BED remains unpadded")
    p.add_argument("--regions", type=Path, help="Optional BED scope for bounded integration runs; intersects autosomes and WES targets")
    p.add_argument("--resume", action="store_true")
    p.add_argument("--dry-run", action="store_true", help="Validate and print the exact plan without writing output")
    return p


def main(argv=None):
    args = parser().parse_args(argv)
    manifest = None
    run_lock = None
    try:
        if args.cpus < 1 or args.memory_gb < 1 or args.wes_padding < 0:
            raise ValidationError("cpus and memory-gb must be positive; wes-padding must be nonnegative")
        if args.queue and args.executor != "slurm":
            raise ValidationError("--queue requires --executor slurm")
        if args.outdir.is_symlink():
            raise ValidationError("Output directory cannot be a symbolic link")
        outdir = args.outdir.expanduser().resolve()
        if not SAFE_PATH.fullmatch(str(outdir)):
            raise ValidationError("Output directory must not contain spaces or shell metacharacters")
        resources = load_resources(args.resources, args.caller)
        rows = load_samples(args.input)
        regions = checked_path(str(args.regions), Path.cwd(), "regions", (".bed",)) if args.regions else None
        groups = make_groups(rows, resources, args.wes_padding, regions)
        if args.sv and not any(g["assay"] == "WGS" for g in groups):
            raise ValidationError("--sv requires at least one WGS group; TIDDIT is disabled for WES")
        config = execution_config(args)
        fingerprint_data = {"sarek_revision": SAREK_REVISION, "nextflow_version": NEXTFLOW_VERSION,
                            "resources": resources, "groups": groups, "caller": args.caller, "sv": args.sv,
                            "wes_padding": args.wes_padding, "regions": str(regions) if regions else None,
                            "regions_sha256": sha256(regions) if regions else None,
                            "launcher_sha256": sha256(Path(__file__)),
                            "sv_workflow_sha256": sha256(ROOT / "configs" / "shortread_sv.nf"),
                            "reads": [{"path": row[key], "size": Path(row[key]).stat().st_size,
                                       "mtime_ns": Path(row[key]).stat().st_mtime_ns} for row in rows for key in ("fastq_1", "fastq_2")]}
        fingerprint = hashlib.sha256(json.dumps(fingerprint_data, sort_keys=True).encode()).hexdigest()
        validate_output(outdir, args.resume, fingerprint)
        plan = {"schema_version": 1, "status": "planned", "fingerprint": fingerprint,
                "workflow": f"nf-core/sarek {SAREK_VERSION}", "provenance": fingerprint_data,
                "nextflow_version": NEXTFLOW_VERSION, "profile": args.profile, "executor": args.executor,
                "cpus": args.cpus, "memory_gb": args.memory_gb,
                "groups": [{"id": g["id"], "assay": g["assay"], "params": make_params(args, resources, g, outdir),
                            "command": nextflow_command(args, g, outdir),
                            "sv_command": nextflow_command(args, g, outdir, True) if args.sv and g["assay"] == "WGS" else None} for g in groups]}
        if args.dry_run:
            print(json.dumps(plan, indent=2))
            return 0
        for executable in ("nextflow", args.profile):
            if not shutil.which(executable):
                raise ValidationError(f"Required executable not found: {executable}")
        if args.executor == "slurm" and not shutil.which("sbatch"):
            raise ValidationError("Slurm executor requires sbatch in PATH")
        env = dict(os.environ, NXF_VER=NEXTFLOW_VERSION)
        version = subprocess.run(["nextflow", "-version"], env=env, capture_output=True, text=True, check=True)
        if not re.search(r"version\s+" + re.escape(NEXTFLOW_VERSION) + r"(?=\s|$)", version.stdout + version.stderr):
            raise ValidationError(f"Nextflow {NEXTFLOW_VERSION} is required")
        outdir.mkdir(parents=True, exist_ok=True)
        run_lock = (outdir / ".run.lock").open("a")
        try:
            fcntl.flock(run_lock, fcntl.LOCK_EX | fcntl.LOCK_NB)
        except BlockingIOError as exc:
            raise ValidationError("Another launcher is already running in this output directory") from exc
        candidate = outdir / "run_manifest.json"
        if candidate.exists() and not args.resume:
            raise ValidationError("A run manifest appeared after validation; refusing to overwrite")
        manifest = candidate
        if args.resume:
            archive_run_files(outdir, ["run_manifest.json", "outputs.json", "calls.csv"])
        plan["status"] = "running"
        plan["started_at"] = datetime.now(timezone.utc).isoformat()
        plan["nextflow_version_output"] = version.stdout.strip()
        write_json(manifest, plan)
        all_calls, outputs = [], []
        for group, item in zip(groups, plan["groups"]):
            base = outdir / "groups" / group["id"]
            inputs = base / "inputs"
            inputs.mkdir(parents=True, exist_ok=True)
            if args.resume:
                archive_run_files(base, ["execution.config", "sarek.params.json", "sv.params.json"] +
                                  [f"{name}.{suffix}" for name in ("sarek", "sv")
                                   for suffix in ("trace.tsv", "report.html", "timeline.html")])
            (base / "execution.config").write_text(config)
            for kind in ("calling", "reporting"):
                (inputs / f"{kind}.bed").write_text("".join(f"{c}\t{s}\t{e}\n" for c, s, e in group[f"{kind}_intervals"]))
            write_csv(inputs / "sarek.csv", ["patient", "sex", "status", "sample", "lane", "fastq_1", "fastq_2"],
                      [{"patient": row["sample"], "sex": "NA", "status": 0,
                        **{k: row[k] for k in ("sample", "lane", "fastq_1", "fastq_2")}} for row in group["samples"]])
            if resources["known_sites_vcf"]:
                known_dir = inputs / "known-sites"
                known_dir.mkdir(exist_ok=True)
                for number, site in enumerate(resources["known_sites_vcf"]):
                    for key, extension in (("path", ".vcf.gz"), ("index", ".vcf.gz.tbi")):
                        link = known_dir / f"known-{number:04d}{extension}"
                        if not link.exists():
                            link.symlink_to(site[key])
            write_json(base / "sarek.params.json", item["params"])
            run_nextflow(item["command"], base, env, base / "sarek.console.log")
            expected = set(row["sample"] for row in group["samples"])
            all_calls.extend(collect_calls(base, expected, args.caller))
            if item["sv_command"]:
                prepare_sv_input(base, expected)
                (base / "execution.config").write_text(execution_config(args, supplemental=True))
                write_json(base / "sv.params.json", {"input": str(inputs / "sv.csv"), "outdir": str(base / "sv"),
                                                     "fasta": resources["reference_fasta"], "fai": resources["reference_fai"],
                                                     "regions": str(inputs / "reporting.bed")})
                run_nextflow(item["sv_command"], base, env, base / "sv.console.log")
                all_calls.extend(collect_sv_calls(base, expected))
            outputs.append({"group": group["id"], "reporting_bed": str(inputs / "reporting.bed"),
                            "calling_bed": str(inputs / "calling.bed"), "sarek_files": inventory_files(base / "sarek"),
                            "sv_files": inventory_files(base / "sv") if item["sv_command"] else []})
        write_csv(outdir / "calls.csv", ["sample", "label", "vcf", "variant_type"], all_calls)
        write_json(outdir / "outputs.json", {"schema_version": 1, "calls": all_calls, "groups": outputs})
        plan["status"] = "complete"
        plan["completed_at"] = datetime.now(timezone.utc).isoformat()
        write_json(manifest, plan)
        print(f"Completed. Callset manifest: {outdir / 'calls.csv'}")
        return 0
    except KeyboardInterrupt:
        if manifest and manifest.exists():
            data = json.loads(manifest.read_text())
            data["status"], data["error"] = "failed", "Interrupted; Nextflow process group was terminated"
            write_json(manifest, data)
        print("ERROR: Interrupted; Nextflow was stopped. Use --resume to retry.", file=sys.stderr)
        return 130
    except (ValidationError, OSError, ValueError, subprocess.CalledProcessError) as exc:
        if manifest and manifest.exists():
            data = json.loads(manifest.read_text())
            data["status"], data["error"] = "failed", str(exc)
            write_json(manifest, data)
        print(f"ERROR: {exc}", file=sys.stderr)
        return 2
    finally:
        if run_lock is not None:
            run_lock.close()


if __name__ == "__main__":
    def interrupt_on_sigterm(signum, frame):
        raise KeyboardInterrupt
    signal.signal(signal.SIGTERM, interrupt_on_sigterm)
    raise SystemExit(main())
