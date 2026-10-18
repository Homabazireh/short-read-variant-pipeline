#!/usr/bin/env python3
"""Acquire bounded public HG003 data from the official DeepVariant case studies.

This is an integration fixture, not a whole-genome accuracy benchmark. Original
source URLs, HTTP metadata, commands, and SHA256 digests are recorded. FASTQ
reconstruction cannot undo preprocessing of the source BAM.
"""
from __future__ import annotations

import argparse
from contextlib import contextmanager
import fcntl
import concurrent.futures
import csv
import gzip
import hashlib
import json
import os
from pathlib import Path
import re
import shutil
import subprocess
import tempfile
import urllib.request

REFERENCE_BASE = "https://ftp.ncbi.nlm.nih.gov/genomes/all/GCA/000/001/405/GCA_000001405.15_GRCh38/seqs_for_alignment_pipelines.ucsc_ids/"
REFERENCE_NAME = "GCA_000001405.15_GRCh38_no_alt_analysis_set.fna.gz"
TRUTH_BASE = "https://ftp-trace.ncbi.nlm.nih.gov/giab/ftp/release/AshkenazimTrio/HG003_NA24149_father/NISTv4.2.1/GRCh38/"
TRUTH_NAME = "HG003_GRCh38_1_22_v4.2.1_benchmark"
READS = {
    "wgs": "https://storage.googleapis.com/deepvariant/case-study-testdata/HG003.novaseq.pcr-free.35x.dedup.grch38_no_alt.chr20.bam",
    "wes": "https://storage.googleapis.com/deepvariant/exome-case-study-testdata/HG003.novaseq.wes_idt.100x.dedup.bam",
    "hifi": "https://storage.googleapis.com/deepvariant/pacbio-case-study-testdata/HG003.SPRQ.pacbio.GRCh38.nov2024.chr20.bam",
    "ont": "https://storage.googleapis.com/deepvariant/ont-case-study-testdata/HG003_R104_sup_merged.80x.chr20.bam",
}
CAPTURE = "https://storage.googleapis.com/deepvariant/exome-case-study-testdata/idt_capture_novogene.grch38.bed"


def digest(path: Path, algorithm: str = "sha256") -> str:
    h = hashlib.new(algorithm)
    with path.open("rb") as handle:
        for block in iter(lambda: handle.read(8 * 1024 * 1024), b""):
            h.update(block)
    return h.hexdigest()


def immutable_json(path: Path, value: dict) -> None:
    if path.exists():
        if json.loads(path.read_text()) != value:
            raise ValueError(f"Existing metadata differs; use a new acquisition directory: {path}")
        return
    path.write_text(json.dumps(value, indent=2) + "\n")


def immutable_text(path: Path, value: str) -> None:
    if path.exists():
        if path.read_bytes() != value.encode():
            raise ValueError(f"Existing derived file differs: {path}")
        return
    path.write_text(value)


@contextmanager
def acquisition_lock(out: Path, marker_name: str, spec: dict):
    """Reject unrelated destinations and serialize all fixture writes."""
    marker = out / marker_name
    if out.exists() and not out.is_dir():
        raise ValueError(f"Acquisition destination is not a directory: {out}")
    if out.exists() and not marker.is_file() and any(p.name != '.acquisition.lock' for p in out.iterdir()):
        raise ValueError(f"Nonempty directory has no matching acquisition marker: {out}")
    out.mkdir(parents=True, exist_ok=True)
    with (out / '.acquisition.lock').open('a+') as handle:
        try:
            fcntl.flock(handle.fileno(), fcntl.LOCK_EX | fcntl.LOCK_NB)
        except BlockingIOError as exc:
            raise ValueError(f"Another acquisition is writing to {out}") from exc
        try:
            if not marker.is_file() and any(p.name != '.acquisition.lock' for p in out.iterdir()):
                raise ValueError(f"Acquisition directory changed before its lock was acquired: {out}")
            immutable_json(marker, spec)
            yield
        finally:
            fcntl.flock(handle.fileno(), fcntl.LOCK_UN)


def validate_derived(marker: Path, source_identity: dict, outputs: list[Path]) -> bool:
    """Resume only complete, checksum-matching derived outputs."""
    if not marker.exists():
        if any(path.exists() for path in outputs):
            raise ValueError(f"Untracked/partial derived output; use a new directory: {marker}")
        return False
    recorded = json.loads(marker.read_text())
    if recorded.get('source_identity') != source_identity:
        raise ValueError(f"Derived input identity changed: {marker}")
    expected = recorded.get('output_sha256', {})
    for path in outputs:
        if not path.is_file() or expected.get(path.name) != digest(path):
            raise ValueError(f"Derived output missing or checksum changed: {path}")
    return True


def record_derived(marker: Path, source_identity: dict, outputs: list[Path], **details) -> None:
    immutable_json(marker, {**details, 'source_identity': source_identity,
                           'output_sha256': {path.name: digest(path) for path in outputs}})


def ensure_reference_fasta(compressed: Path, fasta: Path) -> str:
    """Compare an existing FASTA to decompressed publisher bytes; never rebless it."""
    expected = hashlib.sha256()
    with gzip.open(compressed, 'rb') as source:
        for block in iter(lambda: source.read(8 * 1024 * 1024), b''):
            expected.update(block)
    checksum = expected.hexdigest()
    marker = fasta.with_name(fasta.name + '.derived.json')
    identity = {'compressed_sha256': digest(compressed), 'uncompressed_sha256': checksum}
    if marker.exists() and json.loads(marker.read_text()) != identity:
        raise ValueError(f"Reference derivation identity changed: {marker}")
    if fasta.exists():
        if not fasta.is_file() or digest(fasta) != checksum:
            raise ValueError(f"Existing FASTA differs from the verified compressed reference: {fasta}")
    else:
        with tempfile.TemporaryDirectory(prefix='.reference-', dir=fasta.parent) as directory:
            temporary = Path(directory) / fasta.name
            with gzip.open(compressed, 'rb') as source, temporary.open('wb') as target:
                shutil.copyfileobj(source, target, 8 * 1024 * 1024)
            if digest(temporary) != checksum:
                raise ValueError('Decompressed reference checksum mismatch')
            temporary.replace(fasta)
    immutable_json(marker, identity)
    return checksum


def ensure_fasta_index(fasta: Path, out: Path) -> None:
    index = Path(str(fasta) + '.fai')
    # Generate beside a symlink in a temporary directory, then compare before
    # publishing; this also verifies previously existing untracked FAI files.
    with tempfile.TemporaryDirectory(prefix='.fai-', dir=out) as directory:
        staged = Path(directory) / fasta.name
        staged.symlink_to(fasta)
        run(['samtools', 'faidx', str(staged)], out)
        generated = Path(str(staged) + '.fai')
        if index.exists():
            if digest(index) != digest(generated):
                raise ValueError(f'Existing FASTA index differs from the reference: {index}')
        else:
            generated.replace(index)


def fetch(url: str, target: Path) -> dict:
    metadata_path = target.with_name(target.name + ".source.json")
    if target.exists():
        if not metadata_path.exists():
            raise ValueError(f"Existing download has no provenance: {target}")
        meta = json.loads(metadata_path.read_text())
        if meta["url"] != url or digest(target) != meta["sha256"]:
            raise ValueError(f"Existing download differs from provenance: {target}")
        return meta
    print(f"Downloading {target.name}", flush=True)
    partial = target.with_name(target.name + ".partial")
    with urllib.request.urlopen(url, timeout=180) as response, partial.open("wb") as output:
        meta = {"url": url, "etag": response.headers.get("ETag"),
                "last_modified": response.headers.get("Last-Modified")}
        shutil.copyfileobj(response, output, 8 * 1024 * 1024)
    meta.update(sha256=digest(partial), bytes=partial.stat().st_size)
    if metadata_path.exists():
        previous = json.loads(metadata_path.read_text())
        if previous.get('url') != url or previous.get('sha256') != meta['sha256']:
            raise ValueError(f"Re-downloaded source differs from its original provenance: {target}")
        meta = previous
    partial.replace(target)
    immutable_json(metadata_path, meta)
    return meta


def run(command: list[str], cwd: Path) -> None:
    print("Running: " + " ".join(command), flush=True)
    subprocess.run(command, cwd=cwd, check=True)


def main() -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--outdir", required=True, type=Path)
    parser.add_argument("--region", default="chr20:10000001-10200000",
                        help="1-based inclusive chr20 interval; default is 200 kb")
    parser.add_argument("--platforms", default="wgs,wes,hifi,ont")
    parser.add_argument("--threads", default=4, type=int)
    args = parser.parse_args()
    m = re.fullmatch(r"chr20:(\d+)-(\d+)", args.region)
    if not m or not 1 <= int(m[1]) <= int(m[2]) <= 64444167:
        parser.error("--region must be a valid GRCh38 chr20 interval")
    if not 1 <= args.threads <= 32:
        parser.error("--threads must be between 1 and 32")
    platforms = args.platforms.split(",")
    if any(p not in READS for p in platforms) or len(set(platforms)) != len(platforms):
        parser.error("--platforms must contain unique names from wgs,wes,hifi,ont")
    if not shutil.which("samtools"):
        parser.error("samtools is required (including HTTPS support)")
    if args.outdir.is_symlink():
        parser.error("acquisition output directory cannot be a symbolic link")
    out = args.outdir.expanduser().resolve()
    spec = {"region": args.region, "platforms": platforms, "sample": "HG003"}
    with acquisition_lock(out, 'acquisition.json', spec):
        return acquire(args, out, platforms, m)


def acquire(args, out, platforms, m):
    if shutil.disk_usage(out).free < 8 * 1024**3:
        raise ValueError("at least 8 GiB free disk is required for bounded validation data")
    downloads = [(REFERENCE_BASE + REFERENCE_NAME, out / REFERENCE_NAME),
                 (REFERENCE_BASE + "md5checksums.txt", out / "reference.md5checksums.txt")]
    for suffix in (".vcf.gz", ".vcf.gz.tbi", "_noinconsistent.bed"):
        downloads.append((TRUTH_BASE + TRUTH_NAME + suffix, out / (TRUTH_NAME + suffix)))
    if "wes" in platforms:
        downloads.append((CAPTURE, out / "idt_capture_novogene.grch38.bed"))
    with concurrent.futures.ThreadPoolExecutor(max_workers=3) as pool:
        metadata = list(pool.map(lambda pair: fetch(*pair), downloads))
    checksums = (out / "reference.md5checksums.txt").read_text().splitlines()
    expected_md5 = next(line.split()[0] for line in checksums
                        if line.split()[-1] == "./" + REFERENCE_NAME)
    if digest(out / REFERENCE_NAME, "md5") != expected_md5:
        raise ValueError("Reference failed the NCBI-published MD5 checksum")
    fasta = out / "GRCh38_no_alt_analysis_set.fasta"
    reference_hash = ensure_reference_fasta(out / REFERENCE_NAME, fasta)
    ensure_fasta_index(fasta, out)
    immutable_json(out / "resources.json", {
        "schema_version": 1, "genome_build": "GRCh38",
        "reference_fasta": str(fasta), "reference_sha256": reference_hash,
    })
    immutable_text(out / "evaluation.bed", f"chr20\t{int(m[1])-1}\t{m[2]}\n")
    # Include flanking reads; score only the requested interior interval.
    padded = f"chr20:{max(1, int(m[1])-50000)}-{min(64444167, int(m[2])+50000)}"
    derived = []
    for platform in platforms:
        bam = out / f"HG003.{platform}.subset.bam"
        marker = bam.with_suffix(".complete.json")
        command = ["samtools", "view", "-@", str(args.threads), "-b"]
        if platform in ("wgs", "wes"):
            command += ["--fetch-pairs"]
        command += ["-o", str(bam), READS[platform], padded]
        index = Path(str(bam) + '.bai')
        identity = {'command': command}
        if not validate_derived(marker, identity, [bam, index]):
            with tempfile.TemporaryDirectory(prefix=f'.{platform}-bam-', dir=out) as directory:
                temporary = Path(directory) / bam.name
                temporary_command = list(command)
                temporary_command[temporary_command.index('-o') + 1] = str(temporary)
                run(temporary_command, out)
                run(["samtools", "quickcheck", "-v", str(temporary)], out)
                run(["samtools", "index", "-@", str(args.threads), str(temporary)], out)
                temporary.replace(bam)
                Path(str(temporary) + '.bai').replace(index)
            record_derived(marker, identity, [bam, index], command=command)
        recorded = json.loads(marker.read_text())
        derived.append({"platform": platform, "source_url": READS[platform], **recorded})
        header = subprocess.check_output(["samtools", "view", "-H", str(bam)], text=True)
        immutable_text(out / f"HG003.{platform}.header.sam", header)
        if platform in ("wgs", "wes"):
            fq1, fq2 = out / f"HG003.{platform}.R1.fastq.gz", out / f"HG003.{platform}.R2.fastq.gz"
            singleton = out / f"HG003.{platform}.singletons.fastq.gz"
            fq_marker = out / f"HG003.{platform}.fastq.json"
            fq_identity = {'source_bam_sha256': digest(bam), 'method': 'samtools collate then fastq; primary paired reads; singletons separate'}
            if not validate_derived(fq_marker, fq_identity, [fq1, fq2, singleton]):
                with tempfile.TemporaryDirectory(prefix=f'.{platform}-fastq-', dir=out) as directory:
                    temporary = [Path(directory) / path.name for path in (fq1, fq2, singleton)]
                    collate = subprocess.Popen(["samtools", "collate", "-u", "-O", "-@", str(args.threads), str(bam)], stdout=subprocess.PIPE)
                    try:
                        fastq = subprocess.run(["samtools", "fastq", "-@", str(args.threads), "-n", "-1", str(temporary[0]), "-2", str(temporary[1]), "-0", os.devnull, "-s", str(temporary[2]), "-"], stdin=collate.stdout)
                        collate.stdout.close()
                        if fastq.returncode or collate.wait():
                            raise RuntimeError("Paired FASTQ reconstruction failed")
                    except BaseException:
                        collate.kill()
                        collate.wait()
                        raise
                    for source, destination in zip(temporary, (fq1, fq2, singleton)):
                        source.replace(destination)
                record_derived(fq_marker, fq_identity, [fq1, fq2, singleton])
            import io
            sample_text = io.StringIO(newline='')
            writer = csv.writer(sample_text)
            writer.writerow(["sample", "assay", "lane", "fastq_1", "fastq_2", "target_bed"])
            writer.writerow(["HG003", platform.upper(), "L001", fq1, fq2,
                             str(out / "idt_capture_novogene.grch38.bed") if platform == "wes" else ""])
            immutable_text(out / f"{platform}.samples.csv", sample_text.getvalue())
    # These publisher-specific transformations are explicit and audited; study
    # data are never passed through this fixture-only helper by the pipelines.
    from prepare_public_fixture import prepare_hifi, prepare_wes
    contigs = {line.split('\t')[0]: int(line.split('\t')[1])
               for line in Path(str(fasta) + '.fai').read_text().splitlines()}
    prepare_hifi(out, contigs)
    prepare_wes(out, contigs)
    immutable_json(out / "provenance.json", {
        "purpose": "bounded public-data integration validation; not whole-genome accuracy",
        "sample": "HG003", "region": args.region, "read_fetch_region": padded,
        "reference_sha256": reference_hash, "downloads": metadata, "derived": derived,
        "limitations": ["short-read FASTQs reconstructed from already-preprocessed BAMs",
                        "no independent generalization claim without model training audit"],
    })
    print(f"Prepared public HG003 integration data: {out}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
