# Reference and resource preparation

Keep data, indexes, models, work directories and annotation caches outside the
Git checkout. Resource JSON paths can be absolute or relative to that JSON.
The workflows do not guess an assembly from a filename or automatically lift
coordinates between assemblies.

## Reference

Use a single GRCh38 reference bundle consistently for alignment, calling,
annotation and benchmarking. The public-data integration script uses NCBI's
`GCA_000001405.15_GRCh38_no_alt_analysis_set.fna.gz`, verifies NCBI's published
MD5 checksum, decompresses it, indexes it, and records its SHA256.

```bash
samtools faidx /references/GRCh38.fa
sha256sum /references/GRCh38.fa
```

Set `reference_fasta` and `reference_sha256` in resources JSON. BAM/VCF contig
names and lengths, capture BEDs, truth resources and model metadata must match.
GRCh38 references with additional decoys are distinct bundles: never substitute
one based on its `hg38`/`GRCh38` label alone.

## VEP annotation cache

Use VEP **116**, matching the pinned container, and its GRCh38 indexed cache:

`https://ftp.ensembl.org/pub/release-116/variation/indexed_vep_cache/homo_sapiens_vep_116_GRCh38.tar.gz`

Extract it into a cache root containing `homo_sapiens/116_GRCh38`. This archive
is approximately 28 GB compressed; reserve additional space for extraction.
Record the archive digest and create a SHA256 manifest of the extracted cache:

```bash
cd /references/vep
find homo_sapiens/116_GRCh38 -type f -print0 | sort -z | xargs -0 sha256sum > cache.sha256
sha256sum cache.sha256
```

Put the cache path, manifest path and its checksum into the resource JSON as
described in [annotation](annotation.md). Each file is checked before use. VEP
runs offline against the supplied reference. See the [official cache guide](https://www.ensembl.org/info/docs/tools/vep/script/vep_cache.html).

ClinVar and gnomAD are separate, versioned, indexed VCF resources. Record source
URL, release, build, digest, and selected INFO fields. Make contig naming match
the reference through an explicitly recorded preparation step when needed.
For regional validation, a resource subset is valid only for its declared region;
do not reuse it as a whole-genome annotation database. See the [ClinVar FTP guide](https://www.ncbi.nlm.nih.gov/clinvar/docs/ftp_primer/)
and [gnomAD downloads](https://gnomad.broadinstitute.org/downloads).

## Public integration data

The acquisition script uses the official DeepVariant HG003 WGS, WES, HiFi and
ONT case-study sources and GIAB v4.2.1 truth. It fetches a bounded region with
flanking reads, retains full truth VCF context, records sources and hashes, and
reconstructs paired FASTQs for the short-read test:

```bash
python3 scripts/fetch_public_validation.py \
  --outdir /analysis/public-validation --platforms wgs,wes,hifi,ont
```

The script also runs the audited `prepare_public_fixture.py` helper, which
creates `wes.autosomal.samples.csv` with compatible capture targets and
`hifi.samples.csv` with a reference-compatible BAM and documented sample alias.
Original downloads remain intact and every transformation has a provenance JSON.
Use those derived sample sheets for the corresponding integration runs.

These are **integration fixtures**. Reconstructing FASTQ from a preprocessed BAM
cannot restore removed reads or reverse upstream processing. Additional audited
test-fixture preparation is recorded for the published HiFi decoy dictionary,
sample aliases and capture BED; see the measured validation report. A regional
run does not establish whole-genome accuracy, and a legacy ONT run does not
validate modern ONT chemistry.

Full benchmarks use the same command interfaces with complete public reads,
matched truth and its release-defined assessment regions. NIST's
[GIAB guidance](https://www.nist.gov/programs-projects/faqs-genome-bottle) describes
benchmark scope, exome targeting, and current small/SV releases. Published truth
is not a substitute for independent model-training holdout verification.
