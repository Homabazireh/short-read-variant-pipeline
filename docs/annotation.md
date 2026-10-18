# Offline annotation

`bin/annotate` accepts this CSV and annotates existing callsets, including imported DRAGEN or MegSAP output:

```csv
sample,label,vcf,variant_type
HG002,pipeline,/data/HG002.small.vcf.gz,small
HG002,sniffles,/data/HG002.sv.vcf.gz,sv
```

```bash
bin/annotate --input calls.csv --resources resources.json --outdir results/annotation --profile apptainer --cpus 8 --memory-gb 32
```

Use `--profile docker` locally or add `--executor slurm` on a cluster. Python 3.10+, Nextflow 25.10.4, the chosen container runtime, and shared filesystem access are required. `--resume` reuses a previous run with identical parameters, scientific input fingerprints and implementation hashes. New output directories must be empty; concurrent writers are prevented by an exclusive lock. Previous execution reports and traces are preserved under `execution-history/`. `--dry-run` checks real files, checksums, every REF allele and sample identity; it prints the invocation without running containers or writing results. Paths in CSV/JSON are resolved relative to those files. Do not modify inputs or resources while a run is active.

## Prepare pinned resources

Copy `examples/resources.example.json` and replace every placeholder with your actual resource paths and SHA256 checksums. The reference must be **uncompressed GRCh38**, have a matching `.fai`, and use the same contig naming and sequences as all callsets and annotations. A build label alone cannot establish assembly identity; the launcher also checks all input REF alleles and declared contig lengths. All VCFs must be sorted and contain exactly the named single sample. Split multisample callsets before importing them; do not silently relabel biological samples.

The image is Ensembl VEP **116.0**, pinned by OCI digest. Install the **Ensembl 116 GRCh38 cache** under `vep_cache/homo_sapiens/116_GRCh38`. The cache directory should be read-only. `cache_manifest` is an external text file with one `SHA256  relative/path` entry for every cache file. `cache_manifest_sha256` pins that manifest; preflight verifies its coverage and every listed file. Do not place the manifest inside the release directory. A checksum manifest can be generated with this Python snippet (replace paths):

```python
import hashlib
from pathlib import Path

cache = Path('/resources/vep_cache')
with open('/resources/vep_cache.116.GRCh38.sha256', 'w') as manifest:
    for path in sorted((cache / 'homo_sapiens/116_GRCh38').rglob('*')):
        if path.is_file():
            digest = hashlib.sha256()
            with path.open('rb') as source:
                for block in iter(lambda: source.read(8 * 1024 * 1024), b''):
                    digest.update(block)
            manifest.write(f'{digest.hexdigest()}  {path.relative_to(cache)}\n')
```

Provide dated ClinVar and versioned gnomAD VCF resources, with their actual SHA256, GRCh38 metadata, required INFO field names, and adjacent `.tbi` indexes. These resources must share the reference's contig names. The workflow uses VEP exact-allele matching for these **small-variant** resources. No liftover or contig renaming is performed. The example gnomAD version is 4.1; the version in your manifest must describe the file actually provided. Resource downloads are separate from analysis and may require substantial storage.

An unavailable source may explicitly be `null`. The resulting TSV records `unavailable`, distinct from `no_match` when an installed source has no matching allele. For SVs, these small-variant resources are reported as `not_applicable_small_variant_resource`; missing evidence never becomes a benign classification.

## Behavior and outputs

The workflow retains the original input files. Small variants pass through bcftools normalization and multiallelic splitting with strict REF checking; SVs retain their structural representation. VEP runs offline with gene symbols, HGVS, MANE, canonical flags, biotype, exon/intron numbers, protein identifiers, variant class, allele numbers and preferred-transcript **flags**. It retains all transcript rows rather than selecting one transcript.

Each `callsets/<sample>__<label>__<type>/` directory contains the prepared VCF/index, annotated VCF/index, a TSV with one row per CSQ consequence, and VEP's HTML summary when variants are present. The TSV retains every CSQ field, genotype, filter, full remaining INFO as JSON and source-availability flags. A valid empty VCF produces a header-only TSV and empty indexed VCF; it does not fabricate consequences. SV consequences describe overlap and consequences supported by VEP; this does not imply clinical interpretation, dosage sensitivity or supported interpretation of every breakend/rearrangement. Native warnings remain in Nextflow task logs.

`outputs.json` lists final products and checksums. `provenance.json` pins original callsets, reference, cache manifest, annotation databases, container digests, command and elapsed time. `trace.tsv`, `execution-report.html` and `timeline.html` record execution. Check these before using results; software availability and test success do not constitute clinical validation.

## Primary documentation

- [VEP installation and release 116 cache](https://www.ensembl.org/info/docs/tools/vep/script/vep_download.html)
- [VEP cache content and version matching](https://www.ensembl.org/info/docs/tools/vep/script/vep_cache.html)
- [VEP command options](https://www.ensembl.org/info/docs/tools/vep/script/vep_options.html)
- [VEP custom annotations](https://www.ensembl.org/info/docs/tools/vep/script/vep_custom.html)
- [bcftools normalization](https://samtools.github.io/bcftools/bcftools.html#norm)
