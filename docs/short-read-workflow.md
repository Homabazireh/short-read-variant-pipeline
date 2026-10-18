# Short-read calling workflow

`bin/run` validates an Illumina germline sample sheet and a local resource manifest, then runs nf-core/sarek **3.10.0** at immutable commit `8ccac7ad37b05dd792447763bf9671b719824587` using Nextflow **25.10.4**. Sarek's pinned modules supply FastQC, fastp adapter trimming, BWA-MEM2 alignment, GATK duplicate marking, samtools/mosdepth QC, and MultiQC. DeepVariant chooses WGS or WES models per group; BQSR is skipped for DeepVariant. Software containers are the upstream pinned release versions. No moving institutional configuration is downloaded.

The workflow aligns to the complete supplied reference and calls **autosomes only**. Use an uncompressed GRCh38 FASTA with its adjacent `.fai`; generate it with `samtools faidx /references/GRCh38.fa`. Reference identity is checked against the user-provided SHA-256, and FAI contig names, lengths, byte offsets, and line widths are checked against the FASTA. The manifest's declared build must be GRCh38. A checksum verifies identity, not an independently authenticated assembly label; obtain the reference and its expected digest from a trusted source. Paths must not contain spaces or shell metacharacters, because upstream modules interpolate paths into shell scripts.

## Input contract

Use this exact CSV header:

```csv
sample,assay,lane,fastq_1,fastq_2,target_bed
HG002,WGS,L001,/data/HG002_R1.fastq.gz,/data/HG002_R2.fastq.gz,
HG003,WES,L001,/data/HG003_R1.fastq.gz,/data/HG003_R2.fastq.gz,/references/capture.bed
```

- Each row is one paired-end lane. `sample` identifies the biological sample; repeated samples are merged across unique lanes. IDs may contain letters, digits, dots, hyphens, and underscores and must begin with a letter or digit.
- Assay is `WGS` or `WES`. Every lane of one sample must share its assay and capture BED. Reused FASTQs and duplicate sample/lane pairs are rejected.
- FASTQs must be gzip compressed `.fastq.gz` or `.fq.gz`. The first 100 records are checked for formatting and matched pair identifiers; fastp processes and validates the full stream during execution. The launcher does not scan terabyte-scale FASTQs during preflight.
- `target_bed` is mandatory for WES and empty for WGS. BED coordinates are zero-based, half-open, must match reference contig names, and must lie within the reference. Nonautosomal intervals are excluded. Empty autosomal intersections fail.
- Relative file paths are resolved against the CSV or resource manifest containing them. Example paths and digest placeholders must be replaced before use.

Resources are JSON:

```json
{
  "schema_version": 1,
  "genome_build": "GRCh38",
  "reference_fasta": "/references/GRCh38.fa",
  "reference_sha256": "64_HEXADECIMAL_CHARACTERS_FROM_SHA256SUM"
}
```

Annotation resources use the additional fields described in the annotation documentation and are not downloaded or required by this calling entry point.

## Running

Install Python 3.10+, Nextflow 25.10.4 and its supported Java runtime, and either Apptainer or Docker. Container execution must be permitted on the execution host. Sarek source and container downloads require network access on a first run; resource data are supplied locally. Human BWA-MEM2 indexing and whole-genome calling can need substantially more RAM/disk than small integration runs; set resource limits for your data and retain adequate work storage.

```bash
bin/run --input samples.csv --resources resources.json \
  --outdir /analysis/run01 --profile apptainer \
  --cpus 16 --memory-gb 96 --dry-run

bin/run --input samples.csv --resources resources.json \
  --outdir /analysis/run01 --profile apptainer \
  --cpus 16 --memory-gb 96
```

`--dry-run` validates inputs and prints the exact plan without creating results or requiring a container runtime. CPU and RAM limits cap each process resource request; under the local executor they additionally limit aggregate scheduled requests. They do not enforce actual memory use with local cgroups. The default RAM request ceiling is 96 GB because initial human BWA-MEM2 indexing exceeded 64 GB in integration testing. Sarek process defaults remain in use under those request ceilings. For Slurm add `--executor slurm --queue YOUR_PARTITION`; shared paths must be visible on compute nodes. Docker requires a functioning daemon; the existence of its CLI alone does not establish this.

WGS samples and each distinct WES capture BED run as independent groups. WES calling targets are padded by **100 bp**, clipped to contig boundaries, and merged; set `--wes-padding` to change this. `reporting.bed` always records the original unpadded autosomal capture targets. Final VCFs include the padded calling region. Downstream exome accuracy and target summaries must use the unpadded reporting BED intersected with applicable benchmark confidence regions. Whole-genome alignment metrics and Sarek's padded-interval metrics are kept separately from those downstream target metrics.

`--regions subset.bed` restricts calling and reporting intervals to a bounded autosomal scope, suitable for integration runs and regional comparisons. The complete input reference is still used for alignment. This option does **not** constitute a full-genome benchmark; its path, digest, and interval scope are recorded. For WES it intersects the supplied capture intervals rather than replacing the capture-kit requirement.

### GATK comparison

Select `--caller haplotypecaller` to enable BQSR plus GATK HaplotypeCaller. Add `known_sites_vcf` to resources as a nonempty array of objects containing `path` and `sha256`. Each path must be an indexed `.vcf.gz` with an adjacent `.tbi`, compatible contig names/lengths, and the same GRCh38 reference. Supply established known polymorphism resources, not your test sample's calls or benchmark truth. The wrapper stages these resources and passes their glob as Sarek's `known_indels` BQSR resource input; that upstream name does not restrict the records to indels.

GATK's optional CNN truth-based tranche filtering is explicitly disabled because generic BQSR sites are not an adequate substitute for validated variant-filtering training resources. The resulting GATK comparison callset is the original HaplotypeCaller output; the default DeepVariant output retains caller FILTER values. Joint genotyping and cohort recalibration are outside this initial implementation.

### Structural variants

Add `--sv` to run TIDDIT **3.9.5** for WGS groups. WES-only SV requests fail. Sarek 3.10.0 intentionally supplies no BWA index to its TIDDIT module and therefore adds `--skip_assembly`. This repository instead runs `configs/shortread_sv.nf` after the small-variant workflow, builds a classic BWA index, and runs TIDDIT **with assembly enabled** on the duplicate-marked BAMs identified by Sarek's alignment manifest. The BWA index and raw TIDDIT VCF/ploidy outputs are preserved. The final compressed/indexed SV VCF is restricted to the reporting intervals with bcftools; the complete raw SV output is unvalidated outside the declared autosomal scope. Breakend partners can fall outside the first-breakpoint reporting scope and require explicit benchmark eligibility assessment.

The supplemental workflow uses independently pinned BWA/HTSlib/samtools, TIDDIT, and bcftools container images from the same Sarek release. It records software version outputs and execution traces.

## Outputs, provenance, and recovery

- `calls.csv`: downstream callset interface, with `sample,label,vcf,variant_type`; types are `small` and `sv`. Paths are taken from upstream output manifests and checked before success is reported.
- `outputs.json`: grouped file inventory, original Sarek output paths, and separate calling/reporting BED paths.
- `run_manifest.json`: status, pinned revision, validated resource hashes, sample metadata, input file sizes/mtime, regional scope, generated parameters, exact commands, execution limits, and timing. Raw read content hashes are not computed; read size/mtime detect routine changes for resume but are not content authentication.
- `groups/<group>/sarek`: Sarek's original alignment, variant, QC, version, reference, and restart CSV outputs are preserved. Alignments are BAM with index. Per-group `*.trace.tsv`, execution reports, timelines, and console/Nextflow logs are retained.
- `groups/<group>/inputs/{calling,reporting}.bed`: exact zero-based interval scope. Optional SV outputs are in `groups/<group>/sv`.

Output directories are never silently overwritten. Resume with the same input/resource/settings and `--resume`; retaining the work directories and Nextflow cache is necessary. The launcher rejects changes to scientific inputs or workflow code when resuming. Execution resource limits/profile may be changed on resume; cached task reuse is managed by Nextflow. A lock prevents concurrent launcher writes; interruption terminates and waits for the Nextflow process group before releasing it. Previous run manifests and execution reports/traces are archived under `history/` on resume before current copies are updated. Failed runs retain logs and a failed manifest, and never create a completed callset manifest. Annotation and benchmarking are separate commands, so failed external resource preparation cannot be mistaken for completed variant interpretation.

## Verification and references

Run launcher tests with `python3 -m unittest discover -s tests -p 'test_shortread*.py' -v`. Tests cover checksums, reference mismatch, pairing, WES padding, scoped intervals, unsafe inputs, required BQSR resources, restart safety, revision pinning, and output-manifest completeness. Unit tests and workflow compilation do not measure variant-calling accuracy; see the repository validation status for actual integration and benchmark evidence.

Upstream sources: [Sarek 3.10.0 usage](https://nf-co.re/sarek/3.10.0/docs/usage), [pinned workflow](https://github.com/nf-core/sarek/blob/8ccac7ad37b05dd792447763bf9671b719824587/workflows/sarek.nf), [TIDDIT configuration](https://github.com/nf-core/sarek/blob/8ccac7ad37b05dd792447763bf9671b719824587/conf/modules/tiddit.config), [DeepVariant model selection](https://github.com/nf-core/sarek/blob/8ccac7ad37b05dd792447763bf9671b719824587/conf/modules/deepvariant.config).
