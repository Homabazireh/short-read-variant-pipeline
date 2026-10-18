# Measured validation

Validation date: **2026-10-18**. Real public HG003 regional WGS and WES runs
completed on Linux x86-64 with Apptainer, Nextflow 25.10.4, Java 21, and the
full GRCh38 no-alt analysis reference. The host has an AMD EPYC 9454P and
approximately 251 GiB RAM. The workflow foundation is nf-core/sarek 3.10.0,
commit `8ccac7ad37b05dd792447763bf9671b719824587`.

These are integration and regression measurements. Model-training overlap with
HG003 was not audited. They do not establish whole-genome or whole-exome
accuracy, clinical validity, or equivalence to DRAGEN. No DRAGEN installation
or DRAGEN callset was benchmarked.

## Executed calling workflows

Both runs completed all 17 real Nextflow tasks: reference preparation, FastQC,
fastp, BWA-MEM2 alignment, GATK duplicate marking, coverage/alignment QC,
DeepVariant with the correct WGS/WES model, variant QC, and MultiQC. Caller
VCFs, indexes, alignment files, reports, and `calls.csv` were checked.

| Assay | Input read pairs | Pairs after fastp | Original variant records | Wrapper wall time including index/setup |
| --- | ---: | ---: | ---: | ---: |
| WGS fixture | 40,170 | 39,745 | 641 | 1,166.4 s |
| WES fixture | 1,485 | 1,481 | 5 | 1,057.2 s |

Input FASTQs were reconstructed from the official DeepVariant preprocessed
HG003 BAM examples with flanking reads around `chr20:10000001-10200000`.
Reconstruction does not restore previously removed reads or undo the
publisher's upstream processing. WES capture targets were explicitly restricted
to reference-compatible autosomes; original targets and transformation
provenance were retained. The default 100 bp padding yielded 4,131 calling
bases; **the unpadded reporting intervals contain 2,331 bases**.

Each run built a full-reference BWA-MEM2 index (923 and 917 seconds). Both
reported 69.3 GB peak resident memory in the Nextflow trace. The runs requested
8 CPUs and 64 GB on an adequately sized local host; local resource requests do
not enforce process RSS. The released default is therefore **96 GB**. These
regional timings include initial indexing and image setup, and are not
whole-genome performance estimates. DeepVariant task runtimes were 33.8 s
(WGS) and 30.8 s (WES).

The [per-task measurements](validation/short-read-execution.csv) and
[calling provenance](validation/short-read-calling-provenance.json) retain exact
versions, hashes, assessment intervals, task counts, and original timings.
Post-run changes tightened reference-index verification, interruption cleanup,
version checks, and resume-history preservation; scientific caller parameters
were unchanged. The final launcher correctly rejected resuming the earlier
implementation because its code fingerprint differed, preserving the completed
manifest. These changes and both code hashes are recorded in the evidence.

## Small-variant benchmarks

Both evaluations completed with **GIAB HG003 GRCh38 v4.2.1**, hap.py 0.3.12,
and RTG vcfeval. The full truth VCF was retained. Confidence BEDs were
intersected with the WGS evaluation interval or the WES unpadded reporting
targets before evaluation. hap.py's standard insertion-boundary adjustment
was enabled: the input confidence BED contains 194,051 bases for WGS and
2,331 for WES; the native adjusted confidence sizes are 194,597 and 2,333 bases.
Both sizes and the adjustment policy are retained in the benchmark evidence.

| Assay | Type | Truth TP | Query TP | FN | FP | Precision | Recall | F1 |
| --- | --- | ---: | ---: | ---: | ---: | ---: | ---: | ---: |
| WGS | SNP | 389 | 389 | 0 | 0 | 1.000000 | 1.000000 | 1.000000 |
| WGS | INDEL | 62 | 68 | 0 | 0 | 1.000000 | 1.000000 | 1.000000 |
| WES | SNP | 2 | 2 | 0 | 0 | 1.000000 | 1.000000 | 1.000000 |

These are PASS rows from tiny regional fixtures, not whole-assay accuracy
estimates. The WES benchmark assesses only two truth SNPs; three other query
SNPs are unassessed, and **no WES indel estimate is available**. No genotype
errors were reported for these assessed rows. Truth/query TP may differ for
equivalent haplotypes with different variant representations; precision uses
query TP and recall uses truth TP. Unassessed calls are separate from FP.

See [full-precision metrics](validation/short-read-metrics.csv) and
[benchmark provenance](validation/short-read-provenance.json) for resource
hashes, native confidence sizes, tool/container identities, and scope.

## Annotation and structural-variant checks

Real offline VEP 116 annotation completed on both callsets with the complete
Ensembl 116 GRCh38 cache and actual ClinVar 2026-10-11 / gnomAD 4.1 genomes
resources restricted to the test intervals. Reference/cache/resource hashes,
allele matches, HGVS/MANE fields, and transcript consequences were checked.
The regional annotation resources must not be reused as whole-genome databases.

| Callset | Records after normalization | Consequence rows |
| --- | ---: | ---: |
| WGS DeepVariant | 656 | 14,008 |
| WES DeepVariant | 5 | 124 |

Normalization can split multiallelic records. Match counts in the
[annotation evidence](validation/short-read-annotation.json) count consequence
rows, not unique variants. All transcript consequences remain available;
absent resources and an available resource with no allele match are distinct.

The optional TIDDIT stage completed a separate real-data interface test using
an explicitly prepared chr20-only reference and compatible HG003 BAM. Classic
BWA indexing, TIDDIT 3.9.5 including assembly, and scoped/indexed VCF output
passed; the resulting VCF was valid and empty. All three tasks were cached on
resume. The fixture transformation and output hashes are recorded in the
calling provenance. This test is not a full-reference end-to-end SV benchmark.

The real Truvari workflow also passed intentionally perturbed synthetic tests
(2 SV TP, 1 FN, 1 FP, and 1 genotype error), empty-callset cases, and size-bin
accounting. **No real GIAB SV accuracy result is claimed.** A real SV benchmark
needs release-specific truth regions and an explicit matching policy.

## Reproduction and remaining validation

Install the README prerequisites, then acquire the bounded public fixtures:

```bash
python3 scripts/fetch_public_validation.py --outdir /analysis/hg003 --platforms wgs,wes

bin/run --input /analysis/hg003/wgs.samples.csv \
  --resources /analysis/hg003/resources.json --regions /analysis/hg003/evaluation.bed \
  --outdir /analysis/wgs-run --profile apptainer --cpus 8 --memory-gb 96

bin/run --input /analysis/hg003/wes.autosomal.samples.csv \
  --resources /analysis/hg003/resources.json --regions /analysis/hg003/evaluation.bed \
  --outdir /analysis/wes-run --profile apptainer --cpus 8 --memory-gb 96
```

Use each output `calls.csv` with [annotation](annotation.md). Prepare a
[benchmark CSV](benchmarking.md) containing the complete GIAB truth VCF and
confidence BED. For WGS use the evaluation BED; for WES use the run's exact
`groups/wes-*/inputs/reporting.bed` as the benchmark `--regions` argument,
expanding it to one actual path. Do not benchmark WES over padded targets or
the entire chromosome. Annotation resources require the explicit preparation
in [resources](resources.md) and are stored outside Git.

CI runs validation/report-parser tests and compiles local Nextflow workflows.
Tests cover malformed inputs, reference/index checks, mixed WGS/WES separation,
interval boundaries, output locks, changed-input rejection, interruption
cleanup, and report semantics. Unit tests do not measure biological accuracy.

Pending validation is tracked in [the roadmap](roadmap.md): complete genomes
and exomes with audited training holdouts, actual GIAB SV truth, the GATK/BQSR
comparison profile, Docker and Slurm execution, and complete-dataset resource
measurements. Sex chromosomes, mitochondrial analysis, dedicated CNVs, repeat
expansions, somatic analysis, and joint genotyping are outside this release.
