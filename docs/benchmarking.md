# Reproducible variant benchmarking

Benchmarking accepts pipeline outputs and independent DRAGEN/MegSAP callsets directly. Annotation is not required. `examples/benchmark.example.csv` documents the import format:

```csv
sample,label,vcf,variant_type,truth_vcf,truth_bed,evaluation_regions,truth_sample,truth_release,truth_release_notes,stratifications
HG002,pipeline,/data/query.vcf.gz,small,/truth/truth.vcf.gz,/truth/confident.bed,/data/evaluated.bed,HG002,GIAB_v4.2.1,https://github.com/genome-in-a-bottle/giab_latest_release,
```

The first six columns are required. Optional `truth_sample` defaults to `sample`; an explicit value permits known naming aliases without rewriting VCF headers. `truth_release` and `truth_release_notes` record which benchmark definition applies. Use `variant_type=sv` only with a **GRCh38-compatible SV truth release**; do not pass GRCh37 GIAB SV v0.6 files with a GRCh38 reference or casually lift them over.

```bash
bin/benchmark --input benchmark.csv --resources resources.json --outdir results/benchmark --profile apptainer --cpus 8 --memory-gb 32
```

The common options also support Docker, Slurm, `--resume` and read-only `--dry-run`. This release requires tested Nextflow **25.10.4**. Fresh outputs must be empty, concurrent writers are locked out, and resume checks both scientific input fingerprints and workflow/script hashes. Previous execution traces and reports are archived under `execution-history/` before a resumed run. Benchmarking uses only the top-level reference fields from `resources.json`; no VEP cache or annotation databases are required. Reference checksum, every input REF allele, sample identity, contig lengths, coordinates and BED compatibility are checked before launching tools. Inputs are single-sample sorted VCFs, including valid empty callsets. Checksums of truth VCF, truth BED and evaluation/stratification regions are recorded.

## Define the assessed territory

For partial chromosome tests and WES, provide `evaluation_regions` per CSV row or `--regions /path/to/evaluated.bed` for all rows. **WES must use unpadded capture targets; interval tests must use the actual called intervals.** The launcher intersects these with the benchmark confident regions, merges overlaps, and records the requested input BED base count. When omitted, the whole confident BED is supplied. It never trims the truth VCF to the evaluation BED, preserving flanking context for haplotype comparison.

The hap.py command explicitly retains its standard `--adjust-conf-regions` behavior, which adjusts insertion boundaries so already-included truth insertions can be assessed. Consequently the native confidence size can differ from the requested unpadded BED, including for WES. Reports distinguish `requested_input_bed_bases` from `native_adjusted_confident_bases` (read from hap.py's native extended aggregate). The legacy JSON field `evaluated_bases` means **requested input BED bases only**. Calls outside the tool's resulting assessed territory are unassessed, not false positives. This adjustment does not pre-trim or discard the truth VCF.

Optional `stratifications` is a tab-separated file containing `name<TAB>bed_path`, with paths relative to that file. Each stratum is intersected with the same assessed territory and benchmarked separately. This supports difficult regions such as tandem repeats and segmental duplications. Empty intersections fail clearly. Overlapping strata must not be summed. Small-variant size/subtype detail is retained in hap.py's native extended CSV; the top-level report summarizes SNP/INDEL categories, not a pooled accuracy score. SV reports additionally include INS/DEL size bins of 50–99, 100–499, 500–999, 1,000–9,999 and ≥10,000 bp. These partition completed native matches using each representation’s own type and size; matches crossing a bin boundary retain their matched status. Truth and query denominators remain separate, and size-bin rows must not be added to the aggregate row.

## Algorithms and matching policy

Small variants use **hap.py 0.3.12 with RTG vcfeval**. Query/truth records retain their original representations before hap.py's own preprocessing and haplotype comparison. Native summary, extended, ROC, metrics and annotated comparison outputs are retained. The compact report keeps PASS and ALL rows separately, reports genotype errors when available, and uses `TRUTH.TP/(TRUTH.TP+TRUTH.FN)` for recall and `QUERY.TP/(QUERY.TP+QUERY.FP)` for precision. Those TP counts can differ because representations differ.

SVs use **Truvari 5.3.0** for biallelic, sequence-resolved insertions/deletions of at least 50 bp. Symbolic variants, breakends, CNVs and complex alleles are excluded; `eligibility.tsv` shows input and eligible record counts. This exclusion is explicit and does not mean unsupported events are correct. The pinned policy is `refdist=500`, `pctseq=0.7`, `pctsize=0.7`, `pctovl=0`, `sizemin=50`, `sizefilt=50`, `sizemax=-1`, `pick=single`, PASS calls, no homozygous-reference calls, and no confident-region extension. Type matching remains enabled. SV runs require `--sv-policy policy.json` or an `sv_policy` CSV column. Copy `workflows/benchmark/sv-policy.example.json`, record the exact truth release, primary instructions URL and rationale, and pin the reference/truth VCF/truth BED checksums. The policy explicitly specifies every matching parameter. Similarity and distance thresholds can be adapted to the release instructions; unsupported scopes fail instead of silently changing the assessed events. The policy and its checksum are retained in provenance. Review the truth release instructions before choosing these parameters; the template is not automatically appropriate for every release.

Truvari measures allele matches; genotype disagreement among allele matches is reported separately. Its truth and query TPs retain their respective denominators. Truvari's native summary does not supply an unassessed-call count, so that value is null rather than guessed. Unsupported SV records are never merged into small-variant denominators.

## Reports and interpretation

`metrics.csv`, `metrics.json` and `report.html` contain actual parsed TP/FN/FP counts, precision, recall, F1 and available genotype errors. Undefined ratios remain JSON null / CSV blank / HTML NA. No calls is not treated as perfect performance. All native reports and callsets remain under `callsets/`; `provenance.json`, `trace.tsv`, `execution-report.html` and `timeline.html` preserve tool versions/digests, exact matching parameters, input hashes, runtime and resource use.

Compare equivalent samples, reference bundles, regions and coverage. Track caller model names, training sample overlap, hardware, input coverage and chemistry/basecaller versions with the originating analysis run. HG002 performance is a regression benchmark when the caller's training data include HG002; independent accuracy claims require verified held-out samples. The hand-built fixtures in unit tests carry a GRCh38 label only to exercise validation and are explicitly **synthetic plumbing fixtures**, not GRCh38 biological data or measured GIAB accuracy.

Large public-data benchmarks, real caller end-to-end validation and clinical validation are distinct from these parser/preflight tests. Report only jobs that actually completed, along with the tested territory and limitations. DRAGEN execution requires a licensed installation; importing an existing callset needs no local DRAGEN installation.

## Repeat the container plumbing test

This creates a deliberately random 10 kb sequence, small-variant/SV truth and perturbed/empty queries. The GRCh38 label only exercises the launcher schema. Its expected counts are test assertions, never biological accuracy estimates.

```bash
python3 workflows/benchmark/make_synthetic_fixture.py --outdir /scratch/synthetic-fixture
bin/benchmark --input /scratch/synthetic-fixture/benchmark.csv --resources /scratch/synthetic-fixture/resources.json --outdir /scratch/synthetic-benchmark --profile apptainer --cpus 2 --memory-gb 4
```

The perturbed query must yield SNP TP=1/FN=2/FP=2, indel TP=1/FN=0/FP=0, and SV TP=2/FN=1/FP=1 with one genotype disagreement. A matched insertion crosses the 99/100 bp bin boundary to verify postmatch stratification. Empty queries retain truth false negatives and yield undefined precision. The test runs the real pinned hap.py/RTG and Truvari containers. Adding `--resume` should cache all eight tasks; altering a query, including a change preserving size/mtime, must cause resume rejection.

## Primary documentation

- [GIAB release instructions and release-specific restrictions](https://github.com/genome-in-a-bottle/giab_latest_release)
- [NIST Genome in a Bottle guidance](https://www.nist.gov/programs-projects/faqs-genome-bottle)
- [hap.py output columns and confidence-region behavior](https://github.com/Illumina/hap.py/blob/v0.3.12/doc/happy.md)
- [Truvari benchmark matching and confident regions](https://github.com/ACEnglish/truvari/wiki/bench)
- [Pinned Truvari 5.3.0 benchmark implementation](https://github.com/ACEnglish/truvari/blob/v5.3.0/truvari/bench.py)
- [DeepVariant benchmarking instructions and hap.py container](https://github.com/google/deepvariant/blob/r1.10/docs/deepvariant-quick-start.md)
