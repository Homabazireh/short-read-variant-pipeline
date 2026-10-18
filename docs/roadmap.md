# Implementation and validation roadmap

## Delivered in the initial implementation

- Separate private GitHub repositories for short and long reads, with documented
  contracts, reproducible software versions, container profiles, and CI.
- Human germline GRCh38 autosomal preprocessing, QC, small-variant/SV calling,
  and platform-specific sample/model validation.
- Offline VEP annotation with explicit ClinVar/gnomAD resource provenance.
- hap.py/RTG and Truvari benchmarking with confidence regions, interval scope,
  native outputs, report parsing, and safe restart controls.
- Actual public HG003 regional integration runs and compact measured reports.
- Independent callset interfaces for comparison with exported MegSAP/DRAGEN VCFs.

## Next validation milestones

1. Select complete public WGS/WES/HiFi/ONT datasets and audit each caller model's
   training overlap. Use consistent samples and assessment regions for direct
   comparisons; report WES separately within capture targets.
2. Run whole-genome/whole-exome calls, retain reference/model/resource hashes,
   and report precision, recall, F1, genotype errors, no-calls, unassessed calls,
   and difficult-region strata using the truth release's documentation.
3. Run real GIAB SV truth with release-specific inclusion regions and explicitly
   justified Truvari thresholds. Separate unsupported SV types from measured
   insertions/deletions and retain native match outputs.
4. Execute modern ONT model configurations, the GATK comparison profile in the
   short-read companion, Docker, and Slurm. Record differences from the tested
   Apptainer CPU environment.
5. Measure runtime, peak memory, disk, and cost on complete datasets. Publish
   exact hardware, concurrency, cache state, and interval scope with results.
6. If licensed DRAGEN output is supplied, benchmark it alongside these callers
   using the same sample, reference, truth, and assessment policy.

These milestones are pending; the regional measurements do not fulfill them.
New assay classes or clinical deployment need their own scoped implementation
and analytical validation.
