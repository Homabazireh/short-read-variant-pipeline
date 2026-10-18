# Third-party software and data

The MIT license covers original orchestration, tests and documentation in this
repository. External programs, containers, models, data, annotation resources
and incorporated upstream modules retain their respective licenses and terms.
They are obtained from their publishers and are not relicensed by this project.

This project is inspired by [MegSAP](https://github.com/imgag/megSAP) and the
end-to-end capabilities of Illumina DRAGEN. MegSAP source code is not copied.
DRAGEN is not bundled and its execution requires the user's licensed installation.
Exported DRAGEN and MegSAP VCFs can be annotated and benchmarked through the
documented CSV interfaces. DRAGMAP is a mapper, not the complete DRAGEN platform.

Core upstream projects include [Nextflow](https://github.com/nextflow-io/nextflow),
[nf-core/sarek](https://github.com/nf-core/sarek),
[DeepVariant](https://github.com/google/deepvariant),
[Clair3](https://github.com/HKU-BAL/Clair3),
[Ensembl VEP](https://github.com/Ensembl/ensembl-vep),
[hap.py](https://github.com/Illumina/hap.py), and
[Truvari](https://github.com/ACEnglish/truvari).

Benchmark inputs are public Genome in a Bottle materials. Consult
[NIST's benchmark guidance](https://www.nist.gov/programs-projects/faqs-genome-bottle)
and the individual release README before interpreting measurements. Report
pretrained-model training overlap and distinguish bounded integration checks
from independent whole-genome performance studies.
