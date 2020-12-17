# nf-cleanplex-preprocess.

Nextflow workflow to preprocess and trim primers from Paragon Genomics SARS-CoV-2 CleanPlex Illumina sequence data.

- Adapter and quality trimming with [fastp](https://github.com/OpenGene/fastp)
- Read mapping with [BWA-MEM2](https://github.com/bwa-mem2/bwa-mem2)
- Primer trimming with [fgbio](https://github.com/fulcrumgenomics/fgbio)
- Output to FASTQ for processing with [viralrecon](https://github.com/nf-core/viralrecon)
