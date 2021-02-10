#!/usr/bin/env nextflow

nextflow.enable.dsl=2

// TODO: move params to nextflow.config
params.outdir = "results"
params.reads = "reads/*R{1,2}*.fastq.gz"
params.help = false
params.save_fastp_trimmed = false
params.skip_fastp = false
params.ref_fasta = "$baseDir/data/MN908947.3.fasta"
params.primers_table = null
params.ivar_trim = false
params.run_flex = true


def helpMessage() {
    log.info"""
    ==================================================================
    ${workflow.manifest.name}  ~  version ${workflow.manifest.version}
    ==================================================================

    Git info: $workflow.repository - $workflow.revision [$workflow.commitId]

    Usage:
    The typical command for running the pipeline is as follows:
    nextflow run  --outdir 
    Options:
      --reads                          Help message for "reads" (default: "reads/*R{1,2}*.fastq.gz")
    Other options:
      --outdir                      The output directory where the results will be saved (default: )
      -w/--work-dir                 The temporary directory where intermediate data will be saved
      -profile                      Configuration profile to use. [standard, other_profiles] (default 'standard')
    """.stripIndent()
}

// Show help message if --help specified
if (params.help){
  helpMessage()
  exit 0
}

process FASTP {
  label "process_medium"
  publishDir "${params.outdir}/preprocess/fastp",
             saveAs: { filename ->
                          if (filename.endsWith(".json")) filename
                          else if (filename.endsWith(".fastp.html")) filename
                          else if (filename.endsWith("_fastqc.html")) "fastqc/$filename"
                          else if (filename.endsWith(".zip")) "fastqc/zips/$filename"
                          else if (filename.endsWith(".log")) "log/$filename"
                          else params.save_fastp_trimmed ? filename : null
                    },
             mode: 'copy'

  input:
  tuple val(sample), path(reads)

  output:
  tuple val(sample), path(reads), emit: reads
  path "*.{log,fastp,html}"

  script:
  """
  fastp \\
    --in1 ${reads[0]} \\
    --in2 ${reads[1]} \\
    --out1 ${sample}_1.trim.fastq.gz \\
    --out2 ${sample}_2.trim.fastq.gz \\
    --detect_adapter_for_pe \\
    --cut_front \\
    --cut_tail \\
    --cut_mean_quality 30 \\
    --qualified_quality_phred 30 \\
    --unqualified_percent_limit 10 \\
    --length_required 50 \\
    --trim_poly_x \\
    --thread ${task.cpus} \\
    --json ${sample}.fastp.json \\
    --html ${sample}.fastp.html
  cp .command.log ${sample}.fastp.log
  """
}

process BWA_MEM2_INDEX {
  input:
  path(fasta)

  output:
  tuple path(fasta), path("${fasta}.*")

  script:
  """
  bwa-mem2 index $fasta
  """
}

process BWA_MEM2_MAP {
  label "process_medium"
  publishDir "${params.outdir}/mapping/untrimmed",
             pattern: "*.sorted.{bam,bam.bai}",
             mode: 'copy'

  input:
  tuple val(sample), path(reads), path(ref_fasta), path(ref_index)

  output:
  tuple val(sample), path("*.sorted.{bam,bam.bai}")

  script:
  """
  bwa-mem2 mem -t ${task.cpus} $ref_fasta $reads \\
  | samtools sort -@ ${task.cpus} -T $sample \\
  > ${sample}.sorted.bam
  samtools index ${sample}.sorted.bam
  """
}

process FGBIO_TRIM_PRIMERS {
  label "process_low"
  publishDir "${params.outdir}/mapping/trimmed",
             pattern: "*.trim.{bam,bam.bai}",
             mode: 'copy'

  input:
  tuple val(sample), path(bam), path(primers_table)

  output:
  tuple val(sample), path("*.trim.{bam,bam.bai}")

  script:
  """
  fgbio \\
    --sam-validation-stringency=LENIENT \\
    TrimPrimers \\
      -i ${bam[0]} \\
      -o ${sample}.sorted.fgbio.trim.bam \\
      -p $primers_table \\
      -H true
  samtools index ${sample}.sorted.fgbio.trim.bam
  """
}

process CONVERT_PRIMER_TAB_TO_BED {
  publishDir "${params.outdir}/primers", mode: 'copy'

  input:
  path(primer_info_tab)

  output:
  path(bed)

  script:
  bed = "${file(primer_info_tab).getName()}.bed"
  """
  convert_primer_info_tab_to_bed.py $primer_info_tab $bed
  """
}

process IVAR_TRIM {
  label "process_low"
  publishDir "${params.outdir}/mapping/trimmed",
             pattern: "*.trim.{bam,bam.bai}",
             mode: 'copy'

  input:
  tuple val(sample), path(bam), path(bed)

  output:
  tuple val(sample), path("*.trim.{bam,bam.bai}")

  script:
  prefix = "${sample}.ivar.trim"
  """
  ivar trim \\
    -i ${bam[0]} \\
    -b $bed \\
    -p trim \\
    -q 20 \\
    -m 20 \\
    -s 4
  samtools sort -o ${prefix}.bam trim.bam
  samtools index ${prefix}.bam
  rm trim.bam
  """
}

process PRIMER_TRIMMED_BAM_TO_FASTQ {
  publishDir "${params.outdir}/fastq/primer_trimmed",
             pattern: "*.fastq.gz",
             mode: 'copy'
  input:
  tuple val(sample), path(bam_and_index)

  output:
  tuple val(sample), path("*.fastq.gz")

  script:
  """
  samtools fastq \\
    -1 ${sample}_1.fastq.gz \\
    -2 ${sample}_2.fastq.gz \\
    ${bam_and_index[0]}
  """
}

process CREATE_SAMPLESHEET {
  publishDir "${params.outdir}/for_viralrecon",
             mode: 'copy'
  input:
  val(samples)

  output:
  path "samplesheet.csv"

  script:
  outdir = file(params.outdir)
  sample_names = samples.join(" ")
  """
  echo "sample,fastq_1,fastq_2" > samplesheet.csv
  for sample in $sample_names; do
    r1=\$(realpath "$outdir/fastq/primer_trimmed/\${sample}_1.fastq.gz")
    r2=\$(realpath "$outdir/fastq/primer_trimmed/\${sample}_2.fastq.gz")
    echo "\$sample,\$r1,\$r2" >> samplesheet.csv
  done
  """
}


workflow {
  ch_ref = Channel.fromPath(params.ref_fasta)
  if (params.primers_table) {
    ch_primers_table = Channel.fromPath(params.primers_table)
  } else {
    if (params.run_flex) {
      ch_primers_table = Channel.fromPath("$baseDir/data/SARSCoV2.FLEX.primer_info.tab")
    } else {
      ch_primers_table = Channel.fromPath("$baseDir/data/SARSCoV2.primer_info.tab")
    }
  }

  ch_reads = Channel.fromFilePairs(
        params.reads,
        checkIfExists: true)
      .ifEmpty { exit 1, "No reads specified! Please specify where your reads are, e.g. '--reads \"/path/to/reads/*R{1,2}*.fastq.gz\"' (quotes around reads path required if using `*` and other characters expanded by the shell!)"}
      .map { [ it[0].replaceAll(/_S\d{1,2}_L001/, ""), it[1] ] }

  ch_ref | BWA_MEM2_INDEX

  if (params.skip_fastp) {
    ch_reads_ref = ch_reads | combine(BWA_MEM2_INDEX.out)
  } else {
    FASTP(ch_reads)
    ch_reads_ref = FASTP.out.reads | combine(BWA_MEM2_INDEX.out)
  }
  if (params.ivar_trim) {
    CONVERT_PRIMER_TAB_TO_BED(ch_primers_table)
    ch_reads_ref \
      | BWA_MEM2_MAP \
      | combine(CONVERT_PRIMER_TAB_TO_BED.out) \
      | IVAR_TRIM \
      | PRIMER_TRIMMED_BAM_TO_FASTQ
  } else {
    ch_reads_ref \
      | BWA_MEM2_MAP \
      | combine(ch_primers_table) \
      | FGBIO_TRIM_PRIMERS \
      | PRIMER_TRIMMED_BAM_TO_FASTQ
  }
  
  PRIMER_TRIMMED_BAM_TO_FASTQ.out | map { it[0]} | collect | CREATE_SAMPLESHEET
}
