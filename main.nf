#!/usr/bin/env nextflow
/*
========================================================================================
                  Bacterial genome assembly, pan-genome and mGWAS pipeline
========================================================================================
 Extended from the nf-core bacass pipeline.
 #### Homepage / Documentation
 To be done
----------------------------------------------------------------------------------------
*/


def helpMessage() {
    log.info nfcoreHeader()
    log.info"""
    Usage:
    The typical command for running the pipeline is as follows:
    nextflow run nf-core/bacass --input input.csv --kraken2db 'path-to-kraken2db' -profile docker
    Mandatory arguments:
      -profile                      Configuration profile to use. Can use multiple (comma separated)
                                    Available: ilifu.
      --reads                       The sample sheet containing the paths to the fastq files, as well as sample names.
      --genome                      The reference genome to be used in fasta format. Also acts as an outgroup.
      --gff                         Path to GFF3 file OR (see next arg)
      --gtf                         Path to GTF file
    Pipeline arguments:
      --prokka_args                 Advanced: Extra arguments to Prokka (quote and add leading space)
      --unicycler_args              Advanced: Extra arguments to Unicycler (quote and add leading space)
      --phenotype_info              Dichotomized csv file containing phenotype information for the samples. (Optional, needs to be added below)
    Other options:
      --outdir                      The output directory where the results will be saved
      --email                       Set this parameter to your e-mail address to get a summary e-mail with details of the run sent to you when the workflow exits
      -name                         Name for the pipeline run. If not specified, Nextflow will automatically generate a random mnemonic.

    """.stripIndent()
}

/*
 * SET UP CONFIGURATION VARIABLES
 */

// Show help message
if (params.help){
    helpMessage()
    exit 0
}

// Define regular variables so that they can be overwritten
clip_r1 = params.clip_r1
clip_r2 = params.clip_r2
three_prime_clip_r1 = params.three_prime_clip_r1
three_prime_clip_r2 = params.three_prime_clip_r2
forward_stranded = params.forward_stranded
reverse_stranded = params.reverse_stranded
unstranded = params.unstranded

// Preset trimming options
if (params.pico){
    clip_r1 = 3
    clip_r2 = 0
    three_prime_clip_r1 = 0
    three_prime_clip_r2 = 3
    forward_stranded = true
    reverse_stranded = false
    unstranded = false
}


// Has the run name been specified by the user?
//  this has the bonus effect of catching both -name and --name
custom_runName = params.name
if( !(workflow.runName ==~ /[a-z]+_[a-z]+/) ){
  custom_runName = workflow.runName
}

// Stage config files
ch_multiqc_config = file(params.multiqc_config, checkIfExists: true)


//Validate inputs
if ( params.genome == false ) {
    exit 1, "Must set a reference genome fasta file (--genome)"
}

if ( params.reads == false ) {
    exit 1, "Must set the path to the sample file (--reads) in csv format"
}

// Dealing with GFT and GFF.
if( params.gtf ){
    Channel
        .fromPath(params.gtf)
        .ifEmpty { exit 1, "GTF annotation file not found: ${params.gtf}" }
        .into { gtfFile }
} else if( params.gff ){
    Channel
        .fromPath(params.gff)
        .ifEmpty { exit 1, "GFF annotation file not found: ${params.gff}" }
        .into { gffFile }
} else {
    exit 1, "No GTF or GFF3 annotation specified!"
}


if(params.phenotype_info) {
    Channel
        .fromPath(params.phenotype_info)
        .ifEmpty { exit 1, "Phenotype file not found: ${params.phenotype_info}" }
        .into { phenotype_file }
}


genome_file         = file(params.genome)
sample_sheet        = Channel.fromPath(params.reads)
sample_sheet_QC     = Channel.fromPath(params.reads)
reads_ch            = Channel.fromFilePairs(params.reads)

sample_sheet
  .splitCsv(header:true)
  .map{ row-> tuple(row.number, file(row.R1), file(row.R2)) }
  .set { newSampleChannel }


sample_sheet_QC
  .splitCsv(header:true)
  .map{ row-> tuple(row.number, file(row.R1), file(row.R2)) }
  .set { newSampleChannelFastQC }



/*
 * Parse software version numbers
 */
process get_software_versions {

    output:
      file 'software_versions_mqc.yaml' into software_versions_yaml

    script:
    """
    echo \$(bamtools --version 2>&1) > v_bamtools.txt
    samtools --version > v_samtools.txt
    echo $workflow.manifest.version > v_pipeline.txt
    echo $workflow.nextflow.version > v_nextflow.txt
    preseq &> v_preseq.txt
    multiqc --version > v_multiqc.txt
    trim_galore --version > v_trim_galore.txt
    picard MarkDuplicates --version &> v_picard.txt  || true
    echo \$(bwa 2>&1) > v_bwa.txt
    #scrape_software_versions.py > software_versions_mqc.yaml
    echo "not yet" > software_versions_mqc.yaml
    """
}



/*
 * PREPROCESSING - Convert GFF3 to GTF
 */
if(params.gff){
  process convertGFFtoGTF {
      tag "$gff"

      input:
        file gff from gffFile

      output:
        file "${gff.baseName}.gtf" into gtf_makeSTARindex, gtf_makeBED12, gtf_star, gtf_dupradar, gtf_featureCounts

      script:
      """
      gffread $gff -T -o ${gff.baseName}.gtf
      """
  }
} else {

  process convertGTFtoGFF {

  input:
    file gtf from gtfFile

  output:

    file "${gtf.baseName}.gtf" into gtf_makeSTARindex, gtf_makeBED12, gtf_star, gtf_dupradar, gtf_featureCounts
    file "${gtf.baseName}.gff" into snpeff_gff_old_to_test_deleting

  script:
  """
  gffread $gtf -o ${gtf.baseName}.gff
  """

  }

}



/**********
 * PART 1: Data preparation
 *
 * Process 1A: Create a FASTA genome index (.fai) with samtools for GATK
 */

process prepare_genome_samtools {
  tag "$genome.baseName"

  input:
      file genome from genome_file

  output:
      file "${genome}.fai" into genome_index_ch

  script:
  """
  samtools faidx ${genome}
  """
}


/*
 * Process 1B: Create a FASTA genome sequence dictionary with Picard for GATK
 */

process prepare_genome_picard {
  tag "$genome.baseName"

  input:
      file genome from genome_file
  output:
      file "${genome.baseName}.dict" into genome_dict_ch

  script:
  """
  #picard -XX:ParallelGCThreads=5 -Xmx16G -Xms16G CreateSequenceDictionary R=$genome O=${genome.baseName}.dict
  picard CreateSequenceDictionary R=$genome O=${genome.baseName}.dict
  """
}

/*
 * Process 1C: Create a FASTA genome sequence dictionary for BWA
 */

process prepare_genome_bwa {
  tag "$genome.baseName"

  input:
      file genome from genome_file
  output:
      file "${genome}.amb" into genome_bwa_amb
      file "${genome}.ann" into genome_bwa_ann
      file "${genome}.bwt" into genome_bwa_bwt
      file "${genome}.pac" into genome_bwa_pac
      file "${genome}.sa" into genome_bwa_sa

  script:
  """
  bwa index $genome
  """
}



/*
 * Process 1D: FastQC -  NEED TO EDIT
 */

 process fastqc {
    tag "$name"
    publishDir "${params.outdir}/fastqc", mode: 'copy',
        saveAs: {filename -> filename.indexOf(".zip") > 0 ? "zips/$filename" : "$filename"}

    input:
        set number, file(R1), file(R2) from newSampleChannelFastQC

    output:
        file "*_fastqc.{zip,html}" into fastqc_results
        file "*_fastqc.{zip,html}" into ch_fastqc_results

    script:
    """
    fastqc -q $R1
    fastqc -q $R2
    """
}



/**********
 * PART 2: Data filtering and trimming
 *
 * STEP 1 - Trim Galore!
 *
 * Once again, dealing with different files names and extensions has made things tricky.
 */

process trim_galore {
    label 'low_memory'
    tag "$name"
    publishDir "${params.outdir}/trim_galore", mode: 'copy',
        saveAs: {filename ->
            if (filename.indexOf("_fastqc") > 0) "FastQC/$filename"
            else if (filename.indexOf("trimming_report.txt") > 0) "logs/$filename"
            else if (!params.saveTrimmed && filename == "where_are_my_files.txt") filename
            else if (params.saveTrimmed && filename != "where_are_my_files.txt") filename
            else null
        }

    input:
    set number, file(R1), file(R2) from newSampleChannel

    output:

    file "*_1.fq.gz" into forwardTrimmed
    file "*_2.fq.gz" into reverseTrimmed

    set file("*trimming_report.txt"), file("*_fastqc.{zip,html}") into ch_trimgalore_results_mqc
    val "$number" into sampleNumber
    set number, file("*_1.fq.gz"), file("*_2.fq.gz") into vf_read_pairs
    set number, file("*_1.fq.gz"), file("*_2.fq.gz") into unicycler_read_pairs

    script:
    c_r1 = clip_r1 > 0 ? "--clip_r1 ${clip_r1}" : ''
    c_r2 = clip_r2 > 0 ? "--clip_r2 ${clip_r2}" : ''
    tpc_r1 = three_prime_clip_r1 > 0 ? "--three_prime_clip_r1 ${three_prime_clip_r1}" : ''
    tpc_r2 = three_prime_clip_r2 > 0 ? "--three_prime_clip_r2 ${three_prime_clip_r2}" : ''
    if (params.singleEnd) {
        """
        trim_galore --fastqc --gzip $c_r1 $tpc_r1 $R1 $R2
        """
    } else {
        """
        trim_galore --paired --fastqc --gzip $c_r1 $c_r2 $tpc_r1 $tpc_r2 $R1 $R2

        rename 's/val_1_001/1/' *.fq.gz
        rename 's/val_2_001/2/' *.fq.gz

        rename 's/_val_1//' *.fq.gz
        rename 's/_val_2//' *.fq.gz

        rename 's/_R1/_1/' *.fq.gz
        rename 's/_R2/_2/' *.fq.gz

        """
    }
}



/*
 * Process 2B: Align reads to the reference genome using bwa mem
 */

process read_mapping {
  label 'high_memory'
  input:
    file forwardTrimmed
    file reverseTrimmed
    val sampleNumber
    file genome from genome_file
    file genome_bwa_amb
    file genome_bwa_ann
    file genome_bwa_bwt
    file genome_bwa_pac
    file genome_bwa_sa
  output:
    file "sample_${sampleNumber}_sorted.bam" into bamfiles
    file "sample_${sampleNumber}_sorted.bai" into bamindexfiles
    file "sample_${sampleNumber}_sorted.bam" into bam_rseqc
    file "sample_${sampleNumber}_sorted.bai" into bamindexfiles_rseqc
    file "sample_${sampleNumber}_sorted.bam" into bam_preseq
    file "sample_${sampleNumber}_sorted.bam" into bam_forSubsamp
    file "sample_${sampleNumber}_sorted.bam" into bam_skipSubsamp
    file "sample_${sampleNumber}_sorted.bam" into bam_featurecounts
  script:
  """
  bwa mem $genome $forwardTrimmed $reverseTrimmed | /tools/samtools-1.10/samtools sort -O BAM -o sample_${sampleNumber}_sorted.bam
  /tools/samtools-1.10/samtools index sample_${sampleNumber}_sorted.bam sample_${sampleNumber}_sorted.bai
  """
}



/*
 * Process 2C: Mark duplicate reads with picard - EDIT for QC
 */

process mark_duplicates {
  label 'high_memory'
  tag "${sample_bam.baseName}"

  input:
    file sample_bam from bamfiles
  output:
    file "${sample_bam.baseName}_dedup.bam" into dedup_bamfiles
    file "${sample_bam.baseName}_dedup.bam" into bam_md
    file "${sample_bam.baseName}_dedup.bam.bai"
    file "${sample_bam.baseName}.txt" into picard_results
  script:
    """
    picard MarkDuplicates INPUT=$sample_bam OUTPUT=${sample_bam.baseName}_dedup.bam METRICS_FILE=${sample_bam.baseName}.txt ASSUME_SORTED=true REMOVE_DUPLICATES=false
    samtools index ${sample_bam.baseName}_dedup.bam
    """
}



/**********
 * PART 3: Assembly and annotation
 *
 * STEP 1 - Unicycler (short mode!) for genome assembly
*/

process unicycler {
    label 'high_memory'
    tag "$sample_id"
    publishDir "${params.outdir}/${sample_id}_assembly/unicycler", mode: 'copy'

    input:
        set sample_id, file(fq1), file(fq2) from unicycler_read_pairs

    output:
        set sample_id, file("${sample_id}_assembly.fasta") into (prokka_ch, quast_ch, dfast_ch)
        file("${sample_id}_assembly.fasta") into (ch_assembly_nanopolish_unicycler,ch_assembly_medaka_unicycler)
        file("${sample_id}_assembly.gfa")
        file("${sample_id}_assembly.png")
        file("${sample_id}_unicycler.log")

    script:
    """
    unicycler -1 $fq1 -2 $fq2 --threads ${task.cpus} --pilon_path /usr/local/bin/pilon-1.23.jar --keep 0 -o .
    mv unicycler.log ${sample_id}_unicycler.log
    # rename so that quast can use the name
    mv assembly.gfa ${sample_id}_assembly.gfa
    mv assembly.fasta ${sample_id}_assembly.fasta
    Bandage image ${sample_id}_assembly.gfa ${sample_id}_assembly.png
    """
}


/*
* Assembly qc with quast
*/
process quast {
  tag {"$sample_id"}
  publishDir "${params.outdir}/${sample_id}_assembly/QUAST", mode: 'copy'

  input:
    set sample_id, file(fasta) from quast_ch

  output:
      // multiqc only detects a file called report.tsv. to avoid
      // name clash with other samples we need a directory named by sample
      file("${sample_id}_assembly_QC/")
      file("${sample_id}_assembly_QC/${sample_id}_report.tsv") into quast_logs_ch
      file("v_quast_${sample_id}.txt") into ch_quast_version

  script:
  """
  quast.py -t ${task.cpus} -o ${sample_id}_assembly_QC ${fasta}
  quast.py -v > v_quast_${sample_id}.txt
  mv ${sample_id}_assembly_QC/report.tsv ${sample_id}_assembly_QC/${sample_id}_report.tsv
  """
}



/*
 * Annotation with prokka
 */
process prokka {
    label 'high_memory'
    tag "$sample_id"
    publishDir "${params.outdir}/${sample_id}_assembly/prokka", mode: 'copy'

    input:
        set sample_id, file(fasta) from prokka_ch

    output:
        file("${sample_id}_annotation/sample_${sample_id}.gff") into gff
        // multiqc prokka module is just a stub using txt. see https://github.com/ewels/MultiQC/issues/587
        // also, this only makes sense if we could set genus/species/strain. otherwise all samples
        // are the same
        // file("${sample_id}_annotation/*txt") into prokka_logs_ch

   script:
   """
   prokka --cpus ${task.cpus} --prefix sample_${sample_id} --outdir ${sample_id}_annotation ${fasta}
   """
}

/*
 * Annotation with dfast (if annotation_tool == 'dfast')
 */
process dfast {
    tag "$sample_id"
    publishDir "${params.outdir}/${sample_id}/", mode: 'copy'

    input:
        set sample_id, file(fasta) from dfast_ch
        file (config) from Channel.value(params.dfast_config ? file(params.dfast_config) : "")

    output:
        file("RESULT*")
        file("v_dfast.txt") into ch_dfast_version_for_multiqc

    when: !params.skip_annotation && params.annotation_tool == 'dfast'

    script:
    """
    python3 /usr/local/bin/dfast --genome ${fasta} --config $config
    python3 /usr/local/bin/dfast &> v_dfast.txt 2>&1 || true
    """
}


/**********
 * PART 4: pan-genome analysis and mGWAS
 *
 * STEP 1 - Roary - rapid large-scale prokaryote pan genome analysis
 */

process roary {
    publishDir "${params.outdir}/roary", mode: 'copy'

    input:
        file gff from gff.collect()

    output:
        file("*") into roary
        file("pan_genome_reference.fa") into pan_genome
        file("gene_presence_absence.csv") into mGWAS_GPA

    script:
    """
    roary -e -n -v -r $gff
    """
}


/*
* Create required mGWAS input files using the output of roary
*/

if(params.phenotype_info) {
    process post_roary {

        input:
            file sample_phenotype_info from phenotype_file

        output:
            file("traits.csv") into traits_file

        script:
        """
        # Issue caused by not being able to run binaries.
        #post_roary.py --phenotype_data $sample_phenotype_info
        echo 'this' > traits.csv
        """
    }
}

/*
* Run scoary
*/

if (phenotype_file) {
    process scoary {

    input:
        file traits_file
        file mGWAS_GPA

    output:
        file("*") into scoary

    script:
    """
    scoary -g $mGWAS_GPA -t $traits_file
    """
    }
}



/**********
 * PART 5: Assembly and annotation
 *
 * MultiQC
 */

/*
process multiqc {
    publishDir "${params.outdir}/MultiQC", mode: 'copy'

    input:
        file multiqc_config from ch_multiqc_config
        //file ('software_versions/*') from software_versions_yaml
        file ('quast_logs/*') from quast_logs_ch.collect().ifEmpty([])
        //file ('fastqc/*') from ch_fastqc_results.collect().ifEmpty([])
        //file ('trim_galore/*') from ch_trimgalore_results_mqc.collect()

    output:
        file "*multiqc_report.html" into multiqc_report
        file "*_data"
        file "multiqc_plots"

    script:
    rtitle = custom_runName ? "--title \"$custom_runName\"" : ''
    rfilename = custom_runName ? "--filename " + custom_runName.replaceAll('\\W','_').replaceAll('_+','_') + "_multiqc_report" : ''
    """
    multiqc -f $rtitle $rfilename --config $multiqc_config .
    """
}

*/
