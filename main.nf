#!/usr/bin/env nextflow
/*
========================================================================================
                         nf-core/bacass
========================================================================================
 Bacterial Genome Assembly and Analysis Pipeline.
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
                                    Available: conda, docker, singularity, awsbatch, test and more.
      --input                       The design file used for running the pipeline in TSV format.
    Pipeline arguments:
      --assembler                   Default: "Unicycler", Available: "Canu", "Miniasm", "Unicycler". Short & Hybrid assembly always runs "Unicycler".
      --assembly_type               Default: "Short", Available: "Short"
      --kraken2db                   Path to Kraken2 Database directory
      --prokka_args                 Advanced: Extra arguments to Prokka (quote and add leading space)
      --unicycler_args              Advanced: Extra arguments to Unicycler (quote and add leading space)
    Other options:
      --outdir                      The output directory where the results will be saved
      --email                       Set this parameter to your e-mail address to get a summary e-mail with details of the run sent to you when the workflow exits
      -name                         Name for the pipeline run. If not specified, Nextflow will automatically generate a random mnemonic.

   Skipping options:
      --skip_annotation             Skips the annotation with Prokka
      --skip_kraken2                Skips the read classification with Kraken2
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



// Has the run name been specified by the user?
//  this has the bonus effect of catching both -name and --name
custom_runName = params.name
if( !(workflow.runName ==~ /[a-z]+_[a-z]+/) ){
  custom_runName = workflow.runName
}

// Stage config files
ch_multiqc_config = Channel.fromPath(params.multiqc_config)
ch_output_docs = Channel.fromPath("$baseDir/docs/output.md")



process trim_and_combine {
    label 'medium'

    tag "$sample_id"
    publishDir "${params.outdir}/${sample_id}/trimming/shortreads/", mode: 'copy'

    input:
    set sample_id, file(r1), file(r2) from ch_for_short_trim

    output:
    set sample_id, file("${sample_id}_trm-cmb.R1.fastq.gz"), file("${sample_id}_trm-cmb.R2.fastq.gz") into (ch_short_for_kraken2, ch_short_for_unicycler, ch_short_for_fastqc)
    // not keeping logs for multiqc input. for that to be useful we would need to concat first and then run skewer

    script:
    """
    # loop over readunits in pairs per sample
    pairno=0
    echo "${r1} ${r2}" | xargs -n2 | while read fq1 fq2; do
	skewer --quiet -t ${task.cpus} -m pe -q 3 -n -z \$fq1 \$fq2;
    done
    cat \$(ls *trimmed-pair1.fastq.gz | sort) >> ${sample_id}_trm-cmb.R1.fastq.gz
    cat \$(ls *trimmed-pair2.fastq.gz | sort) >> ${sample_id}_trm-cmb.R2.fastq.gz
    """
}

/*
 * STEP 1 - FastQC FOR SHORT READS
*/
process fastqc {
    label 'small'
    tag "$sample_id"
    publishDir "${params.outdir}/${sample_id}/FastQC", mode: 'copy'

    input:
    set sample_id, file(fq1), file(fq2) from ch_short_for_fastqc

    output:
    file "*_fastqc.{zip,html}" into ch_fastqc_results

    script:
    """
    fastqc -t ${task.cpus} -q ${fq1} ${fq2}
    """
}


/* unicycler (short, long or hybrid mode!)
 */
process unicycler {
    tag "$sample_id"
    publishDir "${params.outdir}/${sample_id}/unicycler", mode: 'copy'

    when: params.assembler == 'unicycler'

    input:
    set sample_id, file(fq1), file(fq2), file(lrfastq) from ch_short_long_joint_unicycler

    output:
    set sample_id, file("${sample_id}_assembly.fasta") into (quast_ch, prokka_ch, dfast_ch)
    set sample_id, file("${sample_id}_assembly.gfa") into bandage_ch
    file("${sample_id}_assembly.fasta") into (ch_assembly_nanopolish_unicycler,ch_assembly_medaka_unicycler)
    file("${sample_id}_assembly.gfa")
    file("${sample_id}_assembly.png")
    file("${sample_id}_unicycler.log")

    script:
    if(params.assembly_type == 'long'){
        data_param = "-l $lrfastq"
    } else if (params.assembly_type == 'short'){
        data_param = "-1 $fq1 -2 $fq2"
    } else if (params.assembly_type == 'hybrid'){
        data_param = "-1 $fq1 -2 $fq2 -l $lrfastq"
    }

    """
    unicycler $data_param --threads ${task.cpus} ${params.unicycler_args} --keep 0 -o .
    mv unicycler.log ${sample_id}_unicycler.log
    # rename so that quast can use the name
    mv assembly.gfa ${sample_id}_assembly.gfa
    mv assembly.fasta ${sample_id}_assembly.fasta
    Bandage image ${sample_id}_assembly.gfa ${sample_id}_assembly.png
    """
}


/* kraken classification: QC for sample purity, only short end reads for now
 */
process kraken2 {
    label 'large'
    tag "$sample_id"
    publishDir "${params.outdir}/${sample_id}/kraken", mode: 'copy'

    input:
    set sample_id, file(fq1), file(fq2) from ch_short_for_kraken2

    output:
    file("${sample_id}_kraken2.report")

    when: !params.skip_kraken2

    script:
	"""
    # stdout reports per read which is not needed. kraken.report can be used with pavian
    # braken would be nice but requires readlength and correspondingly build db
	kraken2 --threads ${task.cpus} --paired --db ${kraken2db} \
		--report ${sample_id}_kraken2.report ${fq1} ${fq2} | gzip > kraken2.out.gz
	"""
}


/* assembly qc with quast
 */
process quast {
  tag {"$sample_id"}
  publishDir "${params.outdir}/${sample_id}/QUAST", mode: 'copy'

  input:
  set sample_id, file(fasta) from quast_ch

  output:
  // multiqc only detects a file called report.tsv. to avoid
  // name clash with other samples we need a directory named by sample
  file("${sample_id}_assembly_QC/")
  file("${sample_id}_assembly_QC/report.tsv") into quast_logs_ch
  file("v_quast.txt") into ch_quast_version

  script:
  """
  quast -t ${task.cpus} -o ${sample_id}_assembly_QC ${fasta}
  quast -v > v_quast.txt
  """
}



/*
 * Annotation with prokka
 */
process prokka {
   label 'large'
   tag "$sample_id"
   publishDir "${params.outdir}/${sample_id}/", mode: 'copy'

   input:
   set sample_id, file(fasta) from prokka_ch

   output:
   file("${sample_id}_annotation/")
   // multiqc prokka module is just a stub using txt. see https://github.com/ewels/MultiQC/issues/587
   // also, this only makes sense if we could set genus/species/strain. otherwise all samples
   // are the same
   // file("${sample_id}_annotation/*txt") into prokka_logs_ch

   when: !params.skip_annotation && params.annotation_tool == 'prokka'

   script:
   """
   prokka --cpus ${task.cpus} --prefix "${sample_id}" --outdir ${sample_id}_annotation ${params.prokka_args} ${fasta}
   """
}

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
   dfast --genome ${fasta} --config $config
   dfast &> v_dfast.txt 2>&1 || true
   """
}



/*
 * Parse software version numbers
 */
process get_software_versions {
    publishDir "${params.outdir}/pipeline_info", mode: 'copy',
    saveAs: {filename ->
        if (filename.indexOf(".csv") > 0) filename
        else null
    }

    input:
    file quast_version from ch_quast_version
    file porechop_version from ch_porechop_version
    file dfast_version from ch_dfast_version_for_multiqc


    output:
    file 'software_versions_mqc.yaml' into software_versions_yaml
    file "software_versions.csv"

    script:
    """
    echo $workflow.manifest.version > v_pipeline.txt
    echo $workflow.nextflow.version > v_nextflow.txt
    fastqc --version > v_fastqc.txt
    multiqc --version > v_multiqc.txt
    prokka -v 2> v_prokka.txt
    skewer -v > v_skewer.txt
    kraken2 -v > v_kraken2.txt
    Bandage -v > v_bandage.txt
    nanopolish --version > v_nanopolish.txt
    miniasm -V > v_miniasm.txt
    racon --version > v_racon.txt
    samtools --version &> v_samtools.txt 2>&1 || true
    minimap2 --version &> v_minimap2.txt
    NanoPlot --version > v_nanoplot.txt
    canu --version > v_canu.txt
    scrape_software_versions.py > software_versions_mqc.yaml
    """
}

/*
 * STEP - MultiQC
 */

process multiqc {
    label 'small'
    publishDir "${params.outdir}/MultiQC", mode: 'copy'

    input:
    file multiqc_config from ch_multiqc_config
    //file prokka_logs from prokka_logs_ch.collect().ifEmpty([])
    file ('quast_logs/*') from quast_logs_ch.collect().ifEmpty([])
    // NOTE unicycler and kraken not supported
    file ('fastqc/*') from ch_fastqc_results.collect().ifEmpty([])
    file ('software_versions/*') from software_versions_yaml.collect()
    file workflow_summary from create_workflow_summary(summary)

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


