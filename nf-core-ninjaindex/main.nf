#!/usr/bin/env nextflow
/*
========================================================================================
                         nf-core/ninjaindex
========================================================================================
 nf-core/ninjaindex Analysis Pipeline.
 #### Homepage / Documentation
 https://github.com/nf-core/ninjaindex
----------------------------------------------------------------------------------------
*/


def helpMessage() {
    // TODO nf-core: Add to this help message with new command line parameters
    log.info nfcoreHeader()
    log.info"""

    Usage:

    The typical command for running the pipeline is as follows:

    nextflow run ninjaindex.nf --genomes 'path/to/*.fna' --outdir ./output -profile docker
    nextflow run main.nf --genomes 's3://path/to/*.fna' --outdir 's3://path/to/' -profile czbiohub_aws
    s3://czbiohub-microbiome/ReferenceDBs/NinjaMap/Narrow/20190911/scv2/reference_fasta/
    nextflow run /path/to/ninjaMap/Nextflow/nf-core-ninjaindex/origin.nf -resume  --genomes 's3://czbiohub-microbiome/Sunit_Jain/Synthetic_Community/ninjaMap/20190720_00_NinjaIndex/uniform10x/setup/reference_fasta/*.fna' --outdir 's3://czbiohub-microbiome/Nextflow/test1' -profile czbiohub_aws >> output


    Mandatory arguments:
      --genomes                     Path to reference genome directory (must be surrounded with quotes)
      -profile                      Configuration profile to use. Can use multiple (comma separated)
                                    Available: conda, docker, singularity, awsbatch, test and more.

    Options:
      --outdir                      The output s3 directory where the results will be saved
      --email                       Set this parameter to your e-mail address to get a summary e-mail with details of the run sent to you when the workflow exits
      -name                         Name for the pipeline run. If not specified, Nextflow will automatically generate a random mnemonic.

    AWSBatch options:
      --awsqueue                    The AWSBatch JobQueue that needs to be set when running on AWSBatch
      --awsregion                   The AWS Region for your AWS Batch job to run on
    """.stripIndent()
}

/*
 * SET UP CONFIGURATION VARIABLES
 */

// Show help emssage
if (params.help){
    helpMessage()
    exit 0
}


//
// NOTE - THIS IS NOT USED IN THIS PIPELINE, EXAMPLE ONLY
// If you want to use the above in a process, define the following:
//   input:
//   file fasta from fasta
//


// Has the run name been specified by the user?
//  this has the bonus effect of catching both -name and --name
custom_runName = params.name
if( !(workflow.runName ==~ /[a-z]+_[a-z]+/) ){
  custom_runName = workflow.runName
}


if( workflow.profile == 'awsbatch'){ // || workflow.profile == 'czbiohub_aws' ) {
  // AWSBatch sanity checking
  if (!params.awsqueue || !params.awsregion) exit 1, "Specify correct --awsqueue and --awsregion parameters on AWSBatch!"
  // Check outdir paths to be S3 buckets if running on AWSBatch
  // related: https://github.com/nextflow-io/nextflow/issues/813
  if (!params.outdir.startsWith('s3:')) exit 1, "Outdir not on S3 - specify S3 Bucket to run on AWSBatch!"
  // Prevent trace files to be stored on S3 since S3 does not support rolling files.
  if (params.tracedir.startsWith('s3:')) exit 1, "Specify a local tracedir or run without trace! S3 cannot be used for tracefiles."
}



// Stage config files
//ch_output_docs = Channel.fromPath("$baseDir/docs/output.md")

/*
 * Create a channel for input genome files
 */
Channel
  	.fromPath(params.genomes)
  	.ifEmpty { exit 1, "Cannot find matching genomes" }
  	//.println()

genome_files = Channel.fromPath(params.genomes)
//split it into three channels
genome_files.into {genomes_ch2; genomes_ch3; genomes_ch4; genomes_ch5}


// Header log info
log.info nfcoreHeader()
def summary = [:]
if(workflow.revision) summary['Pipeline Release'] = workflow.revision
summary['Run Name']         = custom_runName ?: workflow.runName
// TODO nf-core: Report custom parameters here
summary['Max Resources']    = "$params.max_memory memory, $params.max_cpus cpus, $params.max_time time per job"
if(workflow.containerEngine) summary['Container'] = "$workflow.containerEngine - $workflow.container"
summary['Output dir']       = params.outdir
summary['Launch dir']       = workflow.launchDir
summary['Working dir']      = workflow.workDir
summary['Script dir']       = workflow.projectDir
summary['User']             = workflow.userName
summary['Docs dir']         = "$baseDir/docs/output.md"
if(workflow.profile == 'awsbatch'){
   summary['AWS Region']    = params.awsregion
   summary['AWS Queue']     = params.awsqueue
}
summary['Config Profile'] = workflow.profile
if(params.config_profile_description) summary['Config Description'] = params.config_profile_description
if(params.config_profile_contact)     summary['Config Contact']     = params.config_profile_contact
if(params.config_profile_url)         summary['Config URL']         = params.config_profile_url
if(params.email) {
  summary['E-mail Address']  = params.email
}
log.info summary.collect { k,v -> "${k.padRight(18)}: $v" }.join("\n")
log.info "\033[2m----------------------------------------------------\033[0m"

// Check the hostnames against configured profiles
checkHostname()

def create_workflow_summary(summary) {
    def yaml_file = workDir.resolve('workflow_summary_mqc.yaml')
    yaml_file.text  = """
    id: 'nf-core-ninjaindex-summary'
    description: " - this information is collected when the pipeline is started."
    section_name: 'nf-core/ninjaindex Workflow Summary'
    section_href: 'https://github.com/nf-core/ninjaindex'
    plot_type: 'html'
    data: |
        <dl class=\"dl-horizontal\">
${summary.collect { k,v -> "            <dt>$k</dt><dd><samp>${v ?: '<span style=\"color:#999999;\">N/A</a>'}</samp></dd>" }.join("\n")}
        </dl>
    """.stripIndent()

   return yaml_file
}


/*
 * Parse software version numbers
 */
process get_software_versions {
    errorStrategy 'ignore'
    publishDir "${params.outdir}/pipeline_info", mode: 'copy',
    saveAs: {filename ->
        if (filename.indexOf(".csv") > 0) filename
        else null
    }

    output:
    file 'software_versions_ninja.yaml' into software_versions_yaml
    file "software_versions.csv"

    script:
    // TODO nf-core: Get all tools to print their version number here
    """
    echo $workflow.manifest.version > v_pipeline.txt
    echo $workflow.nextflow.version > v_nextflow.txt
    scrape_software_versions.py &> software_versions_ninja.yaml
    """
}



/**
	STEP 1.1
    Concatenate the reference genomes into a single fasta file.
*/

genomes_combined = Channel
    			.fromPath(params.genomes)
    			.collectFile(name: 'all_genomes.fa')

//split it into two channels
genomes_combined.into { genomes_combined1; genomes_combined2 }
/*
*
    STEP 1.1
    Run biogrinder on individual genomes to obtain fastq files with
    predetermined uniform coverage (usually 10x)
		Copy fastq files to S3 bucket
*/

/*process grinder_fastq {

  tag "$fna"
	publishDir "${params.outdir}/biogrinder_fastqs", mode:'copy'

  input:
  file fna from genomes_ch2
  output:
  */
  //set file ("tmp_*/Sync/paired_fastq/${fna.baseName}.R1.fastq.gz"), file ("tmp_*/Sync/paired_fastq/${fna.baseName}.R2.fastq.gz") into zipped_fq
  //file "${fna.baseName}.fq" into fastq_ch
/*
  script:
  """
  run_grinder.sh $fna
  """
}
*/

/*
STEP 1.2
ART generaes synthetic reads instead of grinder
generate reads in fastq with zero-sequencing errors for a paired-end read simulation
The 2nd parameter is used to change the coverage
*/

process art_fastq {

  tag "$fna"
	publishDir "${params.outdir}/art_fastqs", mode:'copy'

  input:
  file fna from genomes_ch2
  output:
  set file ("tmp_*/Sync/paired_fastq/${fna.baseName}.R1.fastq.gz"), file ("tmp_*/Sync/paired_fastq/${fna.baseName}.R2.fastq.gz") into zipped_fq

  script:
  """
  run_art.sh $fna 10
  """
}

//split it into two channels
zipped_fq.into { zipped_fq1; zipped_fq2 }

// Sort files in channels to match fastq vs. genomes
zipped_fq2
         .toSortedList { file -> file[0].name }
				 .flatten()
				 .buffer( size: 2 )
				 .set{ sorted_zipped_fq }

sorted_zipped_fq.into { sorted_zipped_fq1; sorted_zipped_fq2 }

/*
*
   STEP 2
   Generate the bowtie2 index for each filtered genome
	 and align each fastq PE file to the corresponding filtered genomes
     file filtered_fa from sorted_filtered_genome_ch
*/

process bowtie2_mapping {
    cpus 16
    tag "$fq1"
		publishDir "${params.outdir}/bowtie2_mapping", mode:'copy'
    input:
    file all_genome from genomes_combined1.collect()
		set file(fq1), file(fq2) from sorted_zipped_fq1

    output:
    file "tmp_*/Sync/bowtie2/*.name_sorted.markdup.bam" into bam_ch
    file "tmp_*/Sync/bowtie2/*.name_sorted.markdup.bam.bai" into bai_ch

    script:
    """
    run_bowtie2_origin.sh $all_genome $fq1 $fq2 &> bowtie2.log
    """
}

genomes_ch5
				.toSortedList{file -> file.name }
				.flatten()
				.set{sorted_genome_ch}

/*
STEP 3
  		Generate final Ninja Index based on the merged bam file
      NinjaIndex need bam index in the same direcroty even though the script
      doesn't require the index as an input parameter
      This step may need more memory to process
*/
process generate_Ninja_Index {

  memory { 256.GB * task.attempt }
  time { 4.hour * task.attempt }

  errorStrategy { task.exitStatus in 137..140 ? 'retry' : 'terminate' }
  maxRetries 2

  tag "$bam"
	publishDir "${params.outdir}/ninjaIndex", mode:'copy'

  input:
  file 'fasta-dir/*' from genomes_ch4.toSortedList()
  file 'bam-dir/*' from bam_ch.toSortedList()
  file 'bam-dir/*' from bai_ch.toSortedList()

  output:
  file "tmp_*/Sync/ninjaIndex/*.ninjaIndex.binmap.csv"
  file "tmp_*/Sync/ninjaIndex/*.ninjaIndex.fasta"

  script:
  """
	ninjaIndex_multiBams.sh fasta-dir bam-dir "sc"
  """
}


/*
 * STEP 6 - Output Description HTML
 */
/*
if [ -e bamfiles.list ]
then
  touch uniform.merged.bam

fi
process output_documentation {
    publishDir "${params.outdir}/pipeline_info", mode: 'copy'

    input:
    file output_docs from ch_output_docs

    output:
    file "results_description.html"

    script:
    """
    markdown_to_html.r $output_docs results_description.html
    """
}
*/


/*
 * Completion e-mail notification
 */
workflow.onComplete {

    // Set up the e-mail variables
    def subject = "[nf-core/ninjaindex] Successful: $workflow.runName"
    if(!workflow.success){
      subject = "[nf-core/ninjaindex] FAILED: $workflow.runName"
    }
    def email_fields = [:]
    email_fields['version'] = workflow.manifest.version
    email_fields['runName'] = custom_runName ?: workflow.runName
    email_fields['success'] = workflow.success
    email_fields['dateComplete'] = workflow.complete
    email_fields['duration'] = workflow.duration
    email_fields['exitStatus'] = workflow.exitStatus
    email_fields['errorMessage'] = (workflow.errorMessage ?: 'None')
    email_fields['errorReport'] = (workflow.errorReport ?: 'None')
    email_fields['commandLine'] = workflow.commandLine
    email_fields['projectDir'] = workflow.projectDir
    email_fields['summary'] = summary
    email_fields['summary']['Date Started'] = workflow.start
    email_fields['summary']['Date Completed'] = workflow.complete
    email_fields['summary']['Pipeline script file path'] = workflow.scriptFile
    email_fields['summary']['Pipeline script hash ID'] = workflow.scriptId
    if(workflow.repository) email_fields['summary']['Pipeline repository Git URL'] = workflow.repository
    if(workflow.commitId) email_fields['summary']['Pipeline repository Git Commit'] = workflow.commitId
    if(workflow.revision) email_fields['summary']['Pipeline Git branch/tag'] = workflow.revision
    if(workflow.container) email_fields['summary']['Docker image'] = workflow.container
    email_fields['summary']['Nextflow Version'] = workflow.nextflow.version
    email_fields['summary']['Nextflow Build'] = workflow.nextflow.build
    email_fields['summary']['Nextflow Compile Timestamp'] = workflow.nextflow.timestamp

    // Render the TXT template
    def engine = new groovy.text.GStringTemplateEngine()
    def tf = new File("$baseDir/assets/email_template.txt")
    def txt_template = engine.createTemplate(tf).make(email_fields)
    def email_txt = txt_template.toString()

    // Render the HTML template
    def hf = new File("$baseDir/assets/email_template.html")
    def html_template = engine.createTemplate(hf).make(email_fields)
    def email_html = html_template.toString()

    // Render the sendmail template
    def smail_fields = [ email: params.email, subject: subject, email_txt: email_txt, email_html: email_html, baseDir: "$baseDir" ]
    def sf = new File("$baseDir/assets/sendmail_template.txt")
    def sendmail_template = engine.createTemplate(sf).make(smail_fields)
    def sendmail_html = sendmail_template.toString()

    // Send the HTML e-mail
    if (params.email) {
        try {
          if( params.plaintext_email ){ throw GroovyException('Send plaintext e-mail, not HTML') }
          // Try to send HTML e-mail using sendmail
          [ 'sendmail', '-t' ].execute() << sendmail_html
          log.info "[nf-core/ninjaindex] Sent summary e-mail to $params.email (sendmail)"
        } catch (all) {
          // Catch failures and try with plaintext
          [ 'mail', '-s', subject, params.email ].execute() << email_txt
          log.info "[nf-core/ninjaindex] Sent summary e-mail to $params.email (mail)"
        }
    }

    // Write summary e-mail HTML to a file
    /*def output_d = new File( "${params.outdir}/pipeline_info/" )
    if( !output_d.exists() ) {
      output_d.mkdirs()
    }
    def output_hf = new File( output_d, "pipeline_report.html" )
    output_hf.withWriter { w -> w << email_html }
    def output_tf = new File( output_d, "pipeline_report.txt" )
    output_tf.withWriter { w -> w << email_txt }
*/
    c_reset = params.monochrome_logs ? '' : "\033[0m";
    c_purple = params.monochrome_logs ? '' : "\033[0;35m";
    c_green = params.monochrome_logs ? '' : "\033[0;32m";
    c_red = params.monochrome_logs ? '' : "\033[0;31m";

    if (workflow.stats.ignoredCountFmt > 0 && workflow.success) {
      log.info "${c_purple}Warning, pipeline completed, but with errored process(es) ${c_reset}"
      log.info "${c_red}Number of ignored errored process(es) : ${workflow.stats.ignoredCountFmt} ${c_reset}"
      log.info "${c_green}Number of successfully ran process(es) : ${workflow.stats.succeedCountFmt} ${c_reset}"
    }

    if(workflow.success){
        log.info "${c_purple}[nf-core/ninjaindex]${c_green} Pipeline completed successfully${c_reset}"
    } else {
        checkHostname()
        log.info "${c_purple}[nf-core/ninjaindex]${c_red} Pipeline completed with errors${c_reset}"
    }

}


def nfcoreHeader(){
    // Log colors ANSI codes
    c_reset = params.monochrome_logs ? '' : "\033[0m";
    c_dim = params.monochrome_logs ? '' : "\033[2m";
    c_black = params.monochrome_logs ? '' : "\033[0;30m";
    c_green = params.monochrome_logs ? '' : "\033[0;32m";
    c_yellow = params.monochrome_logs ? '' : "\033[0;33m";
    c_blue = params.monochrome_logs ? '' : "\033[0;34m";
    c_purple = params.monochrome_logs ? '' : "\033[0;35m";
    c_cyan = params.monochrome_logs ? '' : "\033[0;36m";
    c_white = params.monochrome_logs ? '' : "\033[0;37m";

    return """    ${c_dim}----------------------------------------------------${c_reset}
                                            ${c_green},--.${c_black}/${c_green},-.${c_reset}
    ${c_blue}        ___     __   __   __   ___     ${c_green}/,-._.--~\'${c_reset}
    ${c_blue}  |\\ | |__  __ /  ` /  \\ |__) |__         ${c_yellow}}  {${c_reset}
    ${c_blue}  | \\| |       \\__, \\__/ |  \\ |___     ${c_green}\\`-._,-`-,${c_reset}
                                            ${c_green}`._,._,\'${c_reset}
    ${c_purple}  nf-core/ninjaindex v${workflow.manifest.version}${c_reset}
    ${c_dim}----------------------------------------------------${c_reset}
    """.stripIndent()
}

def checkHostname(){
    def c_reset = params.monochrome_logs ? '' : "\033[0m"
    def c_white = params.monochrome_logs ? '' : "\033[0;37m"
    def c_red = params.monochrome_logs ? '' : "\033[1;91m"
    def c_yellow_bold = params.monochrome_logs ? '' : "\033[1;93m"
    if(params.hostnames){
        def hostname = "hostname".execute().text.trim()
        params.hostnames.each { prof, hnames ->
            hnames.each { hname ->
                if(hostname.contains(hname) && !workflow.profile.contains(prof)){
                    log.error "====================================================\n" +
                            "  ${c_red}WARNING!${c_reset} You are running with `-profile $workflow.profile`\n" +
                            "  but your machine hostname is ${c_white}'$hostname'${c_reset}\n" +
                            "  ${c_yellow_bold}It's highly recommended that you use `-profile $prof${c_reset}`\n" +
                            "============================================================"
                }
            }
        }
    }
}
