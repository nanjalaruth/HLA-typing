#!/usr/bin/env nextflow
nextflow.enable.dsl=2

def helpMessage() {
    log.info"""
    Usage:
    The typical command for running the pipeline is as follows:
    nextflow run nf-core/hlatyping --input '*_R{1,2}.fastq.gz' -profile docker
    Mandatory arguments:
      --input [file]                  Path to input FastQ or BAM file(s). The path must be enclosed in quotes.
      -profile [str]                  Configuration profile to use. Can use multiple (comma separated).
                                      Options: conda, docker, singularity, test, awsbatch, <institute> and more
    Main options:
      --single_end [bool]             Specifies that the input is single-end reads.
                                      Default: ${params.single_end}
      --bam [bool]                    Specifies that the input is in BAM format.
                                      Default: ${params.bam}
      --seqtype [str]                 Specifies whether the input is DNA or RNA. Options: 'dna', 'rna'
                                      Default: '${params.seqtype}'
    
    Other options:
      --outdir [file]                 The output directory where the results will be saved.
      -name [str]                     Name for the pipeline run. If not specified, Nextflow will automatically generate a random mnemonic.
    """.stripIndent()
}

// Show help message
if (params.help) {
    helpMessage()
    exit 0
}
/*
 * SET UP CONFIGURATION VARIABLES
 */
// Validate inputs
params.input ?: params.input_paths ?: { log.error "No read data provided. Make sure you have used the '--input' option."; exit 1 }()
(params.seqtype == 'rna' || params.seqtype == 'dna') ?: { log.error "No or incorrect sequence type provided, you need to add '--seqtype 'dna'' or '--seqtype 'rna''."; exit 1 }()

/*
 * Create a channel for input read files
 */
if( params.input_paths ){
    if( params.single_end || params.bam) {
        Channel
            .from( params.input_paths )
            .map { row -> [ row[0], [ file( row[1][0], checkIfExists: true ) ] ] }
            .ifEmpty { exit 1, "params.input_paths or params.bams was empty - no input files supplied!" }
            .set { input_data }
    } else {
        Channel
            .from( params.input_paths )
            .map { row -> [ row[0], [ file( row[1][0], checkIfExists: true), file( row[1][1], checkIfExists: true) ] ] }
            .ifEmpty { exit 1, "params.input_paths or params.bams was empty - no input files supplied!" }
            .set { input_data }
        }
} else if (!params.bam){
    Channel
    .fromFilePairs( params.input, size: params.single_end ? 1 : 2 )
    .ifEmpty { exit 1, "Cannot find any reads matching: ${params.input}\nNB: Path needs" +
    "to be enclosed in quotes!\nNB: Path requires at least one * wildcard!\nIf this is single-end data, please specify --single_end on the command line." }
    .set { input_data }
} else {
    Channel
    .fromPath( params.input )
    .map { row -> [ file(row).baseName, [ file( row, checkIfExists: true ) ] ] }
    .ifEmpty { exit 1, "Cannot find any bam file matching: ${params.input}\nNB: Path needs" +
    "to be enclosed in quotes!\n" }
    .dump() //For debugging purposes
    .set { input_data }
}

if( params.bam ) log.info "BAM file format detected."

// Header log info
def summary = [:]
summary['File Type']        = params.bam ? 'BAM' : 'Other (fastq, fastq.gz, ...)'
summary['Seq Type']         = params.seqtype
summary['Input']            = params.input_paths ? params.input_paths : params.input
summary['Data Type']        = params.single_end ? 'Single-End' : 'Paired-End'
summary['Output Dir']       = params.outDir
summary['Launch Dir']       = workflow.launchDir
summary['Working Dir']      = workflow.workDir
summary['Script Dir']       = workflow.projectDir
if(workflow.containerEngine) summary['Container'] = "$workflow.containerEngine - $workflow.container"
if (workflow.profile.contains('awsbatch')) {
    summary['AWS Region']   = params.awsregion
    summary['AWS Queue']    = params.awsqueue
    summary['AWS CLI']      = params.awscli
}
summary['Config Profile'] = workflow.profile
if (params.config_profile_description) summary['Config Profile Description'] = params.config_profile_description
if (params.config_profile_contact)     summary['Config Profile Contact']     = params.config_profile_contact
if (params.config_profile_url)         summary['Config Profile URL']         = params.config_profile_url
summary['Config Files'] = workflow.configFiles.join(', ')
log.info summary.collect { k,v -> "${k.padRight(18)}: $v" }.join("\n")
log.info "-\033[2m--------------------------------------------------\033[0m-"


//  Modules file
include { mapping; sam_bam } from './modules/fastq_preprocessing.nf'
include { extract_hla_reads; index_hla; read_pairs_search; unmapped_reads; 
    combine_reads; map_to_hla_loci; estimate_hla_types; hla_types_out} from './modules/typing.nf'

workflow {
    
    // parameters
    // reads = Channel.from( params.input )
    reads = Channel.fromFilePairs( params.input, checkExists:true )
    ref = Channel.from(params.reference_genome)
    hla_ref = Channel.from(params.hla_ref)
    hla_txt = Channel.from(params.hla_txt_file)

    //////////////////////////////////////////////////////////////////////////////////////////////////////////////
    // fastq preprocessing

    // 1) Map fastq reads to the reference genome
    
    input_data = reads.combine(ref)
    mapping(input_data)

     // 2) Convert SAM file to BAM file and sort the BAM file
    sam_bam(mapping.out)

    ///////////////////////////////////////////////////////////////////////////////////////////////////////////////
    // HLA typing
    // 1) Extracting reads aligned to the HLA loci
    extract_hla_reads(sam_bam.out)

    // 2) Building read name index
    index_hla(sam_bam.out)

    // 3) Convert bam to fastq
    read_pairs_search(extract_hla_reads.out)

    // 4) Extract unmapped reads
    unmapped_reads(sam_bam.out)

    // // 5) Combine reads
    input = read_pairs_search.out
    input.view()
    unmap = unmapped_reads.out
    // combine_reads(input, unmap)

    // // 6) Searching read pairs and their sequences on HLA loci
    // input = combine_reads.out
    // ref = hla_ref
    // // // input.view()
    // map_to_hla_loci(input, ref)

    // // // 7) Estimate hla types
    // input = map_to_hla_loci.out.combine(hla_ref)
    // estimate_hla_types(input)
    // // estimate_hla_types.out.view()
    // hla_types_out(estimate_hla_types.out)

}
