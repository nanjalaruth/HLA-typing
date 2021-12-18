#!/usr/bin/env nextflow
nextflow.enable.dsl=2

include { mapping; sam_bam } from './modules/fastq_preprocessing'
include { extract_hla_reads; index_hla; read_pairs_search; unmapped_reads; 
    combine_reads; map_to_hla_loci; estimate_hla_types} from './modules/typing.nf'

workflow {
    
    // parameters
    reads = Channel.from( params.fastq_input )
    ref = Channel.from(params.reference_genome)
    hla_ref = Channel.from(params.fasta_file)
    hla_txt = Channel.from(params.text_file)

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

    // 5) Combine reads
    input = read_pairs_search.out.combine(unmapped_reads.out)
                .map {dataset, partial_1, partial_2, data, unmapped_1, unmapped_2
                -> [dataset, partial_1, partial_2, unmapped_1, unmapped_2]}
    // input.view()
    combine_reads(input)

    // 6) Searching read pairs and their sequences on HLA loci
    input = combine_reads.out.combine(hla_ref)
    // input.view()
    map_to_hla_loci(input)

    // 7) Estimate hla types
    estimate_hla_types(map_to_hla_loci.out)


}