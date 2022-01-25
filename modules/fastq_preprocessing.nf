
process mapping {
    tag "Mapping ${dataset} to the reference genome"
    publishDir "${params.outDir}/fastq_preprocessing", mode: 'copy', overwrite: false
    label "bigmem"
    
    input:
        tuple val(dataset), path(reads), path(ref)
    output:
        tuple val(dataset), path(mapping_out)
    script:
        mapping_out = "${dataset}.sam"
        if( !params.single_end)
            """
            minimap2 -t 32 -ax sr ${ref} ${reads[0]} ${reads[1]} > ${mapping_out}
            """
        else
            """
            minimap2 -ax sr ${ref} ${reads[0]} > ${mapping_out}  
            """
}


process sam_bam {
    tag "Converting sam to bam and sorting the bam file"
    publishDir "${params.outDir}/fastq_preprocessing", mode: 'copy', overwrite: false
    
    input:
        tuple val(dataset), path(sam_file)
    output:
        tuple val(dataset), path(sorted_bam)
    script:
        bam_file = "${dataset}.bam"
        sorted_bam = "${dataset}.sorted.bam"
        """
        samtools view -@ 20 -Sb ${sam_file} > ${bam_file}
        samtools sort -@ 20 -O bam -o ${sorted_bam} ${bam_file}
        """
}
