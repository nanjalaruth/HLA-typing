process extract_hla_reads {
    tag "Extracting reads aligned to the HLA loci"
    publishDir "${params.outDir}/hla_reads", mode: 'copy', overwrite: false
    
    input:
        tuple val(dataset), path(sorted_bam)
    output:
        tuple val(dataset), path(partial_bam)
    script:
        partial_bam = "${dataset}_partial.bam"
        """
        samtools index ${sorted_bam}
        samtools view -b ${sorted_bam} "6:29000000-34000000" > ${partial_bam}
        """
}

process index_hla {
    tag "Building read name index"
    publishDir "${params.outDir}/typing", mode: 'copy', overwrite: false
    
    input:
        tuple val(dataset), path(sorted_bam)
    script:
        """
        java -jar -Xmx32g -Xms32g /scratch3/users/nanje/HLA-VBSEQ/bamNameIndex.jar index ${sorted_bam}
        """
}

process read_pairs_search {
    tag "Convert bam to fastq"
    publishDir "${params.outDir}/typing", mode: 'copy', overwrite: false
    
    input:
        tuple val(dataset), path(partial_bam)
    output:
        tuple val(dataset), file(fastq_1), file(fastq_2)
    script:
        fastq_1 = "${dataset}_partial_1.fastq"
        fastq_2 = "${dataset}_partial_2.fastq"
        """
        samtools fastq -1 ${fastq_1} -2 ${fastq_2} ${partial_bam}
        """
}

process unmapped_reads {
    tag "Extract unmapped reads"
    publishDir "${params.outDir}/typing", mode: 'copy', overwrite: false
    
    input:
        tuple val(dataset), path(sorted_bam)
    output:
        tuple val(dataset), path(fastq_1), path(fastq_2)
    script:
        fastq_1 = "${dataset}_unmapped_1.fastq"
        fastq_2 = "${dataset}_unmapped_2.fastq"
        """
        samtools view -bh -f 12 ${sorted_bam} > ${dataset}.sorted_unmapped.bam
        java -jar /scratch3/users/nanje/HLA-VBSEQ/picard-tools-1.119/SamToFastq.jar I=${dataset}.sorted_unmapped.bam F=${fastq_1} F2=${fastq_2}
        """
}

process combine_reads {
    tag "Combining reads in FASTQ format"
    publishDir "${params.outDir}/typing", mode: 'copy', overwrite: false
    
    input:
        tuple val(dataset), path(partial_fastq_1), path(partial_fastq_2), path(unmapped_fastq_1), path(unmapped_fastq_2)
    output:
        tuple val(dataset), path(fastq_1), path(fastq_2)
    script:
        fastq_1 = "${dataset}_part_1.fastq"
        fastq_2 = "${dataset}_part_2.fastq"
        """
        cat ${partial_fastq_1} ${unmapped_fastq_1} > ${fastq_1}
        cat ${partial_fastq_2} ${unmapped_fastq_2} > ${fastq_2}
        """
}

process map_to_hla_loci {
    tag "Searching read pairs and their sequences on HLA loci"
    publishDir "${params.outDir}/typing", mode: 'copy', overwrite: false
    
    input:
        tuple val(dataset), path(fastq_1), path(fastq_2), path(ref)
    output:
        tuple val(dataset), path(sam)
    script:
        sam = "${dataset}_part.sam"
        """
        #mapping using -p to map paired end as single-end
        bwa index ${ref}
        bwa mem -t 8 -P -L 10000 -a ${ref} ${fastq_1} ${fastq_2} -p > ${sam}
        """
}

process estimate_hla_types {
    tag "Estimation of HLA types by HLA-VBSeq"
    publishDir "${params.outDir}/typing", mode: 'copy', overwrite: false
    
    input:
        tuple val(dataset), path(part_sam), path(ref)
    output:
        tuple val(dataset), file(hla_txt)
    script:
        hla_txt = "${dataset}_result.txt"
        """
        #alpha_zero is a hyperparameter
        #For paired-end read data:
        #java -jar HLAVBSeq.jar hla_all_v2.fasta NA12878_part.sam NA12878_result.txt --alpha_zero 0.01 --is_paired

        #For single-end read data:
        java -jar /scratch3/users/nanje/HLA-VBSEQ/HLAVBSeq.jar ${ref} ${part_sam} ${hla_txt} --alpha_zero 0.01
        """
}



