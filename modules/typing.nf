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
        java -jar -Xmx32g -Xms32g /usr/local/bin/bamNameIndex.jar index ${sorted_bam}
        """
}

process read_pairs_search {
    tag "Convert bam to fastq"
    publishDir "${params.outDir}/typing", mode: 'copy', overwrite: false
    
    input:
        tuple val(dataset), path(partial_bam)
    output:
        tuple val(dataset), path("${dataset}_unzipped_{1,2}.fastq")
    script:
        if( !params.single_end)
            """
            samtools fastq -1 ${dataset}_unzipped_1.fastq -2 ${dataset}_unzipped_2.fastq ${partial_bam}
            """
        else
            """
            samtools fastq ${partial_bam} > ${dataset}_unzipped_1.fastq 
            """
}

process unmapped_reads {
    tag "Extract unmapped reads"
    publishDir "${params.outDir}/typing", mode: 'copy', overwrite: false
    
    input:
        tuple val(dataset), path(sorted_bam)
    output:
        tuple val(dataset), path("${dataset}_unmapped_{1,2}.fastq")
    script:
        if( !params.single_end)
            """
            samtools view -bh -f 12 ${sorted_bam} > ${dataset}.sorted_unmapped.bam
            java -jar /usr/local/bin/SamToFastq.jar I=${dataset}.sorted_unmapped.bam \
                F=${dataset}_unmapped_1.fastq F2=${dataset}_unmapped_2.fastq
            """
        else
            """
            samtools view -bh -f 12 ${sorted_bam} > ${dataset}.sorted_unmapped.bam
            java -jar /usr/local/bin/SamToFastq.jar I=${dataset}.sorted_unmapped.bam F=${dataset}_unmapped_1.fastq
            """
}
// if( !params.single_end)
process combine_reads {
    tag "Combining reads in FASTQ format"
    publishDir "${params.outDir}/typing", mode: 'copy', overwrite: false
        input:
            tuple val(dataset), path(unzipped)  
            tuple val(dataset), path(unmapped)  
        output:
            tuple val(dataset), path("${dataset}_combined_{1,2}.fastq")
        script:
            if( !params.single_end)
                """
                cat ${unzipped[0]} ${unmapped[0]} > ${dataset}_combined_1.fastq
                cat ${unzipped[1]} ${unmapped[1]} > ${dataset}_combined_2.fastq
                """
            else
                """
                cat ${unzipped[0]} ${unmapped[0]} > ${dataset}_combined_1.fastq
                """
}                


process map_to_hla_loci {
    tag "Searching read pairs and their sequences on HLA loci"
    publishDir "${params.outDir}/typing", mode: 'copy', overwrite: false
    
    input:
        tuple val(dataset), path(fastq), path(ref)
    output:
        tuple val(dataset), path(sam)
    script:
        sam = "${dataset}_part.sam"
        if( !params.single_end)
            """
            #mapping using -p to map paired end as single-end
            #index reference
            bwa index ${ref}
            #sort read pairs
            awk '{printf substr(\$0,1,length-2);getline;printf "\\t"\$0;getline;getline;print "\\t"\$0}' ${fastq[0]} | sort -S 8G -T. > read1.txt
            awk '{printf substr(\$0,1,length-2);getline;printf "\\t"\$0;getline;getline;print "\\t"\$0}' ${fastq[1]} | sort -S 8G -T. > read2.txt
            join read1.txt read2.txt | awk '{print \$1"\\n"\$2"\\n+\\n"\$3 > "r1.fq";print \$1"\\n"\$4"\\n+\\n"\$5 > "r2.fq"}'
            #mapping
            #bwa mem -t 8 -P -L 10000 -a ${ref} r1.fq r2.fq -p > ${sam}
            bwa mem -t 8 -P -L 10000 -a ${ref} r1.fq r2.fq > ${sam}
            """
        else
            """
            #mapping using -p to map paired end as single-end
            #index reference
            bwa index ${ref}

            #mapping
            bwa mem -t 8 -P -L 10000 -a ${ref} ${fastq[0]} > ${sam}
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
        if( !params.single_end)
            """
            #alpha_zero is a hyperparameter
            #For paired-end read data:
            java -jar /usr/local/bin/HLAVBSeq.jar ${ref} ${part_sam} ${hla_txt} --alpha_zero 0.01 --is_paired
            """
        else
            """
            #For single-end read data:
            java -jar /usr/local/bin/HLAVBSeq.jar ${ref} ${part_sam} ${hla_txt} --alpha_zero 0.01
            """
}

process hla_types_out {
    tag "Printing HLA types to a csv file"
    publishDir "${params.outDir}/typing", mode: 'copy', overwrite: false
    
    input:
        tuple val(dataset), path(result_txt)
    output:
        tuple val(dataset), file(hla_types)
    script:
        hla_types = "${dataset}_hlatypes.txt"
        """
        perl /usr/local/bin/parse_result.pl /scratch3/users/nanje/HLA-VBSEQ/Allelelist_v2.txt ${result_txt} | grep "^A\\*" | sort -k2 -n -r  > ${dataset}_HLA_A.txt
        perl /usr/local/bin/parse_result.pl /scratch3/users/nanje/HLA-VBSEQ/Allelelist_v2.txt ${result_txt} | grep "^B\\*" | sort -k2 -n -r  > ${dataset}_HLA_B.txt
        perl /usr/local/bin/parse_result.pl /scratch3/users/nanje/HLA-VBSEQ/Allelelist_v2.txt ${result_txt} | grep "^C\\*" | sort -k2 -n -r  > ${dataset}_HLA_C.txt
        perl /usr/local/bin/parse_result.pl /scratch3/users/nanje/HLA-VBSEQ/Allelelist_v2.txt ${result_txt} | grep "^DRB1\\*" | sort -k2 -n -r | cut -f1 > ${dataset}_HLA_DRB1.txt
        perl /usr/local/bin/parse_result.pl /scratch3/users/nanje/HLA-VBSEQ/Allelelist_v2.txt ${result_txt} | grep "^DQA1\\*" | sort -k2 -n -r | cut -f1 > ${dataset}_HLA_DQA1.txt
        perl /usr/local/bin/parse_result.pl /scratch3/users/nanje/HLA-VBSEQ/Allelelist_v2.txt ${result_txt} | grep "^DQB1\\*" | sort -k2 -n -r | cut -f1 > ${dataset}_HLA_DQB1.txt
        paste ${dataset}_HLA_A.txt ${dataset}_HLA_B.txt ${dataset}_HLA_C.txt ${dataset}_HLA_DRB1.txt ${dataset}_HLA_DQA1.txt ${dataset}_HLA_DQB1.txt > ${dataset}_hla_types
        ( echo -e "HLA_A\tHLA_B\tHLA_C\tHLA_DRB1\tHLA_DQA1\tHLA_DQB1"; cat ${dataset}_hla_types ) > ${hla_types}
        """
}

process hla_4d {
    tag "Printing 4digit HLA types"
    publishDir "${params.outDir}/typing", mode: 'copy', overwrite: false
    
    input:
        tuple val(dataset), path(result_txt), path(hla_txt)
    output:
        tuple val(dataset), file(final_report)
    script:
        final_report = "${dataset}_report.d4.txt"
        if( !params.single_end)
            """
            python /usr/local/bin/call_hla_digits.py -v ${result_txt} -a ${hla_txt} -r 90 -d 4 --ispaired > ${final_report}
            """
        else
            """
            python /usr/local/bin/call_hla_digits.py -v ${result_txt} -a ${hla_txt} -r 90 -d 4 > ${final_report}
            """
}


