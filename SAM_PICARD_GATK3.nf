// il canale per il file del reference genome per indexing
canale_genomeFa_indexing=Channel.fromPath(params.genomeFa)

// il canale per il file del reference genome per dictionary
canale_genomeFa_dictionary=Channel.fromPath(params.genomeFa)

// 


process indexing{
    input:
    path genomeFa_indexing from canale_genomeFa_indexing

    """
    samtools faidx ${genomeFa_indexing}
    """
}

process dictionary{
    input:
    path genomeFa_dictionary from canale_genomeFa_dictionary

    """
    picard CreateSequenceDictionary R=${genomeFa_dictionary}
    """
}
//
process GATK_best_1{
    input:
    path file_bam from out_secondPass

    output:
    path .arrg_s1.bam in GATK_out_1

    """
    picard AddOrReplaceReadGroups I="${runDir}{}Aligned.sortedByCoord.out.bam" O="${picardDir}{}.arrg_s1.bam" SO=coordinate RGID=0 RGPL=ILLUMINA RGLB=lib1 RGPU=group RGSM=cell
    """
}

process GATK_best_2{
    input:
    path .arrg_s1.bam from GATK_out_1

    output:
    path .md_s2.bam in GATK_out_2
    path output_metrics_{}.txt

    """
    picard MarkDuplicates I="${picardDir}{}.arrg_s1.bam" O="${picardDir}{}.md_s2.bam" CREATE_INDEX=true VALIDATION_STRINGENCY=SILENT M="${picardDir}output_metrics_{}.txt"
    """
}
process GATK_best_3{
    input:
    path .md_s2.bam from GATK_out_2

    output:

    """
    samtools index -b "${picardDir}{}.md_s2.bam"
    """
}
process GATK_best_4{
    input:
    path .md_s2.bam from GATK_out_2
    path genomeFa

    output:
    path .snc_s3.bam in GATK_out_3

    """
    gatk3 -Xmx32G -T SplitNCigarReads -R ${genomeFa} -I "${picardDir}{}.md_s2.bam" -o "${picardDir}{}.snc_s3.bam" -rf ReassignOneMappingQuality -RMQF 255 -RMQT 60 -U ALLOW_N_CIGAR_READS
    """
}
process GATK_best_5{
    input:
    path .snc_s3.bam from GATK_out_3
    path genomeFa

    output:
    path .table in GATK_out_4

    """
    gatk3 -Xmx32G -T BaseRecalibrator -R ${genomeFa} -I "${picardDir}{}.snc_s3.bam" -OQ -o "${picardDir}{}.table" -knownSites ${dbSNP_vcf}
    """
}
process GATK_best_6{
    input:
    path .table from GATK_out_4

    output:
    path .bbq_s4.bam in GATK_out_5

    """
    gatk3 -Xmx3G -T PrintReads -kpr -R ${genomeFa} -I "${picardDir}{}.snc_s3.bam" -OQ -BQSR "${picardDir}{}.table" -o "${picardDir}{}.bbq_s4.bam"
    """
}

process haplotypecaller{
    input:
    path .bbq_s4.bam from GATK_out_5
    path genomeFa
    path dbSNP_vcf

    output:
    path .raw.snps.indels.vcf in haplotypecaller_out

    """
    gatk3 -T HaplotypeCaller -R ${genomeFa} -I "${picardDir}{}.bbq_s4.bam" -o "${vcfDir}{}.raw.snps.indels.vcf" -dontUseSoftClippedBases -stand_call_conf 20 --dbsnp ${dbSNP_vcf} 
    """
}

process filtration{
    input:
    path genomeFa
    path .raw.snps.indels.vcf from haplotypecaller_out 

    output:
    path .filtered.vcf in filtration_out

    """
    gatk3 -T VariantFiltration -R ${genomeFa} -V "${vcfDir}{}.raw.snps.indels.vcf" -window 35 -cluster 3 -filterName "FS" -filter \"FS '>' 30.0\" -filterName "QD" -filter \"QD '<' 2.0\" -o "${vcfDir}{}.filtered.vcf"
    """
}


