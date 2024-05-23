#! /usr/bin/env nextflow
nextflow.enable.dsl=2

log.info """\
       ===================================
        Neoantigen nextflow pipeline
       ===================================
        sample: ${params.sample}
        reads: ${params.sample_list}
        outdir: ${params.outdir}
        """
        .stripIndent()

process FASTQC {
    publishDir "${params.outdir}/${params.sample}/01_quality_control", mode:'symlink'
    publishDir "${params.outdir}/${params.sample}/00_log", mode:'symlink',saveAs: {filename -> if (filename =~ /.*.log$/) "${sample_id}_FASTQC.log"}
    tag "${params.sample}"
    //cpus 8
    label "m_cpus"
    //echo true

    input:
      tuple val(sample_id), path(reads)

    output:
      tuple val(sample_id), path("fastqc/*")
      path(".command.{log,err,out}")

    script:
    """
      mkdir fastqc
      fastqc -o fastqc ${reads[0]} ${reads[1]} -t ${task.cpus}
    """
}

process FASTP {
    publishDir "${params.outdir}/${params.sample}/01_quality_control/fastp", mode:'symlink'
    publishDir "${params.outdir}/${params.sample}/00_log", mode:'symlink',saveAs: {filename -> if (filename =~ /.*.log$/) "${sample_id}_FASTP.log"}
    tag "${params.sample}"
    //cpus 8
    label "m_cpus"
    //echo true
 
    input:
      tuple val(sample_id), path(reads)

    output:
      tuple val(sample_id), path("*_{R1,R2}.f*q*")
      path("*{html,json}")
      path(".command.{log,err,out}")

    script:
    """
      #mkdir fastp;cd fastp
      fastp -i ${reads[0]} -I ${reads[1]} -o ${sample_id}_fastp_R1.fq.gz -O ${sample_id}_fastp_R2.fq.gz \
      -w ${task.cpus} --html ${sample_id}.html --json ${sample_id}.json
    """
}

process BWA_INDEX {
    //publishDir "${params.outdir}/${params.sample}/04_rna_expression/index", mode:'symlink'
    //tag "${params.sample}"
    echo true
    beforeScript "export SENTIEON_LICENSE=${params.sentieon_license}"

    script:
    if( file(params.dnaIndexFiles) ) {
      """
      echo 'dna index files for bwa exist'
      """
    }
    else {
      """
      echo 'dna index files for bwa building'
      dir=\$(dirname ${params.ref})
      cd \$dir
      ${params.sentieon_dir}/sentieon bwa index ${params.ref}
      samtools faidx ${params.ref}
      picard CreateSequenceDictionary -R ${params.ref}
      """
    }
}

//Map reads to reference
process SENTIEON_MAP {
    publishDir "${params.outdir}/${params.sample}/02_alignment_variant_calling_phasing", mode:'symlink'
    publishDir "${params.outdir}/${params.sample}/00_log", mode:'symlink',saveAs: {filename -> if (filename =~ /.*.log$/) "${sample_id}_SENTIEON_MAP.log"}
    tag "${params.sample}"
    //cpus 16
    label "l_cpus"
    echo true
    beforeScript "export SENTIEON_LICENSE=${params.sentieon_license}"

    input:
      tuple val(sample_id), path(reads)

    output:
      tuple val(sample_id), path("*.bam*")
      path(".command.{log,err,out}")

    when:
      !(sample_id =~ /RNA|rna|Rna/)

    script:
    """
      (${params.sentieon_dir}/sentieon bwa mem -R "@RG\\tID:${sample_id}\\tSM:${sample_id}\\tPL:${params.platform}" \
      -K 10000000 -t ${task.cpus} ${params.ref} ${reads[0]} ${reads[1]}|| echo -n 'error' ) | \
      ${params.sentieon_dir}/sentieon util sort -r ${params.ref} -o ${sample_id}_sorted.bam -t ${task.cpus} --sam2bam -i -
    """
}

//Calculate data align metrics
process SENTIEON_METRIC {
    publishDir "${params.outdir}/${params.sample}/02_alignment_variant_calling_phasing/data_align_metrics", mode:'symlink'
    publishDir "${params.outdir}/${params.sample}/00_log", mode:'symlink',saveAs: {filename -> if (filename =~ /.*.log$/) "${sample_id}_SENTIEON_METRIC.log"}
    tag "${params.sample}"
    //cpus 8
    label "m_cpus"
    echo true
    beforeScript "export SENTIEON_LICENSE=${params.sentieon_license}"

    input:
      tuple val(sample_id), path(bam)

    output:
      tuple val(sample_id), path("*")
      path(".command.{log,err,out}")

    script:
    """
      ${params.sentieon_dir}/sentieon driver -r ${params.ref} -t ${task.cpus} -i ${bam[0]} \
      --algo GCBias --summary ${sample_id}_gc_summary.txt ${sample_id}_gc_metrics.txt \
      --algo MeanQualityByCycle ${sample_id}_mq_metrics.txt \
      --algo QualDistribution ${sample_id}_qd_metrics.txt \
      --algo AlignmentStat --adapter_seq "" ${sample_id}_aln_metrics.txt \
      --algo InsertSizeMetricAlgo ${sample_id}_is_metrics.txt

      ${params.sentieon_dir}/sentieon plot GCBias -o ${sample_id}_gc_report.pdf ${sample_id}_gc_metrics.txt
      ${params.sentieon_dir}/sentieon plot QualDistribution -o ${sample_id}_qd_report.pdf ${sample_id}_qd_metrics.txt
      ${params.sentieon_dir}/sentieon plot MeanQualityByCycle -o ${sample_id}_mq_report.pdf ${sample_id}_mq_metrics.txt
      ${params.sentieon_dir}/sentieon plot InsertSizeMetricAlgo -o ${sample_id}_is_report.pdf ${sample_id}_is_metrics.txt
    """
}

//Mark duplicates (Remove add –rmdup)
process SENTIEON_DEDUP {
    publishDir "${params.outdir}/${params.sample}/02_alignment_variant_calling_phasing", pattern: "*.bam*", mode:'symlink'
    publishDir "${params.outdir}/${params.sample}/02_alignment_variant_calling_phasing/data_dedup_metrics", pattern: "*.{txt,idx}*", mode:'symlink'
    publishDir "${params.outdir}/${params.sample}/00_log", mode:'symlink',saveAs: {filename -> if (filename =~ /.*.log$/) "${sample_id}_SENTIEON_DEDUP.log"}
    tag "${params.sample}"
    //cpus 8
    label "m_cpus"
    echo true
    beforeScript "export SENTIEON_LICENSE=${params.sentieon_license}"

    input:
      tuple val(sample_id), path(bam)

    output:
      //tuple val(sample_id), path("*")
      tuple val(sample_id), path("*.bam*")
      path("*.{txt,idx}*")
      path(".command.{log,err,out}")

    script:
    """
      ${params.sentieon_dir}/sentieon driver -t ${task.cpus} -i ${bam[0]} \
      --algo LocusCollector --fun score_info ${sample_id}_score.txt

      ${params.sentieon_dir}/sentieon driver -t ${task.cpus} -i ${bam[0]} \
      --algo Dedup --score_info ${sample_id}_score.txt \
      --metrics ${sample_id}_dedup_metrics.txt ${sample_id}_sorted_deduped.bam
    """
}

//Base quality score recalibration (BQSR)
if(params.exome_bed)
    println 'WES data\n'
else
    println 'WGS data or no WES exome bed\n'
process SENTIEON_BQSR {
    publishDir "${params.outdir}/${params.sample}/02_alignment_variant_calling_phasing/data_bqsr_metrics", pattern: "*.{csv,table,pdf}*", mode:'symlink'
    publishDir "${params.outdir}/${params.sample}/00_log", mode:'symlink',saveAs: {filename -> if (filename =~ /.*.log$/) "${sample_id}_SENTIEON_BQSR.log"}
    tag "${params.sample}"
    //cpus 8
    label "m_cpus"
    echo true
    beforeScript "export SENTIEON_LICENSE=${params.sentieon_license}"

    input:
      tuple val(sample_id), path(deduped_bam)

    output:
      //tuple val(sample_id), path("*.bam*")
      path("*.{csv,post,pdf}")
      path("${sample_id}_recal_data.table"), emit: recal_data
      path(".command.{log,err,out}")

    script:
      //def exome_option = "--interval ${params.exome_bed}" ? file(params.exome_bed) : ""
      def exome_option = params.exome_bed ? "--interval ${params.exome_bed}" : " "
    """
      ${params.sentieon_dir}/sentieon driver -r ${params.ref} -t ${task.cpus} -i ${deduped_bam[0]} ${exome_option} \
      --algo QualCal \
      -k ${params.smvdb1} -k ${params.smvdb2} -k ${params.smvdb3} \
      ${sample_id}_recal_data.table

      ${params.sentieon_dir}/sentieon driver -r ${params.ref} -t ${task.cpus} -i ${deduped_bam[0]} ${exome_option} \
      -q ${sample_id}_recal_data.table \
      --algo QualCal \
      -k ${params.smvdb1} -k ${params.smvdb2} -k ${params.smvdb3} \
      ${sample_id}_recal_data.table.post

      ${params.sentieon_dir}/sentieon driver -t ${task.cpus} --algo QualCal --plot \
      --before ${sample_id}_recal_data.table --after ${sample_id}_recal_data.table.post ${sample_id}_recal.csv

      ${params.sentieon_dir}/sentieon plot QualCal -o ${sample_id}_recal_plots.pdf ${sample_id}_recal.csv
    """
}

//HC Variant caller
process SENTIEON_HAPLOTYPER {
    publishDir "${params.outdir}/${params.sample}/02_alignment_variant_calling_phasing", pattern: "*.vcf*", mode:'symlink'
    publishDir "${params.outdir}/${params.sample}/00_log", mode:'symlink',saveAs: {filename -> if (filename =~ /.*.log$/) "${sample_id}_SENTIEON_HAPLOTYPER.log"}
    tag "${params.sample}"
    //cpus 16
    label "l_cpus"
    echo true
    beforeScript "export SENTIEON_LICENSE=${params.sentieon_license}"

    input:
      tuple val(sample_id), path(deduped_bam), path(recal_data_table)
      //path "${sample_id}_recal_data.table"

    output:
      tuple val(sample_id), path("*.vcf*")
      path(".command.{log,err,out}")

    when:
      sample_id =~ /(?i)normal|-n/

    script:
      def exome_option = params.exome_bed ? "--interval ${params.exome_bed}" : " "
    """
      ${params.sentieon_dir}/sentieon driver -r ${params.ref} -t ${task.cpus} ${exome_option} \
      -i ${deduped_bam[0]} -q ${sample_id}_recal_data.table \
      --algo Haplotyper -d ${params.smvdb1} --emit_conf=30 --call_conf=30 ${sample_id}.vcf.gz
    """
      //--emit_conf=30 --call_conf=30 (default)
}

//Somatic variant discovery with matched normal sample
process SENTIEON_TNSCOPE {
    publishDir "${params.outdir}/${params.sample}/02_alignment_variant_calling_phasing", pattern: "*.vcf*", mode:'symlink'
    publishDir "${params.outdir}/${params.sample}/00_log", mode:'symlink',saveAs: {filename -> if (filename =~ /.*.log$/) "${sample_id_t}_SENTIEON_TNSCOPE.log"}
    tag "$sample_id_t"
    //cpus 16
    label "l_cpus"
    echo true
    beforeScript "export SENTIEON_LICENSE=${params.sentieon_license}"

    input:
      tuple val(sample_id_n), path(dedup_bam_n), path(recal_n), val(sample_id_t), path(dedup_bam_t), path(recal_t)
      //tuple val(sample_id_t), path(dedup_bam_t), path(recal_t)

    output:
      tuple val(sample_id_t), path("*.vcf*")
      path(".command.{log,err,out}")

    script:
      def exome_option = params.exome_bed ? "--interval ${params.exome_bed}" : " "
    """
      ${params.sentieon_dir}/sentieon driver -r ${params.ref} -t ${task.cpus} ${exome_option} \
      -i ${dedup_bam_t[0]} -q ${recal_t} \
      -i ${dedup_bam_n[0]} -q ${recal_n} \
      --algo TNscope \
      --tumor_sample ${sample_id_t} --normal_sample ${sample_id_n} \
      -d ${params.smvdb1} ${sample_id_t}.vcf.gz
    """
}

//Phase variants using GATK’s ReadBackedPhasing, not supporting multiple threads
process GATK_PHASING {
    publishDir "${params.outdir}/${params.sample}/02_alignment_variant_calling_phasing", pattern: "*[_phased|_GT].vcf*", mode:'symlink'
    publishDir "${params.outdir}/${params.sample}/00_log", mode:'symlink',saveAs: {filename -> if (filename =~ /.*.log$/) "${params.sample}_GATK_PHASING.log"}
    tag "${params.sample}"
    //cpus 8
    label "m_cpus"
    //echo true

    input:
      tuple val(sample_id_n), path(vcf_n), val(sample_id_t), path(vcf_t), val(sample_id_t), path(dedup_bam_t), path(recal_t)

    output:
      path("*sorted_phased.vcf")
      path("*tumor_normal_GT.vcf*")
      path(".command.{log,err,out}")

    shell:
    '''
      #normal_only_vcf
      sed "s/!{sample_id_n}/!{sample_id_t}/" <(zcat !{vcf_n[0]}|grep -E -v "^chr[0-9a-zA-Z]+_|^HLA|^chrM|chrUn") > !{sample_id_n}_rename.vcf
      #awk 'NR==1' ${sample_id_n}_rename.vcf
      bgzip !{sample_id_n}_rename.vcf
      tabix !{sample_id_n}_rename.vcf.gz

      #tumor_only_vcf
      zcat !{vcf_t[0]}| \
      awk -v FS='\t' -v OFS='\t' '{if($0~/##/){print}else{print $1,$2,$3,$4,$5,$6,$7,$8,$9,$10}}'| \
      awk -v FS='\t' -v OFS='\t' '{if($0~/#/){print}else{split($9,a,":");split($10,b,":");print $1,$2,$3,$4,$5,$6,$7,$8,a[1],b[1]}}'| \
      grep -E -v "^chr[0-9a-zA-Z]+_|^HLA|^chrM|chrUn" > !{sample_id_t}_tumorOnly.vcf
      bgzip !{sample_id_t}_tumorOnly.vcf
      tabix !{sample_id_t}_tumorOnly.vcf.gz

      #CombineVariants
      gatk3 -T CombineVariants -R !{params.ref} -nt !{task.cpus} \
      --variant !{sample_id_n}_rename.vcf.gz \
      --variant !{sample_id_t}_tumorOnly.vcf.gz \
      -o !{params.sample}_germline_plus_somatic.vcf \
      --assumeIdenticalSamples

      #Sort combined VCF using Picard
      picard SortVcf \
      I=!{params.sample}_germline_plus_somatic.vcf \
      O=!{params.sample}_germline_plus_somatic_sorted.vcf \
      SEQUENCE_DICTIONARY=!{params.dict}

      #Phase variants using GATK ReadBackedPhasing
      gatk3 -T ReadBackedPhasing \
      -R !{params.ref} \
      -I !{dedup_bam_t[0]} \
      --variant !{params.sample}_germline_plus_somatic_sorted.vcf \
      -L !{params.sample}_germline_plus_somatic_sorted.vcf \
      -o !{params.sample}_germline_plus_somatic_sorted_phased.vcf

      if [ !{params.AF_method} = "sentieon" ] || [ !{params.AF_method} = "Sentieon" ] || [ !{params.AF_method} = "SENTIEON" ]
      then
          #tumor_normal_vcf, only GT info
          zcat !{vcf_t[0]}| \
          perl -lane 'if(/#/){print}else{@a=split /:/,$F[8];@b=split /:/,$F[9];@c=split /:/,$F[10];print join("\t",@F[0..7],join(":",@a[0,3,2,1]),join(":",@b[0,3,2,1]),join(":",@c[0,3,2,1]))}'| \
          sed 's/AFDP/DP/g'|grep -E -v "^chr[0-9a-zA-Z]+_|^HLA|^chrM|^chrUn"|grep -v "GT::AD" > !{params.sample}_tumor_normal_GT.vcf
          bgzip !{params.sample}_tumor_normal_GT.vcf
          tabix !{params.sample}_tumor_normal_GT.vcf.gz
      else
          #tumor_normal_vcf, only GT info
          zcat !{vcf_t[0]}| \
          perl -lane 'if(/#/){print}else{@GT=split /:/,$F[8];@t=split /:/,$F[9];@n=split /:/,$F[10];print join("\t",@F[0..7],$GT[0],$t[0],$n[0])}'| \
          grep -E -v "^chr[0-9a-zA-Z]+_|^HLA|^chrM|chrUn" > !{params.sample}_tumor_normal_GT.vcf
          bgzip !{params.sample}_tumor_normal_GT.vcf
          tabix !{params.sample}_tumor_normal_GT.vcf.gz
      fi
    '''
}

//VEP for phased_vcf and tumor_normal_vcf
process VEP {
    publishDir "${params.outdir}/${params.sample}/03_variant_annotation/vep", pattern: "*.vcf*", mode:'symlink'
    publishDir "${params.outdir}/${params.sample}/00_log", mode:'symlink',saveAs: {filename -> if (filename =~ /.*.log$/) "${params.sample}_VEP.log"}
    tag "${params.sample}"
    //cpus 1
    echo true
    //conda "/home/dxk/miniconda3/envs/vep"

    input:
      //tuple path(sorted_phased_vcf), path(tumor_normal_GT_vcf), path(tumor_normal_GT_vcf_tbi)
      path(sorted_phased_Tumor_normal_GT_vcf)

    output:
      path("*sorted_phased_vep_anno.vcf.gz*")
      path("*tumor_normal_GT_vep_anno.vcf.gz*")
      path(".command.{log,err,out}")

    script:
      def transcript_option = params.transcript_version ? "--transcript_version" : " "
    """
      #phased_vcf
      vep \
      --input_file ${params.sample}_germline_plus_somatic_sorted_phased.vcf \
      --output_file ${params.sample}_germline_plus_somatic_sorted_phased_vep_anno.vcf \
      --format vcf --vcf --symbol --terms SO --tsl \
      --hgvs --fasta ${params.ref} \
      --offline --cache --dir_cache ${params.cache} \
      --plugin Downstream --plugin Wildtype \
      --dir_plugins ${params.plugin} --pick ${transcript_option}

      bgzip ${params.sample}_germline_plus_somatic_sorted_phased_vep_anno.vcf
      tabix ${params.sample}_germline_plus_somatic_sorted_phased_vep_anno.vcf.gz

      #tumor_normal_vcf, only GT info
      vep \
      --input_file ${params.sample}_tumor_normal_GT.vcf.gz --output_file ${params.sample}_tumor_normal_GT_vep_anno.vcf \
      --format vcf --vcf --symbol --terms SO --tsl \
      --hgvs --fasta ${params.ref} \
      --offline --cache --dir_cache ${params.cache} \
      --plugin Frameshift --plugin Wildtype \
      --dir_plugins ${params.plugin} --pick ${transcript_option}

      bgzip ${params.sample}_tumor_normal_GT_vep_anno.vcf
      tabix ${params.sample}_tumor_normal_GT_vep_anno.vcf.gz
    """
}

process HISAT2_INDEX {
    //publishDir "${params.outdir}/${params.sample}/04_rna_expression/index", mode:'symlink'
    //tag "${params.sample}"
    //cpus 8
    label "m_cpus"
    echo true

    script:
    if( file(params.rnaIndexFiles_h) ) {
      """
      echo 'rna index files for hisat2 exist'
      """
    }
    else {
      """
      echo 'rna index files for hisat2 building'
      dir=\$(dirname ${params.ref})
      hisat2-build -p ${task.cpus} ${params.ref} \$dir/rna_index
      """
    }
}

process HISAT2_MAP {
    publishDir "${params.outdir}/${params.sample}/04_rna_expression", mode:'symlink'
    publishDir "${params.outdir}/${params.sample}/00_log", mode:'symlink',saveAs: {filename -> if (filename =~ /.*.log$/) "${params.sample}_HISAT2_MAP.log"}
    tag "${params.sample}"
    //cpus 16
    label "l_cpus"
    //echo true

    input:
      tuple val(sample_id), path(reads)
      //path index

    output:
      tuple val(sample_id), path("*.bam*")
      path("*flagstat")
      path(".command.{log,err,out}")

    when:
      sample_id =~ /(?i)RNA/
      //hisat2 -p $nt -x $ref -U $fq1 -S $sam
      //index=$(basename ${params.ref})
    script:
    """
      dir=\$(dirname ${params.ref})
      hisat2 -p ${task.cpus} -x \$dir/rna_index -1 ${reads[0]} -2 ${reads[1]} | \
      tee >(samtools flagstat -@ ${task.cpus} - > ${params.sample}_hisat2.flagstat) | \
      samtools sort -@ ${task.cpus} -O BAM | \
      tee ${params.sample}_rna_sorted.bam | \
      samtools index -@ ${task.cpus} - ${params.sample}_rna_sorted.bam.bai
    """
}

/*
process STRINGTIE {
    publishDir "${params.outdir}/${params.sample}/04_rna_expression/stringtie", mode:'symlink'
    publishDir "${params.outdir}/${params.sample}/00_log", mode:'symlink',saveAs: {filename -> if (filename =~ /.*.log$/) "${params.sample}_STRINGTIE.log"}
    tag "${params.sample}"
    //cpus 8
    label "m_cpus"
    //echo true

    input:
      tuple val(sample_id), path(rna_bam)

    output:
      path("*")
      path(".command.{log,err,out}")

    script:
    """
      stringtie ${rna_bam[0]} \
      -p ${task.cpus} \
      -G ${params.gtf} -e -B \
      -o ${sample_id}_stringtie_transcripts.gtf \
      -A ${sample_id}_stringtie_gene_expression.tsv -v
    """
}
*/

process KALLISTO {
    publishDir "${params.outdir}/${params.sample}/04_rna_expression/kallisto", mode:'symlink'
    publishDir "${params.outdir}/${params.sample}/00_log", mode:'symlink',saveAs: {filename -> if (filename =~ /.*.log$/) "${params.sample}_KALLISTO.log"}
    tag "${params.sample}"
    //cpus 8
    label "m_cpus"
    echo true

    input:
      tuple val(sample_id), path(reads)

    output:
      path("*")
      path(".command.{log,err,out}")

    when:
      sample_id =~ /(?i)RNA/

    script:
    if( file(params.rnaIndexFiles_k) ) {
      """
      echo 'rna index files for kallisto exist'
      dir=\$(dirname ${params.ref})
      kallisto quant \
      -i \$dir/rna_index_kallisto.idx \
      -o ./ -t ${task.cpus} -b 100 --fusion \
      ${reads[0]} ${reads[1]}

      paste ./abundance.tsv | cut -f 1,2,5| \
      awk '{if(\$0~/target_id/){print}else{split(\$1,a,"." ); print a[1]"\t"\$2"\t"\$3}}' > transcript_tpm.tsv

      kallisto_gene_matrix.pl --gtf_file=${params.gtf} \
      --kallisto_transcript_matrix_in=transcript_tpm.tsv \
      --kallisto_transcript_matrix_out=gene_tpm.tsv
      """
    }
    else {
      """
      echo 'rna index files for kallisto building'
      dir=\$(dirname ${params.ref})
      kallisto index ${params.cdna} -i \$dir/rna_index_kallisto.idx
      kallisto quant \
      -i \$dir/rna_index_kallisto.idx \
      -o ./ -t ${task.cpus} -b 100 --fusion \
      ${reads[0]} ${reads[1]}

      paste ./abundance.tsv | cut -f 1,2,5| \
      awk '{if(\$0~/target_id/){print}else{split(\$1,a,"." ); print a[1]"\t"\$2"\t"\$3}}' > transcript_tpm.tsv

      kallisto_gene_matrix.pl --gtf_file=${params.gtf} \
      --kallisto_transcript_matrix_in=transcript_tpm.tsv \
      --kallisto_transcript_matrix_out=gene_tpm.tsv
      """
    }
}

/*
process KALLISTO_QUANT {
    publishDir "${params.outdir}/${params.sample}/04_rna_expression/kallisto", mode:'symlink'
    tag "${params.sample}"
    cpus 8
    echo true

    input:
      tuple val(sample_id), path(reads)

    output:
      tuple val(sample_id), path("*")

    when:
      sample_id =~ /(?i)RNA/

    script:
    """
      dir=\$(dirname ${params.ref})
      kallisto quant \
      -i \$dir/rna_index_kallisto.idx \
      -o ./ -t ${task.cpus} -b 100 --fusion \
      ${reads[0]} ${reads[1]}
    """
}
*/

process COVERAGE {
    publishDir "${params.outdir}/${params.sample}/03_variant_annotation/coverage", mode:'symlink'
    publishDir "${params.outdir}/${params.sample}/00_log", mode:'symlink',saveAs: {filename -> if (filename =~ /.*.log$/) "${params.sample}_COVERAGE.log"}
    tag "${params.sample}"
    //cpus 4
    label "s_cpus"
    //echo true

    input:
      //tuple path(tumor_normal_GT_vep_anno_vcf) ,path(tumor_normal_GT_vep_anno_vcf_tbi), path(tumor_normal_GT_vep_anno_vcf_summary), val(sample_id_n), path(dedup_bam_n), path(recal_n), val(sample_id_t), path(dedup_bam_t), path(recal_t), val(sample_id_rna), path(rna_sorted_bam)
      tuple path(tumor_normal_GT_vep_anno_vcf) ,path(tumor_normal_GT_vep_anno_vcf_tbi), val(sample_id_n), path(dedup_bam_n), path(recal_n), val(sample_id_t), path(dedup_bam_t), path(recal_t), val(sample_id_rna), path(rna_sorted_bam)
      //path(coverage_input_files)

    output:
      path("*_anno_decomposed_snv1_indel1_snv2_indel2_RNA_snv_indel.vcf")
      path(".command.{log,err,out}")

    shell:
    '''
      #vt decompose
      gunzip -c !{params.sample}_tumor_normal_GT_vep_anno.vcf.gz > !{params.sample}_tumor_normal_GT_vep_anno.vcf
      vt decompose -s !{params.sample}_tumor_normal_GT_vep_anno.vcf -o !{params.sample}_anno_decomposed.vcf

      grep -v '^#' !{params.sample}_anno_decomposed.vcf| \
      awk -v OFS="\t" '{print $1, $2-1, $2+1}' > !{params.sample}_anno_decomposed.bed
      mkdir split_decomposed_vcf
      for i in {1..22} X Y;do grep -E -w "#|CHROM|chr$i" !{params.sample}_anno_decomposed.vcf > ./split_decomposed_vcf/!{params.sample}_anno_decomposed_chr${i}.vcf;done

      if [ !{params.AF_method} = "sentieon" ] || [ !{params.AF_method} = "Sentieon" ] || [ !{params.AF_method} = "SENTIEON" ]
      then
      #DNA tumor
      SAMPLE1=$(cat !{params.sample}_anno_decomposed.vcf|grep "#CHROM"|cut -f 10)

      #RNA tumor
      samtools view -@ !{task.cpus} -b -L !{params.sample}_anno_decomposed.bed !{params.sample}_rna_sorted.bam > !{params.sample}_rna_sorted_subset.bam
      samtools index -@ !{task.cpus} !{params.sample}_rna_sorted_subset.bam

      mkdir !{params.sample}_RNA_readcount
      for i in {1..22} X Y;do bam_readcount_helper_2019.py \
      ./split_decomposed_vcf/!{params.sample}_anno_decomposed_chr${i}.vcf \
      ${SAMPLE1} !{params.ref} !{params.sample}_rna_sorted_subset.bam chr${i} !{params.sample}_RNA_readcount;done
      cd !{params.sample}_RNA_readcount
      cat chr*_bam_readcount_snv.tsv > !{params.sample}_bam_readcount_snv.tsv
      cat chr*_bam_readcount_indel.tsv > !{params.sample}_bam_readcount_indel.tsv
      #vatools
      cd ..
      vcf-readcount-annotator !{params.sample}_anno_decomposed.vcf \
      ./!{params.sample}_RNA_readcount/!{params.sample}_bam_readcount_snv.tsv \
      RNA -s ${SAMPLE1} -t snv -o !{params.sample}_anno_decomposed_snv1_indel1_snv2_indel2_RNA_snv.vcf

      vcf-readcount-annotator !{params.sample}_anno_decomposed_snv1_indel1_snv2_indel2_RNA_snv.vcf \
      ./!{params.sample}_RNA_readcount/!{params.sample}_bam_readcount_indel.tsv \
      RNA -s ${SAMPLE1} -t indel -o !{params.sample}_anno_decomposed_snv1_indel1_snv2_indel2_RNA_snv_indel.vcf
      else
      #DNA tumor
      #bam_readcount
      SAMPLE1=$(cat !{params.sample}_anno_decomposed.vcf|grep "#CHROM"|cut -f 10)
      samtools view -@ !{task.cpus} -b -L !{params.sample}_anno_decomposed.bed ${SAMPLE1}_sorted_deduped_indelRealigned.bam > ${SAMPLE1}_sorted_deduped_subset.bam
      samtools index -@ !{task.cpus} ${SAMPLE1}_sorted_deduped_subset.bam

      mkdir ${SAMPLE1}_readcount
      for i in {1..22} X Y;do bam_readcount_helper_2019.py \
      ./split_decomposed_vcf/!{params.sample}_anno_decomposed_chr${i}.vcf \
      ${SAMPLE1} !{params.ref} ${SAMPLE1}_sorted_deduped_subset.bam chr${i} ${SAMPLE1}_readcount;done
      cd ${SAMPLE1}_readcount
      cat chr*_bam_readcount_snv.tsv > ${SAMPLE1}_bam_readcount_snv.tsv
      cat chr*_bam_readcount_indel.tsv > ${SAMPLE1}_bam_readcount_indel.tsv
      cd ..
      #vatools
      vcf-readcount-annotator !{params.sample}_anno_decomposed.vcf \
      ./${SAMPLE1}_readcount/${SAMPLE1}_bam_readcount_snv.tsv \
      DNA -s ${SAMPLE1} -t snv -o !{params.sample}_anno_decomposed_snv1.vcf

      vcf-readcount-annotator !{params.sample}_anno_decomposed_snv1.vcf \
      ./${SAMPLE1}_readcount/${SAMPLE1}_bam_readcount_indel.tsv \
      DNA -s ${SAMPLE1} -t indel -o !{params.sample}_anno_decomposed_snv1_indel1.vcf

      #DNA normal
      #bam_readcount
      SAMPLE2=$(cat !{params.sample}_anno_decomposed.vcf|grep "#CHROM"|cut -f 11)
      samtools view -@ !{task.cpus} -b -L !{params.sample}_anno_decomposed.bed ${SAMPLE2}_sorted_deduped_indelRealigned.bam > ${SAMPLE2}_sorted_deduped_subset.bam
      samtools index -@ !{task.cpus} ${SAMPLE2}_sorted_deduped_subset.bam

      mkdir ${SAMPLE2}_readcount
      for i in {1..22} X Y;do bam_readcount_helper_2019.py \
      ./split_decomposed_vcf/!{params.sample}_anno_decomposed_chr${i}.vcf \
      ${SAMPLE2} !{params.ref} ${SAMPLE2}_sorted_deduped_subset.bam chr${i} ${SAMPLE2}_readcount;done
      cd ${SAMPLE2}_readcount
      cat chr*_bam_readcount_snv.tsv > ${SAMPLE2}_bam_readcount_snv.tsv
      cat chr*_bam_readcount_indel.tsv > ${SAMPLE2}_bam_readcount_indel.tsv
      #vatools
      cd ..
      vcf-readcount-annotator !{params.sample}_anno_decomposed_snv1_indel1.vcf \
      ./${SAMPLE2}_readcount/${SAMPLE2}_bam_readcount_snv.tsv \
      DNA -s ${SAMPLE2} -t snv -o !{params.sample}_anno_decomposed_snv1_indel1_snv2.vcf

      vcf-readcount-annotator !{params.sample}_anno_decomposed_snv1_indel1_snv2.vcf \
      ./${SAMPLE2}_readcount/${SAMPLE2}_bam_readcount_indel.tsv \
      DNA -s ${SAMPLE2} -t indel -o !{params.sample}_anno_decomposed_snv1_indel1_snv2_indel2.vcf

      #RNA tumor
      samtools view -@ !{task.cpus} -b -L !{params.sample}_anno_decomposed.bed !{params.sample}_rna_sorted.bam > !{params.sample}_rna_sorted_subset.bam
      samtools index -@ !{task.cpus} !{params.sample}_rna_sorted_subset.bam

      mkdir !{params.sample}_RNA_readcount
      for i in {1..22} X Y;do bam_readcount_helper_2019.py \
      ./split_decomposed_vcf/!{params.sample}_anno_decomposed_chr${i}.vcf \
      ${SAMPLE1} !{params.ref} !{params.sample}_rna_sorted_subset.bam chr${i} !{params.sample}_RNA_readcount;done
      cd !{params.sample}_RNA_readcount
      cat chr*_bam_readcount_snv.tsv > !{params.sample}_bam_readcount_snv.tsv
      cat chr*_bam_readcount_indel.tsv > !{params.sample}_bam_readcount_indel.tsv
      #vatools
      cd ..
      vcf-readcount-annotator !{params.sample}_anno_decomposed_snv1_indel1_snv2_indel2.vcf \
      ./!{params.sample}_RNA_readcount/!{params.sample}_bam_readcount_snv.tsv \
      RNA -s ${SAMPLE1} -t snv -o !{params.sample}_anno_decomposed_snv1_indel1_snv2_indel2_RNA_snv.vcf

      vcf-readcount-annotator !{params.sample}_anno_decomposed_snv1_indel1_snv2_indel2_RNA_snv.vcf \
      ./!{params.sample}_RNA_readcount/!{params.sample}_bam_readcount_indel.tsv \
      RNA -s ${SAMPLE1} -t indel -o !{params.sample}_anno_decomposed_snv1_indel1_snv2_indel2_RNA_snv_indel.vcf
      fi
    '''
}

process EXPRESSION {
    publishDir "${params.outdir}/${params.sample}/03_variant_annotation/expression", mode:'symlink'
    publishDir "${params.outdir}/${params.sample}/00_log", mode:'symlink',saveAs: {filename -> if (filename =~ /.*.log$/) "${params.sample}_EXPRESSION.log"}
    tag "${params.sample}"
    //echo true

    input:
      path(coverage_vcf_Gene_Transcript_expression)

    output:
      path("*_gene_tran.vcf.gz*")
      path(".command.{log,err,out}")

    shell:
      transcript_option = params.transcript_version ? " " : "--ignore-ensembl-id-version"
    '''
      SAMPLE1=$(cat !{params.sample}_anno_decomposed_snv1_indel1_snv2_indel2_RNA_snv_indel.vcf|grep "#CHROM"|cut -f 10)

      awk '{if($0~/target_id/){print "Gene ID""\t""TPM"}else{print $1"\t"$3}}' \
      gene_tpm.tsv > gene_tpm_v.tsv

      #gene
      #kallisto
      exp_gene_tool=stringtie
      vcf-expression-annotator !{params.sample}_anno_decomposed_snv1_indel1_snv2_indel2_RNA_snv_indel.vcf \
      gene_tpm_v.tsv ${exp_gene_tool} gene -s ${SAMPLE1} \
      -o !{params.sample}_anno_decomposed_snv1_indel1_snv2_indel2_RNA_snv_indel_gene.vcf

      #transcript
      #kallisto
      exp_tran_tool=kallisto
      vcf-expression-annotator !{params.sample}_anno_decomposed_snv1_indel1_snv2_indel2_RNA_snv_indel_gene.vcf \
      abundance.tsv ${exp_tran_tool} transcript -s ${SAMPLE1} \
      -o !{params.sample}_anno_decomposed_snv1_indel1_snv2_indel2_RNA_snv_indel_gene_tran.vcf !{transcript_option}

      #GT:DP:AF:AD:RDP:RAF:RAD:GX:TX	0/1:103:0.14563:88,15:34:0.0:34,0:ENSG00000186827|4.011825:ENST00000379236|8.68948

      bgzip !{params.sample}_anno_decomposed_snv1_indel1_snv2_indel2_RNA_snv_indel_gene_tran.vcf
      tabix !{params.sample}_anno_decomposed_snv1_indel1_snv2_indel2_RNA_snv_indel_gene_tran.vcf.gz
      '''
}

process OPTITYPE {
    publishDir "${params.outdir}/${params.sample}/05_hla_i_typing", mode:'symlink'
    publishDir "${params.outdir}/${params.sample}/00_log", mode:'symlink',saveAs: {filename -> if (filename =~ /.*.log$/) "${sample_id}_OPTITYPE.log"}
    tag "${params.sample}"
    //cpus 16
    label "l_cpus"
    //echo true

    input:
      tuple val(sample_id), path(reads)

    output:
      path("*_HLA_I.txt")
      path(".command.{log,err,out}")

    when:
      //!(sample_id =~ /(?i)RNA/)
      //sample_id =~ /(?i)normal|-n/
      //sample_id =~ /(?i)tumor_DNA|-t|(?i)tumor_RNA/
      sample_id =~ /(?i)tumor_DNA|-t/

    //razers3,optitype
    shell:
    '''

      if [[ !{sample_id} =~ "DNA" ]] || [[ !{sample_id} =~ "dna" ]]
      then
      seq_type=dna
      else
      seq_type=rna
      fi

      razers3 -i 95 -m 1 -dr 0 -tc !{task.cpus} \
      -o !{sample_id}_1.bam \
      !{params.opti_hla_ref_path}/hla_reference_${seq_type}.fasta \
      !{reads[0]}
      samtools bam2fq -@ !{task.cpus} !{sample_id}_1.bam > !{sample_id}_1_fished.fastq

      razers3 -i 95 -m 1 -dr 0 -tc !{task.cpus} \
      -o !{sample_id}_2.bam \
      !{params.opti_hla_ref_path}/hla_reference_${seq_type}.fasta \
      !{reads[1]}
      samtools bam2fq -@ !{task.cpus} !{sample_id}_2.bam > !{sample_id}_2_fished.fastq

      singularity run !{params.opti_image} \
      -i !{sample_id}_1_fished.fastq !{sample_id}_2_fished.fastq \
      --${seq_type} -v -p !{sample_id} -o ./

      grep -v "Reads" !{sample_id}_result.tsv|xargs -n 1|grep -E "A|B|C"| \
      awk '{print "HLA-"$0}'|xargs|sed 's/ /,/g' > !{sample_id}_HLA_I.txt
    '''
}

process HLA_LA {
    publishDir "${params.outdir}/${params.sample}/08_hla_ii_typing", mode:'copy'
    publishDir "${params.outdir}/${params.sample}/00_log", mode:'symlink',saveAs: {filename -> if (filename =~ /.*.log$/) "${sample_id}_HLA_LA.log"}
    tag "${params.sample}"
    //cpus 16
    label "l_cpus"
    //conda "/home/dxk/miniconda3/envs/hlala"
    //echo true

    input:
      tuple val(sample_id), path(bam)

    output:
      path("*_HLA_II.txt")
      path("./${sample_id}/hla/R1_bestguess_G.txt")
      path(".command.{log,err,out}")

    when:
      //!(sample_id =~ /(?i)RNA/)
      //sample_id =~ /(?i)normal|-n/
      sample_id =~ /(?i)tumor_DNA|-t/

    shell:
    '''
      HLA-LA.pl --BAM !{bam[0]} --graph PRG_MHC_GRCh38_withIMGT --sampleID !{sample_id} --maxThreads !{task.cpus} --workingDir ./

      cut -f 3,4,6,8-9,11-12 ./!{sample_id}/hla/R1_bestguess_G.txt| \
      awk -v OFS="\t" '{if($2>=0.9 && $5>=0.9){print $1}}'|grep -E -v "^[ABCEFG]"| \
      awk '{if($0~/^D/){print}else{print "HLA-"$0}}'|cut -d : -f 1-2|sort|uniq|xargs|sed 's/ /,/g' \
      > ./!{sample_id}_HLA_II_individual.txt

      #one per line for HLA II
      sed 's/,/ /g' !{sample_id}_HLA_II_individual.txt|xargs -n 1 > !{sample_id}_hla_ii_individual.txt
      #xxx_hla_ii_combine.txt
      comb_hla_ii.py --hla_ii_list !{sample_id}_hla_ii_individual.txt --sample_name !{sample_id} --output_dir ./ 
      cat !{sample_id}_hla_ii_combine.txt|xargs|sed 's/ /,/g' > !{sample_id}_HLA_II_combine.txt

      cat !{sample_id}_HLA_II_individual.txt !{sample_id}_HLA_II_combine.txt|xargs|sed 's/ /,/g' \
      > !{sample_id}_HLA_II.txt
    '''
}

process HLA_I_II {
    publishDir "${params.outdir}/${params.sample}/hla_I_II", mode:'symlink'
    publishDir "${params.outdir}/${params.sample}/00_log", mode:'symlink',saveAs: {filename -> if (filename =~ /.*.log$/) "${sample_id}_HLA_I_II.log"}
    tag "${params.sample}"
    //echo true

    input:
      path(HLA_I_txt_HLA_II_txt)

    output:
      //tuple val(sample_id), path("*_HLA_I_II.txt")
      path("*_HLA_I_II.txt")
      path(".command.{log,err,out}")

    when:
      //!(sample_id =~ /(?i)RNA/)
      //sample_id =~ /(?i)normal|-n/
      sample_id =~ /(?i)tumor_DNA|-t/

    shell:
    '''
      cat !{sample_id}_HLA_I.txt !{sample_id}_HLA_II.txt|xargs|sed 's/ /,/g' > !{sample_id}_HLA_I_II.txt
    '''
}

pass_only_option = params.pass_only ? "--pass-only" : " "
process PVACSEQ_HLA_I {
    publishDir "${params.outdir}/${params.sample}/06_hla_i_neoantigen_detection", pattern: "*", mode:'symlink'
    publishDir "${params.outdir}/${params.sample}/00_log", mode:'symlink',saveAs: {filename -> if (filename =~ /.*.log$/) "${params.sample}_PVACSEQ_HLA_I.log"}
    tag "${params.sample}"
    //cpus 16
    label "l_cpus"
    //conda "/home/dxk/miniconda3/envs/pvactools"

    input:
      //tuple path(anno_cov_exp_vcf), path(anno_cov_exp_vcf_tbi), path(phased_anno_vcf), path(phased_anno_vcf_tbi), path(phased_anno_vcf_summary), path(HLA_I_II_txt)
      path(pvacseq_input_files)

    output:
      //path("**/*.{yml,R,json,tsv,fasta}")
      //path("pvacseq_${params.sample}_*")
      path("*")
      path(".command.{log,err,out}")

    //HLA_I=$(cat ${SAMPLE2}_HLA_I.txt)
    //HLA_II=$(cat ${SAMPLE2}_HLA_II.txt)
    //HLA_II_combine=$(cat ${SAMPLE2}_HLA_II_combine.txt)
    shell:
      //def pass_only_option = params.pass_only ? "--pass-only" : " "
    '''
      bcftools view -v snps !{params.sample}_germline_plus_somatic_sorted_phased_vep_anno.vcf.gz \
      > !{params.sample}_germline_plus_somatic_sorted_phased_vep_anno_snp.vcf
      bgzip !{params.sample}_germline_plus_somatic_sorted_phased_vep_anno_snp.vcf
      tabix !{params.sample}_germline_plus_somatic_sorted_phased_vep_anno_snp.vcf.gz

      bcftools view -i "FORMAT/AF >=!{params.PVACSEQ_AF} && FORMAT/RAF >=!{params.PVACSEQ_RAF}" !{params.sample}_anno_decomposed_snv1_indel1_snv2_indel2_RNA_snv_indel_gene_tran.vcf.gz \
      > !{params.sample}_anno_decomposed_snv1_indel1_snv2_indel2_RNA_snv_indel_gene_tran_AF!{params.PVACSEQ_AF}_RAF!{params.PVACSEQ_RAF}.vcf
      bgzip !{params.sample}_anno_decomposed_snv1_indel1_snv2_indel2_RNA_snv_indel_gene_tran_AF!{params.PVACSEQ_AF}_RAF!{params.PVACSEQ_RAF}.vcf
      tabix !{params.sample}_anno_decomposed_snv1_indel1_snv2_indel2_RNA_snv_indel_gene_tran_AF!{params.PVACSEQ_AF}_RAF!{params.PVACSEQ_RAF}.vcf.gz

      ivcf=!{params.sample}_anno_decomposed_snv1_indel1_snv2_indel2_RNA_snv_indel_gene_tran_AF!{params.PVACSEQ_AF}_RAF!{params.PVACSEQ_RAF}.vcf.gz
      pvcf=!{params.sample}_germline_plus_somatic_sorted_phased_vep_anno_snp.vcf.gz
      inum=$(zgrep -v "#" $ivcf|wc -l)
      pnum=$(zgrep -v "#" $pvcf|wc -l)
      echo "The number of somatic SNVs and indes is: $inum"
      echo "The number of phased SNVs and indes is: $pnum"
      #tumor_DNA
      SAMPLE1=$(zcat $ivcf|grep "#CHROM"|cut -f 10)
      #normal_DNA
      SAMPLE2=$(zcat $ivcf|grep "#CHROM"|cut -f 11)

      HLA=$(cat ${SAMPLE1}_HLA_*.txt)
      vac_method_sh=!{params.vac_method}

      pvacseq run \
      --normal-sample-name ${SAMPLE2} \
      -p $pvcf -t !{task.cpus} !{pass_only_option} \
      --iedb-install-directory !{params.iedb_dir} \
      $ivcf ${SAMPLE1} \
      ${HLA} ${vac_method_sh//-/ } \
      pvacseq_!{params.sample}_${vac_method_sh}
    '''
}

process PVACSEQ_HLA_II {
    publishDir "${params.outdir}/${params.sample}/09_hla_ii_neoantigen_detection", pattern: "*", mode:'symlink'
    publishDir "${params.outdir}/${params.sample}/00_log", mode:'symlink',saveAs: {filename -> if (filename =~ /.*.log$/) "${params.sample}_PVACSEQ_HLA_II.log"}
    tag "${params.sample}"
    //cpus 16
    label "l_cpus"
    //conda "/home/dxk/miniconda3/envs/pvactools"

    input:
      //tuple path(anno_cov_exp_vcf), path(anno_cov_exp_vcf_tbi), path(phased_anno_vcf), path(phased_anno_vcf_tbi), path(phased_anno_vcf_summary), path(HLA_I_II_txt)
      path(pvacseq_input_files)

    output:
      //path("**/*.{yml,R,json,tsv,fasta}")
      //path("pvacseq_${params.sample}_*")
      path("*")
      path(".command.{log,err,out}")

    //HLA_I=$(cat ${SAMPLE2}_HLA_I.txt)
    //HLA_II=$(cat ${SAMPLE2}_HLA_II.txt)
    //HLA_II_combine=$(cat ${SAMPLE2}_HLA_II_combine.txt)
    shell:
      //def pass_only_option = params.pass_only ? "--pass-only" : " "
    '''
      bcftools view -v snps !{params.sample}_germline_plus_somatic_sorted_phased_vep_anno.vcf.gz \
      > !{params.sample}_germline_plus_somatic_sorted_phased_vep_anno_snp.vcf
      bgzip !{params.sample}_germline_plus_somatic_sorted_phased_vep_anno_snp.vcf
      tabix !{params.sample}_germline_plus_somatic_sorted_phased_vep_anno_snp.vcf.gz

      bcftools view -i "FORMAT/AF >=!{params.PVACSEQ_AF} && FORMAT/RAF >=!{params.PVACSEQ_RAF}" !{params.sample}_anno_decomposed_snv1_indel1_snv2_indel2_RNA_snv_indel_gene_tran.vcf.gz \
      > !{params.sample}_anno_decomposed_snv1_indel1_snv2_indel2_RNA_snv_indel_gene_tran_AF!{params.PVACSEQ_AF}_RAF!{params.PVACSEQ_RAF}.vcf
      bgzip !{params.sample}_anno_decomposed_snv1_indel1_snv2_indel2_RNA_snv_indel_gene_tran_AF!{params.PVACSEQ_AF}_RAF!{params.PVACSEQ_RAF}.vcf
      tabix !{params.sample}_anno_decomposed_snv1_indel1_snv2_indel2_RNA_snv_indel_gene_tran_AF!{params.PVACSEQ_AF}_RAF!{params.PVACSEQ_RAF}.vcf.gz

      ivcf=!{params.sample}_anno_decomposed_snv1_indel1_snv2_indel2_RNA_snv_indel_gene_tran_AF!{params.PVACSEQ_AF}_RAF!{params.PVACSEQ_RAF}.vcf.gz
      pvcf=!{params.sample}_germline_plus_somatic_sorted_phased_vep_anno_snp.vcf.gz
      inum=$(zgrep -v "#" $ivcf|wc -l)
      pnum=$(zgrep -v "#" $pvcf|wc -l)
      echo "The number of somatic SNVs and indes is: $inum"
      echo "The number of phased SNVs and indes is: $pnum"
      #tumor_DNA
      SAMPLE1=$(zcat $ivcf|grep "#CHROM"|cut -f 10)
      #normal_DNA
      SAMPLE2=$(zcat $ivcf|grep "#CHROM"|cut -f 11)

      HLA=$(cat ${SAMPLE1}_HLA_*.txt)
      vac_method_sh=!{params.vac_method}

      pvacseq run \
      --normal-sample-name ${SAMPLE2} \
      -p $pvcf -t !{task.cpus} !{pass_only_option} \
      --iedb-install-directory !{params.iedb_dir} \
      $ivcf ${SAMPLE1} \
      ${HLA} ${vac_method_sh//-/ } \
      pvacseq_!{params.sample}_${vac_method_sh}
    '''
}

process PVACSEQ_SPLIT_HLA {
    //publishDir "${params.outdir}/${params.sample}/pvacseq/${params.sample}_\${vac_method_sh}_\${HLA_label}", mode:'symlink'
    publishDir "${params.outdir}/${params.sample}/pvacseq_split_hla", mode:'symlink'
    publishDir "${params.outdir}/${params.sample}/00_log", mode:'symlink',saveAs: {filename -> if (filename =~ /.*.log$/) "${params.sample}_PVACSEQ_SPLIT_HLA.log"}
    tag "${params.sample}"
    //cpus 16
    label "l_cpus"
    //conda "/home/dxk/miniconda3/envs/pvactools"
    //echo true

    input:
      //tuple path(anno_cov_exp_vcf), path(anno_cov_exp_vcf_tbi), path(phased_anno_vcf), path(phased_anno_vcf_tbi), path(phased_anno_vcf_summary), path(HLA_I_II_txt)
      path(pvacseq_input_files)
      each HLA_ele
      //each e1_length

    output:
      //path("**/*.{yml,R,json,tsv,fasta}")
      path("pvacseq_${params.sample}_*")
      path(".command.{log,err,out}")

    //HLA_I=$(cat ${SAMPLE2}_HLA_I.txt)
    //HLA_II=$(cat ${SAMPLE2}_HLA_II.txt)
    //HLA_II_combine=$(cat ${SAMPLE2}_HLA_II_combine.txt)
    shell:
      //def pass_only_option = params.pass_only ? "--pass-only" : " "
    '''
      bcftools view -v snps !{params.sample}_germline_plus_somatic_sorted_phased_vep_anno.vcf.gz \
      > !{params.sample}_germline_plus_somatic_sorted_phased_vep_anno_snp.vcf
      bgzip !{params.sample}_germline_plus_somatic_sorted_phased_vep_anno_snp.vcf
      tabix !{params.sample}_germline_plus_somatic_sorted_phased_vep_anno_snp.vcf.gz

      bcftools view -i "FORMAT/AF >=!{params.PVACSEQ_AF} && FORMAT/RAF >=!{params.PVACSEQ_RAF}" !{params.sample}_anno_decomposed_snv1_indel1_snv2_indel2_RNA_snv_indel_gene_tran.vcf.gz \
      > !{params.sample}_anno_decomposed_snv1_indel1_snv2_indel2_RNA_snv_indel_gene_tran_AF!{params.PVACSEQ_AF}_RAF!{params.PVACSEQ_RAF}.vcf
      bgzip !{params.sample}_anno_decomposed_snv1_indel1_snv2_indel2_RNA_snv_indel_gene_tran_AF!{params.PVACSEQ_AF}_RAF!{params.PVACSEQ_RAF}.vcf
      tabix !{params.sample}_anno_decomposed_snv1_indel1_snv2_indel2_RNA_snv_indel_gene_tran_AF!{params.PVACSEQ_AF}_RAF!{params.PVACSEQ_RAF}.vcf.gz

      ivcf=!{params.sample}_anno_decomposed_snv1_indel1_snv2_indel2_RNA_snv_indel_gene_tran_AF!{params.PVACSEQ_AF}_RAF!{params.PVACSEQ_RAF}.vcf.gz
      pvcf=!{params.sample}_germline_plus_somatic_sorted_phased_vep_anno_snp.vcf.gz
      inum=$(zgrep -v "#" $ivcf|wc -l)
      pnum=$(zgrep -v "#" $pvcf|wc -l)
      echo "The number of somatic SNVs and indes is: $inum"
      echo "The number of phased SNVs and indes is: $pnum"
      #tumor_DNA
      SAMPLE1=$(zcat $ivcf|grep "#CHROM"|cut -f 10)
      #normal_DNA
      SAMPLE2=$(zcat $ivcf|grep "#CHROM"|cut -f 11)

      vac_method_sh=!{params.vac_method}
      HLA_label=$(sed -E -e 's/-/_/g' -e 's/\\*/_/g' -e 's/://g' <(echo !{HLA_ele}))

      pvacseq run \
      --normal-sample-name ${SAMPLE2} \
      -p $pvcf -t !{task.cpus} !{pass_only_option} -s 400 \
      --iedb-install-directory !{params.iedb_dir} \
      $ivcf ${SAMPLE1} \
      !{HLA_ele} ${vac_method_sh//-/ } \
      pvacseq_!{params.sample}_${vac_method_sh}_${HLA_label}
    '''
}
/*
      ivcf=!{params.sample}_anno_decomposed_snv1_indel1_snv2_indel2_RNA_snv_indel_gene_tran.vcf.gz
      pvcf=!{params.sample}_germline_plus_somatic_sorted_phased_vep_anno.vcf.gz
      #tumor_DNA
      SAMPLE1=$(zcat $ivcf|grep "#CHROM"|cut -f 10)
      #normal_DNA
      SAMPLE2=$(zcat $ivcf|grep "#CHROM"|cut -f 11)

      HLA=$(cat ${SAMPLE2}_HLA_I_II.txt)
      vac_method_sh=!{params.vac_method}

      pvacseq run \
      --normal-sample-name ${SAMPLE2} \
      -p $pvcf -t !{task.cpus} \
      --iedb-install-directory !{params.iedb_dir} \
      $ivcf ${SAMPLE1} \
      ${HLA} ${vac_method_sh//-/ } \
      !{params.sample}_${vac_method_sh}
*/

process PEPTIDE_POSTPROCESS_HLA_I_PVACSEQ {
    //publishDir "${params.outdir}/${params.sample}/pvacseq/${params.sample}_\${vac_method_sh}_\${HLA_label}", mode:'symlink'
    publishDir "${params.outdir}/${params.sample}/07_hla_i_neoantigen_filter_rank", pattern: "**/*.{txt,xlsx}", mode:'copy'
    publishDir "${params.outdir}/${params.sample}/00_log", mode:'symlink',saveAs: {filename -> if (filename =~ /.*.log$/) "${params.sample}_PEPTIDE_POSTPROCESS_HLA_I_PVACSEQ.log"}
    tag "${params.sample}"
    //conda "/home/dxk/miniconda3/envs/pvactools"

    input:
      path(pvacseq_all_epitopes)
    output:
      path("${params.sample}_tumor_DNA_HLA_I_all_epitopes/*.{txt,xlsx}")
      path(".command.{log,err,out}")

    //let kmer=!{params.peptide_length}/2
    shell:
    '''
    find ./pvacseq_!{params.sample}_*/MHC_Class_I/ -name "*all_epitopes.tsv"|xargs cat > !{params.sample}_tumor_DNA_HLA_I_all_epitopes_cache
    grep "Chromosome" !{params.sample}_tumor_DNA_HLA_I_all_epitopes_cache|head -n 1 > !{params.sample}_tumor_DNA_HLA_I_header
    grep -v "Chromosome" !{params.sample}_tumor_DNA_HLA_I_all_epitopes_cache > !{params.sample}_tumor_DNA_HLA_I_content
    cat !{params.sample}_tumor_DNA_HLA_I_header !{params.sample}_tumor_DNA_HLA_I_content > !{params.sample}_tumor_DNA_HLA_I_all_epitopes.tsv

    ivcf=!{params.sample}_anno_decomposed_snv1_indel1_snv2_indel2_RNA_snv_indel_gene_tran_AF!{params.PVACSEQ_AF}_RAF!{params.PVACSEQ_RAF}.vcf.gz
    pvcf=!{params.sample}_germline_plus_somatic_sorted_phased_vep_anno_snp.vcf.gz
    tsv=!{params.sample}_tumor_DNA_HLA_I_all_epitopes.tsv
    pvacseq generate_protein_fasta \
    $ivcf !{params.peptide_flank_length} !{params.sample}_flank!{params.peptide_flank_length}_phased.fasta -s !{params.sample}_tumor_DNA \
    --input-tsv $tsv -p $pvcf

    mkdir !{params.sample}_tumor_DNA_HLA_I_all_epitopes
    mv !{params.sample}_tumor_DNA_HLA_I_all_epitopes.tsv ./!{params.sample}_tumor_DNA_HLA_I_all_epitopes

    prepare_input_for_prime_netmhcstabpan.py --in_file_dir ./!{params.sample}_tumor_DNA_HLA_I_all_epitopes \
    --sample_name !{params.sample} --output_dir ./!{params.sample}_tumor_DNA_HLA_I_all_epitopes

    run_prime_netmhcstabpan.py --in_file_dir ./!{params.sample}_tumor_DNA_HLA_I_all_epitopes --sample_name !{params.sample} \
    --prime_path !{params.prime_path} --mixmhcpred_path !{params.mixmhcpred_path} --netmhcstabpan_path !{params.netmhcstabpan_path}

    cut -f 15,19 ./!{params.sample}_tumor_DNA_HLA_I_all_epitopes/!{params.sample}_tumor_DNA_HLA_I_all_epitopes.tsv| \
    sed -E -e "s/HLA Allele/mhc/" -e "s/MT Epitope Seq/pep/" -e "s/\t/,/" > !{params.sample}_tumor_DNA_bigmhc_input.txt

    python !{params.bigmhc_path}/src/predict.py -i=!{params.sample}_tumor_DNA_bigmhc_input.txt \
    -m=!{params.bigmhc_path}/models/bigmhc_!{params.bigmhc_model} -s=!{params.bigmhc_path}/data/pseudoseqs.csv -d cpu -o=!{params.sample}_tumor_DNA_bigmhc_!{params.bigmhc_model}_output.txt

    merge_results_for_pvac_prime_netmhcstabpan_bigmhc.py --sample_name !{params.sample} --path_pvacseq_out ./!{params.sample}_tumor_DNA_HLA_I_all_epitopes \
    --path_prime_netmhcstabpan_out ./!{params.sample}_tumor_DNA_HLA_I_all_epitopes --bigmhc_out !{params.sample}_tumor_DNA_bigmhc_!{params.bigmhc_model}_output.txt \
    --path_merge_result_out ./!{params.sample}_tumor_DNA_HLA_I_all_epitopes

    peptide_postprocessing_for_pvacseq.py --in_file ./!{params.sample}_tumor_DNA_HLA_I_all_epitopes/!{params.sample}_pvacseq_prime_netmhcstabpan_bigmhc.txt \
    --sample_name !{params.sample} --output_dir ./!{params.sample}_tumor_DNA_HLA_I_all_epitopes --tumor_dna_vaf 0.1

    peptide_postprocessing_for_pvacseq.py --in_file ./!{params.sample}_tumor_DNA_HLA_I_all_epitopes/!{params.sample}_pvacseq_prime_netmhcstabpan_bigmhc.txt \
    --sample_name !{params.sample} --output_dir ./!{params.sample}_tumor_DNA_HLA_I_all_epitopes --tumor_dna_vaf 0.02

    cd ./!{params.sample}_tumor_DNA_HLA_I_all_epitopes
    af_tsv=!{params.sample}_peptide_postprocessing_AF0.1.txt
    extend_epitope_dynamic_size_peptide_for_pvacseq.py -e ${af_tsv} -f ../!{params.sample}_*fasta.manufacturability.tsv -t ../pvacseq_!{params.sample}_*/MHC_Class_I/!{params.sample}_tumor_DNA.tsv \
    -s !{params.sample} -l !{params.peptide_flank_length} -x HLA_I_AF0.1

    af_tsv=!{params.sample}_peptide_postprocessing_AF0.02.txt
    extend_epitope_dynamic_size_peptide_for_pvacseq.py -e ${af_tsv} -f ../!{params.sample}_*fasta.manufacturability.tsv -t ../pvacseq_!{params.sample}_*/MHC_Class_I/!{params.sample}_tumor_DNA.tsv \
    -s !{params.sample} -l !{params.peptide_flank_length} -x HLA_I_AF0.02
    '''
}

process PEPTIDE_POSTPROCESS_HLA_II_PVACSEQ {
    //publishDir "${params.outdir}/${params.sample}/pvacseq/${params.sample}_\${vac_method_sh}_\${HLA_label}", mode:'symlink'
    publishDir "${params.outdir}/${params.sample}/10_hla_ii_neoantigen_filter_rank", pattern: "**/*.{txt,xlsx}", mode:'copy'
    publishDir "${params.outdir}/${params.sample}/00_log", mode:'symlink',saveAs: {filename -> if (filename =~ /.*.log$/) "${params.sample}_PEPTIDE_POSTPROCESS_HLA_II_PVACSEQ.log"}
    tag "${params.sample}"
    //conda "/home/dxk/miniconda3/envs/pvactools"

    input:
      path(pvacseq_all_epitopes)
    output:
      path("${params.sample}_tumor_DNA_HLA_II_all_epitopes/*.{txt,xlsx}")
      path(".command.{log,err,out}")

    shell:
    '''
    find ./pvacseq_!{params.sample}_*/MHC_Class_II/ -name "*all_epitopes.tsv"|xargs cat > !{params.sample}_tumor_DNA_HLA_II_all_epitopes_cache
    grep "Chromosome" !{params.sample}_tumor_DNA_HLA_II_all_epitopes_cache|head -n 1 > !{params.sample}_tumor_DNA_HLA_II_header
    grep -v "Chromosome" !{params.sample}_tumor_DNA_HLA_II_all_epitopes_cache > !{params.sample}_tumor_DNA_HLA_II_content
    cat !{params.sample}_tumor_DNA_HLA_II_header !{params.sample}_tumor_DNA_HLA_II_content > !{params.sample}_tumor_DNA_HLA_II_all_epitopes.tsv

    ivcf=!{params.sample}_anno_decomposed_snv1_indel1_snv2_indel2_RNA_snv_indel_gene_tran_AF!{params.PVACSEQ_AF}_RAF!{params.PVACSEQ_RAF}.vcf.gz
    pvcf=!{params.sample}_germline_plus_somatic_sorted_phased_vep_anno_snp.vcf.gz
    tsv=!{params.sample}_tumor_DNA_HLA_II_all_epitopes.tsv

    pvacseq generate_protein_fasta \
    $ivcf !{params.peptide_flank_length} !{params.sample}_flank!{params.peptide_flank_length}_phased.fasta -s !{params.sample}_tumor_DNA \
    --input-tsv $tsv -p $pvcf

    mkdir !{params.sample}_tumor_DNA_HLA_II_all_epitopes
    mv !{params.sample}_tumor_DNA_HLA_II_all_epitopes.tsv ./!{params.sample}_tumor_DNA_HLA_II_all_epitopes

    cd !{params.sample}_tumor_DNA_HLA_II_all_epitopes
    pvacseq binding_filter !{params.sample}_tumor_DNA_HLA_II_all_epitopes.tsv \
    !{params.sample}_tumor_DNA_HLA_II_all_epitopes_bf.tsv -b 500 -c 1 --exclude-NAs

    pvacseq coverage_filter !{params.sample}_tumor_DNA_HLA_II_all_epitopes_bf.tsv \
    !{params.sample}_tumor_DNA_HLA_II_all_epitopes_bf_cf_AF0.1.tsv \
    --normal-cov 5 --tdna-cov 10 --trna-cov 10 --normal-vaf 0.01 --tdna-vaf 0.1 --trna-vaf 0.01 --expn-val 1 --exclude-NAs

    pvacseq coverage_filter !{params.sample}_tumor_DNA_HLA_II_all_epitopes_bf.tsv \
    !{params.sample}_tumor_DNA_HLA_II_all_epitopes_bf_cf_AF0.02.tsv \
    --normal-cov 5 --tdna-cov 10 --trna-cov 10 --normal-vaf 0.01 --tdna-vaf 0.02 --trna-vaf 0.01 --expn-val 1 --exclude-NAs

    pvacseq transcript_support_level_filter !{params.sample}_tumor_DNA_HLA_II_all_epitopes_bf_cf_AF0.1.tsv \
    !{params.sample}_tumor_DNA_HLA_II_all_epitopes_bf_cf_AF0.1_tsl.tsv --maximum-transcript-support-level 1

    pvacseq transcript_support_level_filter !{params.sample}_tumor_DNA_HLA_II_all_epitopes_bf_cf_AF0.02.tsv \
    !{params.sample}_tumor_DNA_HLA_II_all_epitopes_bf_cf_AF0.02_tsl.tsv --maximum-transcript-support-level 1

    pvacseq top_score_filter !{params.sample}_tumor_DNA_HLA_II_all_epitopes_bf_cf_AF0.1_tsl.tsv \
    !{params.sample}_tumor_DNA_HLA_II_all_epitopes_bf_cf_AF0.1_tsl_tsf.tsv -m median

    pvacseq top_score_filter !{params.sample}_tumor_DNA_HLA_II_all_epitopes_bf_cf_AF0.02_tsl.tsv \
    !{params.sample}_tumor_DNA_HLA_II_all_epitopes_bf_cf_AF0.02_tsl_tsf.tsv -m median

    af_tsv=!{params.sample}_tumor_DNA_HLA_II_all_epitopes_bf_cf_AF0.1_tsl_tsf.tsv
    extend_epitope_dynamic_size_peptide_for_pvacseq.py -e ${af_tsv} -f ../!{params.sample}_*fasta.manufacturability.tsv -t ../pvacseq_!{params.sample}_*/MHC_Class_II/!{params.sample}_tumor_DNA.tsv \
    -s !{params.sample} -l !{params.peptide_flank_length} -x HLA_II_AF0.1
    af_tsv=!{params.sample}_tumor_DNA_HLA_II_all_epitopes_bf_cf_AF0.02_tsl_tsf.tsv
    extend_epitope_dynamic_size_peptide_for_pvacseq.py -e ${af_tsv} -f ../!{params.sample}_*fasta.manufacturability.tsv -t ../pvacseq_!{params.sample}_*/MHC_Class_II/!{params.sample}_tumor_DNA.tsv \
    -s !{params.sample} -l !{params.peptide_flank_length} -x HLA_II_AF0.02
    '''
}

process PEPTIDE_POSTPROCESS_HLA_I_II_PVACSEQ {
    //publishDir "${params.outdir}/${params.sample}/pvacseq/${params.sample}_\${vac_method_sh}_\${HLA_label}", mode:'symlink'
    publishDir "${params.outdir}/${params.sample}/11_hla_i_ii_neoantigen_combine", pattern: "*", mode:'copy'
    publishDir "${params.outdir}/${params.sample}/00_log", mode:'symlink',saveAs: {filename -> if (filename =~ /.*.log$/) "${params.sample}_PEPTIDE_POSTPROCESS_HLA_I_II_PVACSEQ.log"}
    tag "${params.sample}"
    //conda "/home/dxk/miniconda3/envs/pvactools"

    input:
      path(pvacseq_HLA_I_epitopes)
      path(pvacseq_HLA_II_epitopes)
    output:
      path("*")
      path(".command.{log,err,out}")

    shell:
    '''
    merge_hla_i_ii_peptide.py \
    -i pvacseq_!{params.sample}_flank!{params.peptide_flank_length}aa_peptide_HLA_I_AF0.1.txt \
    -t pvacseq_!{params.sample}_flank!{params.peptide_flank_length}aa_peptide_HLA_II_AF0.1.txt \
    -s !{params.sample} -l !{params.peptide_flank_length} -x HLA_I_II_AF0.1

    merge_hla_i_ii_peptide.py \
    -i pvacseq_!{params.sample}_flank!{params.peptide_flank_length}aa_peptide_HLA_I_AF0.02.txt \
    -t pvacseq_!{params.sample}_flank!{params.peptide_flank_length}aa_peptide_HLA_II_AF0.02.txt \
    -s !{params.sample} -l !{params.peptide_flank_length} -x HLA_I_II_AF0.02
    '''
}

process STAR_FUSION {
    publishDir "${params.outdir}/${params.sample}/fusion", mode:'symlink'
    publishDir "${params.outdir}/${params.sample}/00_log", mode:'symlink',saveAs: {filename -> if (filename =~ /.*.log$/) "${params.sample}_STAR_FUSION.log"}
    tag "${params.sample}"
    label "m_cpus"
    //echo true
    beforeScript "export SENTIEON_LICENSE=${params.sentieon_license}"

    input:
      tuple val(sample_id), path(reads)

    output:
      path("StarFusionOut/star-fusion_predictions_abridged_agfusion.tsv")
      path(".command.{log,err,out}")

    when:
      sample_id =~ /(?i)RNA/

    script:
      """
      ${params.sentieon_dir}/sentieon STAR --genomeDir ${params.ctat}/ref_genome.fa.star.idx \
      --runThreadN ${task.cpus} \
      --readFilesIn ${reads[0]} ${reads[1]} \
      --outReadsUnmapped None \
      --twopassMode Basic \
      --readFilesCommand "zcat" \
      --outSAMtype BAM Unsorted \
      --outSAMstrandField intronMotif \
      --outSAMunmapped Within \
      --chimSegmentMin 12 \
      --chimJunctionOverhangMin 8 \
      --chimOutJunctionFormat 1 \
      --alignSJDBoverhangMin 10 \
      --alignMatesGapMax 100000 \
      --alignIntronMax 100000 \
      --alignSJstitchMismatchNmax 5 -1 5 5 \
      --outSAMattrRGline ID:GRPundef \
      --chimMultimapScoreRange 3 \
      --chimScoreJunctionNonGTAG -4 \
      --chimMultimapNmax 20 \
      --chimNonchimScoreDropMin 10 \
      --peOverlapNbasesMin 12 \
      --peOverlapMMp 0.1 \
      --alignInsertionFlush Right \
      --alignSplicedMateMapLminOverLmate 0 \
      --alignSplicedMateMapLmin 30 \
      --quantMode GeneCounts

      singularity exec -e -B `pwd` -B ${params.ctat} \
      ${params.starFusion_simg} \
      STAR-Fusion \
      --genome_lib_dir ${params.ctat} \
      -J Chimeric.out.junction \
      -O StarFusionOut

      cd StarFusionOut
      cut -f 1-3,6-30 star-fusion.fusion_predictions.abridged.tsv > star-fusion_predictions_abridged_agfusion.tsv
      """
}
/*
      singularity exec -e -B `pwd` -B !{params.ctat} \
        !{params.starFusion_simg} \
        STAR-Fusion \
        --CPU !{task.cpus} \
        --left_fq !{reads[0]} \
        --right_fq !{reads[1]} \
        --genome_lib_dir !{params.ctat} \
        -O StarFusionOut \
        --FusionInspector validate \
        --examine_coding_effect \
        --denovo_reconstruct
*/

process EASYFUSE {
    publishDir "${params.outdir}/${params.sample}/16_fusion", mode:'copy'
    publishDir "${params.outdir}/${params.sample}/00_log", mode:'symlink',saveAs: {filename -> if (filename =~ /.*.log$/) "${params.sample}_EASYFUSE.log"}
    tag "${params.sample}"
    label "m_cpus"
    //echo true
    //beforeScript "export SENTIEON_LICENSE=${params.sentieon_license}"

    input:
      tuple val(sample_id), path(reads)

    output:
      path("FusionSummary/*_fusRank_1.pred.csv")
      path("FusionSummary/*_fusRank_1.pred.all.csv")
      path(".command.{log,err,out}")

    when:
      sample_id =~ /(?i)RNA/

    script:
      """
      mkdir rna_data
      cp -L ${reads[0]} ./rna_data
      cp -L ${reads[1]} ./rna_data
      singularity exec --containall \
      --bind ${params.easyfuse_ref}:/ref \
      --bind ./rna_data:/data \
      --bind ./:/output \
      ${params.easyfuse_sif} \
      python /code/easyfuse/processing.py -i /data/ -o /output -c /ref/config.ini
      """
}

process AGFUSION {
    publishDir "${params.outdir}/${params.sample}/16_fusion", mode:'symlink'
    publishDir "${params.outdir}/${params.sample}/00_log", mode:'symlink',saveAs: {filename -> if (filename =~ /.*.log$/) "${params.sample}_AGFUSION.log"}
    tag "${params.sample}"
    //cpus 8
    //conda "/home/dxk/miniconda3/envs/agfusion"
    //echo true
    beforeScript "export PYENSEMBL_CACHE_DIR=${params.pyensembl_cache}"

    input:
      path(fusion_result)

    output:
      //path("AGFusion/*/*.{fa,csv}")
      path("AGFusion/")
      path(".command.{log,err,out}")

    shell:
    '''
      agfusion batch \
      -f !{fusion_result} -a easyfuse \
      -db !{params.agfdb} \
      -o ./AGFusion \
      --middlestar --noncanonical
    '''
}

process PVACFUSE_HLA_I {
    publishDir "${params.outdir}/${params.sample}/17_hla_i_neoantigen_detection_fusion", mode:'symlink'
    publishDir "${params.outdir}/${params.sample}/00_log", mode:'symlink',saveAs: {filename -> if (filename =~ /.*.log$/) "${params.sample}_PVACFUSE_HLA_I.log"}
    tag "${params.sample}"
    label "l_cpus"
    //conda "/home/dxk/miniconda3/envs/pvactools"

    input:
      //path(agfusion_dir_HLA_I_II_txt)
      path(agfusion_dir)
      path(HLA_txt)

    output:
      //path("pvacfuse_${params.sample}_*")
      path("*")
      path(".command.{log,err,out}")

    shell:
    '''
      HLA=$(cat *_HLA_*.txt)
      vac_method_sh=!{params.vac_method}

      pvacfuse run \
      ./AGFusion \
      !{params.sample} \
      ${HLA} \
      ${vac_method_sh//-/ } \
      pvacfuse_!{params.sample}_${vac_method_sh} -t !{task.cpus} \
      --iedb-install-directory !{params.iedb_dir}
    '''
}

process PVACFUSE_HLA_II {
    publishDir "${params.outdir}/${params.sample}/19_hla_ii_neoantigen_detection_fusion", mode:'symlink'
    publishDir "${params.outdir}/${params.sample}/00_log", mode:'symlink',saveAs: {filename -> if (filename =~ /.*.log$/) "${params.sample}_PVACFUSE_HLA_II.log"}
    tag "${params.sample}"
    label "l_cpus"
    //conda "/home/dxk/miniconda3/envs/pvactools"

    input:
      //path(agfusion_dir_HLA_I_II_txt)
      path(agfusion_dir)
      path(HLA_txt)

    output:
      //path("pvacfuse_${params.sample}_*")
      path("*")
      path(".command.{log,err,out}")

    shell:
    '''
      HLA=$(cat *_HLA_*.txt)
      vac_method_sh=!{params.vac_method}

      pvacfuse run \
      ./AGFusion \
      !{params.sample} \
      ${HLA} \
      ${vac_method_sh//-/ } \
      pvacfuse_!{params.sample}_${vac_method_sh} -t !{task.cpus} \
      --iedb-install-directory !{params.iedb_dir}
    '''
}

process PEPTIDE_POSTPROCESS_HLA_I_PVACFUSE {
    //publishDir "${params.outdir}/${params.sample}/pvacseq/${params.sample}_\${vac_method_sh}_\${HLA_label}", mode:'symlink'
    publishDir "${params.outdir}/${params.sample}/18_hla_i_neoantigen_filter_rank_fusion", mode:'copy'
    publishDir "${params.outdir}/${params.sample}/00_log", mode:'symlink',saveAs: {filename -> if (filename =~ /.*.log$/) "${params.sample}_PEPTIDE_POSTPROCESS_HLA_I_PVACFUSE.log"}
    tag "${params.sample}"
    //conda "/home/dxk/miniconda3/envs/pvactools"

    input:
      path(agfusion_dir)
      path(pvacfuse_all_epitopes)
      path(easyfuse_all_pred)
      path(easyfuse_pos_pred)
    output:
      path("${params.sample}_tumor_DNA_HLA_I_all_epitopes/*.{tsv,txt,xlsx}")
      path(".command.{log,err,out}")

    shell:
    if ( file("./pvacfuse_${params.sample}_all/MHC_Class_I/${params.sample}.all_epitopes.tsv") ) {
    '''
    find ./pvacfuse_!{params.sample}_*/MHC_Class_I/ -maxdepth 1 -name "*all_epitopes.tsv"|xargs cat > !{params.sample}_tumor_DNA_HLA_I_all_epitopes_cache
    grep "Chromosome" !{params.sample}_tumor_DNA_HLA_I_all_epitopes_cache|head -n 1 > !{params.sample}_tumor_DNA_HLA_I_header
    grep -v "Chromosome" !{params.sample}_tumor_DNA_HLA_I_all_epitopes_cache > !{params.sample}_tumor_DNA_HLA_I_content
    cat !{params.sample}_tumor_DNA_HLA_I_header !{params.sample}_tumor_DNA_HLA_I_content > !{params.sample}_tumor_DNA_HLA_I_all_epitopes.tsv

    agfusion=./AGFusion
    tsv=./!{params.sample}_tumor_DNA_HLA_I_all_epitopes.tsv
    sed 's/Mutation/Index/' $tsv > !{params.sample}_tumor_DNA_HLA_I_all_epitopes_M2I.tsv
    pvacfuse generate_protein_fasta $agfusion !{params.peptide_flank_length} !{params.sample}_flank!{params.peptide_flank_length}_phased.fasta --input-tsv !{params.sample}_tumor_DNA_HLA_I_all_epitopes_M2I.tsv

    mkdir !{params.sample}_tumor_DNA_HLA_I_all_epitopes
    mv !{params.sample}_tumor_DNA_HLA_I_all_epitopes.tsv ./!{params.sample}_tumor_DNA_HLA_I_all_epitopes

    cd !{params.sample}_tumor_DNA_HLA_I_all_epitopes
    pvacfuse binding_filter !{params.sample}_tumor_DNA_HLA_I_all_epitopes.tsv \
    !{params.sample}_tumor_DNA_HLA_I_all_epitopes_bf.tsv -b 500 --exclude-NAs

    pvacfuse top_score_filter !{params.sample}_tumor_DNA_HLA_I_all_epitopes_bf.tsv \
    !{params.sample}_tumor_DNA_HLA_I_all_epitopes_bf_tsl.tsv -m median

    easyfuse_pred_filter.py ../*.pred.csv !{params.sample}_easyfuse_pred_filter.csv
    cat !{params.sample}_easyfuse_pred_filter.csv|perl -lane 'if(/(ENST.*?)_.*?(ENST.*?)\\s/){print "$1_$2"}' > easyfuse_enst_pair.txt
    grep -f <(cat <(echo "Chromosome") easyfuse_enst_pair.txt) !{params.sample}_tumor_DNA_HLA_I_all_epitopes_bf_tsl.tsv > !{params.sample}_tumor_DNA_HLA_I_all_epitopes_bf_tsl_enst.tsv

    af_tsv=!{params.sample}_tumor_DNA_HLA_I_all_epitopes_bf_tsl_enst.tsv
    extend_epitope_dynamic_size_peptide_for_pvacfuse.py -e $af_tsv -f ../!{params.sample}_*fasta.manufacturability.tsv -t ../pvacfuse_!{params.sample}_*/MHC_Class_I/8/!{params.sample}*.fa.tsv \
    -s !{params.sample} -l !{params.peptide_flank_length} -x HLA_I
    '''
    }
     else {
     '''
      mkdir !{params.sample}_tumor_DNA_HLA_I_all_epitopes
      cd !{params.sample}_tumor_DNA_HLA_I_all_epitopes
      touch !{params.sample}_tumor_DNA_HLA_I_all_epitopes.tsv
     '''
     }
}

process PEPTIDE_POSTPROCESS_HLA_II_PVACFUSE {
    //publishDir "${params.outdir}/${params.sample}/pvacseq/${params.sample}_\${vac_method_sh}_\${HLA_label}", mode:'symlink'
    publishDir "${params.outdir}/${params.sample}/20_hla_ii_neoantigen_filter_rank_fusion", mode:'copy'
    publishDir "${params.outdir}/${params.sample}/00_log", mode:'symlink',saveAs: {filename -> if (filename =~ /.*.log$/) "${params.sample}_PEPTIDE_POSTPROCESS_HLA_II_PVACFUSE.log"}
    tag "${params.sample}"
    conda "/home/dxk/miniconda3/envs/pvactools"

    input:
      path(agfusion_dir)
      path(pvacfuse_all_epitopes)
    output:
      path("${params.sample}_tumor_DNA_HLA_II_all_epitopes/*.{tsv,txt,xlsx}")
      path(".command.{log,err,out}")

    shell:
    if ( file("./pvacfuse_${params.sample}_*/MHC_Class_II/*all_epitopes.tsv") ) {
    '''
    find ./pvacfuse_!{params.sample}_*/MHC_Class_II/ -maxdepth 1 -name "*all_epitopes.tsv"|xargs cat > !{params.sample}_tumor_DNA_HLA_II_all_epitopes_cache
    grep "Chromosome" !{params.sample}_tumor_DNA_HLA_II_all_epitopes_cache|head -n 1 > !{params.sample}_tumor_DNA_HLA_II_header
    grep -v "Chromosome" !{params.sample}_tumor_DNA_HLA_II_all_epitopes_cache > !{params.sample}_tumor_DNA_HLA_II_content
    cat !{params.sample}_tumor_DNA_HLA_II_header !{params.sample}_tumor_DNA_HLA_II_content > !{params.sample}_tumor_DNA_HLA_II_all_epitopes.tsv

    agfusion=./AGFusion
    tsv=./!{params.sample}_tumor_DNA_HLA_II_all_epitopes.tsv
    sed 's/Mutation/Index/' $tsv > !{params.sample}_tumor_DNA_HLA_II_all_epitopes_M2I.tsv
    pvacfuse generate_protein_fasta $agfusion !{params.peptide_flank_length} !{params.sample}_flank!{params.peptide_flank_length}_phased.fasta --input-tsv !{params.sample}_tumor_DNA_HLA_II_all_epitopes_M2I.tsv

    mkdir !{params.sample}_tumor_DNA_HLA_II_all_epitopes
    mv !{params.sample}_tumor_DNA_HLA_II_all_epitopes.tsv ./!{params.sample}_tumor_DNA_HLA_II_all_epitopes

    cd !{params.sample}_tumor_DNA_HLA_II_all_epitopes
    pvacfuse binding_filter !{params.sample}_tumor_DNA_HLA_II_all_epitopes.tsv \
    !{params.sample}_tumor_DNA_HLA_II_all_epitopes_bf.tsv -b 500 --exclude-NAs

    pvacfuse top_score_filter !{params.sample}_tumor_DNA_HLA_II_all_epitopes_bf.tsv \
    !{params.sample}_tumor_DNA_HLA_II_all_epitopes_bf_tsl.tsv -m median

    af_tsv=!{params.sample}_tumor_DNA_HLA_II_all_epitopes_bf_tsl.tsv
    extend_epitope_dynamic_size_peptide_for_pvacfuse.py -e $af_tsv -f ../!{params.sample}_*fasta.manufacturability.tsv -t ../pvacfuse_!{params.sample}_*/MHC_Class_II/12/!{params.sample}*.fa.tsv \
    -s !{params.sample} -l !{params.peptide_flank_length} -x HLA_II
    '''
    }
     else {
     '''
      mkdir !{params.sample}_tumor_DNA_HLA_II_all_epitopes
      cd !{params.sample}_tumor_DNA_HLA_II_all_epitopes
      touch !{params.sample}_tumor_DNA_HLA_II_all_epitopes.tsv
     '''
     }
}

process CNVKIT {
    publishDir "${params.outdir}/${params.sample}/12_tp_tmb_tnb/cnvkit", pattern: "*.{cnr,seg,bed}", mode:'symlink'
    publishDir "${params.outdir}/${params.sample}/00_log", mode:'symlink',saveAs: {filename -> if (filename =~ /.*.log$/) "${sample_id_t}_CNVKIT.log"}
    tag "$sample_id_t"
    label "m_cpus"
    //echo true
    //conda "/home/dxk/miniconda3/envs/cnvkit"

    input:
      tuple val(sample_id_n), path(dedup_bam_n), path(recal_n), val(sample_id_t), path(dedup_bam_t), path(recal_t)
      //tuple val(sample_id_t), path(dedup_bam_t), path(recal_t)

    output:
      path("*.{cnr,seg}")
      path("*_baits.bed")
      path(".command.{log,err,out}")

    script:
      //def exome_option_cnvkit = params.exome_bed ? "--interval ${params.exome_bed}" : " "
      if(params.exome_bed) {
        """
        cnvkit.py batch \
        ${dedup_bam_t[0]} --normal ${dedup_bam_n[0]} \
        -m hybrid -p ${task.cpus} \
        --targets ${params.exome_bed} \
        --fasta ${params.ref} --access ${params.access_bed} --annotate ${params.refFlat} \
        --output-reference ${params.sample}_reference.cnn \
        --output-dir ./

        cnvkit.py export seg \
        ${params.sample}_tumor_DNA_sorted_deduped.cns --enumerate-chroms -o ${params.sample}_enumerate_chroms.seg
        echo "OK" > ${params.sample}_hg38_exon_ucsc_baits.bed
        """
      }
      else {
        """
        guess_baits.py \
        ${dedup_bam_t[0]} ${dedup_bam_n[0]} \
        -t ${params.exon_flatten_bed} -p ${task.cpus} -o ${params.sample}_hg38_exon_ucsc_baits.bed

        cnvkit.py batch \
        ${dedup_bam_t[0]} --normal ${dedup_bam_n[0]} \
        -m hybrid -p ${task.cpus} \
        --targets ${params.sample}_hg38_exon_ucsc_baits.bed \
        --fasta ${params.ref} --access ${params.access_bed} --annotate ${params.refFlat} \
        --output-reference ${params.sample}_reference.cnn \
        --output-dir ./

        cnvkit.py export seg \
        ${params.sample}_tumor_DNA_sorted_deduped.cns --enumerate-chroms -o ${params.sample}_enumerate_chroms.seg
        """
      }
}

process CALLSTATE {
    publishDir "${params.outdir}/${params.sample}/12_tp_tmb_tnb/callstate", mode:'symlink'
    publishDir "${params.outdir}/${params.sample}/00_log", mode:'symlink',saveAs: {filename -> if (filename =~ /.*.log$/) "${params.sample}_CALLSTATE.log"}
    tag "${params.sample}"
    label "m_cpus"
    //conda "/home/dxk/miniconda3/envs/callstate"

    input:
      tuple path(cnvkit_baits_bed), val(sample_id_t), path(dedup_bam_t), path(recal_t)
      //tuple val(sample_id_t), path(dedup_bam_t), path(recal_t)

    output:
      path("*_callable.bed")
      path(".command.{log,err,out}")

    shell:
      exome_option_callstate = params.exome_bed ? "${params.exome_bed}" : "${cnvkit_baits_bed}"
    '''
      awk -v OFS='\t' '{if($3-$2 >= 2){print}}' !{exome_option_callstate}|cut -f 1-3 > !{params.sample}_baits_gteq2bp.bed
      callstate -t !{task.cpus} -o !{params.sample}_baits_callstate.bed !{params.sample}_baits_gteq2bp.bed !{dedup_bam_t[0]}
      grep CALLABLE !{params.sample}_baits_callstate.bed|cut -f 1-3|grep -E "chr[0-9]+|chrX|chrY" > !{params.sample}_baits_callable.bed
    '''
}

process PURECN {
    publishDir "${params.outdir}/${params.sample}/12_tp_tmb_tnb/purecn", pattern: "*.{csv,pdf}", mode:'symlink'
    publishDir "${params.outdir}/${params.sample}/00_log", mode:'symlink',saveAs: {filename -> if (filename =~ /.*.log$/) "${params.sample}_PURECN.log"}
    tag "${params.sample}"
    //echo true
    //conda "/home/dxk/miniconda3/envs/purecn"

    input:
      tuple path(cnvkit_seg), path(cnvkit_cnr), path(bam_readcount_vcf), path(bam_readcount_vcf_tbi), val(sample_id_n), path(vcf_n), path(callable_bed)
      //tuple val(sample_id_t), path(dedup_bam_t), path(recal_t)

    output:
      path("*.{csv,pdf}")
      path(".command.{log,err,out}")

    shell:
      //exome_option_purecn = params.exome_bed ? "${params.exome_bed}" : "${callable_bed}"
    '''
      SAMPLE1=$(zcat !{bam_readcount_vcf}|grep "#CHROM"|cut -f 10)
      SAMPLE2=$(zcat !{bam_readcount_vcf}|grep "#CHROM"|cut -f 11)

      bcftools view -s ${SAMPLE1} !{bam_readcount_vcf} > !{params.sample}_tumor_DNA_AF!{params.TMB_AF}.vcf

      zcat !{vcf_n[0]}|grep -E -v "^chr[0-9a-zA-Z]+_|^HLA|^chrM|chrUn"|sed "s/${SAMPLE2}/${SAMPLE1}/" > !{sample_id_n}_rename.vcf

      cat !{params.sample}_tumor_DNA_AF!{params.TMB_AF}.vcf !{sample_id_n}_rename.vcf|grep -v "#"| \
      sort -k1,1V -k2,2n > !{params.sample}_tumor_DNA_AF!{params.TMB_AF}_normal.vcf

      cat <(grep "#" !{params.sample}_tumor_DNA_AF!{params.TMB_AF}.vcf) !{params.sample}_tumor_DNA_AF!{params.TMB_AF}_normal.vcf \
      > !{params.sample}_tumor_DNA_AF!{params.TMB_AF}_normal_h.vcf

      vcf=!{params.sample}_tumor_DNA_AF!{params.TMB_AF}_normal_h.vcf
      Rscript !{params.purecn_path}/PureCN.R --out ./ \
      --sampleid !{params.sample}_tumor_DNA \
      --tumor !{params.sample}_tumor_DNA_sorted_deduped.cnr \
      --seg-file !{params.sample}_enumerate_chroms.seg \
      --vcf $vcf \
      --genome hg38 \
      --fun-segmentation Hclust \
      --force --post-optimize --seed 123 --parallel --min-af !{params.TMB_AF}

      Rscript !{params.purecn_path}/Dx.R --out !{params.sample} \
      --rds !{params.sample}_tumor_DNA.rds \
      --callable !{callable_bed} \
      --exclude !{params.simpleRepeat} \
      --signatures --force
    '''
}

process TP_TMB_TNB {
    publishDir "${params.outdir}/${params.sample}/12_tp_tmb_tnb/tp_tmb_tnb", pattern: "*_TP_TMB_TNB.{txt,xlsx}", mode:'copy'
    publishDir "${params.outdir}/${params.sample}/00_log", mode:'symlink',saveAs: {filename -> if (filename =~ /.*.log$/) "${params.sample}_TP_TMB_TNB.log"}
    tag "${params.sample}"
    //echo true
    //conda "/home/dxk/miniconda3/envs/purecn"

    input:
      path(neoantigen_txt_callable_bed)
      //tuple val(sample_id_t), path(dedup_bam_t), path(recal_t)

    output:
      path("*_TNB.{txt,xlsx}")
      path(".command.{log,err,out}")

    shell:
    '''
      peptide_postprocessing_for_pvacseq.py --in_file ./!{params.sample}_pvacseq_prime_netmhcstabpan*.txt \
      --sample_name !{params.sample} --output_dir ./  --tumor_dna_vaf !{params.TNB_AF}

      output_tp_tmb_tnb.py -p !{params.sample}_tumor_DNA.csv -m *_mutation_burden.csv -n *_peptide_postprocessing_AF!{params.TNB_AF}.txt -s !{params.sample} -o ./
    '''
}

process SEQUENZA_GC_WIGGLE {
    //publishDir "${params.outdir}/${params.sample}/sequenza", mode:'symlink'
    //publishDir "${params.outdir}/${params.sample}/00_log", mode:'symlink',saveAs: {filename -> if (filename =~ /.*.log$/) "${params.sample}_SEQUENZA_BAM2SEQZ_CHR.log"}
    //tag "${params.sample}"
    echo true
    //conda "/home/dxk/miniconda3/envs/sequenza"

    script:
    if( file(params.hg38_gc_wig) ) {
      """
      echo 'hg38_gc_wig file for sequenza exists'
      """
    }
    else {
      """
      echo 'hg38_gc_wig file for sequenza building'
      dir=\$(dirname ${params.ref})
      sequenza-utils gc_wiggle -w 50 --fasta ${params.ref} -o \$dir/sequenza_hg38_gc50bp.wig.gz
      """
    }
}

process SEQUENZA_BAM2SEQZ_CHR {
    publishDir "${params.outdir}/${params.sample}/13_hla_loh/sequenza", pattern: "*.txt", mode:'symlink'
    publishDir "${params.outdir}/${params.sample}/00_log", mode:'symlink',saveAs: {filename -> if (filename =~ /.*.log$/) "${params.sample}_SEQUENZA_BAM2SEQZ_CHR.log"}
    tag "${params.sample}"
    echo true
    //conda "/home/dxk/miniconda3/envs/sequenza"

    input:
      tuple val(sample_id_n), path(dedup_bam_n), path(recal_n), val(sample_id_t), path(dedup_bam_t), path(recal_t)
      each chr
      //tuple val(sample_id_t), path(dedup_bam_t), path(recal_t)

    output:
      path("*seqz.gz")
      path("*.txt")
      path(".command.{log,err,out}")

    script:
    if( file(params.hg38_gc_wig) ) {
      """
      sequenza-utils bam2seqz -n ${dedup_bam_n[0]} -t ${dedup_bam_t[0]} --fasta ${params.ref} -gc ${params.hg38_gc_wig} -C chr${chr} -o ${params.sample}_chr${chr}.seqz.gz
      echo "OK" > ${params.sample}_bam2seqz_chr.txt
      """
    }
    else {
      """
      dir=\$(dirname ${params.ref})
      sequenza-utils bam2seqz -n ${dedup_bam_n[0]} -t ${dedup_bam_t[0]} --fasta ${params.ref} -gc \$dir/sequenza_hg38_gc50bp.wig.gz -C chr${chr} -o ${params.sample}_chr${chr}.seqz.gz
      echo "OK" ${params.sample}_bam2seqz_chr.txt
      """
    }
}

process SEQUENZA {
    publishDir "${params.outdir}/${params.sample}/13_hla_loh/sequenza", pattern: "*.txt", mode:'copy'
    publishDir "${params.outdir}/${params.sample}/00_log", mode:'symlink',saveAs: {filename -> if (filename =~ /.*.log$/) "${params.sample}_SEQUENZA.log"}
    tag "${params.sample}"
    echo true
    //conda "/home/dxk/miniconda3/envs/sequenza"

    input:
      //tuple val(sample_id_n), path(dedup_bam_n), path(recal_n), val(sample_id_t), path(dedup_bam_t), path(recal_t)
      path(seqz)
      //tuple val(sample_id_t), path(dedup_bam_t), path(recal_t)

    output:
      path("*.txt")
      path(".command.{log,err,out}")

    script:
      """
      cat <(zcat ${params.sample}_chr1.seqz.gz|grep "chromosome") <(zcat ${params.sample}_chr*.seqz.gz|grep -v "chromosome"|sort -k1,1V -k2,2n) > ${params.sample}.seqz
      bgzip ${params.sample}.seqz
      tabix -f -s 1 -b 2 -e 2 -S 1 ${params.sample}.seqz.gz
      sequenza-utils seqz_binning --seqz ${params.sample}.seqz.gz -w 50 -o ${params.sample}_small.seqz.gz
      tabix -f -s 1 -b 2 -e 2 -S 1 ${params.sample}_small.seqz.gz
      Rscript ${baseDir}/bin/sequenza_r_part.R ${params.sample}_small.seqz.gz ${params.sample} ./
      """
}

process POLYSOLVER {
    publishDir "${params.outdir}/${params.sample}/13_hla_loh/polysolver", mode:'copy'
    publishDir "${params.outdir}/${params.sample}/00_log", mode:'symlink',saveAs: {filename -> if (filename =~ /.*.log$/) "${params.sample}_POLYSOLVER.log"}
    tag "${params.sample}"
    //echo true
    //conda "/home/dxk/miniconda3/envs/sequenza"

    input:
      tuple val(sample_id_n), path(dedup_bam_n), path(recal_n), val(sample_id_t), path(dedup_bam_t), path(recal_t)
      //tuple val(sample_id_t), path(dedup_bam_t), path(recal_t)

    output:
      path("*.{txt,out,intervals,vcf}")
      path(".command.{log,err,out}")

    shell:
    '''
      bam=!{params.sample}_normal_DNA_sorted_deduped.bam
      race=Unknown
      includeFreq=1
      build=hg38
      format=STDFQ
      insertCalc=0
      outDir=./
      singularity run !{params.poly_image} \
      /home/polysolver/scripts/shell_call_hla_type \
      $bam $race $includeFreq $build $format $insertCalc $outDir

      normal_bam=!{params.sample}_normal_DNA_sorted_deduped.bam
      tumor_bam=!{params.sample}_tumor_DNA_sorted_deduped.bam
      hla=./winners.hla.txt
      outDir=./
      indiv=!{params.sample}
      singularity run !{params.poly_image} \
      /home/polysolver/scripts/shell_call_hla_mutations_from_type \
      $normal_bam $tumor_bam $hla $build $format $outDir $indiv
    '''
}

process DASH {
    publishDir "${params.outdir}/${params.sample}/13_hla_loh/dash", mode:'copy'
    publishDir "${params.outdir}/${params.sample}/00_log", mode:'symlink',saveAs: {filename -> if (filename =~ /.*.log$/) "${params.sample}_DASH.log"}
    tag "${params.sample}"
    //echo true
    label "l_cpus"
    //conda "/home/dxk/miniconda3/envs/dash"
    beforeScript "export SENTIEON_LICENSE=${params.sentieon_license}"

    input:
      tuple val(sample_id_n), path(reads_n), val(sample_id_t), path(reads_t)
      path(sequenza_outfile)
      path(polysolver_outfile)
      //tuple val(sample_id_t), path(dedup_bam_t), path(recal_t)

    output:
      path("DASH/dash_result*")
      path(".command.{log,err,out}")

    shell:
    '''
      purity=$(cat !{params.sample}_purity_ploidy.txt|awk '{print $1}')
      ploidy=$(cat !{params.sample}_purity_ploidy.txt|awk '{print $2}')

      hla=$(cat winners.hla.txt|cut -f 2-3|xargs|tr ' ' ',')

      normal_read=$(zcat !{reads_n[0]}|awk '{s++}END{print s/2}')
      tumor_read=$(zcat !{reads_t[0]}|awk '{s++}END{print s/2}')
      echo "normal_read: $normal_read; tumor_read: $tumor_read"

      python !{params.DASH_dir}/DASH_manuscript_PE_chech_BAM_sentieon.py \
      --purity $purity --ploidy $ploidy --hla_types $hla \
      --normal_fastq_1 !{reads_n[0]} \
      --normal_fastq_2 !{reads_n[1]} \
      --tumor_fastq_1 !{reads_t[0]} \
      --tumor_fastq_2 !{reads_t[1]} \
      --normal_read_count $normal_read --tumor_read_count $tumor_read \
      --hla_somatic_mutations !{params.DASH_dir}/test_data/hla_mutations.merged.vcf \
      --b_allele_flanking_loh 1,1,1 \
      --all_allele_reference !{params.DASH_dir}/test_data/hla_gen.fasta.IMGTv312.txt \
      --model_filename !{params.DASH_dir}/test_data/training.xgboost_model.2021_05_10.p \
      --threads !{task.cpus} \
      --out_dir ./DASH --sentieon_dir !{params.sentieon_dir}
    '''
}

process READS_STAT {
    publishDir "${params.outdir}/${params.sample}/14_reads_stat", mode:'copy'
    publishDir "${params.outdir}/${params.sample}/00_log", mode:'copy',saveAs: {filename -> if (filename =~ /.*.log$/) "${params.sample}_READS_STAT.log"}
    tag "${params.sample}"
    label "m_cpus"
    //echo true
    //conda "/home/dxk/miniconda3/envs/mosdepth"

    input:
      tuple val(sample_id_n), path(dedup_bam_n), path(recal_n), val(sample_id_t), path(dedup_bam_t), path(recal_t), val(sample_id), path(rna_bam), path(callable_bed)

    output:
      path("*txt")
      path(".command.{log,err,out}")

    script:
      def exome_option_reads_stat = params.exome_bed ? "${params.exome_bed}" : "${callable_bed}"
      """
        echo "\n${sample_id_n}" > ${sample_id_n}_reads_stat.txt
        echo "\n${sample_id_t}" > ${sample_id_t}_reads_stat.txt
        echo "\n${sample_id}" > ${sample_id}_reads_stat.txt

        reads_stat.py ${dedup_bam_n[0]} ${task.cpus} >> ${sample_id_n}_reads_stat.txt
        reads_stat.py ${dedup_bam_t[0]} ${task.cpus} >> ${sample_id_t}_reads_stat.txt
        reads_stat.py ${rna_bam[0]} ${task.cpus} >> ${sample_id}_reads_stat.txt

        bed=${exome_option_reads_stat}
        cat \$bed|awk '{SUM += \$3-\$2} END {print "WES bed size: "SUM}' >> ${params.sample}_reads_stat.txt
        num=\$(cat \$bed|awk '{SUM += \$3-\$2} END {print SUM}')

        mosdepth --by \$bed ${sample_id_n} ${dedup_bam_n[0]} -t ${task.cpus}
        mosdepth --by \$bed ${sample_id_t} ${dedup_bam_t[0]} -t ${task.cpus}
        zcat ${sample_id_n}.regions.bed.gz|awk -v bed_size=\$num '{SUM += (\$3-\$2)*\$4}END{print "Mean depth: "SUM/bed_size}' >> ${sample_id_n}_reads_stat.txt
        zcat ${sample_id_t}.regions.bed.gz|awk -v bed_size=\$num  '{SUM += (\$3-\$2)*\$4}END{print "Mean depth: "SUM/bed_size}' >> ${sample_id_t}_reads_stat.txt


        samtools depth -b \$bed ${dedup_bam_n[0]} -@ ${task.cpus} > ${sample_id_n}_WES_bed_base_depth
        awk -v bed_size=\$num '{if(\$3!=0){sum += 1}}END{print "Coverage: "sum/bed_size}' ${sample_id_n}_WES_bed_base_depth \
        >> ${sample_id_n}_reads_stat.txt
        samtools depth -b \$bed ${dedup_bam_t[0]} -@ ${task.cpus} > ${sample_id_t}_WES_bed_base_depth
        awk -v bed_size=\$num '{if(\$3!=0){sum += 1}}END{print "Coverage: "sum/bed_size}' ${sample_id_t}_WES_bed_base_depth \
        >> ${sample_id_t}_reads_stat.txt

        total_base=\$(zcat ${sample_id_n}.per-base.bed.gz| awk '{sum+=(\$3-\$2)*\$4}END{print sum}')
        target_base=\$(zcat ${sample_id_n}.regions.bed.gz| awk '{sum+=(\$3-\$2)*\$4}END{print sum}')
        awk -v total=\$total_base -v target=\$target_base 'BEGIN{print "On target base(%): "target/total*100}' >> ${sample_id_n}_reads_stat.txt

        total_base=\$(zcat ${sample_id_t}.per-base.bed.gz| awk '{sum+=(\$3-\$2)*\$4}END{print sum}')
        target_base=\$(zcat ${sample_id_t}.regions.bed.gz| awk '{sum+=(\$3-\$2)*\$4}END{print sum}')
        awk -v total=\$total_base -v target=\$target_base 'BEGIN{print "On target base(%): "target/total*100}' >> ${sample_id_t}_reads_stat.txt

        echo "${sample_id_n}" > ${params.sample}_HLA_ABC_cov.txt
        samtools coverage -r chr6:29941260-29945884 ${dedup_bam_n[0]} >> ${params.sample}_HLA_ABC_cov.txt
        samtools coverage -r chr6:31353872-31367067 ${dedup_bam_n[0]} >> ${params.sample}_HLA_ABC_cov.txt
        samtools coverage -r chr6:31268749-31272130 ${dedup_bam_n[0]} >> ${params.sample}_HLA_ABC_cov.txt

        echo "\n${sample_id_t}" >> ${params.sample}_HLA_ABC_cov.txt
        samtools coverage -r chr6:29941260-29945884 ${dedup_bam_t[0]} >> ${params.sample}_HLA_ABC_cov.txt
        samtools coverage -r chr6:31353872-31367067 ${dedup_bam_t[0]} >> ${params.sample}_HLA_ABC_cov.txt
        samtools coverage -r chr6:31268749-31272130 ${dedup_bam_t[0]} >> ${params.sample}_HLA_ABC_cov.txt

        cat ${sample_id_n}_reads_stat.txt ${sample_id_t}_reads_stat.txt ${sample_id}_reads_stat.txt >> ${params.sample}_reads_stat.txt
      """
}

process REPORT_FORMAT {
    publishDir "${params.outdir}/${params.sample}/15_report_format", mode:'copy'
    publishDir "${params.outdir}/${params.sample}/00_log", mode:'symlink',saveAs: {filename -> if (filename =~ /.*.log$/) "${params.sample}_REPORT_FORMAT.log"}
    tag "${params.sample}"

    input:
      path(tp_tmb_tnb)
      path(dash)
      path(peptide)

    output:
      path("*")
      path(".command.{log,err,out}")

    shell:
    '''
      echo -e "project\tStemirna_project\npatient\t!{params.sample}" > project.txt

      report.py \
      -t !{params.sample}_TP_TMB_TNB.txt \
      -a ./dash_result_*/DASH.output.txt \
      -n pvacseq_!{params.sample}_flank!{params.peptide_flank_length}aa_peptide_HLA_I_AF0.1.txt \
      -p project.txt -c !{baseDir}/bin/report_cache \
      -o ./HLA_I_AF0.1

      report.py \
      -t !{params.sample}_TP_TMB_TNB.txt \
      -a ./dash_result_*/DASH.output.txt \
      -n pvacseq_!{params.sample}_flank!{params.peptide_flank_length}aa_peptide_HLA_I_AF0.02.txt \
      -p project.txt -c !{baseDir}/bin/report_cache \
      -o ./HLA_I_AF0.02
    '''
}

process FARDEEP {
    publishDir "${params.outdir}/${params.sample}/fardeep", mode:'copy'
    publishDir "${params.outdir}/${params.sample}/00_log", mode:'symlink',saveAs: {filename -> if (filename =~ /.*.log$/) "${params.sample}_FARDEEP.log"}
    tag "${params.sample}"
    //echo true
    //conda "/home/dxk/miniconda3/envs/dash"

    input:
      path(stringtie_outfile)
      //tuple val(sample_id_t), path(dedup_bam_t), path(recal_t)

    output:
      path("*txt")
      path(".command.{log,err,out}")

    shell:
    '''
      file=!{params.sample}_tumor_RNA_stringtie_gene_expression.tsv
      cat $file|awk -v OFS='\t' '$2!~/^-|^\\./{print $2,$NF}'|sed 's/ID\tTPM/Gene\tSample/'|\
      awk -v OFS='\t' '!a[$1]++{print}' > !{params.sample}_fardeep_input.txt

      Rscript !{baseDir}/bin/fardeep_TIL10_LM22_single_sample.R \
      !{params.sample}_fardeep_input.txt !{params.sample}_fardeep_output.txt
    '''
}

//read_pairs_ch = channel.fromFilePairs(params.reads, checkIfExists: true)
//e1_length_ch = channel.of(8, 9, 10, 11)

workflow {
    read_pairs_ch = channel.fromPath(params.sample_list).splitCsv(header:true).map { row-> tuple(row.sample_id, [file(row.read1), file(row.read2)]) }


if(params.run_snv_indel){
    //FASTQC(read_pairs_ch)
    BWA_INDEX()

    if(params.trim_reads){
      FASTP(read_pairs_ch)
      SENTIEON_MAP(FASTP.out[0])
    }else{
      SENTIEON_MAP(read_pairs_ch)
    }

    SENTIEON_METRIC(SENTIEON_MAP.out[0])
    SENTIEON_DEDUP(SENTIEON_MAP.out[0])
    SENTIEON_BQSR(SENTIEON_DEDUP.out[0])

    SENTIEON_DEDUP_ch = SENTIEON_DEDUP.out[0]
    //SENTIEON_DEDUP_REALIGN_ch.view()
    //SENTIEON_DEDUP_REALIGN_ch.flatten().view()
    //SENTIEON_BQSR.out[0].view()
    //SENTIEON_BQSR.out[1].collect().view()

    SENTIEON_BQSR_ch = SENTIEON_BQSR.out[1].map({file -> [file.getName().split('_recal')[0],file]})
    //SENTIEON_BQSR_ch.view()
    //SENTIEON_DEDUP_BQSR_ch = SENTIEON_DEDUP_REALIGN_ch.mix(SENTIEON_BQSR_ch).groupTuple()
    SENTIEON_DEDUP_BQSR_ch = SENTIEON_DEDUP_ch.join(SENTIEON_BQSR_ch)
    //SENTIEON_DEDUP_BQSR_ch.view()
    SENTIEON_HAPLOTYPER(SENTIEON_DEDUP_BQSR_ch)

    //SENTIEON_DEDUP_REALIGN_ch.view()
    //SENTIEON_BQSR_ch.view()
    //SENTIEON_DEDUP_REALIGN_ch.concat(SENTIEON_BQSR_ch).view()
    //SENTIEON_DEDUP_REALIGN_ch.join(SENTIEON_BQSR_ch).view()
    //SENTIEON_DEDUP_REALIGN_ch.concat(SENTIEON_BQSR_ch).collect().view()

    SENTIEON_DEDUP_BQSR_ch
       .branch {
        normal_br: it =~ /(?i)normal|-n/
        tumor_br: it =~ /(?i)tumor|-t/
    }.set {PAIR_ch}
    PAIR_nt_ch = PAIR_ch.normal_br.concat(PAIR_ch.tumor_br).collect()
    SENTIEON_TNSCOPE(PAIR_nt_ch)

    PHASING_ch = SENTIEON_HAPLOTYPER.out[0].concat(SENTIEON_TNSCOPE.out[0],PAIR_ch.tumor_br).collect()
    //PHASING_ch.view()
    GATK_PHASING(PHASING_ch)

    VEP_ch = GATK_PHASING.out[0].concat(GATK_PHASING.out[1]).collect()
    VEP(VEP_ch)

    HISAT2_INDEX()
    if(params.trim_reads){
      HISAT2_MAP(FASTP.out[0])
    }else{
      HISAT2_MAP(read_pairs_ch)
    }

    //STRINGTIE(HISAT2_MAP.out[0])

    if(params.trim_reads){
      KALLISTO(FASTP.out[0])
    }else{
      KALLISTO(read_pairs_ch)
    }

    COVERAGE_ch = VEP.out[1].concat(PAIR_ch.normal_br,PAIR_ch.tumor_br,HISAT2_MAP.out[0]).collect()
    //COVERAGE_ch.view()
    COVERAGE(COVERAGE_ch)

    //EXPRESSION_ch = COVERAGE.out[0].concat(STRINGTIE.out[0],KALLISTO.out[0]).collect()
    EXPRESSION_ch = COVERAGE.out[0].concat(KALLISTO.out[0]).collect()
    EXPRESSION(EXPRESSION_ch)
    //EXPRESSION.out[0].view()
    //VEP.out[0].view()


    if(params.hla_i){
      if(params.trim_reads){
        OPTITYPE(FASTP.out[0])
      }else{
        OPTITYPE(read_pairs_ch)
      }
    }

    if(params.hla_ii){
      HLA_LA(SENTIEON_DEDUP_ch)
    }

    if(params.hla_i_ii){
      if(params.trim_reads){
        OPTITYPE(FASTP.out[0])
      }else{
        OPTITYPE(read_pairs_ch)
      }
      HLA_LA(SENTIEON_DEDUP_ch)
      HLA_I_II_ch = OPTITYPE.out[0].concat(HLA_LA.out[0])
      HLA_I_II(HLA_I_II_ch)
    }


    if(params.hla_i){
      if(params.pvacseq_split_hla){
        PVACSEQ_ch = EXPRESSION.out[0].concat(VEP.out[0]).collect()
        HLA_I_ch = OPTITYPE.out[0].splitCsv().flatten().filter({HLA -> HLA =~/(HLA-A|HLA-B|HLA-C)/})
        PVACSEQ_SPLIT_HLA(PVACSEQ_ch, HLA_I_ch)
        PEPTIDE_POSTPROCESS_HLA_I_ch = PVACSEQ_SPLIT_HLA.out[0].collect()
        //PEPTIDE_POSTPROCESS_HLA_I_ch.view()
        PEPTIDE_POSTPROCESS_HLA_I_PVACSEQ(PEPTIDE_POSTPROCESS_HLA_I_ch)
      }else{
        PVACSEQ_ch = EXPRESSION.out[0].concat(VEP.out[0],OPTITYPE.out[0]).collect()
        PVACSEQ_HLA_I(PVACSEQ_ch)
        PEPTIDE_POSTPROCESS_HLA_I_PVACSEQ(PVACSEQ_HLA_I.out[0])
           }
      }

    if(params.hla_ii){
      if(params.pvacseq_split_hla){
        PVACSEQ_ch = EXPRESSION.out[0].concat(VEP.out[0]).collect()
        HLA_II_ch = HLA_LA.out[0].splitCsv().flatten().flatten().filter({HLA -> HLA =~/(DR|DP|DQ)/})
        PVACSEQ_SPLIT_HLA(PVACSEQ_ch, HLA_II_ch)
      }else{
        PVACSEQ_ch = EXPRESSION.out[0].concat(VEP.out[0],HLA_LA.out[0]).collect()
        PVACSEQ_HLA_II(PVACSEQ_ch)
        PEPTIDE_POSTPROCESS_HLA_II_PVACSEQ(PVACSEQ_HLA_II.out[0])
           }
    }

    if(params.hla_i && params.hla_ii){
       PEPTIDE_POSTPROCESS_HLA_I_II_PVACSEQ(PEPTIDE_POSTPROCESS_HLA_I_PVACSEQ.out[0],PEPTIDE_POSTPROCESS_HLA_II_PVACSEQ.out[0])
    }

    if(params.hla_i_ii){
      if(params.pvacseq_split_hla){
        PVACSEQ_ch = EXPRESSION.out[0].concat(VEP.out[0]).collect()
        HLA_I_II_ch = HLA_I_II.out[0].splitCsv().flatten().filter({HLA -> HLA =~/(HLA-A|HLA-B|HLA-C|DR|DP|DQ)/})
        PVACSEQ_SPLIT_HLA(PVACSEQ_ch, HLA_I_II_ch)
      }else{
        PVACSEQ_ch = EXPRESSION.out[0].concat(VEP.out[0],HLA_I_II.out[0]).collect()
        PVACSEQ(PVACSEQ_ch)
           }
    }
}



if( params.run_indicator){
    CNVKIT(PAIR_nt_ch)
    CALLSTATE_ch = CNVKIT.out[1].concat(PAIR_ch.tumor_br).collect()
    CALLSTATE(CALLSTATE_ch)
    PURECN_ch = CNVKIT.out[0].concat(EXPRESSION.out[0],SENTIEON_HAPLOTYPER.out[0],CALLSTATE.out[0]).collect()
    PURECN(PURECN_ch)

    TP_TMB_TNB_ch = PEPTIDE_POSTPROCESS_HLA_I_PVACSEQ.out[0].concat(PURECN.out[0]).collect()
    TP_TMB_TNB(TP_TMB_TNB_ch)

    PAIR_nt_ch = PAIR_ch.normal_br.concat(PAIR_ch.tumor_br).collect()
    SEQUENZA_GC_WIGGLE()
    chromosome_ch = channel.of(1..22, 'X', 'Y')
    SEQUENZA_BAM2SEQZ_CHR(PAIR_nt_ch,chromosome_ch)
    SEQUENZA(SEQUENZA_BAM2SEQZ_CHR.out[0].collect())
    POLYSOLVER(PAIR_nt_ch)

    if(params.trim_reads){
      FASTP_ch = FASTP.out[0]
      FASTP_ch
         .branch {
          normal_br: it =~ /(?i)normal_DNA/
          tumor_br: it =~ /(?i)tumor_DNA/
      }.set {FASTP_FASTQ_ch}
      FASTP_FASTQ_nt_ch = FASTP_FASTQ_ch.normal_br.concat(FASTP_FASTQ_ch.tumor_br).collect()
      DASH(FASTP_FASTQ_nt_ch,SEQUENZA.out[0],POLYSOLVER.out[0])
    }else{
      read_pairs_ch
               .branch {
          normal_br: it =~ /(?i)normal_DNA/
          tumor_br: it =~ /(?i)tumor_DNA/
      }.set {read_pairs_FASTQ_ch}
      read_pairs_FASTQ_nt_ch = read_pairs_FASTQ_ch.normal_br.concat(read_pairs_FASTQ_ch.tumor_br).collect()
      DASH(read_pairs_FASTQ_nt_ch,SEQUENZA.out[0],POLYSOLVER.out[0])
    }

    READS_STAT_ch = PAIR_nt_ch.concat(HISAT2_MAP.out[0],CALLSTATE.out[0]).collect()
    READS_STAT(READS_STAT_ch)
    //read_pairs_ch.view()
    //FASTP.out[0].collect().view()
    REPORT_FORMAT(TP_TMB_TNB.out[0],DASH.out[0],PEPTIDE_POSTPROCESS_HLA_I_PVACSEQ.out[0])
    //FARDEEP(STRINGTIE.out[0])
}


if(params.run_fusion){
    if(params.trim_reads){
      //STAR_FUSION(FASTP.out[0])
      FASTP(read_pairs_ch)
      EASYFUSE(FASTP.out[0])
    }else{
      //STAR_FUSION(read_pairs_ch)
      EASYFUSE(read_pairs_ch)
    }

    //AGFUSION(STAR_FUSION.out[0])
    AGFUSION(EASYFUSE.out[0])
    //file(AGFUSION.out[0])

    if(params.hla_i){
      PVACFUSE_HLA_I(AGFUSION.out[0],OPTITYPE.out[0])
      PEPTIDE_POSTPROCESS_HLA_I_PVACFUSE(AGFUSION.out[0],PVACFUSE_HLA_I.out[0],EASYFUSE.out[1],EASYFUSE.out[0])
    }

    if(params.hla_ii){
      PVACFUSE_HLA_II(AGFUSION.out[0],HLA_LA.out[0])
      PEPTIDE_POSTPROCESS_HLA_II_PVACFUSE(AGFUSION.out[0],PVACFUSE_HLA_II.out[0])
    }

    if(params.hla_i_ii){
            PVACFUSE(AGFUSION.out[0],HLA_I_II.out[0])
    }
}


/*
dump
    AGFUSION_ch = AGFUSION.out[0]
    .flatten()
    .map({file -> [file.baseName, file]})
    .filter({ fileName, file -> fileName =~/(?i)^[A-Z0-9]/})
    .map({fileName, file -> [file] })
    .collect()
    //.view()

    //PVACFUSE_ch = AGFUSION_ch.concat(HLA_I_II.out[0]).collect()
    //PVACFUSE_ch.view()
    //PVACFUSE(PVACFUSE_ch)
*/
}
