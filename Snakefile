configfile: "config.yaml"

rule all:
    input:
        expand("results/fastqc_report/{sample}_1_fastqc.html", sample=config["samples"]),
        expand("results/fastqc_report/{sample}_2_fastqc.html", sample=config["samples"]),
        expand("results/norm_vcf/{sample}_normalized.vcf.gz", sample=config["samples"]),
        expand("results/norm_vcf/{sample}_normalized.vcf.gz.tbi", sample=config["samples"])

rule Fastqc:
    input:
        "tb_data/{sample}_1.fastq.gz",
        "tb_data/{sample}_2.fastq.gz"
    output:
        "results/fastqc_report/{sample}_1_fastqc.html",
        "results/fastqc_report/{sample}_2_fastqc.html"   
    log:
        "logs/fastqc/{sample}_fastqc.log"
    conda:
        "env.yaml"
    shell:
        """
        fastqc {input} -o results/fastqc_report/ > {log} 2>&1
        """
   

rule Trimming:
    input:
        r1="tb_data/{sample}_1.fastq.gz",
        r2="tb_data/{sample}_2.fastq.gz"
    output:
        r1_paired="results/trimmed_fastq/{sample}_trimmed_1P.fastq.gz",
        r1_unpaired="results/trimmed_fastq/{sample}_trimmed_1U.fastq.gz",
        r2_paired="results/trimmed_fastq/{sample}_trimmed_2P.fastq.gz",
        r2_unpaired="results/trimmed_fastq/{sample}_trimmed_2U.fastq.gz"
    conda:
        "env.yaml"
    log:
        "logs/trimming/{sample}_trimming.log"
    shell:
        """
        trimmomatic PE -phred33 {input.r1} {input.r2} \
        {output.r1_paired} {output.r1_unpaired} {output.r2_paired} {output.r2_unpaired} \
        LEADING:3 TRAILING:3 SLIDINGWINDOW:4:20 MINLEN:36  > {log} 2>&1
        """
  
        
rule Alignment:
    input:
        ref=config["reference"],
        r1_paired="results/trimmed_fastq/{sample}_trimmed_1P.fastq.gz",
        r2_paired="results/trimmed_fastq/{sample}_trimmed_2P.fastq.gz"
    output:
        "results/sam_files/{sample}.sam"
    log:
        "logs/alignment/{sample}_alignment.log"
    params:
        read_group="@RG\\tID:{wildcards.sample}\\tSM:{wildcards.sample}\\tPL:ILLUMINA"
    conda:
        "env.yaml"
    shell:
        """
        bwa mem -K 10000000 -c 100 -R "{params.read_group}" -M -T 50 {input.ref} {input.r1_paired} {input.r2_paired} > {output}  2>> {log}
        """
   

rule Sort_and_markdup:
    input:
        sam="results/sam_files/{sample}.sam"
    output:
        dedup_bam="results/dedup_bam/{sample}_dedup.bam",
        dedup_index="results/dedup_bam/{sample}_dedup.bam.bai",
        sorted_bam=temp("results/sorted_bam/{sample}_sorted.bam"),
        fixmate_bam=temp("results/fixmate_bam/{sample}_fixmate.bam"),
        sorted_fixmate_bam=temp("results/sorted_fixmate_bam/{sample}_sorted_fixmate.bam")
    log:
        "logs/sort_and_markdup/{sample}_sort_and_markdup.log"
    conda:
        "env.yaml"
    shell:
        """
        samtools sort -n -o {output.sorted_bam} {input.sam} 2>> {log}
        samtools fixmate -m {output.sorted_bam} {output.fixmate_bam} 2>> {log}
        samtools sort {output.fixmate_bam} -o {output.sorted_fixmate_bam}  2>> {log}
        samtools markdup {output.sorted_fixmate_bam} {output.dedup_bam} 2>> {log}
        samtools index {output.dedup_bam} 2>> {log}
        """

rule Variant_calling:
    input:
        ref=config["reference"],
        bam="results/dedup_bam/{sample}_dedup.bam",
        bam_index="results/dedup_bam/{sample}_dedup.bam.bai"
    output:
        "results/raw_vcf/{sample}_raw.vcf.gz"
    log:
        "logs/variant_calling/{sample}_variant_calling.log"    
    conda:
        "env.yaml"
    shell:
        """
        gatk HaplotypeCaller -R {input.ref} -I {input.bam} -O {output} -A StrandBiasBySample > {log} 2>&1
        """

rule Normalize_vcf:
    input:
        ref=config["reference"],
        raw_vcf="results/raw_vcf/{sample}_raw.vcf.gz"
    output:
        vcf="results/norm_vcf/{sample}_normalized.vcf.gz",
        tbi="results/norm_vcf/{sample}_normalized.vcf.gz.tbi" 
    log:
        "logs/normalize_vcf/{sample}_normalize_vcf.log"    
    conda:
        "env.yaml"
    shell:
        """
        bcftools norm -f {input.ref} -Oz -o {output.vcf} {input.raw_vcf}  > {log} 2>&1
        tabix -p vcf {output.vcf} >> {log} 2>&1
        """
	
