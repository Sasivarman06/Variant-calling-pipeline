configfile: "config.yaml"

rule all:
    input:
        expand("results/fastqc_report/{sample}_1_fastqc.html",sample=config["samples"]),
        expand("results/fastqc_report/{sample}_2_fastqc.html",sample=config["samples"]),
        expand("results/norm_vcf/{sample}_normalized.vcf.gz", sample=config["samples"]),
        expand("results/norm_vcf/{sample}_normalized.vcf.gz.tbi", sample=config["samples"])

rule fastqc:
    input:
        "tb_data/{sample}_1.fastq.gz",
        "tb_data/{sample}_2.fastq.gz"
    output:
        "results/fastqc_report/{sample}_1_fastqc.html",
        "results/fastqc_report/{sample}_2_fastqc.html"	
    shell:
        """
        fastqc {input} -o results/fastqc_report/
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
    shell:
        """
        trimmomatic PE -phred33 {input.r1} {input.r2} \
        {output.r1_paired} {output.r1_unpaired} {output.r2_paired} {output.r2_unpaired} \
        LEADING:3 TRAILING:3 SLIDINGWINDOW:4:20 MINLEN:36
        """
        
rule Alignment:
    input:
        ref=config["reference"],
        r1_paired="results/trimmed_fastq/{sample}_trimmed_1P.fastq.gz",
        r2_paired="results/trimmed_fastq/{sample}_trimmed_2P.fastq.gz"
    output:
        "results/sam_files/{sample}.sam"
    params:
        read_group="@RG\\tID:{wildcards.sample}\\tSM:{wildcards.sample}\\tPL:ILLUMINA"
    shell:
        """
        bwa mem -K 10000000 -c 100 -R "{params.read_group}" -M -T 50 {input.ref} {input.r1_paired} {input.r2_paired} > {output}
        """

rule sort_and_markdup:
    input:
        sam="results/sam_files/{sample}.sam"
    output:
        dedup_bam="results/dedup_bam/{sample}_dedup.bam",
        dedup_index="results/dedup_bam/{sample}_dedup.bam.bai",
        sorted_bam=temp("results/sorted_bam/{sample}_sorted.bam"),
        fixmate_bam=temp("results/fixmate_bam/{sample}_fixmate.bam"),
        sorted_fixmate_bam=temp("results/sorted_fixmate_bam/{sample}_sorted_fixmate.bam")
    shell:
        """
        samtools sort -n -o {output.sorted_bam} {input.sam}
        samtools fixmate -m {output.sorted_bam} {output.fixmate_bam}
        samtools sort {output.fixmate_bam} -o {output.sorted_fixmate_bam}
        samtools markdup {output.sorted_fixmate_bam} {output.dedup_bam}
        samtools index {output.dedup_bam}
        """

rule Variant_calling:
    input:
        ref=config["reference"],
        bam="results/dedup_bam/{sample}_dedup.bam",
        bam_index="results/dedup_bam/{sample}_dedup.bam.bai"
    output:
        "results/raw_vcf/{sample}_raw.vcf.gz"
    shell:
        """
        gatk HaplotypeCaller -R {input.ref} -I {input.bam} -O {output} -A StrandBiasBySample
        """

rule normalize_vcf:
    input:
        ref=config["reference"],
        raw_vcf="results/raw_vcf/{sample}_raw.vcf.gz"
    output:
        vcf="results/norm_vcf/{sample}_normalized.vcf.gz",
        tbi="results/norm_vcf/{sample}_normalized.vcf.gz.tbi"	
    shell:
        """
        bcftools norm -f {input.ref} -Oz -o {output.vcf} {input.raw_vcf}
        tabix -p vcf {output.vcf}
        """

