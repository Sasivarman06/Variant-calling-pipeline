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
        "../envs/gatk.yaml"
    shell:
        """
        gatk HaplotypeCaller -R {input.ref} -I {input.bam} -O {output} -A StrandBiasBySample > {log} 2>&1
        """
