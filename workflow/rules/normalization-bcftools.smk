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
        "../envs/bcftools.yaml"
    shell:
        """
        bcftools norm -f {input.ref} -Oz -o {output.vcf} {input.raw_vcf}  > {log} 2>&1
        tabix -p vcf {output.vcf} >> {log} 2>&1
        """
	
