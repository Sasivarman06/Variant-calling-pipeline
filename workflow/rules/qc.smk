

rule Fastqc:
    input:
        unpack(get_fastq)
    output:
        "results/fastqc_report/{sample}_1_fastqc.html",
        "results/fastqc_report/{sample}_2_fastqc.html"   
    log:
        "logs/fastqc/{sample}_fastqc.log"
    conda:
        "../envs/fastqc.yaml"
    shell:
        """
        fastqc {input} -o results/fastqc_report/ > {log} 2>&1
        """
