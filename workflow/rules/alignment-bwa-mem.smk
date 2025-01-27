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
        "../envs/bwa.yaml"
    shell:
        """
        bwa mem -K 10000000 -c 100 -R "{params.read_group}" -M -T 50 {input.ref} {input.r1_paired} {input.r2_paired} > {output}  2>> {log}
        """
   
