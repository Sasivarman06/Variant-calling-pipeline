
rule Trimming:
    input:
        unpack(get_fastq)
    output:
        r1_paired="results/trimmed_fastq/{sample}_trimmed_1P.fastq.gz",
        r1_unpaired="results/trimmed_fastq/{sample}_trimmed_1U.fastq.gz",
        r2_paired="results/trimmed_fastq/{sample}_trimmed_2P.fastq.gz",
        r2_unpaired="results/trimmed_fastq/{sample}_trimmed_2U.fastq.gz"
    conda:
        "../envs/trim.yaml"
    log:
        "logs/trimming/{sample}_trimming.log"
    shell:
        """
        trimmomatic PE -phred33 {input.r1} {input.r2} \
        {output.r1_paired} {output.r1_unpaired} {output.r2_paired} {output.r2_unpaired} \
        LEADING:3 TRAILING:3 SLIDINGWINDOW:4:20 MINLEN:36  > {log} 2>&1
        """
