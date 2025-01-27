
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
        "../envs/samtools.yaml"
    shell:
        """
        samtools sort -n -o {output.sorted_bam} {input.sam} 2>> {log}
        samtools fixmate -m {output.sorted_bam} {output.fixmate_bam} 2>> {log}
        samtools sort {output.fixmate_bam} -o {output.sorted_fixmate_bam}  2>> {log}
        samtools markdup {output.sorted_fixmate_bam} {output.dedup_bam} 2>> {log}
        samtools index {output.dedup_bam} 2>> {log}
        """
