import pandas as pd
from snakemake.utils import validate

configfile: "config/config.yaml"
samples = pd.read_csv(config["samples"]).set_index("sample", drop=False)
print(samples)
validate(samples, schema="../schemas/config.schema.yaml")


wildcard_constraints:
    sample="|".join(samples.index)

def get_fastq(wildcards):
    fastqs = samples.loc[(wildcards.sample), ["fq1","fq2"]].dropna()
    if len(fastqs) == 2:
       return {"r1": fastqs.fq1, "r2": fastqs.fq2}
    return {"r1": fastqs.fq1}


