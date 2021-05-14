samples = ["A", "B"]


rule all:
    input:
        "calls/all.vcf",
        "plots/quals.svg"


rule bwa:
    input:
        "data/genome.fa",
        "data/samples/{sample}.fastq"
    output:
        temp("mapped/{sample}.bam")
    spack:
        "envs/spack-mapping.yml"
    threads: 8
    shell:
        "bwa mem -t {threads} {input} | samtools view -Sb - > {output}"


rule sort:
    input:
        "mapped/{sample}.bam"
    output:
        "mapped/{sample}.sorted.bam"
    spack:
        "envs/spack-mapping.yml"
    shell:
        "samtools sort -o {output} {input}"



rule call:
    input:
        fa="data/genome.fa",
        bam=expand("mapped/{sample}.sorted.bam", sample=samples)
    output:
        "calls/all.vcf"
    spack:
        "envs/spack-calling.yml"
    shell:
        "samtools mpileup -g -f {input.fa} {input.bam} | "
        "bcftools call -mv - > {output}"

rule stats:
    input:
        "calls/all.vcf"
    output:
        report("plots/quals.svg", caption="report/calling.rst")
    spack:
        "envs/spack-stats.yml"
    script:
        "scripts/plot-quals.py"
