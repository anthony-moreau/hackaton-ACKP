rule all:
    input: 
        "result/{sample}.counts"

rule feature_count:
    input: 
        "mapping/{sample}.bam",
        "mapping/{sample}.bam.bai",
        "annotation/human_genome_annotation.chr.gtf.gz"
    output: 
        "result/{sample}.counts"
    singularity: 
        "docker://evolbioinfo/subread:v2.0.1"
    shell: 
        "featureCounts -T {threads} -t gene -g gene_id -s 0 -a input.gtf -o {output} mapping/{sample}.bam"

rule uncompress_annotations:
    input:
        "annotation/human_genome_annotation.chr.gtf.gz"
    output:
        "annotation/human_genome_annotation.gtf"
    shell:
        "gunzip -c {input} > {output}"
    

rule get_annotations:
    output: "annotation/human_genome_annotation.chr.gtf.gz"
    shell:
        "wget ftp://ftp.ensembl.org/pub/release-101/gtf/homo_sapiens/Homo_sapiens.GRCh38.101.chr.gtf.gz"

rule samtool_index:
    input: 
        "mapping/{sample}.bam"
    output: 
        "mapping/{sample}.bam.bai"
    singularity: 
        "docker://evolbioinfo/samtools:v1.11" 
    shell:
        "samtools index mapping/{sample}.bam"