chromosomes = ["1", "2","3","4","5","6","7","8","9","10","11","12","13","14","15","16","17","18","19","20","21","22","Mt"]
accessions_numbers = ["SRR628582", "SRR628583", "SRR628584", "SRR628585", "SRR628586", "SRR628587", "SRR628588", "SRR628589"]

rule all:
    input:
        expand("result/{sample}.counts", sample=accessions_numbers)

rule feature_count:
    input: 
        sample_mapping = "mapping/{sample}.bam",
        genome_annotation = "annotation/human_genome_annotation.gtf",
        index = "mapping/{sample}.bam.bai"
    output: 
        "result/{sample}.counts"
    singularity: 
        "docker://evolbioinfo/subread:v2.0.1"
    shell: 
        "featureCounts -T {threads} -t gene -g gene_id -s 0 -a {input.genome_annotation} -o {output} {input.sample_mapping}"

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

rule samtool_index_mapping:
    input: 
        "mapping/{sample}.bam"
    output: 
        "mapping/{sample}.bam.bai"
    singularity: 
        "docker://evolbioinfo/samtools:v1.11"
    shell:
        "samtools index mapping/{sample}.bam"

rule star_mapping:
    input:
        fastq1 = "sequences/{sample}_1.fastq",
        fastq2 = "sequences/{sample}_2.fastq",
        genome_directory = "ref"
    output: 
        "mapping/{sample}.bam"
    singularity:
        "docker://evolbioinfo/star:v2.7.6a"
    shell:
        """
        mkdir mapping
        STAR --outSAMstrandField intronMotif \
        --outFilterMismatchNmax 4 \
        --outFilterMultimapNmax 10 \
        --genomeDir {input.genome_directory} \
        --readFilesIn <(gunzip -c {input.fastq1}) <(gunzip -c {input.fastq2}) \
        --runThreadN {threads} \
        --outSAMunmapped None \
        --outSAMtype BAM SortedByCoordinate \
        --outStd BAM_SortedByCoordinate \
        --genomeLoad NoSharedMemory \
        --limitBAMsortRAM mem \
        > mapping/{sample}.bam
        """

rule star_genome_index:
    input: 
        "sequences/ref.fa"
    output: 
        "ref"
    singularity:
        "docker://evolbioinfo/star:v2.7.6a" 
    shell:
        """
        mkdir ref
        STAR --runThreadN {threads} --runMode genomeGenerate --genomeDir ref/ --genomeFastaFiles {input}
        """

rule create_reference_file:
    input: 
        expand("sequences/{sample}.fa.gz", sample=chromosomes)
    output: 
        "sequences/ref.fa"
    shell:
        "gunzip -c *.fa.gz > ref.fa"

rule download_chromosome_sequence:
    output: 
        "sequences/{chr}.fa.gz"
    shell:
        "wget -o {output} ftp://ftp.ensembl.org/pub/release-101/fasta/homo_sapiens/dna/Homo_sapiens.GRCh38.dna.chromosome.!{input}.fa.gz"

rule download_sample_fastq:
    output: 
        "sequences/{sample}_1.fastq",
        "sequences/{sample}_2.fastq"
    singularity:
        "docker://evolbioinfo/sratoolkit:v2.10.8"
    threads: 8
    shell:
        "faster-qdump {sample} --split-files --include-technical -e {threads} -O sequences"
