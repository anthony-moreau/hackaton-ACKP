chromosomes = ["1", "2","3","4","5","6","7","8","9","10","11","12","13","14","15","16","17","18","19","20","21","22","MT"]
accessions_numbers = ["SRR628582", "SRR628583", "SRR628584", "SRR628585", "SRR628586", "SRR628587", "SRR628588", "SRR628589"]

## Règle all :

rule all:
    input:
        "analysis_R/Analysis.Rdata"
        
# telecharger les données SRA depuis le site internet NCBI

rule get_sra :
    output :
        "SRA_files/{sample}.sra"
    threads : 8
    shell :
        " wget -O {output} https://sra-downloadb.be-md.ncbi.nlm.nih.gov/sos1/sra-pub-run-5/{wildcards.sample}/{wildcards.sample}.1"


# convertir chaque fichier SRA en 2 fichiers fastq

rule conversion_sample_fastq:
    input:
        "SRA_files/{sample}.sra"
    output:
        "fastq_files/{sample}_1.fastq",
        "fastq_files/{sample}_2.fastq"
    singularity:
        "docker://evolbioinfo/sratoolkit:v2.5.7"
    threads:4
    shell:
        " fastq-dump --split-files ./{input} -O fastq_files"

# telecharger les chromosomes du génome humain.
# on choisit de télecharger le génome de référence de réference GRCh38-release-101 (et non pas le release-104)

rule download_chromosome_sequence:
    output:
        "chr_files/{chr}.fa.gz"
    threads: 8
    shell:
        """
        wget -O {output} ftp://ftp.ensembl.org/pub/release-101/fasta/homo_sapiens/dna/Homo_sapiens.GRCh38.dna.chromosome.{wildcards.chr}.fa.gz
        """
# Unziper les fichiers pour chaque chromosome de référence téléchargé et les concaténer en un seul fichier pour obtenir le génome humaine de réference : le fichier ref.fa

rule create_reference_file:
    input:
       expand("chr_files/{chr}.fa.gz", chr= chromosomes)
    output:
        "ref.fa"
    shell:
        "gunzip -c {input} > {output}"

# On index le génome humain de réference on obtient en sortie un dossier ref_index avec l'ensemble des informations nécessaires

rule create_genome_index:
    input:
        "ref.fa"
    output :
        "ref_index/SA",
        "ref_index/SAindex",
    threads : 4
    singularity :
        "docker://evolbioinfo/star:v2.7.6a"
    shell :
        """
        STAR --runThreadN {threads} --runMode genomeGenerate --genomeDir ref_index/ --genomeSAindexNbases 6 --genomeFastaFiles {input}
        """

# On annote le génome de réference

rule get_annotations:
    output:
        genome_annot_zip = "annotations/human_genome_annotation.chr.gtf.gz",
        genome_annot = "annotations/human_genome_annotation.chr.gtf"
    shell:
        """
        wget -O {output.genome_annot_zip} ftp://ftp.ensembl.org/pub/release-101/gtf/homo_sapiens/Homo_sapiens.GRCh38.101.chr.gtf.gz
        gunzip -c {output.genome_annot_zip} > {output.genome_annot}
        """

# mapper le génôme pour obtenir un fichier .bam
# on met les fichiers SA et SAindex en input pour forcer le Snakefile à aller jusqu'au bout et pas lancer le mapping avant que l'indexation soit terminée

rule star_mapping:
    input:
        fastq1 = "fastq_files/{sample}_1.fastq",
        fastq2 = "fastq_files/{sample}_2.fastq",
        genome_SA = "ref_index/SA",
        genome_SAindex = "ref_index/SAindex"
    output:
        "mapping/{sample}.bam"
    threads: 16
    singularity:
        "docker://evolbioinfo/star:v2.7.6a"
    shell:
        """
        STAR --outSAMstrandField intronMotif \
        --outFilterMismatchNmax 4 \
        --outFilterMultimapNmax 10 \
        --genomeDir ref_index/ \
        --readFilesIn {input.fastq1} {input.fastq2} \
        --runThreadN {threads} \
        --outSAMunmapped None \
        --outSAMtype BAM SortedByCoordinate \
        --outStd BAM_SortedByCoordinate \
        --genomeLoad NoSharedMemory \
        --limitBAMsortRAM 50000000000 \
        > {output}
        """

# convertir le fichier .bam en un fichier .bam.bai

rule samtool_index_mapping:
    input:
        "mapping/{sample}.bam"
    output:
        "mapping/{sample}.bam.bai"
    singularity:
        "docker://evolbioinfo/samtools:v1.11"
    shell:
        "samtools index mapping/{wildcards.sample}.bam"

# compter

rule feature_count:
    input:
        sample_mapping = "mapping/{sample}.bam",
        genome_annotation = "annotations/human_genome_annotation.chr.gtf",
        index = "mapping/{sample}.bam.bai"
    output:
        "result/{sample}.counts"
    threads : 4
    singularity:
        "docker://evolbioinfo/subread:v2.0.1"
    shell:
        "featureCounts -T {threads} -t gene -g gene_id -s 0 -a {input.genome_annotation} -o {output} {input.sample_mapping}"
        

      
  
# analyse sur R (et renvoyer une markdown de l'analyse)
        
rule analyse_stat_R:
    input:
        expand("result/{sample}.counts", sample=accessions_numbers)
    output:
        "analysis_R/Analysis.Rdata"
    singularity:"docker://evolbioinfo/deseq2:v1.28.1"
    script:
        "DESeq2.R"
