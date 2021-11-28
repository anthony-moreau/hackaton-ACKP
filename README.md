
### Guide Utilisateur : 
#### Étape 1 : Télécharger le dépot github
Se placer dans un répertoire qui contient suffisamment de mémoire, dans notre cas, on se place dans le répertoire mydatalocal grâce à la commande : 

```shell
cd /mnt/mydatalocal
````


Ensuite télécharger tout le répertoire, pour cela, tapez la commande ci dessous :

```shell
git clone git@github.com:anthony-moreau/hackaton-ACKP
```

![image](https://user-images.githubusercontent.com/90893697/143780706-44e62151-e6d6-4b14-ac81-d2612de44491.png)

#### Étape 2 : Installer conda
(https://www.anaconda.com/products/individual)


#### Étape 3 : Installer via conda snakemake et singularity 
```shell
# Activer conda environment
conda activate 

# Installer snakemake 
# vous pouvez l'installer par Conda/Mamba ou pip : 
# https://snakemake.readthedocs.io/en/stable/getting_started/installation.html#
pip install snakemake

# Activer snakemake environment
conda activate snakemake

# Installer singularity sous snakemake environment. Nous utilisons ici la version 3.6.3.
conda install singularity=3.6.3
```

#### Étape 4 : Lancer le pipepeline
- Le sturcture du répertoire :
```
.
├── Group_4_R_result
│   ├── Analysis.Rdata
│   ├── Analysis.html
│   └── Rplots.pdf
├── Group_4_counts_result
│   ├── SRR628582.counts
│   ├── SRR628582.counts.summary
│   ├── SRR628583.counts
│   ├── ...
│   ├── ...
│   └── SRR628589.counts.summary
├── Groupe_4_Rapport
│   └── Rapport_hackathon_Groupe4
├── DESeq2.R
├── README.md
└── Snakefile
```

Veuillez aller sur le site https://htmlpreview.github.io/ pour afficher le ficher "Analysis.html".  
Avant de lancer le pipeline, assurez-vous que vous êtes dans le dossier avec nos fichiers "DESeq2.R" et "Snakefile".


- Lancer le pipeline 

Se placer dans le dossier du dépôt et lancer le pipeline

```shell
# Vous pouvez préciser le nombre de cores après --cores, et le fichier snakefile en utilisant -s <nom_du_fichier> 
# si vous avez d'autre fichiers snakefile dans le même dossier
snakemake --use-singularity --cores all
```

- Le structure des sorties du Snakefile :
```
.
├── Log.final.out
├── Log.out
├── Log.progress.out
├── Log.std.out
├── SJ.out.tab
├── SRA_files
│   ├── SRR628582.sra
│   ├── SRR628583.sra
│   ├── SRR628584.sra
│   ├── SRR628585.sra
│   ├── SRR628586.sra
│   ├── SRR628587.sra
│   ├── SRR628588.sra
│   └── SRR628589.sra
├── annotations
│   ├──human_genome_annotation.chr.gtf
│   └──human_genome_annotation.chr.gtf.gz
├── chr_files
│   ├── 1.fa.gz
│   ├── 2.fa.gz
│   ├── ...
│   ├── ...
│   ├── 22.fa.gz
│   └── MT.fa.gz
├── fastq_files
│   ├── SRR628582_1.fastq
│   ├── SRR628582_2.fastq
│   ├── SRR628583_1.fastq
│   ├── ...
│   ├── ...
│   └── SRR628589_2.fastq
├── mapping
│   ├── SRR628582.bam
│   ├── SRR628582.bam.bai
│   ├── SRR628583.bam
│   ├── ...
│   ├── ...
│   └── SRR628589.bam.bai
├── ref.fa
├── ref_index
│   ├──Genome
│   ├──Log.out
│   ├──SA
│   ├──SAindex
│   ├──chrLength.txt
│   ├──chrName.txt
│   ├──chrNameLength.txt
│   ├──chrStart.txt
│   └──genomeParameters.txt
├── result
│   ├──SRR628582.counts
│   ├──SRR628582.counts.summary
│   ├──SRR628583.counts
│   ├── ...
│   ├── ...
│   └──SRR628589.counts.summary
├── DESeq2.R
└── Snakefile
```
Voir le rapport pour l'analyse des résultats et des plots obtenus. 

