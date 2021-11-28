Ce dépot contient un pipeline snakemake pour l'analyse de donnée d'expression de l'ADN (RNAseq). Ce pipeline fonctionne sous linux avec la distribution anaconda de python.

Pour le faire tourner sous une machine linux:

- télécharger le dépot github
- installer conda (https://www.anaconda.com/products/individual)
- installer via conda snakemake et singularity dans le même environnement (conda activate puis conda install snakemake/singularity)
- se placer dans le dossier du dépôt et lancer le pipepeline avec la commande snakemake --use-singularity -cores all


### Usage
#### Étape 1 : Télécharger le dépot github
Vous pouvez télécherger tout le répertoire par des méthodes proposées dans le bouton verte "Code" comme ci-dessous :
![image](https://user-images.githubusercontent.com/90893697/143773298-cee3915a-1b96-4367-9b14-662534adaf6c.png)


#### Étape 2 : Installer conda
(https://www.anaconda.com/products/individual)


#### Étape 3 : Installer via conda snakemake et singularity 
```shell
# Activer conda environment
conda activate 

# Installer snakemake 
# vous pouvez l'installer par Conda/Mamba ou pip : https://snakemake.readthedocs.io/en/stable/getting_started/installation.html#

# Activer snakemake environment
conda activate snakemake

# Installer singularity sous snakemake environment. Nous utilisons ici la version 3.6.3.
conda activate singularity=v3.6.3
```

#### Étape 4 : Lancer le pipepeline
- Le sturcture du répertoire :
```
.
├── our_result
│   ├── Analysis.Rdata
│   ├── Rplots.pdf
│   ├── SRR628582.counts
│   ├── SRR628582.counts.summary
│   ├── SRR628583.counts
│   ├── SRR628583.counts.summary
│   ├── ...
│   ├── ...
│   └── SRR628589.counts.summary
├── Analysis.html 
├── DESeq2.R
├── README.md
└── Snakefile
```

- Lancer le pipeline 
```shell
# Avant de lancer le pipeline, assurez-vous que vous êtes dans le dossier avec DESeq2.R et Snakefile
# 
snakemake --use-singularity -cores all
```


