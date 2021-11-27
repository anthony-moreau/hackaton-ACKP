Ce dépot contient un pipeline snakemake pour l'analyse de donnée d'expression de l'ADN (RNAseq). Ce pipeline fonctionne sous linux avec la distribution anaconda de python.

Pour le faire tourner sous une machine linux:

- télécharger le dépot github
- installer conda (https://www.anaconda.com/products/individual)
- installer via conda snakemake et singularity dans le même environnement (conda activate puis conda install snakemake/singularity)
- se placer dans le dossier du dépôt et lancer le pipepeline avec la commande snakemake --use-singularity -cores all
