##############################################################################################
# Fichier DESeq2 permettant l'obtention de : 
#     - un fichier Analysis.Rdata dans le dossier analysis_R
#     - un fichier Rplots.pdf avec les graphes que l'on trace dans ce script R dans le même 
##############################################################################################


## téléchargement des données 

SRR628582 <- read.table("result/SRR628582.counts", head=TRUE)
SRR628583 <- read.table("result/SRR628583.counts", head=TRUE)
SRR628584 <- read.table("result/SRR628584.counts", head=TRUE)
SRR628585 <- read.table("result/SRR628585.counts", head=TRUE)
SRR628586 <- read.table("result/SRR628586.counts", head=TRUE)
SRR628587 <- read.table("result/SRR628587.counts", head=TRUE)
SRR628588 <- read.table("result/SRR628588.counts", head=TRUE)
SRR628589 <- read.table("result/SRR628589.counts", head=TRUE)

###### installation des libraries

# ATTENTION pour installer DESeq2, IL FAUT LA VERSION 4.1 de R d'après internet
# si DESeq2 ne veut toujours pas s'installer, copier coller cela sur votre Console :

# if (!requireNamespace("BiocManager", quietly = TRUE))
#   install.packages("BiocManager")
# BiocManager::install("DESeq2")

library("DESeq2")

### Création des jeu de données et renommer les colonnes 

# On créée dans la variable `data` le jeu de donnée, avec le nom du GeneId en 1ère colonne et les counts des 8 SRR par la suite. 

data <- data.frame(SRR628582$Geneid,SRR628582$mapping.SRR628582.bam,SRR628583$mapping.SRR628583.bam,SRR628584$mapping.SRR628584.bam,SRR628585$mapping.SRR628585.bam,SRR628586$mapping.SRR628586.bam,SRR628587$mapping.SRR628587.bam,SRR628588$mapping.SRR628588.bam,SRR628589$mapping.SRR628589.bam)
colnames(data)[1] <- "Geneid"
colnames(data)[2] <- "SRR628582"
colnames(data)[3] <- "SRR628583"
colnames(data)[4] <- "SRR628584"
colnames(data)[5] <- "SRR628585"
colnames(data)[6] <- "SRR628586"
colnames(data)[7] <- "SRR628587"
colnames(data)[8] <- "SRR628588"
colnames(data)[9] <- "SRR628589"

#On crée les métadonnées: récuperation au format data.frame des labels :
  
samples = c("SRR628582", "SRR628583", "SRR628584", "SRR628585", "SRR628586", "SRR628587", "SRR628588", "SRR628589")
Type_M_WT = c("M", "M", "WT", "WT", "WT", "WT", "WT", "M") # M: mutants, WT: Wild-Type
MetaData = data.frame(samples, Type_M_WT)


#On visualise les premières lignes de la variable `data`

head(data)


#On visualise les données dans la variable `MetaData`

print(MetaData)


### Construction de l'objet DESEQDataSet

dds <- DESeqDataSetFromMatrix(countData=data, colData=MetaData, design=~Type_M_WT, tidy = TRUE)

#Regardons à quoi ressemble le `dds` : 
  
print(dds)

### Lançons la fonction DESeq et Regardons la table des résultats 

dds <- DESeq(dds)
res <- results(dds)
dds$Type_M_WT <- relevel(dds$Type_M_WT, ref ="WT") # on modifie pour que la reference soit le WT et pas le M (qui est pris automatiquement car alphabétiquement il est plus proche)
levels(dds$Type_M_WT) # on vérifie que le WT a bien été pris comme référence 
head(results(dds, tidy=TRUE))

###  Résumé de l'expression différentiel des gènes 

summary(res) 

#On souhaite maintenant comparer avec l'article 2. Pour cela on modifie les paramètres, on prend pvalue = 0,05 et LFC=1,5. 

res2 <- results(dds,alpha = 0.05, lfcThreshold = 1.5, pAdjustMethod = 'none') 
summary(res2)

###  Trier la liste récapitulative par p-value ajusté décroissante

#On affiche la liste des 10 gènes avec les padj les plus faibles donc les 10 gènes avec une Expression Différentielle la plus forte. 

res <- res[order(res$padj),]
head(res, 10) # pour afficher la liste des 10 gènes avec les padj les plus faibles 


#Ainsi on a les identifiants des 10 gènes avec une pdjust la plus faibles, essayons de réaliser des graphiques pour mieux visualiser l'expression différentielle de ces gènes. 

### On réalise des plotCounts pour comparer les comptes normalisés entre les mutants et les WildType pour quelques gènes

#Exemple du gène ENSG00000115524 = SF3B1 dont parlent les articles 

par(mfrow=c(1,1))
plotCounts(dds, gene="ENSG00000115524", intgroup="Type_M_WT")


#Exemple des 9 autres gènes les plus exprimés c'est à dire dont la p-value est la plus faible 

par(mfrow=c(3,3))  # par(mfrow=c(2,3)) si je fais 6 graphes, ça les classes en 2 lignes 3 colonnes par exemple 

plotCounts(dds, gene="ENSG00000164659", intgroup="Type_M_WT")
plotCounts(dds, gene="ENSG00000234261", intgroup="Type_M_WT")
plotCounts(dds, gene="ENSG00000153283", intgroup="Type_M_WT")
plotCounts(dds, gene="ENSG00000214212", intgroup="Type_M_WT")
plotCounts(dds, gene="ENSG00000142621", intgroup="Type_M_WT")
plotCounts(dds, gene="ENSG00000186973", intgroup="Type_M_WT")
plotCounts(dds, gene="ENSG00000286122", intgroup="Type_M_WT")
plotCounts(dds, gene="ENSG00000103449", intgroup="Type_M_WT")
plotCounts(dds, gene="ENSG00000007038", intgroup="Type_M_WT")


### Réalisation du Volcano Plot 

#On réalise ici un Volcano Plot avec la librairie DESeq2, toutefois, sur le deuxième fichier R, on a tracé une version améliorée, le EnhancedVolcano Plot qui permet d'avoir les étiquettes des gènes sur les points. 

#reset par
par(mfrow=c(1,2))

# Faire un volcano plot basique 
with(res, plot(log2FoldChange, -log10(pvalue), pch=20, main="Volcano plot", xlim=c(-10,10)))

# Ajouter des points colorés: padj<.01 et LFC = 4
with(subset(res, padj<.01 ), points(log2FoldChange, -log10(pvalue), pch=20, col="blue"))
with(subset(res, padj<.01 & abs(log2FoldChange)>4), points(log2FoldChange, -log10(pvalue), pch=20, col="red"))
abline(v=0, col="black", lty=3, lwd=1.0)
abline(v=-4, col="black", lty=4, lwd=2.0)
abline(v=4, col="black", lty=4, lwd=2.0)
abline(h=-log10(max(res$pvalue[res$padj<.01], na.rm=TRUE)), col="black", lty=4, lwd=2.0)


# Faire un volcano plot basique 
with(res, plot(log2FoldChange, -log10(pvalue), pch=20, main="Volcano plot", xlim=c(-10,10)))

# Ajouter des points colorés: padj<.002 et LFC = 2
with(subset(res, padj<.0002 ), points(log2FoldChange, -log10(pvalue), pch=20, col="blue"))
with(subset(res, padj<.0002 & abs(log2FoldChange)>2), points(log2FoldChange, -log10(pvalue), pch=20, col="red"))
abline(v=0, col="black", lty=3, lwd=1.0)
abline(v=-2, col="black", lty=4, lwd=2.0)
abline(v=2, col="black", lty=4, lwd=2.0)
abline(h=-log10(max(res$pvalue[res$padj<.0002], na.rm=TRUE)), col="black", lty=4, lwd=2.0)


### ACP avec la librairie DESeq2

#Nous devons d'abord transformer les données bruts de comptage: la fonction vst effectuera une transformation stabilisant la variance.

vsdata <- vst(dds, blind=FALSE)
plotPCA(vsdata, intgroup="Type_M_WT")

save.image(file = "analysis_R/Analysis.Rdata")
