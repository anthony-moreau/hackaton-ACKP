# j'ai d'abord téléchargé tout le github hackahton sur mon bureau 

# téléchargement des données (changer le path vers l'endroit ou se trouve le dossier que vous avez téléchargé )

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

###### création du jeu de données avec uniquement les counts à la fin puis renommer les colonnes 

data <- data.frame(SRR628582$Geneid,SRR628582$mapping.SRR628582.bam,SRR628583$mapping.SRR628583.bam,SRR628584$mapping.SRR628584.bam,SRR628585$mapping.SRR628585.bam,SRR628586$mapping.SRR628586.bam,SRR628587$mapping.SRR628587.bam,SRR628588$mapping.SRR628588.bam,SRR628589$mapping.SRR628589.bam)
data <- rename(data, Geneid = SRR628582.Geneid)
data <- rename(data, SRR628582 = SRR628582.mapping.SRR628582.bam)
data <- rename(data, SRR628583 = SRR628583.mapping.SRR628583.bam)
data <- rename(data, SRR628584 = SRR628584.mapping.SRR628584.bam)
data <- rename(data, SRR628585 = SRR628585.mapping.SRR628585.bam)
data <- rename(data, SRR628586 = SRR628586.mapping.SRR628586.bam)
data <- rename(data, SRR628587 = SRR628587.mapping.SRR628587.bam)
data <- rename(data, SRR628588 = SRR628588.mapping.SRR628588.bam)
data <- rename(data, SRR628589 = SRR628589.mapping.SRR628589.bam)


###### Métadonnées: récuperation au format data.frame des labels :
samples = c("SRR628582", "SRR628583", "SRR628584", "SRR628585", "SRR628586", "SRR628587", "SRR628588", "SRR628589")
Type_M_WT = c("M", "M", "WT", "WT", "WT", "WT", "WT", "M")
MetaData = data.frame(samples, Type_M_WT)


###### Visualiser les données :
head(data)
print(MetaData)


###### Construction de l'objet DESEQDataSet
dds <- DESeqDataSetFromMatrix(countData=data, colData=MetaData, design=~Type_M_WT, tidy = TRUE)

###### Regardons à quoi ressemble de dds 

print(dds)

###### Lancons la fonction DESEQ maintemant 

dds <- DESeq(dds)

###### Regardons la table de résultat

res <- results(dds)
head(results(dds, tidy=TRUE)) #let's look at the results table

######  Summary de l'expression différentiel des gènes 

summary(res) #summary of results

#####  Trier la liste récapitulative par p-value

res <- res[order(res$padj),]
head(res, 11)

##### On réalise des plotCounts pour comparer les comptes normalisés entre les mutants et WT pour quelques gènes

par(mfrow=c(2,3))  # par(mfrow=c(2,3)) si je fais 6 graphes, ça les classes en 2 lignes 3 colonnes par exemple 

# exemple du gène ENSG00000115524 = SF3B1 dont parlent les articles 

plotCounts(dds, gene="ENSG00000115524", intgroup="Type_M_WT")

# exemple des 9 autres gènes les plus exprimés c'est à dire dont la p-value est la plus faible 

plotCounts(dds, gene="ENSG00000164659", intgroup="Type_M_WT")
plotCounts(dds, gene="ENSG00000234261", intgroup="Type_M_WT")
plotCounts(dds, gene="ENSG00000153283", intgroup="Type_M_WT")
plotCounts(dds, gene="ENSG00000214212", intgroup="Type_M_WT")
plotCounts(dds, gene="ENSG00000142621", intgroup="Type_M_WT")
#plotCounts(dds, gene="ENSG00000186973", intgroup="Type_M_WT")
#plotCounts(dds, gene="ENSG00000286122", intgroup="Type_M_WT")
#plotCounts(dds, gene="ENSG00000103449", intgroup="Type_M_WT")
#plotCounts(dds, gene="ENSG00000007038", intgroup="Type_M_WT")

##### Réalisation du Volcano Plot 

#reset par
par(mfrow=c(1,1))
# Faire un volcano plot basique 
with(res, plot(log2FoldChange, -log10(pvalue), pch=20, main="Volcano plot", xlim=c(-10,10)))

# Ajouter des points colorés: bleu si padj<0.01, rouge si log2FC>1 et padj<0.05)
with(subset(res, padj<.01 ), points(log2FoldChange, -log10(pvalue), pch=20, col="blue"))
with(subset(res, padj<.01 & abs(log2FoldChange)>4), points(log2FoldChange, -log10(pvalue), pch=20, col="red"))


####### ACP avec library DESeq2

# Nous devons d'abord transformer les données bruts de comptage
# La fonction vst effectuera une transformation stabilisant la variance.

vsdata <- vst(dds, blind=FALSE)
plotPCA(vsdata, intgroup="Type_M_WT")
