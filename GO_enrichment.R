#!/usr/bin/Rscript

"""
## Date Created: 06-Aug-2019
## Author: N. S. Araujo

## Script to identify GO enriched terms among the DET when compared to the entire transcriptome set.
## Additionally it describes how to use TopGo outuput as GOplot input

"""

####################################
### For T. angustula

## Calling libraries
library(topGO)
library("Rgraphviz")
library("GO.db")

## Creating the object of TopGo

## Custom annotation for BP
geneID2GO <- readMappings(file = "Jt-NURFOR_all_BP_TopGo.map", sep = "\t", IDsep = ",")

## Creating the gene list and indicating the genes of interest
gene_names <- read.table("Jt-NURFOR_all_BP_TopGo.map", sep = "\t")
gene_names <- gene_names$V1
interest_genes <- read.table("G1_vs_G2.DESeq2.DE_results.P1e-3_C2.combined.IDs")
## fix the gene ids (Cluster-1108.0|50|0.969391)
interest_genes$V1 <- gsub("\\|.*","",interest_genes$V1)
interest_genes <- interest_genes$V1

geneList <- factor(as.integer(gene_names %in% interest_genes))
names(geneList) <- gene_names
summary(geneList)

## Combining everything in the TopGo object
GOdata <- new("topGOdata", ontology = "BP", description = "Jt_ForNur", allGenes = geneList, annot = annFUN.gene2GO, gene2GO = geneID2GO)

## Generate GO graph
graph(GOdata) 

## Acessing GOdata
numGenes(GOdata) # annotated genes that will be used
head(geneScore(GOdata, use.names = FALSE))
head(sigGenes(GOdata)) # Significant genes
numSigGenes(GOdata)

## Some statistics for the GO
termStat(GOdata)
allTermStats <- (termStat(GOdata)) # 1766 terms in total
## Genes in each GO
allGOGenes <- genesInTerm(GOdata)
allGOGenes["GO:0007010"] ## with this you can get all genes in the GO term

####
#### Running the enrichment test
######

## Classical p-values (other algorithms and statistics can be used)
resultFisher <- runTest(GOdata, algorithm = "classic", statistic = 'fisher')
resultFisher  ## results

## Table with significantly enriched terms' results
allRes <- GenTable(GOdata, classicFisher = resultFisher, ranksOf = "classicFisher", topNodes= 30)
write.table(allRes, "Jt_ForNur_GOenrich.txt", sep="\t", row.names=FALSE,col.names=TRUE, quote=FALSE)

## Generated subgraph by enriched terms
showSigOfNodes(GOdata, score(resultFisher), firstSigNodes = 30, useInfo = 'np')
printGraph(GOdata, resultFisher, firstSigNodes = 30, useInfo = 'np', pdfSW = TRUE,fn.prefix = "sampleFile")


######
### Preparing files for the GOplot figure
######

## All genes resulting table
## topNodes == Size of allTermStats
allRes <- GenTable(GOdata, classicFisher = resultFisher, ranksOf = "classicFisher", topNodes= 3670)
res <- subset(allRes, select=c("GO.ID", "Term","classicFisher"))
res$Genes <- NA  ## Create new clomun
head(res)

## Now I just have to merge this table with allGOGenes info

for(i in 1:length(allRes$GO.ID)){
  GO <- allRes$GO.ID[i]
  GOinfo <- allGOGenes[GO]
  res$Genes[i] <- GOinfo
}

## unlist out
res$Genes <- vapply(res$Genes, paste, collapse = ", ", character(1L))

## Add new column
res$Category <- "BP" ## all annotation is for biological process in my file
head(res)

## Order table according to GOplot
res.ordered <- res[,c(5,1,2,4,3)]
head(res.ordered)

## Fixing col names according GOplot
colnames(res.ordered) <- c("Category", "ID", "Term", "Genes", "adj_pval")
head(res.ordered)

## Okay, now just saving this file
write.table(res.ordered, "Jt_ForNur_GOstats_4GOPlot.txt", sep="\t", row.names=FALSE,col.names=TRUE, quote=FALSE)


#### For the expression input data

## Get the counts for the transcripts
counts_matrix <- read.table("matrix.counts.matrix.G1_vs_G2.DESeq2.DE_results", sep = "\t", header=T)
head(counts_matrix)

## fix the gene ids (Cluster-1108.0|50|0.969391)
counts_matrix$id <- gsub("\\|.*","",counts_matrix$id)

## Filter the information GOplot needs
gene_expression <- subset(counts_matrix, select=c("id","log2FoldChange","pvalue","padj")) 
head(gene_expression)

## Fixing col names
colnames(gene_expression) <- c("ID", "logFC", "P.Value", "adj.P.Val")
head(gene_expression)

## Okay, now just saving this file
write.table(gene_expression, "Jt_ForNur_GENEstats_4GOPlot.txt", sep="\t", row.names=FALSE,col.names=TRUE, quote=FALSE)


#####################

## For B. terrestris same code was used just ajusting for specific file names and parameters.



