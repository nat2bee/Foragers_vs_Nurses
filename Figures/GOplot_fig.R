#!/usr/bin/Rscript

"""
## Date Created: 6-Aug-2019
## Author: N. S. Araujo

## Script used to create the GOplot figure
"""

## Call libraries
library(GOplot)
library(dplyr)


##################
### For B. terrestris

## Cleaning the environment
rm(list = ls())

## Open datasets
res.ordered <- read.table("Bt_ForNur_GOstats_4GOPlot.txt", sep = "\t", header=T)
gene_expression <- read.table("Bt_ForNur_GENEstats_4GOPlot.txt", sep = "\t", header=T)

# Generate the plotting object
circ <- circle_dat(res.ordered, gene_expression)

## Filter circ for the DET
interest_genes <- read.table("G1_vs_G2.DESeq2.DE_results.P1e-3_C2.combined.IDs")
## fix the gene ids (Cluster-1108.0|50|0.969391)
interest_genes$V1 <- gsub("\\|.*","",interest_genes$V1)
interest_genes <- interest_genes$V1
head(interest_genes)
## all to upper letter
interest_genes <- toupper(interest_genes)
head(interest_genes)

## keep only genes in DET in circ

DEGcirc <- as.data.frame(matrix(,0,length(names(circ))))
names(DEGcirc) <- names(circ)
l <- 0

for(i in 1:length(circ$genes) ){
  g <- circ$genes[i]
  if( g %in% interest_genes){
    l <- l+1
    DEGcirc[l,] <- circ[i,]
  }
}



## Select terms to show in outer circles (third level terms from B. terrestris subgraph)
int_process <- c("cellular metabolic process","primary metabolic process",
                 "nitrogen compound metabolic process","organic substance metabolic process",
                 "transposition")


## Get list of DET that can be plotted
chord_genes <- subset(DEGcirc, select=c("genes","logFC"))
colnames(chord_genes) <- c("ID", "logFC")
chord_genes <- distinct(chord_genes)
head(chord_genes)

## Create chord
chord <- chord_dat(data = DEGcirc, genes = chord_genes, process = int_process)


## Tunning collor to have the same collor for both species
total_terms <- c("cellular metabolic process","primary metabolic process",
                 "nitrogen compound metabolic process","organic substance metabolic process",
                 "transposition",
                 "catabolic process", "oxidation-reduction process")

col.terms <- brewer.pal(n = length(total_terms), name = "Set3")


##### Plot the graphic
pdf.name <- "Bt_GoPlot_DET-L2_loFC.pdf"
pdf (pdf.name, width=14, height=20, compress = F)
GOCluster(DEGcirc, int_process, clust.by = 'logFC', lfc.col = c('yellow', 'black', 'blue'), lfc.min = -6, lfc.max = 6, term.width=2.5, lfc.width=0.7,
          term.col = c(col.terms[1],col.terms[2], col.terms[3],col.terms[4], col.terms[5]))
dev.off()





################################
### For T. diversipes

## Cleaning the environment
rm(list = ls())

## Open datasets
res.ordered <- read.table("Jt_ForNur_GOstats_4GOPlot.txt", sep = "\t", header=T)
gene_expression <- read.table("Jt_ForNur_GENEstats_4GOPlot.txt", sep = "\t", header=T)

# Generate the plotting object
circ <- circle_dat(res.ordered, gene_expression)

## Filter the DET to plot in the graph
interest_genes <- read.table("G1_vs_G2.DESeq2.DE_results.P1e-3_C2.combined.IDs")
## fix the gene ids (Cluster-1108.0|50|0.969391)
interest_genes$V1 <- gsub("\\|.*","",interest_genes$V1)
interest_genes <- interest_genes$V1
head(interest_genes)
## all to upper letter
interest_genes <- toupper(interest_genes)
head(interest_genes)

## keep only genes in DET in circ
DEGcirc <- as.data.frame(matrix(,0,length(names(circ))))
names(DEGcirc) <- names(circ)
l <- 0

for(i in 1:length(circ$genes) ){
  g <- circ$genes[i]
  if( g %in% interest_genes){
    l <- l+1
    DEGcirc[l,] <- circ[i,]
  }
}

## Select terms to show in outer circles (third level terms from T. angustula subgraph)
int_process <- c("cellular metabolic process","primary metabolic process",
                 "nitrogen compound metabolic process","organic substance metabolic process",
                 "catabolic process", "oxidation-reduction process")
                 

## Get list of DET that can be plotted
chord_genes <- subset(DEGcirc, select=c("genes","logFC"))
colnames(chord_genes) <- c("ID", "logFC")
chord_genes <- distinct(chord_genes)
head(chord_genes)


## Create chord
chord <- chord_dat(data = DEGcirc, genes = chord_genes, process = int_process)


## Tunning collor to have the same collor for both species
total_terms <- c("cellular metabolic process","primary metabolic process",
                 "nitrogen compound metabolic process","organic substance metabolic process",
                 "transposition",
                 "catabolic process", "oxidation-reduction process")

col.terms <- brewer.pal(n = length(total_terms), name = "Set3")


#### Generating the graph
pdf.name <- "Jt_GoPlot_DET-L2_logFC.pdf"
pdf (pdf.name, width=14, height=20, compress = F)
GOCluster(DEGcirc, int_process, clust.by = 'logFC', lfc.col = c('yellow', 'black', 'blue'), lfc.min = -6, lfc.max = 6, term.width=2.5, lfc.width=0.7,
term.col = c(col.terms[1],col.terms[2], col.terms[3],col.terms[4], col.terms[6],col.terms[7]))
dev.off()









