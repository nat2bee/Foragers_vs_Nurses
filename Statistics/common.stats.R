#!/usr/bin/Rscript

"""
## Date Created: 9/8/19
## Author: N. S. Araujo

## Teste the significance (p-value) of having the "x" genes in common between two DE datasets using random sampling

## Usage: common.stats(x, sp1, sp2, n1, n2, n = 10000, p = 0.01)

## Inputs: x (number of common genes between the sps, including repetitions), 
## sp1 (R character vector of the gene names from all the annotated transcripts in the transcriptome of sp1 in UPPER CASE, uniq, excluding the umbiguous annotations, e.g. uncharacterized protein, hypothetical protein, hypothetical uncharacterized, putative uncharacterized), 
## sp2 (R character vector of the gene names from all the annotated transcripts in the transcriptome of sp2 in UPPER CASE, uniq,  excluding the umbiguous annotations, e.g. uncharacterized protein, hypothetical protein, hypothetical uncharacterized, putative uncharacterized), 
## n1 (number of annotated DE, or UP expressed, transcripts in sp1 to sample), n2 (number of annotated DE, or UP expressed, transcripts in sp2 to sample).
## s is the number of simulations to run (default is 10000)
## p is the p-value threshold (default is 0.01)

## Output: p-value of having an "x" iqual or greater then the observed "x" by chance
"""


common.stats <- function (x, sp1, sp2, n1, n2, n = 10000, p = 0.01){
  # Open list of transcripts if not an R object
  # sp1 <- readLines(sp1)
  # sp2 <- readLines(sp2)
  
  # Creat necessary variables
  common.total <- c()
  
  # Simulation of the nule hipotesis
  for(i in 1:n){   # From the total transcriptome take radom genes, not related to any biological process
    DE.sp1 <- as.vector(sample(sp1, size = n1, replace = F))
    DE.sp2 <- as.vector(sample(sp2, size = n2, replace = F))
    a <- c()
    b <- c()
    z <- 0
    w <- 0
    
    for(j in 1:length(DE.sp1)){  # Check how many of the random genes are common between the species
      if(sum(grepl(DE.sp1[j], DE.sp2, fixed = T)) > 0){
        z <- z+1
        a[z] <- DE.sp1[j]
      }
    }
    
    for(k in 1:length(DE.sp2)){  
      if(sum(grepl(DE.sp2[k], DE.sp1, fixed = T)) > 0){
        w <- w+1
        b[w] <- DE.sp2[k]
      }
    }
    combined.GOs <- unique(c(a,b))
    common.total[i] <- length(combined.GOs) # save the results of the common genes
  }
  
  h.plot <- hist(common.total, main = "Genes in common from random sampling", xlab = "# of common genes")
  # return(common.total)
  
  # Calculate and return the p-value of getting results equal or greater than the observed x value
  p.value <- sum(common.total >= x)/n
  
  # Returning the result
  if(p.value <= p){
    return(list(h.plot,cat(c("Probability of getting ", x, " genes in common or more randomly from both dataset is significant. \np-value:", p.value))))
  }
  else{
    return(list(h.plot,cat(c("Probability of getting ", x, " genes in common or more randomly from both dataset is not significant. \np-value:", p.value))))
  }
}




## Comparisons between T. angustula and B. terrestris
### All DET
common.stats(x=18, sp1=jt_full.g, sp2=bt_full.g, n1=109, n2=472, n = 10000, p = 0.01)

### Commonly highly expressed in nurser
common.stats(x=9, sp1=jt_full.g, sp2=bt_full.g, n1=87, n2=170, n = 10000, p = 0.01)

### Commonly highly expressed in foragers
common.stats(x=2, sp1=jt_full.g, sp2=bt_full.g, n1=24, n2=317, n = 10000, p = 0.01)
