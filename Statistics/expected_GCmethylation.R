#!/usr/bin/Rscript

"""
## Date Created: 18-Jun-2019
## Author: N. S. Araujo

## Script to test whether the amount of CG or non-CG methylation in nurses transcriptome is different from expected if sites where randonmly slected based on the mean CG content.

### It's expected that random sampling will actually correspond to the % of GC sites in the transcriptome...
### Nevertheless, it is important to remember that in "real life" methylation is not random.
### Therefore results just means that, there are more (or less) CG methylation than would be expected if sites were chosen randomly.
"""

####################
### For T. angustula

## Creat a dataset to represent the C methylation in the transcriptome based on real values; same number of C sites, in the right context.
CG <- c(rep("CG", 404491))
CA <- c(rep("CA", 796461))
CC <- c(rep("CC", 457240))
CT <- c(rep("CT", 961752))
total.Cs <- c(CG,CA,CC,CT)
summary(total.Cs)
head(total.Cs)

## Portion of CG sites in the transcriptome
portionGC <- (length(CG)*100)/length(total.Cs)
# 15.43892

## Amount of CG methylation really found
found.GCmeth <- 0.3590 * 100

## Difference from observed GC and GC sites available
diff.GC <- abs(found.GCmeth - portionGC)

## Define sample size based on the methylation mean
mean.meth <- 0.0124 *100
s.size <- round((length(total.Cs)*mean.meth)/100) ## Sample this number of Cs from the dataset

## nule hipothesis sampling
results <- c()
n <- 10000

for(i in 1:n){   # Get the mean CG methilation from random samples of the total transcritpome
  s1 <- sample(total.Cs, size = s.size, replace = F)
  s.simula1 <- (sum(s1 == "CG")*100)/length(s1)
  results[i] <- s.simula1
}

summary(results)
"""
   Min. 1st Qu.  Median    Mean 3rd Qu.    Max. 
  14.72   15.31   15.44   15.44   15.57   16.18
"""

## Find how many times the same amount of GC sites methylated in real data is methylated in random sampling  
p.value1 <- sum(results >= found.GCmeth)/n 
# 0
## Find how many times the difference observed is greater than the found
p.value2 <- sum(abs(results - portionGC) >= diff.GC)/n 
# 0


####################
### For Bombus

## Creat a dataset to represent the C methylation in the transcriptome based on real values; same number of C sites, in the right context.
CG <- c(rep("CG", 2083666))
CA <- c(rep("CA", 2756292))
CC <- c(rep("CC", 1476880))
CT <- c(rep("CT", 2689407))
total.Cs <- c(CG,CA,CC,CT)
summary(total.Cs)
head(total.Cs)

## Portion of CG sites 
portionGC <- (length(CG)*100)/length(total.Cs)
# 23.13579

## Amount of CG methylation 
found.GCmeth <- 0.5752 * 100
## Define sample size based on the methylation mean
mean.meth <- 0.0066 *100
s.size <- round((length(total.Cs)*mean.meth)/100) ## Sample this number of Cs from the dataset

## nule hipothesis sample
results <- c()
n <- 10000

for(i in 1:n){   # Simule dois conjuntos de dados de todo o transcriptoma para criar a hipotese nula
  s1 <- sample(total.Cs, size = s.size, replace = F)
  s.simula1 <- (sum(s1 == "CG")*100)/length(s1)
  results[i] <- s.simula1
}
summary(results)
"""
   Min. 1st Qu.  Median    Mean 3rd Qu.    Max. 
  22.50   23.02   23.14   23.14   23.26   23.78 
"""

## Find how many times this amount of GC sites is methylated  
p.value <- sum(results >= found.GCmeth)/n 
# 0



