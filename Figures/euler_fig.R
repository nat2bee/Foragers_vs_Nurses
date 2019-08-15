#!/usr/bin/Rscript

"""
## Date Created: 12-Aug-2019
## Author: N. S. Araujo

## Create an euler diagram to illustrate the genes in common among the DET of both species.
"""

## Load library
library(eulerr)

## Sizes based on the number of unique gene annotations after removing non-informative terms (e.g. uncharacterized protein)
## Manual filter has not being performed so some terms might still be repeated at this stage (e.g. "CYTOCHROME C OXIDASE SUBUNIT 1 (FRAGMENT)" and "CYTOCHROME C OXIDASE SUBUNIT 1")
## After manual filter "T.angustula&B.terrestris" = 16 
## But to keep proportions manual filtering was excluded, as in the random testing stats.

fit <- euler(c("T.angustula" = 93, "B.terrestris" = 456, 
               "T.angustula&B.terrestris" = 18 ))

plot(fit, quantities = F, main=list(label="Annotated genes in common", font=2),
     legend = list(labels = c(c(expression(italic("T. angustula")),expression(italic("B. terrestris"))))))


