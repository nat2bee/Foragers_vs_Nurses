#!/usr/bin/Rscript

"""
## Date Created: 14-Aug-2019
## Author: N. S. Araujo###### Other Fig

## Creat the waffle graphs that show the context in which C methylation occurs on different gene sets.
"""

install.packages("waffle")
library(waffle)

######################
### For B. terrestris

## create the pdf to save the figure
pdf ("Bterrestris_MethC.pdf", width=10, height=6)

## Data info to plot
trans.Bombus <- c("CG (57.52%)"=57, "CC (8.67%)" =9, 
                  "CA (18.87%)" =19, "CT (14.89%)"=15)

DEG.Bombus <- c("CG (51.17%)"=51, "CC (11.35%)" =11, 
                "CA (19.87%)" =20, "CT (17.53%)"=18)

UPfor.Bombus <- c("CG (46.59%)"=47, "CC (13.01%)" =13, 
                  "CA (20.98%)" =21, "CT (19.34%)"=19)

UPnur.Bombus <- c("CG (58.21%)"=58, "CC (8.58%)" =9, 
                  "CA (18.45%)" =18, "CT (14.70%)"=15)


## Plot
iron(
  waffle(trans.Bombus, rows=5, size=0.7, colors = c( "darkgrey","lightskyblue4","lightskyblue3",  "lightskyblue")),
  
  waffle(DEG.Bombus, rows=5, size=0.7, colors = c( "darkgrey","lightskyblue4","lightskyblue3",  "lightskyblue")),
  
  waffle(UPfor.Bombus, rows=5, size=0.7, colors = c( "darkgrey","lightskyblue4","lightskyblue3",  "lightskyblue")),
  
  waffle(UPnur.Bombus, rows=5, size=0.7, colors = c( "darkgrey","lightskyblue4","lightskyblue3",  "lightskyblue"))
)

dev.off()





######################
### For T. angustula

## create the pdf to save the figure
pdf ("Tangustula_MethC.pdf", width=10, height=6)

## Data info to plot
trans.Jatai <- c("CG (35.90%)"=36, "CC (16.06%)" =16, 
                 "CA (22.03%)" =22, "CT (25.92%)"=26)

DEG.Jatai <- c("CG (25.49%)"=25, "CC (19.56%)" =20, 
               "CA (27.68%)" =28, "CT (27.17%)"=27)

UPfor.Jatai <- c("CG (27.60%)"=27, "CC (19.00%)" =19, 
                 "CA (27.71%)" =28, "CT (25.67%)"=26)

UPnur.Jatai <- c("CG (24.62%)"=25, "CC (19.35%)" =19, 
                 "CA (28.14%)" =28, "CT (27.74%)"=28)


## Plot
iron(
  waffle(trans.Jatai, rows=5, size=0.7, colors = c("darkgrey", rgb(0.7,0,0,0.6), rgb(0.9,0.2,0,0.5), rgb(1,0,0,0.7))),
  
  waffle(DEG.Jatai, rows=5, size=0.7, colors = c("darkgrey", rgb(0.7,0,0,0.6), rgb(0.9,0.2,0,0.5), rgb(1,0,0,0.7))),
  
  waffle(UPfor.Jatai, rows=5, size=0.7, colors = c("darkgrey", rgb(0.7,0,0,0.6), rgb(0.9,0.2,0,0.5), rgb(1,0,0,0.7))),
  
  waffle(UPnur.Jatai, rows=5, size=0.7, colors = c("darkgrey", rgb(0.7,0,0,0.6), rgb(0.9,0.2,0,0.5), rgb(1,0,0,0.7)))
)


dev.off()
