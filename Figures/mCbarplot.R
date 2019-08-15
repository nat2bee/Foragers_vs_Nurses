#!/usr/bin/Rscript

"""
## Date Created: 18-Jun-2019
## Author: N. S. Araujo

## Barplot figure comparing mC mean methylation among transcripts groups
"""


## Create data frame with mean methylation values to plot
m.data <- data.frame(Barplot1=c(0.66,1.24),Barplot2=c(0.78,1.3),
                    Barplot3=c(0.73,1.57), Barplot4=c(0.95,1.25))

## Plor bars groupped
barplot(as.matrix(m.data), main="", ylab="Mean mC (%)",beside=TRUE,
        cex.axis = 1, space = c(0,0.3), ylim = c(0,2), col=c("gray77","gray37"), border = NA,
        names.arg = c("Transcriptome", "DET", "High foragers", "High nurses")) 

## Add legend
legend(legend = c(expression(italic("B. terrestris")), expression(italic("T. angustula"))), lwd=10,col=c("gray77","gray37"), cex = 1, 
       x = "topright", horiz = T, bty = "n")

## Add values over the bars
text(x=c(c(1:2)-0.2, c(3:4)+0.1,c(5:6)+0.4, c(7:8)+0.7), labels = c("0.66%","1.24%","0.78%*","1.30%","0.73%","1.57%","0.95%*","1.25%"), 
     y = c(0.66,1.24,0.78,1.3,0.73,1.57,0.95,1.25), pos = 3, cex = 0.8, 
     col = "black")

## Add confidence interval bars
arrows(3+0.1, 0.7188236,
       3+0.1, 0.78,angle=90,code=1,length=0.05)

arrows(4+0.1, 1.009685,
       4+0.1, 1.3,angle=90,code=1,length=0.05)

arrows(5+0.4, 0.653384,
       5+0.4, 0.73,angle=90,code=1,length=0.05)

arrows(6+0.4, 1.233139,
       6+0.4, 1.57,angle=90,code=1,length=0.05)

arrows(7+0.7, 0.8483813,
       7+0.7,0.95,angle=90,code=1,length=0.05)

arrows(8+0.7, 0.6776233,
       8+0.7, 1.25,angle=90,code=1,length=0.05)
       
       
