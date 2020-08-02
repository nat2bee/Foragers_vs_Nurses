### Code used for the donuplot showing the proportion of conserved/ taxonomicaly restricted genes
## R version 3.6.3 (2020-02-29) -- "Holding the Windsock"

## Load libraries
library(tidyverse)
library(ggrepel)


## Getting the stats as tibble 
Jt <- read_delim(file = "Jt_sps_orthogroups.txt", delim ="\t")
Bt <- read_delim(file = "Bt_sps_orthogroups.txt", delim ="\t")


## Prepare the tables to be plotted
Jt <- Jt %>%
  mutate(ymax = cumsum(all), # Compute the cumulative percentages (top of each rectangle)
         ymin = c(0, head(ymax, n=-1)), # Compute the bottom of each rectangle
         ymax.DE = cumsum(DE),
         ymin.DE = c(0, head(ymax.DE, n=-1)))


myColors <- c(gray.colors(3, start = 0.7, end = 0.9),"tan1", "tan2","tan3")
names(myColors) <- levels(as.factor(Jt$group))


scale_fill_nat <- function(...){
  ggplot2:::manual_scale(
    'fill', 
    values = myColors, drop=F, 
    ...
  )
}

scale_fill_nat2 <- function(...){
  ggplot2:::manual_scale(
    'fill', 
    values = alpha(myColors,.7), drop=F, 
    ...
  )
}


pdf.name <- "Jt_orthogroups_shared.pdf"
pdf (pdf.name, width=8, height=5, compress = F)

ggplot(Jt) +
  #geom_rect(aes(ymax=ymax.DE, ymin=ymin.DE, xmax=2.9, xmin=0, fill=group )) +
  geom_rect(aes(ymax=ymax, ymin=ymin, xmax=4, xmin=3, fill=group )) +
  xlim(c(0, 4)) + 
  theme(aspect.ratio=1) +
  scale_fill_nat2() +
  labs(fill = "") +
  geom_label_repel(aes(label=paste(all,"%"),x=3.5,y=(ymin+ymax)/2),inherit.aes = TRUE,label.size = NA, fill = NA )+
  coord_polar(theta="y")  + theme_void()

ggplot(Jt) +
  geom_rect(aes(ymax=ymax.DE, ymin=ymin.DE, xmax=2.9, xmin=0, fill=group )) +
  #geom_rect(aes(ymax=ymax, ymin=ymin, xmax=4, xmin=3, fill=group )) +
  xlim(c(0, 4)) + 
  theme(aspect.ratio=1) +
  scale_fill_nat() +
  labs(fill = "") +
  ## remove 0 values
  geom_label_repel(data=subset(Jt,DE != 0),aes(label=paste(DE,"%"),x=1.5,y=(ymin.DE+ymax.DE)/2),inherit.aes = TRUE,label.size = NA, fill = NA )+
  coord_polar(theta="y")  + theme_void()


dev.off()



## For B bombus


## Prepare the tables to be plotted
Bt <- Bt %>%
  mutate(ymax = cumsum(all), # Compute the cumulative percentages (top of each rectangle)
         ymin = c(0, head(ymax, n=-1)), # Compute the bottom of each rectangle
         ymax.DE = cumsum(DE),
         ymin.DE = c(0, head(ymax.DE, n=-1)))


myColors <- c(gray.colors(3, start = 0.7, end = 0.9),"steelblue1", "steelblue3","steelblue")
Bt$group <- factor(Bt$group, levels = c("Apinae", "Corbiculates", "Social corbiculates", "Specie-specific", "Bterrestris (G)", "Bumblebees"))
names(myColors) <- levels(as.factor(Bt$group))


scale_fill_nat <- function(...){
  ggplot2:::manual_scale(
    'fill', 
    values = myColors, drop=F, 
    ...
  )
}

scale_fill_nat2 <- function(...){
  ggplot2:::manual_scale(
    'fill', 
    values = alpha(myColors,.7), drop=F, 
    ...
  )
}


pdf.name <- "Bt_orthogroups_shared.pdf"
pdf (pdf.name, width=8, height=5, compress = F)

ggplot(Bt) +
  #geom_rect(aes(ymax=ymax.DE, ymin=ymin.DE, xmax=2.9, xmin=0, fill=group )) +
  geom_rect(aes(ymax=ymax, ymin=ymin, xmax=4, xmin=3, fill=group )) +
  xlim(c(0, 4)) + 
  theme(aspect.ratio=1) +
  scale_fill_nat2() +
  labs(fill = "") +
  geom_label_repel(aes(label=paste(all,"%"),x=3.5,y=(ymin+ymax)/2),inherit.aes = TRUE,label.size = NA, fill = NA )+
  coord_polar(theta="y")  + theme_void()

ggplot(Bt) +
  geom_rect(aes(ymax=ymax.DE, ymin=ymin.DE, xmax=2.9, xmin=0, fill=group )) +
  #geom_rect(aes(ymax=ymax, ymin=ymin, xmax=4, xmin=3, fill=group )) +
  xlim(c(0, 4)) + 
  theme(aspect.ratio=1) +
  scale_fill_nat() +
  labs(fill = "") +
  ## remove 0 values
  geom_label_repel(data=subset(Bt,DE != 0),aes(label=paste(DE,"%"),x=1.5,y=(ymin.DE+ymax.DE)/2),inherit.aes = TRUE,label.size = NA, fill = NA )+
  coord_polar(theta="y")  + theme_void()


dev.off()
