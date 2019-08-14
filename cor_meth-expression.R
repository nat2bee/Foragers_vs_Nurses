"""
#!/usr/bin/Rscript
## Date Created: 18-Jun-2019
## Author: N. S. Araujo

### Estimate the correlation coeficient between DNA methylation and transcript expression
## Since methylation is based only in nurses bissulfite sequencing, only nurses expression profile will be used.

"""


##################
### For B. terrestris

## Getting total methylation data
meth.data <- read.table("Bt_fornur.mstat.data",header = T)
meth.data <- meth.data[11:length(meth.data$context),]

## Fixing clusters names
meth.data$context <-sub("_", "-", meth.data$context, fixed=TRUE)
meth.data$context <- gsub("_", ".",  meth.data$context)

## Matrix of normalized counts for nurses
trans.data <- read.table("matrix.TPM.not_cross_norm",header = T)

## Match transcripts methylation in counts
sum(meth.data$context %in% row.names(trans.data))
trans.data <- trans.data[row.names(trans.data) %in% meth.data$context,]

## Merge both tables based on the transcript ID
merged.table <- merge(meth.data, trans.data, by.x="context", by.y="row.names", all = TRUE)
merged.table <- na.omit(merged.table) ## exclude missing values

## Plot to visualy verify correlation and estimate Spearman coeficient
## In CG context
plot(merged.table$CG ~ log2(merged.table$Bt1N_RSEM.trinity))
plot( log2(merged.table$Bt1N_RSEM.trinity) ~ merged.table$CG)
cor(x=merged.table$CG, y=merged.table$Bt1N_RSEM.trinity, use = "complete.obs", method = "spearman")
# 0.230354

## In CW context
plot(merged.table$CW ~ log2(merged.table$Bt1N_RSEM.trinity))
cor(x=merged.table$CW, y=merged.table$Bt1N_RSEM.trinity, use = "complete.obs", method = "spearman")
# 0.08365438


##### Same code can be used to different dataset, like only the DET instead of the full transcriptome set, or T. angustula data.

