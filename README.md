# Data and scripts repository from the manuscript:
> **Multiple lineages, same molecular tools: task specialization is commonly regulated across all eusocial bee groups**
*Natalia de Souza Araujo, Yannick Wurm, and Maria Cristina Arias*

## Repository content

### Data files
- **Jt_fornur_dez16_lncCod.fasta.gz**: transcripts assembled based on RNASeq data from foragers and nurses of *T. angustula*. Assembly method described in the manuscript.

- **Bt_fornur.mstat.data**: Estimation of DNA methylation in the transcriptome of *B. terrestris* nurse. Analysis method described in the manuscript.

- **Jt_fornur.mstat.data**: Estimation of DNA methylation in the transcriptome of *T. angustula* nurse. Analysis method described in the manuscript.

 
### Statistics
- **common.stats.R**: R function used to test overlap significance between DET in the two species. Used parameters in the manuscript are included.

- **expected_GCmethylation.R**: R script used to test wether the mean amount of CG methylation observed is greater than expected based on the proportion of CG sites in the transcriptome.

- **methylation_mean_dev.R**: R script used to test whether the mean amount of C methylation observed in a transcripts subset is significantly greater than the general mean.


### Figures
- **mCbarplot.R**: R script to create the DNA methylation barplot
