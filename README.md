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

- **cor_meth-expression.R**: R script to estimate the Spearman' coeficient between gene expression and mC.

- **GO_enrichment.R**: R script to estimate enriched GO terms among the DET using TopGO.


### Figures
- **mCbarplot.R**: R script to create the DNA methylation barplot.

- **mC_waffle.R**: R script to create the DNA methylation waffle plot.

- **GOplot_fig.R**: R script to create the GOplot graph for third level terms from subgraphs of enriched terms.

- **euler_fig.R**: R script to create the Euler diagram of genes in common.


### Others
- **Annocript2TopGo.py**: Python script used to format Annocript output to TopGO input.



'''
## License

This work is distributed under the GPLv3 license. Reuse of code derived from this repository is permitted under two conditions:

Proper attribution (i.e., citation of the associated publication; see CITATION.cff and above).
Publication of reused scripts on an open-access platform, such as Github.




