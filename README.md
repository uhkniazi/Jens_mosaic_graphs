# Jens_mosaic_graphs
Graphs of DE genes from Trevor's data 



# 02_dataset_2.R
The analysis is very similar to the previous script with a few differences. The cutoffs for the reactome term distribution
is done on the negative binomial rather than poisson distribution. The observed to expected ratio distribution appears to 
follow a power law and taking a square root of that and fit a poisson and negative binomial distribution on top - so this
second cutoff is important as previously it was done on the 75% quantile, but not we do it by fitting a distribution on top.
Look at TAG 1, and choose the type of correlation graph desired: i.e. comparisons of flu pos vs flu neg or healthy etc.

