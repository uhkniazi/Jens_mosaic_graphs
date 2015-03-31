# Jens_mosaic_graphs
Graphs of DE genes in the MOSAIC study 

# mosaic_genes_reactome.R
Reads a set of genes and log fold changes from a csv file. Removes any duplicate genes (as gene names are in symbols)
instead of enterez ids. Gets the corresponding enterez ids from org.hs.eg.db and removes gene symbols with different
ids i.e. duplicated gene symbols and duplicated enterez ids. Furthermore, it also removes any genes symbols with no
enterez ids (these happen perhaps due to different versions of annotations). Get the reactome ids from rectome.db using
the enterez ids and remove any NA values i.e. no mapping. Convert this data frame into a bipartite graph. After some
house keeping on the graph, the reactome terms that are too common are removed. This is done by looking at the degree
distribution of the vertices of the second kind and fitting a poisson and negative bionomial distribution on top.
Which ever fits better, we use that (the parameter is the mean of the degree distribution) and anything outside the
0.05 interval is removed i.e. reactome terms that are too common, as this creates too many connections in the graph.
Create the CGraph class object, and project to one dimension, assigning observed to expected probabilities as weights.
Remove low observed to expected probabilities (this parameter can be changed when you want more edges in the graph).
Create a second graph based on the correlations of the genes, and select an appropriate cutoff for the edge criteria
(i.e. a correlation of e.g. -0.7 to 0.7 is removed). Intersect the 2 graphs to get a smaller sized graph for further
analysis.
