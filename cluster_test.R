#NAME: cluster_test.R
#DESC: input the gene lists from the csv file and create reactome based graph
#AUTH: u.niazi@imperial.ac.uk
#Date: 31/03/2015

## source libraries
source('~/Dropbox/Home/Data/R/My_Libraries/NGS_functions.R')
source('~/Dropbox/Home/Data/R/My_Libraries/CGraph.R')
source('~/Dropbox/Home/Data/R/My_Libraries/CGeneAnnotation.R')
library(igraph)
library(RColorBrewer)

# databases 
library(org.Hs.eg.db)
library(GO.db)
library(reactome.db)
#library(KEGG.db)

p.old = par()
## data loading and formatting
dfDat = read.csv(file.choose(), header=T, stringsAsFactors=F)

# get gene symbols
# s = gsub('\'(\\w+)-.+\'', replacement = '\\1', x = dfDat$X)
# dfDat$V1 = s
# remove duplicated genes
f = !duplicated(dfDat$X)
dfDat = dfDat[f,]
# get the enterez ids of the genes
dfEnt = select(org.Hs.eg.db, keys = dfDat$X, columns = 'ENTREZID', keytype = 'SYMBOL')
# remove the duplicated symbols (one instance of it)
f = !duplicated(dfEnt$SYMBOL)
dfEnt = dfEnt[f,]
# match the 2 data frames and get enterez ids
i = match(dfDat$X, dfEnt$SYMBOL)
dfDat$ENTREZID = dfEnt[i, 'ENTREZID']
# some genes symbols do not have entrezids due to old annotations
dfDat = na.omit(dfDat)
# get the respective reactome ids
dfGraph = select(reactome.db, dfDat$ENTREZID, 'REACTOMEID', 'ENTREZID')
dfGraph = na.omit(dfGraph)
# create bipartite graph
oIGbp = graph.data.frame(dfGraph, directed = F)
# set the vertex type variable to make graph bipartite
f = rep(c(T, F), times = c(length(unique(dfGraph$ENTREZID)),length(unique(dfGraph$REACTOMEID))))
V(oIGbp)$type = f
# sanity check - is graph bipartite
is.bipartite(oIGbp)

## graph cleaning
# remove high level reactome terms, as they create too many edges
f = V(oIGbp)$type
ivDegGo = degree(oIGbp, V(oIGbp)[!f])
t = log(ivDegGo)
summary(t)
r = range(t)
r = seq(floor(r[1])-0.5, ceiling(r[2])+0.5, by=1)
# which distribution can approximate the frequency of reactome terms
hist(t, prob=T, main='distribution of reactome terms with number associated genes', breaks=r)
# try negative binomial and poisson distributions
# parameterized on the means
dn = dnbinom(0:5, size = mean(t), mu = mean(t))
dp = dpois(0:5, mean(t))
lines(0:5, dn, col='black', type='b')
lines(0:5, dp, col='red', type='b')
legend('topright', legend =c('nbinom', 'poi'), fill = c('black', 'red'))
# a poisson distribution with mean(t) fits well - use this as cutoff
# however a negative binomial will adjust for overdispertion, try both perhaps
#i = round(exp(qpois(0.05, mean(t), lower.tail = F)))
i = round(exp(qnbinom(0.05, size = mean(t), mu = mean(t), lower.tail = F)))
c = names(which(ivDegGo>i))
v = V(oIGbp)[c]
oIGbp = delete.vertices(oIGbp, v)

# delete any orphan genes left behind, i.e. genes with no reactome terms
d = degree(oIGbp)
oIGbp = delete.vertices(oIGbp, which(d == 0))

# assign associated meta data to the gene vertices
# save the old datafame first for later use
dfDat.bk = dfDat
rownames(dfDat) = as.character(dfDat$ENTREZID)
f = V(oIGbp)$type
n = V(oIGbp)[f]$name
V(oIGbp)[n]$logFC = dfDat[n,'Fold.change']
V(oIGbp)[n]$symbol = dfDat[n,'X']

## graph projection to one dimension
# create the CGraph object and calculate obs to exp weights after projection
obj = CGraph(oIGbp)
# create a projection of the graph 
oIGProj = CGraph.getProjectedGraph(obj)

## some genes are orphans as they don't share reactome terms with others, remove those
d = degree(oIGProj)
oIGProj = delete.vertices(oIGProj, which(d == 0))

## house keeping
# show overexpressed genes as circle and underexpressed as square
V(oIGProj)$shape = ifelse(V(oIGProj)$logFC > 0, 'circle', 'square')
V(oIGProj)$color = ifelse(V(oIGProj)$logFC > 0, 'red', 'blue')
# need this for later plotting using cytoscape to use as a node size parameter
V(oIGProj)$logFC_abs = abs(V(oIGProj)$logFC)

# switch the weights with obs to exp ratio
E(oIGProj)$weight_old = E(oIGProj)$weight
E(oIGProj)$weight = E(oIGProj)$ob_to_ex

## remove low observed to expected probabilities
w = E(oIGProj)$weight
# choose a cutoff by modelling the distribution shape
# it appears that the distribution follows a power law?
# taking square root means we can fit a poisson distribution
w2 = sqrt(w)
r = range(w2)
s = seq(floor(r[1])-0.5, ceiling(r[2])+0.5, by = 1)
hist(w2, prob=T, breaks=s, main='distribution of obs to exp ratios', 
     xlab='square root obs to exp ratio', ylab='')
r = round(r)
dp = dpois(r[1]:r[2], lambda = median(w2))
dn = dnbinom(r[1]:r[2], size = median(w2), mu = median(w2))
lines(r[1]:r[2], dp, col='red', type='b')
lines(r[1]:r[2], dn, col='blue', type='b')
legend('topright', legend = c('poi', 'nbin'), fill = c('red', 'blue'))


# NOTE: this cutoff can be changed, the lower it is the more edges in the graph
# use negative binomial to choose cutoff
c = qnbinom(0.05, size = median(w2), mu=median(w2), lower.tail = F)
f = which(w2 < c)
oIGProj = delete.edges(oIGProj, edges = f)

## do a plot to visualize graph structure
p.old = par(mar=c(1,1,1,1))

plot(oIGProj, vertex.label=NA, vertex.size=1, layout=layout.fruchterman.reingold, vertex.frame.color=NA)
par(p.old)
### Correlations graph step, redo this step for each type of correlation graph required
n = V(oIGProj)$name
dfCounts = dfDat[rownames(dfDat) %in% n,]
c = colnames(dfCounts)
## TAG 1: NOTE: select the type of graph needed
c = grep('Flu\\.Pos', c)
c = grep('Flu\\.Neg', c)
c = grep('Healthy', c)
c = grep('Flu\\.Pos|Flu\\.Neg|Healthy', c)
c = grep('Flu\\.Pos|Flu\\.Neg', c)
dfCounts.fp = dfCounts[,c]

# mCounts = t(dfCounts.fp)
# mDist = as.matrix(dist(t(mCounts)))
# # create graph of this matrix 
# oIGdist = graph.adjacency(mDist, mode='min', weighted=T)
# # remove gene connections that are farther apart form each other
# w = E(oIGdist)$weight
# summary(w)
# par(p.old)
# hist(w)
# # choose a cutoff by modelling the distribution shape
# # it appears that the distribution follows a power law?
# # taking square root means we can fit a poisson distribution
# w2 = sqrt(w)
# r = range(w2)
# s = seq(floor(r[1])-0.5, ceiling(r[2])+0.5, by = 1)
# hist(w2, prob=T, breaks=s, main='distribution of eucledian distances', 
#      xlab='square root of distances', ylab='')
# r = round(r)
# dp = dpois(r[1]:r[2], lambda = median(w2))
# dn = dnbinom(r[1]:r[2], size = median(w2), mu = median(w2))
# lines(r[1]:r[2], dp, col='red', type='b')
# lines(r[1]:r[2], dn, col='blue', type='b')
# legend('topright', legend = c('poi', 'nbin'), fill = c('red', 'blue'))
# 
# # NOTE: this cutoff can be changed, the lower it is the less edges in the graph
# # use poisson/negative binomial to choose cutoff
# #c = qnbinom(0.05, size = median(w2), mu=median(w2), lower.tail = F)
# c = qpois(0.05, median(w2), lower.tail = F)
# f = which(w2 >= c)
# oIGdist = delete.edges(oIGdist, edges = f)
# # get the distance matrix again / adjacency matrix
# mDist = as.matrix(get.adjacency(oIGdist, 'both', attr = 'weight', names = T))
# # invert the weights, things that are farthest will have the smallest number
# # others that are closest to each other will have the most remaining behind
# # e.g. 10 - 5 = 5 no change in weight
# # 10 - 8 = 2 as these 2 are far away, they will have a weight of 2 remaining etc
# mDist = max(mDist) - mDist
# # no connections or weights on diagonals
# diag(mDist) = 0  
# # create the graph of distances/correlations
# oIGdist = graph.adjacency(mDist, mode='min', weighted=T)
# E(oIGdist)$distance = E(oIGdist)$weight

## use correlation matrix to create graph
mCounts = t(dfCounts.fp)
mCor = cor(mCounts)
diag(mCor) = 0

# create the graph of correlations
oIGcor = graph.adjacency(mCor, mode='min', weighted=T)
## house keeping and cleaning graph
c = E(oIGcor)$weight
par(p.old)
hist(c)
summary(c)
E(oIGcor)$cor = E(oIGcor)$weight

# keep only positively correlated genes connected
f = which(c < 0.5)
oIGcor = delete.edges(oIGcor, edges = f)

### graph intersection
l = list(oIGProj, oIGcor)

ig.1 = graph.intersection(l)
# set observed to expected ratio as weight
E(ig.1)$weight = E(ig.1)$ob_to_ex

par(mar=c(1,1,1,1))
plot(ig.1, vertex.label=NA, vertex.size=1, layout=layout.fruchterman.reingold, vertex.frame.color=NA)
# delete orphan nodes
d = degree(ig.1)
ig.2 = delete.vertices(ig.1, which(d == 0))

par(mar=c(1,1,1,1))
plot(ig.2, vertex.label=NA, vertex.size=1, layout=layout.fruchterman.reingold, vertex.frame.color=NA)

# create communities in the graph
com = edge.betweenness.community(ig.2)
dend = as.dendrogram(com)
# plot the heatmap of the communities
col = brewer.pal(9, 'Greens')
n = V(ig.2)$name
mCounts.heat = t(mCounts[,colnames(mCounts) %in% n])
c = colnames(mCounts.heat)
c = gsub('^X.+(Pos|Neg|Healthy).+', '\\1', c)
colnames(mCounts.heat) = c
mCounts.heat = mCounts.heat[,order(colnames(mCounts.heat))]
par(p.old)
heatmap(mCounts.heat, Rowv = dend, scale='row', col=col, cexRow = 0.2, cexCol = 0.5), Colv = NA)

## cut the tree of communities into sub groups and find the reactome 
## pathways associated with those genes. this is a multi step process

# recreate the bipartite graph but only with nodes that are in our final graph
oIGbp = CGraph.getBipartiteGraph(obj)
# get the indices for the vertices of type 1, i.e. genes
f = V(oIGbp)$type
n = V(oIGbp)[f]$name
# get names of genes present in last graph after cutoffs
n2 = V(ig.2)$name
# intersect the names
i = !(n %in% n2)
# choose the names not commong between bipartite graph and the last graph
n = n[i]
oIGbp = delete.vertices(oIGbp, v = n)
d = degree(oIGbp)
# delete orphan nodes
oIGbp = delete.vertices(oIGbp, which(d == 0))
f = V(oIGbp)$type

hc = as.hclust(com)
plot(hc)
memb = membership(com)
memb = cutree(hc, h = 230)
table(memb)
mCent = matrix(NA, nrow = length(unique(memb)), ncol = ncol(mCounts.heat))
rv = rep(NA, length=nrow(mCent))
rv.g = rep(NA, length=nrow(mCounts.heat))
# get the centers and reactome terms for these genes from each cluster
for (i in 1:length(unique(memb))){
  # check if cluster has one member only
  if (sum(memb == i) == 1) {
    mCent[i,] = mCounts.heat[memb == i,]
  } else {
  # else if more than one member, we can use mean 
  mCent[i,] = colMeans(mCounts.heat[memb == i,])}
  rn = rownames(mCounts.heat)[memb == i]
  # get the reactome pathway names
  nei = graph.neighborhood(oIGbp, order = 1, nodes = rn)
  pw = sapply(seq_along(nei), function(x) V(nei[[x]])$name)
  pw = unlist(pw)
  pw = as.data.frame(table(pw))
  rv[i] = as.character(pw[which.max(pw$Freq), 1])
  rv.g[memb == i] = rv[i]
}
rownames(mCent) = rv
mCounts.heat.2 = mCounts.heat
rownames(mCounts.heat.2) = rv.g
dend = hc
dend$labels = rv.g
dend = as.dendrogram(dend)
heatmap(mCounts.heat.2, Rowv = dend, scale='row', col=col, cexRow = 0.1, cexCol = 0.3, Colv = NA)

hc = getHclust(oGr)
l = hc$labels
l2 = getClusterLabels(oGr)
m = mCounts[l,]



# heatmap based on means of each community/cluster
colnames(mCent) = colnames(mCounts.heat)
hc.2 = hclust(dist(mCent), members = table(memb))#, method = 'cent')
heatmap(mCent, Rowv = as.dendrogram(hc.2), scale='row', col=col, cexRow = 0.5, cexCol = 0.3, Colv = NA)
heatmap(log(mCent), Rowv = as.dendrogram(hc.2), scale='row', col=col, cexRow = 0.5, cexCol = 0.5), Colv = NA)

# recreate the bipartite graph but only with nodes that are in our final graph
oIGbp = CGraph.getBipartiteGraph(obj)
# get the indices for the vertices of type 1
f = V(oIGbp)$type
n = V(oIGbp)[f]$name
n2 = V(ig.2)$name
i = !(n %in% n2)
n = n[i]
oIGbp = delete.vertices(oIGbp, v = n)
d = degree(oIGbp)
oIGbp = delete.vertices(oIGbp, which(d == 0))
f = V(oIGbp)$type
V(oIGbp)[f]$color = 'blue'
V(oIGbp)[!f]$color = 'red'
par(mar=c(1,1,1,1))
plot(oIGbp, vertex.label=NA, vertex.size=1, layout=layout.fruchterman.reingold, vertex.frame.color=NA)


### save the graphs for analysis on cytoscape
dir.create('Results',showWarnings = F)
## export the graphs for use in cytoscape
write.graph(oIGProj, 'Results/ds_2_reactome_terms.graphml', format = 'graphml')
write.graph(ig.2, 'Results/ds_2_decor_fp_vs_fn.graphml', format = 'graphml')






