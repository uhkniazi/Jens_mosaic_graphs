#NAME: 02_dataset_2.R
#DESC: input the gene lists from the csv file and create reactome based graph
#AUTH: u.niazi@imperial.ac.uk
#Date: 31/03/2015

## source libraries
source('~/Dropbox/Home/Data/R/My_Libraries/NGS_functions.R')
source('~/Dropbox/Home/Data/R/My_Libraries/CGraph.R')
source('~/Dropbox/Home/Data/R/My_Libraries/CGeneAnnotation.R')
library(igraph)

# databases 
library(org.Hs.eg.db)
library(GO.db)
library(reactome.db)

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
r = seq(r[1]-0.5, r[2]+0.5, by=1)
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
# get the neighbours i.e. reactome ids for each enterez id
for (i in 1:length(n)){
  nv = neighbors(oIGbp, V(oIGbp)[n[i]])
  V(oIGbp)[n[i]]$Reactome = list(V(oIGbp)[nv]$name)
}

## graph projection to one dimension
# create the CGraph object and calculate obs to exp weights after projection
obj = CGraph(oIGbp)
# create a projection of the graph 
oIGProj = CGraph.getProjectedGraph(obj)
rm(obj)

## some genes are orphans as they don't share reactome terms with others, remove those
d = degree(oIGProj)
oIGProj = delete.vertices(oIGProj, which(d == 0))

## house keeping
# show overexpressed genes as circle and underexpressed as square
V(oIGProj)$shape = ifelse(V(oIGProj)$logFC > 0, 'circle', 'square')
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
r = round(range(w2))
s = seq(r[1]-0.5, r[2]+0.5, by = 1)
hist(w2, prob=T, breaks=s, main='distribution of obs to exp ratios', 
     xlab='square root obs to exp ratio', ylab='')
dp = dpois(r[1]:r[2], lambda = median(w2))
dn = dnbinom(r[1]:r[2], size = median(w2), mu = median(w2))
lines(r[1]:r[2], dp, col='red')
lines(r[1]:r[2], dn, col='blue')
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
t = abs(c)
hist(t, prob=T)
# the correlation looks fairly normal
dn = dnorm(seq(0, 1, 0.01), mean(t), sd(t))
lines(seq(0, 1, 0.01), dn)
summary(t)
dn = dnorm(seq(0, 1, 0.01), median(t), mad(t))
lines(seq(0, 1, 0.01), dn, col='red')
# choose appropriate cutoff based on data
f = which((c >= -0.5 & c <= 0.5))
# if we need genes that are not correlated
f = which((c <= -0.5 | c >= 0.5))
oIGcor = delete.edges(oIGcor, edges = f)

### graph intersection
# intersect the 2 graphs
l = list(oIGProj, oIGcor)
ig.1 = graph.intersection(l)
E(ig.1)$weight = E(ig.1)$weight_1
# color the positive and negative edges differently
col = ifelse(E(ig.1)$weight_2 < 0, 'red', 'black')
E(ig.1)$color = col

par(mar=c(1,1,1,1))
plot(ig.1, vertex.label=NA, vertex.size=1, layout=layout.fruchterman.reingold, vertex.frame.color=NA)
# delete orphan nodes
d = degree(ig.1)
ig.2 = delete.vertices(ig.1, which(d == 0))

par(mar=c(1,1,1,1))
plot(ig.2, vertex.label=NA, vertex.size=1, layout=layout.fruchterman.reingold, vertex.frame.color=NA)

### save the graphs for analysis on cytoscape
dir.create('Results',showWarnings = F)
## export the graphs for use in cytoscape
write.graph(oIGProj, 'Results/ds_2_reactome_terms.graphml', format = 'graphml')
write.graph(ig.2, 'Results/ds_2_decor_fp_vs_fn.graphml', format = 'graphml')







