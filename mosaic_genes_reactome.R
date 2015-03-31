#NAME: mosaic_genes_reactome.R
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
dfDat = read.csv(file.choose(), header=FALSE)

# get gene symbols
s = gsub('\'(\\w+)-.+\'', replacement = '\\1', x = dfDat$V1)
dfDat$V1 = s
# remove duplicated genes
f = !duplicated(dfDat$V1)
dfDat = dfDat[f,]
# get the enterez ids of the genes
dfEnt = select(org.Hs.eg.db, keys = dfDat$V1, columns = 'ENTREZID', keytype = 'SYMBOL')
# remove the duplicated symbols (one instance of it)
f = !duplicated(dfEnt$SYMBOL)
dfEnt = dfEnt[f,]
# match the 2 data frames and get enterez ids
i = match(dfDat$V1, dfEnt$SYMBOL)
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
# which distribution can approximate the frequency of reactome terms
hist(t, prob=T, main='distribution of reactome terms with number associated genes', breaks=seq(-0.5, 6.5, by=1))
# try negative binomial and poisson distributions
# parameterized on the means
dn = dnbinom(0:7, size = mean(t), mu = mean(t))
dp = dpois(0:7, mean(t))
lines(0:7, dn, col='black', type='b')
lines(0:7, dp, col='red', type='b')
legend('topright', legend =c('nbinom', 'poi'), fill = c('black', 'red'))
# a poisson distribution with mean(t) fits well - use this as cutoff
i = round(exp(qpois(0.05, mean(t), lower.tail = F)))
c = names(which(ivDegGo>i))
v = V(oIGbp)[c]
oIGbp = delete.vertices(oIGbp, v)

# delete any orphan genes left behind, i.e. genes with no reactome terms
d = degree(oIGbp)
oIGbp = delete.vertices(oIGbp, which(d == 0))

# assign associated meta data to the gene vertices
rownames(dfDat) = as.character(dfDat$ENTREZID)
f = V(oIGbp)$type
n = V(oIGbp)[f]$name
V(oIGbp)[n]$logFC = dfDat[n,'V2']
V(oIGbp)[n]$symbol = dfDat[n,'V1']
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
# switch the weights with obs to exp ratio
E(oIGProj)$weight_old = E(oIGProj)$weight
E(oIGProj)$weight = E(oIGProj)$ob_to_ex

## remove low observed to expected probabilities
x = E(oIGProj)$weight
# NOTE: this cutoff can be changed, the lower it is the more edges in the graph
i = quantile(x, 0.75)
i = which(x <= i)
# remove these edges
oIGProj = delete.edges(oIGProj, edges = i)

## do a plot to visualize graph structure
p.old = par(mar=c(1,1,1,1))

plot(oIGProj, vertex.label=NA, vertex.size=1, layout=layout.fruchterman.reingold, vertex.frame.color=NA)

### Correlations graph step
dfCounts = read.csv(file.choose(), header=F)
# get gene symbols
s = gsub('\'(\\w+)-.+\'', replacement = '\\1', x = dfCounts$V1)
dfCounts$V1 = s
# remove duplicated genes
f = !duplicated(dfCounts$V1)
dfCounts = dfCounts[f,]
# select only those genes that have a matching symbol in other graph
f = dfCounts$V1 %in% V(oIGProj)$symbol
dfCounts = dfCounts[f,]
# get matching enterez ids for the symbols
i = match(dfCounts$V1, V(oIGProj)$symbol)
n = V(oIGProj)$name
n2 = n[i]
mCounts = as.matrix(dfCounts[,-1])
rownames(mCounts) = n2
mCor = cor(t(mCounts))
diag(mCor) = 0

# create the graph of correlations
oIGcor = graph.adjacency(mCor, mode='min', weighted=T)
## house keeping and cleaning graph
c = E(oIGcor)$weight
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
write.graph(oIGProj, 'Results/graph_full.graphml', format = 'graphml')
write.graph(ig.2, 'Results/graph_intersection_0.5.graphml', format = 'graphml')








