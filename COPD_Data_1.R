#NAME: COPD_Data_1.R
#DESC: input the gene lists from the csv file and create reactome based graph and CGraphClust based heatmap
#AUTH: u.niazi@imperial.ac.uk
#Date: 10/06/2015

## source libraries
source('~/Dropbox/Home/Data/R/My_Libraries/NGS_functions.R')
source('../CGraphClust/CGraphClust.R')
library(igraph)

# databases 
library(org.Hs.eg.db)
library(GO.db)
library(reactome.db)

p.old = par()
##################################################################################
## data loading and formatting
dfDat = read.csv(file.choose(), header=T, stringsAsFactors=F)


# get gene symbols as they have an underscore before them
s = gsub('\'(.+)\'', replacement = '\\1', x = dfDat$X)
dfDat$X = s
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
rownames(dfDat) = dfDat$ENTREZID
# get the respective reactome ids
dfGraph = select(reactome.db, dfDat$ENTREZID, 'REACTOMEID', 'ENTREZID')
dfGraph = na.omit(dfGraph)

## create correlation matrix
n = unique(dfGraph$ENTREZID)
# some genes do not have reactome ids
dfCounts = dfDat[rownames(dfDat) %in% n,]
# remove the extra columns
dfCounts.p = dfCounts[,2:34]
mCounts = t(dfCounts.p)
mCor = cor(mCounts)
hist(mCor, prob=T)

# create a factor to subsample the data
c = colnames(dfCounts.p)
fSamples = gsub('(\\w+).*', '\\1', c)
dfCounts.p.o = dfCounts.p[,order(fSamples)]
fSamples = fSamples[order(fSamples)]
fSamples = as.factor(fSamples)
# stabalize the data
mCounts = t(dfCounts.p.o)
mCounts = apply(mCounts, 2, function(x) f_ivStabilizeData(x, fSamples))
rownames(mCounts) = fSamples
mCor = cor(mCounts)
hist(mCor, prob=T)
# # subsample only samples of type 2 and 3
# mCounts = dfCounts.p.o#[,fSamples != '2']
# mCounts = t(mCounts)
# mCor = cor(mCounts)
# hist(mCor)

## create the clustering graph
oGr = CGraphClust(dfGraph, abs(mCor), iCorCut = 0.5)#, iCorCut = 0.7)

### analysis 
# plot the main communities and centrality graphs
fGroups = fSamples
ig = getFinalGraph(oGr)
par(mar=c(1,1,1,1)+0.1)
set.seed(1)
ig = f_igCalculateVertexSizesAndColors(ig, t(mCounts), fGroups, bColor = T, iSize = 20)
plot(getCommunity(oGr), ig, vertex.label=NA, layout=layout_with_fr, 
     vertex.frame.color=NA, mark.groups=NULL, edge.color='lightgrey')
set.seed(1)
ig = getFinalGraph(oGr)
ig = f_igCalculateVertexSizesAndColors(ig, t(mCounts), fGroups, bColor = F, iSize = 30)
plot(getCommunity(oGr), ig, vertex.label=NA, layout=layout_with_fr, 
     vertex.frame.color=NA, edge.color='darkgrey')


# look at the graph centrality properties
set.seed(1)
ig = plot.centrality.graph(oGr)
ig = f_igCalculateVertexSizesAndColors(ig, t(mCounts), fGroups, bColor = F, iSize = 30)
set.seed(1)
plot(ig, vertex.label=NA, layout=layout_with_fr, vertex.frame.color=NA, edge.color='darkgrey')
par(p.old)
plot.centrality.diagnostics(oGr)

# get the centrality parameters
mCent = mPrintCentralitySummary(oGr)
# top  5% of the vertices from each category and largest clique
l = lGetTopVertices(oGr, iQuantile = 0.85)
# top genes based on centrality parameters
dfClique = f_dfGetGeneAnnotation(l$clique)
cvSum = f_csGetGeneSummaryFromGenbank(dfClique$ENTREZID)
cvSum.2 = dfClique$SYMBOL
dfClique$Summary = cvSum[cvSum.2]
write.csv(dfClique, file='Temp/dfClique.csv')

dfDegree = f_dfGetGeneAnnotation(l$degree)
cvSum = f_csGetGeneSummaryFromGenbank(dfDegree$ENTREZID)
cvSum.2 = dfDegree$SYMBOL
dfDegree$Summary = cvSum[cvSum.2]
write.csv(dfDegree, file='Temp/dfDegree.csv')

dfHub = f_dfGetGeneAnnotation(l$hub)
cvSum = f_csGetGeneSummaryFromGenbank(dfHub$ENTREZID)
cvSum.2 = dfHub$SYMBOL
dfHub$Summary = cvSum[cvSum.2]
write.csv(dfHub, file='Temp/dfHub.csv')

dfBetweenness = f_dfGetGeneAnnotation(l$betweenness)
cvSum = f_csGetGeneSummaryFromGenbank(dfBetweenness$ENTREZID)
cvSum.2 = dfBetweenness$SYMBOL
dfBetweenness$Summary = cvSum[cvSum.2]
write.csv(dfBetweenness, file='Temp/dfBetweenness.csv')

dfCloseness = f_dfGetGeneAnnotation(l$closeness)
cvSum = f_csGetGeneSummaryFromGenbank(dfCloseness$ENTREZID)
cvSum.2 = dfCloseness$SYMBOL
dfCloseness$Summary = cvSum[cvSum.2]
write.csv(dfCloseness, file='Temp/dfCloseness.csv')

# plot largest connected graph - clique
set.seed(1)
ig = plot.graph.clique(oGr)
ig = f_igCalculateVertexSizesAndColors(ig, t(mCounts), fGroups, bColor = F, iSize = 20)
par(mar=c(1,1,1,1)+0.1)
set.seed(1)
plot(ig, vertex.label=NA, layout=layout_with_fr, vertex.frame.color=NA, edge.color='lightgrey')

# plot the largest clique
par(mar=c(1,1,1,1)+0.1)
ig = induced_subgraph(getFinalGraph(oGr), vids = unlist(getLargestCliques(oGr)))
ig = f_igCalculateVertexSizesAndColors(ig, t(mCounts), fGroups, bColor = T)
n = V(ig)$name
lab = f_dfGetGeneAnnotation(n)
V(ig)$label = as.character(lab$SYMBOL)
plot(ig, layout=layout_with_fr)
legend('topright', legend = c('Underexpressed', 'Overexpressed'), fill = c('lightblue', 'pink'))
par(p.old)

## we can look at the problem from the other direction and look at clusters instead of genes
# sample plots
# mean expression of groups in every cluster
plot.mean.expressions(oGr, t(mCounts), fGroups, legend.pos = 'bottomleft', main='Total Change in Each Cluster')
# only significant clusters
#par(mar=c(7, 3, 2, 2)+0.1)
plot.significant.expressions(oGr, t(mCounts), fGroups, main='Significant Clusters', lwd=2, bStabalize = F, cex.axis=0.7)
pr.out = plot.components(oGr, t(mCounts), fGroups, bStabalize = F)
par(mar=c(4,2,4,2))
biplot(pr.out, cex=0.8, cex.axis=0.8, arrow.len = 0)

# plot summary heatmaps
# expression in all clusters
plot.heatmap.all(oGr, t(mCounts))
# marginal expression level in each cluster
plot.heatmap.marginal(oGr, t(mCounts))
#plot.heatmap.significant.clusters(oGr, t(mCounts), fGroups, bStabalize = T)
plot.heatmap.significant.clusters(oGr, t(mCounts), fGroups, bStabalize = F)


# ############################# end analysis
# ## plot heatmap of the community and 
# rownames(mCounts) = gsub('(\\w+).*', '\\1', rownames(mCounts))
# plot.heatmap(oGr, t(mCounts))
# plot.heatmap.means(oGr, t(mCounts))
# plot.mean.expressions(oGr, t(mCounts), fGroups = fSamples, main='Mean expression in each cluster', legend.pos = 'bottomright')
# ## ids for our genes of interest
# ig = getFinalGraph(oGr)
# n = V(ig)$name
# V(ig)[n]$symbol = dfDat[n,'X']
# V(ig)[n]$logFC_abs = abs(dfDat[n,'Fold.change'])
# V(ig)[n]$logFC = dfDat[n,'Fold.change']
# V(ig)$color = ifelse(V(ig)$logFC > 0, 'red', 'blue')
# 
# plot(ig, vertex.label=NA, vertex.size=2, layout=layout.fruchterman.reingold, vertex.frame.color=NA)
# 
# dir.create('Results',showWarnings = F)
# ## export the graphs for use in cytoscape
# write.graph(ig, 'Results/COPD_Data_1.graphml', format = 'graphml')
# 
# df = getClusterMapping(oGr)
# colnames(df) = c('gene', 'cluster')
# df = df[order(df$cluster),]
# write.csv(df, file='Results/COPD_Data_1_expressions_clusters.csv')
# 
# pdf(file='Results/ds_3_expressions_clusters.pdf')
# plot(ig, vertex.label=NA, vertex.size=1, layout=layout.fruchterman.reingold, vertex.frame.color=NA)
# plot(getCommunity(oGr),ig, vertex.label=NA, vertex.size=1, layout=layout.fruchterman.reingold, vertex.frame.color=NA)
# plot.mean.expressions(oGr, t(mCounts), fGroups = fSamples, main='Mean expression in each cluster')
# dev.off(dev.cur())
# 



