# Name: Severity_data_V2.R
# Auth: u.niazi@imperial.ac.uk
# Date: 21/07/2015
# Desc: analysis of expression data from mosaic study for trevor


# import the data
source('../CGraphClust/CGraphClust.R')
library(reactome.db)
library(org.Hs.eg.db)

p.old = par()

## data loading and formatting
dfDat = read.csv(file.choose(), header=F, stringsAsFactors=F)

# convert first row to factors
fGroups = dfDat[1,2:(ncol(dfDat))]
fGroups = unlist(fGroups)
names(fGroups) = NULL
fGroups = factor(fGroups)

dfDat = dfDat[-1,]
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
rownames(dfDat) = dfDat$ENTREZID
# get the respective reactome ids
dfGraph = select(reactome.db, dfDat$ENTREZID, 'REACTOMEID', 'ENTREZID')
dfGraph = na.omit(dfGraph)

## create correlation matrix
n = unique(dfGraph$ENTREZID)
# some genes do not have reactome ids
dfCounts = dfDat[rownames(dfDat) %in% n,]
# remove the extra columns
dfCounts = dfCounts[,2:140]
# convert to count matrix
mCounts = t(dfCounts)
# order the count matrix before making heatmaps or plots
rownames(mCounts) = fGroups
mCounts = mCounts[order(fGroups),]
fGroups = fGroups[order(fGroups)]

# try stabalizing the data before proceeding
mCounts = t(mCounts)
mCounts = t(apply(mCounts, 1, function(x) f_ivStabilizeData(x, fGroups)))
colnames(mCounts) = fGroups
mCounts = t(mCounts)
# create correlation matrix
mCor = cor(mCounts)
hist(sample(mCor, 1000, replace = F), prob=T, xlab='', ylab='', main='Correlation of Genes')

# create the graph cluster object
# using absolute correlation vs actual values lead to different clusters
oGr = CGraphClust(dfGraph, abs(mCor), iCorCut = 0.6)
#oGr = CGraphClust(dfGraph, (mCor), iCorCut = 0.7)

# sample plots
# mean expression in every cluster
plot.mean.expressions(oGr, t(mCounts), fGroups, legend.pos = 'topright')
# only significant clusters
plot.significant.expressions(oGr, t(mCounts), fGroups, lwd=2)
# only one cluster
plot.cluster.expressions(oGr, t(mCounts), fGroups, csClustLabel = '1430728', main='1430728 Metabolism - oxidation')
plot.cluster.expressions(oGr, t(mCounts), fGroups, csClustLabel = '109582', main='109582 Hemostasis')
plot.cluster.expressions(oGr, t(mCounts), fGroups, csClustLabel = '168256', main='168256 Immune System')
# without data stabalization
pr.out = plot.components(oGr, t(mCounts), fGroups, bStabalize = F)
biplot(pr.out, cex=0.6, arrow.len=0)

# with some stabalization
pr.out = plot.components(oGr, t(mCounts), fGroups, bStabalize = T)
biplot(pr.out, cex = 0.5, arrow.len = 0)

# graph properties
n = getLargestCliques(oGr)
n = unique(names(unlist(n)))
f_dfGetGeneAnnotation(cvEnterezID = n)
plot.graph.clique(obj = oGr)
# final graph
plot.final.graph(oGr)

# top vertices based on centrality
mCent = mPrintCentralitySummary(oGr)
# top 2% of the vertices from each category and largest clique
l = lGetTopVertices(oGr)
l = unique(unlist(l))
f_dfGetGeneAnnotation(cvEnterezID = l)

# plot summary heatmaps
plot.heatmap.all(oGr, t(mCounts))
plot.heatmap.means(oGr, t(mCounts))

ig = getFinalGraph(oGr)
par(mar=c(1,1,1,1)+0.1)
plot(ig, vertex.label=NA, vertex.size=2, layout=layout_with_fr(ig, weights = E(ig)$ob_to_ex), vertex.frame.color=NA)
plot(getCommunity(oGr), ig, vertex.label=NA, vertex.size=2, layout=layout.fruchterman.reingold, vertex.frame.color=NA)
par(p.old)

df = getClusterMapping(oGr)
colnames(df) = c('gene', 'cluster')
df = df[order(df$cluster),]

write.csv(df, 'Temp/clusters.csv')