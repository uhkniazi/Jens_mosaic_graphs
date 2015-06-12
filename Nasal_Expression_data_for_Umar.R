#NAME: Nasal_Expression_data_for_Umar.R
#DESC: input the gene lists from the csv file and create reactome based graph and CGraphClust based heatmap
#AUTH: u.niazi@imperial.ac.uk
#Date: 11/06/2015

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

# remove duplicated genes
f = !duplicated(dfDat$Names)
dfDat = dfDat[f,]
# get the enterez ids of the genes
dfEnt = select(org.Hs.eg.db, keys = dfDat$Names, columns = 'ENTREZID', keytype = 'SYMBOL')
# remove the duplicated symbols (one instance of it)
f = !duplicated(dfEnt$SYMBOL)
dfEnt = dfEnt[f,]
# match the 2 data frames and get enterez ids
i = match(dfDat$Names, dfEnt$SYMBOL)
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
dfCounts.p = dfCounts[,2:67]
# create count matrix with genes as columns i.e column vectors to create correlation matrix
mCounts = t(dfCounts.p)
mCor = cor(mCounts)
hist(mCor, prob=T, main='correlation matrix')

# order the count matrix according to sequence of samples
c = rownames(mCounts)
fSamples = gsub('(\\w+).*', '\\1', c)
mCounts = mCounts[order(fSamples),]
fSamples = fSamples[order(fSamples)]
rownames(mCounts) = fSamples

## create the clustering graph
oGr = CGraphClust(dfGraph, abs(mCor), iCorCut = 0.2)#, iCorCut = 0.7)

## look at the clusters in the graph first
ig = getFinalGraph(oGr)
vcount(ig)
ecount(ig)
plot(ig, vertex.label=NA, vertex.size=2, layout=layout.fruchterman.reingold, vertex.frame.color=NA)

## plot heatmap of the community and 
#rownames(mCounts) = gsub('(\\w+).*', '\\1', rownames(mCounts))
plot.heatmap(oGr, t(mCounts))
plot.heatmap.means(oGr, t(mCounts))
plot.mean.expressions(oGr, t(mCounts), fGroups = fSamples, main='Mean expression in each cluster', legend.pos = 'topright', lwd=2)
## ids for our genes of interest
ig = getFinalGraph(oGr)
n = V(ig)$name
V(ig)[n]$symbol = dfDat[n,'Names']
V(ig)[n]$logFC_abs = abs(dfDat[n,'Fold.change'])
V(ig)[n]$logFC = dfDat[n,'Fold.change']
V(ig)$color = ifelse(V(ig)$logFC > 0, 'red', 'blue')

par(mar=c(1,1,1,1)+0.1)
plot(ig, vertex.label=NA, vertex.size=2, layout=layout.fruchterman.reingold, vertex.frame.color=NA)
par(p.old)
dir.create('Results',showWarnings = F)
## export the graphs for use in cytoscape
write.graph(ig, 'Results/Nasal_Expression_data_for_Umar.graphml', format = 'graphml')

df = getClusterMapping(oGr)
colnames(df) = c('gene', 'cluster')
df = df[order(df$cluster),]
write.csv(df, file='Results/Nasal_Expression_data_for_Umar.csv')

pdf(file='Results/ds_3_expressions_clusters.pdf')
plot(ig, vertex.label=NA, vertex.size=1, layout=layout.fruchterman.reingold, vertex.frame.color=NA)
plot(getCommunity(oGr),ig, vertex.label=NA, vertex.size=1, layout=layout.fruchterman.reingold, vertex.frame.color=NA)
plot.mean.expressions(oGr, t(mCounts), fGroups = fSamples, main='Mean expression in each cluster')
dev.off(dev.cur())



