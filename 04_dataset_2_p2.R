#NAME: 04_dataset_2_p2.R
#DESC: input the gene lists from the csv file and create reactome based graph
#AUTH: u.niazi@imperial.ac.uk
#Date: 29/04/2015


## source libraries
source('~/Dropbox/Home/Data/R/My_Libraries/NGS_functions.R')
source('CGraphClust.R')
library(igraph)

# databases 
library(org.Hs.eg.db)
library(GO.db)
library(reactome.db)

p.old = par()

## data loading and formatting
dfDat = read.csv(file.choose(), header=T, stringsAsFactors=F)

# get gene symbols
#s = gsub('\'(\\w+)\'', replacement = '\\1', x = dfDat$X)
s = gsub('\'(.+)\'', replacement = '\\1', x = dfDat$X)
dfDat$X = s
s = gsub('(.+)_\\d+', replacement = '\\1', x = dfDat$X)
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
dfCounts.p = dfCounts[,2:140]
mCounts = t(dfCounts.p)
mCor = cor(mCounts)
hist(mCor)
# there is not much happening if we look at all the data, correlations show a unimodal distribution
# create a factor to subsample the data
c = colnames(dfCounts.p)
# fSamples = gsub('X\\.(Flu\\.Pos|Flu\\.Neg|Healthy).+', '\\1', c)
# fSamples = factor(fSamples, levels = c('Healthy', 'Flu.Neg', 'Flu.Pos'))
fSamples = gsub('X(1|0).*', '\\1', c)
fSamples = factor(fSamples, levels = c(1, 0), labels= c('Flu', 'FluLike'))
dfCounts.p.o = dfCounts.p[,order(fSamples)]
fSamples = fSamples[order(fSamples)]
# subsample only samples of type 2 and 3
mCounts = dfCounts.p.o#[,fSamples != '2']
mCounts = t(mCounts)
mCor = cor(mCounts)
hist(mCor)
mCor = abs(mCor)

## create the clustering graph
oGr = CGraphClust(dfGraph, mCor, iCorCut = 0.3)
rownames(mCounts) = fSamples
plot.heatmap(oGr, t(mCounts))
plot.heatmap.means(oGr, t(mCounts))

## add annotation to graph before saving
ig = getFinalGraph(oGr)
n = V(ig)$name
V(ig)[n]$symbol = dfDat[n,'X']
V(ig)[n]$logFC_abs = abs(dfDat[n,'fold.change'])
V(ig)[n]$logFC = dfDat[n,'fold.change']
V(ig)$color = ifelse(V(ig)$logFC > 0, 'red', 'blue')

plot(ig, vertex.label=NA, vertex.size=1, layout=layout.fruchterman.reingold, vertex.frame.color=NA)
plot(getCommunity(oGr), ig, vertex.label=NA, vertex.size=1, layout=layout.fruchterman.reingold, vertex.frame.color=NA)

dir.create('Results/DS_2/New',showWarnings = F)
## export the graphs for use in cytoscape
write.graph(ig, 'Results/DS_2/New/ds_2_new.graphml', format = 'graphml')
