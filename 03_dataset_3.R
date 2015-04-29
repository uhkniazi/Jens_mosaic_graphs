#NAME: 03_dataset_3.R
#DESC: input the gene lists from the csv file and create reactome based graph and CGraphClust based heatmap
#AUTH: u.niazi@imperial.ac.uk
#Date: 27/04/2015

## source libraries
source('~/Dropbox/Home/Data/R/My_Libraries/NGS_functions.R')
source('../CGraphClust/CGraphClust.R')
library(igraph)

# databases 
library(org.Hs.eg.db)
library(GO.db)
library(reactome.db)

p.old = par()
## data loading and formatting
dfDat = read.csv(file.choose(), header=T, stringsAsFactors=F)
## name of PKM2 replaced by PKM
## STIM1 doesn't have reactome term, use STIM2

# get gene symbols as they have an underscore before them
s = gsub('^(\\w+)_\\d+', replacement = '\\1', x = dfDat$X)
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
fSamples = gsub('\\w(\\d).*', '\\1', c)
dfCounts.p.o = dfCounts.p[,order(fSamples)]
fSamples = fSamples[order(fSamples)]
# subsample only samples of type 2 and 3
mCounts = dfCounts.p.o#[,fSamples != '2']
mCounts = t(mCounts)
mCor = cor(mCounts)
hist(mCor)

## create the clustering graph
oGr = CGraphClust(dfGraph, mCor)#, iCorCut = 0.7)

## plot heatmap of the community and 
rownames(mCounts) = gsub('\\w(\\d).*', '\\1', rownames(mCounts))
plot.heatmap(oGr, t(mCounts))
plot.heatmap.means(oGr, t(mCounts))

## ids for our genes of interest
ig = getFinalGraph(oGr)
n = V(ig)$name
V(ig)[n]$symbol = dfDat[n,'X']
V(ig)[n]$logFC_abs = abs(dfDat[n,'Fold.change'])
V(ig)[n]$logFC = dfDat[n,'Fold.change']
V(ig)$color = ifelse(V(ig)$logFC > 0, 'red', 'blue')

pkm = which(V(ig)$symbol == 'PKM')
pkm = V(ig)[pkm]$name

stim1 = which(V(ig)$symbol == 'STIM1')
stim1 = V(ig)[stim1]$name

p = get.all.shortest.paths(ig, pkm, stim1)

ivEdges = numeric()
for (o in 1:length(p$res)){
  iFirst = p$res[[o]][1]
  iLast = p$res[[o]][length(p$res[[o]])]
  # make an edge list
  for (i in 1:(length(p$res[[o]])-1) ) {
    ivE = c(p$res[[o]][i], p$res[[o]][i+1])
    ivEdges = c(ivEdges, ivE)
  } # inner for
} # outer for
# get the ids for these edges
e = get.edge.ids(ig, ivEdges)
E(ig)$color = 'grey'
E(ig)[e]$color = 'red'

plot(ig, vertex.label=NA, vertex.size=1, layout=layout.fruchterman.reingold, vertex.frame.color=NA)

dir.create('Results',showWarnings = F)
## export the graphs for use in cytoscape
write.graph(ig, 'Results/ds_3_pkm_stim.graphml', format = 'graphml')


