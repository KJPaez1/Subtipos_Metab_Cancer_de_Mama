
# CLUSTERING
# Basic Machine Learning
# By Rafael Irizarry and Michael Love

## Load packages
library(tissuesGeneExpression)

## data: The data represent RNA expression levels for eight tissues, each with several individuals
data(tissuesGeneExpression)
View(e)
head(e)
dim(e)

## How many tissue types and how many by each type
table(tissue)

## calculate the distances
d = dist(t(e))

## 
library(rafalib)
mypar()

## perform hierarchical clustering
hc <- hclust(d)
hc
plot(hc,labels=tissue,cex=0.5)

## add color (rafalib package)
myplclust(hc, labels=tissue, lab.col=as.fumeric(tissue), cex=0.5)

## Draw a horizontal line at the height we wish to cut
abline(h = 93)
## Examine how the clusters overlap with the actual tissues and see how many samples there are by cluster
hclusters <- cutree(hc, h=120)
table(true=tissue, cluster=hclusters)

## ask cutree() function to give back a given number of clusters. The function then automatically finds the height that results in the requested number of clusters:
hclusters <- cutree(hc, k=8)
table(true=tissue, cluster=hclusters)

### Selecting the number of clusters is generally a challenging step in practice and an active area of research.


# k-means with the first two genes
km <- kmeans(t(e[1:2,]), centers=7)
names(km)

## left: color represent the actual tissues; right: color represents the clusters that were defined by k-means 
par(mfrow = c(1,2))
plot(e[1,], e[2,], col=as.fumeric(tissue), pch=19)
plot(e[1,], e[2,], col=km$cluster, pch=16)

## see from tabulating the results: so so
table(true=tissue,cluster=km$cluster)

### This is very likely due to the fact the the first two genes are not informative regarding tissue type. 
### We can see this in the first plot above. 

## perform k-means clustering using all of the genes. 
km <- kmeans(t(e), centers=7)

## To visualize this, we can use an MDS plot
mds <- cmdscale(d)
par(mfrow = c(1,2))
plot(mds[,1], mds[,2], pch = 19) 
plot(mds[,1], mds[,2], col=km$cluster, pch=16)

## tabulating the results
table(true=tissue,cluster=km$cluster)

table(tissue)

# HEATMAPS
# Basic Machine Learning
# By Rafael Irizarry and Michael Love

## load package for nice color
library(RColorBrewer)
?RColorBrewer

## defining a color palette
hmcol = colorRampPalette(brewer.pal(9, "GnBu"))(100)

## pick the genes with the top variance over all samples
## load package
library(genefilter)

## calculate the variance and sd
?rowVars
rv = rowVars(e)
View(rv)

## ???
?order
idx = order(-rv)[1:40]
View(idx)

## load gplot package to make heatmap
library(gplots)
library(rafalib)
cols = palette(brewer.pal(8, "Dark2"))[as.fumeric(tissue)]

## samples indicated with one color
head(cbind(colnames(e), cols))

## heatmap
heatmap.2(e[idx,], labCol=tissue,
          trace="none", 
          ColSideColors=cols, 
          col=hmcol

















