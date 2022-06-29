# Clustering Analysis 

## Optimizes graphical parameters and add color by name
library(rafalib)


GSE31448_Subtypes = c(
)
## Optimizes graphical parameters
mypar()

## make plot
GSE58812_Peng_normalized = as.matrix(GSE58812_Peng_normalized) ## 
dGSE58812_Peng_normalized = dist(t(GSE58812_Peng_normalized))
hcGSE58812_Peng_normalized <- hclust(dGSE58812_Peng_normalized)
hcGSE58812_Peng_normalized
plot(hcGSE58812_Peng_normalized, xlab = NULL, cex=0.3, hang = -1)

## Adding labels and color
myplclust(hcGSE58812_Peng_normalized, labels = GSE31448_Subtypes, lab.col=as.fumeric(GSE31448_Subtypes), cex=0.8)

## 
library(dendextend)
dendGSE58812_Peng_normalized = as.dendrogram(hcGSE58812_Peng_normalized)
dendGSE58812_Peng_normalized = color_labels(hcGSE58812_Peng_normalized, 8, col= c("magenta4", "blue4", "firebrick2","hotpink", "yellow4", "darkgoldenrod", "mediumorchid4", "chartreuse4" ), cex = 0)

plot(dendGSE58812_Peng_normalized)
abline(h = 6.6, col = "yellow2", lwd = 1)
#------------------------------------
dendGSE2071_Rosario = color_labels(hc_GSE58812, 1:14)
#------------------------------------
## Adding labels and color
myplclust(hcGSE58812_Peng_normalized, labels=GSE31448_Subtypes, lab.col=as.fumeric(GSE31448_Subtypes), cex=0.9)

## Draw a horizontal line at the height we wish to cut and see the numbers of clusters 
abline(h = 93)
hclusters <- cutree(hc, h=120)
table(true=tissue, cluster=hclusters)

## defining clusters
hclusters <- cutree(hc, k=8)
table(true=tissue, cluster=hclusters)










# ------------------------ NOTES ------------------------------------------



# make density plot
# use file devices to make plots
pdf(file = "myplot.pdf")
with(faithful, plot(eruptions, waiting))
title(main = "Old Faithful Geyser data")
dev.off()
# make cluster and PCA too
# when find diferentially expressed genes, make boxplots



# ------------------------ NOTES ------------------------------------------

