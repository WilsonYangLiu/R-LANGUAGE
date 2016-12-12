rm(list = ls())
#setwd(dir = "/media/wilson/b776f228-366c-4e52-acd6-65df5b458e8c/Project_protiomics/proteinGroup/")
setwd(dir = "C:/Users/Wei-Xin/Downloads/Temp/PCA - R/")
if(!require("FactoMineR")){
  install.packages("FactoMineR")
  library("FactoMineR")
}
if(!require("corrplot")){
  install.packages("corrplot")
  library("corrplot")
}
if(!require("PerformanceAnalytics")){
  install.packages("PerformanceAnalytics")
  library("PerformanceAnalytics")
}
if(!require("factoextra")){
  install.packages("factoextra")
  library("factoextra")
}
if (!require("RColorBrewer")) {
  install.packages("RColorBrewer", dependencies = TRUE)
  library(RColorBrewer)
}

if (!require("gplots")) {
  install.packages("gplots", dependencies = TRUE)
  library(gplots)
}

colName <- c("Protein IDs","Reverse","Potential contaminant","mix_1369B","1369Bp","1322Bp","1589Bp","B5Bg","A9Bp","A13Bp","B13Bg","B15Bg","A19Bp","mix_1369R","1369Rp","1746Rg","B5Rg","B13Rg","A19Rp","1368Bp","A15Bp","na","na","mix_1411B","1411Bp","1617Bg","1820Bg","1679Bg","1385Bg","1746Bg","A5Bp","B9Bg","B19Bg","mix_1411R","1411Rp","1617Rg","1368Rp","A5Rp","A9Rp","B9Rg","B15Rg","B19Rg","na","mix_1820R","1820Rg","1322Rp","1679Rg","1385Rg","1589Rp","A13Rp","A15Rp","na","na")
Data <- read.csv(file = "data.csv", header = TRUE, col.names = colName)
filterData <- Data[which((Data$Reverse != "+")&(Data$Potential.contaminant != "+")), which(colName != "na")]
#write.csv(x = filterData, file = "filterData.csv", row.names = FALSE)

# number of valid data in every lane
# see: http://www.plob.org/2012/08/31/3297.html
tmp <- colName[-c(1:3)]
no0_all <- apply(X = filterData[, -c(1:3)], MARGIN = 2, FUN = function(x) {sum(x != 0)})
bp <- barplot(height = no0_all,width = 1, space = 1, axisnames = FALSE, main = "number of valid data in every lane", xlab = "case", ylab = "number", xpd = TRUE)
axis(1, at = bp, labels = FALSE)
text(x = seq(1.5, 89.5, by = 2), y = -300, srt = 45, adj = 1, labels = tmp[c(tmp != "na")], xpd = TRUE, cex = 0.5)

# proteins that have valid data in all 45 lanes
idx <- c(rep(FALSE, 3), rep(TRUE, 45))
no0 <- apply(X = filterData, MARGIN = 1, FUN = function(x) {all(as.numeric(x[idx]) != 0)})
no0_filterData <- filterData[which(no0), -c(2, 3)]
rowName <- filterData[which(no0), 1]
no0_filterData_rowname <-  `rownames<-`(no0_filterData[, -c(1)], rowName)
#write.csv(x = no0_filterData, file = "no0_filterData.csv", row.names = FALSE)

# divide by mix, and then perform pairwise.t.test, calculate p value
Cancer_1369B <- no0_filterData_rowname[, 1:10] / no0_filterData_rowname$mix_1369B
Cancer_1369R <- no0_filterData_rowname[, 11:18] / no0_filterData_rowname$mix_1369R
Cancer_1411B <- no0_filterData_rowname[, 19:28] / no0_filterData_rowname$mix_1411B
Cancer_1411R <- no0_filterData_rowname[, 29:37] / no0_filterData_rowname$mix_1411R
Cancer_1820R <- no0_filterData_rowname[, 38:45] / no0_filterData_rowname$mix_1820R
# t-test
Cancer <- cbind(Cancer_1369B[-c(1)], Cancer_1369R[-c(1)], Cancer_1411B[-c(1)], Cancer_1411R[-c(1)], Cancer_1820R[-c(1)])
#write.csv(x = Cancer, file = "Cancer.csv")
pvals <- apply(Cancer, 1, function(x) {t.test(x[c(1:9, 15:25)],x[c(10:14, 26:40)],paired = TRUE)$p.value})
#write.table(rowName[which(pvals < 0.05)], file = "sigProtein.txt", quote = FALSE, row.names = FALSE, col.names = FALSE)
adj_pvals <- p.adjust(p = pvals, method = "fdr", n = length(pvals))
#write.table(rowName[which(adj_pvals < 0.05)], file = "sigProtein_fdr.txt", quote = FALSE, row.names = FALSE, col.names = FALSE)

# normalization by row (individuals) (20 to 20, wrong)
#idxB <- c(1:9, 15:25)
#idxR <- c(10:14, 26:40)
#CancerMatrix_pri <- apply(as.matrix(Cancer), MARGIN = 1, FUN = function(test) {c((test[idxB] - mean(test[idxB])) / sd(test[idxB]), (test[idxR] - mean(test[idxR])) / sd(test[idxR]))})
#CancerMatrix <- t(CancerMatrix_pri)

CancerMatrix <- t(scale(as.matrix(Cancer))) # normalization by columns in Cancer
#write.csv(x = CancerMatrix, file = "CancerMatrix.csv")
#CancerMatrix <- as.matrix(read.csv(file = "CancerMatrix.csv", row.names = 1))

# correlation analysis
source(file = "pca_own.R")
# basic info
CancerMatrix_stats <- round(Data_stats(X = CancerMatrix), 2)
# Correlation between variables
CancerMatrix_cor <- round(cor(CancerMatrix), 2)
# Visualize the correlation matrix, only the first 100
corrplot(CancerMatrix_cor[1:100, 1:100], type="upper",tl.pos = "td", tl.cex = 0.2, order="hclust", tl.col="black")
# Make a scatter plot matrix showing the correlation coefficients between variables and the significance levels
# only the first 5
chart.Correlation(CancerMatrix[, 1:5], histogram=TRUE, pch=19)
# Visualize correlation matrix using a heatmap, only the first 100
col<- colorRampPalette(c("blue", "white", "red"))(20)
heatmap(x = CancerMatrix_cor[1:100, 1:100], col = col, symm = TRUE, cexRow = 0.3, cexCol = 0.3)

# perform pca
CancerMatrix.pca <- PCA(CancerMatrix, scale.unit = TRUE, ncp = 5, graph = FALSE)
# Variances of the principal components
eigenvalues <- CancerMatrix.pca$eig
barplot(eigenvalues[, 2], names.arg=1:nrow(eigenvalues), 
        main = "Variances",
        xlab = "Principal Components",
        ylab = "Percentage of variances",
        col ="steelblue")
# Add connected line segments to the plot
lines(x = 1:nrow(eigenvalues), eigenvalues[, 2], 
      type="b", pch=19, col = "red")
# Graph of individus and variables
position <- c(rep(1, 9), rep(0, 5), rep(1, 11), rep(0, 15))
efficacy <- c(rep(1, 3), rep(0, 1), rep(1, 2), rep(0, 2), rep(1, 2), rep(0, 3), rep(1, 4), rep(0, 5), rep(1, 1), rep(0, 2), rep(1, 1), rep(0, 1), rep(1, 3), rep(0, 4), rep(1, 1), rep(0, 2), rep(1, 3))
fviz_pca_ind(X = CancerMatrix.pca, label = "all", habillage = as.matrix(position), addEllipses = TRUE) + theme_minimal()
fviz_pca_ind(X = CancerMatrix.pca, label = "all", habillage = as.matrix(efficacy), addEllipses = TRUE) + theme_minimal()


# heat maps between individual and varialbe
# creates a own color palette from red to green
my_palette <- colorRampPalette(c("green", "yellow", "red"))(n = 29)
# defines the color breaks manually for a "skewed" color transition
# col_breaks = c(seq(-1,0,length=100),               # for red
#                seq(0,0.8,length=100),              # for yellow
#                seq(0.8,1,length=100))              # for green
png("./heatmaps_in_r.png",      
    width = 5*300,        # 5 x 300 pixels
    height = 5*300,
    res = 300,            # 300 pixels per inch
    pointsize = 8)        # smaller font size
heatmap.2(CancerMatrix[,],
          main = "Correlation", # heat map title
          notecol="black",      # change font color of cell labels to black
          density.info="none",  # turns off density plot inside color legend
          trace="none",         # turns off trace lines inside the heat map
          scale = "column",
          margins =c(12,9),     # widens margins around plot
          col=my_palette,       # use on color palette defined earlier
          dendrogram="both",
          hclustfun = hclust,
          labCol = FALSE)    # only draw a row dendrogram
dev.off()
