rm(list = ls())
setwd(dir = "~/Downloads/Temp/PCA/")
#priDatas <- read.csv(file = "4017.csv", header = TRUE)
priDatas <- read.csv(file = "CancerMatrix.csv", row.names = 1)
#Datas <- scale(as.matrix(priDatas))
Datas <- (as.matrix(priDatas))

source(file = "pca_own.R")
newData <- pca.own(X = Datas)

# plot
plot(-60:60, -60:60, type = "n", xlab = "PC1", ylab = "PC2")
points(newData[c(1:9, 15:25), 1], newData[c(1:9, 15:25), 2], col = "red")
points(newData[c(10:14, 26:40), 1], newData[c(10:14, 26:40), 2], col = "blue")
