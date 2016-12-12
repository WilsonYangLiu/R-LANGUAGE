rm(list = ls())
#setwd(dir = "D:/Download/PAM - R/")

if(!require("pamr")){
  install.packages("pamr")
  library("pamr")
}

ProteinGroup <- t(read.csv(file = "ProteinGroup.csv", row.names = 1))

rowName <- rownames(ProteinGroup)
label <- c(rep(1, 3), rep(0, 1), rep(1, 2), rep(0, 2), rep(1, 2), rep(0, 3), rep(1, 4), rep(0, 5), rep(1, 1), rep(0, 2), rep(1, 1), rep(0, 1), rep(1, 3), rep(0, 4), rep(1, 1), rep(0, 2), rep(1, 3))
#label <- c(rep(1, 9), rep(0, 5), rep(1, 11), rep(0, 15))
mydata <- list(x=data.matrix(ProteinGroup), y=factor(label))
khan.train <- pamr.train(mydata)
# Cross-validate the classifier
khan.results<- pamr.cv(khan.train, mydata)
# Plot the cross-validated error curves
png(file = "plotcv.png")
pamr.plotcv(khan.results)
dev.off()
# Plot the cross-validated class probabilities by class
png(file = "plotcvprob.png")
pamr.plotcvprob(khan.results, mydata, threshold=0.8)
dev.off()
# Plot the class centroids
pdf(file = "plotcen.pdf")
mydata_1 <- list(x=data.matrix(ProteinGroup), y=label, genenames=rowName)
pamr.plotcen(khan.train, mydata_1, threshold=0.8)
dev.off()




