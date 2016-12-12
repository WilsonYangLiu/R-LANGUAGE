rm(list = ls())
setwd(dir = "D:/Download/")

if(!require("pamr")){
  install.packages("pamr")
  library("pamr")
}
# data("khan")
# 
# rowName <- as.vector(khan$X1[-1])
# khan <- `rownames<-`(khan, khan$X)
# khan <- khan[-1, -c(1, 2)]
# write.csv(x = khan, file = "D:/Download/khan.csv")
khan <- read.csv(file = "D:/Download/khan.csv", row.names = 1)

label <- c(rep(1, 23), rep(2, 20), rep(3, 12), rep(4, 8))
mydata <- list(x=data.matrix(khan), y=factor(label))
khan.train <- pamr.train(mydata)
# Cross-validate the classifier
khan.results<- pamr.cv(khan.train, mydata)
# Plot the cross-validated error curves
png(file = "plotcv.png")
pamr.plotcv(khan.results)
dev.off()
# Compute the confusion matrix for a particular model (threshold=4.0)
#pamr.confusion(khan.results, threshold=4.0)
# Plot the cross-validated class probabilities by class
png(file = "plotcvprob.png")
pamr.plotcvprob(khan.results, mydata, threshold=4.0)
dev.off()
# Plot the class centroids
pdf(file = "plotcen.pdf")
mydata_1 <- list(x=data.matrix(khan), y=label, genenames=rowName)
pamr.plotcen(khan.train, mydata_1, threshold=4.0)
dev.off()
# Make a gene plot of the most significant genes
mydata_1 <- list(x=data.matrix(khan), y=label)
pamr.geneplot(khan.train, mydata_1, threshold=5.3)
# Estimate false discovery rates and plot them
mydata_1 <- list(x=data.matrix(khan), y=factor(label), geneid=as.character(1:nrow(x)), genenames=paste("g",rowName,sep=""))
khan.train <- pamr.train(mydata_1)
fdr.obj<- pamr.fdr(khan.train, mydata_1)
pamr.plotfdr(fdr.obj)
# List the significant genes
mydata_1 <- list(x=data.matrix(khan), y=factor(label), geneid=as.character(1:nrow(x)), genenames=paste("g",rowName,sep=""))
pamr.listgenes(khan.train, mydata_1, threshold=4.0)
# Try heterogeneity analysis, with class "BL" taken to be the normal group
mydata_1 <- list(x=data.matrix(khan), y=factor(label), geneid=as.character(1:nrow(x)), genenames=paste("g",rowName,sep=""))
khan.train2 <- pamr.train(mydata_1, hetero=4)
khan.results2 <- pamr.cv(khan.train2, mydata_1)
png(file = "plotcv_hetero.png")
pamr.plotcv(khan.results2)
dev.off()
# Look for better threshold scalings
mydata_1 <- list(x=data.matrix(khan), y=factor(label), geneid=as.character(1:nrow(x)), genenames=paste("g",rowName,sep=""))
khan.scales <- pamr.adaptthresh(khan.train)
khan.train3 <- pamr.train(mydata_1, threshold.scale=khan.scales)
khan.results3 <- pamr.cv(khan.train3, mydata_1)
png(file = "plotcv_scale.png")
pamr.plotcv(khan.results3)
dev.off()

