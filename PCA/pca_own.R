if(!require("Hmisc")){
  install.packages("Hmisc")
  library("Hmisc")
}

# ++++++++++++++++++++++++++++
# PCA.own
# ++++++++++++++++++++++++++++
# Performs a principal components analysis on the given data matrix 
# Returns the results of the new Components
# ----------------------------
# X : matrix of oringin data
# center: a logical value indicating whether the variables should be shifted to be zero centered.
# scale.: a logical value indicating whether the variables should be scaled to have unit variance before the analysis takes place.
pca.own <- function(X, center = TRUE, scale. = TRUE){
  # preprocess
  if(scale.){
    Datas <- scale(X, center = TRUE, scale = TRUE)
  } else{
    if(center){
      Datas <- scale(X, center = TRUE, scale = FALSE)
    }
  }
  # Covariance matrix
  Datas.cov <- round(cov(Datas), 6)
  # eigenvectors and the eigenvalues of the covariance matrix
  Datas.eig <- eigen(Datas.cov) # Rate-limiting step
  # the new dataset
  newData <- Datas %*% Datas.eig$vectors[, 1:40]
  colnames(newData) <- paste(c(rep("PC", 40)), as.character(c(1:40)), sep = '')

  return(newData)
}

# ++++++++++++++++++++++++++++
# eig_stats
# ++++++++++++++++++++++++++++
# Return the eingenvalue, variance and cumulative variance.
# Note that the function summary() to extract the eigenvalues and variances from an object of class prcomp
# You can also use the package factoextra to do the same thing.
# ----------------------------
# sdev : sdev <- {prcomp}$sdev. the standard deviations of the principal components comming from the function prcomp or princomp
#       (i.e., the square roots of the eigenvalues of the covariance/correlation matrix, 
#       though the calculation is actually done with the singular values of the data matrix).
eig_stats <- function(sdev){
  # Eigenvalues
  eig <- (sdev)^2
  # Variances in percentage
  variance <- eig*100/sum(eig)
  # Cumulative variances
  cumvar <- cumsum(variance)
  eig_stats <- data.frame(eig = eig, variance = variance,
                        cumvariance = cumvar)
  return(eig_stats)
}

# ++++++++++++++++++++++++++++
# pca_var
# ++++++++++++++++++++++++++++
# Return
# var.cor (var.coord): the coordinates of variables on the principal components.
# var.cos2: quality of representation for variables
# var.contrib: contributions of the variables to the principal components
# ----------------------------
# loading : loadings <- {prcomp}$rotation
# sdev: sdev <- {prcomp}$sdev
pca_var <- function(loading, sdev){
  var.coord <- var.cor <- t(apply(loadings, 1, FUN = function(var.loading, comp.sdev){var.loadings * comp.sdev}, sdev))
  var.cos2 <- var.coord^2
  comp.cos2 <- apply(var.cos2, 2, sum)
  var.contrib <- t(apply(var.cos2, 1, FUN = function(var.cos2, comp.cos2){var.cos2*100 / comp.cos2}, comp.cos2))
  return(list(var.cor, var.cos2, var.contrib))
}

getdistance <- function(ind_row, center, scale){
  return(sum(((ind_row-center)/scale)^2))
}
cos2 <- function(ind.coord, d2){
  return(ind.coord^2/d2)
}
contrib <- function(ind.coord, comp.sdev, n.ind){
  100*(1/n.ind)*ind.coord^2/comp.sdev^2
}
# ++++++++++++++++++++++++++++
# pca_ind
# ++++++++++++++++++++++++++++
# Return
# var.cos2: quality of representation for variables
# var.contrib: contributions of the variables to the principal components
# ----------------------------
# X : oringinal DataSets
# ind.coord: ind.coord <- {prcomp}$x
# sdev: sdev <- {prcomp}$sdev
# center: center <- {prcomp}$center
# scale: scale <- {prcomp}$scale
pac_ind <- function(X, ind.coord, sdev, center, scale){
  d2 <- apply(X, 1, getdistance, center, scale)
  ind.cos2 <- apply(ind.coord, 2, cos2, d2)
  ind.contrib <- t(apply(ind.coord, 1, contrib, sdev, nrow(ind.coord)))
  return(list(ind.cos2, ind.contrib))
}

# ++++++++++++++++++++++++++++
# Data_stats
# ++++++++++++++++++++++++++++
# A simple function for retaining the basic statistics of variables containing 6 columns:
# minimun, first quartile, median, mean, Third quartile and maxium
# ----------------------------
# X : matrix of oringin data
Data_stats <- function(X){
  # the statistic for variables (columns)
  return(data.frame(
    Min = apply(X, 2, min), # minimum
    Q1 = apply(X, 2, quantile, 1/4), # First quartile
    Med = apply(X, 2, median), # median
    Mean = apply(X, 2, mean), # mean
    Q3 = apply(X, 2, quantile, 3/4), # Third quartile
    Max = apply(X, 2, max) # Maximum
  ))
}

# ++++++++++++++++++++++++++++
# flattenCorrMatrix
# ++++++++++++++++++++++++++++
# A simple function for formatting a correlation matrix into a table with 4 columns containing :
# Column 1 : row names (variable 1 for the correlation test)
# Column 2 : column names (variable 2 for the correlation test)
# Column 3 : the correlation coefficients
# Column 4 : the p-values of the correlations
# optional: use Hmisc R package to calculate the correlation and correlation p-values: rcorr()
# ----------------------------
# cormat : matrix of the correlation coefficients
# pmat : matrix of the correlation p-values
flattenCorrMatrix <- function(cormat, pmat) {
  ut <- upper.tri(cormat)
  return(data.frame(
    row = rownames(cormat)[row(cormat)[ut]],
    column = rownames(cormat)[col(cormat)[ut]],
    cor  =(cormat)[ut],
    p = pmat[ut]
  ))
}

# ++++++++++++++++++++++++++++
# cor.mtest
# ++++++++++++++++++++++++++++
# Computing the p-value of correlations (matrix).
# Retuen the matrix of p-value.
# ----------------------------
# mat : is a matrix of data
# ... : further arguments to pass to the native R cor.test function
cor.mtest <- function(mat, ...) {
  mat <- as.matrix(mat)
  n <- ncol(mat)
  p.mat<- matrix(NA, n, n)
  diag(p.mat) <- 0
  for (i in 1:(n - 1)) {
    for (j in (i + 1):n) {
      tmp <- cor.test(mat[, i], mat[, j], ...)
      p.mat[i, j] <- p.mat[j, i] <- tmp$p.value
    }
  }
  colnames(p.mat) <- rownames(p.mat) <- colnames(mat)
  return(p.mat)
}