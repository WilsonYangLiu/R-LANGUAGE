rm(list = ls())
setwd(dir = "C:/Users/Wei-Xin/Downloads/Temp/PCA - R/")
if(!require("FactoMineR")){
  install.packages("FactoMineR")
  library("FactoMineR")
}
if(!require("factoextra")){
  install.packages("factoextra")
  library("factoextra")
}
if(!require("corrplot")){
  install.packages("corrplot")
  library("corrplot")
}
if(!require("PerformanceAnalytics")){
  install.packages("PerformanceAnalytics")
  library("PerformanceAnalytics")
}

data("iris")

#########################
# correlation analysis ##
#########################
source(file = "pca_own.R")
# basic info
iris_stats <- round(Data_stats(X = iris[, -5]), 2)
# Correlation between variables
iris_cor <- round(cor(iris[, -5]), 4)
# Visualize the correlation matrix
corrplot(iris_cor, type="upper",tl.pos = "td", tl.cex = 0.2, order="hclust", tl.col="black")
# Make a scatter plot matrix showing the correlation coefficients between variables and the significance levels
chart.Correlation(iris[, -5], histogram=TRUE, pch=19)
# Visualize correlation matrix using a heatmap
col<- colorRampPalette(c("blue", "white", "red"))(20)
heatmap(x = iris_cor, col = col, symm = TRUE, cexRow = 0.3, cexCol = 0.3)

########
# PCA ##
########
iris.pca <- PCA(X = iris[, -5], graph = FALSE)
# Extract and visualize eigenvalues/variances
get_eig(X = iris.pca)
fviz_eig(X = iris.pca, addlabels = TRUE, hjust = -0.5) + theme_minimal()

# Extract and visualize results for variables
var <- get_pca_var(res.pca = iris.pca); var
fviz_pca_var(X = iris.pca, col.var = "contrib") + scale_color_gradient2(low="white", mid="blue", high="red", midpoint = 96) + theme_minimal()
## Variable contributions to the principal axes
fviz_contrib(X = iris.pca, choice = "var", axes = 1)

# Extract and visualize results for individuals
ind <- get_pca_ind(res.pca = iris.pca); ind
## Use repel = TRUE to avoid overplotting
fviz_pca_ind(X = iris.pca, repel = TRUE, col.ind = "cos2") + scale_color_gradient2(low="blue", mid="white", high="red", midpoint=0.6)+ theme_minimal()
## Color by groups: habillage=iris$Species
p <- fviz_pca_ind(X = iris.pca, geom = "point", habillage = iris$Species, addEllipses = TRUE) + theme_minimal()
print(p)

# Biplot of individuals and variables
fviz_pca_biplot(X = iris.pca, label = "var", habillage = iris$Species, addEllipses = TRUE) + theme_minimal()

############
# cluster ##
############
# k-means
set.seed(123)
# Determine the optimal number of clusters
fviz_nbclust(x = scale(iris[, -5]), FUNcluster = kmeans, method = "gap_stat")
iris.km <- kmeans(x = scale(iris[, -5]), centers = 3, nstart = 25)
fviz_cluster(object = iris.km, data = scale(iris[, -5])) + theme_minimal()+ scale_color_manual(values = c("#00AFBB","#2E9FDF", "#E7B800", "#FC4E07"))+ scale_fill_manual(values = c("#00AFBB","#2E9FDF", "#E7B800", "#FC4E07")) + labs(title= "Partitioning Clustering Plot")
# Hierarchical clustering
iris.hc <- hcut(x = iris[, -5], k = 3, stand = TRUE)
fviz_dend(x = iris.hc, k = 3, rect = TRUE, k_colors = c("#00AFBB","#2E9FDF", "#E7B800", "#FC4E07"), cex = 0.5)


