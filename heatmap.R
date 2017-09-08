# Interactive heatmaps
if (!require("devtools")) install.packages("devtools")
devtools::install_github("rstudio/d3heatmap")

library("d3heatmap")
data("mtcars")
d3heatmap(scale(mtcars), colors = "RdBu",
          k_row = 4, k_col = 2)

# Using `dendextend` to enhance heatmaps
library(dendextend)
# order for rows
Rowv  <- mtcars %>% scale %>% dist %>% hclust %>% as.dendrogram %>%
  set("branches_k_color", k = 3) %>% set("branches_lwd", 1.2) %>%
  ladderize
# Order for columns
# We must transpose the data
Colv  <- mtcars %>% scale %>% t %>% dist %>% hclust %>% as.dendrogram %>%
  set("branches_k_color", k = 2, value = c("orange", "blue")) %>%
  set("branches_lwd", 1.2) %>%
  ladderize

heatmap(scale(mtcars), Rowv = Rowv, Colv = Colv, scale = "none")

library(gplots)
heatmap.2(scale(mtcars), scale = "none", col = bluered(100), 
          Rowv = Rowv, Colv = Colv,
          trace = "none", density.info = "none")

library("d3heatmap")
d3heatmap(scale(mtcars), colors = "RdBu", Rowv = Rowv, Colv = Colv)

# Complex heatmap
if (!require("devtools")) install.packages("devtools")
devtools::install_github("jokergoo/ComplexHeatmap")

library(ComplexHeatmap)
df <- scale(mtcars)
library(circlize)
mycol <- colorRamp2(c(-2, 0, 2), c("blue", "white", "red"))

Heatmap(df, name = "mtcars", col = mycol,
        column_title = "Column title",
        column_title_gp = gpar(fontsize = 14, fontface = "bold"),
        row_title = "Row title",
        row_title_gp = gpar(fontsize = 14, fontface = "bold"))

set.seed(2)
Heatmap(df, name = "mtcars", col = mycol, split = 2) # split into 2 groups

df <- t(df)
# Annotation data frame
annot_df <- data.frame(cyl = mtcars$cyl, am = mtcars$am,  mpg = mtcars$mpg)
# Define colors for each levels of qualitative variables
# Define gradient color for continuous variable (mpg)
col = list(cyl = c("4" = "green", "6" = "gray", "8" = "darkred"),
           am = c("0" = "yellow", "1" = "orange"),
           mpg = colorRamp2(c(17, 25), c("lightblue", "purple"))
)
# Create the heatmap annotation
ha <- HeatmapAnnotation(annot_df, col = col, show_legend = FALSE)
# Combine the heatmap and the annotation
Heatmap(df, name = "mtcars", col = mycol, top_annotation = ha)

library("GetoptLong")
for(an in colnames(annot_df)) {
  seekViewport(qq("annotation_@{an}"))
  grid.text(an, unit(1, "npc") + unit(2, "mm"), 0.5,
            default.units = "npc", just = "left")
}

# Define some graphics to display the distribution of columns
.hist = anno_histogram(df, gp = gpar(fill = "lightblue"))
.density = anno_density(df, type = "line", gp = gpar(col = "blue"))
ha_mix_top = HeatmapAnnotation(hist = .hist, density = .density)
# Define some graphics to display the distribution of rows
.violin = anno_density(df, type = "violin", 
                       gp = gpar(fill = "lightblue"), which = "row")
.boxplot = anno_boxplot(df, which = "row")
ha_mix_right = HeatmapAnnotation(violin = .violin, bxplt = .boxplot,
                                 which = "row", width = unit(4, "cm"))
# Combine annotation with heatmap
Heatmap(df, name = "mtcars", col = mycol,
        column_names_gp = gpar(fontsize = 8),
        top_annotation = ha_mix_top, 
        top_annotation_height = unit(4, "cm")) + ha_mix_right




seekViewport("plot")
grid.rect()
grid.points(x, y)
grid.xaxis()
grid.yaxis()

seekViewport("margin1")
grid.text("Random X", y = unit(1, "lines"))

seekViewport("margin2")
grid.text("Random Y", x = unit(1, "lines"), rot = 90)
