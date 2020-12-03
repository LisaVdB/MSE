.libPaths("C:/Program Files (x86)/R/win-library/3.5")
library(RColorBrewer)
library(gplots)
library(gdata)
library("openxlsx")
library("dynamicTreeCut")
library(reshape2)
library(ggplot2)
library(scales)
library("wesanderson")
library(tidyverse)
library(pheatmap)

setwd("d:/Post-doc/Project_TomatoPepper_Ryan/Tomato Pepper Data/MSE")
setwd("c:/Users/lisav/Downloads")

data = read.xlsx("FPKM_average_nomalized_lowlyexpressedremoved.xlsx", colNames = TRUE, rowNames = TRUE)
head(data)

data_zero = data
data_zero[data_zero == 0] <- .000001
P = data_zero / rowSums(data_zero) # P is the relative expression of each gene in each timepoint
H = rowSums(-P * log2(P))

#data.max = apply(data_zero, 1, FUN = max)

hist(H)
Th = max(H)*0.7 #30% entropy is allowed to be considered an outlier
SE.genes = data_zero[which(H < Th), ]
SE.genes.safe = SE.genes

# Initialize an empty matrix to store the outliers for each gene
outliers.matrix <- matrix(0, nrow(SE.genes), ncol(SE.genes))

for (gene in 1:nrow(SE.genes)) {
  # Normalize the expression pattern of each gene
  E_norm_unsorted = (SE.genes[gene, ] - mean(as.numeric(SE.genes[gene, ]))) / sd(SE.genes[gene, ])
  # Sort the expression of each gene, from the smallest to the largest in each condition
  # Having the expressions sorted faciliates the calculation of the minimum U
  E_norm = unlist(as.matrix(sort.int(as.numeric(E_norm_unsorted), index.return = TRUE))[1])
  # idx is a vector with the indexes of the sorted expression normalized (E_norm)
  idx = unlist(as.matrix(sort.int(as.numeric(E_norm_unsorted), index.return = TRUE))[2])
  
  s.max = 4 # Max num of times that we want to consider as outliers for each gene
  
  #  Initialize an empty matrix to store the outliers for each gene before choosing the ones that lead to the minimum U
  outliers <- matrix(0, (2 * s.max + 1), ncol(SE.genes))
  #  Initialize U
  U <- matrix(0, (2 * s.max + 1), 1)
  
  for (s in 1:s.max) {
    n = length(E_norm) - s
    E_u = E_norm[(s + 1):length(E_norm)]
    outliers[s, idx[1:s]] <- -1
    U[s] = n * log10(sd(E_u)) + sqrt(2) * s * factorial(log10(n)) / n
  }
  for (s in 1:s.max) {
    n = length(E_norm) - s
    E_u = E_norm[1:(length(E_norm) - s)]
    outliers[(s.max + s), idx[(length(E_norm) - s + 1):length(E_norm)]] <- 1
    U[s.max+s] = n * log10(sd(E_u)) + sqrt(2) * s * factorial(log10(n)) / n
  }
  n = length(E_norm) - s.max
  E_u = E_norm[(s.max / 2):(length(E_norm) - s.max / 2)]
  outliers[(2 * s.max + 1), idx[1]] <- -1
  outliers[(2 * s.max + 1), idx[length(E_norm)]] <- 1
  U[2 * s.max + 1] = n * log10(sd(E_u)) + sqrt(2) * s.max * factorial(log10(n)) / n
  outliers.matrix[gene, ] = outliers[which(U == min(U)), ]
}

# Convert matrix of outliers to dataframe and merge with the table 
outliers.matrix = data.frame(outliers.matrix)
colnames(outliers.matrix) = colnames(data)

# Merge 
SE.genes <- cbind(SE.genes, outliers.matrix)
write.table(SE.genes, file = "MSE_output_Normalized.txt", sep = "\t", row.names = TRUE)

# Heatmap
colors = wes_palette("Zissou1", 21, type = "continuous")
group = c("PT_1d", "PT_3d", "PT_5d", "TP_1d","TP_3d","TP_5d", "TT/PP_1d", "TT/PP_3d", "TT/PP_5d")
data_norm = data.frame()
for (i in 1:nrow(SE.genes.safe)) {
  norm = (SE.genes.safe[i,1:9] - min(SE.genes.safe[i,1:9])) / (max(SE.genes.safe[i,1:9]) - min(SE.genes.safe[i,1:9]))
  data_norm = rbind(data_norm, norm)
}
map = pheatmap(data_norm, color = colors, cluster_distance_rows = "correlation", show_rownames = FALSE, 
               angle_col = 90, labels_col = group)
map
ggsave('heatmap_MSE_TomPep.png', plot = map, width = 20, height = 15, units = "cm", limitsize = FALSE, dpi = 600)


#Other heatmap with prior hierarchical clustering
clusters = hclust(dist(data_norm), method = "complete")
clusters_split = cutreeDynamic(clusters, distM = as.matrix(dist(data_norm)), method = "hybrid", deepSplit = 1)
clusters_split = rbind.data.frame(as.matrix(clusters_split))
rownames(clusters_split) = rownames(data_norm)
clusters_split = rownames_to_column(clusters_split, var = "Var1")

df_molten = data_norm[clusters$order,]
df_molten = melt(as.matrix(df_molten))
df_molten = join(df_molten, clusters_split, by = "Var1")

colors2 = wes_palette("Royal1")

# Create theme
newtheme <- theme(axis.ticks = element_blank(), axis.text.x = element_text(size = 10),
                  axis.text.y = element_blank(),
                  legend.text = element_text(size = 10), legend.title = element_text(size = 10, face = "bold"),
                  plot.title = element_text(size = 12, face = "bold", hjust = 0.5),
                  strip.text.x = element_text(size = 6, face = "bold"),
                  strip.text.y = element_text(size = 10, face = "bold"))
# Create plot
p <- ggplot(df_molten, aes(x = Var2, y = Var1, fill = value)) + geom_tile()
p <- p + newtheme + labs(x = "Time points", y = "", fill = "Scaled FPKM")
p <- p + scale_fill_gradientn(colors = colors, na.value = colors2[1]) + ggtitle("MSE selected genes")
p <- p + scale_x_discrete(labels = c(unique(group)))
#p = p + facet_wrap(~ V1, scales = "free")
p
ggsave('heatmap_MSE_TomPep_b.png', plot = p, width = 30, height = 20, units = "cm", limitsize = FALSE, dpi = 600)

#Heatmap grouped according to MSE (6 groups)
data = read.xlsx("MSE_Output_selection_norm.xlsx", colNames = TRUE)
data_norm = data.frame()
for (i in 1:nrow(data)) {
  norm = (data[i,2:10] - min(data[i,2:10])) / (max(data[i,2:10]) - min(data[i,2:10]))
  data_norm = rbind(data_norm, norm)
}

clusters_split = cutreeDynamic(hclust(dist(data_norm[,2:10]), method = "complete"), distM = as.matrix(dist(data_norm[,2:10])), method = "hybrid", deepSplit = 1)
clusters_split = rbind.data.frame(as.matrix(clusters_split))

data_all = cbind(data[,1], data[,20], clusters_split, data_norm)
colnames(data_all)[c(1,2,3)] = c("GeneID", "Group", "Cluster")
data_all = data_all[order(data_all$Cluster),]

kmeans = kmeans(data[,11:19], 20,iter.max = 50, nstart = 20)
data_all2 = cbind(data[,1], data[,20], kmeans$cluster, data_norm)
colnames(data_all2)[c(1,2,3)] = c("GeneID", "Group", "Cluster")
data_all2 = data_all2[order(data_all$Cluster),]

data_all3 = cbind(data[,1], data[,20], kmeans$cluster, data[,11:19])
colnames(data_all3)[c(1,2,3)] = c("GeneID", "Group", "Cluster")
data_all3 = data_all3[order(data_all3$Cluster),]

data_long = gather(data_all2, Sample, Value, -GeneID, -Group, -Cluster)
p <- ggplot(data_long, aes(x = Sample, y = GeneID, fill = Value)) + geom_tile()
p <- p + newtheme + labs(x = "", y = "", fill = "Scaled FPKM")
p <- p + scale_fill_gradientn(colors = colors) + ggtitle("MSE selected genes")
p <- p + scale_x_discrete(labels = c(unique(group)))
p = p + facet_wrap(~ Group, scales = "free")
p
ggsave('heatmap_MSE_TomPep_c.png', plot = p, width = 30, height = 20, units = "cm", limitsize = FALSE, dpi = 600)
