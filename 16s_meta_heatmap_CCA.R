library(pheatmap)
library(openxlsx)
library(vegan)
library(RColorBrewer)
library(dplyr)
library(Hmisc)
library(reshape2)
library(factoextra)
library(readxl)
library(corrplot)
library(igraph)
library(psych)

CorMatrix <- function(cor, p) {
  ut <- upper.tri(cor)
  data.frame(
    row = rownames(cor)[row(cor)[ut]],
    column = rownames(cor)[col(cor)[ut]],
    cor = (cor)[ut],
    p = p[ut]
  )
}

##Genus
genus.table <- read.xlsx("HG_vs_MG.genus_diff.xlsx", sheet = 1, startRow = 1, colNames = T, rowNames = F)
row.names(genus.table) <- genus.table$X1
x1 <- genus.table[genus.table$pValue <= 0.05, -c(1:8)]
## meta
meta.table <- read.xlsx("CG_vs_PG.meta_diff.xlsx", sheet = 1, startRow = 1, colNames = T, rowNames = F)
row.names(meta.table) <- meta.table$name
x2 <- meta.table[, -c(1:2)]

## hclust heatmap
res.hc <- eclust(x1, "hclust", k = 2, graph = TRUE, hc_metric = "spearman", hc_method = "ward.D")
pdf("Genus cluster.pdf", width = 18, height = 10)
fviz_dend(res.hc,
  k = 2,
  cex = 1,
  k_colors = c("#2E9FDF", "#00AFBB", "#E7B800", "#FC4E07"),
  color_labels_by_k = TRUE, horiz = T,
  rect = TRUE
)
dev.off()
# Genus_group <- cbind(res.hc$labels, res.hc$cluster)
# corheatmap
r <- rcorr(x = t(x1), y = t(x2), type = "spearman")
Gmeta_cor <- melt(r$r[-c(1:nrow(x1)), 1:nrow(x1)])
Gmeta_P <- melt(r$P[-c(1:nrow(x1)), 1:nrow(x1)])
cor_value <- cbind(Gmeta_cor, Gmeta_P$value)
colnames(cor_value) <- c("meta", "Genus", "R", "P")
write.csv(cor_value, "Genus Correlation analysis.csv")
cor_Matrices <- as.data.frame(r$r[-c(1:nrow(x1)), 1:nrow(x1)])
P_Matrices <- r$P[-c(1:nrow(x1)), 1:nrow(x1)]
{
  P_Matrices[which(P_Matrices >= 0.05)] <- c("")
  P_Matrices[which(P_Matrices >= 0.01 & P_Matrices < 0.05)] <- c("*")
  P_Matrices[which(P_Matrices >= 0.001 & P_Matrices < 0.01)] <- c("**")
  P_Matrices[which(as.numeric(P_Matrices) < 0.001)] <- c("***")
}
windowsFonts(A = windowsFont("Times New Roman"), B = windowsFont("Arial Black"), C = windowsFont("Comic Sans MS"))
cor_Matrices <- as.data.frame(t(cor_Matrices))
pheatmap(cor_Matrices,
  # kmeans_k = 3,
  clustering_distance_rows = "correlation",
  clustering_method = "ward.D",
  scale = "none",
  cluster_row = TRUE, cluster_col = TRUE,
  # annotation_row=annotation_col,
  legend = TRUE,
  display_numbers = t(P_Matrices), ## 设置显著性符号
  cellwidth = 16, cellheight = 12,
  angle_col = 45,
  # color=colorRampPalette(c("green","black","red"))(1000),
  color = colorRampPalette(c("steelblue", "white", "hotpink"))(1000),
  fontsize = 8, filename = "corrplot.tiff"
)

# CCA
vare.cca <- cca(t(x1), t(x2))
pdf("CCA.pdf",width =10,height= 10)
plot(vare.cca)
dev.off()
summary(vare.cca)

#  corrplot
pdf(
  file = "corrplot.pdf",
  height = unit(min(20, log10(ncol(r$r)) * 12), "cm"),
  width = unit(min(20, log10(ncol(r$r)) * 12), "cm"),
  pointsize = min(25, 25 - log10(ncol(r$r)))
)
corrplot(r$r, type = "lower", order = "hclust", tl.col = "black", tl.srt = 45, tl.cex = 0.5)
dev.off()


if(0) {
# RDA
dune.Manure <- rda(dune ~ Manure, dune.env)
plot(dune.Manure)

# Co-occurrence network

occor = corr.test(t(x1),use="pairwise",method="spearman",adjust="fdr",alpha=.05)
occor.r = r$r
occor.p = r$P
# 确定物种间存在相互作用关系的阈值，将相关性R矩阵内不符合的数据转换为0
occor.r[occor.p>0.05|abs(occor.r)<0.6] = 0 


# 构建igraph对象
igraph = graph_from_adjacency_matrix(occor.r,mode="undirected",weighted=TRUE,diag=FALSE)
igraph
# NOTE:可以设置weighted=NULL,但是此时要注意此函数只能识别相互作用矩阵内正整数，所以应用前请确保矩阵正确。
# 可以按下面命令转换数据
# occor.r[occor.r!=0] = 1
# igraph = graph_from_adjacency_matrix(occor.r,mode="undirected",weighted=NULL,diag=FALSE)

# 是否去掉孤立顶点，根据自己实验而定
# remove isolated nodes，即去掉和所有otu均无相关性的otu 可省略，前期矩阵已处理过
bad.vs = V(igraph)[degree(igraph) == 0]
igraph = delete.vertices(igraph, bad.vs)
igraph
# 将igraph weight属性赋值到igraph.weight
igraph.weight = E(igraph)$weight

# 做图前去掉igraph的weight权重，因为做图时某些layout会受到其影响
E(igraph)$weight = NA
set.seed(123)
plot(igraph, layout=layout_with_kk)
plot(igraph, layout=layout_with_fr)
plot(igraph, layout=layout_as_tree)
plot(igraph,main="Co-occurrence network",vertex.frame.color=NA,vertex.label=NA,edge.width=1,
     vertex.size=5,edge.lty=1,edge.curved=TRUE,margin=c(0,0,0,0))
plot(igraph,main="Co-occurrence network",vertex.frame.color=NA,vertex.label=NA,
     edge.lty=1,edge.curved=TRUE,margin=c(0,0,0,0))
}