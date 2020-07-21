## Heatmap of viral gene vs mock comparison


MatrixA <- c(IAV_cytokines$logFC)
CytokineFC1 <- as.matrix(CytokineFC)
CytokineFC1 [2] <- RSV_cytokines$logFC


CytokineFC1 <- as.data.frame(Cytokines [ ,c(1:2)])
CytokineFC1 [ ,3] <- RSV_cytokines [ ,c(2)]
CytokineFC1 [ ,4] <-  IAV_cytokines [ ,c(2)]
CytokineFC1 [ ,5] <- LungBiopsyCytokines [ ,c(2)]

names(CytokineFC1)[1:5] <- c("Gene", "COV_Vitro", "RSV", "IAV", "COV_Vivo")
CytokineFC1 <- as.data.frame(LogFC_virals [ ,-1])
CytokineFC2 <- CytokineFC1 [ ,-1]

CytokineFC2 <- data.frame(CytokineFC1[,-1], row.names=CytokineFC1[,1])
CytokineFCHM <- as.matrix(CytokineFC2)

library(RColorBrewer)
hmcol<-brewer.pal(11,"RdBu")

heatmap(CytokineFCHM)
coul <- colorRampPalette(brewer.pal(8, "RdBu"))(5000)
heatmap(CytokineFCHM, scale = "column", col = coul,  Colv = NA, na.rm = T, cexRow = 0.5, cexCol = 0.75, xlab = "Lung Tissue Conditions",
        ylab = "Genes", main = "Fold Change")
library(heatmaply)
heatmaply(CytokineFCHM, scale = "column", col = coul, Colv = NA, na.rm = T, cexRow = 0.45, cexCol = 0.75, xlab = "Lung Tissue Conditions",
          ylab = "Genes", main = "Fold Change",  row_dend_left = F)
heatmap(CytokineFCHM, scale = "column", col = coul, Colv = NA, na.rm = T, cexRow = 0.45, cexCol = 0.75, xlab = "Lung Tissue Conditions",
          ylab = "Genes", main = "Fold Change",  row_dend_left = F)
install.packages("viridis")
library("viridis")
heatmap.2(
  CytokineFCHM,
  trace = "none",
  col = viridis(100),
  key = FALSE)
