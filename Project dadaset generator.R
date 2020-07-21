library(tidyverse)
library(edgeR)
library(limma)
library(dbplyr)


## Make counts2 and counts3 (counts & counts1 with gene row)

counts2 <- rownames_to_column(counts1)
names(counts2)[1] <- "Gene"
counts3 <- rownames_to_column(counts)
names(counts3)[1] <- "Gene"



## Make "Cytokines" (Expression data from "process.GSE147507.A549_ACE2") on genes of interest

Cytokines <- subset(process.GSE147507.A549_ACE2.genes, Gene == "CXCL5" | Gene == "CXCL2"
                    | Gene == "CXCL3"
                    | Gene == "CXCL1"
                    | Gene == "CXCL5"
                    | Gene == "CXCR4"
                    | Gene == "CXCL16"
                    | Gene == "CXCL12"
                    | Gene == "IL5"
                    | Gene == "IL6"
                    | Gene == "IL11"
                    | Gene == "IL12A"
                    | Gene == "IL15"
                    | Gene == "IL16"
                    | Gene == "IL21R-AS1"
                    | Gene == "IL23A"
                    | Gene == "IL32"
                    | Gene == "TRAFD1"
                    | Gene == "TRAF4"
                    | Gene == "TRAF6"
                    | Gene == "TRAF1"
                    | Gene == "TRAF3IP2"
                    | Gene == "TRAF7"
                    | Gene == "TRAF3IP2-AS1"
                    | Gene == "TRAF3"
                    | Gene == "TRAF2"
                    | Gene == "TRAF3IP1"
                    | Gene == "TRAF5"
                    | Gene == "TNFAIP3"
                    | Gene == "TNFAIP8"
                    | Gene == "TNF"
                    | Gene == "TNFAIP2"
                    | Gene == "TNFRSF21"
                    | Gene == "TNFRSF12A"
                    | Gene == "TNFAIP6"
                    | Gene == "TNFSF14"
                    | Gene == "TNFSF9"
                    | Gene == "TNFRSF10B"
                    | Gene == "TNFRSF1A"
                    | Gene == "C1QTNF6"
                    | Gene == "TNFRSF1B"
                    | Gene == "C1QTNF1"
                    | Gene == "TNFRSF19"
                    | Gene == "TNFSF15"
                    | Gene == "TNFAIP8L1"
                    | Gene == "TNFRSF10C"
                    | Gene == "TNFRSF25"
                    | Gene == "TNFRSF10A"
                    | Gene == "TNFSF12"
                    | Gene == "TNFAIP1"
                    | Gene == "TNFRSF11A"
                    | Gene == "TNFRSF11B"
                    | Gene == "TNFSF10"
                    | Gene == "TNFAIP8L3"
                    | Gene == "TNFRSF10D"
                    | Gene == "TNFRSF9"
                    | Gene == "IFNA22P"
                    | Gene == "IFNAR1"
                    | Gene == "IFNAR2"
                    | Gene == "IFNE"
                    | Gene == "IFNGR1"
                    | Gene == "IFNGR2"
                    | Gene == "IFNLR1"| Gene == "CCL20"| Gene == "EGR1"| Gene == "CCL4")


## Plot an expression map of Cytokines relevant to "covid1" (covid1 = process gene exp)

plot(Cytokines$AveExpr, Cytokines$LogFC)


##   Make this look pretty (include gene labels on graph for highest and lowest
##   LogFC and AveExpr values)


plot(Covid1$AveExpr, Covid1$logFC)


## Subsetting gene data in "counts" 
## Make rows selectable in "counts" and name as "Gene" 

AllTissuesCytokines <- subset(counts3, Gene == "CXCL5" | Gene == "CXCL2"
                         | Gene == "CXCL3"
                         | Gene == "CXCL1"
                         | Gene == "CXCL5"
                         | Gene == "CXCR4"
                         | Gene == "CXCL16"
                         | Gene == "CXCL12"
                         | Gene == "IL5"
                         | Gene == "IL6"
                         | Gene == "IL11"
                         | Gene == "IL12A"
                         | Gene == "IL15"
                         | Gene == "IL16"
                         | Gene == "IL21R-AS1"
                         | Gene == "IL23A"
                         | Gene == "IL32"
                         | Gene == "TRAFD1"
                         | Gene == "TRAF4"
                         | Gene == "TRAF6"
                         | Gene == "TRAF1"
                         | Gene == "TRAF3IP2"
                         | Gene == "TRAF7"
                         | Gene == "TRAF3IP2-AS1"
                         | Gene == "TRAF3"
                         | Gene == "TRAF2"
                         | Gene == "TRAF3IP1"
                         | Gene == "TRAF5"
                         | Gene == "TNFAIP3"
                         | Gene == "TNFAIP8"
                         | Gene == "TNF"
                         | Gene == "TNFAIP2"
                         | Gene == "TNFRSF21"
                         | Gene == "TNFRSF12A"
                         | Gene == "TNFAIP6"
                         | Gene == "TNFSF14"
                         | Gene == "TNFSF9"
                         | Gene == "TNFRSF10B"
                         | Gene == "TNFRSF1A"
                         | Gene == "C1QTNF6"
                         | Gene == "TNFRSF1B"
                         | Gene == "C1QTNF1"
                         | Gene == "TNFRSF19"
                         | Gene == "TNFSF15"
                         | Gene == "TNFAIP8L1"
                         | Gene == "TNFRSF10C"
                         | Gene == "TNFRSF25"
                         | Gene == "TNFRSF10A"
                         | Gene == "TNFSF12"
                         | Gene == "TNFAIP1"
                         | Gene == "TNFRSF11A"
                         | Gene == "TNFRSF11B"
                         | Gene == "TNFSF10"
                         | Gene == "TNFAIP8L3"
                         | Gene == "TNFRSF10D"
                         | Gene == "TNFRSF9"
                         | Gene == "IFNA22P"
                         | Gene == "IFNAR1"
                         | Gene == "IFNAR2"
                         | Gene == "IFNE"
                         | Gene == "IFNGR1"
                         | Gene == "IFNGR2"
                         | Gene == "IFNLR1"| Gene == "CCL20"| Gene == "EGR1"| Gene == "CCL4")

CovidCytokines <- subset(counts2, Gene == "CXCL5" | Gene == "CXCL2"
                              | Gene == "CXCL3"
                              | Gene == "CXCL1"
                              | Gene == "CXCL5"
                              | Gene == "CXCR4"
                              | Gene == "CXCL16"
                              | Gene == "CXCL12"
                              | Gene == "IL5"
                              | Gene == "IL6"
                              | Gene == "IL11"
                              | Gene == "IL12A"
                              | Gene == "IL15"
                              | Gene == "IL16"
                              | Gene == "IL21R-AS1"
                              | Gene == "IL23A"
                              | Gene == "IL32"
                              | Gene == "TRAFD1"
                              | Gene == "TRAF4"
                              | Gene == "TRAF6"
                              | Gene == "TRAF1"
                              | Gene == "TRAF3IP2"
                              | Gene == "TRAF7"
                              | Gene == "TRAF3IP2-AS1"
                              | Gene == "TRAF3"
                              | Gene == "TRAF2"
                              | Gene == "TRAF3IP1"
                              | Gene == "TRAF5"
                              | Gene == "TNFAIP3"
                              | Gene == "TNFAIP8"
                              | Gene == "TNF"
                              | Gene == "TNFAIP2"
                              | Gene == "TNFRSF21"
                              | Gene == "TNFRSF12A"
                              | Gene == "TNFAIP6"
                              | Gene == "TNFSF14"
                              | Gene == "TNFSF9"
                              | Gene == "TNFRSF10B"
                              | Gene == "TNFRSF1A"
                              | Gene == "C1QTNF6"
                              | Gene == "TNFRSF1B"
                              | Gene == "C1QTNF1"
                              | Gene == "TNFRSF19"
                              | Gene == "TNFSF15"
                              | Gene == "TNFAIP8L1"
                              | Gene == "TNFRSF10C"
                              | Gene == "TNFRSF25"
                              | Gene == "TNFRSF10A"
                              | Gene == "TNFSF12"
                              | Gene == "TNFAIP1"
                              | Gene == "TNFRSF11A"
                              | Gene == "TNFRSF11B"
                              | Gene == "TNFSF10"
                              | Gene == "TNFAIP8L3"
                              | Gene == "TNFRSF10D"
                              | Gene == "TNFRSF9"
                              | Gene == "IFNA22P"
                              | Gene == "IFNAR1"
                              | Gene == "IFNAR2"
                              | Gene == "IFNE"
                              | Gene == "IFNGR1"
                              | Gene == "IFNGR2"
                              | Gene == "IFNLR1"
                              | Gene == "CCL20"| Gene == "EGR1"| Gene == "CCL4")
                         

##
## Mapping cytokines in scatter plot of AveExpr vs LogFC

ggplot(Cytokines, aes(x=AveExpr, y=logFC, colour=t)) + geom_point() + 
  ggtitle("Cytokine Fold Change vs Average Expression in COVID") +
  geom_text(aes(label=ifelse(t>17,as.character(Gene),'')), hjust=1.25, vjust=1)


ggplot(Cytokines, aes(x=AveExpr, y=logFC, color=t)) + scale_color_gradient(low="blue", high="red") + geom_point() + 
  ggtitle("Cytokine Fold Change vs Average Expression in COVID") + 
  geom_text_repel(aes(label=ifelse(logFC>-3,as.character(Gene),'')),segment.size = 0.25, hjust=-0.25, vjust=0.3, size = 2)


## Mapping Total AveExpr vs LogFC

ggplot(Covid1, aes(x=AveExpr, y=logFC, color=t)) + scale_color_gradient(low="black", high="red") + geom_point() + 
  ggtitle("Total Fold Change vs Average Expression in COVID") + 
  geom_text_repel(aes(label=ifelse(logFC>5,as.character(Gene),'')),segment.size = 0.25, hjust=-0.3, vjust=0.3, size = 1.75) 
  
## Heatmapping Cov vs Control
# Fontsize as variable
fontsize <- 0.5


# Plotting 
CytokineHeatmap <- as.matrix(CovidCytokines[ ,c(2:7)])
CytokineHeatmap <- `rownames<-`(CovidCytokines, CovidCytokines$Gene)
CytokineHeatmap <- CytokineHeatmap[,-1]
Geneinfo <- data.frame(Gene = CovidCytokines$Gene)
library(ComplexHeatmap)
Heatmap(FinalCytoHM, cluster_columns = F, 
        row_names_side = "left",
        row_names_gp = gpar(cex = fontsize))

##Include NA values
df3 <- as.data.frame(CytokineHeatmap)
df4 <- df3 %>% replace_with_na_all(condition = ~.x == 0)
CytHM <- as.matrix(df4)
FinalCytoHM <- `rownames<-`(CytHM, CovidCytokines$Gene)


##Plot the final cytokine heatmap

FinalCytoHM <- as.matrix(FinalCytoHM)
coul <- colorRampPalette(brewer.pal(8, "Reds"))(5000)

heatmap(CytokineHeatmap, scale = "column", col = coul,  Colv = NA,  na.rm = T,
        cexRow = 0.5, cexCol = 0.75, xlab = "Conditions", ylab = "Genes", main = "In vitro Covid vs Control")



## Heatmapping Cov vs Control vs Influenza




