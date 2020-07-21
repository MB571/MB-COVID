## Lung tissue analysis
## Subset lung tissue data from "counts"
## Subset cytokine genes from lung tissue data

Lung_tissues <- select(counts,contains("Lung"))
names(Lung_tissues)[1:4] <- c("Healthy Lung Biopsy 1", "Healthy Lung Biopsy 2", "Severe COVID19 Lung Biopsy 1", "Severe COVID19 Lung Biopsy 2")
Lung_tissues1 <- rownames_to_column(Lung_tissues)
names(Lung_tissues1)[1] <- "Gene"
Lung_tissues_cytokines <- subset(Lung_tissues1, Gene == "CXCL5" | Gene == "CXCL2"
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
                             | Gene == "IFNLR1"| Gene == "CCL20"| Gene == "EGR1"| Gene == "CCL4"))

## Return Gene ID's as rownames
## Set 0 to "NA" in data frame
Lung_tissues_cytokines1 <- data.frame(Lung_tissues_cytokines[,-1], row.names=Lung_tissues_cytokines[,1])
Lung_tissues_cytokines2 <- Lung_tissues_cytokines %>% replace_with_na_all(condition = ~.x == 0)
Lung_cytokines <- data.frame(Lung_tissues_cytokines2[,-1], row.names=Lung_tissues_cytokines[,1])
names(Lung_cytokines)[1:4] <- c("Healthy1", "Healthy2", "COVID1", "COVID2")
Lung_cytokines1 <- as.matrix(Lung_cytokines)


##Make heatmap of genes
heatmap(Lung_cytokines1)
coul <- colorRampPalette(brewer.pal(8, "Reds"))(5000)
heatmap(Lung_cytokines1, scale = "column", col = coul,  Colv = NA, na.rm = T, cexRow = 0.5, cexCol = 0.75, xlab = "Lung Tissue Conditions",
        ylab = "Genes", main = "Normal vs COVID19 Lung Biopsy")




## Generating LogFC and Avg Expr data for Healthy and Covid Lung tissues
## Differential Gene expression
names(Lung_tissues) <- c("Healthy1", "Healthy2", "COVID1", "COVID2")
snames <- colnames(Lung_tissues)
group <- substr(snames, 1, nchar(snames) - 2)
View(group)
Fileoutput <- "Process.Lung.Tissues"
LT <- DGEList(Lung_tissues, group=group, genes=GeneID)

# filtering to remove low counts

dim(LT)
cutoff <- 1
drop <- which(apply(cpm(LT), 1, max) < cutoff)
d <- LT[-drop,]
dim(d)
AveLogCPM <- aveLogCPM(d)
hist(AveLogCPM)


LT <- calcNormFactors(LT)
LT <- as.matrix(LT)
design <- model.matrix(~0+group)
design
y <- voom(LT, design, plot = T)
fit <- lmFit(y, design,)
#fit <- eBayes(fit)
head(coef(fit))

COVI_vs_Health <- makeContrasts(groupCOVI - groupHealth, levels = colnames(coef(fit)))
tmp <- contrasts.fit(fit, COVI_vs_Health)
tmp <- eBayes(tmp)
top.table <- topTable(tmp, sort.by = "P", n = Inf)
head(top.table, 10)
length(which(top.table$adj.P.Val < 0.05))
top.table$Gene <- rownames(top.table)
top.table <- top.table[,c("Gene", names(top.table)[1:6])]
write.table(top.table, file = paste(Fileoutput,".genes.txt",sep=""), row.names = F, sep = "\t", quote = F)


## Generate AveExpr vs LogFC graphs for Healthy vs COVID biopsies
## For whole gene list
LungBiopsyData <- top.table
LungBiopsyData

ggplot(LungBiopsyData, aes(x=AveExpr, y=logFC, colour=t)) + geom_point() + 
  ggtitle("COVID infected vs Healthy Lung Biopsy") +
  geom_text(aes(label=ifelse(t>12,as.character(Gene),'')), hjust=1.25, vjust=1)

ggplot(LungBiopsyData, aes(x=AveExpr, y=logFC, color=t)) + scale_color_gradient(low="blue", high="red") + geom_point() + 
  ggtitle("COVID infected vs Healthy Lung Biopsy") + 
  geom_text(aes(label=ifelse(t>12,as.character(Gene),'')), hjust=-0.25, vjust=0.3, size = 2)



## Isolate cytokine list and plot
LungBiopsyCytokines <- subset(LungBiopsyData, Gene == "CXCL5" | Gene == "CXCL2"
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


## Plot Lung biopsy cytokines
library(ggrepel)
ggplot(LungBiopsyCytokines, aes(x=AveExpr, y=logFC, color=adj.P.Val)) + scale_color_gradient(low="blue", high="red") + geom_point() + 
  ggtitle("COVID infected vs Healthy Lung Biopsy") + 
  geom_text_repel(aes(label=ifelse(adj.P.Val>0,as.character(Gene),'')),segment.size = 0.25, hjust=-0.25, vjust=0.3, size = 2)

ggplot(Cytokines, aes(x=AveExpr, y=logFC, color=adj.P.Val)) + scale_color_gradient(low="blue", high="red") + geom_point() + 
  ggtitle("COVID infected vs Healthy Lung Biopsy") + 
  geom_text_repel(aes(label=ifelse(adj.P.Val>0,as.character(Gene),'')),segment.size = 0.25, hjust=-0.25, vjust=0.3, size = 2)

