## Interferon only response BIOPSY DATA

LungBiopsyInterferon <- subset(LungBiopsyData, Gene == "IFNGR1"
                               | Gene == "IFNAR1"
                               | Gene == "IFNLR1"
                               | Gene == "IFNAR2"
                               | Gene == "IFNG"
                               | Gene == "IFNGR2"
                               | Gene == "IFNE"
                               | Gene == "IFNK"
                               | Gene == "IFNA5"
                               | Gene == "IFNL3"
                               | Gene == "IFNA22P"
                               | Gene == "IFNL4"
                               | Gene == "IFNL2"
                               | Gene == "IFNB1"
                               | Gene == "IFNW1"
                               | Gene == "IFNA21"
                               | Gene == "IFNA4"
                               | Gene == "IFNA7"
                               | Gene == "IFNA10"
                               | Gene == "IFNA16"
                               | Gene == "IFNA17"
                               | Gene == "IFNA14"
                               | Gene == "IFNA13"
                               | Gene == "IFNA2"
                               | Gene == "IFNA8"
                               | Gene == "IFNA1"
                               | Gene == "IFNA6"
                               | Gene == "IFNL1"
                               | Gene == "JAKMIP2"
                               | Gene == "JAKMIP1"
                               | Gene == "JAK1"
                               | Gene == "JAK2"
                               | Gene == "JAK3"
                               | Gene == "JAKMIP3"
                               | Gene == "JAKMIP2-AS1"
                               | Gene == "STAT5A"
                               | Gene == "STAT4"
                               | Gene == "STAT6"
                               | Gene == "STAT1"
                               | Gene == "STAT2"
                               | Gene == "STAT5B"
                               | Gene == "STAT3")

install.packages("ggrepel")
library(ggrepel)
ggplot(LungBiopsyInterferon, aes(x=AveExpr, y=logFC, color=t)) + scale_color_gradient(low="blue", high="red") + geom_point() + 
  ggtitle("Type 1 interferon COVID vs Healthy Biopsy data") + 
  geom_text_repel(aes(label=ifelse(t>-5,as.character(Gene),'')),segment.size = 0.25, hjust=-0.25, vjust=0.3,  size = 2)

## Interferon only response $ VIRAL CONDITIONS DATA

InVitroInterferon <- subset(process.GSE147507.A549_ACE2.genes, Gene == "IFNGR1"
                             | Gene == "IFNAR1"
                             | Gene == "IFNLR1"
                             | Gene == "IFNAR2"
                             | Gene == "IFNG"
                             | Gene == "IFNGR2"
                             | Gene == "IFNE"
                             | Gene == "IFNK"
                             | Gene == "IFNA5"
                             | Gene == "IFNL3"
                             | Gene == "IFNA22P"
                             | Gene == "IFNL4"
                             | Gene == "IFNL2"
                             | Gene == "IFNB1"
                             | Gene == "IFNW1"
                             | Gene == "IFNA21"
                             | Gene == "IFNA4"
                             | Gene == "IFNA7"
                             | Gene == "IFNA10"
                             | Gene == "IFNA16"
                             | Gene == "IFNA17"
                             | Gene == "IFNA14"
                             | Gene == "IFNA13"
                             | Gene == "IFNA2"
                             | Gene == "IFNA8"
                             | Gene == "IFNA1"
                             | Gene == "IFNA6"
                             | Gene == "IFNL1"
                             | Gene == "JAKMIP2"
                             | Gene == "JAKMIP1"
                             | Gene == "JAK1"
                             | Gene == "JAK2"
                             | Gene == "JAK3"
                             | Gene == "JAKMIP3"
                             | Gene == "JAKMIP2-AS1"
                             | Gene == "STAT5A"
                             | Gene == "STAT4"
                             | Gene == "STAT6"
                             | Gene == "STAT1"
                             | Gene == "STAT2"
                             | Gene == "STAT5B"
                             | Gene == "STAT3")

ggplot(InVitroInterferon, aes(x=AveExpr, y=logFC, color=t)) + scale_color_gradient(low="blue", high="red") + geom_point() + 
  ggtitle("Type 1 interferon COVID vs Healthy in vitro") + 
  geom_text_repel(aes(label=ifelse(t>-10,as.character(Gene),'')),segment.size = 0.25, hjust=-0.25, vjust=0.3,  size = 2)


## Interferon only response BIOPSY DATA

LungBiopsyCXC <- subset(LungBiopsyData, Gene == "CCL4"
                        | Gene == "CCL8"
                        | Gene == "CCL3"
                        | Gene == "CCL15-CCL14"
                        | Gene == "CCL14"
                        | Gene == "CCL15"
                        | Gene == "CCL3L3"
                        | Gene == "CCL3L1"
                        | Gene == "CCL4L2"
                        | Gene == "CCL4L1"
                        | Gene == "CCL25"
                        | Gene == "CCL11"
                        | Gene == "CCL7"
                        | Gene == "CCL27"
                        | Gene == "CCL19"
                        | Gene == "CCL1"
                        | Gene == "CCL13"
                        | Gene == "CCL18"
                        | Gene == "CCL16"
                        | Gene == "CCL2"
                        | Gene == "CCL20"
                        | Gene == "CCL26"
                        | Gene == "CCL28"
                        | Gene == "CCL24"
                        | Gene == "CCL17"
                        | Gene == "CCL5"
                        | Gene == "CCL22"
                        | Gene == "CCL21"
                        | Gene == "CCL23"
                        | Gene == "CXCL10"
                        | Gene == "CXCL11"
                        | Gene == "CXCL17"
                        | Gene == "CXCL6"
                        | Gene == "CXCL13"
                        | Gene == "CXCL16"
                        | Gene == "CXCL1"
                        | Gene == "CXCL9"
                        | Gene == "CXCL5"
                        | Gene == "CXCL3"
                        | Gene == "CXCL14"
                        | Gene == "CXCL12"
                        | Gene == "CXCL2")

install.packages("ggrepel")
library(ggrepel)
library(ggplot2)
library(gplots)
library(ggrepel)
ggplot(LungBiopsyCXC, aes(x=AveExpr, y=logFC, color=t)) + scale_color_gradient(low="blue", high="red") + geom_point() + 
  ggtitle("CXC COVID vs Healthy Biopsy") + 
  geom_text_repel(aes(label=ifelse(t>-5,as.character(Gene),'')),segment.size = 0.25, hjust=-0.25, vjust=0.3,  size = 2)

## Interferon only response $ VIRAL CONDITIONS DATA

InVitroCXC <- subset(process.GSE147507.A549_ACE2.genes, Gene == "CCL4"
                            | Gene == "CCL8"
                            | Gene == "CCL3"
                            | Gene == "CCL15-CCL14"
                            | Gene == "CCL14"
                            | Gene == "CCL15"
                            | Gene == "CCL3L3"
                            | Gene == "CCL3L1"
                            | Gene == "CCL4L2"
                            | Gene == "CCL4L1"
                            | Gene == "CCL25"
                            | Gene == "CCL11"
                            | Gene == "CCL7"
                            | Gene == "CCL27"
                            | Gene == "CCL19"
                            | Gene == "CCL1"
                            | Gene == "CCL13"
                            | Gene == "CCL18"
                            | Gene == "CCL16"
                            | Gene == "CCL2"
                            | Gene == "CCL20"
                            | Gene == "CCL26"
                            | Gene == "CCL28"
                            | Gene == "CCL24"
                            | Gene == "CCL17"
                            | Gene == "CCL5"
                            | Gene == "CCL22"
                            | Gene == "CCL21"
                            | Gene == "CCL23"
                            | Gene == "CXCL10"
                            | Gene == "CXCL11"
                            | Gene == "CXCL17"
                            | Gene == "CXCL6"
                            | Gene == "CXCL13"
                            | Gene == "CXCL16"
                            | Gene == "CXCL1"
                            | Gene == "CXCL9"
                            | Gene == "CXCL5"
                            | Gene == "CXCL3"
                            | Gene == "CXCL14"
                            | Gene == "CXCL12"
                            | Gene == "CXCL2")

ggplot(InVitroCXC, aes(x=AveExpr, y=logFC, color=t)) + scale_color_gradient(low="blue", high="red") + geom_point() + 
  ggtitle("CXC COVID vs Healthy in vitro") + 
  geom_text_repel(aes(label=ifelse(t>-10,as.character(Gene),'')),segment.size = 0.25, hjust=-0.25, vjust=0.3,  size = 2)

## Gene ontology

go.fisher.1 <- goana(LungBiopsyCXC, species="Hs", geneid = "genes")
head(topGO(go.fisher.1, sort = "up"),10)
head(topGO(go.fisher.1, sort = "down"),10)
write.table(go.fisher.1, file = paste(Fileoutput,".go.txt",sep=""), row.names = F, sep = "\t", quote = F)

