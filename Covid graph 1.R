library(tidyverse)
library(edgeR)
library(limma)
library(dbplyr)


  # Import, read and rename processed file

filename1 <- "process.GSE147507.A549_ACE2.go.txt"
covidfunc <- read.csv(filename1, sep="\t", header=TRUE)
head(covidfunc)

filename <- "process.GSE147507.A549_ACE2.genes.txt"
covid1 <- read.csv(filename, sep="\t", header=TRUE)
head(covid1)


  # Make table of Genes vs log (Fold change) (covid2)

covid2 <- as.data.frame(covid1$Gene) 
covid2$logFC = as.data.frame(covid1$AveExpr)
names(covid2)[1:2] <- c("Gene", "AveExp")
head(covid2)

H <- c(covid2$AveExp [2:8])

(Top_100_Avg_exp_genes_covid)

 ###### Practice above #####

covid <- as.data.frame(covid1)

ggplot(data = covid) + 
  geom_point(mapping = aes(x = AveExpr, y = logFC),
             main="Covid gene expression vs logFC")


## Useful stuff

geom_text(aes(label=Gene), hjust=1.25, vjust=1)
geom_text(aes(label=ifelse(t>26,as.character(Gene),'')))

### practice

ggplot(Cytokines, aes(x=AveExpr, y=LogFC, colour=t)) + geom_point() + 
  geom_text(aes(label=ifelse(t>17,as.character(Gene),'')), hjust=1.25, vjust=1)
  

