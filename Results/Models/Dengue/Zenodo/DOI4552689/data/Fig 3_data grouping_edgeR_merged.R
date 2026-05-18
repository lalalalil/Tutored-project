library(tidyverse)
library(edgeR)
rawdata <- read.csv("/Volumes/Home/sadiehl/MedMyDocs/DENV2-RNASeq/VIGR/2019.11.04 Korin merged/genecounts.csv", check.names=FALSE, stringsAsFactors=FALSE)
colnames(rawdata)

rawdata2 <- rawdata %>% select(Gene.name, Geneid, S004_d180, S004_d188,S004_d208, S006_d180,S006_d188,S006_d208,S010_d180,S010_d188,S010_d208, 
                               S011_d180,S011_d188,S011_d208, S012_d180,S012_d188,S012_d208,S017_d180,S017_d188,S017_d208, 
                               S025_d180,S025_d188,S025_d208, S030_d180,S030_d188,S030_d208, S034_d180,S034_d188,S034_d208, 
                               S036_d180,S036_d188,S036_d208, S046_d180,S046_d188,S046_d208)
#select day 180 and 188
rawdata3 <- rawdata %>% select(1,2,3,4,6,7,9,10,12,13,15,16,18,19,21,22,24,25,27,28,30,31,33,34)

#select day 188 and 208
rawdata2 <- rawdata %>% select(1,2,4,5,7,8,10,11,13,14,16,17,19,20,22,23,25,26,28,29,31,32,34,35)

#select day 180 and 208
rawdata2 <- rawdata %>% select(1,2,3,5,6,8,9,11,12,14,15,17,18,20,21,23,24,26,27,29,30,32,33,35)

head(rawdata3)
y <- DGEList(counts=rawdata3[,3:24], genes=rawdata[,1:2])
o <- order(rowSums(y$counts), decreasing=TRUE)
y <- y[o,]
d <- duplicated(y$genes$Gene.name)
y <- y[!d,]
nrow(y)

y$samples$lib.size <- colSums(y$counts)
y$genes$Geneid <- NULL
y <- calcNormFactors(y)
y$samples


Subject <- factor(c(4,4,6,6,10,10,11,11,12,12,17,17,25,25,30,30,34,34,36,36,46,46))

#select day 180 and 188
Day_0v8 <-factor(c(0,8,0,8,0,8,0,8,0,8,0,8,0,8,0,8,0,8,0,8,0,8))
#select day 188 and 208
Day_8v28 <-factor(c(8,28,8,28,8,28,8,28,8,28,8,28,8,28,8,28,8,28,8,28,8,28))
#select day 180 and 208
Day_0v28 <-factor(c(0,28,0,28,0,28,0,28,0,28,0,28,0,28,0,28,0,28,0,28,0,28))
#select all subjects
Subject <- factor(c(4,4,4,6,6,6,10,10,10,11,11,11,12,12,12,17,17,17,25,25,25,30,30,30,34,34,34,36,36,36,46,46,46))

design <- model.matrix(~Subject+Day_0v8)
rownames(design) <- colnames(y)
design
y <- estimateDisp(y, design, robust=TRUE)
y$common.dispersion
plotBCV(y)

fit <- glmFit(y, design)
lrt <- glmLRT(fit)
tab<-topTags(lrt, n = 1000, sort.by = "logFC")
tab
write.table(tab, file="0v8.txt")
summary(decideTests(lrt,p.value=0.05))
plotMD(lrt,highlight=100,names=lrt$genes$Gene.name)
abline(h=c(-1, 1), col="blue")
go <- goana(lrt)
topGO <- topGO(go, ont="BP", sort="Up", n=30, truncate=30)
write.table(topGO, file="0v8GO.txt")


kegg <- kegga(lrt, geneid = rownames(lrt), FDR = 0.05, trend = FALSE)

topKEGG <- topGO(kegg, ont="BP", sort="Up", n=30, truncate=30)
topKEGG
