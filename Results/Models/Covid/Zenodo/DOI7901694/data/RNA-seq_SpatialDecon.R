library(SpatialDecon) 

biolog<-function(x) { #base-2 biolog
  asinh(x/2)/log(2)
}

####################define LUNG training set########################################################
Human<-download_profile_matrix("Human_Cell_Landscape")
dim(Human) #27341 genes x 102 cell types
#plot(Human[,c(48,89)]) #some Neutrophil types, "Neutrophil.1" "Neutrophil.2"
#abline(0,1)
#Human is in log2-scale already
#rownames(Human) has human gene symbols

colnames(Human)
* [2] "Macrophage"
* [3] "B.cell..Plasmocyte."
* [4] "Fibroblast"
 [6] "T.cell"
 [8] "Endothelial.cell..APC."
 [13] "Monocyte"
*[14] "B.cell..Plasmocyte..1"
*[18] "Fibroblast.1"
 [20] "Endothelial.cell"
 [22] "Dendritic.cell"
*[24] "T.cell.1"
 [27] "Stromal.cell"
*[29] "Endothelial.cell.1"
 [30] "AT2.cell"
 [31] "Monocyte.1"
*[33] "Smooth.muscle.cell"
 [34] "Fibroblast.2"
 [35] "Smooth.muscle.cell.1"
 [36] "Fibroblast.3"
 [37] "B.cell"
 [38] "Epithelial.cell"
 [40] "Neutrophil..RPS.high."
 [41] "Antigen.presenting.cell..RPS.high."
*[42] "Smooth.muscle.cell.2"
 [45] "Macrophage.1"
 [46] "Neutrophil"
*[47] "M2.Macrophage"
 [48] "Neutrophil.1"
 [49] "Epithelial.cell.1"
 [51] "Macrophage.2"
 [52] "Proliferating.T.cell"
 [54] "Goblet.cell"
 [57] "Stratified.epithelial.cell"
 [58] "Stromal.cell.1"
 [59] "Epithelial.cell.2"
 [60] "Epithelial.cell..Intermediated."
 [61] "Sinusoidal.endothelial.cell"
*[66] "Endothelial.cell..endothelial.to.mesenchymal.transition."
 [67] "Basal.cell"
 [69] "Macrophage.3"
 [70] "Fibroblast.4"
 [72] "Stromal.cell.2"
*[76] "Mesothelial.cell"
*[78] "Macrophage.4"
 [79] "Smooth.muscle.cell.3"
 [81] "Epithelial.cell.3"
*[82] "Stratified.epithelial.cell.1"
 [86] "Mast.cell"
 [89] "Neutrophil.2"
 [90] "Goblet.cell.1"
 [91] "Goblet.cell.2"
 [93] "Myeloid.cell"
*[94] "Dendritic.cell.1"
 [95] "Stromal.cell.3"
 [98] "Epithelial.cell.4"
*[100] "B.cell.1"


lungset<-Human[,c(2:4,6,8,13,14,18,20,22,24,27,29:31,,33:38,40:42,45:49,51,52,54,57:61,66,67,69,70,72,76,79,81,82,86,89:91,93:95,98,100)]
colnames(lungset)
 [1] "Macrophage"
 [2] "B.cell..Plasmocyte."
 [3] "Fibroblast"
 [4] "T.cell"
 [5] "Endothelial.cell..APC."
 [6] "Monocyte"
 [7] "B.cell..Plasmocyte..1"
 [8] "Fibroblast.1"
 [9] "Endothelial.cell"
[10] "Dendritic.cell"
[11] "T.cell.1"
[12] "Stromal.cell"
[13] "Endothelial.cell.1"
[14] "AT2.cell"
[15] "Monocyte.1"
[16] "Thyroid.follicular.cell"
[17] "Smooth.muscle.cell"
[18] "Fibroblast.2"
[19] "Smooth.muscle.cell.1"
[20] "Fibroblast.3"
[21] "B.cell"
[22] "Epithelial.cell"
[23] "Neutrophil..RPS.high."
[24] "Antigen.presenting.cell..RPS.high."
[25] "Smooth.muscle.cell.2"
[26] "Macrophage.1"
[27] "Neutrophil"
[28] "M2.Macrophage"
[29] "Neutrophil.1"
[30] "Epithelial.cell.1"
[31] "Macrophage.2"
[32] "Proliferating.T.cell"
[33] "Goblet.cell"
[34] "Stratified.epithelial.cell"
[35] "Stromal.cell.1"
[36] "Epithelial.cell.2"
[37] "Epithelial.cell..Intermediated."
[38] "Sinusoidal.endothelial.cell"
[39] "Endothelial.cell..endothelial.to.mesenchymal.transition."
[40] "Basal.cell"
[41] "Macrophage.3"
[42] "Fibroblast.4"
[43] "Stromal.cell.2"
[44] "Mesothelial.cell"
[45] "Smooth.muscle.cell.3"
[46] "Epithelial.cell.3"
[47] "Stratified.epithelial.cell.1"
[48] "Mast.cell"
[49] "Neutrophil.2"
[50] "Goblet.cell.1"
[51] "Goblet.cell.2"
[52] "Myeloid.cell"
[53] "Dendritic.cell.1"
[54] "Stromal.cell.3"
[55] "Epithelial.cell.4"
[56] "B.cell.1"
################end define LUNG training set######################################################################




###################################read data########################################################################
D<-read.csv("COVID gene expression data.csv",header=TRUE,stringsAsFactors =FALSE)
descr<-D[,2]
data<-data.matrix(D[,-(1:4)])
rownames(data)<-descr

#last 8 samples are COMPLETELY different! They need to be normalized separately
LUNG<-data[,1:17]
BALF<-data[,18:25] #Should be normalized and analyzed separately. Too different.
###############################end read data########################################################################





###############################LUNG analysis########################################################################
#cpm is a kind of rough normalization
data<-LUNG
cpm<-t(t(data)/colSums(data)*1e6)

ldata<-biolog(cpm)
a<-rowMeans(ldata)
hist(a,breaks=1000,ylim=c(0,200))
good<-a>1
gcpm<-cpm[good,]
gdescr<-D[good,1:4]
gldata<-ldata[good,]
ga<-a[good]

#RLE normalization
normfactor<-apply(gldata-ga,2,median)
ngcpm<-t(t(gcpm)/2^normfactor) #normalize all good data
ngldata<-biolog(ngcpm)
#end RLE normalization

ncpm<-t(t(cpm)/2^normfactor) #normalize all data
#ngcpm is the lin-scale data used for deconvolution

sum(rownames(ngcpm) %in% rownames(Human))
11111



#find background
#bg<-derive_GeoMx_background(norm = ncpm, probepool = rep(1, nrow(ncpm)), negnames = descr[!good])

#or set background
bg<-matrix(5,dim(ngcpm)[1],dim(ngcpm)[2]) #bg set to 5 because at lower bg there is an artifact in residuals distr.
rownames(bg)<-rownames(ngcpm)
colnames(bg)<-colnames(ngcpm)
LUNG_Human_res <- spatialdecon(norm = ngcpm, bg = bg, X = lungset, align_genes = TRUE)
#check residuals - looks normal? Not really, but that's all I can do.
hist(LUNG_Human_res$resids,breaks=200)

#plot cellular composition of every sample
pdf("LUNG.cellular.fractions.barplot.pdf",width=14)
par(mar=c(22,3,1,1))
for (i in 1:dim(LUNG_Human_res$prop_of_all)[2]) {
  barplot(LUNG_Human_res$prop_of_all[,i],ylim=c(0,.75),las=2,main=colnames(ngcpm)[i])
  box()
}
dev.off()

#separate COVID from nonCOVID samples, do boxplots
COV<-grep("COV",colnames(ngcpm))
nonCOV<-setdiff(1:dim(ngcpm)[2],COV)
pdf("LUNG.cellular.fractions.boxplots.pdf",width=14,height=10)
par(mar=c(22,3,1,1))
boxplot(t(LUNG_Human_res$prop_of_all[,nonCOV]),las=2,ylim=c(0,.8),main="lung, non-COV")
boxplot(t(LUNG_Human_res$prop_of_all[,COV]),las=2,ylim=c(0,.8),main="lung, COV")
dev.off()
###########################end LUNG analysis########################################################################



############################## BALF analysis########################################################################
#cpm is a kind of rough normalization
data<-BALF
cpm<-t(t(data)/colSums(data)*1e6)

ldata<-biolog(cpm)
a<-rowMeans(ldata)
hist(a,breaks=1000,ylim=c(0,200))
good<-a>3 #roughly log2(10)
gcpm<-cpm[good,]
gdescr<-D[good,1:4]
gldata<-ldata[good,]
ga<-a[good]

#PCA plot:
pdf("BALF_PCA.pdf")
pca.res <- prcomp(t(gldata))
xy<-pca.res$x[,1:2]
rge<-range(xy)
plot(xy,pch=16,cex=2,xlim=rge,ylim=rge,main="principal component analysis",xlab="PC 1",ylab="PC 2")
for (i in 1:dim(gldata)[2]) {
points(xy[i,1],xy[i,2],pch=16,cex=3)
text(xy[i,1],xy[i,2],as.character(i),cex=1,col="white")
}
dev.off()
#remove sample 6 

cpm<-cpm[,-6]
ldata<-ldata[,-6]
a<-rowMeans(ldata)
hist(a,breaks=1000,ylim=c(0,200))
good<-a>3 #roughly log2(10)
gcpm<-cpm[good,]
gdescr<-D[good,1:4]
gldata<-ldata[good,]
ga<-a[good]

#RLE normalization
normfactor<-apply(gldata-ga,2,median)
ngcpm<-t(t(gcpm)/2^normfactor) #normalize all good data
ngldata<-biolog(ngcpm)
#end RLE normalization

ncpm<-t(t(cpm)/2^normfactor) #normalize all data
#ngcpm is the lin-scale data used for deconvolution

sum(rownames(ngcpm) %in% rownames(Human))
4514

#reorder genes in lungset for convenience
olungset<-lungset[,c(1,25,30,40,44,27, 26,28,49,22, 4,11,31, 20,56,7, 2,3,5,6,8,9,10,12,13,14,15,16,17,18,19,21,23,24,29,32,33,34,35,36,37,38,39,41,42,43,45,46,47,48,50,51,52,53,54,55)]

bg<-matrix(8,dim(ngcpm)[1],dim(ngcpm)[2])
rownames(bg)<-rownames(ngcpm)
colnames(bg)<-colnames(ngcpm)
BALF_Human_res <- spatialdecon(norm = ngcpm, bg = bg, X = olungset, align_genes = TRUE)
#check residuals
hist(BALF_Human_res$resids,breaks=200) #not bad

pdf("BALF.cellular.fractions.barplot.pdf",width=14)
par(mar=c(22,3,1,1))
for (i in 1:dim(BALF_Human_res$prop_of_all)[2]) {
  barplot(BALF_Human_res$prop_of_all[,i],ylim=c(0,.75),las=2,main=colnames(ngcpm)[i])
  box()
}
dev.off()

#plot boxplots across all samples
pdf("BALF.boxplots.pdf",width=14)
par(mar=c(22,3,1,1))
boxplot(t(BALF_Human_res$prop_of_all),las=2,ylim=c(0,.7),ylab="cellular fraction",main="BALF")
dev.off()
###########################end BALF analysis########################################################################




##################begin heatmaps of Human matching types####################################
#clumping together related cell types - this defines matching categories
cell.match<-list()
cell.match[["neutrophils"]]<-c("Neutrophil","Neutrophil.1","Neutrophil.2","Neutrophil..RPS.high.")
cell.match[["T cells"]]<-c("T.cell","T.cell.1","Proliferating.T.cell")
cell.match[["dendritic cells"]]<-c("Dendritic.cell","Dendritic.cell.1")
cell.match[["B cells"]]<-c("B.cell","B.cell.1","B.cell..Plasmocyte.","B.cell..Plasmocyte..1")
cell.match[["monocytes"]]<-c("Monocyte","Monocyte.1")
cell.match[["macrophages"]]<-c("Macrophage","Macrophage.1","Macrophage.2","Macrophage.3","Macrophage.4","M2.Macrophage")
cell.match[["antigen-presenting cells, high RPs"]]<-c("Antigen.presenting.cell..RPS.high.")
cell.match[["epithelial cells"]]<-c("Epithelial.cell","Epithelial.cell.1","Epithelial.cell.2","Epithelial.cell.3","Epithelial.cell.4")
cell.match[["smooth muscle cells"]]<-c("Smooth.muscle.cell","Smooth.muscle.cell.1","Smooth.muscle.cell.2","Smooth.muscle.cell.3")


source("heatmap.2_original.R")
for (ctype in c("neutrophils","T cells","dendritic cells","B cells","monocytes","macrophages","epithelial cells","smooth muscle cells")) {
#this section is to show what makes 4 types of neutrophils distinct
ctypedata<-newHuman[,cell.match[[ctype]]]

v<-apply(ctypedata,1,var)
plot(v)
ix<-match(rownames(gdata),rownames(ctypedata))
points(ix,v[ix],col="red")

o<-order(-v)
ctypedata.o<-ctypedata[o[1:120],]

pdf(paste0(ctype,".pdf"),width=12,height=20)
heatmap.2.1(sweep(ctypedata.o, 1, apply(ctypedata.o, 1, max), "/"),key=FALSE,Colv=FALSE,margins = c(25, 10), cexCol = 2,lhei=c(0.05,5),dendrogram="none",trace="none")
dev.off()
}
##################end heatmaps of Human matching types#################################

