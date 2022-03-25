library(Seurat)
library(magrittr)
library(harmony)
library(dplyr)
library(tidyverse)
library(patchwork)
library(ggsci)
library(pheatmap)
setwd("/media/feng/data1/骨肉瘤挖掘/肺转移")
OS10.data <- Read10X(data.dir = "/media/feng/data1/骨肉瘤挖掘/BC10")
OS17.data <- Read10X(data.dir = "/media/feng/data1/骨肉瘤挖掘/BC17")

OS10 <- CreateSeuratObject(counts = OS10.data, project = "patient 10", min.cells = 3, min.features = 200)
OS17 <- CreateSeuratObject(counts = OS17.data, project = "patient 17", min.cells = 3, min.features = 200)

OS.combined <- merge(OS10, y = list(OS17), add.cell.ids = c("P10","P17"), project = "Osteosarcoma")
table(OS.combined$orig.ident)

OS.combined[["percent.mt"]] <- PercentageFeatureSet(OS.combined, pattern = "^MT-")
OS.combined[["percent.rb"]] <- PercentageFeatureSet(OS.combined, pattern = "^RB[SL]")
p <- VlnPlot(OS.combined, features = c("nFeature_RNA", "nCount_RNA", "percent.mt", "percent.rb"), ncol = 4,pt.size=0)
ggsave("vlnplot_before_qc.png",p)
ggsave("vlnplot_before_qc.eps",p)
ggsave("vlnplot_before_qc.pdf",p)


plot1 <- FeatureScatter(OS.combined, feature1 = "nCount_RNA", feature2 = "percent.mt")
plot2 <- FeatureScatter(OS.combined, feature1 = "nCount_RNA", feature2 = "nFeature_RNA")
p <- CombinePlots(plots = list(plot1, plot2))
ggsave("nCount_RNAVSpercent.mt.png",p,width = 13, height = 5.5)
ggsave("nCount_RNAVSpercent.mt.eps",p,width = 13, height = 5.5)
ggsave("nCount_RNAVSpercent.mt.pdf",p,width = 13, height = 5.5)


OS.combined<- subset(OS.combined, subset = nFeature_RNA > 300 & nFeature_RNA < 5000 & percent.mt < 10)#之前是200 6000 25 核糖体无限制
OS.combined <- NormalizeData(OS.combined, normalization.method = "LogNormalize", scale.factor = 10000)
OS.combined <- FindVariableFeatures(OS.combined, selection.method = "vst", nfeatures = 2000)
top10 <- head(VariableFeatures(OS.combined), 10)
OS.combined <- FindVariableFeatures(OS.combined, selection.method = "vst", nfeatures = 2000)
top10 <- head(VariableFeatures(OS.combined), 10)
pdf(file="High variable gene.pdf",width = 13, height = 5.5)
plot1 <- VariableFeaturePlot(OS.combined)
plot2 <- LabelPoints(plot = plot1, points = top10, repel = TRUE)
CombinePlots(plots = list(plot1, plot2))
dev.off()

all.genes <- rownames(OS.combined)
OS.combined <- ScaleData(OS.combined, features = all.genes)

g2m_genes <- cc.genes$g2m.genes
g2m_genes <- CaseMatch(search=g2m_genes, match=rownames(OS.combined))
s_genes <- cc.genes$s.genes    
s_genes <- CaseMatch(search=s_genes, match=rownames(OS.combined))
OS.combined <- CellCycleScoring(OS.combined, g2m.features=g2m_genes, s.features=s_genes)
tmp <- RunPCA(OS.combined, features = c(g2m_genes, s_genes), verbose = F)
p <- DimPlot(tmp, reduction = "pca", group.by = "orig.ident")
ggsave("CellCycle_pca.png", p, width = 8, height = 6)
ggsave("CellCycle_pca.pdf", p, width = 8, height = 6)

mt.genes <- grep("^MT-", rownames(OS.combined), value=T, ignore.case=T)
tmp <- RunPCA(OS.combined, features = mt.genes, verbose = F)
p <- DimPlot(tmp, reduction = "pca", group.by = "orig.ident")
ggsave("mito_pca.png", p, width = 8, height = 6)
ggsave("mito_pca.pdf", p, width = 8, height = 6)
rm(tmp)

rb.genes <- grep("^RP[SL]", rownames(OS.combined), value=T, ignore.case=T)
tmp <- RunPCA(OS.combined, features = rb.genes, verbose = F)
p <- DimPlot(tmp, reduction = "pca", group.by = "orig.ident")
ggsave("ribo_pca.png", p, width = 8, height = 6)
ggsave("mito_pca.pdf", p, width = 8, height = 6)
rm(tmp)

OS.combined <- RunPCA(OS.combined, features = VariableFeatures(object = OS.combined))
ElbowPlot(OS.combined, ndims = 50)

OS.combined <- OS.combined %>% 
  RunHarmony("orig.ident", plot_convergence = TRUE)
harmony_embeddings <- Embeddings(OS.combined, 'harmony')
harmony_embeddings[1:5, 1:5]

OS.combined <- OS.combined %>%
  RunUMAP(reduction="harmony", dims=1:20) %>%
  FindNeighbors(reduction = "harmony", dims = 1:20) %>% 
  FindClusters(resolution = 0.02) %>% 
  identity()

DimPlot(OS.combined, group.by = "orig.ident")
DimPlot(OS.combined)

lung <- readRDS("/media/feng/data1/骨肉瘤挖掘/肺转移/lung_ts.rds")
OS.combined$Donor <- OS.combined$orig.ident
OS.combined$Donor <- ifelse(OS.combined$Donor == "patient 10","BC10","BC17")
OS_lung <- merge(OS.combined, lung, add.cell.ids = c("OSlungmetastasis","Lungtissue"))
table(OS_lung$Donor)


OS_lung <- NormalizeData(OS_lung) %>% FindVariableFeatures() %>%ScaleData(features = rownames(OS_lung))
OS_lung <- RunPCA(OS_lung, npcs=50, verbose=FALSE)
DimPlot(OS_lung, reduction = "pca", group.by = "Donor")

OS_lung <- RunHarmony(OS_lung, group.by.vars="Donor",max.iter.harmony = 15)
ElbowPlot(OS_lung, ndims = 30)

pc.num=1:20  

OS_lung <- OS_lung %>%
  RunUMAP(reduction="harmony", dims=1:20) %>%
  FindNeighbors(reduction = "harmony", dims = 1:20) %>% 
  FindClusters(resolution = 0.5) %>% 
  identity()

my36colors <-c('#E5D2DD', '#53A85F', '#F1BB72', '#F3B1A0', '#D6E7A3', '#57C3F3', '#476D87',
               '#E95C59', '#E59CC4', '#AB3282', '#23452F', '#BD956A', '#8C549C', '#585658',
               '#9FA3A8', '#E0D4CA', '#5F3D69', '#C5DEBA', '#58A4C3', '#E4C755', '#F7F398',
               '#AA9A59', '#E63863', '#E39A35', '#C1E6F3', '#6778AE', '#91D0BE', '#B53E2B',
               '#712820', '#DCC1DD', '#CCE0F5',  '#CCC9E6', '#625D9E', '#68A180', '#3A6963',
               '#968175')

OS_lung@active.ident <- OS_lung$seurat_clusters
new.cluster.ids <- c("T cells", "T cells", "Malignant cells",	"Monocytes",	"Alveolar cells",	"T cells",	"Malignant cells",	"Osteoclasts",	"Endothelial",	"Mast cells",	"Proliferative cells", 	"Macrophages",	"Fibroblasts",	"Fibroblasts",	"Dendritic cells",	"T cells",	"B cells",	"Muscle cells",	"Secretory cells",	"Fibroblasts",	"Ciliated cells",	"Low quality cells",	"Plasma cells",	"T cells",	"Proliferative cells", "T cells"
)
names(new.cluster.ids) <- levels(OS_lung)

OS_lung <- RenameIdents(OS_lung, new.cluster.ids)
p <- DimPlot(OS_lung,cols = my36colors) 
ggsave("allumap.png",p)
ggsave("allumap.eps",p)
ggsave("allumap.pdf",p)
p <- DimPlot(OS_lung,group.by = "Donor") + scale_color_jama()
ggsave("donorumap.png",p)
ggsave("donorumap.eps",p)
ggsave("donorumap.pdf",p)

meta <- as.data.frame(OS_lung@meta.data)
meta$bingli <- meta$Donor
meta$bingli <- gsub("BC10","Conventional",meta$bingli)
meta$bingli <- gsub("BC17","Chondroblastic",meta$bingli)
meta$bingli <- gsub("343B","Normal",meta$bingli)
meta$bingli <- gsub("356C","Normal",meta$bingli)
meta$bingli <- gsub("367C","Normal",meta$bingli)
meta$bingli <- gsub("368C","Normal",meta$bingli)
meta$bingli <- gsub("390C","Normal",meta$bingli)
OS_lung@meta.data$bingli <- meta$bingli

p <- DimPlot(OS_lung,group.by = "bingli") + scale_color_aaas()
ggsave("bingliumap.png",p)
ggsave("bingliumap.eps",p)
ggsave("bingliumap.pdf",p)

Average_OSlung<- AverageExpression(OS_lung, assays = NULL, features = NULL,return.seurat = FALSE, add.ident = NULL, slot = "data",use.scale = FALSE, use.counts = FALSE, verbose = TRUE)
gene<-Average_OSlung$RNA
gene<-as.matrix(gene) 
usegene <- c('CD3E','CD3D', 'CD3G', #NK/T cells
             "SOX9","ALPL","RUNX2",#chondrocyte)
             'CD14', 'VCAN', 'FCN1', #Monocyte
             'SFTPC', 'SLPI', #Type 1 (AGER) and 2 alveolar cells
             "CD163","MRC1","MSR1", #monocyte/macrophage
             "ACP5","ATP6V0D2","CTSK",#Osteoclasts,
             'PECAM1',"EGFL7","PLVAP", #Endothelial
             "KIT" ,"HPGDS","TPSB2", #Mast cell
             'MKI67','CDK1',"TOP2A", #Proliferative cells
             'APOE', 'CD5L', #Macrophage
             'PDGFRA', 'FAP' ,#Fibroblast
             'CLEC9A', 'XCR1', #cDC1
             'CD1C',  'CLEC10A', #cDC2
             'CD19', 'CD79A', 'MS4A1', #B-cell
             'ACTA2', 'TAGLN', #Muscle cells
             'MUC5AC', 'MUC5B', 'SCGB1A1' ,#Secretory
             'EPCAM','PIFO', 'FOXJ1' ,#ciliated
             'XBP1', 'JCHAIN') #Plasma cells

gene <- as.data.frame(gene)
gene <- gene[usegene,]
t <- t(scale(t(gene)))

p = pheatmap(t, cluster_cols = F,cluster_rows = F, color = colorRampPalette(c("navy", "white", "firebrick3"))(50),fontsize = 5,border_color = "white")
ggsave("allmarkerheatmap.png",p,width = 6,height = 6)
ggsave("allmarkerheatmap.eps",p,width = 6,height = 6)
ggsave("allmarkerheatmap.pdf",p,width = 6,height = 6)

saveRDS(OS_lung,"renamesOS_lung.rds")


setwd("/media/feng/data1/骨肉瘤挖掘/肺转移/Tcell")

NKTcell <- subset(OS_lung,idents = c("T cells"))
NKTcell <- NormalizeData(NKTcell) %>% FindVariableFeatures() %>%ScaleData(features = rownames(NKTcell))
NKTcell <- RunPCA(NKTcell, npcs=50, verbose=FALSE)
DimPlot(NKTcell, reduction = "pca", group.by = "Donor")

NKTcell <- RunHarmony(NKTcell, group.by.vars="Donor",max.iter.harmony = 15)
ElbowPlot(NKTcell, ndims = 30)

pc.num=1:10  

NKTcell <- NKTcell %>%
  RunUMAP(reduction = "harmony",dims=pc.num) %>%
  FindNeighbors(reduction = "harmony",dims = pc.num) %>%
  FindClusters(resolution = 0.1) %>%
  identity()

p <- DimPlot(NKTcell, group.by = "Donor") + scale_color_ucscgb()
ggsave("NKTumap_Donor.png",p)
ggsave("NKTumap_Donor.eps",p)
ggsave("NKTumap_Donor.pdf",p)
p <- DimPlot(NKTcell) + scale_color_jama()
ggsave("NKTumap_source.png",p)
ggsave("NKTumap_source.eps",p)
ggsave("NKTumap_source.pdf",p)
p <- DimPlot(NKTcell,group.by = "bingli") + scale_color_aaas()
ggsave("bingliNKTumap.png",p)
ggsave("bingliNKTumap.eps",p)
ggsave("bingliNKTumap.pdf",p)

NKTcell@active.ident <- NKTcell$seurat_clusters
new.cluster.ids <- c("Lung tissue","Lung tissue","Lung tissue","Lung tissue","Lung tissue(Conventional)","Lung tissue(Chondroblastic)")
names(new.cluster.ids) <- levels(NKTcell)
NKTcell <- RenameIdents(NKTcell, new.cluster.ids)
DimPlot(NKTcell, group.by = "Donor")
p <- DimPlot(NKTcell) + scale_color_aaas()
ggsave("Tcellumap.png",p)
ggsave("Tcellumap.eps",p)
ggsave("Tcellumap.pdf",p)

options(stringsAsFactors = F)
library(DESeq2)
library(limma)
library(stringr)
library(ggsci)
library(ggpubr)

NKT_CON <- subset(NKTcell,idents = c("Lung tissue","Lung tissue(Conventional)"))
tmp_singlecell <- as.data.frame(NKT_CON@assays$RNA@counts)
tmp_singlecell <- tmp_singlecell[rowMeans(tmp_singlecell) >0,]
tmp_singlecell <- tmp_singlecell +1
k <- tmp_singlecell[1:10,1:10]
group_singlecell <- as.data.frame(NKT_CON@active.ident)
group_singlecell$ID <- rownames(group_singlecell)

group_list=as.factor(group_singlecell$`NKT_CON@active.ident`)#做好分组因子

colData <- data.frame(row.names=group_singlecell$ID, group_list=group_list)
dds <- DESeqDataSetFromMatrix(countData = ceiling(tmp_singlecell[,rownames(colData)]),
                              colData = colData,
                              design = ~ group_list)
dds2 <- DESeq(dds)  ##第二步，直接用DESeq函数即可
resultsNames(dds2)
res <-  results(dds2, contrast=c("group_list","Lung tissue(Conventional)","Lung tissue"))#1是转移，0是不转移
## 提取你想要的差异分析结果，我们这里是Endothelial cells_lung metastasis组对Endothelial cells_primary组进行比较
resOrdered <- res[order(res$padj),]
resOrdered=as.data.frame(resOrdered)
resOrdered = na.omit(resOrdered)
write.csv(resOrdered,file = "NKT_CONsinglecellecDEseq2DEG.csv")

head(resOrdered)
resOrdered$logP <- -log10(resOrdered$padj)
resOrdered$Group <-"Not-significant"
resOrdered$Group[which((resOrdered$padj <(0.05))&(resOrdered$log2FoldChange >(1)))] <-"Up-regulated"
resOrdered$Group[which((resOrdered$padj <(0.05))&(resOrdered$log2FoldChange <(-1)))] <-"Down-regulated"
table(resOrdered$Group)
resOrdered$Label = ""
resOrdered <- resOrdered[order(resOrdered$padj),]
resOrdered$symbol <- rownames(resOrdered)

up_gene <- head(resOrdered$symbol[which(resOrdered$Group == "Up-regulated")],10)
down_gene <- head(resOrdered$symbol[which(resOrdered$Group == "Down-regulated")],10)
resOrdered_top10gene <- c(as.character(up_gene),as.character(down_gene))
resOrdered$Label[match(resOrdered_top10gene,resOrdered$symbol)] <- resOrdered_top10gene

p <-ggscatter(resOrdered,x = "log2FoldChange",y = "logP",
              color = "Group",
              palette = c("#2f5688","#BBBBBB","#CC0000"),
              size = 1,
              label = resOrdered$Label,
              font.label = 8,
              repel = T,
              xlab = "log2FoldChange",
              ylab = "-log10(P-value)")  + geom_hline(yintercept = 1.3, linetype = "dashed")+
  geom_vline(xintercept = c(-1,1), linetype = "dashed")
ggsave("NKT经典型对比火山图.png",p)
ggsave("NKT经典型对比火山图.pdf",p)
ggsave("NKT经典型对比火山图.eps",p)

ungene1 <- resOrdered[resOrdered$Group == "Up-regulated",]
library("clusterProfiler")
library("org.Hs.eg.db")
library("enrichplot")
library("ggplot2")
genes=as.vector(ungene1[,10])
entrezIDs <- mget(genes, org.Hs.egSYMBOL2EG, ifnotfound=NA)    
entrezIDs <- as.character(entrezIDs)
out=cbind(ungene1$symbol,entrezID=entrezIDs)
dir.create("/media/feng/data1/骨肉瘤挖掘/肺转移/Tcell/NKT_CON_UP/")
setwd("/media/feng/data1/骨肉瘤挖掘/肺转移/Tcell/NKT_CON_UP/")
write.table(out,file="id.txt",sep="\t",quote=F,row.names=F)   
rt=read.table("id.txt",sep="\t",header=T,check.names=F)           
rt=rt[is.na(rt[,"entrezID"])==F,]                                 
gene=rt$entrezID
kk <- enrichGO(gene = gene, OrgDb = org.Hs.eg.db, pvalueCutoff =0.05, qvalueCutoff = 0.05, ont="all", readable =T)
write.table(kk,file="GO.txt",sep="\t",quote=F,row.names = F)                 

p <-barplot(kk, drop = TRUE, showCategory =10,split="ONTOLOGY") + facet_grid(ONTOLOGY~., scale='free')
ggsave("GObarplot.png",p,width = 10,height = 8)
ggsave("GObarplot.eps",p,width = 10,height = 8)
ggsave("GObarplot.pdf",p,width = 10,height = 8)

p <- dotplot(kk,showCategory = 10,split="ONTOLOGY") + facet_grid(ONTOLOGY~., scale='free')
ggsave("GObubble.png",p,width = 10,height = 8)
ggsave("GObubble.eps",p,width = 10,height = 8)
ggsave("GObubble.pdf",p,width = 10,height = 8)

resOrdered <- read.csv("/media/feng/data1/骨肉瘤挖掘/肺转移/Tcell/NKT_CONsinglecellecDEseq2DEG.csv")
head(resOrdered)
resOrdered$logP <- -log10(resOrdered$padj)
resOrdered$Group <-"Not-significant"
resOrdered$Group[which((resOrdered$padj <(0.05))&(resOrdered$log2FoldChange >(1)))] <-"Up-regulated"
resOrdered$Group[which((resOrdered$padj <(0.05))&(resOrdered$log2FoldChange <(-1)))] <-"Down-regulated"
table(resOrdered$Group)
resOrdered$Label = ""
resOrdered <- resOrdered[order(resOrdered$padj),]
resOrdered$symbol <- rownames(resOrdered)

up_gene <- head(resOrdered$symbol[which(resOrdered$Group == "Up-regulated")],10)
down_gene <- head(resOrdered$symbol[which(resOrdered$Group == "Down-regulated")],10)
resOrdered_top10gene <- c(as.character(up_gene),as.character(down_gene))
resOrdered$Label[match(resOrdered_top10gene,resOrdered$symbol)] <- resOrdered_top10gene
ungene1 <- resOrdered[resOrdered$Group == "Up-regulated",]
library("clusterProfiler")
library("org.Hs.eg.db")
library("enrichplot")
library("ggplot2")
setwd("/media/feng/data1/骨肉瘤挖掘/肺转移/Tcell")
dir.create("/media/feng/data1/骨肉瘤挖掘/肺转移/Tcell/NKTCONKEGG_UP")
setwd("/media/feng/data1/骨肉瘤挖掘/肺转移/Tcell/NKTCONKEGG_UP")
genes=as.vector(ungene1$X)
entrezIDs <- mget(genes, org.Hs.egSYMBOL2EG, ifnotfound=NA)    
entrezIDs <- as.character(entrezIDs)
out=cbind(genes,entrezID=entrezIDs)
write.table(out,file="id.txt",sep="\t",quote=F,row.names=F)   
rt=read.table("id.txt",sep="\t",header=T,check.names=F)           
rt=rt[is.na(rt[,"entrezID"])==F,]                                 
gene=rt$entrezID
kk <- enrichKEGG(gene = gene, organism = "hsa", pvalueCutoff =0.05, qvalueCutoff =0.05)
write.table(kk,file="KEGG.txt",sep="\t",quote=F,row.names = F)                 

p <- barplot(kk, drop = TRUE, showCategory = 30)
ggsave("KEGGbarplot.png",p,width = 10,height = 8)
ggsave("KEGGbarplot.eps",p,width = 10,height = 8)
ggsave("KEGGbarplot.pdf",p,width = 10,height = 8)

p <- dotplot(kk, showCategory = 30)
ggsave("KEGGdotplot.png",p,width = 10,height = 8)
ggsave("KEGGdotplot.eps",p,width = 10,height = 8)
ggsave("KEGGdotplot.pdf",p,width = 10,height = 8)

setwd("/media/feng/data1/骨肉瘤挖掘/肺转移/Tcell/")
ungene2 <- resOrdered[resOrdered$Group == "Down-regulated",]
genes=as.vector(ungene2[,10])
entrezIDs <- mget(genes, org.Hs.egSYMBOL2EG, ifnotfound=NA)    
entrezIDs <- as.character(entrezIDs)
out=cbind(ungene1$symbol,entrezID=entrezIDs)
dir.create("/media/feng/data1/骨肉瘤挖掘/肺转移/Tcell/NKT_CON_DOWN/")
setwd("/media/feng/data1/骨肉瘤挖掘/肺转移/Tcell/NKT_CON_DOWN/")
write.table(out,file="id.txt",sep="\t",quote=F,row.names=F)   
rt=read.table("id.txt",sep="\t",header=T,check.names=F)           
rt=rt[is.na(rt[,"entrezID"])==F,]                                 
gene=rt$entrezID
kk <- enrichGO(gene = gene, OrgDb = org.Hs.eg.db, pvalueCutoff =0.05, qvalueCutoff = 0.05, ont="all", readable =T)
write.table(kk,file="GO.txt",sep="\t",quote=F,row.names = F)                 

p <-barplot(kk, drop = TRUE, showCategory =10,split="ONTOLOGY") + facet_grid(ONTOLOGY~., scale='free')
ggsave("GObarplot.png",p,width = 7,height = 8)
ggsave("GObarplot.eps",p,width = 7,height = 8)
ggsave("GObarplot.pdf",p,width = 7,height = 8)

p <- dotplot(kk,showCategory = 10,split="ONTOLOGY") + facet_grid(ONTOLOGY~., scale='free')
ggsave("GObubble.png",p,width = 7,height = 8)
ggsave("GObubble.eps",p,width = 7,height = 8)
ggsave("GObubble.pdf",p,width = 7,height = 8)

ungene2 <- resOrdered[resOrdered$Group == "Down-regulated",]
library("clusterProfiler")
library("org.Hs.eg.db")
library("enrichplot")
library("ggplot2")
setwd("/media/feng/data1/骨肉瘤挖掘/肺转移/Tcell")
dir.create("/media/feng/data1/骨肉瘤挖掘/肺转移/Tcell/NKTCONKEGG_DOWN")
setwd("/media/feng/data1/骨肉瘤挖掘/肺转移/Tcell/NKTCONKEGG_DOWN")
genes=as.vector(ungene2$X)
entrezIDs <- mget(genes, org.Hs.egSYMBOL2EG, ifnotfound=NA)    
entrezIDs <- as.character(entrezIDs)
out=cbind(genes,entrezID=entrezIDs)
write.table(out,file="id.txt",sep="\t",quote=F,row.names=F)   
rt=read.table("id.txt",sep="\t",header=T,check.names=F)           
rt=rt[is.na(rt[,"entrezID"])==F,]                                 
gene=rt$entrezID
kk <- enrichKEGG(gene = gene, organism = "hsa", pvalueCutoff =0.05, qvalueCutoff =0.05)
write.table(kk,file="KEGG.txt",sep="\t",quote=F,row.names = F)                 

p <- barplot(kk, drop = TRUE, showCategory = 30)
ggsave("KEGGbarplot.png",p,width = 10,height = 8)
ggsave("KEGGbarplot.eps",p,width = 10,height = 8)
ggsave("KEGGbarplot.pdf",p,width = 10,height = 8)


p <- dotplot(kk, showCategory = 30)
ggsave("KEGGdotplot.png",p,width = 10,height = 8)
ggsave("KEGGdotplot.eps",p,width = 10,height = 8)
ggsave("KEGGdotplot.pdf",p,width = 10,height = 8)


NKT_CHO <- subset(NKTcell,idents = c("Lung tissue"))

#随机抽取10%的细胞进行差异分析
NKT_CHO1 <- subset(NKTcell,idents = c("Lung tissue(Chondroblastic)"))
NKT_CHO <- merge(NKT_CHO,NKT_CHO1)
tmp_singlecell <- as.data.frame(NKT_CHO@assays$RNA@counts)
tmp_singlecell <- tmp_singlecell[rowMeans(tmp_singlecell) >0,]
tmp_singlecell <- tmp_singlecell +1
k <- tmp_singlecell[1:10,1:10]
group_singlecell <- as.data.frame(NKT_CHO@active.ident)
group_singlecell$ID <- rownames(group_singlecell)

group_list=as.factor(group_singlecell$`NKT_CHO@active.ident`)

colData <- data.frame(row.names=group_singlecell$ID, group_list=group_list)
dds <- DESeqDataSetFromMatrix(countData = ceiling(tmp_singlecell[,rownames(colData)]),
                              colData = colData,
                              design = ~ group_list)
dds2 <- DESeq(dds)  
resultsNames(dds2)
res <-  results(dds2, contrast=c("group_list","Lung tissue(Chondroblastic)","Lung tissue"))

resOrdered <- res[order(res$padj),]
resOrdered=as.data.frame(resOrdered)
resOrdered = na.omit(resOrdered)
write.csv(resOrdered,file = "NKT_CHOsinglecellecDEseq2DEG.csv")

head(resOrdered)
resOrdered$logP <- -log10(resOrdered$padj)

resOrdered$Group <-"Not-significant"

resOrdered$Group[which((resOrdered$padj <(0.05))&(resOrdered$log2FoldChange >(1)))] <-"Up-regulated"
resOrdered$Group[which((resOrdered$padj <(0.05))&(resOrdered$log2FoldChange <(-1)))] <-"Down-regulated"
table(resOrdered$Group)

resOrdered$Label = ""
resOrdered <- resOrdered[order(resOrdered$padj),]
resOrdered$symbol <- rownames(resOrdered)

up_gene <- head(resOrdered$symbol[which(resOrdered$Group == "Up-regulated")],10)
down_gene <- head(resOrdered$symbol[which(resOrdered$Group == "Down-regulated")],10)
resOrdered_top10gene <- c(as.character(up_gene),as.character(down_gene))
resOrdered$Label[match(resOrdered_top10gene,resOrdered$symbol)] <- resOrdered_top10gene

p <-ggscatter(resOrdered,x = "log2FoldChange",y = "logP",
              color = "Group",
              palette = c("#2f5688","#BBBBBB","#CC0000"),
              size = 1,
              label = resOrdered$Label,
              font.label = 8,
              repel = T,
              xlab = "log2FoldChange",
              ylab = "-log10(P-value)")  + geom_hline(yintercept = 1.3, linetype = "dashed")+
  geom_vline(xintercept = c(-1,1), linetype = "dashed")
ggsave("NKT成软骨型对比火山图1.png",p)
ggsave("NKT成软骨型对比火山图1.pdf",p)
ggsave("NKT成软骨型对比火山图1.eps",p)

ungene3 <- resOrdered[resOrdered$Group == "Up-regulated",]
library("clusterProfiler")
library("org.Hs.eg.db")
library("enrichplot")
library("ggplot2")
genes=as.vector(ungene3[,10])
entrezIDs <- mget(genes, org.Hs.egSYMBOL2EG, ifnotfound=NA)    
entrezIDs <- as.character(entrezIDs)
out=cbind(ungene3$symbol,entrezID=entrezIDs)
dir.create("/media/feng/data1/骨肉瘤挖掘/肺转移/Tcell/NKT_CHO_UP/")
setwd("/media/feng/data1/骨肉瘤挖掘/肺转移/Tcell/NKT_CHO_UP/")
write.table(out,file="id.txt",sep="\t",quote=F,row.names=F)   
rt=read.table("id.txt",sep="\t",header=T,check.names=F)           
rt=rt[is.na(rt[,"entrezID"])==F,]                                 
gene=rt$entrezID
kk <- enrichGO(gene = gene, OrgDb = org.Hs.eg.db, pvalueCutoff =0.05, qvalueCutoff = 0.05, ont="all", readable =T)
write.table(kk,file="GO.txt",sep="\t",quote=F,row.names = F)                 

p <-barplot(kk, drop = TRUE, showCategory =10,split="ONTOLOGY") + facet_grid(ONTOLOGY~., scale='free')
ggsave("GObarplot.png",p,width = 11,height = 8)
ggsave("GObarplot.eps",p,width = 11,height = 8)
ggsave("GObarplot.pdf",p,width = 11,height = 8)

p <- dotplot(kk,showCategory = 10,split="ONTOLOGY") + facet_grid(ONTOLOGY~., scale='free')
ggsave("GObubble.png",p,width = 11,height = 8)
ggsave("GObubble.eps",p,width = 11,height = 8)
ggsave("GObubble.pdf",p,width = 11,height = 8)

resOrdered <- read.csv("/media/feng/data1/骨肉瘤挖掘/肺转移/Tcell/NKT_CHOsinglecellecDEseq2DEG.csv")
head(resOrdered)
resOrdered$logP <- -log10(resOrdered$padj)

resOrdered$Group <-"Not-significant"

resOrdered$Group[which((resOrdered$padj <(0.05))&(resOrdered$log2FoldChange >(1)))] <-"Up-regulated"
resOrdered$Group[which((resOrdered$padj <(0.05))&(resOrdered$log2FoldChange <(-1)))] <-"Down-regulated"
table(resOrdered$Group)

resOrdered$Label = ""
resOrdered <- resOrdered[order(resOrdered$padj),]
resOrdered$symbol <- rownames(resOrdered)

up_gene <- head(resOrdered$symbol[which(resOrdered$Group == "Up-regulated")],10)
down_gene <- head(resOrdered$symbol[which(resOrdered$Group == "Down-regulated")],10)
resOrdered_top10gene <- c(as.character(up_gene),as.character(down_gene))
resOrdered$Label[match(resOrdered_top10gene,resOrdered$symbol)] <- resOrdered_top10gene
ungene3 <- resOrdered[resOrdered$Group == "Up-regulated",]
library("clusterProfiler")
library("org.Hs.eg.db")
library("enrichplot")
library("ggplot2")
setwd("/media/feng/data1/骨肉瘤挖掘/肺转移/Tcell")
dir.create("/media/feng/data1/骨肉瘤挖掘/肺转移/Tcell/NKTCHOKEGG_UP")
setwd("/media/feng/data1/骨肉瘤挖掘/肺转移/Tcell/NKTCHOKEGG_UP")
genes=as.vector(ungene3$X)
entrezIDs <- mget(genes, org.Hs.egSYMBOL2EG, ifnotfound=NA)    
entrezIDs <- as.character(entrezIDs)
out=cbind(genes,entrezID=entrezIDs)
write.table(out,file="id.txt",sep="\t",quote=F,row.names=F)   
rt=read.table("id.txt",sep="\t",header=T,check.names=F)           
rt=rt[is.na(rt[,"entrezID"])==F,]                                 
gene=rt$entrezID
kk <- enrichKEGG(gene = gene, organism = "hsa", pvalueCutoff =0.05, qvalueCutoff =0.05)
write.table(kk,file="KEGG.txt",sep="\t",quote=F,row.names = F)                 

p <- barplot(kk, drop = TRUE, showCategory = 30)
ggsave("KEGGbarplot.png",p,width = 10,height = 8)
ggsave("KEGGbarplot.eps",p,width = 10,height = 8)
ggsave("KEGGbarplot.pdf",p,width = 10,height = 8)

p <- dotplot(kk, showCategory = 30)
ggsave("KEGGdotplot.png",p,width = 10,height = 8)
ggsave("KEGGdotplot.eps",p,width = 10,height = 8)
ggsave("KEGGdotplot.pdf",p,width = 10,height = 8)

setwd("/media/feng/data1/骨肉瘤挖掘/肺转移/Tcell/")
ungene4 <- resOrdered[resOrdered$Group == "Down-regulated",]
genes=as.vector(ungene4[,10])
entrezIDs <- mget(genes, org.Hs.egSYMBOL2EG, ifnotfound=NA)    
entrezIDs <- as.character(entrezIDs)
out=cbind(ungene4$symbol,entrezID=entrezIDs)
dir.create("/media/feng/data1/骨肉瘤挖掘/肺转移/Tcell/NKT_CHO_DOWN/")
setwd("/media/feng/data1/骨肉瘤挖掘/肺转移/Tcell/NKT_CHO_DOWN/")
write.table(out,file="id.txt",sep="\t",quote=F,row.names=F)   
rt=read.table("id.txt",sep="\t",header=T,check.names=F)           
rt=rt[is.na(rt[,"entrezID"])==F,]                                 
gene=rt$entrezID
kk <- enrichGO(gene = gene, OrgDb = org.Hs.eg.db, pvalueCutoff =0.05, qvalueCutoff = 0.05, ont="all", readable =T)
write.table(kk,file="GO.txt",sep="\t",quote=F,row.names = F)                 

p <-barplot(kk, drop = TRUE, showCategory =10,split="ONTOLOGY") + facet_grid(ONTOLOGY~., scale='free')
ggsave("GObarplot.png",p,width = 7,height = 8)
ggsave("GObarplot.eps",p,width = 7,height = 8)
ggsave("GObarplot.pdf",p,width = 7,height = 8)

p <- dotplot(kk,showCategory = 10,split="ONTOLOGY") + facet_grid(ONTOLOGY~., scale='free')
ggsave("GObubble.png",p,width = 7,height = 8)
ggsave("GObubble.eps",p,width = 7,height = 8)
ggsave("GObubble.pdf",p,width = 7,height = 8)

ungene4 <- resOrdered[resOrdered$Group == "Down-regulated",]
library("clusterProfiler")
library("org.Hs.eg.db")
library("enrichplot")
library("ggplot2")
setwd("/media/feng/data1/骨肉瘤挖掘/肺转移/Tcell")
dir.create("/media/feng/data1/骨肉瘤挖掘/肺转移/Tcell/NKTCHOKEGG_DOWN")
setwd("/media/feng/data1/骨肉瘤挖掘/肺转移/Tcell/NKTCHOKEGG_DOWN")
genes=as.vector(ungene4$X)
entrezIDs <- mget(genes, org.Hs.egSYMBOL2EG, ifnotfound=NA)    
entrezIDs <- as.character(entrezIDs)
out=cbind(genes,entrezID=entrezIDs)
write.table(out,file="id.txt",sep="\t",quote=F,row.names=F)   
rt=read.table("id.txt",sep="\t",header=T,check.names=F)           
rt=rt[is.na(rt[,"entrezID"])==F,]                                 
gene=rt$entrezID
kk <- enrichKEGG(gene = gene, organism = "hsa", pvalueCutoff =0.05, qvalueCutoff =0.05)
write.table(kk,file="KEGG.txt",sep="\t",quote=F,row.names = F)                 

p <- barplot(kk, drop = TRUE, showCategory = 30)
ggsave("KEGGbarplot.png",p,width = 10,height = 8)
ggsave("KEGGbarplot.eps",p,width = 10,height = 8)
ggsave("KEGGbarplot.pdf",p,width = 10,height = 8)

p <- dotplot(kk, showCategory = 30)
ggsave("KEGGdotplot.png",p,width = 10,height = 8)
ggsave("KEGGdotplot.eps",p,width = 10,height = 8)
ggsave("KEGGdotplot.pdf",p,width = 10,height = 8)

NKTcell@active.ident <- NKTcell$seurat_clusters
new.cluster.ids <- c("C1_Lung tissue","C2_Lung tissue","C3_Lung tissue","C4_Lung tissue","Lung tissue(Conventional)","Lung tissue(Chondroblastic)")
names(new.cluster.ids) <- levels(NKTcell)
NKTcell <- RenameIdents(NKTcell, new.cluster.ids)
usegene <- rbind(ungene1,ungene2,ungene3,ungene4)

Average_NKTcell<- AverageExpression(NKTcell, assays = NULL, features = NULL,return.seurat = FALSE, add.ident = NULL, slot = "data",use.scale = FALSE, use.counts = FALSE, verbose = TRUE)
gene<-Average_NKTcell$RNA
gene<-as.matrix(gene)
usegene <- usegene[,10]
gene <- as.data.frame(gene)
gene <- gene[usegene,]
t <- t(scale(t(gene)))

library(ComplexHeatmap)

heat <- Heatmap(t, 
                col = colorRampPalette(c("navy", "white", "firebrick3"))(50), 
                heatmap_legend_param = list(grid_height = unit(6,'mm')), 
                show_row_names = FALSE,  #不展示基因名称
                top_annotation = HeatmapAnnotation(Group = new.cluster.ids,
                                                   simple_anno_size = unit(2, 'mm'), 
                                                   col = list(Group = c('C1_Lung tissue' = '#E5D2DD', 
                                                                        'C2_Lung tissue' = '#53A85F', 
                                                                        'C3_Lung tissue' = '#F1BB72',
                                                                        "C4_Lung tissue" = '#F3B1A0',
                                                                        "Lung tissue(Conventional)" = '#D6E7A3',
                                                                        "Lung tissue(Chondroblastic)" = '#57C3F3')), 
                                                   show_annotation_name = FALSE), 
                column_names_gp = gpar(fontsize = 10), row_names_gp = gpar(fontsize = 6))

name <- as.data.frame(c(ungene1$symbol[1:7],ungene2$symbol[1:10],ungene3$symbol[1:10],ungene4$symbol[1:10]))
colnames(name) <- c("V1")
name <- as.data.frame(name[!duplicated(name$V1),])
heat + rowAnnotation(link = anno_mark(at = which(rownames(t) %in% name$`name[!duplicated(name$V1), ]`), 
                                      labels = name$`name[!duplicated(name$V1), ]`, labels_gp = gpar(fontsize = 6)))

gene<-as.data.frame(NKTcell@assays$RNA@counts)
usegene <- c("PDCD1","CTLA4","LAG3","TIGIT","HAVCR2","GZMM","GZMK","GZMB","GZMA")
gene <- gene[usegene,]
gene <- t(gene)
gene <- as.data.frame(gene)
gene <- log2(gene + 1)
gene$ID <- rownames(gene)
NKTcell$Celltype <- NKTcell@active.ident
meta <- as.data.frame(NKTcell@meta.data)
meta$ID <- rownames(meta)
meta <- data.frame(meta$Celltype,meta$ID)
colnames(meta) <- c("Celltype","ID")
NKTusgene <- merge(meta,gene,by = "ID")
write.csv(NKTusgene,file = "NKTusgene.csv")

library(dplyr)
library(monocle)
data <- as(as.matrix(NKTcell@assays$RNA@counts), 'sparseMatrix')
pd <- new('AnnotatedDataFrame', data = NKTcell@meta.data)
fData <- data.frame(gene_short_name = row.names(data), row.names = row.names(data))
fd <- new('AnnotatedDataFrame', data = fData)
my_cds <- newCellDataSet(data,
                         phenoData = pd,
                         featureData = fd,
                         lowerDetectionLimit = 0.5,
                         expressionFamily = negbinomial.size())
my_cds <- estimateSizeFactors(my_cds)
my_cds <- estimateDispersions(my_cds)
my_cds <- detectGenes(my_cds, min_expr = 0.1)
head(fData(my_cds))
head(pData(my_cds))
pData(my_cds)$UMI <- Matrix::colSums(exprs(my_cds))
disp_table <- dispersionTable(my_cds)
head(disp_table)
table(disp_table$mean_expression>=0.1)
unsup_clustering_genes <- subset(disp_table, mean_expression >= 0.1)
my_cds <- setOrderingFilter(my_cds, unsup_clustering_genes$gene_id)
plot_ordering_genes(my_cds)
expressed_genes <- row.names(subset(fData(my_cds), num_cells_expressed >= 10))
my_cds_subset <- my_cds[expressed_genes, ]

head(pData(my_cds_subset))
my_cds_subset <- detectGenes(my_cds_subset, min_expr = 0.1)
fData(my_cds_subset)$use_for_ordering <- fData(my_cds_subset)$num_cells_expressed > 0.05 * ncol(my_cds_subset)
table(fData(my_cds_subset)$use_for_ordering)
plot_pc_variance_explained(my_cds_subset, return_all = FALSE)
my_cds_subset <- reduceDimension(my_cds_subset,max_components = 2,norm_method = 'log',num_dim = 10,reduction_method = 'tSNE',verbose = TRUE)
my_cds_subset <- clusterCells(my_cds_subset, verbose = FALSE)
plot_rho_delta(my_cds_subset, rho_threshold = 2, delta_threshold = 10)
my_cds_subset <- clusterCells(my_cds_subset,rho_threshold = 2,delta_threshold = 10,skip_rho_sigma = T,verbose = FALSE)
table(pData(my_cds_subset)$Cluster)
plot_cell_clusters(my_cds_subset)
head(pData(my_cds_subset))
clustering_DEG_genes <- differentialGeneTest(my_cds_subset,fullModelFormulaStr = '~Cluster',cores = 40)
dim(clustering_DEG_genes)

clustering_DEG_genes %>% arrange(qval) %>% head()
my_ordering_genes <- row.names(clustering_DEG_genes)[order(clustering_DEG_genes$qval)][1:1000]
my_cds_subset <- setOrderingFilter(my_cds_subset, ordering_genes = my_ordering_genes)
my_cds_subset <- reduceDimension(my_cds_subset, method = 'DDRTree')
my_cds_subset <- orderCells(my_cds_subset)

setwd("/media/feng/data1/骨肉瘤挖掘/肺转移/Tcell")
my_cds_subset$Celltype <- NKTcell@active.ident
p <- plot_cell_trajectory(my_cds_subset, color_by = "State")
ggsave("monocle_state.png",p)
ggsave("monocle_state.eps",p)
ggsave("monocle_state.pdf",p)

p <- plot_cell_trajectory(my_cds_subset, color_by = "Pseudotime")
ggsave("monocle_Pseudotime.png",p)
ggsave("monocle_Pseudotime.eps",p)
ggsave("monocle_Pseudotime.pdf",p)

p <- plot_cell_trajectory(my_cds_subset, color_by = "Celltype")
ggsave("monocle_Celltype.png",p)
ggsave("monocle_Celltype.eps",p)
ggsave("monocle_Celltype.pdf",p)

p <- plot_cell_trajectory(my_cds_subset, color_by = "Donor")
ggsave("monocle_Donor.png",p)
ggsave("monocle_Donor.eps",p)
ggsave("monocle_Donor.pdf",p)

p <- plot_cell_trajectory(my_cds_subset,markers = c("JUND"),use_color_gradient = TRUE)
ggsave("monocle_JUND.png",p)
ggsave("monocle_JUND.eps",p)
ggsave("monocle_JUND.pdf",p)

p <- plot_cell_trajectory(my_cds_subset,markers = c("GZMK"),use_color_gradient = TRUE)
ggsave("monocle_GZMK.png",p)
ggsave("monocle_GZMK.eps",p)
ggsave("monocle_GZMK.pdf",p)

p <- plot_cell_trajectory(my_cds_subset,markers = c("GZMM"),use_color_gradient = TRUE)
ggsave("monocle_GZMM.png",p)
ggsave("monocle_GZMM.eps",p)
ggsave("monocle_GZMM.pdf",p)

p <- plot_cell_trajectory(my_cds_subset,markers = c("GZMB"),use_color_gradient = TRUE)
ggsave("monocle_GZMB.png",p)
ggsave("monocle_GZMB.eps",p)
ggsave("monocle_GZMB.pdf",p)

p <- plot_cell_trajectory(my_cds_subset,markers = c("GZMA"),use_color_gradient = TRUE)
ggsave("monocle_GZMA.png",p)
ggsave("monocle_GZMA.eps",p)
ggsave("monocle_GZMA.pdf",p)


CD8T <- subset(NKTcell,subset = CD8A > 0 | CD8B > 0 )
cd8t <- as.data.frame(table(CD8T$Donor))
nkt <- as.data.frame(table(NKTcell$Donor))
cd8t$percentage <- cd8t/nkt *100

setwd("/media/feng/data1/骨肉瘤挖掘/肺转移/Mac")
Mac <- subset(OS_lung,idents = c("Macrophages"))
Mac <- NormalizeData(Mac) %>% FindVariableFeatures() %>%ScaleData(features = rownames(Mac))
Mac <- RunPCA(Mac, npcs=50, verbose=FALSE)
DimPlot(Mac, reduction = "pca", group.by = "Donor")

Mac <- RunHarmony(Mac, group.by.vars="Donor",max.iter.harmony = 15)
ElbowPlot(Mac, ndims = 30)

pc.num=1:20  
Mac <- Mac %>%
  RunUMAP(reduction = "harmony",dims=pc.num) %>%
  RunTSNE(reduction = "harmony",dims=pc.num) %>%
  FindNeighbors(dims = pc.num) %>%
  FindClusters(resolution = 0.03) %>%
  identity()

p <- DimPlot(Mac, group.by = "Donor") + scale_color_ucscgb()
ggsave("Macumap_Donor.png",p)
ggsave("Macumap_Donor.eps",p)
ggsave("Macumap_Donor.pdf",p)
p <- DimPlot(Mac) + scale_color_jama()
ggsave("Macumap_source.png",p)
ggsave("Macumap_source.eps",p)
ggsave("Macumap_source.pdf",p)
p <- DimPlot(Mac,group.by = "bingli") + scale_color_jama()
ggsave("Macumap_bingli.png",p)
ggsave("Macumap_bingli.eps",p)
ggsave("Macumap_bingli.pdf",p)

Mac@active.ident <- Mac$seurat_clusters
new.cluster.ids <- c("Lung tissue","Lung tissue","Lung tissue","Lung tissue(Conventional)","Lung tissue","Lung tissue","Lung tissue(Chondroblastic)")
names(new.cluster.ids) <- levels(Mac)
Mac <- RenameIdents(Mac, new.cluster.ids)

p <- DimPlot(Mac) + scale_color_aaas()

ggsave("Macumap.png",p)
ggsave("Macumap.eps",p)
ggsave("Macumap.pdf",p)

p <- TSNEPlot(Mac) + scale_color_aaas()

ggsave("MacTSNE.png",p)
ggsave("MacTSNE.eps",p)
ggsave("MacTSNE.pdf",p)


options(stringsAsFactors = F)
library(DESeq2)
library(limma)
library(stringr)
library(ggsci)
library(ggpubr)
Mac_CON <- subset(Mac,idents = c("Lung tissue","Lung tissue(Conventional)"))
tmp_singlecell <- as.data.frame(Mac_CON@assays$RNA@counts)
tmp_singlecell <- tmp_singlecell[rowMeans(tmp_singlecell) >0,]
tmp_singlecell <- tmp_singlecell +1
k <- tmp_singlecell[1:10,1:10]
group_singlecell <- as.data.frame(Mac_CON@active.ident)
group_singlecell$ID <- rownames(group_singlecell)


group_list=as.factor(group_singlecell$`Mac_CON@active.ident`)

colData <- data.frame(row.names=group_singlecell$ID, group_list=group_list)
dds <- DESeqDataSetFromMatrix(countData = ceiling(tmp_singlecell[,rownames(colData)]),
                              colData = colData,
                              design = ~ group_list)
dds2 <- DESeq(dds)  
resultsNames(dds2)
res <-  results(dds2, contrast=c("group_list","Lung tissue(Conventional)","Lung tissue"))#1是转移，0是不转移
resOrdered <- res[order(res$padj),]
resOrdered=as.data.frame(resOrdered)
resOrdered = na.omit(resOrdered)
write.csv(resOrdered,file = "Mac_CONsinglecellecDEseq2DEG.csv")

head(resOrdered)
resOrdered$logP <- -log10(resOrdered$padj)
resOrdered$Group <-"Not-significant"
resOrdered$Group[which((resOrdered$padj <(0.05))&(resOrdered$log2FoldChange >(1)))] <-"Up-regulated"
resOrdered$Group[which((resOrdered$padj <(0.05))&(resOrdered$log2FoldChange <(-1)))] <-"Down-regulated"
table(resOrdered$Group)
resOrdered$Label = ""
resOrdered <- resOrdered[order(resOrdered$padj),]
resOrdered$symbol <- rownames(resOrdered)

up_gene <- head(resOrdered$symbol[which(resOrdered$Group == "Up-regulated")],10)
down_gene <- head(resOrdered$symbol[which(resOrdered$Group == "Down-regulated")],10)
resOrdered_top10gene <- c(as.character(up_gene),as.character(down_gene))
resOrdered$Label[match(resOrdered_top10gene,resOrdered$symbol)] <- resOrdered_top10gene
p <-ggscatter(resOrdered,x = "log2FoldChange",y = "logP",
              color = "Group",
              palette = c("#2f5688","#BBBBBB","#CC0000"),
              size = 1,
              label = resOrdered$Label,
              font.label = 8,
              repel = T,
              xlab = "log2FoldChange",
              ylab = "-log10(P-value)")  + geom_hline(yintercept = 1.3, linetype = "dashed")+
  geom_vline(xintercept = c(-1,1), linetype = "dashed")
ggsave("Mac经典型对比火山图.png",p)
ggsave("Mac经典型对比火山图.pdf",p)
ggsave("Mac经典型对比火山图.eps",p)


ungene1 <- resOrdered[resOrdered$Group == "Up-regulated",]
library("clusterProfiler")
library("org.Hs.eg.db")
library("enrichplot")
library("ggplot2")
genes=as.vector(ungene1[,10])
entrezIDs <- mget(genes, org.Hs.egSYMBOL2EG, ifnotfound=NA)    
entrezIDs <- as.character(entrezIDs)
out=cbind(ungene1$symbol,entrezID=entrezIDs)
dir.create("/media/feng/data1/骨肉瘤挖掘/肺转移/Mac/Mac_CON_UP/")
setwd("/media/feng/data1/骨肉瘤挖掘/肺转移/Mac/Mac_CON_UP/")
write.table(out,file="id.txt",sep="\t",quote=F,row.names=F)   
rt=read.table("id.txt",sep="\t",header=T,check.names=F)           
rt=rt[is.na(rt[,"entrezID"])==F,]                                 
gene=rt$entrezID
kk <- enrichGO(gene = gene, OrgDb = org.Hs.eg.db, pvalueCutoff =0.05, qvalueCutoff = 0.05, ont="all", readable =T)
write.table(kk,file="GO.txt",sep="\t",quote=F,row.names = F)                 

p <-barplot(kk, drop = TRUE, showCategory =10,split="ONTOLOGY") + facet_grid(ONTOLOGY~., scale='free')
ggsave("GObarplot.png",p,width = 11,height = 8)
ggsave("GObarplot.eps",p,width = 11,height = 8)
ggsave("GObarplot.pdf",p,width = 11,height = 8)

p <- dotplot(kk,showCategory = 10,split="ONTOLOGY") + facet_grid(ONTOLOGY~., scale='free')
ggsave("GObubble.png",p,width = 11,height = 8)
ggsave("GObubble.eps",p,width = 11,height = 8)
ggsave("GObubble.pdf",p,width = 11,height = 8)

ungene1 <- resOrdered[resOrdered$Group == "Up-regulated",]
library("clusterProfiler")
library("org.Hs.eg.db")
library("enrichplot")
library("ggplot2")
setwd("/media/feng/data1/骨肉瘤挖掘/肺转移/Mac")
dir.create("/media/feng/data1/骨肉瘤挖掘/肺转移/Mac/MacCHOKEGG_UP")
setwd("/media/feng/data1/骨肉瘤挖掘/肺转移/Mac/MacCHOKEGG_UP")
ungene1$X <- rownames(ungene1)
genes=as.vector(ungene1$X)
entrezIDs <- mget(genes, org.Hs.egSYMBOL2EG, ifnotfound=NA)    
entrezIDs <- as.character(entrezIDs)
out=cbind(genes,entrezID=entrezIDs)
write.table(out,file="id.txt",sep="\t",quote=F,row.names=F)   
rt=read.table("id.txt",sep="\t",header=T,check.names=F)           
rt=rt[is.na(rt[,"entrezID"])==F,]                                 
gene=rt$entrezID
kk <- enrichKEGG(gene = gene, organism = "hsa", pvalueCutoff =0.05, qvalueCutoff =0.05)
write.table(kk,file="KEGG.txt",sep="\t",quote=F,row.names = F)                 
p <- barplot(kk, drop = TRUE, showCategory = 30)
ggsave("KEGGbarplot.png",p,width = 10,height = 8)
ggsave("KEGGbarplot.eps",p,width = 10,height = 8)
ggsave("KEGGbarplot.pdf",p,width = 10,height = 8)

p <- dotplot(kk, showCategory = 30)
ggsave("KEGGdotplot.png",p,width = 10,height = 8)
ggsave("KEGGdotplot.eps",p,width = 10,height = 8)
ggsave("KEGGdotplot.pdf",p,width = 10,height = 8)

setwd("/media/feng/data1/骨肉瘤挖掘/肺转移/Mac/")
ungene2 <- resOrdered[resOrdered$Group == "Down-regulated",]
genes=as.vector(ungene2[,10])
entrezIDs <- mget(genes, org.Hs.egSYMBOL2EG, ifnotfound=NA)    
entrezIDs <- as.character(entrezIDs)
out=cbind(ungene1$symbol,entrezID=entrezIDs)
dir.create("/media/feng/data1/骨肉瘤挖掘/肺转移/Mac/Mac_CON_DOWN/")
setwd("/media/feng/data1/骨肉瘤挖掘/肺转移/Mac/Mac_CON_DOWN/")
write.table(out,file="id.txt",sep="\t",quote=F,row.names=F)   
rt=read.table("id.txt",sep="\t",header=T,check.names=F)           
rt=rt[is.na(rt[,"entrezID"])==F,]                                 
gene=rt$entrezID
kk <- enrichGO(gene = gene, OrgDb = org.Hs.eg.db, pvalueCutoff =0.05, qvalueCutoff = 0.05, ont="all", readable =T)
write.table(kk,file="GO.txt",sep="\t",quote=F,row.names = F)                 

p <-barplot(kk, drop = TRUE, showCategory =10,split="ONTOLOGY") + facet_grid(ONTOLOGY~., scale='free')
ggsave("GObarplot.png",p,width = 11,height = 8)
ggsave("GObarplot.eps",p,width = 11,height = 8)
ggsave("GObarplot.pdf",p,width = 11,height = 8)

p <- dotplot(kk,showCategory = 10,split="ONTOLOGY") + facet_grid(ONTOLOGY~., scale='free')
ggsave("GObubble.png",p,width = 11,height = 8)
ggsave("GObubble.eps",p,width = 11,height = 8)
ggsave("GObubble.pdf",p,width = 11,height = 8)

ungene2 <- resOrdered[resOrdered$Group == "Down-regulated",]
library("clusterProfiler")
library("org.Hs.eg.db")
library("enrichplot")
library("ggplot2")
setwd("/media/feng/data1/骨肉瘤挖掘/肺转移/Mac")
dir.create("/media/feng/data1/骨肉瘤挖掘/肺转移/Mac/MacCHOKEGG_DOWN")
setwd("/media/feng/data1/骨肉瘤挖掘/肺转移/Mac/MacCHOKEGG_DOWN")
ungene2$X <- rownames(ungene2)
genes=as.vector(ungene2$X)
entrezIDs <- mget(genes, org.Hs.egSYMBOL2EG, ifnotfound=NA)    
entrezIDs <- as.character(entrezIDs)
out=cbind(genes,entrezID=entrezIDs)
write.table(out,file="id.txt",sep="\t",quote=F,row.names=F)   
rt=read.table("id.txt",sep="\t",header=T,check.names=F)           
rt=rt[is.na(rt[,"entrezID"])==F,]                                 
gene=rt$entrezID
kk <- enrichKEGG(gene = gene, organism = "hsa", pvalueCutoff =0.05, qvalueCutoff =0.05)
write.table(kk,file="KEGG.txt",sep="\t",quote=F,row.names = F)                 
p <- barplot(kk, drop = TRUE, showCategory = 30)
ggsave("KEGGbarplot.png",p,width = 10,height = 8)
ggsave("KEGGbarplot.eps",p,width = 10,height = 8)
ggsave("KEGGbarplot.pdf",p,width = 10,height = 8)

p <- dotplot(kk, showCategory = 30)
ggsave("KEGGdotplot.png",p,width = 10,height = 8)
ggsave("KEGGdotplot.eps",p,width = 10,height = 8)
ggsave("KEGGdotplot.pdf",p,width = 10,height = 8)

Mac_CHO <- subset(Mac,idents = c("Lung tissue","Lung tissue(Chondroblastic)"))


tmp_singlecell <- as.data.frame(Mac_CHO@assays$RNA@counts)
tmp_singlecell <- tmp_singlecell[rowMeans(tmp_singlecell) >0,]
tmp_singlecell <- tmp_singlecell +1
k <- tmp_singlecell[1:10,1:10]
group_singlecell <- as.data.frame(Mac_CHO@active.ident)
group_singlecell$ID <- rownames(group_singlecell)

group_list=as.factor(group_singlecell$`Mac_CHO@active.ident`)#做好分组因子

colData <- data.frame(row.names=group_singlecell$ID, group_list=group_list)
dds <- DESeqDataSetFromMatrix(countData = ceiling(tmp_singlecell[,rownames(colData)]),
                              colData = colData,
                              design = ~ group_list)
dds2 <- DESeq(dds)  
resultsNames(dds2)
res <-  results(dds2, contrast=c("group_list","Lung tissue(Chondroblastic)","Lung tissue"))#1是转移，0是不转移
resOrdered <- res[order(res$padj),]
resOrdered=as.data.frame(resOrdered)
resOrdered = na.omit(resOrdered)
write.csv(resOrdered,file = "Mac_CHOsinglecellecDEseq2DEG.csv")

head(resOrdered)
resOrdered$logP <- -log10(resOrdered$padj)
resOrdered$Group <-"Not-significant"
resOrdered$Group[which((resOrdered$padj <(0.05))&(resOrdered$log2FoldChange >(1)))] <-"Up-regulated"
resOrdered$Group[which((resOrdered$padj <(0.05))&(resOrdered$log2FoldChange <(-1)))] <-"Down-regulated"
table(resOrdered$Group)
resOrdered$Label = ""
resOrdered <- resOrdered[order(resOrdered$padj),]
resOrdered$symbol <- rownames(resOrdered)

up_gene <- head(resOrdered$symbol[which(resOrdered$Group == "Up-regulated")],10)
down_gene <- head(resOrdered$symbol[which(resOrdered$Group == "Down-regulated")],10)
resOrdered_top10gene <- c(as.character(up_gene),as.character(down_gene))
resOrdered$Label[match(resOrdered_top10gene,resOrdered$symbol)] <- resOrdered_top10gene
p <-ggscatter(resOrdered,x = "log2FoldChange",y = "logP",
              color = "Group",
              palette = c("#2f5688","#BBBBBB","#CC0000"),
              size = 1,
              label = resOrdered$Label,
              font.label = 8,
              repel = T,
              xlab = "log2FoldChange",
              ylab = "-log10(P-value)")  + geom_hline(yintercept = 1.3, linetype = "dashed")+
  geom_vline(xintercept = c(-1,1), linetype = "dashed")
ggsave("Mac成软骨型对比火山图.png",p)
ggsave("Mac成软骨型对比火山图.pdf",p)
ggsave("Mac成软骨型对比火山图.eps",p)

ungene3 <- resOrdered[resOrdered$Group == "Up-regulated",]
library("clusterProfiler")
library("org.Hs.eg.db")
library("enrichplot")
library("ggplot2")
genes=as.vector(ungene3[,10])
entrezIDs <- mget(genes, org.Hs.egSYMBOL2EG, ifnotfound=NA)    
entrezIDs <- as.character(entrezIDs)
out=cbind(ungene3$symbol,entrezID=entrezIDs)
dir.create("/media/feng/data1/骨肉瘤挖掘/肺转移/Mac/Mac_CHO_UP/")
setwd("/media/feng/data1/骨肉瘤挖掘/肺转移/Mac/Mac_CHO_UP/")
write.table(out,file="id.txt",sep="\t",quote=F,row.names=F)   
rt=read.table("id.txt",sep="\t",header=T,check.names=F)           
rt=rt[is.na(rt[,"entrezID"])==F,]                                 
gene=rt$entrezID
kk <- enrichGO(gene = gene, OrgDb = org.Hs.eg.db, pvalueCutoff =0.05, qvalueCutoff = 0.05, ont="all", readable =T)
write.table(kk,file="GO.txt",sep="\t",quote=F,row.names = F)                 

p <-barplot(kk, drop = TRUE, showCategory =10,split="ONTOLOGY") + facet_grid(ONTOLOGY~., scale='free')
ggsave("GObarplot.png",p,width = 11,height = 8)
ggsave("GObarplot.eps",p,width = 11,height = 8)
ggsave("GObarplot.pdf",p,width = 11,height = 8)

p <- dotplot(kk,showCategory = 10,split="ONTOLOGY") + facet_grid(ONTOLOGY~., scale='free')
ggsave("GObubble.png",p,width = 11,height = 8)
ggsave("GObubble.eps",p,width = 11,height = 8)
ggsave("GObubble.pdf",p,width = 11,height = 8)

ungene3 <- resOrdered[resOrdered$Group == "Up-regulated",]
library("clusterProfiler")
library("org.Hs.eg.db")
library("enrichplot")
library("ggplot2")
setwd("/media/feng/data1/骨肉瘤挖掘/肺转移/Mac")
dir.create("/media/feng/data1/骨肉瘤挖掘/肺转移/Mac/MacCHOKEGG_UP")
setwd("/media/feng/data1/骨肉瘤挖掘/肺转移/Mac/MacCHOKEGG_UP")
ungene3$X <- rownames(ungene3)
genes=as.vector(ungene3$X)
entrezIDs <- mget(genes, org.Hs.egSYMBOL2EG, ifnotfound=NA)    
entrezIDs <- as.character(entrezIDs)
out=cbind(genes,entrezID=entrezIDs)
write.table(out,file="id.txt",sep="\t",quote=F,row.names=F)   
rt=read.table("id.txt",sep="\t",header=T,check.names=F)           
rt=rt[is.na(rt[,"entrezID"])==F,]                                 
gene=rt$entrezID
kk <- enrichKEGG(gene = gene, organism = "hsa", pvalueCutoff =0.05, qvalueCutoff =0.05)
write.table(kk,file="KEGG.txt",sep="\t",quote=F,row.names = F)                 
p <- barplot(kk, drop = TRUE, showCategory = 30)
ggsave("KEGGbarplot.png",p,width = 10,height = 8)
ggsave("KEGGbarplot.eps",p,width = 10,height = 8)
ggsave("KEGGbarplot.pdf",p,width = 10,height = 8)

p <- dotplot(kk, showCategory = 30)
ggsave("KEGGdotplot.png",p,width = 10,height = 8)
ggsave("KEGGdotplot.eps",p,width = 10,height = 8)
ggsave("KEGGdotplot.pdf",p,width = 10,height = 8)

setwd("/media/feng/data1/骨肉瘤挖掘/肺转移/Mac/")
ungene4 <- resOrdered[resOrdered$Group == "Down-regulated",]
genes=as.vector(ungene4[,10])
entrezIDs <- mget(genes, org.Hs.egSYMBOL2EG, ifnotfound=NA)    
entrezIDs <- as.character(entrezIDs)
out=cbind(ungene4$symbol,entrezID=entrezIDs)
dir.create("/media/feng/data1/骨肉瘤挖掘/肺转移/Mac/Mac_CHO_DOWN/")
setwd("/media/feng/data1/骨肉瘤挖掘/肺转移/Mac/Mac_CHO_DOWN/")
write.table(out,file="id.txt",sep="\t",quote=F,row.names=F)   
rt=read.table("id.txt",sep="\t",header=T,check.names=F)           
rt=rt[is.na(rt[,"entrezID"])==F,]                                 
gene=rt$entrezID
kk <- enrichGO(gene = gene, OrgDb = org.Hs.eg.db, pvalueCutoff =0.05, qvalueCutoff = 0.05, ont="all", readable =T)
write.table(kk,file="GO.txt",sep="\t",quote=F,row.names = F)                 

p <-barplot(kk, drop = TRUE, showCategory =10,split="ONTOLOGY") + facet_grid(ONTOLOGY~., scale='free')
ggsave("GObarplot.png",p,width = 10,height = 8)
ggsave("GObarplot.eps",p,width = 10,height = 8)
ggsave("GObarplot.pdf",p,width = 10,height = 8)

p <- dotplot(kk,showCategory = 10,split="ONTOLOGY") + facet_grid(ONTOLOGY~., scale='free')
ggsave("GObubble.png",p,width = 10,height = 8)
ggsave("GObubble.eps",p,width = 10,height = 8)
ggsave("GObubble.pdf",p,width = 10,height = 8)

ungene4 <- resOrdered[resOrdered$Group == "Down-regulated",]
library("clusterProfiler")
library("org.Hs.eg.db")
library("enrichplot")
library("ggplot2")
setwd("/media/feng/data1/骨肉瘤挖掘/肺转移/Mac")
dir.create("/media/feng/data1/骨肉瘤挖掘/肺转移/Mac/MacCHOKEGG_DOWN")
setwd("/media/feng/data1/骨肉瘤挖掘/肺转移/Mac/MacCHOKEGG_DOWN")
ungene4$X <- rownames(ungene4)
genes=as.vector(ungene4$X)
entrezIDs <- mget(genes, org.Hs.egSYMBOL2EG, ifnotfound=NA)    
entrezIDs <- as.character(entrezIDs)
out=cbind(genes,entrezID=entrezIDs)
write.table(out,file="id.txt",sep="\t",quote=F,row.names=F)   
rt=read.table("id.txt",sep="\t",header=T,check.names=F)           
rt=rt[is.na(rt[,"entrezID"])==F,]                                 
gene=rt$entrezID
kk <- enrichKEGG(gene = gene, organism = "hsa", pvalueCutoff =0.05, qvalueCutoff =0.05)
write.table(kk,file="KEGG.txt",sep="\t",quote=F,row.names = F)                 
p <- barplot(kk, drop = TRUE, showCategory = 30)
ggsave("KEGGbarplot.png",p,width = 10,height = 8)
ggsave("KEGGbarplot.eps",p,width = 10,height = 8)
ggsave("KEGGbarplot.pdf",p,width = 10,height = 8)

p <- dotplot(kk, showCategory = 30)
ggsave("KEGGdotplot.png",p,width = 10,height = 8)
ggsave("KEGGdotplot.eps",p,width = 10,height = 8)
ggsave("KEGGdotplot.pdf",p,width = 10,height = 8)

dir.create("/media/feng/data1/骨肉瘤挖掘/肺转移/Mac/GSVA/")
setwd("/media/feng/data1/骨肉瘤挖掘/肺转移/Mac/GSVA/")

library(dplyr)
library(GSVA)
library(GSEABase)
library(limma)
library(Seurat)
library(pheatmap)
library(stringr)
library(ggplot2)
library(psych)
library(Matrix)

genesets <- getGmt("/media/feng/data1/GCTB/GCTB20210406/all/all/MY-he/msigdb.v7.2.symbols.xls") 
Average_MY<- AverageExpression(Mac, assays = NULL, features = NULL,return.seurat = FALSE, add.ident = NULL, slot = "data",use.scale = FALSE, use.counts = FALSE, verbose = TRUE)
gene<-Average_MY$RNA
gene<-as.matrix(gene)
GSVAresult <- gsva(gene, genesets, min.sz=10, max.sz=Inf, tau=1, method="gsva", kcdf="Poisson", mx.diff=TRUE, abs.ranking=FALSE, verbose=TRUE, parallel.sz=10)
t <- t(scale(t(GSVAresult)))
write.csv(t,file = "Mac_GSVAresult.csv")

t <- read.table("MacGSVA.txt",header = T,check.names = F, sep = "\t",row.names = 1)
p = pheatmap(t, cluster_cols = F, color = colorRampPalette(c("navy", "white", "firebrick3"))(50),fontsize = 5,border_color = "white")
ggsave("GSVAmarkerheatmap.png",p,width = 4.5,height = 2.5)
ggsave("GSVAmarkerheatmap.eps",p,width = 4.5,height = 2.5)
ggsave("GSVAmarkerheatmap.pdf",p,width = 4.5,height = 2.5)


usegene <- c("IL1B", 
             "IL1A", 
             "TNF", 
             "IL6", 
             "CXCL9", 
             "CXCL10", 
             "IL12A", 
             "IL12B", 
             "IL23A", 
             "FCGR1A", 
             "FCGR1B", 
             "FCGR1C", 
             "CCR7", 
             "IL8", 
             "CCL5", 
             "HLA-DRA", 
             "IRF5", 
             "IRF1")
M1 <- as.data.frame(usegene)

usegene <- c("IL10", 
             "CD163", 
             "MARCO", 
             "MRC1", 
             "MSR1", 
             "ARG1", 
             "STAB1", 
             "TGM2", 
             "MMP7", 
             "MMP9", 
             "MMP19", 
             "TGFB1", 
             "TGFB2", 
             "TGFB3", 
             "VEGFA", 
             "FN1", 
             "CCL4", 
             "CCL22", 
             "CCL17", 
             "CCL18", 
             "IL4R", 
             "IL7R", 
             "IRF4")
M2 <- as.data.frame(usegene)

Mac <- AddModuleScore(
  object = Mac,
  features = M1,
  ctrl = 5,
  name = 'M1')
Mac <- AddModuleScore(
  object = Mac,
  features = M2,
  ctrl = 5,
  name = 'M2')
Mac$Celltype <- Mac@active.ident
meta <- as.data.frame(Mac@meta.data)
meta$ID <- rownames(meta)
M12scores <- data.frame(meta$ID,meta$M11,meta$M21,meta$Celltype)
colnames(M12scores) <- c("ID","M1_gene_score","M2_gene_score","Source")
write.csv(M12scores,file = "M12scores.csv")

p <- ggplot(M12scores, aes(x = M1_gene_score, y = M2_gene_score, colour = Source)) +geom_point()+theme_bw()+theme(panel.grid=element_blank(),panel.border=element_blank(),axis.line=element_line(size=0.5,colour="black"))+ scale_color_aaas()
ggsave("M12scores.png",p)
ggsave("M12scores.eps",p)
ggsave("M12scores.pdf",p)

Average_MY<- AverageExpression(Mac, assays = NULL, features = NULL,return.seurat = FALSE, add.ident = NULL, slot = "data",use.scale = FALSE, use.counts = FALSE, verbose = TRUE)
gene<-Average_MY$RNA
gene<-as.matrix(gene)
usegene = c("GPNMB","MSR1","AXL","CCL2","MERTK","MRC1","CD163L1","STAB1","F13A1","FOLR2","CD163","CXCL9","CXCL10","TNF","IL10","CCL3","CCL4","CHI3L1","IL1B","IL1RN","CXCL3","CXCL2")
gene = gene[usegene,]
t <- t(scale(t(gene)))
annotation_row = data.frame(
  GeneClass = factor(c("Anti-inflammatory","Anti-inflammatory",
                       "Anti-inflammatory","Pro-inflammatory",
                       "Anti-inflammatory","Anti-inflammatory",
                       "Anti-inflammatory","Anti-inflammatory",
                       "Anti-inflammatory","Anti-inflammatory",
                       "Anti-inflammatory","Pro-inflammatory",
                       "Pro-inflammatory","Pro-inflammatory",
                       "Anti-inflammatory","Pro-inflammatory",
                       "Pro-inflammatory","Anti-inflammatory",
                       "Pro-inflammatory","Anti-inflammatory",
                       "Pro-inflammatory","Pro-inflammatory")))
rownames(annotation_row) <- rownames(gene)
p = pheatmap(t, cluster_cols = F,cluster_rows = F, color = colorRampPalette(c("navy", "white", "firebrick3"))(50),fontsize = 5,border_color = "white",annotation_row = annotation_row)
ggsave("inflammarkerheatmap.png",p,width = 3,height = 4)
ggsave("inflammarkerheatmap.eps",p,width = 3,height = 4)
ggsave("inflammarkerheatmap.pdf",p,width = 3,height = 4)

usegene = c("HLA-DOB","HLA-E","HLA-B","HLA-F","HLA-C","HLA-DQB2","HLA-G","HLA-DRB5","HLA-A","HLA-DMA","HLA-DMB","HLA-DOA","HLA-DQA1","HLA-DQB1","HLA-DQA2","HLA-DPA1","HLA-DPB1","HLA-DRA","HLA-DRB1")
gene = gene[usegene,]
t <- t(scale(t(gene)))
p = pheatmap(t, cluster_cols = F,cluster_rows = F, color = colorRampPalette(c("navy", "white", "firebrick3"))(50),fontsize = 5,border_color = "white")

ggsave("HLAmarkerheatmap.png",p,width = 2,height = 4)
ggsave("HLAmarkerheatmap.eps",p,width = 2,height = 4)
ggsave("HLAmarkerheatmap.pdf",p,width = 2,height = 4)






library(Seurat)
library(tidyverse)
library(patchwork)
setwd("/media/feng/data1/骨肉瘤挖掘/肺转移/")
dir.create("cellphone")
setwd("/media/feng/data1/骨肉瘤挖掘/肺转移/cellphone/")
geneID_10x <- read.table('features.tsv', header=F, sep='\t')
geneID_10x<-geneID_10x[,-3]
names(geneID_10x) <- c('Ensembl','Gene')
OS <- subset(OS_lung,idents = c("Chondroblasts"))
OSNKT <- merge(OS,NKTcell)
sp1<-OSNKT
sp1[["celltype"]]<-sp1@active.ident
sp1_counts <- data.frame(sp1@assays$RNA@data, check.names = F,stringsAsFactors = F)
sp1_counts <- data.frame(Gene=rownames(sp1_counts), sp1_counts,check.names = F,stringsAsFactors = F)
sp1_counts <- inner_join(geneID_10x,sp1_counts)
sp1_counts <- sp1_counts[,-which(names(sp1_counts)=="Gene")]
sp1_meta <- data.frame(Cell=rownames(sp1@meta.data), cell_type=sp1@meta.data$celltype,stringsAsFactors = F)
write.table(sp1_counts, "sp1_counts.txt", sep='\t',row.names = F)
write.table(sp1_meta, "sp1_meta.txt", sep='\t')


Average_T<- AverageExpression(Tcell, assays = NULL, features = NULL,return.seurat = FALSE, add.ident = NULL, slot = "data",use.scale = FALSE, use.counts = FALSE, verbose = TRUE)
gene<-Average_T$RNA
gene<-as.matrix(gene)
usegene <- c("CD8B","CD8A","CD3G","CD3D","CD3E","CD4","CTLA4","TIGIT","LAG3")
gene <- as.data.frame(gene)
gene <- gene[rownames(gene)%in%usegene,]
t <- t(scale(t(gene)))

p = pheatmap(t, cluster_cols = F,cluster_rows = T, color = colorRampPalette(c("navy", "white", "firebrick3"))(50),fontsize = 5,border_color = "white")
ggsave("Tcellmarkerheatmap.png",p,width = 2.5,height = 2)
ggsave("Tcellmarkerheatmap.eps",p,width = 2.5,height = 2)
ggsave("Tcellmarkerheatmap.pdf",p,width = 2.5,height = 2)





library("beyondcell")
library("Seurat")
library("ggplot2")
set.seed(1)
OS_lung@active.ident <- OS_lung$seurat_clusters
new.cluster.ids <- c("T cells", "T cells", "Malignant cells",	"Monocytes",	"Alveolar cells",	"T cells",	"Malignant cells",	"Osteoclasts",	"Endothelial",	"Mast cells",	"Proliferative cells", 	"Macrophages",	"Fibroblasts",	"Fibroblasts",	"Dendritic cells",	"T cells",	"B cells",	"Muscle cells",	"Secretory cells",	"Fibroblasts",	"Ciliated cells",	"Low quality cells",	"Plasma cells",	"T cells",	"Proliferative cells", "T cells"
)
names(new.cluster.ids) <- levels(OS_lung)

OS_lung <- RenameIdents(OS_lung, new.cluster.ids)
MC <- subset(OS_lung,idents = c("Malignant cells"))
sc <- MC
DefaultAssay(sc) <- "RNA"
gs <- GenerateGenesets(PSc)
nopath <- GenerateGenesets(PSc, include.pathways = FALSE)
bc <- bcScore(sc, gs, expr.thres = 0.1) 
bc <- bcUMAP(bc, k.neighbors = 4, res = 0.1, method = "umap-learn")
bc <- bcUMAP(bc, pc = 10, k.neighbors = 4, res = 0.01, method = "umap-learn")
bcClusters(bc, UMAP = "beyondcell", idents = "nFeature_RNA", factor.col = FALSE)
bcClusters(bc, UMAP = "beyondcell", idents = "Phase", factor.col = TRUE)
bc <- bcUMAP(bc, pc = 10, k.neighbors = 20, res = 0.01, method = "umap-learn")
bcClusters(bc, UMAP = "beyondcell", idents = "nFeature_RNA", factor.col = FALSE)
bcClusters(bc, UMAP = "beyondcell", idents = "bc_clusters_res.0.01") 

meta <- bc@meta.data
meta$TC <- meta$bc_clusters_res.0.01
meta$TC <- gsub("0","TC0",meta$TC)
meta$TC <- gsub("1","TC1",meta$TC)
meta$TC <- gsub("2","TC2",meta$TC)
meta$TC <- gsub("3","TC3",meta$TC)
meta$TC <- gsub("4","TC4",meta$TC)
bc@ranks$TC <- meta$TC
bc@meta.data$TC <- meta$TC

bc <- bcRanks(bc)
bc <- bcRanks(bc, idents = "TC")
bc <- bcRanks(bc, idents = "bc_clusters_res.0.01", extended = FALSE)

general <- bc@ranks$general
TC <- bc@ranks$TC

bc4Squares(bc, idents = "TC", lvl = "TC0", top = 5)



FindDrugs(bc, "doxorubicin")

setwd("/media/feng/data1/骨肉瘤挖掘/肺转移/beyondcell")
general <- bc@ranks$general
write.csv(general,file = "beyoudcellgeneral.csv")
TC <- bc@ranks$TC
write.csv(TC,file = "beyoudcellTC.csv")
cluster <-  bc@ranks$bc_clusters_res.0.05
write.csv(cluster,file = "beyoudcellcluster.csv")
p1 <- bcClusters(bc, UMAP = "beyondcell", idents = "bc_clusters_res.0.01", pt.size = 1.5) +  scale_color_futurama()
ggsave("bcClusters.png",p1)
ggsave("bcClusters.pdf",p1)
ggsave("bcClusters.eps",p1)


MC <- readRDS("/media/feng/data1/骨肉瘤挖掘/肺转移/beyondcell/MC.rds")
MC$TC <- meta$TC
MC@active.ident <- as.factor(MC$TC)
class(MC@active.ident)[1]
table(MC@active.ident)
saveRDS(MC,file = "MCTC.rds")

#GSVA
library(dplyr)
library(GSVA)
library(GSEABase)
library(limma)
library(Seurat)
library(pheatmap)
library(stringr)
library(ggplot2)
library(psych)
library(Matrix)

genesets <- getGmt("/media/feng/data1/骨肉瘤挖掘/肺转移/Mac/GSVA/msigdb.v7.2.symbols.xls") 

Average_MY<- AverageExpression(MC, assays = NULL, features = NULL,return.seurat = FALSE, add.ident = NULL, slot = "data",use.scale = FALSE, use.counts = FALSE, verbose = TRUE)
Average_MY <- AverageExpression(MC)
gene<-Average_MY$RNA
gene<-as.matrix(gene)
GSVAresult <- gsva(gene, genesets, min.sz=10, max.sz=Inf, tau=1, method="gsva", kcdf="Poisson", mx.diff=TRUE, abs.ranking=FALSE, verbose=TRUE, parallel.sz=10)
pdf(file = "MY.pdf",width = 4.5,height = 13)
pheatmap(GSVAresult,color = colorRampPalette(c("navy","white","firebrick3"))(50),fontsize_row=5,border_color = "white")
dev.off()


p3 <- bc4Squares(bc, idents = "TC", lvl = "TC0", top = 5)

ggsave("bc4Squarescondition.png",p3)
ggsave("bc4Squarescondition.pdf",p3)
ggsave("bc4Squarescondition.eps",p3)
FindDrugs(bc, "")


bcSignatures(bc, UMAP = "beyondcell", signatures = list(values = "sig_7462"), pt.size = 1.5)


p <- bcSignatures(bc, UMAP = "beyondcell", signatures = list(values = "sig_7462"), pt.size = 1.5)
ggsave("sig_7462.png",p[[1]])
ggsave("sig_7462.pdf",p[[1]])
ggsave("sig_7462.eps",p[[1]])
p <- bcSignatures(bc, UMAP = "beyondcell", signatures = list(values = "sig_20652"), pt.size = 1.5)
ggsave("sig_20652.png",p[[1]])
ggsave("sig_20652.pdf",p[[1]])
ggsave("sig_20652.eps",p[[1]])
p <- bcSignatures(bc, UMAP = "beyondcell", signatures = list(values = "sig_2047"), pt.size = 1.5)
ggsave("sig_2047.png",p[[1]])
ggsave("sig_2047.pdf",p[[1]])
ggsave("sig_2047.eps",p[[1]])
p <- bcSignatures(bc, UMAP = "beyondcell", signatures = list(values = "sig_1257"), pt.size = 1.5)
ggsave("sig_1257.png",p[[1]])
ggsave("sig_1257.pdf",p[[1]])
ggsave("sig_1257.eps",p[[1]])



p <- bcSignatures(bc, UMAP = "beyondcell", signatures = list(values = "sig_8633"), pt.size = 1.5)
ggsave("sig_8633.png",p[[1]])
ggsave("sig_8633.pdf",p[[1]])
ggsave("sig_8633.eps",p[[1]])

p <- bcSignatures(bc, UMAP = "beyondcell", signatures = list(values = "sig_20574"), pt.size = 1.5)
ggsave("sig_20574.png",p[[1]])
ggsave("sig_20574.pdf",p[[1]])
ggsave("sig_20574.eps",p[[1]])

p <- bcSignatures(bc, UMAP = "beyondcell", signatures = list(values = "sig_2649"), pt.size = 1.5)
ggsave("sig_2649.png",p[[1]])
ggsave("sig_2649.pdf",p[[1]])
ggsave("sig_2649.eps",p[[1]])

p <- bcSignatures(bc, UMAP = "beyondcell", signatures = list(values = "sig_1339"), pt.size = 1.5)
ggsave("sig_1339.png",p[[1]])
ggsave("sig_1339.pdf",p[[1]])
ggsave("sig_1339.eps",p[[1]])

p <- bcSignatures(bc, UMAP = "beyondcell", signatures = list(values = "sig_7064"), pt.size = 1.5)
ggsave("sig_7064.png",p[[1]])
ggsave("sig_7064.pdf",p[[1]])
ggsave("sig_7064.eps",p[[1]])



library("patchwork").
bortezomib_IDs <- FindDrugs(bc, "bortezomib")$IDs
bortezomib_IDs1 <- c("sig_1339","sig_8633","sig_896","sig_20574","sig_2649")
bortezomib <- bcSignatures(bc, UMAP = "beyondcell", signatures = list(values = "bortezomib_IDs"), pt.size = 1.5)
wrap_plots(bortezomib, ncol = 3)
bcSignatures(bc, UMAP = "beyondcell", genes = list(values = "PSMA5"))

bcHistogram(bc, signatures = "sig_18868", idents = NULL)
p <- bcHistogram(bc, signatures = "sig_8633", idents = "TC")
ggsave("bcHistogramsig_8633.png",p[[1]])
ggsave("bcHistogramsig_8633.pdf",p[[1]])
ggsave("bcHistogramsig_8633.eps",p[[1]])
p <- bcHistogram(bc, signatures = "sig_1339", idents = "TC")
ggsave("bcHistogramsig_1339.png",p[[1]])
ggsave("bcHistogramsig_1339.pdf",p[[1]])
ggsave("bcHistogramsig_1339.eps",p[[1]])
p <- bcHistogram(bc, signatures = "sig_20574", idents = "TC")
ggsave("bcHistogramsig_20574.png",p[[1]])
ggsave("bcHistogramsig_20574.pdf",p[[1]])
ggsave("bcHistogramsig_20574.eps",p[[1]])
p <- bcHistogram(bc, signatures = "sig_2649", idents = "TC")
ggsave("bcHistogramsig_2649.png",p[[1]])
ggsave("bcHistogramsig_2649.pdf",p[[1]])
ggsave("bcHistogramsig_2649.eps",p[[1]])
p <- bcHistogram(bc, signatures = "sig_7064", idents = "TC")
ggsave("bcHistogramsig_7064.png",p[[1]])
ggsave("bcHistogramsig_7064.pdf",p[[1]])
ggsave("bcHistogramsig_7064.eps",p[[1]])
sig_8633_down <- gs@genelist$sig_8633$down
sig_2649_down <- gs@genelist$sig_2649$down
sig_8633_down_gene <- gene[rownames(gene)%in%sig_8633_down,]
sig_2649_down_gene <- gene[rownames(gene)%in%sig_2649_down,]
sig_8633_down_gene_t <- na.omit(t(scale(t(sig_8633_down_gene))))
sig_2649_down_gene_t <- na.omit(t(scale(t(sig_2649_down_gene))))
sig_8633_down_gene_t <- as.data.frame(sig_8633_down_gene_t)
sub_sig_8633_down_gene_t <- sig_8633_down_gene_t[sig_8633_down_gene_t$TC0>sig_8633_down_gene_t$TC1&sig_8633_down_gene_t$TC0>sig_8633_down_gene_t$TC2,]
write.csv(sub_sig_8633_down_gene_t,file = "sub_sig_8633_down_gene.csv")
p <- pheatmap(sub_sig_8633_down_gene_t,color = colorRampPalette(c("navy","white","firebrick3"))(50),fontsize_row=5,border_color = "white",cluster_cols = F)
ggsave("美伐他汀靶点.png",p,width = 3,height = 12)
ggsave("美伐他汀靶点.pdf",p,width = 3,height = 12)
ggsave("美伐他汀靶点.eps",p,width = 3,height = 12)

sig_2649_down_gene_t <- as.data.frame(sig_2649_down_gene_t)
sub_sig_2649_down_gene_t <- sig_2649_down_gene_t[sig_2649_down_gene_t$TC0>sig_2649_down_gene_t$TC1&sig_2649_down_gene_t$TC0>sig_2649_down_gene_t$TC2,]
write.csv(sub_sig_2649_down_gene_t,file = "sub_sig_2649_down_gene.csv")
p <- pheatmap(sub_sig_2649_down_gene_t,color = colorRampPalette(c("navy","white","firebrick3"))(50),fontsize_row=5,border_color = "white",cluster_cols = F)
ggsave("奥芬达唑靶点.png",p,width = 3,height = 12)
ggsave("奥芬达唑靶点.pdf",p,width = 3,height = 12)
ggsave("奥芬达唑靶点.eps",p,width = 3,height = 12)
ungene3 <- as.data.frame(rownames(sub_sig_8633_down_gene_t))
library("clusterProfiler")
library("org.Hs.eg.db")
library("enrichplot")
library("ggplot2")
genes=as.vector(rownames(sub_sig_8633_down_gene_t))
entrezIDs <- mget(genes, org.Hs.egSYMBOL2EG, ifnotfound=NA)    
entrezIDs <- as.character(entrezIDs)
out=cbind(ungene3$`rownames(sub_sig_8633_down_gene_t)`,entrezID=entrezIDs)
dir.create("/media/feng/data1/骨肉瘤挖掘/肺转移/beyondcell/美伐他丁/GODOWN/")
setwd("/media/feng/data1/骨肉瘤挖掘/肺转移/beyondcell/美伐他丁/GODOWN/")
write.table(out,file="id.txt",sep="\t",quote=F,row.names=F)   
rt=read.table("id.txt",sep="\t",header=T,check.names=F)           
rt=rt[is.na(rt[,"entrezID"])==F,]                                 
gene=rt$entrezID
kk <- enrichGO(gene = gene, OrgDb = org.Hs.eg.db, pvalueCutoff =0.05, qvalueCutoff = 0.05, ont="all", readable =T)
write.table(kk,file="GO.txt",sep="\t",quote=F,row.names = F)                 

p <-barplot(kk, drop = TRUE, showCategory =10,split="ONTOLOGY") + facet_grid(ONTOLOGY~., scale='free')
ggsave("GObarplot.png",p,width = 11,height = 8)
ggsave("GObarplot.eps",p,width = 11,height = 8)
ggsave("GObarplot.pdf",p,width = 11,height = 8)

p <- dotplot(kk,showCategory = 10,split="ONTOLOGY") + facet_grid(ONTOLOGY~., scale='free')
ggsave("GObubble.png",p,width = 11,height = 8)
ggsave("GObubble.eps",p,width = 11,height = 8)
ggsave("GObubble.pdf",p,width = 11,height = 8)

dir.create("/media/feng/data1/骨肉瘤挖掘/肺转移/beyondcell/美伐他丁/KEGGDOWN/")
setwd("/media/feng/data1/骨肉瘤挖掘/肺转移/beyondcell/美伐他丁/KEGGDOWN/")
entrezIDs <- mget(genes, org.Hs.egSYMBOL2EG, ifnotfound=NA)    
entrezIDs <- as.character(entrezIDs)
out=cbind(genes,entrezID=entrezIDs)
write.table(out,file="id.txt",sep="\t",quote=F,row.names=F)   
rt=read.table("id.txt",sep="\t",header=T,check.names=F)           
rt=rt[is.na(rt[,"entrezID"])==F,]                                 
gene=rt$entrezID
kk <- enrichKEGG(gene = gene, organism = "hsa", pvalueCutoff =0.05, qvalueCutoff =0.05)
write.table(kk,file="KEGG.txt",sep="\t",quote=F,row.names = F)                 
p <- barplot(kk, drop = TRUE, showCategory = 30)
ggsave("KEGGbarplot.png",p,width = 10,height = 8)
ggsave("KEGGbarplot.eps",p,width = 10,height = 8)
ggsave("KEGGbarplot.pdf",p,width = 10,height = 8)

p <- dotplot(kk, showCategory = 30)
ggsave("KEGGdotplot.png",p,width = 10,height = 8)
ggsave("KEGGdotplot.eps",p,width = 10,height = 8)
ggsave("KEGGdotplot.pdf",p,width = 10,height = 8)

ungene3 <- as.data.frame(rownames(sub_sig_2649_down_gene_t))
library("clusterProfiler")
library("org.Hs.eg.db")
library("enrichplot")
library("ggplot2")
genes=as.vector(rownames(sub_sig_2649_down_gene_t))
entrezIDs <- mget(genes, org.Hs.egSYMBOL2EG, ifnotfound=NA)    
entrezIDs <- as.character(entrezIDs)
out=cbind(ungene3$`rownames(sub_sig_2649_down_gene_t)`,entrezID=entrezIDs)
dir.create("/media/feng/data1/骨肉瘤挖掘/肺转移/beyondcell/奥芬达唑/GODOWN/")
setwd("/media/feng/data1/骨肉瘤挖掘/肺转移/beyondcell/奥芬达唑/GODOWN/")
write.table(out,file="id.txt",sep="\t",quote=F,row.names=F)   
rt=read.table("id.txt",sep="\t",header=T,check.names=F)           
rt=rt[is.na(rt[,"entrezID"])==F,]                                 
gene=rt$entrezID
kk <- enrichGO(gene = gene, OrgDb = org.Hs.eg.db, pvalueCutoff =0.05, qvalueCutoff = 0.05, ont="all", readable =T)
write.table(kk,file="GO.txt",sep="\t",quote=F,row.names = F)                 

p <-barplot(kk, drop = TRUE, showCategory =10,split="ONTOLOGY") + facet_grid(ONTOLOGY~., scale='free')
ggsave("GObarplot.png",p,width = 11,height = 8)
ggsave("GObarplot.eps",p,width = 11,height = 8)
ggsave("GObarplot.pdf",p,width = 11,height = 8)

p <- dotplot(kk,showCategory = 10,split="ONTOLOGY") + facet_grid(ONTOLOGY~., scale='free')
ggsave("GObubble.png",p,width = 11,height = 8)
ggsave("GObubble.eps",p,width = 11,height = 8)
ggsave("GObubble.pdf",p,width = 11,height = 8)

dir.create("/media/feng/data1/骨肉瘤挖掘/肺转移/beyondcell/奥芬达唑/KEGGDOWN/")
setwd("/media/feng/data1/骨肉瘤挖掘/肺转移/beyondcell/奥芬达唑/KEGGDOWN/")
entrezIDs <- mget(genes, org.Hs.egSYMBOL2EG, ifnotfound=NA)    
entrezIDs <- as.character(entrezIDs)
out=cbind(genes,entrezID=entrezIDs)
write.table(out,file="id.txt",sep="\t",quote=F,row.names=F)   
rt=read.table("id.txt",sep="\t",header=T,check.names=F)           
rt=rt[is.na(rt[,"entrezID"])==F,]                                 
gene=rt$entrezID
kk <- enrichKEGG(gene = gene, organism = "hsa", pvalueCutoff =0.05, qvalueCutoff =0.05)
write.table(kk,file="KEGG.txt",sep="\t",quote=F,row.names = F)                 
p <- barplot(kk, drop = TRUE, showCategory = 30)
ggsave("KEGGbarplot.png",p,width = 10,height = 8)
ggsave("KEGGbarplot.eps",p,width = 10,height = 8)
ggsave("KEGGbarplot.pdf",p,width = 10,height = 8)

p <- dotplot(kk, showCategory = 30)
ggsave("KEGGdotplot.png",p,width = 10,height = 8)
ggsave("KEGGdotplot.eps",p,width = 10,height = 8)
ggsave("KEGGdotplot.pdf",p,width = 10,height = 8)

pheatmap(sub_sig_2649_down_gene_t,color = colorRampPalette(c("navy","white","firebrick3"))(50),fontsize_row=5,border_color = "white",cluster_cols = F)

bcCellCycle(bc, signatures = "sig_2649")


##CSOmap##
library(Seurat)
library(CSOmapR)
library(CSOmapR.demo)
invisible(TPM)
invisible(LR)
invisible(labelData)
setwd("/media/feng/data1/骨肉瘤挖掘/肺转移/Tcell/ALLCSOmap")
test <- OS_lung
new.cluster.ids <- c("T cells","Malignant cells","Other cells","Other cells","Other cells","Endothelial",
                     "Other cells","Other cells","Macrophages","Other cells","Other cells","Other cells",
                     "Other cells","Other cells","Other cells","Other cells","Other cells")
names(new.cluster.ids) <- levels(test)
test <- RenameIdents(test, new.cluster.ids)
test$Celltype <- test@active.ident
Idents(test) <- "Donor"
OScon <- subset(test,idents = c("BC10"))
OScho <- subset(test,idents = c("BC17"))
saveRDS(OScon,file= "/media/feng/data1/骨肉瘤挖掘/肺转移/Tcell/allOScon.rds")
saveRDS(OScho,file= "/media/feng/data1/骨肉瘤挖掘/肺转移/Tcell/allOScho.rds")
tmp <- sample(colnames(OScon),1564) %>% sort()
test5257 <- OScon[,tmp]
count5257 <- as.data.frame(test5257@assays$RNA@counts)

genelength <- read.table("/media/feng/data1/骨肉瘤挖掘/肺转移/gene.length.txt",header = F,sep = "\t")
colnames(genelength) <- c("ID","length")
count5257 <- count5257[rownames(count5257)%in%genelength$ID,]
count5257$ID <- rownames(count5257)
count5257 <- merge(count5257,genelength,by = "ID")
rownames(count5257) <- count5257[,1]
count5257 <- count5257[,-1]
count5257 <- count5257/count5257$length
shendu <- apply(count5257,2,sum)
TMP5257 <- count5257/shendu
TMP5257 <- TMP5257*10^6
write.table(TMP5257,file = "CON_CSOmap1564TMP.txt")
Celltype <- as.data.frame(test5257$Celltype)
Celltype$cells <- rownames(Celltype)
colnames(Celltype) <- c("labels","cells")
Celltype$Celltype <- Celltype$labels
Celltype <- Celltype[,-1]
TMP5257 <- TMP5257[,-1565]
write.table(Celltype,file = "CON_CSOmap1564celltype.txt")


affinityMat = getAffinityMat(TMP5257, LR, verbose = T)

coords_res = runExactTSNE_R(
  X = affinityMat,
  no_dims = 3,
  max_iter = 1000,
  verbose = T
)

saveRDS(coords_res,"/media/feng/data1/骨肉瘤挖掘/肺转移/Tcell/1564CON_CSOmap_coords_res.rds")
coords = coords_res$Y
rownames(coords) <- colnames(TMP5257)
colnames(coords) <- c('x', 'y', 'z')

require(dplyr)
coords_tbl = bind_cols(cellName = rownames(coords), as.data.frame(coords))

join_vec = setNames(colnames(Celltype)[1], nm = colnames(coords_tbl)[1])
cellinfo_tbl = left_join(coords_tbl, Celltype, by = join_vec)

density_obj = getDensity3D(cellinfo_tbl$x, cellinfo_tbl$y, cellinfo_tbl$z)
cellinfo_tbl = cellinfo_tbl %>% mutate(density = density_obj)

p_3Ddensity_Celltype= plot3D(cellinfo_tbl, color_by = "Celltype", title = "3D density")

p_3Ddensity_density = plot3D(cellinfo_tbl, color_by = "density", title = "3D density")
saveRDS(cellinfo_tbl,"/media/feng/data1/骨肉瘤挖掘/肺转移/Tcell/1564CON_CSOmap_cellinfo_tbl.rds")


xy <- coords[,1:2]
xy <- as.data.frame(xy)
colnames(xy) <- c("x","y")
xy$cells <- rownames(xy)
xy <- merge(xy ,Celltype ,by = "cells")
cols <- c("#80C6B1","#F39E7F","#A1AFD0","#E09EC4","#B1D475")
p <- ggplot(xy, aes(x = x, y = y, colour = Celltype, size=0.5)) +geom_point()+theme_bw()+theme(panel.grid=element_blank(),panel.border=element_blank(),axis.line=element_line(size=0.5,colour="black"))+ scale_color_manual(name="",labels = c("T cells", "Malignant cells", "Other cells", "Endothelial","Macrophages" ), 
                                                                                                                                                                                                                            values =c("#80C6B1","#F39E7F","#A1AFD0","#E09EC4","#B1D475"))

ggsave("CONCSOmap_xy.png",p)
ggsave("CONCSOmap_xy.eps",p)
ggsave("CONCSOmap_xy.pdf",p)

distance <- as.data.frame(coords)
distance$cells <- rownames(distance)
distance <- merge(distance,Celltype,by = "cells")
distance$distance <- sqrt((distance$x)^2+(distance$y)^2+(distance$z)^2)
write.csv(distance,file = "1564CONdistance.csv")
signif_results = getSignificance(coords, labels = cellinfo_tbl$Celltype, verbose = T)
contribution_list = getContribution(TMP5257, LR, signif_results$detailed_connections)
saveRDS(signif_results,file = "1564CON_CSOmap_signif_results.rds")
saveRDS(contribution_list,file = "1564CON_CSOmap_contribution_list.rds")
a <- as.data.frame(contribution_list$`Malignant cells---Malignant cells`)
write.csv(a,file = "1564CON_CSOmap_genecontribution.csv")
CD63up <- TMP5257
CD63up <- t(CD63up)
CD63up <- as.data.frame(CD63up)
CD63up$cells <- rownames(CD63up)
allup <- merge(CD63up,Celltype,by = "cells")
CD63up <- allup[allup$Celltype == "Malignant cells",]
CD63up <- as.data.frame(CD63up)
sum(CD63up$CD63)/1011
CD63up$CD63 <- 1032161
other <- allup[!allup$Celltype == "Malignant cells",]
CD63up <- rbind(CD63up,other)
CD63up <- CD63up[,-24192]
rownames(CD63up) <- CD63up$cells
CD63up <- CD63up[,-1]
CD63up <- t(CD63up)

affinityMat = getAffinityMat(CD63up, LR, verbose = T)

coords_res = runExactTSNE_R(
  X = affinityMat,
  no_dims = 3,
  max_iter = 1000,
  verbose = T
)

saveRDS(coords_res,"/media/feng/data1/骨肉瘤挖掘/肺转移/Tcell/CD63up_1564CON_CSOmap_coords_res.rds")
coords = coords_res$Y
rownames(coords) <- colnames(CD63up)
colnames(coords) <- c('x', 'y', 'z')

require(dplyr)
coords_tbl = bind_cols(cellName = rownames(coords), as.data.frame(coords))

join_vec = setNames(colnames(Celltype)[1], nm = colnames(coords_tbl)[1])
cellinfo_tbl = left_join(coords_tbl, Celltype, by = join_vec)

density_obj = getDensity3D(cellinfo_tbl$x, cellinfo_tbl$y, cellinfo_tbl$z)
cellinfo_tbl = cellinfo_tbl %>% mutate(density = density_obj)

p_3Ddensity_Celltype_CD63up= plot3D(cellinfo_tbl, color_by = "Celltype", title = "3D density")

p_3Ddensity_density_CD63up = plot3D(cellinfo_tbl, color_by = "density", title = "3D density")
saveRDS(cellinfo_tbl,"/media/feng/data1/骨肉瘤挖掘/肺转移/Tcell/CD63up_1564CON_CSOmap_cellinfo_tbl.rds")

xy <- coords[,1:2]
xy <- as.data.frame(xy)
colnames(xy) <- c("x","y")
xy$cells <- rownames(xy)
xy <- merge(xy ,Celltype ,by = "cells")

cols <- c("#80C6B1","#F39E7F","#A1AFD0","#E09EC4","#B1D475")
p <- ggplot(xy, aes(x = x, y = y, colour = Celltype, size=0.5)) +geom_point()+theme_bw()+theme(panel.grid=element_blank(),panel.border=element_blank(),axis.line=element_line(size=0.5,colour="black"))+ scale_color_manual(name="",labels = c("T cells", "Malignant cells", "Other cells", "Endothelial","Macrophages" ), 
                                                                                                                                                                                                                            values =c("#80C6B1","#F39E7F","#A1AFD0","#E09EC4","#B1D475"))


ggsave("CD63up_CONCSOmap_xy.png",p)
ggsave("CD63up_CONCSOmap_xy.eps",p)
ggsave("CD63up_CONCSOmap_xy.pdf",p)

distance <- as.data.frame(coords)
distance$cells <- rownames(distance)
distance <- merge(distance,Celltype,by = "cells")
distance$distance <- sqrt((distance$x)^2+(distance$y)^2+(distance$z)^2)
write.csv(distance,file = "CD63up_CONdistance.csv")

CD63down <- TMP5257
CD63down <- t(CD63down)
CD63down <- as.data.frame(CD63down)
CD63down$cells <- rownames(CD63down)
allup <- merge(CD63down,Celltype,by = "cells")
CD63down <- allup[allup$Celltype == "Malignant cells",]
CD63down <- as.data.frame(CD63down)
sum(CD63down$CD63)/1011
CD63down$CD63 <- 0
other <- allup[!allup$Celltype == "Malignant cells",]
CD63down <- rbind(CD63down,other)
CD63down <- CD63down[,-24192]
rownames(CD63down) <- CD63down$cells
CD63down <- CD63down[,-1]
CD63down <- t(CD63down)

affinityMat = getAffinityMat(CD63down, LR, verbose = T)

coords_res = runExactTSNE_R(
  X = affinityMat,
  no_dims = 3,
  max_iter = 1000,
  verbose = T
)

saveRDS(coords_res,"/media/feng/data1/骨肉瘤挖掘/肺转移/Tcell/CD63down_1564CON_CSOmap_coords_res.rds")
coords = coords_res$Y
rownames(coords) <- colnames(CD63down)
colnames(coords) <- c('x', 'y', 'z')

require(dplyr)
coords_tbl = bind_cols(cellName = rownames(coords), as.data.frame(coords))

join_vec = setNames(colnames(Celltype)[1], nm = colnames(coords_tbl)[1])
cellinfo_tbl = left_join(coords_tbl, Celltype, by = join_vec)

density_obj = getDensity3D(cellinfo_tbl$x, cellinfo_tbl$y, cellinfo_tbl$z)
cellinfo_tbl = cellinfo_tbl %>% mutate(density = density_obj)

p_3Ddensity_Celltype_CD63down = plot3D(cellinfo_tbl, color_by = "Celltype", title = "3D density")

p_3Ddensity_density_CD63down = plot3D(cellinfo_tbl, color_by = "density", title = "3D density")
saveRDS(cellinfo_tbl,"/media/feng/data1/骨肉瘤挖掘/肺转移/Tcell/CD63down_1564CON_CSOmap_cellinfo_tbl.rds")


##################Z等于0二维平面展示##################

xy <- coords[,1:2]
xy <- as.data.frame(xy)
colnames(xy) <- c("x","y")
xy$cells <- rownames(xy)
xy <- merge(xy ,Celltype ,by = "cells")

cols <- c("#80C6B1","#F39E7F","#A1AFD0","#E09EC4","#B1D475")

p <- ggplot(xy, aes(x = x, y = y, colour = Celltype, size=0.5)) +geom_point()+theme_bw()+theme(panel.grid=element_blank(),panel.border=element_blank(),axis.line=element_line(size=0.5,colour="black"))+ scale_color_manual(name="",labels = c("T cells", "Malignant cells", "Other cells", "Endothelial","Macrophages" ), 
                                                                                                                                                                                                                            values =c("#80C6B1","#F39E7F","#A1AFD0","#E09EC4","#B1D475"))
ggsave("CD63down_CONCSOmap_xy.png",p)
ggsave("CD63down_CONCSOmap_xy.eps",p)
ggsave("CD63down_CONCSOmap_xy.pdf",p)

distance <- as.data.frame(coords)
distance$cells <- rownames(distance)
distance <- merge(distance,Celltype,by = "cells")
distance$distance <- sqrt((distance$x)^2+(distance$y)^2+(distance$z)^2)
write.csv(distance,file = "CD63down_CONdistance.csv")
signif_results = getSignificance(coords, labels = cellinfo_tbl$Celltype, verbose = T)
contribution_list = getContribution(TMP5257, LR, signif_results$detailed_connections)

allCD63down <- TMP5257
allCD63down <- t(allCD63down)
allCD63down <- as.data.frame(allCD63down)
rownames(allCD63down) <- gsub("[.]","-",rownames(allCD63down))
allCD63down$cells <- rownames(allCD63down)
allCD63down <- merge(allCD63down,Celltype,by = "cells")

sum(allCD63down$CD63)/1564
allCD63down$CD63 <- 0
allCD63down <- allCD63down[,-24192]
rownames(allCD63down) <- allCD63down$cells
allCD63down <- allCD63down[,-1]
allCD63down <- t(allCD63down)

affinityMat = getAffinityMat(allCD63down, LR, verbose = T)

coords_res = runExactTSNE_R(
  X = affinityMat,
  no_dims = 3,
  max_iter = 1000,
  verbose = T
)

saveRDS(coords_res,"/media/feng/data1/骨肉瘤挖掘/肺转移/Tcell/allCD63down_1564CON_CSOmap_coords_res.rds")
coords = coords_res$Y
rownames(coords) <- colnames(allCD63down)
colnames(coords) <- c('x', 'y', 'z')

require(dplyr)
coords_tbl = bind_cols(cellName = rownames(coords), as.data.frame(coords))

join_vec = setNames(colnames(Celltype)[1], nm = colnames(coords_tbl)[1])
cellinfo_tbl = left_join(coords_tbl, Celltype, by = join_vec)

density_obj = getDensity3D(cellinfo_tbl$x, cellinfo_tbl$y, cellinfo_tbl$z)
cellinfo_tbl = cellinfo_tbl %>% mutate(density = density_obj)

p_3Ddensity_Celltype_allCD63down = plot3D(cellinfo_tbl, color_by = "Celltype", title = "3D density")

p_3Ddensity_density_allCD63down = plot3D(cellinfo_tbl, color_by = "density", title = "3D density")
saveRDS(cellinfo_tbl,"/media/feng/data1/骨肉瘤挖掘/肺转移/Tcell/allCD63down_1564CON_CSOmap_cellinfo_tbl.rds")

xy <- coords[,1:2]
xy <- as.data.frame(xy)
colnames(xy) <- c("x","y")
xy$cells <- rownames(xy)
xy <- merge(xy ,Celltype ,by = "cells")

cols <- c("#80C6B1","#F39E7F","#A1AFD0","#E09EC4","#B1D475")

p <- ggplot(xy, aes(x = x, y = y, colour = Celltype, size=0.5)) +geom_point()+theme_bw()+theme(panel.grid=element_blank(),panel.border=element_blank(),axis.line=element_line(size=0.5,colour="black"))+ scale_color_manual(name="",labels = c("T cells", "Malignant cells", "Other cells", "Endothelial","Macrophages" ), 
                                                                                                                                                                                                                            values =c("#80C6B1","#F39E7F","#A1AFD0","#E09EC4","#B1D475"))
ggsave("CD63down_CONCSOmap_xy.png",p)
ggsave("CD63down_CONCSOmap_xy.eps",p)
ggsave("CD63down_CONCSOmap_xy.pdf",p)

distance <- as.data.frame(coords)
distance$cells <- rownames(distance)
distance <- merge(distance,Celltype,by = "cells")
distance$distance <- sqrt((distance$x)^2+(distance$y)^2+(distance$z)^2)
write.csv(distance,file = "allCD63down_CONdistance.csv")
#Get significance
signif_results = getSignificance(coords, labels = cellinfo_tbl$Celltype, verbose = T)
contribution_list = getContribution(TMP5257, LR, signif_results$detailed_connections)

count5257 <- as.data.frame(OScho@assays$RNA@counts)

genelength <- read.table("/media/feng/data1/骨肉瘤挖掘/肺转移/gene.length.txt",header = F,sep = "\t")
colnames(genelength) <- c("ID","length")
count5257 <- count5257[rownames(count5257)%in%genelength$ID,]
count5257$ID <- rownames(count5257)
count5257 <- merge(count5257,genelength,by = "ID")
rownames(count5257) <- count5257[,1]
count5257 <- count5257[,-1]
count5257 <- count5257/count5257$length
shendu <- apply(count5257,2,sum)
TMP5257 <- count5257/shendu
TMP5257 <- TMP5257*10^6
write.table(TMP5257,file = "CHO_CSOmap2659TMP.txt")
Celltype <- as.data.frame(OScho$Celltype)
Celltype$cells <- rownames(Celltype)
colnames(Celltype) <- c("labels","cells")
Celltype$Celltype <- Celltype$labels
Celltype <- Celltype[,-1]
TMP5257 <- TMP5257[,-2660]
write.table(Celltype,file = "CHO_CSOmap2659celltype.txt")

affinityMat = getAffinityMat(TMP5257, LR, verbose = T)

coords_res = runExactTSNE_R(
  X = affinityMat,
  no_dims = 3,
  max_iter = 1000,
  verbose = T
)

saveRDS(coords_res,"/media/feng/data1/骨肉瘤挖掘/肺转移/Tcell/2659CHO_CHOmap_coords_res.rds")
coords = coords_res$Y
rownames(coords) <- colnames(TMP5257)
colnames(coords) <- c('x', 'y', 'z')

require(dplyr)
coords_tbl = bind_cols(cellName = rownames(coords), as.data.frame(coords))

join_vec = setNames(colnames(Celltype)[1], nm = colnames(coords_tbl)[1])
cellinfo_tbl = left_join(coords_tbl, Celltype, by = join_vec)

density_obj = getDensity3D(cellinfo_tbl$x, cellinfo_tbl$y, cellinfo_tbl$z)
cellinfo_tbl = cellinfo_tbl %>% mutate(density = density_obj)

p_3Ddensity_Celltype = plot3D(cellinfo_tbl, color_by = "Celltype", title = "3D density")

p_3Ddensity_density = plot3D(cellinfo_tbl, color_by = "density", title = "3D density")
saveRDS(cellinfo_tbl,"/media/feng/data1/骨肉瘤挖掘/肺转移/Tcell/2659CHO_CSOmap_cellinfo_tbl.rds")

xy <- coords[,1:2]
xy <- as.data.frame(xy)
colnames(xy) <- c("x","y")
xy$cells <- rownames(xy)
xy <- merge(xy ,Celltype ,by = "cells")

cols <- c("#80C6B1","#F39E7F","#A1AFD0","#E09EC4","#B1D475")

p <- ggplot(xy, aes(x = x, y = y, colour = Celltype, size=0.5)) +geom_point()+theme_bw()+theme(panel.grid=element_blank(),panel.border=element_blank(),axis.line=element_line(size=0.5,colour="black"))+ scale_color_manual(name="",labels = c("T cells", "Malignant cells", "Other cells", "Endothelial","Macrophages" ), 
                                                                                                                                                                                                                            values =c("#80C6B1","#F39E7F","#A1AFD0","#E09EC4","#B1D475"))
ggsave("CHOCSOmap_xy.png",p)
ggsave("CHOCSOmap_xy.eps",p)
ggsave("CHOCSOmap_xy.pdf",p)

distance <- as.data.frame(coords)
distance$cells <- rownames(distance)
distance <- merge(distance,Celltype,by = "cells")
distance$distance <- sqrt((distance$x)^2+(distance$y)^2+(distance$z)^2)
write.csv(distance,file = "2659CHOdistance.csv")
signif_results = getSignificance(coords, labels = cellinfo_tbl$Celltype, verbose = T)
contribution_list = getContribution(TMP5257, LR, signif_results$detailed_connections)
saveRDS(signif_results,file = "2659CHO_CSOmap_signif_results.rds")
saveRDS(contribution_list,file = "2659CHO_CSOmap_contribution_list.rds")
a <- as.data.frame(contribution_list$`Malignant cells---Malignant cells`)
write.csv(a,file = "1564CON_CSOmap_genecontribution.csv")


###CNV###
library(Seurat)
library(infercnv)
library(tidyverse)
setwd("/media/feng/data1/骨肉瘤挖掘/肺转移/CNV/")
OS_lung <- readRDS("/media/feng/data1/骨肉瘤挖掘/肺转移/renamesOS_lung.rds")
OS <- subset(OS_lung,idents = c("Chondroblasts"))
ref <- subset(OS_lung,idents = c("NK/T cells","Monocytes","Mast cells","Macrophages","DCs","B cells","Plasma cells"))
tset1 <-merge(OS,y= ref)

tset1[["celltype"]]<-tset1@active.ident

cellAnnota <- subset(tset1@meta.data, select='celltype')
exprMatrix <- as.matrix(GetAssayData(tset1, slot='counts'))
infercnv_obj = CreateInfercnvObject(delim = '\t',
                                    raw_counts_matrix = exprMatrix,
                                    annotations_file = cellAnnota,
                                    gene_order_file = '/media/feng/data1/骨肉瘤挖掘/肺转移/CNV/geneLocate.txt',
                                    ref_group_names = c("NK/T cells","Monocytes","Mast cells","Macrophages","DCs","B cells","Plasma cells"))
dir.create("/media/feng/data1/骨肉瘤挖掘/肺转移/CNV/malignant")
infercnv_obj = infercnv::run(infercnv_obj,
                             cutoff=0.1, 
                             out_dir='/media/feng/data1/骨肉瘤挖掘/肺转移/CNV/malignant', 
                             cluster_by_groups=TRUE, 
                             denoise=TRUE,
                             HMM=TRUE)
