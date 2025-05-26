rm(list = ls())
gc()
options(stringsAsFactors = F)
set.seed(123)
library(Seurat)
library(SeuratObject)
library(ggplot2)
library(tidyverse)
library(patchwork)
library(RColorBrewer)
library(harmony)
library(DoubletFinder)
library(dplyr)
library(cowplot)
library(clusterProfiler)
library(org.Hs.eg.db)
setwd("/Lustre03/data/luotingting/pro/WAT/sn31.2205/result333/")
assays <- dir("/Lustre03/data/WeiRanlei/shell/count_mtx/")
dir <- paste0("/Lustre03/data/WeiRanlei/shell/count_mtx/",assays,"/",assays,"_matrix_10X")
names(dir) = c(assays)
snRNAlist <- list()

for(i in 1:length(dir)){
  counts <- Read10X(data.dir = dir[i])
  snRNAlist[[i]] <- CreateSeuratObject(counts, min.cells=3, min.features = 200)
  snRNAlist[[i]]$sample<-assays[[i]]
  #double
  dim.num = 1:30
  pc.num = 1:30
  snRNAlist[[i]] <- SCTransform(snRNAlist[[i]])
  snRNAlist[[i]] <- RunPCA(snRNAlist[[i]])
  snRNAlist[[i]] <- RunUMAP(snRNAlist[[i]], dims = dim.num)
  sweep.res.list <- paramSweep_v3(snRNAlist[[i]], PCs = pc.num, sct = T)
  sweep.stats <- summarizeSweep(sweep.res.list, GT = FALSE)  
  bcmvn <- find.pK(sweep.stats)
  pK_bcmvn <- bcmvn$pK[which.max(bcmvn$BCmetric)] %>% as.character() %>% as.numeric()
  DoubletRate = 0.04
  nExp_poi <- round(DoubletRate*nrow(snRNAlist[[i]]@meta.data))
  snRNAlist[[i]] <- doubletFinder_v3(snRNAlist[[i]], PCs = pc.num, pN = 0.25, pK = pK_bcmvn, 
                                     nExp = nExp_poi, reuse.pANN = F, sct = T)
  colnames(snRNAlist[[i]]@meta.data)<-c("orig.ident", "nCount_RNA","nFeature_RNA","sample","nCount_SCT","nFeature_SCT","pANN","DF.classifications")
}
snRNA <- merge(snRNAlist[[1]], snRNAlist[2:length(snRNAlist)])
saveRDS(snRNA, "rds/snRNA_merge_dbfinder.rds")

snRNA <-subset(snRNA,subset= DF.classifications %in% "Singlet")
saveRDS(snRNA, "rds/snRNA_merge_delet.dbfinder.rds")

MT.genes <- c("ND1","ND2","COX1","COX2","ATP8","ATP6","COX32","ND3","ND4L","ND4","ND5","ND6","CYTB")
MT.genes <- CaseMatch(MT.genes, rownames(snRNA))
snRNA[["percent.mt"]] <- PercentageFeatureSet(snRNA, features=MT.genes)
VlnPlot(snRNA, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), ncol = 3)
ggsave(filename = "fliter_vlnplot.pdf",height = 10,width = 18)

table(snRNA$sample)
minGene=500
maxGene=5000
maxUMI=15000
minUMI=1000
pctMT=15
snRNA <- subset(snRNA, subset = nCount_RNA < maxUMI  &  nCount_RNA > minUMI & 
                  nFeature_RNA > minGene & nFeature_RNA < maxGene & 
                  percent.mt < pctMT )
snRNA

MT.genes <- c("ND1","ND2","COX1","COX2","ATP8","ATP6","COX32","ND3","ND4L","ND4","ND5","ND6","CYTB")
MT.genes <- CaseMatch(MT.genes, rownames(snRNA))
keep_gene<-setdiff(rownames(snRNA),MT.genes)
snRNA<-snRNA[keep_gene,]

gene2type<- read.table("/Lustre03/data/luotingting/soft/ensembl/original/SW_gene2types.txt", header = F)
protein<-gene2type[gene2type$V2 == "protein_coding",]
protein_gene<-protein$V1
protein_keep_gene<-intersect(rownames(snRNA),protein_gene)
snRNA<-snRNA[protein_keep_gene,]
saveRDS(snRNA,"rds/snRNA_merge_dbfinder_fliter.mt.lnc.gene.rds")

snRNA <- SCTransform(snRNA, assay = "RNA",vars.to.regress ="nCount_RNA",verbose = FALSE)
snRNA <- RunPCA(snRNA,verbose = FALSE)
snRNA <- RunHarmony(object = snRNA, group.by.vars="sample", assay.use="SCT",plot_convergence = TRUE)
saveRDS(snRNA,"rds/snRNA_merge_dbfinder_fliter.mt.lnc.gene.sct.rds")

pc.num=1:50
snRNA <- snRNA %>%
  RunUMAP(reduction = "harmony", dims = pc.num, min.dist=0.01) %>%
  FindNeighbors(reduction = "harmony", dims = pc.num) %>%
  FindClusters(verbose = FALSE) %>%
  identity()
saveRDS(snRNA, "rds/snRNA_merge_dbfinder_fliter.mt.lnc.gene.sct.cluster.rds")

dir.create("dimplot")
p2<-DimPlot(snRNA, reduction = "umap",label=T)
path2<-c("dimplot/cluster_umap.pdf")
ggsave(filename =path2, width = 8, height = 6)

dir.create("dotplot")
p1<-DotPlot(object =snRNA,group.by = "seurat_clusters",features =c("ADIPOR2","ADIPOQ","PLIN4","LIPE","FBN1","DCN","COL6A3","LPAR1","PECAM1","FLT1","EMCN","NOS3","LDB2","PTPRB","MYH11","MYO1B","TPM1","SYNPO2","ACTA2","HEYL","MCAM","KRT14","KRT7","ELF3","ADAM28","EZR","TM4SF1","C4A","FN1","VIM","CD163","RBPJ","ITGAM","BANK1","DOCK2","CD3D","CD3E","CD5","CLIP1","ADGRE5","ITGAL","ADGRE1"))+ RotatedAxis()
path1<-c("dotplot/maker_cluster.pdf")
ggsave(filename = path1,height = 10,width = 15)

dir.create("table")
cluster_table<-table(snRNA$sample,snRNA$seurat_clusters)
path3<-c("table/cluster_sample_table.csv")
write.csv(cluster_table,path3)
cluster_table_prop<-prop.table(cluster_table,1)
path4<-c("table/cluster_sample_table_prop.csv")
write.csv(cluster_table_prop,path4)

snRNA <-FindClusters(snRNA,resolution = 0.3,verbose = FALSE)
snRNA$cell_type<-NA
for (i in c("1","17","12","20")){
  snRNA@meta.data[which(snRNA@meta.data$seurat_clusters==i),"cell_type"]<-"Adipocytes"
}
for (i in c("2","13","23")){
  snRNA@meta.data[which(snRNA@meta.data$seurat_clusters==i),"cell_type"]<-"FAPs"
}
for (i in c("4","8","11")){
  snRNA@meta.data[which(snRNA@meta.data$seurat_clusters==i),"cell_type"]<-"Endothelial"
}
for (i in c("6")){
  snRNA@meta.data[which(snRNA@meta.data$seurat_clusters==i),"cell_type"]<-"Smooth_muscle"
}
for (i in c("3","18")){
  snRNA@meta.data[which(snRNA@meta.data$seurat_clusters==i),"cell_type"]<-"Mesothelium"
}
for (i in c("14","15","22","25","26")){
  snRNA@meta.data[which(snRNA@meta.data$seurat_clusters==i),"cell_type"]<-"Fibroblast"
}
for (i in c("5","21","24")){
  snRNA@meta.data[which(snRNA@meta.data$seurat_clusters==i),"cell_type"]<-"Macrophage"
}
for (i in c("16")){
  snRNA@meta.data[which(snRNA@meta.data$seurat_clusters==i),"cell_type"]<-"B_cell"
}
for (i in c("7")){
  snRNA@meta.data[which(snRNA@meta.data$seurat_clusters==i),"cell_type"]<-"T_cell"
}
for (i in c("9")){
  snRNA@meta.data[which(snRNA@meta.data$seurat_clusters==i),"cell_type"]<-"Dendritic"
}
for (i in c("0","10","27","19")){
  snRNA@meta.data[which(snRNA@meta.data$seurat_clusters==i),"cell_type"]<-"Pre_adipocytes"
}

keep_clusters <- c(18, 8, 0, 3, 19, 7, 1, 5, 4, 11, 14, 6, 21, 20, 16, 9, 24, 15, 2, 13, 12, 22, 17)
snRNA <- subset(snRNA, idents = keep_clusters, ident.use = "SCT_snn_res.0.3")

Idents(snRNA)="cell_type"
DefaultAssay(snRNA) <- "RNA"  
snRNA <- NormalizeData(snRNA)
all.genes <- rownames(snRNA)
snRNA <- ScaleData(snRNA, features = all.genes)
snRNA.findallmarkers <- FindAllMarkers(snRNA,only.pos = TRUE, logfc.threshold = 0.25)
write.csv(snRNA.findallmarkers,"table/snRNA.findallmarkers.csv")




