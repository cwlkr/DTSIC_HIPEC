## install dwls algorithm ##
library(devtools)
library(data.table)
library(tidyverse)
install_bitbucket("yuanlab/dwls")
install_github('RGLab/MAST')

ov.rnaseq <- fread('aligned_bulk/readcounts_realigment_starH38104_tpm.txt',header=T) %>% rename(external_gene_id='name')

ov.rnaseq %>% filter(!(external_gene_id  %like% '\\.\\d')) %>%
  filter(!(external_gene_id  %like% 'LIN.*')) %>%
  filter(!(external_gene_id  %like% '.*-.*')) -> target

genelist = target$external_gene_id
target %>% select(-gene_id, -external_gene_id) -> target
target %>% as.matrix() -> target
target %>% dim
rownames(target) <- genelist
nCores <- 16

pv <- 0.01
lfc <- 0.5
path <- 'dwls_results.csv'

diff.frame <- readRDS('seurat_output/markers.rds')

# extract the marker genes <- all genes that are significant in a vector
diff.frame %>% filter(p_val_adj < pv & avg_logFC > lfc) %>% group_by(cluster) %>% arrange(cluster, desc(avg_logFC)) %>% pull(gene)-> signif_genes

markerGenes <- sapply(signif_genes %>% unique,function(x) strsplit(x,"\\--")[[1]][1]) %>% unname
summary(diff.frame$p_val_adj)

#only take marker genes that are present in you target data set
markerGenes <- intersect(markerGenes, rownames(target))


new.cluster.ids <- c('Immune_cells_1_makrophages',
                     'Ovarian_cells',
                     'Tumor_cells_1',
                     'Tumor_cells_2',
                     'Tumor_cells_3',
                     'Fibroblasts_1',
                     'Tumor_cells_4',
                     'Tumor_cells_5',
                     'Smooth_muscle_1',
                     'Fibroblasts_2',
                     'Fibroblasts_3',
                     'Fibroblasts_4',
                     'Immune_cells_2_cd8t_cell',
                     'Endothelial_cells',
                     'Mesothelial_cells',
                     'Tumor_cells_6',
                     'Immune_cells_3_b_cells',
                     'Immune_cells_4_makrophages',
                     'Smooth_muscle_2',
                     'Tumor_cells_7',
                     'Fallopian_tube_epithelial_1',
                     'Immune_cells_5_b_cells',
                     'Tumor_cells_8',
                     'Fibroblasts_5',
                     'Fibroblasts_6')

pbmc.s <- readRDS('seurat_output/seurat_object.rds')
library(Seurat)
DimPlot(pbmc.s, reduction = "umap")

Idents(pbmc.s) <- pbmc.s$seurat_clusters
sc = pbmc.s

sc@meta.data$celltype <- NA
for ( i in 0:max(as.numeric(as.character(sc@meta.data$seurat_clusters)))){
  sc@meta.data$celltype[sc@meta.data$seurat_clusters == i] <- new.cluster.ids[i+1]
}
names(new.cluster.ids) <- levels(sc)
sc <- RenameIdents(sc, new.cluster.ids)
refCellTypes <- unique(sc@meta.data$celltype)

## construct reference profiles by averaging all cells belonging to the same cell type ##
refProfiles <- matrix(ncol=length(refCellTypes),nrow=nrow(sc))
scData <- GetAssayData(sc,slot = "counts")
scData <- apply(scData,2,function(x) (x/sum(x))*1000000)

for ( i in 1:ncol(refProfiles)){
  refProfiles[,i] <- apply(scData[,rownames(sc@meta.data)[sc@meta.data$celltype == refCellTypes[i]]],1,function(x) mean(x))

}
colnames(refProfiles) <- refCellTypes
rownames(refProfiles) <- sapply( rownames(sc),function(x) strsplit(x,"\\--")[[1]][1]) %>% unname


## target is the matrix with the real tumor samples ##
## refProfiles is the signature matrix ##
## markerGenes is the vector containing all the marker genes ##

## normalize reference profiles and target dataset ##
rt <- as.data.frame(apply(target,2,function(x) x/sum(x))[markerGenes,])
sig <- as.matrix(apply(refProfiles,2,function(x) x/sum(x))[markerGenes,])
require(foreach)
## run paralellized dwls on nCores number of cores ##
require(parallel)
cl <- makeCluster(nCores)
require(doSNOW)
registerDoSNOW(cl)
pb <- txtProgressBar(max = ncol(rt), style = 3)
progressBar <- function(n) setTxtProgressBar(pb, n)
opts <- list(progress=progressBar)

results <- foreach(x = c(1:ncol(rt)), .combine=rbind, .options.snow=opts) %dopar% {
  library(DWLS)
  solveDampenedWLS(sig,rt[,x])
}
stopCluster(cl)
rownames(results) <- colnames(rt)
results <- as.data.frame(results)

## calculate some statistics on the deconvolution results ##
results$pearsonCor <- sapply(c(1:nrow(results)),function(x) cor(sig %*% as.numeric(results[x,c(1:ncol(sig))]),rt[,x]))
results$pearsonCorPval <- sapply(c(1:nrow(results)),function(x) cor.test(sig %*% as.numeric(results[x,c(1:ncol(sig))]),rt[,x])$p.value)
results$spearmanCor <- sapply(c(1:nrow(results)),function(x) cor(sig %*% as.numeric(results[x,c(1:ncol(sig))]),rt[,x],method = "spearman"))
results$spearmanCorPval <- sapply(c(1:nrow(results)),function(x) cor.test(sig %*% as.numeric(results[x,c(1:ncol(sig))]),rt[,x],method = "spearman")$p.value)
results$rmse <- sapply(c(1:nrow(results)),function(x) sqrt(mean((sig %*% as.numeric(results[x,c(1:ncol(sig))])-rt[,x])**2)))

## calculate p-value for deconvolution results, based on random tumor samples (same procedure as used in cibersort) ##
## number of permutations ##

perm <- 1000
cl <- makeCluster(nCores)
registerDoSNOW(cl)
pb <- txtProgressBar(max = perm, style = 3)
progressBar <- function(n) setTxtProgressBar(pb, n)
opts <- list(progress=progressBar)
corDist <- foreach ( i = c(1:perm), .combine=c , .options.snow=opts) %dopar% {
  library(DWLS)
  set.seed(i)
  rtlist <- as.list(data.matrix(rt))
  randomMixture <- as.numeric(rtlist[sample(length(rtlist),nrow(sig))])
  randomMixture <- randomMixture/sum(randomMixture)
  cf <- solveDampenedWLS(sig,randomMixture)
  return(cor(sig %*% cf,randomMixture))
}

results$dwlsPval <- sapply(c(1:nrow(results)),function(x) sum(corDist > results$pearsonCor[x]))

## save results ##
write.table(results, path)
