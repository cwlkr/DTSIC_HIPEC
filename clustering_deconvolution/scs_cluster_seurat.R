# partly adapted fom Olbrecht et al., 2021, Genome Medicine

require(tidyverse)

if (!require("BiocManager", quietly = TRUE))
  install.packages("BiocManager")

require(future)
require('Seurat')
require("rtracklayer")

# data availability see Olbrecht et al., 2021, Genome Medicine
count_matrix <- readRDS("olbrecht_data/2095-Olbrecht_counts_matrix.rds")

meta <- read_csv("olbrecht_data/2093-Olbrecht_metadata.csv")

cell_names <- data.frame(cell_label = colnames(count_matrix),
                         sample_name = sub(".*_(.*)",
                                           "\\1",
                                           colnames(count_matrix)),
                         stringsAsFactors = FALSE,
                         row.names = colnames(count_matrix))

meta_seurat <- left_join(x = cell_names,
                         y = as.data.frame(meta),
                         by = "sample_name")
rownames(meta_seurat) <- meta_seurat$cell_label

print(paste0('nr genes: ',dim(count_matrix)[1], ' nr cells: ', dim(count_matrix)[2]))

pbmc <- CreateSeuratObject(count_matrix, min.cells = 3, min.features = 150, meta.data = meta_seurat)
colnames(count_matrix)
rownames(count_matrix)

pbmc[["percent.mt"]] <- PercentageFeatureSet(pbmc, pattern = "^MT-")

data <- GetAssayData(pbmc,slot = "counts")
gtf <- rtracklayer::import("ref_gen/Homo_sapiens.GRCh38.104.gtf")
gtf_df <- as.data.frame(gtf)

geneBiotypes <- c("IG_D_gene","IG_J_gene","IG_C_gene","IG_V_gene","protein_coding","TR_C_gene","TR_D_gene","TR_J_gene","TR_V_gene")
proteinCodingGenes <- unique(gtf_df$gene_name[gtf_df$gene_biotype %in% geneBiotypes])
pbmc[["percent.proteinCoding"]] <- (apply(data[proteinCodingGenes[proteinCodingGenes %in% rownames(data)],],2,sum)/apply(data,2,sum))*100


# QC plots of number of genes detected, number of UMIs and percentage of mitochondrial reads
VlnPlot(object = pbmc, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), ncol = 3)

FeatureScatter(object = pbmc, feature1 = "percent.proteinCoding", feature2 = "percent.mt")

plot1 <- FeatureScatter(pbmc, feature1 = "nCount_RNA", feature2 = "percent.mt")
plot2 <- FeatureScatter(pbmc, feature1 = "nCount_RNA", feature2 = "nFeature_RNA")
CombinePlots(plots = list(plot1, plot2))


pbmc.s <- subset(pbmc, subset = nCount_RNA > 1000 & nCount_RNA < 100000 & percent.mt < 50)
print(pbmc.s)
pbmc.s <- subset(x = pbmc.s, features = proteinCodingGenes[proteinCodingGenes %in% rownames(pbmc.s)] )


# pbmc.s
pbmc.s <- FindVariableFeatures(object = pbmc.s)

plot1 <- VariableFeaturePlot(pbmc.s)
top10 <- head(VariableFeatures(pbmc.s), 10)
plot2 <- LabelPoints(plot = plot1, points = top10, repel = TRUE)

CombinePlots(plots = list(plot1, plot2))

pbmc.s <- NormalizeData(pbmc.s, normalization.method = "LogNormalize", scale.factor = 10000)


cc.genes <- readLines(con = "olbrecht_data/regev_lab_cell_cycle_genes.txt")
s.genes <- cc.genes[1:43]
g2m.genes <- cc.genes[44:97]

pbmc.s <- CellCycleScoring(object = pbmc.s,
                      s.features = s.genes,
                      g2m.features = g2m.genes,
                      set.ident = TRUE)


all.genes <- rownames(pbmc.s)

plan(strategy = "multicore", workers = 18)
options(future.globals.maxSize= (500 * 1024 ^ 2)*10)

pbmc.s <- ScaleData(pbmc.s, features = all.genes, vars.to.regress = c("nCount_RNA",
                                                                     "percent.mito",
                                                                     "S.Score",
                                                                     "G2M.Score"))
saveRDS(pbmc.s, 'seurat_output/scaled_out_countpercsg2m.rds')

# library(Seurat)
# pbmc.s = readRDS('scaled_out_countpercsg2m.rds')

pbmc.s <- RunPCA(pbmc.s, features = VariableFeatures(object = pbmc.s))
p1 <- ElbowPlot(pbmc.s)
print(p1)
# find neighbors how many?
pbmc.s <- FindNeighbors(pbmc.s, dims = 1:20)
pbmc.s <- FindClusters(pbmc.s, resolution = 0.5)

pbmc.s <- RunUMAP(pbmc.s, dims = 1:20)
#pdf(file = 'umap_custering_plbrecht.pdf')
DimPlot(pbmc.s, reduction = "umap")
#dev.off()

original_cluster_idents <- Idents(pbmc.s)
Idents(pbmc.s) <- pbmc.s$sample_name

#pdf(file = 'umap_custering_olbrecht_samplename.pdf')
DimPlot(pbmc.s, reduction = "umap")
#dev.pff()
Idents(pbmc.s) <- original_cluster_idents


#pdf(file = 'gene_names_olbrecht.pdf')
DotPlot(object = pbmc.s,features = rev(c("ESR1","AR","PGR","ERBB2","KRT15","KRT23","KRT81","CDH1","EPCAM","CD24","CD44","MKI67","KRT5","KRT14","KRT17","MYH11","MCAM","MYL9","ACTA2","SOX18","COL1A1","CXCL12","SPARC","VIM","POSTN","MMP11","FN1","VWF","CDH5","KDR","PTPRB","PTPRC","CD14","CD74","CD68","CD163","CSF1R","CD3E","CD3D","CD4","CD8A","CD79A","SSR4","IGLL5","IGLL1","CTSG","CPA3"))) + RotatedAxis()
#dev.off()

plan(strategy = "multicore", workers = 18)
diff_expt = FindAllMarkers(pbmc.s)
summary(diff_expt)

saveRDS(pbmc.s, 'seurat_output/seurat_object.rds')
saveRDS(diff.frame, 'seurat_output/markers.rds')


# write the differently expressed genes into excel sheets
# library(openxlsx)
#wb <- createWorkbook("deg_per_c_all_olbrecht_n.xlsx")

#for(i in 1:(max(as.numeric(diff.expr.signif$cluster)))){
#  addWorksheet(wb, paste("cluster", (i -1)))
#}
#for(i in 1:(max(as.numeric(diff.expr.signif$cluster)))){
#  writeData(wb, i, diff.expr.signif%>%filter(cluster==(i-1)))
#}
#saveWorkbook(wb, file = "cleaned_rerun/deg_per_c_all_olbrecht_n.xlsx", overwrite = TRUE)
