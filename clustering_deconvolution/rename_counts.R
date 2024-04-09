require(BiocManager)
#BiocManager::install("AnnotationDbi")
library(EnsDb.Hsapiens.v75)
library(AnnotationDbi)
library(tidyverse)


edb <- EnsDb.Hsapiens.v75

## Get all transcripts defined in Ensembl (version 75):
tx <- genes(edb, columns=c("tx_id", "gene_id", "gene_name"))
mapping <- data_frame(gene_id=tx$gene_id, name=tx$gene_name)

file_list  <- list.files('aligned_rsem/', pattern = 'gz.ge', full.names = T)

## combine everything together into a list of gene to TPM per patient
name <- file_list[1]
n <- str_extract(name, "\\d.*(?=_R)")

g1 <- read_tsv(name, col_types = cols())
mapped_genes = inner_join(g1, mapping, copy=T, by=c('gene_id')) %>% distinct

mapped_genes %>% select(gene_id, expected_count) -> mapped_genes
colnames(mapped_genes) = c('gene_id', n)
#mapped_genes  %>% inner_join(., mapping,by=c('gene_id')) %>% distinct -> out_to_cor


for  (name in file_list[-1]){
  n <- str_extract(name, "\\d.*(?=_R)")

  g2 <- read_tsv(name, col_types = cols())
  mapped_genes2 = inner_join(g2, mapping, copy=T, by=c('gene_id')) %>% distinct

  mapped_genes2 %>% select(gene_id, expected_count) -> mapped_genes2
  colnames(mapped_genes2) <- c('gene_id', n)
  #mapped_genes  %>% inner_join(., mapping,by=c('gene_id')) %>% distinct -> out_to_cor2
  mapped_genes %>% inner_join(., mapped_genes2, by=c('gene_id')) -> mapped_genes
}


mapped_genes  %>% inner_join(., mapping,by=c('gene_id')) %>% distinct -> out_to_cor

## write list out!
write_tsv(out_to_cor, 'readcounts_realigment_starH38104_expected_counts.txt')

name <- file_list[1]
n <- str_extract(name, "\\d.*(?=_R)")

g1 <- read_tsv(name, col_types = cols())
mapped_genes <- inner_join(g1, mapping, copy=T, by=c('gene_id')) %>% distinct

mapped_genes %>% select(gene_id, expected_count) -> mapped_genes
colnames(mapped_genes) = c('gene_id', n)


for  (name in file_list[-1]){
  n <- str_extract(name, "\\d.*(?=_R)")

  g2 <- read_tsv(name, col_types = cols())
  mapped_genes2 <- inner_join(g2, mapping, copy=T, by=c('gene_id')) %>% distinct

  mapped_genes2 %>% select(gene_id, TPM) -> mapped_genes2
  colnames(mapped_genes2) <- c('gene_id', n)
  #mapped_genes  %>% inner_join(., mapping,by=c('gene_id')) %>% distinct -> out_to_cor2
  mapped_genes %>% inner_join(., mapped_genes2, by=c('gene_id')) -> mapped_genes
}

mapped_genes  %>% inner_join(., mapping,by=c('gene_id')) %>% distinct -> out_to_cor

## write list out!
write_tsv(out_to_cor, 'readcounts_realigment_starH38104_tpm.txt')
