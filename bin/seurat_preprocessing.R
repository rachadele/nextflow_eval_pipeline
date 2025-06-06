
library(Seurat)
library(reticulate)
use_condaenv("/home/rschwartz/anaconda3/envs/r4.3/")
library(sceasy)
library(argparse)
library(tidyr)
source("/space/grp/rschwartz/rschwartz/nextflow_eval_pipeline/bin/seurat_functions.R")
options(future.globals.maxSize = 5 * 1024^3)  # 5 GB

parser = argparse::ArgumentParser(description = "Convert H5AD to H5Seurat.")
parser$add_argument("--h5ad_file", type="character", help="Path to H5AD file.", default = "/space/grp/rschwartz/rschwartz/nextflow_eval_pipeline/sctransform_sample_subsets_dataset_id_only_mmus/ba/63263d2b41cd71c11ee2be5a43941d/GSE214244.1_1161367_processed.h5ad")
parser$add_argument("--normalization_method", type="character", help="Normalization method", default="SCT")
parser$add_argument("--dims", type="integer", help="Number of dimensions", default=50)
parser$add_argument("--nfeatures", type="integer", help="Number of variable features to use for dim reduction", default=2000)

args = parser$parse_args()
h5ad_file = args$h5ad_file
normalization_method = args$normalization_method
dims = args$dims
nfeatures = args$nfeatures

sceasy_seurat <- sceasy::convertFormat(h5ad_file, from="anndata", to="seurat")

if ("feature_id" %in% colnames(sceasy_seurat@assays$RNA[[]])) {
  sceasy_seurat <- rename_features(sceasy_seurat, column_name="feature_id")
}

if (normalization_method=="LogNormalize") {
  sceasy_seurat <- sceasy_seurat %>% NormalizeData(normalization.method=normalization_method) %>% 
  FindVariableFeatures(nfeatures = nfeatures) %>% 
  ScaleData() %>% RunPCA(npcs=dims)
} else if (normalization_method=="SCT") {
  sceasy_seurat <- sceasy_seurat %>%
    SCTransform(verbose=TRUE, variable.features.n=nfeatures) %>%
    RunPCA(npcs=dims, assay="SCT")
} else {
  stop("Normalization method not recognized.")
}

saveRDS(sceasy_seurat, file = gsub(".h5ad",".rds",h5ad_file))
message(paste("Converted H5AD to RDS and saved to", gsub(".h5ad",".rds",h5ad_file)))