
library(Seurat)
library(reticulate)
use_condaenv("/home/rschwartz/anaconda3/envs/r4.3/")
library(sceasy)
library(argparse)
library(tidyr)
source("/space/grp/rschwartz/rschwartz/nextflow_eval_pipeline/bin/seurat_functions.R")
options(future.globals.maxSize = 3 * 1024^3)  # 2 GB

parser = argparse::ArgumentParser(description = "Convert H5AD to H5Seurat.")
parser$add_argument("--h5ad_file", type="character", help="Path to H5AD file.", default = "/space/grp/rschwartz/rschwartz/biof501_proj/refs/whole_cortex.h5ad")parser$add_argument("--normalization_method", type="character", help="Normalization method", default="LogNormalize")
parser$add_argument("--dims", type="integer", help="Number of dimensions", default=50)
parser$add_argument("--batch_key", type="character", help="Batch key", default="dataset_title")
parser$add_argument("--nfeatures", type="integer", help="Number of variable features to use for dim reduction", default=2000)

args = parser$parse_args()
normalization_method = args$normalization_method
dims = args$dims
batch_key = args$batch_key
n_features = args$nfeatures
h5ad_file = args$h5ad_file

sceasy_seurat <- sceasy::convertFormat(h5ad_file, from="anndata", to="seurat")

batch_key_value <- eval(parse(text = paste0("seurat_obj$", batch_key)))
seurat_obj[[assay]] <- split(seurat_obj[[assay]], f = batch_key_value)

if (normalization_method == "SCT") {
	assay="SCT"
	seurat_obj <- SCTransform(seurat_obj)

} 
else if (normalization_method == "LogNormalize") {
	assay="RNA"
	seurat_obj <- NormalizeData(seurat_obj)
	seurat_obj <- FindVariableFeatures(seurat_obj)
	seurat_obj <- ScaleData(seurat_obj)
	seurat_obj <- RunPCA(seurat_obj)
}


seurat_obj <- IntegrateLayers(object = seurat_obj, method = HarmonyIntegration, 
				orig.reduction = "pca", new.reduction = "pca", assay=assay,
				normalization.method=normalization_method)


saveRDS(sceasy_seurat, file = gsub(".h5ad",".rds",h5ad_file))
message(paste("Converted H5AD to RDS and saved to", gsub(".h5ad",".rds",h5ad_file)))