
library(Seurat)
library(reticulate)
use_condaenv("/home/rschwartz/anaconda3/envs/r4.3/")
library(sceasy)
library(argparse)
library(tidyr)
source("/space/grp/rschwartz/rschwartz/nextflow_eval_pipeline/bin/seurat_functions.R")
options(future.globals.maxSize = 3 * 1024^3)  # 2 GB

parser = argparse::ArgumentParser(description = "Convert H5AD to H5Seurat.")
parser$add_argument("--seurat_obj", type="character", help="Path to rds file.", default = "/space/grp/rschwartz/rschwartz/single_cell_stuff/homo_sapiens/rds/refs/whole_cortex.subsample.500.rds")
parser$add_argument("--normalization_method", type="character", help="Normalization method", default="LogNormalize")
parser$add_argument("--dims", type="integer", help="Number of dimensions", default=50)
parser$add_argument("--batch_key", type="character", help="Batch key", default="dataset_title")
args = parser$parse_args()
seurat_obj = args$seurat_obj
normalization_method = args$normalization_method
dims = args$dims
batch_key = args$batch_key

seurat_obj <- readRDS(seurat_obj)

batch_key_value <- eval(parse(text = paste0("seurat_obj$", batch_key)))
seurat_obj[["RNA"]] <- split(seurat_obj[["RNA"]], f = batch_key_value)

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

#seurat_obj <- RunUMAP(seurat_obj, reduction = "pca", dims = 1:20)
#p <- DimPlot(seurat_obj, group.by = batch_key)

saveRDS(seurat_obj, basename(args$seurat_obj))

