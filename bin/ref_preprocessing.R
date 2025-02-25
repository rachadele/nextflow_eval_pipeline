
library(Seurat)
library(reticulate)
use_condaenv("/home/rschwartz/anaconda3/envs/r4.3/")
library(sceasy)
library(argparse)
library(tidyr)
source("/space/grp/rschwartz/rschwartz/nextflow_eval_pipeline/bin/seurat_functions.R")
options(future.globals.maxSize = 5 * 1024^3)  # 5 GB

parser = argparse::ArgumentParser(description = "Convert H5AD to H5Seurat.")
parser$add_argument("--h5ad_file", type="character", help="Path to H5AD file.", default = "/space/grp/rschwartz/rschwartz/hs_nf_results/c6/10de8aba87d83363dcbde1602eb72c/whole_cortex.h5ad")
parser$add_argument("--normalization_method", type="character", help="Normalization method", default="LogNormalize")
parser$add_argument("--dims", type="integer", help="Number of dimensions", default=30)
parser$add_argument("--batch_key", type="character", help="Batch key", default="dataset_title")
parser$add_argument("--nfeatures", type="integer", help="Number of variable features to use for dim reduction", default=2000)
parser$add_argument("--batch-correct", action="store_true", help="Batch correct", default=TRUE)
parser$add_argument("--k.anchor", type="integer", help="Number of anchors", default=5)
parser$add_argument("--k.weight", type="integer", help="k.weight", default=30)
parser$add_argument("--k.score", type="integer", help="k.score", default=30)
args = parser$parse_args()
normalization_method = args$normalization_method
dims = args$dims
batch_key = args$batch_key
n_features = args$nfeatures
h5ad_file = args$h5ad_file
integration_method = args$integration_method
batch_correct = args$batch_correct
sceasy_seurat <- sceasy::convertFormat(h5ad_file, from="anndata", to="seurat")
k.anchor = args$k.anchor
k.weight = args$k.weight
k.score = args$k.score

batch_key_value <- eval(parse(text = paste0("sceasy_seurat$", batch_key)))

if ("feature_id" %in% colnames(sceasy_seurat@assays$RNA[[]])) {
  sceasy_seurat <- rename_features(sceasy_seurat, column_name="feature_id")
}

# If integration method is not null
if (batch_correct) {

    # Split by batch for integration
    sceasy_list <- SplitObject(sceasy_seurat, split.by = batch_key)
	if (length(sceasy_list) > 1) {
		sceasy_list <- Filter(function(x) ncol(x) > dims, sceasy_list)

		if (normalization_method == "SCT") {
			assay <- "SCT"
			sceasy_list <- lapply(sceasy_list, function(x) SCTransform(x, 
										verbose = FALSE, 
										variable.features.n = n_features))

			# Select integration features
			features <- SelectIntegrationFeatures(sceasy_list)		
			# If SCT integration, use PrepSCTIntegration
			sceasy_list <- PrepSCTIntegration(
								sceasy_list,
								anchor.features = features
								)

		} else if (normalization_method == "LogNormalize") {
			assay <- "RNA"
			sceasy_list <- lapply(sceasy_list, function(x) {
				x <- NormalizeData(x)
				x <- FindVariableFeatures(x)
				x <- ScaleData(x)
				return(x)
			})
			features <- SelectIntegrationFeatures(sceasy_list)
		}
		# Find integration anchors
		anchors <- FindIntegrationAnchors(object.list = sceasy_list, 
					dims = 1:dims, k.anchor=k.anchor, k.score=k.score,
					normalization.method = normalization_method, 
					anchor.features = features)


		anchors_per_pair <- table(anchors@anchors$dataset1, anchors@anchors$dataset2)		
		min_anchors <- min(anchors_per_pair[anchors_per_pair > 0])  # Ignore zero values
		k.anchor <- min(min_anchors, k.anchor)
		# Integrate data
		sceasy_seurat <- IntegrateData(anchorset = anchors, 
				dims = 1:dims, k.weight=k.weight, 
				normalization.method = normalization_method, 
				new.assay.name = "integrated")
		# Run PCA on integrated data
		sceasy_seurat <- sceasy_seurat %>% ScaleData() %>%
			RunPCA(npcs = dims, assay = "integrated")
	} else {
		message("Only one batch found. Skipping integration.")
	}
} else {

    # Normalize and scale data without integration
    if (normalization_method == "SCT") {
        sceasy_seurat <- sceasy_seurat %>%
            SCTransform(verbose = FALSE, variable.features.n = n_features) %>%
            RunPCA(npcs = dims, assay = "SCT")
    } else if (normalization_method == "LogNormalize") {
        sceasy_seurat <- sceasy_seurat %>%
            NormalizeData(normalization.method = normalization_method) %>%
            FindVariableFeatures(nfeatures = n_features) %>%
            ScaleData() %>%
            RunPCA(npcs = dims)
    }
}

sceasy_seurat <- RunUMAP(sceasy_seurat, reduction = "pca", dims = 1:30)
p <- DimPlot(sceasy_seurat, group.by = batch_key, label = TRUE)
p

saveRDS(sceasy_seurat, file = gsub(".h5ad", ".rds", h5ad_file))
message(paste("Converted H5AD to RDS and saved to", gsub(".h5ad", ".rds", h5ad_file)))
