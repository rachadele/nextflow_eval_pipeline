
library(Seurat)
library(reticulate)
use_condaenv("/home/rschwartz/anaconda3/envs/r4.3/")
library(sceasy)

parser = argparse::ArgumentParser(description = "Convert H5AD to H5Seurat.")
parser$add_argument("--h5ad_file", type="character", help="Path to H5AD file.", default = "/space/grp/rschwartz/rschwartz/biof501_proj/refs/whole_cortex.h5ad")
parser$add_argument("--normalization_method", type="character", help="Normalization method", default="LogNormalize")
parser$add_argument("--dims", type="integer", help="Number of dimensions", default=50)
parser$add_argument("--nfeatures", type="integer", help="Number of variable features to use for dim reduction", default=2000)

args = parser$parse_args()
h5ad_file = args$h5ad_file
normalization_method = args$normalization_method
dims = args$dims
n_features = args$nfeatures

sceasy_seurat <- sceasy::convertFormat(h5ad_file, from="anndata", to="seurat")
sceasy_seurat <- sceasy_seurat %>% NormalizeData(normalization.method=normalization_method) %>% FindVariableFeatures(nfeatures = n_features) %>% ScaleData() %>% RunPCA(npcs=dims)

saveRDS(sceasy_seurat, file = gsub(".h5ad",".rds",h5ad_file))
message(paste("Converted H5AD to RDS and saved to", gsub(".h5ad",".rds",h5ad_file)))