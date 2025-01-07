
library(Seurat)
library(reticulate)
use_condaenv("/home/rschwartz/anaconda3/envs/r4.3/")
library(sceasy)

parser = argparse::ArgumentParser(description = "Convert H5AD to H5Seurat.")
parser$add_argument("--h5ad_file", type="character", help="Path to H5AD file.", default = "/space/grp/rschwartz/rschwartz/biof501_proj/refs/whole_cortex.h5ad")
args = parser$parse_args()
h5ad_file = args$h5ad_file

sceasy_seurat <- sceasy::convertFormat(h5ad_file, from="anndata", to="seurat")
saveRDS(sceasy_seurat, file = gsub(".h5ad",".rds",h5ad_file))
message(paste("Converted H5AD to RDS and saved to", gsub(".h5ad",".rds",h5ad_file)))