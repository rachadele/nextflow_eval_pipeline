#list files
#read in seurat file, subsample to 10k cells per sample
#write to original path gsub "rds" for "subsample.rds"


source("/space/grp/rschwartz/rschwartz/nextflow_eval_pipeline/bin/seurat_functions.R")
library(reticulate)
use_condaenv("/home/rschwartz/anaconda3/envs/r4.3/")
#load_required_packages()
library(argparse)
library(gplots)
#set.seed(123)
library(dplyr)
library(ggplot2)
parser <- ArgumentParser(description = "Subsample seurat.")
parser$add_argument("--subsample_query", type="character", help="number of cells to subsample", default=100)
parser$add_argument("--query_path", type="character", help="query_path", required=TRUE)
parser$add_argument("--relabel_path", type="character", help="relabel_path", required=TRUE)
parser$add_argument("--subset_columns", type="character", nargs='+', help="Columns in meta.data to split by", required=FALSE)
parser$add_argument("--subset_values", type="character", nargs='+', help="Values in subset_columns to filter by", required=FALSE)
parser$add_argument("--ref_keys", type="character", nargs='+', help="Reference keys to use for relabeling", required=TRUE)

args <- parser$parse_args()
subsample_query=args$subsample_query
query_path=args$query_path
relabel_path=args$relabel_path
subset_columns=args$subset_columns
subset_values=args$subset_values
ref_keys=args$ref_keys

seurat_obj <- readRDS(query_path)

# Split the Seurat object if subset_columns and subset_values are provided
if (!is.null(subset_columns) && !is.null(subset_values)) {
    for (i in seq_along(subset_columns)) {
        seurat_obj <- seurat_obj[ ,seurat_obj@meta.data[[subset_columns]][i] == subset_values[i]]
        #seurat_obj <- subset(seurat_obj, subset = get(subset_columns[i]) == subset_values[i])
    }
}


seurat_mapped <- relabel_wrapper(seurat_obj, relabel_path=relabel_path)
Idents(seurat_mapped) <- seurat_mapped@meta.data[[ref_keys[1]]]
#for subsampling by cell type, see get_census.R
seurat_mapped <- seurat_mapped[ ,seurat_mapped@meta.data[[ref_keys[1]]] != "unknown"]
#seurat_subsampled <- seurat_mapped[, sample(colnames(seurat_mapped), size = subsample_query, replace=F)]

output_file <- gsub(".rds","_processed.rds",basename(query_path))

# Save the subsampled Seurat object
saveRDS(seurat_subsampled, file = output_file)
message(paste("Subsampled Seurat object saved to", output_file))