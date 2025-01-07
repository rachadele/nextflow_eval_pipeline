#source("/space/grp/rschwartz/rschwartz/nextflow_eval_pipeline/bin/get_census.R")
source("/space/grp/rschwartz/rschwartz/nextflow_eval_pipeline/bin/seurat_functions.R")
library(reticulate)
use_condaenv("~/anaconda3/envs/r4.3/")
#load_required_packages()
library(Seurat)
library(dplyr)
library(gplots)
library(argparse)
library(cellxgene.census)
set.seed(123)

parser <- ArgumentParser(description = "Get .rds file references from cellxgene census")
parser$add_argument("--organism", type="character", help="Organism", default="homo_sapiens")
parser$add_argument("--subsample_ref", type="numeric", help="Number of cells to subsample from each cell type in reference", default=5)
parser$add_argument("--split_column", type="character", help="How to split census references", default="dataset_id")
parser$add_argument('--ref_collections', type="character", nargs='+', help="List of reference collections", 
                        default = c("Transcriptomic cytoarchitecture reveals principles of human neocortex organization", "SEA-AD: Seattle Alzheimer’s Disease Brain Cell Atlas"))
parser$add_argument("--relabel_path", type="character", help="Path to relabel file", default="/space/grp/rschwartz/rschwartz/nextflow_eval_pipeline/meta/census_map_human.tsv")
parser$add_argument("--census_version", type="character", help="Census version", default = "2024-07-01")

args <- parser$parse_args()

# Parse the comma separated list into a vector
ref_collections <- args$ref_collections
print(ref_collections)
organism = args$organism
subsample_ref = args$subsample_ref
split_column = args$split_column
relabel_r=args$relabel_path
census_version = args$census_version

get_census(organism=organism, census_version=census_version,subsample=subsample_ref, 
    split_column=split_column, relabel_path=relabel_r,
    ref_collections=ref_collections)