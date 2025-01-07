#source("/space/grp/rschwartz/rschwartz/nextflow_eval_pipeline/bin/get_census.R")
source("/space/grp/rschwartz/rschwartz/nextflow_eval_pipeline/bin/seurat_functions.R")
library(reticulate)
use_condaenv("~/anaconda3/envs/r4.3/")
library(Seurat)
library(dplyr)
library(gplots)
library(argparse)
library(cellxgene.census)
library(data.tree)
set.seed(123)


parser <- ArgumentParser(description = "Process Seurat objects and transfer labels.")
#parser$add_argument("--batch_key", type="character", help="Key for integrating dataset", default="sample")
#parser$add_argument("--organism", type="character", help="Organism", default="human")
parser$add_argument("--integration_method", type="character", help="Integration method of query and reference", default="pcaproject")
parser$add_argument("--ref_keys", type="character", nargs="*", help="List of reference keys to pass to query_transfer", default=c("rachel_subclass", "rachel_class","rachel_family"))
parser$add_argument("--dims", type="integer", help="Number of dimensions", default=20)
parser$add_argument("--max.features", type="integer", help="Maximum number of features", default=200)
parser$add_argument("--k.anchor", type="integer", help="Number of anchors", default=5)
parser$add_argument("--k.score", type="integer", help="?", default=30)
#parser$add_argument("--project.query", type="logical", help="whether to project query on to reference PCA (if unset, determined automatically from size of datasets)")
parser$add_argument("--cutoff", type="numeric", help="Cutoff threshold for label transfer prediction scores", default=0)
parser$add_argument("--ref_path", type="character", help="path to references")
parser$add_argument("--query_path", type="character", help="path to query")
parser$add_argument("--k.weight", type="integer", help="k.weight", default=20)

args <- parser$parse_args()
# Extract arguments from the parsed list    
#batch_key <- args$batch_key
#organism <- args$organism
integration_method <- args$integration_method
ref_keys <- args$ref_keys
dims <- args$dims
max.features <- args$max.features
k.anchor <- args$k.anchor
k.score <- args$k.score
project.query <- NULL
k.weight <- args$k.weight
ref_path <- args$ref_path
query_path <- args$query_path

ref = readRDS(ref_path)
query = readRDS(query_path)

#potentially need to add in options
#default = lognormalize, 2000 variable features
#using 50 PCs in nextflow

#could also bach correct reference and set harmony embedding to be named "pca"
#then do pcaproject on harmony embeddings
# or do the same for SCVI_ embeddings from census

#rename features in reference to match ensembl ID
ref <- rename_features(ref, column_name="feature_id")

ref <- ref %>% NormalizeData() %>% FindVariableFeatures() %>% ScaleData() %>% RunPCA(npcs=dims)
query <- query %>% NormalizeData() %>% FindVariableFeatures() %>% ScaleData() %>% RunPCA(npcs=dims)

prediction_scores <- transfer_multiple_labels(
        query=query, reference=ref, reduction=integration_method, 
        ref_keys=ref_keys, dims=dims, max.features=max.features, 
        k.anchor=k.anchor, k.score=k.score, project.query=project.query, k.weight=20)



query_name = basename(query_path) %>% gsub("_processed.rds", "", .)
ref_name = basename(ref_path) %>% gsub(".rds", "", .)

prediction_scores <- prediction_scores %>% as.data.frame() %>%
        select(-c("predicted.id","prediction.score.max")) %>%
        rename_all(~gsub("prediction.score.", "", .)) %>% 
        rename_all(~gsub("\\.", " ", .)) %>%  # Replace dots with spaces
        rename("L2/3-6 IT" = "L2 3 6 IT")  # Correctly rename the column


write.table(prediction_scores, file=paste0(query_name,"_",ref_name,"_prediction_scores_seurat.tsv"), sep="\t", row.names=FALSE)
write.table(query@meta.data, file=paste0(query_name,".obs.relabel.tsv"), row.names=FALSE, sep= "\t")