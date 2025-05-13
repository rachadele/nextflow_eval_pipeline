#source("/space/grp/rschwartz/rschwartz/nextflow_eval_pipeline/bin/seurat_functions.R")
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
parser$add_argument("--integration_method", type="character", help="Integration method of query and reference", default="pcaproject")
parser$add_argument("--ref_keys", type="character", nargs="*", help="List of reference keys to pass to query_transfer", default=c("subclass", "class","family","global"))
parser$add_argument("--dims", type="integer", help="Number of dimensions", default=50)
parser$add_argument("--max.features", type="integer", help="Maximum number of features", default=200)
parser$add_argument("--k.anchor", type="integer", help="Number of anchors", default=5)
parser$add_argument("--k.score", type="integer", help="?", default=30)
parser$add_argument("--cutoff", type="numeric", help="Cutoff threshold for label transfer prediction scores", default=0)
parser$add_argument("--ref_path", type="character", help="path to references", default="/space/grp/rschwartz/rschwartz/hs_nf_results/3f/360d86ce9fcedb9ffe518a0a1fb90d/whole_cortex.rds")
parser$add_argument("--query_path", type="character", help="path to query", default = "/space/grp/rschwartz/rschwartz/hs_nf_results/3a/48b09b36aece477cd81d80559c3932/rosmap_R7944883_processed.rds")
parser$add_argument("--k.weight", type="integer", help="k.weight", default=50)
parser$add_argument("--normalization_method", type="character", help="Normalization method", default="LogNormalize")
parser$add_argument("--nfeatures", type="integer", help="Number of variable features to use for dim reduction", default=2000)


args <- parser$parse_args()
# Extract arguments from the parsed list    
#batch_key <- args$batch_key
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
normalization_method <- args$normalization_method
n_features <- args$nfeatures

ref = readRDS(ref_path)
query = readRDS(query_path)

transfer_multiple_labels <- function(
            query, reference, reduction, normalization_method="SCT",
            ref_keys, dims, max.features, 
            k.anchor, k.score, k.weight, project.query=NULL) {
    
    #threshold_dims <- 2000
    if (is.null(project.query)) {
        if ((ncol(query) > ncol(reference))) {
            project.query = TRUE
        } else {
            project.query = FALSE
        }
    }


    anchors <- FindTransferAnchors(reference=reference, 
        normalization.method=normalization_method, 
        query=query, 
        # assay should default to whatever normalization was applied
        npcs=dims, dims=1:dims, 
        reduction = reduction, # had to remove project query for SCTransform to work
        max.features=max.features, k.anchor=k.anchor, k.score=k.score,
        reference.reduction="pca") # use precomputed PCA from seurat_processing step


    k.weight = min(k.weight, floor(nrow(anchors@anchors) / k.score ))
    key = ref_keys[1] # assumes keys are ordered
    #change the k.weight back to 50 or dynamically set ?
    predictions <- TransferData(anchorset = anchors, refdata=key, reference=reference, weight.reduction=reduction, k.weight = k.weight)
    return(predictions)

}


prediction_scores <- transfer_multiple_labels(normalization_method=normalization_method,
        query=query, reference=ref, reduction=integration_method, 
        ref_keys=ref_keys, dims=dims, 
        max.features=max.features, 
        k.anchor=k.anchor, k.score=k.score, 
        project.query=project.query, k.weight=k.weight)



query_name = basename(query_path) %>% gsub("_processed.rds", "", .)
ref_name = basename(ref_path) %>% gsub(".rds", "", .)

prediction_scores <- prediction_scores %>%
  as.data.frame() %>%
  select(-any_of(c("predicted.id", "prediction.score.max"))) %>%
  rename_with(~ gsub("prediction.score.", "", .)) %>%
  rename_with(~ gsub("\\.", " ", .)) %>%
  rename_with(~ gsub("L([0-9]+) ([0-9]+) (.*)", "L\\1/\\2 \\3", .)) %>%
  rename_with(~ gsub("L2/3 6 IT", "L2/3-6 IT", .)) %>%
  rename_with(~ gsub("CA1 ProS", "CA1-ProS", .)) %>%
  rename_with(~ gsub("L4 RSP ACA", "L4 RSP-ACA", .)) %>%
  rename_with(~ gsub("CA2 IG FC", "CA2-IG-FC", .)) %>%
  rename_with(~ gsub("SUB ProS", "SUB-ProS", .)) %>%
  rename_with(~ gsub("L4 IT ET", "L4 IT/ET", .)) %>%
  rename_with(~ gsub("Cajal Retzius", "Cajal-Retzius", .))


write.table(prediction_scores, file=paste0(query_name,"_",ref_name,"_prediction_scores_seurat.tsv"), sep="\t", row.names=FALSE)
write.table(query@meta.data, file=paste0(query_name,".obs.relabel.tsv"), row.names=FALSE, sep= "\t")