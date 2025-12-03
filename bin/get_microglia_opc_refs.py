
#!/user/bin/python3

import warnings
warnings.filterwarnings("ignore")
import os
import random
import argparse
import yaml
import scanpy as sc
import numpy as np
import pandas as pd
import anndata as ad
import cellxgene_census
import cellxgene_census.experimental
import scvi
import utils
from utils import *


# Function to parse command line arguments
def parse_arguments():
    parser = argparse.ArgumentParser(description="Download model file based on organism, census version, and tree file.")
    parser.add_argument('--organism', type=str, default='mus_musculus', help='Organism name (e.g., homo_sapiens)')
    parser.add_argument('--census_version', type=str, default='2024-07-01', help='Census version (e.g., 2024-07-01)')
    parser.add_argument('--subsample_ref', type=int, default=50)
    parser.add_argument('--relabel_path', type=str, default="/space/grp/rschwartz/rschwartz/nextflow_eval_pipeline/meta/census_map_mouse_author.tsv")
    parser.add_argument('--ref_collections', type=str, nargs = '+', default = [
        "A taxonomy of transcriptomic cell types across the isocortex and hippocampal formation",
        "An integrated transcriptomic and epigenomic atlas of mouse primary motor cortex cell types",
        "Tabula Muris Senis"
    ]) 
    parser.add_argument('--split_column', type=str, default="dataset_id")
    parser.add_argument('--ref_keys', type=str, nargs="+", default=["subclass","class","family","global"])
    parser.add_argument('--seed', type=int, default=42)
    parser.add_argument('--ref_cell_type', type=str, default="CNS macrophage", help="Cell type to filter for reference (e.g., 'CNS macrophage')")
    if __name__ == "__main__":
        known_args, _ = parser.parse_known_args()
        return known_args
        

def get_aggregate_census(census_version="2024-07-01", organism="mus_musculus", subsample=50, split_column="dataset_id", dims=50, organ="brain",
               ref_collections=None,
               relabel_path="/space/grp/rschwartz/rschwartz/nextflow_eval_pipeline/meta/census_map_mouse_author.tsv", 
               seed=42, 
               ref_keys=["subclass","class","family","global"],
               original_celltypes=None):

    census = cellxgene_census.open_soma(census_version=census_version)
    dataset_info = census.get("census_info").get("datasets").read().concat().to_pandas()
    cellxgene_obs = get_filtered_obs(census, organism, organ=organ, is_primary=True, disease="normal")
    
    cellxgene_obs = cellxgene_obs.merge(dataset_info, on="dataset_id", suffixes=(None,"_y"))
    cellxgene_obs.drop(columns=['soma_joinid_y'], inplace=True)
    cellxgene_obs_filtered = cellxgene_obs[cellxgene_obs['collection_name'].isin(ref_collections)] 

    # Adjust organism naming for compatibility
    organism_name_mapping = {
        "homo_sapiens": "Homo sapiens",
        "mus_musculus": "Mus musculus"
    }
    organism = organism_name_mapping.get(organism, organism)

    cell_columns = [
        "assay", "cell_type", "tissue",
        "tissue_general", "suspension_type",
        "disease", "dataset_id", "development_stage",
        "soma_joinid","observation_joinid"
    ]
    
    # add author_cell_type to obs
    # this will enable proper relabeling and subsampling
    # need to add it back in after getting ids
    if isinstance(original_celltypes, pd.DataFrame) and not original_celltypes.empty:
        cellxgene_obs_filtered = map_author_labels(cellxgene_obs_filtered, original_celltypes)
   
    filtered_ids = cellxgene_obs_filtered['soma_joinid'].values
    adata = extract_data(
        cellxgene_obs_filtered, filtered_ids,
        subsample=subsample, organism=organism,
        census=census,
        cell_columns=cell_columns, 
        dataset_info=dataset_info, dims=dims,
        relabel_path=relabel_path, 
        ref_keys=ref_keys, seed = seed, 
        original_celltypes=original_celltypes
    )
    # Clean up categorical columns
    for col in adata.obs.columns:
        if adata.obs[col].dtype.name == 'category':
            adata.obs[col] = pd.Categorical(adata.obs[col].cat.remove_unused_categories())
    return adata
         
def main():
    # Parse command line arguments
    args = parse_arguments()

    # Set organism and census_version from arguments
    organism = args.organism
    census_version = args.census_version
    subsample_ref = args.subsample_ref
    relabel_path = args.relabel_path
    ref_collections = args.ref_collections
    split_column = args.split_column
    ref_keys = args.ref_keys
    SEED = args.seed
    ref_cell_type = args.ref_cell_type

    if organism == "mus_musculus":
        original_celltypes = get_original_celltypes(columns_file=f"/space/grp/rschwartz/rschwartz/nextflow_eval_pipeline/meta/author_cell_annotations/{census_version}/original_celltype_columns.tsv",
                                                    author_annotations_path=f"/space/grp/rschwartz/rschwartz/nextflow_eval_pipeline/meta/author_cell_annotations/{census_version}") 
    else:
        original_celltypes = None


    adata = get_aggregate_census(
        census_version=census_version,
        organism=organism,
        subsample=subsample_ref,
        split_column=split_column,
        relabel_path=relabel_path,
        ref_collections=ref_collections,
        seed=SEED,
        ref_keys=ref_keys,
        original_celltypes=original_celltypes
    )

    # Filter to CNS macrophage after author label mapping
    if ref_cell_type == "CNS macrophage":
        if "family" in adata.obs.columns:
            ref_ct_subset = adata[adata.obs["family"] == ref_cell_type].copy()
        else:
            raise ValueError("'family' column not found in AnnData obs after label mapping.")
    elif ref_cell_type == "OPC":
        ref_ct_subset = adata[adata.obs["class"] == ref_cell_type].copy()  
        

    # make naming dynamic
    
    
    outdir = "refs"
    os.makedirs(outdir, exist_ok=True)
    ref_ct_subset.write(os.path.join(outdir, f"{ref_cell_type.lower()}_reference.h5ad"))
    ref_ct_subset.obs.to_csv(os.path.join(outdir, f"{ref_cell_type.lower()}_reference.obs.tsv"), sep="\t")

      
if __name__ == "__main__":
    main()