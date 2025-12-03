
#!/user/bin/python3

from pathlib import Path
import os
import sys
import scanpy as sc
import numpy as np
import pandas as pd
import anndata as ad
import cellxgene_census
import scvi
from scipy.sparse import csr_matrix

import argparse
import anndata as ad
import os

def parse_arguments():
    parser = argparse.ArgumentParser(description="Subset a query h5ad file to a given cell type.")
    parser.add_argument('--query_path', type=str, default="/space/grp/rschwartz/rschwartz/nextflow_eval_pipeline/mmus/55/827a840f83dad8889cca0217fd1b86/GSE247339.2_1051974_GSM7887403_processed.h5ad", help='Path to input query h5ad file')
    parser.add_argument('--cell_type', type=str, default="OPC", help='Cell type to subset')
    parser.add_argument('--obs_column', type=str, default='family', help='Column in .obs to filter on (default: cell_type)')
    if __name__ == "__main__":
        known_args, _ = parser.parse_known_args()
        return known_args
        

def main():
    args = parse_arguments()
    adata = ad.read_h5ad(args.query_path)
    query_name = os.path.basename(args.query_path).replace("_processed.h5ad", "")
    cell_type = args.cell_type
    if args.obs_column not in adata.obs.columns:
        raise ValueError(f"Column '{args.obs_column}' not found in AnnData.obs.")
    
    if cell_type == "CNS macrophage":
        subset = adata[adata.obs[args.obs_column] == cell_type].copy()
    elif cell_type == "OPC":
        subset = adata[adata.obs[args.obs_column].isin(["OPC","Oligodendrocyte", "Neural stem cell"])].copy()
        
    ct_name = args.cell_type.replace(" ", "_")
    out_path = f"{query_name}_{ct_name}_processed.h5ad"
    subset.write_h5ad(out_path)
    print(f"Subsetted AnnData saved to {out_path}")

if __name__ == "__main__":
    main()