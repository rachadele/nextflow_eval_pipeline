
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
import warnings
import cellxgene_census
import cellxgene_census.experimental
import scvi
from sklearn.ensemble import RandomForestClassifier
import adata_functions
from adata_functions import *
from pathlib import Path
from sklearn.ensemble import RandomForestClassifier
from sklearn.model_selection import train_test_split
from sklearn.metrics import accuracy_score, classification_report, confusion_matrix, roc_curve, auc
from sklearn.preprocessing import label_binarize
import matplotlib.pyplot as plt
import seaborn as sns
import json
import argparse
import os
import json
from types import SimpleNamespace
import yaml

# Function to parse command line arguments
def parse_arguments():
    parser = argparse.ArgumentParser(description="Classify cells given 1 ref and 1 query")
  #  parser.add_argument('--census_version', type=str, default='2024-07-01', help='Census version (e.g., 2024-07-01)')
    parser.add_argument('--tree_file', type=str, default="/space/grp/rschwartz/rschwartz/nextflow_eval_pipeline/meta/master_hierarchy.json")
    parser.add_argument('--query_path', type=str, default="/space/grp/rschwartz/rschwartz/nextflow_eval_pipeline/work/60/ed91620a1659a23ba5b68f8c21028c/lim_Cingulate.obs.relabel.tsv")
    parser.add_argument('--ref_name', type=str, default="Dissection:_Anterior_cingulate_cortex_ACC.h5ad") #nargs ="+")
    parser.add_argument('--ref_keys', type=str, nargs='+', default=["rachel_subclass", "rachel_class", "rachel_family"])
    parser.add_argument('--cutoff', type=float, default=0, help = "Cutoff threshold for positive classification")
    parser.add_argument('--probs', type=str, default="/space/grp/rschwartz/rschwartz/nextflow_eval_pipeline/work/4a/164f6559104b8e872e88bc617411a2/probs/lim_Cingulate_processed_Dissection:_Anterior_cingulate_cortex_ACC.prob.df.tsv")
    parser.add_argument('--mapping_file', type=str, default="/space/grp/rschwartz/rschwartz/nextflow_eval_pipeline/meta/census_map_human.tsv")
    parser.add_argument('--query_tissue', type=str, default="anterior cingulate cortex")
    parser.add_argument('--ref_tissue_mapping', type=str)
    
    if __name__ == "__main__":
        known_args, _ = parser.parse_known_args()
        return known_args

    
    
def main():
    SEED = 42
    random.seed(SEED)         # For `random`
    np.random.seed(SEED)      # For `numpy`
    # For `torch`'
    scvi.settings.seed = SEED # For `scvi`
    # Parse command line arguments
    args = parse_arguments()
    # query_name = args.query_name
   # ref_name = args.ref_name
    tree_file = args.tree_file
    query_path = args.query_path
    ref_name = args.ref_name
    ref_keys = args.ref_keys
    cutoff = args.cutoff
    query_tissue = args.query_tissue.replace("_", " ")
    ref_tissue_mapping = args.ref_tissue_mapping
    
    # Load data
    ref_tissue_mapping = yaml.load(open(ref_tissue_mapping), Loader=yaml.FullLoader)
    ref_tissue=ref_tissue_mapping[ref_name]
    
    prob_df = pd.read_csv(args.probs, sep="\t")
    mapping_df = pd.read_csv(args.mapping_file, sep="\t")
    query_name = os.path.basename(query_path).replace(".obs.relabel.tsv", "")
    query = pd.read_csv(query_path, sep="\t")
    
    #ref_name = os.path.basename(ref_path).replace(".h5ad", "").replace(".rds", "")
    # Read the JSON tree file
    with open(tree_file, 'r') as file:
        tree = json.load(file)
        
    #rocs = roc_analysis(probabilities=prob_df, query=query, key=ref_keys[0])
    #outdir = os.path.join("roc", query_name, ref_name)
  #  os.makedirs(outdir, exist_ok=True)
   # plot_roc_curves(metrics=rocs, title=f"{query_name} vs {ref_name}", save_path=os.path.join(outdir, "roc_results.png"))
    #roc_df = process_roc(rocs, ref_name=ref_name, query_name=query_name)
    #roc_df.to_csv(os.path.join(outdir, f"{query_name}_{ref_name}.roc.df.tsv"),sep="\t", index=False)
        
    # Classify cells and evaluate
    query = classify_cells(query, ref_keys, cutoff=cutoff, probabilities=prob_df, tree=tree)
    outdir = "predicted_meta"
    os.makedirs(outdir, exist_ok=True)
    query.to_csv(os.path.join(outdir,f"{query_name}_{ref_name}.predictions.{cutoff}.tsv"), index=False, sep="\t")
    class_metrics = eval(query, ref_keys, mapping_df)
    class_metrics = update_classification_report(class_metrics, ref_keys)

    # Plot confusion matrices
    for key in ref_keys:
        outdir = os.path.join("confusion", query_name, ref_name)
        os.makedirs(outdir, exist_ok=True)
        plot_confusion_matrix(query_name, ref_name, key, class_metrics[key]["confusion"], output_dir=outdir)

    # Collect F1 scores
    f1_data = []
    for key in ref_keys:
        classification_report = class_metrics[key]["classification_report"]
        for label, metrics in classification_report.items():
            if label not in ["macro avg", "micro avg", "weighted avg", "accuracy"]:
                f1_data.append({
                    'query': query_name,
                    'reference': ref_name,
                    'label': label,
                    'f1_score': metrics['f1-score'],
                    'macro_f1': classification_report.get('macro avg', {}).get('f1-score', None),
                    'micro_f1': classification_report.get('micro avg', {}).get('f1-score', None),
                    'weighted_f1': classification_report.get('weighted avg', {}).get('f1-score', None),
                    'key': key,
                    'cutoff': cutoff,
                    'query_tissue': query_tissue,
                    'ref_tissue': ref_tissue
                })

    # Save F1 scores to a file
    df = pd.DataFrame(f1_data)
    outdir = "f1_results"
    os.makedirs(outdir, exist_ok=True)
    df.to_csv(os.path.join(outdir, f"{query_name}_{ref_name}.f1.scores.tsv"), sep="\t", index=False)
    
if __name__ == "__main__":
    main()