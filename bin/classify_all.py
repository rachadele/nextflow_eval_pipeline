
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
from sklearn.metrics import precision_recall_curve, average_precision_score, PrecisionRecallDisplay

# Function to parse command line arguments
def parse_arguments():
    parser = argparse.ArgumentParser(description="Classify cells given 1 ref and 1 query")
  #  parser.add_argument('--census_version', type=str, default='2024-07-01', help='Census version (e.g., 2024-07-01)')
    parser.add_argument('--query_path', type=str, default="/space/grp/rschwartz/rschwartz/nextflow_eval_pipeline/work/60/ed91620a1659a23ba5b68f8c21028c/lim_Cingulate.obs.relabel.tsv")
    parser.add_argument('--ref_name', type=str, default="Dissection:_Anterior_cingulate_cortex_ACC.h5ad") #nargs ="+")
    parser.add_argument('--ref_keys', type=str, nargs='+', default=["subclass", "class", "family"])
    parser.add_argument('--cutoff', type=float, default=0, help = "Cutoff threshold for positive classification")
    parser.add_argument('--probs', type=str, default="/space/grp/rschwartz/rschwartz/nextflow_eval_pipeline/work/4a/164f6559104b8e872e88bc617411a2/probs/lim_Cingulate_processed_Dissection:_Anterior_cingulate_cortex_ACC.prob.df.tsv")
    parser.add_argument('--mapping_file', type=str, default="/space/grp/rschwartz/rschwartz/nextflow_eval_pipeline/meta/census_map_human.tsv")
    parser.add_argument('--ref_region_mapping', type=str)

    
    if __name__ == "__main__":
        known_args, _ = parser.parse_known_args()
        return known_args

def get_unique_value(df, column, default=None):
    return df[column].unique()[0] if column in df.columns else default

    
def main():
    SEED = 42
    random.seed(SEED)         # For `random`
    np.random.seed(SEED)      # For `numpy`
    # For `torch`'
    scvi.settings.seed = SEED # For `scvi`
    # Parse command line arguments
    args = parse_arguments()
    query_path = args.query_path
    ref_name = args.ref_name
    ref_keys = args.ref_keys
    cutoff = args.cutoff
    #query_region = args.query_region.replace("_", " ")
    ref_region_mapping = args.ref_region_mapping
    #disease = args.query_disease
    #sex = args.query_sex
    #dev_stage = args.query_dev_stage
    
     
    # Load data
    ref_region_mapping = yaml.load(open(ref_region_mapping), Loader=yaml.FullLoader)
    ref_region=ref_region_mapping[ref_name]
    
    prob_df = pd.read_csv(args.probs, sep="\t")
    mapping_df = pd.read_csv(args.mapping_file, sep="\t")
    query_name = os.path.basename(query_path).replace(".obs.relabel.tsv", "")
    query = pd.read_csv(query_path, sep="\t")
    #for factor in factors:
    

    query_region = get_unique_value(query, 'region')
    disease = get_unique_value(query, 'disease')
    sex = get_unique_value(query, 'sex')
    dev_stage = get_unique_value(query, 'dev_stage')
    treatment = get_unique_value(query, 'treatment')
    genotype = get_unique_value(query, 'genotype')
    strain = get_unique_value(query, 'strain')
    age = get_unique_value(query, 'age')


    pr_metrics = pr_analysis(prob_df, query, key=ref_keys[0])
            plt.figure(figsize=(7, 8))
    
    class_labels = prob_df.columns.values
    # make colors for each class
    colors = sns.color_palette("husl", len(class_labels))
    for i, color in zip(range(len(class_labels)), colors):
        display = PrecisionRecallDisplay(recall=pr_metrics["recall"][i], 
                                         precision=pr_metrics["precision"][i], 
                                         average_precision=pr_metrics["average_precision"][i], 
                                         color=color)
        display.plot(ax=ax, name=f"Precision-recall for class {i}", color=color, despine=True)
        plt.annotate(f"AP={pr_metrics['average_precision'][i]:.2f}",
                       "optimal threshold"={pr_metrics['optimal_threshold'][i]:.2f},)
   
    # Add the legend for the iso-F1 curves
    handles, labels = display.ax_.get_legend_handles_labels()

    # Set the legend and the axes
    ax.legend(handles=handles, labels=labels, loc="best")
    ax.set_title("Precision-Recall curve for f{query_name} vs {ref_name}")

    # Show the plot
    plt.savefig(f"{query_name}_{ref_name}_pr_curve.png")
    plt.close()
    
    
    
    # Classify cells and evaluate
    query = classify_cells(query, ref_keys, cutoff=cutoff, probabilities=prob_df, mapping_df=mapping_df)
    outdir = "predicted_meta"
    os.makedirs(outdir, exist_ok=True)
    query.to_csv(os.path.join(outdir,f"{query_name}_{ref_name}.predictions.{cutoff}.tsv"), index=False, sep="\t")
    
    # map valid labels for given query granularity and evaluate
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
        accuracy = class_metrics[key]["accuracy"]
        label_accuracies = class_metrics[key]["label_accuracies"]
        print(label_accuracies)
        for label, metrics in classification_report.items():
            if label not in ["macro avg", "micro avg", "weighted avg", "accuracy"]:
                f1_data.append({
                    'query': query_name,
                    'reference': ref_name,
                    'label': label,
                    'f1_score': metrics['f1-score'],
                    'weighted_f1': classification_report.get('weighted avg', {}).get('f1-score', None),
                    'key': key,
                    'cutoff': cutoff,
                    'ref_region': ref_region,
                    'accuracy': accuracy,
                    'label_accuracy': label_accuracies.get(label, None),
                }
        )
                


    # Save F1 scores to a file
    df = pd.DataFrame(f1_data)
    
    fields_dict = {'disease': disease, 'sex': sex, 'dev_stage': dev_stage, 
                   'query_region': query_region, 
                   'treatment': treatment, 'genotype': genotype, 'strain': strain, 'age': age}   
    for field, value in fields_dict.items():
        df[field] = value if value is not None else np.nan


    outdir = "f1_results"
    os.makedirs(outdir, exist_ok=True)
    df.to_csv(os.path.join(outdir, f"{query_name}_{ref_name}.f1.scores.tsv"), sep="\t", index=False)
    
if __name__ == "__main__":
    main()
    

