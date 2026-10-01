
#!/user/bin/python3

import os
import signal
signal.signal(signal.SIGPIPE, signal.SIG_DFL)
import numpy as np
import pandas as pd
import scvi
from utils import *
import argparse
import yaml

# Function to parse command line arguments
def parse_arguments():
    parser = argparse.ArgumentParser(description="Classify cells given 1 ref and 1 query")
  #  parser.add_argument('--census_version', type=str, default='2024-07-01', help='Census version (e.g., 2024-07-01)')
    parser.add_argument('--query_path', type=str, default="/space/grp/rschwartz/rschwartz/nextflow_eval_pipeline/2025-01-30/mus_musculus/100/dataset_id/SCT/gap_false/ref_50_cutoff_0/scvi/GSE152715.2/whole_cortex/GSE152715.2_1052248_GSM4624685/probs/GSE152715.2_1052248_GSM4624685.obs.relabel.tsv")
    parser.add_argument('--ref_name', type=str, default="whole_cortex") #nargs ="+")
    parser.add_argument('--ref_keys', type=str, nargs='+', default=["subclass", "class", "family","global"])
    parser.add_argument('--cutoff', type=float, default=0, help = "Cutoff threshold for positive classification")
    parser.add_argument('--probs', type=str, default="/space/grp/rschwartz/rschwartz/nextflow_eval_pipeline/2025-01-30/mus_musculus/100/dataset_id/SCT/gap_false/ref_50_cutoff_0/scvi/GSE152715.2/whole_cortex/GSE152715.2_1052248_GSM4624685/probs/probs/GSE152715.2_1052248_GSM4624685_whole_cortex.prob.df.tsv")
    parser.add_argument('--mapping_file', type=str, default="/space/grp/rschwartz/rschwartz/nextflow_eval_pipeline/meta/census_map_human.tsv")
    parser.add_argument('--ref_region_mapping', type=str, default="")
    parser.add_argument('--study_name', type=str, default="GSE152715.2")
    parser.add_argument('--ref_counts', type=str, default=None, help="TSV with ref label counts per key (columns: key, label, ref_support)")
    parser.add_argument('--use_gap', action='store_true', help="Use gap analysis for classification")
    parser.add_argument('--method', type=str, default=None, help="Classification method (scvi_rf, scvi_knn, seurat)")
    
    if __name__ == "__main__":
        known_args, _ = parser.parse_known_args()
        return known_args

def get_unique_value(df, column, default=None):
    if column in df.columns:
    # check how many unique values there are
        if len(df[column].unique()) == 1:
            return df[column].unique()[0] 
        else:
            return default


def write_unlabeled_unknown_rates(query, is_unlabeled, key, query_name, study_name, ref_name, method, cutoff, use_gap):
    """Per-cell predictions and QC for labeled vs author_unlabeled cells, and the fraction of each called "unknown"."""
    outdir = "unlabeled_qc"
    os.makedirs(outdir, exist_ok=True)
    qc_cols = ["sample_id", "cell_id", "cell_type", key, f"predicted_{key}", "confidence",
               "total_counts", "n_genes_by_counts", "nCount_RNA", "nFeature_RNA", "pct_counts_mito",
               "umi_outlier", "genes_outlier", "counts_outlier", "predicted_doublet", "total_outlier"]
    cells = query[[c for c in qc_cols if c in query.columns]].copy()
    cells.insert(0, "author_unlabeled", is_unlabeled.values)
    cells["predicted_unknown"] = query[f"predicted_{key}"].astype(str) == "unknown"
    cells.to_csv(os.path.join(outdir, f"{query_name}_{ref_name}.unlabeled_cells.{cutoff}.tsv.gz"), sep="\t", index=False, compression="gzip")

    records = []
    for group, grp in cells.groupby("author_unlabeled"):
        records.append({
            'query': query_name, 'study': study_name, 'reference': ref_name, 'method': method,
            'cutoff': cutoff, 'use_gap': use_gap, 'key': key,
            'group': "author_unlabeled" if group else "labeled",
            'n_cells': len(grp),
            'n_unknown': int(grp["predicted_unknown"].sum()),
            'frac_unknown': grp["predicted_unknown"].mean(),
            'frac_mad_outlier': grp["total_outlier"].astype(str).str.lower().eq("true").mean() if "total_outlier" in grp else np.nan,
            'median_total_counts': grp["total_counts"].median() if "total_counts" in grp else np.nan,
            'median_n_genes': grp["n_genes_by_counts"].median() if "n_genes_by_counts" in grp else np.nan,
        })
    pd.DataFrame(records).to_csv(os.path.join(outdir, f"{query_name}_{ref_name}.unlabeled_unknown_rate.{cutoff}.tsv.gz"), sep="\t", index=False, compression="gzip")


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
    ref_region_mapping = args.ref_region_mapping
    study_name = args.study_name
    if args.use_gap:
        use_gap = True
    else:
        use_gap = False 
    
    # Load ref_counts for ref_support lookup (key -> label -> count)
    if args.ref_counts:
        ref_counts_df = pd.read_csv(args.ref_counts, sep="\t")
        ref_counts_lookup = {
            key: grp.set_index('label')['ref_support'].to_dict()
            for key, grp in ref_counts_df.groupby('key')
        }
    else:
        ref_counts_lookup = {}
    flat_ref_counts = flatten_ref_counts(ref_counts_lookup, ref_keys)

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

    os.makedirs("pr_curves", exist_ok=True) 
    
    # Classify cells and evaluate
    query = classify_cells(query=query, ref_keys=ref_keys, cutoff=cutoff, probabilities=prob_df, mapping_df=mapping_df, use_gap=use_gap)

    # author_unlabeled cells are scored separately and kept out of all metrics and predictions below
    if "author_unlabeled" in query.columns:
        is_unlabeled = query["author_unlabeled"].astype(str).str.lower() == "true"
        write_unlabeled_unknown_rates(query, is_unlabeled, ref_keys[0], query_name, study_name, ref_name, args.method, cutoff, use_gap)
        query = query[~is_unlabeled].reset_index(drop=True)

    outdir = os.path.join("predicted_meta")
    os.makedirs(outdir, exist_ok=True)

    # map valid labels for given query granularity and evaluate
    query = map_valid_labels(query, ref_keys, mapping_df)  
    class_metrics = evaluate_sample_predictions(query, ref_keys, mapping_df, ref_counts_lookup=ref_counts_lookup)
    
    query.to_csv(os.path.join(outdir,f"{query_name}_{ref_name}.predictions.{cutoff}.tsv.gz"), index=False, sep="\t", compression='gzip')

    # Plot confusion matrices
    for key in ref_keys:
        outdir = os.path.join("confusion")
        os.makedirs(outdir, exist_ok=True)
        plot_confusion_matrix(query_name, ref_name, key, class_metrics[key]["confusion"], output_dir=outdir)

    # get total cell count for the sample
    total_cell_count = query.shape[0]
    
    # Collect per-label classification metrics
    metrics_records = []
    for key in ref_keys:
        label_metrics = class_metrics[key]["label_metrics"]
        weighted_metrics = class_metrics[key]["weighted_metrics"]
        macro_metrics = class_metrics[key]["macro_metrics"]
        nmi = class_metrics[key]["nmi"]
        ari = class_metrics[key]["ari"]
        overall_accuracy = class_metrics[key]["overall_accuracy"]
        predicted_counts = query[f"predicted_{key}"].value_counts().to_dict()
            # get total cell count for the sample
        for label, metrics in label_metrics.items():
            #if label not in ["macro avg", "micro avg", "weighted avg", "accuracy"]:
                metrics_records.append({
                    'query': query_name,
                    'study': study_name,
                    'reference': ref_name,
                    'label': label,
                    'f1_score': metrics['f1_score'],
                    'accuracy': metrics['accuracy'],
                    'precision': metrics['precision'],
                    'recall': metrics['recall'],
                    'support': metrics['support'],
                    'predicted_support': predicted_counts.get(label, 0),
                    'ref_support': flat_ref_counts.get(label, 0),
                    'weighted_f1': weighted_metrics.get('f1_score', None),
                    'weighted_precision': weighted_metrics.get('precision', None),
                    'weighted_recall': weighted_metrics.get('recall', None),
                    # add macro averages
                    'macro_f1': macro_metrics.get('f1_score', None),
                    'macro_precision': macro_metrics.get('precision', None),
                    'macro_recall': macro_metrics.get('recall', None),
                    # add micro averages
                    'micro_f1': class_metrics[key]["micro_metrics"].get('f1_score', None),
                    'micro_precision': class_metrics[key]["micro_metrics"].get('precision', None),
                    'micro_recall': class_metrics[key]["micro_metrics"].get('recall', None),
                    'nmi': nmi,
                    'ari': ari,
                    'overall_accuracy': overall_accuracy,
                    'key': key,
                    'cutoff': cutoff,
                    'ref_region': ref_region,
                    'total_cell_count': total_cell_count
                }
        )

    # Save classification metrics to a file
    df = pd.DataFrame(metrics_records)
    
    fields_dict = {'disease': disease, 'sex': sex, 'dev_stage': dev_stage, 
                   'query_region': query_region, 
                   'treatment': treatment, 'genotype': genotype, 'strain': strain, 'age': age}   
    for field, value in fields_dict.items():
        df[field] = value if value is not None else np.nan


    outdir = "label_transfer_metrics"
    os.makedirs(outdir, exist_ok=True)
    df.to_csv(os.path.join(outdir, f"{query_name}_{ref_name}.summary.scores.tsv.gz"), sep="\t", index=False, compression='gzip')
    
if __name__ == "__main__":
    main()
    

