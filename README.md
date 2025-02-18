# nextflow_eval_pipeline

## Table of Contents
- [Description](#description)
- [Methods](#methods)
- [Reference data](#reference)
- [Test data](#test)
- [Cell type taxonomy](#taxonomy)
- [Requirements](#requirements)
- [Installation](#installation)
- [Usage](#usage)
- [Input](#input)
- [Output](#output)
- [Workflow Description](#workflow_description)

---

# Description

Single-cell RNA sequencing (scRNA-seq) provides crucial insights into cell-type-specific gene expression, particularly in the context of disease. However, meta-analysis of scRNA-seq datasets remains challenging due to inconsistent and often absent cell-type annotations across publicly available repositories, such as the Gene Expression Omnibus (GEO) [cite GO]. Automated re-annotation of author cell type labels using classifiers trained on high-quality “reference” data  is a common approach to single-cell meta-analysis, enabling investigators to group cell types across datasets. [Pasquini][Lotfollahi].

The selection of an appropriate cell type classifier and related parameters is the subject of much benchmarking, particularly for brain datasets with highly differentiated cell types [Pasquini][Lotfollahi][Leucken][Lopez][Hao]. However, typical benchmarking of automated cell-type annotation (also known as "label transfer") evaluates performance on a hold-out set of the same data or a closely related dataset of a different modality. There is little investigation of how well a joint embedding of reference and test dataset captures biological information "in the wild", that is on public datasets comprising diverse brain regions, conditions and drug treatments. Furthermore, there exist several BRAIN initiative single-cell datasets which are frequently used as “references” for annotation of human and mouse neocortex [Jorstad][Gabitto][Bakken][BRAIN initiative]. Benchmarking of both reference datasets and label transfer methods on public data instead of holdhout sets is therefore essential to our re-annotation pipeline. Here we present an initial evaluation of two prominent cell type annotation strategies on scRNAseq human neocortex: a Random Forest classifier trained on cell embeddings using single-cell variational inference (SCVI) [Lopez], and A Gaussian kernel trained on PCA projection of query onto reference datasets using the Seurat package [Lopez].

Reference data comprises 2 studies, 10 brain regions, 12 individual dissections, and 3 levels of cell type granularity. Preliminary test data includes 6 studies comprising 3 brain regions, 7 diseases and 3 developmental stages. We separated each study into its respective samples to extract condition, sex and developmental stage factors. We assessed the impact of choice of reference, annotation method, reference subsampling, confidence filtering, and experimental factors of the test data on annotation performance. We first performed ANOVA on weighted F1 scores (aggregated across all cell types) to assess influence of these parameters, and found that thresholding had the most influence on F1 scores, which decrease as threshold is increased. Disease, development stage, and sex did not impact performance. We thn then fit a linear model to weighted F1 scores without thresholding and regressed out the effects of individual studies. Re-annotation method was found to interact with reference dataset, with Seurat outperforming SCVI for individual reference datasets, and SCVI outperforming Seurat when reference data is aggregated. This is to be expected, as Seurat PCA projection does not account for batch effects in the reference, while it is an explicit feature of SCVI. Finally, we found that SCVI scales more efficiently to large datasets than Seurat with respect to computation time and CPU usage.


## Methods:
A Gaussian kernel trained on Dual PCA projection of reference and query datasets using the Seurat package [cite Seurat]
A random forest classifier trained on latent embedding for reference and query data using a pre-trained single-cell variational autoencoder (scvi) model [Lopez]

## Reference data:


## Test data:

A manually curated "ground truth" reference, based on the BRAIN initiative cell type taxonomy, was used to assess model performance:

## Cell type taxonomy:




# References
Puntambekar, S., et al. "Cell-level metadata are indispensable for documenting single-cell sequencing datasets." PLoS Biol., 2021.
Pasquini, G., et al. "Automated methods for cell type annotation on scRNA-seq data." Computational and Structural Biotechnology Journal, 2021.
Lotfollahi, M., et al. "The Future of Rapid and Automated Single-Cell Data Analysis Using Reference Mapping." Cell, 2024.
Jorstad, N.L., et al. "Transcriptomic Cytoarchitecture Reveals Principles of Human Neocortex Organization." Science, 2023.
Abdulla, S., et al. "CZ CELL×GENE Discover: A Single-Cell Data Platform." bioRxiv, 2023.
Lopez, R., et al. "Deep Generative Modeling for Single-Cell Transcriptomics." Nature Methods, 2018.





![workflow DAG](dag.png)
