# scGPT: Generative AI Analysis of Single Cell Sequencing Data
### Safa Saeed, Honieh Hemati, and Saara Kriplani


* [Background](#background)
* [What is scGPT?](#what-is-scGPT)
* [How does scGPT work?](#How-does-scGPT-work?)
* [Application I: Gene Perturbations](#applications-i-gene-perturbations)
* [Application II: Multi-Batch Integration](#applications-ii-multi-batch-integration)
* [Tutorial - Cell Annotation App](#tutorial-cell-annotation-app)
* [Sources](#sources)

---

## Background
Single cell RNA sequencing (scRNA-seq) is a technique used to extract and analyze RNA from an individual cell. Usually, with bulk RNA sequencing, we take RNA from a large population of cells and average gene expression values, which can lose information about the heterogeneity of individual cells. scRNA-seq provides a solution to this issue, allowing us to focus on the functions and expression patterns of distinct cell types.

The main output of scRNA-seq is a gene expression matrix, where we get values corresponding to the expression levels of specific genes within the cell of interest. We can use this expression matrix to generate an expression profile for the cell, allowing different cells to be compared and clustered together. 

---

## What is scGPT?

Single cell Generative Pre-Training Transformer, or scGPT, is a tool introduced by the University of Toronto's Wang Lab published earlier this year. This tool is an exciting and significant advancement in single cell analysis. 

Single-cell sequencing produces a wide variety of information that must be analyzed effectively to draw conclusions about things like gene expression and gene interactions. Current ML methods of analyzing single cell data are specific to a certain task, and are thus limited in scope. In comparison, scGPT uses a “Pretraining universally, Fine-tuning on Demand” approach which enables a diverse range of downstream tasks. 

In other words, scGPT is an all-in-one tool that can be easily customized for any single-cell omics use. 

#### scGPT can be fine-tuned for a variety of downstream analyses, including...
 
- Clustering
- Batch Correction
- Cell Type Annotation 
- Gene network inference
- Multi-omics integration 
- Perturbation prediction


## How does scGPT work?
 
scGPT utilizes advanced techniques in machine learning to analyze and predict gene interactions, leveraging a vast and diverse dataset of single-cell RNA sequencing (scRNA-seq) data. To begin, researchers gathered data on 33 million cells from a cell atlas to pretain the model, including information on cell type, gene interactions, and more. 

The three main inputs taken from the scRNA-seq data are 

- Gene name tokens
- Gene expression value
- Condition tokens (e.g. whether the cell came from a cancer patient)

These values are mapped into representative numerical vectors which are calculated by PyTorch embedding layers. Next, there is a self-attention transformer that takes in embedding vectors and calculates attention scores, which are higher for genes that are more likely to interact.

After the pretraining phase, the embedding and attention process can be repeated on a smaller, focused dataset for fine-tuning for specific purposes. The final attention-weighted embedding vectors for each gene now contain information that can be used for downstream analysis

 The core feature of scGPT is the **self-attention transformer**, which incorporates context, i.e. how is one gene influenced by the other genes.



## Application I: Gene Perturbations

One application of scGPT is the prediction of gene perturbation effects. In other words, it can determine how altering a gene can impact expression levels. 

Based on Section A (top left) in the figure below, fine-tuned scGPT (seen in pink) performs better than two leading methods, GEARS and Linear, at predicting gene expression changes compared to actual experimental results. 

![Figure2](./Figures/gene_perturbations_cropped.png)

In Section B (top right), we can see change in expression values compared to control levels for two perturbed genes, DAD1 and KCTD16. Evidently, the value predicted by scGPT (in pink) falls within the true (in blue) expression values for each gene. 

In Section D and E, we can see clusters labeled by perturbed genes, and the predicted gene expression profiles cluster distinctly based on which gene is modified. 


## Application II: Multi-Batch Integration 

Another useful application is multi batch and multi omic integration. It is often difficult to draw conclusions from data collected from different cell batches or different forms of expression data, and scGPT excels at identifying cell types even with this added complexity.

The datasets in this figure are multi-batch and the bottom two rows are multi-omic so like CHIp seq and ATAC seq in addition to RNAseq. We can see that scGPT outperforms other models like Seurat and scGLUE, with the highest AvgBIO score for all datasets. Also, scGPT is better at finding distinct clusters in this multi-batch and -omic data. Looking at the second row, you can see in the other models this brown cluster is not distinctly identified wheras in scGPT it is as clearly separated from other cell types. 




## Tutorial - Cell Annotation App 

## Sources 