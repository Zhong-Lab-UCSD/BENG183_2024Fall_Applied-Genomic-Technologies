# ScRNA Sequence - Seurat 
### By Nichole Mora, Ehsun Yazdani, Nicholas Hubbard
* [What is Seurat](#what-is-seurat)
* [Why Should You Use Seurat](#why-should-you-use-seurat)
* [What is Clustering](#what-is-clustering)
* [Case Study](#case-study)<be>
  * [Loading Data](#loading-data)
  * [Perform Clustering with Seurat Commands](#perform-clustering-with-seurat-commands)
* [Different Applications Performance Compared to Seurat](#different-applications-performance-compared-to-seurat)
  * [Scanpy](#scanpy)
  * [Monocle](#monocle)
  * [CellRanger](#cell-ranger)
  * [Harmony](#harmony)
* [Selected methods comparison](#selected-methods-comparison)




# What is Seurat 

Seurat is a powerful open-source R package specifically designed for the analysis and interpretation of single-cell RNA sequencing (scRNA-seq) data. Developed by the Satija Lab at the New York Genome Center, Seurat has rapidly become a prominent tool in the field of single-cell genomics. Its primary purpose is to facilitate the comprehensive exploration of cellular heterogeneity by providing an integrated workflow that encompasses data preprocessing, normalization, dimensionality reduction, clustering, and visualization.

One of the standout features of Seurat is its ability to handle complex datasets with high dimensionality, making it particularly well-suited for studies involving thousands to millions of individual cells. By leveraging advanced statistical methods and algorithms, Seurat enables researchers to uncover subtle differences in gene expression that may correspond to distinct cell types or states within a heterogeneous population.


## Why Should You Use Seurat

Seurat is widely favored in the single-cell genomics community for several compelling reasons. Firstly, it offers a comprehensive and streamlined workflow that covers all essential steps of single-cell data analysis. From initial quality control and data normalization to more sophisticated analyses like clustering and differential expression, Seurat provides a cohesive set of tools that simplify the entire analytical process.

Another significant advantage of Seurat is its scalability. As single-cell studies continue to grow in size and complexity, the ability to efficiently process large datasets becomes crucial. Seurat is designed to handle extensive datasets without compromising on performance, making it an ideal choice for large-scale studies that demand robust computational capabilities.

Community support is another area where Seurat excels. Being one of the most popular tools in its domain, Seurat has extensive documentation, numerous tutorials, and an active user community. This wide varierty of resources facilitates the learning curve for new users and provides ample support for troubleshooting and optimizing analyses.

Furthermore, Seurat's flexibility allows it to integrate seamlessly with other bioinformatics tools and pipelines. This interoperability enables researchers to customize their analyses and incorporate additional methods as needed, enhancing the overall versatility of their workflows.


## What is Clustering

Clustering is a fundamental unsupervised machine learning technique used to group similar data points based on their inherent features. In the realm of single-cell RNA sequencing, clustering plays a role in identifying distinct cell populations within a heterogeneous sample. By analyzing gene expression profiles, clustering algorithms can discern groups of cells that exhibit similar patterns, which often correspond to specific cell types or functional states.

Clustering serves as a crucial step in unraveling the complexity of cellular ecosystems, allowing researchers to map out the diverse landscape of cell types present in a given tissue or organism. This, in turn, facilitates a deeper understanding of biological processes such as development, differentiation, and disease progression.

Seurat these clustering methods to ensure accurate and meaningful groupings of cells. The process begins with dimensionality reduction, where high-dimensional gene expression data is transformed into a lower-dimensional space using techniques like Principal Component Analysis (PCA) or Uniform Manifold Approximation and Projection (UMAP). This reduction simplifies the data while preserving its essential structure, making it more amenable to clustering.

Once the data is in a reduced dimensional space, Seurat constructs a nearest-neighbor graph that captures the relationships between cells based on their similarity. Clustering algorithms, such as the Louvain or Leiden methods, are applied to this graph to identify clusters of closely related cells. These clusters are subsequently visualized in the reduced dimensional space, often with distinct colors representing different groups, allowing for a intuitive interpretation of the cellular landscape.

Evaluating the quality of clusters is an essential aspect of the analysis. Researchers assess whether the identified clusters correspond to known cell types or reveal novel populations by examining marker genes—genes that are uniquely expressed in each cluster. Additionally, the stability of clusters is evaluated by testing their consistency across different subsamples or varying analysis parameters, ensuring that the groupings are robust and biologically meaningful.


# Case Study
Now that we know what Seurat is and what it is used for, we will demonstrate how to use Seurat to perform clustering and the result of the clustering through a visualization. 

The following Case Study focuses on mouse sc-RNA sequence data and demonstrates clustering.


## Loading Data
To begin analyzing data to cluster, we must first load the data and create a Seurat object. 

In this study, we are given the expression data, metadata, PCA coordinates, and t-SNE coordinates which we will load in the following code.

```bash
#load the expression matrix
expression_data <- read.table("~/SCP404/expression/outputs_5ca76079328cee0c8dad60c0_66021478-5f80-4c3e-80b4-c2fbca32fce6/mouse_E18_nuclei_analysis.scp.expr.txt", 
                              header = TRUE, 
                              row.names = 1, 
                              sep = "\t")

#create Seurat object
seurat_obj <- CreateSeuratObject(counts = expression_data)


#load the data
metadata <- read.table("~/SCP404/metadata/outputs_5ca76079328cee0c8dad60c0_66021478-5f80-4c3e-80b4-c2fbca32fce6/mouse_E18_nuclei_analysis.scp.metadata.txt", 
                       header = TRUE, 
                       row.names = 1, 
                       sep = "\t")
pca_coords <- read.table("~/SCP404/cluster/outputs_5ca76079328cee0c8dad60c0_66021478-5f80-4c3e-80b4-c2fbca32fce6/mouse_E18_nuclei_analysis.scp.X_diffmap_pca.coords.txt", 
                         header = TRUE, 
                         row.names = 1, 
                         sep = "\t")

tsne_coords <- read.table("~/SCP404/cluster/outputs_5ca76079328cee0c8dad60c0_66021478-5f80-4c3e-80b4-c2fbca32fce6/mouse_E18_nuclei_analysis.scp.X_tsne.coords.txt", 
                          header = TRUE, 
                          row.names = 1, 
                          sep = "\t")'''
```
## Perform Clustering With Seurat Commands
In order to begin clustering, it is important to understand the commands required to perform the clustering process and create a visualization of the results. Below are all the Seurat commands with a description of each command you can use with the Seurat package.

```bash
NormalizeData
```
* Normalizes the gene expression data for each cell to consider differences in sequencing depth or library size to make it comparable across cells.
* Output: Normalized expression values which demonstrates the differences due to biological variation.

```bash
FindVariableFeatures
```
* Identifies genes with high variability across cells by looking at ranked genes.
* Output: A list of genes with high variability which is useful for dimensionality reduction and clustering.
```bash
ScaleData 
```
* Centers and scales the expression data for each gene in order to make it simple to compare across genes
* Output: Scaled data showing that all genes contribute equally to PCA.
```bash
RunPCA
```
* Performs principal component analysis to reduce the dimensionality of the dataset. PCA looks at axes that capture the maximum variance in data. 
* Output: A reduced dimensional dataset with cells represented in few dimensions, useful for clustering and visualization.
```bash
RunTSNE
```
* Reduction method for visualizing high-dimensional data into a two-dimensional or three-dimensional space for visualization.
* Output: a 2D or 3D plot demonstrating cells grouped depending on their similarity in gene expression. 
```bash
DimPlot
```
* Plots the cells in a 2D space and groups them by metadata column and color codes cells based on the metadata provided.
* Output: A visualization of the clustering.


### The Code to Perform Clustering
```bash
seurat_obj <- FindNeighbors(seurat_obj, dims = 1:10) 
seurat_obj <- FindClusters(seurat_obj, resolution = 0.5) 

seurat_obj <- NormalizeData(seurat_obj)

seurat_obj <- FindVariableFeatures(seurat_obj)

seurat_obj <- ScaleData(seurat_obj)

#required before t-SNE
seurat_obj <- RunPCA(seurat_obj, verbose = FALSE)

#run t-SNE
seurat_obj <- RunTSNE(seurat_obj, dims = 1:10) 

DimPlot(seurat_obj, reduction = "tsne", group.by = "seurat_clusters")
```
## Clustering Plot
The image below represents the clustering visualization of cells grouped depending on their similarity in gene expression of the mouse single-cell RNA sequencing data for this case study.
#### Key Aspects:
* Color code is used to indicate cluster identity.
* Clusters can reveal different cell types, subtypes within a cell type, and responses to certain treatments or conditions.
* Clustering plot helps us understand the diversity of cell populations and identify patterns in the data.
![](https://github.com/nmora2/ScRNA-Seq---Seurat/blob/c5fd38651e2a375fc65739eff41a0a31cdb1584a/Screenshot%202024-12-07%20at%204.59.19%20PM.png)


# Different Applications Performance Compared to Seurat
Now that we have looked at Seurat itself, we can see how it compares to alternative software and even how to optimize a workflow by combining these different choices.

## Scanpy
Scanpy is a toolkit developed to analyze single-cell expression data. It is designed for preprocessing, visualization, clustering, and more. Built with sparse matrix support and optimized algorithms, this alternative shows improved scalability compared to Seurat. It is also ideal for Python users due to its streamlined integration with other Python tools. Scanpy is often more challenging for new users due to its modular structure and lack of documentation. Considering all those advantages and disadvantages, it is optimal to use Scanpy for its quick preprocessing, and then switch to Seurat for improved visualization and heightened analysis.

## Monocle
Another similar tool to Seurat and Scanpy is Monocle. Built on the R platform, this software boasts support for the analysis of pseudotime and emphasizes trajectory analysis. Given its focus on these two forms of analysis, Monocle doesn't have as extensive capabilities for comprehensive clustering or integration, nor is it optimized for larger datasets. Thus, it makes sense to use Monocle following Scanpy or Seurat's preprocessing and clustering when interested in pseudotime and trajectory analysis.

## Cell Ranger
Cell Ranger is a standalone software developed by 10x Genomics for a multitude of purposes, including barcode processing, single-cell gene counting, sample demultiplexing, among many others. Its greatest advantages are apparent when used in conjunction with other hardware produced by the company. Due to its focus on end-to-end processing, Cell Ranger doesn't allow for significant user input during intermediary steps or parameter switching. Users will often generate a gene-cell matrix using Cell Ranger, and then use Seurat to filter out low-quality cells and genes.

## Harmony
Harmony is a toolkit developed for the R platform that focuses on batch correction and dataset integration  while avoiding overcorrecting and erasing small variations. It is capable of handling large datasets given its utilization of parallel processing and sparse matrix support, similar to Scanpy. On the other hand, Harmony is difficult for new users when used as a standalone software implementation and has limited visualization support. Thus, it is most often used to correct batch effects before further analysis can be performed with Seurat.







# Selected Methods Comparison
<table>
 <tbody>
    <tr>
        <th>Toolkit</td>
        <th>Advantages</td>
        <th>Disadvantages</td>
        <th>Workflow</td>
    </tr>
    <tr>
        <td>Seurat</td>
        <td><ul><li>Dimensionality Reeduction</li><li>Clustering</li><li>Visualization</li><li>Integration</li></li></ul></td>
        <td><ul><li>Struggles with differentiation of similar cell types</li><li>Lack of capability for predicting rare cell populations</li></ul></td>
        <td>N/A</td>
    </tr>
    <tr>
        <td>Scanpy</td>
        <td><ul><li>Improved Scalability</li><li>Python integration</li><li>Conservation of machine memory</li></ul></td>
        <td><ul><li>More challenging for new users</li><li>Less emphasis on visualization</li><li>Less documentation</li></ul></td>
        <td>Use scanpy for preprocessor for its speed and memory capability, then switch to Seurat for visuaization and heightened analysis</td>
    </tr>
    <tr>
        <td>Monocle</td>
        <td><ul><li>Dynamic gene expression anaysis</li><li>Support for multi-modal datasets</li><li>Distinct visualizations for trajectories and pseudotime</li></ul></td>
        <td><ul><li>No focus on comprehensive clustering or integration</li><li>Not optimized for large datasets</li><li>LLess versatility in plotting options</li></ul></td>
        <td>Monocle is used following Seurat's preprocessing for pseudotime and trajectory analysis</td>
    </tr>
    <tr>
        <td>Cell Ranger</td>
        <td><ul><li>Ease of use</li><li>Optimized for high performance, esepecially in conjunction with othere 10x Genomics products</li><li>Spatial transcriptomics</li></ul></td>
        <td><ul><li>Limited flexibility in modifying intermeediary steps</li><li>Requires significant computational resources</li></ul></td>
        <td>Generate a gene-cell matrix with Cell Ranger, then use Seurat to filter out low-quality cells and genes.</td>
    </tr>
    <tr>
        <td>Harmony</td>
        <td><ul><li>Batch correection and dataset integration</li><li>Avoids overcorrecting and erasing small biological variation</li><li>Uses parallel processing and sparse matrices for large dataset support</li></ul></td>
        <td><ul><li>Difficult for new users</li><li>Limited visualization support</li></ul></td>
        <td>Harmony is useed to correct batch effects before further analysis and/or visualization</td>
    </tr>
 </tbody>
</table>

















# References
* Seurat Official Website: https://satijalab.org/seurat
* Seurat GitHub Repository: https://github.com/satijalab/seurat
* Seurat Tutorials: https://satijalab.org/seurat/articles/pbmc3k_tutorial.html
* Comprehensive Analysis Workflow: https://satijalab.org/seurat/articles/essential_features.html
* Seurat Integration Tutorial: https://satijalab.org/seurat/articles/integration_introduction.html
* Mouse Clustering Case study: https://singlecell.broadinstitute.org/single_cell/study/SCP404/demo-cellranger-sccloud

Academic Publications:
* Original Seurat Paper (2015):
Satija, R., et al. (2015). "Spatial reconstruction of single-cell gene expression data." Nature Biotechnology.
https://doi.org/10.1038/nbt.3192
* Seurat v3 Paper (2019):
Stuart, T., et al. (2019). "Comprehensive integration of single-cell data." Cell.
https://doi.org/10.1016/j.cell.2019.05.031
* Seurat v4 Paper (2021):
Hao, Y., et al. (2021). "Integrated analysis of multimodal single-cell data." Cell.
https://doi.org/10.1016/j.cell.2021.04.048
