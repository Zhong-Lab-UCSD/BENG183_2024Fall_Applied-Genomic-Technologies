# ScRNA Sequence - Seurat 
### By Nichole Mora, Ehsun Yazdani, Nicholas Hubbard
* [What is Seurat](#what-is-seurat)
* [Why Should You Use Seurat](#why-should-you-use-seurat)
* [What is Clustering](#what-is-clustering)
* [Case Study](#case-study)<be>
  * [Loading Data](#loading-data)
  * [Perform Clustering with Seurat Commands](#perform-clustering-with-seurat-commands)
* [Different Applications Performance Compared to Seurat](#different-applications-performance-compared-to-seurat)
  * [Insert Different Application1](#insert-different-application1)
  * [Insert Different Application2](#insert-different-application2)
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
The Case Study focuses on mouse sc-RNA sequence data and demonstrates clustering using SCP's analysis pipeline.


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


## Insert Different Application1


## Insert Different Application2








# Selected Methods Comparison
## AN EXAMPLE FROM TEMPLATE but use this to compare applications and include seurat
<table>
 <tbody>
    <tr>
        <th>Method</td>
        <th>Targets</td>
        <th>Resolution</td>
        <th>Notes</td>
    </tr>
    <tr>
        <td>3C <a href="http://refhub.elsevier.com/S2001-0370(17)30093-4/rf0535">[3]</a></td>
        <td>one-vs-one</td>
        <td>~1–10 kb<br></td>
        <td><ul><li>Sequence of bait locus must be known</li><li>Easy data analysis</li><li>Low throughput</li></ul></td>
    </tr>
    <tr>
    <td>4C <a href="http://refhub.elsevier.com/S2001-0370(17)30093-4/rf0545">[4]</a></td>
    <td>one-vs-all</td>
    <td>~2 kb</td>
    <td><ul><li>Sequence of bait locus must be known</li><li>Detects novel contacts</li><li>Long-range contacts</li></ul></td>
    </tr>
    <tr>
    <td>5C <a href="http://refhub.elsevier.com/S2001-0370(17)30093-4/rf0550">[5]</a></td>
    <td>many-vs-many</td>
    <td>~1 kb</td>
    <td><ul><li>High dynamic range</li><li>Complete contact map of a locus</li><li>3C with ligation-mediated amplification (LMA) of a ‘carbon copy’ library of oligos designed across restriction fragment junctions of interest
3C</li></ul></td>
    </tr>
    <tr>
    <td>Hi-C <a href="http://refhub.elsevier.com/S2001-0370(17)30093-4/rf0300">[6]</a></td>
    <td>all-vs-all</td>
    <td>0.1–1 Mb</td>
    <td><ul><li>Genome-wide nucleosome core positioning</li><li>Relative low resolution</li><li>High cost</li></ul></td>
    </tr>
    <tr>
    <td>ChIA-PET <a href="http://refhub.elsevier.com/S0168-9525(15)00063-3/sbref1405">[7]</a></td>
    <td>Interaction of whole genome mediated by protein</td>
    <td>Depends on read depth and the size of the genome region bound by the protein of interest</td>
    <td><ul><li>Lower noise with ChIP</li><li>Biased method since selected protein</li></ul></td>
    </tr>
 </tbody>
</table>

















# References
* Seurat Official Website: https://satijalab.org/seurat
* Seurat GitHub Repository: https://github.com/satijalab/seurat
* Seurat Tutorials: https://satijalab.org/seurat/articles/pbmc3k_tutorial.html
* Comprehensive Analysis Workflow: https://satijalab.org/seurat/articles/essential_features.html
* Seurat Integration Tutorial: https://satijalab.org/seurat/articles/integration_introduction.html

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
