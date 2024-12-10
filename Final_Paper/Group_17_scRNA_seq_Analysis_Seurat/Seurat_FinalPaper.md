# Single-Cell RNA Sequencing Analysis Using [Seurat](https://satijalab.org/seurat/)

By: Joseph Hwang, Andrew Quach, Alisa Vu 

## **Table of Contents**
1. [Introduction](#Introduction)
2. [Learning Objective](#LearningObjective)
3. [Workflow](#Workflow)
4. [Real World Applications](#RealWorld)
5. [Strengths](#Strengths)
6. [Limitations](#Limitations)
7. [Conclusion](#Conclusion)

# **Introduction<a name="Introduction"></a>**

Imagine standing inches away from a large painting where all you see are thousands of individually-colored dots. At first glance, you might think the dots are random and chaotic. Stepping back, however, you realize the tiny dots actually form a coherent picture. This style of painting was created using "Pointillism", a technique partially developed by French Artist Georges Pierre Seurat. 

![point](https://principlearttalk.com/wp-content/uploads/2015/05/seurat-detail.jpg)

Now, imagine Pointilism but with a biological twist. Instead of painted dots, you have single cells, and instead of a brush, you have the bioinformatic tool called Seurat. Seurat is a software for single-cell RNA sequencing (scRNA-seq) analysis that, similar to the painter Georges Seurat, transforms individual data points into biological insights. Seurat gives researchers the ability to paint a clearer picture among complex biological data. With this connection between art and science, let's dive into why Seurat is so powerful in uncovering the hidden beauty and complexity of scRNA-seq data. 


**What is Single-Cell RNA Sequencing (scRNA-seq)?**

Single-cell RNA sequencing is a relatively new sequencing technology that was first discovered in 2009. It is widely used in research, because it allows scientists to analyze gene expression levels at the level of individual cells. In contrast to bulk RNA sequencing, which measures the average gene expression across a given sample of cells, scRNA-seq captures the expression of individual cells. This enables the discovery of cell heterogeneity and identification of distinct cell types, states, and functions of a given sample. As a result, this technology has many applications in research to understand tissues, developmental processes, and diseases.


**Why Seurat?**

Seurat is a very popular tool to use when working with scRNA-seq data, because it provides numerous utilities and an end-to-end workflow. Seurat can be used to do data preprocessing, dimensionality reduction, clustering, differential expression analysis, data integration, and visualization. As a result, it is a popular choice of tool when trying to extract biological meaningful insights from scRNA-seq data, especially when working with hetergeneous samples and high-dimensional data. While Seurat is an R package, there are alternatives for other programming languages, such as Scanpy for Python and Harmony for Python/R.

# **Learning Objective**<a name="LearningObjective"></a>

In this chapter, we will analyze scRNA-seq data of Peripheral Blood Mononuclear Cells (PMBC) from 10X Genomics. By analyzing single-cell data of PBMCs, we can reveal cell heterogeneity by identifying distinct immune cell types based on their gene expression profiles. Note that for our analysis, we are not running CellRanger on the raw FASTQ file (raw form of scRNA-seq data), but instead starting from the counts matrix. The counts matrix can be found [here](https://cf.10xgenomics.com/samples/cell/pbmc3k/pbmc3k_filtered_gene_bc_matrices.tar.gz). 


# **Workflow<a name="Workflow"></a>**

- Loading single-cell data into Seurat
- Quality Control and Filtering
- Normalization and Scaling
- Using dimensionality reduction techniques
- Clustering cells
- Visualizing and identifying cell type clusters

### 0. Setup

First, we'll need to install the Seurat R package. To install, we will need R version 4.0 or greater. Seurat is available on CRAN:
```
install.packages('Seurat')
library(Seurat)
```
Although not required, we also recommend installing the following packages to enhance speed, performance, and functionality. 
```
# For better speed and performance
setRepositories(ind = 1:3, addURLs = c('https://satijalab.r-universe.dev', 'https://bnprks.r-universe.dev/'))
install.packages(c("BPCells", "presto", "glmGamPoi"))
if (!requireNamespace("remotes", quietly = TRUE)) {
  install.packages("remotes")
}

# Added functionalities
install.packages('Signac')
remotes::install_github("satijalab/seurat-data", quiet = TRUE)
remotes::install_github("satijalab/azimuth", quiet = TRUE)
remotes::install_github("satijalab/seurat-wrappers", quiet = TRUE)
```


### 1. Loading the Data

Seurat stores data in a Seurat Object. We can load in our counts matrix data using `Read10X()`. This function returns a unique molecular identified (UMI) count matrix, where the values in the matrix represent the number of molecules for each gene (row) that are detected in each cell (column). Next, we use the count matrix to create a Seurat Object.
```
library(dplyr)
library(Seurat)
library(patchwork)

# Load the PBMC dataset
pbmc.data <- Read10X(data.dir = "/brahms/mollag/practice/filtered_gene_bc_matrices/hg19/")
# Initialize the Seurat object with the raw (non-normalized data).
pbmc <- CreateSeuratObject(counts = pbmc.data, project = "pbmc3k", min.cells = 3, min.features = 200)
pbmc
```
A small snippet of what the counts matrix looks like:
```
3 x 30 sparse Matrix of class "dgCMatrix"
                                                                   
CD3D  4 . 10 . . 1 2 3 1 . . 2 7 1 . . 1 3 . 2  3 . . . . . 3 4 1 5
TCL1A . .  . . . . . . 1 . . . . . . . . . . .  . 1 . . . . . . . .
MS4A1 . 6  . . . . . . 1 1 1 . . . . . . . . . 36 1 2 . . 2 . . . .
```
Note that `.` in the matrix represent 0s.

### 2. Quality Control and Data Filtering

Now we need to preprocess and filter out poor quality data. Generally, there are three levels of filtering: 
- **Filtering cells**: We filter out low-quality cells which often have very few genes expressed. In addition, cell doublets or multiplets may have very high number of genes expressed as well, which we filter as well. 
- **Filtering genes**: We filter out genes that are not expressed or expressed at very low levels.
- **Filtering cells with high mitrochondrial gene expression**: Low-quality, damaged, and dying cells often are indicated by high numbers of mitochondrial transcripts. Since we do not want to cluster our cells based on cell stress levels, we filter cells with high percentage of reads coming from mitochondrial genes. 

```
# Calculate and stash mitochondrial counts in Seurat Object
pbmc[["percent.mt"]] <- PercentageFeatureSet(pbmc, pattern = "^MT-")
```
We can visualize feature-feature relationships to get an understanding of our data using Seurat's `VlnPlot()` and `FeatureScatter()` functions. This can help us determine how we want to filter and preprocess our data. 

```
# QC metrics as a violin plot
VlnPlot(pbmc, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), ncol = 3)
```
![qc2-1](https://github.com/user-attachments/assets/8fa9b3e8-8578-41b0-8ff6-3963ec947cf0)

```
# feature-feature relationships 
plot1 <- FeatureScatter(pbmc, feature1 = "nCount_RNA", feature2 = "percent.mt")
plot2 <- FeatureScatter(pbmc, feature1 = "nCount_RNA", feature2 = "nFeature_RNA")
plot1 + plot2
```
![qc2-2](https://github.com/user-attachments/assets/e765260a-79d8-4d10-af87-82808ba70401)

In our case, we filter cells that have over 2,500 or less than 200 unique genes. In addition, we filter cells that have greater than 5% mitochondrial counts. 
```
pbmc <- subset(pbmc, subset = nFeature_RNA > 200 & nFeature_RNA < 2500 & percent.mt < 5)
```
### 3. Normalization and Scaling
Finally, we want to normalize the data so that we can compare gene expressions across different cells. One way of doing this is by employing global-scaling normalization or log normalization. This method normalizes the expression value of a given gene i in cell j using the formula:

![normalized_expression_formula_corrected](https://github.com/user-attachments/assets/d2362295-054c-4718-8f5a-ed77ab9aa968)

In Seurat, we can use `NormalizeData()` to normalize our values.
```
pbmc <- NormalizeData(pbmc, normalization.method = "LogNormalize", scale.factor = 10000)
```
In addition to normalization, we need to scale the data. This is a standard preprocessing step before applying dimensional reduction techniques. When we scale the data, we shift expression of each gene so that the mean expression across cells is 0 and variance is 1. This allows us to give equal weight for all genes without having highly-expressed genes dominating analyses. 
```
all.genes <- rownames(pbmc)
pbmc <- ScaleData(pbmc, features = all.genes)
```

### 4. Dimensionality Reduction
Now we want to determine how many principal components (PCs) we want to include in our analyses. Seurat reduces the dataset to a smaller set of PCs as a way to summarize key variations in the data while minimizing noise. Each PC can capture biological variations such as variations in cell types or states. PCs are derived from principal component analysis (PCA), which is a dimensionality reduction technique.
```
pbmc <- RunPCA(pbmc, features = VariableFeatures(object = pbmc))
```
We can then use an elbow plot, which displays a ranking of principal components based on the percentage of variance, to determine how many PCs to retain for downstream analysis. We want to identify the elbow in the plot, or the point where the curve flattens, indicating that any additional PCs will have a minimal effect to the analysis.
```
ElbowPlot(pbmc)
```
![elbow](https://github.com/user-attachments/assets/0cb2f7db-6407-425a-9946-0f543230fbf4)


In this example, the elbow exists between PC 1-10, suggesting that most true signal is captured here and therefore we only want to perform downstream analysis on the first 10 PCs.

### 5. Clustering Cells
Seurat uses `FindNeighbors()` and `FindClusters()` to group cells with similar expression patterns.
```
pbmc <- FindNeighbors(pbmc, dims = 1:10)
pbmc <- FindClusters(pbmc, resolution = 0.5)
```
The `FindNeighbors()` function constructs a K-nearest neighbor (KNN) graph based on the euclidean distance in PCA space and uses Jaccard similarity to refine edge weights between each pair of cells.
The `FindClusters()` function applies the Louvain algorithm, which is an algorithm that can maximize the modularity of the clustering and how well the graph is divided into clusters. The clusters are formed by grouping cells that have higher edge weights. The resolution parameter controls how the cells are divided:
- a higher resolution means more clusters with smaller, finer groups
- a lower resolution means fewer clusters with larger, broader groups


### 6. Visualization
Seurat allows us to use non-linear dimensional reduction techniques like UMAP to visualize the datasets. It can help us place similar cells together in low-dimensional space. 
```
pbmc <- RunUMAP(pbmc, dims = 1:10)
DimPlot(pbmc, reduction = "umap")
```

![umapplot-1](https://github.com/user-attachments/assets/d11b154c-4754-4792-a7b8-9b91e175dc9f)

Each dot in the graph represents a cell.

### 7. Identifying Marker Genes and Cell Type Annotation
Finally, we need to assign cell types to our clusters. One way to do this is to use marker genes that are known to be expressed in certain cell types. Conveniently, Seurat has a function called `FindAllMarkers()` that identifies positive and negative markers for all clusters. For example, MS4A1 is a marker for B cells, CD8A is a marker for CD8+ T, etc.
```
pbmc.markers <- FindAllMarkers(pbmc, only.pos = TRUE)
pbmc.markers %>%
    group_by(cluster) %>%
    dplyr::filter(avg_log2FC > 1)
```

For our dataset, we can use canonical markers to identify our clustering to known cell types:
![Screenshot 2024-12-10 111306](https://github.com/user-attachments/assets/933bb0e2-c7d7-45e2-9626-5aeae0c499a4)

After using expression patterns of marker genes to assign our clusters to cell types, we can use the command shown below to label each cluster on the plot. 
```
new.cluster.ids <- c("Naive CD4 T", "CD14+ Mono", "Memory CD4 T", "B", "CD8 T", "FCGR3A+ Mono",
    "NK", "DC", "Platelet")
names(new.cluster.ids) <- levels(pbmc)
pbmc <- RenameIdents(pbmc, new.cluster.ids)
DimPlot(pbmc, reduction = "umap", label = TRUE, pt.size = 0.5) + NoLegend()
```

![labelplot-1](https://github.com/user-attachments/assets/9d79e5cf-046f-40af-91b0-5c6c40495e64)

# Real World Applications of Seurat<a name="RealWorld"></a>
-  Cancer Research
  - Seurat can be used to understand tumor heterogeneity by identifying distinct cell populations within a tumor. For example, researchers can target rare cancer stem cells, assess immune infiltration, and map gene expression between healthy and diseased cells.
- Developmental Biology:
  - Cell differentation over time during development can also be studied using Seurat. Through analysis of gene expression dynamics, lineage trajectories and roles of regulators in embryonic and organ development can be revealed.
- Immune Profiling:
  - Seurat can also be used in characterizing immune cell diversity in health and disease.
- Integration Across Datasets:
  - Additionally, Seurat can be used to combine datasets from different experiments such as combining RNA and protein data to improve analyses through a more comprehensive appraoch.

# **Strengths of Seurat<a name="Strengths"></a>**
- Scalability
  - Handles large datasets with thousands of cells and hundreds of genes.

- Visualization Power
  - Provides methods such as Uniform Manifold Approfiximation and Projection (UMAP) and t-distributed Stochastic Neighbor Embedding (t-SNE) for visualizing data. This allows for researchers to understand relationships between cell populations and detect patterns.

- Flexibility
  - Supports integration of other R packages and accomodates multimodal data.
  - Parts of Seurat can also be installed as the program is modular allowing for just the necessities to be installed, saving space and time.

- Extensive Documentation
  - As an R program, there is a lot of online resources and a large community to ask questions to. This makes Serurat very accessible to larger amount of researchers with varying expertise levels.


# **Limitations of Seurat<a name="Limitations"></a>**
- Technical Expertise:
  - Due to Seurat being an R program, users need to be proficient in R in order to fully utilize Seurat. However, there is a lot of resources to get help.

- Computational Demands:
  - Although Seurat can handle large datasets there will be siginficant memory and processing power usage because of it.

- Interpretation Challenges:
  - Clustering results are sensitive to normalization and parameters that are set when running Seurat. Users have to know what they are looking for in the data.

# **Conclusion<a name="Conclusion"></a>**
Despite some limitations, Seurat is indispensible for uncovering biological insights from scRNA-seq data. Its flexibility, scalability, and visualization capabilites make it a valuable tool for uncovering biological insights in various fields from cancer to developmental biology. With continued support from the community and its creators, Seurat will remain a tool that all bioinformaticists need in their toolbelt.

## **References**
1. Satija Lab. "Guided Clustering Tutorial." Seurat v5, Satija Lab, https://satijalab.org/seurat/articles/pbmc3k_tutorial.html
2. Satija Lab. "Installation Guide for Seurat v5." Seurat v5, Satija Lab, Accessed 10 Dec. 2024, https://satijalab.org/seurat/articles/install_v5
3. Principle Gallery. "Technique Tuesday: Pointillism Take Two." Principle Gallery, https://www.principlegallery.com/technique-tuesday-pointillism-take-two/
4. Wikipedia contributors. "Georges Seurat." Wikipedia, Wikimedia Foundation, https://en.wikipedia.org/wiki/Georges_Seurat
5. Jovic, Dragomirka et al. “Single-cell RNA sequencing technologies and applications: A brief overview.” Clinical and translational medicine vol. 12,3 (2022). https://doi.org/10.1002/ctm2.694
7. 10x Genomics. "3k PBMCs from a Healthy Donor (v1.1)." 10x Genomics, https://www.10xgenomics.com/datasets/3-k-pbm-cs-from-a-healthy-donor-1-standard-1-1-0


