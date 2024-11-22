# Advanced Machine Learning: Applications in Spatial Transcriptomics
### By Dani Rahman, Leo Joseph, Yashwin Madakamutil

* [What is Spatial Transcriptomics?](#what-is-spatial-transcriptomics)
* [Relevance to Cancer](#relevance-to-cancer)
* [Applications I: Clustering](#applications-i-clustering)
* [Applications II: Deconvolution](#applications-ii-deconvolution)
* [Machine Learning Models for Spatial Transcriptomics](#machine-learning-models-for-spatial-transcriptomics)
* [HistoGene: An In-Depth Look](#histogene-an-in-depth-look)
* [Sources](#sources)

---

## What is Spatial Transcriptomics?

Spatial transcriptomics (ST) is a cutting-edge approach that combines genomics and spatial information to localize gene expression within tissue samples. It bridges the gap left by single-cell RNA sequencing, which lacks the ability to maintain spatial context due to tissue dissociation.

### Key Techniques
- **In Situ Capture**: Tags tissue coordinates with oligonucleotide barcodes, mapping extracted mRNA back to specific tissue regions.
- **In Situ Hybridization**: Hybridizes oligos to RNA at known positions and gleans spatial data from that.
-**In Situ Sequencing**: Uses a circular probe arm and an RCA amplification to tag known RNA

### Advantages
- Provides a detailed map of cellular activity within a tissue's structural context.
- Facilitates studies of complex systems like tumor microenvironments or brain tissue architecture.

- In Situ Hyb and Seq are both more cost efficient than In Situ Capture, and provide more accurate results. Their main trade off is that they are reliant on targeting known genes.

- In Situ Capture has the advantage of being able to identify novel spatial information for genes

<!-- ![Spatial Transcriptomics Example](img/spatial_example.png) -->
![](https://github.com/drahmanucsd/beng183proj/blob/main/img2.png)
![](https://github.com/drahmanucsd/beng183proj/blob/main/img1.png)



---

## Relevance to Cancer

Tumor heterogeneity, which refers to genetic variability and differing gene expression within a tumor, is a major challenge in oncology. 

Understanding this variability is important, because this way researchers can tailor personalized treatments on a patient by patient basis.

- By overlaying spatial info with sequencing data, researchers have been able to identify tissues that show evidence of infiltration by cancer cells, which generally has a better prognosis and patient outcome if identified in the early onset stages. 

### Importance of Spatial Data
- **Identifies Cellular Subpopulations**: Based on the spatial data of differing gene expression, spatial data is able to identify cell populations.
- **Tailors Therapies**: Enables development of targeted treatments for specific tumor regions based on the knowledge of differing gene expression at these locations.
- **Reveals Microenvironments**: Investigates the interaction between tumor cells and their surroundings.

<!-- ![Tumor Heterogeneity](img/tumor_heterogeneity.png) -->

![](https://github.com/drahmanucsd/beng183proj/blob/main/img3.png)

---

## Applications I: Clustering

Clustering techniques are invaluable for analyzing spatial transcriptomics data. By grouping data based on similarity, clustering reveals patterns in tissue organization and gene expression.

## Dimensionality Reduction

### **Principal Component Analysis (PCA)**  
PCA is a linear dimensionality reduction technique designed to reduce the number of features in a dataset while retaining as much variance as possible. This is achieved by transforming the original data into a new coordinate system defined by **principal components**, which are linear combinations of the original features.  

#### **Key Steps**:  
- **Standardization**: Normalize data to ensure each feature contributes equally to variance.  
- **Covariance Matrix**: Compute the covariance matrix of the standardized data to measure feature relationships.  
- **Eigenvalues and Eigenvectors**: Identify the eigenvalues (variance) and eigenvectors (directions) of the covariance matrix.  
- **Projection**: Project the data onto a smaller subset of principal components ranked by eigenvalue magnitude.  

#### **Advantages**:  
- Effective for high-dimensional datasets where features are highly correlated.  
- Enables noise reduction and visualization of multi-dimensional data.  
- Provides interpretable axes representing variance in the data.  

#### **Limitations**:  
- Assumes linear relationships between features, which may not always hold.  
- Does not preserve local or global geometric structures in non-linear datasets.  

---
![](https://github.com/drahmanucsd/beng183proj/blob/main/img4.png)
![](https://github.com/drahmanucsd/beng183proj/blob/main/img5.png)
### **Uniform Manifold Approximation and Projection (UMAP)**  
UMAP is a non-linear dimensionality reduction technique that emphasizes preserving the local and global structure of data when projecting it into lower dimensions (e.g., 2D or 3D). It is particularly suited for visualizing clusters and patterns in complex datasets.  

#### **Key Concepts**:  
- **Graph Construction**: Creates a weighted graph of the data, where nodes represent data points and edge weights represent their similarities (based on nearest neighbors).  
- **Manifold Approximation**: Models the high-dimensional data as a manifold, aiming to preserve the relationships between points.  
- **Optimization**: Uses stochastic gradient descent to optimize the layout of points in lower dimensions while minimizing distortion of the graph structure.  

#### **Advantages**:  
- Captures both local relationships (e.g., nearest neighbors) and global structure (e.g., clusters).  
- Efficient and faster than similar methods like t-SNE.  
- Ideal for exploratory data analysis and visualizing patterns in large datasets.  

#### **Limitations**:  
- Less interpretable than PCA due to its non-linear nature.  
- Sensitive to hyperparameters such as the number of neighbors and minimum distance between points.  
- Requires careful tuning to balance local and global preservation.  

---

### **Comparison**:  
| Feature                | PCA                               | UMAP                              |
|------------------------|-----------------------------------|-----------------------------------|
| **Type**              | Linear                           | Non-linear                       |
| **Preserves**         | Global variance                  | Local and global relationships   |
| **Speed**             | Faster for small datasets        | Scales well to large datasets    |
| **Interpretability**  | High (eigenvectors as axes)      | Low (manifold projection)        |
| **Applications**      | Feature reduction, preprocessing | Visualization, clustering         |

Both techniques are valuable tools depending on the dataset's nature and the analysis goals. PCA is best for variance retention and interpretability, while UMAP excels in visualizing clusters in complex, non-linear datasets.

### Applications
- **Inter-Tumor Heterogeneity**: Distinguishes genetic differences across tumor sections.
- **Mapping Subpopulations**: Identifies unique cell types and their spatial distributions.

<!-- ![Clustering Example](img/clustering.png) -->
![](https://github.com/drahmanucsd/beng183proj/blob/main/img6.png)

---

## Applications II: Deconvolution

Deconvolution predicts the proportions of various cell types in a tissue, enhancing the resolution of spatial transcriptomics data.

### Redeconve
- Utilizes linear regression to estimate cell type proportions.
- Efficiently detects differences in cell populations, such as cytotoxic T cells in immune responses.

Deconvolution reduces costs and computational complexity, making it accessible for large-scale studies.

<!-- ![Deconvolution Example](img/deconvolution.png) -->
![](https://github.com/drahmanucsd/beng183proj/blob/main/img7.png)

---

## Machine Learning Models for Spatial Transcriptomics

Machine learning models are essential for processing and interpreting spatial transcriptomics data:

1. **Convolutional Neural Networks (CNNs) + Multi-Layer Perceptrons (MLPs)**:
   - Encode spatial relationships and classify gene expression patterns.

2. **Vision Transformers (ViT)**:
   - Segment images into tokens and learn their interrelations.
   - Effective for long-range spatial dependency modeling.

3. **Graph Neural Networks (GNNs)**:
   - Represent tissue regions as nodes in a graph.
   - Capture spatial dependencies through node connectivity.

### Impacts
- **Pathology Advances**: Provides high-resolution diagnostic tools.
- **Cost-Effective Precision Medicine**: Accelerates and reduces the cost of personalized therapies.

<!-- ![Machine Learning Models](img/ml_models.png) -->
![](https://github.com/drahmanucsd/beng183proj/blob/main/img8.jpeg)
<!-- ![](https://github.com/drahmanucsd/beng183proj/blob/main/img9.png) -->


---

## HistoGene: An In-Depth Look

HistoGene is an advanced spatial transcriptomics platform that integrates histological imaging with gene expression data, providing a comprehensive view of tissue architecture and molecular function.

### How HistoGene Works
1. **Sample Preparation**:
   - Tissue samples are placed on specialized slides with capture probes that bind to mRNA molecules. 
   - The sample is stained for histological imaging, allowing visualization of cellular structures.

2. **RNA Capture and Sequencing**:
   - mRNA from each tissue region hybridizes to spatially barcoded probes.
   - Extracted RNA is sequenced, and each read is mapped back to its spatial coordinate.

3. **Data Integration**:
   - High-resolution histological images are overlaid with gene expression maps.
   - Machine learning models analyze spatial relationships between histological features and gene activity.

4. **Machine Learning Enhancements**:
   - **Clustering**: Groups regions with similar expression profiles to reveal functional areas.
   - **Predictive Models**: Identifies potential biomarkers by correlating spatial gene expression with disease phenotypes.

### Key Applications
- **Tumor Analysis**:
  - Identifies regions of high genetic variability and immune cell infiltration.
  - Pinpoints tumor-stroma interactions critical for metastasis.
- **Neurological Research**:
  - Maps gene expression in brain tissues, revealing functional regions and pathways.
- **Drug Development**:
  - Predicts drug efficacy by analyzing molecular profiles within heterogeneous tissues.

### Advantages
- Combines histological and transcriptomic insights for holistic tissue analysis.
- Improves diagnostic accuracy and enhances drug target discovery.

<!-- ![HistoGene Workflow](img/histogene.png) -->
![](https://github.com/drahmanucsd/beng183proj/blob/main/img10.gif)

---

## Sources
1. [Visium HD Spatial Gene Expression - 10x Genomics](https://www.10xgenomics.com/products/visium-hd-spatial-gene-expression)
2. [Overview of In-Situ Capturing Technologies - ResearchGate](https://www.researchgate.net/figure/Overview-of-in-situ-capturing-technologies-ST-10X-Visium-Microscopic-glass-slides_fig3_341130420)
3. [Spatial Transcriptomics and In Situ Analysis - Nature](https://www.nature.com/articles/s41592-022-01409-2)
4. [Definition of Tumor Heterogeneity - National Cancer Institute](https://www.cancer.gov/publications/dictionaries/cancer-terms/def/tumor-heterogeneity)
5. [Spatial Heterogeneity of Tumors: ST Techniques - ResearchGate](https://www.researchgate.net/figure/Spatial-heterogeneity-of-the-tumor-A-ST-techniques-have-been-used-to-characterize-the_fig1_364518905)
6. [Spatial Assignment of Cell Types: Clustering and UMAP Visualization - ResearchGate](https://www.researchgate.net/figure/Spatial-assignment-of-cell-types-a-Clustering-and-UMAP-visualization-of-cells-based-on_fig5_359945870)
7. [Advances in Spatial Transcriptomics - Nature Communications](https://www.nature.com/articles/s41467-023-43600-9#citeas)
8. [Emerging Methods in Spatial Transcriptomics - Oxford Academic](https://academic.oup.com/bib/article/25/1/bbad464/7494746)

---
`


