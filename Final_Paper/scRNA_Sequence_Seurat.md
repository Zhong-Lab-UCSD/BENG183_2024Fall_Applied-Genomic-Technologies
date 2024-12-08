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

The fundamental object of 3C(Chromosome Conformation Capture) techniques and 3C-derived methods is to understand the physical wiring diagram of the genome by identifying the physical interaction between chromosomes. 

To capture the interaction (crosslink between strings), there are few steps in general:
- Take a snapshot of the flowing cells - **Crosslink** with fixative agent (formaldehyde)
- Zoom in on crosslinked part and exclude untangled parts - **Digested** with a restriction enzyme
- Analyze the components come from the same chromatin - **Reverse crosslink** and **sequence**
- Finish the jigsaw puzzle and get the results - **Align** the reads and **summarize** the contacts

> Based on these general ideas, then we'll dive deeper by walking through two of the most popular  techniques and then briefly introduce some other methods. 

## Why Should You Use Seurat

![](/assets/1-s2.0-S1360138518300827-gr1b2_lrg.jpg)
[Figure1](https://doi.org/10.1016/j.tplants.2018.03.014). Schematic Representation of Chromosome Conformation Capture (3C) and 3C-Derived Methods. These methods help to elucidate nuclear organization by detecting physical interactions between genetic elements located throughout the genome. Abbreviations: IP, immunoprecipitation; RE, restriction enzyme. **Figure by Sotelo-Silveira, Mariana, et al. Trends in Plant Science (2018).**

To better understand the difference between these methods, I'd like to distingush them between the following couple of aspects:

 1) Specificity - What does _one, all, many_ mean<a name="2321"></a>
‘1’, ‘Many’ and ‘All’ indicate how many loci are interrogated in a given experiment. For example, ‘1 versus All’ indicates that the experiment probes the interaction profile between 1 locus and all other potential loci in the genome. ‘All versus All’ means that one can detect the interaction profiles of all loci, genome-wide, and their interactions with all other genomic loci [1].

These kind of specificity is determined by the primer when people use **specific primers** before PCR. 

# What is Clustering









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
## Perform Clustering With Clustering Commands
In order to begin clustering, it is important to know which commands to use to create a clustering plot. Below are all the Seurat commands with a description of each command you can use with the Seurat package.

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

![](https://github.com/nmora2/ScRNA-Seq---Seurat/blob/c5fd38651e2a375fc65739eff41a0a31cdb1584a/Screenshot%202024-12-07%20at%204.59.19%20PM.png)


# Different Applications Performance Compared to Seurat
- Hi-C original: [Lieberman-Aiden et al., Science 2010](doi: 10.1126/science.1181369)
- Hi-C 1.0: [Belton-JM et al., Methods 2012](doi: 10.1016/j.ymeth.2012.05.001)
- In situ Hi-C: [Rao et al., Cell 2014](doi: 10.1016/j.cell.2014.11.021)
- Single cell Hi-C: [Nagano et al., Genome Biology 2015](https://doi.org/10.1186/s13059-015-0753-7)
- DNase Hi-C [Ma, Wenxiu, Methods et al](https://www.ncbi.nlm.nih.gov/pubmed/25437436)
- Hi-C 2.0: [Belaghzal et al., Methods 2017](https://www.ncbi.nlm.nih.gov/pubmed/28435001)
- DLO-Hi-C: [Lin et al., Nature Genetics 2018](https://doi.org/10.1038/s41588-018-0111-2)
- Hi-C improving: [Golloshi et al., Methods 2018](https://www.biorxiv.org/content/biorxiv/early/2018/02/13/264515.full.pdf)
- Arima 1-day Hi-C: [Ghurye et al., BioRxiv 2018](https://www.biorxiv.org/content/early/2018/02/07/261149)

## Insert Different Application1
ChIA-PET is another method that combines ChIP and pair-end sequencing to analysis the chromtin interaction. It allows for targeted binding factors such as: estrogen receptor alpha, CTCF-mediated loops, RNA polymerase II, and a combination of key architectural factors. on the one hand, it has the benefit of achieving a higher resolution compared to Hi-C, as only ligation products involving the immunoprecipitated molecule are sequenced, on the other hand, ChIA-PET has systematic biases due to ChIP process:
- Only one type of binding factor selected
- Different antibodies
- ChIP conditions

## Insert Different Application2








# Selected Methods Comparison
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
