# Visualizing Differential Gene Expression with Volcano Plots
### By Ishvari Desai, Emily Wei, Soyeon Lee(Amy)

* [Introduction](#introduction)
* [Volcano Plot](#volcano-plot)
* [Generating Volcano Plots](#generating-volcano-plots)
* [Case Study: Lung Cancer](#case-study-lung-cancer)
* [Case Study: Results](#case-study-results)

## Introduction

Gene expression studies often result in large datasets that contain detailed measurements of RNA or protein levels across various experimental conditions. To derive meaningful insights from these datasets, effective visualization methods are crucial. These visualizations provide a clear picture of how genes are expressed differently under specific conditions.

Visualizing gene expression allows researchers to understand and interpret differences in gene activity between various biological conditions. For example, visualizations enable researchers to pinpoint specific genes of interest, interpret complex datasets, and present results in an intuitive and comprehensible format. They play a critical role in hypothesis generation, pathway analysis, and identifying candidate biomarkers or therapeutic targets.

### Differential Gene Expression
Differential gene expression (DGE) refers to the process of identifying genes whose expression levels vary significantly between two or more experimental conditions. These changes may indicate biological responses to treatments, diseases, or environmental factors.

### Key Tools for DGE Analysis
Commonly used tools for DGE include:

* DESeq2: R package for robust statistical testing of differential expression.

* edgeR: Another R package widely used for count-based expression data.

These tools typically output statistical metrics like fold change and adjusted p-values, which are essential for volcano plot visualization.


## Volcano Plot
A volcano plot is a scatter plot designed to visualize gene expression changes between conditions. It combines the magnitude of expression changes (fold change) with their statistical significance (p-value) to highlight the most impactful genes.

<img width="400" alt="Screen Shot 2024-12-09 at 3 45 37 PM" src="https://github.com/user-attachments/assets/63c7e94c-2cad-4bf7-bddb-15bee2d9b90c">

The provided image depicts a volcano plot, a scatter plot used to visualize differential gene expression. The x-axis represents the log2 fold change, which indicates the magnitude of change in gene expression between two conditions; values further from zero signify greater changes. The y-axis displays the -log10 p-value, representing the statistical significance of these changes, with higher values indicating stronger confidence in the results. 

Up-regulated Genes:

* Highlighted in red, on the right side of the plot.
* Represent genes with significantly increased expression in the condition of interest.

Down-regulated Genes:
* Highlighted in blue, on the left side of the plot.
* Represent genes with significantly decreased expression.

Non-significant Genes:
* Plotted in grey, representing genes with no statistically significant expression changes.

The horizontal dashed line corresponds to the significance threshold, typically set at a p-value of 0.05; genes above this line are considered statistically significant. The fold change threshold is also used to classify genes based on their biological significance, making the volcano plot an effective tool for identifying key genes for further investigation.


## Generating Volcano Plots

Volcano plots can be generated using tools like R and Python.

* R: Libraries such as ggplot2 and EnhancedVolcano are widely used.

* Python: Libraries like Matplotlib and Seaborn offer similar functionality

Using the ggplot2 library, we map gene expression changes and highlight regions of significance.

**Code Example in R:**
```
library(ggplot2)

ggplot(data, aes(x=log2FoldChange, y=-log10(pvalue))) +
  geom_point(aes(color=significance)) +
  geom_hline(yintercept=-log10(0.05), linetype="dashed") +
  geom_vline(xintercept=c(-1, 1), linetype="dashed") +
  theme_minimal()

```
*Key Features:*

* geom_point: Plots each gene as a point.

* Color Coding: Points are colored based on significance.

* Dashed Lines: Represent p-value and fold-change thresholds.

*Output:*
A clear volcano plot highlighting up-regulated and down-regulated genes.

**Code Example in Python:**
```
import matplotlib.pyplot as plt
import numpy as np

# Example data
data = {'log2FoldChange': np.random.normal(size=1000), 
        'pvalue': np.random.uniform(0.001, 0.1, size=1000)}

# Add significance column
data['significance'] = ['Significant' if abs(x) > 1 and y < 0.05 else 'Not Significant' 
                        for x, y in zip(data['log2FoldChange'], data['pvalue'])]

# Extract values
x = data['log2FoldChange']
y = -np.log10(data['pvalue'])
colors = ['red' if sig == 'Significant' else 'gray' for sig in data['significance']]

# Plot
plt.figure(figsize=(8, 6))
plt.scatter(x, y, c=colors, alpha=0.7)
plt.axhline(-np.log10(0.05), color='blue', linestyle='dashed', linewidth=1)
plt.axvline(-1, color='green', linestyle='dashed', linewidth=1)
plt.axvline(1, color='green', linestyle='dashed', linewidth=1)
plt.title('Volcano Plot')
plt.xlabel('log2 Fold Change')
plt.ylabel('-log10(P-value)')
plt.show()
```
*Key Features:*

* plt.scatter: Plots genes as points, colored based on significance.

* Dashed Lines: Highlight thresholds for significance and fold-change.

* Customization: Matplotlib allows extensive styling and customization.

*Output:*
A volcano plot similar to R’s visualization, emphasizing significant gene expression changes.

Both R and Python excel at generating volcano plots. The choice of tool depends on the user’s familiarity and the specific requirements of the analysis.

## Case Study: Lung Cancer

<img width="720" alt="research_paper_title" src="https://github.com/user-attachments/assets/3d5f3e26-6d84-41de-9976-493b86c4b6b5">

### Differential Expression in Cancer vs. Normal Tissue

Researchers often use volcano plots to visualize how certain genes might be expressed differentially in cancerous tissue versus normal tissue. Such insight provides a clearer target for gene therapies and other methodologies in the specific form of cancer they are researching. 

### Focus of the Paper

PAICS is a metabolic enzyme that has been identified as an oncogene, which is a mutated gene that often causes cancer in many tumor types, such as breast cancer and prostate cancer. In general, cancer is the result of uncontrolled cellular proliferation of mutated cells. The specific molecular mechanism that causes such proliferation tends to differ among different cancer types. Whether it be the mutation of the p53 tumor suppressor gene, preventing necessary apoptosis of mutated cells, or the overactivation of molecular components that promote cellular proliferation, there are many possible ways in which cancer can be caused. Because different types of cancer often have different causes, the role of PAICS in lung adenocarcinoma, or LADC, is unknown and is the focus of this paper's research. 

The authors of this paper highlighted the importance of their research by explaining that lung cancer is the leading cause of cancer-related mortality. About 78.1 of every 10,000 people in China are reported to have suffered from LADC specifically. Due to the fact that metabolic reprogramming has been the focus of cancer research recently, understanding the relationship between the dysregulated metabolism of enzymes and metastasis is necessary to producing effective lung cancer treatments. 

## Case Study: Results

### Plot Criteria

Utilizing the genomics data available on the database cBioPortal to obtain 558 genes that were significantly coexpressed with PAICS in LADC samples, the researchers created a volcano plot with **threshold criteria** of a *log2 fold change of +/- 1* and a *p-values less than 0.05*. 

**The resulting volcano plot looked like:**

<img width="392" alt="Screen Shot 2024-12-09 at 12 55 12 AM" src="https://github.com/user-attachments/assets/b7c99009-821a-448c-b076-3d9c78846797">

*Figure A's Description:* "The expression differences of coexpression genes, obtained from cBioPortal, between altered and unaltered PAICS expression group, are shown in a volcano plot.

Zhou, Shuyi et al. (2019) PAICS is hypomethylated and highly expressed in LADC. ScienceDirect, Gene, Vol. 692, 1-8.

### Plot Results/Analysis

Based on the threshold, **215 genes** were identified as PAICS-associated codifferentially expressed genes. They were thus likely involved in biological pathways related to PAICS expression, helping researchers to better understand potential therapeutic targets. 

*After further analysis, the researchers discovered an important association between the coexpressed genes and tumor development:*

The discovery they made was that a majority of the genes co-expressed with PAICS were found to be involved in cell-proliferation related processes like cellular division and mitotic nuclear division. Uncontrolled cellular proliferation is what results in metastatic tumors (cancer). Thus, these results suggest that PAICS-related genes might be playing a role in promoting this process. 

Since many of these co-expressed genes were associated with cellular proliferation, this volcano plot indicates that PAICS itself likely plays a role in these cancer-causing processes as well. 

## Summary

In summary, volcano plots serve as a very important tool in research that looks at differential gene expression, such as cancer research. There are many possible tools that are available to generate a volcano plot in R or python. These tools make it easy for researchers to generate a volcano plot based on their collected data and use that plot to determine which genes were differentially expressed. From there, they can look at the specific functionality of those genes to try and understand the molecular mechanism for the type of disease or cancer they are researching. 

## Sources 
* [Hemtools](https://hemtools.readthedocs.io/en/latest/content/Bioinformatics_Core_Competencies/Volcanoplot.html)
* [Erikduan](https://erikaduan.github.io/posts/2021-01-02-volcano-plots-with-ggplot2/)
* [Understanding Volcano Plots](https://www.htgmolecular.com/blog/2022-08-25/understanding-volcano-plots)
* [Data Analysis and Visualization](https://mkempenaar.github.io/gene_expression_analysis/chapter-5.html)
* ["Roles of highly expressed PAICS in lung adenocarcinoma"](https://www.sciencedirect.com/science/article/pii/S0378111919300228?via%3Dihub)
