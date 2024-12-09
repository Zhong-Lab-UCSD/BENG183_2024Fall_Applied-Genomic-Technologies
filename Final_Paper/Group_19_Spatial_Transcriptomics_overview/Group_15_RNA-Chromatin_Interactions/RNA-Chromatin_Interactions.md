# RNA-Chromatin Interactions: Applications of Genomic Technologies

### By Amick Licup, Chloe Keggen, Jerold Reputana

---

## Overview 
The human body is one of the most complex mysteries scientists across disciplines and time have worked tirelessly to unravel. Even with bleeding-edge technological advancements, we are still only beginning to understand the vast & complex systems/interactions that come from and affect our bodies. While that may be the case, that hasn’t stopped us from developing and utilizing genomic technologies to tackle specific biological questions– one at a time in a series of logical steps that reveal more of the truth.

In particular, diabetes is one of the most studied health issues and affects a majority of people worldwide. While there are a multitudinous amount of studies on the topic, we will be focusing on a particular research paper that delves into the role of chromatin-associated RNA (caRNA) in endothelial cell dysfunction– a critical factor in diseases like diabetes.

Specifically, what questions did the scientists ask? What genomic technologies did they then use to answer those questions? What insights did that provide them and what were their next logical steps? These are the questions we will be asking ourselves in discovering how these scientists solved this specific aspect of a larger problem with genomic technologies.

## Endothelial Cells & Dysfunction

## What are caRNA & RNA-chromatin interactions?

## Single-cell RNA-seq
In the process of focusing on RNA-chromatin interactions, there were several steps involving applications of different genomic technologies to solve a series of questions. One of the first problems the scientists faced was on how to identify the effects of diabetes stress conditions on epithelial cells– specifically HUVECs (Human Umbilical Vein Endothelial Cells). Essentially, the solution to this problem was to measure & identify the gene expression of endothelial cells (EC) via single-cell RNA-seq, or scRNA-seq. 

As a reminder, scRNA-seq provides a higher resolution and insight, compared to the average of many cells in traditional RNA-seq by focusing on individual cells. The RNA of these cells is then converted to cDNA and amplified for sequencing/mapping and determining gene expression patterns.

![](https://github.com/alicup29/BENG183_2024Fall_Applied-Genomic-Technologies/Final_Paper/Group_15_RNA-Chromatin_Interactions/img1.jpg)

In setting up the in-vitro experiment, HUVECs were subjected to three temporally different concentrations of “...high glucose and TNFα [that]…” to mimic the diabetic hyperglycemia and chronic inflammation that induces EC dysfunction (Calandrelli, et al. 2020). Day 0 corresponds to normal glucose conditions, Day 3 to high glucose/TNFα conditions, and Day 7 had the same conditions but for 7 instead of 3 days.

![](https://github.com/alicup29/BENG183_2024Fall_Applied-Genomic-Technologies/Final_Paper/Group_15_RNA-Chromatin_Interactions/img2.jpg)

Each group was subjected to scRNA-seq resulting in various graphs that reveal “...a profound gene expression and phenotypic change…” (Calandrelli, et al. 2020).

![](https://github.com/alicup29/BENG183_2024Fall_Applied-Genomic-Technologies/Final_Paper/Group_15_RNA-Chromatin_Interactions/img3.jpg)

Following scRNA-seq, one such graph is this expression heatmap that displays the top differentially expressed genes of the treated cells, grouped by function and ordered specifically by a SERPINE1 biomarker, which is indicative of stress-induced dysfunction. Additionally, the heatmap is z-scaled, representing the standard deviations away from the mean SERPINE1 expression. After classifying the top differentially expressed genes, there is a significant increase in expression of inflammatory/immune response genes, which are characteristic of EC reactions to diabetes stress. In terms of the z-scale, there is a majority shift from very low expression (blue in day 0) to extremely high expression (yellow in day 3 and 7).

## Hi-C & iMARGI

As a result of scRNA-seq, the scientists were able to find which differentially expressed genes should be investigated, leading to the next question– how is this EC stress reflected in genomic interactions? 

In order to answer this question, they first utilized Hi-C (Hierarchical-Chromatin organization), a genomics method used to study the 3D organization of the genome and identify physical interactions between different regions of DNA within the cell nucleus. By cross-linking interacting DNA and associated proteins, then fragmenting and sequencing the DNA strands, Hi-C provides important insight into chromatin organization and gene regulation. For example, topologically associated domains (TADs) which are contiguous and long regions subset of chromosomes that are spatially clustered together and frequently interact. 

![](https://github.com/alicup29/BENG183_2024Fall_Applied-Genomic-Technologies/Final_Paper/Group_15_RNA-Chromatin_Interactions/img4.jpg)

As a result of the in situ Hi-C experiments, there were no significant differences between the proportions of read pairs classified as intrachromosomal or genomic-level TAD/A & B compartment changes across all three days. Essentially, the stress conditions “...did not significantly perturb the major 3D genome features…” of the ECs (Calandrelli, et al. 2020).

![](https://github.com/alicup29/BENG183_2024Fall_Applied-Genomic-Technologies/Final_Paper/Group_15_RNA-Chromatin_Interactions/img5.jpg)

If a majority of the 3D genome organization in ECs was not affected, then the next question to be asked was– what about changes to RNA-chromatin interactions? In order to tackle this problem, in-situ MApping of RNA-Genome Interaction (iMARGI) was used. iMARGI is used to analyze how RNA molecules physically associate w/ specific regions of DNA at specific genomic loci by crosslinking RNA and local chromatin (DNA), then using NGS to align and determine which RNA interacts with specific chromatin regions.

![](https://github.com/alicup29/BENG183_2024Fall_Applied-Genomic-Technologies/Final_Paper/Group_15_RNA-Chromatin_Interactions/img6.jpg)

This process produces RNA-DNA hybrids, or “read pairs”. Between day 0 and day 3 & 7 sample ECs, there was a significant difference in the unique iMARGI read pairs classified as interchromosomal, jumping from 34.4% (day 0) to 62.7% (days 3 & 7)-- indicating that diabetes stress conditions “...induced interchromosomal RNA–chromatin interactions in ECs…” (Calandrelli, et al. 2020).

From here, it becomes a matter of asking what type of interchromosomal interaction causes this EC response. Looking further into the obtained iMARGI read pairs (majority RNA ends over DNA ends), the specific locations are a majority RNA-chromatin interactions at super enhancer regions. Specifically, long-intergenic non-coding RNA (a type of chromatin-associated RNA) from the super enhancer “Linc607SE” are the main culprit, and the start of the super enhancer interaction cascade originating from cell stress and ultimately causing EC dysfunction.

![](https://github.com/alicup29/BENG183_2024Fall_Applied-Genomic-Technologies/Final_Paper/Group_15_RNA-Chromatin_Interactions/img7.jpg)

Overall, it was through a mixture of several key genomic technologies (scRNA-seq, Hi-C, & iMARGI) and analyses of their output data that allowed the scientists to unravel the cascade of RNA-chromatin interactions, with just EC dysfunction.

## Key Findings

This study explores how stressed-induced RNA-chromatin interactions promote endothelial dysfunction, specifically a condition linked to diseases like diabetes. The researchers used different advanced genomic technologies, including single-cell RNA sequencing (scRNA-seq), Hi-C, and iMARGI sequencing, in order to analyze endothelial cells under stress caused by high glucose and inflammation. 

Their findings revealed three critical points:
1. RNA-chromatin interactions change during stress. By using iMARGI sequencing, the researchers discovered that stress caused an increase in interchromosomal RNA-chromatin interactions, especially among regions called super enhancers. Those super enhancers are key regulators of gene expressions. 

2. The role of LINC00607. One RNA that the researchers thought was significant was LINC00607 which emerged as a major promoter of endothelial dysfunction by interacting with chromatin regions near SERPINE1 (a pro-inflammatory and pro-fibrotic gene). By disrupting the interaction has led to reduced dysfunction-related gene expression and improved cell behavior. 

3. The stable 3D genome structure. When the researchers performed Hi-C sequencing, the results showed that the overall 3D structure of the genome remained stable. This suggests that RNA-chromatin interactions can primarily promote dysfunction independently of major structural changes. 

## Summary

This research ties closely to the concepts we have learned in class about applied genomic technologies and their real work applications. For instance, sc-RNA sequencing can be used to reveal cell to cell variability in gene expression during dysfunction, while iMARGI sequencing can be used to directly map how RNA interacts with chromatin. These are the tools we have learned in class, where applying these techniques were used to answer real world questions about endothelial dysfunction. Additionally, the process of discovering therapeutic implications. The research suggests that targeting RNA-chromatin interactions, like the one involving LINC00607, could lead to new treatments for diseases driven by endothelial dysfunction, such as diabetes or cardiovascular conditions. 

Overall, this study highlights the power of genomic tools in investigating how molecular interactions shape cell behavior under stress. It not only advances our understanding of endothelial dysfunction but also demonstrates how applied genomic technologies can bridge fundamental research and therapeutic discoveries. This research article is a clear example of how the concepts we studied in class directly inform the proper workflow in research and discover future treatments. 
