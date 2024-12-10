# RNA-Chromatin Interactions: Applications of Genomic Technologies

### By Amick Licup, Chloe Keggen, Jerold Reputana

* [Overview](#overview)
* [Endothelial Cells & Dysfunction](#endothelial-cells--dysfunction)
* [What are caRNA & RNA-chromatin interactions?](#what-are-carna--rna-chromatin-interactions)
* [Single-cell RNA-seq](#single-cell-rna-seq)
* [Hi-C & iMARGI](#hi-c--imargi)
* [Key Findings](#key-findings)
* [Summary](#summary)
* [Sources](#sources)

---

## Overview 
The human body is one of the most complex mysteries scientists across disciplines and time have worked tirelessly to unravel. Even with bleeding-edge technological advancements, we are still only beginning to understand the vast & complex systems/interactions that come from and affect our bodies. While that may be the case, that hasn’t stopped us from developing and utilizing genomic technologies to tackle specific biological questions– one at a time in a series of logical steps that reveal more of the truth.

In particular, diabetes is one of the most studied health issues and affects a majority of people worldwide. While there are a multitudinous amount of studies on the topic, we will be focusing on a particular research paper that delves into the role of chromatin-associated RNA (caRNA) in endothelial cell dysfunction– a critical factor in diseases like diabetes.

Specifically, what questions did the scientists ask? What genomic technologies did they then use to answer those questions? What insights did that provide them and what were their next logical steps? These are the questions we will be asking ourselves in discovering how these scientists solved this specific aspect of a larger problem with genomic technologies.

## Endothelial Cells & Dysfunction
Endothelial cells line all blood vessels, located throughout the entire body, branching from the heart to tissues and organs. These cells are essential in vascular health and homeostasis, including regulating blood flow, blood vessel permeability, and angiogenesis. Endothelial dysfunction is when the endothelium fails to perform its critical physiological functions, making the body defenseless and at risk to certain diseases and complications. For example, in sepsis or severe infections, vascular leakage may occur where albumin, a necessary protein, leaks out of blood vessels because the endothelial cells fail to line the vessel. (Kim, Cardiology University of Washington) Similarly, endothelial cells are necessary for angiogenesis, the creation of new blood vessels, but when this dysfunction occurs, the production of certain growth factors can cause tumor growth. This disease is especially pertinent for those with diabetes, whose high glucose levels cause endothelial damage. Nitric oxide in vascular health is important for vasodilation, relaxing the muscle cells of vessels, therefore lowering blood pressure. NO reduces the risk of blood clot formation, which reduces heart attacks and strokes. In dysfunction, endothelial cells produce less NO due to a downregulation of endothelial nitric oxide synthase. Reduced NO levels lead to impaired vasodilation, causing hypertension and increased stress on the heart. Endothelial dysfunction is thus a key factor in cardiovascular diseases, diabetes, and cancers; understanding the mechanisms behind cell damage is critical for developing therapies. 

![](https://github.com/alicup29/beng183proj/blob/main/img8.png)

## What are caRNA & RNA-chromatin interactions?

Chromatin-associated RNA, caRNAs, are RNA molecules that physically associate with chromatin to regulate gene expression. They include long non-coding RNAs, small nuclear RNAs, and enhancer RNAs. These RNA molecules directly or indirectly interact with chromatin, a “mixture of DNA and proteins that form the chromosomes found in the cells of organisms” (NHGRI). caRNAs are transcribed from one gene, then travel and attach to a genomic sequence, most often at super enhancers or transcription start sites, on another gene. Super enhancers are a cluster of transcriptional enhancers across a long range of genomic DNA. In diabetic patients’ endothelial cells specifically, super enhancer activity increases in stress. Non-coding RNAs are transcribed at an elevated rate, including caRNAs. The caRNAs then travel to a different enhancer site, a mechanism called ‘inter-super-enhancer communication’.

![](https://github.com/alicup29/beng183proj/blob/main/img11.png)
![](https://github.com/alicup29/beng183proj/blob/main/img9.png)

At the new gene enhancer site, caRNA’s can do multiple things. Focusing on two major factors: first, they can recruit chromatin-modifying enzymes such as histone acetyltransferases and methyltransferases, which modify histones at enhancer sites. This recruitment promotes a change in chromatin structure, which can either suppress or promote transcription. However, during stress and this increase in RNA transcription, caRNAs may recruit modifiers to regions where they are not needed. 
Second, caRNAs may prevent or form chromatin loops by tethering distant genomic regions together. Abnormal loops can disrupt the communication between the enhancers and promoters. For example, caRNAs can stabilize chromatin loops at loci that are essential in immune responses, leading to a higher expression of inflammatory genes. 
In diabetic patients, caRNAs can recruit chromatin modifiers that suppress transcription of endothelial nitric oxide synthase, reducing NO levels and therefore disrupting vasodilation. This complication results in hypertension and an increased risk of cardiovascular diseases. Similarly, as a response to hyperglycemia, the stabilization of chromatin loops causes an upregulation in fibrotic gene programs. There are several more side effects of abnormal caRNA transcription in endothelial cells in diabetic stress, such as an enhanced expression of oxidative stress, disruption of angiogenesis, and dysregulation of the lipid metabolism– all of which have major implications on the body. 

![](https://github.com/alicup29/beng183proj/blob/main/img10.png)

## Single-cell RNA-seq
In the process of focusing on RNA-chromatin interactions, there were several steps involving applications of different genomic technologies to solve a series of questions. One of the first problems the scientists faced was on how to identify the effects of diabetes stress conditions on epithelial cells– specifically HUVECs (Human Umbilical Vein Endothelial Cells). Essentially, the solution to this problem was to measure & identify the gene expression of endothelial cells (EC) via single-cell RNA-seq, or scRNA-seq. 

As a reminder, scRNA-seq provides a higher resolution and insight, compared to the average of many cells in traditional RNA-seq by focusing on individual cells. The RNA of these cells is then converted to cDNA and amplified for sequencing/mapping and determining gene expression patterns.

![](https://github.com/alicup29/beng183proj/blob/main/img1.jpg)

In setting up the in-vitro experiment, HUVECs were subjected to three temporally different concentrations of “...high glucose and TNFα [that]…” to mimic the diabetic hyperglycemia and chronic inflammation that induces EC dysfunction (Calandrelli, et al. 2020). Day 0 corresponds to normal glucose conditions, Day 3 to high glucose/TNFα conditions, and Day 7 had the same conditions but for 7 instead of 3 days.

![](https://github.com/alicup29/beng183proj/blob/main/img2.png)

Each group was subjected to scRNA-seq resulting in various graphs that reveal “...a profound gene expression and phenotypic change…” (Calandrelli, et al. 2020).

![](https://github.com/alicup29/beng183proj/blob/main/img3.png)

Following scRNA-seq, one such graph is this expression heatmap that displays the top differentially expressed genes of the treated cells, grouped by function and ordered specifically by a SERPINE1 biomarker, which is indicative of stress-induced dysfunction. Additionally, the heatmap is z-scaled, representing the standard deviations away from the mean SERPINE1 expression. After classifying the top differentially expressed genes, there is a significant increase in expression of inflammatory/immune response genes, which are characteristic of EC reactions to diabetes stress. In terms of the z-scale, there is a majority shift from very low expression (blue in day 0) to extremely high expression (yellow in day 3 and 7).

## Hi-C & iMARGI

As a result of scRNA-seq, the scientists were able to find which differentially expressed genes should be investigated, leading to the next question– how is this EC stress reflected in genomic interactions? 

In order to answer this question, they first utilized Hi-C (Hierarchical-Chromatin organization), a genomics method used to study the 3D organization of the genome and identify physical interactions between different regions of DNA within the cell nucleus. By cross-linking interacting DNA and associated proteins, then fragmenting and sequencing the DNA strands, Hi-C provides important insight into chromatin organization and gene regulation. For example, topologically associated domains (TADs) which are contiguous and long regions subset of chromosomes that are spatially clustered together and frequently interact. 

![](https://github.com/alicup29/beng183proj/blob/main/img4.png)

As a result of the in situ Hi-C experiments, there were no significant differences between the proportions of read pairs classified as intrachromosomal or genomic-level TAD/A & B compartment changes across all three days. Essentially, the stress conditions “...did not significantly perturb the major 3D genome features…” of the ECs (Calandrelli, et al. 2020).

![](https://github.com/alicup29/beng183proj/blob/main/img5.png)

If a majority of the 3D genome organization in ECs was not affected, then the next question to be asked was– what about changes to RNA-chromatin interactions? In order to tackle this problem, in-situ MApping of RNA-Genome Interaction (iMARGI) was used. iMARGI is used to analyze how RNA molecules physically associate w/ specific regions of DNA at specific genomic loci by crosslinking RNA and local chromatin (DNA), then using NGS to align and determine which RNA interacts with specific chromatin regions.

![](https://github.com/alicup29/beng183proj/blob/main/img6.png)

This process produces RNA-DNA hybrids, or “read pairs”. Between day 0 and day 3 & 7 sample ECs, there was a significant difference in the unique iMARGI read pairs classified as interchromosomal, jumping from 34.4% (day 0) to 62.7% (days 3 & 7)-- indicating that diabetes stress conditions “...induced interchromosomal RNA–chromatin interactions in ECs…” (Calandrelli, et al. 2020).

From here, it becomes a matter of asking what type of interchromosomal interaction causes this EC response. Looking further into the obtained iMARGI read pairs (majority RNA ends over DNA ends), the specific locations are a majority RNA-chromatin interactions at super enhancer regions. Specifically, long-intergenic non-coding RNA (a type of chromatin-associated RNA) from the super enhancer “Linc607SE” are the main culprit, and the start of the super enhancer interaction cascade originating from cell stress and ultimately causing EC dysfunction.

![](https://github.com/alicup29/beng183proj/blob/main/img7.png)

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

## Sources
1. [CHROMATIN - National Human Genome Research Institute](https://www.genome.gov/genetics-glossary/Chromatin#:~:text=00%3A00,fit%20in%20the%20cell%20nucleus)
2. [Targeting Super-Enhancers for Disease Treatment and Diagnosis - Science Direct](https://www.sciencedirect.com/science/article/pii/S1016847823005198)
3. [Stress-induced RNA–chromatin interactions promote endothelial dysfunction - PubMed Central](https://pmc.ncbi.nlm.nih.gov/articles/PMC7566596/)
