# MERFISH: Multiplexed error-robust fluorescence in situ hybridization
### By Nick Monell, Benjamin Esser and Tusha Karnani

## Contents
* [What is MERFISH?](#what-is-merfish)
* [Procedure](#procedure)
* [Computational Decoding](#computational-decoding)
* [Output](#output)
* [Applications of MERFISH](#applications-of-merfish)
* [Conclusion](#conclusion)
* [References](#references)

## What is MERFISH?
#### FISH
#### MERFISH Overview
![](https://github.com/tkarnani/BENG183_2024Fall_Applied-Genomic-Technologies/blob/main/Final_Paper/Group_9_MERFISH/Images/merfish1.png)

## Computational Decoding
MERFISH utilizes a high-throughput encoding system to differentiate between hundreds to thousands of RNA species in a single experiment. The encoding system is as follows:

* We will assign an 'n-bit' binary string to each RNA species where n is the number of imaging rounds in our procedure. Each round will correspond to either a 1 (fluorescence detected) or 0 (no fluorescence) on an RNA based on if a fluorescently labeled probe is attached to that RNA. For example, if we wanted to encode the gene ACTC1 and used 12 imaging rounds, we would use a string of twelve 1's and 0's that may look like `010010101000`, indicating that we expect a fluorescent readout in imaging rounds 2, 5, 7, and 9 at every spatial position where the ACTC1 gene is located
* Important to note, RNA may be scattered everywhere on the surface we are imaging. If the goal is to determine both the count and spatial location of many different types of RNA, we must first locate where each RNA is (generally) and then identify which species is in each position. Locating the general RNA position is as easy as determining wherever there is a fluorescent peak in any imaging round. Determining which RNA it is depends on what readout we see over all imaging rounds for that position.
  
<div align="center">
<img src="https://raw.githubusercontent.com/tkarnani/BENG183_2024Fall_Applied-Genomic-Technologies/main/Final_Paper/Group_9_MERFISH/Images/encoding.jpg" width="45%" style="display: block; margin: auto;"/>

Figure 2: Demonstration of detected RNA spots (general) and how to decode each spot into an RNA species based on fluorescent readout over 16 hybridization rounds. Kok Hao Chen et al.,Spatially resolved, highly multiplexed RNA profiling in single cells. <i>Science</i> <b>348</b>, aaa6090 (2015). <a href="https://doi.org/10.1126/science.aaa6090">DOI:10.1126/science.aaa6090</a>
</div>

* Considering all possible n-length binary strings, we can encode $2^n - 1$ gene species using this system (removing 1 for the string containing only 0's as this would be an undetectable position). This makes it highly scalable for multiplexing. This means the minimum number of imaging rounds we must perform, n, is $\log_2(x + 1)$ where x is the number of gene species we are interested in.
* For our encoding system, we will apply some "rule" to regulate which binary strings we can use to encode RNA species. Remember, MERFISH is Error Robust meaning if a probe doesn't bind properly in a particular hybridization round or fluorescent readout is unresolvable, we should still be able to determine which RNA species was there based on the rest of our code. One such system is Modified Hamming Distance 4 (MDH4), consisting of n-bit binary strings with a maximum of 4 'fluorescence-on' (1) bits that are at least different from any other string by 4-bit positions. Valid codes would be `000000001111` or `100100100100` while an invalid code may be `110001101010`. Because we know that each code is supposed to have 4 on-bits if there is any single-bit error we can correct it by referencing the readout with our codebook, finding which string position the error is in, and decoding which RNA species it is. 
The benefit of encoding and decoding in MERFISH is the redundancy and error correction, allowing us to resolve the count and position of every RNA species even in extensive datasets.
## Validation
The results of decoding should be validated to ensure we are confident in the RNA counts for each species.
#### Control Threshold
Within the encoding scheme designed for the experiment, encode several 'control words', or n-bit strings that follow the same encoding rules (MDH4, Reed-Muller, BCH, etc.) but do not encode for any RNA species. This means that we do not have any designed probe sequence for these codes, so their respective fluorescent readout should not appear in the data or we have made a misidentification. Generate a confidence score for each decoded readout, the maximum confidence for our control words is the confidence threshold we will use to flag other RNA detection. Any spots that have a confidence lower than this threshold are unreliable as they achieve the same confidence as known errors. 

<div align="center">
<img src="https://raw.githubusercontent.com/tkarnani/BENG183_2024Fall_Applied-Genomic-Technologies/main/Final_Paper/Group_9_MERFISH/Images/validation.jpg" width="30%" style="display: block; margin: auto;"/>

Figure 3: Comparison of confidence ratio between detected RNA (blue) and control words (red) with dashed line representing maximum confidence ratio of control words. Kok Hao Chen et al.,Spatially resolved, highly multiplexed RNA profiling in single cells. <i>Science</i> <b>348</b>, aaa6090 (2015). <a href="https://doi.org/10.1126/science.aaa6090">DOI:10.1126/science.aaa6090</a>
</div>

#### RNA-Seq Reference
Use the counts for each RNA species as a form of validation. Compare the actual RNA counts for each species with a known method of RNA-seq. If copy numbers show a strong correlation, we have strong evidence of the method's accuracy.

## Output
#### Expression Noise
#### Expression Covariation
#### Spatial Distribution
Some RNA transcripts enriched in the perinuclear region, some enriched in the cell periphery, and some scattered throughout the cell.
It determines the correlation coefficients for the spatial density profiles of all pairs of RNA species and organized these RNAs according to the pairwise correlations again using a hierarchical clustering approach.
The spatial pattern that observed reflects their cotranslational enrichment at the ER since they pass through the same/similar secretion pathways.
![](https://github.com/tkarnani/BENG183_2024Fall_Applied-Genomic-Technologies/tree/main/Final_Paper/Group_9_MERFISH/Images/spatial.jpeg)


## Applications of MERFISH
MERFISH (Multiplexed Error-Robust Fluorescence In Situ Hybridization) has revolutionized spatial transcriptomics by enabling the high-throughput and spatially resolved analysis of gene expression. Its ability to detect thousands of RNA species while preserving spatial context has found applications across various fields of biology and medicine.
#### Developmental Biology
- **Tissue Morphogenesis**: MERFISH reveals how gene expression drives the formation and differentiation of tissues during embryonic development.
- **Cell Lineage Tracing**: It enables the study of how single cells contribute to tissue formation, providing a spatial view of developmental trajectories.
#### Cancer Research
- **Tumor Microenvironment**: MERFISH maps gene expression in tumor and surrounding stromal cells, providing a spatial understanding of the tumor microenvironment.
It identifies interactions between cancer cells, immune cells, and other stromal components.
- **Understanding Metastasis**: By spatially resolving gene expression in metastatic tumors, MERFISH reveals pathways that cancer cells use to invade new tissues.
- **Therapeutic Targeting**: It helps pinpoint specific cell types or molecular pathways for targeted therapies, such as immune checkpoint inhibitors.
#### Infectious Diseases and Immunology
- **Host-Pathogen Interactions**: It helps study how pathogens interact with host cells, revealing spatial patterns of immune activation or suppression.
- **Response to Vaccines**: It helps assess immune responses to vaccines in different tissues by tracking spatial gene expression changes.

## Conclusion
MERFISH exemplifies the evolution of biological research tools, moving beyond standard sequencing to address not just the "what" but also the "where" of gene expression. This spatial resolution allows researchers to explore how cells interact within their native environments and how these interactions change during development, disease, or treatment.

## References
1. [Kok Hao Chen et al., Spatially resolved, highly multiplexed RNA profiling in single cells, Science](https://www.science.org/doi/10.1126/science.aaa6090)
