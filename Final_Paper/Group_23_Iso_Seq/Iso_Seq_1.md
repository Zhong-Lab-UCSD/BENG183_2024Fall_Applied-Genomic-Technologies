# Isoforms Classification and Quantification 
**Author: Jiarun Liu, Sicheng Jing, Zhijun Qian**

---
## What are RNA Isoforms?
RNA isoforms are variant forms of messenger RNA (mRNA) transcribed from a single gene through mechanisms such as alternative splicing, the use of alternative promoters, or alternative polyadenylation. These different isoforms can ultimately give rise to a diversity of protein products with different functions or regulatory properties, which in turn can further increase the complexity and variety of cellular functions. Generating different RNA isoforms from the same gene can allow organisms to increase their proteomic diversity and regulate gene expression without an increase in genes.

<figure style="text-align: center;">
  <img src="https://raw.githubusercontent.com/njdjyxz/BENG_183_Final/main/Iso_Seq_1/mRNA_isoforms.png" alt="mRNA Isoforms" width="500"/>
  <figcaption>
    <strong>Fig. 1: Alternative Splicing Mechanisms</strong>
    <br/>
    A single gene may be transcribed into several distinct mRNA variants called isoforms through alternative splicing mechanisms. This figure shows six common types of splicing events. (https://www.nature.com/articles/s41467-018-03402-w)
  </figcaption>
</figure>

## Iso-Seq: An Advantage Method in Isoform Detaction
### What is Iso-Seq?
Isoform Sequencing, also called Iso-Seq, is a cutting-edge sequencing technology able to perform a thorough analysis of full-length RNA transcripts with no need for assembly. This technology, developed by Pacific Biosciences, utilizes long-read sequencing so as to efficiently obtain complete isoforms in single reads; thereby providing profound insight into alternative splicing events, transcript diversity, and gene structure. This approach surpasses conventional short-read sequencing techniques in reducing uncertainties in the reconstruction of transcripts, therefore facilitating the recognition of novel isoforms and improving genome annotation. Iso-Seq is critical to advancing our understanding of the regulation and expression of genes, with huge implications in domains like genomics, transcriptomics, and biomedical research.
### Advanced Workflow
<figure style="text-align: center;">
  <img src="https://raw.githubusercontent.com/njdjyxz/BENG_183_Final/main/Iso_Seq_1/Iso_Seq_workflow.jpg" alt="mRNA Isoforms" width="700"/>
  <figcaption>
    <strong>Fig. 2: Iso-Seq Workflow Demo</strong>
    <br/>
    The figure demonstrates the whole complete workflow this chapter will discuss. (Isoform characteristics and Bulk data will not be mentioned)
  </figcaption>
</figure>

#### Experiment Prepration
Similar to standard RNA sequencing protocols, the Iso-Seq workflow begins with the reverse transcription of mRNA samples to generate a complementary DNA (cDNA) library. This essential step converts messenger RNA molecules into stable cDNA, preserving the original sequence information required for accurate downstream analysis. To enhance throughput and cost-effectiveness, especially for single-cell isoform sequencing, the Multiplexed Arrays Sequencing (MAS-Seq) method is employed. MAS-Seq concatenates individual cDNA molecules into longer, ordered arrays, enabling the generation of PacBio Hi-Fi long-read sequences from these concatenated fragments. This approach increases the number of cDNA sequences that can be processed simultaneously, eliminating the need for orthogonal short-read sequencing and thereby streamlining the workflow for high-throughput single-cell isoform analysis.

<figure style="text-align: center;">
  <img src="https://raw.githubusercontent.com/njdjyxz/BENG_183_Final/main/Iso_Seq_1/Mas-Seq.png" alt="mRNA Isoforms" width="500"/>
  <figcaption>
    <strong>Fig. 3: MAS-Seq Library Prep</strong>
  </figcaption>
</figure>

#### Read Segmentation (Skera)
Once the MAS-Seq library is prepared, it undergoes sequencing using PacBio long-read sequencing technologies—the Sequel II/IIe or Revio systems. The concatenated cDNA arrays are sequenced to generate high-fidelity (Hi-Fi) reads, which then undergo bioinformatic processing to recover the original single-cell cDNA sequences. The Skera package is then used to synthesize the long consensus reads into a number of unique sequences in order to assemble the different RNA isoforms present in the sample efficiently. This assembly methodology contributes substantially to the identification and characterization of the different transcript variants, hence providing a broad insight into gene expression and regulation in the biological system being studied. Combining MAS-Seq with PacBio Hi-Fi sequencing increases high-throughput and accurate analysis of single-cell isoforms, expanding our ability to study transcriptomic complexity efficiently and in a cost-effective manner.
<figure style="text-align: center;">
  <img src="https://raw.githubusercontent.com/njdjyxz/BENG_183_Final/main/Iso_Seq_1/Skera.png" alt="mRNA Isoforms" width="500"/>
  <figcaption>
    <strong>Fig. 4: "Skera" Reads Long Concatenated Read Back Into Transcripts</strong>
  </figcaption>
</figure>

Example Code:
~~~Bash
# download HiFi reads for MAS-Seq PBMCs run on Sequel IIe
wget https://downloads.pacbcloud.com/public/dataset/MAS-Seq/DATA-SQ2-PBMC_5kcells/0-CCS/m64476e_220618_014917.hifi_reads.bam

# download MAS adapter fasta
wget https://downloads.pacbcloud.com/public/dataset/MAS-Seq/REF-MAS_adapters/MAS-Seq_Adapter_v1/mas16_primers.fasta

# run skera split to generate segmented reads
skera split m64476e_220618_014917.hifi_reads.bam mas16_primers.fasta segmented.bam
~~~

#### Remove Primers(Lima), polyA tails, CCS, Extract UMI/BCs, Deduplication (IsoSeq3)

<figure style="text-align: center;">
  <img src="https://raw.githubusercontent.com/njdjyxz/BENG_183_Final/main/Iso_Seq_1/Schematic_Workflow.png" alt="mRNA Isoforms" width="600"/>
  <figcaption>
    <strong>Fig. 5: IsoSeq3 Schematic Workflow</strong>
  </figcaption>
</figure>

1. **Primer Removal**
   
   Primers and barcodes are removed using the `lima` tool in `--isoseq` mode. This step ensures that unwanted template-switching oligos (TSOs) are eliminated and that sequences are correctly oriented from 5’ to 3’, regardless of whether the sample is barcoded.

2. **Tagging**
   
   Unique Molecular Identifiers (UMIs) and cell barcodes are clipped from the reads and associated with them using the `isoseq tag` command. This tagging is essential for later deduplication and accurate assignment of reads to their corresponding cells.

3. **Refinement**
   
   The tagged reads undergo refinement, which includes trimming poly(A) tails and removing unintended concatemers. This results in full-length non-concatemer (FLNC) reads, ensuring high-quality data for downstream analysis.

4. **Merge SMRT Cells**
   
   If multiple SMRT cells are used in the sequencing run, their respective `fltnc.bam` files are merged. This consolidation streamlines the data processing and ensures a comprehensive dataset for subsequent steps.

5. **Cell Barcode Correction and Real Cell Identification**
   
   Barcode errors are identified and corrected using a predefined whitelist, ensuring accurate assignment of reads to their respective cells. This step also involves distinguishing real cells from erroneous barcodes, enhancing the reliability of the single-cell data.

6. **Deduplication**
   
   PCR duplicates are removed by clustering reads based on UMIs and corrected cell barcodes using the `isoseq groupdedup` tool. This deduplication process generates consensus sequences for each unique transcript, improving data accuracy and reducing computational requirements.

<figure style="text-align: center;">
  <img src="https://raw.githubusercontent.com/njdjyxz/BENG_183_Final/main/Iso_Seq_1/High-level_Workflow.png" alt="mRNA Isoforms" width="800"/>
  <figcaption>
    <strong>Fig. 6: IsoSeq3 High-Level Workflow with File Type</strong>
  </figcaption>
</figure>

#### Align to Genome (pbmm2/minimap2), Collapse/Merge GFFs(IsoSeq3)
The consensus reads, after deduplication, are aligned to the reference genome using pbmm2, a robust aligner optimized for PacBio HiFi data. This is a very important step in correctly identifying and verifying the presence of various RNA isoforms within the sample. After alignment, Iso-Seq3 was used to assess the collapse of redundant reads to ensure that each identified isoform was uniquely represented. This confirmation assures the integrity and diversity of the transcriptomic data for further downstream analyses, enabling the comprehensive understanding of gene expression patterns.
<figure style="text-align: center;">
  <img src="https://raw.githubusercontent.com/njdjyxz/BENG_183_Final/main/Iso_Seq_1/Collapse.png" alt="mRNA Isoforms" width="600"/>
  <figcaption>
    <strong>Fig. 6: IsoSeq3 High-Level Workflow with File Type</strong>
  </figcaption>
</figure>

Sample Code:

~~~Bash
# Map reads using pbmm2 before collapsing
pbmm2 align --preset ISOSEQ --sort <input.bam> <ref.fa> <mapped.bam>

# Collapse mapped reads into unique isoforms using isoseq collapse.
isoseq collapse <mapped.bam> <collapse.gff>

#################################################################################
# Main Output File: collapse.gff contains the collapsed isoforms in gff format. #
#################################################################################
~~~

