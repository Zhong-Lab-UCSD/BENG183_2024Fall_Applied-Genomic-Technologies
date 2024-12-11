# Isoforms Classification and Quantification 
**Author: Jiarun Liu, Sicheng Jing, Zhijun Qian**

---
## Table of Contents

1. [What are RNA Isoforms?](#what-are-rna-isoforms)
2. [Iso-Seq: An Advantage Method in Isoform Detaction](#iso-seq-an-advantage-method-in-isoform-detaction)
   1. [What is Iso-Seq?](#what-is-iso-seq)
   2. [Advanced Workflow](#advanced-workflow)
      1. [Experiment Preparation](#experiment-preparation)
      2. [Read Segmentation (Skera)](#read-segmentation-skera)
      3. [Remove Primers(Lima), polyA tails, CCS, Extract UMI/BCs, Deduplication (IsoSeq3)](#remove-primerslima-polya-tails-ccs-extract-umibcs-deduplication-isoseq3)
      4. [Align to Genome (pbmm2/minimap2), Collapse/Merge GFFs(IsoSeq3)](#align-to-genome-pbmm2minimap2-collapsemerge-gffsisoseq3)
      5. [Isoform Classification and Quantifying Isoform-Specific Expression (SQANTI3/IsoQuant)](#isoform-classification-and-quantifying-isoform-specific-expression-sqanti3isoquant)
          1. [Why We Need to Classify Isoforms?](#why-we-need-to-classify-isoforms)
          2. [SQANTI3](#sqanti3)
          3. [IsoQuant](#isoquant)
      4. [Finding Special Genes (IGV)](#finding-special-genes-igv)
      5. [Alternative Polyadenylation (APA) Analysis](#alternative-polyadenylation-apa-analysis)
6. [Application of Isoseq on Micronuclei](#application-of-isoseq-on-micronuclei)
    1. [Understanding Micronuclei](#understanding-micronuclei)
    2. [Mechanisms of Micronuclei Formation](#mechanisms-of-micronuclei-formation)
    3. [Micronuclei and Cancer](#micronuclei-and-cancer)
10. [Source](#source)

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
  <img src="https://raw.githubusercontent.com/njdjyxz/BENG_183_Final/main/Iso_Seq_1/Iso_Seq_workflow.png" alt="mRNA Isoforms" width="700"/>
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

---

<!-- ## Isoform Analysis

After obtaining isoform information from the data pre-processing stage, we can proceed with various downstream analyses to uncover hidden insights and underlying relationships. -->

### Isoform Classification and Quantifying Isoform-Specific Expression (SQANTI3/IsoQuant)

#### Why we need to classify isoforms?

Classifying isoforms is essential in genomics and transcriptomics because it helps unravel the complexity of gene expression and its regulation. Different isoforms may have varying levels of expression across samples or conditions. Classifying isoforms enables the quantification of their expression, which is crucial for studying differential expression and identifying biomarkers.

#### SQANTI3

**SQANTI3** (Structural and Quality Annotation of Novel Transcript Isoforms 3) is a comprehensive tool designed for the annotation and quality assessment of transcript isoforms obtained from long-read sequencing technologies, such as PacBio Iso-Seq or Oxford Nanopore. It is particularly useful for analyzing full-length transcripts and assessing their structural accuracy and functional relevance.

<div style="display: flex; flex-direction: column; align-items: center; text-align: center;">
  <img src="https://raw.githubusercontent.com/njdjyxz/BENG_183_Final/main/Iso_Seq_2/SQANTI3.png" width="600"/>
  <figcaption><strong>Fig. 7: SQANTI3 overview</strong></figcaption>
</div>

To get the classification information, running sqanti3_qc.py is enough. Here is the breifly instruction:

```shell
python sqanti3_qc.py  sample1.collapsed.sorted.gff \
		gencode.v39.annotation.sorted.gtf \
		human_GRCh38_no_alt_analysis_set.fasta \
		-t 20 -fl sample1.collapsed.abundance.txt -d $outdir --CAGE_peak \
		refTSS_v3.3_human_coordinate.hg38.sorted.bed --polyA_motif_list \
		polyA.list.txt --report both --isoAnnotLite                    
```
Positional arguments:
1. Isoforms (FASTA/FASTQ) or GTF format. It is recommended to provide them in GTF format, but if it is needed to map the sequences to the genome use a FASTA/FASTQ file with the --fasta option. The obtained `*collapsed.gff`is the output from the `isoseq collapse`.

2. Reference annotation file (GTF format)

3. Reference genome (Fasta format)

As a result of running SQANTI3 QC, the tool will create a series of output files in the specified directory (`--dir` or `-d` flag in the QC script), and the classification results are contained in `_classification.txt`. It is a tab-separated file where transcripts are rows and QC attributes are columns. Isoforms are identified by their ID in the input long-read GTF file (`*collapsed.gff`). 


<figure style="text-align: center;">
  <img src="https://raw.githubusercontent.com/njdjyxz/BENG_183_Final/main/Iso_Seq_2/classification_file.png" width="1000"/>
  <figcaption>
    <strong>Fig. 8: Sample of clssification output file</strong>
  </figcaption>
</figure>

**Columns 1**:`isoform`: the isoform ID. Usually in `PB.X.Y` format. `PB.X`refers to the gene ID, and `.Y`means the number of transcript belongs to this gene.

Columns 2:`chrom`: chromosome.

Columns 3:`strand`: strand.

Columns 4:`length`: isoform length.

**Columns 6**:`structural_category`: one of the categories ["full-splice_match", "incomplete-splice_match", "novel_in_catalog", "novel_not_in_catalog", "genic", "antisense", "fusion", "intergenic", "genic_intron"]

**Columns 15**:`subcategory`: additional splicing categorization, separated by semi-colons.

We need to refer to the official definitions of these categories and subcategories to gain a complete understanding of the results:

<figure style="text-align: center;">
  <img src="https://raw.githubusercontent.com/njdjyxz/BENG_183_Final/main/Iso_Seq_2/figuras_class_SQ3.png" width="600"/>
  <figcaption>
    <strong>Fig. 9: Isoform categories</strong>
  </figcaption>
</figure>

- **FSM (Full Splice Match)**: meaning the reference and query isoform have the same number of exons and each internal junction agree. The exact 5' start and 3' end can differ by any amount.
- **ISM (Incomplete Splice Match)**: the query isoform has fewer 5' exons than the reference, but each internal junction agree. The exact 5' start and 3' end can differ by any amount.
- **NIC (Novel In Catalog)**: the query isoform does not have a FSM or ISM match, but is using a combination of known donor/acceptor sites.
- **NNC (Novel Not in Catalog)**: the query isoform does not have a FSM or ISM match, and has at least one donor or acceptor site that is not annotated.
- **Antisense**: the query isoform does not have overlap a same-strand reference gene but is anti-sense to an annotated gene.
- **Genic Intron**: the query isoform is completely contained within an annotated intron.
- **Genic Genomic**: the query isoform overlaps with introns and exons.
- **Intergenic**: the query isoform is in the intergenic region.

Some categories also include subcategories, as outlined below:

<figure style="text-align: center;">
  <img src="https://raw.githubusercontent.com/njdjyxz/BENG_183_Final/main/Iso_Seq_2/figure_fsm_subcat_SQ3.png" width="600"/>
  <figcaption>
    <strong>Fig. 10: Isoform subcategories(FSM)</strong>
  </figcaption>
</figure>

<figure style="text-align: center;">
  <img src="https://raw.githubusercontent.com/njdjyxz/BENG_183_Final/main/Iso_Seq_2/figure_ism_subcat_SQ3.png" width="600"/>
  <figcaption>
    <strong>Fig. 11: Isoform subcategories(ISM)</strong>
  </figcaption>
</figure>

<figure style="text-align: center;">
  <img src="https://raw.githubusercontent.com/njdjyxz/BENG_183_Final/main/Iso_Seq_2/figure_nic_nnc_subcat_SQ3.png" width="600"/>
  <figcaption>
    <strong>Fig. 12: Isoform subcategories(NIC and NNC)</strong>
  </figcaption>
</figure>

After classifying each transcript, the next step is to quantify the number of reads aligned to each transcript for each sample (assuming we have four samples, each containing a different number of cells). Ultimately, this process provides the distribution of reads across various categories (and subcategories) within each sample. It is important to note that categories or subcategories such as **Genic Intron**, **Intergenic**, **mono-exon**, etc., are not considered transcripts. Therefore, it is advisable to remove rows belonging to these categories before proceeding with further analysis, especially if the file size is very large.

#### IsoQuant

**IsoQuant** is a software tool designed for the quantification, error correction, and functional annotation of transcript isoforms derived from long-read sequencing platforms such as PacBio Iso-Seq or Oxford Nanopore. It is widely used in transcriptomics studies to analyze complex transcriptomes, resolve isoform diversity, and enhance downstream functional insights. 

<figure style="text-align: center;">
  <img src="https://raw.githubusercontent.com/njdjyxz/BENG_183_Final/main/Iso_Seq_2/IsoQuant.png" width="600"/>
  <figcaption>
    <strong>Fig. 13: IsoQuant pipline</strong>
  </figcaption>
</figure>

 Here is the breifly instruction:

```shell
python isoquant.py --bam $fileBam \
                  --reference $reffile \
                  --genedb $gtffile \
                  -d pacbio_ccs \
                  --read_group tag:CB \
                  --output $output_dir \
                  --threads 20  \
                  --check_canonical \
                  --count_exons \
                  --prefix ${prefix} \
                  --no_model_construction
```

Positional arguments:

1. Sorted and indexed BAM file
2. Reference sequence in FASTA format (can be gzipped);
3. *Optionally*, you may provide a reference gene annotation in GTF/GFF format (can be gzipped).

The process generates numerous output files, many of which depend on the options specified by the user in the command line. For transcript data, the final results are stored in `SAMPLE_ID.transcript_counts.tsv`, a TSV file containing raw read counts for each reference transcript.

<figure style="text-align: center;">
  <img src="https://raw.githubusercontent.com/njdjyxz/BENG_183_Final/main/Iso_Seq_2/transcriptCount.png" width="600"/>
  <figcaption>
    <strong>Fig. 14: Sample of trascript count tsv file</strong>
  </figcaption>
</figure>

Ultimately, we can combine the two output files (`_classification.txt` and `SAMPLE_ID.transcript_counts.tsv`) to generate the final statistics. The following code demonstrates how to achieve this:

```shell
#--------------------Sample 1-----------------------
# Calculate the sum of each row and remove rows where the result is 0.
# Here, the sum represents the number of reads aligned to a specific transcript in sample 1.
# If the sum is 0, it means that no reads in sample 1 are associated with this transcript.
awk 'NR==1 {next} {sum=0; for(i=2; i<=NF; i++) sum+=$i; if(sum != 0) print $1, sum}' SAMPLE_ID.transcript_counts.tsv > sample1_readIsoClass.tsv

# For each retained row, use the PB number to look up the corresponding classification in the `_classification.txt` file.
# Include the PB number, the sum from the previous step, and the associated category in the new output file.
awk 'FNR==NR {b[$1]=$6; next} {
    if ($1 in b) {
        new_col = b[$1];
        printf "%s\t%s\t%s", $1, $2, new_col;
        printf "\n";
    } 
}' _classification.txt sample1_readIsoClass.tsv > sample1_readIsoClass_withClass.tsv

# Combine the sums for rows that belong to the same category.
awk '{sum[$3]+=$2} END {for (category in sum) print category, sum[category]}' sample1_readIsoClass_withClass.tsv > sample1_classNum.tsv
```

Assuming we have 4 samples, the same actions need to be performed on each sample. Below is the final output along with the corresponding visualization plot.

<figure style="text-align: center;">
  <img src="https://raw.githubusercontent.com/njdjyxz/BENG_183_Final/main/Iso_Seq_2/result.png" width="1000"/>
  <figcaption>
    <strong>Fig. 15: Results of isoforms classification and quantification</strong>
  </figcaption>
</figure>

Based on these results, we can draw meaningful conclusions, such as sample 1 having more NNC than the other samples. With this, we have successfully completed a comprehensive workflow for isoform classification and quantification!

### Finding special genes (IGV)

Many genes produce multiple isoforms through alternative splicing or other regulatory mechanisms. Examining these isoforms helps identify functional differences, as different isoforms may play distinct roles in cellular processes. One of the most straightforward ways to explore isoform differences between samples is by using IGV (Integrative Genomics Viewer). Below is an example involving the gene *SRSF11*, which exhibits unique isoforms in sample 1 that are absent in the other three samples. By identifying specific genes of interest, we can conduct further analyses to uncover underlying mechanisms and gain deeper insights into their biological significance.

<figure style="text-align: center;">
  <img src="https://raw.githubusercontent.com/njdjyxz/BENG_183_Final/main/Iso_Seq_2/igv.png" width="600"/>
  <figcaption>
    <strong>Fig. 16: Example of gene transcript that differs between samples through igv</strong>
  </figcaption>
</figure>

### Alternative Polyadenylation (APA) Analysis

Isoforms generated through APA events can influence transcript stability, localization, and translation. Classifying these isoforms helps in understanding post-transcriptional regulation.

---
## Application of Isoseq on Micronuclei
From the above discussion we know that IsoSeq enables the study of full-length transcripts, providing valuable insights into RNA splicing and transcriptional disruptions.      
Due to the specific traits of micronuclei, IsoSeq becomes a perfect tool to study micronuclei 

### Understanding Micronuclei
**What are Micronuclei?**   
Micronuclei have a named as "the Lost Chromosomes". It was first discovered in 19th centure by von Hansemann.   
Micronuclei are cytoplasmic structures that contain entire chromosomes or chromosomal fragments. They usually result from errors in segregation during mitosis.   

**Why Study Micronuclei?**   
Micronuclei are key indicators in cancer biology research.     
The exposure of micronuclear chromosome to the cytosolic milieu, through rupture and collapse of the micronuclear envelope during interphase, can raise several pathologic responses which drive tumor progression.

### Mechanisms of Micronuclei Formation
**How did micronulei being produced?**   
Micronuclei arise from mis-segregating events where an entire or a part of a chromosome is left outside of the primary nucleus. And then, the micronucleus ends randomly in one of the daughter cells.  

**The Various Fate of a Micronucleus**
![Various fate of a micronucleus](https://raw.githubusercontent.com/njdjyxz/BENG_183_Final/main/Iso_Seq_3/fate%20of%20micronuclei.png)
**Fig. 8: Alternative splicing mechanisms**  

a. The micronucleus be reincorporated into the primary nucleus through microtubules correct attachment in mitosis, giving rise to two daughter cells, one aneuploid and one diploid.   
b. The micronucleus  persist in a micronuclear state by mis-segregating again in the subsequent mitosis, forming a diploid daughter cell but micronucleated and an aneuploid daughter cell.   
c. The micronucleus undergo replication and persist in a micronuclear state through mis-segregation, forming two micronucleated, diploid daughter cells.   
d. The micronucleus be extruded from the cell, forming after mitotic division two aneuploid daughter cells that are both missing the micronucleated chromosome.  
[Source](https://pubmed.ncbi.nlm.nih.gov/38197599/#&gid=article-figures&pid=figure-1-uid-0)   

**Micronuclei functions as the markers of DNA damage. Compared to mitoric figures, its presence  is easier to detect in tissues. As a result, microneclei is a best case to study ongoing mis-segregation. And work as a reliable marker of chromosomal instability.**

### Micronuclei and Cancer
**Micronuceli Rupture**  
Micronuclear rupture is a fundamental step in cancer progression. It can activate immune signaling pathways.    

**Fig. 9: Genomic consequences of micronuclear rupture**  
![result of Micronuceli rupture](https://raw.githubusercontent.com/njdjyxz/BENG_183_Final/refs/heads/main/Iso_Seq_3/rupture.png)
Upon nuclear envelope collapse, the genetic material contained in the micronucleus undergo profound changes and rearrangements.
   
       
          
**Fig. 10: Inflammatory consequences of micronuclear rupture**   
![cGas-STING](https://raw.githubusercontent.com/njdjyxz/BENG_183_Final/refs/heads/main/Iso_Seq_3/cGas-STING.png)
Due to micronuclear envelope collapse, cytosolic exposure of genetic material initiate the nucleic-acid recognition, and the innate immune cGAS–STING pathway and downstream inflammatory signaling.  
[Source](https://pubmed.ncbi.nlm.nih.gov/38197599/#&gid=article-figures&pid=figure-3-uid-2)

**Therefore, micronuclei have emerged as a prognostic marker in patients with cancer and as a potential predictor of therapeutic outcome.**  
**By identifying unique splicing events and RNA modifications, IsoSeq helps uncover the molecular mechanisms driven by micronuclei.**

## Source   
1. https://github.com/Magdoll/cDNA_Cupcake/wiki/Iso-Seq-Single-Cell-Analysis:-Recommended-Analysis-Guidelines
2. https://pacbio.cn/wp-content/uploads/Technical-overview-MAS-Seq-library-preparation-using-MAS-Seq-for-10x-Single-Cell-3-kit.pdf
3. https://www.pacb.com/wp-content/uploads/Application-note-MAS-Seq-for-single-cell-isoform-sequencing.pdf
4. https://skera.how/
5. https://isoseq.how/
6. https://www.nature.com/articles/s41467-018-03402-w
7. https://github.com/ConesaLab/SQANTI3
8. https://github.com/ablab/IsoQuant?tab=readme-ov-file
9. Di Bona M, Bakhoum SF. *Micronuclei and Cancer.* Cancer Discov. 2024 Feb 8;14(2):214-226. doi: [10.1158/2159-8290.CD-23-1073](https://doi.org/10.1158/2159-8290.CD-23-1073). PMID: 38197599; PMCID: PMC11265298.
10. https://pubmed.ncbi.nlm.nih.gov/38197599/#&gid=article-figures&pid=figure-1-uid-0
11. https://pubmed.ncbi.nlm.nih.gov/38197599/#&gid=article-figures&pid=figure-3-uid-2
12. https://pubmed.ncbi.nlm.nih.gov/38197599/#&gid=article-figures&pid=figure-4-uid-3