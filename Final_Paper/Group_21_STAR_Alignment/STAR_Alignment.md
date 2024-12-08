# **STAR: A Brief Introduction to Its Algorithm and Usage**
### Name: Yang Han, Zhaogu Sun, Simona Wu (Team 21)
## **Background and Introduction:** 
RNA sequencing (RNA-Seq) is prevalently used for analyzing gene expression and transcriptomic profiles. It involves several key computational steps to process raw sequencing data, with read alignment to a reference genome being a critical step. STAR, which stands for "Spliced Transcripts Alignment to a Reference", is a commonly used tool for this RNA-seq. STAR is known for its ability to handle large RNA-Seq datasets.

STAR operates on Linux via the command line, aligning trimmed .fastq files to an indexed reference genome. The tool produces output files in .sam or .bam formats, which contain the aligned reads necessary for downstream analyses, such as gene quantification and splice junction identification.

The RNA-Seq workflow begins with biological sample preparation and sequencing, followed by quality control steps like adapter trimming and sequence evaluation using tools such as FASTQC. STAR facilitates the alignment process, accurately mapping reads to the genome while accounting for splicing events. Once the reads are aligned, they are quantified to associate them with genes, enabling statistical analysis to identify differentially expressed genes.

STAR’s alignment capability supports essential steps in RNA-Seq workflows, ensuring the accurate representation of transcriptomic data for downstream interpretation.
![](https://github.com/TonyYangHan/BENG183_2024Fall_Applied-Genomic-Technologies/blob/main/Final_Paper/Group_21_STAR_Alignment/Graphs/RNAseq_Workflow.png)

## **STAR Alignment Overview**
### **Seed-and-Extension Strategy in STAR Alignment**
STAR employs a "seed-and-extension" strategy to efficiently map sequencing reads to a reference genome. This approach balances computational efficiency with alignment accuracy, making it powerful for RNA-Seq data analysis.

### **Seed Phase**
The first phase of the seed-and-extension strategy involves identifying short sequences, called "seeds," within the sequencing reads. These seeds are chosen because they exactly match subsequences in the reference genome. By focusing on exact matches, STAR can utilize efficient string search algorithms such as the Burrows-Wheeler Transform (BWT) or Suffix Arrays to rapidly locate potential regions of alignment. This step significantly reduces the search space, narrowing the possible locations where the full read might align.

#### **Advantages of Seeds:**
- Efficiency: Exact matching enables the use of precomputed indices like the Suffix Array or BWT, speeding up the search process.

- Scalability: This approach allows the alignment algorithm to handle large genomes and datasets.

#### **Disadvantages of Seeds:**
- Memory Intensity: Building and storing indices for the reference genome can require substantial computational resources.

- High Computational Cost: Despite the efficiency of exact match algorithms, the sheer volume of reads in high-throughput sequencing datasets can make this step computationally expensive.

### **Extension Phase**
Once potential seed locations are identified, the algorithm extends these matches to align the entire read. This involves dealing with mismatches, insertions, deletions, and splicing events that may occur within the read. STAR is especially good at this phase by:

- Allowing for gaps in alignments to accommodate introns, which is particularly important for RNA-Seq data.

- Using scoring systems to evaluate alignment quality and select the best match for each read.

## **STAR Algorithm Illustration**
### Step1: Finding the longest prefixes in reads that exactly matches the reference genome
![](https://github.com/TonyYangHan/BENG183_2024Fall_Applied-Genomic-Technologies/blob/main/Final_Paper/Group_21_STAR_Alignment/Graphs/Step1.png)

STAR starts by finding the longest prefix in reads that exactly match some (>= 1) locations in the reference genome. The current longest prefix in a read is called a Maximal Mapping Prefix. As stated previously, STAR uses uncompressed suffix arrays to conduct string search in the reference genome. This reduces the runtime of aligning a sequence from O(n^2) to at most O(n*log2(n)) and greatly reduces the memory usage from O(n^2) to O(n). By statistical probability, longer sequences tend to map to fewer locations in the reference genome, which offers less genome coordinates to look at at later steps.  Additionally, by doing this, STAR minimizes the unmapped regions in reads, which saves some effort for later alignment tasks.

### Step 2: Mapping the the unmapped regions to the reference genome
#### Scenario 1: Able to find exact match in the reference genome for the unmapped portion: 
![](https://github.com/TonyYangHan/BENG183_2024Fall_Applied-Genomic-Technologies/blob/main/Final_Paper/Group_21_STAR_Alignment/Graphs/Step2_Scenario1.png)

- If STAR is able to find some exact match for the entire unmapped portion at some genome locations after the first MMP, then we are done. However, in many cases, we are not so lucky.

#### Scenario 2: Minor mismatch between unmapped region and the reference genome
![](https://github.com/TonyYangHan/BENG183_2024Fall_Applied-Genomic-Technologies/blob/main/Final_Paper/Group_21_STAR_Alignment/Graphs/Step2_Scenario2.png)

- As STAR is mapping the second MMP (second longest exact match prefix) to some locations in the reference genome, STAR encountered some mismatch that prevented us from further matching to extend the MMP. In this case, STAR will count how many mismatches STAR have seen while STAR continues our effort to match the MMP to the reference genome. If the mismatch is short and STAR reaches the next exact matching region, STAR will extend the second MMP to include the mismatch and the next MMP. Such mismatch will be reported in the final output of STAR.

#### Scenario 3: Major mismatch occurs
![](https://github.com/TonyYangHan/BENG183_2024Fall_Applied-Genomic-Technologies/blob/main/Final_Paper/Group_21_STAR_Alignment/Graphs/Step2_Scenario3.png)

- If STAR encounters a significantly mismatched region or regions with poor sequencing quality that decreases the alignment quality considerably, then such region is soft-clipped and will not be included in the aligned reads. However, the location and length of the soft-clipped region is still reported in the final output. 

Due to length limitations, this paper will not demonstrate the details of the scoring system of STAR during alignment.

Step 2 will be repeated until STAR has matched all portions in all mappable reads to the reference genome. STAR will then proceed to the next stage to stitch together mapped portions of reads.

## Clustering, Stitching, Scoring:
![](https://github.com/TonyYangHan/BENG183_2024Fall_Applied-Genomic-Technologies/blob/main/Final_Paper/Group_21_STAR_Alignment/Graphs/Step3_StitchingReads.png)

- After completing the seed searching steps, STAR proceeds with the clustering, stitching, and scoring processes to reconstruct full, accurate reads from fragmented or partially aligned sequence segments. The pipeline begins by identifying uniquely mapped seeds, which are designated as anchor points for subsequent alignment. These anchor points serve as central hubs around which neighboring seeds are clustered, enabling the identification of fragments that likely originate from the same read. 

- Once the clustering is complete, the next step involves stitching the clustered seeds to reconstruct a continuous read. This is achieved by evaluating potential alignments based on a scoring system that considers various factors, including insertions and deletions (indels), skipped regions, soft-clipped regions, matches, and mismatches. The scoring process aims to identify the alignment that most accurately represents the original sequencing read while minimizing errors.

## Usage
### 1. Building Index for Reference Genome:
```
$STAR --runThreadN 16 --runMode genomeGenerate --genomeDir Your_directory/ --genomeFastaFiles reference.fa --sjdbGTFfile gtf_file.gtf --sjdbOverhang 100
```

#### Explanation: 
- --runMode genomeGenerate：Specifictly indicate to STAR that the current goal is to do genome index generation
- --genomeFastaFiles：Point to the Fasta file that contains the reference genome
- --sjdbGTFfile: point out the path to the GTF file, which contains annotated genes and their locations in the genome. This annotation helps STAR recognize known splice junctions.
- --sjdbOverhang: Set the length of the overhang for splice junctions. Here we use 100 as an example.

### 2. Aligning Reads to the Reference Genome Using STAR:
```
$STAR --runThreadN 16 --genomeDir chrX_STAR_index/ --readFilesIn inputFile.fastq/ --outSAMtype BAM SortedByCoordinate --outFileNamePrefix ./STAR
```

#### Explanation: 
- --genomeDir: Specifies the directory containing the indexed reference genome
- --outSAMtype: Tell STAR to output the alignment result in which format

## Expected Output
1. Prefix_Aligned.sortedByCoord.out.bam: This file contains the aligned reads in BAM format, organized by their genomic coordinates. It serves as the primary output of the alignment.

2. Prefix_Log.out: This log file provides general details about the alignment process, such as the total number of reads processed, the count of uniquely aligned reads, and the number of reads mapped to multiple locations.

3. Prefix_Log.final.out: This file gives a summary of the alignment statistics, providing a quick snapshot of the alignment overall success.

4. Prefix_Log.progress.out: This file tracks the alignment progress in real-time, which is especially helpful for monitoring longer runs.

5. Prefix_SJ.out.tab: This file lists the splice junctions identified during the alignment, which can be useful for identifying novel splicing events or for comparing with known annotations.

## Sources
1. [STAR: ultrafast universal RNA-seq aligner](https://pmc.ncbi.nlm.nih.gov/articles/PMC3530905/pdf/bts635.pdf)
2. [STAR Manual](https://github.com/alexdobin/STAR/blob/master/doc/STARmanual.pdf)
3. [Introduction to RNA-Seq using high-performance computing](https://hbctraining.github.io/Intro-to-rnaseq-hpc-O2/lessons/03_alignment.html)












