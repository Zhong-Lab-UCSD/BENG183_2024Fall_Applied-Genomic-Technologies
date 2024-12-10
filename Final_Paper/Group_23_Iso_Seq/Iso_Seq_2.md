## Isoform Analysis

After obtaining isoform information from the data pre-processing stage, we can proceed with various downstream analyses to uncover hidden insights and underlying relationships.

### 1.Isoform Classification and Quantifying Isoform-Specific Expression

#### 1.1 Why we need to classify isoforms?

Classifying isoforms is essential in genomics and transcriptomics because it helps unravel the complexity of gene expression and its regulation. Different isoforms may have varying levels of expression across samples or conditions. Classifying isoforms enables the quantification of their expression, which is crucial for studying differential expression and identifying biomarkers.

#### 1.2 SQANTI3

**SQANTI3** (Structural and Quality Annotation of Novel Transcript Isoforms 3) is a comprehensive tool designed for the annotation and quality assessment of transcript isoforms obtained from long-read sequencing technologies, such as PacBio Iso-Seq or Oxford Nanopore. It is particularly useful for analyzing full-length transcripts and assessing their structural accuracy and functional relevance.

<figure style="text-align: center;">
  <img src="https://raw.githubusercontent.com/njdjyxz/BENG_183_Final/main/Iso_Seq_2/SQANTI3.png" width="600"/>
  <figcaption>
    <strong>Fig. 7: SQANTI3 overview</strong>
  </figcaption>
</figure>


To get the classification information, running sqanti3_qc.py is enough. Here is the breifly instruction:

```shell
python sqanti3_qc.py  sample1.collapsed.sorted.gff \
											gencode.v39.annotation.sorted.gtf \
											human_GRCh38_no_alt_analysis_set.fasta \
							-t 20 -fl sample1.collapsed.abundance.txt -d $outdir --CAGE_peak
							refTSS_v3.3_human_coordinate.hg38.sorted.bed --polyA_motif_list
							polyA.list.txt --report both --isoAnnotLite                    
```

Positional arguments:

1. isoforms
   Isoforms (FASTA/FASTQ) or GTF format. It is recommended to provide them in GTF format, but if it is needed to map the sequences to the genome use a FASTA/FASTQ file with the --fasta option. The obtained `*collapsed.gff`is the output from the `isoseq collapse`.

2. annotation 
   Reference annotation file (GTF format)

3. genome

   Reference genome (Fasta format)

As a result of running SQANTI3 QC, the tool will create a series of output files in the specified directory (`--dir` or `-d` flag in the QC script), and the classification results are contained in `_classification.txt`. It is a tab-separated file where transcripts are rows and QC attributes are columns. Isoforms are identified by their ID in the input long-read GTF file (`*collapsed.gff`). 

<figure style="text-align: center;">
  <img src="https://raw.githubusercontent.com/njdjyxz/BENG_183_Final/main/Iso_Seq_2/classification_file.png" width="600"/>
  <figcaption>
    <strong>Fig. 8: Sample of clssification output file</strong>
  </figcaption>
</figure>

**Columns 1**:`isoform`: the isoform ID. Usually in `PB.X.Y` format. `PB.X`refers to the gene ID, and `.Y`means the number of transcript belongs to this gene.

Columns 2:`chrom`: chromosome.

Columns 3:``strand`: strand.

Columns 4:``length`: isoform length.

**Columns 6**:``structural_category`: one of the categories ["full-splice_match", "incomplete-splice_match", "novel_in_catalog", "novel_not_in_catalog", "genic", "antisense", "fusion", "intergenic", "genic_intron"]

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

#### 1.3 IsoQuant

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
  <img src="https://raw.githubusercontent.com/njdjyxz/BENG_183_Final/main/Iso_Seq_2/result.png" width="600"/>
  <figcaption>
    <strong>Fig. 15: Results of isoforms classification and quantification</strong>
  </figcaption>
</figure>

Based on these results, we can draw meaningful conclusions, such as sample 1 having more NNC than the other samples. With this, we have successfully completed a comprehensive workflow for isoform classification and quantification!

### 2. Finding special genes

Many genes produce multiple isoforms through alternative splicing or other regulatory mechanisms. Examining these isoforms helps identify functional differences, as different isoforms may play distinct roles in cellular processes. One of the most straightforward ways to explore isoform differences between samples is by using IGV (Integrative Genomics Viewer). Below is an example involving the gene *SRSF11*, which exhibits unique isoforms in sample 1 that are absent in the other three samples. By identifying specific genes of interest, we can conduct further analyses to uncover underlying mechanisms and gain deeper insights into their biological significance.

<figure style="text-align: center;">
  <img src="https://raw.githubusercontent.com/njdjyxz/BENG_183_Final/main/Iso_Seq_2/igv.png" width="600"/>
  <figcaption>
    <strong>Fig. 16: Example of gene transcript that differs between samples through igv</strong>
  </figcaption>
</figure>

### 3. Alternative Polyadenylation (APA) Analysis

Isoforms generated through APA events can influence transcript stability, localization, and translation. Classifying these isoforms helps in understanding post-transcriptional regulation.