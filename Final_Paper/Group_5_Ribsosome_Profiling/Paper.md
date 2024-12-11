# Ribosome Profiling

## Procedures
### DNA library Preparation     
Similar to RNA-seq, the preparation of ribosome profiling libraries requires cell lysis, mRNA purification, and reverse transcription. Here's an overview of the steps:
1) We first need to freeze the translation in cells before lysing them by treating them with elongation inhibitors, such as cycloheximide. To interrogate the initiation site of translation, the cells should be treated with harringtonine, which would cause the ribosomes to accumulate at initiation codons. 
2) After appropriate treatments, the cell should now be ready to be lysed by lysis reagents. For animal-derived samples, such as tissues that need to be physically disrupted, it is a robust approach to cryogenically pulverize the tissues and then thaw the tissues in the translation inhibitors. Note that these first two steps may vary depending on the types of samples and species.
3) Nuclease would then be used to digest regions on mRNA that are not protected by ribosomes. 
4) Next, an ultracentrifuge needs to be performed to spin down and separate the mRNA-ribosome complex with the naked mRNA. The ribosomes would then be washed off the mRNA.
The mRNA sequences are now obtained. The next steps are to sequence them.
5) Linker sequences would first be ligated to the 3’ end of the mRNA fragments.
6) Reverse transcriptase and primers would then be added to perform a reverse transcription.
7) Among the cDNA reverse transcription products, many are rRNA contamination that needed to be eliminated by hybridizing them with biotinylated sense strand oligonucleotides and then treating with streptavidin which would bind to the biotin. Note that contamination may be still present after this step, further elimination of this contamination would need to be done in the later bioinformatics data analysis.
8) The rest of the cDNA would next be PCR amplified to create a ribosome footprint library.

9) ### Data Analysis  
The bioinformatics data analysis of ribosome profiling has many overlaps with that of RNA-seq, and many pipelines and packages have been established for different analysis purposes. 
The graph below shows an overview of a STAR Protocol of ribosome profiling data analysis pipeline specialized in ribosome pausing analysis ([link to this paper](https://star-protocols.cell.com/protocols/1899)):
1) Preparation: raw read fastq files and reference genome of the species need to be downloaded (previously published datasets can be downloaded from SRA). 
2) Preprocessing: similar to RNA-seq, before the alignment we also need to trim off the adapter sequence from the reads and do quality control. One extra step is that we need to discard reads that are too short because they might be ncRNA contamination. These can be done in one step using fastx_toolkit using fastx_clipper command with various options:
   ```
   fastx_clipper -Q<phred_score> -a <adapter_sequence> -l <minimum_read_length_to_keep> -i <input.fastq> -o <output.fastq>
   ```
4) Read alignment: to further eliminate contamination from ncRNA, the raw reads need to be first aligned to the ncRNA reference genome, and discard the aligned reads 
The creation of ncRNA indices and the alignment to the indices can be done using bowtie2:
   ```
   bowtie2-build <reference_ncRNA_sequences.fa> <prefix_of_output_files>
   bowtie2 -L <seed_substring_length> --un=<path_where_unmatched_reads_would_be_saved> -x <reference_ncRNA_indices> > <output.sam>
   ```
   We next use the rest of the reads to align to the reference genome. We can use STAR for the creation of reference genome indices and the alignment with a GTF annotation file:
   ```
   STAR --runMode genomeGenerate --genomeDir <reference_genome_directory> --genomeFastaFiles <reference_genome.fa> --sjdbGTFfile     <annotation_file.gtf>
   STAR --genomeDir <path_to_reference_genome_indices> --readFilesIn <input_reads.fastq> --outFileNamePrefix           <prefix_of_output_file> --outSAMtype BAM SortedByCoordinate
   ```
4) The steps after this vary depending on the purpose of the analysis, for example: <br>To detect the change in translation efficiency, both the reads from ribosome profiling and mRNA seq would be mapped to CDS regions and exons using featureCounts, and by using a package called RiboDiff, the ratio between the counts of  ribosome profiling reads and RNA-seq reads of each CDS or exon would be calculated and compared between different treatments ([manual of RiboDiff](https://github.com/ratschlab/RiboDiff)). 
<br>To determine ribosome distribution, a count array (lists that record the number of reads that are mapped to each codon of a transcript), would be generated with the bam file resulting from the genome alignment along with GTF genome annotation file using a package called Plastid. By using LOESS smoothing with the count array, we can then calculate cumulative ribosome distribution on transcripts ([manual of Plastid](https://plastid.readthedocs.io/en/latest/generated/plastid.html); [manual of LOESS](https://pypi.org/project/loess/#documentation))

| Package Name  |                 Application               |  Organism  |
|:-------------:|:-----------------------------------------:|:----------:|
| HRIBO      | Bacterial Ribosome Profiling data analysis | Prokaryotes |
| ORFik      | Translation complex profiling, ribosome complex profiling, gene expression analysis      |   Prokaryotes and eukaryotes |
| RiboA |  Calculation of accurate A-site offset values   |   Prokaryotes and eukaryotes |


Ever since the introduction of ribosome profiling by Ingolia et al. in 2009, the application of ribosome profiling data has been continually growing, to better integrate the data, RiboSeq.Org was recently introduced. This web browser portal curates tens of thousands of datasets from many studies and provides various data analysis and visualization tools ([link to the portal](https://rdp.ucc.ie/home)). 

