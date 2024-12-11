# AlphaFold Protein Structure Database

### By Serena Chuang, Jessica Wang, Annapurna Saladi (Group 36)

* [Protein Structure](#protein-structure)
* [What is AlphaFold?](#what-is-alphafold)
* [How to use?](#how-to-use)
* [Achievements and Limitations](#achievements-and-limitations)
* [Sources](#sources)

---

## Protein Structure

Proteins are commonly known as the fundamental building blocks for all living organisms. They carry out important tasks in cell, from regulation, structure, and functions of the organ's they belong to. When proteins fail to do their function or hinders the body's function, we see these larger organs fail, in the event of diseases, like in Alzheimer's[1] or ALS[2]. In order to understand how these proteins function, and what happens when they do not function, we can look to studying their specialized structure. Proteins are constructed from amino acid sequences that fold in a certain way, resulting in a mechanism that is specially designed to handling a certain task. For example, a hemoglobin protein that is used for oxygen transport will have specialized binding sites for oxygen molecules, and otpimize the binding and transport of those molecules so they can be delivered across the body[3]. From knowing what shape a protein has, we can understand how it interactes with other molecules in the body, and what role they play in pathways.

Over the years, we have tried to look at protein structure using a variety of methods, each having their own positives and negatives. Three of the main methods we have today are NMR spectroscopy, x-ray crystallography, and cryogenic electron microscopy. One thing these methods have in common is that they are experimental, meaning they rely on collecting real data of the proteins, and 

Because of these 

Proteins fold in very complex ways that rely more factors than just the peptides that make up its sequence:

structure affected by number of internal and external forces, including hydrogen bonding, hydrophobic interaction, dimers, electrostatic interactions, van der Waals forces, environment, etc.

Protein can fold in many different ways along an energy landscape, when created they follow the lowest energy landscape, but this can be hard to find

Proteins can take on a number of different conformations depending on their role, predicting this dynamic behavior can be difficult 



---

## What is AlphaFold?
AlphaFold is an AI system developed by Google DeepMind to predict a protein’s 3D structure using its amino acid sequence. Currently, the AlphaFold protein structure database contains over 200 million predicted protein structures, including the complete human proteome. It recently won the 2024 Nobel Prize in Chemistry. 

### Importance of AlphaFold
AlphaFold has solved the “protein folding problem”, can predict protein structures with high accuracy, and accelerates research in biotechnology and drug development fields as it contributes to targeted drug discovery. 

### Algorithms
There are 3 distinct AlphaFold algorithms: 1, 2, and 3. 

* Alphafold 1 was released in 2018 and introduced a novel approach to protein folding using machine learning. 

* Alphafold 2 incorporated an attention-based neural network, Evoformer, and calculated a “pair representation” for each residue pair to increase its overall accuracy [4].  

* Alphafold 3, released earlier this year in May, is able to model the structures and interactions of proteins with DNA, RNA, ligands, and ions. By modeling protein-ligand interactions, AlphaFold 3 is able to model complexes rather than single chains to contribute to drug development [5]. 


---

## How to use?



### Usage

There are two ways to use AlphaFold. First, there is the [AlphaFold Protein Structure Database](https://alphafold.ebi.ac.uk/), which provides free accessible protein structure predictions. This database contains over 200 million entries, which covers a substantial amount of Uniprot, the database that contains protein sequences and annotations. 

![AlphaFold Database Image](https://github.com/serrachow/beng183images/blob/main/alphafold_db.png)

For example, let's use `PTBP1` as an example gene. We can type in `PTBP1` into the search bar, and select the first option, which shows `Polypyrimidine tract-binding protein 1`. 

![AlphaFold Database Example](https://github.com/serrachow/beng183images/blob/main/alphafold_db_example.png)

A list of what the AlphaFold Protein Structure Database offers:

* `PDB file` - The Protein Data Bank file contains atomic-level coordinates of the protein structure, including atoms, residues, chains, and structural information. This can be downloaded and viewed using tools, such as Chimera.
* `mmCIF file` - The maromolecular Crystallographic Information File is an alternative to the PDB file and is better suited for larger structures. This file also contains metadata about the protein structure and can also be viewed using tools, such as Chimera. 
* `Predicted aligned error map` - The PAE map shows the confidence of the AlphaFold prediction at specific areas in the structure. Values closer to zero indicate a higher confidence, and this map can be visualized using ChimeraX.
* `Structure viewer` - The 3D protein structure viewer features the 3D model, with users being able to rotate, zoom, and pan the structure. Underneath, there are colors showing how confident the structure is, and there is an interactive panel to select different residues and chains.
* `AlphaMissense Pathogenicity Heatmap` - This heatmap shows how likely a missense mutation at a specific residue position is pathogenic or benign, with red being the most likely pathogenic to green, which is likely benign.

Alternatively, if you wanted to run your own simulations and predictions, there is the [AlphaFold v2 pipeline Github Repository](https://github.com/google-deepmind/alphafold), which is an open source package, that you can download and run on a Linux machine. 

![AlphaFold Github](https://github.com/serrachow/beng183images/blob/main/alphafold_github.png)

The exact instructions are listed in the `README.md` document in the Github repository, but we will go over the general idea of how AlphaFold predictions work. 

AlphaFold takes in two types of input data: first a protein sequence of amino acids and second, multiple sequence alignments (MSAs), which are homologous sequences that provide extra evolutionary information. These should be both in FASTA format. MSA help identify co-evolving residues in proteins and provide more information for spatial relationships. 

The main inputs AlphaFold take in are:

* `--fasta_paths` - This takes in the full path and filename to your test data.
* `--output_dir` - This takes in the full path the desired output directory.
* `--model_preset` - This controls which AlphaFold model to run. There are two main types of models, monomer, which is used for simpler protein sequences, and multimer, which is used for a multi-sequence FASTA file to generate a multi-chained protein. 
* `--max_template_date` - The max template date can ensure predictions stay realistic by setting cut off dates for evolutionary information.

The AlphaFold outputs include the computed MSAs, unrelaxed structures (PDB file), relaxed structures (PDB file) that include an extra relaxation step to improve goemetry, ranked structures (PDB file) that contains structures ordered by model confidence, raw model outputs, prediction metadata, and section timings. 

### Framework

AlphaFold is built on a complex neural network focused on end-to-end prediction with a core component being the Evoformer. This is the central block of AlphaFold and is responsible for processing the MSAs and their pair features, which capture residue-residue relationships including the distance and orientation of the protein. In addition, another core components is the structure module that refines the 3D atomic coordinates and ensures that physical realism is kept, such as using geometric constraints like bond lengths and angles. Lastly, structures are repeatedly fed back into the network to refine the structure, resolve inconsistencies, and improve overall accuracy.

The attention mechanism in neural networks allow models to focus on specific relevant parts of input data when making predictions. This helps models efficiently capture relationships and dependencies across multiple elements and improve overall accuracy. AlphaFold uses attention in multiple ways: self-attention in MSAs capture co-evolutionary signals among amino acids, self-attention in 3D structures help predict spatial proximity and structural interactions, and cross-attention between representations integrates MSAs and intermediate structures together.

AlphaFold was first trained on supervised learning on PDB data to generate a new dataset of predicted structures. Then it was later trained on the same architecture from scratch using a mixture of PDB data and the new dataset of predicted structures. This self-distillation procedure improves the accuracy of the resulting network.

---

## Achievements and Limitations



---

## Sources

[1] Penke B, Bogár F, Paragi G, Gera J, Fülöp L. Key Peptides and Proteins in Alzheimer's Disease. Curr Protein Pept Sci. 2019;20(6):577-599. doi: 10.2174/1389203720666190103123434. PMID: 30605056.

[2] Maguire G. Amyotrophic lateral sclerosis as a protein level, non-genomic disease: Therapy with S2RM exosome released molecules. World J Stem Cells. 2017 Nov 26;9(11):187-202. doi: 10.4252/wjsc.v9.i11.187. PMID: 29312526; PMCID: PMC5745587.

[3] Pittman RN. Regulation of Tissue Oxygenation. San Rafael (CA): Morgan & Claypool Life Sciences; 2011. Chapter 4, Oxygen Transport. Available from: https://www.ncbi.nlm.nih.gov/books/NBK54103/

[4] Jumper, J., Evans, R., Pritzel, A. et al. Highly accurate protein structure prediction with AlphaFold. Nature 596, 583–589 (2021). https://doi.org/10.1038/s41586-021-03819-2

[5] Abramson, J., Adler, J., Dunger, J. et al. Accurate structure prediction of biomolecular interactions with AlphaFold 3. Nature 630, 493–500 (2024). https://doi.org/10.1038/s41586-024-07487-w
---
