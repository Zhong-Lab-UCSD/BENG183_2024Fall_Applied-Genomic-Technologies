# AlphaFold Protein Structure Database

### By Serena Chuang, Jessica Wang, Annapurna Saladi (Group 36)

* [Protein Structure](#protein-structure)
* [What is AlphaFold?](#what-is-alphafold)
* [How to use?](#how-to-use)
* [Applications](#applications)
* [Ahicevements and Limitations](#achievements-and-limitations)
* [Sources](#sources)

---

## Protein Structure

---

## What is AlphaFold?

---

## How to use?

# Usage

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

# Framework

AlphaFold is built on a complex neural network focused on end-to-end prediction with a core component being the Evoformer. This is the central block of AlphaFold and is responsible for processing the MSAs and their pair features, which capture residue-residue relationships including the distance and orientation of the protein. In addition, another core components is the structure module that refines the 3D atomic coordinates and ensures that physical realism is kept, such as using geometric constraints like bond lengths and angles. Lastly, structures are repeatedly fed back into the network to refine the structure, resolve inconsistencies, and improve overall accuracy.

The attention mechanism in neural networks allow models to focus on specific relevant parts of input data when making predictions. This helps models efficiently capture relationships and dependencies across multiple elements and improve overall accuracy. AlphaFold uses attention in multiple ways: self-attention in MSAs capture co-evolutionary signals among amino acids, self-attention in 3D structures help predict spatial proximity and structural interactions, and cross-attention between representations integrates MSAs and intermediate structures together.

AlphaFold was first trained on supervised learning on PDB data to generate a new dataset of predicted structures. Then it was later trained on the same architecture from scratch using a mixture of PDB data and the new dataset of predicted structures. This self-distillation procedure improves the accuracy of the resulting network.

---

## Applications

---

## Achievements and Limitations

---

## Sources

---
