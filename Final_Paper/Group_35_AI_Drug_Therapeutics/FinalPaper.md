# Revolutionizing Drug Discovery: Harnessing AlphaFold and AI to Develop Liver Cancer Therapeutics
#### Team 35: Kim, Camilla, and Jenny

## Table of Contents
* [Introduction](#introduction)
* [Key Tools and Technologies](#key-tools-and-technologies)
* [Workflow](#workflow)
* [Results and Summary](#results-and-summary)
* [Input and Output Data](#input-and-output-data)
* [Potential Applications](#potential-applications)
* [References](#references)

## Introduction
Traditional drug discovery relies on slow, expensive trial-and-error chemistry, limiting the scope of exploration and innovation. As the COVID-19 pandemic has shown, the public increasingly expects the rapid development of effective treatments. Artificial intelligence has the potential to transform the speed and efficiency of drug discovery, as it has done for numerous scientific and engineering disciplines over the past decade. AlphaFold, PandaOmics, and Chemistry42 have helped identify a novel hit molecule against a novel target without an experimental structure, starting from target selection towards hit identification, in a cost-and-time efficient manner. 

### Hepatocellular Carcinoma
The focus of the applications of these AI tools will be on Hepatocellular Carcinoma (HCC) due to its high prevalence in liver cancers and lack of effective treatments.HCC is the most common type of primary liver cancer accounting for about 75-85% of all liver cancers. This specific case study is the first reported example that successfully utilized AlphaFold-predicted protein structures to identify a confirmed hit for a novel target in early drug discovery. 

## Key Tools and Technologies

### AlphaFold

AlphaFold is an AI system developed by Google DeepMind that predicts a protein's 3D structure from its amino acid sequence with high accuracy. It provides critical insights into protein functions and binding sites, bypassing the need for experimental structural data and is a pivotal tool in accelerating the drug discovery process by enabling structure based design. AlphaFold has revolutionized structural biology, by reducing time and resources needed to determine protein structures, which traditionally required labor intensive techniques like X-ray crystallography or NMR spectroscopy. 

> AlphaFold Protein Structure Database: https://alphafold.ebi.ac.uk/

### PandaOmics

PandaOmics is an AI-driven platform developed by Insilico Medicine for therapeutic target and biomarker discovery. By integrating diverse datasets such as gene expression, proteomics, and transcriptomics with scientific literature, it identifies and ranks potential drug targets. It then filters these targets based on novelty, druggability, safety, and disease relevance, streamlining research focus. 

> PandaOmics: https://pharma.ai/pandaomics

### Chemistry42

Chemistry42 is a generative chemistry platform developed by Insilico Medicine that designs and optimizes novel drug-like molecules. It uses machine learning to evaluate compound interactions with target proteins, prioritizing the best candidates for testing and synthesis. 
Some key features of Chemistry 42 include: 
- Generative Chemistry - designs new small molecules optimized for specific targets, and facilitating de novo design
- ADMET Profiling - predicts and refines absorption, distribution, metabolism, excretion, and toxicity profiles of molecules
- Alchemistry - estimates relative binding free energies
- Golden Cubes - predicts kinome activity

> Chemistry42: https://pharma.ai/chemistry42


## Workflow
We will now dive into the general workflow for this project using hepatocellular carcinoma (HCC). 
This will demonstrate the usefulness of various AI tools for commercial purposes such as drug discovery.

### Target Identification
PandaOmics was used to analyze OMICs and textual data that associate genes with the disease of interest from 10 different datasets related to HCC and returned a ranked list of the top 20 
targets which could be sorted by novelty, safety, tissue specificity, etc.
The text data selects genes that are frequently mentioned in various scientific text sources. On the other hand, OMICs scores explore the molecular connection of the genes with disease 
based on gene variants, knockout/overexpression studies, differential expression, etc.
This tool allowed the researchers to narrow down the scope of the study's therapeutics targets using multidimensional filtering criteria. Of these 20 targets, cyclin-dependent 
kinase 20 (CDK20) was selected as the initial target due to its strong association with HCC as it is overexpressed in many tumor cell lines, limited information about its structure and lack of approved drugs or 
clinical compounds targeting it. Now that an initial target has been determined, it is time to move on to generation of the compound and testing.

### Compound Generation and Testing
Alphafold was used to predict the 3D structure of CDK2. The surface of the protein is covered in methyl probes whose bonding interaction energy is calculated and then clustered and scored to produce a list of binding sites using Chemistry42. Using this approach, Chemistry42 was then able to use the Alphafold-generated 3D model to produce 8,918 potential compounds. Out of these compounds, after molecular docking and clustering, 7 were selected for 
synthesis and biological testing. Using kinase binding assays, the compound ISM042-2-001 showed a binding affinity (K_d) of 9.2 μM. These results can then be used for the final step.

### Optimization
Using the initial results, a second round of compound generation using Chemistry42 was conducted which generated 16 molecules with the hope of improving binding affinity. Of these 16 molecules, 6 were selected for synthesis and testing was conducted to produce more refined candidates. Of these candidates, the molecule ISM042-2-048 
was discovered with nanomolar efficacy, representing a significant improvement in potency. Now you should be able to explain the basics of the general workflow of how AI tools, such as
Alphafold, can be utilized for drug discovery. 

Below you will find a table that briefly summarizes the workflow steps and the AI tools used in each step. Test your knowledge by removing pieces of the table and try to fill it out
yourself!

| Step | Tool/program | Describe basic usage and function of tool in this workflow                |
|-----------|:----------------:|----:              |
| Target Identification   | PandaOmics           | Analyze data and return potential targets         |
| Compound Generation and Testing  | Alphafold, Chemistry42        | Predict 3D structure of CDK2, generate potential compounds |
| Optimization    | All              | Produce more refined candidates |

## Citations
* Ren, F., Ding, X., Zheng, M., Korzinkin, M., Cai, X., Zhu, W., Mantsyzov, A., Aliper, A., Aladinskiy, V., Cao, Z., Kong, S., Long, X., Liu, B. H. M., Liu, Y., Naumov, V., Shneyderman, A., Ozerov, I. V., Wang, J., Pun, F. W., … Zhavoronkov, A. (2023, January 10). Alphafold accelerates artificial intelligence powered drug discovery: Efficient discovery of a novel cdk20 small molecule inhibitor. Chemical Science. https://pubs.rsc.org/en/Content/ArticleLanding/2023/SC/D2SC05709C 
* Warner, E. (2024, February 12). New Study uses alphafold and AI to accelerate design of novel drug for liver cancer. UToronto - Faculty of Arts & Science. https://www.artsci.utoronto.ca/news/new-study-uses-alphafold-and-ai-accelerate-design-novel-drug-liver-cancer 

