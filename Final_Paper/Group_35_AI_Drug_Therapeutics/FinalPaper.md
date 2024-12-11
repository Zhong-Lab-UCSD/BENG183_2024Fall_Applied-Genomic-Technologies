# Revolutionizing Drug Discovery: Harnessing AlphaFold and AI to Develop Liver Cancer Therapeutics
#### Team 35: Kimberly, Camilla, and Jenny

## Table of Contents
* [Introduction](#introduction)
* [Key Tools and Technologies](#key-tools-and-technologies)
* [Workflow](#workflow)
* [Results and Summary](#results-and-summary)
* [Input and Output Data](#input-and-output-data)
* [Potential Applications](#potential-applications)
* [Citations](#citations)

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
Using the initial results, a second round of compound generation using Chemistry42 was conducted which generated 16 molecules with the hope of improving binding affinity. Of these 16 molecules, 6 were selected for synthesis and testing to produce more refined candidates. Of these candidates, the molecule ISM042-2-048 
was discovered with nanomolar efficacy, representing a significant improvement in potency. Now you should be able to explain the basics of the general workflow of how AI tools, such as
Alphafold, can be utilized for drug discovery. 

Below you will find a table that briefly summarizes the workflow steps and the AI tools used in each step. Test your knowledge by removing pieces of the table and try to fill it out
yourself!

| Step | Tool/program | Describe basic usage and function of tool in this workflow                |
|-----------|:----------------:|----:              |
| Target Identification   | PandaOmics           | Analyze data and return potential targets         |
| Compound Generation and Testing  | Alphafold, Chemistry42        | Predict 3D structure of CDK2, generate potential compounds |
| Optimization    | All              | Produce more refined candidates |

### Results and Summary 

  The discovery of ISM042-2-001 as an initial hit molecule and its optimization to ISM042-2-048 within an unprecedented 30-day timeframe marks a groundbreaking achievement in drug discovery. This feat is a demonstration of the transformative potential of AI-driven methodologies in accelerating and enhancing the drug development process.

  The rapid timeline for discovering ISM042-2-048 illustrates the power of AI in streamlining complex processes. The entire process, from data analysis to compound synthesis and testing, was completed in just 30 days. The integration of AI platforms significantly compressed the timeline, demonstrating a new standard for efficiency in pharmaceutical research. In comparison, traditional methods of drug discovery require years of development. Target Identification and Validation takes approximately 1-2 years, a lengthy and expensive process in which scientists identify a biological target (ex. protein, receptor) involved in a disease and validate its relevance. Then, hit identification and lead discovery takes another 2-3 years, as large compound libraries must be screened to find “hit” molecules that interact with the target. Hits are then refined into “lead” compounds through iterative testing and optimization to improve potency, selectivity, and pharmacokinetics (Hughes et al.)

![image](https://github.com/user-attachments/assets/9f79beb4-4308-43c7-94de-1df1a8b91b02)

Figure by Hughes et al., Br J Pharmacol. 2011 Mar;162(6):1239–1249.


  In this specific case, for the study of HCC related drug development, of the 8,918 compounds generated, only seven were synthesized and tested in the initial round. This remarkable efficiency reflects Chemistry42's predictive capabilities, which prioritized candidates with the highest potential for success. The AI’s ability to minimize unnecessary synthesis and testing resulted in significant time and resource savings. AlphaFold’s predictive capabilities played a critical role by generating precise 3D models of CDK20, a key target in hepatocellular carcinoma (HCC). By eliminating the need for labor-intensive and costly experimental structural determination, researchers could focus on designing molecules with high specificity and binding efficacy. The results speak volumes about the efficacy of the AI-driven approach. ISM042-2-001, the initial compound, demonstrated a promising binding affinity (K_d = 9.2 μM), validating the AI’s ability to identify viable therapeutic candidates. Subsequent optimization rounds led to ISM042-2-048, a molecule with nanomolar potency. This dramatic improvement highlights the capability of AI to deliver optimized drug candidates with clinical potential in record time (Ren et al.).

![image](https://github.com/user-attachments/assets/c360c87f-f752-4e41-84dd-e6306bc3f189)

Figure by Ren et al., Chem. Sci., 2023,14, 1443-1452

(A) Representative Binding Affinity Curve for ISM042-2-048 in CDK20 Kinase Binding Assay

A binding affinity curve demonstrates the strength of interaction between a ligand (ISM042-2-048) and a target protein (CDK20).
The curve is generated based on experimental data where increasing concentrations of the ligand are tested to determine how well it binds to the protein. The binding affinity is quantified by the dissociation constant (K_d), which reflects the concentration of ligand required for half of the protein's binding sites to be occupied. Data points in the curve represent the mean binding measurements from duplicate wells (replicates) in one experiment, which ensures consistency and reliability. This curve helps assess how effectively ISM042-2-048 binds to CDK20, which is crucial in evaluating its potential as a drug candidate.

(B) Predicted Binding Pose for ISM042-2-048 in CDK20

The predicted binding pose is a computational model that shows how the molecule (ISM042-2-048) is expected to interact with the binding pocket of the CDK20 protein at the molecular level. It illustrates the specific interactions between ISM042-2-048 and the CDK20 protein, such as hydrogen bonds, hydrophobic interactions, and ionic interactions, which are critical in stabilizing the binding and influencing the drug's effectiveness. This provides insights into the molecular mechanism of action for ISM042-2-048 and helps identify key residues in CDK20 that are essential for binding, which can guide further optimization of the compound for improved efficacy.

### Input and Output Data

#### Input Data

1. **Omics Data**  
   The platform integrates large-scale, high-dimensional datasets such as gene expression profiles and proteomics. These data sources provide an in-depth view of biological processes, enabling the identification of critical pathways and potential therapeutic targets. This comprehensive analysis forms the foundation for precise and informed intervention strategies. 

3. **Textual Data**  
   The system also incorporates textual data from scientific literature and curated databases. This information is crucial for understanding disease associations, predicting the relevance of potential targets, and prioritizing compounds for development. By synthesizing these insights, the platform establishes meaningful connections between molecular mechanisms and clinical outcomes, guiding its discovery processes with greater accuracy.

4. **AlphaFold-Predicted Structures**  
   Additionally, the platform utilizes AlphaFold-predicted structures to enhance its capabilities. For instance, the computationally modeled 3D structure of CDK20 provides valuable insights into the structural basis of protein function. Even in the absence of experimentally determined structures, these predictions enable the exploration of protein interactions and functional mechanisms. 

##### Output Data

1. **Therapeutic Targets**  
   The workflow provides a prioritized list of therapeutic targets, tailored to specific disease contexts. For hepatocellular carcinoma (HCC), 20 potential candidates were ranked based on their relevance to disease mechanisms and druggability. This ensures a strategic emphasis on high-impact, actionable targets that align with therapeutic priorities.

2. **Compound Library**  
   An extensive virtual compound library is generated through data-driven screening processes. This facilitates the identification of molecules with predicted activity against selected targets, enabling a refined approach to compound selection. For HCC, 8,918 compounds were initially screened, and seven were selected for synthesis based on their predicted binding affinities and pharmacological profiles. This targeted selection minimizes resource expenditure while maintaining precision and maximizing potential therapeutic effectiveness.

3. **Validation Outcomes**  
   The workflow delivers two validated compounds as key outcomes. The first is an initial lead molecule exhibiting micromolar binding affinity, demonstrating therapeutic feasibility. The second is an optimized molecule with nanomolar potency, highlighting the system’s ability to refine and enhance drug candidates for greater efficacy and precision.

### Potential Applications

1. **Broader Disease Applications**  
   The methodology can be extended to other cancers and diseases, enabling the identification of novel therapeutic targets and tailored treatments across various conditions.

2. **Personalized Medicine**  
   By utilizing patient-specific omics data, the platform can design treatments tailored to individual genetic and molecular profiles, enhancing efficacy and safety.

3. **Expanding Druggable Targets**  
   AlphaFold enables the exploration of previously undruggable proteins by predicting their 3D structures, unlocking new possibilities for therapeutic intervention.


## Citations
* Ren, F., Ding, X., Zheng, M., Korzinkin, M., Cai, X., Zhu, W., Mantsyzov, A., Aliper, A., Aladinskiy, V., Cao, Z., Kong, S., Long, X., Liu, B. H. M., Liu, Y., Naumov, V., Shneyderman, A., Ozerov, I. V., Wang, J., Pun, F. W., … Zhavoronkov, A. (2023, January 10). Alphafold accelerates artificial intelligence powered drug discovery: Efficient discovery of a novel cdk20 small molecule inhibitor. Chemical Science. https://pubs.rsc.org/en/Content/ArticleLanding/2023/SC/D2SC05709C 
* Warner, E. (2024, February 12). New Study uses alphafold and AI to accelerate design of novel drug for liver cancer. UToronto - Faculty of Arts & Science. https://www.artsci.utoronto.ca/news/new-study-uses-alphafold-and-ai-accelerate-design-novel-drug-liver-cancer
* Hughes JP, Rees S, Kalindjian SB, Philpott KL. Principles of early drug discovery. Br J Pharmacol. 2011 Mar;162(6):1239-49. doi: 10.1111/j.1476-5381.2010.01127.x. PMID: 21091654; PMCID: PMC3058157.
* Callaway, Ewen. "Major AlphaFold Upgrade Offers Boost for Drug Discovery: Latest Version of the AI Models How Proteins Interact with Other Molecules — but DeepMind Restricts Access to the Tool." Nature, 8 May 2024, https://doi.org/10.1038/d41586-024-01383-z.

