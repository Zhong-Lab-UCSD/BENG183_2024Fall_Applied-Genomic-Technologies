# Group_16_DeepAEG

## Introduction

### Background
Cancer is a disease that affects millions of people worldwide, and there is no single cure for cancer due to the heterogeneity of patients, the diversity of the disease itself, and the uncertainty of
drug efficacy. Because of this, it is incredibly difficult to treat resulting in losses worldwide. A hopeful solution to this issue is the idea of individualized cancer treatments, in which we can utilize
a patient's genotype and their unique drug reactions in order to improve therapeutic effect while reducing drug usage side effects. In order to efficiently narrow down treatments for patients and predict
drug responses, we can utilize machine learning. 
### Limitations of Past Models
Several models have been proposed, such as simple logistic regression models, SVMs, and more as well as more complicated models like tCNNs(maybe embed
links) and DeepTTA. However, these models face limitations that cause the model to lose vital information during feature learning. One example of this is the loss of important molecular structural information such
as charge and stereoisomerism due to the nature of the molecular fingerprinting/SMILES that is used to linearly encode drug structure. 

Due to this, graph based approaches gained popularity because the natural structure of drugs can easily be represented using a graph as well as the power of GCNs that can integrate multiple omics features.
DeepCDR is a popular graph-based model that achieves excellent results compared to traditional models, but its performance is still impeded due to the loss of molecular bond information due to the 
intricate nature of molecular edge characterization and the constraints that come with updating edges in graph-based systems. Molecular bonds are as vital as molecular structure, so the loss of this information
results in the results of the model taking a huge hit. Because of this, it is absolutely necessary to optimize and improve the edge updating fusion algorithm.

### The Advantages of DeepAEG 
To overcome the limitations described above, the Zhejiang lab proposed a novel multi-source heterogeneous graph convolutional neural network, also known as DeepAEG. 