# Expression Quantitative Trait Loci (eQTL): Bridging Genotype and Phenotype

## 1. Introduction to eQTLs

With the emergence of Genome Wide Association Studies (GWAS) in recent years, scientists have identified thousands of single nucleotide polymorphisms (SNPs) associated with complex diseases [1]. Despite this method being able to pinpoint potential risk loci that could explain these traits, the underlying mechanisms remain unavailable due to the majority of effects residing in non-coding regions of the genome [2], as well as the presence of pleiotropy and linkage [3]. A way to target the problem of discovery of causal variants behind these traits lies in leveraging Expression Quantitative Trait Loci (eQTLs). Briefly, eQTLs are genomic loci where variants influence the expression levels of one or more genes, providing a direct connection between genetic variation and transcriptional regulation. By directly linking genetic variants to changes in gene expression levels, eQTLs represent a powerful approach for understanding how genetic variation influences gene expression. In the era of genomics, eQTLs serve as crucial tools for bridging the gap between genetic variants and their functional impacts on cellular processes. 

### 1.1 What are eQTLs?

eQTLs are genomic regions that contain sequence variants (typically SNPs) that influence the expression level of one or more genes. These regulatory relationships help us understand:
- How genetic variation affects gene expression
- The molecular mechanisms underlying disease associations
- The functional significance of variants identified in genome-wide association studies (GWAS)

### 1.2 Linkage Disequilibrium

#### 1.2.1 Linkage & Recombination

Recombination is a process in which DNA is broken and recombined to form a haplotype structure that could not exist otherwise. In eukaryotic cells, recombination usually takes place during meiosis [4]. 

Linkage is the idea that sites that are closer together in the genome typically are inherited together [5]. Thus, analyzing a set of haplotypes unveils a combination of alleles that tend to appear together.

The important thing to understand is that linkage propagates linkage disequilibrium (LD), that is, the nonrandom combination of alleles. On the other hand, recombination breaks down LD.

#### 1.2.2 Linkage Disequilibrium (LD)

How do we measure LD? 

LD can be thought of as the correlation between various SNPs in a sample. Thus, we can convert our alleles to a numerical encoding, where 0 indicates there are no copies of the alternate allele, 1 indicates one copy of the alternate alleles, and 2 indicates two copies of the alternate alleles. Thus, across many samples we can create a SNP by SNP matrix and measure LD by computing the Pearson correlation coefficient [6]. A coefficient of 1 indicates that SNPs are in perfect LD, that is, they are always inherited together. A coefficient of 0 indicates that SNPs are completely independent of each other. 

![unnamed](https://github.com/user-attachments/assets/4c5e7154-d4ee-4a6d-89e0-9824852143ba)

Figure 1. LD matrix demonstration [7]

The above image demonstrates how an LD matrix might look like. You will notice that the diagonal has a correlation coefficient of 1. This makes sense because a SNP must be in perfect correlation with itself. As we move farther from the diagonal, the LD breaks down because the SNPs are farther away from one another and thus less likely to be inherited together.

An important property of LD is that it is typically population dependent. Different populations will have different linkage patterns and thus different LD structure. The below image highlights this.

![unnamed (1)](https://github.com/user-attachments/assets/0a053dc6-31b7-4a34-8e05-31953c18ce31)

Figure 2. LD matrix for different cohorts. [8]

YRI is Yoruba, CEU is Europeans, and CHB + JPT is Chinese and Japanese respectively. 

### 2.1 Generating eQTLS

#### 2.1.1 Data

To generate an eQTL, you will need expression data in the form of normalized RNA sequencing data. You will also need genotypes either from whole genome sequencing or genotyping arrays. There exist consortiums that have massive amounts of available data like the UK biobank. 

Once you have acquired the data, you will need to standardize both the expression and genotyping data to mean zero variance 1. You will also need to regress out covariates in your data. This usually comes from associated metadata and the top 10 principal components in both expression and genotyping data.

In R, regressing out the covariates might look something like this:

```r 
train_gene_Y_resid <- resid(lm(train_gene_Y ~ 
training_gt_pca$rotation[,1:10] + 
training_ge_pca$rotation[,1:10]))
 ```

#### 2.1.2 Regression

The next step is to regress the expression levels on genotyping data. The regression is a very standard ordinary least squares. You have to compute a regression between every SNP per gene. So, if you have n genes and m SNPs, you will effectively compute nm regressions, each with their unique effect size i.e. the slope.

In R, this might look something like this:

```r 
eqtl_analysis <- function(i) {
    mod <- lm(train_Y_resid ~ training_full_snp_matrix_std[, i]) 
    summary(mod)$coefficients[2, c(1, 2, 4)] 
}
```

You would apply this function to every gene and it iterates through the SNPs. 


### 2.2 Statistical fine-mapping (Causal Inference)

#### 2.2.1 Fine-mapping: overview
To address the limitations of GWAS in identifying causal variants, particularly due to linkage disequilibrium (LD), fine-mapping methods have been developed. Earlier methods include CAVIAR [9] and FINEMAP [10] that provided substantial advances in pinpointing potentially causal variants. However, these methods face challenges with scalability: as more genes/SNPs appear in the analysis, it becomes much more computationally costly to perform fine-mapping.

#### 2.2.2 Sum of Single Effects model
In recent years, SuSiE (Sum Of Single Effects) [11] has emerged as a powerful alternative. SuSiE addresses these limitation with an Iterative Bayesian Stepwise Selection (IBSS) algorithm, based on the assumption that phenotype effect can be explained by a sum of single effects:

![unnamed (2)](https://github.com/user-attachments/assets/d49db34b-090e-4e59-bc76-8460ced5e029)

Figure 3. Iterative Bayesian Stepwise Selection

Here we show how IBSS works. The central function behind the algorithm is Single Effect Regression (SER), which linearly fits y = Xb + e, where X – genotypes matrix, b – effect vector with 1 SNP, e – is an N-vector of error terms, y – phenotypes of N individuals. The SER function calculates α – posterior inclusion probability (PIP), indicating how likely it is that this single-effect component is nonzero; μ – posterior mean of the effect size; σ – posterior standard deviation of the effect size. Briefly, SER computes PIP by calculating Bayes Factor, which quantifies the strength of evidence in favor of a hypothesis that SNP is causal compared to a null hypothesis (b=0), and then normalizes it across all effects. More detail on how it is achieved is described in the SuSiE’s paper supplementary material [11].
We now describe how the algorithm works line by line.
- Line 1, we assume no effects are present. 
- Line 2, repeat means that we are entering a while loop that will not be finished until the convergence criterion is satisfied (Line 7). SuSiE assumes that IBSS is a Evidence Lower Bound (ELBO) algorithm, meaning that convergence is satisfied if there exists no increase in ELBO (i.e. in variational inference this means our posterior distribution approximation is similar to the predicted). Alternatively, if SuSiE doesn’t converge within 100 or a customly specified number of iterations, the algorithm stops and returns a warning to the user.
- Line 3 starts a for loop within the L SNPs of the gene.
- Line 4, compute the “expected residuals” after removing the contribution of the l-th effect. This means we remove all the effects explained by other SNPs and only keep the remaining L effect.
- Line 5, we fit the SER function that returns us an estimate of how strongly the l-th effect contributes to the model.
- Line 6, integrates the effect size estimate and probability of this effect being causal.
- Line 7, convergence criteria (mentioned in line 2).
Once the algorithm converges, we get a distribution of PIPs, posterior effects, and errors. This enables us to do further analysis on which SNPs are causal. It is important to note that the IBSS algorithm works for individual level data. This means one needs to provide Individual Genotypes and Phenotypes data to perform fine-mapping. To perform fine-mapping on eQTLs, SuSiE developers introduce “sufficient statistics” that are derived from eQTLs. These Sufficient Statistics in combination with the LD matrix of a cohort can be used to perform fine-mapping on eQTLs. The modified algorithm is called IBSS-ss:

![unnamed (3)](https://github.com/user-attachments/assets/16efe23d-c393-42da-a1a4-b2770f58fb02)

Figure 4. IBSS-ss algorithm.

Here XTX – correlation from LD matrix, XTy – SNP-phenotype association (sufficient statistics) that is derived from eQTLs. Conceptually, the algorithm is very similar to IBSS with the only difference being that we compute residuals in ρ that are derived from sufficient statistics rather than individual level data. More on the argument on how sufficient statistics are derived in the SuSiE-RSS [12].

#### 2.2.3 Fine-mapping in practice
![unnamed (4)](https://github.com/user-attachments/assets/65575067-e088-4958-a9fb-d5cd127ba20b)

Figure 5. Example of how SuSiE fine-mapping works.

To understand fine-mapping in practice, we refer to a figure from the SuSiE-RSS paper [12] that illustrates how a refinement step can improve the performance of SuSiE-RSS. Although we do not detail this optional refinement procedure in our algorithm description, it can enhance results when the initial solution is not globally optimal. Regardless, the figure provides an excellent example of how fine-mapping works. In panel 5A, we see a zoomed-in view of a simulated GWAS locus from the UK Biobank dataset. Initially, numerous SNPs are present in the region, but after applying SuSiE-RSS, they are consolidated into a smaller set of putative causal variants (SMA, SNP1, and SNP2). Panels 5B and 5C then show the fine-mapping results in terms of Posterior Inclusion Probabilities (PIPs).

### 2.3 Polygenic Risk Score

#### 2.3.1 Introduction to Polygenic Risk Scores
Polygenic Risk Scores (PRS) represent a powerful approach for translating genetic findings into clinically relevant predictions. While individual genetic variants typically have small effects on complex traits, PRS aggregate these weak signals across many genomic regions to generate meaningful predictions about phenotypes.

#### 2.3.2 Implementation Example
Recent work [13] has demonstrated the effectiveness of multiple PRS approaches. The GWS PRS, using genome-wide significant SNPs, provides a baseline approach utilizing 65 key genetic markers. The eQTL PRS expands this by incorporating 961 expression-related SNPs, offering additional predictive power through the integration of functional genetic information.
Combined approaches have shown stronger results accounting for overlapping SNPs within genetic loci and linkage disequilibrium.

#### 2.3.3 Performance in Clinical Settings

![unnamed](https://github.com/user-attachments/assets/edc0bcfe-e7cc-4d5a-b354-08106dc54930)

Figure 6. PRS construction. [13]

Studies have shown promising results in PRS applications across different contexts. The GWS PRS achieves an odds ratio of 1.27 (CI: 1.20-1.35, p < 0.001), while the eQTL PRS shows comparable performance with an OR of 1.24 (CI: 1.17-1.32, p < 0.001). Most notably, the combined PRS approach demonstrates enhanced predictive power with an OR of 1.37 (CI: 1.29-1.45, p < 0.001).

These scores have proven particularly valuable in disease risk prediction, ranging from cancer susceptibility to diabetes risk assessment. Their utility extends to identifying high-risk individuals who might benefit from enhanced screening protocols.

#### 2.3.4 Current Limitations in PRS

In the context of genetic prediction, raw GWAS results alone provide limited utility for predicting complex traits. A PRS addresses this by aggregating weak signals spread across many genomic regions. At its core, a PRS can be viewed as a machine learning task: taking an input x (one's genotype and potentially other variables) and producing an output y (the individual's predicted trait value).

A recent paper [14] outlines two main approaches to training PRS. The first involves using a sufficiently large cohort of genotyped and phenotyped individuals to train a PRS from scratch using standard machine-learning algorithms on individual-level data. The second, more common approach involves meta-analyzing summary statistics from published GWAS results into a linear model.

Despite over a decade of refinement with exponentially increasing cohort sizes, the predictive power of most PRS remains limited. While there have been successes with specific phenotypes like height, most diseases and clinically relevant traits fall short of genetics-based risk assessment potential. Current models explain only a small fraction of trait heritability and often don't reach clinical relevance. This happens mainly due to PRS being able to explain heritability only based on cis-genetic expression, while not accounting environment variance and other source of heritability (e.g. nonlinear effects).

The review discusses attempts to improve PRS using nonlinear models, including support vector machines, random forests, and deep neural networks. However, these attempts have generally failed to outperform simple linear models.

## Conclusion

eQTLs have enabled us in a better way to understand the underlying mechanisms behind complex traits by combining wide genotype information with expression data, enabling per variant analysis. This type of analysis empowers scientists to design further downstream wet lab experiments. Thus it is important that these methods are continued to be optimized to best represent biology.

eQTLs enable anonymous faster fine-mapping via sufficient statistics that can be utilized to find causal SNPs behind the genes. Current implementation of SuSiE has led to development of multi-ancestry fine-mapping methods such as SuSiEx [15], SuShiE [16] to improve fine-mapping for limited data in non-European populations. Additionally, recently a new method RSparsePro [17] emerged to tackle the limitation of SuSiE-RSS associated with allele flips in LD matrices.

## Citations

[1] Henriksson, J., Chen, X., Gomes, T., Ullah, U., Meyer, K. B., Miragaia, R., Duddy, G., Pramanik, J., Yusa, K., Lahesmaa, R., & Teichmann, S. A. (2019). Genome-wide CRISPR Screens in T Helper Cells Reveal Pervasive Crosstalk between Activation and Differentiation. Cell, 176(4), 882-896.e18. https://doi.org/10.1016/j.cell.2018.11.044

[2] Gallagher, M. D., & Chen-Plotkin, A. S. (2018). The Post-GWAS Era: From association to function. The American Journal of Human Genetics, 102(5), 717–730. https://doi.org/10.1016/j.ajhg.2018.04.002

[3] Visscher, P. M., Wray, N. R., Zhang, Q., Sklar, P., McCarthy, M. I., Brown, M. A., & Yang, J. (2017). 10 years of GWAS Discovery: Biology, Function, and Translation. The American Journal of Human Genetics, 101(1), 5–22. https://doi.org/10.1016/j.ajhg.2017.06.005

[4] Halldorsson, B. V., Palsson, G., Stefansson, O. A., Jonsson, H., Hardarson, M. T., Eggertsson, H. P., & Stefansson, K. (2019). Characterizing mutagenic effects of recombination through a sequence-level genetic map. Science, 363(6425), eaau1043. https://doi.org/10.1126/science.aau1043

[5] McVean, G. A. (2002). A genealogical interpretation of linkage disequilibrium. Genetics, 162(2), 987-991. https://doi.org/10.1093/genetics/162.2.987

[6] Song, Y. S. (2016). Na Li and Matthew Stephens on modeling linkage disequilibrium. Genetics, 203(3), 1005-1006. https://doi.org/10.1534/genetics.116.187955

[7] Chen, H., Majumdar, A., Wang, L., Kar, S., Brown, K., Feng, H., Turman, C., Dennis, J., Easton, D., Michailidou, K., Simard, J., Bishop, D., Cheng, I., Huyghe, J., Schmit, S., O'Mara, T., Spurdle, A., Gharahkhani, P., Schumacher, J., & Lindstrom, S. (2021). Large-scale cross-cancer fine-mapping of the 5p15.33 region reveals multiple independent signals. Human Genetics and Genomics Advances, 2, 100041. https://doi.org/10.1016/j.xhgg.2021.100041

[8] The International HapMap Consortium. (2005). A haplotype map of the human genome. Nature, 437, 1299-1320. https://doi.org/10.1038/nature04226

[9] Hormozdiari, F., Kostem, E., Kang, E. Y., Pasaniuc, B., & Eskin, E. (2014). Identifying Causal Variants at Loci with Multiple Signals of Association. Genetics, 198(2), 497–508. https://doi.org/10.1534/genetics.114.167908

[10] Benner, C., Spencer, C. C., Havulinna, A. S., Salomaa, V., Ripatti, S., & Pirinen, M. (2016). FINEMAP: efficient variable selection using summary data from genome-wide association studies. Bioinformatics, 32(10), 1493–1501. https://doi.org/10.1093/bioinformatics/btw018

[11] Wang, G., Sarkar, A., Carbonetto, P., & Stephens, M. (2020). A Simple New Approach to Variable Selection in Regression, with Application to Genetic Fine Mapping. Journal of the Royal Statistical Society Series B (Statistical Methodology), 82(5), 1273–1300. https://doi.org/10.1111/rssb.12388

[12] Zou, Y., Carbonetto, P., Wang, G., & Stephens, M. (2022). Fine-mapping from summary data with the "Sum of Single Effects" model. PLoS Genetics, 18(7), e1010299. https://doi.org/10.1371/journal.pgen.1010299

[13] Gabriel, A. a. G., Atkins, J. R., Penha, R. C. C., Smith-Byrne, K., Gaborieau, V., Voegele, C., Abedi-Ardekani, B., Milojevic, M., Olaso, R., Meyer, V., Boland, A., Deleuze, J. F., Zaridze, D., Mukeriya, A., Swiatkowska, B., Janout, V., Schejbalová, M., Mates, D., Stojšić, J., . . . McKay, J. D. (2022). Genetic analysis of lung cancer and the germline impact on somatic mutation burden. JNCI Journal of the National Cancer Institute, 114(8), 1159–1166. https://doi.org/10.1093/jnci/djac087

[14] Brandes, N., Weissbrod, O., & Linial, M. (2022). Open problems in human trait genetics. Genome Biology, 23(1). https://doi.org/10.1186/s13059-022-02697-9 

[15] Yuan, K., Longchamps, R. J., Pardiñas, A. F., Yu, M., Chen, T., Lin, S., Chen, Y., Lam, M., Liu, R., Xia, Y., Guo, Z., Shi, W., Shen, C., Daly, M. J., Neale, B. M., Feng, Y. A., Lin, Y., Chen, C., O'Donovan, M. C., . . . Huang, H. (2024). Fine-mapping across diverse ancestries drives the discovery of putative causal variants underlying human complex traits and diseases. Nature Genetics, 56(9), 1841–1850. https://doi.org/10.1038/s41588-024-01870-z

[16] Lu, Z., Wang, X., Carr, M., Kim, A., Gazal, S., Mohammadi, P., Wu, L., Gusev, A., Pirruccello, J., Kachuri, L., & Mancuso, N. (2024). Improved multi-ancestry fine-mapping identifies cis-regulatory variants underlying molecular traits and disease risk. medRxiv (Cold Spring Harbor Laboratory). https://doi.org/10.1101/2024.04.15.24305836 

[17] Zhang, W., Lu, T., Sladek, R., Dupuis, J., & Lettre, G. (2024). Robust fine-mapping in the presence of linkage disequilibrium mismatch. bioRxiv (Cold Spring Harbor Laboratory). https://doi.org/10.1101/2024.10.29.620968
