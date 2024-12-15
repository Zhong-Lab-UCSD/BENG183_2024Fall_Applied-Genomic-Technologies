Group 41

Shaheer Imran, Corey Nguyen, Somtochukwu Ikeanyi

# Controlling False Discovery Rate in RNA-Seq with Benjamini-Hochberg (BH)

1.  The Challenge of Multiple Hypothesis Testing in RNA-Seq

2.  What is False Discovery Rate (FDR)?

3.  Why Control FDR in RNA-Seq?

4.  Benjamini-Hochberg Procedure for FDR Control

5.  Applying Benjamini-Hochberg to RNA-Seq Data

6.  Example Use Case: Genome-Wide Association Study

7.  Limitations of Benjamini-Hochberg

8.  Conclusion


## 1. The Challenge of Multiple Hypothesis Testing in RNA-Seq

RNA-Seq is widely used to analyze gene expression by comparing transcript levels between different conditions, like treatment and control groups. While it’s incredibly useful, the method has a major statistical challenge: testing thousands of genes at the same time for differential expression increases the chance of Type I errors (false positives). In other words, the more tests you run, the higher the likelihood of finding results that look significant but aren’t actually meaningful.

The main issue comes down to finding the right balance between two goals:

-   Minimizing False Positives: Making sure that the results represent real biological differences and not just random noise.

-   Maintaining Statistical Power: Still being able to detect small but real changes in gene expression.


If this balance isn’t managed properly, it can lead to unreliable conclusions and wasted efforts chasing results that don’t hold up. Controlling false positives is key to ensuring that the insights gained from RNA-Seq are both accurate and useful.

## 2. What is False Discovery Rate (FDR)?

False Discovery Rate (FDR) is a statistical method used to tackle the problem of multiple hypothesis testing. It represents the expected proportion of false positives out of all the significant results: FDR=FP/(FP+TP). Compared to the Family-Wise Error Rate (FWER), which focuses on completely avoiding false positives, FDR takes a different approach by allowing a small, controlled proportion of false positives. This makes FDR more practical for large datasets like RNA-Seq, where strict FWER control could result in missing many true findings.

Using FDR ensures that the results are statistically reliable while still capturing meaningful discoveries. This balance is especially important in high-throughput experiments like RNA-Seq, where it’s crucial to detect true biological signals without getting bogged down by random noise or overly conservative methods.

## 3. Why Control FDR in RNA-Seq?

1.  Minimizing False Positives
	- 1.1 FWER vs BH
	- 1.2 BH vs Alpha
	- 1.3 Large-Scale Datasets
2. Preserving Statistical Power
	- 2.1 Still Detecting True Positives
	- 2.2 Finding Balance
3.  Practical Implications

Understanding the importance of controlling the False Discovery Rate (FDR) is crucial in research, especially in fields involving large-scale data analysis like genomics or other -omics technologies. FDR is a statistical method used to reduce the risk of type I errors (false positives) in multiple hypothesis testing.

### 1. Minimizing False Positives

#### 1.1 FWER vs BH (Family-Wise Error Rate vs Benjamini-Hochberg)

In the context of RNA-Seq data analysis, controlling the error rate is crucial to ensure the reliability of the conclusions drawn from the data. The Family-Wise Error Rate (FWER) and the Benjamini Hochberg (BH) procedure address this issue differently. FWER controls the probability of making even one type I error across all tests, making it quite stringent. This approach is typically used when the consequence of a single false discovery is particularly high. However, its stringency can lead to a substantial loss in statistical power, particularly problematic in high throughput settings like RNA-Seq.



In contrast, the Benjamini-Hochberg procedure controls the False Discovery Rate (FDR), the expected proportion of false discoveries among the rejected hypotheses. This is less conservative than FWER, allowing for more discoveries while still controlling the rate of type I errors. In RNA-Seq, where thousands of genes are tested simultaneously, BH is preferable as it balances error control and discovery rate.
![](https://lh7-rt.googleusercontent.com/docsz/AD_4nXfIxUg5h2qEZpKyJziL6-RFJ64qIBEkTmKaN_UFgAgAMEbwzU7oPnWJH1kp4JVxseoU73GR0gvawcvW61V1Ip_vs5hbEKP54WahOWC4AMZPuUxzO1-FdCuy9uV9Ooew8UEiw6Gk?key=9oICL3L3ve8OraBENXK-GbGB)

Figure: Diagram illustrating the concept of FDR versus traditional error measures (e.g., family-wise error rate). This diagram helps visualize the difference in error accumulation between methods.

#### 1.2 BH vs Alpha

The choice of alpha, the significance level in the BH procedure directly affects the number of hypotheses rejected. Setting a lower alpha can reduce the FDR, ensuring a more conservative approach and minimizing false positives. However, this can also diminish the number of true positives detected. RNA-Seq studies often involve multiple comparisons across a vast array of genes and a balance must be struck between controlling false positives and retaining enough power to detect true differentially expressed genes. The BH method provides a mechanism to adjust p–values based on the rank and the total number of tests, thus tailoring the alpha to the dataset’s specific characteristics.



![](https://lh7-rt.googleusercontent.com/docsz/AD_4nXd-K-Kz3qCq0LrJ0paHQubGzJzwLPFf0Unb0--GODjIdTEQP6aUcYYUILx8mv_WOhXXyUn-_iNyuQC0RuTdCmD4w_XDgjF5RwU6vpAnk1SfQpq9OKjfBdgV-TouZ3ojMlA3exU_NQ?key=9oICL3L3ve8OraBENXK-GbGB)![](https://lh7-rt.googleusercontent.com/docsz/AD_4nXdRJ3HNQIkFePfACJeP-pC7bJKZUhJlzy_-nxvOYQnygdL1STd59HHkPZSCl-w8CX-TyDENmeXapgFQaidRp-oBEwMgg-XztQ6KDmXfPMA581kqOjMIAaJfWT6YOJar2rzVe0to_A?key=9oICL3L3ve8OraBENXK-GbGB)

Figure: FDR control in in silico experiments and simulations. Observed FDR (y-axis) for various α-level cutoffs (x-axis). Demonstrates FDR decreases with alpha level.


#### 1.3 Large-Scale Datasets

RNA-Seq is typically characterized by large-scale datasets, often involving comparisons of thousands of genes across multiple conditions or time points. The sheer volume of simultaneous tests introduces a high chance of encountering false positives. The BH procedure is particularly well-suited for these scenarios because it efficiently controls the FDR while accommodating the large scale of the data. This control is crucial for making valid inferences in studies where the cost of follow-up experiments is high, ensuring that the findings are both statistically and scientifically robust.

### 2. Preserving Statistical Power

#### 2.1 Still Detecting True Positives

One of the primary benefits of the BH procedure over more stringent methods like Bonferroni is its ability to preserve statistical power. By controlling the FDR rather than the FWER, it avoids overly penalizing the dataset for the number of comparisons. This is particularly important in RNA-Seq analyses, where the biological significance often hinges on the ability to detect subtle but true changes in gene expression across conditions. The BH method ensures that researchers can still identify these true positives without an overwhelming risk of false positives.

#### 2.2 Finding Balance

The key advantage of the BH procedure in RNA-Seq data analysis lies in its ability to find a balance between controlling false positives and maintaining the ability to detect true positives. By adjusting p-values based on their rank and the total number of tests, the BH method can adapt to the data's characteristics, providing a flexible yet robust tool for hypothesis testing. This balance is crucial for advancing biological and medical research where confirming new or unexpected results can lead to significant advancements.

![](https://lh7-rt.googleusercontent.com/docsz/AD_4nXcFDAkqDc15G0UeBdgTuwIrAHKHYAYG0c50j9SSps19Gc8bzfvB_KFAefwj4ei-wdSjTbiwsiAaNbulWmr2GrS8AHlIt-PGxK3gz6TD7gwwK_qdukpYXZe2is1w_MxiRGtYJi6ZxQ?key=9oICL3L3ve8OraBENXK-GbGB)

Figure: Graph plotting p-values against their rank. Top dashed line indicates alpha. Bottom dashed line indicates FWER. Solid line indicates BH. Demonstrates both statistical power and control of FDR by Benjamini Hochberg.

### 3. Practical Implications

In practical terms, applying the BH procedure in RNA-Seq studies means researchers can trust the results to reflect a reasonable compromise between sensitivity (detecting true positives) and specificity (not claiming false positives). This trust is crucial for subsequent stages of research, including validation studies, biological interpretation, and therapeutic targeting. By providing a reliable method of controlling the error rate, the BH procedure helps pave the way for these applications, ensuring that findings from RNA-Seq data are both scientifically valid and actionable.

Reference: Joshua N Sampson, Simina M Boca, Steven C Moore, Ruth Heller, FWER and FDR control when testing multiple mediators, Bioinformatics, Volume 34, Issue 14, July 2018, Pages 2418–2424, https://doi.org/10.1093/bioinformatics/bty064

## 4. Benjamini-Hochberg Procedure for FDR Control

This section breaks down the Benjamini-Hochberg procedure into detailed steps, providing a clear guide for its application.

1.  Ranking P-values

2.  Calculating Adjusted P-values

3.  Adjusting P-values Sequentially

4.  Rejecting Null Hypotheses

### 1. Ranking P-values

The initial step in the Benjamini-Hochberg procedure is to rank all obtained p-values in ascending order. This ranking is pivotal because the correction process depends on the relative positioning of each p-value within the dataset. By ordering them from the smallest to the largest, the procedure can prioritize the assessment of findings most likely to be significant, ensuring that the most compelling results are evaluated with the most stringent criteria. This step is fundamental to controlling the FDR because it sets the foundation for sequential adjustments that account for the multiplicity of the tests. Moreover, this method acknowledges that not all p-values contribute equally to the family-wise error rate, giving precedence to those most likely to represent true discoveries.

### 2. Calculating Adjusted P-values

Following the ranking of p-values, the next step is to calculate the adjusted p-values using the formula (m/i)×p(i), where m is the total number of tests, and ii is the rank of the individual p-value. This formula adjusts each p-value by scaling it according to its rank relative to the total number of tests. The multiplication of the original p-value by the factor (m/i) serves to inflate the p-value of each test based on its position in the ordered list. This inflation is crucial as it compensates for the multiple testing problem by controlling the expected proportion of incorrectly rejected null hypotheses, thereby maintaining the false discovery rate at or below a desired level, typically set at 0.05 or 0.01.

### 3. Adjusting P-values Sequentially

The sequential adjustment of p-values is a distinctive feature of the Benjamini-Hochberg procedure. Starting from the highest ranked p-value (smallest p-value), each subsequent p-value is adjusted to be no smaller than any that have come before it. This ensures that the adjustment process is monotonically increasing, which is crucial to maintaining the integrity of the FDR control across all tests. This method guards against the possibility that a less significant result (a higher p-value) receives a disproportionately low adjusted p-value due to the sequential nature of the adjustment process.

### 4. Rejecting Null Hypotheses

The final decision in the Benjamini-Hochberg procedure involves rejecting null hypotheses for all tests where the adjusted p-values are less than or equal to a predefined threshold α, typically set at 0.05. This threshold represents the maximum proportion of false discoveries (false positives) that the researcher is willing to tolerate. By setting this criterion, the procedure ensures that the overall rate of type I errors (false positives) among all significant tests is controlled, making the findings statistically robust.

![](https://lh7-rt.googleusercontent.com/docsz/AD_4nXdOqmjqlu1ZdqQTTf3aQGaiPEImoBwFTkuJvJCRug0i81drOQBTgpXmySYwCRi11r1aGwrYNO-NaXQE6KmVjiiX2zn9TfYZB_gNfVHjYRFa8rJll3u3qIuJ1Jfg21Av3b3TC8t0?key=9oICL3L3ve8OraBENXK-GbGB)

Figure: Flowchart illustrating the sequential steps of the Benjamini-Hochberg procedure.

Reference: Benjamini, Y., & Hochberg, Y. (1995). Controlling the false discovery rate: a practical and powerful approach to multiple testing. Journal of the Royal Statistical Society. Series B (Methodological), 289-300.

## 5. Applying Benjamini-Hochberg to RNA-Seq Data

This section explains the role of the Benjamini-Hochberg procedure in RNA-Seq studies. It highlights the importance of controlling the False Discovery Rate (FDR) when dealing with large datasets typical in RNA-Seq to ensure the accuracy and reliability of the gene expression data being analyzed.

1.  Extracting P-values from RNA-Seq Data

2.  Adjusting P-values Using the BH Method

3.  Interpreting Adjusted Results

4.  Consequences of Applying FDR Control in RNA-Seq


### 1. Extracting P-values from RNA-Seq Data

In RNA-Seq analysis, extracting p-values typically involves using statistical tools like DESeq2 or edgeR. These tools analyze differences in gene expression between conditions or groups by modeling the count data, thus providing p-values for each gene.

### 2. Adjusting P-values Using the Benjamini-Hochberg Method

The procedure starts by ranking all the obtained p-values in ascending order. Each p-value is then adjusted according to its rank to control the FDR, ensuring that the proportion of false positives remains below a pre-specified threshold, commonly set at 0.05.

### 3. Interpreting Adjusted Results

After applying the Benjamini-Hochberg adjustment, the next step is to interpret these adjusted p-values. This involves identifying which genes show statistically significant differences in expression. This selection helps direct further research and validation efforts towards those genes most likely to be biologically significant.

### 4. Consequences of Applying FDR Control in RNA-Seq

By controlling the FDR, researchers can improve the reliability of their findings. This minimizes the risk of pursuing false leads based on erroneous data, thus conserving resources and focusing efforts on the most promising hypotheses.  Applying the Benjamini-Hochberg procedure can substantially affect the conclusions drawn from RNA-Seq data. By rigorously controlling the rate of false discoveries, researchers can ensure that their results are robust and replicable, enhancing the scientific value of their work.

Reference: Anders, Simon, and Wolfgang Huber. "Differential expression of RNA-Seq data at the gene level–the DESeq package." Heidelberg, Germany: European Molecular Biology Laboratory (EMBL) 10 (2012): f1000research.

## 6. Alternative Scenario: Genome-Wide Association Study

The Benjamini-Hochberg procedure is also widely applied in other high-dimensional genomics studies, such as Genome-Wide Association Studies (GWAS). In a typical GWAS, millions of single nucleotide polymorphisms (SNPs) are tested for association with a disease or phenotype. The steps are similar to RNA-Seq applications:

1.  Compute the p-value for each SNP.

2.  Sort the SNPs by p-value in ascending order.

3.  Calculate the adjusted p-values using the BH procedure.

4.  Determine which SNPs have adjusted p-values below a predefined significance threshold (e.g., 0.01).

By controlling the FDR, the researcher can identify SNPs that are most likely to be truly associated with the disease, while minimizing the inclusion of false positive findings. Suppose one finds that 200 SNPs meet the threshold; An FDR control of 1% ensures only 2 are false positives. It cannot be overstated how important ensuring the validity of one’s results are in scientific experimentation and Benjamini-Hochberg is a tool which can serve this purpose in many different projects the bioinformatician works on.

## 7. Limitations of Benjamini-Hochberg

While the Benjamini-Hochberg procedure is widely used, it does have limitations. The method assumes that the tests are independent or positively dependent. In RNA-Seq, however, many genes are co-expressed, meaning that their expression levels are correlated. In such cases, the BH procedure may be too lenient, as the correlation between tests violates the assumption of independence.

Additionally, while the BH procedure controls FDR, it does not eliminate false positives entirely. Researchers should therefore interpret the results with caution, considering the possibility of false discoveries even after FDR control. In some situations, alternative methods, such as Storey’s q-value approach, may be more appropriate, particularly when the tests exhibit complex dependencies.

## 8. Conclusion

Controlling the False Discovery Rate (FDR) is essential in RNA-Seq experiments to ensure that the results reflect genuine biological findings, rather than artifacts of multiple hypothesis testing. The Benjamini-Hochberg procedure offers an effective and widely used method for controlling FDR while maintaining statistical power. By adjusting p-values and limiting false discoveries, the BH procedure enables researchers to identify significant genes without over-interpreting insignificant ones. Despite its limitations, it remains an invaluable tool for RNA-Seq and other high-throughput genomic studies, helping to guide researchers toward meaningful, reproducible conclusions.
