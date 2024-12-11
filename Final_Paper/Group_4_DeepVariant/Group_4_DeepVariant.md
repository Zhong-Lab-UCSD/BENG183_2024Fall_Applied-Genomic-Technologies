# Handout: DeepVariant - Transforming Genomics with AI

![DeepVariant Logo](Figures/dv_logo.png)

---

## 1. Introduction
### 1.1 Background
### 1.2 Importance of Variant Calling

---

## 2. The Challenge of Variant Calling
### 2.1 Sources of Difficulty
### 2.2 DeepVariant’s Solution

---

## 3. Deep Learning and CNNs
### 3.1 Deep Learning Overview
### 3.2 Convolutional Neural Networks (CNNs)

---

## 4. Performance
### 4.1 Key Metrics
- **Primary Training and Testing**:
  - Trained on the CEPH NA12878 sample and tested on the unseen Ashkenazi male NA24385 sample.
  - Achieved:
    - **SNP F1 Score**: 99.95%
    - **Indel F1 Score**: 98.98%
- **Error Reduction**:
  - Demonstrated over 50% fewer errors per genome (4,652 errors) compared to GATK, FreeBayes, SAMtools, 16GT, and Strelka, where the next-best method had 9,531 errors.

### 4.2 Robustness Across Applications
- **Genome Builds**:
  - Nearly identical performance across GRCh37 and GRCh38:
    - GRCh37 → GRCh38: F1 = 99.45%
    - GRCh38 → GRCh38: F1 = 99.53%
  - Highlights generalizability across different genome references.
- **Cross-Species Generalization**:
  - Human-trained model achieved F1 = 98.29% on mouse data, outperforming mouse-trained models (F1 = 97.84%).
  - Demonstrates adaptability to species differences and varying sequencing parameters.
- **Sequencing Platforms**:
  - Robust to sequencing technologies:
    - Illumina TruSeq (50× 2 × 148 bp reads on HiSeq 2500).
    - Custom mouse sequencing (27× 2 × 100 bp reads on Genome Analyzer II).

---

---

## 5. Challenges and Improvements
### 5.1 Challenges
- **Exome Datasets**:
  - Lower performance due to:
    - Fewer variants (20k–30k vs. 4–5M in whole genomes).
    - Non-uniform coverage and errors from exome capture or amplification technologies.
  - High false positive rates initially, reducing the model’s accuracy.

### 5.2 Improvements
- **Whole-Genome Sequencing (WGS)**:
  - Retrained models improved precision (PPV):
    - **SOLiD**: PPV = 99.0% (candidates 82.5%, final 76.6% sensitivity).
    - **PacBio**: PPV = 97.3% (candidates 93.4%, final 88.5% sensitivity).
- **Exome Performance**:
  - Significant PPV improvements after retraining:
    - **Ion Ampliseq**: Increased from 8.1% to 99.7%.
    - **TruSeq**: Increased from 65.3% to 99.3%.
  - Sensitivity reductions were small but manageable:
    - Ion Ampliseq: 91.9% → 89.3%.
    - TruSeq: 94.0% → 92.6%.
- **Preprocessing Enhancements**:
  - Additional preprocessing steps contributed to accuracy improvements, reducing false positives.


---

## 6. How DeepVariant Works
### 6.1 Analytical Pipeline
#### Step 1: Candidate Variant Identification
#### Step 2: Encode Variants as Pileup Images
#### Step 3: CNN Processing
#### Step 4: Variant Calling
#### Step 5: Accessing DeepVariant
- The DeepVariant pipeline is open-source and available on GitHub: [DeepVariant Repository](https://github.com/google/deepvariant).

---

## 7. Advantages of DeepVariant
- **Accuracy**:
  - Near-perfect F1 scores for SNPs and indels.
  - 50% fewer errors compared to other tools.
- **Generalizability**:
  - Consistent performance across genome builds, species, and sequencing parameters.
- **Adaptability**:
  - Demonstrates robust cross-species application (e.g., human to mouse).
- **Versatility**:
  - Compatible with various sequencing platforms and experimental protocols.
- **Efficiency**:
  - Reduces reliance on handcrafted statistical models through deep learning.

---

## 8. Conclusion

---

## References