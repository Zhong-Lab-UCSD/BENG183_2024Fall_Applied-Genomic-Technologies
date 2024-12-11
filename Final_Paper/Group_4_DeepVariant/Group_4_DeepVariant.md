# DeepVariant: Transforming Genomics with AI

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

## 3. How DeepVariant Works
### 3.1 Analytical Pipeline
#### Step 1: Candidate Variant Identification
#### Step 2: Encode Variants as Pileup Images
#### Step 3: CNN Processing
#### Step 4: Variant Calling
#### Step 5: Accessing DeepVariant

---

## 4. Performance
### 4.1 Key Metrics
- **Primary Training and Testing**:
  - Trained on CEPH NA12878 and tested on unseen Ashkenazi male NA24385.
  - Achieved:
    - **SNP F1 Score**: 99.95%.
    - **Indel F1 Score**: 98.98%.
- **Error Reduction**:
  - Over 50% fewer errors compared to GATK, FreeBayes, SAMtools, 16GT, and Strelka.

### 4.2 Robustness Across Applications
- **Genome Builds**:
  - GRCh37 → GRCh38: F1 = 99.45%.
  - GRCh38 → GRCh38: F1 = 99.53%.
  - Demonstrates generalizability across references.
- **Cross-Species Generalization**:
  - Human-trained model achieved F1 = 98.29% on mouse data.
  - Outperformed mouse-trained models (F1 = 97.84%).
- **Sequencing Platforms**:
  - Performs well with:
    - Illumina TruSeq (50× 2 × 148 bp reads on HiSeq 2500).
    - Custom mouse sequencing (27× 2 × 100 bp reads on Genome Analyzer II).

---

## 5. Challenges and Improvements
### 5.1 Challenges
- Exome datasets present difficulties:
  - Fewer variants (~20k–30k vs. ~4–5M in whole genomes).
  - Non-uniform coverage and sequencing errors increase false positives.

### 5.2 Improvements
- **Whole-Genome Sequencing (WGS)**:
  - Retrained models achieved:
    - **SOLiD**: PPV = 99.0% (candidates 82.5%, final 76.6% sensitivity).
    - **PacBio**: PPV = 97.3% (candidates 93.4%, final 88.5% sensitivity).
- **Exome Performance**:
  - Ion Ampliseq PPV improved from 8.1% to 99.7%.
  - TruSeq PPV improved from 65.3% to 99.3%.
  - Small sensitivity reductions:
    - Ion Ampliseq: 91.9% → 89.3%.
    - TruSeq: 94.0% → 92.6%.
- **Preprocessing**:
  - Preprocessing steps contributed to significant accuracy improvements.

---

## 6. Advantages of DeepVariant
- **Accuracy**:
  - Near-perfect F1 scores for SNPs and indels.
  - 50% fewer errors than other tools.
- **Generalizability**:
  - Consistent performance across genome builds, species, and platforms.
- **Adaptability**:
  - Effective cross-species applications, e.g., human to mouse.
- **Versatility**:
  - Compatible with various sequencing platforms and protocols.
- **Efficiency**:
  - Eliminates reliance on handcrafted statistical models through deep learning.

---

## 7. Future Directions

---

## 8. Conclusion

---

## References
- DeepVariant GitHub: [https://github.com/google/deepvariant](https://github.com/google/deepvariant)
- Nature Biotechnology Paper: [DeepVariant Paper](https://www.nature.com/articles/nbt.4235)