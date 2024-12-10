# 1.1 MEME (Multiple Expectation maximizations for Motif Elicitation) 
1. [Introduction](#111)<br>
    1.1. [About MEME Suite](#1111)<br>
    1.2. [Importance of Motif-based Analysis](#1112)
2. [Overview of MEME Parameters](#112)<br>
    2.1. [Motif Discovery Modes](#1121)<br>
3. [Algorithm](#113)
5. [Example](#114)
6. [Limitations](115)

## 1.1.1 Introduction<a name="111"></a>

#### About MEME Suite<a name="1111"></a>
MEME Suite is an extensive collection of motif-based sequence analysis tools, widely used in the scientific community for identifying and analyzing motifs found in biological sciences. The basis of algorithms behind the MEME Suite tools were first published in 1994 by Timothy L. Bailey and Charles Elkan from UC San Diego's Computer Science and Engineering department. Its tools are specialized for motif discoveries, enrichments, scannings, comparisons, and more. This chapter will be focused on MEME, which is a motif discovery tool within MEME Suite. 

#### Importance of Motif-based Analysis<a name="1112"></a>
Motifs are short repetitive patterns in DNA, RNA, and protein sequences that are biologically significant. They play a vital role in understanding biological processes because they often correspond to functional or regulatory elements within the genome. Two types of motifs are: gapped motifs, which are recurring, variable length motifs, and ungapped motifs, which are recurring, fixed-length motifs. For example, the TATA box, an important promoter sequence, is an ungapped motif of fixed length, whereas the bacterial σ70 promoter: TTGACA– (16-19 bp gap) –TATAAT, is a gapped motif of variable length.

## 1.1.2 Overview of MEME Parameters<a name="112"></a>

The MEME tool’s basic function is to discover ungapped motifs in inputted group(s) of unaligned sequences through either the Classic or the Discriminative/Differential mode. It also requires the selection of specific parameters to improve the quality of the motif search. In general, all modes will require and generate the following: 
- Input: Group of related sequences (DNA, RNA, or protein) in FastA/BED format
- Output: As many motifs as requested and graphed in the pictogram

#### Motif Discovery Modes<a name="1121"></a>
<table>
 <tbody>
    <tr>
        <th>Motif Discovery Mode</td>
        <th>Description</td>
        <th>Basic Usage</td>
    </tr>
    <tr>
        <td>Classic Mode</td>
        <td>Site distribution: Informs MEME of expected motif distribution to improve quality of motif search</td>
        <td><ul><li>Input primary sequences: Provide one set of sequences and MEME discovers motifs enriched in this set. Enrichment is measured relative to a (higher order) random model based on frequencies of the letters in your sequences</li></ul></td>
    </tr>
    <tr>
    <td>Discriminative Mode</td>
    <td><ul><li>Searching for motifs that appear more in the target than the control set</li></ul><ul><li>Position-specific priors (PSPs) are generated for the target set, which are used to focus MEME on patterns in target set that don’t appear in control</li></ul></td>
    <td><ul><li>Provide two sets of sequences: primary (target), secondary (control)</li></ul></td>
    </tr>
    <tr>
    <td>Differential Mode</td>
    <td><ul><li>Searching for motifs that appear more in the target than the control set</li></ul><ul><li>Simply counts how often motifs appear in the target set compared to the control. It gives higher score to motifs that appear more often appear in the target set</li></ul></td>
    <td><ul><li>Provide two sets of sequences: primary (target), secondary (control)</li></ul></td>
    </tr>
 </tbody>
</table>

## 1.1.3 MEME Algorithm<a name="113"></a>

## 1.1.4 Example of Usage<a name="114"></a>

## 1.1.5 Limitations<a name="115"></a>

# References
