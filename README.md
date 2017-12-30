# motifOverlapR README 

##Contents
1. Introduction
2. What is motifOverlapR ?
3. How to Use Motif OverlapR ?
4. Cite us


## 1 Introudction

Transcription factors (TFs) are proteins that play pivotal roles in cell maintence, development, differeintation by regulation gene transcription in a sequence dependent manner. Typically, to identify TF binding sites one would conduct chromatin immunopreciptation followed by sequencing (ChIP-seq) however, due to the numerous TF's and cell types in which each TF is expressed this is not feasible. Additionally, many TF's are expressed lowly hence in various cell types thus making it difficult to study them using ChIP protocols due to lack of highly specific antibody. Therefore, due to lack of antibody and numeorus cell type-TF combinations it is not always possible to conduct ChIP-seq to identify TFBS. 


## 2 What is motifOverlapR
MotifOverlapR is a computational bioinformatics tool developed that integrates genomic, epigneomic and transcriptomic data to predict TFBS in specific cell types. Due to the 


## 3 How to Use Motif OverlapR
1. Select the transcription factor of interest's motif
  a. From the drop down menu of motifs on the Jaspar Database
  b. Upload a position weight matrix in the same format as the Jaspar database # ** see link **
  c. Type in a DNA sequence that represents the motif of the TF of interest
2. Select how accurately/stringently you wish to identify DNA sequences that match the selected motif
3. Increasing stringency by identifying motifs in conserved promoter regions
4. Select the cell type/Tissue in which the to computationally predict the TFBS
5. Upload a list of genes that showed differential expression in the Knockout/Knockdown of the Transcription Factor
6. Interrogate using the genome browser or export the results to a bedfile

## 4 Cite Us
