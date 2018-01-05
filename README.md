# motifOverlapR README 

## Contents
1. Introduction
2. What is motifOverlapR ?
3. How to Use Motif OverlapR ?
4. Data Sources
5. Cite us


## 1 Introudction

Transcription factors (TFs) are proteins that play pivotal roles in cellular development, differeintation and homeostatsis by binding to specific genomic sites and regulating nearby gene expression. Typically, to identify Transcription Factor Binding Sites (TFBS) one would conduct chromatin immunopreciptation followed by sequencing (ChIP-seq) which can be used in conjunction with transcriptomic analysis of gene expression to identify genes regulated by the transcription factor. However, this approach is not feasible for all TF's because many TFs are lack a sufficently high quality antibody that can detect endogenous expression. Consequently, there is a need for computational identification of TFBS that can be used in conjunction with the transcriptomic data to identify genes regulated by specific TFs in specific cell types. Here we present, motifOverlapR, a genomic tool that can predict TFBS through the integration of multiple genomic and epigenomic data types.  


## 2 What is motifOverlapR
MotifOverlapR is a computational bioinformatics tool which predicts TFBS by integration of genomic,  cell type specific epigneomic and transcriptic omic data types. motifOverlapR firstly identifies all DNA sequences that match the Position weight matrix of the Transcription Factor of interest, therefore identifying all possible motifs. once, identified, motifOverlapR identifies which of these motifs fall in regions of the genome that are found to play integral roles in gene regulation, such as promoter regions of enhancer regions. Consequence, of many motifs that play important roles in gene regulation are frequently conserved, motifOverlapR can also identify motifs in promoter sequences that fall into conserved regions (regions with a PhyloP score >0). Typically, enhancers evolve rapidly hence, the TFBS located in them are subject to higher turnover rates hence this conserved option is only available for promoter regions. To identify which sites are accessible to the TF of interest, motifOverlapR will identify which of the regulatory motif instances are located in active chromatin states in the cell type in which the TF is expressed. These sites are considered the unbaised predicted TFBS as they integrate genome wide data. An additional feature of motifOverlapR is to integrate transcriptomic data by uploading a list of genes that showed differential expression in the absence or mutation of the TF. motifOverlapR will return predicted TFBS that would regulate these genes expression. 


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

## 4 Data Sources

 Genomic Sequences were attained from the UCSC BS.genome Package
 Conservation data was downloaded from the UCSC FTP site
 Epigneomics data was downloaded off of the epigenomics roadmap using the annotationHub package
 Enhancer targets were assigned using the CAGE expression off of the fantom5 database
 
 


## 5 Cite Us

Awais Choudhry
