
Welcome to the motifOverlapR wiki!
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


## 3 How to Use motifOverlapR

#### The Contents
motifOverlapR contains 3 tabs as illustrated in the figure below
1. The computational Pipeline
2. The genome browser and motif Table
3. Gene Ontology
![](https://raw.githubusercontent.com/awaisc/motifOverlapR/master/ImagesForReadMe/motifOverlapRLayout.jpg)

#### The Computational Pipeline
In the first tab is the computational pipeline and two data tables containing the results of the computational pipeline.

There are 5 total inputs
Input 1 - the Position Weight Matrix to match to the genome. There are 3 seperate options for this.

  1. Type in the DNA sequence (No spaces)
  2. Upload a custom position weight matrix in the JASPAR format (including the header) in a csv format
  3. Select an existing PWM matrix from the JASPAR database

Input 2. - How accuretly you want the DNA sequence to match the position weight matrix selected
This is a numeric value from 1-100 whereby 100 will identify motifs that match the position weight matrix perfectly.

Input 3. - Identify motifs that are located in conserved promoter regions between vertebrate genomes.
This is a checkbox that when selected will identify motifs that are in promoter regions that are falling into genomic regions that are completely conserved. 

Please be aware, this will download a 9GB file from UCSC that contains the PhyloP scores for the entire hg19 genome. Hence, this option may take hours to load the first time.

Input 4 - Select the cell type/Tissue type to predict the sites within
This input will determine the epigenomic environment for each cell type hence allowing us to identify which motifs are located in active chromatin regions and therefore accessible to the TF. Typically, this will be the cell type in which the TF is expressed.

Input 5 - Optional: Upload a list of genes that showed differential expression in the absence of the TF. 
This is an optional input whereby uploading a list of genes (in the gens symbols format) will limit predicted TFBS to regulatory modules (ie promoters and enhancers) that regulate genes that showed differential expression.  

![](https://raw.githubusercontent.com/awaisc/motifOverlapR/master/ImagesForReadMe/ComputationalTabResults.jpg)


#### Results Tables

The raw results are illustrated as two tables based on the regulatory module they're found in: either promoters or enhancers.

Both tables contain the following information
1. Genomic coordinates of the motif
2. The score of the motif according to the position weight matrix.
3. The DNA String
4. The Gene regulated (as both a gene symbol and a UCSC transcript ID). 


![](https://raw.githubusercontent.com/awaisc/motifOverlapR/master/ImagesForReadMe/ComputationalTabResults.jpg)



These Tables are filterable by each of these columns to select for specific motifs. 

![](https://raw.githubusercontent.com/awaisc/motifOverlapR/master/ImagesForReadMe/ComputationalTabResultsFilterAbility.jpg)


#### Genome Browser Tab/ Table

 motifOverlapR contains an inbuilt genome browser that contains the results of the computational pipeline (therefore the two data tables) to interrogate the results further and more freely (containing a raw motif instance track as well). 

The inputs for the genome browser are the genomic coordinates you wish to view and clicking go. Please be aware, this is computationally intensive hence requires a lot of ram. Thus, if you are running into issues, please try reducing the size of the region you're interrogating.

![](https://raw.githubusercontent.com/awaisc/motifOverlapR/master/ImagesForReadMe/GenomeBrowserTabInputs.jpg)

The predicted motif track correlates with the datatable below showing the information regarding these motifs.

![](https://raw.githubusercontent.com/awaisc/motifOverlapR/master/ImagesForReadMe/GenomeBrowserTabPredictedSitesTable.jpg)

The colour of the Selected cell type track correlates with the chromHMM legend, informing you of the epigenomic state for each genomic region. This is illustrated in the figure where the dark green region correlates with active transcription and the red region is an active transcription start site.  

![](https://raw.githubusercontent.com/awaisc/motifOverlapR/master/ImagesForReadMe/GenomeBrowserTabChromHMMLegend.jpg)


#### Gene Ontology Tab

motifOverlapR also contains a gene ontology tab that conducts an enrichment based on the UCSC gene transcripts that are predicted to be regulated by the predicted TFBS.

Here, we can see many activating ontology terms for the Arnt transcription factor.

![](https://raw.githubusercontent.com/awaisc/motifOverlapR/master/ImagesForReadMe/GeneOntologyTabResults.jpg)



## 4 Data Sources

 Genomic Sequences were attained from the UCSC BS.genome Package
 Conservation data was downloaded from the UCSC FTP site
 Epigneomics data was downloaded off of the epigenomics roadmap using the annotationHub package
 Enhancer targets were assigned using the CAGE expression off of the fantom5 database
 
 


## 5 Cite Us

Awais Choudhry
