#
# This is the server logic of a Shiny web application. You can run the 
# application by clicking 'Run App' above.
#
# Find out more about building applications with Shiny here:
# 
#    http://shiny.rstudio.com/
#

require(shiny)
require(Gviz)
require(GenomicInteractions)
require(magrittr)
require(TxDb.Hsapiens.UCSC.hg19.knownGene)
require(org.Hs.eg.db)
require(BSgenome.Hsapiens.UCSC.hg19)
require(readxl)
require(dplyr)
require(tidyr)
require(JASPAR2016)
require(TFBSTools)
require(readr)
require(AnnotationHub)
require(seqLogo)
require(readr)
require(limma)
require(biomaRt)


shinyServer(function(input, output) {


  #######################################################
  ### Identifying which are located in enhancer regions and which are located in promoter regions
  ### And Identifying where they are located in human genome
  #######################################################
  ##############################################
  #Render Human GVIZ plot 1
  #############################################33
  
  JASPARR2016Matrices <- getMatrixSet(JASPAR2016, list() )
  namedJasparDataBase<-lapply(1:length(JASPARR2016Matrices), function(x){
    names(JASPARR2016Matrices[x])<-JASPARR2016Matrices[[x]]@name  
  })
  
  names(JASPARR2016Matrices)<-namedJasparDataBase
  
  ## Visualize your motif
  
  output$VisualizeTFMotif<-renderPlot({
    
    if(input$CustomTFMatrix==TRUE){
      
      PathToCustomMatrix<-input$JasparCustomMatrix
    
      #Import and convert to custom Matrix
      TextTranscriptionFactor <- read_delim(PathToCustomMatrix$datapath, 
                                            "\t", escape_double = FALSE, col_names = FALSE, 
                                            trim_ws = TRUE, skip = 1)%>%as.matrix()
      
      
      #Sum the first column and divide all numbers by that to normalize and make it a frequency matrix.
      
      PWMTextTranscriptionFactor<-(TextTranscriptionFactor/
                                     TextTranscriptionFactor%>%
                                     colSums())%>%as.matrix()
      
      assign("matrixForMatching", TextTranscriptionFactor, .GlobalEnv)
      
      seqLogo(PWMTextTranscriptionFactor)
      
      
  } else if(input$TypeInSequence==TRUE){
     
    
    TranscriptionFactor<-round(PWM(input$CustomDNASequence)*dim((PWM(input$CustomDNASequence)))[2]) %T>%seqLogo()
  
    assign("matrixForMatching", TranscriptionFactor, .GlobalEnv)
    
    } else {
   
    ## Rewrite this because some of these matrices dont have even numbers of nucelotides for eahc position #WTF 
      # OKay so the issue appears to be with repeat numbers Eg 0.333333333333333333333. NOt sure how to deal with this. 
      
    (JASPARR2016Matrices[[paste0(input$TranscriptionFactorPWM)]]@profileMatrix/
       colSums(JASPARR2016Matrices[[paste0(input$TranscriptionFactorPWM)]]@profileMatrix))%>%
        as.matrix()%>%
        seqLogo()
    
      JasparMatrixForMatching<-JASPARR2016Matrices[[paste0(input$TranscriptionFactorPWM)]]@profileMatrix
      assign("matrixForMatching", JasparMatrixForMatching, .GlobalEnv)
  }
    })
  
  
  
  
  output$HumangvizPlot <- renderDataTable({
    
    input$action
    
    isolate({
     
     ########################################
     ###Data Inputs
     #######################################
     
     ##Importing Bed File with Granges
     promoterTracks <- read.delim("../DataFiles/Gene Tracks/Human/hg19bedWithNames.bed.gz")%>%
       makeGRangesFromDataFrame(
         keep.extra.columns=TRUE,
         ignore.strand=FALSE,
         seqinfo=NULL,
         seqnames.field="hg19.knownGene.chrom",
         start.field= "hg19.knownGene.txStart",
         end.field="hg19.knownGene.txEnd",
         strand.field="hg19.knownGene.strand",
         starts.in.df.are.0based=FALSE)%>%promoters()
     assign("promoterTracks", promoterTracks, .GlobalEnv)
     
     
     
     # ########################################################
     ### Using the TXDB object to make a Promoter Grange with Gene symbols instead of importing.
     #Limited because it does not get every gene symbol name (inparticular those in the scaffolds)
     ###########################################################
     # GENEIDToSymbol <- cbind.data.frame(x[mappedkeys(org.Hs.egSYMBOL)])
     # 
     # GeneIDDataFrame<-as.data.frame(Knownegne)
     # 
     # GeneIDDataFrame$GENEID<-GeneIDDataFrame$GENEID%>%as.character()
     # 
     # promoterTracks<-left_join(GeneIDDataFrame, GENEIDToSymbol,
     #                                by= c("GENEID" ="gene_id"))%>%
     #   makeGRangesFromDataFrame(., keep.extra.columns = TRUE)%>%promoters()
     # mcols(promoterTracks) <- cbind.data.frame( "hg19.kgXref.geneSymbol" =mcols(promoterTracks)$symbol)
     # 
     ### Pipe line for identifying potential TFBS
     
     
     # Epigenomic Chromatin state 
     
     ah<-AnnotationHub()
     
     assign("ah", ah, .GlobalEnv)
     ###################################
     ###Promoter Enhancer assoication data table Input
     ###################################
         #Read in the table
     EnhancerPromoterAssoications <-read_delim("../DataFiles/Enhancers Track/Human/hg19_enhancer_promoter_correlations_distances_cell_type.txt", 
                                               "\t", escape_double = FALSE, trim_ws = TRUE)
     assign("EnhancerPromoterAssoications", EnhancerPromoterAssoications, .GlobalEnv)
     
     
     
     PromotersAssoicatedWithEnhancers<-separate(EnhancerPromoterAssoications, col=promoter, into= c("promoter", "strand"), sep= ',')%>%
       separate(., col= promoter, into= c("chr", "start", "end"))%>%makeGRangesFromDataFrame(keep.extra.columns = TRUE)
     
     
     
     
     
     EnhancerGrangeWithTargets<-EnhancerPromoterAssoications$enhancer%>%GRanges()
     
     mcols(EnhancerGrangeWithTargets)<-PromotersAssoicatedWithEnhancers
     
     
     EnhancerGrangeWithTargetsSiginificant <- subset(EnhancerGrangeWithTargets, `FDR`!=0 & `FDR`<0.05)
     
     
     IntersectBetweenOverlappingPromoterRanges<-findOverlaps(EnhancerGrangeWithTargetsSiginificant$X, promoterTracks)
     
     GenesRegulatedByEnhancers<-as.data.frame(mcols(promoterTracks))[IntersectBetweenOverlappingPromoterRanges%>%subjectHits(),]
     
     EnhancersWithGeneTargetsGrange<-EnhancerGrangeWithTargetsSiginificant[IntersectBetweenOverlappingPromoterRanges%>%queryHits()]
     
     mcols(EnhancersWithGeneTargetsGrange) <- cbind.data.frame(mcols(EnhancersWithGeneTargetsGrange),GenesRegulatedByEnhancers)
     mcols(EnhancersWithGeneTargetsGrange) <- mcols(EnhancersWithGeneTargetsGrange)[12:20] 
     
     #We cn use this to identify targets of motifs in enhancer regions
     assign("EnhancersWithGeneTargetsGrange", EnhancersWithGeneTargetsGrange, .GlobalEnv)
     
     
     
     #ChromHMM Function from the coMet package (Edited for speed)
     ##ChromHMM Track Generator specifically for humans
     chromHMMTrackGenerator<-function (gen = "hg19",
                                       chr, 
                                       from, 
                                       to,
                                       bedFile, 
                                       featureDisplay = featureDisplay, 
                                       colorcase = "roadmap15") 
     {
       desiredRegion <- subset(bedFile, end > from & 
                                 start < to & seqnames == chr)
       
       track <- AnnotationTrack(desiredRegion, 
                                stacking = "dense",
                                col.line="black",
                                feature = (mcols(desiredRegion))$abbr,
                                genome = "hg19",
                                strand= "*",
                                name = "Cell Type Selected")
       
       displayPars(track) <- list(`1_TssA` = "#FF0000", `2_TssAFlnk` = "#FF6E00", 
                                  `3_TxFlnk` = "#32CD32", `4_Tx` = "#008000", `5_TxWk` = "#006400", 
                                  `6_EnhG` = "#C2E105", `7_Enh` = "#FFFF00", `8_ZNF/Rpts` = "#66CDAA", 
                                  `9_Het` = "#8A91D0", `10_TssBiv` = "#CD5C5C", `11_BivFlnk` = "#E9967A", 
                                  `12_EnhBiv` = "#BDB76B", `13_ReprPC` = "#3A3838", 
                                  `14_ReprPCWk` = "#808080", `15_Quies` = "#DCDCDC", 
                                  Empty = "#ffffff")
       return(track)
     }
     
     assign("chromHMMTrackGenerator", chromHMMTrackGenerator, .GlobalEnv)
     
     ###############################################
     ###Interaction Track Based on the CAGE Expression Enhancer-promoter assoication
     #################################################
     
     CageExpressionGenomicIntearctions <- readRDS(file = "../DataFiles/Interactions/Human/EnhancerPromoterAssoicationRObject")
     
     IntearctionTrack<-CageExpressionGenomicIntearctions%>%InteractionTrack()
     
     assign("IntearctionTrack", IntearctionTrack, .GlobalEnv)
     
     
     
     ### Pipe line for identifying potential TFBS
     
     ### Enhancer sites
     #To do this, we took the CAGE TSS-enhancer assoication off of Fantom5 Database to identify what genes are 
     #Regulated by these motifs in enhancer regions
     
     
     
     #####################################################################
     ###Start of the pipe Line
     ###################################################################
     
     
     ### Input Transcription Factor Matrix from the seqlogo function
     
     assign("TranscriptionFactorPWM",matrixForMatching, .GlobalEnv )
     
     
       genomicLocationOfMotifs<-matchPWM(matrixForMatching, 
                                         BSgenome.Hsapiens.UCSC.hg19,
                                         input$MatchPercentage)
       
       assign("genomicLocationOfMotifs", genomicLocationOfMotifs, .GlobalEnv)
     
    
     
     ### Identify which of these motifs are located in enhancer regions
     MotifsInEnhancers<-subsetByOverlaps(genomicLocationOfMotifs ,
                                         EnhancersWithGeneTargetsGrange)
     
     if(input$Conserved==TRUE){
       
       ConservedValue<-input$Conserved
       
       assign("ConservedValue", ConservedValue , .GlobalEnv)
       
       if(!exists("ConservedRegionsInPromoters")){
         
         # Import the Conserved region and run the following code
         
         ConservedRegionsInPromoters<-readRDS("../DataFiles/Conserved Region/Human/PromoterConservedRegions")
         assign("ConservedRegionsInPromoters", ConservedRegionsInPromoters, .GlobalEnv)
         
         MotifsInConservedPromoterRegions<-subset(genomicLocationOfMotifs, 
                                                  countOverlaps(genomicLocationOfMotifs,
                                                                ConservedRegionsInPromoters)>=genomicLocationOfMotifs$string[[1]]%>%length)
       } else{
         # If it does exist it will run this and and overlap it to get motifs that are fully conserved
         
         MotifsInConservedPromoterRegions<-subset(genomicLocationOfMotifs, 
                                                  countOverlaps(genomicLocationOfMotifs,
                                                                ConservedRegionsInPromoters)>=genomicLocationOfMotifs$string[[1]]%>%length)
       }
       
     } else{
       ConservedValue<-input$Conserved
       
       assign("ConservedValue", ConservedValue , .GlobalEnv)
       
       # If we are not looking for conserved promoter regions
       MotifsInConservedPromoterRegions<-subsetByOverlaps(genomicLocationOfMotifs, promoterTracks)
       
     }
     
     #### Combing the promoter and enhancer motifs into a single Grange List object
     
     MotifsInPromotersAndEnhancers<-list("Promoters"=MotifsInConservedPromoterRegions,
                                         "Enhancers"=MotifsInEnhancers)
     
     assign("MotifsInPromotersAndEnhancers", MotifsInPromotersAndEnhancers, .GlobalEnv)
     
     
     
     
     
     ################################################
     #### Cell Type specific epigneonic Analysis
     ##################################################

     CellTypeToPredict<-input$CellTypeToPredict
     
     assign("CellTypeToPredict", CellTypeToPredict, .GlobalEnv )

     ### Download the Chromatin state file From the annotation hub

     epiFiles <- query(ah, c(paste0(CellTypeToPredict,"_15_coreMarks_mnemonics"), "EpigenomeRoadMap") )
     
     chromatinState<-epiFiles[[paste0("AH", epiFiles@.db_uid)]]
     
     
     #Assign the bedfile to the global environment for analysis later on downstream
     
     assign("chromatinState", chromatinState, .GlobalEnv)
     
     ## Subsetting the Chromatin states for active states to identify motifs in these regions
     ActiveChromatinStates<-c("10_TssBiv",
                              "7_Enh",
                              "1_TssA",
                              "11_BivFlnk",
                              "2_TssAFlnk", 
                              "5_TxWk",
                              "4_Tx",
                              "8_ZNF/Rpts",
                              "6_EnhG", 
                              "12_EnhBiv",
                              "3_TxFlnk")
     
     
     ActiveChromatinRegions<-chromatinState[chromatinState$abbr %in% ActiveChromatinStates,]
     
     #Overlap active reginos with motifs in CRMs giving us a completely unbiased set of results
     UnbiasedPredictedMotifs<-lapply(MotifsInPromotersAndEnhancers,function(x){subsetByOverlaps( x,
                                                                                                 ActiveChromatinRegions)})
     
     assign("UnbiasedPredictedMotifs", UnbiasedPredictedMotifs, .GlobalEnv )

     ###################################################
     ###Assigning Gene targets to eac motif based on the assoication of the CRM it is located within. 
     ####################################################
     
     ### First we identify the genes that are regulated by our unbiased predicted results using genomic annotation 
     
     
     #Overlap the promoter regions of genes with unbiased motifs returning promoters with a predicted TFBS
     OverlappingRangeOfMOtifsInPromoters<-findOverlaps(promoterTracks, UnbiasedPredictedMotifs$Promoters)
     
     unbiasedPromoterMotifs<-UnbiasedPredictedMotifs$Promoters[OverlappingRangeOfMOtifsInPromoters%>%subjectHits()]
     
     
     mcols(unbiasedPromoterMotifs)<- cbind.data.frame(mcols(unbiasedPromoterMotifs),
                                                      "Genes Regulated" = promoterTracks[OverlappingRangeOfMOtifsInPromoters%>%queryHits()]$hg19.kgXref.geneSymbol,
                                                      "UCSC Transcript ID" = promoterTracks[OverlappingRangeOfMOtifsInPromoters%>%queryHits()]$hg19.kgXref.kgID)
   
     
     ## Now lets get the promoters of genes regulated by motifs in enhancres
     OverlappingRangeEnhancersMotifs<-findOverlaps(EnhancersWithGeneTargetsGrange, UnbiasedPredictedMotifs$Enhancers)
     
     EnhancerTargets<-cbind.data.frame(EnhancersWithGeneTargetsGrange$hg19.kgXref.geneSymbol,EnhancersWithGeneTargetsGrange$X.hg19.knownGene.name)[OverlappingRangeEnhancersMotifs%>%queryHits(),]
     
     colnames(EnhancerTargets)<- c("Genes Regulated" , "UCSC Transcript ID")
     
     MotifsInEnhancers<- UnbiasedPredictedMotifs$Enhancers[OverlappingRangeEnhancersMotifs%>%subjectHits()]
     
     mcols(MotifsInEnhancers)<- cbind.data.frame(mcols(MotifsInEnhancers), EnhancerTargets)
     
     
     UnbiasedMotifsPredicted<-c("Promoters With Gene Targets" = unbiasedPromoterMotifs,
                                "Enhancers With Gene Targets" = MotifsInEnhancers)
     
     assign("UnbiasedMotifsPredicted", UnbiasedMotifsPredicted, .GlobalEnv)
     
     
  
     
     
     #######################################################
     ## Left Joining with differentially expressed gene list
     #########################################################
     
     if(input$DifferentialExpressedGenes==TRUE){
       
       
       differenitallyExpressedGenesList<-input$differenitallyExpressedGenesList
       
       assign("differenitallyExpressedGenesList", differenitallyExpressedGenesList, .GlobalEnv)
       
       ## Genes who showed differenital expression with an enhancer that was correlated with its expression
       enhancerTargetsOfTF<-subset(UnbiasedMotifsPredicted$`Enhancers With Gene Targets`,  
                                   `Genes Regulated` %in% differenitallyExpressedGenesList) 
       
       
       ##Genes who showed differenital expression with promoter targets
       promoterTargetsOfTF<-subset(UnbiasedMotifsPredicted$`Promoters With Gene Targets` ,
                                   `Genes Regulated` %in% differenitallyExpressedGenesList) 
       
       
       
       ##Predicted Sites in the Regulatory Elements of these genes
       
       
       returnObjectDifferentialSites<-c("Promoter Predicted Sites" = promoterTargetsOfTF,
                       "Enhancer Predicted Sites" = enhancerTargetsOfTF)%>%unlist()
                      
       GenomeBrowserBiasedSites<-c(returnObjectDifferentialSites$`Promoter Predicted Sites`,
                                   returnObjectDifferentialSites$`Enhancer Predicted Sites`)%>%unlist()
       
       assign("PredictedTFBS", returnObjectDifferentialSites, .GlobalEnv)
       assign("returnObjectDifferentialSites", returnObjectDifferentialSites, .GlobalEnv)
     
     }
     else {
       
       ## Returning unbiased results without transcriptomic data
       
       
       returnObjectUnbaised<-c(
         "Promoter Predicted Sites"= UnbiasedMotifsPredicted$`Enhancers With Gene Targets`,
         "Enhancer Predicted Sites"= UnbiasedMotifsPredicted$`Promoters With Gene Targets` )%>%unlist()
       
       assign("returnObjectUnbaised", returnObjectUnbaised, .GlobalEnv)
        
       GenomeBrowserUnbiasedSites<-c(returnObjectUnbaised$`Promoter Predicted Sites`,
                                     returnObjectUnbaised$`Enhancer Predicted Sites`)%>%unlist()
       
       assign("PredictedTFBS", GenomeBrowserUnbiasedSites, .GlobalEnv)
     }
     
PromoterPredictedSites <- UnbiasedMotifsPredicted$`Promoters With Gene Targets`%>%as.data.frame()
  })
    } )
  

  
  output$EnhancerPredictedSites<-renderDataTable({
    EnhancerPredictedSites <- UnbiasedMotifsPredicted$`Enhancers With Gene Targets`%>%as.data.frame()
  })


















output$GenomeBrowser<-renderPlot({ 
  
  input$GenomeBrowserAction
  isolate({
  
  if(!exists("chrM")){
    ###############################################################################################
    ###Genome browser part for when you change chromosomes
    ##############################################################################################
    chrM<-input$chrM
    assign("chrM", chrM, .GlobalEnv)
    
    # Ideogram Track
    humanIdeogramTrack<-IdeogramTrack(chromosome = input$chrM, genome="hg19",name= "Ideogram")
    
    #Genome Axis Track for apprxoimate location
    gHumanTrack<-GenomeAxisTrack(name= "Axis")
    
    assign("humanIdeogramTrack", humanIdeogramTrack, .GlobalEnv)
    assign("gHumanTrack", gHumanTrack, .GlobalEnv)
    
    
    ###############################################
    ####Identifying motifs in CRM regions
    ################################################
    
    
    
    # Gene Track with symbols :D
    knownGenes <- GeneRegionTrack(TxDb.Hsapiens.UCSC.hg19.knownGene, 
                                  genome="hg19", 
                                  chromosome=input$chrM, 
                                  showId=TRUE,
                                  geneSymbol=TRUE, 
                                  name="UCSC")
    
    symbols <- unlist(mapIds(org.Hs.eg.db, gene(knownGenes),
                             "SYMBOL", "ENTREZID", 
                             multiVals = "first"))
    
    symbol(knownGenes) <- symbols[gene(knownGenes)]
    
    assign("knownGenes", knownGenes, .GlobalEnv)
    
    
    
    #Promoter and Enhancer Tracks for each chormosome Track
    promotertrackChromosomeSpecific <- promoterTracks%>%subset(. , 
                                                               seqnames==input$chrM)%>%AnnotationTrack(., name= "PromoterTrack", 
                                                                                                       genome="hg19")
    geneTrackChromosomeSpecific <- knownGenes
    
    
    
    EnhancersHumanChromosomeSpecific <- EnhancersWithGeneTargetsGrange%>%subset(. ,
                                                                seqnames==input$chrM)%>%AnnotationTrack(., name = "Enhancers",
                                                                                                        genome = "hg19")
    
    assign("EnhancersHumanChromosomeSpecific", EnhancersHumanChromosomeSpecific, .GlobalEnv)
    assign("promotertrackChromosomeSpecific", promotertrackChromosomeSpecific, .GlobalEnv)
    assign("geneTrackChromosomeSpecific", geneTrackChromosomeSpecific, .GlobalEnv)
    
    #Chromosome Specific Predicted Motifs
    PredictedTFBSTrack<-PredictedTFBS%>%subset(seqnames==input$chrM)%>%AnnotationTrack(genome = "hg19", 
                                                                                       stacking = "dense",
                                                                                       strand= "*",
                                                                                       col.line="black",
                                                                                       name="Predicted TFBS")
    
    
    
    assign("PredictedTFBSTrack", PredictedTFBSTrack, .GlobalEnv)
    ########################################################
    ## Re render each time anything changes
    ###########################################################
    # Raw Motif Instances
    RawMotifInstancesTrack<-subset(genomicLocationOfMotifs, 
                                   seqnames==input$chrM & start > input$fromM & end< input$toM)%>%AnnotationTrack(.,
                                                                                                                  genome = "hg19",
                                                                                                                  stacking = "dense", 
                                                                                                                  col.line="black",
                                                                                                                  name="All Motif Instances")
    assign("RawMotifInstancesTrack", RawMotifInstancesTrack, .GlobalEnv)
    
    
    # Chromosome For Predicted Motifs
    chromatinStatesTrack<-chromHMMTrackGenerator(gen="hg19", 
                                                 chr= input$chrM, 
                                                 from  = input$fromM,
                                                 to = input$toM,
                                                 bedFile = chromatinState,
                                                 featureDisplay = "all",
                                                 colorcase='roadmap15')
    
    assign("chromatinStatesTrack", chromatinStatesTrack, .GlobalEnv)
    
    
    
    plotTracks(trackList =c(humanIdeogramTrack,
                            gHumanTrack, 
                            IntearctionTrack,
                            EnhancersHumanChromosomeSpecific,
                            PredictedTFBSTrack,
                            RawMotifInstancesTrack, 
                            promotertrackChromosomeSpecific, 
                            geneTrackChromosomeSpecific,
                            chromatinStatesTrack), 
               sizes= c(1,1,3,1,1,1,1,3,3),
               from =input$fromM, 
               to= input$toM,
               chromosome= input$chrM,
               cex.title = 0.72, 
               rotation.title = 0, 
               showAxis = FALSE, 
               background.title = "white",
               lwd.title = 2, 
               title.width = 2, 
               cex.main = 5, 
               col = NULL, 
               fontcolor.title = "black")
    
   
  } else if(!chrM==input$chrM){ 
        
        ###############################################################################################
        ###Genome browser part for when you change chromosomes
        ##############################################################################################
        chrM<-input$chrM
        assign("chrM", chrM, .GlobalEnv)
        
        humanIdeogramTrack<-IdeogramTrack(chromosome = input$chrM, genome="hg19",name= "Ideogram")
        gHumanTrack<-GenomeAxisTrack(name= "Axis")
        
        assign("humanIdeogramTrack", humanIdeogramTrack, .GlobalEnv)
        assign("gHumanTrack", gHumanTrack, .GlobalEnv)
        
        
        ###############################################
        ####Identifying motifs in CRM regions
        ################################################
        
        
        
        # Gene Track with symbols :D
        knownGenes <- GeneRegionTrack(TxDb.Hsapiens.UCSC.hg19.knownGene, 
                                      genome="hg19", 
                                      chromosome=input$chrM, 
                                      showId=TRUE,
                                      geneSymbol=TRUE, 
                                      name="UCSC")
        
        symbols <- unlist(mapIds(org.Hs.eg.db, gene(knownGenes),
                                 "SYMBOL", "ENTREZID", 
                                 multiVals = "first"))
        
        symbol(knownGenes) <- symbols[gene(knownGenes)]
        
        assign("knownGenes", knownGenes, .GlobalEnv)
        
        
        
        #Promoter and Enhancer Tracks for each chormosome Track
        promotertrackChromosomeSpecific <- promoterTracks%>%subset(. , 
                                                                   seqnames==input$chrM)%>%AnnotationTrack(., name= "PromoterTrack", 
                                                                                                           genome="hg19")
        geneTrackChromosomeSpecific <- knownGenes
        
        
        
        EnhancersHumanChromosomeSpecific <- EnhancersWithGeneTargetsGrange%>%subset(. ,
                                                                    seqnames==input$chrM)%>%AnnotationTrack(., name = "Enhancers",
                                                                                                            genome = "hg19")
        
        assign("EnhancersHumanChromosomeSpecific", EnhancersHumanChromosomeSpecific, .GlobalEnv)
        assign("promotertrackChromosomeSpecific", promotertrackChromosomeSpecific, .GlobalEnv)
        assign("geneTrackChromosomeSpecific", geneTrackChromosomeSpecific, .GlobalEnv)
        
        #Chromosome Specific Predicted Motifs
        PredictedTFBSTrack<-PredictedTFBS%>%subset(seqnames==input$chrM)%>%AnnotationTrack(genome = "hg19", 
                                                                                           stacking = "dense",
                                                                                           strand= "*",
                                                                                           col.line="black",
                                                                                           name="Predicted TFBS")
        
        
        
        assign("PredictedTFBSTrack", PredictedTFBSTrack, .GlobalEnv)
        ########################################################
        ## Re render each time anything changes
        ###########################################################
        # Raw Motif Instances
        RawMotifInstancesTrack<-subset(genomicLocationOfMotifs, 
                                       seqnames==input$chrM & start > input$fromM & end< input$toM)%>%AnnotationTrack(.,
                                                                                                                      genome = "hg19",
                                                                                                                      stacking = "dense", 
                                                                                                                      col.line="black",
                                                                                                                      name="All Motif Instances")
        assign("RawMotifInstancesTrack", RawMotifInstancesTrack, .GlobalEnv)
        
        
        # Chromosome For Predicted Motifs
        
        # Chromosome For Predicted Motifs
        chromatinStatesTrack<-chromHMMTrackGenerator(gen="hg19", 
                                                     chr= input$chrM, 
                                                     from  = input$fromM,
                                                     to = input$toM,
                                                     bedFile = chromatinState,
                                                     featureDisplay = "all",
                                                     colorcase='roadmap15')
        
        
        assign("chromatinStatesTrack", chromatinStatesTrack, .GlobalEnv)
        
        
        
        plotTracks(trackList =c(humanIdeogramTrack,
                                gHumanTrack, 
                                IntearctionTrack,
                                EnhancersHumanChromosomeSpecific,
                                PredictedTFBSTrack,
                                RawMotifInstancesTrack, 
                                promotertrackChromosomeSpecific, 
                                geneTrackChromosomeSpecific,
                                chromatinStatesTrack), 
                   sizes= c(1,1,3,1,1,1,1,3,3),
                   from =input$fromM, 
                   to= input$toM,
                   chromosome= input$chrM,
                   cex.title = 0.72, 
                   rotation.title = 0, 
                   showAxis = FALSE, 
                   background.title = "white",
                   lwd.title = 2, 
                   title.width = 2, 
                   cex.main = 5, 
                   col = NULL, 
                   fontcolor.title = "black")
        
        } else {
          
          
          ########################################################
          ## Re render each time anything changes
          ###########################################################
          # Raw Motif Instances
          RawMotifInstancesTrack<-subset(genomicLocationOfMotifs, 
                                         seqnames==input$chrM & start > input$fromM & end< input$toM)%>%AnnotationTrack(.,
                                                                                                                        genome = "hg19",
                                                                                                                        stacking = "dense", 
                                                                                                                        col.line="black",
                                                                                                                        name="All Motif Instances")
          assign("RawMotifInstancesTrack", RawMotifInstancesTrack, .GlobalEnv)
          
          
          # Chromosome For Predicted Motifs
          
          # Chromosome For Predicted Motifs
          chromatinStatesTrack<-chromHMMTrackGenerator(gen="hg19", 
                                                       chr= input$chrM, 
                                                       from  = input$fromM,
                                                       to = input$toM,
                                                       bedFile = chromatinState,
                                                       featureDisplay = "all",
                                                       colorcase='roadmap15')
          
          assign("chromatinStatesTrack", chromatinStatesTrack, .GlobalEnv)
          
          
          
          plotTracks(trackList =c(humanIdeogramTrack,
                                  gHumanTrack, 
                                  IntearctionTrack,
                                  EnhancersHumanChromosomeSpecific,
                                  PredictedTFBSTrack,
                                  RawMotifInstancesTrack, 
                                  promotertrackChromosomeSpecific, 
                                  geneTrackChromosomeSpecific,
                                  chromatinStatesTrack), 
                     sizes= c(1,1,3,1,1,1,1,3,3),
                     from =input$fromM, 
                     to= input$toM,
                     chromosome= input$chrM,
                     cex.title = 0.72, 
                     rotation.title = 0, 
                     showAxis = FALSE, 
                     background.title = "white",
                     lwd.title = 2, 
                     title.width = 2, 
                     cex.main = 5, 
                     col = NULL, 
                     fontcolor.title = "black")
        }

  })
})
  
 output$DataTablePredictedSites<- renderDataTable({

   RenderDataFrame <- subset(PredictedTFBS, start>=input$fromM & end<= input$toM & seqnames== input$chrM)%>%as.data.frame()
   }
  )

  
 
 
  
  
  
  
  output$LegendsPlot<- renderImage({
    
    list(
      src = "www/EpigenomicsRoadMapLegendHMM.jpeg",
      contentType = "image/jpeg",
      alt = "Human/Epigenomics Road Map Legend"
    )
    
  }, deleteFile = FALSE)
  
  
  
  
  output$GeneOntologyResults <- renderDataTable({
    
    data <- unique(c(as.character(UnbiasedMotifsPredicted$`Promoters With Gene Targets`$`UCSC Transcript ID`), 
                     as.character(UnbiasedMotifsPredicted$`Enhancers With Gene Targets`$`UCSC Transcript ID`)))
    
    
    human_mart <- biomaRt::useMart("ENSEMBL_MART_ENSEMBL", dataset = "hsapiens_gene_ensembl")
    
    ucscToEntrez <- getBM(attributes = c("ucsc", "entrezgene", "external_gene_name"),
                          filters = "ucsc",
                          values = list(data), mart = human_mart)
    
    uniservser <- getBM(attributes = c("ucsc", "entrezgene", "external_gene_name"), mart = human_mart)
    
    geneOntology<-goana(de = ucscToEntrez$entrezgene, universe = uniservser$entrezgene, FDR=0.05, species = "Hs")
    
    
    GeneOntologyResultsSOrted <- topGO(geneOntology)

    
  })
  

  
  
  
}
)
