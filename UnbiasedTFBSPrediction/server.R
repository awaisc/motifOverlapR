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



shinyServer(function(input, output) {


  #######################################################
  ### Identifying which are located in enhancer regions and which are located in promoter regions
  ### And Identifying where they are located in human genome
  #######################################################
  ##############################################
  #Render Human GVIZ plot 1
  #############################################33
  output$HumangvizPlot <- renderPlot({
   if(!exists("chrM")){
     
     ########################################
     ###Data Inputs
     #######################################
     #Enhancers
     EnhancersHuman<-import("DataFiles/Enhancers Track/Human/human_permissive_enhancers_phase_1_and_2.bed.gz")
     
     assign("EnhancersHuman", EnhancersHuman, .GlobalEnv)
     
     ConservedRegionsInPromoters<-readRDS("DataFiles/Conserved Region/Human/PromoterConservedRegions")
     assign("ConservedRegionsInPromoters", ConservedRegionsInPromoters, .GlobalEnv)
     
     ##Importing Bed File with Granges
     promoterTracks <- read.delim("/media/awais/NewDrivewho/Downloads/hg19bedWithNames.bed.gz")%>%
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
     
     
     # Epigenomic Chromatin state 
     
     ah<-AnnotationHub()
     ###################################
     ###Promoter Enhancer assoication data table Input
     ###################################
         #Read in the table
     EnhancerPromoterAssoications <-read_delim("/media/awais/NewDrivewho/UnbiasedTFBSPrediction/Unbiased TFBS Prediction/DataFiles/Enhancers Track/Human/hg19_enhancer_promoter_correlations_distances_cell_type.txt", 
                                               "\t", escape_double = FALSE, trim_ws = TRUE)
     assign("EnhancerPromoterAssoications", EnhancerPromoterAssoications, .GlobalEnv)
     
      
     ### Pipe line for identifying potential TFBS
     
     ### Enhancer sites
     #To do this, we took the CAGE TSS-enhancer assoication off of Fantom5 Database to identify what genes are 
     #Regulated by these motifs in enhancer regions
     

     ### Extract PWM and give them a name from the JASPAR database
     PFMatrixList <- getMatrixSet(JASPAR2016, list() )
     namedJasparDataBase<-lapply(1:length(PFMatrixList), function(x){
       names(PFMatrixList[x])<-PFMatrixList[[x]]@name  
     })
     
     names(PFMatrixList)<-namedJasparDataBase
     
     
     
     ### Select a PWM Matrix and match to the genome
     
     TranscriptionFactorPWM<-input$TranscriptionFactorPWM
     
     assign("TranscriptionFactorPWM", TranscriptionFactorPWM, .GlobalEnv)
     
     
     matrix<-PFMatrixList[[paste0(TranscriptionFactorPWM)]]@profileMatrix
     
     genomicLocationOfMotifs<-matchPWM(matrix, 
                                BSgenome.Hsapiens.UCSC.hg19,
                                input$MatchPercentage)
     
     assign("genomicLocationOfMotifs", genomicLocationOfMotifs, .GlobalEnv)
     ### Identify which of these motifs are located in enhancer regions
     MotifsInEnhancers<-subsetByOverlaps(genomicLocationOfMotifs ,
                                         EnhancersHuman)
     
     if(input$Conserved==TRUE){
     MotifsInConservedPromoterRegions<-subset(genomicLocationOfMotifs, 
                                              countOverlaps(genomicLocationOfMotifs,
                                                            ConservedRegionsInPromoters)>=genomicLocationOfMotifs$string[[1]]%>%length)
     } else(
       MotifsInConservedPromoterRegions<-subsetByOverlaps(genomicLocationOfMotifs, promoterTracks)
       
     )
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
                                                      "Genes Regulated" = promoterTracks[OverlappingRangeOfMOtifsInPromoters%>%queryHits()]$hg19.kgXref.geneSymbol)
   
     
     ## Now lets get the promoters of genes regulated by motifs in enhancres
     
     
     PromotersAssoicatedWithEnhancers<-separate(EnhancerPromoterAssoications, col=promoter, into= c("promoter", "strand"), sep= ',')%>%
     separate(., col= promoter, into= c("chr", "start", "end"))%>%makeGRangesFromDataFrame(keep.extra.columns = FALSE)
     
     
     
     
     
     EnhancerGrangeWithTargets<-EnhancerPromoterAssoications$enhancer%>%GRanges()
     
     mcols(EnhancerGrangeWithTargets)<-cbind.data.frame("Promoters"=PromotersAssoicatedWithEnhancers , EnhancerPromoterAssoications[3:6]) 
     
     
     
     EnhancerGrangeWithTargets <- subset(EnhancerGrangeWithTargets, `FDR`!=0 & `FDR`<0.05)
     
     #Split the colulmn with the enhancer information and make it into a GRange object that we overlap with our unbaised enhancer results
     

     assign("PredictedEnhancerTargets", PredictedEnhancerTargets, .GlobalEnv)
     
     #Gene the gene list from this GRANGE
     GenesWithAMotifInTheEnhancer <- PredictedEnhancerTargets[!duplicated(PredictedEnhancerTargets$gene)]
     
     
     #Combine them so we can overlap with Differenitally expressed genes
     DirectTargets<-list("Promoter Direct Targets" =GenesWithAMotifInThePromoter, 
                         "Enhancer Direct Targets"= GenesWithAMotifInTheEnhancer)
     
     
     assign("DirectTargets", DirectTargets, .GlobalEnv)
     
     
     #######################################################
     ## Left Joining with differentially expressed gene list
     #########################################################
     
     if(input$DifferentialExpressedGenes==TRUE){
       
       
       differenitallyExpressedGenesList<-input$differenitallyExpressedGenesList
       
       assign("differenitallyExpressedGenesList", differenitallyExpressedGenesList, .GlobalEnv)
       
       ## Genes who showed differenital expression with enhancer targets
       enhancerTargetsOfTF<-subset(DirectTargets$`Enhancer Direct Targets`, 
                                   gene %in% differenitallyExpressedGenesList) 
       
       
       ##Genes who showed differenital expression with promoter targets
       promoterTargetsOfTF<-subset(DirectTargets$`Promoter Direct Targets`,
                                   hg19.kgXref.geneSymbol %in% differenitallyExpressedGenesList) 
       
       
       
       ##Predicted Sites in the Regulatory Elements of these genes
       
       
       returnObjectDifferentialSites<-c("Promoter Predicted Sites" = subsetByOverlaps(UnbiasedPredictedMotifs$Promoters, promoterTargetsOfTF),
                       "Enhancer Predicted Sites" = subsetByOverlaps(UnbiasedPredictedMotifs$Enhancers, enhancerTargetsOfTF))%>%unlist()
                      
       GenomeBrowserBiasedSites<-c(returnObjectDifferentialSites$`Promoter Predicted Sites`,
                                   returnObjectDifferentialSites$`Enhancer Predicted Sites`)%>%unlist()
       
       assign("PredictedTFBS", returnObjectDifferentialSites, .GlobalEnv)
       assign("returnObjectDifferentialSites", returnObjectDifferentialSites, .GlobalEnv)
     
     }
     else {
       
       ## Returning unbiased results without transcriptomic data
       
       
       returnObjectUnbaised<-c(
         "Promoter Predicted Sites"= UnbiasedPredictedMotifs$Promoters,
         "Enhancer Predicted Sites"= UnbiasedPredictedMotifs$Enhancers)%>%unlist()
       
       assign("returnObjectUnbaised", returnObjectUnbaised, .GlobalEnv)
        
       GenomeBrowserUnbiasedSites<-c(returnObjectUnbaised$`Promoter Predicted Sites`,
                                     returnObjectUnbaised$`Enhancer Predicted Sites`)%>%unlist()
       
       assign("PredictedTFBS", GenomeBrowserUnbiasedSites, .GlobalEnv)
     }
   
     
     ###############################################################################################
     ###Genome browser part   
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
      
      
    
      PredictedTFBSTrack<-PredictedTFBS%>%AnnotationTrack(genome = "hg19", 
                                                          stacking = "dense",
                                                          strand= "*",
                                                          col.line="black",
                                                          name="Predicted TFBS")
     
      
      
      assign("PredictedTFBSTrack", PredictedTFBSTrack, .GlobalEnv)
    
      
      
      
      ##ChromHMM Track Generator specifically for humans
      chromHMMTrackGenerator<-function (gen = "hg19",
                                        chr, 
                                        from, 
                                        to,
                                        bedFile, 
                                        featureDisplay = featureDisplay, 
                                        colorcase = "roadmap15") 
      {
        desiredRegion <- subset(get(bedFile), end > from & 
                                  start < to & seqnames == chr)
        track <- AnnotationTrack(desiredRegion, 
                                 stacking = "dense",
                                 col.line="black",
                                 feature = (mcols(desiredRegion))$name,
                                 genome = "hg19",
                                 strand= "*",
                                 name = paste(bedFile))
        
        displayPars(track) <- list(`1_TssA` = "#FF0000", `2_TssAFlnk` = "#FF6E00", 
                                   `3_TxFlnk` = "#32CD32", `4_Tx` = "#008000", `5_TxWk` = "#006400", 
                                   `6_EnhG` = "#C2E105", `7_Enh` = "#FFFF00", `8_ZNF/Rpts` = "#66CDAA", 
                                   `9_Het` = "#8A91D0", `10_TssBiv` = "#CD5C5C", `11_BivFlnk` = "#E9967A", 
                                   `12_EnhBiv` = "#BDB76B", `13_ReprPC` = "#3A3838", 
                                   `14_ReprPCWk` = "#808080", `15_Quies` = "#DCDCDC", 
                                   Empty = "#ffffff")
        return(track)
      }
      
      chromatinStatesTrack<-chromHMMTrackGenerator(gen="hg19", 
                                                  chr= input$chrM, 
                                                  from  = input$fromM,
                                                  to = input$toM,
                                                  bedFile = chromatinState,
                                                  featureDisplay = "all",
                                                  colorcase='roadmap15')
      
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
      
      #Promoter and Enhancer Tracks for each chormosome Track
      promotertrackChromosomeSpecific<-promoterTracks%>%subset(. , 
                                                               seqnames==input$chrM)%>%AnnotationTrack(., name= "PromoterTrack", 
                                                                                                       genome="hg19")
      geneTrackChromosomeSpecific<-knownGenes
      EnhancersHumanChromosomeSpecific<-EnhancersHuman%>%subset(. ,
                                                                seqnames==input$chrM)%>%AnnotationTrack(., name = "Enhancers",
                                                                                                        genome = "hg19")
      
      # Raw Motif Instances
      RawMotifInstancesTrack<-subset(genomicPositions, 
                                seqnames==input$chrM & start > input$fromM & end< input$toM)%>%AnnotationTrack(.,genome = "hg19",
                                                                                                               stacking = "dense", 
                                                                                                               col.line="black",
                                                                                                               name="All Motifs")
      

      
      assign("chromatinStatesTrack", chromatinStatesTrack, .GlobalEnv)
      assign("EnhancersHumanChromosomeSpecific", EnhancersHumanChromosomeSpecific, .GlobalEnv)
      assign("promotertrackChromosomeSpecific", promotertrackChromosomeSpecific, .GlobalEnv)
      assign("geneTrackChromosomeSpecific", geneTrackChromosomeSpecific, .GlobalEnv)
      assign("knownGenes", knownGenes, .GlobalEnv)
      assign("chromHMMTrackGenerator", chromHMMTrackGenerator, .GlobalEnv)
      
      
      
      
        plotTracks(trackList =c(humanIdeogramTrack,
                                gHumanTrack, 
                                EnhancersHumanChromosomeSpecific,
                                PredictedTFBSTrack,
                                RawMotifInstancesTrack, 
                                promotertrackChromosomeSpecific, 
                                geneTrackChromosomeSpecific,
                                chromatinStatesTrack), 
                   sizes= c(1,1,1,1,1,1,3),
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
        
      } else if(!TranscriptionFactorPWM==input$TranscriptionFactorPWM){
        
        TranscriptionFactorPWM<-input$TranscriptionFactorPWM
        
        assign("TranscriptionFactorPWM", TranscriptionFactorPWM, .GlobalEnv)
        
        matrix<-PFMatrixList[[paste0(TranscriptionFactorPWM)]]@profileMatrix
        
        genomicPositions<-matchPWM(matrix, 
                                   BSgenome.Hsapiens.UCSC.hg19,
                                   input$MatchPercentage)
        
        assign("genomicPositions", genomicPositions, .GlobalEnv)
        ### Identify which of these motifs are located in enhancer regions
        MotifsInEnhancers<-subsetByOverlaps(genomicPositions ,
                                            EnhancersHuman)
        
        
        MotifsInConservedPromoterRegions<-subset(genomicPositions, 
                                                 countOverlaps(genomicPositions,
                                                               ConservedRegionsInPromoters)>=genomicPositions$string[[1]]%>%length)
        
        #### Combing the promoter and enhancer motifs into a single Grange List object
        
        MotifsInPromotersAndEnhancers<-list("Promoters"=MotifsInConservedPromoterRegions,
                                            "Enhancers"=MotifsInEnhancers)
        
        assign("MotifsInPromotersAndEnhancers", MotifsInPromotersAndEnhancers, .GlobalEnv)
        
        #### Reading in the data required for proper naming of these chromatin states bed files and deduplicating it
        
        # TableS1 <- read_excel("~/Downloads/nature14248-s2/TableS1.xlsx", 
        #                       sheet = "AdditionalQCScores")
        # # Deduplicating the list based on the epigenomic id
        # UniqueTable<-TableS1[!duplicated(TableS1$EID),]
        # saveRDS(UniqueTable,"DataFiles/ChromHMM/human/EpigneomicEnvironmentToCellTypeConverter")
        
        CellTypeToPredict<-input$CellTypeToPredict
        
        chromatinState<-import(paste0( "~/Downloads/all.dense.browserFiles/", 
                                       dir("~/Downloads/all.dense.browserFiles/", 
                                           pattern = UniqueTable$EID[UniqueTable$`Standardised epigenome name`== CellTypeToPredict]
                                       )
        )
        )
        
        ## Subsetting the Chromatin states for active states
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
        
        
        ActiveChromatinRegions<-chromatinState[chromatinState$name %in% ActiveChromatinStates,]
        
        #Overlap active reginos with motifs in CRMs giving us a completely unbiased set of results
        UnbiasedPredictedMotifs<-lapply(MotifsInPromotersAndEnhancers,function(x){subsetByOverlaps( x,
                                                                                                    ActiveChromatinRegions)})
        
        
        ############################################
        #### Writing the unbiased predicted Sites to a bedfile for exporting to the genome browser/bed file
        ############################################
        
        ##Promoter and enhancer sites unbiasedly predicted Sites
        AllAccessibleTFBSUnbiased<-rbind.data.frame(
          promoterDataFrameBedFile <- data.frame(seqnames=seqnames(UnbiasedPredictedMotifs$Promoters),
                                                 starts=start(UnbiasedPredictedMotifs$Promoters),
                                                 ends=end(UnbiasedPredictedMotifs$Promoters),
                                                 names=c(rep(".", length(UnbiasedPredictedMotifs$Promoters))),
                                                 scores=c(rep(".", length(UnbiasedPredictedMotifs$Promoters))),
                                                 strands=strand(UnbiasedPredictedMotifs$Promoters)),
          
          enhancerDataFrameBedFile <- data.frame(seqnames=seqnames(UnbiasedPredictedMotifs$Enhancers),
                                                 starts=start(UnbiasedPredictedMotifs$Enhancers),
                                                 ends=end(UnbiasedPredictedMotifs$Enhancers),
                                                 names=c(rep(".", length(UnbiasedPredictedMotifs$Enhancers))),
                                                 scores=c(rep(".", length(UnbiasedPredictedMotifs$Enhancers))),
                                                 strands=strand(UnbiasedPredictedMotifs$Enhancers)))%T>%
          write.table(.,
                      file="unbiasedPredictedTFBS.bed",
                      quote=FALSE, 
                      sep="\t", 
                      row.names=FALSE, 
                      col.names=FALSE)
        
        
        ######################################
        ###Left Join with differentialy expressed gene list
        #######################################
        
        ### First we identify the genes that are regulated by our unbiased predicted results
        
        
        #Overlap the promoter regions of genes with unbiased motifs returning promoters with a predicted TFBS
        promoterTranscriptionFactors<-subsetByOverlaps(promoterTracks,
                                                       UnbiasedPredictedMotifs$Promoters)
        GenesWithAMotifInThePromoter<-promoterTranscriptionFactors[!duplicated(promoterTranscriptionFactors$hg19.kgXref.geneSymbol)]
        
        
        
        
        ### Enhancer sites
        #To do this, we took the CAGE TSS-enhancer assoication off of Fantom5 Database to identify what genes are 
        #Regulated by these motifs in enhancer regions
        
        #Read in the table
        EnhancerPromoterAssoications <- read.delim("~/Downloads/human.associations.hdr.txt")
        
        
        #Split the colulmn with the enhancer information and make it into a GRange object that we overlap with our unbaised enhancer results
        
        PredictedTargets<-separate(EnhancerPromoterAssoications,
                                   col =promoter, 
                                   into = c("p", "gene"))%>%select(.,c("enhancer", "gene","cor", "fdr","distance"))%>%
          separate(., col=enhancer, into= c("chromosome", 
                                            "start",
                                            "end"))%>%
          makeGRangesFromDataFrame(.,
                                   keep.extra.columns=TRUE,
                                   ignore.strand=TRUE,
                                   seqinfo=NULL,
                                   seqnames.field="chromosome",
                                   start.field="start",
                                   end.field="end",
                                   starts.in.df.are.0based=FALSE)%>%subset(.,
                                                                           fdr<=0.05 & fdr!=0.000000e+00 )%>%subsetByOverlaps(. ,
                                                                                                                              UnbiasedPredictedMotifs$Enhancers)
        
        #Gene the gene list from this GRANGE
        GenesWithAMotifInTheEnhancer <- PredictedTargets[!duplicated(PredictedTargets$gene)]
        
        
        #Combine them so we can overlap with Differenitally expressed genes
        DirectTargets<-list("Promoter Direct Targets" =GenesWithAMotifInThePromoter, 
                            "Enhancer Direct Targets"= GenesWithAMotifInTheEnhancer)
        
        ## Genes with an promoter and enhacer motif
        # PromotersWithEnhancerSites<-subset(DirectTargets[[1]], 
        #                                    hg19.kgXref.geneSymbol %in% DirectTargets[[2]]$gene)
        # 
        # 
        
        
        #######################################################
        ## Left Joining with differentially expressed gene list
        #########################################################
        
        if(exists("differenitallyExpressedGenesList")){
          
          
          differenitallyExpressedGenesList<-input$differenitallyExpressedGenesList
          assign("differenitallyExpressedGenesList", differenitallyExpressedGenesList, .GlobalEnv)
          
          ## Genes who showed differenital expression with enhancer targets
          enhancerTargetsOfTF<-subset(DirectTargets[[2]], 
                                      gene %in% differenitallyExpressedGenesList) 
          
          
          ##Genes who showed differenital expression with promoter targets
          promoterTargetsOfTF<-subset(DirectTargets[[1]],
                                      hg19.kgXref.geneSymbol %in% differenitallyExpressedGenesList) 
          
          
          
          ##Predicted Sites in the Regulatory Elements of these genes
          
          
          returnObjectDifferentialSites<-c("Promoter Predicted Sites" = subsetByOverlaps(UnbiasedPredictedMotifs$Promoters, promoterTargetsOfTF),
                                           "Enhancer Predicted Sites" = subsetByOverlaps(UnbiasedPredictedMotifs$Enhancers, enhancerTargetsOfTF))%>%unlist()
          
          
          assign("PredictedTFBS", returnObjectDifferentialSites, .GlobalEnv)
          
        }
        else {
          
          ## Returning unbiased results without transcriptomic data
          
          
          returnObjectUnbaised<-c(
            "Promoter Predicted Sites"= subsetByOverlaps(UnbiasedPredictedMotifs$Promoters, promoterTargetsOfTF),
            "Promoter Predicted Sites"= subsetByOverlaps(UnbiasedPredictedMotifs$Enhancers, promoterTargetsOfTF)
          )%>%unlist()
          
          
          assign("PredictedTFBS", returnObjectUnbaised, .GlobalEnv)
        }
        
        
        ################################
        ###Genome browser part
        ###############################
        chrM<-input$chrM
        assign("chrM", chrM, .GlobalEnv)
        
        humanIdeogramTrack<-IdeogramTrack(chromosome = input$chrM, genome="hg19",name= "Ideogram")
        gHumanTrack<-GenomeAxisTrack(name= "Axis")
        
        assign("humanIdeogramTrack", humanIdeogramTrack, .GlobalEnv)
        assign("gHumanTrack", gHumanTrack, .GlobalEnv)
        
        
        ###############################################
        ####Identifying motifs in CRM regions
        ################################################
        
        
        
        PredictedTFBSTrack<-PredictedTFBS%>%AnnotationTrack(genome = "hg19", 
                                                            stacking = "dense",
                                                            strand= "*",
                                                            col.line="black",
                                                            name="Predicted TFBS")
        
        
        
        assign("PredictedTFBSTrack", PredictedTFBSTrack, .GlobalEnv)
        assign("EnhancersHuman", EnhancersHuman, .GlobalEnv)
        
        
        
        
        ##ChromHMM Track Generator specifically for humans
        chromHMMTrackGenerator<-function (gen = "hg19",
                                          chr, 
                                          from, 
                                          to,
                                          bedFile, 
                                          featureDisplay = featureDisplay, 
                                          colorcase = "roadmap15") 
        {
          desiredRegion <- subset(get(bedFile), end > from & 
                                    start < to & seqnames == chr)
          track <- AnnotationTrack(desiredRegion, 
                                   stacking = "dense",
                                   col.line="black",
                                   feature = (mcols(desiredRegion))$name,
                                   genome = "hg19",
                                   strand= "*",
                                   name = paste(bedFile))
          
          displayPars(track) <- list(`1_TssA` = "#FF0000", `2_TssAFlnk` = "#FF6E00", 
                                     `3_TxFlnk` = "#32CD32", `4_Tx` = "#008000", `5_TxWk` = "#006400", 
                                     `6_EnhG` = "#C2E105", `7_Enh` = "#FFFF00", `8_ZNF/Rpts` = "#66CDAA", 
                                     `9_Het` = "#8A91D0", `10_TssBiv` = "#CD5C5C", `11_BivFlnk` = "#E9967A", 
                                     `12_EnhBiv` = "#BDB76B", `13_ReprPC` = "#3A3838", 
                                     `14_ReprPCWk` = "#808080", `15_Quies` = "#DCDCDC", 
                                     Empty = "#ffffff")
          return(track)
        }
        
        chromatinStatesTrack<-chromHMMTrackGenerator(gen="hg19", 
                                                     chr= input$chrM, 
                                                     from  = input$fromM,
                                                     to = input$toM,
                                                     bedFile = chromatinState,
                                                     featureDisplay = "all",
                                                     colorcase='roadmap15')
        
        # Gene Track with symbols :D
        knownGenes <- GeneRegionTrack(TxDb.Hsapiens.UCSC.hg19.knownGene, 
                                      genome="hg19", 
                                      chromosome="chrX", 
                                      showId=TRUE,
                                      geneSymbol=TRUE, 
                                      name="UCSC")
        
        symbols <- unlist(mapIds(org.Hs.eg.db, gene(knownGenes),
                                 "SYMBOL", "ENTREZID", 
                                 multiVals = "first"))
        symbol(knownGenes) <- symbols[gene(knownGenes)]
        
        #Promoter and Enhancer Tracks for each chormosome Track
        promotertrackChromosomeSpecific<-promoterTracks%>%subset(. , 
                                                                 seqnames==input$chrM)%>%AnnotationTrack(., name= "PromoterTrack", 
                                                                                                         genome="hg19")
        geneTrackChromosomeSpecific<-knownGenes
        EnhancersHumanChromosomeSpecific<-EnhancersHuman%>%subset(. ,
                                                                  seqnames==input$chrM)%>%AnnotationTrack(., name = "Enhancers",
                                                                                                          genome = "hg19")
        
        # Raw Motif Instances
        RawMotifInstancesTrack<-subset(genomicPositions, 
                                       seqnames==input$chrM & start > input$fromM & end< input$toM)%>%AnnotationTrack(.,genome = "hg19",
                                                                                                                      stacking = "dense", 
                                                                                                                      col.line="black",
                                                                                                                      name="All Motifs")
        
        
        
        assign("chromatinStatesTrack", chromatinStatesTrack, .GlobalEnv)
        assign("RawMotifInstancesTrack", RawMotifInstancesTrack, .GlobalEnv)
        assign("EnhancersHumanChromosomeSpecific", EnhancersHumanChromosomeSpecific, .GlobalEnv)
        assign("promotertrackChromosomeSpecific", promotertrackChromosomeSpecific, .GlobalEnv)
        assign("geneTrackChromosomeSpecific", geneTrackChromosomeSpecific, .GlobalEnv)
        assign("knownGenes", knownGenes, .GlobalEnv)
        assign("chromHMMTrackGenerator", chromHMMTrackGenerator, .GlobalEnv)
        
        
        
        
        plotTracks(trackList =c(humanIdeogramTrack,
                                gHumanTrack, 
                                EnhancersHumanChromosomeSpecific,
                                PredictedTFBSTrack,
                                RawMotifInstancesTrack, 
                                promotertrackChromosomeSpecific, 
                                geneTrackChromosomeSpecific,
                                chromatinStatesTrack), 
                   sizes= c(1,1,1,1,1,1,3),
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
        
    } else if(!CellTypeToPredict==input$CellTypeToPredict){
        CellTypeToPredict<-input$CellTypeToPredict
      
        assign("CellTypeToPredict", CellTypeToPredict, .GlobalEnv)
        
        chromatinState<-import(paste0( "~/Downloads/all.dense.browserFiles/", 
                                       dir("~/Downloads/all.dense.browserFiles/", 
                                           pattern = UniqueTable$EID[UniqueTable$`Standardised epigenome name`== CellTypeToPredict]
                                       )
        )
        )
        
        ## Subsetting the Chromatin states for active states
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
        
        
        ActiveChromatinRegions<-chromatinState[chromatinState$name %in% ActiveChromatinStates,]
        
        #Overlap active reginos with motifs in CRMs giving us a completely unbiased set of results
        UnbiasedPredictedMotifs<-lapply(MotifsInPromotersAndEnhancers,function(x){subsetByOverlaps( x,
                                                                                                    ActiveChromatinRegions)})
        
        
        ############################################
        #### Writing the unbiased predicted Sites to a bedfile for exporting to the genome browser/bed file
        ############################################
        
        # ##Promoter and enhancer sites unbiasedly predicted Sites
        # AllAccessibleTFBSUnbiased<-rbind.data.frame(
        #   promoterDataFrameBedFile <- data.frame(seqnames=seqnames(UnbiasedPredictedMotifs$Promoters),
        #                                          starts=start(UnbiasedPredictedMotifs$Promoters),
        #                                          ends=end(UnbiasedPredictedMotifs$Promoters),
        #                                          names=c(rep(".", length(UnbiasedPredictedMotifs$Promoters))),
        #                                          scores=c(rep(".", length(UnbiasedPredictedMotifs$Promoters))),
        #                                          strands=strand(UnbiasedPredictedMotifs$Promoters)),
        #   
        #   enhancerDataFrameBedFile <- data.frame(seqnames=seqnames(UnbiasedPredictedMotifs$Enhancers),
        #                                          starts=start(UnbiasedPredictedMotifs$Enhancers),
        #                                          ends=end(UnbiasedPredictedMotifs$Enhancers),
        #                                          names=c(rep(".", length(UnbiasedPredictedMotifs$Enhancers))),
        #                                          scores=c(rep(".", length(UnbiasedPredictedMotifs$Enhancers))),
        #                                          strands=strand(UnbiasedPredictedMotifs$Enhancers)))%T>%
        #   write.table(.,
        #               file="unbiasedPredictedTFBS.bed",
        #               quote=FALSE, 
        #               sep="\t", 
        #               row.names=FALSE, 
        #               col.names=FALSE)
        
        
        ######################################
        ###Left Join with differentialy expressed gene list
        #######################################
        
        ### First we identify the genes that are regulated by our unbiased predicted results
        
        
        #Overlap the promoter regions of genes with unbiased motifs returning promoters with a predicted TFBS
        promoterTranscriptionFactors<-subsetByOverlaps(promoterTracks,
                                                       UnbiasedPredictedMotifs$Promoters)
        GenesWithAMotifInThePromoter<-promoterTranscriptionFactors[!duplicated(promoterTranscriptionFactors$hg19.kgXref.geneSymbol)]
        
        
        
        
        ### Enhancer sites
        #To do this, we took the CAGE TSS-enhancer assoication off of Fantom5 Database to identify what genes are 
        #Regulated by these motifs in enhancer regions
        
        #Read in the table
        EnhancerPromoterAssoications <- read.delim("~/Downloads/human.associations.hdr.txt")

        #Split the colulmn with the enhancer information and make it into a GRange object that we overlap with our unbaised enhancer results
        
        PredictedTargets<-separate(EnhancerPromoterAssoications,
                                   col =promoter, 
                                   into = c("p", "gene"))%>%select(.,c("enhancer", "gene","cor", "fdr","distance"))%>%
          separate(., col=enhancer, into= c("chromosome", 
                                            "start",
                                            "end"))%>%
          makeGRangesFromDataFrame(.,
                                   keep.extra.columns=TRUE,
                                   ignore.strand=TRUE,
                                   seqinfo=NULL,
                                   seqnames.field="chromosome",
                                   start.field="start",
                                   end.field="end",
                                   starts.in.df.are.0based=FALSE)%>%subset(.,
                                                                           fdr<=0.05 & fdr!=0.000000e+00 )%>%subsetByOverlaps(. ,
                                                                                                                              UnbiasedPredictedMotifs$Enhancers)
        
        #Gene the gene list from this GRANGE
        GenesWithAMotifInTheEnhancer <- PredictedTargets[!duplicated(PredictedTargets$gene)]
        
        
        #Combine them so we can overlap with Differenitally expressed genes
        DirectTargets<-list("Promoter Direct Targets" =GenesWithAMotifInThePromoter, 
                            "Enhancer Direct Targets"= GenesWithAMotifInTheEnhancer)
        
        ## Genes with an promoter and enhacer motif
        # PromotersWithEnhancerSites<-subset(DirectTargets[[1]], 
        #                                    hg19.kgXref.geneSymbol %in% DirectTargets[[2]]$gene)
        # 
        # 
        
        
        #######################################################
        ## Left Joining with differentially expressed gene list
        #########################################################
        
        if(exists("differenitallyExpressedGenesList")){
          
          
          differenitallyExpressedGenesList<-input$differenitallyExpressedGenesList
          assign("differenitallyExpressedGenesList", differenitallyExpressedGenesList, .GlobalEnv)
          
          ## Genes who showed differenital expression with enhancer targets
          enhancerTargetsOfTF<-subset(DirectTargets[[2]], 
                                      gene %in% differenitallyExpressedGenesList) 
          
          
          ##Genes who showed differenital expression with promoter targets
          promoterTargetsOfTF<-subset(DirectTargets[[1]],
                                      hg19.kgXref.geneSymbol %in% differenitallyExpressedGenesList) 
          
          
          
          ##Predicted Sites in the Regulatory Elements of these genes
          
          
          returnObjectDifferentialSites<-c("Promoter Predicted Sites" = subsetByOverlaps(UnbiasedPredictedMotifs$Promoters, promoterTargetsOfTF),
                                           "Enhancer Predicted Sites" = subsetByOverlaps(UnbiasedPredictedMotifs$Enhancers, enhancerTargetsOfTF))%>%unlist()
          
          
          assign("PredictedTFBS", returnObjectDifferentialSites, .GlobalEnv)
          
        }
        else {
          
          ## Returning unbiased results without transcriptomic data
          
          
          returnObjectUnbaised<-c(
            "Promoter Predicted Sites"= subsetByOverlaps(UnbiasedPredictedMotifs$Promoters, promoterTargetsOfTF),
            "Promoter Predicted Sites"= subsetByOverlaps(UnbiasedPredictedMotifs$Enhancers, promoterTargetsOfTF)
          )%>%unlist()
          
          
          assign("PredictedTFBS", returnObjectUnbaised, .GlobalEnv)
        }
        
        
        ################################
        ###Genome browser part
        ###############################
        chrM<-input$chrM
        assign("chrM", chrM, .GlobalEnv)
        
        humanIdeogramTrack<-IdeogramTrack(chromosome = input$chrM, genome="hg19",name= "Ideogram")
        gHumanTrack<-GenomeAxisTrack(name= "Axis")
        
        assign("humanIdeogramTrack", humanIdeogramTrack, .GlobalEnv)
        assign("gHumanTrack", gHumanTrack, .GlobalEnv)
        
        
        ###############################################
        ####Identifying motifs in CRM regions
        ################################################
        
        
        
        PredictedTFBSTrack<-PredictedTFBS%>%AnnotationTrack(genome = "hg19", 
                                                            stacking = "dense",
                                                            strand= "*",
                                                            col.line="black",
                                                            name="Predicted TFBS")
        
        
        
        assign("PredictedTFBSTrack", PredictedTFBSTrack, .GlobalEnv)
        
        ##ChromHMM Track Generator specifically for humans
       
        
        chromatinStatesTrack<-chromHMMTrackGenerator(gen="hg19", 
                                                     chr= input$chrM, 
                                                     from  = input$fromM,
                                                     to = input$toM,
                                                     bedFile = chromatinState,
                                                     featureDisplay = "all",
                                                     colorcase='roadmap15')
        
        # Gene Track with symbols :D
        knownGenes <- GeneRegionTrack(TxDb.Hsapiens.UCSC.hg19.knownGene, 
                                      genome="hg19", 
                                      chromosome="chrX", 
                                      showId=TRUE,
                                      geneSymbol=TRUE, 
                                      name="UCSC")
        
        symbols <- unlist(mapIds(org.Hs.eg.db, gene(knownGenes),
                                 "SYMBOL", "ENTREZID", 
                                 multiVals = "first"))
        symbol(knownGenes) <- symbols[gene(knownGenes)]
        
        #Promoter and Enhancer Tracks for each chormosome Track
        promotertrackChromosomeSpecific<-promoterTracks%>%subset(. , 
                                                                 seqnames==input$chrM)%>%AnnotationTrack(., name= "PromoterTrack", 
                                                                                                         genome="hg19")
        geneTrackChromosomeSpecific<-knownGenes
        EnhancersHumanChromosomeSpecific<-EnhancersHuman%>%subset(. ,
                                                                  seqnames==input$chrM)%>%AnnotationTrack(., name = "Enhancers",
                                                                                                          genome = "hg19")
        
        # Raw Motif Instances
        RawMotifInstancesTrack<-subset(genomicPositions, 
                                       seqnames==input$chrM & start > input$fromM & end< input$toM)%>%AnnotationTrack(.,genome = "hg19",
                                                                                                                      stacking = "dense", 
                                                                                                                      col.line="black",
                                                                                                                      name="All Motifs")
        
        
        
        assign("chromatinStatesTrack", chromatinStatesTrack, .GlobalEnv)
        assign("RawMotifInstancesTrack", RawMotifInstancesTrack, .GlobalEnv)
        assign("EnhancersHumanChromosomeSpecific", EnhancersHumanChromosomeSpecific, .GlobalEnv)
        assign("promotertrackChromosomeSpecific", promotertrackChromosomeSpecific, .GlobalEnv)
        assign("geneTrackChromosomeSpecific", geneTrackChromosomeSpecific, .GlobalEnv)
        assign("knownGenes", knownGenes, .GlobalEnv)
        
        
        
        
        plotTracks(trackList =c(humanIdeogramTrack,
                                gHumanTrack, 
                                EnhancersHumanChromosomeSpecific,
                                PredictedTFBSTrack,
                                RawMotifInstancesTrack, 
                                promotertrackChromosomeSpecific, 
                                geneTrackChromosomeSpecific,
                                chromatinStatesTrack), 
                   sizes= c(1,1,1,1,1,1,3),
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
        
        
      }else if(!chrM==input$chrM){ 
        chrM<-input$chrM
        assign("chrM", chrM, .GlobalEnv)
        humanIdeogramTrack<-IdeogramTrack(chromosome = input$chrM, genome="hg19",name= "Ideogram")
        
        
        chromHMM_RoadMapAll<-lapply(c("Pancreas",
                                      "PancreasIslets",
                                      "fetalBrainFemale",
                                      "fetalBrainMale",
                                      "H9NeuronCells",
                                      "H9NeuronProgenitorCells"), function(x){chromHMMTrackGenerator(gen="hg19",
                                                                                                     chr=input$chrM,
                                                                                                     from = input$fromM,
                                                                                                     to = input$toM,
                                                                                                     bedFile = x,
                                                                                                     featureDisplay = "all", 
                                                                                                     colorcase='roadmap15' )})
        
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
        
        #Promoter and Motif Track
        promotertrackChromosomeSpecific<-promoterTracks%>%subset(. , 
                                                                 seqnames==input$chrM)%>%AnnotationTrack(., name= "PromoterTrack", 
                                                                                                         genome="hg19")
        geneTrackChromosomeSpecific<-knownGenes
        EnhancersHumanChromosomeSpecific<-EnhancersHuman%>%subset(. ,
                                                                  seqnames==input$chrM)%>%AnnotationTrack(., name = "Enhancers",
                                                                                                          genome = "hg19")
        
        Arx6merHumanTrack<-subset(arxMotifsHumanRaw, 
                                  seqnames==input$chrM & start > input$fromM & end< input$toM)%>%AnnotationTrack(.,genome = "hg19",
                                                                                                                 stacking = "dense", 
                                                                                                                 col.line="black",
                                                                                                                 name="ALL ARX Motifs",
                                                                                                                 feature= (mcols(.))$Model)
        displayPars(Arx6merHumanTrack) <- list(`6mer` = "#e6194b", 
                                               `T6mer2` = "#3cb44b", 
                                               `Jolma` = "#0082c8", 
                                               `P6mer4` = "#008080")
        
        
        
        
        assign("humanIdeogramTrack", humanIdeogramTrack, .GlobalEnv)
        assign("chromHMM_RoadMapAll", chromHMM_RoadMapAll, .GlobalEnv)
        assign("Arx6merHumanTrack", Arx6merHumanTrack, .GlobalEnv)
        assign("EnhancersHumanChromosomeSpecific", EnhancersHumanChromosomeSpecific, .GlobalEnv)
        assign("promotertrackChromosomeSpecific", promotertrackChromosomeSpecific, .GlobalEnv)
        assign("geneTrackChromosomeSpecific", geneTrackChromosomeSpecific, .GlobalEnv)
        assign("knownGenes", knownGenes, .GlobalEnv)
        
        
        
        
        if(input$contactProbabilities==TRUE){
          plotTracks(trackList = c(humanIdeogramTrack,
                                   gHumanTrack,
                                   contactProbabilities, 
                                   EnhancersHumanChromosomeSpecific,
                                   ARXEnhancerMotifs,
                                   Arx6merHumanTrack, 
                                   promotertrackChromosomeSpecific, 
                                   geneTrackChromosomeSpecific,
                                   chromHMM_RoadMapAll), 
                     sizes= c(1,1,2,1,1,1,1,3,rep(1,6)),
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
        } else{
          
          plotTracks(trackList =c(humanIdeogramTrack,
                                  gHumanTrack,
                                  interactionsHumanBrain, 
                                  EnhancersHumanChromosomeSpecific,
                                  ARXEnhancerMotifs,
                                  Arx6merHumanTrack, 
                                  promotertrackChromosomeSpecific, 
                                  geneTrackChromosomeSpecific,
                                  chromHMM_RoadMapAll), 
                     sizes= c(1,1,2,1,1,1,1,3,rep(1,6)),
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
        }} else if(input$contactProbabilities==TRUE) {
          
          
          chromHMM_RoadMapAll<-lapply(c("Pancreas",
                                        "PancreasIslets",
                                        "fetalBrainFemale",
                                        "fetalBrainMale",
                                        "H9NeuronCells",
                                        "H9NeuronProgenitorCells"), function(x){chromHMMTrackGenerator(gen="hg19",
                                                                                                       chr=input$chrM,
                                                                                                       from = input$fromM,
                                                                                                       to = input$toM,
                                                                                                       bedFile = x,
                                                                                                       featureDisplay = "all", 
                                                                                                       colorcase='roadmap15' )})
          # Arx All motifs Track
          
          Arx6merHumanTrack<-subset(arxMotifsHumanRaw, 
                                    seqnames==input$chrM & start > input$fromM & end< input$toM)%>%AnnotationTrack(.,genome = "hg19",
                                                                                                                   stacking = "dense", 
                                                                                                                   col.line="black",
                                                                                                                   name="ALL ARX Motifs",
                                                                                                                   feature= (mcols(.))$Model)
          
          displayPars(Arx6merHumanTrack) <- list(`6mer` = "#e6194b", 
                                                 `T6mer2` = "#3cb44b", 
                                                 `Jolma` = "#0082c8", 
                                                 `P6mer4` = "#008080")
          
          
          
          assign("chromHMM_RoadMapAll", chromHMM_RoadMapAll, .GlobalEnv)
          assign("Arx6merHumanTrack", Arx6merHumanTrack, .GlobalEnv)
          assign("EnhancersHumanChromosomeSpecific", EnhancersHumanChromosomeSpecific, .GlobalEnv)
          assign("promotertrackChromosomeSpecific", promotertrackChromosomeSpecific, .GlobalEnv)
          assign("knownGenes", knownGenes, .GlobalEnv)
          
          
          
          plotTracks(trackList = c(humanIdeogramTrack,
                                   gHumanTrack,
                                   contactProbabilities, 
                                   EnhancersHumanChromosomeSpecific,
                                   ARXEnhancerMotifs,
                                   Arx6merHumanTrack, 
                                   promotertrackChromosomeSpecific, 
                                   geneTrackChromosomeSpecific,
                                   chromHMM_RoadMapAll), 
                     sizes= c(1,1,2,1,1,1,1,3,rep(1,6)),
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
          
        }else {
          
          chromHMM_RoadMapAll<-lapply(c("Pancreas",
                                        "PancreasIslets",
                                        "fetalBrainFemale",
                                        "fetalBrainMale",
                                        "H9NeuronCells",
                                        "H9NeuronProgenitorCells"), function(x){chromHMMTrackGenerator(gen="hg19",
                                                                                                       chr=input$chrM,
                                                                                                       from = input$fromM,
                                                                                                       to = input$toM,
                                                                                                       bedFile = x,
                                                                                                       featureDisplay = "all", 
                                                                                                       colorcase='roadmap15' )})
          # Arx All motifs Track
          
          Arx6merHumanTrack<-subset(arxMotifsHumanRaw, 
                                    seqnames==input$chrM & start > input$fromM & end< input$toM)%>%AnnotationTrack(.,genome = "hg19",
                                                                                                                   stacking = "dense", 
                                                                                                                   col.line="black",
                                                                                                                   name="ALL ARX Motifs",
                                                                                                                   feature= (mcols(.))$Model)
          
          displayPars(Arx6merHumanTrack) <- list(`6mer` = "#e6194b", 
                                                 `T6mer2` = "#3cb44b", 
                                                 `Jolma` = "#0082c8", 
                                                 `P6mer4` = "#008080")
          
          
          
          plotTracks(trackList = c(humanIdeogramTrack,
                                   gHumanTrack,
                                   interactionsHumanBrain, 
                                   EnhancersHumanChromosomeSpecific,
                                   ARXEnhancerMotifs,
                                   Arx6merHumanTrack, 
                                   promotertrackChromosomeSpecific, 
                                   geneTrackChromosomeSpecific,
                                   chromHMM_RoadMapAll), 
                     sizes= c(1,1,2,1,1,1,1,3,rep(1,6)),
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
  
 output$DataTablePredictedSites<- renderDataTable(
    subset(PredictedTFBSTrack, start>=input$fromM & end<= input$toM, seqnames>= input$chrM)
  )
  
  
  
  
  
  
  
  
  
  
  
  
  
  
  
  
  
  
  
  
  
  
  
  
  
  
  
  
  
  
  
  
  output$LegendsPlot<- renderImage({
    
    list(
      src = "www/EpigenomicsRoadMapLegendHMM.jpeg",
      contentType = "image/jpeg",
      alt = "Human/Epigenomics Road Map Legend"
    )
    
  }, deleteFile = FALSE)
  
  
  output$LegendsPlotMouse<- renderImage({
    
    list(
      src = "www/mm9ChromHMMStates.jpeg",
      contentType = "image/jpeg",
      alt = "Human/Epigenomics Road Map Legend"
    )
    
  }, deleteFile = FALSE)
  
  
})
