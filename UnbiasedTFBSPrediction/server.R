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
require(tibble)


shinyServer(function(input, output) {


  ######################################################
  ## Querying the JASPAR data Base for the matrices and names of the TFs
  #######################################################

  JASPARR2016Matrices <- getMatrixSet(JASPAR2016, list() )
  namedJasparDataBase<-lapply(1:length(JASPARR2016Matrices), function(x){
    names(JASPARR2016Matrices[x])<-JASPARR2016Matrices[[x]]@name
  })

  names(JASPARR2016Matrices) <- namedJasparDataBase



  ##########################################################
  ## Visualize the PWM of the TF
  ##########################################################

  output$VisualizeTFMotif<-renderPlot({

    if(input$CustomTFMatrix==TRUE){

      PathToCustomMatrix <- input$JasparCustomMatrix

      #Import and convert to custom Matrix
      TextTranscriptionFactor <- read_delim(PathToCustomMatrix$datapath,
                                            "\t", escape_double = FALSE, col_names = FALSE,
                                            trim_ws = TRUE, skip = 1)%>%as.matrix()


      #Sum the first column and divide all numbers by that to normalize and make it a frequency matrix.

      PWMTextTranscriptionFactor <- (TextTranscriptionFactor/
                                       TextTranscriptionFactor%>%
                                       colSums())%>%as.matrix()

      assign("matrixForMatching", TextTranscriptionFactor, .GlobalEnv)

      seqLogo(PWMTextTranscriptionFactor)


    } else if(input$TypeInSequence==TRUE){


      TranscriptionFactor <- round(PWM(input$CustomDNASequence)*dim((PWM(input$CustomDNASequence)))[2]) %T>%seqLogo()

      assign("matrixForMatching", TranscriptionFactor, .GlobalEnv)

    } else {

      ## Rewrite this because some of these matrices dont have even numbers of nucelotides for eahc position #WTF
      # OKay so the issue appears to be with repeat numbers Eg 0.333333333333333333333. NOt sure how to deal with this.


      (JASPARR2016Matrices[[paste0(input$TranscriptionFactorPWM)]]@profileMatrix %>%
         apply(MARGIN = 2, function(x){x / sum(x)}))%>%seqLogo()


      JasparMatrixForMatching <- JASPARR2016Matrices[[paste0(input$TranscriptionFactorPWM)]]@profileMatrix
      assign("matrixForMatching", JasparMatrixForMatching, .GlobalEnv)
    }
  })

  ###########################################################
  #### Data Inputs
  ###########################################################

  ##Importing Bed File for the promoter regions with gene symbol names
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


  #########################################################
  # Initiating the annotation hub object to get the  Chromatin state GRANGE off of it
  ##############################################################
  ah <- AnnotationHub()

  assign("ah", ah, .GlobalEnv)

  EpigenoomicsConverter <- cbind.data.frame("H1 Cell Line"="E003",
                                            "H1 BMP4 Derived Mesendoderm Cultured Cells"="E004",
                                            "H1 BMP4 Derived Trophoblast Cultured Cells"="E005",
                                            "H1 Derived Mesenchymal Stem Cells"="E006",
                                            "H1 Derived Neuronal Progenitor Cultured Cells"="E007",
                                            "H9 Cell Line"="E008",
                                            "IMR90 fetal lung fibroblasts Cell Line"="E017",
                                            "iPS DF 6.9 Cell Line"="E021",
                                            "iPS DF 19.11 Cell Line"="E022",
                                            "Breast variant Human Mammary Epithelial Cells (vHMEC)"="E028",
                                            "Primary monocytes from peripheral blood"="E029",
                                            "Primary B cells from peripheral blood"="E032",
                                            "Primary T cells from cord blood"="E033",
                                            "Primary T cells from peripheral blood"="E034",
                                            "Primary Natural Killer cells from peripheral blood"="E046",
                                            "Primary hematopoietic stem cells G-CSF-mobilized Female"="E050",
                                            "Primary hematopoietic stem cells G-CSF-mobilized Male"="E051",
                                            "Foreskin Fibroblast Primary Cells skin01"="E055",
                                            "Foreskin Fibroblast Primary Cells skin02"="E056",
                                            "Foreskin Keratinocyte Primary Cells skin02"="E057",
                                            "Foreskin Melanocyte Primary Cells skin01"="E059",
                                            "Fetal Adrenal Gland"="E080",
                                            "Fetal Brain Male"="E081",
                                            "Fetal Brain Female"="E082",
                                            "Fetal Heart"="E083",
                                            "Fetal Intestine Large"="E084",
                                            "Fetal Intestine Small"="E085",
                                            "Fetal Kidney"="E086",
                                            "Fetal Lung"="E088",
                                            "Fetal Muscle Trunk"="E089",
                                            "Fetal Muscle Leg"="E090",
                                            "Placenta"="E091",
                                            "Fetal Stomach"="E092",
                                            "Fetal Thymus"="E093",
                                            "Gastric"="E094",
                                            "Ovary"="E097",
                                            "Pancreas"="E098",
                                            "Psoas Muscle"="E100",
                                            "Small Intestine"="E109",
                                            "A549 EtOH 0.02pct Lung Carcinoma Cell Line"="E114",
                                            "GM12878 Lymphoblastoid Cell Line"="E116",
                                            "HeLa-S3 Cervical Carcinoma Cell Line"="E117",
                                            "HepG2 Hepatocellular Carcinoma Cell Line"="E118",
                                            "HMEC Mammary Epithelial Primary Cells"="E119",
                                            "HSMM Skeletal Muscle Myoblasts Cell Line"="E120",
                                            "HSMM cell derived Skeletal Muscle Myotubes Cell Line"="E121",
                                            "HUVEC Umbilical Vein Endothelial Cells Cell Line"="E122",
                                            "K562 Leukemia Cell Line"="E123",
                                            "Monocytes-CD14+ RO01746 Cell Line"="E124",
                                            "NH-A Astrocytes Cell Line"="E125",
                                            "NHDF-Ad Adult Dermal Fibroblast Primary Cells"="E126",
                                            "NHEK-Epidermal Keratinocyte Primary Cells"="E127",
                                            "NHLF Lung Fibroblast Primary Cells"="E128",
                                            "H9 Derived Neuronal Progenitor Cultured Cells"="E009",
                                            "H9 Derived Neuron Cultured Cells"="E010",
                                            "Dnd41 TCell Leukemia Cell Line"="E115",
                                            "Osteoblast Primary Cells"="E129",
                                            "hESC Derived CD184+ Endoderm Cultured Cells"="E011",
                                            "hESC Derived CD56+ Ectoderm Cultured Cells"="E012",
                                            "hESC Derived CD56+ Mesoderm Cultured Cells"="E013",
                                            "HUES48 Cell Line"="E014",
                                            "HUES6 Cell Line"="E015",
                                            "HUES64 Cell Line"="E016",
                                            "iPS-18 Cell Line"="E019",
                                            "iPS-20b Cell Line"="E020",
                                            "Bone Marrow Derived Cultured Mesenchymal Stem Cells"="E026",
                                            "Primary T helper memory cells from peripheral blood 2"="E037",
                                            "Primary T helper naive cells from peripheral blood"="E038",
                                            "Primary T helper naive cells from peripheral blood"="E039",
                                            "Primary T helper memory cells from peripheral blood 1"="E040",
                                            "Primary T helper cells PMA-I stimulated"="E041",
                                            "Primary T helper 17 cells PMA-I stimulated"="E042",
                                            "Primary T helper cells from peripheral blood"="E043",
                                            "Primary T regulatory cells from peripheral blood"="E044",
                                            "Primary T cells effector/memory enriched from peripheral blood"="E045",
                                            "Primary T killer naive cells from peripheral blood"="E047",
                                            "Primary T killer memory cells from peripheral blood"="E048",
                                            "Mesenchymal Stem Cell Derived Chondrocyte Cultured Cells"="E049",
                                            "Foreskin Keratinocyte Primary Cells skin03"="E058",
                                            "Foreskin Melanocyte Primary Cells skin03"="E061",
                                            "Primary mononuclear cells from peripheral blood"="E062",
                                            "Adipose Nuclei"="E063",
                                            "Aorta"="E065",
                                            "Liver"="E066",
                                            "Brain Angular Gyrus"="E067",
                                            "Brain Anterior Caudate"="E068",
                                            "Brain Cingulate Gyrus"="E069",
                                            "Brain Hippocampus Middle"="E071",
                                            "Brain Inferior Temporal Lobe"="E072",
                                            "Brain_Dorsolateral_Prefrontal_Cortex"="E073",
                                            "Brain Substantia Nigra"="E074",
                                            "Colonic Mucosa"="E075",
                                            "Colon Smooth Muscle"="E076",
                                            "Duodenum Smooth Muscle"="E078",
                                            "Esophagus"="E079",
                                            "Pancreatic Islets"="E087",
                                            "Left Ventricle"="E095",
                                            "Lung"="E096",
                                            "Placenta Amnion"="E099",
                                            "Rectal Mucosa Donor 29"="E101",
                                            "Rectal Mucosa Donor 31"="E102",
                                            "Rectal Smooth Muscle"="E103",
                                            "Right Atrium"="E104",
                                            "Right Ventricle"="E105",
                                            "Sigmoid Colon"="E106",
                                            "Skeletal Muscle Female"="E108",
                                            "Stomach Smooth Muscle"="E111",
                                            "Thymus"="E112",
                                            "Spleen"="E113",
                                            "ES-I3 Cell Line"="E001",
                                            "ES-WA7 Cell Line"="E002",
                                            "iPS-15b Cell Line"="E018",
                                            "Mesenchymal Stem Cell Derived Adipocyte Cultured Cells"="E023",
                                            "ES-UCSF4  Cell Line"="E024",
                                            "Adipose Derived Mesenchymal Stem Cell Cultured Cells"="E025",
                                            "Breast Myoepithelial Primary Cells"="E027",
                                            "Primary neutrophils from peripheral blood"="E030",
                                            "Primary B cells from cord blood"="E031",
                                            "Primary hematopoietic stem cells"="E035",
                                            "Primary hematopoietic stem cells short term culture"="E036",
                                            "Muscle Satellite Cultured Cells"="E052",
                                            "Cortex derived primary cultured neurospheres"="E053",
                                            "Ganglion Eminence derived primary cultured neurospheres"="E054",
                                            "Brain Germinal Matrix"="E070",
                                            "Duodenum Mucosa"="E077",
                                            "Skeletal Muscle Male"="E107",
                                            "Stomach Mucosa"="E110")%>%t()%>%set_colnames("Epigenomic Road Map")%>%as.data.frame()%>%rownames_to_column(var= "Name Of Cell")


  ###############################################################################################
  ###Promoter Enhancer assoication data table to assign gene targets to enhancers using CAGE expression
  ###############################################################################################

  #Permissive Enhancer
  permissiveEnhancer <- import.bed("../DataFiles/Enhancers Track/Human/human_permissive_enhancers_phase_1_and_2.bed.gz")

  assign("permissiveEnhancer", permissiveEnhancer, .GlobalEnv)

  #Read in the table
  EnhancerPromoterAssoications <- read_delim("../DataFiles/Enhancers Track/Human/hg19_enhancer_promoter_correlations_distances_cell_type.txt.gz",
                                            "\t", escape_double = FALSE, trim_ws = TRUE)

  assign("EnhancerPromoterAssoications", EnhancerPromoterAssoications, .GlobalEnv)



  PromotersAssoicatedWithEnhancers<-separate(EnhancerPromoterAssoications, col=promoter, into= c("promoter", "strand"), sep= ',')%>%
    separate(., col= promoter, into= c("chr", "start", "end"))%>%makeGRangesFromDataFrame(keep.extra.columns = TRUE)




  ## Convert Enhancer To Granges
  EnhancerGrangeWithTargets <- EnhancerPromoterAssoications$enhancer%>%GRanges()

  # Add everything to the metadata columns
  mcols(EnhancerGrangeWithTargets) <- PromotersAssoicatedWithEnhancers


  # Select for significant values only, - inlcuding taking out the 0.0000 as im not confident they're legit.
  EnhancerGrangeWithTargetsSiginificant <- subset(EnhancerGrangeWithTargets, `FDR`!=0 & `FDR`<0.05)

  # Generate a vector for targets of enhancers by overlapping promoter regions with promoter GRanges
  IntersectBetweenOverlappingPromoterRanges <- findOverlaps(EnhancerGrangeWithTargetsSiginificant$X, promoterTracks)

  # Subset for Grange for promoter tarcks by the hits from the findOverlaps vector above (this will also order it)
  GenesRegulatedByEnhancers <- as.data.frame(mcols(promoterTracks))[IntersectBetweenOverlappingPromoterRanges%>%subjectHits(),]

  # Subset the EnhancerGrange by the same vector (except the vector for this Grange) hence ordering them both in the same order
  EnhancersWithGeneTargetsGrange <- EnhancerGrangeWithTargetsSiginificant[IntersectBetweenOverlappingPromoterRanges%>%queryHits()]

  ## Combine them
  mcols(EnhancersWithGeneTargetsGrange) <- cbind.data.frame(mcols(EnhancersWithGeneTargetsGrange),GenesRegulatedByEnhancers)

  ## Removing the redudent data! (I had to use this due to dplyr not enjoying columns having the same name)
  mcols(EnhancersWithGeneTargetsGrange) <- mcols(EnhancersWithGeneTargetsGrange)[12:20]


  #We can use this to identify gene targets of motifs in enhancer regions
  assign("EnhancersWithGeneTargetsGrange", EnhancersWithGeneTargetsGrange, .GlobalEnv)



  #ChromHMM Function from the coMet package (Edited for speed and so it can use the Grange input over a dataframe)
  ##ChromHMM Track Generator specifically for human

  chromHMMTrackGenerator<-function (gen = "hg19",
                                    chr,
                                    from,
                                    to,
                                    bedFile,
                                    featureDisplay = featureDisplay,
                                    colorcase = "roadmap15",
                                    name )
  {
    desiredRegion <- subset(bedFile, end > from &
                              start < to & seqnames == chr)

    track <- AnnotationTrack(desiredRegion,
                             stacking = "dense",
                             col.line="black",
                             feature = (mcols(desiredRegion))$abbr,
                             genome = "hg19",
                             strand= "*",
                             name = paste(name, "chromHMM"))

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

  IntearctionTrack <- CageExpressionGenomicIntearctions%>%InteractionTrack()
  
  
  displayPars(IntearctionTrack) = list(col.interactions="red", 
                                        col.anchors.fill ="blue",
                                        col.anchors.line = "black",
                                        interaction.dimension="height", 
                                        interaction.measure ="counts",
                                        plot.trans=FALSE,
                                        plot.outside = TRUE, 
                                        col.outside="lightblue", 
                                        anchor.height = 0.1)
  

  assign("IntearctionTrack", IntearctionTrack, .GlobalEnv)




  output$EnhancerTable <- renderDataTable({

    input$ComputeTranscriptionFactorSites

    isolate({


      #####################################################################
      ###Start of the pipe Line
      ###################################################################


      ### Input Transcription Factor Matrix from the seqlogo function

      assign("TranscriptionFactorPWM", matrixForMatching, .GlobalEnv )


      genomicLocationOfMotifs <- matchPWM(matrixForMatching,
                                        BSgenome.Hsapiens.UCSC.hg19,
                                        paste0(input$MatchPercentage,'%'))

      assign("genomicLocationOfMotifs", genomicLocationOfMotifs, .GlobalEnv)



      ### Identify which of these motifs are located in enhancer regions
      MotifsInEnhancers <- subsetByOverlaps(genomicLocationOfMotifs ,
                                            EnhancersWithGeneTargetsGrange)

      MotifsInPermissiveEnhancers <- subsetByOverlaps(genomicLocationOfMotifs,
                                                      permissiveEnhancer)

      MotifsInPromoters <- subsetByOverlaps(genomicLocationOfMotifs,
                                            promoterTracks)

      if(input$Conserved == TRUE){

        ConservedValue <- input$Conserved

        assign("ConservedValue", ConservedValue , .GlobalEnv)

        if(!exists("ConservedRegionsInPromoters")){

          # Import the Conserved region and run the following code
          if(!file.exists("../DataFiles/Conserved Region/Human/phyloP100WayUCSCDataTrack.bw")){



            dir.create("../DataFiles/Conserved Region")


            dir.create("../DataFiles/Conserved Region/Human")

            #Download it
            download.file(url = "ftp://hgdownload.soe.ucsc.edu/goldenPath/hg19/phyloP100way/hg19.100way.phyloP100way.bw",
                          destfile = "../DataFiles/Conserved Region/Human/phyloP100WayUCSCDataTrack.bw",
                          method = "curl")

            # Read in the relative promoter regions and subset for regions with a score greater than 0 (or the neutral drift)

            ConservedRegionsInPromoters <-
              import("../DataFiles/Conserved Region/Human/phyloP100WayUCSCDataTrack.bw", which= MotifsInPromoters)%>%subset(. , score>0)


            assign("ConservedRegionsInPromoters", ConservedRegionsInPromoters, .GlobalEnv)

            # Identify motifs whoses entire motif is in this region.

            MotifsInConservedPromoterRegions <- subset(MotifsInPromoters,
                                                       countOverlaps(MotifsInPromoters,
                                                                     ConservedRegionsInPromoters)>=genomicLocationOfMotifs$string[[1]]%>%length)


          } else {

            ConservedRegionsInPromoters <-
              import("../DataFiles/Conserved Region/Human/phyloP100WayUCSCDataTrack.bw", which= MotifsInPromoters)%>%subset(. , score>0)


            assign("ConservedRegionsInPromoters", ConservedRegionsInPromoters, .GlobalEnv)

            # Identify motifs whoses entire motif is in this region.

            MotifsInConservedPromoterRegions<-subset(MotifsInPromoters,
                                                     countOverlaps(MotifsInPromoters,
                                                                   ConservedRegionsInPromoters)>=genomicLocationOfMotifs$string[[1]]%>%length)
          }

        } else {

          ConservedRegionsInPromoters <-
            import("../DataFiles/Conserved Region/Human/phslopy100WayUCSCDataTrack.bw", which= MotifsInPromoters)%>%subset(. , score>0)


          assign("ConservedRegionsInPromoters", ConservedRegionsInPromoters, .GlobalEnv)

          # Identify motifs whoses entire motif is in this region.

          MotifsInConservedPromoterRegions <- subset(MotifsInPromoters,
                                                     countOverlaps(MotifsInPromoters,
                                                                   ConservedRegionsInPromoters)>=genomicLocationOfMotifs$string[[1]]%>%length)


        }

      } else {

        ConservedValue <- input$Conserved

        assign("ConservedValue", ConservedValue , .GlobalEnv)

        # If we are not looking for conserved promoter regions
        MotifsInConservedPromoterRegions <- subsetByOverlaps(genomicLocationOfMotifs, promoterTracks)

      }

      #### Combing the promoter and enhancer motifs into a single Grange List object

      MotifsInPromotersAndEnhancers <- list("Promoters" = MotifsInConservedPromoterRegions,
                                            "Enhancers" = MotifsInEnhancers,
                                            "Permissive Enhancers" = MotifsInPermissiveEnhancers)

      assign("MotifsInPromotersAndEnhancers", MotifsInPromotersAndEnhancers, .GlobalEnv)





      ################################################
      #### Cell Type specific epigneonic Analysis
      ##################################################

      CellTypeToPredict <- input$CellTypeToPredict

      assign("CellTypeToPredict", CellTypeToPredict, .GlobalEnv )

      ### Download the Chromatin state file From the annotation hub

      epiFiles <- query(ah, c(paste0(CellTypeToPredict,"_15_coreMarks_mnemonics"), "EpigenomeRoadMap") )

      chromatinState <- epiFiles[[paste0("AH", epiFiles@.db_uid)]]


      #Assign the bedfile to the global environment for analysis later on downstream

      assign("chromatinState", chromatinState, .GlobalEnv)


      ## Name for the chromHMM Track Down Stream

      NameOfCell <- subset(EpigenoomicsConverter, `Epigenomic Road Map`==CellTypeToPredict)$`Name Of Cell`

      assign("NameOfCell", NameOfCell, .GlobalEnv)


      ## Subsetting the Chromatin states for active states to identify motifs in these regions
      ActiveChromatinStates <- c("10_TssBiv",
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


      ActiveChromatinRegions <- chromatinState[chromatinState$abbr %in% ActiveChromatinStates,]

      #Overlap active reginos with motifs in CRMs giving us a completely unbiased set of results
      UnbiasedPredictedMotifs <- lapply(MotifsInPromotersAndEnhancers,function(x){


        OverlapBetweenRangesObject <- findOverlaps( x , ActiveChromatinRegions)

        EpigenomicLocation <- ActiveChromatinRegions$name[OverlapBetweenRangesObject%>%subjectHits()]

        OrderingQueryGRange <- x [OverlapBetweenRangesObject%>%queryHits()]

        mcols(OrderingQueryGRange) <- cbind.data.frame(mcols(OrderingQueryGRange),
                                                       "Epigenomic Location" = EpigenomicLocation)

        OrderingQueryGRange

      })

      assign("UnbiasedPredictedMotifs", UnbiasedPredictedMotifs, .GlobalEnv )

      ###################################################
      ###Assigning Gene targets to eac motif based on the assoication of the CRM it is located within.
      ####################################################

      ### First we identify the genes that are regulated by our unbiased predicted results using genomic annotation


      #Overlap the promoter regions of genes with unbiased motifs returning promoters with a predicted TFBS
      OverlappingRangeOfMOtifsInPromoters <- findOverlaps(promoterTracks, UnbiasedPredictedMotifs$Promoters)

      unbiasedPromoterMotifs <- UnbiasedPredictedMotifs$Promoters[OverlappingRangeOfMOtifsInPromoters%>%subjectHits()]


      mcols(unbiasedPromoterMotifs) <- cbind.data.frame(mcols(unbiasedPromoterMotifs),
                                                        "Genes Regulated" = promoterTracks[OverlappingRangeOfMOtifsInPromoters%>%queryHits()]$hg19.kgXref.geneSymbol,
                                                        "UCSC Transcript ID" = promoterTracks[OverlappingRangeOfMOtifsInPromoters%>%queryHits()]$hg19.kgXref.kgID)


      ## Now lets get the promoters of genes regulated by motifs in enhancres
      OverlappingRangeEnhancersMotifs <- findOverlaps(EnhancersWithGeneTargetsGrange, UnbiasedPredictedMotifs$Enhancers)

      EnhancerTargets <- cbind.data.frame(EnhancersWithGeneTargetsGrange$hg19.kgXref.geneSymbol,
                                          EnhancersWithGeneTargetsGrange$X.hg19.knownGene.name)[OverlappingRangeEnhancersMotifs%>%queryHits(),]

      colnames(EnhancerTargets) <- c("Genes Regulated" , "UCSC Transcript ID")

      MotifsInEnhancers <- UnbiasedPredictedMotifs$Enhancers[OverlappingRangeEnhancersMotifs%>%subjectHits()]

      mcols(MotifsInEnhancers) <- cbind.data.frame(mcols(MotifsInEnhancers), EnhancerTargets)

      EnhancerMotifsNotTargeted <- UnbiasedPredictedMotifs$`Permissive Enhancers`[!UnbiasedPredictedMotifs$`Permissive Enhancers` %in% MotifsInEnhancers ]

      mcols(EnhancerMotifsNotTargeted) <- cbind.data.frame(mcols(EnhancerMotifsNotTargeted),
                                                           "Genes Regulated" = "NA",
                                                           "UCSC Transcript ID" = "NA")

      UnbiasedMotifsPredicted <- c("Promoters With Gene Targets" = unbiasedPromoterMotifs,
                                   "Enhancers With Gene Targets" = sort(c(MotifsInEnhancers, EnhancerMotifsNotTargeted))%>%unlist())



      assign("UnbiasedMotifsPredicted", UnbiasedMotifsPredicted, .GlobalEnv)





      #######################################################
      ## Left Joining with differentially expressed gene list
      #########################################################

      if(input$DifferentialExpressedGenes == TRUE){


        differenitallyExpressedGenesPath <- input$differenitallyExpressedGenesList

        differenitallyExpressedGenesList <-
          read.table(differenitallyExpressedGenesPath$datapath , quote="\"", stringsAsFactors = FALSE)

        assign("differenitallyExpressedGenesList", differenitallyExpressedGenesList, .GlobalEnv)

        ## Genes who showed differenital expression with an enhancer that was correlated with its expression
        enhancerTargetsOfTF <- subset(UnbiasedMotifsPredicted$`Enhancers With Gene Targets`,
                                      `Genes Regulated` %in% differenitallyExpressedGenesList[,1])


        ##Genes who showed differenital expression with promoter targets
        promoterTargetsOfTF <- subset(UnbiasedMotifsPredicted$`Promoters With Gene Targets`,
                                      `Genes Regulated` %in% differenitallyExpressedGenesList[,1])



        ##Predicted Sites in the Regulatory Elements of these genes


        returnObjectDifferentialSites <- c("Promoter Predicted Sites" = promoterTargetsOfTF,
                                           "Enhancer Predicted Sites" = enhancerTargetsOfTF)%>%unlist()

        mcols(returnObjectDifferentialSites$`Promoter Predicted Sites`) <- cbind.data.frame(mcols(returnObjectDifferentialSites$`Promoter Predicted Sites`),
                                                                                   "Regulatory Module" = rep("Promoter",
                                                                                                             length(returnObjectDifferentialSites$`Promoter Predicted Sites`)))

        mcols(returnObjectDifferentialSites$`Enhancer Predicted Sites`) <- cbind.data.frame(mcols(returnObjectDifferentialSites$`Enhancer Predicted Sites`),
                                                                                   "Regulatory Module" = rep("Enhancer",
                                                                                                             length(returnObjectDifferentialSites$`Enhancer Predicted Sites`)))

        GenomeBrowserUnbiasedSites <- c(returnObjectDifferentialSites$`Promoter Predicted Sites`,
                                        returnObjectDifferentialSites$`Enhancer Predicted Sites`)%>%unlist()

        assign("PredictedTFBS", GenomeBrowserUnbiasedSites, .GlobalEnv)

        assign("returnObjectDifferentialSites", returnObjectDifferentialSites, .GlobalEnv)

      }
      else {

        ## Returning unbiased results without transcriptomic data


        returnObjectUnbaised<-c(
          "Promoter Predicted Sites"= UnbiasedMotifsPredicted$`Enhancers With Gene Targets`,
          "Enhancer Predicted Sites"= UnbiasedMotifsPredicted$`Promoters With Gene Targets` )%>%unlist()

        assign("returnObjectUnbaised", returnObjectUnbaised, .GlobalEnv)

        mcols(returnObjectUnbaised$`Promoter Predicted Sites`) <- cbind.data.frame(mcols(returnObjectUnbaised$`Promoter Predicted Sites`),
                                                                                   "Regulatory Module" = rep("Promoter",
                                                                                                             length(returnObjectUnbaised$`Promoter Predicted Sites`)))

        mcols(returnObjectUnbaised$`Enhancer Predicted Sites`) <- cbind.data.frame(mcols(returnObjectUnbaised$`Enhancer Predicted Sites`),
                                                                                   "Regulatory Module" = rep("Enhancer",
                                                                                                             length(returnObjectUnbaised$`Enhancer Predicted Sites`)))

        GenomeBrowserUnbiasedSites <- c(returnObjectUnbaised$`Promoter Predicted Sites`,
                                        returnObjectUnbaised$`Enhancer Predicted Sites`)%>%unlist()

        assign("PredictedTFBS", GenomeBrowserUnbiasedSites, .GlobalEnv)
      }

      
      PromoterPredictedSites <-   GenomeBrowserUnbiasedSites%>%as.data.frame()
      
      assign("PromoterPredictedSites", PromoterPredictedSites, .GlobalEnv )
      
      
      
      
     
    
      
    })
  } )


 
      ##################################################################
      #### Download Exports
      ##################################################################

  output$rawMotifPositions <- downloadHandler('rawMotifGenomicPositions.bed', content = function(file) {
    write.table(as.data.frame(genomicLocationOfMotifs), file  ,
                sep="\t",
                row.names = FALSE,
                col.names = FALSE,
                quote = FALSE,
                append = FALSE)
  })    
  
  
  output$returnObjectUnbaised <- downloadHandler('UnbiasedGenomicPositions.bed', content = function(file) {
    write.table( PromoterPredictedSites, file  ,
                 sep="\t",
                 row.names = FALSE,
                 col.names = FALSE,
                 quote = FALSE,
                 append = FALSE)
  })
  
  
  output$ProcessedMotifPositions <- downloadHandler('ProcessedGenomicPositions.bed', content = function(file) {
    
    write.table( PromoterPredictedSites, file  ,
                sep="\t",
                row.names = FALSE,
                col.names = FALSE,
                quote = FALSE,
                append = FALSE)
    })
    
  
  output$RegualtoryModuleMotifs <- downloadHandler('ProcessedGenomicPositions.bed', content = function(file) {
    
    RegulatoryModuleDataFrame <- rbind.data.frame(as.data.frame(MotifsInPromotersAndEnhancers$Promoters),
    as.data.frame(MotifsInPromotersAndEnhancers$`Permissive Enhancers`))
    
    write.table( RegulatoryModuleDataFrame, file  ,
                 sep="\t",
                 row.names = FALSE,
                 col.names = FALSE,
                 quote = FALSE,
                 append = FALSE)
  })













  output$GenomeBrowser <- renderPlot({

    input$GenomeBrowserAction
    isolate({

      if(!exists("chrM")){


        GenomeBrowserTFMatrix <- matrixForMatching

        assign("GenomeBrowserTFMatrix", GenomeBrowserTFMatrix, .GlobalEnv)


        ###############################################################################################
        ###Genome browser part for when you change chromosomes
        ##############################################################################################

        chr <- input$chrM

        # Ideogram Track
        humanIdeogramTrack <- IdeogramTrack(chromosome = chr,
                                            genome="hg19",
                                            name= "Ideogram")

        #Genome Axis Track for apprxoimate location
        gHumanTrack <- GenomeAxisTrack(name= "Genomic Axis Track")

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
                                                                   seqnames==input$chrM)%>%AnnotationTrack(.,
                                                                                                           name= "PromoterTrack",
                                                                                                           genome="hg19")
        
        displayPars(promotertrackChromosomeSpecific) <- 
          list(fill = "red")
        geneTrackChromosomeSpecific <- knownGenes



        EnhancersHumanChromosomeSpecific <- reduce(c(reduce(permissiveEnhancer), reduce(EnhancersWithGeneTargetsGrange)))%>%
          subset(. ,
                 seqnames==input$chrM)%>%AnnotationTrack(.,
                                                         name = "Enhancers",
                                                         genome = "hg19")
        displayPars(EnhancersHumanChromosomeSpecific) <- 
          list(fill = "#911eb4")

        assign("EnhancersHumanChromosomeSpecific", EnhancersHumanChromosomeSpecific, .GlobalEnv)
        assign("promotertrackChromosomeSpecific", promotertrackChromosomeSpecific, .GlobalEnv)
        assign("geneTrackChromosomeSpecific", geneTrackChromosomeSpecific, .GlobalEnv)

        #Chromosome Specific Predicted Motifs
        PredictedTFBSTrack <- PredictedTFBS%>%subset(seqnames==input$chrM)%>%AnnotationTrack(genome = "hg19",
                                                                                           stacking = "dense",
                                                                                           strand= "*",
                                                                                           col.line="black",
                                                                                           name="Predicted TFBS")
        
        displayPars(PredictedTFBSTrack) <-  list(fill = "Black")
        



        assign("PredictedTFBSTrack", PredictedTFBSTrack, .GlobalEnv)
        ########################################################
        ## Re render each time anything changes
        ###########################################################
        # Raw Motif Instances
        RawMotifInstancesTrack <- subset(genomicLocationOfMotifs,
                                         seqnames== input$chrM & start > input$fromM & end < input$toM)%>%AnnotationTrack(.,
                                                                                                                          genome = "hg19",
                                                                                                                          stacking = "dense",
                                                                                                                          col.line="black",
                                                                                                                          name="All Motif Instances")
        displayPars(RawMotifInstancesTrack) <-  list(fill = "grey")
        assign("RawMotifInstancesTrack", RawMotifInstancesTrack, .GlobalEnv)


        # Chromosome For Predicted Motifs
        chromatinStatesTrack<-chromHMMTrackGenerator(gen="hg19",
                                                     chr= input$chrM,
                                                     from  = input$fromM,
                                                     to = input$toM,
                                                     bedFile = chromatinState,
                                                     featureDisplay = "all",
                                                     colorcase='roadmap15',
                                                     name = NameOfCell)

        assign("chromatinStatesTrack", chromatinStatesTrack, .GlobalEnv)



        plotTracks(trackList =c(humanIdeogramTrack,
                                gHumanTrack,
                                IntearctionTrack,
                                EnhancersHumanChromosomeSpecific,
                                promotertrackChromosomeSpecific,
                                PredictedTFBSTrack,
                                RawMotifInstancesTrack,
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
                   fontcolor.title = "black",
                   legend=TRUE)


        ### Putting it last so that if something goes wrong, it'll re run the code when you click refresh
        chrM <- input$chrM
        assign("chrM", chrM, .GlobalEnv)



      } else if(!identical(matrixForMatching,GenomeBrowserTFMatrix)) {



        GenomeBrowserTFMatrix <- matrixForMatching

        assign("GenomeBrowserTFMatrix", GenomeBrowserTFMatrix, .GlobalEnv)


        ###############################################################################################
        ###Genome browser part for when you change chromosomes
        ##############################################################################################

        chr <- input$chrM

        # Ideogram Track
        humanIdeogramTrack <- IdeogramTrack(chromosome = chr,
                                            genome="hg19",
                                            name= "Ideogram")

        #Genome Axis Track for apprxoimate location
        gHumanTrack <- GenomeAxisTrack(name= "Genomic Axis Track")

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
                                                                   seqnames==input$chrM)%>%AnnotationTrack(.,
                                                                                                           name= "PromoterTrack",
                                                                                                           genome="hg19")
        
        displayPars(promotertrackChromosomeSpecific) <- 
          list(fill = "red")
        geneTrackChromosomeSpecific <- knownGenes
        
        
        
        EnhancersHumanChromosomeSpecific <- reduce(c(reduce(permissiveEnhancer), reduce(EnhancersWithGeneTargetsGrange)))%>%
          subset(. ,
                 seqnames==input$chrM)%>%AnnotationTrack(.,
                                                         name = "Enhancers",
                                                         genome = "hg19")
        displayPars(EnhancersHumanChromosomeSpecific) <- 
          list(fill = "#911eb4")
        
        assign("EnhancersHumanChromosomeSpecific", EnhancersHumanChromosomeSpecific, .GlobalEnv)
        assign("promotertrackChromosomeSpecific", promotertrackChromosomeSpecific, .GlobalEnv)
        assign("geneTrackChromosomeSpecific", geneTrackChromosomeSpecific, .GlobalEnv)
        
        #Chromosome Specific Predicted Motifs
        PredictedTFBSTrack <- PredictedTFBS%>%subset(seqnames==input$chrM)%>%AnnotationTrack(genome = "hg19",
                                                                                             stacking = "dense",
                                                                                             strand= "*",
                                                                                             col.line="black",
                                                                                             name="Predicted TFBS")
        
        displayPars(PredictedTFBSTrack) <-  list(fill = "Black")
        
        
        
        
        assign("PredictedTFBSTrack", PredictedTFBSTrack, .GlobalEnv)
        ########################################################
        ## Re render each time anything changes
        ###########################################################
        # Raw Motif Instances
        RawMotifInstancesTrack <- subset(genomicLocationOfMotifs,
                                         seqnames== input$chrM & start > input$fromM & end < input$toM)%>%AnnotationTrack(.,
                                                                                                                          genome = "hg19",
                                                                                                                          stacking = "dense",
                                                                                                                          col.line="black",
                                                                                                                          name="All Motif Instances")
        displayPars(RawMotifInstancesTrack) <-  list(fill = "grey")
        assign("RawMotifInstancesTrack", RawMotifInstancesTrack, .GlobalEnv)

        # Chromosome For Predicted Motifs
        chromatinStatesTrack<-chromHMMTrackGenerator(gen="hg19",
                                                     chr= input$chrM,
                                                     from  = input$fromM,
                                                     to = input$toM,
                                                     bedFile = chromatinState,
                                                     featureDisplay = "all",
                                                     colorcase='roadmap15',
                                                     name = NameOfCell)

        assign("chromatinStatesTrack", chromatinStatesTrack, .GlobalEnv)



        plotTracks(trackList =c(humanIdeogramTrack,
                                gHumanTrack,
                                IntearctionTrack,
                                EnhancersHumanChromosomeSpecific,
                                promotertrackChromosomeSpecific,
                                PredictedTFBSTrack,
                                RawMotifInstancesTrack,
                                geneTrackChromosomeSpecific,
                                chromatinStatesTrack),
                   sizes= c(1,1,3,1,1,1,1,3,3),
                   from = input$fromM,
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
                   fontcolor.title = "black",
                   legend=TRUE)

        ### Putting it last so that if something goes wrong, it'll re run the code when you click refresh
        chrM <- input$chrM
        assign("chrM", chrM, .GlobalEnv)

      }else if(!chrM == input$chrM){

        ###############################################################################################
        ###Genome browser part for when you change chromosomes
        ##############################################################################################

        chr <- input$chrM

        # Ideogram Track
        humanIdeogramTrack <- IdeogramTrack(chromosome = chr,
                                            genome="hg19",
                                            name= "Ideogram")

        gHumanTrack <- GenomeAxisTrack(name= "Axis")

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
                                                                   seqnames==input$chrM)%>%AnnotationTrack(.,
                                                                                                           name= "PromoterTrack",
                                                                                                           genome="hg19")
        
        
        displayPars(promotertrackChromosomeSpecific) <- 
          list(fill = "red")
        geneTrackChromosomeSpecific <- knownGenes
        
        
        
        EnhancersHumanChromosomeSpecific <- reduce(c(reduce(permissiveEnhancer), reduce(EnhancersWithGeneTargetsGrange)))%>%
          subset(. ,
                 seqnames==input$chrM)%>%AnnotationTrack(.,
                                                         name = "Enhancers",
                                                         genome = "hg19")
        displayPars(EnhancersHumanChromosomeSpecific) <- 
          list(fill = "#911eb4")
        
        assign("EnhancersHumanChromosomeSpecific", EnhancersHumanChromosomeSpecific, .GlobalEnv)
        assign("promotertrackChromosomeSpecific", promotertrackChromosomeSpecific, .GlobalEnv)
        assign("geneTrackChromosomeSpecific", geneTrackChromosomeSpecific, .GlobalEnv)
        
        #Chromosome Specific Predicted Motifs
        PredictedTFBSTrack <- PredictedTFBS%>%subset(seqnames==input$chrM)%>%AnnotationTrack(genome = "hg19",
                                                                                             stacking = "dense",
                                                                                             strand= "*",
                                                                                             col.line="black",
                                                                                             name="Predicted TFBS")
        
        displayPars(PredictedTFBSTrack) <-  list(fill = "Black")
        
        
        
        
        assign("PredictedTFBSTrack", PredictedTFBSTrack, .GlobalEnv)
        
        
        ########################################################
        ## Re render each time anything changes
        ###########################################################
        # Raw Motif Instances
        RawMotifInstancesTrack <- subset(genomicLocationOfMotifs,
                                         seqnames== input$chrM & start > input$fromM & end < input$toM)%>%AnnotationTrack(.,
                                                                                                                          genome = "hg19",
                                                                                                                          stacking = "dense",
                                                                                                                          col.line="black",
                                                                                                                          name="All Motif Instances")
        displayPars(RawMotifInstancesTrack) <-  list(fill = "grey")
        assign("RawMotifInstancesTrack", RawMotifInstancesTrack, .GlobalEnv)


        # Chromosome For Predicted Motifs
        chromatinStatesTrack<-chromHMMTrackGenerator(gen="hg19",
                                                     chr= input$chrM,
                                                     from  = input$fromM,
                                                     to = input$toM,
                                                     bedFile = chromatinState,
                                                     featureDisplay = "all",
                                                     colorcase='roadmap15',
                                                     name = NameOfCell)


        assign("chromatinStatesTrack", chromatinStatesTrack, .GlobalEnv)



        plotTracks(trackList =c(humanIdeogramTrack,
                                gHumanTrack,
                                IntearctionTrack,
                                EnhancersHumanChromosomeSpecific,
                                promotertrackChromosomeSpecific,
                                PredictedTFBSTrack,
                                RawMotifInstancesTrack,
                                geneTrackChromosomeSpecific,
                                chromatinStatesTrack),
                   sizes= c(1,1,3,1,1,1,1,3,3),
                   from =24500000,
                   to= 25000000,
                   chromosome= "chrX",
                   cex.title = 0.72,
                   rotation.title = 0,
                   showAxis = FALSE,
                   background.title = "white",
                   lwd.title = 2,
                   title.width = 2,
                   cex.main = 5,
                   col = NULL,
                   fontcolor.title = "black",
                   legend=TRUE)
        
        
        chrM <- input$chrM
        assign("chrM", chrM, .GlobalEnv)


      } else {


        ########################################################
        ## Re render each time anything changes
        ###########################################################
        # Raw Motif Instances
        RawMotifInstancesTrack <- subset(genomicLocationOfMotifs,
                                         seqnames== input$chrM & start > input$fromM & end < input$toM)%>%AnnotationTrack(.,
                                                                                                                          genome = "hg19",
                                                                                                                          stacking = "dense",
                                                                                                                          col.line="black",
                                                                                                                          name="All Motif Instances")
        displayPars(RawMotifInstancesTrack) <-  list(fill = "grey")
        assign("RawMotifInstancesTrack", RawMotifInstancesTrack, .GlobalEnv)

        # Chromosome For Predicted Motifs

        chromatinStatesTrack<-chromHMMTrackGenerator(gen="hg19",
                                                     chr= input$chrM,
                                                     from  = input$fromM,
                                                     to = input$toM,
                                                     bedFile = chromatinState,
                                                     featureDisplay = "all",
                                                     colorcase='roadmap15',
                                                     name = NameOfCell)

        assign("chromatinStatesTrack", chromatinStatesTrack, .GlobalEnv)



        plotTracks(trackList =c(humanIdeogramTrack,
                                gHumanTrack,
                                IntearctionTrack,
                                EnhancersHumanChromosomeSpecific,
                                promotertrackChromosomeSpecific,
                                PredictedTFBSTrack,
                                RawMotifInstancesTrack,
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
                   fontcolor.title = "black",
                   legend=TRUE)
      }

    })
  }, height = 625)

  output$DataTablePredictedSites<- renderDataTable({

    RenderDataFrame <- subset(PredictedTFBS, start >= input$fromM & end <= input$toM & seqnames == input$chrM)%>%as.data.frame()
  }
  )








  output$LegendsPlot<- renderImage({

    list(
      src = "www/ChromatinStateLegend.png",
      contentType = "image/jpeg",
      alt = "Human/Epigenomics Road Map Legend",
      width= 200,
      height= 400
    )

  }, deleteFile = FALSE)




  
  
  
  
  
  
  
  ####################################### 
  ### Gene Ontology Setup the biomart Object & Univserser
  #######################################


  ## Biomart Object
     human_mart <- biomaRt::useMart("ENSEMBL_MART_ENSEMBL",
                                     dataset = "hsapiens_gene_ensembl")

      assign("human_mart", human_mart, .GlobalEnv)

      
      ## Get the gene Ids for all Symbols 
            uniservser <- getBM(attributes = c("ucsc", "entrezgene", "external_gene_name"),
                          mart = human_mart)

      assign("uniservser", uniservser, .GlobalEnv)

      
      ##############################################
      # Data Table setup Click if Tf changes
      ##############################################
      
  GeneOntology <- observe({
    
    

    input$GeneOntology
    

    isolate({
      
      

      ## Gene List
      geneListToConvert <- unique(c(as.character(UnbiasedMotifsPredicted$`Promoters With Gene Targets`$`UCSC Transcript ID`),
                       as.character(UnbiasedMotifsPredicted$`Enhancers With Gene Targets`$`UCSC Transcript ID`)))

      assign("geneListToConvert", geneListToConvert, .GlobalEnv)
      
      
      
      
      ## Entrez Ids for uniservers  
      ucscToEntrez <- uniservser [uniservser$ucsc %in%  geneListToConvert,]

      assign("ucscToEntrez", ucscToEntrez, .GlobalEnv)



# Conduct the Gene Ontology
      geneOntology <- goana(de = ucscToEntrez$entrezgene,
                            universe = uniservser$entrezgene,
                            FDR=0.05,
                            species = "Hs"
                            )



      GeneOntologyResultsSorted <- topGO(geneOntology,
                                         number = dim(geneOntology)[1],
                                         truncate.term = NULL)%>%rownames_to_column(var = "GO ID")%>%
        set_colnames(value = c("GO ID",
                               "GO Term",
                               "Ontology Type",
                               "Number of Genes In ID",
                               "Number of DE Genes In Term",
                               "P Value"))


      quickGoHtmlButton <- lapply( GeneOntologyResultsSorted$`GO ID`, function(val) {
        
        paste0('<a href="https://www.ebi.ac.uk/QuickGO/term/', val,'"','target="_blank" class="btn btn-primary">', val, '</a>')
        
      })
      
      GeneOntologyResultsSorted$`GO ID` <- quickGoHtmlButton
      
      
########## Bon Fonerani P Value Adjustment
      
      GeneOntologyResultsSorted$`Adj P Value` <- GeneOntologyResultsSorted$`P Value`*dim(filter(GeneOntologyResultsSorted, `Number of DE Genes In Term`>0 ))[1]
      
      GeneOntologyResultsSorted$`Adj P Value`[GeneOntologyResultsSorted$`Adj P Value`>=1] <-1
      
      
      assign("GeneOntologyResultsSorted", GeneOntologyResultsSorted, .GlobalEnv)

      



    })
  })
      
      


      output$GeneOntologyTableResults <- renderDataTable({
  
  GeneTable <- filter(GeneOntologyResultsSorted, `Number of DE Genes In Term` >= input$GeneOntologyNumeric &
           `Ontology Type` %in% input$OntologyType & `Adj P Value` <= input$GeneOntologyPValueCutOff & `P Value` <= input$GeneOntologyRawPvalue
         )
         
  GeneTable
})


      
      
      
      
      })
