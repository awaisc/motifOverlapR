#
# This is the server logic of a Shiny web application. You can run the
# application by clicking 'Run App' above.
#
# Find out more about building applications with Shiny here:
#
#    http://shiny.rstudio.com/
#
options(shiny.maxRequestSize = 200*1024^2, repos = BiocInstaller::biocinstallRepos())
getOption("repos")
options("stringsAsFactors" = FALSE)
require(shiny)
require(Gviz)
require(GenomicInteractions)
require(magrittr)
require(TxDb.Hsapiens.UCSC.hg19.knownGene)
require(Homo.sapiens)
require(BSgenome.Hsapiens.UCSC.hg19)
require(readxl)
require(dplyr)
require(tidyr)
require(JASPAR2016)
require(TFBSTools)
require(readr)
library(AnnotationHub)
require(seqLogo)
require(readr)
require(limma)
require(biomaRt)
require(tibble)

shinyServer(function(input, output, session) {
  
  
  
  USER <- reactiveValues(Logged = FALSE , session = session$user) 
  
  
  source("Login.R",  local = TRUE)
  

  

  ## Querying the JASPAR data Base for the matrices and names of the TFs

  JASPARR2016Matrices <- getMatrixSet(JASPAR2016, list() )
  
  namedJasparDataBase <- lapply(1:length(JASPARR2016Matrices), function(x){
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

      seqLogo::seqLogo(PWMTextTranscriptionFactor)


    } else if(input$TypeInSequence==TRUE){


      TranscriptionFactor <- round(PWM(input$CustomDNASequence)*dim((PWM(input$CustomDNASequence)))[2]) %T>%seqLogo::seqLogo()

      assign("matrixForMatching", TranscriptionFactor, .GlobalEnv)


    } else {

      ## Rewrite this because some of these matrices dont have even numbers of nucelotides for eahc position #WTF
      # OKay so the issue appears to be with repeat numbers Eg 0.333333333333333333333. NOt sure how to deal with this.


      (JASPARR2016Matrices[[paste0(input$TranscriptionFactorPWM)]]@profileMatrix %>%
         apply(MARGIN = 2, function(x){x / sum(x)}))%>%seqLogo::seqLogo()


      JasparMatrixForMatching <- JASPARR2016Matrices[[paste0(input$TranscriptionFactorPWM)]]@profileMatrix
      assign("matrixForMatching", JasparMatrixForMatching, .GlobalEnv)
    }
  })

  ###########################################################
  #### Data Inputs
  ###########################################################
geneGrange <- readr::read_delim("../DataFiles/Gene Tracks/Human/hg19GeneLevelFile.txt", delim = "\t") %>%
    makeGRangesFromDataFrame(
      keep.extra.columns=TRUE,
      ignore.strand=FALSE,
      seqinfo=NULL,
      seqnames.field="chrom",
      start.field= "start",
      end.field="end",
      strand.field="strand")
  

  #########################################################
  # Initiating the annotation hub object to get the  Chromatin state GRANGE off of it
  ##############################################################
  ah <- AnnotationHub(localHub = F)

  assign("ah", ah, .GlobalEnv)

  EpigenomicsConverter <- cbind.data.frame("H1 Cell Line"="E003",
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
                                            "Primary T helper naive cells from peripheral blood "="E039",
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
                                            "Stomach Mucosa"="E110")%>%
    t()%>%
    set_colnames("Epigenomic Road Map")%>%
    as.data.frame()%>%
    rownames_to_column(var= "Name Of Cell")

  output$UIOfMotifOverlapRWithPW <- renderUI({
    #Make sure they are logged in!
    if(USER$Logged == TRUE){
      
      tabsetPanel(
        
        
        tabPanel("Computational Prediction Componment", 
                 tabsetPanel(
                   tabPanel("Filter TFBS",
                            
                            box(
                              h4("Motif/TFBS Input"),
                              p("Please only select one at a time"),
                              checkboxInput("TypeInSequence", "Type In DNA Sequence", value =  FALSE),
                              conditionalPanel(
                                condition = "input.TypeInSequence == true && input.CustomTFMatrix == false && input.JasparDataBase == false && input.CustomPredictedSites == false" ,
                                textInput("CustomDNASequence", "Type in DNA Sequence", value = "TAATTA" )
                              ),
                              
                              checkboxInput("CustomTFMatrix", "Upload Custom TF Position Weight Matrix", value =  FALSE),
                              conditionalPanel(
                                condition = "input.TypeInSequence == false && input.CustomTFMatrix == true && input.JasparDataBase == false && input.CustomPredictedSites == false" ,
                                fileInput("JasparCustomMatrix", "Upload Custom Jaspar Formated Position Weight Matrix", multiple = FALSE)
                              ),
                              
                              checkboxInput("JasparDataBase", "Select Position Weight Matrix from Jaspar Database",value =  FALSE),
                              conditionalPanel(
                                condition = "input.TypeInSequence == false && input.CustomTFMatrix == false && input.JasparDataBase == true && input.CustomPredictedSites == false" ,
                                fluidRow(
                                  column(4,selectizeInput(inputId= "TranscriptionFactorPWM", 
                                                          label = "Transcription Factor", 
                                                          choices = unlist(namedJasparDataBase), 
                                                          selected = "FOXP3")),
                                  column(4, numericInput(inputId= "MatchPercentage", label = "% Match of motif", value = 90, min = 0, max= 100)))),
                              
                              
                              checkboxInput("CustomPredictedSites", label = "Upload pre-identifed TFBS", value = TRUE),
                              conditionalPanel(
                                condition = "input.TypeInSequence == false && input.CustomTFMatrix == false && input.JasparDataBase == false && input.CustomPredictedSites == true" ,
                                p("Please use the Upload custom TFBS tab!")
                              ),
                              
                              
                              h4("Cis Regulatory module Options"),
                              checkboxInput("CustomPromoterLengths", label = "Advanced Cis Regulatory Module Options", value = FALSE),
                              conditionalPanel(
                                condition = "input.CustomPromoterLengths == true",
                                fluidRow(
                                  #Promoter Option
                                  column(2, numericInput(inputId= "PromoterStart", label = "Promoter Start/
                                                         Upstream of TSS", value = 5000, min = 0, width = 100)),
                                  
                                  column(2,numericInput(inputId= "PromoterFinish", label = "Promoter End/
                                                        DownStream of TSS", value = 5000, min = 0, width= 100)),
                                  # Enhancer Options
                                  column(6, 
                                         sliderInput("correlationCutOff", "Minimum enhancer-promoter correlation score", 
                                                     min = 0, max= 1, value = 0.7))
                                  )
                                  ),
                              
                              
                              
                              h4("Epigenomic Environment"),
                              selectizeInput("CellTypeToPredict", label = "Cell Type/ Epigenomic Environment",
                                             choices = set_names(as.character(EpigenomicsConverter$`Epigenomic Road Map`),
                                                                 as.character(EpigenomicsConverter$`Name Of Cell`))%>%
                                               c(),
                                             multiple= FALSE, selected= "E033"),
                              
                              # Transcriptomic List
                              h4("Transcriptomic Input"),
                              checkboxInput("DifferentialExpressedGenes", label = "Differentially expressed genes", value = FALSE),
                              conditionalPanel(
                                condition = "input.DifferentialExpressedGenes == true",
                                p("Please upload your diferentially expressed
                                  gene symbol list in the: Upload DE Genes tab")),
                              
                              # Action button
                              actionButton("ComputeTranscriptionFactorSites", "Computationally Predict Sites"), 
                              width = 4, height = 900),
                            
                            box(title = "Position Weight Matrix", 
                                fluidRow(withSpinner(plotOutput("VisualizeTFMotif", width="90%", height = 600),
                                                     type = getOption("spinner.type", default = 3),
                                                     color = getOption("spinner.color", default = "#0275D8"),
                                                     color.background = getOption("spinner.color.background", default = "#FFFFFF"))),
                                br(),
                                h4("Downloads"),
                                fluidRow(
                                  #Download buttons
                                  column( 3 ,
                                          downloadButton('rawMotifPositions', 'Download all
                                                         TFBS'),
                                          
                                          downloadButton('RegualtoryModuleMotifs', 'Download TFBS 
                                                         in CRMs'),
                                          
                                          downloadButton('returnObjectUnbaised', 'Download TFBS 
                                                         in active regions'), 
                                          
                                          downloadButton('TranscriptomicFiltered', 'Download the 
                                                         Transcriptomic TFBS')
                                          ),
                                  
                                  column( 8 ,
                                          p("Download the unfiltered motifs or uploaded TFBS"), 
                                          
                                          p("Download the motifs/TFBS in CRMs"),
                                          
                                          p("Download the motifs/TFBS filtered for CRM and active chromatin states"), 
                                          
                                          p("Download the motifs/TFBS in CRM and active chromatin states regulating DE genes"))
                                  
                                  
                                          ),
                                collapsible = TRUE, width= 8, height = 900
                                  ),
                            
                            
                            
                            box(title = "Predicted TFBS ", width= 12, collapsible = TRUE,
                                withSpinner(dataTableOutput("EnhancerTable"),
                                            type = getOption("spinner.type", default = 3),
                                            color = getOption("spinner.color", default = "#0275D8"),
                                            color.background = getOption("spinner.color.background", default = "#FFFFFF")))),
                   
                   
                   tabPanel("Upload custom TFBS",
                            box( title = "Import spreadsheet options",
                                 
                                 # Input: Select a file ----
                                 fileInput("CustomPredictedSitesGenomicSites", "Upload bed/csv file containing TFBS",
                                           multiple = FALSE,
                                           accept = c("text/csv",
                                                      "text/comma-separated-values,text/plain",
                                                      ".csv")),
                                 p("The genomic coordinates must be aligned to the hg19 genome"),
                                 
                                 
                                 # Input: Checkbox if file has header ----
                                 checkboxInput("header", "Header", TRUE),
                                 
                                 
                                 # Input: Checkbox if user wants to keep extra columns ----
                                 checkboxInput("customKeepExtraColumns", "Keep columns beyond genomic coordinates", FALSE),
                                 
                                 # Input: Select separator ----
                                 selectizeInput("sep", "Separator",
                                                choices = c(Comma = ",",
                                                            Semicolon = ";",
                                                            Tab = "\t"),
                                                selected = ","),
                                 
                                 # Input: Select quotes ----
                                 selectizeInput("quote", "Quote",
                                                choices = c(None = "",
                                                            "Double Quote" = '"',
                                                            "Single Quote" = "'"),
                                                selected = '"'),
                                 numericInput("skip", "Skip", value = 0), width = 4, height = 500
                            ),
                            
                            
                            # Output: Data file ----
                            box(
                              withSpinner( dataTableOutput("PreviewOfTheCustomTFBS"),
                                           type = getOption("spinner.type", default = 3),
                                           color = getOption("spinner.color", default = "#0275D8"),
                                           color.background = getOption("spinner.color.background", default = "#FFFFFF")),
                              width = 8
                              
                              
                              
                            ),
                            box(
                              title = "Please input the names of the columns with the corresponding data",
                              withSpinner(uiOutput("NamesChromosome"),
                                          type = getOption("spinner.type", default = 3),
                                          color = getOption("spinner.color", default = "#0275D8"),
                                          color.background = getOption("spinner.color.background", default = "#FFFFFF")),
                              withSpinner(uiOutput("NamesStart"),
                                          type = getOption("spinner.type", default = 3),
                                          color = getOption("spinner.color", default = "#0275D8"),
                                          color.background = getOption("spinner.color.background", default = "#FFFFFF")),
                              withSpinner(uiOutput("NamesEnd"),
                                          type = getOption("spinner.type", default = 3),
                                          color = getOption("spinner.color", default = "#0275D8"),
                                          color.background = getOption("spinner.color.background", default = "#FFFFFF")),
                              withSpinner(uiOutput("NameStrand"),
                                          type = getOption("spinner.type", default = 3),
                                          color = getOption("spinner.color", default = "#0275D8"),
                                          color.background = getOption("spinner.color.background", default = "#FFFFFF")),
                              
                              
                              
                              width = 12)
                            
                            
                            
                   ),
                   tabPanel("Upload your DE Genes",
                            box( title = "Import spreadsheet options",
                                 
                                 # Input: Select a file ----
                                 fileInput("differenitallyExpressedGenesList",
                                           "Upload Differentially Expressed Gene Symbol list",
                                           multiple = FALSE,
                                           accept = c("text/csv",
                                                      "text/comma-separated-values,text/plain",
                                                      ".csv")),
                                 
                                 
                                 
                                 # Input: Checkbox if file has header ----
                                 checkboxInput("headerDE", "Header", TRUE),
                                 
                                 # Input: Select separator ----
                                 selectizeInput("sepDE", "Separator",
                                                choices = c(Comma = ",",
                                                            Semicolon = ";",
                                                            Tab = "\t"),
                                                selected = ","),
                                 
                                 # Input: Select quotes ----
                                 selectizeInput("quoteDE", "Quote",
                                                choices = c(None = "",
                                                            "Double Quote" = '"',
                                                            "Single Quote" = "'"),
                                                selected = '"'),
                                 numericInput("skipDE", "Skip", value = 0), width = 4, height = 500
                            ),
                            
                            # Main panel for displaying outputs ----
                            box(
                              
                              # Output: Data file ----
                              dataTableOutput("PreviewDEGenes"), width = 8
                              
                            ),
                            
                            box(
                              title = "Please input the names of the columns with the corresponding data",
                              
                              withSpinner(uiOutput("DEGeneColumn"),
                                          type = getOption("spinner.type", default = 3),
                                          color = getOption("spinner.color", default = "#0275D8"),
                                          color.background = getOption("spinner.color.background", default = "#FFFFFF")),
                              # 
                              # textInput("customDeLfcColumn", label = "Name of the LFC Column"),
                              
                              
                              width = 12)
                   ))
),

tabPanel("Genome Browser and Data Table",
         box(title= "Genome Browser Inputs",
             
             checkboxInput("GeneSearching", "Search by genes", value = FALSE),
             conditionalPanel("input.GeneSearching==true", 
                              
                              selectizeInput("gene", 
                                             label = "Search by gene symbols", 
                                             choices = set_names(geneGrange$gene, geneGrange$gene),
                                             selected = "FOXP3"),
                              numericInput("zoomOutby", label = "# of bases to zoom out around the gene 
                                           Negative values will zoom in", value = 50000)
                              ),
             conditionalPanel("input.GeneSearching==false", 
                              numericInput("fromM", "Starting Base",value = 49006897),
                              
                              numericInput("toM", "Finishing Base", value = 49321288),
                              
                              selectInput("chrM", label = h4("Chromosome"),
                                          
                                          choices = list("Chromosome 1" = "chr1",
                                                         "Chromosome 2" = "chr2",
                                                         "Chromosome 3" = "chr3",
                                                         "Chromosome 4" = "chr4",
                                                         "Chromosome 5" = "chr5",
                                                         "Chromosome 6" = "chr6",
                                                         "Chromosome 7" = "chr7",
                                                         "Chromosome 8" = "chr8",
                                                         "Chromosome 9" = "chr9",
                                                         "Chromosome 10" = "chr10",
                                                         "Chromosome 11" = "chr11",
                                                         "Chromosome 12" = "chr12",
                                                         "Chromosome 13" = "chr13",
                                                         "Chromosome 14" = "chr14",
                                                         "Chromosome 15" = "chr15",
                                                         "Chromosome 16" = "chr16",
                                                         "Chromosome 17" = "chr17",
                                                         "Chromosome 18" = "chr18",
                                                         "Chromosome 19" = "chr19",
                                                         "Chromosome 20" = "chr20",
                                                         "Chromosome 21" = "chr21",
                                                         "Chromosome 22" = "chr22",
                                                         "Chromosome X" = "chrX",
                                                         "Chromosome Y" = "chrY"), selected = "chrX")
                              ),
             actionButton("GenomeBrowserAction", "Refresh Genome Browser"),
             br(),
             
             # activate me later whne u get the server side running
             downloadButton(outputId = 'SaveHighResGVIZImage', label = 'High-Res 
                                  Genome Browser Screen shot'),
             
             width = 3, height = 700),
         box( title = "Genome Browser", width = 9, collapsible = TRUE, height = 700,
              withSpinner(plotOutput("GenomeBrowser"),
                          type = getOption("spinner.type", default = 3),
                          color = getOption("spinner.color", default = "#0275D8"),
                          color.background = getOption("spinner.color.background", default = "#FFFFFF"))),
         
         box(title = "Legend for ChromHMM", collapsible = TRUE, width=3,
             withSpinner(plotOutput("LegendsPlot"),
                         type = getOption("spinner.type", default = 3),
                         color = getOption("spinner.color", default = "#0275D8"),
                         color.background = getOption("spinner.color.background", default = "#FFFFFF"))),
         box( title = "TFBS shown in browser", collapsible = TRUE, width = 9,
              tabsetPanel(
                tabPanel("Multi-omic filtered TFBS",
                         withSpinner(dataTableOutput("DataTablePredictedSites"),
                                     type = getOption("spinner.type", default = 3),
                                     color = getOption("spinner.color", default = "#0275D8"),
                                     color.background = getOption("spinner.color.background", default = "#FFFFFF"))),
                tabPanel("Raw/ Non filtered TFBS",
                         withSpinner(dataTableOutput("DataTableRawSites"),
                                     type = getOption("spinner.type", default = 3),
                                     color = getOption("spinner.color", default = "#0275D8"),
                                     color.background = getOption("spinner.color.background", default = "#FFFFFF")))))
         
         ),



tabPanel("Gene Ontology",
         box( 
           
           
           actionButton("GeneOntology", "Re-compute Gene Ontology Results"),
           h4("Click me if the transcription Factor or Cell type changes."),
           selectInput("OntologyType", "Gene Ontology Type", list("Molecular Function" = "MF",
                                                                  "Biological Processes" = "BP",
                                                                  "Cellular Component" = "CC"), 
                       selected = "BP", 
                       multiple = TRUE,
                       selectize = TRUE),
           sliderInput("GeneOntologyPValueCutOff", "Adjusted P value Cut off", min = 0, max= 1, value = 0.05),
           sliderInput("GeneOntologyRawPvalue", "P value Cut off", min = 0, max= 1, value = 0.05),
           numericInput("GeneOntologyNumeric", label = "Minimum Number of DE genes in Term", value = 0),
           
           width= 3, height= 600),
         
         
         box( title = "Gene Ontology Results", width = 9,
              withSpinner(dataTableOutput("GeneOntologyTableResults"),
                          type = getOption("spinner.type", default = 3),
                          color = getOption("spinner.color", default = "#0275D8"),
                          color.background = getOption("spinner.color.background", default = "#FFFFFF"))), height= 600)


)}
    
    
    
  })
  
  ###############################################################################################
  ###Promoter Enhancer assoication data table to assign gene targets to enhancers using CAGE expression
  ###############################################################################################

  #Permissive Enhancer
  permissiveEnhancer <- import.bed("../DataFiles/Enhancers Track/Human/human_permissive_enhancers_phase_1_and_2.bed.gz")

  assign("permissiveEnhancer", permissiveEnhancer, .GlobalEnv)

#   #Read in the table
#   EnhancerPromoterAssoications <- read_delim("../DataFiles/Enhancers Track/Human/hg19_enhancer_promoter_correlations_distances_cell_type.txt.gz",
#                                             "\t", escape_double = FALSE, trim_ws = TRUE)
# 
#   assign("EnhancerPromoterAssoications", EnhancerPromoterAssoications, .GlobalEnv)
# 
# 
# 
#   PromotersAssoicatedWithEnhancers<-separate(EnhancerPromoterAssoications, col=promoter, into= c("promoter", "strand"), sep= ',')%>%
#     separate(., col= promoter, into= c("chr", "start", "end"))%>%makeGRangesFromDataFrame(keep.extra.columns = TRUE)
# 
# 
# 
# 
#   ## Convert Enhancer To Granges
#   EnhancerGrangeWithTargets <- EnhancerPromoterAssoications$enhancer%>%GRanges()
# 
#   # Add everything to the metadata columns
#   mcols(EnhancerGrangeWithTargets) <- PromotersAssoicatedWithEnhancers
# 
# 
#   # Select for significant values only,
#   EnhancerGrangeWithTargetsSiginificant <- subset(EnhancerGrangeWithTargets, `FDR`<0.05)
# 
#   # Generate a vector for targets of enhancers by overlapping promoter regions with promoter GRanges
#   IntersectBetweenOverlappingPromoterRanges <- findOverlaps(EnhancerGrangeWithTargetsSiginificant$X, promoterTracks)
# 
#   # Subset for Grange for promoter tarcks by the hits from the findOverlaps vector above (this will also order it)
#   GenesRegulatedByEnhancers <- as.data.frame(mcols(promoterTracks))[IntersectBetweenOverlappingPromoterRanges%>%subjectHits(),]
# 
#   # Subset the EnhancerGrange by the same vector (except the vector for this Grange) hence ordering them both in the same order
#   EnhancersWithGeneTargetsGrange <- EnhancerGrangeWithTargetsSiginificant[IntersectBetweenOverlappingPromoterRanges%>%queryHits()]
# 
#   ## Combine them
#   mcols(EnhancersWithGeneTargetsGrange) <- cbind.data.frame(mcols(EnhancersWithGeneTargetsGrange),GenesRegulatedByEnhancers)
# 
#   ## Removing the redudent data! (I had to use this due to dplyr not enjoying columns having the same name)
#   mcols(EnhancersWithGeneTargetsGrange) <- mcols(EnhancersWithGeneTargetsGrange)[12:20]
# 
# saveRDS(EnhancersWithGeneTargetsGrange, "../DataFiles/Interactions/Human/EnhancerPromoterInteractionDataFrame.rds")


EnhancersWithGeneTargetsGrange <- readRDS("../DataFiles/Interactions/Human/EnhancerPromoterInteractionDataFrame.rds")

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
    # edited coMet function to use the Genomic Ranges instead of a dataframe
    desiredRegion <- subset(bedFile, end > from &
                              start < to & seqnames == chr)
    #Adapted to use the GRange instead of the dataframe and name the track according to the epigenomic environment number assigned
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

  genomicRangesBlackListedNames <- c("seqnames", "ranges", "strand", "seqlevels", "seqlengths", "isCircular", "start", "end", "width", "element")
  assign("genomicRangesBlackListedNames", genomicRangesBlackListedNames, .GlobalEnv)
  
  ###############################################
  ###Interaction Track Based on the CAGE Expression Enhancer-promoter assoication
  #################################################

  CageExpressionGenomicIntearctions <- readRDS(file = "../DataFiles/Interactions/Human/EnhancerPromoterAssoicationRObject")

  IntearctionTrack <- CageExpressionGenomicIntearctions%>%InteractionTrack(name = "Fantom5 Enhancer-
                                                                           Promoter Assoication
                                                                           Track")
  
  
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


      #Data Pipeline
      ####### Custom Data import so users can adjust the promoter regions as they please
      
      #########Importing Bed File for the promoter regions with gene symbol names
      promoterTracks <- geneGrange%>%
        promoters(upstream = input$PromoterStart,
                  downstream = input$PromoterFinish )
      
      
      assign("promoterTracks", promoterTracks, .GlobalEnv)
      
      ###### If the user uploads a custom set of TFBS
  if(input$CustomPredictedSites==TRUE) {
    
    if(input$CustomPredictedSitesGenomicSites$datapath %>% is.null() ){
        
       # Example DataFile!!
      CustomPredictedSitesPath <- "../ExampleFiles/ChIP_Data_Hg19.csv"
        
        GenomicPositions <- readr::read_delim(file = CustomPredictedSitesPath,
                                                      ",", escape_double = FALSE,
                                             col_names = T, 
                                                      trim_ws = TRUE)%>%
          dplyr::select(c("seqnames","start", "end"))%>% # Remove the irrelevant columns just for the sake of the cleaning up the final datatable
          makeGRangesFromDataFrame(.,
                                   keep.extra.columns=FALSE, # Keep extra columns for additional info but there wont be colnames
                                   ignore.strand=TRUE,
                                   seqinfo=NULL,
                                   seqnames.field= "seqnames",
                                   start.field="start",
                                   end.field= "end",
                                   starts.in.df.are.0based=FALSE)%>%
          sort()
        
        seqlevelsStyle(GenomicPositions) <- "UCSC"
        
        # Assigning to motif positions
        assign("genomicLocationOfMotifs", GenomicPositions, .GlobalEnv)
        
    } else {
      
      # Extract the file path for the uploaded object
      GenomicPositionsdf<- read_delim( input$CustomPredictedSitesGenomicSites$datapath,
                                    col_names = input$header,
                                    delim = input$sep,
                                    quote = input$quote,
                                    skip = input$skip)
      
      # Check Colnames
      # Check Strand
      if(input$customStrand != ""){
        if(! input$customStrand %in% colnames(GenomicPositionsdf)){
          print("Strand column not found in your data frame, check your custom column names!")
          stop()
        }
      }
      
      #Check Colnames 
      if(input$customKeepExtraColumns){
        if(input$customStrand!=""){
          if( isTRUE( 
            colnames(GenomicPositionsdf) %in% 
                      genomicRangesBlackListedNames[!genomicRangesBlackListedNames %in% 
                                                                                      c(input$customChromosome, 
                                                                                        input$customStart,
                                                                                        input$customEnd,
                                                                                        input$customStrand)]
            )
         ){ print("Error, cannot have columns of the following names in your metadata columns!,
              ranges", "seqlevels", "seqlengths", "isCircular", "width", "element")
            stop() }
        } else {
          if( isTRUE( 
            colnames(GenomicPositionsdf) %in% 
            genomicRangesBlackListedNames[!genomicRangesBlackListedNames %in% 
                                          c(input$customChromosome, 
                                            input$customStart,
                                            input$customEnd)]
          )
          ){ print("Error, cannot have columns of the following names in your metadata columns:
              ranges, seqlevels, seqlengths, isCircular, width, element")
            stop() }
            
          } 
      }
      
      GenomicPositionsdf <-
        subset(GenomicPositionsdf,
        # Remove any rows with missing values in the chrom, strand or ranges columns
        (dplyr::select(GenomicPositionsdf, 
                       c(input$customChromosome,input$customStartRange, input$customEndRange ,input$customStartRange))%>%
           is.na(.)%>%
           apply(. ,MARGIN = 1, FUN = sum)==0) 
        )
      
      # Generate a new column called width
      GenomicPositionsdf$width <- subtract(GenomicPositionsdf[input$customEndRange],
                                           GenomicPositionsdf[input$customStartRange])[[1]] # Vectorizing
    
      # Ensure all ranges are positive ranges and swap the start and stops for negative ranges
      #Filter positive ranges and put them aside
      PositiveRanges <-   filter(GenomicPositionsdf, width>=0)
      
      #Extract negative ranges
      NegativeRanges <- filter(GenomicPositionsdf, width < 0)
      #Save the start column of negative ranges
      NegativeStartX2 <- dplyr::select(NegativeRanges,(input$customStartRange))[[1]]
      
      
      #Replace start column of Negative range with end column
      NegativeRanges[input$customStartRange] <- NegativeRanges[input$customEndRange]
      #Replace end column of negative range with start column thus resolving this issue
      NegativeRanges[input$customEndRange] <- NegativeStartX2
      
      if( input$customStrand != ""){
      # Generate the GRanges
      GenomicPositions <- bind_rows(
        PositiveRanges,
        NegativeRanges)%>%
        dplyr::select( -c("width"))%>%
        dplyr::select(c(   input$customChromosome, # Format based on UCSC bedfile names
                           input$customStartRange,
                           input$customEndRange,
                           input$customStrand))%>%# Remove the irreleant columns just for the sake of the data
        makeGRangesFromDataFrame(.,
                                 keep.extra.columns=input$customKeepExtraColumns, # Keep extra columns for additional info but there wont be colnames as bed files dont have col names
                                 ignore.strand=FALSE,
                                 seqinfo=NULL,
                                 seqnames.field= input$customChromosome, # Format based on UCSC bedfile names
                                 start.field=input$customStartRange,
                                 end.field= input$customEndRange,
                                 strand.field = input$customStrand,
                                 starts.in.df.are.0based=FALSE)
      
      seqlevelsStyle(GenomicPositions) <- "UCSC"
      
      
      # Assigning to motif positions
      assign("genomicLocationOfMotifs", GenomicPositions, .GlobalEnv)
      
      } else {
        # Generate the GRanges
        GenomicPositions <- bind_rows(
          PositiveRanges,
          NegativeRanges)%>%
          dplyr::select( -c("width"))%>%
          dplyr::select(c(   input$customChromosome, # Format based on UCSC bedfile names
                              input$customStartRange,
                             input$customEndRange))%>%#Remove the irreleant columns just for the sake of the data
          makeGRangesFromDataFrame(.,
                                   keep.extra.columns= input$customKeepExtraColumns, # Keep extra columns for additional info but there wont be colnames as bed files dont have col names
                                   ignore.strand=TRUE,
                                   seqinfo=NULL,
                                   seqnames.field= input$customChromosome, # Format based on UCSC bedfile names
                                   start.field=input$customStartRange,
                                   end.field= input$customEndRange,
                                   starts.in.df.are.0based=FALSE)
        
        seqlevelsStyle(GenomicPositions) <- "UCSC"
        
        # Assigning to motif positions
        assign("genomicLocationOfMotifs", GenomicPositions, .GlobalEnv)
       
      }
    }
        
      }  else {
        
      ### Input Transcription Factor Matrix from the seqlogo function

      assign("TranscriptionFactorPWM", matrixForMatching, .GlobalEnv )


      genomicLocationOfMotifs <- matchPWM(matrixForMatching,
                                        BSgenome.Hsapiens.UCSC.hg19,
                                        paste0(input$MatchPercentage,'%'))

      assign("genomicLocationOfMotifs", genomicLocationOfMotifs, .GlobalEnv)
  
}


      EnhancersWithGeneTargetsGrange <- subset(EnhancersWithGeneTargetsGrange, 
                                               correlation >= input$correlationCutOff )
      
      ### Identify which of these motifs are located in enhancer regions
      MotifsInEnhancers <- subsetByOverlaps(genomicLocationOfMotifs ,
                                            EnhancersWithGeneTargetsGrange)

      mcols(MotifsInEnhancers) <- cbind.data.frame(mcols(MotifsInEnhancers), "CRM"="Enhancer")
      
      MotifsInPermissiveEnhancers <- subsetByOverlaps(genomicLocationOfMotifs,
                                                      permissiveEnhancer)
      mcols(MotifsInPermissiveEnhancers) <- cbind.data.frame(mcols(MotifsInPermissiveEnhancers), "CRM"="Enhancer")

      # Used to identify which regions of the conserved bigwig file to import
      MotifsInPromoters <- subsetByOverlaps(genomicLocationOfMotifs,
                                            promoterTracks)%>%reduce()
      



        # If we are not looking for conserved promoter regions
        MotifsInConservedPromoterRegions <- subsetByOverlaps(genomicLocationOfMotifs, promoterTracks)
        mcols(MotifsInConservedPromoterRegions) <- cbind.data.frame(mcols(MotifsInConservedPromoterRegions), "CRM"="Promoter")


      #### Combing the promoter and enhancer motifs into a single Grange List object

      MotifsInPromotersAndEnhancers <- list("Promoters" = MotifsInConservedPromoterRegions,
                                            "Enhancers" = MotifsInEnhancers,
                                            "Permissive Enhancers" = MotifsInPermissiveEnhancers)

      assign("MotifsInPromotersAndEnhancers", MotifsInPromotersAndEnhancers, .GlobalEnv)





      ################################################
      #### Cell Type specific epigneonic Analysis
      ##################################################


      ### Download the Chromatin state file From the annotation hub

      chromatinState <- query(ah, c(paste0(input$CellTypeToPredict,"_15_coreMarks_mnemonics"), "EpigenomeRoadMap") )[[1]]


      #Assign the bedfile to the global environment for analysis later on downstream

      assign("chromatinState", chromatinState, .GlobalEnv)


      ## Name for the chromHMM Track Down Stream
      Test <- input$CellTypeToPredict  
      assign("Test", Test, .GlobalEnv)
      NameOfCell <- subset(EpigenomicsConverter, `Epigenomic Road Map`== Test)$`Name Of Cell`

      assign("NameOfCell", NameOfCell, .GlobalEnv)


      ## Subsetting the Chromatin states for active states to identify motifs in these regions


      ActiveChromatinRegions <- chromatinState[chromatinState$abbr %in% c("10_TssBiv",
                                                                          "7_Enh",
                                                                          "1_TssA",
                                                                          "11_BivFlnk",
                                                                          "2_TssAFlnk",
                                                                          "5_TxWk",
                                                                          "4_Tx",
                                                                          "6_EnhG",
                                                                          "12_EnhBiv",
                                                                          "3_TxFlnk"),]

      #Overlap active regions with motifs in CRMs giving us a completely unbiased set of results
      UnbiasedPredictedMotifs <- lapply(MotifsInPromotersAndEnhancers,function(x){


        OverlapBetweenRangesObject <- findOverlaps( x, ActiveChromatinRegions)

        EpigenomicLocation <- ActiveChromatinRegions$name[OverlapBetweenRangesObject%>%subjectHits()]

        OrderingQueryGRange <- x [OverlapBetweenRangesObject%>%queryHits()]

        mcols(OrderingQueryGRange) <- cbind.data.frame(mcols(OrderingQueryGRange),
                                                       "Epigenomic Location" = EpigenomicLocation)

        return(OrderingQueryGRange)

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
                                                        "Genes Regulated" = promoterTracks[OverlappingRangeOfMOtifsInPromoters%>%queryHits()]$gene)


      ## Now lets get targets to genes predicted to be regulated by TFBS in enhancers
      
      ## Enhancers with motifs
      #Get Overlapping ranges
      OverlappingRangeEnhancersMotifs <- findOverlaps(EnhancersWithGeneTargetsGrange, UnbiasedPredictedMotifs$Enhancers)
      
      #Get gene targets vector in the order of Ovrelapping Ranges above
      EnhancerTargets <- EnhancersWithGeneTargetsGrange$hg19.kgXref.geneSymbol[OverlappingRangeEnhancersMotifs%>%queryHits()]
      
      #Get the enhancer motifs into the same order
      MotifsInEnhancers <- UnbiasedPredictedMotifs$Enhancers[OverlappingRangeEnhancersMotifs%>%subjectHits()]
      
      #Paste the genes regulated into the metadata column of the enhancers
      mcols(MotifsInEnhancers) <- cbind.data.frame(mcols(MotifsInEnhancers), 
                                                   "Genes Regulated" = EnhancerTargets)
      
      #Enhancers TFBS without a target (therefore not in the previous enhancer dataframe)
      
      EnhancerMotifsNotTargeted <- UnbiasedPredictedMotifs$`Permissive Enhancers`[!UnbiasedPredictedMotifs$`Permissive Enhancers`
                                                                                  %in% MotifsInEnhancers ]
      
      #Paste NA for these enhancers.
      mcols(EnhancerMotifsNotTargeted) <- cbind.data.frame(mcols(EnhancerMotifsNotTargeted),
                                                           "Genes Regulated" = "NA")
      
      #Combine everything into a singe list!
      UnbiasedMotifsPredicted <- GenomicRangesList("Promoters With Gene Targets" = unbiasedPromoterMotifs,
                                   "Enhancers With Gene Targets" = GRangesList(MotifsInEnhancers,
                                                                               EnhancerMotifsNotTargeted)%>%
                                     unlist())

      assign("UnbiasedMotifsPredicted", UnbiasedMotifsPredicted, .GlobalEnv)


      #######################################################
      ## Left Joining with differentially expressed gene list
      #########################################################

      if(input$DifferentialExpressedGenes == TRUE){


        differenitallyExpressedGenesList <- read_delim(input$differenitallyExpressedGenesList$datapath,
                           col_names = input$headerDE,
                           delim = input$sepDE,
                           quote = input$quoteDE, 
                           skip = input$skipDE)
        

        assign("differenitallyExpressedGenesList", differenitallyExpressedGenesList, .GlobalEnv)

        ## Genes who showed differenital expression with an enhancer that was correlated with its expression
        enhancerTargetsOfTF <- subset(UnbiasedMotifsPredicted$`Enhancers With Gene Targets`,
                                      `Genes Regulated` %in% differenitallyExpressedGenesList[[input$customDEGeneColumn]])


        ##Genes who showed differenital expression with promoter targets
        promoterTargetsOfTF <- subset(UnbiasedMotifsPredicted$`Promoters With Gene Targets`,
                                      `Genes Regulated` %in% differenitallyExpressedGenesList[[input$customDEGeneColumn]])



        ##Predicted Sites in the Regulatory Elements of these genes


        returnObjectDifferentialSites <- GRangesList("Promoter Predicted Sites" = promoterTargetsOfTF,
                                           "Enhancer Predicted Sites" = enhancerTargetsOfTF)%>%unlist()
        names(returnObjectDifferentialSites) <- NULL

        assign("PredictedTFBS", returnObjectDifferentialSites, .GlobalEnv)
        

      }
      else {

        ## Returning unbiased results without transcriptomic data
        returnObjectUnbaised <- GRangesList(
          "Promoter Predicted Sites"= UnbiasedMotifsPredicted$`Enhancers With Gene Targets`,
          "Enhancer Predicted Sites"= UnbiasedMotifsPredicted$`Promoters With Gene Targets` )%>%unlist()

        names(returnObjectUnbaised) <- NULL

        assign("PredictedTFBS", returnObjectUnbaised, .GlobalEnv)
      }

#Printing the table 
      names(PredictedTFBS) <- NULL
      PredictedTFBSTable <- as.data.frame(PredictedTFBS)
      PredictedTFBSTable
    
      
    })
  }, options = list(scrollX = TRUE) )


 
      ##################################################################
      #### Download Exports
      #################################################################

  output$rawMotifPositions <- downloadHandler('rawMotifGenomicPositions.bed', content = function(file) {
    write.table(as.data.frame(genomicLocationOfMotifs), file  ,
                sep="\t",
                row.names = FALSE,
                col.names = FALSE,
                quote = FALSE,
                append = FALSE)
  })    
  
  output$RegualtoryModuleMotifs <- downloadHandler('CRMMotifs.bed', content = function(file) {
    
    RegulatoryModuleDataFrame <- rbind.data.frame(
      # Only using the Permissive enhancer motifs as the CAGE enhancer-promoter motifs are duplicates of this track and there are not needed
      as.data.frame(MotifsInPromotersAndEnhancers$Promoters)%>%mutate("CRM"= "Promoter"),
      as.data.frame(MotifsInPromotersAndEnhancers$`Permissive Enhancers`)%>%mutate("CRM"= "Enhancer")
    )
    write.table( RegulatoryModuleDataFrame, file  ,
                 sep="\t",
                 row.names = FALSE,
                 col.names = FALSE,
                 quote = FALSE,
                 append = FALSE)
  })
  
  output$returnObjectUnbaised <- downloadHandler('EpigneomicFiltered.bed', content = function(file) {
    write.table( rbind.data.frame(as.data.frame(UnbiasedMotifsPredicted$`Promoters With Gene Targets`),
                                  as.data.frame(UnbiasedMotifsPredicted$`Enhancers With Gene Targets`)), file  ,
                 sep="\t",
                 row.names = FALSE,
                 col.names = FALSE,
                 quote = FALSE,
                 append = FALSE)
  })
  
  
  output$TranscriptomicFiltered <- downloadHandler('TranscriptomicFiltered.bed', content = function(file) {
    
    write.table( as.data.frame(PredictedTFBS), file  ,
                sep="\t",
                row.names = FALSE,
                col.names = FALSE,
                quote = FALSE,
                append = FALSE)
    })
    


  output$ReducedTFBS <- downloadHandler('DeduplicatedFullyFilteredTFBS.bed',
                                        content = function(file) {
    
    write.table( PredictedTFBS%>%reduce()%>%as.data.frame(), file  ,
                 sep="\t",
                 row.names = FALSE,
                 col.names = FALSE,
                 quote = FALSE,
                 append = FALSE)
  })
  






  ###############################################################################################
  ###Genome browser Reactive function
  ##############################################################################################
  
  
  output$GenomeBrowser <- renderPlot({
    
    input$GenomeBrowserAction
    isolate({
      
      if(input$GeneSearching==TRUE){
        tempGeneHolder <- subset(geneGrange, gene == input$gene )%>%as.data.frame()
        
        toM <- sum(tempGeneHolder$end, input$zoomOutby)
        fromM <- subtract(tempGeneHolder$start, input$zoomOutby)
        chrM <- tempGeneHolder$seqnames%>%as.character
        
        assign("chrM", chrM, .GlobalEnv)
        assign("fromM", fromM, .GlobalEnv)
        assign("toM", toM, .GlobalEnv)
      } else {
        toM <- input$toM
        fromM <- input$fromM
        chrM <- input$chrM
        assign("chrM", chrM, .GlobalEnv)
        assign("fromM", fromM, .GlobalEnv)
        assign("toM", toM, .GlobalEnv)
      }
      
      
      
      if(!exists("chr")){
        
          
        #Assign current TF matrix to global 
          assign("GenomeBrowserTFMatrix", matrixForMatching, .GlobalEnv)
        
        
        ###############################################################################################
        ###Re run this code if the chromosome changes
        ##############################################################################################
        
        # Ideogram Track
        humanIdeogramTrack <- IdeogramTrack(chromosome = chrM,
                                            genome="hg19",
                                            name= "Ideogram")
        displayPars(humanIdeogramTrack) <- list(fontsize= 20)
        #Genome Axis Track for apprxoimate location
        gHumanTrack <- GenomeAxisTrack(name= "Genomic Axis Track")
        displayPars(gHumanTrack) <- list(fontsize = 20)
        assign("humanIdeogramTrack", humanIdeogramTrack, .GlobalEnv)
        assign("gHumanTrack", gHumanTrack, .GlobalEnv)
        
        
        ###############################################
        ####Identifying motifs in CRM regions
        ################################################
        
        
        
        # Gene Track with symbols :D
        knownGenes <- GeneRegionTrack(TxDb.Hsapiens.UCSC.hg19.knownGene,
                                      genome="hg19",
                                      chromosome=chrM,
                                      showId=TRUE,
                                      geneSymbol=TRUE,
                                      stacking = "full",
                                      name="UCSC 
                                      Known
                                      Genes")
        # Map the ids between symbol and entrezid
        symbols <- unlist(mapIds(Homo.sapiens, gene(knownGenes),
                                 "SYMBOL", "ENTREZID",
                                 multiVals = "first"))
        #Map the ids
        symbol(knownGenes) <- symbols[gene(knownGenes)]
        geneTrackChromosomeSpecific <- knownGenes
        #Increase font size
        displayPars(geneTrackChromosomeSpecific) <- list(fontsize.group = 20)
        
        assign("knownGenes", knownGenes, .GlobalEnv)
        assign("geneTrackChromosomeSpecific", geneTrackChromosomeSpecific, .GlobalEnv)
        
        
        #Promoter and Enhancer Tracks for each chormosome Track
        promotertrackChromosomeSpecific <- promoterTracks%>%subset(. ,
                                                                   seqnames==chrM)%>%AnnotationTrack(. ,
                                                                                                     name = "Promoters",
                                                                                                     genome = "hg19")
        
        displayPars(promotertrackChromosomeSpecific) <- 
          list(fill = "red")
        
        
        EnhancersHumanChromosomeSpecific <- reduce(GRangesList(reduce(permissiveEnhancer), 
                                                               reduce(EnhancersWithGeneTargetsGrange))%>%unlist())%>%
          subset(. ,
                 seqnames==chrM)%>%AnnotationTrack(.,
                                                   name = "Enhancers",
                                                   genome = "hg19")
        displayPars(EnhancersHumanChromosomeSpecific) <- 
          list(fill = "#911eb4")
        
        assign("EnhancersHumanChromosomeSpecific", EnhancersHumanChromosomeSpecific, .GlobalEnv)
        assign("promotertrackChromosomeSpecific", promotertrackChromosomeSpecific, .GlobalEnv)
        
        #Chromosome Specific Predicted Motifs
        PredictedTFBSTrack <- PredictedTFBS%>%subset(seqnames==chrM)%>%AnnotationTrack(genome = "hg19",
                                                                                       stacking = "dense",
                                                                                       strand= "*",
                                                                                       col.line="black",
                                                                                       name="Multi-omic
                                                                                       filtered TFBS")
        
        displayPars(PredictedTFBSTrack) <-  list(fill = "Black")
        
        
        
        
        assign("PredictedTFBSTrack", PredictedTFBSTrack, .GlobalEnv)
        ########################################################
        ## Re render each time anything changes
        ###########################################################
        # Raw Motif Instances
        RawMotifInstancesTrack <- subset(genomicLocationOfMotifs,
                                         seqnames== chrM & start > fromM & end < toM)%>%AnnotationTrack(.,
                                                                                                        genome = "hg19",
                                                                                                        stacking = "dense",
                                                                                                        col.line="black",
                                                                                                        name="All TFBS")
        displayPars(RawMotifInstancesTrack) <-  list(fill = "grey")
        assign("RawMotifInstancesTrack", RawMotifInstancesTrack, .GlobalEnv)
        
        
        # Chromosome For Predicted Motifs
        chromatinStatesTrack<-chromHMMTrackGenerator(gen="hg19",
                                                     chr= chrM,
                                                     from  = fromM,
                                                     to = toM,
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
                   sizes= c(1,1.5,3,1,1,1,1,3,3),
                   from =fromM,
                   to= toM,
                   chromosome= chrM,
                   cex.title = 1,
                   rotation.title = 0,
                   showAxis = FALSE,
                   background.title = "white",
                   lwd.title = 2,
                   title.width = 3,
                   cex.main = 5,
                   col = NULL,
                   fontcolor.title = "black",
                   legend=TRUE)
        
        
        
        ### Putting it last so that if something goes wrong, it'll re run the code when you click refresh
        assign("chr", chrM, .GlobalEnv)
        
        
        
        
      } else if(!identical(matrixForMatching,GenomeBrowserTFMatrix)) {
        
        
          # assign new tf matrix over the current to the global
          assign("GenomeBrowserTFMatrix", matrixForMatching, .GlobalEnv)
        
        ###############################################################################################
        ###Genome browser part for when you change chromosomes
        ##############################################################################################

        
        # Ideogram Track
        humanIdeogramTrack <- IdeogramTrack(chromosome = chrM,
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
                                      chromosome=chrM,
                                      showId=TRUE,
                                      geneSymbol=TRUE,
                                      name="UCSC
                                      Genes")
        
        symbols <- unlist(mapIds(org.Hs.eg.db, gene(knownGenes),
                                 "SYMBOL", "ENTREZID",
                                 multiVals = "first"))
        
        symbol(knownGenes) <- symbols[gene(knownGenes)]
        
        assign("knownGenes", knownGenes, .GlobalEnv)
        
        
        #Promoter and Enhancer Tracks for each chormosome Track
        promotertrackChromosomeSpecific <- promoterTracks%>%subset(. ,
                                                                   seqnames==chrM)%>%AnnotationTrack(.,
                                                                                                     name= "Promoters",
                                                                                                     genome="hg19")
        
        displayPars(promotertrackChromosomeSpecific) <- 
          list(fill = "red")
        geneTrackChromosomeSpecific <- knownGenes
        
        
        
        EnhancersHumanChromosomeSpecific <- reduce(GRangesList(reduce(permissiveEnhancer), 
                                                               reduce(EnhancersWithGeneTargetsGrange))%>%unlist())%>%
          subset(. ,
                 seqnames==chrM)%>%AnnotationTrack(.,
                                                   name = "Enhancers",
                                                   genome = "hg19")
        displayPars(EnhancersHumanChromosomeSpecific) <- 
          list(fill = "#911eb4")
        
        assign("EnhancersHumanChromosomeSpecific", EnhancersHumanChromosomeSpecific, .GlobalEnv)
        assign("promotertrackChromosomeSpecific", promotertrackChromosomeSpecific, .GlobalEnv)
        assign("geneTrackChromosomeSpecific", geneTrackChromosomeSpecific, .GlobalEnv)
        
        #Chromosome Specific Predicted Motifs
        PredictedTFBSTrack <- PredictedTFBS%>%subset(seqnames==chrM)%>%AnnotationTrack(genome = "hg19",
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
                                         seqnames== chrM & start > fromM & end < toM)%>%AnnotationTrack(.,
                                                                                                        genome = "hg19",
                                                                                                        stacking = "dense",
                                                                                                        col.line="black",
                                                                                                        name="Unfiltered TFBS")
        displayPars(RawMotifInstancesTrack) <-  list(fill = "grey")
        assign("RawMotifInstancesTrack", RawMotifInstancesTrack, .GlobalEnv)
        
        # Chromosome For Predicted Motifs
        chromatinStatesTrack<-chromHMMTrackGenerator(gen="hg19",
                                                     chr= chrM,
                                                     from  = fromM,
                                                     to = toM,
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
                   sizes= c(1,1.5,3,1,1,1,1,3,3),
                   from = fromM,
                   to= toM,
                   chromosome= chrM,
                   cex.title = 1,
                   rotation.title = 0,
                   showAxis = FALSE,
                   background.title = "white",
                   lwd.title = 2,
                   title.width = 3,
                   cex.main = 5,
                   col = NULL,
                   fontcolor.title = "black",
                   legend=TRUE)
        
        
        ### Putting it last so that if something goes wrong, it'll re run the code when you click refresh
        assign("chr", chrM, .GlobalEnv)
        
      }else if(!chr == chrM){
        
        ###############################################################################################
        ###Genome browser part for when you change chromosomes
        ##############################################################################################
        
        # Ideogram Track
        humanIdeogramTrack <- IdeogramTrack(chromosome = chrM,
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
                                      chromosome=chrM,
                                      showId=TRUE,
                                      geneSymbol=TRUE,
                                      name="UCSC
                                      Known
                                      Genes")
        
        symbols <- unlist(mapIds(org.Hs.eg.db, gene(knownGenes),
                                 "SYMBOL", "ENTREZID",
                                 multiVals = "first"))
        
        symbol(knownGenes) <- symbols[gene(knownGenes)]
        
        assign("knownGenes", knownGenes, .GlobalEnv)
        
        
        #Promoter and Enhancer Tracks for each chormosome Track
        promotertrackChromosomeSpecific <- promoterTracks%>%subset(. ,
                                                                   seqnames==chrM)%>%AnnotationTrack(.,
                                                                                                     name= "Promoter",
                                                                                                     genome="hg19")
        
        
        displayPars(promotertrackChromosomeSpecific) <- 
          list(fill = "red")
        geneTrackChromosomeSpecific <- knownGenes
        
        
        
        EnhancersHumanChromosomeSpecific <- reduce(GRangesList(reduce(permissiveEnhancer), 
                                                               reduce(EnhancersWithGeneTargetsGrange))%>%unlist())%>%
          subset(. ,
                 seqnames==chrM)%>%AnnotationTrack(.,
                                                   name = "Enhancers",
                                                   genome = "hg19")
        
        displayPars(EnhancersHumanChromosomeSpecific) <- 
          list(fill = "#911eb4")
        
        assign("EnhancersHumanChromosomeSpecific", EnhancersHumanChromosomeSpecific, .GlobalEnv)
        assign("promotertrackChromosomeSpecific", promotertrackChromosomeSpecific, .GlobalEnv)
        assign("geneTrackChromosomeSpecific", geneTrackChromosomeSpecific, .GlobalEnv)
        
        #Chromosome Specific Predicted Motifs
        PredictedTFBSTrack <- PredictedTFBS%>%subset(seqnames==chrM)%>%AnnotationTrack(genome = "hg19",
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
                                         seqnames== chrM & start > fromM & end < toM)%>%AnnotationTrack(.,
                                                                                                        genome = "hg19",
                                                                                                        stacking = "dense",
                                                                                                        col.line="black",
                                                                                                        name="Unfiltered TFBS")
        displayPars(RawMotifInstancesTrack) <-  list(fill = "grey")
        assign("RawMotifInstancesTrack", RawMotifInstancesTrack, .GlobalEnv)
        
        
        # Chromosome For Predicted Motifs
        chromatinStatesTrack<-chromHMMTrackGenerator(gen="hg19",
                                                     chr= chrM,
                                                     from  = fromM,
                                                     to = toM,
                                                     bedFile = chromatinState,
                                                     featureDisplay = "all",
                                                     colorcase='roadmap15',
                                                     name = NameOfCell)
        
        
        assign("chromatinStatesTrack", chromatinStatesTrack, .GlobalEnv)
        
        
        plotTracks(trackList = c(humanIdeogramTrack,
                                gHumanTrack,
                                IntearctionTrack,
                                EnhancersHumanChromosomeSpecific,
                                promotertrackChromosomeSpecific,
                                PredictedTFBSTrack,
                                RawMotifInstancesTrack,
                                geneTrackChromosomeSpecific,
                                chromatinStatesTrack),
                   sizes= c(1,1.5,3,1,1,1,1,3,3),
                   from = fromM,
                   to = toM,
                   chromosome= chrM,
                   cex.title = 1,
                   rotation.title = 0,
                   showAxis = FALSE,
                   background.title = "white",
                   lwd.title = 2,
                   title.width = 3,
                   cex.main = 5,
                   col = NULL,
                   fontcolor.title = "black",
                   legend=TRUE)
        
        
        ### Putting it last so that if something goes wrong, it'll re run the code when you click refresh
        assign("chr", chrM, .GlobalEnv)
        
        
      } else {
        
        
        ########################################################
        ## Re render each time anything changes
        ###########################################################
        # Raw Motif Instances
        RawMotifInstancesTrack <- subset(genomicLocationOfMotifs,
                                         seqnames== chrM & start > fromM & end < toM)%>%
          AnnotationTrack(.,
                          genome = "hg19",
                          stacking = "dense",
                          col.line="black",
                          name="Unfiltered TFBS")
        displayPars(RawMotifInstancesTrack) <-  list(fill = "grey")
        assign("RawMotifInstancesTrack", RawMotifInstancesTrack, .GlobalEnv)
        
        # Chromosome For Predicted Motifs
        
        chromatinStatesTrack<-chromHMMTrackGenerator(gen="hg19",
                                                     chr= chrM,
                                                     from  = fromM,
                                                     to = toM,
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
                   sizes= c(1,1.5,3,1,1,1,1,3,3),
                   from =fromM,
                   to= toM,
                   chromosome= chrM,
                   cex.title = 1,
                   rotation.title = 0,
                   showAxis = FALSE,
                   background.title = "white",
                   lwd.title = 2,
                   title.width = 3,
                   cex.main = 5,
                   col = NULL,
                   fontcolor.title = "black",
                   legend=TRUE)
        
      }
    })
  }, height = 625 )
  
  
  
  
  ###############################################################################################
  ###Table below Genome browser for motif Information
  ##############################################################################################
  
  output$DataTablePredictedSites<- renderDataTable({
    
    if(input$GeneSearching==TRUE){
      tempGeneHolder <- subset(geneGrange, gene == input$gene )%>%as.data.frame()
      
      toM <- sum(tempGeneHolder$end, input$zoomOutby)
      fromM <- subtract(tempGeneHolder$start, input$zoomOutby)
      chrM <- tempGeneHolder$seqnames%>%as.character
      
      assign("chrM", chrM, .GlobalEnv)
      assign("fromM", fromM, .GlobalEnv)
      assign("toM", toM, .GlobalEnv)
    } else {
      toM <- input$toM
      fromM <- input$fromM
      chrM <- input$chrM
      assign("chrM", chrM, .GlobalEnv)
      assign("fromM", fromM, .GlobalEnv)
      assign("toM", toM, .GlobalEnv)
    }
    
    RenderDataFrame <- subset(PredictedTFBS, start >= fromM & end <= toM & seqnames == chrM)%>%as.data.frame()
  }, options = list(scrollX = TRUE)  )
  
  output$DataTableRawSites<- renderDataTable({
    
    if(input$GeneSearching==TRUE){
      tempGeneHolder <- subset(geneGrange, gene == input$gene )%>%as.data.frame()
      
      toM <- sum(tempGeneHolder$end, input$zoomOutby)
      fromM <- subtract(tempGeneHolder$start, input$zoomOutby)
      chrM <- tempGeneHolder$seqnames%>%as.character
      
      assign("chrM", chrM, .GlobalEnv)
      assign("fromM", fromM, .GlobalEnv)
      assign("toM", toM, .GlobalEnv)
    } else {
      toM <- input$toM
      fromM <- input$fromM
      chrM <- input$chrM
      assign("chrM", chrM, .GlobalEnv)
      assign("fromM", fromM, .GlobalEnv)
      assign("toM", toM, .GlobalEnv)
    }
    
    
    RenderDataFrame <- subset(genomicLocationOfMotifs, start >= fromM & end <= toM & seqnames == chrM)%>%
      as.data.frame()
  }, options = list(scrollX = TRUE) )
  
  
  output$LegendsPlot<- renderImage({

    list(
      src = "www/ChromatinStateLegend.png",
      contentType = "image/jpeg",
      alt = "Human/Epigenomics Road Map Legend",
      width= 300,
      height= 400
    )

  }, deleteFile = FALSE)
  
  ###############################################################################################
  ### Preview of Data! 
  ##############################################################################################
  ### TFBS
  output$PreviewOfTheCustomTFBS <- renderDataTable({
    
    # input$file1 will be NULL initially. After the user selects
    # and uploads a file, head of that data file by default,
    # or all rows if selected, will be shown.
    
    req(input$CustomPredictedSitesGenomicSites)
    
    df <- read_delim(input$CustomPredictedSitesGenomicSites$datapath,
                     col_names = input$header,
                     delim  = input$sep,
                     quote = input$quote,
                     skip = input$skip)%>%head(.)
    
    assign( "PreviewOfTheCustomTFBS",  colnames(df), .GlobalEnv)
    
    return( head(df, n=10))
    
  }, options = list(scrollX = TRUE))
  
  output$NamesChromosome <- renderUI({
    req(input$CustomPredictedSitesGenomicSites)
    input$header
    input$sep
    input$quote
    input$skip
    
    ChromColDector <-  grep(
      paste(
        c("seq", "chr"), 
        collapse = "|"), 
      PreviewOfTheCustomTFBS,ignore.case=TRUE,value=TRUE) [1]
    
    
    selectizeInput("customChromosome", label = "Chromosome Information", choices = PreviewOfTheCustomTFBS,
                   selected = ChromColDector)
  })
  
  output$NamesStart <- renderUI({
    req(input$CustomPredictedSitesGenomicSites)
    input$header
    input$sep
    input$quote
    input$skip
    
    StartColDector <-  grep(
      paste(
        c("start", "begin"), 
        collapse = "|"), 
      PreviewOfTheCustomTFBS,ignore.case=TRUE,value=TRUE) [1]
    
    selectizeInput("customStartRange", label = "Start Column Name", choices = PreviewOfTheCustomTFBS,
                   selected = StartColDector)
  })
  
  output$NamesEnd <- renderUI({
    req(input$CustomPredictedSitesGenomicSites)
    input$header
    input$sep
    input$quote
    input$skip
    
    EndColDector <-  grep(
      paste(
        c("end", "stop" , "finish"), 
        collapse = "|"), 
      PreviewOfTheCustomTFBS,ignore.case=TRUE,value=TRUE) [1]
    
    selectizeInput("customEndRange", label = "End Column Name", choices = PreviewOfTheCustomTFBS, 
                   selected = EndColDector)
  })
  
  
  output$NameStrand <- renderUI({
    req(input$CustomPredictedSitesGenomicSites)
    input$header
    input$sep
    input$quote
    input$skip
    
    
    StrandColDector <-  grep(
      paste(
        c("strand"), 
        collapse = "|"), 
      PreviewOfTheCustomTFBS,ignore.case=TRUE,value=TRUE) [1]
    
    if(is.na(StrandColDector)){

    textInput("customStrand",
                label = "No strand column detected, please type it in or leave blank")
    } else {
      selectInput("customStrand", label = "Strand Column Name",
                  selected = StrandColDector, choices = c(PreviewOfTheCustomTFBS, ""), selectize = TRUE)
    }
  })
  
  ###DE gene list
  
  output$PreviewDEGenes <- renderDataTable({
    
    # input$file1 will be NULL initially. After the user selects
    # and uploads a file, head of that data file by default,
    # or all rows if selected, will be shown.
    
    req(input$differenitallyExpressedGenesList)
    
    dfDE <- read_delim(file = input$differenitallyExpressedGenesList$datapath,
                     col_names = input$headerDE,
                     delim  = input$sepDE,
                     quote = input$quoteDE,
                     skip = input$skipDE)%>%head(n=10)
    
    assign("DEdataFrameColnames", colnames(dfDE), .GlobalEnv)
    
    return( head(dfDE, n=10))
    
  }, options = list(scrollX = TRUE))

  output$DEGeneColumn<- renderUI({
    
    
    req(input$differenitallyExpressedGenesList)
    input$headerDE
    input$sepDE
    input$quoteDE
    input$skipDE
    
    geneSymbolColDector <-  grep(
      paste(
        c("symbol", "external", "HGNC"), 
        collapse = "|"), 
      DEdataFrameColnames,ignore.case=TRUE,value=TRUE) [1]
    
    if(is.na(geneSymbolColDector)){
      
      textInput("customDEGeneColumn",
                label = "Failed to Detect DE column, please type it in here!")
    } else {
      selectInput("customDEGeneColumn", label = "Name of the DE gene column",
                  selected = geneSymbolColDector, choices = DEdataFrameColnames, selectize = TRUE)
      
      }
    })




# export Gviz "Screen shot" that can be downloaded
   output$SaveHighResGVIZImage <- downloadHandler('GenomeBrowserImage',
                                                content = function(file) {
                                                  
                                                  tiff(filename = paste0(file, ".tiff"),
                                                       width = 3400, height = 1600,
                                                       units = "px", pointsize = 12,
                                                       compression = c("none"),
                                                       bg = "white", res = 150)
                                
                                                  #Gviz plot code
                                                  plotTracks(trackList =c(humanIdeogramTrack,
                                                                          gHumanTrack,
                                                                          IntearctionTrack,
                                                                          EnhancersHumanChromosomeSpecific,
                                                                          promotertrackChromosomeSpecific,
                                                                          PredictedTFBSTrack,
                                                                          RawMotifInstancesTrack,
                                                                          geneTrackChromosomeSpecific,
                                                                          chromatinStatesTrack),
                                                             sizes= c(1,1.5,3,1,1,1,1,3,3),
                                                             from =input$fromM,
                                                             to= input$toM,
                                                             chromosome= input$chrM,
                                                             cex.title = 1,
                                                             rotation.title = 0,
                                                             showAxis = FALSE,
                                                             background.title = "white",
                                                             lwd.title = 1.1,
                                                             title.width = 3,
                                                             cex.main = 5,
                                                             col = NULL,
                                                             fontcolor.title = "black",
                                                             legend=TRUE)
                                                  
                                                  dev.off()
                                                     })



  
  
  
  
  
  ####################################### 
  ### Gene Ontology Setup the biomart Object & Univserser
  #######################################


  # ## Biomart Object
  #    human_mart <- biomaRt::useMart("ENSEMBL_MART_ENSEMBL",
  #                                    dataset = "hsapiens_gene_ensembl")
  # 
  #     assign("human_mart", human_mart, .GlobalEnv)
  # 
  # 
  #     ## Get the gene Ids for all Symbols
  #           uniservser <- getBM(attributes = c("ucsc", "entrezgene", "external_gene_name"),
  #                         mart = human_mart)
  #           saveRDS(universer, "../DataFiles/Gene Tracks/Human/bioMartObject.rds")
            uniservser <- readRDS("../DataFiles/Gene Tracks/Human/bioMartObject.rds")
      assign("uniservser", uniservser, .GlobalEnv)

      
      ##############################################
      # Data Table setup Click if Tf changes
      ##############################################
      
  GeneOntology <- observe({
    
    
# if the user clicks the button on the TF page
    input$GeneOntology

    # If the user clicks the button on the home page of motifOverlapR
    input$CustomPredictedSitesGenomicSites
    

    isolate({
      if(exists("PredictedTFBS")){
        
      ## Get the gene symbol list from the df
      geneListToConvert <- unique(c(as.character(PredictedTFBS$`Genes Regulated`)))
      

      
      ## Entrez Ids for uniservers  
      ucscToEntrez <- uniservser [uniservser$external_gene_name %in%  geneListToConvert,]

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
      
      # Remove terms without a go id
      GeneOntologyResultsSorted <- filter(GeneOntologyResultsSorted, `Number of DE Genes In Term` > 0)
      
########## Bon Fonerani P Value Adjustment
      
      GeneOntologyResultsSorted$`Adjusted P Value (FDR)` <- p.adjust(GeneOntologyResultsSorted$`P Value`, "fdr")
      
      
      assign("GeneOntologyResultsSorted", GeneOntologyResultsSorted, .GlobalEnv)

      



   } }) # close the isolate brackets
  }) # close the reactive table brackets
      
      


      output$GeneOntologyTableResults <- renderDataTable({
        if(exists("GeneOntologyResultsSorted")){
        input$CustomPredictedSitesGenomicSites
        input$GeneOntology
        
  GeneTable <- filter(GeneOntologyResultsSorted, `Number of DE Genes In Term` >= input$GeneOntologyNumeric &
           `Ontology Type` %in% input$OntologyType & `Adjusted P Value (FDR)` <= input$GeneOntologyPValueCutOff & `P Value` <= input$GeneOntologyRawPvalue
         )
        } else {
          print("Please click Compute TFBS before running this tab!")
        }
  return(GeneTable)
  })


      
      
      
      
      })
