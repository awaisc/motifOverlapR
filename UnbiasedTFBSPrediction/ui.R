
options(shiny.maxRequestSize = 200*1024^2, repos = BiocInstaller::biocinstallRepos())
getOption("repos")

## Install all packages needed if packages not installed!!!!
list.of.packages <- c("shiny", 
                      "shinycssloaders",
                      "shinydashboard", 
                      "gridExtra", 
                      "Gviz",  
                      "GenomicInteractions",
                      "GenomicAlignments",
                      "rtracklayer", 
                      "magrittr",
                      "parallel", 
                      "TxDb.Hsapiens.UCSC.hg19.knownGene",
                      "BSgenome.Hsapiens.UCSC.hg19",
                      "org.Hs.eg.db",
                      "TFBSTools",
                      "JASPAR2016", 
                      "seqLogo",
                      "limma",
                      "biomaRt",
                      "AnnotationHub",
                      "readr",
                      "dplyr",
                      "tidyr",
                      "readxl",
                      "tibble",
                      "Homo.sapiens")

new.packages <- list.of.packages[!(list.of.packages %in% installed.packages()[,"Package"])]
## try http:// if https:// URLs are not supported
source("https://bioconductor.org/biocLite.R")
if(length(new.packages)) BiocInstaller::biocLite(new.packages, suppressUpdates = TRUE)
#################################
######Code that needs to be run once
#######################################

library(shiny)
library(shinycssloaders)
library(shinydashboard)
library(gridExtra)

dashboardPage(
  dashboardHeader(title = "motifOverlapR"),
  dashboardSidebar(disable = TRUE),
  
    
  
  dashboardBody(
      div(class = "login",
          uiOutput("uiLogin"),
          textOutput("pass")),
    uiOutput("UIOfMotifOverlapRWithPW")

    
    )
)
