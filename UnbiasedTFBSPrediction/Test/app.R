library(shinydashboard)
library(DT)
library(shiny)

ui <- dashboardPage(
  dashboardHeader(title = "Info boxes"),
  dashboardSidebar(),
  dashboardBody(
    fluidRow(DT::dataTableOutput('data')),
    fluidRow(p(class = 'text-center', downloadButton('x3', 'Download Filtered Data')))
  )
)

server <- function(input, output) {
  df <- reactiveValues(data = data.frame(
    Value1 = 1:10,
    Value2 = c("A", "B", "C", "D", "E"),
    stringsAsFactors = FALSE,
    row.names = 1:10
  ))
  
  output$data <- DT::renderDataTable(
    df$data, server =  TRUE, filter = 'top', escape = FALSE, selection = 'none')
  
  # download the filtered data
  output$x3 = downloadHandler('rawMotifGenomicPositions.bed', content = function(file) {
    write.table(as.data.frame(genomicLocationOfMotifs), file  ,
                sep="\t",
                row.names = FALSE,
                col.names = FALSE,
                quote = FALSE,
                append = FALSE)
  })    
}
shinyApp(ui = ui, server = server)
