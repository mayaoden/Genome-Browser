library(shiny)
library(shinyjs)
library(plotly)

ui <- fluidPage(
  titlePanel("Genome Browser - NacreMITFIntegrated_allcell_GFP_positive"),
  hr(), 
  
  sidebarLayout(
    sidebarPanel(
      uiOutput("subsetGeneSelector"),
      uiOutput("geneExpressionSelector")
    ),
    
    mainPanel(
      tabsetPanel(
        tabPanel("Gene Expression Plot",
                 plotOutput("geneExpressionPlot", width = "800px", height = "900px"),
                 style = "display: flex; justify-content: center; align-items: center;"
        ),
        tabPanel("Another Tab",
                 # Add content for the second tab here
                 p("Content for the second tab goes here.")
        )
        # You can add more tabs as needed
      )
    )
  )
)