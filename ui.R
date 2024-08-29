library(shiny)
library(shinyjs)
library(plotly)

ui <- fluidPage(
  titlePanel("Genome Browser - NacreMITFIntegrated_allcell_GFP_positive"),
  hr(),
  
  tabsetPanel(
    tabPanel("Gene Expression Subsetting",
             sidebarLayout(
               sidebarPanel(
                 uiOutput("subsetGeneSelector")
               ),
               mainPanel(
                 mainPanel(
                   div(
                     plotOutput("geneSubsetPlot", height = "750px", width = "920px")
                   )
                 )
               )
             )
    ),
    tabPanel("Gene Expression Plot",
             sidebarLayout(
               sidebarPanel(
                 uiOutput("geneExpressionSelector")
               ),
               mainPanel(
                 div(
                   plotOutput("geneExpressionPlot", height = "800px")
                 )
               )
             )
    )
  )
)
