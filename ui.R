library(shiny)
library(shinyjs)
library(plotly)

ui <- fluidPage(
  titlePanel("Genome Browser - NacreMITFIntegrated_allcell_GFP_positive"),
  hr(),
  
  tabsetPanel(
    id = "tabs",
    tabPanel("Gene Expression Subsetting",
             sidebarLayout(
               sidebarPanel(
                 uiOutput("subsetGeneSelector"),
                 downloadButton("downloadDEGsTable", "Download DEGs Table")
               ),
               mainPanel(
                 div(
                   plotOutput("geneSubsetPlot", height = "700px", width = "910px")
                 )
               )
             )
    ),
    tabPanel("Gene Expression Plot",
             sidebarLayout(
               sidebarPanel(
                 uiOutput("geneExpressionSelector")  # UI output for gene selection
               ),
               mainPanel(
                 div(
                   plotOutput("geneExpressionPlot", height = "700px")  # Fixed dimensions
                 )
               )
             )
    ),
    tabPanel("Violin Plot",
             sidebarLayout(
               sidebarPanel(
                 uiOutput("violinPlotGeneSelector")
               ),
               mainPanel(
                 div(
                   plotOutput("violinPlot", height = "700px")
                 )
               )
             )
    ),
    tabPanel("Monocle3 / Psuedotime",
             sidebarLayout(
               sidebarPanel(
                 uiOutput("rootNodeSelector") 
               ),
               mainPanel(
                 div(
                   plotOutput("rootNodePlots", height = "600px")
                 )
               )
             )
    )
  )
)
