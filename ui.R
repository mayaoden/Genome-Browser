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
                 div(
                   plotOutput("geneSubsetPlot", height = "740px", width = "910px")
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
                   plotOutput("geneExpressionPlot", height = "800px")  # Fixed dimensions
                 )
               )
             )
    ),
    tabPanel("Violin Plot",
             sidebarLayout(
               sidebarPanel(
                 uiOutput("violinPlotGeneSelector")  # UI output for gene selection
               ),
               mainPanel(
                 div(
                   plotOutput("violinPlot", height = "740px")
                 )
               )
             )
    )
  )
)
