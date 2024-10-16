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
                   plotOutput("geneSubsetPlot", height = "600px", width = "910px")
                 )
               )
             )
    ),
    tabPanel("Cluster Composition",
      fluidRow(
        column(width = 8, offset = 2,  
        plotOutput("clusterCompositionPlot", width = "100%", height = "500px")
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
    ),
    tabPanel("Violin Plot",
             sidebarLayout(
               sidebarPanel(
                 uiOutput("violinPlotGeneSelector")
               ),
               mainPanel(
                 div(
                   plotOutput("violinPlot", height = "600px")
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
                   plotOutput("rootNodePlots") #height = "300px"
                 )
               )
             )
    )
  )
)
