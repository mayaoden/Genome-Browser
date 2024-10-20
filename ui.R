library(shiny)
library(shinyjs)
library(plotly)

ui <- fluidPage(
  titlePanel("Genome Browser - NacreMITFIntegrated_allcell_GFP_positive"),
  hr(),
  
  tabsetPanel(
    id = "tabs",
    tabPanel("Select Dataset",
             sidebarLayout(
               sidebarPanel(
                 uiOutput("datasetSelector"),
               ),
               mainPanel(
                 div(
                 )
               )
             )
    ),
    tabPanel("Gene Subsetting",
             sidebarLayout(
               sidebarPanel(
                 uiOutput("subsetGeneSelector"),
               ),
               mainPanel(
                 div(
                   plotOutput("geneSubsetPlot", height = "600px", width = "910px")
                 )
               )
             )
    ),
    tabPanel("DEGs",
             sidebarLayout(
               sidebarPanel(
                 uiOutput("clusterSelector"),
                 downloadButton("downloadDEGsTable", "Download DEGs Table"),
                 p(),
                 p("When no cluster is selected, the DEGs between each cluster
                    and all other cells are calculated. If a cluster is 
                    selected, the DEGs between WT and Nacre cells are 
                    calculated. Positive log fold change and pct.1 correspond 
                    with WT, while negative log fold change and pct.2 correspond
                    with Nacre.")
               ),
               mainPanel(
                 div(
                   plotOutput("geneSubsetPlot1", height = "600px", width = "910px")
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
