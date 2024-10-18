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
                 p("If a cluster is selected, the DE genes between WT and Nacre cells are calculated. 
            Positive log fold-change indicates that the gene is more highly expressed in WT,
            whereas negative values correspond to Nacre. Similarly, pct.1 corresponds with WT 
            and pct.2 with Nacre.")
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
