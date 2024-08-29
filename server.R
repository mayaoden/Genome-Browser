library(shiny)
library(Seurat)
library(patchwork)
library(ggplot2)
library(glmGamPoi)

# Load your dataset
NacreMITFIntegrated_allcell_GFP_positive <- readRDS("./NacreMITFIntegrated_allcell_GFP_positive.rds")

# Getting Expressed Genes
gene_expression_data <- GetAssayData(NacreMITFIntegrated_allcell_GFP_positive, slot = "data")
expressed_genes <- rownames(gene_expression_data)[Matrix::rowSums(gene_expression_data > 0) > 0]

server <- function(input, output) {

  output$subsetGeneSelector <- renderUI({
    selectizeInput("subsetGene", "Select Gene to Subset By", 
                   choices = c("", expressed_genes), 
                   selected = "")
  })
  
  filtered_data <- reactive({
    if (is.null(input$subsetGene) || input$subsetGene == "") {
      NacreMITFIntegrated_allcell_GFP_positive
    } else {
      subset_gene <- input$subsetGene
      gene_name <- sym(subset_gene)
      subsetted_data <- NacreMITFIntegrated_allcell_GFP_positive %>% subset(!!gene_name > 0)
      subsetted_data <- SCTransform(subsetted_data, conserve.memory=TRUE)
      subsetted_data <- RunPCA(subsetted_data)
      subsetted_data  <- RunUMAP(subsetted_data,reduction = "pca", dims = 1:30)
      subsetted_data <- FindNeighbors(subsetted_data, reduction = "pca", dims = 1:30)
      subsetted_data <- FindClusters(subsetted_data, resolution = 1.2)
    }
  })
  
  output$geneExpressionSelector <- renderUI({
    selectizeInput("expressionGene", "Select Gene of Interest", 
                   choices = c("", expressed_genes), 
                   selected = "")
  })
  
  output$geneExpressionPlot <- renderPlot({
    
    data_to_plot <- filtered_data()

    if (is.null(input$expressionGene) || input$expressionGene == "") {
      dummy_feature <- rep(0, ncol(data_to_plot))
      names(dummy_feature) <- colnames(data_to_plot)
      data_to_plot[["."]] <- dummy_feature
      
      p1 <- FeaturePlot(data_to_plot, features = ".", label = TRUE, cols = c("gray", "gray")) + NoLegend() + ggtitle(" ")
      p2 <- FeaturePlot(data_to_plot, features = ".", label = TRUE, cols = c("gray", "gray"), split.by = "Type") + NoLegend() +
        patchwork::plot_layout(ncol = 2, nrow = 2)
    } else {
      p1 <- FeaturePlot(data_to_plot, features = input$expressionGene, label = TRUE, cols = c("gray", "red")) + NoLegend()
      p2 <- FeaturePlot(data_to_plot, features = input$expressionGene, label = TRUE, cols = c("gray", "red"), split.by = "Type") + NoLegend() +
        patchwork::plot_layout(ncol = 2, nrow = 2)
    }
    (p1 / p2 + plot_layout(heights = c(3, 6)))
  })
  
  output$geneSubsetPlot <- renderPlot({
    
    data_to_plot <- filtered_data()
    
    p1 <- DimPlot(data_to_plot, label = TRUE)
    p2 <- DimPlot(data_to_plot, group.by = "Type") + 
          (DimPlot(data_to_plot, split.by = "Type", label = TRUE) + NoLegend()) + 
          plot_layout(ncol = 2, widths = c(1, 2))
    
    plot_title <- ifelse(is.null(input$subsetGene) || input$subsetGene == "", 
                         "NacreMITFIntegrated_allcell_GFP_positive", 
                         paste("NacreMITFIntegrated_allcell_GFP_positive_", input$subsetGene, "_subset", sep = ""))
      
    (p1 / p2 ) + plot_annotation(plot_title, theme = theme(plot.title = element_text(face = "bold", hjust = 0.5, size = 15)))
  })
}
