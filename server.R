library(shiny)
library(Seurat)
library(patchwork)
library(ggplot2)
library(glmGamPoi)
library(openxlsx)
library(DT)
library(SeuratWrappers)
library(monocle3)
library(cowplot)

# Load your dataset
NacreMITFIntegrated_allcell_GFP_positive <- readRDS("./NacreMITFIntegrated_allcell_GFP_positive.rds")

# Getting Expressed Genes
gene_expression_data <- GetAssayData(NacreMITFIntegrated_allcell_GFP_positive, layer = "data")
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
      subsetted_data  <- RunUMAP(subsetted_data, reduction = "pca", dims = 1:30)
      subsetted_data <- FindNeighbors(subsetted_data, reduction = "pca", dims = 1:30)
      subsetted_data <- FindClusters(subsetted_data, resolution = 1.2)
    }
  })
  
  output$geneExpressionSelector <- renderUI({
    selectizeInput("expressionGene", "Select Gene of Interest", 
                   choices = c("", expressed_genes), 
                   selected = "")
  })
  
  output$violinPlotGeneSelector <- renderUI({
    selectizeInput("violinPlotGenes", "Select Gene of Interest", 
                   choices = c("", expressed_genes), 
                   selected = "",
                   multiple = TRUE)
  })
  
  output$geneExpressionPlot <- renderPlot({
    data_to_plot <- filtered_data()
    
    if (is.null(input$expressionGene) || input$expressionGene == "") {
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
  
  output$violinPlot <- renderPlot({
    data_to_plot <- filtered_data()
    
    selected_genes <- input$violinPlotGenes
    
    if (is.null(selected_genes) || length(selected_genes) == 0) {
      VlnPlot(data_to_plot, features = ".", group.by = "seurat_clusters", split.by = "Type")
    } else if (length(selected_genes) > 1){    
      VlnPlot(data_to_plot, features = selected_genes, group.by = "seurat_clusters", split.by = "Type", stack = TRUE, sort = FALSE, flip = TRUE)
    } else {
      VlnPlot(data_to_plot, features = selected_genes, group.by = "seurat_clusters", split.by = "Type")
    }  
  })
  
  output$downloadDEGsTable <- downloadHandler(
    filename = function() {
      filename <- ifelse(is.null(input$subsetGene) || input$subsetGene == "", 
                         "NacreMITFIntegrated_allcell_GFP_positive_DEGs.xlsx", 
                         paste("NacreMITFIntegrated_allcell_GFP_positive_", input$subsetGene, "_subset_DEGs.xlsx", sep = ""))
      return(filename)
    },
    content = function(file) {
      data_to_plot <- filtered_data()
      clusters <- levels(data_to_plot$seurat_clusters)
      wb <- createWorkbook()
      for (cluster_id in clusters) {
        tryCatch({
          cell_indices <- WhichCells(data_to_plot, idents = cluster_id)
          if (length(cell_indices) == 0) {
            warning(paste("No cells found for cluster", cluster_id))
            next
          }
          markers <- FindMarkers(data_to_plot, ident.1 = cluster_id,
                                 only.pos = TRUE, min.pct = 0.25, logfc.threshold = 0.25)
          markers_df <- as.data.frame(markers)
          gene_names <- rownames(markers_df)
          
          sheet_name <- paste("Cluster_", cluster_id)
          addWorksheet(wb, sheetName = sheet_name)
          writeData(wb, sheet = sheet_name, cbind("Gene" = gene_names, markers_df))
          
          print(paste("Cluster", cluster_id, "done."))
        }, error = function(e) {
          warning(paste("Error processing cluster", cluster_id, ": ", conditionMessage(e)))
        })
      }
      saveWorkbook(wb, file, overwrite = TRUE)
    }
  )
  
  output$rootNodeSelector <- renderUI({
    cluster_numbers <- levels(filtered_data()$seurat_clusters)
    selectizeInput("selectedRootNodes", "Select Clusters", 
                   choices = c("", cluster_numbers), 
                   selected = "", 
                   multiple = TRUE)
  })
  

  output$rootNodePlots <- renderPlot({
    
    data_to_plot <- filtered_data()
    selected_root_nodes <- input$selectedRootNodes
    print(selected_root_nodes)
    
    getCDSObject <- function(seurat_object) {
      
      seurat_object_name <- deparse(substitute(seurat_object))
      cds <- as.cell_data_set(seurat_object)
      
      cds <- tryCatch({
        cluster_cells(cds)
      }, error = function(e) {
        message("Default clustering method failed: ", e$message)
        message("Attempting Louvain clustering method.")
        cluster_cells(cds, cluster_method = "louvain")
      })
      
      original_umap_clusters <- seurat_object$seurat_clusters
      colData(cds)$orig_umap_cluster <- original_umap_clusters
      
      cds <- learn_graph(cds, use_partition = TRUE, verbose = FALSE)
    }
    
    cds <- getCDSObject(data_to_plot)
    cds_WT <- getCDSObject(subset(data_to_plot, Type == "WT"))
    cds_Nacre <- getCDSObject(subset(data_to_plot, Type == "Nacre"))

    plot_title <- ifelse(is.null(input$subsetGene) || input$subsetGene == "", 
                         "NacreMITFIntegrated_allcell_GFP_positive", 
                         paste("NacreMITFIntegrated_allcell_GFP_positive_", input$subsetGene, "_subset", sep = ""))
    
    if (is.null(selected_root_nodes) || length(selected_root_nodes) == 0) {
      p1 <- plot_cells(cds,
                 color_cells_by = "orig_umap_cluster",
                 label_groups_by_cluster = TRUE,
                 label_leaves = FALSE,
                 label_branch_points = FALSE, 
                 group_label_size = 5)
      p2 <- plot_cells(cds_WT,
                       color_cells_by = "orig_umap_cluster",
                       label_groups_by_cluster = TRUE,
                       label_leaves = FALSE,
                       label_branch_points = FALSE, 
                       group_label_size = 5) +
            ggtitle("WT") + theme(plot.title = element_text(hjust = 0.5))
      p3 <- plot_cells(cds_Nacre,
                       color_cells_by = "orig_umap_cluster",
                       label_groups_by_cluster = TRUE,
                       label_leaves = FALSE,
                       label_branch_points = FALSE, 
                       group_label_size = 5) +
            ggtitle("Nacre") + theme(plot.title = element_text(hjust = 0.5))
      
      p1 + p2 + p3 + plot_annotation(plot_title, theme = theme(plot.title = element_text(face = "bold", hjust = 0.5, size = 15)))
    } else {
      
      getPsuedotime <- function(cds_object, root_node) {
        root_cells <- colnames(cds)[cds$orig_umap_cluster == root_node]
        cds_pseudotime <- order_cells(cds, root_cells = root_cells)
      }
      
      getPsuedotimeGraph <- function(cds_pseudotime, cds_WT_pseudotime, cds_Nacre_pseudotime, root_node) {
        p1 <- plot_cells(cds_pseudotime,
                         color_cells_by = "pseudotime",
                         group_cells_by = "orig_umap_cluster",
                         label_cell_groups = FALSE,
                         label_groups_by_cluster = FALSE,
                         label_leaves = FALSE,
                         label_branch_points = FALSE,
                         label_roots = FALSE,
                         trajectory_graph_color = "grey60") +
          ggtitle(paste("Root Node", root_node)) + theme(plot.title = element_text(hjust = 0.5))
        
        p2 <- plot_cells(cds_WT_pseudotime,
                         color_cells_by = "pseudotime",
                         group_cells_by = "orig_umap_cluster",
                         label_cell_groups = FALSE,
                         label_groups_by_cluster = FALSE,
                         label_leaves = FALSE,
                         label_branch_points = FALSE,
                         label_roots = FALSE,
                         trajectory_graph_color = "grey60") +
          ggtitle("WT") + theme(plot.title = element_text(hjust = 0.5))
        
        p3 <- plot_cells(cds_Nacre_pseudotime,
                         color_cells_by = "pseudotime",
                         group_cells_by = "orig_umap_cluster",
                         label_cell_groups = FALSE,
                         label_groups_by_cluster = FALSE,
                         label_leaves = FALSE,
                         label_branch_points = FALSE,
                         label_roots = FALSE,
                         trajectory_graph_color = "grey60") +
          ggtitle("Nacre") + theme(plot.title = element_text(hjust = 0.5))
        
        p <- p1 + p2 + p3
      }
      
      plot_list <- list()
      
      for (root_node in selected_root_nodes) {
        
        cds <- getPsuedotime(cds, root_node)
        cds_WT <- getPsuedotime(cds_WT, root_node)
        cds_Nacre <- getPsuedotime(cds_Nacre, root_node)
        
        plot_list[[length(plot_list) + 1]] <- getPsuedotimeGraph(cds, cds_WT, cds_Nacre, root_node)
      }
      
      plot_grid(plotlist = plot_list, ncol = 1) + 
        plot_annotation(title = plot_title, 
                        theme = theme(plot.title = 
                                        element_text(hjust = 0.5, face = "bold")))
    }  
  })
}
