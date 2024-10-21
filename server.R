library(shiny)
library(Seurat)
library(patchwork)
library(ggplot2)
library(glmGamPoi)
library(openxlsx)
library(DT)
library(R.utils)
library(SeuratWrappers)
library(monocle3)
library(cowplot)

server <- function(input, output) {
  
  selected_dataset <- reactive({
    req(input$selectedDataset)
    readRDS(paste0("./", input$selectedDataset, ".rds"))
  })
  
  expressed_genes <- reactive({
    req(selected_dataset())
    gene_expression_data <- GetAssayData(selected_dataset(), layer = "data")
    rownames(gene_expression_data)[Matrix::rowSums(gene_expression_data > 0) > 0]
  })
  
  output$subsetGeneSelector <- renderUI({
    selectizeInput("subsetGene", "Select Gene to Subset By", 
                   choices = c("", expressed_genes()), 
                   selected = "")
  })
  
  filtered_data <- reactive({
    if (is.null(input$subsetGene) || input$subsetGene == "") {
      selected_dataset()
    } else {
      subset_gene <- input$subsetGene
      gene_name <- sym(subset_gene)
      subsetted_data <- selected_dataset() %>% subset(!!gene_name > 0)
      subsetted_data <- SCTransform(subsetted_data, conserve.memory=TRUE)
      subsetted_data <- RunPCA(subsetted_data)
      subsetted_data  <- RunUMAP(subsetted_data, reduction = "pca", dims = 1:30)
      subsetted_data <- FindNeighbors(subsetted_data, reduction = "pca", dims = 1:30)
      subsetted_data <- FindClusters(subsetted_data, resolution = 1.2)
    }
  })
  
  plot_title <- reactive({
    if (is.null(input$subsetGene) || input$subsetGene == "") {
      input$selectedDataset
    } else {
      paste(input$selectedDataset, "_", input$subsetGene, "_subset", sep = "")
    }
  })
  
  output$geneExpressionSelector <- renderUI({
    selectizeInput("expressionGene", "Select Gene of Interest", 
                   choices = c("", expressed_genes()), 
                   selected = "")
  })
  
  output$violinPlotGeneSelector <- renderUI({
    selectizeInput("violinPlotGenes", "Select Gene of Interest", 
                   choices = c("", expressed_genes()), 
                   selected = "",
                   multiple = TRUE)
  })
  
  output$geneExpressionPlot <- renderPlot({

    data_to_plot <- filtered_data();
    
    if (is.null(input$expressionGene) || input$expressionGene == "") {
      p1 <- FeaturePlot(data_to_plot, features = ".", label = TRUE, label.size = 5, cols = c("gray", "gray")) + NoLegend() + ggtitle(" ")
      p2 <- FeaturePlot(data_to_plot, features = ".", label = TRUE, label.size = 5, cols = c("gray", "gray"), split.by = "Type") + NoLegend() +
        patchwork::plot_layout(ncol = 2, nrow = 2)
    } else {
      p1 <- FeaturePlot(data_to_plot, features = input$expressionGene, label = TRUE, label.size = 5, cols = c("gray", "red")) + NoLegend()  + 
        theme(plot.title = element_text(face = "plain"))
      p2 <- FeaturePlot(data_to_plot, features = input$expressionGene, label = TRUE, label.size = 5, cols = c("gray", "red"), split.by = "Type") + NoLegend() +
        patchwork::plot_layout(ncol = 2, nrow = 2)
    }
    (p1 / p2 + plot_layout(heights = c(3, 6))) + 
      plot_annotation(plot_title(), theme = theme(plot.title = element_text(face = "bold", hjust = 0.5, size = 15)))
  })
  
  output$rootNodeSelector <- renderUI({
    cluster_numbers <- levels(filtered_data()$seurat_clusters)
    selectizeInput("selectedRootNode", "Select Root Node Clusters", 
                   choices = c("", cluster_numbers), 
                   selected = "")
  })
  
  output$geneSubsetPlot <- renderPlot({
    data_to_plot <- filtered_data()
    
    p1 <- DimPlot(data_to_plot, label = TRUE, label.size = 5)
    p2 <- DimPlot(data_to_plot, group.by = "Type") + 
      (DimPlot(data_to_plot, split.by = "Type", label = TRUE, label.size = 5) + NoLegend()) + 
      plot_layout(ncol = 2, widths = c(1, 2))
    
    (p1 / p2 ) + plot_annotation(plot_title(), theme = theme(plot.title = element_text(face = "bold", hjust = 0.5, size = 15)))
  })
  
  output$geneSubsetPlot1 <- renderPlot({
    data_to_plot <- filtered_data()
    
    p1 <- DimPlot(data_to_plot, label = TRUE, label.size = 5)
    p2 <- DimPlot(data_to_plot, group.by = "Type") + 
      (DimPlot(data_to_plot, split.by = "Type", label = TRUE, label.size = 5) + NoLegend()) + 
      plot_layout(ncol = 2, widths = c(1, 2))
    
    (p1 / p2 ) + plot_annotation(plot_title(), theme = theme(plot.title = element_text(face = "bold", hjust = 0.5, size = 15)))
  })
  
  output$violinPlot <- renderPlot({
    data_to_plot <- filtered_data()
    
    selected_genes <- input$violinPlotGenes
    
    if (is.null(selected_genes) || length(selected_genes) == 0) {
      VlnPlot(data_to_plot, features = ".", group.by = "seurat_clusters", split.by = "Type") + 
        theme(plot.title = element_text(face = "plain")) + 
        plot_annotation(plot_title(), theme = theme(plot.title = element_text(face = "bold", hjust = 0.5, size = 15)))
    } else if (length(selected_genes) > 1){    
      VlnPlot(data_to_plot, features = selected_genes, group.by = "seurat_clusters", split.by = "Type", stack = TRUE, sort = FALSE, flip = TRUE) + 
        theme(plot.title = element_text(face = "plain")) + 
        plot_annotation(plot_title(), theme = theme(plot.title = element_text(face = "bold", hjust = 0.5, size = 15)))
    } else {
      VlnPlot(data_to_plot, features = selected_genes, group.by = "seurat_clusters", split.by = "Type") + 
        theme(plot.title = element_text(face = "plain")) + 
        plot_annotation(plot_title(), theme = theme(plot.title = element_text(face = "bold", hjust = 0.5, size = 15)))
    }  
  })
  
  output$downloadDEGsTable <- downloadHandler(
    filename = function() {
      plot_title <- plot_title()
      if (input$cluster == "All") {
        plot_title <- paste0(plot_title, "_DEGS.xlsx") 
      } else {
        plot_title <- paste0(plot_title, "_cluster_", input$cluster, "_DEGs.xlsx")
      }
      return(plot_title)
    },
    content = function(file) {
      data_to_plot <- filtered_data()
      wb <- createWorkbook()
      
      if (input$cluster == "All") {
        clusters <- levels(data_to_plot$seurat_clusters)
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
      } else {
        cluster <- input$cluster
        data_to_plot <- subset(data_to_plot, seurat_clusters == cluster)
        Idents(data_to_plot) <- data_to_plot$Type
        markers <- FindMarkers(data_to_plot, ident.1 = "WT", ident.2 = "Nacre")
        markers_df <- as.data.frame(markers)
        gene_names <- rownames(markers_df)
        addWorksheet(wb, "Sheet1")
        writeData(wb, "Sheet1", cbind("Gene" = gene_names, markers_df))
      }
      saveWorkbook(wb, file, overwrite = TRUE)
    }
  )
  
  cds_objects <- reactive({
    data_to_plot <- filtered_data()
    
    getCDSObject <- function(seurat_object) {
      cds <- as.cell_data_set(seurat_object)
      cds <- tryCatch({
        cluster_cells(cds)
      }, error = function(e) {
        message("Default clustering method failed: ", e$message)
        message("Attempting Louvain clustering method.")
        cluster_cells(cds, cluster_method = "louvain")
      })
      colData(cds)$orig_umap_cluster <- seurat_object$seurat_clusters
      cds <- learn_graph(cds, use_partition = TRUE, verbose = FALSE)
      return(cds)
    }
    
    list(
      cds_WT = getCDSObject(subset(data_to_plot, Type == "WT")),
      cds_Nacre = getCDSObject(subset(data_to_plot, Type == "Nacre"))
    )
  })
  
  output$rootNodePlots <- renderPlot({
    selected_root_node <- input$selectedRootNode
    cds_list <- cds_objects()
    
    p1 <- plot_cells(cds_list$cds_WT, 
                     color_cells_by = "orig_umap_cluster", 
                     group_cells_by = "orig_umap_cluster", 
                     label_groups_by_cluster = TRUE, 
                     label_leaves = FALSE, 
                     label_branch_points = FALSE, 
                     group_label_size = 5) +
      ggtitle("WT") + theme(plot.title = element_text(hjust = 0.5))
    
    p2 <- plot_cells(cds_list$cds_Nacre,
                     color_cells_by = "orig_umap_cluster",
                     group_cells_by = "orig_umap_cluster", 
                     label_groups_by_cluster = TRUE,
                     label_leaves = FALSE,
                     label_branch_points = FALSE,
                     group_label_size = 5) +
      ggtitle("Nacre") + theme(plot.title = element_text(hjust = 0.5))
    
    p_monocle3 <- p1 + p2
    
    if (is.null(selected_root_node) || selected_root_node == "") {
      return(p_monocle3 + 
              plot_annotation(
              plot_title(),
              theme = theme(plot.title = element_text(face = "bold", hjust = 0.5, size = 15))
      ))
    }
    
    getPsuedotime <- function(cds_object, root_node) {
      root_group <- colnames(cds_object)[cds_object$orig_umap_cluster == root_node]
      order_cells(cds_object, root_cells = root_group)
    }
    
    getPsuedotimeGraph <- function(cds_WT_pseudotime, cds_Nacre_pseudotime, root_node) {
      p1 <- plot_cells(cds_WT_pseudotime,
                       color_cells_by = "pseudotime",
                       group_cells_by = "orig_umap_cluster",
                       label_cell_groups = FALSE,
                       label_groups_by_cluster = TRUE,
                       label_leaves = FALSE,
                       label_branch_points = FALSE,
                       label_roots = FALSE,
                       trajectory_graph_color = "grey60") +
        ggtitle("WT") + theme(plot.title = element_text(hjust = 0.5))
      
      p2 <- plot_cells(cds_Nacre_pseudotime,
                       color_cells_by = "pseudotime",
                       group_cells_by = "orig_umap_cluster",
                       label_cell_groups = FALSE,
                       label_groups_by_cluster = TRUE,
                       label_leaves = FALSE,
                       label_branch_points = FALSE,
                       label_roots = FALSE,
                       trajectory_graph_color = "grey60") +
        ggtitle("Nacre") + theme(plot.title = element_text(hjust = 0.5))
      
      p1 + p2
    }
    
    cds_WT <- getPsuedotime(cds_list$cds_WT, selected_root_node)
    cds_Nacre <- getPsuedotime(cds_list$cds_Nacre, selected_root_node)
    
    p_pseudotime <- getPsuedotimeGraph(cds_WT, cds_Nacre, selected_root_node) +
      plot_annotation(
        subtitle = paste0("Root Node: ", selected_root_node),
        theme = theme(plot.subtitle = element_text(hjust = 0.5, size = 12))
      )
    p_combinded <- p_monocle3 / p_pseudotime + plot_annotation(
      plot_title(),
      subtitle = paste0("Root Node: ", selected_root_node),
      theme = theme(plot.title = element_text(face = "bold", hjust = 0.5, size = 15),
                    plot.subtitle = element_text(hjust = 0.5, size = 12))
    )
    p_combinded
  }, height = function() {
    if (is.null(input$selectedRootNode) || input$selectedRootNode == "") {
      return(300) 
    } else {
      return(600)
    }
  })
  
  output$clusterCompositionPlot <- renderPlot({ 
    data_to_plot <- filtered_data()
    total_cells_per_type <- table(data_to_plot$Type)
    total_cells_per_type_df <- as.data.frame(total_cells_per_type)
    colnames(total_cells_per_type_df) <- c("Type", "TotalCells")
    total_cells_per_type_df$Percentage <- total_cells_per_type_df$TotalCells / 
      sum(total_cells_per_type_df$TotalCells) * 100
    pie_chart <- ggplot(total_cells_per_type_df, aes(x = "", y = TotalCells, 
                                                     fill = Type)) + geom_bar(stat = "identity") + geom_text(aes(
                                                       label = paste0(round(Percentage), "%")), position = position_stack(
                                                         vjust = 0.5)) + coord_polar("y", start = 0) + scale_fill_discrete(name = "Type") + theme_void()
    cell_counts <- table(data_to_plot$seurat_clusters, 
                         data_to_plot$Type)
    cell_counts_df <- as.data.frame(cell_counts)
    colnames(cell_counts_df) <- c("Cluster", "Type", "CellCount")
    total_cells_per_type <- tapply(data_to_plot$Type, 
                                   data_to_plot$Type, length)
    cell_counts_df$Proportion <- cell_counts_df$CellCount / 
      total_cells_per_type[cell_counts_df$Type]
    
    bar_graph <- ggplot(cell_counts_df, aes(x = Cluster, y = Proportion, 
                                            fill = Type)) +
      geom_bar(stat = "identity", position = "dodge") +
      labs(x = "Cluster", y = "Proportion of Cells") + 
      theme_minimal() +
      theme(axis.text.x = element_text(angle = 45, hjust = 1))
    
    (bar_graph + NoLegend()) + pie_chart + plot_layout(ncol = 2, widths = c(4, 1)) + 
      plot_annotation(plot_title(), theme = theme(plot.title = element_text(face = "bold", hjust = 0.5, size = 15)))
  })
  
  output$clusterSelector <- renderUI({
    cluster_numbers <- levels(filtered_data()$seurat_clusters)
    selectizeInput("cluster", "Select Clusters", 
                   choices = c("All", cluster_numbers), 
                   selected = "All")
  })
  
  output$datasetSelector <- renderUI({
    selectizeInput("selectedDataset", "Select Dataset", 
                   choices = c("", 
                               "NacreMITFIntegrated_allcell_GFP_positive",
                               "NacreMITFIntegrated_allcell_GFP_negative"), 
                   selected = "")
  })
}