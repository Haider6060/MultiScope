# ============================================================
# 🖌 EVOscope Multiome Priming Visualization — Server
# ============================================================

multiome_visualization_server <- function(id, priming_results) {
  moduleServer(id, function(input, output, session) {
    
    plot_obj <- reactiveVal(NULL)
    plot_status <- reactiveVal("Run Multiome Priming first, then generate plots.")
    
    clean_theme <- function() {
      theme_classic(base_size = 15) +
        theme(
          panel.background = element_rect(fill = "white", color = NA),
          plot.background = element_rect(fill = "white", color = NA),
          legend.background = element_rect(fill = "white", color = NA),
          axis.line = element_line(color = "black", linewidth = 0.7),
          axis.text = element_text(color = "black"),
          axis.title = element_text(color = "black"),
          plot.title = element_text(face = "bold", hjust = 0.5, size = 18),
          plot.subtitle = element_text(hjust = 0.5, size = 12),
          legend.position = "right"
        )
    }
    
    output$plot_status <- renderText({
      plot_status()
    })
    
    observeEvent(input$generate_plot, {
      
      req(priming_results$cell_scores)
      req(priming_results$gene_scores)
      
      cell_scores <- priming_results$cell_scores
      gene_scores <- priming_results$gene_scores
      
      # ======================================================
      # 0. Cell Priming UMAP
      # ======================================================
      if (input$plot_type == "priming_umap") {
        
        if (!all(c("UMAP1", "UMAP2") %in% colnames(cell_scores))) {
          
          set.seed(123)
          
          umap_input <- cell_scores[, c("RNA_ATAC_Correlation", "Regulatory_Difference_Score")]
          umap_input <- as.data.frame(lapply(umap_input, as.numeric))
          umap_input <- umap_input[complete.cases(umap_input), , drop = FALSE]
          
          if (nrow(umap_input) < 10) {
            showNotification("❌ Not enough cells to generate UMAP.", type = "error", duration = 6)
            plot_status("❌ Not enough cells to generate UMAP.")
            return(NULL)
          }
          
          um <- uwot::umap(
            umap_input,
            n_neighbors = 15,
            min_dist = 0.30,
            metric = "euclidean"
          )
          
          cell_scores <- cell_scores[complete.cases(
            cell_scores[, c("RNA_ATAC_Correlation", "Regulatory_Difference_Score")]
          ), , drop = FALSE]
          
          cell_scores$UMAP1 <- um[, 1]
          cell_scores$UMAP2 <- um[, 2]
        }
        
        p <- ggplot(
          cell_scores,
          aes(
            x = UMAP1,
            y = UMAP2,
            color = Regulatory_Difference_Score
          )
        ) +
          geom_point(size = 1.6, alpha = 0.85) +
          viridis::scale_color_viridis(
            option = "plasma",
            name = "Priming Score"
          ) +
          clean_theme() +
          theme(
            axis.text = element_blank(),
            axis.ticks = element_blank(),
            axis.line = element_blank()
          ) +
          labs(
            title = "Cell Priming UMAP",
            subtitle = "Cells colored by RNA–ATAC regulatory difference score",
            x = "UMAP 1",
            y = "UMAP 2"
          )
        
        plot_obj(p)
        plot_status("✅ Cell priming UMAP generated.")
      }
      
      # ======================================================
      # 1. Cell Priming Score Distribution
      # ======================================================
      if (input$plot_type == "cell_distribution") {
        
        p <- ggplot(cell_scores, aes(x = Regulatory_Difference_Score)) +
          geom_histogram(
            aes(y = after_stat(density)),
            bins = 40,
            fill = "#2E86DE",
            color = "white",
            alpha = 0.85
          ) +
          geom_density(
            color = "#E74C3C",
            linewidth = 1.3
          ) +
          clean_theme() +
          labs(
            title = "Cell Priming Score Distribution",
            subtitle = "Distribution of RNA–ATAC regulatory difference scores across cells",
            x = "Regulatory Difference Score",
            y = "Density"
          )
        
        plot_obj(p)
        plot_status("✅ Cell priming distribution plot generated.")
      }
      
      # ======================================================
      # 2. Top 20 High-Priming Cells — Lollipop
      # ======================================================
      if (input$plot_type == "top_cells") {
        
        df <- head(
          cell_scores[order(cell_scores$Regulatory_Difference_Score, decreasing = TRUE), ],
          20
        )
        
        df$Cell <- factor(df$Cell, levels = rev(df$Cell))
        
        p <- ggplot(df, aes(x = Cell, y = Regulatory_Difference_Score)) +
          geom_segment(
            aes(
              x = Cell,
              xend = Cell,
              y = min(Regulatory_Difference_Score) * 0.95,
              yend = Regulatory_Difference_Score
            ),
            linewidth = 0.8,
            color = "#E74C3C"
          ) +
          geom_point(
            size = 3.5,
            color = "#C0392B"
          ) +
          coord_flip() +
          clean_theme() +
          labs(
            title = "Top 20 High-Priming Cells",
            subtitle = "Cells ranked by highest regulatory difference score",
            x = "Cell",
            y = "Regulatory Difference Score"
          )
        
        plot_obj(p)
        plot_status("✅ Top 20 high-priming cells plot generated.")
      }
      
      # ======================================================
      # 3. Top 20 Primed Genes — Line-Dot Plot
      # ======================================================
      if (input$plot_type == "top_genes") {
        
        df <- head(
          gene_scores[order(gene_scores$Regulatory_Difference_Score, decreasing = TRUE), ],
          20
        )
        
        df$Gene <- factor(df$Gene, levels = df$Gene)
        
        p <- ggplot(df, aes(x = Gene, y = Regulatory_Difference_Score, group = 1)) +
          geom_line(
            color = "#27AE60",
            linewidth = 1.2
          ) +
          geom_point(
            aes(size = Regulatory_Difference_Score),
            color = "#E74C3C",
            alpha = 0.95
          ) +
          scale_size_continuous(range = c(3, 7)) +
          clean_theme() +
          theme(
            axis.text.x = element_text(angle = 45, hjust = 1),
            legend.position = "none"
          ) +
          labs(
            title = "Top 20 Primed Genes",
            subtitle = "Line-dot profile of top RNA–ATAC regulatory difference genes",
            x = "Gene",
            y = "Regulatory Difference Score"
          )
        
        plot_obj(p)
        plot_status("✅ Top 20 primed genes line-dot plot generated.")
      }
      
      # ======================================================
      # 4. High vs Low Priming Boxplot
      # ======================================================
      if (input$plot_type == "high_low_boxplot") {
        
        df <- cell_scores[cell_scores$Priming_Group %in% c("High_Priming", "Low_Priming"), ]
        df$Priming_Group <- factor(df$Priming_Group, levels = c("Low_Priming", "High_Priming"))
        
        p <- ggplot(df, aes(x = Priming_Group, y = Regulatory_Difference_Score, fill = Priming_Group)) +
          geom_boxplot(
            color = "black",
            linewidth = 0.7,
            outlier.size = 1.7,
            alpha = 0.9
          ) +
          scale_fill_manual(
            values = c(
              "Low_Priming" = "#3498DB",
              "High_Priming" = "#E74C3C"
            )
          ) +
          clean_theme() +
          labs(
            title = "High vs Low Priming Cells",
            subtitle = "Comparison of regulatory difference scores between priming states",
            x = "Priming Group",
            y = "Regulatory Difference Score",
            fill = "Priming Group"
          )
        
        plot_obj(p)
        plot_status("✅ High vs Low priming boxplot generated.")
      }
      
      # ======================================================
      # 5. Priming Group Counts
      # ======================================================
      if (input$plot_type == "group_counts") {
        
        df <- as.data.frame(table(cell_scores$Priming_Group))
        colnames(df) <- c("Priming_Group", "Count")
        df$Priming_Group <- factor(df$Priming_Group, levels = c("Low_Priming", "Middle", "High_Priming"))
        
        p <- ggplot(df, aes(x = Priming_Group, y = Count, fill = Priming_Group)) +
          geom_col(
            color = "black",
            width = 0.65,
            alpha = 0.9
          ) +
          scale_fill_manual(
            values = c(
              "Low_Priming" = "#3498DB",
              "Middle" = "#95A5A6",
              "High_Priming" = "#E74C3C"
            )
          ) +
          clean_theme() +
          labs(
            title = "Priming Group Counts",
            subtitle = "Number of cells assigned to each priming state",
            x = "Priming Group",
            y = "Number of Cells",
            fill = "Priming Group"
          )
        
        plot_obj(p)
        plot_status("✅ Priming group count plot generated.")
      }
      
      # ======================================================
      # 6. Top Pathway Enrichment
      # ======================================================
      if (input$plot_type == "primed_pathways") {
        
        if (!exists("gene_sets") || is.null(gene_sets)) {
          showNotification("❌ Hallmark gene sets not loaded.", type = "error", duration = 6)
          plot_status("❌ Hallmark gene sets not loaded.")
          return(NULL)
        }
        
        pathway_list <- gene_sets
        
        if (inherits(pathway_list, "GeneSetCollection")) {
          pathway_list <- lapply(pathway_list, function(x) GSEABase::geneIds(x))
          names(pathway_list) <- sapply(gene_sets, function(x) GSEABase::setName(x))
        }
        
        gene_scores_df <- as.data.frame(gene_scores)
        gene_scores_df <- gene_scores_df[!is.na(gene_scores_df[["Regulatory_Difference_Score"]]), ]
        gene_scores_df <- gene_scores_df[!is.na(gene_scores_df[["Gene"]]), ]
        
        stats_vec <- setNames(
          as.numeric(gene_scores_df[["Regulatory_Difference_Score"]]),
          as.character(gene_scores_df[["Gene"]])
        )
        
        stats_vec <- stats_vec[!is.na(stats_vec)]
        stats_vec <- sort(stats_vec, decreasing = TRUE)
        
        if (length(stats_vec) < 20) {
          showNotification("❌ Not enough ranked genes for pathway enrichment.", type = "error", duration = 6)
          plot_status("❌ Not enough ranked genes for pathway enrichment.")
          return(NULL)
        }
        
        fgsea_res <- fgsea::fgseaMultilevel(
          pathways = pathway_list,
          stats = stats_vec,
          minSize = 5,
          maxSize = 500,
          scoreType = "pos"
        )
        
        fgsea_res <- as.data.frame(fgsea_res)
        fgsea_res <- fgsea_res[!is.na(fgsea_res$NES), ]
        fgsea_res <- fgsea_res[order(fgsea_res$NES, decreasing = TRUE), ]
        top_pathways <- head(fgsea_res, 10)
        
        if (nrow(top_pathways) == 0) {
          showNotification("❌ No enriched pathways found.", type = "warning", duration = 6)
          plot_status("❌ No enriched pathways found.")
          return(NULL)
        }
        
        top_pathways$pathway <- gsub("HALLMARK_", "", top_pathways$pathway)
        top_pathways$pathway <- gsub("_", " ", top_pathways$pathway)
        top_pathways$pathway <- stringr::str_to_title(top_pathways$pathway)
        top_pathways$pathway <- factor(top_pathways$pathway, levels = rev(top_pathways$pathway))
        
        p <- ggplot(top_pathways, aes(x = pathway, y = NES, fill = NES)) +
          geom_col(color = "black", alpha = 0.9, width = 0.72) +
          geom_point(aes(y = NES), size = 3.2, color = "black") +
          coord_flip() +
          viridis::scale_fill_viridis(option = "plasma", name = "NES") +
          clean_theme() +
          labs(
            title = "Top Pathway Enrichment",
            subtitle = "Hallmark pathways enriched among high regulatory-priming genes",
            x = "Pathway",
            y = "Normalized Enrichment Score"
          )
        
        plot_obj(p)
        plot_status("✅ Top pathway enrichment plot generated.")
      }
    })
    
    output$priming_plot <- renderPlot({
      req(plot_obj())
      plot_obj()
    })
    
    output$download_plot <- downloadHandler(
      filename = function() {
        paste0(input$plot_type, "_1000dpi_", Sys.Date(), ".png")
      },
      content = function(file) {
        req(plot_obj())
        ggsave(
          filename = file,
          plot = plot_obj(),
          width = 10,
          height = 7,
          dpi = 1000,
          bg = "white"
        )
      }
    )
  })
}