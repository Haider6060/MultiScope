# ============================================================
# 🧬 EVOscope Multiome Priming Module — Server
# ============================================================

multiome_priming_server <- function(id, multiome_obj) {
  moduleServer(id, function(input, output, session) {
    
    priming_results <- reactiveValues(
      cell_scores = NULL,
      gene_scores = NULL,
      status_msg = "Module ready. Upload a valid multiome object and click Run."
    )
    
    output$priming_status <- renderText({
      priming_results$status_msg
    })
    
    output$progress_text <- renderUI({
      HTML(gsub("\n", "<br>", priming_results$status_msg))
    })
    
    observeEvent(input$run_priming, {
      
      obj <- multiome_obj()
      req(obj)
      
      if (!is.list(obj) || !all(c("RNA", "ATAC") %in% names(obj))) {
        showNotification("❌ Please upload a valid multiome .rds file containing RNA and ATAC matrices.",
                         type = "error", duration = 6)
        priming_results$status_msg <- "❌ Invalid multiome object. Required fields: RNA and ATAC."
        return(NULL)
      }
      
      withProgress(message = "Running Multiome Priming Analysis", value = 0, {
        
        incProgress(0.05, detail = "Loading object...")
        priming_results$status_msg <- "🔄 Step 1/6: Loading multiome object..."
        
        rna_mat <- obj$RNA
        atac_mat <- obj$ATAC
        
        if (is.null(rna_mat) || is.null(atac_mat)) {
          priming_results$status_msg <- "❌ RNA or ATAC matrix is missing."
          showNotification("❌ RNA or ATAC matrix is missing.", type = "error", duration = 6)
          return(NULL)
        }
        
        incProgress(0.10, detail = "Preparing matrix names...")
        priming_results$status_msg <- "🔄 Step 2/6: Preparing matrix names..."
        
        if (!is.null(obj$genes)) {
          rownames(rna_mat) <- obj$genes
          rownames(atac_mat) <- obj$genes
        }
        
        if (!is.null(obj$cells)) {
          colnames(rna_mat) <- obj$cells
          colnames(atac_mat) <- obj$cells
        }
        
        if (is.null(rownames(rna_mat)) || is.null(rownames(atac_mat)) ||
            is.null(colnames(rna_mat)) || is.null(colnames(atac_mat))) {
          priming_results$status_msg <- "❌ RNA/ATAC matrices must have rownames and colnames."
          showNotification("❌ Missing gene or cell names.", type = "error", duration = 6)
          return(NULL)
        }
        
        incProgress(0.20, detail = "Matching genes and cells...")
        priming_results$status_msg <- "🔄 Step 3/6: Matching genes and cells..."
        
        common_cells <- intersect(colnames(rna_mat), colnames(atac_mat))
        common_genes <- intersect(rownames(rna_mat), rownames(atac_mat))
        
        if (length(common_cells) == 0 || length(common_genes) == 0) {
          priming_results$status_msg <- "❌ No matched cells/genes found."
          showNotification("❌ No matched cells or genes found.", type = "error", duration = 6)
          return(NULL)
        }
        
        rna_mat <- rna_mat[common_genes, common_cells, drop = FALSE]
        atac_mat <- atac_mat[common_genes, common_cells, drop = FALSE]
        
        incProgress(0.35, detail = "Computing cell-level RNA–ATAC correlations...")
        priming_results$status_msg <- paste0(
          "✅ Matched ", length(common_genes), " genes and ", length(common_cells), " cells.\n",
          "🔄 Step 4/6: Computing cell-level RNA–ATAC correlations..."
        )
        
        cell_cor <- sapply(seq_len(ncol(rna_mat)), function(i) {
          suppressWarnings(
            cor(
              as.numeric(rna_mat[, i]),
              as.numeric(atac_mat[, i]),
              method = input$cor_method,
              use = "complete.obs"
            )
          )
        })
        
        names(cell_cor) <- colnames(rna_mat)
        
        decoupling_score <- 1 - cell_cor
        names(decoupling_score) <- colnames(rna_mat)
        
        incProgress(0.60, detail = "Classifying priming groups...")
        priming_results$status_msg <- "🔄 Step 5/6: Classifying High/Low Priming cells..."
        
        high_thresh <- quantile(decoupling_score, probs = input$high_cutoff / 100, na.rm = TRUE)
        low_thresh  <- quantile(decoupling_score, probs = input$low_cutoff / 100, na.rm = TRUE)
        
        priming_group <- rep("Middle", length(decoupling_score))
        priming_group[decoupling_score >= high_thresh] <- "High_Priming"
        priming_group[decoupling_score <= low_thresh] <- "Low_Priming"
        
        cell_scores <- data.frame(
          Cell = names(decoupling_score),
          RNA_ATAC_Correlation = as.numeric(cell_cor),
          Regulatory_Difference_Score = as.numeric(decoupling_score),
          Priming_Group = priming_group,
          stringsAsFactors = FALSE
        )
        
        incProgress(0.80, detail = "Computing gene-level RNA–ATAC correlations...")
        priming_results$status_msg <- "🔄 Step 6/6: Computing gene-level RNA–ATAC correlations..."
        
        gene_cor <- sapply(seq_len(nrow(rna_mat)), function(i) {
          suppressWarnings(
            cor(
              as.numeric(rna_mat[i, ]),
              as.numeric(atac_mat[i, ]),
              method = "pearson",
              use = "complete.obs"
            )
          )
        })
        
        names(gene_cor) <- rownames(rna_mat)
        
        gene_decoupling <- 1 - gene_cor
        names(gene_decoupling) <- rownames(rna_mat)
        
        gene_scores <- data.frame(
          Gene = names(gene_decoupling),
          RNA_ATAC_Correlation = as.numeric(gene_cor),
          Regulatory_Difference_Score = as.numeric(gene_decoupling),
          stringsAsFactors = FALSE
        )
        
        gene_scores <- gene_scores[!is.na(gene_scores$RNA_ATAC_Correlation), ]
        
        incProgress(1.0, detail = "Finalizing results...")
        
        priming_results$cell_scores <- cell_scores
        priming_results$gene_scores <- gene_scores
        
        priming_results$status_msg <- paste0(
          "✅ Multiome Priming Analysis completed.\n\n",
          "Cells analyzed: ", nrow(cell_scores), "\n",
          "Genes analyzed: ", nrow(gene_scores), "\n",
          "High-Priming cells: ", sum(cell_scores$Priming_Group == "High_Priming"), "\n",
          "Low-Priming cells: ", sum(cell_scores$Priming_Group == "Low_Priming"), "\n\n",
          "You can now download tables or generate plots."
        )
        
        showNotification("✅ Multiome Priming Analysis completed.", type = "message", duration = 5)
      })
    })
    
    output$download_cell_scores <- downloadHandler(
      filename = function() {
        paste0("Multiome_Cell_Priming_Scores_", Sys.Date(), ".csv")
      },
      content = function(file) {
        req(priming_results$cell_scores)
        write.csv(priming_results$cell_scores, file, row.names = FALSE)
      }
    )
    
    output$download_gene_scores <- downloadHandler(
      filename = function() {
        paste0("Multiome_Gene_Priming_Scores_", Sys.Date(), ".csv")
      },
      content = function(file) {
        req(priming_results$gene_scores)
        write.csv(priming_results$gene_scores, file, row.names = FALSE)
      }
    )
    
    return(priming_results)
  })
}