# ============================================================
# 🧠 Multiome PrimingMain Server (Multiome Priming + Visualization Added)
# ============================================================

server <- function(input, output, session) {
  
  dataset <- reactiveVal(NULL)
  
  # ============================================================
  # 📂 Load and Prepare Uploaded Object
  # ============================================================
  observeEvent(input$datafile, {
    req(input$datafile)
    file_path <- input$datafile$datapath
    
    tryCatch({
      obj <- readRDS(file_path)
      
      if (inherits(obj, "Seurat")) {
        obj <- prepare_seurat_object(obj)
        
      } else if (is.list(obj) && all(c("RNA", "ATAC", "genes", "cells") %in% names(obj))) {
        message("✅ Multiome object detected.")
        
      } else {
        showNotification(
          "❌ Uploaded file is neither a Seurat object nor a valid Multiome object.",
          type = "error",
          duration = 6
        )
        return(NULL)
      }
      
      dataset(obj)
      assign("global_dataset", obj, envir = .GlobalEnv)
      
      showNotification("✅ Dataset loaded successfully!", type = "message", duration = 5)
      
    }, error = function(e) {
      showNotification(
        paste("❌ Error loading file:", e$message),
        type = "error",
        duration = 6
      )
    })
  })
  
  # ============================================================
  # 📊 Dataset Summary Panel
  # ============================================================
  output$data_summary <- renderPrint({
    obj <- dataset()
    req(obj)
    
    cat("📁 File Name:", input$datafile$name, "\n")
    cat("🧬 Object Class:", class(obj), "\n")
    
    if (inherits(obj, "Seurat")) {
      cat("🔹 Number of Cells:", ncol(obj), "\n")
      cat("🔹 Number of Genes:", nrow(obj), "\n")
      cat("⭐ Default Assay:", DefaultAssay(obj), "\n")
      
    } else if (is.list(obj) && all(c("RNA", "ATAC", "genes", "cells") %in% names(obj))) {
      cat("🔹 Number of Cells:", length(obj$cells), "\n")
      cat("🔹 Number of Genes:", length(obj$genes), "\n")
      cat("⭐ Contains RNA & ATAC matrices\n")
    }
  })
  
  # ============================================================
  # ⚙️ Load Existing Module Servers
  # ============================================================
  clustering_server("cluster_module")
  entropy_server(input, output, session)
  dispersion_server("dispersion_module")
  pathway_server("pathway_module")
  eps_server("eps_module")
  eps_visualization_server("visual_module")
  benchmark_server("benchmark_module")
  
  # ============================================================
  # 🧬 Multiome Priming Module
  # ============================================================
  multiome_results <- multiome_priming_server(
    id = "multiome_priming_module",
    multiome_obj = dataset
  )
  
  # ============================================================
  # 📊 Multiome Visualization Module
  # ============================================================
  multiome_visualization_server(
    id = "multiome_vis_module",
    priming_results = multiome_results
  )
  
  # ============================================================
  # ✅ Status Messages
  # ============================================================
  cat("✅ Multiome Primingserver initialized.\n")
  cat("📦 Modules loaded: clustering, entropy, dispersion, pathway, EPS, visualization, benchmarking, multiome priming, multiome visualization.\n")
}