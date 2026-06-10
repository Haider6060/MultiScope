# ============================================================
# 🧬  EVOscope Explorer
# ============================================================

library(shiny)
library(shinydashboard)

# ============================================================
# 📦 Source UI Module Files
# ============================================================
source("multiome_priming_ui.R")
source("multiome_visualization_ui.R")  # New visualization module

ui <- dashboardPage(
  
  # ---------- HEADER ----------
  dashboardHeader(
    title = span("🧬 EVOscope", style = "font-weight:600; color:white;"),
    titleWidth = 280
  ),
  
  # ---------- SIDEBAR ----------
  dashboardSidebar(
    width = 280,
    tags$head(
      tags$style(HTML("
        .main-sidebar {
          background: linear-gradient(180deg,#6a11cb 0%,#2575fc 100%);
          color: white;
        }
        .sidebar-menu > li > a {
          color: white;
          font-weight: 500;
          font-size: 15px;
        }
        .sidebar-menu > li > a:hover {
          background-color: rgba(255,255,255,0.15);
        }
        .sidebar-menu .active > a {
          background-color: rgba(255,255,255,0.25)!important;
        }
      "))
    ),
    sidebarMenu(
      menuItem("📁 Data Upload",         tabName = "upload",           icon = icon("folder-open")),
      menuItem("🧩 Clustering",          tabName = "clustering",       icon = icon("object-group")),
      menuItem("📊 Entropy Calculation", tabName = "entropy",          icon = icon("chart-area")),
      menuItem("🌐 Dispersion Modeling", tabName = "dispersion",       icon = icon("project-diagram")),
      menuItem("🧬 Pathway Diversity",   tabName = "pathway",          icon = icon("dna")),
      menuItem("⚖️ EPS Integration",     tabName = "eps",              icon = icon("balance-scale")),
      menuItem("🎨 EPS Visualization",   tabName = "visualization",    icon = icon("palette")),
      menuItem("📈 Benchmarking",        tabName = "benchmark",        icon = icon("chart-line")),
      menuItem("🧬 Multiome Priming",    tabName = "multiome_priming", icon = icon("vial")),
      menuItem("📊 Multiome Visualization", tabName = "multiome_visualization", icon = icon("chart-bar"))
    )
  ),
  
  # ---------- BODY ----------
  dashboardBody(
    tags$head(
      tags$style(HTML("
        body {background-color:#f7f8fa; font-family:'Segoe UI',sans-serif;}
        .content-wrapper {background-color:#f7f8fa;}
        .main-card {
          background-color:#ffffff;
          border-radius:10px;
          box-shadow:0 0 8px rgba(0,0,0,0.08);
          padding:25px; margin:20px;
        }
        .summary-box {
          background-color:#f0f4ff;
          border-left:5px solid #2575fc;
          border-radius:6px;
          padding:12px 15px;
          margin-top:10px;
          font-size:14px;
        }
        .btn-primary {
          background-color:white!important;
          color:#2575fc!important;
          border-color:#2575fc!important;
          font-weight:500;
          border-radius:6px;
        }
        h2,h3,h4 {font-weight:600;}
      "))
    ),
    
    tabItems(
      # ========== TAB 1: DATA UPLOAD ==========
      tabItem(
        tabName = "upload",
        fluidRow(
          box(
            title = "📁 Upload Your Data",
            width = 4,
            solidHeader = TRUE,
            status = "primary",
            
            fileInput(
              "datafile",
              "Upload Input Data (.rds)",
              accept = ".rds"
            ),
            
            tags$div(
              style = "font-size: 13px; color: #555; margin-top: 8px;",
              tags$p(tags$b("EPS modules:"), " Seurat object containing RNA or SCT assay."),
              tags$p(tags$b("Multiome Priming:"), " Prepared RNA–ATAC multiome object with matched genes and cells.")
            ),
            
            tags$div(
              class = "help-text",
              "⚠️ Upload either a preprocessed Seurat object for RNA modules or a prepared multiome object for Multiome Priming.",
              tags$ul(
                tags$li("RNA modules: Seurat object with RNA or SCT assay"),
                tags$li("Multiome module: list object with RNA and ATAC matrices"),
                tags$li("Matched cell and gene names are required for multiome analysis")
              )
            ),
            
            br(),
            actionButton(
              "load_btn",
              "🚀 Load Dataset",
              class = "btn-primary"
            )
          ),
          
          box(
            title = "Welcome to Evoscope!",
            width = 8,
            solidHeader = TRUE,
            status = "info",
            
            p("EVOscope is an integrated framework for quantifying transcriptional adaptability and regulatory priming at single-cell resolution. The Evolutionary Potential Score (EPS) module analyzes scRNA-seq data to measure cellular adaptability, whereas the Multiome Priming module analyzes paired RNA and chromatin accessibility profiles to identify regulatory discordance and primed cellular states."),
            
            br(),
            h4("📊 Dataset Summary"),
            
            tags$div(
              class = "summary-box",
              verbatimTextOutput("data_summary")
            ),
            
            br(),
            
            tags$div(
              style = "color:#888;font-size:13px;",
              "Note: RNA modules require Seurat objects. Multiome Priming requires a prepared list containing RNA and ATAC matrices."
            )
          )
        )
      ),
      
      # ========== TAB 2: CLUSTERING ==========
      tabItem(
        tabName = "clustering",
        clustering_ui("cluster_module")
      ),
      
      # ========== TAB 3: ENTROPY ==========
      tabItem(
        tabName = "entropy",
        entropy_ui()
      ),
      
      # ========== TAB 4: DISPERSION ==========
      tabItem(
        tabName = "dispersion",
        dispersion_ui("dispersion_module")
      ),
      
      # ========== TAB 5: PATHWAY DIVERSITY ==========
      tabItem(
        tabName = "pathway",
        pathway_ui("pathway_module")
      ),
      
      # ========== TAB 6: EPS INTEGRATION ==========
      tabItem(
        tabName = "eps",
        div(
          class = "main-card",
          eps_ui("eps_module")
        )
      ),
      
      # ========== TAB 7: EPS VISUALIZATION ==========
      tabItem(
        tabName = "visualization",
        div(
          class = "main-card",
          eps_visualization_ui("visual_module")
        )
      ),
      
      # ========== TAB 8: BENCHMARKING ==========
      tabItem(
        tabName = "benchmark",
        div(
          class = "main-card",
          benchmark_ui("benchmark_module")
        )
      ),
      
      # ========== TAB 9: MULTIOME PRIMING ==========
      tabItem(
        tabName = "multiome_priming",
        div(
          class = "main-card",
          multiome_priming_ui("multiome_priming_module")
        )
      ),
      
      # ========== TAB 10: MULTIOME VISUALIZATION ==========
      tabItem(
        tabName = "multiome_visualization",
        div(
          class = "main-card",
          multiome_visualization_ui("multiome_vis_module")
        )
      )
    )
  )
)