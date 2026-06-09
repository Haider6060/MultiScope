# ============================================================
# 🧬 EVOscope Multiome Priming Module — User Interface
# ============================================================

multiome_priming_ui <- function(id) {
  ns <- NS(id)
  
  tagList(
    fluidRow(
      
      # ---------- Left Column: Input & Settings ----------
      column(
        6,
        box(
          title = "🧬 Multiome Priming Analysis",
          status = "primary",
          solidHeader = TRUE,
          width = 12,
          
          h4("Input Requirement"),
          p("This module uses a prepared multiome object containing:"),
          tags$ul(
            tags$li("RNA matrix"),
            tags$li("ATAC gene activity matrix"),
            tags$li("Matched cell names"),
            tags$li("Matched gene names")
          ),
          
          tags$hr(),
          
          h4("Analysis Settings"),
          
          selectInput(
            ns("cor_method"),
            "Correlation Method",
            choices = c("Spearman" = "spearman", "Pearson" = "pearson"),
            selected = "spearman"
          ),
          
          numericInput(
            ns("high_cutoff"),
            "High-Priming Cutoff (%)",
            value = 90,
            min = 50,
            max = 99,
            step = 1
          ),
          
          numericInput(
            ns("low_cutoff"),
            "Low-Priming Cutoff (%)",
            value = 10,
            min = 1,
            max = 50,
            step = 1
          ),
          
          tags$hr(),
          
          actionButton(
            ns("run_priming"),
            "Run Multiome Priming Analysis",
            style = "background-color:white; color:black; border:1px solid #ccc;"
          )
        )
      ),
      
      # ---------- Right Column: Status & Downloads ----------
      column(
        6,
        box(
          title = "📌 Analysis Status",
          status = "info",
          solidHeader = TRUE,
          width = 12,
          
          h4("Module Output"),
          p("After running, this module will calculate:"),
          tags$ul(
            tags$li("Cell Priming Score"),
            tags$li("Gene Priming Score"),
            tags$li("High-Priming cells"),
            tags$li("Low-Priming cells"),
            tags$li("Genes enriched in High-Priming cells")
          ),
          
          tags$hr(),
          
          h4("Status / Progress"),
          verbatimTextOutput(ns("priming_status")),  # Live status
          uiOutput(ns("progress_text")),            # Optional, for HTML styled messages
          
          tags$hr(),
          
          downloadButton(
            ns("download_cell_scores"),
            "Download Cell Priming Scores",
            style = "background-color:white; color:black; border:1px solid #ccc;"
          ),
          
          tags$br(),
          tags$br(),
          
          downloadButton(
            ns("download_gene_scores"),
            "Download Gene Priming Scores",
            style = "background-color:white; color:black; border:1px solid #ccc;"
          )
        )
      )
    )
  )
}