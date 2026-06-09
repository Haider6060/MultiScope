# ============================================================
# 🖌 EVOscope Multiome Priming Visualization — UI
# ============================================================

multiome_visualization_ui <- function(id) {
  ns <- NS(id)
  
  tagList(
    fluidRow(
      
      column(
        4,
        box(
          title = "📊 Visualization Settings",
          status = "primary",
          solidHeader = TRUE,
          width = 12,
          
          selectInput(
            ns("plot_type"),
            "Select Plot Type",
            choices = c(
              "Cell Priming UMAP" = "priming_umap",
              "Cell Priming Score Distribution" = "cell_distribution",
              "Top 20 Primed Cells" = "top_cells",
              "Top 20 Primed Genes" = "top_genes",
              "High vs Low Priming Boxplot" = "high_low_boxplot",
              "Priming Group Counts" = "group_counts",
              "Top Pathway Enrichment" = "primed_pathways"
            ),
            selected = "priming_umap"
          ),
          
          tags$hr(),
          
          actionButton(
            ns("generate_plot"),
            "Generate Plot",
            style = "background-color:white; color:black; border:1px solid #ccc;"
          ),
          
          tags$br(),
          tags$br(),
          
          downloadButton(
            ns("download_plot"),
            "Download Plot PNG 1000 dpi",
            style = "background-color:white; color:black; border:1px solid #ccc;"
          ),
          
          tags$hr(),
          
          h4("Status"),
          verbatimTextOutput(ns("plot_status"))
        )
      ),
      
      column(
        8,
        box(
          title = "📈 Plot Preview",
          status = "info",
          solidHeader = TRUE,
          width = 12,
          plotOutput(ns("priming_plot"), height = "550px")
        )
      )
    )
  )
}