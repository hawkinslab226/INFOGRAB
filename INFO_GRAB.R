
#setwd("~/Documents/Internship/InfoGrab")

required_packages <- c(
  "shiny", "shinyBS", "shinythemes", "shinycustomloader", "shinycssloaders", "ggplot2", "pheatmap", "dplyr", "readr", "readxl", 
  "pheatmap", "DT", "stringr", "purrr", "tibble", "reshape2", 
  "RColorBrewer", "gplots", "scales", "data.table", "gridExtra", 
  "plotly", "lubridate", "shinyjs", "shinydashboard", 
  "shinycssloaders", "shinyWidgets", "htmltools", "tools", 
  "htmlwidgets", "DESeq2", "ggplotify", "grid", "graphics", "igraph", "seriation", "GGally"
  )

install_if_missing <- function(pkg) {
  if (!require(pkg, character.only = TRUE)) {
    install.packages(pkg, dependencies = TRUE)
    }
  }

lapply(required_packages, install_if_missing)

if (!require(DESeq2)) {
  BiocManager::install("DESeq2", force = TRUE)
}

if (!require(pcaExplorer)) {
  BiocManager::install("pcaExplorer", force = TRUE)
}

lapply(required_packages, library, character.only = TRUE)

library(shiny)
library(shinyBS)
library(shinythemes)
library(shinyjs)
library(shinydashboard)
library(shinycssloaders)
library(shinycustomloader)
library(dplyr)
library(ggplot2)
library(ggrepel)
library(pheatmap)
library(DT)
library(stringr)
library(viridis)
library(readr)
library(tidyr)
library(tibble)
library(purrr)
library(magrittr)
library(htmltools)
library(scales)
library(plotly)
library(DESeq2)
library(pcaExplorer)
library(ggplotify)
library(htmlwidgets)
library(tools)
library(grid)
library(graphics)
library(igraph)
library(GGally)

if (!exists("comparison_results")) {
  comparison_results <- readRDS("data/comparison_results.rds")
}

if (!exists("RPKM_data")) {
  RPKM_data <- read.delim("data/logRPKM0625_filter.txt")
}

if (!exists("phenodata")) {
  phenodata <- read.csv("data/Pheno_data071523.csv")
}

if (!exists("cnt_data")) {
  cnt_data <- read.delim("data/all_cnt_0625.txt")
}

mycolors1 <- c(
  `Kidney`="#43009A",
  `Trachea`="#990099",
  `BCell`="#0DFAFA",
  `TCell(spleen)`="#13B7B7",
  `Bursa`="#004C99",
  `Thymus`="#2685E4",
  `Macrophage(D0)`="#F02B6D",
  `Macrophage(D3)`="#FF0091",
  `Macrophage(D6)`="#F572BC",
  `Macrophage(lung)`="#D06AAA",
  `Monocyte(blood)`="#91155B",
  `Ileum`="#98D55A",
  `Jejunum`="#4C9900",
  `Proximal.Cecum`="#CCFF99",
  `Iliotibial.major`="#A1122A",
  `Pectoralis.major`="#FFCCCC",
  `Isthmus`="#FF8000",
  `Magnum`="#DE5100",
  `Shell.Gland`="#FFC78E",
  `Ovary`="#CCCC00"
)

############################################################

# INFO GRAB

# Cite BioRender
# Cite pcaExplorer

tissue_names <- unique(phenodata$Tissue)

systems <- unique(phenodata$System)
tissue_by_system <- lapply(systems, function(system) {
  unique(phenodata$Tissue[phenodata$System == system])
})
names(tissue_by_system) <- systems


ui <- navbarPage(
  title = tags$a(
    tags$div(
      tags$img(src = "UW_logo.png", height = "50px", style = "padding: 5px; margin-top: -17px;"),
      tags$img(src = "WesternU_logo.png", height = "60px", style = "padding: 5px; margin-top: -19px;")
    )
  ),
  windowTitle = "INFO GRAB",
  theme = shinytheme("paper"),
  id = "navbar",
  useShinyjs(),
  tags$head(
    tags$link(rel = "stylesheet", href = "https://fonts.googleapis.com/css2?family=Sen:wght@400;700&display=swap"),
    tags$link(rel = "stylesheet", href = "https://fonts.googleapis.com/css2?family=Ubuntu:wght@400;700&display=swap"),
    tags$script(HTML("
      $(document).on('shiny:connected', function(event) {
        Shiny.setInputValue('startAnalysis', true);
      });
    ")),
    tags$style(HTML("
      .navbar-nav > li > a, .navbar-brand {
        font-family: 'Sen', sans-serif !important;
        font-size: 15px;
      }
      h1 {
        font-family: 'Ubuntu', sans-serif !important;
        font-size: 40px !important;
        font-weight: bold !important;
      }
      h2 {
        font-family: 'Ubuntu', sans-serif !important;
        font-size: 34px !important;
        font-weight: bold !important;
      }
      h3 {
        font-family: 'Ubuntu', sans-serif !important;
        font-size: 28px !important;
        font-weight: bold !important;
      }
      h4 {
        font-family: 'Sen', sans-serif !important;
        font-size: 20px !important;
      }
      .form-group > label {
        font-family: 'Sen', sans-serif !important;
      }
      .navbar-nav > li > a, .navbar-brand {
        font-family: 'Sen', sans-serif !important;
        font-size: 15px;
        color: white !important;
      }
      .navbar-default {
        background-color: #4B0082 !important;
        border-color: #4B0082 !important;
      }
      .navbar-default .navbar-nav > .active > a, 
      .navbar-default .navbar-nav > .active > a:focus, 
      .navbar-default .navbar-nav > .active > a:hover {
        color: white !important;
        background-color: grey !important;
      }
      .btn {
        background-color: #4B0082 !important;
        color: white !important;
      }
      .btn:hover {
        background-color: #4B0082 !important;
        color: white !important;
      }
      #info-grab-link:hover {
        background-color: grey !important;
        color: white !important;
      }
      #progress-container {
        position: fixed;
        bottom: 10px;
        right: 10px;
        width: 300px;
        z-index: 1050;
      }
    "))
  ),
  tabPanel("INFO GRAB",
           fluidPage(
             div(style = "text-align: center;",
                 h2("INFO GRAB"),
                 h3("INteractive Functional Ontology Genomic Regulatory Element Analysis and Browsing"),
                 tags$hr(style = "height:1px; border:none; color:#300; background-color:#300;"),
                 h5("Welcome to INFO GRAB. This app allows you to perform differential expression analysis and 
                   visualize the results through various plots and tables, looking at tissues in the chicken (Gallus gallus)."),
                 tags$hr(style = "height:1px; border:none; color:#300; background-color:#300;"),
                 img(src = "chicken.png", id = "chicken-image", style = "width: 40%; max-width: 100%; height: auto;")
             ),
             uiOutput("tissue_display")
           )
  ),
  
  tabPanel("Analysis",
           tabsetPanel(
             tabPanel("Home",
                      sidebarLayout(
                        sidebarPanel(
                          selectInput("tissue_select1", "Tissue 1", choices = tissue_names),
                          selectInput("tissue_select2", "Tissue 2", choices = tissue_names[-1]),
                          tags$hr(style = "height:1px; border:none; color:#300; background-color:#300;"),
                          selectInput("selected_fdr", "FDR", 
                                      choices = c("No filter" = NA, "0.01" = 0.01, "0.05" = 0.05, "0.1" = 0.1), 
                                      selected = "0.05"),
                          sliderInput("selected_fc", "FC (log2 fold change)", min = 0, max = 20, step = 1, value = 0),
                          tags$hr(style = "height:1px; border:none; color:#300; background-color:#300;"),
                          sliderInput("num_of_digits", "Number of digits", value = 5, min = 2, max = 20),
                          tags$hr(style = "height:1px; border:none; color:#300; background-color:#300;"),
                          downloadButton("download_DEA_data", "Download CSV"),
                          tags$hr(style = "height:1px; border:none; color:#300; background-color:#300;"),
                          tags$p("FDR: The expected proportion of false discoveries among the rejected hypotheses."),
                          tags$p("FC: A measure describing how much a quantity changes between an original and a subsequent measurement.")
                        ),
                        mainPanel(
                          uiOutput("title"),
                          
                          withLoader(tableOutput("summary_info"), type = "html", loader = "dnaspin"),
                          br(),
                          textInput("searched_gene", "Search for genes (comma-separated)", width = '600px'),
                          actionButton("search_gene", "Search"),
                          actionButton("random_gene", "Search 10 Random Genes"),
                          actionButton("clear_search_main", "Clear"),
                          DTOutput("searched_gene_info"),
                          br(),
                          uiOutput("collapsePanels")
                        )
                      )
             ),
             
             tabPanel("Correlation Plot",
                      sidebarLayout(
                        sidebarPanel(
                          radioButtons("correlation_mode", "Select Mode:",
                                       choices = c("Selected Tissues" = "selected", "Tissues by System" = "system"),
                                       selected = "selected"),
                          conditionalPanel(
                            condition = "input.correlation_mode == 'system'",
                            selectInput("correlation_system", "Select System:", choices = systems)
                          )
                        ),
                        mainPanel(
                          withLoader(plotOutput("correlation_plot", width = "90%", height = "750px"), type = "html", loader = "dnaspin")
                        )
                      )
             ),
             
             tabPanel("Volcano Plot",
                      sidebarLayout(
                        sidebarPanel(
                          textInput("searched_gene_volcano", "Search for gene:", width = '600px'),
                          actionButton("search_gene_volcano", "Search"),
                          actionButton("clear_search_gene_volcano", "Clear"),
                          tags$hr(style = "height:1px; border:none; color:#300; background-color:#300;"),
                          sliderInput("n_labels", "Number of Labels", min = 0, max = 100, value = 10),
                          tags$hr(style = "height:1px; border:none; color:#300; background-color:#300;"),
                          selectInput("volcano_background_col", "Background Color", 
                                      choices = c("White" = "White","Grey" = "Grey"), selected = "White"),
                          tags$hr(style = "height:1px; border:none; color:#300; background-color:#300;"),
                          downloadButton("download_volcano_plot", "Download PNG"),
                          tags$hr(style = "height:1px; border:none; color:#300; background-color:#300;"),
                          tags$p("Each point represents a gene, plotted by its log2 fold change and the negative 
                                 logarithm of its adjusted p-value. Genes with significant changes (high fold change 
                                 and low p-value) are highlighted."),
                          br(),
                          tags$p("Genes farthest from the origin along the x-axis (log2 fold change) have the most 
                                 significant changes in expression. The y-axis shows the significance level, with points 
                                 higher up indicating more significant changes.")
                        ),
                        mainPanel(
                          uiOutput("title_volcano_plot"),
                          withLoader(plotOutput("volcano_plot", width = "90%", height = "750px"), type = "html", loader = "dnaspin")
                        )
                      )
             ),
             tabPanel("Heatmap",
                      sidebarLayout(
                        sidebarPanel(
                          selectInput("sort_by", "Sort by:", choices = c("FDR" = "FDR", "FC" = "FC")),
                          selectInput("top_n", "Number of Top Genes:", 
                                      choices = c("5" = 5, "10" = 10, "20" = 20, "50" = 50, "100" = 100, "500" = 500, "1000" = 1000, "All DE Genes" = Inf), 
                                      selected = "10"),
                          tags$hr(style = "height:1px; border:none; color:#300; background-color:#300;"),
                          textInput("search_genes_heatmap", "Search for Specific Genes (comma-separated):"),
                          actionButton("search_gene_heatmap", "Search"),
                          actionButton("random_gene_heatmap", "Search 10 Random Genes"),
                          actionButton("clear_search", "Clear"),
                          tags$hr(style = "height:1px; border:none; color:#300; background-color:#300;"),
                          selectInput("color_scale", "Select Color Scale", 
                                      choices = c("Blue-White-Red" = "blue_white_red", 
                                                  "Green-Black-Red" = "green_black_red",
                                                  "Purple-White-Green" = "purple_white_green",
                                                  "Viridis" = "viridis",
                                                  "Plasma" = "plasma",
                                                  "Cividis" = "cividis",
                                                  "Inferno" = "inferno")),
                          tags$hr(style = "height:1px; border:none; color:#300; background-color:#300;"),
                          downloadButton("download_heatmap_plot", "Download PNG"),
                          tags$hr(style = "height:1px; border:none; color:#300; background-color:#300;"),
                          tags$p("Each cell represents the expression level of a gene in a sample, with colors indicating 
                                 the level of expression."),
                          br(),
                          tags$p("The color scale is based on z-scores, which represent the number of standard deviations 
                                 away from the mean expression level. High positive z-scores indicate higher expression levels, 
                                 while high negative z-scores indicate lower expression levels.")
                        ),
                        mainPanel(
                          uiOutput("title_heatmap"),
                          withLoader(plotOutput("heatmap_plot", width = "80%", height = "750px"), type = "html", loader = "dnaspin")
                        )
                      )
             ),
             
             tabPanel("PCA",
                      sidebarLayout(
                        sidebarPanel(),
                        mainPanel()
                      )),
             
             tabPanel("Boxplot",
                      sidebarLayout(
                        sidebarPanel(
                          textInput("boxplot_searched_gene", "Search for genes (comma-separated)", width = '600px'),
                        ),
                        mainPanel()
                      )),
             
             tabPanel("Help",
                      mainPanel(
                        h2("Analysis"),
                        h3("Home"),
                        p("In the 'Home' tab, you can select two tissues to compare and adjust the analysis parameters."),
                        tags$ul(
                          tags$li(tags$b("Tissue 1 and Tissue 2:"), " Select the tissues you want to compare."),
                          tags$li(tags$b("FDR:"), " Choose the false discovery rate threshold for filtering the results."),
                          tags$li(tags$b("FC (log2 fold change):"), " Set the fold change threshold for filtering the results."),
                          tags$li(tags$b("Number of digits:"), " Specify the number of decimal places for displaying values."),
                          tags$li(tags$b("Download CSV:"), " Download the filtered data as a CSV file."),
                          tags$li(tags$b("Search for genes:"), " Enter genes to search in the results."),
                          tags$li(tags$b("Search 10 Random Genes:"), " Display information for 10 randomly genes.")
                        ),
                        h3("Volcano Plot"),
                        p("The 'Volcano Plot' tab allows you to visualize the results of the differential expression analysis."),
                        tags$ul(
                          tags$li(tags$b("Search for gene:"), " Enter a gene name to highlight it on the plot."),
                          tags$li(tags$b("Number of Labels:"), " Specify the number of top genes to label on the plot."),
                          tags$li(tags$b("Download PNG:"), " Download the volcano plot as a PNG file.")
                        ),
                        h3("Heatmap"),
                        p("The 'Heatmap' tab provides a heatmap visualization of the expression levels of genes."),
                        tags$ul(
                          tags$li(tags$b("Sort by:"), " Choose to sort the genes by FDR or fold change. 
                                  Based on the selection, the top genes that value will be displayed."),
                          tags$li(tags$b("Number of Top Genes:"), " Select the number of top genes to include in the heatmap."),
                          tags$li(tags$b("Search for Specific Genes:"), " Enter gene names (comma-separated) to highlight on the heatmap."),
                          tags$li(tags$b("Select Color Scale:"), " Choose the color scale for the heatmap."),
                          tags$li(tags$b("Download PNG:"), " Download the heatmap as a PNG file.")
                        ),
                        h3("Parameters Explanation"),
                        tags$p("False Discovery Rate (FDR): The expected proportion of false discoveries among the rejected hypotheses."),
                        tags$p("Fold Change (FC): A measure describing how much a quantity changes between an original and a subsequent measurement.")
                      )
             )
           )
  )
)

server <- function(input, output, session) {
  
  if (!exists("comparison_results")) {
    comparison_results <- readRDS("data/comparison_results.rds")
  }
  
  if (!exists("RPKM_data")) {
    RPKM_data <- read.delim("data/logRPKM0625_filter.txt")
  }
  
  if (!exists("phenodata")) {
    phenodata <- read.csv("data/Pheno_data071523.csv")
  }
  
  if (!exists("cnt_data")) {
    cnt_data <- read.delim("data/all_cnt_0625.txt")
  }
  
  observeEvent(input$startAnalysis, {
    runjs('$("#progress-container").show();')
    runjs('$("#progress-bar .progress-bar").css("width", "0%");')
    runjs('$("#progress-text").text("Initializing analysis...");')
    
    withProgress(message = 'Loading...', value = 0, {
      incProgress(0.1, detail = "Updating tissue select inputs...")
      runjs('$("#progress-bar .progress-bar").css("width", "10%");')
      updateSelectInput(session, "tissue_select2", 
                        choices = tissue_names[tissue_names != input$tissue_select1], 
                        selected = isolate(input$tissue_select2))
      
      updateSelectInput(session, "tissue_select1", 
                        choices = tissue_names[tissue_names != input$tissue_select2], 
                        selected = isolate(input$tissue_select1))
      
      incProgress(0.5, detail = "Processing data...")
      runjs('$("#progress-bar .progress-bar").css("width", "30%");')
      filtered_data_combined()
      expressed_genes_by_tissue()
      
      
      incProgress(0.1, detail = "Finalizing...")
      runjs('$("#progress-bar .progress-bar").css("width", "100%");')
      runjs('$("#progress-container").hide();')
    })
  })
  
  # INFO GRAB Tab
  
  output$tissue_display <- renderUI({
    lapply(names(tissue_by_system), function(system) {
      list(
        h3(paste(system, "System")),
        div(style = "display: flex; flex-wrap: wrap; gap: 15px; align-items: center;",
            lapply(tissue_by_system[[system]], function(tissue) {
              div(style = "width: calc(33.33% - 10px); display: flex; align-items: center;",
                  {
                    img_path <- paste0(gsub(" ", "_", tissue), ".png")
                    img_src <- file.path("www", img_path)
                    if (file.exists(img_src)) {
                      img(src = img_path, style = "width: 30px; height: auto; margin-right: 5px;")
                    }
                  },
                  actionLink(inputId = paste0("tissue_", gsub(" ", "_", tissue)), 
                             label = tissue, 
                             style = "font-size: 18px; text-align: left;")
              )
            })
        ),
        tags$hr(style = "height:1px; border:none; color:#300; background-color:#300;")
      )
    })
  })
  
  observe({
    lapply(names(tissue_by_system), function(system) {
      lapply(tissue_by_system[[system]], function(tissue) {
        observeEvent(input[[paste0("tissue_", gsub(" ", "_", tissue))]], {
          samples <- phenodata %>% filter(Tissue == tissue) %>% pull(Sample_name)
          samples_display <- phenodata %>% filter(Tissue == tissue) %>% pull(Sample_name2)
          expressed_genes_count <- expressed_genes_by_tissue() %>%
            filter(Tissue == tissue) %>%
            pull(Expressed_Genes)
          
          tissue_data <- RPKM_data %>% select(samples)
          
          gene_expression_sums <- tissue_data %>% 
            rowSums(na.rm = TRUE) %>%
            as.data.frame() %>%
            setNames("Expression_Sum") %>%
            rownames_to_column(var = "Gene") %>%
            arrange(desc(Expression_Sum))
          
          top_genes <- gene_expression_sums %>% head(5)
          
          showModal(modalDialog(
            title = paste(tissue, "Information"),
            HTML(
              paste(
                "Samples: ", paste(samples_display, collapse = ", "), "<br><hr>",
                "Total Expressed Genes: ", expressed_genes_count, "<br><hr>",
                "Top 5 Most Expressed Genes (RPKM):<br>",
                paste(top_genes$Gene, ": ", round(top_genes$Expression_Sum, 2), sep = "", collapse = "<br>")
              )
            ),
            easyClose = TRUE,
            footer = NULL
          ))
        }, ignoreInit = TRUE)
      })
    })
  })
  
  # Data
  
  observe({
    selected_tissue1 <- input$tissue_select1
    updateSelectInput(session, "tissue_select2", 
                      choices = tissue_names[tissue_names != selected_tissue1], 
                      selected = isolate(input$tissue_select2))
  })
  
  observe({
    selected_tissue2 <- input$tissue_select2
    updateSelectInput(session, "tissue_select1", 
                      choices = tissue_names[tissue_names != selected_tissue2], 
                      selected = isolate(input$tissue_select1))
  })
  
  order_val <- reactiveVal(0)
  
  unfiltered_data <- reactive({
    req(input$tissue_select1, input$tissue_select2)
    
    if (paste("Tissue", input$tissue_select1, "vs", input$tissue_select2, sep = "_") %in% names(comparison_results)) {
      comparison_list <- comparison_results[[paste("Tissue", input$tissue_select1, "vs", input$tissue_select2, sep = "_")]]
      order_val(0)
    } 
    else if (paste("Tissue", input$tissue_select2, "vs", input$tissue_select1, sep = "_") %in% names(comparison_results)) {
      comparison_list <- comparison_results[[paste("Tissue", input$tissue_select2, "vs", input$tissue_select1, sep = "_")]]
      order_val(1)
    }
    
    comparison_list <- as.data.frame(comparison_list)
    comparison_list$Gene <- rownames(comparison_list)
    
    
    comparison_list <- comparison_list %>%
      mutate(padj = case_when(
        padj == 0 ~ .Machine$double.xmin,
        TRUE ~ padj
      ))
    
    comparison_list
    
  })
  
  filtered_data_fdr <- reactive({
    data <- unfiltered_data()
    req(data)
    
    if (!is.na(as.numeric(input$selected_fdr))) {
      data <- data %>% filter(padj <= as.numeric(input$selected_fdr))
    }
    
    data
  })
  
  filtered_data_fc <- reactive({
    data <- unfiltered_data()
    req(data)
    
    if (!is.na(as.numeric(input$selected_fc))) {
      data <- data %>% filter(abs(log2FoldChange) >= as.numeric(input$selected_fc))
    }
    
    data
  })
  
  filtered_data_combined <- reactive({
    
    data <- filtered_data_fdr()
    req(data)
    
    if (as.numeric(input$selected_fc) != 0) {
      data <- data %>% filter(abs(log2FoldChange) >= as.numeric(input$selected_fc))
    }
    
    data
  })
  
  
  output$title <- renderUI({
    req(input$tissue_select1, input$tissue_select2)
    order <- order_val()
    
    if (order == 0) {
      return(h3(paste(input$tissue_select1, " (+) vs. ", input$tissue_select2, " (-)", sep = "")))
    }
    
    if (order == 1) {
      return(h3(paste(input$tissue_select2, " (+) vs. ", input$tissue_select1, " (-)", sep = "")))
    }
    
  })
  
  # Expressed genes counts
  expressed_genes_by_tissue <- reactive({
    # May need to change (or make user able to modify)
    threshold <- 1  
    expressed_genes <- RPKM_data %>%
      rowwise() %>%
      mutate(across(everything(), ~ . > threshold)) %>%
      summarise(across(everything(), sum))
    
    tissue_counts <- phenodata %>%
      filter(Sample_name %in% colnames(RPKM_data)) %>%
      group_by(Tissue) %>%
      summarise(Expressed_Genes = sum(colnames(RPKM_data) %in% Sample_name * expressed_genes))
    
    tissue_counts
  })
  
  
  # Summary of data
  output$summary_info <- renderUI({
    data <- filtered_data_combined()
    order <- order_val()
    req(order)
    data_padj <- filtered_data_fdr()
    data_fc <- filtered_data_fc()
    unfiltered <- unfiltered_data()
    req(data, unfiltered)
    
    total_genes <- nrow(unfiltered)
    total_significant_genes <- nrow(data_padj)
    total_fc_genes <- nrow(data_fc)
    total_combined_genes <- nrow(data)
    
    if (order == 0) {
      upregulated_genes_tissue_1 <- sum(data$log2FoldChange > 0)
      upregulated_genes_tissue_2 <- sum(data$log2FoldChange < 0)
    } else {
      upregulated_genes_tissue_1 <- sum(data$log2FoldChange < 0)
      upregulated_genes_tissue_2 <- sum(data$log2FoldChange > 0)
    }
    
    mean_log2fc <- mean(data$log2FoldChange)
    mean_adj_pvalue <- format(mean(data$padj, na.rm = TRUE), scientific = T, digits = input$num_of_digits)
    
    html_output <- HTML(paste(
      "<table style='width: 50%;'>",
      "<tr><th></th><th></th></tr>",
      "<tr><td>Total genes</td><td>", total_genes, "</td></tr>",
      "<tr><td>Total significant genes</td><td>", total_significant_genes, "</td></tr>",
      "<tr><td>FC > parameter</td><td>", total_fc_genes, "</td></tr>",
      "<tr><td>Significant and FC > parameter</td><td>", total_combined_genes, "</td></tr>",
      "<tr><td colspan='2'><hr></td></tr>",
      "<tr><td>Upregulated genes in ", input$tissue_select1, "</td><td>", upregulated_genes_tissue_1, "</td></tr>",
      "<tr><td>Upregulated genes in ", input$tissue_select2, "</td><td>", upregulated_genes_tissue_2, "</td></tr>",
      "<tr><td colspan='2'><hr></td></tr>",
      "<tr><td>Mean log2 Fold Change</td><td>", round(mean_log2fc, 2), "</td></tr>",
      "<tr><td>Mean Adjusted p-value</td><td>", formatC(mean_adj_pvalue, format = "f", digits = 5), "</td></tr>",
      "<tr><td colspan='2'><hr></td></tr>",
      "</table>"
    ))
    
    return(html_output)
  })
  
  # Search for genes
  
  gene_search_result <- reactiveVal()
  
  observeEvent(input$search_gene, {
    data <- unfiltered_data()
    req(data)
    
    searched_genes <- strsplit(input$searched_gene, ",\\s*")[[1]]
    
    if (length(searched_genes) > 0) {
      gene_data <- data %>% filter(toupper(Gene) %in% toupper(searched_genes))
      not_found_genes <- setdiff(toupper(searched_genes), toupper(gene_data$Gene))
      
      if (length(not_found_genes) > 0) {
        showNotification(paste("The following genes were not found:", paste(not_found_genes, collapse = ", ")), type = "warning")
      }
      
      if (nrow(gene_data) == 0) {
        gene_search_result(data.frame(Statistic = "Genes not found.", Value = ""))
      } else {
        gene_search_result(gene_data)
      }
    }
  })
  
  observeEvent(input$clear_search_main, {
    gene_search_result(NULL)
    updateTextInput(session, "searched_gene", value = "")
  })
  
  
  observeEvent(input$random_gene, {
    data <- unfiltered_data()
    req(data)
    
    random_genes <- sample(data$Gene, 10)
    updateTextInput(session, "searched_gene", value = paste(random_genes, collapse = ", "))
    gene_data <- data %>% filter(Gene %in% random_genes)
    gene_search_result(gene_data)
  })
  
  output$searched_gene_info <- renderDT({
    gene_data <- gene_search_result()
    
    if (is.null(gene_data)) {
      return(NULL)
    } else if ("Statistic" %in% colnames(gene_data) && gene_data$Statistic == "Genes not found.") {
      return(datatable(gene_data, 
                       options = list(dom = 't', paging = FALSE, ordering = FALSE), rownames = FALSE))
    } else {
      gene_data <- gene_data %>%
        arrange(padj) %>%
        mutate(padj = format(padj, scientific = T, digits = input$num_of_digits), 
               log2FoldChange = format(log2FoldChange, scientific = T, digits = input$num_of_digits)) %>%
        select(Gene, log2FoldChange, padj)
      
      return(datatable(gene_data, 
                       options = list(dom = 't', paging = FALSE, ordering = FALSE), rownames = FALSE))
    }
  })
  
  # Tables
  
  tissue_titles <- reactive({
    if (order_val() == 0) {
      list(
        upregulated = paste("Top 10 Most Upregulated Genes in", input$tissue_select1),
        downregulated = paste("Top 10 Most Upregulated Genes in", input$tissue_select2)
      )
    }
    else {
      list(
        upregulated = paste("Top 10 Most Upregulated Genes in", input$tissue_select2),
        downregulated = paste("Top 10 Most Upregulated Genes in", input$tissue_select1)
      )
    }
  })
  
  
  output$collapsePanels <- renderUI({
    bsCollapse(id = "collapseTables", open = "",
               shinyBS::bsCollapsePanel("Ten Smallest p-values", DTOutput("table_smallest_pvalues")),
               shinyBS::bsCollapsePanel("Ten Largest Absolute Fold Changes", DTOutput("table_largest_fc")),
               shinyBS::bsCollapsePanel(title = tissue_titles()$upregulated, DTOutput("table_upregulated")),
               shinyBS::bsCollapsePanel(title = tissue_titles()$downregulated, DTOutput("table_downregulated"))
    )
  })
  
  output$table_smallest_pvalues <- renderDT({
    data <- filtered_data_combined()
    req(data)
    
    if (nrow(data) == 0) {
      return(datatable(data.frame(Statistic = "No results found for the selected criteria.", Value = "")))
    }
    
    smallest_pvalues <- data %>% 
      arrange(padj) %>%
      mutate(padj = format(padj, scientific = T, digits = input$num_of_digits), 
             log2FoldChange = format(log2FoldChange, scientific = T, digits = input$num_of_digits)) %>%
      head(10) %>%
      select(Gene, log2FoldChange, padj)
    
    datatable(smallest_pvalues, 
              options = list(dom = 't', paging = FALSE, ordering = FALSE), rownames = FALSE)
  })
  
  output$table_largest_fc <- renderDT({
    data <- filtered_data_combined()
    req(data)
    
    if (nrow(data) == 0) {
      return(datatable(data.frame(Statistic = "No results found for the selected criteria.", Value = "")))
    }
    
    largest_fc <- data %>% 
      arrange(desc(abs(log2FoldChange))) %>% 
      mutate(padj = format(padj, scientific = T, digits = input$num_of_digits), 
             log2FoldChange = format(log2FoldChange, scientific = T, digits = input$num_of_digits)) %>%
      head(10) %>%
      select(Gene, log2FoldChange, padj)
    
    datatable(largest_fc, options = list(dom = 't', paging = FALSE, ordering = FALSE), rownames = FALSE)
  })
  
  output$table_upregulated <- renderDT({
    data <- filtered_data_combined()
    req(data)
    
    if (nrow(data) == 0) {
      return(datatable(data.frame(Statistic = "No results found for the selected criteria.", Value = "")))
    }
    
    upregulated <- data %>% 
      arrange(desc(log2FoldChange)) %>% 
      mutate(padj = format(padj, scientific = T, digits = input$num_of_digits), 
             log2FoldChange = format(log2FoldChange, scientific = T, digits = input$num_of_digits)) %>%
      head(10) %>%
      select(Gene, log2FoldChange, padj)
    
    datatable(upregulated, options = list(dom = 't', paging = FALSE, ordering = FALSE), rownames = FALSE)
  })
  
  output$table_downregulated <- renderDT({
    data <- filtered_data_combined()
    req(data)
    
    if (nrow(data) == 0) {
      return(datatable(data.frame(Statistic = "No results found for the selected criteria.", Value = "")))
    }
    
    downregulated <- data %>% 
      arrange(log2FoldChange) %>% 
      mutate(padj = format(padj, scientific = T, digits = input$num_of_digits), 
             log2FoldChange = format(log2FoldChange, scientific = T, digits = input$num_of_digits)) %>%
      head(10) %>%
      select(Gene, log2FoldChange, padj)
    
    datatable(downregulated, options = list(dom = 't', paging = FALSE, ordering = FALSE), rownames = FALSE)
  })
  
  output$table_high_expression <- renderDT({
    data <- filtered_data_combined()
    req(data)
    
    if (nrow(data) == 0) {
      return(datatable(data.frame(Statistic = "No results found for the selected criteria.", Value = "")))
    }
    
    high_expression <- data %>% 
      arrange(desc(baseMean)) %>% 
      head(10) %>%
      select(Gene, baseMean, log2FoldChange, padj)
    
    datatable(high_expression, options = list(dom = 't', paging = FALSE, ordering = FALSE), rownames = FALSE)
  })
  
  output$table_smallest_pvalues <- renderDT({
    data <- filtered_data_combined()
    req(data)
    
    if (nrow(data) == 0) {
      return(datatable(data.frame(Statistic = "No results found for the selected criteria.", Value = "")))
    }
    
    smallest_pvalues <- data %>% 
      arrange(padj) %>% 
      mutate(padj = format(padj, scientific = T, digits = input$num_of_digits), 
             log2FoldChange = format(log2FoldChange, scientific = T, digits = input$num_of_digits)) %>%
      head(10) %>%
      select(Gene, log2FoldChange, padj)
    
    datatable(smallest_pvalues, options = list(dom = 't', paging = FALSE, ordering = FALSE), rownames = FALSE)
  })
  
  output$download_DEA_data <- downloadHandler(
    
    filename = function(){"DEA_data.csv"}, 
    content = function(fname){
      data <- filtered_data_combined()
      data <- data %>%
        select(-Gene)
      write.csv(data, fname)
    }
  )
  
  ### Volcano Plot
  
  output$title_volcano_plot <- renderUI({
    req(input$tissue_select1, input$tissue_select2)
    h3(paste("Volcano Plot:", input$tissue_select1, "vs.", input$tissue_select2))
  })
  
  search_gene_volcano <- reactiveVal(NULL)
  
  observeEvent(input$search_gene_volcano, {
    search_gene <- input$searched_gene_volcano
    search_gene_volcano(search_gene)
  })
  
  observeEvent(input$clear_search_gene_volcano, {
    search_gene_volcano(NULL)
    updateTextInput(session, "searched_gene_volcano", value = "")
  })
  
  
  generate_volcano_plot <- function(data) {
    selected_fdr <- as.numeric(input$selected_fdr)
    selected_fc <- as.numeric(input$selected_fc)
    order <- order_val()
    searched_gene <- search_gene_volcano()
    
    data <- data %>%
      mutate(Gene = toupper(Gene),
             color = case_when(
               (is.na(selected_fdr) | padj < selected_fdr) & log2FoldChange > selected_fc & order == 0 ~ "Upregulated in Tissue 1",
               (is.na(selected_fdr) | padj < selected_fdr) & log2FoldChange < -selected_fc & order == 0 ~ "Downregulated in Tissue 2",
               (is.na(selected_fdr) | padj < selected_fdr) & log2FoldChange > selected_fc & order == 1 ~ "Upregulated in Tissue 2",
               (is.na(selected_fdr) | padj < selected_fdr) & log2FoldChange < -selected_fc & order == 1 ~ "Downregulated in Tissue 1",
               TRUE ~ "Non-significant"
             ))
    
    if (order == 0) {
      volcano_plot_labels <- c(
        "Upregulated in Tissue 1" = paste("Upregulated in", input$tissue_select1),
        "Downregulated in Tissue 2" = paste("Upregulated in", input$tissue_select2),
        "Non-significant" = "Non-significant"
      )
      label_1 <- input$tissue_select1
      label_2 <- input$tissue_select2
    } else {
      volcano_plot_labels <- c(
        "Upregulated in Tissue 2" = paste("Upregulated in", input$tissue_select2),
        "Downregulated in Tissue 1" = paste("Upregulated in", input$tissue_select1),
        "Non-significant" = "Non-significant"
      )
      label_1 <- input$tissue_select2
      label_2 <- input$tissue_select1
    }
    
    top_genes <- data %>%
      filter((is.na(selected_fdr) | padj < selected_fdr) & abs(log2FoldChange) >= selected_fc) %>%
      arrange(padj) %>%
      head(input$n_labels)
    
    if (input$volcano_background_col == "White") {
      background <- element_blank()
      label_color <- "black"
      searched_gene_color <- "#000000"
      arrow_color <- "black"
      vline_color <- "darkblue"
      hline_color <- "darkred"
    } else { # "Grey"
      background <- element_rect(fill = "grey37",
                                 colour = "grey37",
                                 size = 0.5, linetype = "solid")
      label_color <- "white"
      searched_gene_color <- "#FFFFFF"
      arrow_color <- "white"
      vline_color <- "skyblue"
      hline_color <- "pink"
    }
    
    p <- ggplot(data, aes(x = log2FoldChange, y = -log10(padj))) +
      geom_point(aes(color = color), alpha = 0.5) +
      scale_color_manual(values = c(
        "Upregulated in Tissue 1" = mycolors1[[input$tissue_select1]],
        "Downregulated in Tissue 2" = mycolors1[[input$tissue_select2]],
        "Upregulated in Tissue 2" = mycolors1[[input$tissue_select2]],
        "Downregulated in Tissue 1" = mycolors1[[input$tissue_select1]],
        "Non-significant" = "grey"
      ), labels = volcano_plot_labels) +
      labs(x = "Log2 Fold Change", y = "-log10 Adjusted P-Value", color = "Expression: ") +
      theme_minimal() +
      theme(
        text = element_text(family = "Sen"),
        legend.position = "bottom",
        legend.box = "horizontal",
        legend.box.just = "center",
        legend.text = element_text(size = 14), 
        legend.title = element_text(size = 16),
        legend.background = element_rect(fill = "#f0f0f0", color = "black"),
        legend.key = element_rect(fill = "white", color = "black"),
        legend.key.size = unit(1.5, "lines"), 
        axis.title = element_text(size = 16),
        axis.text = element_text(size = 14),
        panel.grid.major = element_blank(), 
        panel.grid.minor = element_blank(),
        panel.background = background,
        axis.line = element_line(color = "black"),
      ) +
      guides(color = guide_legend(override.aes = list(shape = 15, size = 5))) +
      geom_text_repel(data = top_genes, aes(label = Gene), max.overlaps = 1000, color = label_color) +
      geom_hline(yintercept = -log10(selected_fdr), linetype = "dashed", color = hline_color, alpha = 1/3) +
      geom_vline(xintercept = c(-selected_fc, selected_fc), linetype = "dashed", color = vline_color, alpha = 1/3) +
      annotate("text", x = 0 + max(data$log2FoldChange) * 0.5, 
               y = max(-log10(data$padj), na.rm = TRUE) * 1.2, 
               label = label_1, size = 5, hjust = 0.5, color = label_color) +
      annotate("text", x = 0 - max(data$log2FoldChange) * 0.5, 
               y = max(-log10(data$padj), na.rm = TRUE) * 1.2, 
               label = label_2, size = 5, hjust = 0.5, color = label_color) +
      annotate("segment", x = 0 + max(data$log2FoldChange) * 0.4, 
               y = max(-log10(data$padj), na.rm = TRUE) * 1.175, 
               xend = 0 + max(data$log2FoldChange) * 0.6, 
               yend = max(-log10(data$padj), na.rm = TRUE) * 1.175,
               arrow = arrow(length = unit(0.3, "cm")), color = arrow_color) +
      annotate("segment", x = 0 - max(data$log2FoldChange) * 0.4, 
               y = max(-log10(data$padj), na.rm = TRUE) * 1.175, 
               xend = 0 - max(data$log2FoldChange) * 0.6, 
               yend = max(-log10(data$padj), na.rm = TRUE) * 1.175,
               arrow = arrow(length = unit(0.3, "cm")), color = arrow_color)
    
    if (!is.null(searched_gene) && toupper(searched_gene) %in% data$Gene) {
      p <- p + geom_point(data = data %>% filter(Gene == toupper(searched_gene)), 
                          aes(x = log2FoldChange, y = -log10(padj)), color = searched_gene_color, size = 3)
    } else if (!is.null(searched_gene) && !(toupper(searched_gene) %in% data$Gene)) {
      showNotification(paste("Gene", searched_gene, "not found."), type = "warning")
    }
    
    return(p)
  }
  
  output$volcano_plot <- renderPlot({
    data <- unfiltered_data()
    req(data)
    
    p <- generate_volcano_plot(data)
    print(p)
  })
  
  output$download_volcano_plot <- downloadHandler(
    filename = function() { paste('volcano.png') },
    content = function(file) {
      data <- unfiltered_data()
      req(data)
      
      p <- generate_volcano_plot(data)
      ggsave(file, plot = p, width = 20, height = 13, units = "in", dpi = 300, bg = "white")
    }
  )
  
  ### Heatmap
  
  output$title_heatmap <- renderUI({
    req(input$tissue_select1, input$tissue_select2)
    h3(paste("Heatmap:", input$tissue_select1, "vs.", input$tissue_select2))
  })
  
  search_genes_heatmap <- reactiveVal(NULL)
  
  observeEvent(input$search_gene_heatmap, {
    data <- unfiltered_data()
    req(data)
    
    searched_genes <- strsplit(input$search_genes_heatmap, ",\\s*")[[1]]
    
    if (length(searched_genes) > 0) {
      gene_data <- data %>% filter(Gene %in% searched_genes)
      not_found_genes <- setdiff(searched_genes, gene_data$Gene)
      
      if (length(not_found_genes) > 0) {
        showNotification(paste("The following genes were not found:", paste(not_found_genes, collapse = ", ")), type = "warning")
      }
      
      if (nrow(gene_data) > 0) {
        search_genes_heatmap(gene_data$Gene)
      }
    }
  })
  
  observeEvent(input$random_gene_heatmap, {
    data <- filtered_data_combined()
    req(data)
    
    random_genes <- sample(data$Gene, 10)
    search_genes_heatmap(random_genes)
    updateTextInput(session, "search_genes_heatmap", value = paste(random_genes, collapse = ", "))
  })
  
  observeEvent(input$clear_search, {
    search_genes_heatmap(NULL)
    updateTextInput(session, "search_genes_heatmap", value = "")
  })
  
  generate_plot <- function(data, search_genes) {
    if (input$sort_by == "FDR") {
      data <- data %>% arrange(padj)
    } else {
      data <- data %>% arrange(desc(abs(log2FoldChange)))
    }
    
    search_genes <- search_genes_heatmap()
    
    if (!is.null(search_genes) && length(search_genes) > 0) {
      top_genes <- search_genes
    } else {
      top_genes <- data %>% 
        head(as.numeric(input$top_n)) %>% 
        pull(Gene)
    }
    
    selected_tissues <- c(input$tissue_select1, input$tissue_select2)
    
    samples <- phenodata %>%
      filter(Tissue %in% selected_tissues) %>%
      pull(Sample_name)
    
    rpkm_filtered <- RPKM_data %>%
      filter(rownames(.) %in% top_genes)
    
    samples <- samples[samples %in% colnames(rpkm_filtered)]
    
    rpkm_filtered <- rpkm_filtered[, samples, drop = FALSE]
    
    colnames(rpkm_filtered) <- gsub("Sample_", "", colnames(rpkm_filtered))
    
    annotation <- phenodata %>%
      filter(Sample_name %in% samples) %>%
      select(Sample_name, Tissue) %>%
      mutate(Sample_name = gsub("Sample_", "", Sample_name)) %>%
      column_to_rownames(var = "Sample_name")
    
    present_tissues <- unique(annotation$Tissue)
    annotation_colors <- list(Tissue = mycolors1[present_tissues])
    
    n_genes <- nrow(rpkm_filtered)
    fontsize_row <- ifelse(n_genes > 50, 7, 14)
    show_rownames <- n_genes <= 100
    
    color_scale <- switch(input$color_scale,
                          "blue_white_red" = c("royalblue", "white", "firebrick3"),
                          "green_black_red" = c("springgreen2", "black", "firebrick2"),
                          "purple_white_green" = c("purple", "white", "springgreen4"),
                          "viridis" = viridis(100),
                          "plasma" = plasma(100),
                          "cividis" = cividis(100),
                          "inferno" = inferno(100))
    
    heatmap <- pheatmap(
      rpkm_filtered, 
      scale = "row",
      cluster_rows = TRUE, 
      cluster_cols = TRUE, 
      show_rownames = show_rownames, 
      show_colnames = TRUE,
      annotation_col = annotation,
      annotation_colors = annotation_colors,
      color = colorRampPalette(color_scale)(100),
      cutree_cols = 2,
      cutree_rows = 2,
      angle_col = 45,
      fontsize = 13,
      fontsize_col = 14,
      fontsize_row = fontsize_row,
      width = 30
    )
    
    ggheatmap <- as.ggplot(heatmap$gtable) + 
      theme(text = element_text(family = "Sen"), 
            legend.position = "right", 
            legend.justification = "center")
    
    ggheatmap
  }
  
  output$heatmap_plot <- renderPlot({
    data <- filtered_data_combined()
    req(data)
    search_genes <- search_genes_heatmap()
    generate_plot(data, search_genes)
  })
  
  output$download_heatmap_plot <- downloadHandler(
    filename = function() {paste('heatmap.png')},
    content = function(file) {
      data <- filtered_data_combined()
      req(data)
      
      png(file, width = 18, height = 13, units = "in", res = 300)
      generate_plot(data, search_genes_heatmap())
      dev.off()
    }
  )
  
  ### Correlation Plot
  
  # Use RPKM data
  
  # Either all tissues, selected tissues, or by system
  
  # Use pcaExplorer
  
  observe({
    selected_samples <- NULL
    
    if (input$correlation_mode == "selected") {
      selected_tissues <- c(input$tissue_select1, input$tissue_select2)
      selected_samples <- phenodata %>%
        filter(Tissue %in% selected_tissues) %>%
        arrange(Tissue) %>%
        pull(Sample_name2)
    } else if (input$correlation_mode == "system") {
      selected_system <- input$correlation_system
      selected_tissues <- phenodata %>%
        filter(System == selected_system) %>%
        pull(Tissue) %>%
        unique()
      
      selected_samples <- phenodata %>%
        filter(Tissue %in% selected_tissues) %>%
        arrange(Tissue) %>%
        pull(Sample_name2)
    }
    
    req(selected_samples)
    
    N <- 1000  # Number of top variable genes to select
    gene_variability <- apply(RPKM_data, 1, var)
    top_genes <- names(sort(gene_variability, decreasing = TRUE)[1:N])
    
    sample_names1 <- phenodata$Sample_name[match(selected_samples, phenodata$Sample_name2)]
    selected_data <- RPKM_data[top_genes, sample_names1, drop = FALSE]
    
    colnames(selected_data) <- selected_samples
    
    sample_tissues <- phenodata$Tissue[match(selected_samples, phenodata$Sample_name2)]
    sample_colors <- mycolors1[sample_tissues]
    
    if (ncol(selected_data) > 1) {
      output$correlation_plot <- renderPlot({
        pairs(selected_data, 
              upper.panel = function(x, y) {
                points(x, y)
              },
              lower.panel = function(x, y) {
                points(x, y)
              })
      }, height = 800, width = 1200)
    } else {
      output$correlation_plot <- renderPlot({
        plot.new()
        text(0.5, 0.5, "Not enough data to compute correlations", cex = 1.5)
      })
    }
  })
}

### PCA

#Should be able to select colors and symbols

### Boxplots

# Search for gene

shinyApp(ui = ui, server = server)
