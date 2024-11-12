
#setwd("~/Documents/Internship/InfoGrab")

required_packages <- c(
  "shiny", "shinyBS", "shinythemes", "shinycustomloader", "shinycssloaders", 
  "ggplot2", "pheatmap", "dplyr", "readr", "readxl", "DT", "stringr", 
  "purrr", "tibble", "reshape2", 
  "RColorBrewer", "gplots", "scales", "data.table", "gridExtra", 
  "lubridate", "shinyjs", "shinydashboard", 
  "shinycssloaders", "shinyWidgets", "htmltools", "tools", 
  "htmlwidgets", "DESeq2", "ggplotify", "grid", "graphics", 
  "igraph", "seriation", "GGally", "knitr", "devtools", "stringr",
  "BiocManager", "VariantAnnotation", "bslib"
)

install_if_missing <- function(pkg) {
  if (!require(pkg, character.only = TRUE)) {
    install.packages(pkg, dependencies = TRUE)
  }
}

lapply(required_packages, install_if_missing)

if (!requireNamespace("BiocManager", quietly = TRUE)) {
  install.packages("BiocManager")
  BiocManager::install("VariantAnnotation", force = TRUE)
  BiocManager::install("preprocessCore")
}

if (!require(DESeq2)) {
  BiocManager::install("DESeq2", force = TRUE)
}

if (!require(pcaExplorer)) {
  BiocManager::install("pcaExplorer", force = TRUE)
}

if (!require(rtracklayer)) {
  install.packages("https://cran.r-project.org/src/contrib/Archive/rjson/rjson_0.2.20.tar.gz", repos = NULL, type = "source")
  BiocManager::install("rtracklayer", force = TRUE)
}

if (!require(tispec)) {
 library(devtools)
 library(remotes)
 remotes::install_github('roonysgalbi/tispec')
}

if (!require(igvShiny)) {
  library(devtools)
  install_github("paul-shannon/igvShiny")
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
library(DESeq2)
library(pcaExplorer)
library(ggplotify)
library(htmlwidgets)
library(tools)
library(grid)
library(graphics)
library(igraph)
library(GGally)
library(gridExtra)
library(tispec)
library(knitr)
library(stringr)
library(JBrowseR)
library(bslib)
library(igvShiny)
library(rtracklayer)
library(VariantAnnotation)
library(Rsamtools)
library(preprocessCore)


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
  `Iliotibialis.lateralis`="#A1122A",
  `Pectoralis.major`="#FFCCCC",
  `Isthmus`="#FF8000",
  `Magnum`="#DE5100",
  `Shell.Gland`="#FFC78E",
  `Ovary`="#CCCC00"
)

mycolors2 <- c(
  `Immune` = "#31E1F7", 
  `Respiratory` = "#400D51", 
  `Excretory` ="#FF7777",
  `Muscular` = "#D800A6",
  `Intestine` ="#6499E9", 
  `Reproductive` = "#836FFF")

############################################################

# INFO GRAB

# Cite BioRender
# Cite pcaExplorer
# Cite IGV (https://igv.org/doc/webapp/##citing-igv)

# Need to update how the data loads

chrom_map <- read_tsv("browser_data_for_app/chrm-ncbi.txt", col_names = c("Chrm_ncbi", "Chromosome"))

# Load the GFF file and fix chromosome names
gff_file <- "browser_data_for_app/genomeAnnoatation_gallus.updated.gff"
gff_data <- import(gff_file, format = "gff")
gff_df <- as.data.frame(gff_data)

# Rename the seqnames column to 'chr' for compatibility
colnames(gff_df)[colnames(gff_df) == "seqnames"] <- "chr"

# Merge chromosome mapping with GFF data
gff_df <- gff_df %>%
  left_join(chrom_map, by = c("chr" = "Chrm_ncbi")) %>%
  mutate(chr = Chromosome) %>%
  dplyr::select(-Chromosome)

tissue_names <- unique(phenodata$Tissue)

systems <- unique(phenodata$System)

tissue_by_system <- lapply(systems, function(system) {
  unique(phenodata$Tissue[phenodata$System == system])
})

names(tissue_by_system) <- systems

##################

ui <- navbarPage(
  title = tags$a(
    tags$div(
      tags$a(
        href = "https://www.washington.edu",
        target = "_blank",
        tags$img(src = "UW_logo.png", height = "50px", style = "padding: 5px; margin-top: -17px;")
      ),
      tags$a(
        href = "https://www.westernu.edu",
        target = "_blank",
        tags$img(src = "WesternU_logo.png", height = "60px", style = "padding: 5px; margin-top: -19px;")
      ),
      tags$a(
        href = "https://www.hawkinslab.org",
        target = "_blank",
        tags$img(src = "hawkins_lab.png", height = "60px", style = "padding: 5px; margin-top: -19px;")
      )
      
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
     #tissue-selection-menu {
        position: fixed;
        bottom: 20px;
        right: 20px;
        width: 250px;
        z-index: 1000;
        background-color: white;
        border: 1px solid #ccc;
        border-radius: 8px;
        box-shadow: 2px 2px 5px rgba(0,0,0,0.3);
      }
      #tissue-selection-menu h4 {
        margin: 10px;
        font-size: 16px;
      }
      #tissue-selection-menu .close-btn {
        position: absolute;
        top: 5px;
        right: 10px;
        cursor: pointer;
        font-size: 18px;
      }
      #tissue-selection-menu .content {
        padding: 10px;
      }
      #tissue-selection-menu input, #tissue-selection-menu select, #tissue-selection-menu button {
        width: 100%;
        margin-bottom: 10px;
      }
    "))
  ),
  tabPanel("INFO GRAB",
           fluidPage(
             div(style = "text-align: center;",
                 h2("INFO GRAB"),
                 h3("INteractive Functional Ontology for Gene Regulation Analysis and Browsing"),
                 tags$hr(style = "height:1px; border:none; color:#300; background-color:#300;"),
                 h5("Welcome to INFO GRAB. This app allows you to perform differential expression analysis and 
                   visualize the results through various plots, looking at tissues in the chicken (Gallus gallus)."),
                 tags$hr(style = "height:1px; border:none; color:#300; background-color:#300;"),
                 img(src = "chicken.png", id = "chicken-image", style = "width: 40%; max-width: 100%; height: auto;")
             ),
             uiOutput("tissue_display")
           )
  ),
  
  tabPanel("Differential Analysis",
           tabsetPanel(
             id = "differential_analysis_tab",
             tabPanel("DE Home",
                      sidebarLayout(
                        sidebarPanel(
                          selectInput("tissue_select1", "Tissue 1", choices = tissue_names),
                          selectInput("tissue_select2", "Tissue 2", choices = tissue_names[-1]),
                          tags$hr(style = "height:1px; border:none; color:#300; background-color:#300;"),
                          actionButton("add_all_de_cart", "Grab All DE Genes"),
                          uiOutput("de_genes_tissue1_button"),
                          uiOutput("de_genes_tissue2_button"),
                          tags$hr(style = "height:1px; border:none; color:#300; background-color:#300;"),
                          selectInput("selected_fdr", "FDR", 
                                      choices = c("No filter" = NA, "0.01" = 0.01, "0.05" = 0.05, "0.1" = 0.1), 
                                      selected = "0.05"),
                          sliderInput("selected_fc", "FC (log2 fold change)", min = 0, max = 20, step = 1, value = 0),
                          tags$hr(style = "height:1px; border:none; color:#300; background-color:#300;"),
                          sliderInput("num_of_digits", "Number of digits", value = 3, min = 2, max = 20),
                          tags$hr(style = "height:1px; border:none; color:#300; background-color:#300;"),
                          downloadButton("download_DEA_data", "Download CSV"),
                          tags$hr(style = "height:1px; border:none; color:#300; background-color:#300;"),
                          tags$p("FDR: The expected proportion of false discoveries among the rejected hypotheses."),
                          tags$p("FC: A measure describing how much a quantity changes between an original and a subsequent measurement."),
                          uiOutput("fold_change_equation")
                        ),
                        mainPanel(
                          uiOutput("title"),
                          
                          withLoader(tableOutput("summary_info"), type = "html", loader = "dnaspin"),
                          br(),
                          textInput("searched_gene", "Search for genes (comma-separated)", width = '600px', value = "", placeholder = "Search for genes here"),
                          actionButton("search_gene", "Search"),
                          actionButton("random_gene", "Search 10 Random Genes"),
                          actionButton("clear_search_main", "Clear"),
                          actionButton("add_to_cart_home", "Grab Genes"),
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
                            condition = "input.correlation_mode == 'selected'",
                            checkboxInput("use_filtered_data", "Use Filtered Data", value = TRUE)
                          ),
                          conditionalPanel(
                            condition = "input.correlation_mode == 'system'",
                            selectInput("correlation_system", "Select System:", choices = systems)
                          ),
                          tags$hr(style = "height:1px; border:none; color:#300; background-color:#300;"),
                          downloadButton("download_correlation_plot", "Download PNG"),
                          tags$hr(style = "height:1px; border:none; color:#300; background-color:#300;"),
                          tags$p("This section allows you to visualize the correlation between the expression 
                                 levels of selected genes across different tissues or systems. Correlation plots are useful 
                                 for understanding how gene expression patterns are related between different conditions or tissue types."),
                          tags$p("The diagonal panels show density plots of the gene expression 
                                 distribution for each sample, while the upper panels display the correlation 
                                 coefficients between pairs of samples.")
                        ),
                        mainPanel(
                          uiOutput("title_cor_plot"),
                          withLoader(plotOutput("correlation_plot", width = "90%", height = "750px"), type = "html", loader = "dnaspin")
                        )
                      )
             ),
             
             
             tabPanel("Volcano Plot",
                      sidebarLayout(
                        sidebarPanel(
                          textInput("searched_gene_volcano", "Search for gene(s): (comma-separated):", width = '600px', value = "", placeholder = "Search for genes here"),
                          actionButton("search_gene_volcano", "Search"),
                          actionButton("clear_search_gene_volcano", "Clear"),

                          tags$hr(style = "height:1px; border:none; color:#300; background-color:#300;"),
                          
                            radioButtons("labeling_metric", "Label genes by:",
                                         choices = list("Smallest FDR" = "fdr", "Largest Absolute Fold Change" = "fc"),
                                         selected = "fdr"),
                          
                          
                          radioButtons("labeling_option", "Labeling Options:",
                                       choices = list("Basic labels" = "basic", "Advanced labels" = "advanced"),
                                       selected = "basic"),
                          
                          conditionalPanel(
                            condition = "input.labeling_option == 'basic'",
                            sliderInput("n_labels", "Number of Labels", min = 0, max = 100, value = 10)
                          ),
                          
                          conditionalPanel(
                            condition = "input.labeling_option == 'advanced'",
                            sliderInput("n_labels_left", "Number of Labels (Left)", min = 0, max = 50, value = 5),
                            sliderInput("n_labels_right", "Number of Labels (Right)", min = 0, max = 50, value = 5)
                          ),

                          actionButton("add_labeled_genes", "Grab Labeled Genes"),
                          
                          tags$hr(style = "height:1px; border:none; color:#300; background-color:#300;"),
                          downloadButton("download_volcano_plot", "Download PNG"),
                          tags$hr(style = "height:1px; border:none; color:#300; background-color:#300;"),
                          tags$p("Each point represents a gene, plotted by its log2 fold change and the negative 
                          logarithm of its adjusted p-value. Genes with significant changes (high fold change 
                                 and low p-value) are highlighted."),
                          tags$p("Genes with the largest log2 fold changes (farthest from the center on the x-axis) 
                                 show the most substantial differences in expression between the two tissues. 
                                 The y-axis represents the significance level, with higher points indicating more 
                                 statistically significant changes.")
                        ),
                        mainPanel(
                          uiOutput("title_volcano_plot"),
                          withLoader(plotOutput("volcano_plot", width = "90%", height = "750px"), 
                                     type = "html", loader = "dnaspin")
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
                          textInput("search_genes_heatmap", "Search for Specific Genes (comma-separated):", value = "", placeholder = "Search for genes here"),
                          actionButton("search_gene_heatmap", "Search"),
                          actionButton("random_gene_heatmap", "Search 10 Random Genes"),
                          actionButton("clear_search", "Clear"),
                          tags$hr(style = "height:1px; border:none; color:#300; background-color:#300;"),
                          actionButton("add_to_cart_heatmap", "Grab Displayed Genes"),
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
                          tags$p("Each cell represents the expression level of a gene in a sample, with colors indicating the level of expression."),
                          tags$p("The color scale is based on z-scores, which represent the number of standard deviations away from the mean expression level. High positive z-scores indicate higher expression levels, while high negative z-scores indicate lower expression levels.")
                        ),
                        mainPanel(
                          uiOutput("title_heatmap"),
                          withLoader(plotOutput("heatmap_plot", width = "85%", height = "750px",
                                                brush = brushOpts(id = "heatmap_brush", resetOnNew = TRUE)),
                                     type = "html", loader = "dnaspin"),
                          tags$hr(),
                          tags$p(tags$b("Brushing: Drag to select genes for further analysis."))
                        )
                      )
             ),
             
             tabPanel("PCA",
                      sidebarLayout(
                        sidebarPanel(
                          radioButtons("pca_mode", "Select Mode:",
                                       choices = c("Selected Tissues" = "selected", "Tissues by System" = "system", "All Tissues" = "all"),
                                       selected = "selected"),
                          conditionalPanel(
                            condition = "input.pca_mode == 'system'",
                            selectInput("pca_system_select", "Select System:", choices = unique(phenodata$System))
                          ),
                          tags$hr(style = "height:1px; border:none; color:#300; background-color:#300;"),
                          downloadButton("download_pca_plot", "Download PNG"),
                          tags$hr(style = "height:1px; border:none; color:#300; background-color:#300;"),
                          tags$p("This section allows you to perform Principal Component Analysis (PCA) on your data. 
                                 PCA is a dimensionality reduction technique that projects high-dimensional data into a lower-dimensional space, 
                                 making it easier to visualize patterns and relationships in the data."),
                          tags$p("The first two principal components (PC1 and PC2) are plotted, showing the most variance in the data. 
                                 Each point represents a sample, colored by tissue type and shaped by system.")
                        ),
                        mainPanel(
                          uiOutput("title_pca_plot"),
                          withLoader(plotOutput("pca_plot", height = "600px"), type = "html", loader = "dnaspin")
                        )
                      )
             ),
             
             conditionalPanel(
               condition = "input.differential_analysis_tab != 'DE Home'",
               div(
                 id = "tissue-selection-menu",
                 wellPanel(
                   actionButton("toggle_tissue_menu", "Select Tissues"),
                   hidden(
                     div(
                       id = "tissue_menu_content",
                       selectInput("tissue_select1_menu", "Tissue 1", choices = tissue_names),
                       selectInput("tissue_select2_menu", "Tissue 2", choices = tissue_names),
                       actionButton("apply_tissue_selection_menu", "Apply")
                     )
                   )
                 )
               )
             )
           )
  ),
  
  tabPanel("Gene Expression Visualization",
           sidebarLayout(
             sidebarPanel(
               textInput("barplot_searched_gene", "Search for genes (comma-separated)", width = '600px', value = "", placeholder = "Search for genes here"),
               actionButton("search_gene_barplot", "Search"),
               actionButton("random_gene_barplot", "Search 5 Random Genes"),
               actionButton("clear_search_gene_barplot", "Clear"),
               actionButton("add_to_cart_barplot", "Grab Genes"),
               actionButton("transfer_to_heatmap", "Show Heatmap of Searched Genes"),
               tags$hr(style = "height:1px; border:none; color:#300; background-color:#300;"),
               downloadButton("download_barplot", "Download PNG"),
               tags$hr(style = "height:1px; border:none; color:#300; background-color:#300;"),
               tags$p("This section provides a barplot of the average expression levels of selected genes across 
                      different tissues. This visualization is helpful for comparing how genes are expressed in different tissue types."),
               tags$p("You can search for specific genes and display their average expression 
                      levels across the selected tissues. The plot allows for a quick comparison of gene expression patterns and identification of tissue-specific expression.")
             ),
             mainPanel(
               uiOutput("dynamic_title"),
               uiOutput("dynamic_barplot_output")
             )
           )
  ),
  
  tabPanel("Tissue-Specific Analysis",
           tabsetPanel(
             tabPanel("Heatmap",
                      sidebarLayout(
                        sidebarPanel(
                          selectInput("n_variable_genes", "Number of Genes (Ordered by Variance)", 
                                      choices = c(50, 100, 500, 1000), selected = 500),
                          tags$hr(style = "height:1px; border:none; color:#300; background-color:#300;"),
                          checkboxGroupInput("system_choice", "Select System(s):", 
                                             choices = c("All Systems", unique(phenodata$System)),
                                             selected = "All Systems"),
                          uiOutput("tissue_selection_ui"),
                          tags$hr(style = "height:1px; border:none; color:#300; background-color:#300;"),
                          selectInput("color_scale_tissue_specific", "Select Color Scale", 
                                      choices = c("Purple-White-Green" = "purple_white_green",
                                                  "Blue-White-Red" = "blue_white_red",
                                                  "Green-Black-Red" = "green_black_red",
                                                  "Blue-Yellow-Purple" = "cyan_yellow_purple",
                                                  "Viridis" = "viridis",
                                                  "Plasma" = "plasma",
                                                  "Cividis" = "cividis",
                                                  "Inferno" = "inferno"),
                                      selected = "purple_white_green"),
                          tags$hr(style = "height:1px; border:none; color:#300; background-color:#300;"),
                          actionButton("add_to_cart_tissue_specific_heatmap", "Grab Displayed Genes"),
                          tags$hr(style = "height:1px; border:none; color:#300; background-color:#300;"),
                          downloadButton("download_tissue_specific_heatmap", "Download PNG"),
                          tags$hr(style = "height:1px; border:none; color:#300; background-color:#300;"),
                          tags$p("This section displays a tissue-specific heatmap, which focuses on the expression of genes that are most variable across the selected tissues."),
                          tags$p("The heatmap is color-coded to represent different levels of gene expression, allowing you to easily identify patterns and differences in gene expression between tissues. This is particularly useful for identifying tissue-specific genes and understanding how gene expression varies across different biological systems.")
                        ),
                        mainPanel(
                          uiOutput("title_tissue_specific"),
                          withLoader(plotOutput("tissue_specific_heatmap", width = "100%", height = "750px",
                                                click = "heatmap_click", 
                                                brush = brushOpts(id = "heatmap_tissue_specific_brush", resetOnNew = TRUE)), 
                                     type = "html", loader = "dnaspin"),
                          tags$hr(),
                          tags$p(tags$b("Brushing: Drag to select genes for further analysis."))
                        )
                      )
             )
           )
  ),
  
  
  # Genome Browser
  
  tabPanel("Genome Viewer",
           fluidPage(
             tags$head(
               tags$style(HTML("
        .shiny-input-container {
          margin-bottom: 0px; /* Remove bottom margin */
        }
        .btn-primary {
          margin-top: 24px; /* Align button with input fields */
          width: 40%;
        }
        #loadTrack, #load_geneASE_track, #load_snp, #load_gff {
          width: 40%;
        }
        .action-button {
          margin-top: 5px;
        }
      "))
             ),
             fluidRow(
               column(4,
                      selectInput("genome_select", "Select a Genome",
                                  choices = c("galGal6"), selected = "galGal6"),
                      br(),
                      textInput("genome_browser_search", "Search For a Locus or a Gene"),
                      actionButton("search_button", "Search", class = "btn-primary")
               ),
               column(2,
                      radioButtons("input_type", "Add a Track (GFF3, BAM, BED, or VCF):",
                                   choices = c("Local File" = "file", "URL" = "url"),
                                   selected = "file")
               ),
               column(4,
                      conditionalPanel(
                        condition = "input.input_type == 'file'",
                        fileInput("file", "Track File", width = '100%')
                      ),
                      conditionalPanel(
                        condition = "input.input_type == 'url'",
                        textInput("url", "Track URL"),
                        textInput("index", "Index URL")
                      ),
                      actionButton("loadTrack", "Load Custom Track", class = "btn-primary"),
                      actionButton("load_gff", "Load Annotation Track"),
                      actionButton("load_geneASE_track", "Load Gene ASE Track", class = "action-button"),
                      actionButton("load_snp", "Load SNP ASE Track", class = "action-button"),
                      actionButton("grab_all_ase_genes", "Grab All ASE Genes", class = "btn-primary") # <-- Added here
               )
             ),
             br(),
             fluidRow(
               column(12,
                      igvShinyOutput("igvShiny", height = "1000px")
               )
             )
           )
  ),
  
  
  
  tabPanel("Genes",
           sidebarLayout(
             sidebarPanel(
               tags$textarea(id = "gene_text_area", style = "position: absolute; left: -9999px;"),
               actionButton("copy_cart_genes", "Copy Genes to Clipboard"),
               actionButton("clear_cart", "Clear Genes"),
               tags$hr(style = "height:1px; border:none; color:#300; background-color:#300;"),
               textInput("add_genes_input", "Add Genes (comma separated):", value = "", placeholder = "Search for genes here"),
               actionButton("add_genes_button", "Add Genes"),
               tags$hr(style = "height:1px; border:none; color:#300; background-color:#300;"),
               sliderInput("num_of_digits_tau", "Number of digits", value = 3, min = 2, max = 20),
               tags$hr(style = "height:1px; border:none; color:#300; background-color:#300;"),
               
               # Summary Statistics
               tags$h3("Summary Statistics"),
               tableOutput("summary_stats_table"),
               
               tags$h3("Tau Value Distribution"),
               plotOutput("tau_histogram", height = "350px"),
               
               tags$h3("Allele-Specific Genes Count"),
               plotOutput("allele_specific_bar_chart", height = "350px"),
               
               tags$hr(style = "height:1px; border:none; color:#300; background-color:#300;"),
               
               tags$h3("Most Frequent Tissues"),
               plotOutput("top_max_tissues_table", height = "350px"),
               
               
               tags$hr(style = "height:1px; border:none; color:#300; background-color:#300;"),
               
               tags$h3("Useful Links"),
               div(style = "font-size: 15px;",
                   tags$ul(
                     tags$li(a(href = "http://bioinformatics.sdstate.edu/go/", "ShinyGO: Gene Ontology Enrichment Analysis", target = "_blank")),
                     tags$li(a(href = "https://davidbioinformatics.nih.gov/tools.jsp", "DAVID: Database for Annotation, Visualization and Integrated Discovery", target = "_blank")),
                     tags$li(a(href = "https://pantherdb.org", "PantherDB: Panther Classification System", target = "_blank"))
                   )
               ),
               tags$hr(style = "height:1px; border:none; color:#300; background-color:#300;"),
               downloadButton("download_results", "Download Genes (CSV)")
             ),
             mainPanel(
               withLoader(DTOutput("cart_gene_table"), loader = "dnaspin")
             )
           )
  ),
  
  tabPanel("Help",
           fluidPage(
             h2("App Version and Information"),
             
             tabsetPanel(
               
               # Differential Analysis Help Tab
               tabPanel("Differential Analysis Help",
                        mainPanel(
                          h2("Analysis Help"),
                          h3("DE Home"),
                          p("In the 'DE Home' tab, you can select two tissues to compare and adjust the analysis parameters."),
                          tags$ul(
                            tags$li(tags$b("Tissue 1 and Tissue 2:"), " Select the tissues you want to compare."),
                            tags$li(tags$b("FDR:"), " Choose the false discovery rate threshold for filtering the results."),
                            tags$li(tags$b("FC (log2 fold change):"), " Set the fold change threshold for filtering the results."),
                            tags$li(tags$b("Number of digits:"), " Specify the number of decimal places for displaying values."),
                            tags$li(tags$b("Download CSV:"), " Download the filtered data as a CSV file."),
                            tags$li(tags$b("Search for genes:"), " Enter genes to search in the results."),
                            tags$li(tags$b("Search 10 Random Genes:"), " Display information for 10 randomly selected genes."),
                            tags$li(tags$b("Clear:"), " Clear the gene search input."),
                            tags$li(tags$b("Fold Change Equation:"), " Displays the log2 fold change equation for the selected tissues.")
                          ),
                          h3("Correlation Plot"),
                          p("The 'Correlation Plot' tab allows you to visualize the correlation between expression levels of genes across different tissues."),
                          tags$ul(
                            tags$li(tags$b("Select Mode:"), " Choose between 'Selected Tissues' and 'Tissues by System'."),
                            tags$li(tags$b("Use Filtered Data:"), " Option to use filtered data based on selected FDR and FC thresholds."),
                            tags$li(tags$b("Select System:"), " Select a specific system to view correlation within that system."),
                            tags$li(tags$b("Correlation Plot Legend:"), " The plot includes a legend that matches the tissue colors with their corresponding tissues.")
                          ),
                          h3("Volcano Plot"),
                          p("The 'Volcano Plot' tab allows you to visualize the results of the differential expression analysis."),
                          tags$ul(
                            tags$li(tags$b("Search for gene:"), " Enter a gene name to highlight it on the plot."),
                            tags$li(tags$b("Labeling Options:"), " Choose between basic or advanced labeling options for the plot."),
                            tags$li(tags$b("Background:"), " Choose the background color for the plot."),
                            tags$li(tags$b("Download PNG:"), " Download the volcano plot as a PNG file."),
                            tags$li(tags$b("Advanced Labeling Options:"), " Allows the user to specify the number of labels on each side of the volcano plot separately.")
                          ),
                          h3("Heatmap"),
                          p("The 'Heatmap' tab provides a heatmap visualization of the expression levels of genes."),
                          tags$ul(
                            tags$li(tags$b("Sort by:"), " Choose to sort the genes by FDR or fold change."),
                            tags$li(tags$b("Number of Top Genes:"), " Select the number of top genes to include in the heatmap."),
                            tags$li(tags$b("Search for Specific Genes:"), " Enter gene names (comma-separated) to highlight on the heatmap."),
                            tags$li(tags$b("Random Gene Search:"), " Search for 10 random genes and highlight them on the heatmap."),
                            tags$li(tags$b("Select Color Scale:"), " Choose the color scale for the heatmap."),
                            tags$li(tags$b("Download PNG:"), " Download the heatmap as a PNG file.")
                          ),
                          h3("PCA"),
                          p("The 'PCA' tab allows you to perform Principal Component Analysis on the selected tissues."),
                          tags$ul(
                            tags$li(tags$b("Select Mode:"), " Choose between 'Selected Tissues', 'Tissues by System', and 'All Tissues'."),
                            tags$li(tags$b("Select System:"), " Select a specific system for PCA analysis when 'Tissues by System' mode is selected."),
                            tags$li(tags$b("Download PCA Plot:"), " Download the PCA plot as a PNG file.")
                          ),
                          h3("Parameters Explanation"),
                          tags$p("False Discovery Rate (FDR): The expected proportion of false discoveries among the rejected hypotheses."),
                          tags$p("Fold Change (FC): A measure describing how much a quantity changes between an original and a subsequent measurement.")
                        )
               ),
               
               # Gene Expression Visualization Help Tab
               tabPanel("Gene Expression Visualization Help",
                        mainPanel(
                          h2("Gene Expression Visualization Help"),
                          h3("Overview"),
                          p("The 'Gene Expression' tab allows users to visualize the average expression levels of selected genes across different tissues using bar plots."),
                          h3("Selecting Tissues and Genes"),
                          tags$ul(
                            tags$li(tags$b("Tissue 1 and Tissue 2:"), " Select two tissues for comparison."),
                            tags$li(tags$b("Top Genes:"), " The top 5 commonly expressed genes are automatically selected and displayed."),
                            tags$li(tags$b("Search for Genes:"), " You can search for specific genes by entering their names into the search box."),
                            tags$li(tags$b("Search 5 Random Genes:"), " Display 5 random genes from the dataset."),
                            tags$li(tags$b("Grab Genes:"), " Grab the genes displayed in the bar plot to the list of grabbed gene."),
                            tags$li(tags$b("Transfer to Heatmap:"), " Transfer searched genes to the heatmap for further exploration.")
                          ),
                          h3("Bar Plot"),
                          p("The bar plot shows the average expression levels of selected genes across tissues, allowing for easy comparison."),
                          h3("Transferring to Heatmap"),
                          p("Selected genes can be transferred to the heatmap for detailed visualization of expression patterns."),
                          tags$ul(
                            tags$li(tags$b("Heatmap Modal:"), " The heatmap is displayed in a modal dialog with clustering and color-coded annotations."),
                            tags$li(tags$b("Download Heatmap:"), " Download the heatmap as a PNG file.")
                          )
                        )
               ),
               
               # Tissue-Specific Analysis Help Tab
               tabPanel("Tissue-Specific Analysis Help",
                        mainPanel(
                          h2("Tissue-Specific Analysis Help"),
                          h3("Overview"),
                          p("The 'Tissue-Specific Analysis' tab visualizes tissue-specific genes across tissues and systems using heatmaps."),
                          h3("Selecting Systems and Tissues"),
                          tags$ul(
                            tags$li(tags$b("System Choice:"), " Select one or more biological systems to filter the tissues."),
                            tags$li(tags$b("All Systems:"), " Select 'All Systems' to include all available tissues."),
                            tags$li(tags$b("Tissue Choice:"), " Choose tissues to visualize."),
                            tags$li(tags$b("Number of Variable Genes:"), " Specify the number of most variable genes to display.")
                          ),
                          h3("Heatmap Visualization"),
                          p("Visualizes the most variable tissue-specific genes across selected tissues."),
                          tags$ul(
                            tags$li(tags$b("Tau Scores:"), " Only genes with Tau scores of ", tags$b("0.85"), " or higher are included."),
                            tags$li(tags$b("Clustering:"), " Genes are clustered by correlation and tissues are ordered by system.")
                          ),
                          h3("Interactive Heatmap Features"),
                          p("Explore tissue-specific gene expression interactively by brushing genes."),
                          h3("Grabbed Genes"),
                          p("Selected genes can be added to the genes of interest list for further analysis.")
                        )
               ),
               
               # Genome Browser Help Tab
               tabPanel("Genome Browser Help",
                        mainPanel(
                          h2("Genome Browser Help"),
                          h3("Overview"),
                          p("The Genome Browser allows you to visualize genomic data by uploading or linking to file formats like GFF3, BAM, BED, and VCF."),
                          h3("Supported File Types"),
                          tags$ul(
                            tags$li(tags$b("GFF3/GFF:"), " Upload for gene annotations."),
                            tags$li(tags$b("BAM:"), " Upload for aligned sequencing reads."),
                            tags$li(tags$b("BED:"), " Upload for defining genomic regions."),
                            tags$li(tags$b("VCF:"), " Upload for storing genetic variations.")
                          ),
                          h3("File Upload"),
                          p("Upload files from your local machine. Supported formats include GFF3, BAM, BED, and VCF."),
                          h3("Loading Tracks from URLs"),
                          p("Load tracks directly from URLs."),
                          h3("Navigating the Genome"),
                          p("Enter a genomic region in the search field to navigate to specific regions."),
                          h3("Common Errors"),
                          tags$ul(
                            tags$li(tags$b("Unsupported File Type:"), " Ensure that uploaded files are in the supported formats."),
                            tags$li(tags$b("Chromosome Mapping Errors:"), " If the chromosome mapping file is missing or incorrectly formatted, an error will be shown.")
                          )
                        )
               ),
               
               # Gene Cart Help Tab
               tabPanel("Genes Help",
                        mainPanel(
                          h2("Grabbed Gene and Allele-Specific Expression (ASE) Help"),
                          h3("Genes List Functionality"),
                          tags$ul(
                            tags$li(tags$b("Grab Genes:"), " Add genes from various sections of the app for further analysis."),
                            tags$li(tags$b("Clearing the List:"), " Clear all genes from the cart."),
                            tags$li(tags$b("Copying Genes from List:"), " Copy the gene list for use in other applications."),
                            tags$li(tags$b("Downloading the List:"), "Download a CSV file of your list of genes, associated Tau scores, and maximum expressing tissue."),
                          ),
                          h3("Tau Score Calculation"),
                          p("The Tau score is a measure of tissue specificity, ranging from 0 to 1, with 1 being more tissue specific.", href = "https://www.washington.edu"),
                          h3("Displaying and Interacting with Gene Data"),
                          p("Gene expression data is aggregated by tissue, with Tau scores calculated for each gene."),
                        )
               ),
               
               # Version Tab
               tabPanel("Version",
                        br(),
                        p("Version: 1.0.0"),
                        p("Developed by: Oliver Brown (ombrown@uw.edu)"),
                        tags$a(href = "https://www.hawkinslab.org", "Hawkins Lab", target = "_blank"),
               ),
               
               # Credits Tab
               tabPanel("Credits",
                        h4("Acknowledgments:"),
                        br(),
                        h4("Citations:"),
                        br(),
                        h6("IGV:"),
                        tags$ul(
                          tags$li("Robinson JT, Thorvaldsdóttir H, Winckler W, et al. Integrative Genomics Viewer. *Nature Biotechnology* 29, 24–26 (2011)."),
                          tags$li("Thorvaldsdóttir H, Robinson JT, Mesirov JP. Integrative Genomics Viewer (IGV): high-performance genomics data visualization. *Briefings in Bioinformatics* 14, 178–192 (2013)."),
                          tags$li("Robinson JT, Thorvaldsdóttir H, Wenger AM, et al. Variant Review with the Integrative Genomics Viewer (IGV). *Cancer Research* 77(21):31-34 (2017)."),
                          tags$li("Robinson JT, Thorvaldsdóttir H, Turner D, Mesirov JP. igv.js: an embeddable JavaScript implementation of the Integrative Genomics Viewer (IGV). *Bioinformatics* 39(1), btac830 (2023).")
                        ),
                        br(),
                        h6("BioRender:"),
                        br(),
                        h6("pcaExplorer:")
               )
             )
           )
  )
  
)

##########

server <- function(input, output, session) {
  
  cart_genes <- reactiveVal(character(0))
  
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
                      img(src = img_path, style = "width: 80px; height: auto; margin-right: 13px;")
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
          
          tissue_data <- RPKM_data %>% 
            dplyr::select(samples)
          
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
                paste(top_genes$Gene, sep = "", collapse = "<br>")
              )
            ),
            easyClose = TRUE,
            footer = NULL
          ))
        }, ignoreInit = TRUE)
      })
    })
  })
  
  # Create data
  
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
      return(h3(paste(input$tissue_select1, " vs. ", input$tissue_select2, sep = "")))
    }
    
    if (order == 1) {
      return(h3(paste(input$tissue_select2, " vs. ", input$tissue_select1, sep = "")))
    }
    
  })
  
  # Expressed genes counts
  expressed_genes_by_tissue <- reactive({
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
  
  output$fold_change_equation <- renderUI({
    req(input$tissue_select1, input$tissue_select2)
    order <- order_val()
    
    if (order == 0) {
      equation <- paste0("<b>log<sub>2</sub> fold change = log<sub>2</sub>(", input$tissue_select1, "/", input$tissue_select2, ")</b>")
    } else {
      equation <- paste0("<b>log<sub>2</sub> fold change = log<sub>2</sub>(", input$tissue_select2, "/", input$tissue_select1, ")</b>")
    }
    
    HTML(equation)
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
        mutate(FalseDiscoveryRate = format(padj, scientific = T, digits = input$num_of_digits), 
               log2FoldChange = format(log2FoldChange, scientific = T, digits = input$num_of_digits)) %>%
        dplyr::select(Gene, log2FoldChange, FalseDiscoveryRate)
      
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
  
  output$table_smallest_pvalues <- renderDT({
    data <- filtered_data_combined()
    req(data)
    
    if (nrow(data) == 0) {
      return(datatable(data.frame(Statistic = "No results found for the selected criteria.", Value = "")))
    }
    
    smallest_pvalues <- data %>% 
      arrange(padj) %>%
      mutate(FalseDiscoveryRate = format(padj, scientific = T, digits = input$num_of_digits), 
             log2FoldChange = format(log2FoldChange, scientific = T, digits = input$num_of_digits)) %>%
      head(10) %>%
      dplyr::select(Gene, log2FoldChange, FalseDiscoveryRate)
    
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
      mutate(FalseDiscoveryRate = format(padj, scientific = T, digits = input$num_of_digits), 
             log2FoldChange = format(log2FoldChange, scientific = T, digits = input$num_of_digits)) %>%
      head(10) %>%
      dplyr::select(Gene, log2FoldChange, FalseDiscoveryRate)
    
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
      mutate(FalseDiscoveryRate = format(padj, scientific = T, digits = input$num_of_digits), 
             log2FoldChange = format(log2FoldChange, scientific = T, digits = input$num_of_digits)) %>%
      head(10) %>%
      dplyr::select(Gene, log2FoldChange, FalseDiscoveryRate)
    
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
      mutate(FalseDiscoveryRate = format(padj, scientific = T, digits = input$num_of_digits), 
             log2FoldChange = format(log2FoldChange, scientific = T, digits = input$num_of_digits)) %>%
      head(10) %>%
      dplyr::select(Gene, log2FoldChange, FalseDiscoveryRate)
    
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
      dplyr::select(Gene, baseMean, log2FoldChange, padj)
    
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
      mutate(FalseDiscoveryRate = format(padj, scientific = T, digits = input$num_of_digits), 
             log2FoldChange = format(log2FoldChange, scientific = T, digits = input$num_of_digits)) %>%
      head(10) %>%
      dplyr::select(Gene, log2FoldChange, FalseDiscoveryRate)
    
    datatable(smallest_pvalues, options = list(dom = 't', paging = FALSE, ordering = FALSE), rownames = FALSE)
  })
  
  output$download_DEA_data <- downloadHandler(
    
    filename = function(){"DEA_data.csv"}, 
    content = function(fname){
      data <- filtered_data_combined()
      data <- data %>%
        dplyr::select(-Gene)
      write.csv(data, fname)
    }
  )
  
  observeEvent(input$add_all_de_cart, {
    data <- filtered_data_combined()
    req(data)
    
    all_de_genes <- data$Gene
    add_to_cart(all_de_genes)
    
    showNotification(paste("Genes grabbed."), type = "message")
  })
  
  
  observeEvent(input$add_to_cart_home, {
    genes <- gene_search_result()$Gene
    add_to_cart(genes)
    
    showNotification("Genes grabbed.", type = "message")
  })
  
  observeEvent(input$add_smallest_pvalues_cart, {
    data <- filtered_data_combined()
    req(data)
    
    smallest_pvalues_genes <- data %>% 
      arrange(padj) %>% 
      head(10) %>%
      pull(Gene)
    
    add_to_cart(smallest_pvalues_genes)
    
    showNotification("Genes grabbed.", type = "message")
  })
  
  observeEvent(input$add_largest_fc_cart, {
    data <- filtered_data_combined()
    req(data)
    
    largest_fc_genes <- data %>% 
      arrange(desc(abs(log2FoldChange))) %>% 
      head(10) %>%
      pull(Gene)
    
    add_to_cart(largest_fc_genes)
    
    showNotification(paste("Genes grabbed."), type = "message")
  })
  
  observeEvent(input$add_upregulated_cart, {
    data <- filtered_data_combined()
    req(data)
    
    upregulated_genes <- data %>% 
      arrange(desc(log2FoldChange)) %>% 
      head(10) %>%
      pull(Gene)
    
    add_to_cart(upregulated_genes)
    
    showNotification(paste("Genes grabbed."), type = "message")
  })
  
  observeEvent(input$add_downregulated_cart, {
    data <- filtered_data_combined()
    req(data)
    
    downregulated_genes <- data %>% 
      arrange(log2FoldChange) %>% 
      head(10) %>%
      pull(Gene)
    
    add_to_cart(downregulated_genes)
    
    showNotification(paste("Genes grabbed."), type = "message")
  })
  
  output$collapsePanels <- renderUI({
    bsCollapse(id = "collapseTables", open = "",
               shinyBS::bsCollapsePanel("Ten Smallest p-values", 
                                        DTOutput("table_smallest_pvalues"),
                                        br(),
                                        actionButton("add_smallest_pvalues_cart", "Grab Genes")
               ),
               shinyBS::bsCollapsePanel("Ten Largest Absolute Fold Changes", 
                                        DTOutput("table_largest_fc"),
                                        br(),
                                        actionButton("add_largest_fc_cart", "Grab Genes")
               ),
               shinyBS::bsCollapsePanel(title = tissue_titles()$upregulated, 
                                        DTOutput("table_upregulated"),
                                        br(),
                                        actionButton("add_upregulated_cart", "Grab Genes")
               ),
               shinyBS::bsCollapsePanel(title = tissue_titles()$downregulated, 
                                        DTOutput("table_downregulated"),
                                        br(),
                                        actionButton("add_downregulated_cart", "Grab Genes")
               )
    )
  })
  
  # Popup menu

  observe({
    shinyjs::useShinyjs()  # Enable shinyjs
  })
  
  tissue_menu_visible <- reactiveVal(FALSE) 
  
  observeEvent(input$toggle_tissue_menu, {
    if (tissue_menu_visible()) {
      shinyjs::runjs('$("#tissue-selection-menu").css("height", "75px");')
      tissue_menu_visible(FALSE) 
    } else {

      shinyjs::runjs('$("#tissue-selection-menu").css("height", "400px");')  
      tissue_menu_visible(TRUE)  
    }
    toggle("tissue_menu_content") 
  })
  
  shinyjs::runjs('$("#tissue-selection-menu").css("height", "100px");')
  
  observe({
    updateSelectInput(session, "tissue_select2_menu",
                      choices = tissue_names[tissue_names != input$tissue_select1_menu],
                      selected = isolate(input$tissue_select2_menu))
  })
  
  observe({
    updateSelectInput(session, "tissue_select1_menu",
                      choices = tissue_names[tissue_names != input$tissue_select2_menu],
                      selected = isolate(input$tissue_select1_menu))
  })
  
  # Handle the "Apply" button click
  # Handle the "Apply" button click
  observeEvent(input$apply_tissue_selection_menu, {
    req(input$tissue_select1_menu, input$tissue_select2_menu)
    
    # Update the main tissue selection inputs
    updateSelectInput(session, "tissue_select1", selected = input$tissue_select1_menu)
    updateSelectInput(session, "tissue_select2", selected = input$tissue_select2_menu)
    
    # Close the tissue selection menu after applying
    shinyjs::runjs('$("#tissue-selection-menu").css("height", "75px");')  # Set shorter height when closed
    hide("tissue_menu_content")  # Hide the menu content
    tissue_menu_visible(FALSE)  # Update state to closed
    
    # Optionally, you can trigger any reactivity needed, e.g., re-run analyses or update plots
    # For example, you might invalidate reactive expressions or force re-rendering of plots
  })
  
  
  # Ensure that when tissues are changed on the DE Home tab, the menu reflects those changes
  observe({
    updateSelectInput(session, "tissue_select1_menu", selected = input$tissue_select1)
    updateSelectInput(session, "tissue_select2_menu", selected = input$tissue_select2)
  })
  
  
  ### Volcano Plot
  
  brushed_genes <- reactiveVal(NULL)
  
  
  output$title_volcano_plot <- renderUI({
    req(input$tissue_select1, input$tissue_select2)
    h3(paste("Volcano Plot:", input$tissue_select1, "vs.", input$tissue_select2))
  })
  
  search_gene_volcano <- reactiveVal(NULL)
  
  observeEvent(input$search_gene_volcano, {
    search_genes <- strsplit(toupper(input$searched_gene_volcano), "[, ]+")[[1]]  # Split by comma or space
    search_gene_volcano(search_genes)
  })
  
  observeEvent(input$clear_search_gene_volcano, {
    search_gene_volcano(NULL)
    updateTextInput(session, "searched_gene_volcano", value = "")  # Clear the search bar
  })
  
  
  generate_volcano_plot <- function(data) {
    selected_fdr <- as.numeric(input$selected_fdr)
    selected_fc <- as.numeric(input$selected_fc)
    order <- order_val()
    searched_genes <- search_gene_volcano()  # Use search for multiple genes
    
    # Dynamically determine contrasting text color for searched genes based on tissue colors
    searched_gene_color <- "black"  # Default color for highlighting searched genes
    
    # Add a column to indicate gene labeling logic
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
    
    # Define the plot with adjusted margins
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
      ggplot2::theme_minimal() +
      ggplot2::theme(
        text = ggplot2::element_text(family = "Sen"),
        legend.position = "bottom",
        legend.text = ggplot2::element_text(size = 14),
        legend.title = ggplot2::element_text(size = 16),
        legend.background = ggplot2::element_rect(fill = "#f0f0f0", color = "black"),
        legend.key = ggplot2::element_rect(fill = "white", color = "black"),
        legend.key.size = ggplot2::unit(1.5, "lines"),
        axis.title = ggplot2::element_text(size = 16),
        axis.text = ggplot2::element_text(size = 14),
        plot.margin = ggplot2::margin(40, 40, 40, 40),
        panel.grid.major = ggplot2::element_blank(),
        panel.grid.minor = ggplot2::element_blank(),
        panel.background = ggplot2::element_rect(fill = "white", color = NA),
        axis.line = ggplot2::element_line(color = "black")
      ) +
      ggplot2::guides(color = ggplot2::guide_legend(override.aes = list(shape = 15, size = 5))) +
      geom_hline(yintercept = -log10(selected_fdr), linetype = "dashed", color = "red", alpha = 1/3) +
      geom_vline(xintercept = c(-selected_fc, selected_fc), linetype = "dashed", color = "blue", alpha = 1/3) + 
      annotate("text", x = 0 + max(data$log2FoldChange) * 0.5, 
               y = max(-log10(data$padj), na.rm = TRUE) * 1.2, 
               label = label_1, size = 5, hjust = 0.5) +
      annotate("text", x = 0 - max(data$log2FoldChange) * 0.5, 
               y = max(-log10(data$padj), na.rm = TRUE) * 1.2, 
               label = label_2, size = 5, hjust = 0.5) +
      annotate("segment", x = 0 + max(data$log2FoldChange) * 0.4, 
               y = max(-log10(data$padj), na.rm = TRUE) * 1.175, 
               xend = 0 + max(data$log2FoldChange) * 0.6, 
               yend = max(-log10(data$padj), na.rm = TRUE) * 1.175,
               arrow = arrow(length = unit(0.3, "cm"))) +
      annotate("segment", x = 0 - max(data$log2FoldChange) * 0.4, 
               y = max(-log10(data$padj), na.rm = TRUE) * 1.175, 
               xend = 0 - max(data$log2FoldChange) * 0.6, 
               yend = max(-log10(data$padj), na.rm = TRUE) * 1.175,
               arrow = arrow(length = unit(0.3, "cm")))
    
    # If searched genes exist, only label searched genes
    if (!is.null(searched_genes) && length(searched_genes) > 0) {
      p <- p + geom_point(data = data %>% filter(Gene %in% searched_genes),
                          aes(x = log2FoldChange, y = -log10(padj)),
                          color = searched_gene_color, size = 3) +
        geom_text_repel(data = data %>% filter(Gene %in% searched_genes),
                        aes(label = Gene),
                        max.overlaps = Inf,  # Ensure all searched genes get labeled
                        box.padding = 0.5,   # Add space between labels and points
                        point.padding = 0.5, # Add space between labels and the points
                        color = searched_gene_color, 
                        size = 5)            # Label size and color for searched genes
    } else {
      # If no searched genes, label top genes based on FDR or FC
      if (input$labeling_option == "basic") {
        if (input$labeling_metric == "fdr") {
          # Label by FDR (smallest p-values)
          top_genes <- data %>%
            filter((is.na(selected_fdr) | padj < selected_fdr) & abs(log2FoldChange) >= selected_fc) %>%
            arrange(padj) %>%
            head(input$n_labels)
        } else {
          # Label by Fold Change (largest absolute log2 fold changes)
          top_genes <- data %>%
            filter((is.na(selected_fdr) | padj < selected_fdr)) %>%
            arrange(desc(abs(log2FoldChange))) %>%
            head(input$n_labels)
        }
      } else {
        if (input$labeling_metric == "fdr") {
          # Label by FDR (smallest p-values) for advanced labeling
          top_genes_left <- data %>%
            filter(log2FoldChange < -selected_fc & (is.na(selected_fdr) | padj < selected_fdr)) %>%
            arrange(padj) %>%
            head(input$n_labels_left)
          
          top_genes_right <- data %>%
            filter(log2FoldChange > selected_fc & (is.na(selected_fdr) | padj < selected_fdr)) %>%
            arrange(padj) %>%
            head(input$n_labels_right)
        } else {
          # Label by Fold Change (largest absolute log2 fold changes) for advanced labeling
          top_genes_left <- data %>%
            filter(log2FoldChange < -selected_fc & (is.na(selected_fdr) | padj < selected_fdr)) %>%
            arrange(desc(abs(log2FoldChange))) %>%
            head(input$n_labels_left)
          
          top_genes_right <- data %>%
            filter(log2FoldChange > selected_fc & (is.na(selected_fdr) | padj < selected_fdr)) %>%
            arrange(desc(abs(log2FoldChange))) %>%
            head(input$n_labels_right)
        }
        
        # Combine left and right labeled genes
        top_genes <- bind_rows(top_genes_left, top_genes_right)
      }
      
      # Add the general top gene labels
      p <- p + geom_text_repel(data = top_genes, aes(label = Gene),
                               max.overlaps = 1000, color = "black")
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
  
  observeEvent(input$add_de_genes_tissue1, {
    data <- unfiltered_data() %>%
      mutate(Gene = toupper(Gene)) 
    
    selected_fdr <- as.numeric(input$selected_fdr)
    selected_fc <- as.numeric(input$selected_fc)
    order <- order_val()
    
    if (order == 0) {
      de_genes_tissue1 <- data %>%
        filter(padj < selected_fdr & log2FoldChange > selected_fc) %>%
        pull(Gene)
    } else {
      de_genes_tissue1 <- data %>%
        filter(padj < selected_fdr & log2FoldChange < -selected_fc) %>%
        pull(Gene)
    }
    
    if (length(de_genes_tissue1) > 0) {
      add_to_cart(de_genes_tissue1)
      showNotification(paste("Grabbed", length(de_genes_tissue1), "DE genes upregulated in", input$tissue_select1, "."), type = "message")
    } else {
      showNotification(paste("No DE genes found upregulated in", input$tissue_select1), type = "warning")
    }
  })
  
  output$de_genes_tissue1_button <- renderUI({
    req(input$tissue_select1)
    actionButton("add_de_genes_tissue1", paste("Grab all ppregulated genes from", input$tissue_select1))
  })
  
  output$de_genes_tissue2_button <- renderUI({
    req(input$tissue_select2)
    actionButton("add_de_genes_tissue2", paste("Grab all upregulated genes from", input$tissue_select2))
  })
  
  

  observeEvent(input$add_de_genes_tissue2, {
    data <- unfiltered_data() %>%
      mutate(Gene = toupper(Gene))  # Ensure genes are uppercase
    
    selected_fdr <- as.numeric(input$selected_fdr)
    selected_fc <- as.numeric(input$selected_fc)
    order <- order_val()
    
    # Determine if we are adding genes upregulated in Tissue 2
    if (order == 0) {
      # Upregulated in Tissue 2: log2FoldChange < -selected_fc
      de_genes_tissue2 <- data %>%
        filter(padj < selected_fdr & log2FoldChange < -selected_fc) %>%
        pull(Gene)
    } else {
      # Upregulated in Tissue 2 when order == 1 (log2FoldChange > selected_fc)
      de_genes_tissue2 <- data %>%
        filter(padj < selected_fdr & log2FoldChange > selected_fc) %>%
        pull(Gene)
    }
    
    if (length(de_genes_tissue2) > 0) {
      add_to_cart(de_genes_tissue2)
      showNotification(paste("Grabbed", length(de_genes_tissue2), "DE genes upregulated in", input$tissue_select2), type = "message")
    } else {
      showNotification(paste("No DE genes found upregulated in", input$tissue_select2), type = "warning")
    }
  })
  
  # Event to add labeled genes to the gene cart
  observeEvent(input$add_labeled_genes, {
    data <- unfiltered_data() %>%
      mutate(Gene = toupper(Gene))  # Ensure genes are uppercase
    
    selected_fdr <- as.numeric(input$selected_fdr)
    selected_fc <- as.numeric(input$selected_fc)
    
    # Logic to identify labeled genes based on labeling option (basic or advanced)
    if (input$labeling_option == "basic") {
      # For basic labeling, select top genes based on the selected metric (FDR or FC)
      if (input$labeling_metric == "fdr") {
        # Label by FDR (smallest p-values)
        labeled_genes <- data %>%
          filter((padj < selected_fdr) & abs(log2FoldChange) >= selected_fc) %>%
          arrange(padj) %>%
          head(input$n_labels) %>%
          pull(Gene)
      } else {
        # Label by Fold Change (largest absolute log2 fold changes)
        labeled_genes <- data %>%
          filter((padj < selected_fdr)) %>%
          arrange(desc(abs(log2FoldChange))) %>%
          head(input$n_labels) %>%
          pull(Gene)
      }
    } else {
      # For advanced labeling, select top genes separately for left and right of the volcano plot
      if (input$labeling_metric == "fdr") {
        # Label by FDR (smallest p-values) for advanced labeling
        labeled_genes_left <- data %>%
          filter(log2FoldChange < -selected_fc & padj < selected_fdr) %>%
          arrange(padj) %>%
          head(input$n_labels_left) %>%
          pull(Gene)
        
        labeled_genes_right <- data %>%
          filter(log2FoldChange > selected_fc & padj < selected_fdr) %>%
          arrange(padj) %>%
          head(input$n_labels_right) %>%
          pull(Gene)
      } else {
        # Label by Fold Change (largest absolute log2 fold changes) for advanced labeling
        labeled_genes_left <- data %>%
          filter(log2FoldChange < -selected_fc & padj < selected_fdr) %>%
          arrange(desc(abs(log2FoldChange))) %>%
          head(input$n_labels_left) %>%
          pull(Gene)
        
        labeled_genes_right <- data %>%
          filter(log2FoldChange > selected_fc & padj < selected_fdr) %>%
          arrange(desc(abs(log2FoldChange))) %>%
          head(input$n_labels_right) %>%
          pull(Gene)
      }
      
      # Combine left and right labeled genes
      labeled_genes <- unique(c(labeled_genes_left, labeled_genes_right))
    }
    
    # Add the labeled genes to the gene cart
    if (length(labeled_genes) > 0) {
      add_to_cart(labeled_genes)
      showNotification(paste("Grabbed", length(labeled_genes), "labeled genes."), type = "message")
    } else {
      showNotification("No labeled genes found.", type = "warning")
    }
  })
  
  
  
  
  ### Heatmap
  
  rpkm_filtered <- reactiveVal()
  ordered_genes <- reactiveVal(NULL) 
  
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
      gene_data <- data %>% filter(toupper(Gene) %in% toupper(searched_genes))
      not_found_genes <- setdiff(toupper(searched_genes), toupper(gene_data$Gene))
      
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
      dplyr::select(Sample_name, Tissue) %>%
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
    
    rpkm_filtered(rpkm_filtered)
    
    heatmap <- pheatmap(
      rpkm_filtered, 
      scale = "row",
      cluster_rows = (n_genes > 1), 
      cluster_cols = TRUE, 
      show_rownames = show_rownames, 
      show_colnames = TRUE,
      annotation_col = annotation,
      annotation_colors = annotation_colors,
      color = colorRampPalette(color_scale)(100),
      legend_labels = c("low", "medium", "high"),
      cutree_cols = 2,
      cutree_rows = 1,
      angle_col = 45,
      fontsize = 13,
      fontsize_col = 14,
      fontsize_row = fontsize_row,
      width = 30,
      treeheight_row = 0,
      treeheight_col = 0
    )
    
    ordered_genes(heatmap$tree_row$labels[heatmap$tree_row$order])
    
    rpkm_filtered(rpkm_filtered)
    
    heatmap
  }
  
  output$heatmap_plot <- renderPlot({
    data <- filtered_data_combined()
    req(data)
    search_genes <- search_genes_heatmap()
    generate_plot(data, search_genes)
  })
  
  observeEvent(input$heatmap_brush, {
    brush_data <- input$heatmap_brush
    
    if (!is.null(brush_data)) {
      heatmap_rows <- ordered_genes()
      total_height <- length(heatmap_rows)
      
      start_row_index <- total_height - ceiling(brush_data$ymin / (1 / total_height)) + 1
      end_row_index <- total_height - ceiling(brush_data$ymax / (1 / total_height)) + 1
      
      start_row_index <- max(min(start_row_index, total_height), 1)
      end_row_index <- max(min(end_row_index, total_height), 1)
      
      selected_rows <- heatmap_rows[end_row_index:start_row_index]
      
      selected_genes(selected_rows)

      showModal(modalDialog(
        title = "Selected Genes",
        paste(paste(selected_rows, collapse = ", ")),
        footer = tagList(
          actionButton("add_to_cart_tissue_specific_modal", "Grab Genes"),
          modalButton("Close")
        ),
        easyClose = TRUE
      ))
    }
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
  
  observeEvent(input$add_to_cart_heatmap, {
    
    add_to_cart(ordered_genes())
    
    showNotification("Genes grabbed.", type = "message")
  })
  
  
  ### Correlation Plot
  
  debounced_correlation_mode <- debounce(reactive(input$correlation_mode), 1500)
  debounced_correlation_system <- debounce(reactive(input$correlation_system), 3000)
  
  output$title_cor_plot <- renderUI({
    req(input$tissue_select1, input$tissue_select2)
    
    if (debounced_correlation_mode() == "selected") {
      h3(paste("Correlation Plot:", input$tissue_select1, "and", input$tissue_select2))
    } else if (debounced_correlation_mode() == "system") {
      h3(paste("Correlation Plot:", input$correlation_system, "system"))
    }
  })
  
  generate_correlation_plot <- reactive({
    req(debounced_correlation_mode())
    
    selected_samples <- NULL
    if (debounced_correlation_mode() == "selected") {
      req(input$tissue_select1, input$tissue_select2)
      selected_tissues <- c(input$tissue_select1, input$tissue_select2)
      selected_samples <- phenodata %>%
        filter(Tissue %in% selected_tissues) %>%
        arrange(Tissue) %>%
        pull(Sample_name2)
    } else if (debounced_correlation_mode() == "system") {
      req(input$correlation_system)
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
    
    if (input$correlation_mode == "selected" && isTRUE(input$use_filtered_data)) {
      req(filtered_data_combined())
      filtered_genes <- filtered_data_combined()$Gene
      if (length(filtered_genes) == 0) {
        showNotification("No genes available after filtering. Using top variable genes instead.", type = "warning")
        N <- 1000
        gene_variability <- apply(RPKM_data, 1, var)
        top_genes <- names(sort(gene_variability, decreasing = TRUE)[1:N])
      } else {
        top_genes <- filtered_genes
      }
    } else {
      N <- 1000
      gene_variability <- apply(RPKM_data, 1, var)
      top_genes <- names(sort(gene_variability, decreasing = TRUE)[1:N])
    }
    
    sample_names1 <- phenodata$Sample_name[match(selected_samples, phenodata$Sample_name2)]
    
    filtered_RPKM <- RPKM_data[rownames(RPKM_data) %in% top_genes, , drop = FALSE]
    selected_data <- filtered_RPKM[, sample_names1, drop = FALSE]
    colnames(selected_data) <- selected_samples
    
    create_diag <- function(data, mapping, ...) {
      varname <- as_label(mapping$x)
      tissue <- phenodata$Tissue[match(varname, phenodata$Sample_name2)]
      color <- mycolors1[tissue]
      ggplot(data = data, mapping = mapping) +
        geom_density() +
        theme_void() +
        ggplot2::theme(panel.background = element_rect(fill = color, colour = "black", size = 2))
    }
    
    custom_ggally_cor <- function(data, mapping, ...) {
      x <- eval_data_col(data, mapping$x)
      y <- eval_data_col(data, mapping$y)
      cor_value <- cor(x, y, use = "complete.obs")
      label <- format(cor_value, digits = 2)
      
      tile_data <- data.frame(cor_value = cor_value)
      
      ggplot(data = tile_data) +
        geom_tile(aes(x = 0.5, y = 0.5, fill = cor_value), width = -Inf, height = Inf, alpha = 1) +
        geom_text(aes(x = 0.5, y = 0.5, label = label), size = 5.5, hjust = 0.5, vjust = 0.5) +
        scale_fill_gradient2(low = "blue", mid = "white", high = "firebrick3", limits = c(-1, 1), na.value = NA) +
        theme_void() +
        ggplot2::theme(legend.position = "none")
    }
    
    create_upper <- function(data, mapping, ...) {
      custom_ggally_cor(data = data, mapping = mapping) +
        theme_minimal() +
        ggplot2::theme(plot.title = element_blank(), axis.title = element_blank())
    }
    
    ggpairs(
      data = as.data.frame(selected_data),
      upper = list(continuous = create_upper),
      diag = list(continuous = create_diag),
      lower = list(continuous = wrap("points", color = "black", alpha = 0.5))
    ) + ggplot2::theme(
      legend.position = "none",
      legend.box = "horizontal",
      legend.box.just = "center",
      legend.text = element_text(size = 14),
      legend.title = element_text(size = 16),
      legend.background = element_rect(fill = "#f0f0f0", color = "black"),
      legend.key = element_rect(fill = "white", color = "black"),
      legend.key.size = unit(1.5, "lines")
    )
  })
  
  
  output$correlation_plot <- renderPlot({
    req(generate_correlation_plot())
    print(generate_correlation_plot())
  })
  
  output$download_correlation_plot <- downloadHandler(
    filename = function() {
      paste("Correlation", ".png", sep = "")
    },
    content = function(file) {
      ggsave(file, plot = generate_correlation_plot(), device = "png", width = 14, height = 10)
    }
  )
  
  
  
  ### PCA
  
  output$title_pca_plot <- renderUI({
    req(input$tissue_select1, input$tissue_select2)
    
    if (input$pca_mode == "selected") {
      h3(paste("PCA:", input$tissue_select1, "and", input$tissue_select2))
    } else if (input$pca_mode == "system") {
      h3(paste("PCA:", input$pca_system_select, "System"))
    }
    else if (input$pca_mode == "all") {
      h3(paste("PCA: All Tissues"))
    }
  })
  
  observe({
    selected_samples <- NULL
    
    if (input$pca_mode == "selected") {
      selected_tissues <- c(input$tissue_select1, input$tissue_select2)
      selected_samples <- phenodata %>%
        filter(Tissue %in% selected_tissues) %>%
        arrange(Tissue) %>%
        pull(Sample_name2)
      
    } else if (input$pca_mode == "system") {
      selected_system <- input$pca_system_select
      selected_tissues <- phenodata %>%
        filter(System == selected_system) %>%
        pull(Tissue) %>%
        unique()
      
      selected_samples <- phenodata %>%
        filter(Tissue %in% selected_tissues) %>%
        arrange(Tissue) %>%
        pull(Sample_name2)
      
    } else if (input$pca_mode == "all") {
      selected_samples <- phenodata %>%
        pull(Sample_name2)
    }
    
    req(selected_samples)
    
    N <- 1000
    gene_variability <- apply(RPKM_data, 1, var)
    top_genes <- names(sort(gene_variability, decreasing = TRUE)[1:N])
    
    sample_names1 <- phenodata$Sample_name[match(selected_samples, phenodata$Sample_name2)]
    
    filtered_RPKM <- RPKM_data[rownames(RPKM_data) %in% top_genes, , drop = FALSE]
    
    selected_data <- filtered_RPKM[, sample_names1, drop = FALSE]
    
    colnames(selected_data) <- selected_samples
    
    sample_tissues <- phenodata$Tissue[match(selected_samples, phenodata$Sample_name2)]
    
    output$pca_plot <- renderPlot({
      
      pca_data <- prcomp(t(selected_data), scale. = F)
      pca_df <- as.data.frame(pca_data$x)
      pca_df$Sample <- rownames(pca_df)
      pca_df <- merge(pca_df, phenodata, by.x = "Sample", by.y = "Sample_name2")
      
      ggplot_pca <- ggplot(pca_df, aes(x = PC1, y = PC2, color = Tissue, shape = System)) +
        geom_point(size = 4) +
        geom_polygon(data = pca_df, aes(x = PC1, y = PC2, group = Tissue, fill = Tissue), alpha = 0.2) +
        scale_color_manual(values = mycolors1) +
        labs(x = "Principal Component 1", y = "Principal Component 2") +
        theme_minimal()
      
      print(ggplot_pca)
      
      output$download_pca_plot <- downloadHandler(
        filename = function() {
          paste("PCA", ".png", sep = "")
        },
        content = function(file) {
          ggsave(file, plot = ggplot_pca, device = "png", width = 10, height = 8, bg = "white")
        }
      )
    })
  })
  
  ### Gene Expression
  
  options(ragg.max_dim = 10000000)
  
  barplot_data <- reactiveVal(NULL)
  current_plot <- reactiveVal(NULL)
  
  # Add a reactive value to keep track of the current title
  current_title <- reactiveVal("")
  
  common_genes <- function(tissue1, tissue2, RPKM_data, phenodata) {
    # Filter phenotype data to get the samples belonging to the selected tissues
    tissue1_samples <- phenodata %>%
      filter(Tissue == tissue1) %>%
      pull(Sample_name)
    
    tissue2_samples <- phenodata %>%
      filter(Tissue == tissue2) %>%
      pull(Sample_name)
    
    # Select the columns in RPKM_data that match the filtered sample names
    selected_samples <- RPKM_data %>%
      dplyr::select(all_of(c(tissue1_samples, tissue2_samples)))
    
    # Calculate the average expression across the selected tissues for each gene
    top_genes <- selected_samples %>%
      rowMeans(na.rm = TRUE) %>%
      sort(decreasing = TRUE) %>%
      head(5) %>%
      names()
    
    return(top_genes)
  }
  
  observeEvent(input$differential_analysis_tab, {
    
    if (input$differential_analysis_tab == "Gene Expression") {
      # Ensure data is available before proceeding
      req(input$tissue_select1, input$tissue_select2, RPKM_data, phenodata)
      
      top_genes <- common_genes(input$tissue_select1, input$tissue_select2, RPKM_data, phenodata)
      
      if (length(top_genes) > 0) {
        updateTextInput(session, "barplot_searched_gene", value = paste(top_genes, collapse = ", "))
        
        # Add a check to force render if needed
        if (!is.null(top_genes) && length(top_genes) > 0) {
          updateBarplotData(top_genes)
        }
        
        current_title(paste("Showing Top 5 Commonly Expressed Genes in", input$tissue_select1, "and", input$tissue_select2))
      } else {
        showNotification("No common genes found between selected tissues.", type = "error")
      }
    }
  })
  
  
  
  
  observeEvent(input$search_gene_barplot, {
    req(input$barplot_searched_gene)
    
    searched_genes <- strsplit(input$barplot_searched_gene, ",\\s*")[[1]]
    searched_genes <- toupper(searched_genes)
    
    found_genes <- intersect(toupper(rownames(RPKM_data)), searched_genes)
    
    if (length(found_genes) == 0) {
      showNotification("No genes found.", type = "warning")
      barplot_data(NULL)
      current_plot(NULL)
      return(NULL)
    }
    
    updateBarplotData(found_genes)
    
    # Update the title
    current_title("Showing Searched Genes")
  })
  
  observeEvent(input$random_gene_barplot, {
    random_genes <- sample(rownames(RPKM_data), 5)
    updateTextInput(session, "barplot_searched_gene", value = paste(random_genes, collapse = ", "))
    updateBarplotData(random_genes)
    
    # Update the title
    current_title("Showing searched genes")
  })
  
  observeEvent(input$clear_search_gene_barplot, {
    updateTextInput(session, "barplot_searched_gene", value = "")
    barplot_data(NULL)
    current_plot(NULL)
  })
  
  output$dynamic_title <- renderUI({
    h2(current_title(), style = "text-align: center; margin-top: 20px;")
  })
  
  
  updateBarplotData <- function(genes) {
    if (is.null(genes) || length(genes) == 0) {
      showNotification("No valid genes to plot.", type = "error")
      output$dynamic_barplot_output <- renderUI({
        h4("No data available for the selected genes and tissues.", style = "text-align: center; color: red;")
      })
      return(NULL)
    }
    
    selected_data <- RPKM_data[genes, , drop = FALSE]
    
    if (nrow(selected_data) == 0 || ncol(selected_data) == 0) {
      showNotification("No valid data available for the selected genes and tissues.", type = "error")
      output$dynamic_barplot_output <- renderUI({
        h4("No data available for the selected genes and tissues.", style = "text-align: center; color: red;")
      })
      return(NULL)
    }
    
    selected_data$Gene <- rownames(selected_data)
    melted_data <- reshape2::melt(selected_data)
    tissue_info <- phenodata[match(melted_data$variable, phenodata$Sample_name), "Tissue"]
    melted_data$Tissue <- tissue_info
    
    if (nrow(melted_data) == 0 || all(is.na(melted_data$Tissue))) {
      showNotification("No valid data available for the selected genes and tissues.", type = "error")
      output$dynamic_barplot_output <- renderUI({
        h4("No data available for the selected genes and tissues.", style = "text-align: center; color: red;")
      })
      return(NULL)
    }
    
    average_expression_per_tissue <- melted_data %>%
      group_by(Gene, Tissue) %>%
      summarize(average_value = mean(value, na.rm = TRUE), .groups = 'drop') %>%
      ungroup()
    
    barplot_data(average_expression_per_tissue)
    
    output$dynamic_barplot_output <- renderUI({
      plot_height <- 300 + 300 * length(genes)
      plotOutput("barplot_plot", height = paste0(plot_height, "px"))
    })
    
    output$barplot_plot <- renderPlot({
      p <- ggplot(average_expression_per_tissue, aes(x = Tissue, y = average_value, fill = Tissue)) +
        geom_bar(stat = "identity", color = "black", width = 0.8) +
        facet_wrap(~ Gene, scales = "free", ncol = 1) +
        theme_minimal() +
        ggplot2::theme(
          axis.text.x = element_text(angle = 45, hjust = 1, size = 12),
          panel.spacing = unit(0.5, "lines"),
          panel.border = element_rect(color = "black", fill = NA, size = 0.5),
          strip.background = element_blank(),
          strip.placement = "outside",
          panel.background = element_blank(),
          strip.text = element_text(size = 14, face = "bold"),
          legend.position = "none"
        ) +
        scale_fill_manual(values = mycolors1) +
        labs(y = "Average Expression (Avg. RPKM)", x = "Tissue")
      
      current_plot(p)
      print(p)
    })
  }
  
  output$download_barplot <- downloadHandler(
    filename = function() {
      paste("Barplot", ".png", sep = "")
    },
    content = function(file) {
      ggsave(file, plot = current_plot(), device = "png", width = 10, height = length(barplot_data()) * 4, bg = "white")
    }
  )
  
  observeEvent(input$add_to_cart_barplot, {
    genes <- unique(barplot_data()$Gene)
    add_to_cart(genes)
    
    showNotification("Genes grabbed.", type = "message")
  })

  observeEvent(input$transfer_to_heatmap, {

    transferred_genes <- strsplit(input$barplot_searched_gene, ",\\s*")[[1]]
    transferred_genes <- toupper(transferred_genes)
    
    rpkm_filtered <- RPKM_data %>%
      filter(rownames(.) %in% transferred_genes)
    
    samples <- phenodata %>%
      pull(Sample_name)
    
    samples <- samples[samples %in% colnames(rpkm_filtered)]
    
    rpkm_filtered <- rpkm_filtered[, samples, drop = FALSE]
    
    colnames(rpkm_filtered) <- gsub("Sample_", "", colnames(rpkm_filtered))
    
    annotation <- phenodata %>%
      filter(Sample_name %in% samples) %>%
      dplyr::select(Sample_name, Tissue) %>%
      mutate(Sample_name = gsub("Sample_", "", Sample_name)) %>%
      column_to_rownames(var = "Sample_name")
    
    present_tissues <- unique(annotation$Tissue)
    annotation_colors <- list(Tissue = mycolors1[present_tissues])
    
    sorted_samples <- annotation %>%
      arrange(Tissue) %>%
      rownames()
    
    rpkm_filtered <- rpkm_filtered[, sorted_samples]
    
    num_genes <- length(transferred_genes)
    cell_height <- ifelse(num_genes > 50, 5, ifelse(num_genes > 30, 8, ifelse(num_genes > 20, 15, ifelse(num_genes > 10, 25, ifelse(num_genes > 5, 50, 100)))))
    
    show_genes <- num_genes <= 30
    
    # Heatmap generation
    showModal(
      modalDialog(
        title = "Gene Expression Heatmap for Searched Genes",
        size = "l",  # This will make the modal larger
        plotOutput("heatmap_in_modal", width = "100%", height = "750px"),
        downloadButton("download_modal_heatmap", "Download Heatmap"),
        easyClose = TRUE,
        footer = modalButton("Close"),
        tags$style(".modal-dialog { width: 80%; }")
      )
    )
    
    output$heatmap_in_modal <- renderPlot({

      pheatmap(
        rpkm_filtered,
        scale = "row",
        cluster_rows = TRUE,
        cluster_cols = FALSE,
        show_rownames = show_genes,
        show_colnames = TRUE,
        annotation_col = annotation,
        annotation_colors = annotation_colors,
        color = colorRampPalette(c("royalblue", "white", "firebrick3"))(100),  # Example of color scale
        legend_labels = c("low", "medium", "high"),
        fontsize_col = 11,
        fontsize_row = 12,
        angle_col = 45,
        treeheight_row = 0,
        treeheight_col = 0,
        cellwidth = 23,
      )
    })
    
    # Add the download handler for the heatmap
    output$download_modal_heatmap <- downloadHandler(
      filename = function() {
        paste("heatmap", ".png", sep = "")
      },
      content = function(file) {
        # Save the heatmap as a PNG file
        png(file, width = 15, height = 10, units = "in", res = 300)
        pheatmap(
          rpkm_filtered,
          scale = "row",
          cluster_rows = TRUE,
          cluster_cols = FALSE,  # Prevent clustering of columns to keep pairs together
          show_rownames = TRUE,
          show_colnames = TRUE,
          annotation_col = annotation,
          annotation_colors = annotation_colors,
          color = colorRampPalette(c("royalblue", "white", "firebrick3"))(100),  # Example of color scale
          legend_labels = c("low", "medium", "high"),
          fontsize = 10,
          angle_col = 45,
          treeheight_row = 0,
          treeheight_col = 0,
          cellwidth = 17,
          cellheight = cell_height
        )
        dev.off()
      }
    )
  })
  
  ##########################
  
  # Tissue specific analysis
  
  selected_genes <- reactiveVal()
  rpkm_heatmap_ordered <- reactiveVal()
  ordered_tissue_specific_genes <- reactiveVal(NULL)
  
  output$title_tissue_specific <- renderUI({
    h3("Tissue-Specific Genes")
  })
  
  
  observeEvent(input$system_choice, {
    all_systems_selected <- setdiff(unique(phenodata$System), "All Systems")
    
    if ("All Systems" %in% input$system_choice && length(input$system_choice) > 1) {
      updateCheckboxGroupInput(session, "system_choice", 
                               selected = setdiff(input$system_choice, "All Systems"))
    }
    
    if (all(all_systems_selected %in% input$system_choice)) {
      updateCheckboxGroupInput(session, "system_choice", 
                               selected = "All Systems")
    }
  })
  
  output$tissue_selection_ui <- renderUI({
    req(input$system_choice)
    
    if ("All Systems" %in% input$system_choice) {
      tissues <- unique(phenodata$Tissue)
    } else {
      tissues <- unique(phenodata$Tissue[phenodata$System %in% input$system_choice])
    }
    
    checkboxGroupInput("tissue_choice", "Select Tissue(s):", choices = tissues, selected = tissues)
  })
  
  output$tissue_specific_heatmap <- renderPlot({
    req(input$tissue_choice)
    current_plot(NULL)
    
    selected_tissues <- input$tissue_choice
    
    samples <- phenodata %>%
      filter(Tissue %in% selected_tissues) %>%
      pull(Sample_name)
    
    samples <- samples[samples %in% colnames(RPKM_data)]
    
    rpkm_filtered <- RPKM_data[, samples, drop = FALSE]
    
    non_zero_genes <- apply(rpkm_filtered, 1, function(row) !all(row == 0))
    rpkm_filtered <- rpkm_filtered[non_zero_genes, ]
    
    colnames(rpkm_filtered) <- gsub("Sample_", "", colnames(rpkm_filtered))
    
    tau_scores <- calcTau(rpkm_filtered)
    
    filtered_tau <- tau_scores[tau_scores$tau >= 0.85, ]
    rpkm_filtered <- rpkm_filtered[rownames(rpkm_filtered) %in% rownames(filtered_tau), ]
    
    rpkm_tau <- as.matrix(rpkm_filtered)
    
    variances <- apply(t(rpkm_tau), 2, var)

      top_genes <- order(variances, decreasing = TRUE)[1:input$n_variable_genes]
      rpkm_heatmap <- rpkm_tau[top_genes, ]
    
    
    annotation <- phenodata %>%
      filter(Sample_name %in% samples) %>%
      dplyr::select(Sample_name, Tissue, System) %>%
      mutate(Sample_name = gsub("Sample_", "", Sample_name)) %>%
      column_to_rownames(var = "Sample_name")
    
    present_tissues <- unique(annotation$Tissue)
    
    annotation_colors <- list(
      Tissue = mycolors1[unique(annotation$Tissue)],
      System = mycolors2[unique(annotation$System)]
      )
    
    n_genes <- nrow(rpkm_tau)
    
    color_scale_tissue_specific <- switch(input$color_scale_tissue_specific,
                                          "blue_white_red" = c("royalblue", "white", "firebrick3"),
                                          "green_black_red" = c("springgreen2", "black", "firebrick2"),
                                          "purple_white_green" = c("purple", "white", "springgreen4"),
                                          "cyan_yellow_purple" = c("purple", "lightyellow", "blue"),
                                          "viridis" = viridis(100),
                                          "plasma" = plasma(100),
                                          "cividis" = cividis(100),
                                          "inferno" = inferno(100))
    
    annotation_order <- annotation %>%
      arrange(System, Tissue) %>%
      rownames()
    
    rpkm_heatmap_ordered <- rpkm_heatmap[, annotation_order]
    
    rpkm_heatmap_ordered(rpkm_heatmap_ordered)
    
    show_row_names <- nrow(rpkm_heatmap_ordered) <= 50
    
    heatmap <- pheatmap(
      rpkm_heatmap_ordered, 
      annotation_col = annotation[annotation_order, ],
      annotation_colors = annotation_colors,
      color = colorRampPalette(color_scale_tissue_specific)(100),
      clustering_distance_rows = "correlation",
      clustering_distance_cols = "correlation",
      clustering_method = "complete",
      scale = "row",
      drop_levels = TRUE,
      cluster_rows = TRUE,
      cluster_cols = FALSE,
      show_rownames = show_row_names,
      show_colnames = TRUE,
      legend_labels = c("low", "medium", "high"),
      legend = TRUE,
      angle_col = 45,
      fontsize = 12,
      fontsize_col = 9,
      treeheight_col = 0
    )
    
    ordered_tissue_specific_genes(heatmap$tree_row$labels[heatmap$tree_row$order])
    
    current_plot(heatmap)
    
    print(heatmap)
  })
  
  observeEvent(input$heatmap_tissue_specific_brush, {
    brush_data <- input$heatmap_tissue_specific_brush
    
    if (!is.null(brush_data)) {
      heatmap_rows <- ordered_tissue_specific_genes()
      total_height <- length(heatmap_rows)
      
      start_row_index <- total_height - ceiling(brush_data$ymin * total_height) + 1
      end_row_index <- total_height - ceiling(brush_data$ymax * total_height) + 1
      
      start_row_index <- max(min(start_row_index, total_height), 1)
      end_row_index <- max(min(end_row_index, total_height), 1)

      selected_rows <- heatmap_rows[end_row_index:start_row_index]
      
      selected_genes(selected_rows)
      
      showModal(modalDialog(
        title = "Selected Genes",
        paste(paste(selected_rows, collapse = ", ")),
        footer = tagList(
          actionButton("add_to_cart_tissue_specific_modal", "Grab Genes"),
          modalButton("Close")
        ),
        easyClose = TRUE
      ))
    }
  })
  
  output$download_tissue_specific_heatmap <- downloadHandler(
    filename = function() {
      paste("Tissue_specific_heatmap", ".png", sep = "")
    },
    content = function(file) {
      ggsave(file, plot = current_plot(), device = "png", width = 10, height = 11, bg = "white")
    }
  )
  
  observeEvent(input$add_to_cart_tissue_specific_modal, {
    genes_to_add <- selected_genes()
    
    if (!is.null(genes_to_add)) {
      add_to_cart(genes_to_add)
      
      showNotification("Genes grabbed.", type = "message")
      
      removeModal()
    }
  })
  
  
  observeEvent(input$add_to_cart_tissue_specific_heatmap, {
    
    add_to_cart(ordered_tissue_specific_genes())
    
    showNotification("Genes grabbed.", type = "message")
  })
  
  ##############
  
  # Genome Browser

  options(shiny.maxRequestSize = 1024 * 1024^2) 
  
  convert_lists_to_chars <- function(df) {
    df[] <- lapply(df, function(x) if (is.list(x)) sapply(x, toString) else x)
    return(df)
  }
  
  observeEvent(input$genome_select, {
    genome <- input$genome_select
    options <- parseAndValidateGenomeSpec(genomeName = genome)
    
    output$igvShiny <- renderIgvShiny({
      igvShiny(options)
    })
  })
  
  observeEvent(input$search_button, {
    req(input$genome_browser_search)
    showGenomicRegion(session = session, id = "igvShiny", region = input$genome_browser_search)
  })
  
  chrom_map_path <- "browser_data_for_app/chrm-ncbi.txt"
  chrom_map <- NULL
  
  if (file.exists(chrom_map_path)) {
    chrom_map <- read_tsv(chrom_map_path, col_names = c("Chrm_ncbi", "Chromosome"))
  } else {
    showNotification("Chromosome mapping file not found.", type = "error")
  }
  
  observeEvent(input$loadTrack, {
    colorTable <- list("gene" = "blue")
    colorByAttribute <- "type"
    displayMode <- "EXPANDED"
    visibilityWindow <- 1000000
    
    withProgress(message = 'Loading track...', value = 0, {
      incProgress(0.25)
      
      if (input$input_type == "file" && !is.null(input$file)) {
        file_path <- input$file$datapath
        file_type <- tools::file_ext(input$file$name)
        
        tryCatch({
          if (file_type == "gff3" || file_type == "gff") {
            tbl.gff <- rtracklayer::import(file_path, format = file_type)
            tbl.gff <- convert_lists_to_chars(as.data.frame(tbl.gff))
            loadGFF3TrackFromLocalData(
              session = session,
              id = "igvShiny",
              trackName = input$file$name,
              tbl.gff3 = tbl.gff,
              color = "gray",
              colorTable = colorTable,
              colorByAttribute = colorByAttribute,
              displayMode = displayMode,
              trackHeight = 50,
              visibilityWindow = visibilityWindow,
              deleteTracksOfSameName = TRUE
            )
            showNotification("GFF/GFF3 Track loaded successfully", type = "message")
            
          } else if (file_type == "bam") {
            bamFile <- Rsamtools::BamFile(file_path)
            loadBamTrackFromLocalData(
              session = session,
              id = "igvShiny",
              trackName = "Uploaded BAM Track",
              bamFile = bamFile,
              color = "grey",
              trackHeight = 50,
              displayMode = displayMode,
              deleteTracksOfSameName = TRUE
            )
            showNotification("BAM Track loaded successfully", type = "message")
            
          } else if (file_type == "bed") {
            tbl.bed <- rtracklayer::import(file_path, format = "bed")
            tbl.bed <- as.data.frame(tbl.bed)
            
            if (!is.character(tbl.bed$seqnames)) {
              tbl.bed$seqnames <- as.character(tbl.bed$seqnames)
            }
            
            colnames(tbl.bed)[colnames(tbl.bed) == "seqnames"] <- "chr"
            
            loadBedTrack(
              session = session,
              id = "igvShiny",
              trackName = input$file$name,
              tbl = tbl.bed,
              color = "gray",
              trackHeight = 50,
              deleteTracksOfSameName = TRUE
            )
            showNotification("BED Track loaded successfully", type = "message")
            
          } else if (file_type == "vcf") {
            vcf_data <- VariantAnnotation::readVcf(file_path)
            loadVcfTrack(
              session = session,
              id = "igvShiny",
              trackName = input$file$name,
              vcfData = vcf_data,
              deleteTracksOfSameName = TRUE
            )
            showNotification("VCF Track loaded successfully", type = "message")
            
          } else {
            showNotification("Unsupported file type", type = "error")
          }
          incProgress(0.75)
          
        }, error = function(e) {
          showNotification(paste("Error loading track:", e$message), type = "error")
        })
        
      } else if (input$input_type == "url" && !is.null(input$url)) {
        url <- input$url
        file_type <- tools::file_ext(url)
        
        tryCatch({
          if (file_type == "gff3" || file_type == "gff") {
            loadGFF3TrackFromURL(
              session = session,
              id = "igvShiny",
              trackName = "URL GFF/GFF3 Track",
              gff3URL = url,
              color = "gray",
              colorTable = colorTable,
              colorByAttribute = colorByAttribute,
              displayMode = displayMode,
              trackHeight = 50,
              visibilityWindow = visibilityWindow,
              deleteTracksOfSameName = TRUE
            )
            showNotification("GFF/GFF3 Track from URL loaded successfully", type = "message")
            
          } else if (file_type == "bam") {
            loadBamTrackFromURL(
              session = session,
              id = "igvShiny",
              trackName = "URL BAM Track",
              bamURL = url,
              displayMode = displayMode,
              deleteTracksOfSameName = TRUE
            )
            showNotification("BAM Track from URL loaded successfully", type = "message")
            
          } else if (file_type == "bed") {
            loadBedTrackFromURL(
              session = session,
              id = "igvShiny",
              trackName = "URL BED Track",
              url = url,
              color = "gray",
              trackHeight = 50,
              deleteTracksOfSameName = TRUE
            )
            showNotification("BED Track from URL loaded successfully", type = "message")
            
          } else if (file_type == "vcf") {
            vcf_data <- VariantAnnotation::readVcf(url)
            loadVcfTrack(
              session = session,
              id = "igvShiny",
              trackName = "URL VCF Track",
              vcfData = vcf_data,
              deleteTracksOfSameName = TRUE
            )
            showNotification("VCF Track from URL loaded successfully", type = "message")
            
          } else {
            showNotification("Unsupported URL file type", type = "error")
          }
          incProgress(0.75)
          
        }, error = function(e) {
          showNotification(paste("Error loading track from URL:", e$message), type = "error")
        })
        
      } else {
        showNotification("No file or URL provided", type = "warning")
      }
      
      incProgress(1)
    })
  })

  observeEvent(input$load_geneASE_track, {
    gene_ase_path <- "browser_data_for_app/gene_ase_converted_cleaned_with_strand.bed"
    
    if (file.exists(gene_ase_path) && !is.null(chrom_map)) {
      withProgress(message = 'Loading Gene ASE Track...', value = 0, {
        tryCatch({
          # Import the BED file
          gene_ase_data <- rtracklayer::import(gene_ase_path, format = "bed")
          gene_ase_df <- as.data.frame(gene_ase_data)
          
          # Rename seqnames column to chr
          colnames(gene_ase_df)[colnames(gene_ase_df) == "seqnames"] <- "chr"
          
          # Merge chromosome mapping and adjust strand
          gene_ase_df <- gene_ase_df %>%
            left_join(chrom_map, by = c("chr" = "Chrm_ncbi")) %>%
            mutate(chr = Chromosome) %>%
            dplyr::select(-Chromosome) %>%
            mutate(strand = ifelse(strand == "*", "+", strand))
          
          # Rearrange columns: place width after start and end
          gene_ase_df <- gene_ase_df %>%
            dplyr::select(chr, start, end, name, score, strand)  # Rearrange columns as per BED format
          
          # Load BED track
          loadBedTrack(
            session = session,
            id = "igvShiny",
            trackName = "Gene ASE Track",
            tbl = gene_ase_df,
            color = "blue",
            trackHeight = 70,
            deleteTracksOfSameName = TRUE
          )
          
          incProgress(0.5)
        }, error = function(e) {
          showNotification(paste("Error loading Gene ASE Track:", e$message), type = "error")
        })
      })
    }
  })
  
  observeEvent(input$load_snp, {
    showNotification("Loading may take a few minutes", type = "warning", duration = 10)

    withProgress(message = 'Loading VCF file...', value = 0, {
      vcf_file <- "browser_data_for_app/ASE.SNP.gallus.updated02.vcf"

      if (file.exists(vcf_file)) {
        showNotification("VCF file found.", type = "message")

        tryCatch({
          vcf_data <- VariantAnnotation::readVcf(vcf_file, "galGal6")

          snp_data <- rowRanges(vcf_data)
          fixed_fields <- as.data.frame(fixed(vcf_data))
          snp_df <- data.frame(
            chr = as.character(seqnames(snp_data)),
            start = start(snp_data),
            end = end(snp_data),
            REF = fixed_fields$REF,
            ALT = sapply(fixed_fields$ALT, function(x) paste(x, collapse = ","))
          )

          loadBedTrack(
            session = session,
            id = "igvShiny",
            trackName = "SNP ASE Track",
            tbl = snp_df,
            color = "darkgreen",
            trackHeight = 70,
            deleteTracksOfSameName = TRUE
          )

          showNotification("Track loaded successfully.", type = "message")
        }, error = function(e) {
          showNotification(paste("Error loading track into IGV:", e$message), type = "error")
        })
      } else {
        showNotification("VCF file not found.", type = "error")
      }
    })
  })

  observeEvent(input$load_gff, {
    gff_file <- "browser_data_for_app/genomeAnnoatation_gallus.updated.gff"
    chrom_map_path <- "browser_data_for_app/chrm-ncbi.txt"  # Assuming this is the chromosome mapping file path
    
    withProgress(message = 'Loading GFF Track...', value = 0, {
      tryCatch({
        if (file.exists(gff_file)) {
          # Import the GFF file
          gff_data <- rtracklayer::import(gff_file, format = "gff")
          gff_df <- as.data.frame(gff_data)
          
          # Load the chromosome mapping file
          chrom_map <- readr::read_tsv(chrom_map_path, col_names = c("Chrm_ncbi", "Chromosome"))
          
          # Debugging: print chromosome map to ensure it is loaded correctly
          print("Chromosome mapping:")
          print(head(chrom_map))
          
          # Rename seqnames column to chr to match expected column names
          colnames(gff_df)[colnames(gff_df) == "seqnames"] <- "chr"
          
          # Convert chr column from factor to character
          gff_df$chr <- as.character(gff_df$chr)
          
          # Merge the chromosome mapping with GFF dataframe to replace NCBI codes with chromosome names
          gff_df <- gff_df %>%
            left_join(chrom_map, by = c("chr" = "Chrm_ncbi")) %>%
            mutate(chr = Chromosome) %>%  # Replace chr with mapped chromosome names
            dplyr::select(-Chromosome)    # Drop the Chromosome column after mapping
          
          # Debugging: print head of dataframe after chromosome mapping
          print("GFF dataframe after applying chromosome mapping:")
          print(head(gff_df))
          
          # Use existing gene_name column and fill any NA values
          gff_df$gene_name[is.na(gff_df$gene_name)] <- "NA"
          
          # Debugging: print head of dataframe after handling gene_name
          print("GFF dataframe after handling gene_name column:")
          print(head(gff_df))
          
          # Ensure the correct column order for IGV
          gff_df <- gff_df %>%
            dplyr::select(chr, start, end, gene_name, score, strand)  # Reorder columns
          
          # Debugging: print head of dataframe after reordering columns
          print("GFF dataframe after reordering columns:")
          print(head(gff_df))
          
          # Load the track in IGV
          loadBedTrack(
            session = session,
            id = "igvShiny",
            trackName = "Annotation Track",
            tbl = gff_df,
            color = "darkred",
            trackHeight = 70,
            deleteTracksOfSameName = TRUE
          )
          
          incProgress(1)
        } else {
          showNotification("GFF file not found.", type = "error")
        }
      }, error = function(e) {
        print("Error loading GFF Track:")
        print(e$message)
        showNotification(paste("Error loading GFF Track:", e$message), type = "error")
      })
    })
  })
  
  observeEvent(input$grab_all_ase_genes, {
    gene_ase_path <- "browser_data_for_app/gene_ase_converted_cleaned_with_strand.bed"
    
    if (file.exists(gene_ase_path)) {
      withProgress(message = 'Grabbing All ASE Genes...', value = 0, {
        tryCatch({

          gene_ase_data <- rtracklayer::import(gene_ase_path, format = "bed")
          gene_ase_df <- as.data.frame(gene_ase_data)
          
          gene_ase_df <- gene_ase_df %>%
            mutate(chr = as.character(seqnames)) %>%
            dplyr::select(chr, start, end, name, score, strand)
          
          all_ase_genes <- toupper(gene_ase_df$name)
          
          add_to_cart(all_ase_genes)
          
          showNotification("All ASE genes have been grabbed successfully.", type = "message")
          
          incProgress(1)
        }, error = function(e) {
          showNotification(paste("Error grabbing ASE genes:", e$message), type = "error")
        })
      })
    } else {
      showNotification("Gene ASE data not found.", type = "error")
    }
  })
  
  
  ##############
  
  # Genes
  
  add_to_cart <- function(new_genes) {
    current_cart <- cart_genes()
    updated_cart <- unique(c(current_cart, new_genes))
    cart_genes(updated_cart)
  }
  
  observeEvent(input$clear_cart, {
    showModal(
      modalDialog(
        title = "Confirm Clear Cart",
        "Are you sure you want to clear all genes from the cart?",
        easyClose = FALSE,
        footer = tagList(
          modalButton("Cancel"),
          actionButton("confirm_clear_cart", "Yes, Clear Cart")
        )
      )
    )
  })
  
  observeEvent(input$confirm_clear_cart, {
    cart_genes(character(0))
    removeModal()
    showNotification("Gene cart has been cleared.", type = "message")
  })
  
  calcTau_custom <- function(expression_matrix) {
    if (is.vector(expression_matrix)) {
      expression_matrix <- matrix(expression_matrix, nrow = 1)
      rownames(expression_matrix) <- genes
    }
    tau_values <- apply(expression_matrix, 1, function(x) {
      if (all(is.na(x)) || sum(x, na.rm = TRUE) == 0) {
        return(NA)
      }
      x_norm <- x / max(x, na.rm = TRUE)
      tau <- sum(1 - x_norm, na.rm = TRUE) / (length(x) - 1)
      return(tau)
    })
    return(data.frame(gene = rownames(expression_matrix), tau = tau_values))
  }
  
  allele_specific_genes <- reactiveVal()
  
  observe({
    gene_ase_path <- "browser_data_for_app/gene_ase_converted_cleaned_with_strand.bed"
    gene_ase_data <- rtracklayer::import(gene_ase_path, format = "bed")
    gene_ase_df <- as.data.frame(gene_ase_data)
    allele_specific_genes(unique(gene_ase_df$name))
  })
  
  gene_expression_data_tissue_reactive <- reactive({
    genes <- cart_genes()
    if (length(genes) == 0) {
      return(NULL)
    }
    
    gene_expression_data <- RPKM_data[rownames(RPKM_data) %in% genes, , drop = FALSE]
    
    sample_to_tissue <- phenodata %>% dplyr::select(Sample_name, Tissue)
    
    common_samples <- intersect(colnames(gene_expression_data), sample_to_tissue$Sample_name)
    
    gene_expression_data <- gene_expression_data[, common_samples, drop = FALSE]
    sample_to_tissue <- sample_to_tissue[sample_to_tissue$Sample_name %in% common_samples, ]
    
    sample_to_tissue <- sample_to_tissue[match(common_samples, sample_to_tissue$Sample_name), ]
    
    sample_tissue_map <- setNames(sample_to_tissue$Tissue, sample_to_tissue$Sample_name)
    
    gene_expression_data_tissue_list <- lapply(unique(sample_tissue_map), function(tissue) {
      samples_in_tissue <- names(sample_tissue_map)[sample_tissue_map == tissue]
      if (length(samples_in_tissue) == 1) {
        as.matrix(gene_expression_data[, samples_in_tissue, drop = FALSE])
      } else {
        rowMeans(gene_expression_data[, samples_in_tissue, drop = FALSE], na.rm = TRUE)
      }
    })
    
    gene_expression_data_tissue <- do.call(cbind, gene_expression_data_tissue_list)
    colnames(gene_expression_data_tissue) <- unique(sample_tissue_map)
    rownames(gene_expression_data_tissue) <- rownames(gene_expression_data)
    
    gene_expression_data_tissue <- as.matrix(gene_expression_data_tissue)
    
    gene_expression_data_tissue[is.na(gene_expression_data_tissue)] <- 0
    
    return(gene_expression_data_tissue)
  })
  
  allele_specific_genes_reactive <- reactive({
    gene_ase_path <- "browser_data_for_app/gene_ase_converted_cleaned_with_strand.bed"
    gene_ase_data <- rtracklayer::import(gene_ase_path, format = "bed")
    gene_ase_df <- as.data.frame(gene_ase_data)
    
    return(unique(gene_ase_df$name))
  })
  
  output$cart_gene_table <- renderDT({
    gene_expression_data_tissue <- gene_expression_data_tissue_reactive()
    if (is.null(gene_expression_data_tissue)) {
      return(datatable(data.frame(Gene = "Gene list is empty.")))
    }
    
    tau_scores <- calcTau_custom(gene_expression_data_tissue)
    
    tau_values <- tau_scores$tau
    names(tau_values) <- tau_scores$gene
    
    max_tissues <- apply(gene_expression_data_tissue, 1, function(x) {
      colnames(gene_expression_data_tissue)[which.max(x)]
    })
    
    allele_specific_genes <- allele_specific_genes_reactive()
    
    result_df <- data.frame(
      Gene = rownames(gene_expression_data_tissue),
      Tau = format(tau_values[rownames(gene_expression_data_tissue)], scientific = F, digits = input$num_of_digits_tau),
      MaxTissue = max_tissues,
      AlleleSpecific = ifelse(toupper(rownames(gene_expression_data_tissue)) %in% toupper(allele_specific_genes), "Yes", "No"),
      stringsAsFactors = FALSE,
      row.names = NULL
    )
    
    return(datatable(result_df, options = list(pageLength = 10)))
  })
  
  observeEvent(input$copy_cart_genes, {
    genes <- cart_genes()
    allele_specific_genes <- allele_specific_genes_reactive()
    
    if (length(genes) > 0) {
      showModal(modalDialog(
        title = "Copy Grabbed Gene",
        radioButtons("copy_gene_type", "Select Genes to Copy:",
                     choices = c("All Genes" = "all",
                                 "Allele-Specific Genes" = "allele_specific",
                                 "Non-Allele-Specific Genes" = "non_allele_specific"),
                     selected = "all"),
        radioButtons("separator_choice", "Separator:",
                     choices = c("Comma" = ", ", "New Line" = "\n"),
                     selected = ", "),
        textAreaInput("gene_list_display", "Grabbed Genes (Copy with CTRL + A then CTRL + C)", "", rows = 15, cols = 100),
        footer = modalButton("Close")
      ))
      
      observe({
        separator <- input$separator_choice
        gene_type <- input$copy_gene_type
        
        if (!is.null(gene_type)) {
          if (gene_type == "allele_specific") {
            genes_to_copy <- genes[toupper(genes) %in% toupper(allele_specific_genes)]
          } else if (gene_type == "non_allele_specific") {
            genes_to_copy <- genes[!toupper(genes) %in% toupper(allele_specific_genes)]
          } else {
            genes_to_copy <- genes
          }
          
          genes_string <- paste(genes_to_copy, collapse = separator)
          updateTextAreaInput(session, "gene_list_display", value = genes_string)
        }
      })
    } else {
      showNotification("Cart is empty.", type = "warning")
    }
  })
  
  
  
  output$download_results <- downloadHandler(
    filename = function() {
      paste("gene_cart_results.csv", sep = "")
    },
    content = function(file) {
      gene_expression_data_tissue <- gene_expression_data_tissue_reactive()
      if (!is.null(gene_expression_data_tissue)) {
        
        tau_scores <- calcTau_custom(gene_expression_data_tissue)
        
        tau_values <- tau_scores$tau
        names(tau_values) <- tau_scores$gene
        
        max_tissues <- apply(gene_expression_data_tissue, 1, function(x) {
          colnames(gene_expression_data_tissue)[which.max(x)]
        })
        
        allele_specific_genes <- allele_specific_genes_reactive()
        
        result_df <- data.frame(
          Gene = rownames(gene_expression_data_tissue),
          Tau = format(tau_values[rownames(gene_expression_data_tissue)], scientific = F, digits = input$num_of_digits_tau),
          MaxTissue = max_tissues,
          AlleleSpecific = ifelse(rownames(gene_expression_data_tissue) %in% allele_specific_genes, "Yes", "No"),
          stringsAsFactors = FALSE,
          row.names = NULL
        )
        
        write.csv(result_df, file, row.names = FALSE)
      }
    }
  )
  
  output$summary_stats_table <- renderTable({
    gene_expression_data_tissue <- gene_expression_data_tissue_reactive()
    if (is.null(gene_expression_data_tissue)) {
      return(NULL)
    }
    
    tau_scores <- calcTau_custom(gene_expression_data_tissue)
    
    summary_stats <- data.frame(
      Metric = c("Number of Genes", "Average Tau", "Median Tau", "Tau Variance", "Min Tau", "Max Tau"),
      Value = c(
        nrow(gene_expression_data_tissue),
        round(mean(tau_scores$tau, na.rm = TRUE), input$num_of_digits_tau),
        round(median(tau_scores$tau, na.rm = TRUE), input$num_of_digits_tau),
        round(var(tau_scores$tau, na.rm = TRUE), input$num_of_digits_tau),
        round(min(tau_scores$tau, na.rm = TRUE), input$num_of_digits_tau),
        round(max(tau_scores$tau, na.rm = TRUE), input$num_of_digits_tau)
      )
    )
    
    transposed_summary_stats <- t(summary_stats)
    
    colnames(transposed_summary_stats) <- NULL
    
    return(transposed_summary_stats)
  }, rownames = FALSE, colnames = FALSE)
  
  output$tau_histogram <- renderPlot({
    gene_expression_data_tissue <- gene_expression_data_tissue_reactive()
    if (is.null(gene_expression_data_tissue)) {
      return(NULL)
    }
    
    tau_scores <- calcTau_custom(gene_expression_data_tissue)
    tau_values <- tau_scores$tau
    
    hist(
      tau_values,
      breaks = 30,
      col = "darkgreen",
      border = "white",
      main = "Distribution of Tau Values",
      xlab = "Tau",
      ylab = "Frequency"
    )
  })
  
  output$allele_specific_bar_chart <- renderPlot({
    allele_specific_genes <- allele_specific_genes_reactive()
    gene_expression_data_tissue <- gene_expression_data_tissue_reactive()
    if (is.null(gene_expression_data_tissue)) {
      return(NULL)
    }
    
    allele_specific_count <- sum(toupper(rownames(gene_expression_data_tissue)) %in% toupper(allele_specific_genes))
    non_allele_specific_count <- nrow(gene_expression_data_tissue) - allele_specific_count
    
    bar_data <- data.frame(
      Category = c("Allele-Specific", "Non-Allele-Specific"),
      Count = c(allele_specific_count, non_allele_specific_count)
    )
    
    barplot(bar_data$Count, names.arg = bar_data$Category, col = c("blue", "gray"),
            main = "Number of Allele-Specific Genes", ylab = "Count")
  })
  
  output$top_max_tissues_table <- renderPlot({
    gene_expression_data_tissue <- gene_expression_data_tissue_reactive()
    if (is.null(gene_expression_data_tissue)) {
      return(NULL)
    }
    
    max_tissues <- apply(gene_expression_data_tissue, 1, function(x) {
      colnames(gene_expression_data_tissue)[which.max(x)]
    })
    
    tissue_frequency <- table(max_tissues)
    
    top_tissues <- sort(tissue_frequency, decreasing = TRUE)
    top_tissues_df <- data.frame(Tissue = names(top_tissues), Frequency = as.vector(top_tissues))
    
    par(mar = c(10, 4, 4, 2)) 
    
    bar_pos <- barplot(
      top_tissues_df$Frequency[1:10],
      names.arg = NA,  
      col = "darkorange",
      main = "Top 10 Frequent Tissues",
      ylab = "Frequency",
      las = 1,  
      cex.names = 0.8  
    )
    
  
    text(x = bar_pos, 
         y = par("usr")[3] - 0.1 * max(top_tissues_df$Frequency),  
         labels = top_tissues_df$Tissue[1:10], 
         srt = 45, 
         adj = 1, 
         xpd = TRUE,  
         cex = 0.8)  
  })
  
  observeEvent(input$add_genes_button, {
    new_genes <- unlist(strsplit(input$add_genes_input, ",\\s*"))
    
    new_genes <- new_genes[new_genes != ""]
    
    if (length(new_genes) > 0) {
      add_to_cart(toupper(new_genes))
      
      updateTextInput(session, "add_genes_input", value = "")
      
      showNotification("Genes have been added.", type = "message")
    } else {
      showNotification("Please enter valid gene names.", type = "warning")
    }
  })
  
  
}

shinyApp(ui = ui, server = server)
