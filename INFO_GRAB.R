
###
# INFO GRAB App
# Author: Oliver Brown
# Email: ombrown@uw.edu
###

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
  "BiocManager", "VariantAnnotation", "bslib", "ggdist", "ggiraph"
)

install_if_missing <- function(pkg) {
  if (!require(pkg, character.only = TRUE)) {
    install.packages(pkg, dependencies = TRUE)
  }
}

lapply(required_packages, install_if_missing)

if (!requireNamespace("BiocManager", quietly = TRUE)) {
  install.packages("BiocManager")
  BiocManager::install(c("Rhtslib", "genefilter", "rtracklayer", "DESeq2", 
                         "VariantAnnotation", "Rsamtools", 
                         "GenomicAlignments", "BSgenome", 
                         "preprocessCore", "GenomicFeatures", "edgeR"), force = T)
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
library(ggdist)
library(ggiraph)


# load required package
library(stringr)

# 1) Load data (if not already in memory)
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

# 2) Fix sample names for macrophages
phenodata$Sample_name2 <- ifelse(
  grepl("macrophage", phenodata$Tissue, ignore.case = TRUE) &
    !grepl("^macrophage_", phenodata$Sample_name2, ignore.case = TRUE),
  paste0("macrophage_", phenodata$Sample_name2),
  phenodata$Sample_name2
)

orig_names <- phenodata$Sample_name
new_names  <- orig_names
mask <- grepl("macrophage", phenodata$Tissue, ignore.case = TRUE) &
  !grepl("^Sample_macrophage_", new_names, ignore.case = TRUE)
if (any(mask)) {
  new_names[mask] <- sub("^(?i)Sample[_ ]+", "Sample_macrophage_", new_names[mask], perl = TRUE)
}
phenodata$Sample_name <- new_names
name_mapping <- setNames(new_names, orig_names)

update_sample_names <- function(data_obj) {
  cols <- colnames(data_obj)
  colnames(data_obj) <- sapply(cols, function(x) if (x %in% names(name_mapping)) name_mapping[[x]] else x)
  data_obj
}

RPKM_data <- update_sample_names(RPKM_data)
cnt_data  <- update_sample_names(cnt_data)

# 3) Define your color palettes
mycolors1 <- c(
  `Kidney`="#43009A",
  `Trachea`="#990099",
  `B cells`="#0DFAFA",
  `T cell (Spleen)`="#13B7B7",
  `Bursa`="#004C99",
  `Thymus`="#2685E4",
  `Macrophage at Day 0 Differentiation (D0)`="#F02B6D",
  `Macrophage at Day 3 Differentiation (D3)`="#FF0091",
  `Macrophage at Day 6 Differentiation (D6)`="#F572BC",
  `Macrophage (From Lung)`="#D06AAA",
  `Monocyte (Blood)`="#91155B",
  `Ileum`="#98D55A",
  `Jejunum`="#4C9900",
  `Proximal Cecum`="#CCFF99",
  `Iliotibialis Lateralis (M.Tight)`="#A1122A",
  `Pectoralis Major (M.Breast)`="#FFCCCC",
  `Isthmus`="#FF8000",
  `Magnum`="#DE5100",
  `Shell Gland`="#FFC78E",
  `Ovary`="#CCCC00"
)
mycolors2 <- c(
  `Immune System (Tissues and Cells)` = "#31E1F7",
  `Respiratory System`               = "#400D51",
  `Excretory System`                 = "#FF7777",
  `Muscular System (Chicken Meat)`   = "#D800A6",
  `Intestinal System`                = "#6499E9",
  `Reproductive Female Duct System`  = "#836FFF"
)

# 4) Define the typo‐corrections
repl <- c(
  "Iliotibialis Lateralis \\(M\\.Tight\\)" = "Iliotibialis Lateralis (Thigh)",
  "Pectoralis Major \\(M\\.Breast\\)"      = "Pectoralis Major (Breast)"
)

# 5) Apply corrections everywhere
phenodata$Tissue    <- str_replace_all(phenodata$Tissue,    repl)
names(mycolors1)    <- str_replace_all(names(mycolors1),    repl)
names(mycolors2)    <- str_replace_all(names(mycolors2),    repl)
names(comparison_results) <- str_replace_all(names(comparison_results), repl)

# 6) Recompute any derived tissue lists
tissue_names     <- unique(phenodata$Tissue)
systems          <- unique(phenodata$System)
tissue_by_system <- setNames(
  lapply(systems, function(s) unique(phenodata$Tissue[phenodata$System == s])),
  systems
)


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
      @keyframes blink {
      0% { background-color: yellow; }
      50% { background-color: orange; }
      100% { background-color: yellow; }
      }
      .blink {
        animation: blink 1s infinite;
      }
    .navbar {
      position: fixed;
      top: 0;
      width: 100%;
      z-index: 1000;
    }
    /* Add top padding to the body so content isn't hidden behind the navbar */
    body {
      padding-top: 70px;  /* Adjust this value based on your navbar's height */
    }
    
    /* Wrap plot outputs with extra space at bottom */
    .plot-container {
      padding-bottom: 350px;  /* Adjust this value as needed */
    }
    "))
  ),
  tabPanel("INFO GRAB", value = "de_tab",
           fluidPage(
             div(style = "text-align: center;",
                 h2("INFO GRAB"),
                 h3("INteractive Functional Ontology for Gene Regulation Analysis and Browsing"),
                 tags$hr(style = "height:1px; border:none; color:#300; background-color:#300;"),
                 h5(tags$b("Welcome to INFO GRAB."), " This Shiny Application allows you to perform differential expression analysis and 
                   data visualization of the transcriptome, looking at tissues in the chicken (Gallus gallus)."),
                 tags$hr(style = "height:1px; border:none; color:#300; background-color:#300;"),
                 img(src = "chicken.png", id = "chicken-image", style = "width: 40%; max-width: 100%; height: auto;")
             ),
             uiOutput("tissue_display")
           )
  ),
  
  # Add labels
  
  tabPanel("Differential Analysis", value = "de_tab",
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
                          conditionalPanel(
                            condition = "input.correlation_mode == 'system'",
                            actionButton("confirm_cor_plot", "Confirm Plot")
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
                          div(class = "plot-container",
                              withLoader(plotOutput("correlation_plot", width = "90%", height = "750px"), type = "html", loader = "dnaspin")
                          )
                        )
                      )
             ),
             
             
             tabPanel("Volcano Plot",
                      sidebarLayout(
                        sidebarPanel(
                          textInput("searched_gene_volcano", "Label plot: (comma-separated) e.g. EGR1, CTCF", width = '600px', value = "", placeholder = "Search for genes here"),
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
                          div(class = "plot-container",
                              withLoader(plotOutput("volcano_plot", width = "90%", height = "750px"), 
                                         type = "html", loader = "dnaspin")
                          )
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
                          tags$p(
                            "The color scale is based on ", tags$b("z-scores"),
                            ", which represent the number of standard deviations away from the mean expression level. High positive ", tags$b("z-scores"),
                            " indicate higher expression levels, while high negative ", tags$b("z-scores"),
                            " indicate lower expression levels."
                          )
                        ),
                        mainPanel(
                          uiOutput("title_heatmap"),
                          withLoader(plotOutput("heatmap_plot", width = "85%", height = "750px",
                                                brush = brushOpts(id = "heatmap_brush", resetOnNew = TRUE)),
                                     type = "html", loader = "dnaspin"),
                          tags$hr(),
                          tags$p(tags$b("Brushing: hover your mouse over the plot, click, then drag to select genes."))
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
                          div(class = "plot-container",
                              withLoader(plotOutput("pca_plot", height = "600px"), type = "html", loader = "dnaspin")
                          )
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
  
  tabPanel("Expression Profiles", value = "tissue_tab",
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
  
  tabPanel("Tissue-Specific Analysis", value = "tissue_tab",
           
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
                          tags$p("This section displays a tissue-specific heatmap, which focuses on the expression of genes that are most tissue specific across the selected tissues."),
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
  
  shiny::tabPanel(
    title = "Allelic Expression",
    value = "ase_tab",
    shiny::sidebarLayout(
      
      # ─ Sidebar with controls ────────────────────────────────────────────────────
      shiny::sidebarPanel(
        # Data type selector
        shiny::radioButtons(
          inputId  = "raincloud_data_type",
          label    = "Data Type:",
          choices  = c("SNP ASE" = "snp", "Gene ASE" = "gene"),
          selected = "snp"
        ),
        
        # How many top features?
        shiny::numericInput(
          inputId = "top_n_raincloud",
          label   = "Number of Top Features (Raincloud):",
          value   = 10, min = 1, max = 10000, step = 1
        ),
        
        # FDR cutoff
        shiny::sliderInput(
          inputId = "pval_cutoff",
          label   = "P-value Cutoff (FDR):",
          min     = 0, max = 0.05,
          value   = 0.05, step = 0.005
        ),
        
        # Tissue selection
        shiny::selectInput(
          inputId  = "tissue_select_raincloud",
          label    = "Select Tissues:",
          choices  = NULL,
          selected = NULL,
          multiple = TRUE
        ),
        
        # Search box only when "All Tissues" **and** gene mode
        shiny::conditionalPanel(
          condition = "input.raincloud_data_type=='gene' && input.tissue_select_raincloud.includes('All Tissues')",
          shiny::textInput(
            inputId = "raincloud_searched_gene",
            label   = "Search for Genes (comma-separated):",
            value   = ""
          )
        ),
        
        # Search box only when "All Tissues" **and** SNP mode
        shiny::conditionalPanel(
          condition = "input.raincloud_data_type=='snp' && input.tissue_select_raincloud.includes('All Tissues')",
          shiny::textInput(
            inputId = "raincloud_searched_snp",
            label   = "Search for SNP IDs (comma-separated):",
            value   = ""
          )
        ),
        
        # Random / clear buttons
        shiny::actionButton("random_genes_raincloud", "10 Random Genes"),
        shiny::actionButton("random_snps_raincloud",  "10 Random SNPs"),
        shiny::actionButton("clear_search_raincloud", "Clear Search"),
        
        # Plot toggles
        shiny::checkboxInput("show_points",    "Show Individual Points", TRUE),
        shiny::checkboxInput("show_mean_line", "Show Mean Line",        TRUE),
        
        # Download + grab buttons
        shiny::downloadButton("download_raincloud_plot", "Download PNG"),
        
        shiny::conditionalPanel(
          condition = "input.raincloud_data_type=='gene'",
          shiny::actionButton("grab_selected_raincloud", "Grab Selected"),
          shiny::actionButton("grab_displayed_raincloud", "Grab All Displayed")
        ),
        
        # Info panel
        ## in your UI, sidebarPanel, instead of the single info div:
        
        # Gene info (only when “Gene ASE” is selected)
        shiny::conditionalPanel(
          condition = "input.raincloud_data_type == 'gene'",
          tags$div(
            style = "margin:15px 0; padding:10px; background:#f9f9f9; border:1px solid #ddd;",
            tags$strong("Mean effect:"), "the average of sample allelic imbalance estimates, 
            reflecting the overall magnitude of allele-specific expression in a gene."
          )
        ),
        
        # SNP info (exactly your original AE text)
        shiny::conditionalPanel(
          condition = "input.raincloud_data_type == 'snp'",
          tags$div(
            style = "margin:15px 0; padding:10px; background:#f9f9f9; border:1px solid #ddd;",
            tags$strong("Allelic Expression (AE):"), 
            "measures how far the reference-allele fraction (Rratio) is from 0.5:",
            tags$br(),
            HTML("&nbsp;&nbsp;<em>AE</em> = |0.5 – R<sub>ratio</sub>|, where R<sub>ratio</sub> = Ref reads / Total reads"),
            tags$br(), tags$br(),
            "Values of AE near 0 indicate nearly equal expression from both alleles.  ",
            "Values of AE approaching 0.5 indicate strong bias toward one allele.",
            tags$br(), tags$strong("Common threshold:"), 
            "AE ≥ 0.3 is often used to call a gene “allele-specific,” but you may tighten or loosen this cut-off."
          )
        )
      ),
      
      # ─ Main panel with the raincloud plot ───────────────────────────────────────
      shiny::mainPanel(
        shiny::h3(shiny::textOutput("raincloud_title")),
        withLoader(
          ggiraph::girafeOutput("raincloud_plot", width = "100%", height = "800px"),
          type   = "html",
          loader = "dnaspin"
        )
      )
    )
  ),
  
  # Genome Browser
  
  tabPanel(
    "Genome Viewer",
    value = "ase_tab",
    fluidPage(
      tags$head(
        tags$style(
          HTML("
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
          ")
        )
      ),
      
      fluidRow(
        column(
          width = 4,
          textInput("genome_browser_search", "Search For a Locus or a Gene"),
          actionButton("search_button", "Search", class = "btn-primary")
        ),
        column(
          width = 2,
          radioButtons(
            inputId = "input_type",
            label   = "Add a Track (GFF3, BAM, BED, or VCF):",
            choices = c("Local File" = "file", "URL" = "url"),
            selected = "file"
          )
        ),
        column(
          width = 4,
          conditionalPanel(
            condition = "input.input_type == 'file'",
            fileInput("file", "Track File", width = '100%')
          ),
          conditionalPanel(
            condition = "input.input_type == 'url'",
            textInput("url", "Track URL"),
            textInput("index", "Index URL")
          ),
          actionButton("loadTrack",          "Load Custom Track", class = "btn-primary"),
          actionButton("load_gff",           "Load Annotation Track"),
          actionButton("load_geneASE_track", "Load Gene ASE Track", class = "action-button"),
          actionButton("load_snp",           "Load SNP ASE Track", class = "action-button"),
          actionButton("grab_all_ase_genes", "Grab All ASE Genes", class = "btn-primary")
        )
      ),
      
      br(),
      
      ## ROW 2a: Tissues checkbox for ASE, including "All Tissues"
      fluidRow(
        column(
          width = 4,
          checkboxGroupInput(
            inputId = "ase_tissue_filter",
            label   = "Filter ASE by Tissue (select one or more):",
            choices = character(0),  # Populated by server
            selected = NULL
          )
        )
      ),
      
      br(),
      
      ## ROW 2b: Single vs Multiple track radio button
      fluidRow(
        column(
          width = 4,
          radioButtons(
            inputId = "track_mode",
            label   = "Track Mode:",
            choices = c(
              "Single Track (Combine All Selected Tissues)" = "single",
              "Multiple Tracks (One per Tissue)"            = "multiple"
            ),
            selected = "multiple"
          )
        )
      ),
      
      br(),
      
      ## ROW 3: IGV Output
      fluidRow(
        column(
          width = 12,
          igvShinyOutput("igvShiny", height = "1000px")
        )
      )
    )
  ),
  
  # In your UI definition (within the "Genes" tabPanel), 
  # find the section where you have the current layout for the genes cart.
  # Below is an example of how to integrate the mode selection input.
  
  tabPanel("Genes", value = "genes_tab",
           sidebarLayout(
             sidebarPanel(
               tags$textarea(id = "gene_text_area", style = "position: absolute; left: -9999px;"),
               actionButton("copy_cart_genes", "Copy Genes to Clipboard"),
               selectInput("download_type", "Download Data:", 
                           choices = c("DE CSV" = "de",
                                       "Tissue CSV" = "tissue",
                                       "ASE CSV" = "ase"),
                           selected = "combined"),
               downloadButton("download_results", "Download Selected CSV"),
               
               actionButton("delete_selected_genes", "Delete Selected Genes"),
               actionButton("clear_cart", "Clear Genes"),
               tags$hr(style = "height:1px; border:none; color:#300; background-color:#300;"),
               textInput("add_genes_input", "Add Genes (comma separated):", value = "", placeholder = "Search for genes here"),
               actionButton("add_genes_button", "Add Genes"),
               tags$hr(style = "height:1px; border:none; color:#300; background-color:#300;"),
               
               # NEW: Mode selection input
               selectInput("gene_mode", "Select Mode:",
                           choices = c("Differential Expression" = "de",
                                       "Tissue-Specific Expression" = "tissue",
                                       "ASE Expression" = "ase"),
                           selected = "tissue"), 
               
               tags$h3("Summary Statistics"),
               tableOutput("summary_stats_table"),
               
               tags$h3("Tau Value Distribution"),
               plotOutput("tau_histogram", height = "350px"),
               
               tags$h3("Allele-Specific Genes Count"),
               plotOutput("allele_specific_bar_chart", height = "350px"),
               
               tags$hr(style = "height:1px; border:none; color:#300; background-color:#300;"),
               
               h3("Mode-Specific Visualization"),
               plotOutput("mode_specific_plot", height = "300px"), # Adjust height as needed
               
               tags$h3("Useful Links"),
               div(style = "font-size: 15px;",
                   tags$ul(
                     tags$li(a(href = "http://bioinformatics.sdstate.edu/go/", 
                               "ShinyGO: Gene Ontology Enrichment Analysis", target = "_blank")),
                     tags$li(a(href = "https://davidbioinformatics.nih.gov/tools.jsp", 
                               "DAVID: Database for Annotation, Visualization and Integrated Discovery", target = "_blank")),
                     tags$li(a(href = "https://pantherdb.org", 
                               "PantherDB: Panther Classification System", target = "_blank"))
                   )
               )
             ),
             mainPanel(
               withLoader(DTOutput("cart_gene_table"), loader = "dnaspin")
             )
           )
  ),
  
  tabPanel(
    "Help",
    fluidPage(
      h2("App Version and Information"),
      tabsetPanel(
        
        # 1) Differential Expression (DE) Analysis
        tabPanel("Differential Expression (DE) Analysis Help",
                 mainPanel(
                   h2("Differential Expression Analysis"),
                   h3("Overview"),
                   p("INFO GRAB’s DE module utilizes RNA-seq counts and phenotype metadata (via DESeq2) to pinpoint genes that are significantly differentially expressed between two tissues."),
                   
                   h3("Features and Parameters"),
                   tags$ul(
                     tags$li(tags$b("Tissue Selection:"), 
                             " Pick any two tissues to calculate log2 fold changes."),
                     tags$li(tags$b("FDR Control:"), 
                             " Apply an adjusted p-value threshold to control false discoveries."),
                     tags$li(tags$b("Fold Change (FC) Threshold:"), 
                             " Focus on genes exceeding a chosen log2 fold-change limit."),
                     tags$li(tags$b("Digits Precision:"), 
                             " Select the number of decimal places for numeric outputs."),
                     tags$li(tags$b("Search & Random Genes:"), 
                             " Locate specific genes or display a random set for easy exploration."),
                     tags$li(tags$b("Download CSV:"), 
                             " Export the current DE results for offline analysis."),
                     tags$li(tags$b("Fold Change Equation:"), 
                             " View the formula for log2 fold change between the selected tissues.")
                   ),
                   
                   h3("Visual Summaries"),
                   p("Interactive tables allow sorting and filtering by fold change, FDR, and p-values to streamline data inspection."),
                   
                   h3("Technical Process"),
                   p("The module uses DESeq2 for normalization and statistical testing, ensuring compatibility with FAANG data standards.")
                 )
        ),
        
        # 2) Correlation Plot
        tabPanel("Correlation Plot Help",
                 mainPanel(
                   h2("Correlation Plot"),
                   h3("Overview"),
                   p("This plot visualizes pairwise correlation of gene expression among tissues or grouped systems, providing insight into expression similarities and differences."),
                   
                   h3("Features and Parameters"),
                   tags$ul(
                     tags$li(tags$b("Select Mode:"), 
                             " Either compare specific tissues or group multiple tissues by their biological system."),
                     tags$li(tags$b("Use Filtered Data:"), 
                             " Optionally apply FDR and FC thresholds; otherwise, the top 1,000 most tissue specific genes are used."),
                     tags$li(tags$b("System Selection:"), 
                             " Narrow down to tissues within a single system, e.g., immune or respiratory."),
                     tags$li(tags$b("Legend:"), 
                             " Each tissue is assigned a distinct color for easy identification.")
                   ),
                   
                   h3("Plot Details"),
                   p("Diagonal panels show density plots of single-tissue expression distributions, while upper panels display correlation coefficients between tissue pairs."),
                   
                   h3("Download Option"),
                   p("A PNG of the correlation matrix can be saved for external reference or presentation.")
                 )
        ),
        
        # 3) Volcano Plot
        tabPanel("Volcano Plot Help",
                 mainPanel(
                   h2("Volcano Plot"),
                   h3("Overview"),
                   p("Plots each gene’s log2 fold change against –log10(FDR), emphasizing highly significant and strongly regulated genes."),
                   
                   h3("Features and Customization"),
                   tags$ul(
                     tags$li(tags$b("Search Gene:"), 
                             " Label specific genes by name in the plot."),
                     tags$li(tags$b("Basic or Advanced Labeling:"), 
                             " Control how many genes are labeled on each side of the plot."),
                     tags$li(tags$b("Background Color:"), 
                             " Adjust the color scheme for clarity."),
                     tags$li(tags$b("Download PNG:"), 
                             " Save the volcano plot as a PNG for easy sharing or publication.")
                   ),
                   
                   h3("Advanced Labeling"),
                   p("Advanced labeling separates the labeling on left vs. right to highlight top genes by fold change or FDR.")
                 )
        ),
        
        # 4) Heatmap
        tabPanel("Heatmap Help",
                 mainPanel(
                   h2("Heatmap"),
                   h3("Overview"),
                   p("Displays z-score transformed gene expression across selected tissues, unveiling patterns in expression intensities and co-expression."),
                   
                   h3("Parameters and Customization"),
                   tags$ul(
                     tags$li(tags$b("Sort by FDR or Fold Change:"), 
                             " Emphasize either statistical significance or expression magnitude."),
                     tags$li(tags$b("Top Genes:"), 
                             " Specify how many top-ranked genes to include."),
                     tags$li(tags$b("Search & Random Genes:"), 
                             " Pinpoint certain genes or explore a random sample."),
                     tags$li(tags$b("Color Scale:"), 
                             " Customize the gradient palette to best represent expression levels.")
                   ),
                   
                   h3("How to Use Brushing"),
                   helpText(
                     "“Brushing” is how you interactively draw a box around points (in the volcano or correlation plots) or rows (in the heatmaps) to select them:",
                     tags$ol(
                       tags$li("Click and hold your left mouse button anywhere in the plot or heatmap."),
                       tags$li("Drag diagonally to form a rectangle over the points or rows you’re interested in."),
                       tags$li("Release the mouse button: anything inside that rectangle becomes your “selection.”"),
                       tags$li("A pop‑up modal will appear showing the genes you brushed. From there you can “Grab Genes” to add them to your cart.")
                     ),
                   )
                 )
        ),
        
        # 5) PCA
        tabPanel("PCA Help",
                 mainPanel(
                   h2("Principal Component Analysis (PCA)"),
                   h3("Overview"),
                   p("Reduces dimensionality by projecting gene expression data onto the principal components, revealing major sources of variance across tissues."),
                   
                   h3("PCA Modes and Options"),
                   tags$ul(
                     tags$li(tags$b("Selected Tissues, System, or All:"), 
                             " Restrict PCA to particular tissues, a biological system, or include every tissue."),
                     tags$li(tags$b("System Selection:"), 
                             " When a single system is chosen, focus solely on that system’s tissues."),
                     tags$li(tags$b("Top Tissue Specific Genes:"), 
                             " Up to 1,000 of the most tissue specific genes are used for clarity."),
                     tags$li(tags$b("Download Plot:"), 
                             " Save the PCA visualization as a PNG file.")
                   ),
                   
                   h3("Visualization Details"),
                   p("Each sample is colored by tissue and shaped by system, helping to identify clusters and outliers efficiently.")
                 )
        ),
        
        # 6) Raincloud Plot
        tabPanel("Raincloud Plot Help",
                 mainPanel(
                   h2("Raincloud Plot"),
                   h3("Overview"),
                   p("Combines a violin or half-violin display with box plots and individual data points in a single figure. INFO GRAB often uses raincloud plots to present allele-specific expression (ASE) distributions."),
                   
                   h3("Key Features"),
                   tags$ul(
                     tags$li(tags$b("Facet by Tissue:"), 
                             " Compare multiple tissues side by side or view them combined."),
                     tags$li(tags$b("Random Gene Selection:"), 
                             " Load random genes quickly to observe distribution patterns."),
                     tags$li(tags$b("Mean Line & Points:"), 
                             " Optionally highlight the mean or show jittered data points."),
                     tags$li(tags$b("Downloadable Plot:"), 
                             " Export the raincloud figure as a PNG.")
                   ),
                   
                   h3("Use Cases"),
                   p("This format is well-suited for ASE or other metrics, illustrating both the overarching distribution and each individual value in one plot.")
                 )
        ),
        
        # 7) Tissue-Specific Analysis
        tabPanel("Tissue-Specific Analysis Help",
                 mainPanel(
                   h2("Tissue-Specific Expression Analysis"),
                   h3("Overview"),
                   p("Identifies genes with tissue-focused expression using a tau threshold (e.g., tau ≥ 0.85) and highlights their expression across chosen systems or tissues."),
                   
                   h3("Features"),
                   tags$ul(
                     tags$li(tags$b("System & Tissue Selection:"), 
                             " Target the system(s) or tissues you wish to explore."),
                     tags$li(tags$b("Tau Threshold:"), 
                             " Set at 0.85 by default to capture strongly tissue-specific genes."),
                     tags$li(tags$b("Top Tissue Specific Genes:"), 
                             " Narrow the heatmap to genes with the greatest variance.")
                   ),
                   
                   h3("How to Use Brushing"),
                   helpText(
                     "“Brushing” is how you interactively draw a box around points (in the volcano or correlation plots) or rows (in the heatmaps) to select them:",
                     tags$ol(
                       tags$li("Click and hold your left mouse button anywhere in the plot or heatmap."),
                       tags$li("Drag diagonally to form a rectangle over the points or rows you’re interested in."),
                       tags$li("Release the mouse button: anything inside that rectangle becomes your “selection.”"),
                       tags$li("A pop‑up modal will appear showing the genes you brushed. From there you can “Grab Genes” to add them to your cart.")
                     ),
                   )
                 )
        ),
        
        # 8) Genome Viewer
        tabPanel("Genome Viewer Help",
                 mainPanel(
                   h2("Genome Browser"),
                   h3("Overview"),
                   p("Allows visualization of GFF3, BAM, BED, and VCF files in the chicken genome. Preloaded annotation and ASE tracks expedite initial inspection of genes and variants."),
                   
                   h3("Supported File Types"),
                   tags$ul(
                     tags$li(tags$b("GFF3:"), " Stores gene feature annotations."),
                     tags$li(tags$b("BAM:"), " Displays read alignments in genomic context."),
                     tags$li(tags$b("BED:"), " Highlights particular regions of interest."),
                     tags$li(tags$b("VCF:"), " Contains variant calls.")
                   ),
                   
                   h3("Preloaded Tracks"),
                   p("SNP ASE tracks and a basic gene annotation track are included, giving immediate insight into regulatory variants and transcription units."),
                   
                   h3("Common Issues"),
                   tags$ul(
                     tags$li(tags$b("File Compatibility:"), 
                             " Check that file formats match one of the accepted types."),
                     tags$li(tags$b("Chromosome Mapping:"), 
                             " Confirm your file references the same genome assembly naming system.")
                   )
                 )
        ),
        
        # 9) Genes
        tabPanel("Genes Help",
                 mainPanel(
                   h2("Genes Management and Visualization"),
                   h3("Overview"),
                   p("Facilitates adding, organizing, and examining gene lists for different analyses (Differential Expression, Tissue-Specific Expression, ASE). This section also displays summaries such as tau distributions and ASE frequencies."),
                   
                   h3("Key Features"),
                   tags$ul(
                     tags$li(tags$b("Add Genes:"), 
                             " Enter or search gene names to include them in your working list."),
                     tags$li(tags$b("Clear & Copy:"), 
                             " Remove genes from the list or copy them to your clipboard for use in external tools."),
                     tags$li(tags$b("Download Data:"), 
                             " Export your selected genes in DE, Tissue-Specific, ASE, or Combined CSV files."),
                     tags$li(tags$b("Switch Mode:"), 
                             " Choose how your gene list is viewed: show DE results, tissue-based stats, or ASE data."),
                     tags$li(tags$b("Summary Statistics:"), 
                             " Inspect total genes, average tau values, or ASE gene counts, depending on the active mode."),
                     tags$li(tags$b("Tau Value Distribution:"), 
                             " Visualize tissue-specificity with a histogram of tau values for your selected genes."),
                     tags$li(tags$b("Allele-Specific Genes Count:"), 
                             " Bar chart indicating how many listed genes exhibit ASE."),
                     tags$li(tags$b("Mode-Specific Plot:"), 
                             " Quickly generate mini plots (e.g., Volcano, Bar, Box) tailored to the current mode."),
                     tags$li(tags$b("Useful Links:"), 
                             " Access external tools such as ShinyGO, DAVID, or PantherDB for advanced gene enrichment analyses.")
                   ),
                   
                   h3("Genes Table"),
                   p("A dynamic table displays your selected genes, updating columns based on the chosen mode (DE, tissue, or ASE).")
                 )
        ),
        
        tabPanel("Version",
                 br(),
                 p("Version: 1.0.0"),
                 p("Developed by: Oliver Brown (ombrown@uw.edu), Andressa Oliveira de Lima, and R. David Hawkins")
        ),
        
        tabPanel("Credits",
                 h4("Acknowledgments"),
                 p(),
                 br(),
                 h4("Citations"),
                 tags$ul(
                   tags$li("de Lima, A. O., Ng, T. T., Sparling, B., Griggs, L. M., Lai, K., Drechsler, Y., & Hawkins, R. D. (2024). An updated Gallus gallus genome annotation through multi-tissue transcriptome analysis. Under Review."),
                   tags$li("Di Tommaso, P., Chatzou, M., Floden, E. W., Barja, P. P., Palumbo, E., & Notredame, C. (2017). Nextflow enables reproducible computational workflows. Nature Biotechnology, 35(4), 316–319. https://doi.org/10.1038/nbt.3820"),
                   tags$li("Hu, Z.-L., Park, C. A., & Reecy, J. M. (2019). Building a livestock genetic and genomic information knowledgebase through integrative developments of Animal QTLdb and CorrDB. Nucleic Acids Research, 47(D1), D701–D710. https://doi.org/10.1093/nar/gky1084"),
                   tags$li("Hu, Z.-L., Park, C. A., & Reecy, J. M. (2020). AgAnimalGenomes: Improving the accessibility of livestock genomic resources. Frontiers in Genetics, 11, 615. https://doi.org/10.3389/fgene.2020.00615"),
                   tags$li("Overbey, E.G., McCarthy, F., Mason, A.S., Fleming, D., Luna, A.M., Dickel, D., & Schmidt, C.J. (2021). Transcriptomes of an array of chicken ovary, intestinal, and immune cells and tissues. Frontiers in Genetics, 12, 664424. https://doi.org/10.3389/fgene.2021.664424"),
                   tags$li("Lüleci, H.B., & Yılmaz, A. (2022). Robust and rigorous identification of tissue-specific genes by statistically extending tau score. BioData Mining, 15(31). https://doi.org/10.1186/s13040-022-00315-9"),
                   tags$li("Ge, S.X., Jung, D., & Yao, R. (2020). ShinyGO: A graphical gene-set enrichment tool for animals and plants. Bioinformatics, 36(8), 2628–2629. https://doi.org/10.1093/bioinformatics/btz931"),
                   tags$li("Huang, D.W., Sherman, B.T., & Lempicki, R.A. (2009). Systematic and integrative analysis of large gene lists using DAVID bioinformatics resources. Nature Protocols, 4, 44–57. https://doi.org/10.1038/nprot.2008.211"),
                   tags$li("Mi, H., Muruganujan, A., Ebert, D., Huang, X., & Thomas, P.D. (2023). PANTHER version 17: expanded annotation data from Gene Ontology and Reactome pathways, and data analysis tool enhancements. Nucleic Acids Research, 51(D1), D541–D549. https://doi.org/10.1093/nar/gkad103"),
                   tags$li("Love, M. I., Huber, W., & Anders, S. (2014). Moderated estimation of fold change and dispersion for RNA-seq data with DESeq2. Genome Biology, 15(12), 550."),
                   tags$li("BioRender.com. BioRender software was used to create figures. BioRender, 2024. Available at https://biorender.com."),
                   tags$li("Morgan, M., & Pagès, H. (2024). Rhtslib: HTSlib high-throughput sequencing library as an R package (Version 3.2.0) [R package]. Bioconductor. Available at https://bioconductor.org/packages/Rhtslib."),
                   tags$li("Gentleman, R., Carey, V., Bates, D., Bolstad, B., Dettling, M., Dudoit, S., ... & Zhang, J. (2024). genefilter: Methods for filtering genes from high-throughput experiments (Version 1.80.0) [R package]. Bioconductor. https://bioconductor.org/packages/genefilter."),
                   tags$li("Lawrence, M., Gentleman, R., & Carey, V. (2024). rtracklayer: R interface to genome annotation files (Version 1.60.0) [R package]. Bioconductor. Available at https://bioconductor.org/packages/rtracklayer."),
                   tags$li("Love, M. I., Huber, W., & Anders, S. (2014). Moderated estimation of fold change and dispersion for RNA-seq data with DESeq2. Genome Biology, 15(12), 550."),
                   tags$li("Obenchain, V., Lawrence, M., & Carey, V. (2024). VariantAnnotation: Annotation of Genetic Variants (Version 1.48.0) [R package]. Bioconductor. Available at https://bioconductor.org/packages/VariantAnnotation."),
                   tags$li("Morgan, M., & Pagès, H. (2024). Rsamtools: Binary alignment (BAM), FASTA, variant call (BCF), and tabix file import (Version 2.16.0) [R package]. Bioconductor. https://bioconductor.org/packages/Rsamtools."),
                   tags$li("Pagès, H., Obenchain, V., & Morgan, M. (2024). GenomicAlignments: Representation and manipulation of short genomic alignments (Version 1.42.0) [R package]. Bioconductor. Available at https://bioconductor.org/packages/GenomicAlignments."),
                   tags$li("Pagès, H. (2024). BSgenome: Software infrastructure for efficient representation of full genomes and their SNPs (Version 1.68.0) [R package]. Bioconductor. Available at https://bioconductor.org/packages/BSgenome."),
                   tags$li("Bolstad, B. M. (2024). preprocessCore: A collection of pre-processing functions (Version 1.62.0) [R package]. Bioconductor. Available at https://bioconductor.org/packages/preprocessCore."),
                   tags$li("Carlson, M., & Pagès, H. (2024). GenomicFeatures: Tools for making and manipulating transcript-centric annotations (Version 1.58.0) [R package]. Bioconductor. https://bioconductor.org/packages/GenomicFeatures."),
                   tags$li("Robinson, M. D., McCarthy, D. J., & Smyth, G. K. (2010). edgeR: A Bioconductor package for differential expression analysis of digital gene expression data. Bioinformatics, 26(1), 139–140."),
                   tags$li("Marini, F., & Binder, H. (2019). pcaExplorer: An R/Bioconductor package for interacting with RNA-seq principal components. BMC Bioinformatics, 20, 331."),
                   tags$li("Shannon, P. (2024). igvShiny: Integrating Interactive Genomics Viewer (IGV) with Shiny applications (Version 1.0.0) [R package]. GitHub. Available at https://github.com/paul-shannon/igvShiny."),
                   tags$li("Chang, W., Cheng, J., Allaire, J., Sievert, C., & McPherson, J. (2024). shiny: Web Application Framework for R (Version 1.7.4) [R package]. RStudio. Available at https://shiny.rstudio.com."),
                   tags$li("Bailey, E. (2024). shinyBS: Twitter Bootstrap Components for Shiny (Version 0.61) [R package]. https://CRAN.R-project.org/package=shinyBS."),
                   tags$li("Chang, W. (2024). shinythemes: Themes for Shiny (Version 1.2.0) [R package]. RStudio. https://CRAN.R-project.org/package=shinythemes."),
                   tags$li("Attali, D. (2024). shinyjs: Easily Improve the User Experience of Your Shiny Apps in Seconds (Version 2.1.0) [R package]. https://CRAN.R-project.org/package=shinyjs."),
                   tags$li("Chang, W., & Borges, B. (2024). shinydashboard: Create Dashboards with 'Shiny' (Version 0.7.2) [R package]. https://CRAN.R-project.org/package=shinydashboard."),
                   tags$li("Sali, A. (2024). shinycssloaders: Add Loading Animations to a 'shiny' Output While It's Recalculating (Version 1.0.0) [R package]. https://CRAN.R-project.org/package=shinycssloaders."),
                   tags$li("Sali, A. (2024). shinycustomloader: Add Custom Loading Animations to 'shiny' Outputs (Version 0.9.0) [R package]. https://CRAN.R-project.org/package=shinycustomloader."),
                   tags$li("Wickham, H., François, R., Henry, L., Müller, K., & Vaughan, D. (2023). dplyr: A Grammar of Data Manipulation (Version 1.1.4) [R package]. https://dplyr.tidyverse.org."),
                   tags$li("Wickham, H. (2016). ggplot2: Elegant Graphics for Data Analysis. Springer-Verlag New York. https://ggplot2.tidyverse.org."),
                   tags$li("Slowikowski, K. (2024). ggrepel: Automatically Position Non-Overlapping Text Labels with 'ggplot2' (Version 0.9.2) [R package]. https://CRAN.R-project.org/package=ggrepel."),
                   tags$li("Kolde, R. (2019). pheatmap: Pretty Heatmaps (Version 1.0.12) [R package]. https://CRAN.R-project.org/package=pheatmap."),
                   tags$li("Xie, Y., Cheng, J., & Tan, X. (2024). DT: A Wrapper of the JavaScript Library 'DataTables' (Version 0.28) [R package]. https://CRAN.R-project.org/package=DT."),
                   tags$li("Wickham, H. (2019). stringr: Simple, Consistent Wrappers for Common String Operations (Version 1.4.0) [R package]. https://CRAN.R-project.org/package=stringr."),
                   tags$li("Garnier, S. (2018). viridis: Default Color Maps from 'matplotlib' (Version 0.5.1) [R package]. https://CRAN.R-project.org/package=viridis."),
                   tags$li("Wickham, H., Hester, J., & Bryan, J. (2024). readr: Read Rectangular Text Data (Version 2.1.4) [R package]. https://CRAN.R-project.org/package=readr."),
                   tags$li("Wickham, H., & Girlich, M. (2024). tidyr: Tidy Messy Data (Version 1.2.1) [R package]. https://CRAN.R-project.org/package=tidyr."),
                   tags$li("Müller, K., & Wickham, H. (2024). tibble: Simple Data Frames (Version 3.1.8) [R package]. https://CRAN.R-project.org/package=tibble."),
                   tags$li("Henry, L., & Wickham, H. (2020). purrr: Functional Programming Tools (Version 0.3.4) [R package]. https://CRAN.R-project.org/package=purrr."),
                   tags$li("Bache, S. M., & Wickham, H. (2014). magrittr: A Forward-Pipe Operator for R (Version 1.5) [R package]. https://CRAN.R-project.org/package=magrittr."),
                   tags$li("Cheng, J., Sievert, C., & Xie, Y. (2024). htmltools: Tools for HTML (Version 0.5.4) [R package]. https://CRAN.R-project.org/package=htmltools."),
                   tags$li("Wickham, H., & Seidel, D. (2020). scales: Scale Functions for Visualization (Version 1.1.1) [R package]. https://CRAN.R-project.org/package=scales."),
                   tags$li("Yu, G. (2020). ggplotify: Convert Plot to 'grob' or 'ggplot' Object (Version 0.0.5) [R package]. https://CRAN.R-project.org/package=ggplotify."),
                   tags$li("Vaidyanathan, R., Allaire, J., & Cheng, J. (2024). htmlwidgets: HTML Widgets for R (Version 1.6.2) [R package]. https://CRAN.R-project.org/package=htmlwidgets."),
                   tags$li("R Core Team. (2024). tools: Tools for Package Development [R package]. R Foundation. https://CRAN.R-project.org/package=tools."),
                   tags$li("R Core Team. (2024). grid: The Grid Graphics Package [R package]. R Foundation. https://CRAN.R-project.org/package=grid."),
                   tags$li("R Core Team. (2024). graphics: The R Graphics Package [R package]. R Foundation. https://CRAN.R-project.org/package=graphics."),
                   tags$li("Csardi, G., & Nepusz, T. (2006). The igraph software package for complex network research. InterJournal, Complex Systems, 1695. https://igraph.org."),
                   tags$li("Schloerke, B., Crowley, J., Cook, D., Hofmann, H., & Wickham, H. (2024). GGally: Extension to 'ggplot2' (Version 2.1.2) [R package]. https://CRAN.R-project.org/package=GGally."),
                   tags$li("Auguie, B. (2017). gridExtra: Miscellaneous Functions for 'Grid' Graphics (Version 2.3) [R package]. https://CRAN.R-project.org/package=gridExtra."),
                   tags$li("Xie, Y. (2015). Dynamic Documents with R and knitr (2nd ed.). Chapman and Hall/CRC. https://yihui.org/knitr/."),
                   tags$li("Buels, R., & Yao, E. (2024). JBrowseR: R Interface to JBrowse (Version 1.0.0) [R package]. https://CRAN.R-project.org/package=JBrowseR."),
                   tags$li("Sass, J., & Schloerke, B. (2024). bslib: Custom 'Bootstrap' 'Sass' Themes for 'shiny' and 'rmarkdown' (Version 0.4.2) [R package]. https://CRAN.R-project.org/package=bslib."),
                   tags$li("Kay, M. (2024). ggdist: Visualizations of Distributions and Uncertainty (Version 3.0.1) [R package]. https://CRAN.R-project.org/package=ggdist."),
                   tags$li("Gohel, D. (2024). ggiraph: Make 'ggplot2' Graphics Interactive (Version 0.8.3) [R package]. https://CRAN.R-project.org/package=ggiraph.")
                 )
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
    tagList(
      # Instruction header with added help information.
      div(style = "text-align: center; margin-bottom: 30px;",
          h1("Welcome to INFO GRAB", style = "font-size: 3em; margin-bottom: 10px;"),
          p("Click on any tissue below to view detailed sample and gene information.", 
            style = "font-size: 1.2em; color: #555;"),
          p("Please visit the Help page/tab to get an overview and explanation of the functions offered.", 
            style = "font-size: 1.1em; color: #555;")
      ),
      # Loop over each system and its tissues.
      lapply(names(tissue_by_system), function(system) {
        tagList(
          h2(system, style = "text-align: center; margin-bottom: 20px;"),
          div(style = "display: flex; flex-wrap: wrap; justify-content: center; gap: 20px;",
              lapply(tissue_by_system[[system]], function(tissue) {
                local_tissue <- tissue  # capture the tissue name locally
                div(class = "tissue-card", 
                    style = "border: 1px solid #ccc; border-radius: 5px; padding: 10px; width: 200px; text-align: center;",
                    {
                      img_path <- paste0(gsub(" ", "_", local_tissue), ".png")
                      if (file.exists(file.path("www", img_path))) {
                        tags$img(src = img_path, 
                                 style = "width: 100px; height: auto; margin-bottom: 10px;")
                      } else {
                        tags$div("No image", style = "height: 100px; margin-bottom: 10px; color: #aaa;")
                      }
                    },
                    actionLink(inputId = paste0("tissue_", gsub(" ", "_", local_tissue)),
                               label = local_tissue,
                               style = "font-size: 1.2em; font-weight: bold; color: #337ab7; text-decoration: none;"),
                    p("Click to view details", style = "font-size: 0.9em; color: #666;")
                )
              })
          ),
          tags$hr(style = "margin-top:30px; margin-bottom:30px;")
        )
      })
    )
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
          
          tissue_data <- RPKM_data %>% dplyr::select(all_of(samples))
          
          # Top 5 Most Expressed Genes (RPKM)
          gene_expression_sums <- tissue_data %>% 
            rowSums(na.rm = TRUE) %>% 
            as.data.frame() %>% 
            setNames("Expression_Sum") %>% 
            rownames_to_column(var = "Gene") %>% 
            arrange(desc(Expression_Sum))
          top_rpkm_genes <- gene_expression_sums %>% head(5)
          top_rpkm_list <- paste(top_rpkm_genes$Gene, collapse = "<br>")
          
          # Top 5 ASE Genes (by Score)
          ase_data <- gene_ase_df()
          print(ase_data)
          top_ase_genes <- ase_data %>% 
            filter(Tissue == tissue) %>% 
            arrange(desc(mean.s)) %>%  # Using mean.s as score
            head(5)
          if (nrow(top_ase_genes) > 0) {
            top_ase_list <- paste(top_ase_genes$Gene, collapse = "<br>")
          } else {
            top_ase_list <- "No ASE genes found for this tissue"
          }
          
          # Top 5 Tissue-Specific (TS) Genes using tau score
          ts_data <- calcTau_custom(tissue_data)
          top_ts_genes <- ts_data %>% arrange(desc(tau)) %>% head(5)
          top_ts_list <- paste(top_ts_genes$gene, collapse = "<br>")
          
          showModal(modalDialog(
            title = paste(tissue),
            HTML(
              paste0(
                "Names: ", paste(samples_display, collapse = ", "), "<br><hr>",
                "Total Expressed Genes: ", expressed_genes_count, "<br><hr>",
                "<b>Top 5 Most Expressed Genes (RPKM):</b><br>", top_rpkm_list, "<br><hr>",
                "<b>Top 5 Tissue-Specific Genes (TS):</b><br>", top_ts_list, "<br><hr>",
                "<b>Top 5 ASE Genes (by Score):</b><br>", top_ase_list, "<br><hr>"
              )
            ),
            easyClose = TRUE,
            footer = NULL
          ))
        }, ignoreInit = TRUE)
      })
    })
  })
  
  
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
    filename = function() { 
      "DEA_data.csv" 
    },
    content = function(fname) {
      data <- filtered_data_combined()
      if ("log2FoldChange" %in% names(data)) {
        if (order_val() == 0) {
          names(data)[names(data) == "log2FoldChange"] <- 
            paste0("log2FC (", input$tissue_select1, " up / ", input$tissue_select2, " down)")
        } else {
          names(data)[names(data) == "log2FoldChange"] <- 
            paste0("log2FC (", input$tissue_select2, " up / ", input$tissue_select1, " down)")
        }
      }
      write.csv(data, fname, row.names = FALSE)
    }
  )
  
  
  observeEvent(input$add_all_de_cart, {
    data <- filtered_data_combined()
    req(data)
    
    all_de_genes <- data$Gene
    request_add_genes(all_de_genes)
    
    showNotification(paste("Genes grabbed."), type = "message")
  })
  
  
  observeEvent(input$add_to_cart_home, {
    genes <- gene_search_result()$Gene
    request_add_genes(genes)
    
    showNotification("Genes grabbed.", type = "message")
  })
  
  observeEvent(input$add_smallest_pvalues_cart, {
    data <- filtered_data_combined()
    req(data)
    
    smallest_pvalues_genes <- data %>% 
      arrange(padj) %>% 
      head(10) %>%
      pull(Gene)
    
    request_add_genes(smallest_pvalues_genes)
    
    showNotification("Genes grabbed.", type = "message")
  })
  
  observeEvent(input$add_largest_fc_cart, {
    data <- filtered_data_combined()
    req(data)
    
    largest_fc_genes <- data %>% 
      arrange(desc(abs(log2FoldChange))) %>% 
      head(10) %>%
      pull(Gene)
    
    request_add_genes(largest_fc_genes)
    
    showNotification(paste("Genes grabbed."), type = "message")
  })
  
  observeEvent(input$add_upregulated_cart, {
    data <- filtered_data_combined()
    req(data)
    
    upregulated_genes <- data %>% 
      arrange(desc(log2FoldChange)) %>% 
      head(10) %>%
      pull(Gene)
    
    request_add_genes(upregulated_genes)
    
    showNotification(paste("Genes grabbed."), type = "message")
  })
  
  observeEvent(input$add_downregulated_cart, {
    data <- filtered_data_combined()
    req(data)
    
    downregulated_genes <- data %>% 
      arrange(log2FoldChange) %>% 
      head(10) %>%
      pull(Gene)
    
    request_add_genes(downregulated_genes)
    
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
    data <- unfiltered_data() 
    
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
      request_add_genes(de_genes_tissue1)
      showNotification(paste("Grabbed", length(de_genes_tissue1), "DE genes upregulated in", input$tissue_select1, "."), type = "message")
    } else {
      showNotification(paste("No DE genes found upregulated in", input$tissue_select1), type = "warning")
    }
  })
  
  output$de_genes_tissue1_button <- renderUI({
    req(input$tissue_select1)
    actionButton("add_de_genes_tissue1", paste("Grab all upregulated genes from", input$tissue_select1))
  })
  
  output$de_genes_tissue2_button <- renderUI({
    req(input$tissue_select2)
    actionButton("add_de_genes_tissue2", paste("Grab all upregulated genes from", input$tissue_select2))
  })
  
  
  
  observeEvent(input$add_de_genes_tissue2, {
    data <- unfiltered_data()
    
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
      request_add_genes(de_genes_tissue2)
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
      request_add_genes(labeled_genes)
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
    
    annotation$Tissue <- as.character(annotation$Tissue)
    present_tissues <- unique(annotation$Tissue)
    
    missing_colors <- setdiff(present_tissues, names(mycolors1))
    if (length(missing_colors) > 0) {
      stop(paste("Missing colors for tissues:", paste(missing_colors, collapse = ", ")))
    }
    
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
    
    # Draw the heatmap
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
      treeheight_col = 0
    )
    
    # Add legend title
    #grid::grid.text(
    #label = "Z-Score",
    #x = .94,          
    #y = 0.786,       
    #gp = grid::gpar(fontsize = 12, fontface = "bold")
    #)
    
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
  
  # Declare a dedicated reactive value for the differential heatmap selection
  diff_selected_genes <- reactiveVal(NULL)
  
  observeEvent(input$heatmap_brush, {
    brush_data <- input$heatmap_brush
    if (!is.null(brush_data)) {
      heatmap_rows <- isolate(ordered_genes())
      total_height <- length(heatmap_rows)
      
      # Define a vertical offset to correct for the misalignment.
      # Adjust this value as needed (e.g., 0.1 for 10% of the plot height).
      y_offset <- 0.1
      
      # Adjust the brush y-values by the offset.
      adj_ymin <- brush_data$ymin - y_offset
      adj_ymax <- brush_data$ymax - y_offset
      
      # Ensure the adjusted values remain within [0,1].
      adj_ymin <- max(min(adj_ymin, 1), 0)
      adj_ymax <- max(min(adj_ymax, 1), 0)
      
      # Convert the adjusted normalized y-values to row indices.
      start_idx <- ceiling((1 - adj_ymax) * total_height)
      end_idx   <- floor((1 - adj_ymin) * total_height)
      
      # Bound the indices within the valid range.
      start_idx <- max(min(start_idx, total_height), 1)
      end_idx   <- max(min(end_idx, total_height), 1)
      
      # Swap indices if needed.
      if (start_idx > end_idx) {
        tmp <- start_idx
        start_idx <- end_idx
        end_idx <- tmp
      }
      
      sel_genes <- heatmap_rows[start_idx:end_idx]
      diff_selected_genes(sel_genes)
      
      showModal(modalDialog(
        title = "Selected Genes",
        HTML(paste(sel_genes, collapse = ", ")),
        footer = tagList(
          actionButton("add_to_cart_diff_modal", "Grab Genes"),
          modalButton("Close")
        ),
        easyClose = TRUE
      ))
    }
  })
  
  
  observeEvent(input$add_to_cart_diff_modal, {
    # Use isolate() to lock in the current selection
    genes_to_add <- isolate(diff_selected_genes())
    if (!is.null(genes_to_add) && length(genes_to_add) > 0) {
      request_add_genes(genes_to_add)  # Your function to add genes to the cart
      showNotification("Genes grabbed.", type = "message")
    } else {
      showNotification("No genes selected.", type = "warning")
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
    
    request_add_genes(ordered_genes())
    
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
        showNotification("No genes available after filtering. Using top tissue specific genes instead.", type = "warning")
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
      varname    <- rlang::as_label(mapping$x)
      tissue     <- phenodata$Tissue[match(varname, phenodata$Sample_name2)]
      color_val  <- mycolors1[[tissue]]
      
      ggplot2::ggplot(data = data, mapping = mapping) +
        ggplot2::geom_density(fill = color_val, alpha = 0.4) +
        ggplot2::theme_void() +
        ggplot2::theme(
          panel.background = ggplot2::element_rect(
            fill   = color_val,
            colour = "black",
            size   = 2
          )
        )
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
    
    ggp <- ggpairs(
      data = as.data.frame(selected_data),
      upper = list(continuous = create_upper),
      diag  = list(continuous = create_diag),
      lower = list(continuous = wrap("points", color = "black", alpha = 0.5))
    )
    
    n_tiss <- if (debounced_correlation_mode() == "selected") {
      length(unique(c(input$tissue_select1, input$tissue_select2)))
    } else {
      length(unique(phenodata$Tissue[phenodata$System == input$correlation_system]))
    }
    
    # pick a font size that scales down when there are many tissues
    strip_font <- if (n_tiss <= 5) {
      14
    } else if (n_tiss <= 10) {
      10
    } else {
      8
    }
    
    ggp +
      ggplot2::theme_minimal() +
      # force strips outside the plot panel
      ggplot2::theme(
        strip.placement    = "outside",
        strip.background   = ggplot2::element_blank(),
        # move top strip labels up a bit less aggressively
        strip.text.x = ggplot2::element_text(
          angle  = 0,
          hjust  = 0.5,
          vjust  = 0.5,
          size   = strip_font
        ),
        # make sure nothing gets clipped
        plot.margin     = ggplot2::margin(10, 10, 10, 10),
        panel.spacing.x = grid::unit(0.5, "lines"),
        # right strip labels horizontally, nudged left
        strip.text.y       = ggplot2::element_text(
          angle   = 0,
          hjust   = 1,
          vjust   = 0.5,
          size    = strip_font,
          margin  = ggplot2::margin(l = 6)  # add left margin
        ),
        panel.spacing.y    = grid::unit(0.5, "lines"),
        # legend down below, horizontal
        legend.position    = "bottom",
        legend.direction   = "horizontal",
        legend.title       = ggplot2::element_blank(),
        legend.text        = ggplot2::element_text(size = strip_font),
        legend.key.height  = grid::unit(0.4, "cm"),
        legend.key.width   = grid::unit(0.4, "cm")
      )
    
    
  })
  
  system_corr_plot <- eventReactive(input$confirm_cor_plot, {
    generate_correlation_plot()
  })
  
  output$correlation_plot <- renderPlot({
    shiny::withProgress(message = "Generating correlation plot...", value = 0, {
      shiny::incProgress(0.2)
      if (input$correlation_mode == "system") {
        req(system_corr_plot())
        p <- system_corr_plot()
      } else {
        p <- generate_correlation_plot()
      }
      shiny::incProgress(0.7)
      print(p)
      shiny::incProgress(1)
    })
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
        theme_minimal() +
        ggplot2::theme(
          axis.text = element_text(size = 14),  # Font size for axis tick labels
          axis.title = element_text(size = 16), # Font size for axis titles
          legend.text = element_text(size = 12), # Optional: Font size for legend text
          legend.title = element_text(size = 14) # Optional: Font size for legend titles
        )
      
      
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
      req(input$tissue_select1, input$tissue_select2, RPKM_data, phenodata)
      
      top_genes <- common_genes(input$tissue_select1, input$tissue_select2, RPKM_data, phenodata)
      
      if (length(top_genes) > 0) {
        updateTextInput(session, "barplot_searched_gene", value = paste(top_genes, collapse = ", "))
        
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
          axis.text.y = element_text(size = 18),
          axis.title.y = element_text(size = 24),
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
    request_add_genes(genes)
    
    showNotification("Genes grabbed.", type = "message")
  })
  
  observeEvent(input$transfer_to_heatmap, {
    transferred_genes <- strsplit(input$barplot_searched_gene, ",\\s*")[[1]]
    transferred_genes <- toupper(transferred_genes)
    
    if (length(transferred_genes) < 2) {
      showNotification(
        "Please select at least two genes to generate a heatmap.",
        type = "warning"
      )
      return()
    }
    
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
    
    print("Annotation")
    print(annotation)
    
    annotation$Tissue <- as.character(annotation$Tissue)
    present_tissues <- unique(annotation$Tissue)
    
    print("Present tissues")
    print(present_tissues)
    
    missing_colors <- setdiff(present_tissues, names(mycolors1))
    if (length(missing_colors) > 0) {
      stop(paste("Missing colors for tissues:", paste(missing_colors, collapse = ", ")))
    }
    
    annotation_colors <- list(Tissue = mycolors1[present_tissues])
    
    print("annotation colors")
    print(annotation_colors)
    
    sorted_samples <- annotation %>%
      arrange(Tissue) %>%
      rownames()
    
    rpkm_filtered <- rpkm_filtered[, sorted_samples]
    
    num_genes <- length(transferred_genes)
    cell_height <- ifelse(num_genes > 50, 5,
                          ifelse(num_genes > 30, 8,
                                 ifelse(num_genes > 20, 15,
                                        ifelse(num_genes > 10, 25,
                                               ifelse(num_genes > 5, 50, 100)))))
    
    show_genes <- num_genes <= 30
    
    showModal(
      modalDialog(
        title = "Gene Expression Heatmap",
        size = "l",
        plotOutput("heatmap_in_modal", width = "95%", height = "750px"),
        selectInput("color_scale_modal", "Select Color Scale",
                    choices = c("Blue-White-Red" = "blue_white_red",
                                "Green-Black-Red" = "green_black_red",
                                "Purple-White-Green" = "purple_white_green",
                                "Viridis" = "viridis",
                                "Plasma" = "plasma",
                                "Cividis" = "cividis",
                                "Inferno" = "inferno")),
        br(),
        downloadButton("download_modal_heatmap", "Download Heatmap"),
        easyClose = TRUE,
        footer = modalButton("Close"),
        tags$style(".modal-dialog { width: 90%; }")
      )
    )
    
    # Reactive expression to get the color scale based on user selection
    color_scale_reactive <- reactive({
      switch(input$color_scale_modal,
             "blue_white_red" = colorRampPalette(c("royalblue", "white", "firebrick3"))(100),
             "green_black_red" = colorRampPalette(c("springgreen2", "black", "firebrick2"))(100),
             "purple_white_green" = colorRampPalette(c("purple", "white", "springgreen4"))(100),
             "viridis" = viridis(100),
             "plasma" = plasma(100),
             "cividis" = cividis(100),
             "inferno" = inferno(100))
    })
    
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
        color = color_scale_reactive(),
        legend_labels = c("low", "medium", "high"),
        fontsize_col = 11,
        fontsize_row = 12,
        angle_col = 45,
        cellwidth = 23
      )
      
      grid::grid.text(
        label = "Z-Score",  
        x = .74,           
        y = 0.78,            
        gp = grid::gpar(fontsize = 10, fontface = "bold")  
      )
    })
    
    # Add the download handler for the heatmap
    output$download_modal_heatmap <- downloadHandler(
      filename = function() {
        paste("heatmap", ".png", sep = "")
      },
      content = function(file) {
        png(file, width = 15, height = 10, units = "in", res = 300)
        pheatmap(
          rpkm_filtered,
          scale = "row",
          cluster_rows = TRUE,
          cluster_cols = FALSE,
          show_rownames = TRUE,
          show_colnames = TRUE,
          annotation_col = annotation,
          annotation_colors = annotation_colors,
          color = color_scale_reactive(),
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
    
    
    # Process annotation data
    annotation <- phenodata %>%
      filter(Sample_name %in% samples) %>%
      dplyr::select(Sample_name, Tissue, System) %>%
      mutate(Sample_name = gsub("Sample_", "", Sample_name)) %>%
      column_to_rownames(var = "Sample_name")
    
    # Ensure Tissue and System columns are characters
    annotation$Tissue <- as.character(annotation$Tissue)
    annotation$System <- as.character(annotation$System)
    
    # Get unique tissues and systems
    present_tissues <- unique(annotation$Tissue)
    present_systems <- unique(annotation$System)
    
    print("Present tissues:")
    print(present_tissues)
    
    print("Present systems:")
    print(present_systems)
    
    # Check if all tissues and systems have colors in mycolors1 and mycolors2
    missing_tissue_colors <- setdiff(present_tissues, names(mycolors1))
    missing_system_colors <- setdiff(present_systems, names(mycolors2))
    
    # Stop with an error if there are missing colors
    if (length(missing_tissue_colors) > 0) {
      stop(paste("Missing colors for tissues:", paste(missing_tissue_colors, collapse = ", ")))
    }
    if (length(missing_system_colors) > 0) {
      stop(paste("Missing colors for systems:", paste(missing_system_colors, collapse = ", ")))
    }
    
    # Define annotation colors
    annotation_colors <- list(
      Tissue = mycolors1[present_tissues],
      System = mycolors2[present_systems]
    )
    
    print("Annotation colors:")
    print(annotation_colors)
    
    
    print(annotation_colors)
    
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
  
  # Create a separate reactive value for tissue-specific selected genes
  tissue_specific_selected_genes <- reactiveVal(NULL)
  
  # Create a debounced reactive for the brush input (500ms delay)
  tissue_specific_brush <- debounce(reactive(input$heatmap_tissue_specific_brush), 500)
  
  observeEvent(tissue_specific_brush(), {
    brush_data <- tissue_specific_brush()
    if (!is.null(brush_data)) {
      heatmap_rows <- ordered_tissue_specific_genes()
      total_height <- length(heatmap_rows)
      
      start_row_index <- total_height - ceiling(brush_data$ymin * total_height) + 1
      end_row_index   <- total_height - ceiling(brush_data$ymax * total_height) + 1
      start_row_index <- max(min(start_row_index, total_height), 1)
      end_row_index   <- max(min(end_row_index, total_height), 1)
      
      selected_rows <- heatmap_rows[end_row_index:start_row_index]
      if(length(selected_rows) == 0){
        showNotification("No genes selected from brush.", type = "warning")
        return()
      }
      tissue_specific_selected_genes(selected_rows)
      
      showModal(modalDialog(
        title = "Selected Genes",
        paste(selected_rows, collapse = ", "),
        footer = tagList(
          actionButton("add_to_cart_tissue_specific_modal", "Grab Genes"),
          modalButton("Close")
        ),
        easyClose = FALSE  # Prevent auto-closing when clicking outside
      ))
    }
  })
  
  
  observeEvent(input$add_to_cart_tissue_specific_modal, {
    genes_to_add <- tissue_specific_selected_genes()
    if (!is.null(genes_to_add) && length(genes_to_add) > 0) {
      request_add_genes(genes_to_add)
      showNotification("Genes grabbed.", type = "message")
    } else {
      showNotification("No genes selected.", type = "warning")
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
  
  
  observeEvent(input$add_to_cart_tissue_specific_heatmap, {
    selected_genes <- isolate(ordered_tissue_specific_genes())
    if(length(selected_genes) == 0){
      showNotification("No genes are currently displayed.", type = "warning")
    } else {
      request_add_genes(selected_genes)
      showNotification("Genes grabbed.", type = "message")
    }
  })
  
  
  ##############
  
  # Raincloud Plot
  
  # 1) Load ASE data --------------------------------------------------------
  
  # SNP ASE
  # 2) Load SNP ASE data
  snp_ase_path <- "data/all.snpAse_0906.txt"
  snp_ase_df   <- shiny::reactiveVal(data.frame())
  shiny::observe({
    req(file.exists(snp_ase_path))
    withProgress(message = "Loading SNP ASE...", value = 0, {
      df <- read.delim(snp_ase_path, stringsAsFactors = FALSE)
      if (!"AE" %in% colnames(df)) {
        df <- df %>% mutate(
          Rratio = REF_COUNT/(REF_COUNT+ALT_COUNT),
          AE     = abs(0.5 - Rratio)
        )
      }
      
      # Recode here
      df$Tissue <- str_replace_all(df$Tissue, repl)
      
      if (!"Tissue" %in% colnames(df)) stop("Missing Tissue column")
      snp_ase_df(df)
      
      updateSelectInput(session, "tissue_select_raincloud",
                        choices  = c("All Tissues", sort(unique(df$Tissue))),
                        selected = "All Tissues"
      )
      
      incProgress(1)
    })
  })
  
  
  # 2) Search / random helpers -----------------------------------------------
  
  ## 3) Helpers for clearing & random selection ####
  selected_genes <- shiny::reactiveVal(NULL)
  selected_snps  <- shiny::reactiveVal(NULL)
  
  # clear both gene & SNP selections
  # in your server, replace the existing clear_search_raincloud observer
  # with this one (or just add the two extra lines):
  
  # in server.R, after you’ve declared selected_genes/selected_snps:
  shiny::observeEvent(input$clear_search_raincloud, {
    # 1) Clear out whatever was typed in the gene & SNP boxes:
    shiny::updateTextInput(session,
                           inputId = "raincloud_searched_gene",
                           value   = "")
    shiny::updateTextInput(session,
                           inputId = "raincloud_searched_snp",
                           value   = "")
    
    # 2) Reset the tissue picker back to ALL:
    #    (must match exactly the “All Tissues” choice in your selectInput)
    shiny::updateSelectInput(session,
                             inputId  = "tissue_select_raincloud",
                             selected = "All Tissues")
    
    # 3) Clear any internal reactiveVals just in case you’re also
    #    using those to drive logic (won’t hurt to reset them):
    selected_genes(NULL)
    selected_snps(NULL)
  })
  
  # random genes
  # random genes
  shiny::observeEvent(input$random_genes_raincloud, {
    df <- gene_ase_df()
    if (nrow(df) >= input$top_n_raincloud) {
      genes <- sample(unique(df$Gene), size = input$top_n_raincloud)
      selected_genes(genes)
      shiny::updateTextInput(
        session,
        inputId = "raincloud_searched_gene",
        value   = paste(genes, collapse = ",")
      )
    }
  })
  
  # random SNPs
  shiny::observeEvent(input$random_snps_raincloud, {
    df <- snp_ase_df()
    if (nrow(df) >= input$top_n_raincloud) {
      ids <- sample(unique(df$ID), size = input$top_n_raincloud)
      selected_snps(ids)
      shiny::updateTextInput(
        session,
        inputId = "raincloud_searched_snp",
        value   = paste(ids, collapse = ",")
      )
    }
  })
  
  ## ——— 1) Load Gene ASE data ——————————————————————————————————————
  # — 1) Load Gene ASE data, coerce score unambiguously —
  gene_ase_path <- "data/all.geneAse.info827.txt"
  gene_ase_df   <- shiny::reactiveVal(data.frame())
  shiny::observe({
    if (!file.exists(gene_ase_path)) return()
    withProgress(message = "Loading Gene ASE...", value = 0, {
      df <- read.delim(gene_ase_path, stringsAsFactors = FALSE)
      df$mean.s <- as.numeric(df$mean.s)
      df <- df[, c("Gene","n.vars","mean.s","fdr","Tissue")]
      
      # Recode here
      df$Tissue <- str_replace_all(df$Tissue, repl)
      
      gene_ase_df(df)
      updateSelectInput(session, "tissue_select_raincloud",
                        choices  = c("All Tissues", sort(unique(df$Tissue))),
                        selected = "All Tissues"
      )
      incProgress(1)
    })
  })
  
  
  ## ——— 3) The unified raincloud‐plot function ————————————————————————
  make_raincloud_plot <- function(data,
                                  yvar,
                                  show_mean_line,
                                  show_points,
                                  facet_by_tissue = FALSE) {
    # Deduplicate
    data <- dplyr::distinct(data)
    
    # Base ggplot
    p <- ggplot2::ggplot(
      data,
      ggplot2::aes_string(
        x = if (facet_by_tissue) "Tissue" else "1",
        y = yvar
      )
    )
    
    # Add interactive points (horizontal jitter small, vertical jitter tiny)
    if (show_points) {
      idvar   <- if (yvar == "mean.s") "Gene" else "ID"
      tooltip <- if (yvar == "mean.s") {
        "paste0('Gene: ', Gene, '<br>Mean effect: ', round(mean.s, 3))"
      } else {
        "paste0('SNP: ', ID, '<br>AE: ', round(AE, 3))"
      }
      p <- p + ggiraph::geom_point_interactive(
        ggplot2::aes_string(
          x       = "0.87",
          color   = idvar,
          tooltip = tooltip,
          data_id = idvar
        ),
        position = ggplot2::position_jitter(width = 0.02, height = 0.03),
        shape    = 16,
        alpha    = 0.7,
        size     = 2
      )
    }
    
    # Half-violin (if enough points)
    if (nrow(data) > 2) {
      p <- p + gghalves::geom_half_violin(
        side   = "r",
        trim   = FALSE,
        fill   = "grey",
        alpha  = 0.6,
        color  = "black"
      )
    }
    
    # Boxplot
    p <- p + ggplot2::geom_boxplot(
      width         = 0.1,
      outlier.shape = NA,
      alpha         = 0.6,
      color         = "black"
    )
    
    # Optional mean line
    if (show_mean_line) {
      p <- p + ggplot2::stat_summary(
        fun   = mean,
        geom  = "crossbar",
        width = 0.075,
        color = "red"
      )
    }
    
    # Flip coordinates, theme, labels
    p <- p +
      ggplot2::coord_flip() +
      ggplot2::theme_minimal(base_size = 14) +
      ggplot2::labs(
        x     = NULL,
        y     = if (yvar == "mean.s") "Mean Effect" else "AE",
        title = if (yvar == "mean.s")
          "Raincloud Plot: Gene ASE (mean effect)"
        else
          "Raincloud Plot: SNP ASE (AE)"
      ) +
      ggplot2::theme(
        axis.ticks.y    = ggplot2::element_blank(),
        plot.title      = ggplot2::element_text(hjust = 0.5),
        legend.position = "none"
      )
    
    # Facet by tissue if requested
    if (facet_by_tissue) {
      p <- p + ggplot2::facet_wrap(~ Tissue, scales = "free", ncol = 1)
    }
    
    # Return interactive ggiraph object
    ggiraph::girafe(
      ggobj     = p,
      width_svg = 8,
      height_svg= 6
    )
  }
  
  # 4) Reactive selection of gene vs snp data --------------------------------
  
  ## 4) Reactive dataset with proper Tissue & search filtering ####
  filtered_data_raincloud <- shiny::reactive({
    # pick the right dataframe
    df <- if (input$raincloud_data_type == "gene") {
      gene_ase_df() %>% dplyr::filter(fdr <= input$pval_cutoff)
    } else {
      snp_ase_df()  %>% dplyr::filter(FDR <= input$pval_cutoff)
    }
    
    # tissue filter
    sel_t <- input$tissue_select_raincloud
    if (!is.null(sel_t) && !"All Tissues" %in% sel_t) {
      df <- df %>% dplyr::filter(Tissue %in% sel_t)
    }
    
    # SEARCH: priority to random selection, else textInput
    if (input$raincloud_data_type == "gene") {
      
      if (!is.null(selected_genes())) {
        df <- df %>% dplyr::filter(Gene %in% selected_genes())
      } else {
        genes <- strsplit(input$raincloud_searched_gene, ",\\s*")[[1]]
        if (length(genes) && nzchar(genes[1])) {
          df <- df %>% dplyr::filter(Gene %in% genes)
        }
      }
      # top N genes by mean effect
      topN <- df %>%
        dplyr::group_by(Gene) %>%
        dplyr::summarise(m = mean(score, na.rm = TRUE), .groups = "drop") %>%
        dplyr::arrange(dplyr::desc(m)) %>%
        head(input$top_n_raincloud) %>%
        dplyr::pull(Gene)
      df <- df %>% dplyr::filter(Gene %in% topN)
      
    } else {
      
      if (!is.null(selected_snps())) {
        df <- df %>% dplyr::filter(ID %in% selected_snps())
      } else {
        ids <- strsplit(input$raincloud_searched_snp, ",\\s*")[[1]]
        if (length(ids) && nzchar(ids[1])) {
          df <- df %>% dplyr::filter(ID %in% ids)
        }
      }
      # top N SNPs by AE
      topN <- df %>%
        dplyr::group_by(ID) %>%
        dplyr::summarise(ae = mean(AE, na.rm = TRUE), .groups = "drop") %>%
        dplyr::arrange(dplyr::desc(ae)) %>%
        head(input$top_n_raincloud) %>%
        dplyr::pull(ID)
      df <- df %>% dplyr::filter(ID %in% topN)
      
    }
    df
  })
  
  # 5) Render, download & grab selections -----------------------------------
  
  output$raincloud_plot <- ggiraph::renderGirafe({
    df    <- filtered_data_raincloud()
    shiny::req(nrow(df) > 0)
    facet <- !"All Tissues" %in% input$tissue_select_raincloud
    yvar  <- if (input$raincloud_data_type == "gene") "mean.s" else "AE"
    
    make_raincloud_plot(
      data            = df,
      yvar            = yvar,
      show_mean_line  = input$show_mean_line,
      show_points     = input$show_points,
      facet_by_tissue = facet
    )
  })
  
  output$download_raincloud_plot <- shiny::downloadHandler(
    filename = function()   "raincloud.png",
    content  = function(file) {
      df    <- filtered_data_raincloud()
      facet <- !("All Tissues" %in% input$tissue_select_raincloud)
      yvar  <- if (input$raincloud_data_type=="gene") "score" else "AE"
      p     <- make_static_raincloud_plot(df, yvar, input$show_mean_line, input$show_points, facet)
      ggplot2::ggsave(file, plot = p, width = 10, height = 6)
    }
  )
  
  ## only grab when in gene‐ASE mode:
  shiny::observeEvent(input$grab_selected_raincloud, {
    req(input$raincloud_data_type == "gene")
    sel <- input$raincloud_plot_selected
    if (length(sel)) {
      request_add_genes(sel)
      shiny::showNotification(paste(length(sel), "genes grabbed"), type="message")
    }
  })
  
  ## and same for “Grab All Displayed”:
  shiny::observeEvent(input$grab_displayed_raincloud, {
    req(input$raincloud_data_type == "gene")
    df  <- filtered_data_raincloud()
    ids <- unique(df$Gene)
    if (length(ids)) {
      request_add_genes(ids)
      shiny::showNotification(paste(length(ids), "genes added"), type="message")
    }
  })
  
  
  
  #############
  
  # Genome Browser
  
  ###########################################
  # PART A: Tissue-Gene info from the TXT
  ###########################################
  
  gene_ase_df <- reactiveVal(data.frame())
  
  observe({
    info_path <- "browser_data_for_app/all.geneAse.info827.txt"
    req(file.exists(info_path))
    df_info <- read.delim(info_path, stringsAsFactors = FALSE)
    
    # Recode here
    df_info$Tissue <- str_replace_all(df_info$Tissue, repl)
    
    gene_ase_df(df_info)
  })
  
  
  # Add "All Tissues" to the checkbox
  observe({
    df <- gene_ase_df()
    req(df)
    
    all_tissues <- sort(unique(df$Tissue))
    # Insert "All Tissues" at the beginning
    all_tissues <- c("All Tissues", all_tissues)
    
    updateCheckboxGroupInput(
      session,
      "ase_tissue_filter",
      choices  = all_tissues,
      selected = NULL
    )
  })
  
  
  ###########################################
  # PART B: Initialize IGV + Searching
  ###########################################
  
  observeEvent("GCF_016699485.2", {
    genome  <- "GCF_016699485.2"
    igv_opts <- parseAndValidateGenomeSpec(genomeName = genome)
    
    output$igvShiny <- renderIgvShiny({
      igvShiny(igv_opts)
    })
  })
  
  observeEvent(input$search_button, {
    req(input$genome_browser_search)
    showGenomicRegion(session, "igvShiny", region=input$genome_browser_search)
  })
  
  
  ###########################################
  # CHROMOSOME MAPPING (Optional)
  ###########################################
  chrom_map_path <- "browser_data_for_app/chrm-ncbi.txt"
  chrom_map <- NULL
  if (file.exists(chrom_map_path)) {
    chrom_map <- readr::read_tsv(chrom_map_path, col_names=c("Chrm_ncbi","Chromosome"))
  } else {
    showNotification("chrm-ncbi.txt not found. Chrom mapping disabled.", type="error")
  }
  
  
  ###########################################
  # CUSTOM TRACKS (GFF, BAM, BED, VCF)
  ###########################################
  
  convert_lists_to_chars <- function(df) {
    df[] <- lapply(df, function(x) if (is.list(x)) sapply(x, toString) else x)
    return(df)
  }
  
  observeEvent(input$loadTrack, {
    withProgress(message='Loading track...', value=0, {
      incProgress(0.2)
      
      if (input$input_type == "file" && !is.null(input$file)) {
        file_path <- input$file$datapath
        file_type <- tools::file_ext(input$file$name)
        tryCatch({
          if (file_type %in% c("gff3","gff")) {
            tbl.gff <- rtracklayer::import(file_path, format=file_type)
            tbl.gff <- convert_lists_to_chars(as.data.frame(tbl.gff))
            loadGFF3TrackFromLocalData(
              session=session,
              id="igvShiny",
              trackName=input$file$name,
              tbl.gff3=tbl.gff,
              color="gray",
              deleteTracksOfSameName=TRUE
            )
            showNotification("GFF/GFF3 track loaded successfully.", type="message")
            
          } else if (file_type == "bam") {
            bamFile <- Rsamtools::BamFile(file_path)
            loadBamTrackFromLocalData(
              session=session,
              id="igvShiny",
              trackName="Uploaded BAM Track",
              bamFile=bamFile,
              color="grey",
              deleteTracksOfSameName=TRUE
            )
            showNotification("BAM Track loaded successfully.", type="message")
            
          } else if (file_type == "bed") {
            tbl.bed <- rtracklayer::import(file_path, format="bed")
            tbl.bed <- as.data.frame(tbl.bed)
            if ("seqnames" %in% colnames(tbl.bed)) {
              colnames(tbl.bed)[colnames(tbl.bed)=="seqnames"] <- "chr"
            }
            loadBedTrack(
              session=session,
              id="igvShiny",
              trackName=input$file$name,
              tbl=tbl.bed,
              color="gray",
              deleteTracksOfSameName=TRUE
            )
            showNotification("BED Track loaded successfully.", type="message")
            
          } else if (file_type == "vcf") {
            vcf_data <- VariantAnnotation::readVcf(file_path)
            loadVcfTrack(
              session=session,
              id="igvShiny",
              trackName=input$file$name,
              vcfData=vcf_data,
              deleteTracksOfSameName=TRUE
            )
            showNotification("VCF Track loaded successfully.", type="message")
            
          } else {
            showNotification("Unsupported file type.", type="error")
          }
          incProgress(0.8)
        }, error=function(e) {
          showNotification(paste("Error loading track:", e$message), type="error")
        })
        
      } else if (input$input_type == "url" && !is.null(input$url)) {
        showNotification("Loading track from URL not fully implemented here.", type="message")
        incProgress(0.8)
      } else {
        showNotification("No file or URL provided.", type="warning")
      }
      
      incProgress(1)
    })
  })
  
  
  ###########################################
  # LOAD ANNOTATION (GFF)
  ###########################################
  observeEvent(input$load_gff, {
    gff_file <- "browser_data_for_app/genomeAnnoatation_gallus.updated.gff"
    if (!file.exists(gff_file)) {
      showNotification("GFF file not found.", type="error")
      return(NULL)
    }
    
    withProgress(message='Loading GFF Track...', value=0, {
      tryCatch({
        gff_data <- rtracklayer::import(gff_file, format="gff")
        gff_df   <- as.data.frame(gff_data)
        
        if (!is.null(chrom_map)) {
          if ("seqnames" %in% colnames(gff_df)) {
            colnames(gff_df)[colnames(gff_df)=="seqnames"] <- "chr"
          }
          gff_df$chr <- as.character(gff_df$chr)
          gff_df <- gff_df %>%
            left_join(chrom_map, by=c("chr"="Chrm_ncbi")) %>%
            mutate(chr=Chromosome) %>%
            dplyr::select(-Chromosome)
        }
        
        if (!"gene_name" %in% colnames(gff_df)) {
          gff_df$gene_name <- "NA"
        } else {
          gff_df$gene_name[is.na(gff_df$gene_name)] <- "NA"
        }
        
        gff_df <- gff_df %>%
          dplyr::select(chr, start, end, gene_name, score, strand) %>%
          mutate(
            chr   = as.character(chr),
            start = as.integer(start),
            end   = as.integer(end),
            score = as.numeric(score),
            strand= as.character(strand)
          )
        
        loadBedTrack(
          session=session,
          id="igvShiny",
          trackName="Annotation Track",
          tbl=gff_df,
          color="darkred",
          trackHeight=70,
          deleteTracksOfSameName=TRUE
        )
      }, error=function(e) {
        showNotification(paste("Error loading GFF Track:", e$message), type="error")
      })
    })
  })
  
  
  ###########################################
  # LOAD GENE ASE TRACKS (SINGLE or MULTIPLE)
  ###########################################
  observeEvent(input$load_geneASE_track, {
    # 1) Check Tissue selection
    if (is.null(input$ase_tissue_filter) || length(input$ase_tissue_filter) == 0) {
      showNotification("No tissues selected! Please pick at least one tissue.", type = "warning")
      return(NULL)
    }
    
    df_info <- gene_ase_df()
    required_cols_df <- c("Tissue", "Gene")
    if (!all(required_cols_df %in% names(df_info))) {
      showNotification(
        paste("df_info is missing:", 
              paste(setdiff(required_cols_df, names(df_info)), collapse = ", ")),
        type = "error"
      )
      return(NULL)
    }
    
    if (nrow(df_info) == 0) {
      showNotification("No Tissue info loaded from TXT.", type = "error")
      return(NULL)
    }
    
    chosen_tissues <- input$ase_tissue_filter
    if ("All Tissues" %in% chosen_tissues) {
      chosen_tissues <- sort(unique(df_info$Tissue))
    }
    
    # 2) Load the BED
    bed_path <- "browser_data_for_app/gene_ase_converted_cleaned_with_strand.bed"
    if (!file.exists(bed_path)) {
      showNotification("Gene ASE BED file not found.", type = "error")
      return(NULL)
    }
    
    withProgress(message = 'Loading Gene ASE Track(s)...', value = 0, {
      incProgress(0.2)
      tryCatch({
        ase_bed <- rtracklayer::import(bed_path, format = "bed")
        bed_df  <- as.data.frame(ase_bed)
        if (!is.data.frame(bed_df)) {
          showNotification("Could not convert BED to a data frame.", type = "error")
          return(NULL)
        }
        
        if (nrow(bed_df) == 0) {
          showNotification("Imported BED data is empty.", type = "error")
          return(NULL)
        }
        
        # Rename columns
        if ("seqnames" %in% names(bed_df)) {
          names(bed_df)[names(bed_df) == "seqnames"] <- "chr"
        }
        if ("name" %in% names(bed_df)) {
          names(bed_df)[names(bed_df) == "name"] <- "Gene"
        }
        if (!"Gene" %in% names(bed_df)) {
          showNotification("No 'Gene' column in BED.", type = "error")
          return(NULL)
        }
        
        # Re-check data frame status
        if (!is.data.frame(bed_df)) {
          showNotification("bed_df is no longer a data frame after renaming.", type = "error")
          return(NULL)
        }
        
        needed_cols_bed <- c("chr", "start", "end", "Gene")
        if (!all(needed_cols_bed %in% names(bed_df))) {
          showNotification(
            paste("BED missing:", 
                  paste(setdiff(needed_cols_bed, names(bed_df)), collapse = ", ")),
            type = "error"
          )
          return(NULL)
        }
        
        bed_df$Gene <- toupper(bed_df$Gene)
        if (!is.data.frame(bed_df)) {
          showNotification("bed_df is no longer a data frame after Gene conversion.", type = "error")
          return(NULL)
        }
        
        color_map <- setNames(rainbow(length(chosen_tissues)), chosen_tissues)
        track_mode <- input$track_mode
        if (!track_mode %in% c("single", "multiple")) {
          showNotification("Invalid track mode.", type = "error")
          return(NULL)
        }
        
        # SINGLE TRACK
        if (track_mode == "single") {
          combined_genes <- df_info %>%
            dplyr::filter(Tissue %in% chosen_tissues) %>%
            dplyr::pull(Gene) %>%
            toupper() %>%
            unique()
          
          if (length(combined_genes) == 0) {
            showNotification("No Genes found for these Tissues.", type = "warning")
            return(NULL)
          }
          
          single_bed <- bed_df %>% dplyr::filter(Gene %in% combined_genes)
          if (nrow(single_bed) == 0) {
            showNotification("No BED entries match these Genes.", type = "warning")
            return(NULL)
          }
          
          final_bed <- single_bed %>%
            dplyr::rename(name = Gene) %>%
            dplyr::select(chr, start, end, name, score, strand) %>%
            dplyr::mutate(
              chr   = as.character(chr),
              start = as.integer(start),
              end   = as.integer(end),
              name  = as.character(name),
              score = as.numeric(score),
              strand= as.character(strand)
            )
          
          loadBedTrack(
            session  = session,
            id       = "igvShiny",
            trackName= "Gene ASE: All Selected Tissues (Single Track)",
            tbl      = final_bed,
            color    = "purple",
            trackHeight = 70,
            deleteTracksOfSameName = FALSE
          )
          
          # MULTIPLE TRACKS
        } else {
          for (tissue in chosen_tissues) {
            tissue_genes <- df_info %>%
              dplyr::filter(Tissue == tissue) %>%
              dplyr::pull(Gene) %>%
              toupper() %>%
              unique()
            
            if (length(tissue_genes) == 0) {
              showNotification(paste("No Genes for Tissue:", tissue), type = "warning")
              next
            }
            
            subset_bed <- bed_df %>% dplyr::filter(Gene %in% tissue_genes)
            if (nrow(subset_bed) == 0) {
              showNotification(paste("No BED match for Tissue:", tissue), type = "warning")
              next
            }
            
            final_bed <- subset_bed %>%
              dplyr::rename(name = Gene) %>%
              dplyr::select(chr, start, end, name, score, strand) %>%
              dplyr::mutate(
                chr   = as.character(chr),
                start = as.integer(start),
                end   = as.integer(end),
                name  = as.character(name),
                score = as.numeric(score),
                strand= as.character(strand)
              )
            
            loadBedTrack(
              session  = session,
              id       = "igvShiny",
              trackName= paste("Gene ASE:", tissue),
              tbl      = final_bed,
              color    = color_map[[tissue]],
              trackHeight = 70,
              deleteTracksOfSameName = FALSE
            )
          }
        }
        
        incProgress(1)
        showNotification(paste("Loaded Gene ASE track(s). Mode:", track_mode), type = "message")
      }, error = function(e) {
        showNotification(paste("Error loading Gene ASE tracks:", e$message), type = "error")
      })
    })
  })
  
  
  
  
  ###########################################
  # LOAD SNP ASE TRACK (Ignoring Tissue)
  ###########################################
  observeEvent(input$load_snp, {
    showNotification("Loading SNP ASE track. This may take a few minutes...", type="message")
    
    vcf_path <- "browser_data_for_app/ASE.SNP.gallus.updated02.vcf"
    if (!file.exists(vcf_path)) {
      showNotification("SNP VCF file not found.", type="error")
      return(NULL)
    }
    
    withProgress(message='Loading SNP ASE VCF...', value=0, {
      incProgress(0.2)
      tryCatch({
        vcf_data <- readVcf(vcf_path, "GCF_016699485.2")
        snp_data <- rowRanges(vcf_data)
        ff       <- as.data.frame(fixed(vcf_data))
        
        snp_df <- data.frame(
          chr   = as.character(seqnames(snp_data)),
          start = start(snp_data),
          end   = end(snp_data),
          name  = paste0(ff$REF, "/", sapply(ff$ALT, function(x) paste(x, collapse=","))),
          score = 0,
          strand= "+"
        )
        
        snp_df <- snp_df %>%
          mutate(
            chr   = as.character(chr),
            start = as.integer(start),
            end   = as.integer(end),
            name  = as.character(name),
            score = as.numeric(score),
            strand= as.character(strand)
          )
        
        loadBedTrack(
          session=session,
          id="igvShiny",
          trackName="SNP ASE Track",
          tbl=snp_df,
          color="darkgreen",
          trackHeight=70,
          deleteTracksOfSameName=TRUE
        )
        
        incProgress(1)
        showNotification("SNP ASE track loaded successfully.", type="message")
        
      }, error=function(e) {
        showNotification(paste("Error loading SNP ASE track:", e$message), type="error")
      })
    })
  })
  
  
  ###########################################
  # GRAB ALL ASE GENES
  ###########################################
  observeEvent(input$grab_all_ase_genes, {
    bed_path <- "browser_data_for_app/gene_ase_converted_cleaned_with_strand.bed"
    if (!file.exists(bed_path)) {
      showNotification("Gene ASE bed file not found for 'Grab All ASE Genes'.", type = "error")
      return(NULL)
    }
    
    withProgress(message = 'Grabbing All ASE Genes...', value = 0, {
      tryCatch({
        bed_data <- rtracklayer::import(bed_path, format = "bed")
        bed_df <- as.data.frame(bed_data)
        all_ase_genes <- character(0)
        
        if ("name" %in% colnames(bed_df)) {
          all_ase_genes <- toupper(bed_df$name)
        }
        
        if (length(all_ase_genes) == 0) {
          showNotification("No ASE genes found in the bed file.", type = "error")
          return(NULL)
        }
        
        # Add all ASE genes to the gene cart
        request_add_genes(all_ase_genes)
        
        showNotification(paste("All ASE Genes grabbed:", length(all_ase_genes)),
                         type = "message")
        incProgress(1)
      }, error = function(e) {
        showNotification(paste("Error grabbing ASE genes:", e$message),
                         type = "error")
      })
    })
  })
  
  
  ##############
  
  # Genes
  
  previous_tab <- reactiveVal("de_tab")  # Initial assumption
  
  # Observe navigation changes
  observeEvent(input$navbar, {
    current_tab <- input$navbar
    
    # If the user navigates to the Genes tab
    if (current_tab == "genes_tab") {
      # Determine which mode to select based on the previous tab
      prev <- previous_tab()
      
      if (prev == "de_tab") {
        updateSelectInput(session, "gene_mode", selected = "de")
      } else if (prev == "tissue_tab") {
        updateSelectInput(session, "gene_mode", selected = "tissue")
      } else if (prev == "ase_tab") {
        updateSelectInput(session, "gene_mode", selected = "ase")
      } else {
        # Default mode if coming from some other tab or unknown state
        updateSelectInput(session, "gene_mode", selected = "tissue")
      }
    }
    
    # Update previous_tab for the next navigation
    previous_tab(current_tab)
  })
  
  new_genes_to_add <- reactiveVal(NULL)
  
  # Low-level function: adds genes without prompting
  add_to_cart <- function(new_genes) {
    current_cart <- cart_genes()
    updated_cart <- unique(c(current_cart, new_genes))
    cart_genes(updated_cart)
  }
  
  # Function to request adding genes, showing a modal
  # Updated request_add_genes function
  request_add_genes <- function(new_genes) {
    new_genes <- toupper(new_genes)
    new_genes_to_add(new_genes)
    current_count <- length(cart_genes())
    new_count <- length(new_genes)
    showModal(modalDialog(
      title = "Add Genes to Cart",
      HTML(paste0(
        "Currently, you have <b>", current_count, "</b> gene(s) in the cart.<br>",
        "You are about to add <b>", new_count, "</b> new gene(s).<br><br>",
        "What would you like to do?"
      )),
      footer = tagList(
        actionButton("add_only_btn", "Add to Existing"),
        actionButton("clear_and_add_btn", "Clear and Add"),
        modalButton("Cancel")
      ),
      easyClose = FALSE
    ))
  }
  
  # When the user confirms gene addition via "Add to Existing"
  observeEvent(input$add_only_btn, {
    req(new_genes_to_add())
    add_to_cart(new_genes_to_add())
    new_genes_to_add(NULL)
    showNotification("Genes have been added.", type = "message")
    shinyjs::runjs("
      var geneTab = $('a[data-value=\"genes_tab\"]');
      geneTab.addClass('blink');
      setTimeout(function(){ geneTab.removeClass('blink'); }, 2000);
  ")
    removeModal()
  })
  
  # When the user confirms gene addition via "Clear and Add"
  observeEvent(input$clear_and_add_btn, {
    req(new_genes_to_add())
    cart_genes(character(0))  # Clear cart
    add_to_cart(new_genes_to_add())
    new_genes_to_add(NULL)
    showNotification("Cart cleared, new genes added.", type = "message")
    shinyjs::runjs("
      var geneTab = $('a[data-value=\"genes_tab\"]');
      geneTab.addClass('blink');
      setTimeout(function(){ geneTab.removeClass('blink'); }, 2000);
  ")
    removeModal()
  })
  
  
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
    gene_ase_path <- "data/gene_ase_converted_cleaned_with_strand.bed"
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
    gene_ase_path <- "data/gene_ase_converted_cleaned_with_strand.bed"
    gene_ase_data <- rtracklayer::import(gene_ase_path, format = "bed")
    gene_ase_df <- as.data.frame(gene_ase_data)
    
    return(unique(gene_ase_df$name))
  })
  
  # In the server function, after you've defined all required reactive values and before output$cart_gene_table:
  
  # Reactive expression to track current mode
  current_gene_mode <- reactive({
    req(input$gene_mode)
    input$gene_mode
  })
  
  # Reactive expression to generate the table data based on the selected mode
  mode_based_gene_table <- reactive({
    mode <- current_gene_mode()
    genes <- cart_genes()
    
    # If cart is empty
    if (length(genes) == 0) {
      return(data.frame(Message = "No genes in cart"))
    }
    
    if (mode == "de") {
      # — your existing DE branch unchanged —
      de_data <- unfiltered_data()
      if (is.null(de_data) || nrow(de_data) == 0) {
        return(data.frame(Message = "No DE data available"))
      }
      gene_df <- data.frame(Gene = genes, stringsAsFactors = FALSE)
      if (!"Gene" %in% colnames(de_data)) {
        de_data <- tibble::rownames_to_column(de_data, var = "Gene")
      }
      merged <- dplyr::left_join(gene_df, de_data, by = "Gene")
      # rename log2FoldChange column
      if ("log2FoldChange" %in% names(merged)) {
        if (order_val() == 0) {
          names(merged)[names(merged) == "log2FoldChange"] <-
            paste0("log2FC (", input$tissue_select1, " up / ", input$tissue_select2, " down)")
        } else {
          names(merged)[names(merged) == "log2FoldChange"] <-
            paste0("log2FC (", input$tissue_select2, " up / ", input$tissue_select1, " down)")
        }
      }
      return(merged)
      
    } else if (mode == "tissue") {
      
      # 1) full tissue × gene matrix
      full_mat <- gene_expression_data_tissue_reactive()
      if (is.null(full_mat) || nrow(full_mat) == 0) {
        return(data.frame(Message = "No tissue-specific data available for these genes"))
      }
      
      # 2) subset to cart genes
      cart   <- as.character(cart_genes())
      expr_mat <- full_mat[rownames(full_mat) %in% cart, , drop = FALSE]
      if (nrow(expr_mat) == 0) {
        return(data.frame(Message = "No tissue-specific data available for these genes"))
      }
      
      # 3) coerce to plain numeric matrix
      expr_mat <- as.matrix(expr_mat)
      storage.mode(expr_mat) <- "numeric"
      
      # 4) compute Tau
      tau_df <- calcTau_custom(expr_mat)
      
      # 5) compute mean RPKM
      mean_expr <- rowMeans(expr_mat, na.rm = TRUE)
      
      # 6) find tissue with max expression
      top_tissue <- colnames(expr_mat)[max.col(expr_mat, ties.method = "first")]
      
      # 7) assemble results
      result_df <- data.frame(
        Gene      = rownames(expr_mat),
        Tau       = tau_df$tau[match(rownames(expr_mat), tau_df$gene)],
        Mean_RPKM = mean_expr,
        Top_Tissue= top_tissue,
        stringsAsFactors = FALSE
      )
      
      # 8) preserve cart order and fill NAs
      gene_df   <- data.frame(Gene = cart, stringsAsFactors = FALSE)
      joined_df <- dplyr::left_join(gene_df, result_df, by = "Gene")
      return(joined_df)
    }
    
    
    else if (mode == "ase") {
      # your original ASE branch unchanged
      ase_data <- gene_ase_df()
      if (is.null(ase_data) || nrow(ase_data) == 0) {
        return(data.frame(Message = "No ASE data available"))
      }
      ase_data$Gene <- toupper(ase_data$Gene)
      if ("mean.s" %in% colnames(ase_data) && !("score" %in% colnames(ase_data))) {
        ase_data <- ase_data %>% dplyr::rename(score = mean.s)
      }
      if (!("score" %in% colnames(ase_data))) {
        ase_data$score <- NA
      }
      gene_df <- data.frame(Gene = toupper(genes), stringsAsFactors = FALSE)
      result_df <- dplyr::left_join(gene_df, ase_data, by = "Gene")
      if (nrow(result_df) > 0 && (is.null(result_df$score) || length(result_df$score) == 0)) {
        result_df$score <- rep(NA, nrow(result_df))
      }
      result_df$ASE_Status <- ifelse(is.na(result_df$score), "Not ASE", "ASE")
      return(result_df)
      
    } else {
      return(data.frame(Message = "Mode not recognized"))
    }
  })
  
  
  
  # Note: We are not yet using mode_based_gene_table() in output$cart_gene_table. 
  # We will do that in the next chunks once we populate these placeholders 
  # with the actual data retrieval and formatting code.
  
  
  # Render the gene cart table with full statistics and multiple row selection enabled:
  output$cart_gene_table <- renderDT({
    # Assume mode_based_gene_table() returns your full gene statistics
    df <- mode_based_gene_table()
    
    # If the table is empty or has only a message, pass through unchanged.
    if (is.null(df) || (is.data.frame(df) && nrow(df) == 0) || 
        (is.character(df) && length(df) == 0)) {
      return(datatable(
        data.frame(Message = "No genes in cart"),
        options = list(dom = 't', paging = FALSE),
        rownames = FALSE
      ))
    }
    
    # If we're in DE mode and the fold change column exists,
    # rename it to reflect the current tissue pair.
    if (current_gene_mode() == "de" && "log2FoldChange" %in% names(df)) {
      if (order_val() == 0) {
        names(df)[names(df) == "log2FoldChange"] <- paste("log2FC (", 
                                                          input$tissue_select1, " up / ", 
                                                          input$tissue_select2, " down)", sep = "")
      } else {
        names(df)[names(df) == "log2FoldChange"] <- paste("log2FC (", 
                                                          input$tissue_select2, " up / ", 
                                                          input$tissue_select1, " down)", sep = "")
      }
    }
    
    datatable(df,
              selection = "multiple",
              options = list(pageLength = 10),
              rownames = FALSE)
  })
  
  
  # Observer to delete selected rows from the gene cart.
  observeEvent(input$delete_selected_genes, {
    selectedRows <- input$cart_gene_table_rows_selected
    if (!is.null(selectedRows) && length(selectedRows) > 0) {
      currentCart <- cart_genes()
      
      # If the cart is stored as a data frame:
      if (is.data.frame(currentCart)) {
        if (length(selectedRows) == nrow(currentCart)) {
          # Remove all rows; return an empty data frame with the same columns.
          updatedCart <- currentCart[FALSE, , drop = FALSE]
        } else {
          updatedCart <- currentCart[-selectedRows, , drop = FALSE]
        }
      } else if (is.character(currentCart)) {
        # If stored as a character vector:
        if (length(selectedRows) == length(currentCart)) {
          updatedCart <- character(0)
        } else {
          updatedCart <- currentCart[-selectedRows]
        }
      } else {
        # Fall back to unchanged if the type is unexpected.
        updatedCart <- currentCart
      }
      
      cart_genes(updatedCart)
      showNotification("Selected genes removed", type = "message")
    } else {
      showNotification("No genes selected", type = "warning")
    }
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
  
  # Helper function to get data for a specified mode
  get_data_for_mode <- function(mode, genes, unfiltered_data, gene_expression_data_tissue_reactive, gene_ase_df, calcTau_custom) {
    if (mode == "de") {
      de_data <- unfiltered_data()
      if (is.null(de_data) || nrow(de_data) == 0) return(NULL)
      gene_df <- data.frame(Gene = genes, stringsAsFactors = FALSE)
      if (!"Gene" %in% colnames(de_data)) {
        de_data <- tibble::rownames_to_column(de_data, var = "Gene")
      }
      merged_data <- dplyr::left_join(gene_df, de_data, by = "Gene")
      if (ncol(merged_data) == 1) return(NULL) # No data found
      return(merged_data)
      
    } else if (mode == "tissue") {
      
      # 1. Get the per‐tissue RPKM matrix
      expr_mat <- gene_expression_data_tissue_reactive()
      if (is.null(expr_mat) || nrow(expr_mat)==0) {
        return(data.frame(Message="No tissue‐specific data available for these genes"))
      }
      
      # 2. Compute Tau
      tau_vals <- calcTau_custom(expr_mat)$tau
      names(tau_vals) <- rownames(expr_mat)
      
      # 3. For each gene, pick the tissue with max expression
      max_tissue <- apply(expr_mat, 1, function(x) {
        colnames(expr_mat)[which.max(x)]
      })
      
      # 4. Assemble into one wide data.frame
      tissue_wide <- as.data.frame(expr_mat, stringsAsFactors=FALSE)
      tissue_wide$Gene    <- rownames(tissue_wide)
      tissue_wide$Tau     <- tau_vals[ tissue_wide$Gene ]
      tissue_wide$Tissue  <- max_tissue[ tissue_wide$Gene ]
      
      # reorder columns: Gene | Tau | Tissue | <tissue1> | <tissue2> | ...
      result_df <- tissue_wide %>%
        dplyr::select(Gene, Tau, Tissue, dplyr::everything())
      
      # make sure you preserve the cart order (optional)
      cart_df <- data.frame(Gene=genes, stringsAsFactors=FALSE)
      result_df <- dplyr::left_join(cart_df, result_df, by="Gene")
      
      return(result_df)
    }
    else if (mode == "ase") {
      ase_data <- gene_ase_df()
      if (is.null(ase_data) || nrow(ase_data) == 0) {
        return(data.frame(Message = "No ASE data available"))
      }
      # Ensure ase_data has a 'score' column
      if (!("score" %in% colnames(ase_data))) {
        ase_data$score <- NA
      }
      gene_df <- data.frame(Gene = toupper(cart_genes()), stringsAsFactors = FALSE)
      result_df <- dplyr::left_join(gene_df, ase_data, by = "Gene")
      
      if (length(result_df$score) == 0) {
        result_df$score <- rep(NA, nrow(result_df))
      }
      result_df$ASE_Status <- ifelse(is.na(result_df$score), "Not ASE", "ASE")
      return(result_df)
    }
    
    
    
    return(NULL)
  }
  
  # Adjusting the download_results handler
  output$download_results <- downloadHandler(
    filename = function() {
      paste0("gene_cart_results_", input$download_type, ".csv")
    },
    content = function(file) {
      mode_chosen <- input$download_type
      genes <- cart_genes()
      
      # If cart is empty
      if (length(genes) == 0) {
        write.csv(data.frame(Message = "No genes in cart"), file, row.names = FALSE)
        return()
      }
      
      if (mode_chosen == "de") {
        # — Differential Expression export (unchanged) —
        data <- get_data_for_mode("de", genes, unfiltered_data, gene_expression_data_tissue_reactive, gene_ase_df, calcTau_custom)
        if (is.null(data)) {
          write.csv(data.frame(Message = "No DE data available"), file, row.names = FALSE)
        } else {
          # Rename log2FoldChange column
          if ("log2FoldChange" %in% names(data)) {
            if (order_val() == 0) {
              names(data)[names(data) == "log2FoldChange"] <-
                paste0("log2FC (", input$tissue_select1, " up / ", input$tissue_select2, " down)")
            } else {
              names(data)[names(data) == "log2FoldChange"] <-
                paste0("log2FC (", input$tissue_select2, " up / ", input$tissue_select1, " down)")
            }
          }
          write.csv(data, file, row.names = FALSE)
        }
        
      } else if (mode_chosen == "tissue") {
        # — Improved Tissue-Specific export —
        # 1) Which tissues were selected in the heatmap?
        tissues <- isolate(input$tissue_choice)
        # 2) Pull the matching samples
        samples <- phenodata %>%
          filter(Tissue %in% tissues) %>%
          pull(Sample_name)
        # 3) Subset the RPKM matrix
        expr <- RPKM_data[genes, samples, drop = FALSE]
        # 4) Compute per-tissue averages
        expr_by_tissue <- sapply(tissues, function(t) {
          cols <- phenodata$Sample_name[phenodata$Tissue == t]
          rowMeans(expr[, cols, drop = FALSE], na.rm = TRUE)
        }, USE.NAMES = TRUE)
        colnames(expr_by_tissue) <- tissues
        # 5) Compute Tau
        tau_df <- calcTau_custom(expr_by_tissue)
        # 6) Build the output data.frame
        df_out <- data.frame(
          Gene   = rownames(expr_by_tissue),
          Tau    = tau_df$tau,
          Tissue = tissues[apply(expr_by_tissue, 1, which.max)],
          expr_by_tissue,
          stringsAsFactors = FALSE
        )
        # 7) Write CSV
        write.csv(df_out, file, row.names = FALSE)
        
      } else if (mode_chosen == "ase") {
        # — ASE export (unchanged) —
        data <- get_data_for_mode("ase", genes, unfiltered_data, gene_expression_data_tissue_reactive, gene_ase_df, calcTau_custom)
        if (is.null(data)) {
          write.csv(data.frame(Message = "No ASE data available"), file, row.names = FALSE)
        } else {
          write.csv(data, file, row.names = FALSE)
        }
        
      } else if (mode_chosen == "combined") {
        # — Combined export (unchanged) —
        de_data     <- get_data_for_mode("de",     genes, unfiltered_data, gene_expression_data_tissue_reactive, gene_ase_df, calcTau_custom)
        tissue_data <- get_data_for_mode("tissue", genes, unfiltered_data, gene_expression_data_tissue_reactive, gene_ase_df, calcTau_custom)
        ase_data    <- get_data_for_mode("ase",    genes, unfiltered_data, gene_expression_data_tissue_reactive, gene_ase_df, calcTau_custom)
        combined <- de_data
        if (!is.null(tissue_data)) combined <- full_join(combined, tissue_data, by = "Gene")
        if (!is.null(ase_data   )) combined <- full_join(combined, ase_data,    by = "Gene")
        if ("log2FoldChange" %in% names(combined)) {
          if (order_val() == 0) {
            names(combined)[names(combined) == "log2FoldChange"] <-
              paste0("log2FC (", input$tissue_select1, " up / ", input$tissue_select2, " down)")
          } else {
            names(combined)[names(combined) == "log2FoldChange"] <-
              paste0("log2FC (", input$tissue_select2, " up / ", input$tissue_select1, " down)")
          }
        }
        write.csv(combined, file, row.names = FALSE)
      }
    }
  )
  
  
  output$summary_stats_table <- renderTable({
    mode <- current_gene_mode()
    df <- mode_based_gene_table()
    
    # If the data frame only has a Message column, return NULL
    if (ncol(df) == 1 && "Message" %in% colnames(df)) {
      return(NULL)
    }
    
    if (mode == "de") {
      # Differential Expression Mode Stats
      # Assume df has columns: Gene, log2FoldChange, padj
      # Filter out rows with missing values as needed
      if (!"log2FoldChange" %in% colnames(df)) {
        return(NULL) # Can't compute DE stats if columns are missing
      }
      
      num_genes <- nrow(df)
      avg_lfc <- mean(df$log2FoldChange, na.rm = TRUE)
      median_lfc <- median(df$log2FoldChange, na.rm = TRUE)
      
      # Fraction passing FDR = 0.05 if padj is present
      if ("padj" %in% colnames(df)) {
        passing_fdr <- sum(df$padj < 0.05, na.rm = TRUE)
        fdr_fraction <- passing_fdr / num_genes
      } else {
        fdr_fraction <- NA
      }
      
      summary_stats <- data.frame(
        Metric = c("Number of Genes", "Average log2FC", "Median log2FC", "Fraction padj<0.05"),
        Value = c(num_genes, round(avg_lfc, 3), round(median_lfc, 3), round(fdr_fraction, 3))
      )
      
      return(t(summary_stats))
      
    } else if (mode == "tissue") {
      # Tissue-Specific Mode: Show Tau stats as before
      gene_expression_data_tissue <- gene_expression_data_tissue_reactive()
      if (is.null(gene_expression_data_tissue)) {
        return(NULL)
      }
      
      tau_scores <- calcTau_custom(gene_expression_data_tissue)
      tau_values <- tau_scores$tau
      
      summary_stats <- data.frame(
        Metric = c("Number of Genes", "Average Tau", "Median Tau", "Tau Variance", "Min Tau", "Max Tau"),
        Value = c(
          nrow(gene_expression_data_tissue),
          round(mean(tau_values, na.rm = TRUE), 4),
          round(median(tau_values, na.rm = TRUE), 4),
          round(var(tau_values, na.rm = TRUE), 4),
          round(min(tau_values, na.rm = TRUE), 4),
          round(max(tau_values, na.rm = TRUE), 4)
        )
      )
      
      return(t(summary_stats))
      
    } else if (mode == "ase") {
      # ASE Mode: Show stats about ASE vs Non-ASE genes
      # df should have ASE_Status and maybe score
      if (!"ASE_Status" %in% colnames(df)) {
        return(NULL)
      }
      
      num_genes <- nrow(df)
      ase_count <- sum(df$ASE_Status == "ASE", na.rm = TRUE)
      not_ase_count <- sum(df$ASE_Status == "Not ASE", na.rm = TRUE)
      
      # If score column exists, compute average ASE score among ASE genes
      if ("score" %in% colnames(df)) {
        ase_scores <- df$score[df$ASE_Status == "ASE" & !is.na(df$score)]
        avg_ase_score <- if(length(ase_scores) > 0) mean(ase_scores, na.rm = TRUE) else NA
      } else {
        avg_ase_score <- NA
      }
      
      summary_stats <- data.frame(
        Metric = c("Number of Genes", "ASE Genes", "Not ASE Genes", "Average ASE Score (ASE genes)"),
        Value = c(num_genes, ase_count, not_ase_count, round(avg_ase_score, 3))
      )
      
      return(t(summary_stats))
      
    } else {
      # If mode not recognized or no data
      return(NULL)
    }
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
    mode <- current_gene_mode()
    genes <- cart_genes()
    
    if (length(genes) == 0) {
      # No genes in cart
      plot.new()
      text(0.5, 0.5, "No genes selected.", cex = 1)
      return(NULL)
    }
    
    if (mode == "de") {
      # DE mode: show ASE distribution across the two selected tissues
      req(input$tissue_select1, input$tissue_select2)
      
      # Fetch ASE data
      ase_data <- gene_ase_df()
      if (is.null(ase_data) || nrow(ase_data) == 0) {
        plot.new()
        text(0.5, 0.5, "No ASE data available.", cex = 1)
        return(NULL)
      }
      
      Tissue1 <- input$tissue_select1
      Tissue2 <- input$tissue_select2
      
      # Find ASE genes in each tissue (score non-NA indicates ASE)
      ase_genes_t1 <- unique(ase_data %>% 
                               filter(Gene %in% genes, Tissue == Tissue1, !is.na(score)) %>% 
                               pull(Gene))
      ase_genes_t2 <- unique(ase_data %>% 
                               filter(Gene %in% genes, Tissue == Tissue2, !is.na(score)) %>% 
                               pull(Gene))
      
      # Genes ASE in both tissues
      both <- intersect(ase_genes_t1, ase_genes_t2)
      # Genes ASE only in Tissue 1
      t1_only <- setdiff(ase_genes_t1, both)
      # Genes ASE only in Tissue 2
      t2_only <- setdiff(ase_genes_t2, both)
      
      bar_data <- data.frame(
        Category = c(paste("ASE in", Tissue1), paste("ASE in", Tissue2), "ASE in Both"),
        Count = c(length(t1_only), length(t2_only), length(both))
      )
      
      barplot(bar_data$Count, names.arg = bar_data$Category, col = c("blue", "orange", "green"),
              main = "ASE Genes Distribution", ylab = "Count")
      
    } else {
      # Non-DE mode: original behavior
      allele_specific_genes_list <- allele_specific_genes_reactive()
      gene_expression_data_tissue <- gene_expression_data_tissue_reactive()
      if (is.null(gene_expression_data_tissue)) {
        return(NULL)
      }
      
      allele_specific_count <- sum(toupper(rownames(gene_expression_data_tissue)) %in% toupper(allele_specific_genes_list))
      non_allele_specific_count <- nrow(gene_expression_data_tissue) - allele_specific_count
      
      bar_data <- data.frame(
        Category = c("Allele-Specific", "Non-Allele-Specific"),
        Count = c(allele_specific_count, non_allele_specific_count)
      )
      
      barplot(bar_data$Count, names.arg = bar_data$Category, col = c("blue", "gray"),
              main = "Number of Allele-Specific Genes", ylab = "Count")
    }
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
  
  output$mode_specific_plot <- renderPlot({
    mode <- current_gene_mode()
    genes <- cart_genes()
    
    if (mode == "de") {
      data <- unfiltered_data()
      if (is.null(data) || nrow(data) == 0) {
        plot.new()
        text(0.5, 0.5, "No DE Data", cex = 1)
        return()
      }
      if (!"Gene" %in% colnames(data)) {
        data <- tibble::rownames_to_column(data, var = "Gene")
      }
      
      # Replace zero padj to avoid -Inf and calculate -log10(padj)
      data$padj <- ifelse(data$padj == 0, 1e-300, data$padj)
      data$negLogP <- -log10(data$padj)
      
      # Build comparison label based on order
      comp_label <- if (order_val() == 0) {
        paste(input$tissue_select1, "vs", input$tissue_select2)
      } else {
        paste(input$tissue_select2, "vs", input$tissue_select1)
      }
      
      # Plot the volcano plot with base R plotting:
      plot(data$log2FoldChange, data$negLogP, pch = 19, cex = 0.5,
           xlab = paste0("log2FC (", comp_label, ")"),
           ylab = "-log10 P", main = "Volcano Plot",
           cex.main = 1, cex.lab = 0.9, cex.axis = 0.9)
      
      # Highlight genes that are in the gene cart
      if (!is.null(genes) && length(genes) > 0) {
        # If the cart is a data frame with a "Source Comparison" column, filter by the current comparison.
        if (is.data.frame(genes) && "Source Comparison" %in% names(genes)) {
          current_genes <- genes$Gene[genes$`Source Comparison` == comp_label]
        } else {
          current_genes <- genes
        }
        
        highlighted <- data[data$Gene %in% current_genes, ]
        if (nrow(highlighted) > 0) {
          # Simply replot the highlighted points at the same coordinates in red.
          points(highlighted$log2FoldChange, highlighted$negLogP, pch = 19, col = "red")
          # Optionally, add text labels at the exact coordinates.
          text(highlighted$log2FoldChange, highlighted$negLogP, labels = highlighted$Gene,
               pos = 3, cex = 0.7, col = "red")
        }
      }
      
      # Add tissue labels to indicate upregulation direction.
      if (!is.null(input$tissue_select1) && !is.null(input$tissue_select2)) {
        x_range <- range(data$log2FoldChange, na.rm = TRUE)
        y_range <- range(data$negLogP, na.rm = TRUE)
        # Left side: negative values indicate upregulation in tissue_select2.
        text(x = x_range[1], y = y_range[2] * 0.8, 
             labels = paste("Up in", input$tissue_select2), pos = 4, cex = 0.8, col = "blue")
        # Right side: positive values indicate upregulation in tissue_select1.
        text(x = x_range[2], y = y_range[2] * 0.8, 
             labels = paste("Up in", input$tissue_select1), pos = 2, cex = 0.8, col = "blue")
      }
      
    } else if (mode == "tissue") {
      # Tissue mode code remains unchanged.
      table_data <- mode_based_gene_table()
      if ("Mean_RPKM" %in% names(table_data) && !all(is.na(table_data$Mean_RPKM))) {
        table_data <- table_data[order(table_data$Mean_RPKM, decreasing = TRUE), ]
        par(mar = c(7, 4, 2, 1))
        barplot(table_data$Mean_RPKM, names.arg = table_data$Gene, las = 2, cex.names = 0.7,
                main = "Avg Expression", ylab = "Mean RPKM", cex.main = 0.9, cex.lab = 0.8, cex.axis = 0.8)
      } else {
        plot.new()
        text(0.5, 0.5, "No Tissue Data", cex = 1)
      }
      
    } else if (mode == "ase") {
      table_data <- mode_based_gene_table()
      if ("score" %in% names(table_data)) {
        ase_genes <- table_data[table_data$ASE_Status == "ASE" & !is.na(table_data$score), ]
        not_ase_genes <- table_data[table_data$ASE_Status == "Not ASE", ]
        if (nrow(ase_genes) > 0) {
          boxplot(ase_genes$score, main = "ASE Scores", ylab = "score",
                  cex.main = 0.9, cex.lab = 0.8, cex.axis = 0.8)
          mtext(paste("ASE:", nrow(ase_genes), "Not ASE:", nrow(not_ase_genes)),
                side = 3, cex = 0.8)
        } else {
          plot.new()
          text(0.5, 0.5, "No ASE Genes", cex = 1)
        }
      } else {
        plot.new()
        text(0.5, 0.5, "No ASE Data", cex = 1)
      }
      
    } else {
      plot.new()
      text(0.5, 0.5, "Select a mode or add genes", cex = 1)
    }
  })
  
}

shinyApp(ui = ui, server = server)
