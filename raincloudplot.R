library(shiny)
library(ggplot2)
library(gghalves)
library(dplyr)
library(ggiraph)

# Raincloud plot function
make_raincloud_plot <- function(data, show_mean_line, show_points, facet_by_tissue = FALSE) {
  data <- dplyr::distinct(data)
  
  gg <- ggplot2::ggplot(data, ggplot2::aes(y = score, x = if (facet_by_tissue) Tissue else 1)) +
    {if (show_points) ggiraph::geom_point_interactive(
      ggplot2::aes(
        x = 0.87,
        color = Gene,
        tooltip = paste("Gene:", Gene, "<br>Score:", round(score, 4)),
        data_id = Gene
      ),
      position = ggplot2::position_jitter(width = 0.06, height = 0),
      shape = 16,
      alpha = 0.7,
      size = 2
    ) else NULL} +
    {if (nrow(data) > 2) gghalves::geom_half_violin(
      side = "r",
      trim = FALSE,
      fill = "grey",
      alpha = 0.6,
      color = "black"
    ) else NULL} +
    ggplot2::geom_boxplot(
      width = 0.1,
      outlier.shape = NA,
      alpha = 0.6,
      color = "black"
    ) +
    {if (show_mean_line) ggplot2::stat_summary(fun = mean, geom = "crossbar", width = 0.075, color = "red") else NULL} +
    ggplot2::coord_flip() +
    ggplot2::theme_minimal(base_size = 14) +
    ggplot2::scale_x_discrete(breaks = NULL) +
    ggplot2::scale_color_viridis_d(option = "D", end = 0.9) +
    ggplot2::labs(
      x = NULL,
      y = "Score",
      title = "Raincloud Plot"
    ) +
    ggplot2::theme(
      axis.ticks.y = ggplot2::element_blank(),
      plot.title = ggplot2::element_text(hjust = 0.5),
      legend.position = "none"
    )
  
  if (facet_by_tissue) {
    gg <- gg + ggplot2::facet_wrap(~Tissue, scales = "free", ncol = 2)
  }
  
  ggiraph::girafe(ggobj = gg, width_svg = 8, height_svg = 6)
}

make_static_raincloud_plot <- function(data, show_mean_line, show_points, facet_by_tissue = FALSE) {
  data <- dplyr::distinct(data)
  
  gg <- ggplot2::ggplot(data, ggplot2::aes(y = score, x = if (facet_by_tissue) Tissue else 1)) +
    {if (show_points) ggplot2::geom_point(
      ggplot2::aes(
        x = 0.87,
        color = Gene
      ),
      position = ggplot2::position_jitter(width = 0.06, height = 0),
      shape = 16,
      alpha = 0.7,
      size = 2
    ) else NULL} +
    {if (nrow(data) > 2) gghalves::geom_half_violin(
      side = "r",
      trim = FALSE,
      fill = "grey",
      alpha = 0.6,
      color = "black"
    ) else NULL} +
    ggplot2::geom_boxplot(
      width = 0.1,
      outlier.shape = NA,
      alpha = 0.6,
      color = "black"
    ) +
    {if (show_mean_line) ggplot2::stat_summary(fun = mean, geom = "crossbar", width = 0.075, color = "red") else NULL} +
    ggplot2::coord_flip() +
    ggplot2::theme_minimal(base_size = 14) +
    ggplot2::scale_x_discrete(breaks = NULL) +
    ggplot2::scale_color_viridis_d(option = "D", end = 0.9) +
    ggplot2::labs(
      x = NULL,
      y = "Score",
      title = "Raincloud Plot"
    ) +
    ggplot2::theme(
      axis.ticks.y = ggplot2::element_blank(),
      plot.title = ggplot2::element_text(hjust = 0.5),
      legend.position = "none"
    )
  
  if (facet_by_tissue) {
    gg <- gg + ggplot2::facet_wrap(~Tissue, scales = "free", ncol = 2)
  }
  
  return(gg)
}

# UI
ui <- shiny::navbarPage(
  "INFO GRAB",
  shiny::tabPanel(
    "Allele Expression",
    shiny::sidebarLayout(
      shiny::sidebarPanel(
        shiny::numericInput(
          "top_n_raincloud",
          "Number of Top Genes (Raincloud):",
          value = 10,
          min = 1,
          max = 10000,
          step = 1
        ),
        shiny::sliderInput("pval_cutoff", "P-value Cutoff (FDR):", min = 0, max = 0.05, value = 0.05, step = 0.005),
        shiny::selectInput("tissue_select_raincloud", "Select Tissues:", choices = NULL, selected = NULL, multiple = TRUE),
        shiny::conditionalPanel(
          condition = "input.tissue_select_raincloud.includes('All Tissues')",
          shiny::selectizeInput(
            "gene_search_select",
            "Search for Genes:",
            choices = NULL,
            selected = NULL,
            multiple = TRUE,
            options = list(placeholder = "Type to search genes...")
          )
        ),
        shiny::actionButton("random_genes", "10 Random Genes"),
        shiny::actionButton("reset_top_genes", "Reset to Top Genes"), # New Button
        shiny::checkboxInput("show_points", "Show Individual Points", TRUE),
        shiny::checkboxInput("show_mean_line", "Show Mean Line", TRUE),
        shiny::downloadButton("download_raincloud_plot", "Download PNG"),
        shiny::actionButton("grab_all_ase_genes", "Grab All ASE Genes")
      ),
      shiny::mainPanel(
        shiny::h3("Raincloud Plot by Gene"),
        ggiraph::girafeOutput("raincloud_plot", width = "100%", height = "1200px")
      )
    )
  )
)

# Server
server <- function(input, output, session) {
  gene_ase_path <- "browser_data_for_app/all.geneAse.info827.txt"
  
  gene_ase_df <- shiny::reactiveVal(data.frame(Gene = character(), n.vars = numeric(), score = numeric(), fdr = numeric(), Tissue = character()))
  selected_genes <- shiny::reactiveVal(NULL)
  
  shiny::observe({
    if (file.exists(gene_ase_path)) {
      shiny::withProgress(message = 'Loading Gene ASE Data...', value = 0, {
        tryCatch({
          df <- utils::read.delim(gene_ase_path, header = TRUE, stringsAsFactors = FALSE, sep = "\t")
          df <- df %>% dplyr::rename(score = mean.s) %>% dplyr::select(Gene, n.vars, score, fdr, Tissue)
          gene_ase_df(df)
          tissues <- sort(unique(df$Tissue))
          shiny::updateSelectInput(session, "tissue_select_raincloud", choices = c("All Tissues", tissues), selected = "All Tissues")
          shiny::incProgress(1)
        }, error = function(e) {
        })
      })
    } else {
    }
  })
  
  shiny::observeEvent(input$tissue_select_raincloud, {
    df <- gene_ase_df()
    shiny::req(df)
    
    if ("All Tissues" %in% input$tissue_select_raincloud) {
      shiny::updateSelectizeInput(session, "gene_search_select", choices = sort(unique(df$Gene)), selected = NULL)
    } else {
      shiny::updateSelectizeInput(session, "gene_search_select", choices = NULL, selected = NULL)
    }
  })
  
  shiny::observeEvent(input$gene_search_select, {
    if (!is.null(input$gene_search_select) && length(input$gene_search_select) > 0) {
      # Clear random genes when a gene is searched
      selected_genes(NULL)
      shiny::showNotification("Random genes list cleared. Showing selected genes.", type = "message")
    }
  })
  
  shiny::observeEvent(input$reset_top_genes, {
    # Clear random genes
    selected_genes(NULL)
    
    # Clear the search bar
    shiny::updateSelectizeInput(session, "gene_search_select", selected = NULL)
    
    # Notify the user
    shiny::showNotification("Reset to Top Genes. Showing the top selected genes.", type = "message")
  })
  
  
  
  shiny::observeEvent(input$random_genes, {
    df <- gene_ase_df()
    shiny::req(df)
    
    # Clear searched genes in the UI and backend
    shiny::updateSelectizeInput(session, "gene_search_select", selected = NULL)
    session$sendInputMessage("gene_search_select", list(value = ""))
    
    if (is.null(input$tissue_select_raincloud) || length(input$tissue_select_raincloud) == 0) {
      return()
    }
    
    random_genes_list <- list()
    
    if ("All Tissues" %in% input$tissue_select_raincloud) {
      all_genes <- unique(df$Gene)
      if (length(all_genes) < 10) {
        shiny::showNotification("Not enough genes in the dataset for random selection.", type = "warning")
        selected_genes(NULL)
        return()
      }
      random_genes <- sample(all_genes, size = 10, replace = FALSE)
      random_genes_list[["All Tissues"]] <- random_genes
    } else {
      for (tissue in input$tissue_select_raincloud) {
        tissue_df <- df %>% dplyr::filter(Tissue == tissue)
        if (dplyr::n_distinct(tissue_df$Gene) < 10) {
          next
        }
        random_genes <- sample(unique(tissue_df$Gene), size = 10, replace = FALSE)
        random_genes_list[[tissue]] <- random_genes
      }
    }
    
    if (length(random_genes_list) == 0) {
      selected_genes(NULL)
      return()
    }
    
    selected_genes(random_genes_list)
    shiny::showNotification("Random genes selected successfully!", type = "message")
  })
  
  
  
  filtered_data <- shiny::reactive({
    df <- gene_ase_df()
    shiny::req(df)
    
    # Filter by FDR cutoff
    df <- df %>% dplyr::filter(fdr <= input$pval_cutoff)
    
    # If 'All Tissues' is selected and specific genes are searched
    if ("All Tissues" %in% input$tissue_select_raincloud && !is.null(input$gene_search_select)) {
      df <- df %>% dplyr::filter(Gene %in% input$gene_search_select)
    }
    
    # Handle random genes selection
    random_genes_list <- selected_genes()
    if (!is.null(random_genes_list)) {
      if ("All Tissues" %in% names(random_genes_list)) {
        df <- df %>% dplyr::filter(Gene %in% random_genes_list[["All Tissues"]])
      } else {
        df <- df %>% dplyr::filter(Tissue %in% names(random_genes_list) & Gene %in% unlist(random_genes_list))
      }
    } else {
      # Filter for top N genes
      if ("All Tissues" %in% input$tissue_select_raincloud) {
        top_genes <- df %>%
          dplyr::group_by(Gene) %>%
          dplyr::summarise(mean_score = mean(score, na.rm = TRUE), .groups = "drop") %>%
          dplyr::arrange(dplyr::desc(mean_score)) %>%
          dplyr::slice_head(n = input$top_n_raincloud)
        df <- df %>% dplyr::filter(Gene %in% top_genes$Gene)
      } else {
        df <- df %>% dplyr::filter(Tissue %in% input$tissue_select_raincloud)
        top_genes <- df %>%
          dplyr::group_by(Tissue, Gene) %>%
          dplyr::summarise(mean_score = mean(score, na.rm = TRUE), .groups = "drop") %>%
          dplyr::arrange(Tissue, dplyr::desc(mean_score)) %>%
          dplyr::group_by(Tissue) %>%
          dplyr::slice_head(n = input$top_n_raincloud)
        df <- df %>% dplyr::filter(Gene %in% top_genes$Gene)
      }
    }
    
    if (nrow(df) == 0) {
      shiny::req(FALSE)
    }
    df
  })
  
  
  output$raincloud_plot <- ggiraph::renderGirafe({
    data <- filtered_data()
    shiny::req(data)
    facet_by_tissue <- !("All Tissues" %in% input$tissue_select_raincloud)
    make_raincloud_plot(data, input$show_mean_line, input$show_points, facet_by_tissue)
  })
  
  output$download_raincloud_plot <- shiny::downloadHandler(
    filename = function() { "raincloud_plot.png" },
    content = function(file) {
      data <- filtered_data()
      shiny::req(data)
      facet_by_tissue <- !("All Tissues" %in% input$tissue_select_raincloud)
      plot <- make_static_raincloud_plot(data, input$show_mean_line, input$show_points, facet_by_tissue)
      ggplot2::ggsave(file, plot, width = 10, height = 6)
    }
  )
  
}

shiny::shinyApp(ui = ui, server = server)
