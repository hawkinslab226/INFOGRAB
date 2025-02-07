#######################################################
## app.R -- Single-file Shiny app with:
##  1) "All Tissues" in the checkbox
##  2) Single-track vs multiple-track radio button
##  3) Tissue-based Gene ASE loading
#######################################################

library(shiny)
library(igvShiny)
library(rtracklayer)
library(dplyr)
library(readr)
library(VariantAnnotation)  # readVcf for VCF

#setwd("~/Documents/Internship/InfoGrab")

ui <- fluidPage(
  
  # Single "Genome Viewer" tabPanel
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
      
      ## ROW 1: Search region + track input type
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
  )
)


server <- function(input, output, session) {
  
  ###########################################
  # PART A: Tissue-Gene info from the TXT
  ###########################################
  
  gene_ase_df <- reactiveVal(data.frame())
  
  observe({
    info_path <- "browser_data_for_app/all.geneAse.info827.txt"
    if (!file.exists(info_path)) {
      showNotification("Cannot find all.geneAse.info827.txt. No Tissue list will be available.", type="error")
      return(NULL)
    }
    
    df_info <- read.delim(info_path, header=TRUE, stringsAsFactors=FALSE)
    if (!all(c("Gene", "Tissue") %in% colnames(df_info))) {
      showNotification("Missing 'Gene' or 'Tissue' columns in the TXT file.", type="error")
      return(NULL)
    }
    
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
  # GRAB ALL ASE GENES (Example)
  ###########################################
  observeEvent(input$grab_all_ase_genes, {
    bed_path <- "browser_data_for_app/gene_ase_converted_cleaned_with_strand.bed"
    if (!file.exists(bed_path)) {
      showNotification("Gene ASE bed file not found for 'Grab All ASE Genes'.", type="error")
      return(NULL)
    }
    
    withProgress(message='Grabbing All ASE Genes...', value=0, {
      tryCatch({
        bed_data <- rtracklayer::import(bed_path, format="bed")
        bed_df   <- as.data.frame(bed_data)
        
        # If bed uses "name" as the gene
        all_ase_genes <- character(0)
        if ("name" %in% colnames(bed_df)) {
          all_ase_genes <- toupper(bed_df$name)
        }
        
        showNotification(
          paste("All ASE Genes grabbed:", length(all_ase_genes)),
          type="message"
        )
        
        incProgress(1)
      }, error=function(e) {
        showNotification(
          paste("Error grabbing ASE genes:", e$message),
          type="error"
        )
      })
    })
  })
}

shinyApp(ui=ui, server=server)
