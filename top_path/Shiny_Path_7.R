



library(shiny)
library(DT)
library(monocle3)
library(Seurat)
library(DEsingle)
library(ggplot2)
library(ggforce)
library(png)
library(grid)
library(talklr)
library(stats)
library(edgeR)
library(parallel)
library(MAST)
library(shinycssloaders)
library(AUCell)
library(digest)
library(Matrix)
library(igraph)
# install.packages("ggrepel")  # if needed
library(ggrepel)

library(text2vec)    # fallback uses GloVe if word2vec pkg missing
# optional: if installed, we'll prefer it
# install.packages("word2vec")  # only if you don't have it
# install.packages("text2vec")
library(treemap)
library(RColorBrewer)

# BiocManager::install("AUCell")
# devtools::install_github("aertslab/SCENIC") 


# Pathways file
#pathway_choices <- readLines("top_path/pathways.txt")
#pathway_choices_with_all <- c("All", pathway_choices)

meta.file <- tempdir()

CellChatDB <- readRDS("database/CellChatDB.rds")
receptor_ligand <- read.csv("database/receptor_ligand.csv")
#receptor_ligand <- receptor_ligand[ , -14]


# Source your function files
source("top_path/top_path_function_2.0.R")
source("top_path/Run_scMLnet.R")
source("top_path/Draw_MLnet.R")
source("top_path/plot_lr_wiring.R")
source("top_path/convert_h5ad_to_cds.R")
source("top_path/util.R")




# ---- EMPTY-SAFE replacement for .mode_weighted_magnitude ----
mode_weighted_magnitude_tbl  <- function(df, key_cols, wcol, scol) {
  df <- as.data.frame(df)
  need <- c(key_cols, wcol, scol)
  miss <- setdiff(need, names(df))
  if (length(miss)) stop("Missing columns: ", paste(miss, collapse=", "))
  
  # keep complete rows on weight & score
  df <- df[complete.cases(df[, c(wcol, scol), drop = FALSE]), , drop = FALSE]
  
  # if nothing to aggregate, return a correctly-shaped empty frame
  if (!nrow(df)) {
    out <- data.frame(
      setNames(replicate(length(key_cols), character(0), simplify = FALSE), key_cols),
      aggregate = numeric(0),
      score.perc = numeric(0),
      stringsAsFactors = FALSE
    )
    return(out)
  }
  
  # weighted mean of 'scol' by 'wcol' over key_cols
  ws <- df[[wcol]]
  sc <- df[[scol]]
  grp <- df[key_cols]
  
  num <- stats::aggregate(sc * ws, by = grp, FUN = sum)
  den <- stats::aggregate(ws,       by = grp, FUN = sum)
  names(num)[ncol(num)] <- "aggregate"
  names(den)[ncol(den)] <- "w"
  
  out <- merge(num, den, by = key_cols, all = TRUE)
  out$aggregate <- ifelse(is.na(out$w) | out$w == 0, NA_real_, out$aggregate / out$w)
  
  mx <- suppressWarnings(max(out$aggregate, na.rm = TRUE))
  out$score.perc <- if (is.finite(mx) && mx > 0) 100 * out$aggregate / mx else NA_real_
  
  out$w <- NULL
  out[order(-out$aggregate), , drop = FALSE]
}


library(base64enc)
logo_uri <- dataURI(file = "www/toppath_logo.png", mime = "image/png")


ui <- fluidPage(
  tags$head(tags$style(HTML("
  .tab-content { padding-top: 8px; }
  .dataTables_wrapper { overflow-x: auto; }
  .status-msg { margin-bottom: 14px; display: block; }   /* add spacing under status line */
"))),
  # div(style="text-align:center;", #10px 0 6px
  #     img(src = logo_uri, height = 80, alt = "TopPath 2.0")),
  
  #light blue center
  # tags$head(tags$style(HTML("
  #   .brand-hero {
  #     background: linear-gradient(135deg,#f7fbff 0%, #eef5ff 100%);
  #     border-radius: 14px;
  #     box-shadow: 0 8px 22px rgba(0,0,0,0.06);
  #     padding: 0px 0px 0px;
  #     margin: 0px 0px 0px;
  #     text-align: center;
  #   }
  #   .brand-hero img {
  #     max-height: 200px;
  #     height: auto;
  #     width: auto;
  #     vertical-align: middle;
  #     filter: drop-shadow(0 2px 4px rgba(0,0,0,.08));
  #   }
  #   .brand-title {
  #     margin-top: 8px;
  #     font-weight: 700;
  #     font-size: 22px;
  #     color: #1f2d3d;
  #     letter-spacing: .2px;
  #   }
  #   .brand-sub {
  #     margin-top: 2px;
  #     color: #5c6773;
  #     font-size: 13px;
  #   }
  # "))),
  # div(class = "brand-hero",
  #     img(src = logo_uri, alt = "TopPath 2.0")
  #     #div(class = "brand-title", "TopPath 2.0"),
  #     #div(class = "brand-sub", "Top Pathway 2.0 Analysis:")
  # ),
  
  
  # thin side
  tags$head(tags$style(HTML("
    .brand-bar {
      position: sticky; top: 0; z-index: 1000;
      background: #ffffffcc;
      -webkit-backdrop-filter: blur(6px);
      backdrop-filter: blur(6px);
      border-bottom: 1px solid #eef1f5;
      padding: 8px 14px;
      display: flex; align-items: center; gap: 12px;
    }
    .brand-bar img { max-height: 80px; filter: drop-shadow(0 1px 2px rgba(0,0,0,.05)); }
    .brand-name { font-weight: 700; font-size: 18px; color: #1f2d3d; }
    .brand-pill {
      margin-left: auto;
      background:#eef5ff; color:#3a6ee8;
      padding: 4px 10px; border-radius: 9999px; font-size: 12px; font-weight:600;
    }
  "))),
  div(class = "brand-bar",
      img(src = logo_uri, alt = "TopPath 2.0"),
      span(class = "brand-name", "Top Pathway Analysis:"),
      span(class = "brand-pill", "v2.0 | Ruohola-Baker Lab")
      #span(class = "brand-pill", "Ruohola-Baker Lab")
      
  ),
  
  # titlePanel("Top Pathway 2.0 Analysis:"),
  sidebarLayout(
    sidebarPanel(
      fileInput("cdsfile", "Upload Monocle CDS (.rds) or Scanpy (.h5ad)", accept = c(".rds", ".h5ad")),
      uiOutput("filetype_status"),
      selectizeInput("clustering_col", "Choose clustering column", choices = NULL, selected = NULL, options = list(placeholder = "Select clustering column...")),
      selectizeInput("receiver_cell", "Starting Cell", choices = NULL, selected = NULL, options = list(placeholder = "Select receiver cell...")),
      selectizeInput("diff_target_cell", "Diff Target Cell", choices = NULL, selected = NULL, options = list(placeholder = "Select diff target cell...")),
      selectizeInput("sender_cells", "Sender Cells (multiple allowed)", choices = NULL, selected = NULL, multiple = TRUE, options = list(placeholder = "Select one or more sender cells...")),
      # ---- Pathway filters + selector ---------------------------------------------
      checkboxGroupInput(
        "lr_classes", "Interaction types",
        choices  = c("Secreted Signaling",
                     "ECM-Receptor",
                     "Cell-Cell Contact",
                     "Non-protein Signaling"),
        selected = "Secreted Signaling"   # <- default
      ),
      
      selectizeInput(
        "path_included", "Pathways Included",
        choices  = NULL,                   # filled by server
        selected = "All",
        multiple = TRUE,
        options  = list(placeholder = "Select pathways or All...",
                        plugins = list("remove_button"))
      )
      
      ,
      selectInput("de_method", "Differential Expression Method",
                  choices = c("DEsingle (cell-level)" = "DEsingle", "edgeR (Pseudo-bulk)" = "edgeR", "MAST (cell-level)" = "MAST"),
                  selected = "MAST"
      ),
      ## ---- SIDEBAR: learning + SCENIC (replace the three checkboxes + sliders) ----
      checkboxInput("simple_learning", "Simple Learning (Embeding and Clustering)", value = FALSE),
      # --- Simple Learning options (shown only when Simple Learning is ON) ---
      conditionalPanel(
        condition = "input.simple_learning",
        div(style = "margin-left: 8px; border-left: 3px solid #eee; padding-left: 10px;",
            # --- Simple-learning pruning controls ---
            sliderInput("k_per_cluster", "Max pathways per cluster (pruned)",
                        min = 1, max = 5, value = 2, step = 1),
            conditionalPanel(
              condition = "false",
              selectInput("rank_method", "Prune ranking",
                          choices = c("Pareto front" = "pareto",
                                      "Weighted (Active/Unique)" = "weighted"),
                          selected = "pareto")
            )
            ,
            
        )
      ),
      
      # Merge old 'deep_learning' and 'true_gnn' into one:
      checkboxInput("deep_learning", "Deep Learning (Graph Nerual Network)", value = FALSE),
      
      # SCENIC toggle, with its params hidden unless enabled
      checkboxInput("use_scenic", "Use SCENIC TF activity weighting (AUCell)", FALSE),
      
      # Show the sliders only when SCENIC is ON
      conditionalPanel(
        condition = "input.use_scenic",
        div(style = "margin-left: 8px; border-left: 3px solid #eee; padding-left: 10px;",
            sliderInput("scenic_w_range", "Weight cap [min, max]", min = 0.25, max = 4, value = c(0.5, 2), step = 0.25),
            sliderInput("scenic_beta", "Sensitivity (β)", min = 0.5, max = 5, value = 2, step = 0.1)
        )
      ),
      
      
      div(
        style = "text-align:center;margin-top:15px;margin-bottom:5px;",
        actionButton("run", "Run Analysis"),
        actionButton("show_example", "Example"),
        actionButton("cancel_analysis", "Cancel", style = "color:white;background-color:#b20000;")
      )
    ),
    mainPanel(
      uiOutput("abbrev_legend"),
      plotOutput("selection_plot", height = "600px"),
      uiOutput("analysis_status"),
      
      # One tabset for everything below
      tabsetPanel(
        id = "main_tabs",
        tabPanel(
          "Treemaps",
          tabsetPanel(
            id = "treemap_tabs",   # <— add this
            tabPanel("Pathway",           imageOutput("pathwayPlot",      height = "750px")),
            tabPanel("Ligand by Pathway", imageOutput("ligandPlot",       height = "750px")),
            tabPanel("Ligand Only",       imageOutput("ligandOnlyPlot",   height = "750px")),
            tabPanel("Receptor Only",     imageOutput("receptorOnlyPlot", height = "750px")),
            #tabPanel("Interaction by Pathway", imageOutput("interactionPlot",     height = "750px")),
            tabPanel("Interaction Only",       imageOutput("interactionOnlyPlot", height = "750px"))
            
          )
          
        ),
        tabPanel(
          "Tables",
          div(style = "margin: 12px 0;",
              selectInput(
                "table_choice", "Result Table",
                choices = c("Pathway (summary + p)" = "pathway_n",
                            "Ligand (summary + p)"  = "lig_rank_all",
                            "Receptor (summary + p)"= "rec_rank_all",
                            "Interaction (summary+p)" = "int_rank_all"),   # <--- NEW
                selected = "pathway_n"
              )
              
          ),
          DTOutput("result_table")
        ),
        tabPanel(
          "Simple Learning",
          tabsetPanel(
            tabPanel("Embedding", imageOutput("simple_scatter", height = "520px")),
            tabPanel("Clustered Treemap",          imageOutput("simple_treemap", height = "740px")),
            tabPanel("Pruned Treemap",             imageOutput("simple_pruned_treemap", height = "740px")),
            tabPanel("Explanations", uiOutput("explain_gallery")),
            tabPanel("Active/Unique",
                     fluidRow(
                       column(12, imageOutput("simple_active_unique_img",   height = "600px")
                       )
                     )
            ),
            tabPanel("Pruned LR Treemap", imageOutput("simple_lr_pruned_treemap", height = "740px"))
          )
          
        ),
        tabPanel(
          "Deep Learning",
          tabsetPanel(
            tabPanel("GNN Pathways",       imageOutput("true_gnn_treemap", height = "740px")),
            tabPanel("Embedding",            imageOutput("true_simple_scatter", height = "520px")),
            tabPanel("Clustered Treemap",  imageOutput("true_simple_treemap", height = "740px")),
            tabPanel("Pruned Treemap",     imageOutput("true_simple_pruned_treemap", height = "740px")),
            tabPanel("Explanations",       uiOutput("true_explain_gallery", height = "100%")),
            tabPanel(
              "Active/Unique",
              fluidRow(
                column(12,
                       imageOutput("deep_active_unique_img",   height = "600px"),
                       br(),
                       imageOutput("deep_activity_unique_img", height = "600px")
                )
              )
            )
            
            
          )
        ),
        
        tabPanel(
          "talklr",
          fluidPage(
            fluidRow(
              column(
                8,
                selectizeInput(
                  "talklr_interaction",
                  "Choose LR interaction",
                  choices = NULL,
                  options = list(placeholder = "type to search…"),
                  width = "100%"
                )
              ),
              column(
                4,
                tags$div(style = "display:none;",
                         sliderInput(
                           "talklr_thresh",
                           "Expression threshold",
                           min = 0, max = 5,
                           value = if (exists("thresh_exprs", inherits = TRUE)) get("thresh_exprs") else 0.5,
                           step = 0.05
                         )
                ))
            ),
            hr(),
            plotOutput("talklr_plot", height = "520px")
            
            
          )
        ),
        tabPanel(
          "talklr pathway",
          fluidPage(
            fluidRow(
              column(
                8,
                selectizeInput(
                  "talklr_pathway",
                  "Choose pathway",
                  choices = NULL,
                  options = list(placeholder = "type to search…"),
                  width = "100%"
                )
              ),
              column(
                4,
                helpText("Aggregates only interactions that are currently selectable in the talklr tab.")
              )
            ),
            hr(),
            plotOutput("talklr_pathway_plot", height = "520px")
          )
        )
        
        
        
      )
    )
    
    
    
    
    
  )
)

server <- function(input, output, session) {
  cds_data <- reactiveVal(NULL)
  abbrev_lookup <- reactiveVal(NULL)
  label_for <- function(full) {
    tbl <- abbrev_lookup()
    if (is.null(tbl)) return(full)
    i <- match(full, tbl$full)
    ifelse(is.na(i), full, tbl$abbr[i])
  }
  # file-safe short label (uses your abbrev_lookup via label_for)
  short_file_label <- function(full) {
    lbl <- label_for(full)              # your short version
    if (is.null(lbl) || is.na(lbl) || !nzchar(lbl)) lbl <- full
    lbl <- gsub("\\s+", "", lbl)        # drop spaces
    lbl <- gsub("[^A-Za-z0-9+_.-]", "_", lbl)  # keep nice filename chars
    lbl
  }
  
  talklr_store <- reactiveValues(lr = NULL, glom = NULL)
  meta_data      <- reactiveVal(NULL) 
  analysis_result <- reactiveVal()
  std_tables_rv   <- reactiveVal(NULL)   # <-- add this line
  last_params <- reactiveValues(
    receiver_cell = NULL, 
    diff_target_cell = NULL, 
    normalization_label = NULL, 
    output_dir = NULL
  )
  
  learn_paths <- reactiveValues(
    simple_scatter = NULL,
    simple_cluster_treemap = NULL,
    simple_pruned_treemap = NULL,
    gnn_treemap = NULL,
    explain_pngs = character(0),
    true_gnn_treemap = NULL,
    true_explain_pngs = character(0)
  )
  learn_tables <- reactiveValues()
  learn_misc   <- reactiveValues()
  
  
  # Where the 4 text files live (change if needed)
  # ---- where the four files live ------------------------------------------------
  path_dir <- "top_path"
  
  category_files <- c(
    "Secreted Signaling"     = file.path(path_dir, "Secreted_Signaling_pathways.txt"),
    "ECM-Receptor"           = file.path(path_dir, "ECM_Receptor_pathways.txt"),
    "Cell-Cell Contact"      = file.path(path_dir, "Cell_Cell_Contact_pathways.txt"),
    "Non-protein Signaling"  = file.path(path_dir, "Non_protein_Signaling_pathways.txt")
  )
  
  .read_path_list <- function(types) {
    types <- intersect(types, names(category_files))
    if (!length(types)) return(character(0))
    out <- unlist(lapply(types, function(t) {
      fp <- category_files[[t]]
      if (!file.exists(fp)) {
        warning(sprintf("Missing pathway list for '%s': %s", t, fp))
        return(character(0))
      }
      x <- readLines(fp, warn = FALSE)
      x <- trimws(x); x <- x[nzchar(x)]
      x
    }), use.names = FALSE)
    sort(unique(out))
  }
  
  
  
  
  
  
  # ---- keep the pathway dropdown in sync ---------------------------------------
  
  
  
  DEFAULT_LR <- "Secreted Signaling"
  
  observeEvent(input$lr_classes, {
    cur <- input$lr_classes
    
    # If user deselects everything, snap back to DEFAULT_LR
    if (is.null(cur) || !length(cur)) {
      cur <- DEFAULT_LR
      updateCheckboxGroupInput(session, "lr_classes", selected = cur)
    }
    
    # Recompute pathways for the currently ticked types
    paths <- .read_path_list(cur)
    choices_now <- c("All", paths)
    
    # DO NOT preserve the previous selection: we want a clean default
    updateSelectizeInput(
      session, "path_included",
      choices  = choices_now,
      selected = "All",
      server   = TRUE
    )
  }, ignoreInit = FALSE)
  
  
  
  paths_current <- reactive({
    cur <- input$lr_classes
    if (is.null(cur) || !length(cur)) cur <- DEFAULT_LR
    .read_path_list(cur)
  })
  
  # ---- downstream helper: the effective set you should use ---------------------
  selected_paths <- reactive({
    cur <- input$lr_classes
    if (is.null(cur) || !length(cur)) cur <- DEFAULT_LR
    paths <- .read_path_list(cur)
    
    sel <- input$path_included %||% "All"
    if ("All" %in% sel) paths else intersect(sel, paths)
  })
  
  
  
  
  # same CSV used by the converter
  ortholog_file <- "database/mouse_to_human_orthologs.csv"
  
  load_ortholog_map <- function() {
    if (file.exists(ortholog_file)) {
      tryCatch(read.csv(ortholog_file, stringsAsFactors = FALSE), error=function(e) NULL)
    } else NULL
  }
  
  auto_humanize_genes <- function(symbols, cds=NULL) {
    # if converter already humanized, this is a cheap no-op
    info <- if (!is.null(cds)) attr(cds, "import_info") else NULL
    sp   <- if (!is.null(info)) info$species_guess else NA_character_
    syms <- as.character(symbols)
    if (identical(sp, "mouse")) {
      ortho <- load_ortholog_map()
      syms <- if (is.null(ortho) || nrow(ortho)==0) toupper(syms) else {
        mouse_uc <- toupper(ortho$mouse_symbol); human_uc <- toupper(ortho$human_symbol)
        idx <- match(toupper(syms), mouse_uc)
        out <- human_uc[idx]; out[is.na(out) | out==""] <- toupper(syms[is.na(out) | out==""])
        out
      }
    }
    toupper(syms)
  }
  
  
  # where we will serve images from (www is auto-served by Shiny)
  plots_www_dir <- normalizePath(file.path("www", "plots"), mustWork = FALSE)
  if (!dir.exists(plots_www_dir)) dir.create(plots_www_dir, recursive = TRUE)
  
  # keep final, absolute PNG paths here (so renderImage can use them)
  plot_paths <- reactiveValues(
    path = NULL, ligand = NULL, ligand_only = NULL, receptor_only = NULL,
    interaction = NULL, interaction_only = NULL
  )
  
  
  # convenience: build the 6 treemap names
  png_names <- function(receiver_label, target_label, norm) {
    c(
      sprintf("Tree_path_%s_to_%s_%s.png",             receiver_label, target_label, norm),
      sprintf("Tree_ligand_%s_to_%s_%s.png",           receiver_label, target_label, norm),
      sprintf("Tree_ligand_only_%s_to_%s_%s.png",      receiver_label, target_label, norm),
      sprintf("Tree_receptor_only_%s_to_%s_%s.png",    receiver_label, target_label, norm),
      sprintf("Tree_interaction_%s_to_%s_%s.png",      receiver_label, target_label, norm),
      sprintf("Tree_interaction_only_%s_to_%s_%s.png", receiver_label, target_label, norm)  # ensure suffix
    )
  }
  
  
  
  
  
  
  
  observeEvent(input$cdsfile, {
    req(input$cdsfile)
    ext <- tools::file_ext(input$cdsfile$name)
    
    note_id <- showNotification(
      if (ext == "h5ad") "Converting .h5ad → Monocle3 CDS…" else "Reading .rds…",
      type = "message", duration = NULL, closeButton = FALSE
    )
    on.exit(removeNotification(note_id), add = TRUE)
    
    withProgress(message = if (ext == "h5ad") "Converting .h5ad → CDS" else "Loading .rds", value = 0, {
      incProgress(0.2, detail = "Reading file")
      cds <- switch(
        ext,
        "h5ad" = convert_h5ad_to_cds(input$cdsfile$datapath),
        "rds"  = {
          obj <- readRDS(input$cdsfile$datapath)
          cls <- class(obj)
          if (inherits(obj, "cell_data_set")) {
            output$filetype_status <- renderUI(htmltools::div(class="status-msg",
                                                              HTML("<b>Detected Monocle CDS (.rds).</b> Loading…")))
            obj
          } else if ("Seurat" %in% cls) {
            output$filetype_status <- renderUI(htmltools::div(class="status-msg",
                                                              HTML("<b>Detected Seurat object (.rds).</b> Converting to Monocle3 CDS…")))
            convert_seurat_to_cds(obj)
          } else {
            stop("Unsupported .rds object (class: ", paste(cls, collapse=", "), "). Provide a Monocle3 CDS or a Seurat object.")
          }
        },
        stop("Unsupported file type: ", ext)
      )
      
      incProgress(0.6, detail = "Normalizing gene symbols")
      cds <- normalize_cds_symbols(cds)
      
      cds_data(cds)
      meta <- as.data.frame(colData(cds)); meta_data(meta)
      
      col_sizes <- sapply(meta, function(x) length(unique(x)))
      keep_cols <- names(col_sizes)[col_sizes > 1 & col_sizes < 200]
      updateSelectizeInput(session, "clustering_col",
                           choices = keep_cols,
                           selected = character(0),
                           options  = list(placeholder = "Select clustering column..."))
      
      # Reset pickers
      abbrev_lookup(NULL)
      updateSelectizeInput(session, "receiver_cell", choices = NULL)
      updateSelectizeInput(session, "diff_target_cell", choices = NULL)
      updateSelectizeInput(session, "sender_cells", choices = NULL)
      
      info <- attr(cds, "import_info")
      src  <- if (!is.null(info)) info$source_assay else "unknown"
      rnd  <- if (!is.null(info) && isTRUE(info$rounded_counts)) " (rounded to integers)" else ""
      origin_txt <- if (ext == "h5ad") "Scanpy (.h5ad) — converted to Monocle3 CDS."
      else if (!is.null(info) && startsWith(src, "Seurat:")) "Seurat (.rds) — converted to Monocle3 CDS."
      else "Monocle CDS (.rds)."
      
      output$filetype_status <- renderUI({
        htmltools::div(class = "status-msg",
                       HTML(sprintf("<b>Detected %s</b> Using assay: %s%s. Symbols in <code>gene_short_name</code>.",
                                    origin_txt, src, rnd)))
      })
    })
    
  })
  
  
  
  
  
  output$filetype_status <- renderUI({
    req(input$cdsfile)
    filetype <- tools::file_ext(input$cdsfile$name)
    cds <- cds_data()
    if (is.null(cds)) {
      msg <- if (filetype == "h5ad") "Reading .h5ad…" else "Reading .rds…"
    } else {
      msg <- if (filetype == "h5ad")
        "Detected Scanpy (.h5ad) — converted to Monocle3 CDS; gene symbols in gene_short_name."
      else
        "Detected Monocle CDS (.rds)."
    }
    HTML(paste0( msg, "<br>"))
  })
  
  
  
  
  
  
  # Helper to get cluster choices safely
  get_clusters <- function() {
    req(cds_data(), input$clustering_col)
    meta <- as.data.frame(colData(cds_data()))
    unique(as.character(meta[[input$clustering_col]]))
  }
  
  
  output$receiver_cell_ui <- renderUI({
    req(input$clustering_col, meta_data())
    clusters <- unique(as.character(meta_data()[[input$clustering_col]]))
    selectizeInput("receiver_cell", "Receiver Cell", choices = clusters, selected = NULL)
  })
  # Repeat for diff_target_cell_ui and sender_cell_ui
  
  
  output$diff_target_cell_ui <- renderUI({
    req(input$clustering_col, meta_data())
    clusters <- unique(as.character(meta_data()[[input$clustering_col]]))
    if (is.null(clusters)) return(NULL)
    selectizeInput(
      "diff_target_cell", "Diff Target Cell",
      choices = c("", clusters),
      selected = "",
      multiple = FALSE,
      options = list(placeholder = "Select diff target cell...", server = TRUE)
    )
  })
  
  output$sender_cell_ui <- renderUI({
    req(input$clustering_col, meta_data())
    clusters <- unique(as.character(meta_data()[[input$clustering_col]]))
    if (is.null(clusters)) return(NULL)
    selectizeInput(
      "sender_cells", "Sender Cells (multiple allowed)",
      choices = clusters,
      selected = NULL,
      multiple = TRUE,
      options = list(placeholder = "Select one or more sender cells...", server = TRUE)
    )
  })
  
  # Abbreviation legend: show ONLY if abbreviations are used
  output$abbrev_legend <- renderUI({
    map <- abbrev_lookup(); if (is.null(map)) return(NULL)
    
    used <- c(
      input$receiver_cell,
      if (!is.null(input$diff_target_cell) && nzchar(input$diff_target_cell)) input$diff_target_cell,
      input$sender_cells
    )
    used <- unique(used[!is.na(used) & nzchar(used)])
    if (!length(used)) return(NULL)
    
    # Keep only rows for the used labels (and in the same order)
    sub <- merge(
      data.frame(full = used, .order = seq_along(used), stringsAsFactors = FALSE),
      map,
      by = "full",
      all.x = TRUE, sort = FALSE
    )
    sub <- sub[!is.na(sub$abbr), c("abbr","full",".order")]
    # Only when abbr actually differs (case-insensitive)
    sub <- sub[tolower(sub$abbr) != tolower(sub$full), , drop = FALSE]
    if (!nrow(sub)) return(NULL)
    
    sub <- sub[order(sub$.order), c("abbr","full")]
    
    HTML(paste0(
      "<b>Legend:</b><br>",
      paste0("<span style='font-family:monospace'>", sub$abbr, " = ", sub$full, "</span>", collapse = "<br>")
    ))
  })
  
  
  
  
  
  
  
  output$selection_plot <- renderPlot({
    receiver <- input$receiver_cell
    diff_target <- input$diff_target_cell
    senders <- input$sender_cells
    
    show_treemap <- FALSE
    if (!is.null(input$run) && input$run > 0 && !is.null(receiver) && receiver != "") {
      pdf_file <- file.path(tempdir(), paste0("Tree_path_", receiver, ".pdf"))
      show_treemap <- file.exists(pdf_file)
    }
    
    center_x <- 0
    center_y <- 0
    base_offset_x <- 4.5
    
    senders_n <- if (!is.null(senders)) length(senders) else 0
    shrink_factor <- if (senders_n > 3) 0.90 ^ pmax(0, senders_n - 3) else 1
    r_main <- 1.0 * shrink_factor
    r_sender <- 0.7 * shrink_factor + 0.05
    offset_x <- base_offset_x * shrink_factor + 0.1
    sender_y_offset <- -3.5 * shrink_factor
    sender_x_pad <- (2.2 * shrink_factor) + 0.25
    spread <- if (senders_n > 1) (senders_n - 1) * sender_x_pad else 0
    xmin <- min(-2.5, -spread/2 - r_sender - 1)
    xmax <- max(7, offset_x + r_main + 2 + spread/2)
    ymin <- min(-4.5, center_y + sender_y_offset - r_sender - 1.2)
    ymax <- max(3.5, center_y + r_main + 2)
    
    gg <- ggplot() +
      theme_void() +
      coord_fixed(xlim = c(xmin, xmax), ylim = c(ymin, ymax), clip = "off")
    
    if (!is.null(receiver) && receiver != "") {
      gg <- gg +
        ggforce::geom_circle(aes(x0 = center_x, y0 = center_y, r = r_main),
                             fill = "#97c9ea", color = "black", linewidth = 1) +
        annotate("text", x = center_x, y = center_y,
                 label = label_for(receiver), size = 7 * shrink_factor, fontface = "bold")
    }
    
    if (!is.null(diff_target) && diff_target != "") {
      gg <- gg +
        ggforce::geom_circle(aes(x0 = center_x + offset_x, y0 = center_y, r = r_main),
                             fill = "#cc6b37", color = "black", linewidth = 1) +
        annotate("text", x = center_x + offset_x, y = center_y,
                 label = label_for(diff_target), size = 7 * shrink_factor, fontface = "bold") +
        geom_segment(aes(x = center_x + r_main, y = center_y,
                         xend = center_x + offset_x - r_main, yend = center_y),
                     arrow = arrow(type = "closed", length = unit(0.28, "inches")),
                     linewidth = 1.4)
    }
    
    if (!is.null(senders) && senders_n > 0) {
      n <- senders_n
      total_width <- if (n > 1) (n - 1) * sender_x_pad else 0
      sx <- center_x - total_width/2 + (0:(n-1)) * sender_x_pad
      sy <- rep(center_y + sender_y_offset, n)
      if (n == 3) {
        sy[1] <- sy[1] - 0.1 * shrink_factor
        sy[3] <- sy[3] - 0.1 * shrink_factor
      }
      circle_df <- data.frame(sx = sx, sy = sy, label = senders, stringsAsFactors = FALSE)
      circle_df$disp <- vapply(circle_df$label, label_for, character(1))
      
      
      
      
      arc_spread <- pi/2
      center_angle <- -pi/2
      arrow_thetas <- if (n == 1) center_angle else
        seq(center_angle - arc_spread/2, center_angle + arc_spread/2, length.out = n)
      arc_width <- 0.35
      arc_radius <- r_main + 0.03 * shrink_factor
      dot_radius <- 0.09 * shrink_factor + 0.01
      dot_dist <- r_main + 0.01 + dot_radius + 0.13 * shrink_factor
      dot_x <- center_x + dot_dist * cos(arrow_thetas)
      dot_y <- center_y + dot_dist * sin(arrow_thetas)
      
      arc_df <- data.frame()
      for (i in 1:n) {
        t0 <- arrow_thetas[i] - arc_width/2
        t1 <- arrow_thetas[i] + arc_width/2
        arc_points <- seq(t0, t1, length.out = 40)
        arc_x <- center_x + arc_radius * cos(arc_points)
        arc_y <- center_y + arc_radius * sin(arc_points)
        arc_df <- rbind(arc_df, data.frame(x = arc_x, y = arc_y, group = i))
      }
      
      gg <- gg +
        ggforce::geom_circle(data = circle_df, aes(x0 = sx, y0 = sy, r = r_sender),
                             fill = "white", color = "black", linewidth = 1) +
        # use the abbreviated label here:
        geom_text(data = circle_df, aes(x = sx, y = sy, label = disp),
                  size = 5 * shrink_factor + 0.7, fontface = "bold") +
        geom_segment(aes(x = sx, y = sy + r_sender, xend = dot_x, yend = dot_y),
                     data = cbind(circle_df, dot_x, dot_y), linewidth = 0.9) +
        geom_point(aes(x = dot_x, y = dot_y), size = 7.2 * shrink_factor + 0.1,
                   color = "black", shape = 16) +
        geom_path(data = arc_df, aes(x = x, y = y, group = group),
                  color = "#455a64", linewidth = 5 * shrink_factor + 0.2, lineend = "round")
    }
    
    
    
    
    
    gg
  })
  
  
  
  # 2) helper
  get_tbl <- function(name) {
    ar <- analysis_result(); if (is.null(ar)) return(NULL)
    out <- switch(name,
                  "pathway_n"     = ar$pathway_n,
                  "lig_rank_all"  = ar$lig_rank_all,
                  "rec_rank_all"  = ar$rec_rank_all,     # <--- NEW
                  "pval_pathway"  = ar$pval_pathway,
                  "pval_ligand"   = ar$pval_ligand,
                  "pval_receptor" = ar$pval_receptor,    # <--- NEW
                  NULL
    )
    if (is.null(out)) return(NULL)
    as.data.frame(out)
  }
  
  # --- once, when server starts
  # one-time DT init with canonical schema
  output$result_table <- DT::renderDT(
    DT::datatable(
      data.frame(
        kind = character(), entity = character(), pathway = character(),
        score = numeric(), score.perc = character(),
        p.value = numeric(), significant = character()
      ),
      options = list(pageLength = 10, scrollX = TRUE, deferRender = TRUE),
      rownames = FALSE
    ),
    server = TRUE
  )
  
  proxy <- DT::dataTableProxy("result_table")
  outputOptions(output, "result_table", suspendWhenHidden = FALSE)
  
  
  get_std_tbl <- function(which) {
    x <- std_tables_rv(); if (is.null(x)) return(NULL)
    x[[which]]
  }
  
  observeEvent(list(input$table_choice, std_tables_rv()), {
    req(std_tables_rv())
    which <- switch(input$table_choice,
                    "pathway_n"    = "pathway",
                    "lig_rank_all" = "ligand",
                    "rec_rank_all" = "receptor",
                    "int_rank_all" = "interaction",  # NEW
                    "pathway"
    )
    
    tbl <- get_std_tbl(which); req(tbl)
    DT::replaceData(proxy, tbl, resetPaging = FALSE, rownames = FALSE)
  }, ignoreInit = FALSE)
  
  
  
  
  
  # Abbreviate: initials for multi-token names; drop "cell(s)";
  # preserve short markers, numbers, and special tokens.
  abbreviate_clusters <- function(clusters) {
    clusters <- as.character(clusters)
    
    # split on spaces/hyphens/slashes; keep '_' and '+' inside tokens
    split_tokens <- function(s) {
      s1 <- gsub("([a-z])([A-Z])", "\\1 \\2", s, perl = TRUE)  # camelCase breaker
      s2 <- gsub("[\\-/]+", " ", s1, perl = TRUE)              # hyphen/slash -> space
      toks <- unlist(strsplit(s2, "\\s+"))
      toks[nzchar(toks)]
    }
    
    is_cells_word <- function(tok) grepl("^(?i)cells?$", tok, perl = TRUE)
    
    cap2 <- function(s) {
      if (!nzchar(s)) return("")
      if (nchar(s) == 1) toupper(s) else paste0(toupper(substr(s,1,1)), tolower(substr(s,2,2)))
    }
    
    process_token <- function(tok, force_initials = FALSE) {
      if (is_cells_word(tok)) return(NA_character_)
      
      # (A) letters ':' digits  -> compress left + strip leading zeros on number
      # e.g. "MmusDv:0000036" -> "MDv36"
      tok_trim <- trimws(tok)
      if (grepl("^[A-Za-z]+:[0-9]+$", tok_trim, perl = TRUE)) {
        left  <- sub(":[0-9]+$", "", tok_trim, perl = TRUE)
        right <- sub("^.*:",     "", tok_trim, perl = TRUE)
        right <- sub("^0+", "", right); if (right == "") right <- "0"
        
        # compress 'left': FIRST uppercase + LAST uppercase "hump"
        ups <- gregexpr("[A-Z]", left, perl = TRUE)[[1]]
        if (length(ups) >= 2 && ups[1] != -1) {
          first_uc <- substr(left, ups[1], ups[1])
          last_uc  <- ups[length(ups)]
          last_hump <- substr(left, last_uc, nchar(left))  # e.g., "Dv"
          left_abbr <- paste0(first_uc, last_hump)
        } else {
          left_abbr <- toupper(substr(left, 1, 1))
        }
        return(paste0(left_abbr, right))
      }
      
      # (B) letters '_' digits  -> two-letter titlecase of letters + digits (underscore removed)
      # e.g. "embryo_74" -> "Em74"
      if (grepl("^[A-Za-z]+_\\d+$", tok)) {
        letters <- sub("_\\d+$", "", tok)
        digits  <- sub("^.*_",   "", tok)
        return(paste0(cap2(letters), digits))
      }
      
      # (C) Preserve tokens with '+' or '_' as-is (e.g., "SP6+", "09_11w")
      if (grepl("[+_]", tok)) return(tok)
      
      # (D) Mixed alpha+digits like CD4, p53 -> keep whole (uppercased)
      if (grepl("^[A-Za-z]+\\d+$", tok) || grepl("^\\d+[A-Za-z]+$", tok)) return(toupper(tok))
      
      # (E) Pure numbers → keep
      if (grepl("^\\d+$", tok)) return(tok)
      
      # (F) Pure letters
      if (grepl("^[A-Za-z]+$", tok)) {
        if (force_initials) {
          # in multi-token names: keep very short markers (<=3) intact; else initial
          if (nchar(tok) <= 3) return(toupper(tok)) else return(toupper(substr(tok, 1, 1)))
        } else {
          # in 1–2 token names: keep short tokens (<=5) intact; else initial
          if (nchar(tok) <= 5) return(tok) else return(toupper(substr(tok, 1, 1)))
        }
      }
      
      # (G) Fallback: compose from parts (letters -> initial, digits kept, '+'/'_' kept)
      parts <- regmatches(tok, gregexpr("([A-Za-z]+)|(\\d+)|([+_])", tok, perl = TRUE))[[1]]
      if (length(parts)) {
        out <- vapply(parts, function(p) {
          if (grepl("^[A-Za-z]+$", p)) toupper(substr(p, 1, 1))
          else if (grepl("^\\d+$", p)) p
          else if (p %in% c("+","_")) p
          else ""
        }, character(1))
        paste0(out, collapse = "")
      } else {
        toupper(tok)
      }
    }
    
    to_abbrev <- function(s) {
      if (is.na(s) || s == "") return("")
      toks_all <- split_tokens(s)
      toks     <- toks_all[!vapply(toks_all, is_cells_word, logical(1))]
      if (!length(toks)) return(toupper(s))
      
      # Use initials if there are 3+ tokens after dropping "cell(s)"
      force_initials <- length(toks) >= 3
      parts <- Filter(Negate(is.na), lapply(toks, process_token, force_initials = force_initials))
      ab <- paste0(unlist(parts), collapse = "")
      if (ab == "") toupper(s) else ab
    }
    
    abbr <- vapply(clusters, to_abbrev, character(1))
    abbr <- make.unique(abbr, sep = "_")
    data.frame(full = clusters, abbr = abbr, stringsAsFactors = FALSE)
  }
  
  
  
  
  
  
  
  
  
  
  # When the clustering column is chosen, update abbreviations table ONCE
  observeEvent(input$clustering_col, {
    req(meta_data(), input$clustering_col)
    clusters <- unique(as.character(meta_data()[[input$clustering_col]]))
    table <- abbreviate_clusters(clusters)
    abbrev_lookup(table)  # store mapping for label_for()
    
    updateSelectizeInput(session, "receiver_cell",
                         choices  = clusters,
                         selected = character(0),
                         options  = list(placeholder = "Select receiver cell...")
    )
    
    updateSelectizeInput(session, "diff_target_cell",
                         choices  = clusters,
                         selected = character(0),
                         options  = list(placeholder = "Select diff target cell...")
    )
    
    updateSelectizeInput(session, "sender_cells",
                         choices  = clusters,
                         selected = character(0),
                         options  = list(placeholder = "Select one or more sender cells...")
    )
    
    
  })
  
  
  
  # Example button
  observeEvent(input$show_example, {
    demo_table <- readRDS("www/example_result_table.rds")
    # If demo_table doesn't match the canonical schema, coerce or just replace raw:
    DT::replaceData(proxy, as.data.frame(demo_table), resetPaging = TRUE, rownames = FALSE)
    
    output$analysis_status <- renderUI({
      HTML("<span style='color: #388e3c;'><b>Showing example output.</b></span>")
    })
  })
  
  
  # Cancel button logic: only a dummy message here
  observeEvent(input$cancel_analysis, {
    output$analysis_status <- renderUI({
      HTML("<span style='color: #b20000;'><b>Cancel only works for background jobs. This analysis cannot be canceled here.</b></span>")
    })
  })
  
  observeEvent(input$run, {
    req(cds_data())  # <-- make sure a CDS exists
    req(
      input$clustering_col != "",
      input$receiver_cell != "",
      input$diff_target_cell != "",
      length(input$sender_cells) > 0,
      length(input$path_included) > 0
    )
    output$analysis_status <- renderUI({HTML("<b>Analysis running ...</b>")})
    withProgress(message = "Running Analysis", value = 0, {
      receiver_cell    <- input$receiver_cell
      diff_target_cell <- input$diff_target_cell
      in_groups        <- input$sender_cells
      path_included <- selected_paths()
      req(length(path_included) > 0)  # guard against empty selection
      
      pair <- c(receiver_cell, diff_target_cell)
      
      
      ## === DEG with method in filename ===
      deg_file <- paste0(meta.file, "/deg_0.05_", paste(pair, collapse = '_'), "_", input$de_method, ".rds")
      if (file.exists(deg_file)) {
        deg <- readRDS(deg_file)
        message("Loaded existing DEG file: ", deg_file)
      } else {
        if (input$de_method == "DEsingle") {
          incProgress(0.1, detail = "Differential Expression")
          temp_cds <- cds_data()[, colData(cds_data())[[input$clustering_col]] %in% pair]
          group <- factor(colData(temp_cds)[[input$clustering_col]])
          suppressWarnings({
            results <- DEsingle(counts = counts(temp_cds), group = group, parallel = TRUE)
          })
          results.classified <- DEtype(results = results, threshold = 0.05)
          results.sig <- results.classified[results.classified$pvalue < 0.05, ]
          deg <- results.sig[,c(12,13,11,20,21,23,24)]
          deg$Gene <- rowData(temp_cds[rownames(deg),])$gene_short_name
          deg <- deg[,c(8,1,2,3,4,5,6,7)]
          deg$norm_total_mean_1 <- deg$norm_total_mean_1 * 1000
          deg$norm_total_mean_2 <- deg$norm_total_mean_2 * 1000
          colnames(deg) <- c("Gene","Cluster1_total_mean", "Cluster2_total_mean", "foldChange" , "pvalue", "pvalue.adj.FDR","Type", "State")
          saveRDS(deg, deg_file)
        } else if (input$de_method == "edgeR") {
          incProgress(0.1, detail = "Differential Expression (edgeR)")
          
          temp_cds <- cds_data()[, colData(cds_data())[[input$clustering_col]] %in% pair]
          
          # --- 1) Get a counts-like matrix & sanitize ---
          expr_mat <- as.matrix(counts(temp_cds))
          
          # (optional) if you really want to collapse duplicate rownames, do it here
          if (any(duplicated(rownames(expr_mat)))) {
            expr_mat <- rowsum(expr_mat, group = rownames(expr_mat))
          }
          
          # Heuristic check for non-integer-looking values
          nz <- expr_mat[expr_mat > 0]
          if (length(nz) == 0L) stop("All-zero matrix after subsetting; edgeR needs raw counts.")
          if (median(nz) < 1 || max(nz) < 20) {
            message("Counts look non-integer/normalized; rounding to nearest integer.")
          }
          
          expr_mat[is.na(expr_mat)] <- 0
          expr_mat[expr_mat < 0] <- 0
          if (any(abs(expr_mat - round(expr_mat)) > 1e-8)) expr_mat <- round(expr_mat)
          
          # Drop all-zero rows/cols
          keep_genes  <- rowSums(expr_mat) > 0
          expr_mat    <- expr_mat[keep_genes, , drop = FALSE]
          keep_cells  <- colSums(expr_mat) > 0
          expr_mat    <- expr_mat[, keep_cells, drop = FALSE]
          if (ncol(expr_mat) == 0 || nrow(expr_mat) == 0) stop("No non-zero counts remain; cannot run edgeR.")
          
          # --- 2) Align metadata & ensure two groups exist ---
          meta <- as.data.frame(colData(temp_cds))
          meta <- meta[keep_cells, , drop = FALSE]
          meta$cluster <- as.character(meta[[input$clustering_col]])
          meta <- meta[meta$cluster %in% pair, , drop = FALSE]
          expr_mat <- expr_mat[, rownames(meta), drop = FALSE]
          if (length(unique(meta$cluster)) != 2) {
            stop("edgeR needs exactly 2 groups; got: ", paste(unique(meta$cluster), collapse = ", "))
          }
          
          # --- 3) Build robust pseudo-replicates per group ---
          # Use up to 3 reps but never more than the number of cells in that group.
          make_reps <- function(idx, n_max = 3) {
            n <- length(idx)
            if (n <= 3) return(setNames(rep(1, n), idx))  # single rep if tiny group
            k <- min(n_max, n)                            # number of reps
            # split as evenly as possible
            bins <- cut(seq_len(n), breaks = k, labels = FALSE)
            setNames(bins, idx)
          }
          
          meta$pseudorep <- NA_character_
          for (cl in unique(meta$cluster)) {
            idx <- which(meta$cluster == cl)
            bins <- make_reps(idx, n_max = 3)
            meta$pseudorep[idx] <- paste0(cl, "_rep", bins)
          }
          
          # Aggregate counts to pseudo-bulk
          pseudo_bulk <- rowsum(t(expr_mat), meta$pseudorep)  # reps x genes
          pseudo_bulk <- t(pseudo_bulk)                       # genes x reps
          
          # Ensure both groups represented
          rep_groups <- factor(sub("_rep\\d+$", "", colnames(pseudo_bulk)))
          if (length(levels(rep_groups)) != 2) {
            stop("Pseudo-bulk construction failed to create two groups.")
          }
          
          # --- 4) edgeR DE ---
          suppressPackageStartupMessages(library(edgeR))
          dge <- DGEList(counts = pseudo_bulk, group = rep_groups)
          dge <- calcNormFactors(dge)
          design <- model.matrix(~ rep_groups)
          dge <- estimateDisp(dge, design)
          fit <- glmQLFit(dge, design)
          qlf <- glmQLFTest(fit, coef = 2)
          res <- topTags(qlf, n = Inf)$table
          
          # --- 5) Means per cluster (average over that cluster’s reps) ---
          cl1 <- levels(rep_groups)[1]
          cl2 <- levels(rep_groups)[2]
          cl1_cols <- grep(paste0("^", cl1, "_rep"), colnames(pseudo_bulk))
          cl2_cols <- grep(paste0("^", cl2, "_rep"), colnames(pseudo_bulk))
          Cluster1_total_mean <- rowMeans(pseudo_bulk[, cl1_cols, drop = FALSE])
          Cluster2_total_mean <- rowMeans(pseudo_bulk[, cl2_cols, drop = FALSE])
          
          # Map Ensembl -> symbols from the CDS we actually used (post filtering)
          # --- 5) Map Ensembl -> symbols from the CDS actually used ---
          gene_symbols <- rowData(temp_cds)$gene_short_name
          names(gene_symbols) <- rownames(temp_cds)
          
          # Keep only the symbols for the genes in your filtered expr_mat
          gene_symbols <- gene_symbols[rownames(expr_mat)]
          
          
          deg <- data.frame(
            Gene                 = gene_symbols[rownames(res)],
            Ensembl              = rownames(res),
            Cluster1_total_mean  = Cluster1_total_mean[rownames(res)],
            Cluster2_total_mean  = Cluster2_total_mean[rownames(res)],
            foldChange           = res$logFC,
            pvalue               = res$PValue,
            pvalue.adj.FDR       = res$FDR,
            Type                 = ifelse(res$logFC > 0, "up", "down"),
            State                = ifelse(res$logFC > 0, "up", "down"),
            stringsAsFactors = FALSE
          )
          
          saveRDS(deg, deg_file)
        }
        else if (input$de_method == "MAST") {
          incProgress(0.1, detail = "Differential Expression")
          temp_cds <- cds_data()[, colData(cds_data())[[input$clustering_col]] %in% pair]
          expr_mat <- as.matrix(counts(temp_cds))
          cell_meta <- data.frame(group = factor(colData(temp_cds)[[input$clustering_col]]))
          sca <- FromMatrix(log2(expr_mat + 1), cell_meta, data.frame(gene=rownames(expr_mat)))
          
          # Fit model
          zlmCond <- zlm(~group, sca)
          group_levels <- levels(cell_meta$group)
          contrast_name <- paste0('group', group_levels[2])
          
          summaryCond <- summary(zlmCond, doLRT=contrast_name)
          summaryDt <- summaryCond$datatable
          # Extract logFC as a NUMERIC VECTOR
          logFC_table <- summaryDt[component == "logFC" & contrast == contrast_name]
          if ("coef" %in% names(logFC_table)) {
            logFC_vec <- as.numeric(logFC_table$coef)
          } else if ("logFC" %in% names(logFC_table)) {
            logFC_vec <- as.numeric(logFC_table$logFC)
          } else {
            stop("Cannot find logFC/coef column in summaryDt.")
          }
          names(logFC_vec) <- logFC_table$primerid
          # P-values from the Hurdle test (component == 'H')
          hurdle_table <- summaryDt[component == "H" & contrast == contrast_name]
          pval <- as.numeric(hurdle_table$`Pr(>Chisq)`)
          adj_pval <- p.adjust(pval, method = "fdr")
          names(pval) <- hurdle_table$primerid
          names(adj_pval) <- hurdle_table$primerid
          # Combine into DEG table
          common_genes <- intersect(names(logFC_vec), names(pval))
          deg <- data.frame(
            Gene = rowData(temp_cds)[common_genes, ]$gene_short_name
            # Then reorder columns as needed
            ,
            Cluster1_total_mean = NA,
            Cluster2_total_mean = NA,
            foldChange = logFC_vec[common_genes],
            pvalue = pval[common_genes],
            pvalue.adj.FDR = adj_pval[common_genes],
            Type = ifelse(logFC_vec[common_genes] > 0, "up", "down"),
            State = ifelse(logFC_vec[common_genes] > 0, "up", "down"),
            stringsAsFactors = FALSE
          )
          deg <- deg[, c("Gene", "Cluster1_total_mean", "Cluster2_total_mean", "foldChange", "pvalue", "pvalue.adj.FDR", "Type", "State")]
          saveRDS(deg, deg_file)
        }
      }
      
      ## === Enrichment + talklr (stable on-disk cache; no tempdir, no reloads) ===
      
      # 0) Stable cache inside the project (survives restarts; portable)
      cache_dir <- "cache"
      if (!dir.exists(cache_dir)) dir.create(cache_dir, recursive = TRUE)
      
      # --- Enrichment (receiver & target) ---
      
      # Common objects for enrichment
      temp_cds_r <- cds_data()[, colData(cds_data())[[input$clustering_col]] %in% in_groups]
      colData(temp_cds_r)[[input$clustering_col]] <- droplevels(colData(temp_cds_r)[[input$clustering_col]])
      
      GCMat <- counts(temp_cds_r)
      rownames(GCMat) <- rowData(temp_cds_r)$gene_short_name
      colnames(GCMat) <- as.character(colnames(temp_cds_r))
      
      # Write clustering file into cache (NOT database/)
      barfile <- file.path(cache_dir, "barcodetype.txt")
      clustering <- data.frame(
        Barcode = colnames(temp_cds_r),
        Cluster = colData(temp_cds_r)[[input$clustering_col]],
        check.names = FALSE
      )
      write.table(clustering, barfile, sep = "\t", row.names = FALSE)
      BarCluTable <- read.table(barfile, sep = "\t", header = TRUE, stringsAsFactors = FALSE)
      
      pval  <- 0.15
      logfc <- 0.15
      cores <- max(1, parallel::detectCores() - 1)
      
      # Receiver enrichment
      RecClu  <- receiver_cell
      LigClu  <- setdiff(in_groups, c(RecClu, diff_target_cell))
      RecClus <- getHighExpGene(GCMat, BarCluTable, RecClu, LigClu, pval, logfc, cores)
      
      # Target enrichment
      RecClu_t  <- diff_target_cell
      LigClu_t  <- setdiff(in_groups, RecClu_t)
      RecClus_t <- getHighExpGene(GCMat, BarCluTable, RecClu_t, LigClu_t, pval, logfc, cores)
      
      # Canonical names used later; keep in memory
      enrich_1 <- RecClus
      enrich_2 <- RecClus_t
      
      # (Optional) Save snapshots for debugging (we don't read them back)
      saveRDS(enrich_1, file.path(cache_dir, paste0(receiver_cell, "_enriched.rds")))
      saveRDS(enrich_2, file.path(cache_dir, paste0(diff_target_cell, "_enriched.rds")))
      
      # --- talklr: always recompute; do NOT load old .Rdata ---
      
      talklr_file <- file.path(cache_dir, paste0("talklr_", receiver_cell, ".Rdata"))
      
      temp_cds2 <- cds_data()[, colData(cds_data())[[input$clustering_col]] %in% c(receiver_cell, in_groups)]
      colData(temp_cds2)[[input$clustering_col]] <- droplevels(colData(temp_cds2)[[input$clustering_col]])
      
      cell_group_df <- data.frame(
        cell  = colnames(temp_cds2),
        group = colData(temp_cds2)[[input$clustering_col]],
        check.names = FALSE
      )
      
      glom_normal <- aggregate_gene_expression(temp_cds2, cell_group_df = cell_group_df, norm_method = "size_only")
      rownames(glom_normal) <- rowData(temp_cds2)$gene_short_name
      glom_normal <- glom_normal / as.vector(table(colData(temp_cds2)[[input$clustering_col]]))
      glom_normal <- as.data.frame(glom_normal)
      glom_normal$genes <- toupper(as.character(rownames(glom_normal)))
      glom_normal <- glom_normal[, c(ncol(glom_normal), 1:(ncol(glom_normal)-1))]
      glom_normal[, 2:ncol(glom_normal)] <- lapply(glom_normal[, 2:ncol(glom_normal)], as.numeric)
      
      min_exprs    <- 1.713682e-11
      thresh_exprs <- 0.909563e-08
      assign("thresh_exprs", thresh_exprs, inherits = TRUE)
      updateSliderInput(
        session, "talklr_thresh",
        value = thresh_exprs, min = 0,
        max = max(thresh_exprs * 50, 1e-3),
        step = signif(thresh_exprs / 5, 2)
      )
      
      path_included <- selected_paths(); req(length(path_included) > 0)
      receptor_ligand_sub <- receptor_ligand[receptor_ligand$pathway_name %in% path_included, ]
      
      lr_glom_normal <- make_expressed_net(
        glom_normal,
        expressed_thresh = as.numeric(thresh_exprs),
        receptor_ligand_sub,
        KL_method    = "product",
        pseudo_count = as.numeric(min_exprs)
      )
      lr_glom_normal <- dplyr::arrange(lr_glom_normal, dplyr::desc(KL))
      lr_glom_normal <- lr_glom_normal[lr_glom_normal$Pair.Evidence == "literature supported", ]
      lr_glom_normal <- coerce_lr_numeric_blocks(lr_glom_normal, glom_normal)
      
      lr_glom_normal <- restrict_lr_to_receiver(
        lr = lr_glom_normal, gl = glom_normal,
        receiver_id      = input$receiver_cell,
        expressed_thresh = as.numeric(thresh_exprs),
        min_edge_frac    = 0.00,
        allow_self       = TRUE
      )
      
      # Save for debugging (optional); main objects stay in memory
      save(list = c("lr_glom_normal", "glom_normal", "deg"), file = talklr_file)
      
      talklr_store$lr   <- lr_glom_normal
      talklr_store$glom <- glom_normal
      
      # --- Downstream section using deg / enrich_* ---
      
      temp_cds2 <- cds_data()[, colData(cds_data())[[input$clustering_col]] %in% c(receiver_cell, diff_target_cell, in_groups)]
      colData(temp_cds2)[[input$clustering_col]] <- droplevels(colData(temp_cds2)[[input$clustering_col]])
      
      GCMat <- counts(temp_cds2)
      rownames(GCMat) <- rowData(temp_cds2)$gene_short_name
      colnames(GCMat) <- as.character(colnames(temp_cds2))
      
      clustering <- data.frame(
        Barcode = colnames(temp_cds2),
        Cluster = colData(temp_cds2)[[input$clustering_col]],
        check.names = FALSE
      )
      barfile <- file.path(cache_dir, "barcodetype.txt")
      write.table(clustering, barfile, sep = "\t", row.names = FALSE)
      BarCluFile  <- barfile
      BarCluTable <- read.table(BarCluFile, sep = "\t", header = TRUE, stringsAsFactors = FALSE)
      
      # Ensure 'deg' exists; read if you have a file, else fallback to empty
      if (!exists("deg", inherits = TRUE)) {
        deg <- tryCatch(
          readRDS(deg_file),
          error = function(e) data.frame(Gene = character(), State = character(), stringsAsFactors = FALSE)
        )
      }
      
      filtered <- c(
        setdiff(intersect(c(enrich_2), c(subset(deg, deg$State %in% c("down", "up"))$Gene)), enrich_1),
        setdiff(subset(deg, deg$State %in% c("up"))$Gene, c(enrich_1, enrich_2)),
        intersect(intersect(enrich_1, enrich_2), subset(deg, deg$State %in% c("down", "up"))$Gene)
      )
      
      # Pathway analysis
      RecClu <- receiver_cell
      LigClu <- setdiff(in_groups, c(RecClu, diff_target_cell))
      pval <- 0.05
      logfc <- 0.15
      LigRecLib <- "database/LigRec.txt"
      TFTarLib <- "database/TFTargetGene.txt"
      RecTFLib <- "database/RecTF.txt"
      # Always set BarCluFile
      BarCluFile <- "database/barcodetype.txt"
      #clustering <- data.frame(Barcode = colnames(temp_cds2), Cluster = colData(temp_cds2)[[input$clustering_col]])
      #write.table(clustering, BarCluFile, sep = "\t")
      
      netList <- RunMLnet(GCMat, BarCluFile, RecClu, LigClu,
                          pval, logfc,
                          LigRecLib, TFTarLib, RecTFLib, 
                          lr_glom_normal = lr_glom_normal, 
                          filtered = filtered,
                          RecClusGenes = enrich_1,
                          LigClusGenes = enrich_2)
      
      # before scoring/plotting:
      # If you want to insist targets be DEGs, pass deg$Gene and set require_deg_targets=TRUE
      has_edges <- function(nl) {
        nrow(as.data.frame(nl$LigRec)) + nrow(as.data.frame(nl$RecTF)) + nrow(as.data.frame(nl$TFTar)) > 0
      }
      
      if (!has_edges(netList)) {
        message("GNN skip: graph has no usable edges after pruning.")
        # emit empty canonical tables and blank plots, then return
        analysis_result(list(
          pathway_n=data.frame(),
          lig_rank_all=data.frame(),
          rec_rank_all=data.frame(),
          int_rank_all=data.frame()
        ))
        return(invisible(NULL))
      }
      
      netList <- prune_net_to_full_paths(
        netList,
        deg_genes = if (exists("deg")) as.character(deg$Gene) else NULL,
        require_deg_targets = FALSE
      )
      
      
      workdir <- paste(user.file,"metadata", sep="/")
      PyHome <- "python3"
      DrawMLnet(netList, LigClu, RecClu, workdir, PyHome, plotMLnet = FALSE)
      saveRDS(netList, paste0(out.file, "/", workdir, "/netList.rds"))
      
      # --- Optional SCENIC weighting (AUCell) ---
      tf_w <- if (isTRUE(input$use_scenic)) {
        compute_scenic_tf_weights(
          cds = cds_data(),
          scMlnet_results = netList,
          cluster_col = input$clustering_col,
          receiver_cell = input$receiver_cell,
          target_cell   = input$diff_target_cell,
          beta = 2, w_min = 0.5, w_max = 2  # your current settings
        )
      } else NULL
      
      
      # Downstream reporting, treemap and table
      # --- Downstream reporting: write PNGs to a temp folder, then copy into www/plots
      normalization_label <- "db"
      
      # Let top_pathway write to a temp output dir
      tmp_outdir <- tempfile("treemaps_")
      dir.create(tmp_outdir, recursive = TRUE, showWarnings = FALSE)
      
      # inside observeEvent(input$run, { ... })
      map <- abbrev_lookup()  # current full<->abbr table
      receiver_short <- short_file_label_map(input$receiver_cell,   map)
      target_short   <- short_file_label_map(input$diff_target_cell, map)
      
      
      
      result <- top_pathway(
        scMlnet_results = netList,
        deg             = deg,
        receiver_cell   = receiver_cell,      # full (inputs unchanged)
        target_cell     = diff_target_cell,   # full
        lr_glom_normal  = lr_glom_normal,
        method          = "KL_rec",
        output_dir      = tmp_outdir,
        tf_weight_vec   = tf_w,
        scenic_alpha    = 1.0,
        scenic_gamma    = 1.5,
        scenic_clip     = c(0.5, 4),
        scenic_center   = TRUE,
        ## NEW: what to write on the plots & filenames
        receiver_label  = receiver_short,
        target_label    = target_short
      )
      
      
      
      incProgress(0.15, detail = "Preparing outputs")
      
      # Store tables (prefer result$db if present)
      res_norm <- if (!is.null(result$db)) result$db else result
      analysis_result(res_norm)
      
      message("Tables available: ", paste(names(analysis_result()), collapse=", "))
      message("pathway_n class: ", paste(class(analysis_result()$pathway_n), collapse=", "))
      message("pathway_n dim: ", tryCatch(paste(dim(analysis_result()$pathway_n), collapse=" x "), error=function(e) "no dim"))
      
      cat("#rows pathway_n:", nrow(result$db$pathway_n), "\n")
      cat("#rows int_rank_all:", nrow(result$db$int_rank_all), "\n")
      cat("#rows rec_rank_all:", nrow(result$db$rec_rank_all), "\n")
      cat("Has RecTF?", !is.null(netList$RecTF), "\n")
      cat("Has TFTar?", !is.null(netList$TFTar), "\n")
      
      
      
      # ------- Learning modules (optional) -------
      # 1) Prepare inputs both modules will reuse
      res_norm <- if (!is.null(result$db)) result$db else result
      pw_tbl   <- as.data.frame(res_norm$pathway_n)  # pathway, score, etc.
      
      # 'expres' magnitudes already computed inside top_pathway() before normalization
      # If you need it here, recompute the (weighted, positive magnitude) vector again,
      # or pass it back from top_pathway. We'll reconstruct quickly from 'deg' as below:
      expres_mag <- {
        v <- as.numeric(deg$foldChange); names(v) <- as.character(deg$Gene)
        v <- v[!is.na(v) & !is.nan(v)]
        v[v ==  Inf] <- max(v[is.finite(v)], na.rm = TRUE)
        v[v == -Inf] <- min(v[is.finite(v)], na.rm = TRUE)
        # make positive magnitude
        if (any(v < 0, na.rm = TRUE)) v <- abs(v) else { idx <- v < 1 & v > 0; v[idx] <- 1/v[idx] }
        v[is.na(v)] <- 0; v
      }
      
      # Where to save images
      learn_outdir <- file.path(plots_www_dir, "learn")
      if (!dir.exists(learn_outdir)) dir.create(learn_outdir, recursive = TRUE)
      
      if (!dir.exists(plots_www_dir)) dir.create(plots_www_dir, recursive = TRUE)
      
      # SIMPLE (existing)
      if (isTRUE(input$simple_learning)) {
        # Simple: Active vs Unique and Activity vs Unique
        # SIMPLE
        pw_cov_simple <- pathway_coverage_tables(
          scMlnet_results = list(LigRec = netList$LigRec, RecTF = netList$RecTF, TFTar = netList$TFTar),
          receptor_ligand = receptor_ligand,
          deg = deg,
          tf_weight_vec = tf_w,
          limit_pathways = unique(result$db$pathway_n$pathway)     # <— NEW
        )
        
        
        # 1) Active (coverage frac) vs Unique
        simp_au_png <- file.path(plots_www_dir, sprintf(
          "Simple_active_unique_%s_to_%s.png", receiver_short, target_short))
        
        png(simp_au_png, width = 900, height = 600)
        df <- data.frame(
          x = pw_cov_simple$targ_coverage_frac,
          y = pw_cov_simple$targ_unique_sum,
          label = pw_cov_simple$pathway
        )
        df <- df[is.finite(df$x) & is.finite(df$y), ]
        
        p <- ggplot(df, aes(x, y)) +
          geom_point(size = 3) +
          ggrepel::geom_text_repel(aes(label = label), max.overlaps = Inf, size = 3.2) +
          labs(
            x = "Active (fraction of target signal covered)",
            y = "Unique target signal (sum)",
            title = "Active vs Unique — pathways (Simple)"
          ) +
          theme_minimal(base_size = 13)
        print(p)
        dev.off()
        
        
        # 2) Activity score (from pathway_n) vs Unique
        simp_avsu_png <- file.path(plots_www_dir, sprintf(
          "Simple_activity_unique_%s_to_%s.png", receiver_short, target_short))
        
        png(simp_avsu_png, width = 900, height = 600)
        act <- setNames(result$db$pathway_n$score, result$db$pathway_n$pathway)
        df <- data.frame(
          x = act[pw_cov_simple$pathway],
          y = pw_cov_simple$targ_unique_sum,
          label = pw_cov_simple$pathway
        )
        df <- df[is.finite(df$x) & is.finite(df$y), ]
        
        
        
        
        p <- ggplot(df, aes(x, y)) +
          geom_point(size = 3) +
          ggrepel::geom_text_repel(aes(label = label), max.overlaps = Inf, size = 3.2) +
          labs(
            x = "Activity score (pathway_n$score)",
            y = "Unique target signal (sum)",
            title = "Activity vs Unique — pathways (Simple)"
          ) +
          theme_minimal(base_size = 13)
        print(p)
        dev.off()
        
        
        
        
        learn_paths$simple_active_unique   <- simp_au_png
        learn_paths$simple_activity_unique <- simp_avsu_png
        
        if (is.null(result$db$int_rank_all) && !is.null(result$db$lr_rank_all)) {
          result$db$int_rank_all <- result$db$lr_rank_all
        }
        
        sl <- simple_learning_run(
          lig_rank_all     = result$db$lig_rank_all,
          pathway_n        = result$db$pathway_n,
          int_rank_all     = result$db$int_rank_all,      # used for LR-pruned
          rec_rank_all     = result$db$rec_rank_all,
          scMlnet_results  = list(LigRec = netList$LigRec,
                                  RecTF  = netList$RecTF,
                                  TFTar  = netList$TFTar),
          receptor_ligand  = receptor_ligand,
          deg              = deg, 
          tf_weight_vec    = tf_w,
          out_dir          = plots_www_dir,
          receiver         = receiver_short,
          target           = target_short,
          norm_label       = "db",
          mode             = "tf_centered",       # ensure SIMPLE uses ligand-weighted, tf_centered "ligand_weighted","ligand_tf_weighted","interaction_weighted","ligand_set",
          # "tf_centered","tftarget_centered","tf_centered_unweighted","tftarget_centered_unweighted"
          mirror.y         = FALSE,
          rank_method      = input$rank_method,
          k_per_cluster    = input$k_per_cluster,
          w_active         = 0.5, 
          w_unique         = 0.5,
          embed_method     ="umap", #umap/ mds/ pca/ auto
          cluster_method   ="kmeans", #kmeans hdbscan    dbscan
          n_clusters       = 4,
          umap_args = list(min_dist = 0.15, n_neighbors = 6)  ,
          # list(
          #   metric = "cosine",
          #   n_neighbors = k,
          #   min_dist = md,
          #   init = "spectral",
          #   pca = pmin(50, ncol(M), n - 1),         # light PCA cap
          #   nn_method = "brute",                    # exact kNN (N is small)
          #   dens_frac = 0.3, dens_lambda = 2.0,     # densMAP to damp “tear-apart”
          #   ret_model = FALSE, scale = FALSE,
          #   seed = seed, verbose = FALSE
          # )
          
        )
        
        # publish all paths (note: bubble is used for BOTH simple_* plots to “merge” them)
        learn_paths$simple_scatter            <- sl$paths$scatter
        learn_paths$simple_cluster_treemap    <- sl$paths$cluster_treemap
        learn_paths$simple_pruned_treemap     <- sl$paths$pruned_treemap
        learn_paths$simple_lr_pruned_treemap  <- sl$paths$lr_pruned_treemap
        learn_paths$simple_active_unique      <- sl$paths$au_bubble   # merged bubble
        learn_paths$simple_activity_unique    <- sl$paths$au_bubble   # merged bubble
        
        # NEW: guard + fallback for the LR-pruned treemap
        if (!is.null(sl$paths$lr_pruned_treemap) && nzchar(sl$paths$lr_pruned_treemap) &&
            file.exists(sl$paths$lr_pruned_treemap)) {
          learn_paths$simple_lr_pruned_treemap <- sl$paths$lr_pruned_treemap
        } else {
          message("[Simple] LR pruned treemap missing; creating placeholder.")
          png_fallback <- file.path(plots_www_dir, sprintf(
            "Simple_lr_pruned_treemap_%s_to_%s_placeholder.png", receiver_short, target_short))
          png(png_fallback, 900, 700); plot.new(); text(.5, .5, "No LR tiles after pruning", cex = 1.4); dev.off()
          learn_paths$simple_lr_pruned_treemap <- png_fallback
        }
        
        message("int_rank_all rows: ", nrow(result$db$int_rank_all))
        message("LR pruned path: ", sl$paths$lr_pruned_treemap)
        message("Exists? ", file.exists(sl$paths$lr_pruned_treemap))
        
        
      }
      
      # Also produce the "Deep Explanations (GNN-like)" here, but skip its treemap
      dl <- deep_learning_run(
        pathway_n    = result$db$pathway_n,
        int_rank_all = result$db$int_rank_all,
        rec_rank_all = result$db$rec_rank_all,
        scMlnet_results = list(
          LigRec = netList$LigRec,
          RecTF  = netList$RecTF,
          TFTar  = netList$TFTar
        ),
        deg           = deg,
        tf_weight_vec = tf_w,
        out_dir       = plots_www_dir,
        receiver      = receiver_short,
        target        = target_short,
        norm_label    = "db",
        topK = 3, edges_per_path = 20,
        k_rt = 40, k_tft = 80, q_rt = 0.80
      )
      if (!is.null(dl)) {
        # don't surface dl$paths$gnn anywhere (that’s the redundant treemap)
        learn_paths$explain_pngs <- dl$paths$explain
      }
      
      
      # TRUE (Deep Learning) — single checkbox now
      if (isTRUE(input$deep_learning)) {
        tg <- true_gnn_run(
          scMlnet_results = list(LigRec = netList$LigRec,
                                 RecTF  = netList$RecTF,
                                 TFTar  = netList$TFTar),
          deg         = deg,
          out_dir     = plots_www_dir,
          receiver    = receiver_short,
          target      = target_short,
          norm_label  = "db",
          tf_weight_vec = tf_w,
          epochs = 150, d_model = 32, learn_rate = 1e-2
        )
        learn_paths$true_gnn_treemap  <- tg$paths$gnn
        learn_paths$true_explain_pngs <- tg$paths$explain
        
        # unify/derive tables (handles older true_gnn_run)
        tbl <- ensure_true_tables(tg)
        
        if (!is.null(tbl$pathway_n) && nrow(tbl$pathway_n) &&
            !is.null(tbl$lig_rank_all) && nrow(tbl$lig_rank_all)) {
          
          # --- True → Simple-style, but based on interaction-weighted similarity ---
          sl_true <- simple_learning_run(
            lig_rank_all = tbl$lig_rank_all,
            pathway_n    = tbl$pathway_n,
            int_rank_all = tbl$int_rank_all,
            scMlnet_results = list(LigRec = netList$LigRec, RecTF = netList$RecTF, TFTar = netList$TFTar), # we reuse your pruned net
            deg           = deg, receptor_ligand = receptor_ligand, tf_weight_vec = tf_w,
            out_dir      = file.path(plots_www_dir, "true"),
            receiver     = receiver_short, target = target_short, norm_label = "db",
            mode         = "tftarget_weighted",
            mirror.y     = FALSE
          )
          learn_paths$true_simple_scatter         <- sl_true$paths$scatter
          learn_paths$true_simple_cluster_treemap <- sl_true$paths$cluster_treemap
          learn_paths$true_simple_pruned_treemap  <- sl_true$paths$pruned_treemap
          
          
        } else {
          learn_paths$true_simple_scatter <- learn_paths$true_simple_cluster_treemap <- learn_paths$true_simple_pruned_treemap <- NULL
          message("TRUE: missing tables after GNN; skipping simple-style plots.")
        }
        
        # ---------- Deep tab analytics: Active vs Unique + minimal cover ----------
        
        
        # compute & stash
        # DEEP
        pw_cov <- pathway_coverage_tables(
          scMlnet_results = list(LigRec = netList$LigRec, RecTF = netList$RecTF, TFTar = netList$TFTar),
          receptor_ligand = receptor_ligand, deg = deg, tf_weight_vec = tf_w,
          limit_pathways = unique(tbl$pathway_n$pathway)           # <— NEW
        )
        
        learn_tables$deep_active_unique <- pw_cov
        
        # 1) Active vs Unique
        au_png <- file.path(plots_www_dir, sprintf(
          "Deep_active_unique_%s_to_%s.png", receiver_short, target_short))
        
        png(au_png, width = 900, height = 600)
        df <- data.frame(
          x = pw_cov$targ_coverage_frac,
          y = pw_cov$targ_unique_sum,
          label = pw_cov$pathway
        )
        df <- df[is.finite(df$x) & is.finite(df$y), ]
        
        p <- ggplot(df, aes(x, y)) +
          geom_point(size = 3) +
          ggrepel::geom_text_repel(aes(label = label), max.overlaps = Inf, size = 3.2) +
          labs(
            x = "Active (fraction of target signal covered)",
            y = "Unique target signal (sum)",
            title = "Active vs Unique — pathways (Deep)"
          ) +
          theme_minimal(base_size = 13)
        print(p)
        dev.off()
        
        
        learn_paths$deep_active_unique <- au_png
        
        # 2) Activity (gnn pathway score if present, else classic) vs Unique
        deep_scores <- tbl$pathway_n
        act <- setNames(as.numeric(deep_scores$score %||% deep_scores$aggregate), as.character(deep_scores$pathway))
        axu_png <- file.path(plots_www_dir, sprintf(
          "Deep_activity_unique_%s_to_%s.png", receiver_short, target_short))
        
        png(axu_png, width = 900, height = 600)
        deep_scores <- tbl$pathway_n
        act <- setNames(as.numeric(deep_scores$score %||% deep_scores$aggregate),
                        as.character(deep_scores$pathway))
        df <- data.frame(
          x = act[pw_cov$pathway],
          y = pw_cov$targ_unique_sum,
          label = pw_cov$pathway
        )
        df <- df[is.finite(df$x) & is.finite(df$y), ]
        
        p <- ggplot(df, aes(x, y)) +
          geom_point(size = 3) +
          ggrepel::geom_text_repel(aes(label = label), max.overlaps = Inf, size = 3.2) +
          labs(
            x = "Activity score (Deep)",
            y = "Unique target signal (sum)",
            title = "Activity vs Unique — pathways (Deep)"
          ) +
          theme_minimal(base_size = 13)
        print(p)
        dev.off()
        
        
        learn_paths$deep_activity_unique <- axu_png
        
        # Minimal cover chips (~80%)
        Rsets <- build_reachability(list(LigRec = netList$LigRec, RecTF = netList$RecTF, TFTar = netList$TFTar), receptor_ligand)
        w_tg  <- .target_mag_from_deg(deg)
        min_cover <- greedy_cover(Rsets$targets_by_pathway, w_tg, threshold = 0.8)
        
        
      }
      
      
      
      
      
      
      
      
      # ---- Publish std tables from top_pathway (with fallback) ----
      
      # Helper: return an empty canonical table
      .empty_canon <- function() {
        data.frame(
          kind = character(), entity = character(), pathway = character(),
          score = numeric(), score.perc = character(),
          p.value = numeric(), significant = character(),
          stringsAsFactors = FALSE
        )
      }
      
      # Helper: coerce any df to canonical columns/types and decorate % text
      .coerce_to_canon <- function(df, kind_name) {
        CANON <- c("kind","entity","pathway","score","score.perc","p.value","significant")
        if (is.null(df) || !nrow(df)) {
          out <- .empty_canon(); out$kind <- kind_name; return(out)
        }
        # If already canonical, just coerce types and format
        nm <- names(df)
        # Try to infer "entity" if it isn't present
        if (!"entity" %in% nm) {
          if ("ligand"   %in% nm) df$entity <- df$ligand
          else if ("receptor" %in% nm) df$entity <- df$receptor
          else if ("pathway"  %in% nm) df$entity <- df$pathway
          else if ("name"     %in% nm) df$entity <- df$name
          else df$entity <- NA_character_
        }
        if (!"score" %in% nm && "aggregate" %in% nm) df$score <- df$aggregate
        if (!"score.perc" %in% nm && "score_perc" %in% nm) df$score.perc <- df$score_perc
        if (!"p.value" %in% nm && "pvalue" %in% nm) df$p.value <- df$pvalue
        if (!"significant" %in% nm && "sig" %in% nm) df$significant <- df$sig
        
        df$kind        <- kind_name
        df$entity      <- as.character(df$entity)
        if (!"pathway" %in% names(df)) df$pathway <- NA_character_
        df$pathway     <- as.character(df$pathway)
        df$score       <- suppressWarnings(as.numeric(df$score))
        # Accept either raw % (0–100) or fraction (0–1); leave as-is if NA
        if ("score.perc" %in% names(df)) {
          sp <- suppressWarnings(as.numeric(df$score.perc))
          # Heuristic: if all non-NA values are <= 1, treat as fraction
          frac_like <- all(is.na(sp) | sp <= 1)
          if (frac_like) sp <- sp * 100
          df$score.perc <- sp
        } else {
          df$score.perc <- NA_real_
        }
        if ("p.value"   %in% names(df)) df$p.value   <- suppressWarnings(as.numeric(df$p.value)) else df$p.value <- NA_real_
        if (!"significant" %in% names(df)) df$significant <- NA_character_
        
        out <- df[, intersect(CANON, names(df)), drop = FALSE]
        for (nm in setdiff(CANON, names(out))) {
          out[[nm]] <- if (nm %in% c("score","score.perc","p.value")) NA_real_ else NA_character_
        }
        out <- out[CANON]
        # Pretty percent text
        out$score.perc <- ifelse(is.na(out$score.perc), NA, sprintf("%.1f%%", as.numeric(out$score.perc)))
        # Order by score desc when present
        suppressWarnings({
          ord <- order(-as.numeric(out$score))
          if (length(ord) && all(is.finite(ord))) out <- out[ord, , drop = FALSE]
        })
        out
      }
      
      # If top_pathway() emitted unified std tables, use them directly
      if (!is.null(res_norm$std) && is.list(res_norm$std)) {
        std_in <- res_norm$std
        std_tables_rv(list(
          pathway     = .coerce_to_canon(std_in$pathway,     "pathway"),
          ligand      = .coerce_to_canon(std_in$ligand,      "ligand"),
          receptor    = .coerce_to_canon(std_in$receptor,    "receptor"),
          interaction = .coerce_to_canon(std_in$interaction, "interaction")
        ))
      } else {
        # -------- Fallback: build canonical tables from legacy objects --------
        safe_df <- function(x) {
          if (is.null(x)) return(data.frame())
          y <- try(as.data.frame(x), silent = TRUE)
          if (inherits(y, "try-error")) data.frame() else y
        }
        
        # PATHWAY (summary + p)
        path_sum <- safe_df(res_norm$pathway_n)             # pathway, score, score.perc, score.perc.txt
        path_pv  <- safe_df(res_norm$pval_pathway)          # name, p.value, significant
        if (nrow(path_sum)) {
          path_tbl <- merge(
            transform(path_sum, pathway = as.character(pathway)),
            setNames(path_pv[c("name","p.value","significant")], c("pathway","p.value","significant")),
            by = "pathway", all.x = TRUE
          )
        } else path_tbl <- data.frame()
        
        # LIGAND (summary + p)
        lig_sum <- safe_df(res_norm$lig_rank_all)           # ligand, pathway, score, score.perc, score.perc.txt
        lig_pv  <- safe_df(res_norm$pval_ligand)            # name, p.value, significant
        if (nrow(lig_sum)) {
          lig_tbl <- merge(
            transform(lig_sum, ligand = as.character(ligand), pathway = as.character(pathway)),
            setNames(lig_pv[c("name","p.value","significant")], c("ligand","p.value","significant")),
            by = "ligand", all.x = TRUE
          )
        } else lig_tbl <- data.frame()
        
        # RECEPTOR (summary + p)
        rec_sum <- safe_df(res_norm$rec_rank_all)           # receptor, pathway, score, score.perc, score.perc.txt
        rec_pv  <- safe_df(res_norm$pval_receptor)          # name, p.value, significant
        if (nrow(rec_sum)) {
          # If receptor names got lost as rownames, restore
          if (!"receptor" %in% names(rec_sum)) {
            rn <- rownames(rec_sum)
            rec_sum$receptor <- if (!is.null(rn) && any(nzchar(rn))) rn else NA_character_
          }
          right <- if (nrow(rec_pv)) setNames(rec_pv[c("name","p.value","significant")], c("receptor","p.value","significant"))
          else data.frame(receptor=character(), p.value=numeric(), significant=character(), stringsAsFactors=FALSE)
          rec_tbl <- merge(
            transform(rec_sum, receptor = as.character(receptor), pathway = as.character(pathway)),
            right, by = "receptor", all.x = TRUE
          )
        } else rec_tbl <- data.frame()
        
        # Build canonical outputs (interaction not available in legacy path; emit empty)
        std_tables_rv(list(
          pathway     = .coerce_to_canon(path_tbl, "pathway"),
          ligand      = .coerce_to_canon(lig_tbl,  "ligand"),
          receptor    = .coerce_to_canon(rec_tbl,  "receptor"),
          interaction = .empty_canon()
        ))
      }
      
      
      
      
      
      
      # Copy the three expected PNGs to www/plots and update reactive paths
      normalization_label <- "db"
      
      need  <- png_names(receiver_short, target_short, normalization_label)
      srcs  <- file.path(tmp_outdir, need)
      dests <- file.path(plots_www_dir, need)
      
      to_copy <- file.exists(srcs)
      if (any(to_copy)) file.copy(srcs[to_copy], dests[to_copy], overwrite = TRUE)
      
      plot_paths$path             <- if (file.exists(dests[1])) dests[1] else NULL
      plot_paths$ligand           <- if (file.exists(dests[2])) dests[2] else NULL
      plot_paths$ligand_only      <- if (file.exists(dests[3])) dests[3] else NULL
      plot_paths$receptor_only    <- if (file.exists(dests[4])) dests[4] else NULL
      plot_paths$interaction      <- if (file.exists(dests[5])) dests[5] else NULL
      plot_paths$interaction_only <- if (file.exists(dests[6])) dests[6] else NULL
      
      
      
      # remember for any other modules that still look at last_params
      last_params$receiver_cell       <- receiver_cell
      last_params$diff_target_cell    <- diff_target_cell
      last_params$normalization_label <- normalization_label
      last_params$output_dir          <- plots_www_dir  # now we serve from www/plots
      
      # Debug info in the UI
      #existing_pngs <- list.files(plots_www_dir, pattern="\\.png$", full.names = TRUE)
      output$analysis_status <- renderUI({
        htmltools::tagList(
          HTML("<span style='color:#388e3c'><b>Analysis finished.</b></span>"))
        
      })
      
      
      
      
    }) # end withProgress
    
  })
  
  
  
  
  observeEvent(input$treemap_tabs, {
    req(analysis_result())  # make sure we actually have tables
    selected <- switch(input$treemap_tabs,
                       "Pathway"               = "pathway_n",
                       "Ligand by Pathway"     = "lig_rank_all",
                       "Ligand Only"           = "lig_rank_all",
                       "Receptor Only"         = "rec_rank_all",
                       "Interaction by Pathway"= "int_rank_all",   # NEW
                       "Interaction Only"      = "int_rank_all",   # NEW
                       "pathway_n"
    )
    updateSelectInput(session, "table_choice", selected = selected)
    
  }, ignoreInit = TRUE)
  
  
  output$pathwayPlot <- renderImage({
    req(plot_paths$path)
    list(src = plot_paths$path, contentType = "image/png", width = 600, height = 750)
  }, deleteFile = FALSE)
  
  output$ligandPlot <- renderImage({
    req(plot_paths$ligand)
    list(src = plot_paths$ligand, contentType = "image/png", width = 600, height = 750)
  }, deleteFile = FALSE)
  
  output$ligandOnlyPlot <- renderImage({
    req(plot_paths$ligand_only)
    list(src = plot_paths$ligand_only, contentType = "image/png", width = 600, height = 750)
  }, deleteFile = FALSE)
  
  output$receptorOnlyPlot <- renderImage({
    req(plot_paths$receptor_only)
    list(src = plot_paths$receptor_only, contentType = "image/png", width = 600, height = 750)
  }, deleteFile = FALSE)
  
  output$interactionPlot <- renderImage({
    req(plot_paths$interaction)
    list(src = plot_paths$interaction, contentType = "image/png", width = 600, height = 750)
  }, deleteFile = FALSE)
  
  output$interactionOnlyPlot <- renderImage({
    req(plot_paths$interaction_only)
    list(src = plot_paths$interaction_only, contentType = "image/png", width = 600, height = 750)
  }, deleteFile = FALSE)
  
  
  output$simple_scatter <- renderImage({
    req(learn_paths$simple_scatter)
    list(src = learn_paths$simple_scatter, contentType = "image/png", width = 800, height = 650)
  }, deleteFile = FALSE)
  
  output$simple_treemap <- renderImage({
    req(learn_paths$simple_cluster_treemap)
    list(src = learn_paths$simple_cluster_treemap, contentType = "image/png", width = 600, height = 750)
  }, deleteFile = FALSE)
  
  output$simple_pruned_treemap <- renderImage({
    req(learn_paths$simple_pruned_treemap)
    list(src = learn_paths$simple_pruned_treemap, contentType = "image/png", width = 600, height = 750)
  }, deleteFile = FALSE)
  
  # Simple tab Active/Unique
  output$simple_active_unique_img <- renderImage({
    req(learn_paths$simple_active_unique)
    list(src = learn_paths$simple_active_unique, contentType = "image/png")
  }, deleteFile = FALSE)
  
  output$simple_activity_unique_img <- renderImage({
    list(src = NULL)
  }, deleteFile = FALSE)
  
  
  # Deep tab Active/Unique
  output$deep_active_unique_img <- renderImage({
    req(learn_paths$deep_active_unique)
    list(src = learn_paths$deep_active_unique, contentType = "image/png")
  }, deleteFile = FALSE)
  
  
  
  
  # output$gnn_treemap <- renderImage({
  #   req(learn_paths$gnn_treemap)
  #   list(src = learn_paths$gnn_treemap, contentType = "image/png", width = 600, height = 750)
  # }, deleteFile = FALSE)
  
  output$true_simple_scatter <- renderImage({
    req(learn_paths$true_simple_scatter)
    list(src = learn_paths$true_simple_scatter, contentType = "image/png", width = 800, height = 650)
  }, deleteFile = FALSE)
  
  output$true_simple_treemap <- renderImage({
    req(learn_paths$true_simple_cluster_treemap)
    list(src = learn_paths$true_simple_cluster_treemap, contentType = "image/png", width = 600, height = 750)
  }, deleteFile = FALSE)
  
  output$true_simple_pruned_treemap <- renderImage({
    req(learn_paths$true_simple_pruned_treemap)
    list(src = learn_paths$true_simple_pruned_treemap, contentType = "image/png", width = 600, height = 750)
  }, deleteFile = FALSE)
  
  
  
  # For the 3 explanation subgraphs
  # --- helper: turn absolute disk path into a URL under /www ---
  # UI: a container that will hold N imageOutputs
  output$explain_gallery <- renderUI({
    paths <- learn_paths$explain_pngs
    if (length(paths) == 0L) return(tags$i("No explanation subgraphs available."))
    ids <- paste0("explain_img_", seq_along(paths))
    tagList(lapply(ids, function(id) imageOutput(id, height = "auto")))
  })
  
  # Server: create one renderImage per file (absolute paths work here)
  observe({
    paths <- learn_paths$explain_pngs
    if (length(paths) == 0L) return()
    ids <- paste0("explain_img_", seq_along(paths))
    lapply(seq_along(paths), function(i) {
      local({
        mypath <- paths[i]
        myid   <- ids[i]
        output[[myid]] <- renderImage({
          list(src = mypath, contentType = "image/png", alt = basename(mypath))
        }, deleteFile = FALSE)
      })
    })
  })
  
  
  output$simple_lr_pruned_treemap <- renderImage({
    req(learn_paths$simple_lr_pruned_treemap)
    list(
      src = learn_paths$simple_lr_pruned_treemap,
      contentType = "image/png",
      width = 600, height = 740
    )
  }, deleteFile = FALSE)
  
  output$deep_active_unique_img <- renderImage({
    list(src = learn_paths$deep_active_unique, contentType = "image/png")
  }, deleteFile = FALSE)
  
  output$deep_min_cover_ui <- renderUI({
    sel <- learn_misc$deep_min_cover_paths %||% character(0)
    if (!length(sel)) return(NULL)
    tags$div(
      tags$span("Minimal cover (≈80% targets): "),
      lapply(sel, function(p)
        tags$span(class = "badge bg-primary", style="margin-right:6px;", p))
    )
  })
  
  output$deep_active_unique_table <- renderTable({
    if (isTRUE(input$show_au_table)) learn_tables$deep_active_unique
  })
  
  output$true_gnn_treemap <- renderImage({
    req(learn_paths$true_gnn_treemap)
    list(src = learn_paths$true_gnn_treemap, contentType = "image/png", width = 600, height = 750)
  }, deleteFile = FALSE)
  
  output$true_explain_gallery <- renderUI({
    paths <- learn_paths$true_explain_pngs
    if (length(paths) == 0L) return(tags$i("No explanation subgraphs available from the True GNN."))
    ids <- paste0("true_explain_img_", seq_along(paths))
    tagList(lapply(ids, function(id) imageOutput(id, height = "auto")))
  })
  observe({
    paths <- learn_paths$true_explain_pngs
    if (length(paths) == 0L) return()
    ids <- paste0("true_explain_img_", seq_along(paths))
    lapply(seq_along(paths), function(i) {
      local({
        mypath <- paths[i]; myid <- ids[i]
        output[[myid]] <- renderImage({ list(src = mypath, contentType = "image/png") }, deleteFile = FALSE)
      })
    })
  })
  
  
  
  
  
  
  # ---- helpers ---------------------------------------------------------------
  
  # hardcode thresholds here
  LR_THRESH        <- 0
  MIN_EDGE_FRAC    <- 0          # you can raise later (e.g., 0.02)
  ALLOW_SELF_EDGES <- TRUE      # set TRUE if you want autocrine pairs included
  EXCLUDE_PURE_SELF_ONLY <- FALSE
  
  
  
  # angles for senders (all below the receiver)
  .sender_angles <- function(n_senders, arc = c(220, 320)) {
    if (n_senders <= 0) return(numeric(0))
    if (n_senders == 1) return(3*pi/2)                         # bottom
    if (n_senders == 3) return(c(240, 270, 300) * pi/180)      # keep well below
    seq(arc[1], arc[2], length.out = n_senders) * pi/180       # general case
  }
  
  # Grab int_rank_all from wherever you keep it (reactive or global)
  # ---- helpers: interaction table for talklr -----------------------------------
  .int_tbl <- reactive({
    # 1) your real results live here
    ar <- analysis_result()
    if (!is.null(ar)) {
      if (!is.null(ar$int_rank_all) && nrow(ar$int_rank_all)) {
        return(as.data.frame(ar$int_rank_all, stringsAsFactors = FALSE))
      }
      if (!is.null(ar$db) && !is.null(ar$db$int_rank_all) && nrow(ar$db$int_rank_all)) {
        return(as.data.frame(ar$db$int_rank_all, stringsAsFactors = FALSE))
      }
    }
    
    # 2) legacy globals (if present)
    if (exists("int_rank_all", inherits = TRUE)) {
      return(as.data.frame(get("int_rank_all", inherits = TRUE), stringsAsFactors = FALSE))
    }
    
    # 3) last-ditch: build choices from lr_glom_normal
    if (exists("lr_glom_normal", inherits = TRUE)) {
      lg <- as.data.frame(get("lr_glom_normal", inherits = TRUE), stringsAsFactors = FALSE)
      lc <- intersect(names(lg), c("interaction","ligand","Ligand","Ligand.ApprovedSymbol"))[1]
      rc <- intersect(names(lg), c("receptor","Receptor","Receptor.ApprovedSymbol"))[1]
      if (!is.na(lc) && lc != "interaction" && !is.na(rc)) {
        lg$interaction <- paste(lg[[lc]], lg[[rc]], sep = "_")
      }
      if ("interaction" %in% names(lg)) return(lg["interaction"])
    }
    
    NULL
  })
  
  
  # Make a robust "interaction" column (ligand_receptor) if needed
  .mk_interaction_col <- function(df) {
    df <- as.data.frame(df, stringsAsFactors = FALSE)
    if (!"interaction" %in% names(df)) {
      lc <- intersect(names(df), c("ligand","Ligand","Ligand.ApprovedSymbol"))[1]
      rc <- intersect(names(df), c("receptor","Receptor","Receptor.ApprovedSymbol"))[1]
      if (!is.na(lc) && !is.na(rc)) {
        df$interaction <- paste(df[[lc]], df[[rc]], sep = "_")
      }
    }
    df
  }
  
  # Find the row index in lr_glom_normal for a selected "LIGAND_RECEPTOR"
  .find_lr_row <- function(sel, lr_glom_normal) {
    # 1) direct "interaction" column
    if ("interaction" %in% names(lr_glom_normal)) {
      i <- match(sel, as.character(lr_glom_normal$interaction))
      if (!is.na(i)) return(i)
    }
    # 2) ligand + receptor columns
    lc <- intersect(names(lr_glom_normal), c("ligand","Ligand","Ligand.ApprovedSymbol"))[1]
    rc <- intersect(names(lr_glom_normal), c("receptor","Receptor","Receptor.ApprovedSymbol"))[1]
    if (!is.na(lc) && !is.na(rc)) {
      parts <- strsplit(sel, "_", fixed = TRUE)[[1]]
      if (length(parts) >= 2) {
        li <- toupper(parts[1]); ri <- toupper(paste(parts[-1], collapse = "_"))
        i <- which(toupper(lr_glom_normal[[lc]]) == li & toupper(lr_glom_normal[[rc]]) == ri)
        if (length(i)) return(i[1])
      }
    }
    # 3) rownames as "L_R"
    if (!is.null(rownames(lr_glom_normal))) {
      i <- match(sel, rownames(lr_glom_normal))
      if (!is.na(i)) return(i)
    }
    NA_integer_
  }
  
  # ---- populate the dropdown -------------------------------------------------
  # observeEvent(.int_tbl(), {
  #   df <- .int_tbl(); req(df)
  #   df <- .ensure_interaction_id(df)
  #   choices <- sort(unique(df$interaction_id[!is.na(df$interaction_id) & nzchar(df$interaction_id)]))
  #   updateSelectizeInput(session, "talklr_interaction", choices = choices, server = TRUE)
  # }, ignoreInit = FALSE)
  
  
  # ---- NEW: auto-detect ligand/receptor per-cell blocks in lr_glom_normal ----
  # Robustly locate the ligand/receptor per-cell blocks inside lr_glom_normal
  detect_lr_blocks <- function(lr, gl) {
    lr <- as.data.frame(lr, stringsAsFactors = FALSE)
    ncell <- ncol(gl) - 1L
    if (ncell <= 0) return(list(lig = integer(0), rec = integer(0)))
    
    cells <- colnames(gl)[-1]
    nm <- colnames(lr)
    
    # ---- A) Name-based: each cell name should occur twice (lig & rec), order = cells
    pos_list <- lapply(cells, function(s) which(toupper(nm) == toupper(s)))
    have_first  <- all(lengths(pos_list) >= 1)
    have_second <- all(lengths(pos_list) >= 2)
    
    if (have_first && have_second) {
      lig <- vapply(pos_list, function(ix) ix[1], integer(1))
      rec <- vapply(pos_list, function(ix) ix[2], integer(1))
      return(list(lig = lig, rec = rec))
    }
    
    # ---- B) Numeric fallback: take the *last* 2*ncell numeric columns (not necessarily contiguous)
    num_idx <- which(vapply(lr, function(x) is.numeric(x) || is.integer(x), logical(1)))
    if (length(num_idx) >= 2L * ncell) {
      last <- tail(num_idx, 2L * ncell)
      lig <- last[seq_len(ncell)]
      rec <- last[(ncell + 1):(2L * ncell)]
      return(list(lig = lig, rec = rec))
    }
    
    # ---- C) Legacy 17-based guess (clipped)
    lig <- 17:(17 + ncell - 1L); lig <- lig[lig <= ncol(lr)]
    rec <- (17 + ncell):(17 + 2L*ncell - 1L); rec <- rec[rec <= ncol(lr)]
    list(lig = lig, rec = rec)
  }
  
  
  coerce_lr_numeric_blocks <- function(lr, gl) {
    lr <- as.data.frame(lr, stringsAsFactors = FALSE)
    blk <- detect_lr_blocks(lr, gl)
    if (length(blk$lig) && length(blk$rec)) {
      lr[, blk$lig] <- lapply(lr[, blk$lig, drop = FALSE], function(x) as.numeric(as.character(x)))
      lr[, blk$rec] <- lapply(lr[, blk$rec, drop = FALSE], function(x) as.numeric(as.character(x)))
    }
    lr
  }
  
  
  restrict_lr_to_receiver <- function(
    lr, gl, receiver_id,
    expressed_thresh,
    min_edge_frac = 0.00,
    allow_self    = TRUE
  ) {
    lr  <- as.data.frame(lr, stringsAsFactors = FALSE)
    blk <- detect_lr_blocks(lr, gl)
    lig_cols <- blk$lig; rec_cols <- blk$rec
    
    if (!length(lig_cols) || !length(rec_cols))
      stop("Could not locate ligand/receptor per-cell blocks in lr_glom_normal.")
    
    cells <- colnames(gl)[-1]
    rcol  <- match(receiver_id, cells)
    if (is.na(rcol)) stop("Receiver '", receiver_id, "' not found in glom_normal.")
    
    # --- Coerce to numeric and zero-fill bad values (CRITICAL) ---
    L <- as.matrix(sapply(lr[, lig_cols, drop = FALSE], function(x) as.numeric(as.character(x))))
    R <- as.matrix(sapply(lr[, rec_cols, drop = FALSE], function(x) as.numeric(as.character(x))))
    L[!is.finite(L)] <- 0
    R[!is.finite(R)] <- 0
    
    # Threshold like talklr
    L[L < expressed_thresh] <- 0
    R[R < expressed_thresh] <- 0
    
    keep <- logical(nrow(lr))
    for (i in seq_len(nrow(lr))) {
      lig <- L[i, ]; rec <- R[i, ]
      if (!any(lig > 0) || !any(rec > 0)) next
      
      M <- outer(lig, rec, `*`)
      if (!allow_self) diag(M) <- 0
      
      # use na.rm=TRUE so we never feed NA to `if`
      s <- sum(M, na.rm = TRUE)
      if (s > 0 && min_edge_frac > 0) {
        M <- M / s
        M[M < min_edge_frac] <- 0
      }
      
      senders_to_r <- which(M[, rcol] > 0)
      if (length(senders_to_r)) {
        keep[i] <- TRUE
        if (length(senders_to_r) < length(lig)) L[i, -senders_to_r] <- 0
        if (ncol(R) > 1) R[i, -rcol] <- 0
      }
    }
    
    lr2 <- lr[keep, , drop = FALSE]
    if (nrow(lr2)) {
      lr2[, lig_cols] <- L[keep, , drop = FALSE]
      lr2[, rec_cols] <- R[keep, , drop = FALSE]
    }
    lr2
  }
  
  
  
  
  # ---- helper: build clean choices for the LR selector ----
  # Build the list of selectable interactions, keeping only those
  # that will actually render at least one edge under the same gating
  # used for plotting.
  make_interaction_choices <- function(
    lr, gl,
    thresh = 0,
    min_edge_frac = 0,
    allow_self = FALSE,
    exclude_pure_self_only = TRUE,
    fallback_if_empty = FALSE,
    receiver_id = NULL,                 # <--- NEW
    require_to_receiver = TRUE          # <--- NEW
  ) {
    lr <- as.data.frame(lr, stringsAsFactors = FALSE)
    
    # 0) interaction id if missing
    if (!"interaction" %in% names(lr)) {
      Lc <- intersect(names(lr), c("ligand","Ligand","Ligand.ApprovedSymbol"))[1]
      Rc <- intersect(names(lr), c("receptor","Receptor","Receptor.ApprovedSymbol"))[1]
      if (!is.na(Lc) && !is.na(Rc)) {
        lr$interaction <- paste(lr[[Lc]], lr[[Rc]], sep = "_")
      }
    }
    lr$interaction_id <- toupper(as.character(lr$interaction))
    lr$interaction_id[!nzchar(lr$interaction_id)] <- NA
    
    # 1) detect ligand/receptor per-cell blocks (robust, with legacy fallback)
    .detect_lr_blocks <- function(df, gl) {
      ncell <- ncol(gl) - 1L
      if (ncell <= 0) return(list(lig = integer(0), rec = integer(0)))
      num_idx <- which(vapply(df, function(x) is.numeric(x) || is.integer(x), logical(1)))
      if (length(num_idx) >= 2L * ncell) {
        for (start in rev(seq_len(ncol(df) - (2L*ncell) + 1L))) {
          window <- start:(start + 2L*ncell - 1L)
          if (all(window %in% num_idx)) {
            lig <- window[seq_len(ncell)]
            rec <- window[(ncell + 1):(2L*ncell)]
            return(list(lig = lig, rec = rec))
          }
        }
      }
      # legacy 17-based guess (clipped)
      lig <- 17:(17 + ncell - 1L); lig <- lig[lig <= ncol(df)]
      rec <- (17 + ncell):(17 + 2L*ncell - 1L); rec <- rec[rec <= ncol(df)]
      list(lig = lig, rec = rec)
    }
    
    idx <- .detect_lr_blocks(lr, gl)
    lig_col <- idx$lig; rec_col <- idx$rec
    if (!length(lig_col) || !length(rec_col)) {
      ch <- sort(unique(lr$interaction_id[!is.na(lr$interaction_id)]))
      return(ch)
    }
    
    cell_names <- colnames(gl)[-1]
    recv_col <- if (!is.null(receiver_id) && nzchar(receiver_id))
      match(receiver_id, cell_names) else NA_integer_
    
    # 2) keep only rows that will have at least one visible edge
    has_edges <- vapply(seq_len(nrow(lr)), function(i) {
      lig <- suppressWarnings(as.numeric(lr[i, lig_col])); lig[!is.finite(lig)] <- 0
      rec <- suppressWarnings(as.numeric(lr[i, rec_col])); rec[!is.finite(rec)] <- 0
      
      lig[lig < thresh] <- 0
      rec[rec < thresh] <- 0
      if (!any(lig > 0) || !any(rec > 0)) return(FALSE)
      
      m <- outer(lig, rec, `*`)
      if (!allow_self) diag(m) <- 0
      s <- sum(m); if (s <= 0) return(FALSE)
      
      if (min_edge_frac > 0) {
        m <- m / s
        m[m < min_edge_frac] <- 0
      }
      
      if (allow_self && exclude_pure_self_only) {
        off <- m; diag(off) <- 0
        if (!any(off > 0)) return(FALSE)
      }
      
      # <--- require at least one edge into the UI-selected receiver
      if (isTRUE(require_to_receiver) && !is.na(recv_col)) {
        if (!any(m[, recv_col] > 0)) return(FALSE)
      }
      any(m > 0)
    }, logical(1))
    
    choices <- sort(unique(lr$interaction_id[has_edges & !is.na(lr$interaction_id)]))
    
    if (!length(choices) && isTRUE(fallback_if_empty)) {
      choices <- sort(unique(lr$interaction_id[!is.na(lr$interaction_id)]))
    }
    choices
  }
  
  
  
  
  
  # ---- observer that feeds the selectize (drop this where you currently update it) ----
  observe({
    req(talklr_store$lr, talklr_store$glom)
    thr <- isolate(input$talklr_thresh %||% LR_THRESH)
    choices <- make_interaction_choices(
      talklr_store$lr, talklr_store$glom,
      thresh = input$talklr_thresh %||% 0,
      min_edge_frac = 0.00,
      allow_self = TRUE,
      exclude_pure_self_only = FALSE,
      receiver_id = isolate(input$receiver_cell),      # <--- NEW
      require_to_receiver = TRUE                       # <--- NEW
    )
    
    prev <- isolate(input$talklr_interaction)
    sel  <- if (!is.null(prev) && prev %in% choices) prev else (if (length(choices)) choices[1] else NULL)
    updateSelectizeInput(session, "talklr_interaction",
                         choices = choices, selected = sel, server = TRUE)
  })
  
  
  
  
  
  
  # --- talklr-lite: construct an LR graph safely ---
  make_lr_graph <- function(ligand_exprs,
                            receptor_exprs,
                            cell_labels,
                            thresh = 0,
                            min_edge_frac = 0.02,   # 2% of total mass (was 0.1 = 10%)
                            allow_self = FALSE) {
    if (!requireNamespace("igraph", quietly = TRUE)) {
      stop("Package 'igraph' is required.")
    }
    
    ligand_exprs   <- as.numeric(ligand_exprs)
    receptor_exprs <- as.numeric(receptor_exprs)
    stopifnot(length(ligand_exprs) == length(receptor_exprs),
              length(cell_labels)  == length(ligand_exprs))
    
    # threshold raw expressions
    ligand_exprs[ligand_exprs     < thresh] <- 0
    receptor_exprs[receptor_exprs < thresh] <- 0
    
    n <- length(cell_labels)
    # outer product: sender i -> receiver j
    mat <- outer(ligand_exprs, receptor_exprs, `*`)
    tot <- sum(mat)
    if (tot > 0) {
      mat <- mat / tot
      mat[mat < min_edge_frac] <- 0
    }
    if (!allow_self) diag(mat) <- 0
    
    dimnames(mat) <- list(cell_labels, cell_labels)
    
    g <- igraph::graph_from_adjacency_matrix(
      mat, mode = "directed", weighted = TRUE, diag = allow_self
    )
    # drop zero-weight edges (igraph keeps them)
    if (igraph::ecount(g)) {
      g <- igraph::delete_edges(g, which(igraph::E(g)$weight <= 0))
    }
    
    # set widths / labels (scale widths 1..8)
    if (igraph::ecount(g)) {
      w <- igraph::E(g)$weight
      igraph::E(g)$width <- 1 + 7 * (w / max(w))
      igraph::E(g)$label <- sprintf("%.2f", w)
    }
    g
  }
  
  # --- draw function (keeps par() local to this device) ---
  # 1) Make a single, robust interaction ID column for lr_glom_normal
  # make a single, correct interaction_id column
  .ensure_interaction_id <- function(df) {
    df <- as.data.frame(df, stringsAsFactors = FALSE)
    
    # if the data already has an 'interaction' column, trust it as-is
    if ("interaction" %in% names(df)) {
      df$interaction_id <- toupper(as.character(df$interaction))
      return(df)
    }
    
    # otherwise construct from ligand + receptor
    pick_lig <- intersect(names(df), c("ligand","Ligand","Ligand.ApprovedSymbol"))[1]
    pick_rec <- intersect(names(df), c("receptor","Receptor","Receptor.ApprovedSymbol"))[1]
    
    if (!is.na(pick_lig) && !is.na(pick_rec)) {
      L <- toupper(as.character(df[[pick_lig]]))
      R <- toupper(as.character(df[[pick_rec]]))
      df$interaction_id <- paste(L, R, sep = "_")
    } else {
      df$interaction_id <- NA_character_
    }
    df
  }
  
  
  # 2) Find the row for a selected ID
  .find_lr_row <- function(sel, df) {
    i <- match(toupper(sel), toupper(df$interaction_id))
    if (is.na(i)) NA_integer_ else i
  }
  
  # 3) Build an igraph from ligand/receiver vectors + threshold (base plotting)
  .plot_lr_graph <- function(lig, rec, labels, thresh, title_txt = NULL,
                             size_base = 24, arrow_base = 0.45, width_max = 10) {
    
    stopifnot(length(lig) == length(rec), length(labels) == length(lig))
    
    # If the threshold kills everything, relax it to the largest value
    # that still yields ≥1 edge (same logic as talklr: li>=t & rj>=t)
    mat_ok <- outer(lig >= thresh, rec >= thresh, "&")
    if (!any(mat_ok, na.rm = TRUE)) {
      # smallest max that keeps at least one TRUE
      t_l <- max(lig[is.finite(lig)], na.rm = TRUE)
      t_r <- max(rec[is.finite(rec)], na.rm = TRUE)
      t_auto <- min(t_l, t_r)
      if (is.finite(t_auto) && t_auto > 0) {
        thresh <- t_auto * 0.999  # just under the max so ≥1 edge survives
        mat_ok <- outer(lig >= thresh, rec >= thresh, "&")
      }
    }
    
    # adjacency by product, gated by threshold (identical to your talklr logic)
    adj <- outer(lig, rec, "*")
    adj[!mat_ok] <- 0
    
    # normalize weights for display scaling only
    w <- adj
    s <- sum(w)
    if (s > 0) w <- w / s
    
    dimnames(w) <- list(labels, labels)
    
    requireNamespace("igraph", quietly = TRUE)
    g <- igraph::graph_from_adjacency_matrix(w, mode = "directed", weighted = TRUE)
    
    # scales: keep thick edges but not comically so
    wv <- igraph::E(g)$weight
    if (length(wv) && any(wv > 0, na.rm = TRUE)) {
      wnorm <- wv / max(wv, na.rm = TRUE)
      igraph::E(g)$width <- pmax(1, width_max * sqrt(wnorm))
    } else {
      igraph::E(g)$width <- 1
    }
    
    co <- igraph::layout_in_circle(g)
    
    oldpar <- graphics::par(no.readonly = TRUE)
    on.exit(graphics::par(oldpar), add = TRUE)
    graphics::par(mar = c(0.3, 0.3, 2.5, 0.3))  # top room for title
    
    # IMPORTANT: use graphics::plot() (or igraph::plot.igraph), not igraph::plot
    graphics::plot(
      g,
      layout             = co,
      edge.curved        = TRUE,
      edge.arrow.size    = arrow_base,
      edge.color         = "black",
      edge.label         = if (length(wv)) round(wv, 2) else NA,
      edge.label.cex     = 1.2,
      edge.label.font    = 2,
      vertex.size        = size_base,
      vertex.color       = "white",
      vertex.label       = labels,
      vertex.label.cex   = 1.5,
      vertex.label.font  = 2,
      vertex.label.degree= -pi/2,
      vertex.label.dist  = 3
    )
    
    if (!is.null(title_txt) && nzchar(title_txt)) {
      graphics::title(main = title_txt, cex.main = 1.4, line = 0.5)
    }
    
    invisible(g)
  }
  
  
  #option B
  library(ggplot2)
  
  # ---- wiring plot (ggplot) ----
  
  
  library(ggplot2)
  library(grid)   # unit()
  # ---- option B (clean ggplot wiring) -----------------------------------------
  suppressPackageStartupMessages({
    library(ggplot2)
  })
  
  # Draw a clean LR wiring diagram with curved edges and safe self-loops.
  lr_wiring_ggplot <- function(
    ligand_exprs, receptor_exprs, cell_labels,
    thresh = 0, min_edge_frac = 0.02, allow_self = FALSE,
    receiver_col = "#97c9ea", sender_col = "white",
    r_receiver = 0.20, r_sender = 0.15,
    edge_max = 4, base_curv = 0.35, title_txt = NULL,
    receiver_id = NULL
  ) {
    if (!requireNamespace("ggplot2", quietly = TRUE)) stop("ggplot2 required")
    if (!requireNamespace("grid", quietly = TRUE)) stop("grid required")
    have_ggforce <- requireNamespace("ggforce", quietly = TRUE)
    
    ligand_exprs   <- as.numeric(ligand_exprs)
    receptor_exprs <- as.numeric(receptor_exprs)
    stopifnot(length(ligand_exprs) == length(receptor_exprs),
              length(cell_labels)  == length(ligand_exprs))
    
    ## threshold + product (talklr-like gating)
    lig <- ligand_exprs; lig[lig < thresh] <- 0
    rec <- receptor_exprs; rec[rec < thresh] <- 0
    mat <- outer(lig, rec, `*`)
    
    tot <- sum(mat)
    if (tot > 0) {
      mat <- mat / tot
      if (min_edge_frac > 0) mat[mat < min_edge_frac] <- 0
    }
    if (!allow_self) diag(mat) <- 0
    
    n <- length(cell_labels)
    if (n == 0) {
      return(list(plot = ggplot2::ggplot() + ggplot2::theme_void(), has_edges = FALSE))
    }
    
    ## receiver at top, others around the *lower* arc
    recv_idx <- if (!is.null(receiver_id) && nzchar(receiver_id))
      match(toupper(receiver_id), toupper(cell_labels)) else NA_integer_
    if (is.na(recv_idx)) recv_idx <- 1L
    
    n <- length(cell_labels)
    theta <- numeric(n)
    theta[recv_idx] <- pi/2
    others <- setdiff(seq_len(n), recv_idx)
    theta[others] <- .sender_angles(length(others))   # <— new
    
    nodes <- data.frame(
      id   = cell_labels,
      x0   = cos(theta),
      y0   = sin(theta),
      r    = ifelse(seq_len(n) == recv_idx, r_receiver, r_sender),
      fill = ifelse(seq_len(n) == recv_idx, receiver_col, sender_col),
      stringsAsFactors = FALSE
    )
    
    ## push senders a bit lower to avoid the receiver (tweak 0.00–0.15 as you like)
    if (length(others)) nodes$y0[others] <- nodes$y0[others] - 0.08
    
    
    ## collect edges
    keep <- which(mat > 0, arr.ind = TRUE)
    if (!nrow(keep)) {
      p <- ggplot2::ggplot() +
        ggplot2::coord_equal(xlim = c(-1.6, 1.6), ylim = c(-1.3, 1.3), expand = FALSE) +
        ggplot2::theme_void() + ggplot2::ggtitle(title_txt)
      return(list(plot = p, has_edges = FALSE))
    }
    
    edges <- data.frame(
      i = keep[,1], j = keep[,2],
      from = cell_labels[keep[,1]],
      to   = cell_labels[keep[,2]],
      w    = as.numeric(mat[keep]),
      stringsAsFactors = FALSE
    )
    edges$size <- 0.5 + (edge_max - 0.5) * edges$w / max(edges$w)
    
    edges_self <- edges[edges$i == edges$j, , drop = FALSE]
    edges_xy   <- edges[edges$i != edges$j, , drop = FALSE]
    
    ## --- geometry for non-self edges (no overlap; default inward) --------------
    if (nrow(edges_xy)) {
      sx <- nodes$x0[edges_xy$i]; sy <- nodes$y0[edges_xy$i]
      tx <- nodes$x0[edges_xy$j]; ty <- nodes$y0[edges_xy$j]
      rs <- nodes$r[edges_xy$i];  rt <- nodes$r[edges_xy$j]
      
      dx <- tx - sx; dy <- ty - sy
      d  <- sqrt(dx*dx + dy*dy); d[d == 0] <- 1e-9
      ux <- dx/d;  uy <- dy/d
      
      tip_pad <- 0.02
      edges_xy$xs <- sx + ux * rs
      edges_xy$ys <- sy + uy * rs
      edges_xy$xe <- tx - ux * (rt + tip_pad)
      edges_xy$ye <- ty - uy * (rt + tip_pad)
      
      ## decide "inward" curvature sign per directed edge
      ## positive curvature in geom_curve means "bend to the left"
      inward_sign_for <- function(ii, jj) {
        sx <- nodes$x0[ii]; sy <- nodes$y0[ii]
        tx <- nodes$x0[jj]; ty <- nodes$y0[jj]
        mx <- (sx + tx)/2; my <- (sy + ty)/2          # midpoint vector
        vx <- (tx - sx);  vy <- (ty - sy)             # chord vector
        nLx <- -vy;      nLy <-  vx                   # left-hand normal
        # if left points outward (dot>0), inward is to the right => -1
        if (mx * nLx + my * nLy > 0) -1 else 1
      }
      edges_xy$in_sign <- mapply(inward_sign_for, edges_xy$i, edges_xy$j)
      
      ## group by unordered pairs
      pair_key <- paste(pmin(edges_xy$i, edges_xy$j),
                        pmax(edges_xy$i, edges_xy$j), sep = "_")
      edges_xy$curv_sign <- 0
      
      for (k in unique(pair_key)) {
        idx <- which(pair_key == k)
        if (length(idx) == 1) {
          ## single direction -> bend inward by default
          edges_xy$curv_sign[idx] <- edges_xy$in_sign[idx]
        } else {
          ## two directions -> thicker inward, thinner outward
          ord <- idx[order(edges_xy$size[idx], decreasing = TRUE)]
          hi  <- ord[1]; lo <- ord[-1]
          edges_xy$curv_sign[hi] <- edges_xy$in_sign[hi]      # inward
          edges_xy$curv_sign[lo] <- -edges_xy$in_sign[lo]     # outward (opposite side)
        }
      }
    }
    
    ## --- geometry for self-loops -----------------------------------------------
    if (nrow(edges_self)) {
      i  <- edges_self$i
      ph <- atan2(nodes$y0[i], nodes$x0[i])
      loop_r <- nodes$r[i] + 0.08
      delta  <- 0.60
      
      xs <- nodes$x0[i] + loop_r * cos(ph - delta)
      ys <- nodes$y0[i] + loop_r * sin(ph - delta)
      xe <- nodes$x0[i] + loop_r * cos(ph + delta)
      ye <- nodes$y0[i] + loop_r * sin(ph + delta)
      
      edges_self$xs <- xs; edges_self$ys <- ys; edges_self$xe <- xe; edges_self$ye <- ye
    }
    
    .drop_identical <- function(df) {
      if (!nrow(df)) return(df)
      df[(abs(df$xs - df$xe) + abs(df$ys - df$ye)) > 1e-12, , drop = FALSE]
    }
    edges_xy   <- .drop_identical(edges_xy)
    edges_self <- .drop_identical(edges_self)
    
    ## --- plotting ---------------------------------------------------------------
    p <- ggplot2::ggplot()
    
    if (have_ggforce) {
      p <- p +
        ggforce::geom_circle(
          data = nodes, ggplot2::aes(x0 = x0, y0 = y0, r = r, fill = fill),
          color = "black", linewidth = 1
        ) +
        ggplot2::scale_fill_identity()
    } else {
      p <- p +
        ggplot2::geom_point(
          data = nodes, ggplot2::aes(x0, y0, fill = I(fill)),
          shape = 21, size = 95 * nodes$r, color = "black", stroke = 1
        )
    }
    
    p <- p + ggplot2::geom_text(data = nodes, ggplot2::aes(x0, y0, label = id),
                                fontface = 2, size = 6)
    
    if (nrow(edges_xy)) {
      e_pos <- edges_xy[edges_xy$curv_sign > 0, , drop = FALSE]
      e_neg <- edges_xy[edges_xy$curv_sign < 0, , drop = FALSE]
      
      if (nrow(e_pos)) {
        p <- p + ggplot2::geom_curve(
          data = e_pos,
          ggplot2::aes(x = xs, y = ys, xend = xe, yend = ye, linewidth = size),
          curvature =  base_curv,
          arrow = grid::arrow(type = "closed", length = grid::unit(0.18, "in")),
          lineend = "round", colour = "black",
          show.legend = FALSE
        )
      }
      if (nrow(e_neg)) {
        p <- p + ggplot2::geom_curve(
          data = e_neg,
          ggplot2::aes(x = xs, y = ys, xend = xe, yend = ye, linewidth = size),
          curvature = -base_curv,
          arrow = grid::arrow(type = "closed", length = grid::unit(0.18, "in")),
          lineend = "round", colour = "black",
          show.legend = FALSE
        )
      }
    }
    
    if (nrow(edges_self)) {
      p <- p + ggplot2::geom_curve(
        data = edges_self,
        ggplot2::aes(x = xs, y = ys, xend = xe, yend = ye, linewidth = size),
        curvature =  base_curv,
        arrow = grid::arrow(type = "closed", length = grid::unit(0.18, "in")),
        lineend = "round", colour = "black",
        show.legend = FALSE
      )
    }
    
    p <- p +
      ggplot2::scale_linewidth_identity() +
      ggplot2::coord_equal(xlim = c(-1.6, 1.6), ylim = c(-1.3, 1.38), expand = FALSE) +
      ggplot2::theme_void() +
      ggplot2::theme(
        legend.position = "none",
        plot.title.position = "plot",
        plot.title = ggplot2::element_text(hjust = 0.5, face = "bold", size = 20,
                                           margin = ggplot2::margin(b = 10))
      ) +
      ggplot2::ggtitle(title_txt)
    
    
    
    list(plot = p, has_edges = TRUE)
  }
  
  
  
  
  
  
  
  
  
  
  
  
  
  
  
  
  
  
  #---- plot wiring on selection ----------------------------------------------
  # # hardcode thresholds here
  # LR_THRESH      <- 0     # raw gating on ligand/receptor
  # MIN_EDGE_FRAC  <- 0     # keep all non-zero edges
  # 
  output$talklr_plot <- renderPlot({
    req(input$talklr_interaction, talklr_store$lr, talklr_store$glom)
    lr_glom_normal <- talklr_store$lr
    glom_normal    <- talklr_store$glom
    cols <- detect_lr_blocks(lr_glom_normal, glom_normal)
    lig_col <- cols$lig; rec_col <- cols$rec
    
    lr <- as.data.frame(lr_glom_normal)
    if (!"interaction" %in% names(lr)) {
      Lc <- intersect(names(lr), c("ligand","Ligand","Ligand.ApprovedSymbol"))[1]
      Rc <- intersect(names(lr), c("receptor","Receptor","Receptor.ApprovedSymbol"))[1]
      lr$interaction <- paste(lr[[Lc]], lr[[Rc]], sep = "_")
    }
    lr$interaction_id <- toupper(as.character(lr$interaction))
    idx <- match(toupper(input$talklr_interaction), lr$interaction_id)
    validate(need(!is.na(idx), sprintf("Interaction '%s' not found.", input$talklr_interaction)))
    
    lig_vec <- suppressWarnings(as.numeric(lr[idx, lig_col, drop = TRUE])); lig_vec[!is.finite(lig_vec)] <- 0
    rec_vec <- suppressWarnings(as.numeric(lr[idx, rec_col, drop = TRUE])); rec_vec[!is.finite(rec_vec)] <- 0
    labs    <- colnames(glom_normal)[-1]
    
    # Receiver highlight comes ONLY from the UI selection:
    recv_lbl <- label_for(input$receiver_cell)
    if (!recv_lbl %in% labs) recv_lbl <- input$receiver_cell  # fall back to full name if needed
    
    res <- lr_wiring_ggplot(
      ligand_exprs   = lig_vec,
      receptor_exprs = rec_vec,
      cell_labels    = labs,
      thresh         = isolate(input$talklr_thresh %||% LR_THRESH),
      min_edge_frac  = MIN_EDGE_FRAC,
      allow_self     = ALLOW_SELF_EDGES,
      title_txt      = input$talklr_interaction,
      receiver_id    = recv_lbl   # <-- highlight exactly the receiver chosen in the UI
      #default_clockwise = FALSE   # <-- restores your previous default bend
    )
    validate(need(res$has_edges, NULL))
    print(res$plot)
  })
  
  
  
  outputOptions(output, "talklr_plot", suspendWhenHidden = FALSE)
  
  
  
  ########## talklr pathway
  
  
  
  # ensure we have: interaction_id (L_R) and pathway_std (string)
  .standardize_lr <- function(lr) {
    lr <- as.data.frame(lr, stringsAsFactors = FALSE)
    
    # interaction_id
    if (!"interaction" %in% names(lr)) {
      Lc <- intersect(names(lr), c("ligand","Ligand","Ligand.ApprovedSymbol"))[1]
      Rc <- intersect(names(lr), c("receptor","Receptor","Receptor.ApprovedSymbol"))[1]
      if (!is.na(Lc) && !is.na(Rc)) lr$interaction <- paste(lr[[Lc]], lr[[Rc]], sep = "_")
    }
    lr$interaction_id <- toupper(as.character(lr$interaction))
    
    # pathway
    pw_col <- intersect(names(lr), c("pathway_name","pathway","Pathway","pathwayName"))[1]
    if (!is.na(pw_col)) {
      lr$pathway_std <- as.character(lr[[pw_col]])
    } else {
      # fall back to receptor_ligand join
      Lc <- intersect(names(lr), c("ligand","Ligand","Ligand.ApprovedSymbol"))[1]
      Rc <- intersect(names(lr), c("receptor","Receptor","Receptor.ApprovedSymbol"))[1]
      if (!is.na(Lc) && !is.na(Rc) && exists("receptor_ligand", inherits = TRUE)) {
        rl <- get("receptor_ligand", inherits = TRUE)
        key_lr <- paste(toupper(lr[[Lc]]), toupper(lr[[Rc]]), sep = "_")
        rl$key <- paste(toupper(rl$Ligand.ApprovedSymbol), toupper(rl$Receptor.ApprovedSymbol), sep = "_")
        m <- match(key_lr, rl$key)
        lr$pathway_std <- as.character(rl$pathway_name[m])
      } else {
        lr$pathway_std <- NA_character_
      }
    }
    lr
  }
  
  # the set of interactions that would draw *right now* in the talklr tab
  talklr_allowed_interactions <- reactive({
    req(talklr_store$lr, talklr_store$glom)
    make_interaction_choices(
      talklr_store$lr, talklr_store$glom,
      thresh = input$talklr_thresh %||% 0,
      min_edge_frac = 0.00,
      allow_self = TRUE,
      exclude_pure_self_only = FALSE,
      receiver_id = isolate(input$receiver_cell),
      require_to_receiver = TRUE
    )
  })
  
  
  # draw-from-matrix version of the wiring plot (same geometry/logic)
  wiring_from_matrix <- function(
    mat, cell_labels,
    receiver_col = "#97c9ea", sender_col = "white",
    r_receiver = 0.20, r_sender = 0.15,
    edge_max = 3.5, base_curv = 0.35, title_txt = NULL,
    receiver_id = NULL
  ) {
    stopifnot(length(cell_labels) == nrow(mat), ncol(mat) == length(cell_labels))
    # reuse the body of lr_wiring_ggplot *starting* from where `mat` exists.
    # The code below mirrors the second half of lr_wiring_ggplot verbatim.
    
    n <- length(cell_labels)
    recv_idx <- if (!is.null(receiver_id) && nzchar(receiver_id))
      match(toupper(receiver_id), toupper(cell_labels)) else NA_integer_
    if (is.na(recv_idx)) recv_idx <- 1L
    
    n <- length(cell_labels)
    theta <- numeric(n)
    theta[recv_idx] <- pi/2
    others <- setdiff(seq_len(n), recv_idx)
    theta[others] <- .sender_angles(length(others))   # <— new
    
    nodes <- data.frame(
      id   = cell_labels,
      x0   = cos(theta),
      y0   = sin(theta),
      r    = ifelse(seq_len(n) == recv_idx, r_receiver, r_sender),
      fill = ifelse(seq_len(n) == recv_idx, receiver_col, sender_col),
      stringsAsFactors = FALSE
    )
    
    if (length(others)) nodes$y0[others] <- nodes$y0[others] - 0.08
    
    
    keep <- which(mat > 0, arr.ind = TRUE)
    if (!nrow(keep)) {
      p <- ggplot2::ggplot() +
        ggplot2::coord_equal(xlim = c(-1.6, 1.6), ylim = c(-1.3, 1.38), expand = FALSE) +
        ggplot2::theme_void() + ggplot2::ggtitle(title_txt)
      return(list(plot = p, has_edges = FALSE))
    }
    
    edges <- data.frame(
      i = keep[,1], j = keep[,2],
      from = cell_labels[keep[,1]],
      to   = cell_labels[keep[,2]],
      w    = as.numeric(mat[keep]),
      stringsAsFactors = FALSE
    )
    edges$size <- 0.5 + (edge_max - 0.5) * edges$w / max(edges$w)
    
    edges_self <- edges[edges$i == edges$j, , drop = FALSE]
    edges_xy   <- edges[edges$i != edges$j, , drop = FALSE]
    
    if (nrow(edges_xy)) {
      sx <- nodes$x0[edges_xy$i]; sy <- nodes$y0[edges_xy$i]
      tx <- nodes$x0[edges_xy$j]; ty <- nodes$y0[edges_xy$j]
      rs <- nodes$r[edges_xy$i];  rt <- nodes$r[edges_xy$j]
      
      dx <- tx - sx; dy <- ty - sy
      d  <- sqrt(dx*dx + dy*dy); d[d == 0] <- 1e-9
      ux <- dx/d;  uy <- dy/d
      
      tip_pad <- 0.02
      edges_xy$xs <- sx + ux * rs
      edges_xy$ys <- sy + uy * rs
      edges_xy$xe <- tx - ux * (rt + tip_pad)
      edges_xy$ye <- ty - uy * (rt + tip_pad)
      
      inward_sign_for <- function(ii, jj) {
        sx <- nodes$x0[ii]; sy <- nodes$y0[ii]
        tx <- nodes$x0[jj]; ty <- nodes$y0[jj]
        mx <- (sx + tx)/2; my <- (sy + ty)/2
        vx <- (tx - sx);  vy <- (ty - sy)
        nLx <- -vy;       nLy <-  vx
        if (mx * nLx + my * nLy > 0) -1 else 1
      }
      edges_xy$in_sign <- mapply(inward_sign_for, edges_xy$i, edges_xy$j)
      
      pair_key <- paste(pmin(edges_xy$i, edges_xy$j),
                        pmax(edges_xy$i, edges_xy$j), sep = "_")
      edges_xy$curv_sign <- 0
      for (k in unique(pair_key)) {
        idx <- which(pair_key == k)
        if (length(idx) == 1) {
          edges_xy$curv_sign[idx] <- edges_xy$in_sign[idx]
        } else {
          ord <- idx[order(edges_xy$size[idx], decreasing = TRUE)]
          hi  <- ord[1]; lo <- ord[-1]
          edges_xy$curv_sign[hi] <-  edges_xy$in_sign[hi]
          edges_xy$curv_sign[lo] <- -edges_xy$in_sign[lo]
        }
      }
    }
    
    if (nrow(edges_self)) {
      i  <- edges_self$i
      ph <- atan2(nodes$y0[i], nodes$x0[i])
      loop_r <- nodes$r[i] + 0.08
      delta  <- 0.60
      xs <- nodes$x0[i] + loop_r * cos(ph - delta)
      ys <- nodes$y0[i] + loop_r * sin(ph - delta)
      xe <- nodes$x0[i] + loop_r * cos(ph + delta)
      ye <- nodes$y0[i] + loop_r * sin(ph + delta)
      edges_self$xs <- xs; edges_self$ys <- ys; edges_self$xe <- xe; edges_self$ye <- ye
    }
    
    .drop_identical <- function(df) {
      if (!nrow(df)) return(df)
      df[(abs(df$xs - df$xe) + abs(df$ys - df$ye)) > 1e-12, , drop = FALSE]
    }
    edges_xy   <- .drop_identical(edges_xy)
    edges_self <- .drop_identical(edges_self)
    
    p <- ggplot2::ggplot()
    
    if (requireNamespace("ggforce", quietly = TRUE)) {
      p <- p +
        ggforce::geom_circle(
          data = nodes, ggplot2::aes(x0 = x0, y0 = y0, r = r, fill = fill),
          color = "black", linewidth = 1
        ) +
        ggplot2::scale_fill_identity()
    } else {
      p <- p +
        ggplot2::geom_point(
          data = nodes, ggplot2::aes(x0, y0, fill = I(fill)),
          shape = 21, size = 95 * nodes$r, color = "black", stroke = 1
        )
    }
    
    p <- p + ggplot2::geom_text(data = nodes, ggplot2::aes(x0, y0, label = id),
                                fontface = 2, size = 6)
    
    if (nrow(edges_xy)) {
      e_pos <- edges_xy[edges_xy$curv_sign > 0, , drop = FALSE]
      e_neg <- edges_xy[edges_xy$curv_sign < 0, , drop = FALSE]
      if (nrow(e_pos)) {
        p <- p + ggplot2::geom_curve(
          data = e_pos,
          ggplot2::aes(x = xs, y = ys, xend = xe, yend = ye, linewidth = size),
          curvature =  base_curv,
          arrow = grid::arrow(type = "closed", length = grid::unit(0.18, "in")),
          lineend = "round", colour = "black",
          show.legend = FALSE
        )
      }
      if (nrow(e_neg)) {
        p <- p + ggplot2::geom_curve(
          data = e_neg,
          ggplot2::aes(x = xs, y = ys, xend = xe, yend = ye, linewidth = size),
          curvature = -base_curv,
          arrow = grid::arrow(type = "closed", length = grid::unit(0.18, "in")),
          lineend = "round", colour = "black",
          show.legend = FALSE
        )
      }
    }
    
    if (nrow(edges_self)) {
      p <- p + ggplot2::geom_curve(
        data = edges_self,
        ggplot2::aes(x = xs, y = ys, xend = xe, yend = ye, linewidth = size),
        curvature =  base_curv,
        arrow = grid::arrow(type = "closed", length = grid::unit(0.18, "in")),
        lineend = "round", colour = "black",
        show.legend = FALSE
      )
    }
    
    p <- p +
      ggplot2::scale_size_identity() +
      ggplot2::coord_equal(xlim = c(-1.6, 1.6), ylim = c(-1.3, 1.38), expand = FALSE) +
      ggplot2::theme_void() +
      ggplot2::theme(
        plot.title.position = "plot",
        plot.title = ggplot2::element_text(
          hjust = 0.5, face = "bold", size = 20,
          margin = ggplot2::margin(b = 10)
        )
      ) +
      ggplot2::ggtitle(title_txt)
    
    list(plot = p, has_edges = TRUE)
  }
  
  # Sum edge matrices across the interactions belonging to one pathway.
  aggregate_pathway_matrix <- function(
    lr, gl, pathway, thresh = 0, min_edge_frac = 0,
    allow_self = TRUE, receiver_id = NULL, allowed_interactions = NULL
  ) {
    lr <- .mk_interaction_col(lr)
    idx <- detect_lr_blocks(lr, gl)
    lig_col <- idx$lig; rec_col <- idx$rec
    if (!length(lig_col) || !length(rec_col)) return(NULL)
    
    # pathway column (robust pick)
    pw_col <- intersect(names(lr), c("pathway_name","pathway","Pathway","pathwayName"))[1]
    if (is.na(pw_col)) return(NULL)
    
    # keep rows in selected pathway + restricted to allowed interactions (if provided)
    sub <- lr[ toupper(lr[[pw_col]]) == toupper(pathway), , drop = FALSE]
    if (!is.null(allowed_interactions)) {
      sub <- sub[toupper(sub$interaction) %in% toupper(allowed_interactions) |
                   toupper(sub$interaction_id) %in% toupper(allowed_interactions), , drop = FALSE]
    }
    if (!nrow(sub)) return(NULL)
    
    cells <- colnames(gl)[-1]
    M <- matrix(0, nrow = length(cells), ncol = length(cells))
    for (i in seq_len(nrow(sub))) {
      lig <- suppressWarnings(as.numeric(sub[i, lig_col])); lig[!is.finite(lig)] <- 0
      rec <- suppressWarnings(as.numeric(sub[i, rec_col])); rec[!is.finite(rec)] <- 0
      lig[lig < thresh] <- 0
      rec[rec < thresh] <- 0
      if (any(lig > 0) && any(rec > 0)) {
        M <- M + outer(lig, rec, `*`)
      }
    }
    if (!allow_self) diag(M) <- 0
    
    # normalize & prune tiny edges relative to the aggregated total
    tot <- sum(M)
    if (tot > 0) {
      M <- M / tot
      if (min_edge_frac > 0) M[M < min_edge_frac] <- 0
    }
    dimnames(M) <- list(cells, cells)
    M
  }
  
  
  observe({
    req(talklr_store$lr)
    
    lr <- .standardize_lr(talklr_store$lr)
    allowed_ids <- talklr_allowed_interactions()
    
    keep_rows <- which(toupper(lr$interaction_id) %in% toupper(allowed_ids))
    paths <- sort(unique(lr$pathway_std[keep_rows]))
    paths <- paths[nzchar(paths)]
    
    # also respect the top-level "Pathways Included" if present
    sel_paths <- tryCatch(selected_paths(), error = function(e) character(0))
    if (length(sel_paths)) paths <- intersect(paths, sel_paths)
    
    updateSelectizeInput(
      session, "talklr_pathway",
      choices = paths,
      selected = if (length(paths)) paths[1] else NULL,
      server = TRUE
    )
  })
  
  
  
  
  output$talklr_pathway_plot <- renderPlot({
    req(input$talklr_pathway, talklr_store$lr, talklr_store$glom)
    lr  <- .standardize_lr(talklr_store$lr)
    gl  <- talklr_store$glom
    ids <- talklr_allowed_interactions()
    validate(need(length(ids) > 0, "No selectable interactions under current settings."))
    
    # limit rows: selected pathway + currently-allowed interactions
    rows <- which(toupper(lr$pathway_std) == toupper(input$talklr_pathway) &
                    toupper(lr$interaction_id) %in% toupper(ids))
    validate(need(length(rows) > 0, "No interactions in this pathway with the current settings."))
    
    # detect ligand/receptor blocks once
    idx <- detect_lr_blocks(lr, gl)
    lig_col <- idx$lig; rec_col <- idx$rec
    validate(need(length(lig_col) && length(rec_col), "Could not locate ligand/receptor columns."))
    
    cells <- colnames(gl)[-1]
    M <- matrix(0, length(cells), length(cells), dimnames = list(cells, cells))
    
    thr <- input$talklr_thresh %||% 0
    for (i in rows) {
      lig <- suppressWarnings(as.numeric(lr[i, lig_col])); lig[!is.finite(lig)] <- 0
      rec <- suppressWarnings(as.numeric(lr[i, rec_col])); rec[!is.finite(rec)] <- 0
      lig[lig < thr] <- 0; rec[rec < thr] <- 0
      if (any(lig > 0) && any(rec > 0)) M <- M + outer(lig, rec, `*`)
    }
    
    # optional: keep self-edges; if you want to drop them, uncomment next line
    # diag(M) <- 0
    
    validate(need(sum(M) > 0, "No edges available for this pathway with the current settings."))
    
    # normalize only for display thickness (don’t zero out small ones)
    M <- M / sum(M)
    
    labs <- cells
    recv_lbl <- label_for(input$receiver_cell); if (!recv_lbl %in% labs) recv_lbl <- input$receiver_cell
    
    res <- wiring_from_matrix(
      mat = M, cell_labels = labs,
      title_txt  = input$talklr_pathway,
      receiver_id= recv_lbl
    )
    print(res$plot)
  })
  
  
  
  
  
  
  
  
  
  
  
  
  
  
  
}


shinyApp(ui = ui, server = server)

