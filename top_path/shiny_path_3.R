options(shiny.maxRequestSize = 500*1024^2)
setwd("C:/Users/afrit/Documents/R_work")

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

# Pathways file
pathway_choices <- readLines("top_path/pathways.txt")
pathway_choices_with_all <- c("All", pathway_choices)

out.file <- "top_path_output"
user.file <- "shiny_path"
meta.file <- paste0(out.file, "/", user.file, "/metadata")
CellChatDB <- readRDS("database/CellChatDB.rds")
receptor_ligand <- read.csv("database/receptor_ligand.csv")
receptor_ligand <- receptor_ligand[ , -14]


# Source your function files
source("top_path/top_path_function_2.0.R")
source("top_path/Run_scMLnet.R")
source("top_path/Draw_MLnet.R")
source("top_path/plot_lr_wiring.R")
source("top_path/convert_h5ad_to_cds.R")

ui <- fluidPage(
  tags$head(tags$style(HTML("
  .tab-content { padding-top: 8px; }
  .dataTables_wrapper { overflow-x: auto; }
"))),
  
  titlePanel("Top Pathway Analysis:"),
  sidebarLayout(
    sidebarPanel(
      fileInput("cdsfile", "Upload Monocle CDS (.rds) or Scanpy (.h5ad)", accept = c(".rds", ".h5ad")),
      uiOutput("filetype_status"),
      selectizeInput("clustering_col", "Choose clustering column", choices = NULL, selected = NULL, options = list(placeholder = "Select clustering column...")),
      selectizeInput("receiver_cell", "Receiver Cell", choices = NULL, selected = NULL, options = list(placeholder = "Select receiver cell...")),
      selectizeInput("diff_target_cell", "Diff Target Cell", choices = NULL, selected = NULL, options = list(placeholder = "Select diff target cell...")),
      selectizeInput("sender_cells", "Sender Cells (multiple allowed)", choices = NULL, selected = NULL, multiple = TRUE, options = list(placeholder = "Select one or more sender cells...")),
      selectizeInput("path_included", "Pathways Included",
                     choices = pathway_choices_with_all, multiple = TRUE,
                     options = list(server = TRUE, placeholder = "Select pathways or All...")),
      selectInput("de_method", "Differential Expression Method",
                  choices = c("DEsingle (cell-level)" = "DEsingle", "edgeR (Pseudo-bulk)" = "edgeR", "MAST (cell-level)" = "MAST"),
                  selected = "DEsingle"
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
            tabPanel("Ligand Only",       imageOutput("ligandOnlyPlot",   height = "750px"))
            
          )
          
        ),
        tabPanel(
          "Tables",
          div(style = "margin: 12px 0;",
              selectInput(
                "table_choice", "Result Table",
                choices = c(
                  "Pathway summary" = "pathway_n",
                  "Ligand summary"  = "lig_rank_all",
                  "Pathway P-values"= "pval_pathway",
                  "Ligand P-values" = "pval_ligand"
                ),
                selected = "pathway_n"
              )
          ),
          DTOutput("result_table")
        )
        
      )
    )
    
    
    
    
    
  )
)

server <- function(input, output, session) {
  abbrev_lookup <- reactiveVal(NULL)
  label_for <- function(full) {
    tbl <- abbrev_lookup()
    if (is.null(tbl)) return(full)
    i <- match(full, tbl$full)
    ifelse(is.na(i), full, tbl$abbr[i])
  }
  
  
  
  meta_data <- reactiveVal(NULL)
  analysis_result <- reactiveVal()
  last_params <- reactiveValues(
    receiver_cell = NULL, 
    diff_target_cell = NULL, 
    normalization_label = NULL, 
    output_dir = NULL
  )
  
  # where we will serve images from (www is auto-served by Shiny)
  plots_www_dir <- normalizePath(file.path("www", "plots"), mustWork = FALSE)
  if (!dir.exists(plots_www_dir)) dir.create(plots_www_dir, recursive = TRUE)
  
  # keep final, absolute PNG paths here (so renderImage can use them)
  plot_paths <- reactiveValues(path = NULL, ligand = NULL, ligand_only = NULL)
  
  # convenience: build the 3 treemap names
  png_names <- function(receiver, target, norm) {
    c(
      sprintf("Tree_path_%s_to_%s_%s.png",          receiver, target, norm),
      sprintf("Tree_ligand_%s_to_%s_%s.png",        receiver, target, norm),
      sprintf("Tree_ligand_only_%s_to_%s_%s.png",   receiver, target, norm)
    )
  }
  
  
  
  
  observeEvent(input$cdsfile, {
    req(input$cdsfile)
    # Use your actual data load:
    meta <- as.data.frame(pData(readRDS(input$cdsfile$datapath)))
    meta_data(meta)
    
    # Dynamically get columns with multiple unique values
    col_sizes <- sapply(meta, function(x) length(unique(x)))
    keep_cols <- names(col_sizes)[col_sizes > 1 & col_sizes < 200]
    updateSelectizeInput(session, "clustering_col", choices = keep_cols, selected = keep_cols[1])
    
    # Reset legend and cell pickers
    abbrev_lookup(NULL)
    updateSelectizeInput(session, "receiver_cell", choices = NULL)
    updateSelectizeInput(session, "diff_target_cell", choices = NULL)
    updateSelectizeInput(session, "sender_cells", choices = NULL)
  })
  
  
  cds_data <- reactive({
    req(input$cdsfile)
    filetype <- tools::file_ext(input$cdsfile$name)
    if (filetype == "rds") {
      readRDS(input$cdsfile$datapath)
    } else if (filetype == "h5ad") {
      convert_h5ad_to_cds(input$cdsfile$datapath)
    } else {
      stop("Unsupported file type.")
    }
  })
  
  output$filetype_status <- renderUI({
    req(input$cdsfile)
    filetype <- tools::file_ext(input$cdsfile$name)
    msg <- switch(filetype,
                  "rds" = "Detected Monocle CDS (.rds)",
                  "h5ad" = "Detected Scanpy (.h5ad) — converting to CDS.",
                  "Unknown or unsupported file type.")
    HTML(paste0("<b>", msg, "</b>"))
  })
  
  
  
  
  
  # Helper to get cluster choices safely
  get_clusters <- function() {
    if (is.null(input$clustering_col) || input$clustering_col == "") return(NULL)
    meta <- as.data.frame(pData(cds_data()))
    clusters <- unique(as.character(meta[[input$clustering_col]]))
    clusters
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
    table <- abbrev_lookup()
    if (is.null(table)) return(NULL)
    # Only show legend if abbr actually differs from full (case insensitive)
    show_legend <- any(tolower(table$abbr) != tolower(table$full))
    if (!show_legend) return(NULL)
    HTML(paste0("<b>Legend:</b><br>",
                paste0("<span style='font-family:monospace'>",
                       table$abbr, " = ", table$full,
                       collapse = "<br>")))
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
                 label = receiver, size = 7 * shrink_factor, fontface = "bold")
    }
    
    if (!is.null(diff_target) && diff_target != "") {
      gg <- gg +
        ggforce::geom_circle(aes(x0 = center_x + offset_x, y0 = center_y, r = r_main),
                             fill = "#cc6b37", color = "black", linewidth = 1) +
        annotate("text", x = center_x + offset_x, y = center_y,
                 label = diff_target, size = 7 * shrink_factor, fontface = "bold") +
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
      circle_df <- data.frame(sx = sx, sy = sy, label = senders)
      
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
        geom_text(data = circle_df, aes(x = sx, y = sy, label = label),
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
  
  
 
  
  # render once
  output$result_table <- DT::renderDT(
    DT::datatable(
      data.frame(),  # placeholder
      options = list(pageLength = 10, scrollX = TRUE, deferRender = TRUE),
      rownames = FALSE
    ),
    server = TRUE
  )
  
  # 1) keep this so the table renders even when the tab is hidden
  outputOptions(output, "result_table", suspendWhenHidden = FALSE)
  
  # 2) helper
  get_tbl <- function(name) {
    ar <- analysis_result(); if (is.null(ar)) return(NULL)
    out <- switch(name,
                  "pathway_n"    = ar$pathway_n,
                  "lig_rank_all" = ar$lig_rank_all,
                  "pval_pathway" = ar$pval_pathway,
                  "pval_ligand"  = ar$pval_ligand,
                  NULL
    )
    if (is.null(out)) return(NULL)
    as.data.frame(out)
  }
  
  # 3) initial render (empty)
  output$result_table <- DT::renderDT(
    DT::datatable(data.frame(), options = list(pageLength = 10, scrollX = TRUE, deferRender = TRUE), rownames = FALSE),
    server = TRUE
  )
  proxy <- DT::dataTableProxy("result_table")
  
  # 4) smart updater: re-render when schema changes, otherwise replaceData
  last_cols <- reactiveVal(NULL)
  
  observeEvent(list(input$table_choice, analysis_result()), {
    req(analysis_result())
    name <- if (isTruthy(input$table_choice)) input$table_choice else "pathway_n"
    tbl  <- get_tbl(name); req(tbl)
    
    cols <- names(tbl)
    if (is.null(last_cols()) || !identical(cols, last_cols())) {
      # schema changed → full render
      output$result_table <- DT::renderDT(
        DT::datatable(tbl, options = list(pageLength = 10, scrollX = TRUE, deferRender = TRUE), rownames = FALSE),
        server = TRUE
      )
      last_cols(cols)
    } else {
      # same schema → cheap swap
      DT::replaceData(proxy, tbl, resetPaging = FALSE, rownames = FALSE)
    }
  }, ignoreInit = FALSE)
  
  
  proxy <- DT::dataTableProxy("result_table")
  
  get_tbl <- function(name) {
    ar <- analysis_result(); if (is.null(ar)) return(NULL)
    out <- switch(name,
                  "pathway_n"    = ar$pathway_n,
                  "lig_rank_all" = ar$lig_rank_all,
                  "pval_pathway" = ar$pval_pathway,
                  "pval_ligand"  = ar$pval_ligand,
                  NULL
    )
    if (is.null(out)) return(NULL)
    as.data.frame(out)
  }
  
  observeEvent(list(input$table_choice, analysis_result()), {
    req(analysis_result())
    name <- if (isTruthy(input$table_choice)) input$table_choice else "pathway_n"
    tbl  <- get_tbl(name); req(tbl)
    DT::replaceData(proxy, tbl, resetPaging = FALSE, rownames = FALSE)
  }, ignoreInit = FALSE)
  
  
  
  
  # Helper: Generate abbreviations only if needed
  abbreviate_clusters <- function(clusters) {
    print("abbreviate_clusters called with:")
    print(clusters)
    clusters <- as.character(clusters)
    abbr <- ifelse(nchar(clusters) <= 3, clusters, toupper(substring(gsub("[^A-Za-z]", "", clusters), 1, 1:2, 1:2)))
    abbr <- ifelse(nchar(clusters) <= 3, clusters, make.unique(abbr, sep = "_"))
    data.frame(full = clusters, abbr = abbr, stringsAsFactors = FALSE)
  }
  
  output$treemapPlot <- renderUI({
    img_src <- file.path("plots", paste0("Tree_path_", input$receiver_cell, "_to_", input$diff_target_cell, "_", input$normalization_label, ".png"))
    tags$img(src = img_src, width = "600px", height = "750px")
  })
  
  
  # When the clustering column is chosen, update abbreviations table ONCE
  observeEvent(input$clustering_col, {
    req(meta_data(), input$clustering_col)
    clusters <- unique(as.character(meta_data()[[input$clustering_col]]))
    # Only abbreviate if any are long
    need_abbr <- any(nchar(clusters) > 3)
    if (need_abbr) {
      abbr <- ifelse(nchar(clusters) <= 3, clusters, toupper(substring(gsub("[^A-Za-z]", "", clusters), 1, 2)))
      abbr <- ifelse(nchar(clusters) <= 3, clusters, make.unique(abbr, sep = "_"))
      table <- data.frame(full = clusters, abbr = abbr, stringsAsFactors = FALSE)
      abbrev_lookup(table)
    } else {
      abbrev_lookup(NULL)
    }
    updateSelectizeInput(session, "receiver_cell", choices = clusters, selected = NULL)
    updateSelectizeInput(session, "diff_target_cell", choices = clusters, selected = NULL)
    updateSelectizeInput(session, "sender_cells", choices = clusters, selected = NULL)
  })
  
  
  # Example button
  observeEvent(input$show_example, {
    demo_table <- readRDS("www/example_result_table.rds")
    output$treemap <- renderPlot({
      img <- png::readPNG("www/example_treemap.png")
      grid::grid.raster(img)
    })
    output$result_table <- renderDT({
      as.data.frame(demo_table)
    })
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
      path_included    <- if ("All" %in% input$path_included) pathway_choices else input$path_included
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
          incProgress(0.1, detail = "Differential Expression")
          temp_cds <- cds_data()[, colData(cds_data())[[input$clustering_col]] %in% pair]
          expr_mat <- as.matrix(counts(temp_cds))
          meta <- as.data.frame(colData(temp_cds))
          meta$cluster <- as.character(meta[[input$clustering_col]])
          meta <- meta[meta$cluster %in% pair, ]
          expr_mat <- expr_mat[, rownames(meta)] # Make sure columns match meta
          
          n_rep <- 3  # Set number of pseudo-replicates
          set.seed(42) # For reproducibility
          
          meta$pseudorep <- NA
          
          for (cl in unique(meta$cluster)) {
            idx <- which(meta$cluster == cl)
            # Shuffle and split into n_rep groups
            splits <- split(idx, cut(seq_along(idx), breaks = n_rep, labels = FALSE))
            for (i in seq_along(splits)) {
              meta$pseudorep[splits[[i]]] <- paste0(cl, "_rep", i)
            }
          }
          
          # Aggregate counts by pseudo-replicate
          pseudo_bulk <- rowsum(t(expr_mat), meta$pseudorep)
          pseudo_bulk <- t(pseudo_bulk) # genes x pseudo-replicates
          
          # Assign cluster label to each pseudo-replicate
          pseudo_groups <- factor(sapply(strsplit(colnames(pseudo_bulk), "_"), `[`, 1))
          
          dge <- DGEList(counts = pseudo_bulk, group = pseudo_groups)
          dge <- calcNormFactors(dge)
          design <- model.matrix(~pseudo_groups)
          dge <- estimateDisp(dge, design)
          fit <- glmQLFit(dge, design)
          qlf <- glmQLFTest(fit, coef=2)
          res <- topTags(qlf, n=Inf)$table
          
          # Calculate means for each cluster (averaging pseudo-reps)
          cl1 <- levels(pseudo_groups)[1]
          cl2 <- levels(pseudo_groups)[2]
          cl1_cols <- grep(paste0("^", cl1, "_rep"), colnames(pseudo_bulk))
          cl2_cols <- grep(paste0("^", cl2, "_rep"), colnames(pseudo_bulk))
          Cluster1_total_mean <- rowMeans(pseudo_bulk[, cl1_cols, drop=FALSE])
          Cluster2_total_mean <- rowMeans(pseudo_bulk[, cl2_cols, drop=FALSE])
          
          gene_symbols <- rowData(temp_cds)$gene_short_name
          names(gene_symbols) <- rownames(temp_cds)
          
          # Now, for each gene (rownames(res)), get the symbol
          deg <- data.frame(
            Gene = gene_symbols[rownames(res)],
            Ensembl = rownames(res),  # Optional: keep original ID as well
            Cluster1_total_mean = Cluster1_total_mean[rownames(res)],
            Cluster2_total_mean = Cluster2_total_mean[rownames(res)],
            foldChange = res$logFC,
            pvalue = res$PValue,
            pvalue.adj.FDR = res$FDR,
            Type = ifelse(res$logFC > 0, "up", "down"),
            State = ifelse(res$logFC > 0, "up", "down")
          )
          
          saveRDS(deg, deg_file)
        } else if (input$de_method == "MAST") {
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
      
      ## === Enrichment checks (receiver & target) ===
      enrich_1_file <- paste0(meta.file, "/", receiver_cell, "_enriched.rds")
      enrich_2_file <- paste0(meta.file, "/", diff_target_cell, "_enriched.rds")
      if (file.exists(enrich_1_file)) {
        RecClus <- readRDS(enrich_1_file)
        message("Loaded existing enrichment for: ", receiver_cell)
      } else {
        temp_cds_r <- cds_data()[, colData(cds_data())[[input$clustering_col]] %in% in_groups]
        colData(temp_cds_r)[[input$clustering_col]] <- droplevels(colData(temp_cds_r)[[input$clustering_col]])
        GCMat <- counts(temp_cds_r)
        rownames(GCMat) <- rowData(temp_cds_r)$gene_short_name
        colnames(GCMat) <- as.character(colnames(temp_cds_r))
        clustering <- data.frame(Barcode = colnames(temp_cds_r), Cluster = colData(temp_cds_r)[[input$clustering_col]])
        write.table(clustering, "database/barcodetype.txt", sep = "\t")
        BarCluFile <- "database/barcodetype.txt"
        BarCluTable <- read.table(BarCluFile, sep = "\t", header = TRUE, stringsAsFactors = FALSE)
        RecClu <- receiver_cell
        LigClu <- setdiff(in_groups, c(RecClu, diff_target_cell))
        pval <- 0.15
        logfc <- 0.15
        cores_num <- parallel::detectCores()
        cores <- ifelse(cores_num > 1, cores_num - 1, 1)
        RecClus <- getHighExpGene(GCMat, BarCluTable, RecClu, LigClu, pval, logfc, cores)
        saveRDS(RecClus, enrich_1_file)
      }
      if (file.exists(enrich_2_file)) {
        RecClus_t <- readRDS(enrich_2_file)
        message("Loaded existing enrichment for: ", diff_target_cell)
      } else {
        temp_cds_r <- cds_data()[, colData(cds_data())[[input$clustering_col]] %in% in_groups]
        colData(temp_cds_r)[[input$clustering_col]] <- droplevels(colData(temp_cds_r)[[input$clustering_col]])
        GCMat <- counts(temp_cds_r)
        rownames(GCMat) <- rowData(temp_cds_r)$gene_short_name
        colnames(GCMat) <- as.character(colnames(temp_cds_r))
        clustering <- data.frame(Barcode = colnames(temp_cds_r), Cluster = colData(temp_cds_r)[[input$clustering_col]])
        write.table(clustering, "database/barcodetype.txt", sep = "\t")
        BarCluFile <- "database/barcodetype.txt"
        BarCluTable <- read.table(BarCluFile, sep = "\t", header = TRUE, stringsAsFactors = FALSE)
        RecClu_t <- diff_target_cell
        LigClu_t <- setdiff(in_groups, RecClu_t)
        pval <- 0.15
        logfc <- 0.15
        cores_num <- parallel::detectCores()
        cores <- ifelse(cores_num > 1, cores_num - 1, 1)
        RecClus_t <- getHighExpGene(GCMat, BarCluTable, RecClu_t, LigClu_t, pval, logfc, cores)
        saveRDS(RecClus_t, enrich_2_file)
      }
      
      ## === talklr Rdata check ===
      talklr_file <- paste0(meta.file, "/talklr_", receiver_cell, ".Rdata")
      if (file.exists(talklr_file)) {
        load(talklr_file)
        message("Loaded existing talklr Rdata: ", talklr_file)
      } else {
        temp_cds2 <- cds_data()[, colData(cds_data())[[input$clustering_col]] %in% c(receiver_cell, diff_target_cell, in_groups)]
        colData(temp_cds2)[[input$clustering_col]] <- droplevels(colData(temp_cds2)[[input$clustering_col]])
        cell_group_df <- data.frame(cell = colnames(temp_cds2), group = colData(temp_cds2)[[input$clustering_col]])
        glom_normal <- aggregate_gene_expression(temp_cds2, cell_group_df = cell_group_df, norm_method = "size_only")
        rownames(glom_normal) <- rowData(temp_cds2)$gene_short_name
        glom_normal <- glom_normal / as.vector(table(colData(temp_cds2)[[input$clustering_col]]))
        glom_normal <- as.data.frame(glom_normal)
        glom_normal$genes <- rownames(glom_normal)
        glom_normal <- glom_normal[, c(ncol(glom_normal), 1:(ncol(glom_normal)-1))]
        print("glom_normal:")
        print(glom_normal)
        min_exprs <- 1.713682e-11
        thresh_exprs <- .909563e-08
        receptor_ligand_sub <- receptor_ligand[receptor_ligand$pathway_name %in% path_included, ]
        print("path_included:")
        print(path_included)
        print("receptor_ligand_sub:")
        print(receptor_ligand_sub)
        # Make sure glom_normal$genes is character and upper-case
        glom_normal$genes <- toupper(as.character(glom_normal$genes))
        
        cat("First 10 glom_normal$genes:\n")
        print(head(glom_normal$genes, 10))
        
        cat("First 10 Ligand.ApprovedSymbol:\n")
        print(head(receptor_ligand_sub$Ligand.ApprovedSymbol, 10))
        
        cat("First 10 Receptor.ApprovedSymbol:\n")
        print(head(receptor_ligand_sub$Receptor.ApprovedSymbol, 10))
        
        cat("Number of genes in both glom_normal$genes and ligand symbols:\n")
        print(length(intersect(glom_normal$genes, receptor_ligand_sub$Ligand.ApprovedSymbol)))
        
        cat("Number of genes in both glom_normal$genes and receptor symbols:\n")
        print(length(intersect(glom_normal$genes, receptor_ligand_sub$Receptor.ApprovedSymbol)))
        
        cat("Number of ligand-receptor pairs with both genes expressed:\n")
        n_expressed_pairs <- nrow(receptor_ligand_sub[
          receptor_ligand_sub$Ligand.ApprovedSymbol %in% glom_normal$genes &
            receptor_ligand_sub$Receptor.ApprovedSymbol %in% glom_normal$genes, ])
        print(n_expressed_pairs)
        cat("Column classes of glom_normal before passing to make_expressed_net:\n")
        print(sapply(glom_normal, class))
        
        # Convert all expression columns to numeric (except first column, 'genes')
        glom_normal[, 2:ncol(glom_normal)] <- lapply(glom_normal[, 2:ncol(glom_normal)], as.numeric)
        
        # Optional: check for NAs (indicating conversion failed on some values)
        cat("Any NA in glom_normal?\n")
        print(any(is.na(as.matrix(glom_normal[, 2:ncol(glom_normal)]))))
        
        
        lr_glom_normal <- make_expressed_net(glom_normal, expressed_thresh = as.numeric(thresh_exprs),
                                             receptor_ligand_sub, KL_method = 'product',
                                             pseudo_count = as.numeric(min_exprs))
        lr_glom_normal <- dplyr::arrange(lr_glom_normal, desc(KL))
        lr_glom_normal <- lr_glom_normal[lr_glom_normal$Pair.Evidence == "literature supported", ]
        lig_col <- 17:(17 + ncol(glom_normal) - 2)
        rec_col <- (17 + ncol(glom_normal) - 1):((17 + ncol(glom_normal) - 1) + ncol(glom_normal) - 2)
        lig_mat <- lr_glom_normal[, lig_col]
        colnames(lig_mat) <- colnames(glom_normal)[-1]
        recived_in_sig <- c()
        edgelist_all <- matrix(, nrow = 0, ncol = 2)
        edgelist_rec_only <- matrix(, nrow = 0, ncol = 2)
        for (i in seq_len(nrow(lr_glom_normal))) {
          net <- plot_lr_wiring(as.numeric(lr_glom_normal[i, lig_col]),
                                as.numeric(lr_glom_normal[i, rec_col]),
                                cell_labels = colnames(glom_normal)[-1],
                                thresh = thresh_exprs)
          if (!is.null(net) && inherits(net, "igraph")) {
            edgelist <- igraph::as_edgelist(net)
            edgelist_all <- rbind(edgelist_all, edgelist)
            if (receiver_cell %in% edgelist[,2]) {
              recived_in_sig <- c(recived_in_sig, i)
              senders <- edgelist[which(edgelist[,2] %in% receiver_cell), , drop = FALSE][,1]
              lig_mat[i, !colnames(lig_mat) %in% senders] <- 0
              edgelist_rec_only <- rbind(edgelist_rec_only, edgelist)
            }
          }
        }
        lr_glom_normal <- lr_glom_normal[recived_in_sig, ]
        lig_mat <- lig_mat[recived_in_sig, ]
        save(list = c("lr_glom_normal", "glom_normal", "lig_mat", "edgelist_rec_only", "deg"),
             file = talklr_file)
      }
      print("lr_glom_normal:")
      print(lr_glom_normal)
      temp_cds2 <- cds_data()[, colData(cds_data())[[input$clustering_col]] %in% c(receiver_cell, diff_target_cell, in_groups)]
      colData(temp_cds2)[[input$clustering_col]] <- droplevels(colData(temp_cds2)[[input$clustering_col]])
      GCMat <- counts(temp_cds2)
      rownames(GCMat) <- rowData(temp_cds2)$gene_short_name
      colnames(GCMat) <- as.character(colnames(temp_cds2))
      clustering <- data.frame(Barcode = colnames(temp_cds2), Cluster = colData(temp_cds2)[[input$clustering_col]])
      write.table(clustering, "database/barcodetype.txt", sep = "\t")
      BarCluFile <- "database/barcodetype.txt"
      BarCluTable <- read.table(BarCluFile, sep = "\t", header = TRUE, stringsAsFactors = FALSE)
      enrich_1 <- readRDS(paste0(meta.file,"/", receiver_cell, "_enriched.rds"))
      enrich_2 <- readRDS(paste0(meta.file,"/", diff_target_cell, "_enriched.rds"))
      deg <- readRDS(deg_file)
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
      workdir <- paste(user.file,"metadata", sep="/")
      PyHome <- "python3"
      DrawMLnet(netList, LigClu, RecClu, workdir, PyHome, plotMLnet = TRUE)
      saveRDS(netList, paste0(out.file, "/", workdir, "/netList.rds"))
      
      # Downstream reporting, treemap and table
      # --- Downstream reporting: write PNGs to a temp folder, then copy into www/plots
      normalization_label <- "db"
      
      # Let top_pathway write to a temp output dir
      tmp_outdir <- tempfile("treemaps_")
      dir.create(tmp_outdir, recursive = TRUE, showWarnings = FALSE)
      
      result <- top_pathway(
        scMlnet_results = netList,
        deg             = deg,
        receiver_cell   = receiver_cell,
        target_cell     = diff_target_cell,
        lr_glom_normal  = lr_glom_normal,
        method          = "KL_rec",
        output_dir      = tmp_outdir
      )
      
      incProgress(0.15, detail = "Preparing outputs")
      
      # Store tables (prefer result$db if present)
      res_norm <- if (!is.null(result$db)) result$db else result
      analysis_result(res_norm)
      message("Tables available: ", paste(names(analysis_result()), collapse=", "))
      message("pathway_n class: ", paste(class(analysis_result()$pathway_n), collapse=", "))
      message("pathway_n dim: ", tryCatch(paste(dim(analysis_result()$pathway_n), collapse=" x "), error=function(e) "no dim"))
      
      
      # Copy the three expected PNGs to www/plots and update reactive paths
      need <- png_names(receiver_cell, diff_target_cell, normalization_label)
      srcs <- file.path(tmp_outdir, need)
      dests <- file.path(plots_www_dir, need)
      
      # copy if they exist
      to_copy <- file.exists(srcs)
      if (any(to_copy)) file.copy(srcs[to_copy], dests[to_copy], overwrite = TRUE)
      
      # point image outputs at the absolute files under www/plots
      plot_paths$path        <- if (file.exists(dests[1])) dests[1] else NULL
      plot_paths$ligand      <- if (file.exists(dests[2])) dests[2] else NULL
      plot_paths$ligand_only <- if (file.exists(dests[3])) dests[3] else NULL
      
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
  
  observeEvent(analysis_result(), {
    req(analysis_result())
    # render once with the current selected table (or pathway_n)
    name <- if (isTruthy(input$table_choice)) input$table_choice else "pathway_n"
    tbl  <- get_tbl(name); req(tbl)
    
    output$result_table <- DT::renderDT(
      DT::datatable(
        tbl,
        options = list(pageLength = 10, scrollX = TRUE, deferRender = TRUE),
        rownames = FALSE
      ),
      server = TRUE
    )
  }, ignoreInit = FALSE)
  
  
  observeEvent(input$treemap_tabs, {
    req(analysis_result())  # make sure we actually have tables
    selected <- switch(input$treemap_tabs,
                       "Pathway"            = "pathway_n",
                       "Ligand by Pathway"  = "lig_rank_all",
                       "Ligand Only"        = "lig_rank_all",
                       "pathway_n"  # fallback
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
  
  
  
}


shinyApp(ui = ui, server = server)

