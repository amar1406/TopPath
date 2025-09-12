options(shiny.maxRequestSize = 500*1024^2)
setwd("C:/Users/afrit/Documents/R_work")

library(future)
plan(multisession) # use background processes

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
library(promises)
library(magick)


# Pathways file
pathway_choices <- readLines("top_path/pathways.txt")
pathway_choices_with_all <- c("All", pathway_choices)

out.file <- "top_path_output"
user.file <- "shiny_path"
meta.file <- paste0(out.file, "/", user.file, "/metadata")
CellChatDB <- readRDS("database/CellChatDB.rds")

# Source your function files
source("top_path/top_path_function.R")
source("top_path/Run_scMLnet.R")
source("top_path/Draw_MLnet.R")
source("top_path/plot_lr_wiring.R")

ui <- fluidPage(
  titlePanel("Top Pathway Analysis:"),
  sidebarLayout(
    sidebarPanel(
  fileInput("cdsfile", "Upload Monocle CDS (.rds)", accept = ".rds"),
  uiOutput("cluster_col_ui"),
  uiOutput("receiver_cell_ui"),
  uiOutput("diff_target_cell_ui"),
  uiOutput("sender_cell_ui"),
  selectizeInput("path_included", "Pathways Included",
                 choices = pathway_choices_with_all, multiple = TRUE,
                 options = list(server = TRUE, placeholder = "Select pathways or All...")),
  # BUTTONS GROUPED AND ALIGNED
  div(
    style = "text-align:center;",
    actionButton("run", "Run Analysis"),
    actionButton("show_example", "Example"),
    actionButton("cancel_analysis", "Cancel", style = "color:white;background-color:#b20000;")
  )
  
)
,
    mainPanel(
      plotOutput("selection_plot", height = "600px"),
      shinycssloaders::withSpinner(plotOutput("treemap")),
      DTOutput("result_table"),
      uiOutput("analysis_status")
    )
    
  )
)

server <- function(input, output, session) {
  cds_data <- reactive({
    req(input$cdsfile)
    readRDS(input$cdsfile$datapath)
  })
  
  # Cancel token for background
  cancel_env <- reactiveValues(future = NULL, cancel_now = FALSE)
  
  # --- UI rendering (unchanged) ---
  output$cluster_col_ui <- renderUI({
    req(cds_data())
    meta <- as.data.frame(pData(cds_data()))
    col_sizes <- sapply(meta, function(x) length(unique(x)))
    keep_cols <- names(col_sizes)[col_sizes > 1 & col_sizes < 200]
    if(length(keep_cols) == 0) {
      helpText("No suitable clustering columns found.")
    } else {
      selectizeInput("clustering_col", "Choose clustering column",
                     choices = c("", keep_cols), selected = "",
                     options = list(
                       placeholder = "Select clustering column...",
                       server = TRUE
                     ))
    }
  })
  output$receiver_cell_ui <- renderUI({
    if(is.null(input$clustering_col) || input$clustering_col == "") {
      selectizeInput("receiver_cell", "Receiver Cell",
                     choices = "", selected = "", multiple = FALSE,
                     options = list(placeholder = "Select receiver cell...", server = TRUE))
    } else {
      meta <- as.data.frame(pData(cds_data()))
      clusters <- unique(as.character(meta[[input$clustering_col]]))
      selectizeInput("receiver_cell", "Receiver Cell",
                     choices = c("", clusters), selected = "",
                     multiple = FALSE,
                     options = list(placeholder = "Select receiver cell...", server = TRUE))
    }
  })
  output$diff_target_cell_ui <- renderUI({
    if(is.null(input$clustering_col) || input$clustering_col == "") {
      selectizeInput("diff_target_cell", "Diff Target Cell",
                     choices = "", selected = "", multiple = FALSE,
                     options = list(placeholder = "Select diff target cell...", server = TRUE))
    } else {
      meta <- as.data.frame(pData(cds_data()))
      clusters <- unique(as.character(meta[[input$clustering_col]]))
      selectizeInput("diff_target_cell", "Diff Target Cell",
                     choices = c("", clusters), selected = "",
                     multiple = FALSE,
                     options = list(placeholder = "Select diff target cell...", server = TRUE))
    }
  })
  output$sender_cell_ui <- renderUI({
    if(is.null(input$clustering_col) || input$clustering_col == "") {
      selectizeInput("sender_cells", "Sender Cells (multiple allowed)",
                     choices = character(0), selected = NULL,
                     multiple = TRUE,
                     options = list(placeholder = "Select one or more sender cells...", server = TRUE))
    } else {
      meta <- as.data.frame(pData(cds_data()))
      clusters <- unique(as.character(meta[[input$clustering_col]]))
      selectizeInput("sender_cells", "Sender Cells (multiple allowed)",
                     choices = clusters, selected = NULL,
                     multiple = TRUE,
                     options = list(placeholder = "Select one or more sender cells...", server = TRUE))
    }
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
    shrink_factor <- if (senders_n > 3) 0.90 ^pmax(0, senders_n - 3) else 1
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
        geom_circle(aes(x0 = center_x, y0 = center_y, r = r_main), fill = "#97c9ea", color = "black", linewidth = 1) +
        annotate("text", x = center_x, y = center_y, label = receiver, size = 7 * shrink_factor, fontface = "bold")
    }
    if (!is.null(diff_target) && diff_target != "") {
      gg <- gg +
        geom_circle(aes(x0 = center_x + offset_x, y0 = center_y, r = r_main), fill = "#cc6b37", color = "black", linewidth = 1) +
        annotate("text", x = center_x + offset_x, y = center_y, label = diff_target, size = 7 * shrink_factor, fontface = "bold") +
        geom_segment(aes(x = center_x + r_main, y = center_y, xend = center_x + offset_x - r_main, yend = center_y),
                     arrow = arrow(type = "closed", length = unit(0.28, "inches")), linewidth = 1.4)
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
      arrow_thetas <- if (n == 1) center_angle else seq(center_angle - arc_spread/2, center_angle + arc_spread/2, length.out = n)
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
        geom_text(data = circle_df, aes(x = sx, y = sy, label = label), size = 5 * shrink_factor + 0.7, fontface = "bold") +
        geom_segment(
          aes(x = sx, y = sy + r_sender, xend = dot_x, yend = dot_y),
          data = cbind(circle_df, dot_x, dot_y),
          linewidth = 0.9
        ) +
        geom_point(aes(x = dot_x, y = dot_y), size = 7.2 * shrink_factor + 0.1, color = "black", shape = 16) +
        geom_path(data = arc_df, aes(x = x, y = y, group = group), color = "#455a64", linewidth = 5 * shrink_factor + 0.2, lineend = "round")
    }
    if (show_treemap) {
      pdf_file <- file.path(tempdir(), paste0("Tree_path_", receiver, ".pdf"))
      png_file <- tempfile(fileext = ".png")
      system(sprintf("convert -density 150 %s %s", shQuote(pdf_file), shQuote(png_file)))
      treemap_img <- png::readPNG(png_file)
      gg <- gg +
        annotation_raster(treemap_img, xmin = center_x - 3, xmax = center_x + offset_x + 3, ymin = center_y + 2, ymax = center_y + 7)
    }
    gg
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
  
  # Cancel button logic
  observeEvent(input$cancel_analysis, {
    cancel_env$cancel_now <- TRUE
    if (!is.null(cancel_env$future)) {
      cancel(cancel_env$future)
      cancel_env$future <- NULL
    }
    output$analysis_status <- renderUI({
      HTML("<span style='color: #b20000;'><b>Analysis canceled.</b></span>")
    })
  })
  
  # The main analysis, runs in background
  observeEvent(input$run, {
    req(
      input$clustering_col != "",
      input$receiver_cell != "", 
      input$diff_target_cell != "",
      length(input$sender_cells) > 0,
      length(input$path_included) > 0
    )
    output$analysis_status <- renderUI({HTML("<b>Analysis running in background ...</b>")})
    output$treemap <- renderPlot(NULL)
    output$result_table <- renderDT(NULL)
    cancel_env$cancel_now <- FALSE
    
    # Everything is passed as value, not reactive, to the future!
    cdsfile <- input$cdsfile$datapath
    clustering_col <- input$clustering_col
    receiver_cell <- input$receiver_cell
    diff_target_cell <- input$diff_target_cell
    in_groups <- input$sender_cells
    path_included <- if ("All" %in% input$path_included) pathway_choices else input$path_included
    CellChatDB2 <- CellChatDB
    meta.file2 <- meta.file
    user.file2 <- user.file
    out.file2 <- out.file
    
    cancel_env$future <- future({
      # Everything here runs in the background!
      progress <- function(detail=NULL) NULL
      library(monocle3)
      library(DEsingle)
      library(talklr)
      library(stats)
      library(dplyr)
      cds <- readRDS(cdsfile)
      pair <- c(receiver_cell, diff_target_cell)
      
      # Differential Expression step
      progress("DE step")
      temp_cds <- cds[, colData(cds)[[clustering_col]] %in% pair]
      group <- factor(colData(temp_cds)[[clustering_col]])
      results <- DEsingle(counts = counts(temp_cds), group = group, parallel = TRUE)
      results.classified <- DEtype(results = results, threshold = 0.05)
      results.sig <- results.classified[results.classified$pvalue < 0.05, ]
      deg <- results.sig[,c(12,13,11,20,21,23,24)]
      deg$Gene <- rowData(temp_cds[rownames(deg),])$gene_short_name
      deg <- deg[,c(8,1,2,3,4,5,6,7)]
      deg$norm_total_mean_1 <- deg$norm_total_mean_1 * 1000
      deg$norm_total_mean_2 <- deg$norm_total_mean_2 * 1000
      colnames(deg) <- c("Gene","Cluster1_total_mean", "Cluster2_total_mean", "foldChange" , "pvalue", "pvalue.adj.FDR","Type", "State")
      saveRDS(deg, paste0(meta.file2,"/deg_0.05_",paste(pair,collapse = '_'), ".rds"))
      # Gene enrichment (receiver)
      temp_cds_r <- cds[, colData(cds)[[clustering_col]] %in% in_groups]
      colData(temp_cds_r)[[clustering_col]] <- droplevels(colData(temp_cds_r)[[clustering_col]])
      GCMat <- counts(temp_cds_r)
      rownames(GCMat) <- rowData(temp_cds_r)$gene_short_name
      colnames(GCMat) <- as.character(colnames(temp_cds_r))
      clustering <- data.frame(Barcode = colnames(temp_cds_r), Cluster = colData(temp_cds_r)[[clustering_col]])
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
      saveRDS(RecClus, paste0(meta.file2,"/",RecClu,"_enriched.rds"))
      # Gene enrichment (target)
      RecClu_t <- diff_target_cell
      LigClu_t <- setdiff(in_groups, RecClu_t)
      RecClus_t <- getHighExpGene(GCMat, BarCluTable, RecClu_t, LigClu_t, pval, logfc, cores)
      saveRDS(RecClus_t, paste0(meta.file2,"/",RecClu_t,"_enriched.rds"))
      # talklr input
      temp_cds2 <- cds[, colData(cds)[[clustering_col]] %in% c(receiver_cell, diff_target_cell, in_groups)]
      colData(temp_cds2)[[clustering_col]] <- droplevels(colData(temp_cds2)[[clustering_col]])
      cell_group_df <- data.frame(cell = colnames(temp_cds2), group = colData(temp_cds2)[[clustering_col]])
      glom_normal <- aggregate_gene_expression(temp_cds2, cell_group_df = cell_group_df, norm_method = "size_only")
      rownames(glom_normal) <- rowData(temp_cds2)$gene_short_name
      glom_normal <- glom_normal / as.vector(table(colData(temp_cds2)[[clustering_col]]))
      glom_normal <- as.data.frame(glom_normal)
      glom_normal$genes <- rownames(glom_normal)
      glom_normal <- glom_normal[, c(ncol(glom_normal), 1:(ncol(glom_normal)-1))]
      min_exprs <- 1.713682e-11
      thresh_exprs <- .909563e-08
      receptor_ligand <- readRDS("database/receptor_ligand.rds")
      receptor_ligand_sub <- receptor_ligand[receptor_ligand$Ligand.ApprovedSymbol %in%
                                               CellChatDB2[CellChatDB2$pathway_name %in% path_included, ]$ligand, ]
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
      GCMat <- counts(temp_cds2)
      rownames(GCMat) <- rowData(temp_cds2)$gene_short_name
      colnames(GCMat) <- as.character(colnames(temp_cds2))
      clustering <- data.frame(Barcode = colnames(temp_cds2), Cluster = colData(temp_cds2)[[clustering_col]])
      write.table(clustering, "database/barcodetype.txt", sep = "\t")
      BarCluFile <- "database/barcodetype.txt"
      BarCluTable <- read.table(BarCluFile, sep = "\t", header = TRUE, stringsAsFactors = FALSE)
      enrich_1 <- readRDS(paste0(meta.file2,"/", receiver_cell, "_enriched.rds"))
      enrich_2 <- readRDS(paste0(meta.file2,"/", diff_target_cell, "_enriched.rds"))
      deg <- readRDS(paste0(meta.file2,"/deg_0.05_", paste(pair, collapse = '_'), ".rds"))
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
      save(list = c("lr_glom_normal", "glom_normal", "lig_mat", "edgelist_rec_only", "deg"),
           file = paste0(meta.file2,"/talklr_", receiver_cell, ".Rdata"))
      netList <- RunMLnet(GCMat, BarCluFile, RecClu, LigClu,
                          pval, logfc,
                          LigRecLib, TFTarLib, RecTFLib, lr_glom_normal=lr_glom_normal)
      workdir <- paste(user.file2,"metadata", sep="/")
      PyHome <- "python3"
      DrawMLnet(netList, LigClu, RecClu, workdir, PyHome, plotMLnet = TRUE)
      saveRDS(netList, paste0(out.file2, "/", workdir, "/netList.rds"))
      result <- top_pathway(
        scMlnet_results   = netList,
        deg              = deg,
        receiver_cell    = receiver_cell,
        lr_glom_normal   = lr_glom_normal,
        show.sub.pathway = FALSE,
        method           = "KL_rec",
        output_dir       = tempdir()
      )
      list(
        pdf = file.path(tempdir(), paste0("Tree_path_", receiver_cell, ".pdf")),
        result_table = result$pathway_n
      )
    })
    
    # Handle future result (or error)
    cancel_env$future %...>% (function(res) {
      if (is.null(res)) return(NULL)
      # Render outputs if not canceled!
      output$treemap <- renderPlot({
        if (file.exists(res$pdf)) {
          png_file <- tempfile(fileext = ".png")
          system(sprintf("convert -density 150 %s %s", shQuote(res$pdf), shQuote(png_file)))
          img <- png::readPNG(png_file)
          grid::grid.raster(img)
        }
      })
      output$result_table <- renderDT({
        as.data.frame(res$result_table)
      })
      output$analysis_status <- renderUI({
        HTML("<span style='color: #388e3c;'><b>Analysis finished.</b></span>")
      })
    }) %...!% (function(e){
      output$analysis_status <- renderUI({
        HTML(sprintf("<span style='color: #b20000;'><b>Analysis canceled or error:<br/>%s</b></span>", as.character(e)))
      })
      output$treemap <- renderPlot(NULL)
      output$result_table <- renderDT(NULL)
    })
  })
}

shinyApp(ui = ui, server = server)

