top_pathway <- function(
    scMlnet_results,
    deg,
    receiver_cell,
    target_cell,
    TF_fun       = "sum",
    main_fun     = "sum",
    path_fun     = "sum",
    lr_glom_normal,
    show.sub.pathway = FALSE,
    method       = "KL_rec",
    proximity_order = NULL,
    # -- Named color palette. Adjust as needed for your own pathways.
    #   The names must match EXACTLY how your data labels each pathway.
    #   If you have more or fewer pathways, add/remove them here.
    palette = c(
      NOTCH   = "#9AA556",
      WNT     = "#E98FC5",
      ncWNT   = "#64B1DF",
      EGF     = "#66E1E1",
      IGF     = "#D5AB90",
      NRG     = "#8B9D63",
      BMP     = "#E88E50",
      ACTIVIN = "#C2E579",
      TGFb    = "#DBC6DE",
      NT      = "#DE8D9F",
      VEGF    = "#AAD8E1",
      HH      = "#6EE07F",
      FGF     = "#8580DF",
      HGF     = "#D9E83D",
      PDGF    = "#62E4B7",
      GDNF    = "#E2E2A1"
    ),
    mirror.y = FALSE,
    output_dir = "."
) {
  suppressWarnings({
    # ---------------------------------------------------------
    # 1)  Read any reference data (if needed in your environment)
    # ---------------------------------------------------------
    tf_ori_data  <- read.delim("./database/TFTargetGene.txt", header = TRUE)
    rec_ori_data <- read.delim("./database/RecTF.txt",        header = TRUE)
    lig_ori_data <- read.delim("./database/LigRec.txt",       header = TRUE)
    
    tf_size_factor  <- table(tf_ori_data$TF)
    rec_size_factor <- table(rec_ori_data$Receptor)
    lig_size_factor <- table(lig_ori_data$Ligand)
    
    # Example: If you rely on a CellChat-like DB:
    # Make sure "CellChatDB" is loaded in your environment
    # It should have columns "ligand", "pathway_name", etc.
    
    # ---------------------------------------------------------
    # 2)  Prepare expression data
    # ---------------------------------------------------------
    expres <- deg$foldChange
    names(expres) <- deg$Gene
    
    # Replace ±Inf with max/min finite
    expres[expres ==  Inf] <- max(expres[is.finite(expres)], na.rm = TRUE)
    expres[expres == -Inf] <- min(expres[is.finite(expres)] & expres != 0, na.rm = TRUE)
    # Convert values <1 by reciprocal
    expres[expres < 1] <- 1 / expres[expres < 1]
    expres[expres ==  Inf] <- max(expres[is.finite(expres)], na.rm = TRUE)
    expres[expres == -Inf] <- min(expres[is.finite(expres)] & expres != 0, na.rm = TRUE)
    
    # Basic receptor expression from lr_glom_normal
    rec_exp <- setNames(
      lr_glom_normal[!duplicated(lr_glom_normal$Receptor.ApprovedSymbol),
                     paste0("receptor_", receiver_cell)],
      unique(lr_glom_normal$Receptor.ApprovedSymbol)
    )
    
    KL_DF <- lr_glom_normal[, c("Ligand.ApprovedSymbol","Receptor.ApprovedSymbol","KL")]
    
    # Split the columns in scMlnet_results:
    #   scMlnet_results$LigRec, scMlnet_results$RecTF, scMlnet_results$TFTar
    scMlnet_results$LigRec <- stringr::str_split(
      scMlnet_results$LigRec, "_", simplify = TRUE
    )
    scMlnet_results$RecTF <- stringr::str_split(
      scMlnet_results$RecTF, "_", simplify = TRUE
    )
    scMlnet_results$TFTar <- stringr::str_split(
      scMlnet_results$TFTar, "_", simplify = TRUE
    )
    
    # ---------------------------------------------------------
    # 3)  Compute TF -> target aggregates
    # ---------------------------------------------------------
    scMlnet_results$TFTar <- as.data.frame(scMlnet_results$TFTar)
    scMlnet_results$TFTar$V3 <- expres[scMlnet_results$TFTar$V2]
    
    tf_target <- scMlnet_results$TFTar %>%
      dplyr::group_by(V1) %>%
      dplyr::summarise(aggregate = get(TF_fun)(V3), .groups = "drop")
    
    tf_n <- setNames(tf_target$aggregate, tf_target$V1)
    # Normalize by TF size factor
    tf_n <- tf_n / tf_size_factor[names(tf_n)]
    
    # ---------------------------------------------------------
    # 4)  Compute Receptor-level aggregates
    # ---------------------------------------------------------
    scMlnet_results$RecTF <- as.data.frame(scMlnet_results$RecTF)
    scMlnet_results$RecTF$tf_n <- 0
    
    # Safely assign tf_n for each row matching a TF name
    names_tf <- names(tf_n)
    for (tf_i in seq_along(tf_n)) {
      idx <- scMlnet_results$RecTF$V2 %in% names_tf[tf_i]
      if (any(idx)) {
        scMlnet_results$RecTF[idx, "tf_n"] <- tf_n[tf_i]
      }
    }
    
    # Summation (or mean, etc.)
    rec_n <- scMlnet_results$RecTF %>%
      dplyr::group_by(V1) %>%
      dplyr::summarise(aggregate = get(main_fun)(tf_n), .groups = "drop")
    rec_n <- setNames(rec_n$aggregate, rec_n$V1)
    
    # ---------------------------------------------------------
    # 5)  Method-specific ways to incorporate receptor expression
    # ---------------------------------------------------------
    if (method == "KL_rec") {
      rec_exp <- setNames(
        lr_glom_normal[!duplicated(lr_glom_normal$Receptor.ApprovedSymbol), "KL"],
        unique(lr_glom_normal$Receptor.ApprovedSymbol)
      )
      # Rescale rec_exp to same range as rec_n
      rec_exp <- scales::rescale(
        rec_exp,
        from = range(rec_exp, na.rm = TRUE),
        to   = range(rec_n,   na.rm = TRUE)
      )
      rec_n <- rec_n + rec_exp[names(rec_n)]
    } else if (method == "min_max") {
      rec_exp <- scales::rescale(
        rec_exp,
        from = range(rec_exp, na.rm = TRUE),
        to   = range(rec_n,   na.rm = TRUE)
      )
      rec_n <- rec_n + rec_exp[names(rec_n)]
    } else if (method == "scale_no_center") {
      rec_n <- scale(rec_n, center = FALSE, scale = TRUE) +
        scale(rec_exp[names(rec_n)], center = FALSE, scale = TRUE)
      rec_n <- setNames(rec_n[,1], rownames(rec_n))
    } else if (method == "scale_sd") {
      rec_n <- scale(
        rec_n, center = FALSE, 
        scale  = sd(rec_n, na.rm = TRUE)
      ) + scale(
        rec_exp[names(rec_n)], center = FALSE, 
        scale = sd(rec_exp[names(rec_n)], na.rm = TRUE)
      )
      rec_n <- setNames(rec_n[,1], rownames(rec_n))
    } else if (method == "rescale") {
      rec_n <- scales::rescale(as.vector(rec_n)) +
        scales::rescale(rec_exp[names(rec_n)])
    } else if (method == "KL_rec_lig") {
      # etc. (Handle KL merging for both ligand & receptor).
      # ...
    } else if (method == "NULL") {
      # If you'd like to do nothing:
      # message("No additional method scaling applied.")
    }
    
    # Adjust for receptor size factor
    rec_n <- rec_n / rec_size_factor[names(rec_n)]
    
    # ---------------------------------------------------------
    # 6)  Combine receptor aggregates -> ligand
    # ---------------------------------------------------------
    scMlnet_results$LigRec <- as.data.frame(scMlnet_results$LigRec)
    scMlnet_results$LigRec$rec_n <- 0
    
    names_rec <- names(rec_n)
    for (r in seq_along(rec_n)) {
      idx <- scMlnet_results$LigRec$V2 %in% names_rec[r]
      if (any(idx)) {
        scMlnet_results$LigRec[idx, "rec_n"] <- rec_n[r]
      }
    }
    
    lig_n <- scMlnet_results$LigRec %>%
      dplyr::group_by(V1) %>%
      dplyr::summarise(score = get(main_fun)(rec_n), .groups = "drop")
    lig_n <- setNames(lig_n$score, lig_n$V1)
    
    # ---------------------------------------------------------
    # 7)  Summarize by pathway & prepare for treemap
    # ---------------------------------------------------------
    # Create rank table
    lig_rank_all <- lig_n[order(lig_n, decreasing = TRUE)]
    lig_rank_all <- as.data.frame(lig_rank_all)
    colnames(lig_rank_all) <- "score"
    lig_rank_all$pathway <- ""
    
    # Suppose you have a data frame named CellChatDB that
    # has columns 'ligand' and 'pathway_name'.
    #
    # We match each ligand name to its pathway:
    for (p in rownames(lig_rank_all)) {
      match_idx <- CellChatDB$ligand == p
      if (any(match_idx)) {
        lig_rank_all[p, "pathway"] <- CellChatDB$pathway_name[match_idx][1]
      }
    }
    
    pathway_n <- lig_rank_all %>%
      dplyr::group_by(pathway) %>%
      dplyr::summarise(score = get(path_fun)(score), .groups = "drop")
    
    # Order by descending
    pathway_n <- pathway_n[order(pathway_n$score, decreasing = TRUE), ]
    # Convert raw scores into a percentage
    total_score <- sum(pathway_n$score)
    pathway_n$score <- (pathway_n$score / total_score) * 100
    pathway_n$score.perc <- scales::percent(pathway_n$score / 100, accuracy = 0.1)
    
    # ---------------------------------------------------------
    # 8)  Factor the pathway column to fix color assignment
    # ---------------------------------------------------------
    # Gather all *possible* pathway names from your palette:
    all_pathways <- names(palette)
    
    # Force the factor level so treemap won't reorder or skip missing ones
    pathway_n$pathway      <- factor(pathway_n$pathway,      levels = all_pathways)
    lig_rank_all$pathway   <- factor(lig_rank_all$pathway,   levels = all_pathways)
    
    # ---------------------------------------------------------
    # 9)  Create the treemap(s)
    # ---------------------------------------------------------
    dir.create(output_dir, showWarnings = FALSE, recursive = TRUE)
    
    # Main plot
    library(treemap)
    pdf(file.path(output_dir, paste0("Tree_path_", receiver_cell,"_to_",target_cell, ".pdf")), width = 6, height = 7.5)
    if (!show.sub.pathway) {
      treemap::treemap(
        pathway_n,
        index   = c("pathway", "score.perc"),
        vSize   = "score",
        type    = "categorical",
        vColor  = "pathway",
        palette = palette,  # use our named palette
        title   = paste(receiver_cell,"-->",target_cell),
        
        border.col  = c("black", "black"),
        border.lwds = c(3, 1),
        mirror.y    = mirror.y,
        
        fontsize.labels  = c(8, 8),
        fontcolor.labels = c("black", "black"),
        fontface.labels  = c(2, 1),
        bg.labels        = c("transparent"),
        align.labels     = list(c("left", "top"), c("right", "bottom")),
        overlap.labels   = 0,
        inflate.labels   = TRUE,
        position.legend  = "none"
      )
      dev.off()
      show.sub.pathway <- TRUE
      saveRDS(pathway_n,file.path(output_dir, paste0("Tree_path_", receiver_cell,"_to_",target_cell, ".rds")))
    }
    
    # Sub-pathway plot
    if (show.sub.pathway) {
      pdf(file.path(output_dir, paste0("Tree_path_", receiver_cell,"_to_",target_cell, "_sub.pdf")), width = 6, height = 7.5)
      # Combine pathway_n columns with lig_rank_all
      # so each "lig" has a sub region in the treemap
      lig_rank_all$lig <- rownames(lig_rank_all)
      merged_df <- merge(lig_rank_all, pathway_n, by.x = "pathway", by.y = "pathway")
      
      saveRDS(lig_rank_all,file.path(output_dir, paste0("lig_rank_all_", receiver_cell,"_to_",target_cell, ".rds")))
      saveRDS(merged_df,file.path(output_dir, paste0("merged_df_", receiver_cell,"_to_",target_cell, ".rds")))
      treemap::treemap(
        merged_df,
        index   = c("pathway", "score.perc", "lig"),
        vSize   = "score.x",       # 'score.x' is from lig_rank_all
        type    = "categorical",
        vColor  = "pathway",
        palette = palette,
        title   = paste(receiver_cell,"-->",target_cell, "_ sub pathways"),
        
        border.col  = c("black", "black", "white"),
        border.lwds = c(3, 1, 3),
        mirror.y    = mirror.y,
        
        fontsize.labels  = c(8, 0, 5),
        fontcolor.labels = c("black", "black", "white"),
        fontface.labels  = 2,
        bg.labels        = c("transparent"),
        align.labels     = list(
          c("left", "top"),
          c("right", "top"),
          c("center", "bottom")
        ),
        overlap.labels   = 0,
        inflate.labels   = TRUE,
        position.legend  = "none"
      )
      dev.off()
    }
    
    ### test for pvalue
    source("Pathways_pvalue_function.R")
    result_with_pval <- empirical_pathway_pvals(pathway_n, lig_rank_all, n_perm=1000)
    print(result_with_pval)
    
    
  })
}