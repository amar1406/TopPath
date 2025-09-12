top_pathway <- function(
    scMlnet_results,
    deg,
    receiver_cell,
    target_cell,
    TF_fun       = "sum",
    main_fun     = "sum",
    path_fun     = "sum",
    lr_glom_normal,
    method       = NULL,
    proximity_order = NULL,
    normalization_mode = c("db", "detected"),
    pval_cutoff = 0.05,
    palette = c(
      NOTCH="#9AA556", WNT="#E98FC5", ncWNT="#64B1DF", EGF="#66E1E1",
      IGF="#D5AB90", NRG="#8B9D63", BMP="#E88E50", ACTIVIN="#C2E579",
      TGFb="#DBC6DE", NT="#DE8D9F", VEGF="#AAD8E1", HH="#6EE07F",
      FGF="#8580DF", HGF="#D9E83D", PDGF="#62E4B7", GDNF="#E2E2A1"
    ),
    receiver_label = NULL,
    target_label   = NULL,
    mirror.y = FALSE,
    output_dir = ".",
    # ---- NEW: SCENIC weight controls ----
    tf_weight_vec = NULL,       # named numeric vector: names=TFs, values=weights
    scenic_alpha  = 1.0,        # 0=no effect, 1=full effect, 0.5=half
    scenic_gamma  = 1.0,        # >1 amplifies differences (e.g. 1.5–2)
    scenic_clip   = c(0.5, 4),  # clip range after transforms
    scenic_center = TRUE        # re-center to mean 1
) {
  
  suppressWarnings({
    normalization_mode <- match.arg(normalization_mode)
    message("=== Normalization mode: ", normalization_mode, " ===")
    dir.create(output_dir, showWarnings = FALSE, recursive = TRUE)
    
    # right after arg parsing:
    if (is.null(receiver_label)) receiver_label <- receiver_cell
    if (is.null(target_label))   target_label   <- target_cell
    
    # ---- Load reference data (with mode-aware TF->Target) ----
    # New file includes Mode: activation / repression / unknown
    tf_withmode  <- read.delim("./database/TFTargetGene_withmode.txt", header = TRUE, stringsAsFactors = FALSE)
    rec_ori_data <- read.delim("./database/RecTF.txt",                  header = TRUE, stringsAsFactors = FALSE)
    lig_ori_data <- read.delim("./database/LigRec.txt",                 header = TRUE, stringsAsFactors = FALSE)
    
    # Make sure the TF/Target/Mode columns exist & normalized
    # Accept common variants of the target column name
    tcol <- c("Target","Targets","target","targets")
    tcol <- tcol[tcol %in% names(tf_withmode)]
    if (length(tcol) == 0L) stop("TFTargetGene_withmode.txt must have a Target(s) column.")
    names(tf_withmode)[match(tcol[1], names(tf_withmode))] <- "Target"
    if (!"TF"   %in% names(tf_withmode)) stop("TFTargetGene_withmode.txt must have a TF column.")
    if (!"Mode" %in% names(tf_withmode)) tf_withmode$Mode <- "unknown"
    
    # Normalize case & allowed values
    tf_withmode$TF     <- toupper(trimws(tf_withmode$TF))
    tf_withmode$Target <- toupper(trimws(tf_withmode$Target))
    tf_withmode$Mode   <- tolower(trimws(tf_withmode$Mode))
    tf_withmode$Mode[!tf_withmode$Mode %in% c("activation","repression","unknown")] <- "unknown"
    
    # Size factors from the (updated) TF table
    tf_size_factor_db  <- table(tf_withmode$TF)
    rec_size_factor_db <- table(rec_ori_data$Receptor)
    lig_size_factor_db <- table(lig_ori_data$Ligand)
    
    # ---- Build raw fold-change vector (keep sign for weighting) ----
    expres <- as.numeric(deg$foldChange)
    names(expres) <- as.character(deg$Gene)
    # Drop NA/NaN early
    keep <- !is.na(expres) & !is.nan(expres)
    expres <- expres[keep]; gene_names <- names(expres)
    
    # Replace +/-Inf with finite extrema to avoid math issues
    finite_vals <- expres[is.finite(expres)]
    if (length(finite_vals)) {
      maxv <- max(finite_vals, na.rm = TRUE)
      minv <- min(finite_vals, na.rm = TRUE)
      expres[expres ==  Inf] <- maxv
      expres[expres == -Inf] <- minv
    }
    
    # ---- Compute majority Mode PER TARGET GENE restricted to this network ----
    # Use scMlnet_results$TFTar pairs for relevance
    TFTar_pairs <- as.data.frame(stringr::str_split(scMlnet_results$TFTar, "_", simplify = TRUE), stringsAsFactors = FALSE)
    if (ncol(TFTar_pairs) != 2L) stop("Unexpected format in scMlnet_results$TFTar (expected 'TF_Target').")
    colnames(TFTar_pairs) <- c("TF","Target")
    TFTar_pairs$TF     <- toupper(trimws(TFTar_pairs$TF))
    TFTar_pairs$Target <- toupper(trimws(TFTar_pairs$Target))
    
    # Join with mode table on (TF,Target)
    pairs_mode <- merge(
      TFTar_pairs,
      tf_withmode[, c("TF","Target","Mode")],
      by = c("TF","Target"),
      all.x = TRUE
    )
    pairs_mode$Mode[is.na(pairs_mode$Mode)] <- "unknown"
    
    # Majority rule per Target: count activation vs repression; ties -> unknown
    major_mode_by_target <- aggregate(Mode ~ Target, data = pairs_mode, FUN = function(v) {
      a <- sum(v == "activation"); r <- sum(v == "repression")
      if (a > r) "activation" else if (r > a) "repression" else "unknown"
    })
    mode_lookup <- setNames(major_mode_by_target$Mode, major_mode_by_target$Target)
    
    # Map a Mode to every gene in 'expres' (by uppercase name); default unknown
    genes_uc <- toupper(gene_names)
    gene_modes <- mode_lookup[genes_uc]
    gene_modes[is.na(gene_modes)] <- "unknown"
    
    # ---- Mode-aware weighting BEFORE taking magnitude/abs ----
    # We support two input styles:
    #  (A) signed log fold-change (has negatives)  -> weight then abs()
    #  (B) ratio fold-change (>0)                  -> weight, then invert <1 to its reciprocal
    
    is_signed <- any(expres < 0, na.rm = TRUE)  # heuristic: negatives present
    weights <- rep(1, length(expres))
    
    if (is_signed) {
      pos <- expres > 0
      neg <- expres < 0
      # weights by sign x mode
      weights[pos & gene_modes == "activation"] <- 2
      weights[pos & gene_modes == "repression"] <- 0.5
      weights[neg & gene_modes == "activation"] <- 0.5
      weights[neg & gene_modes == "repression"] <- 2
      # Unknown -> *1 (already)
      expres <- abs(expres * weights)
    } else {
      # Treat <1 as "down" direction
      up    <- expres >= 1
      down  <- expres > 0 & expres < 1
      weights[up   & gene_modes == "activation"] <- 2
      weights[up   & gene_modes == "repression"] <- 0.5
      weights[down & gene_modes == "activation"] <- 0.5
      weights[down & gene_modes == "repression"] <- 2
      wexp <- expres * weights
      # Convert to magnitude (>1): invert any values still <1
      idx <- is.finite(wexp) & wexp > 0 & wexp < 1
      wexp[idx] <- 1 / wexp[idx]
      expres <- wexp
    }
    
    # Final cleanup for safety (all positive magnitudes at this point)
    expres <- as.numeric(expres)
    names(expres) <- gene_names
    expres <- expres[is.finite(expres)]
    
    
    # Detected-based normalization
    LigRec_split <- stringr::str_split(scMlnet_results$LigRec, "_", simplify = TRUE)
    RecTF_split  <- stringr::str_split(scMlnet_results$RecTF,  "_", simplify = TRUE)
    TFTar_split  <- stringr::str_split(scMlnet_results$TFTar,  "_", simplify = TRUE)
    detected_size_factor <- function(x) table(factor(x, levels = unique(x)))
    
    compute_pathway_network <- function(
    tf_size_factor, rec_size_factor, lig_size_factor, normalization_label
    ) {
      # --- REDO split objects to not get clobbered by dplyr
      scMlnet_results$LigRec <- stringr::str_split(scMlnet_results$LigRec, "_", simplify = TRUE)
      scMlnet_results$RecTF  <- stringr::str_split(scMlnet_results$RecTF,  "_", simplify = TRUE)
      scMlnet_results$TFTar  <- stringr::str_split(scMlnet_results$TFTar,  "_", simplify = TRUE)
      
      # TF->target
      # TF->target
      scMlnet_results$TFTar <- as.data.frame(scMlnet_results$TFTar)
      scMlnet_results$TFTar$V3 <- expres[scMlnet_results$TFTar$V2]
      
      tf_target <- scMlnet_results$TFTar %>%
        dplyr::group_by(V1) %>%  # V1 = TF
        dplyr::summarise(aggregate = get(TF_fun)(V3), .groups = "drop")
      
      # ---- SCENIC TF weights (optional) BEFORE TF-size normalization ----
      if (!is.null(tf_weight_vec) && length(tf_weight_vec)) {
        w <- tf_weight_vec[as.character(tf_target$V1)]
        w[is.na(w)] <- 1
        if (isTRUE(scenic_center) && any(is.finite(w))) {
          m <- mean(w[is.finite(w)], na.rm = TRUE)
          if (is.finite(m) && m > 0) w <- w / m
        }
        if (is.finite(scenic_gamma) && scenic_gamma != 1) w <- w ^ scenic_gamma
        if (is.finite(scenic_alpha) && scenic_alpha != 1) w <- 1 + scenic_alpha * (w - 1)
        if (length(scenic_clip) == 2 && all(is.finite(scenic_clip))) {
          w <- pmin(pmax(w, scenic_clip[1]), scenic_clip[2])
        }
        tf_target$aggregate <- tf_target$aggregate * w
      }
      
      tf_n <- setNames(tf_target$aggregate, tf_target$V1)
      tf_n <- tf_n / tf_size_factor[names(tf_n)]
      
      # ---------------- Receptor scores (from RecTF) ----------------
      scMlnet_results$RecTF <- as.data.frame(scMlnet_results$RecTF)
      # Map TF scores to each RecTF edge
      if (length(tf_n)) {
        tf_map <- tf_n
        idx <- match(scMlnet_results$RecTF$V2, names(tf_map))  # V2 is TF
        scMlnet_results$RecTF$tf_n <- ifelse(is.na(idx), 0, tf_map[idx])
      } else {
        scMlnet_results$RecTF$tf_n <- 0
      }
      
      rec_df <- scMlnet_results$RecTF %>%
        dplyr::group_by(V1) %>%             # V1 = Receptor
        dplyr::summarise(aggregate = get(main_fun)(.data$tf_n), .groups = "drop")
      rec_n <- setNames(rec_df$aggregate, rec_df$V1)
      rec_n <- rec_n / rec_size_factor[names(rec_n)]
      
      # ---------------- Ligand scores (from LigRec) -----------------
      scMlnet_results$LigRec <- as.data.frame(scMlnet_results$LigRec)
      # Map receptor scores to each LigRec edge
      if (length(rec_n)) {
        rec_map <- rec_n
        idx <- match(scMlnet_results$LigRec$V2, names(rec_map))  # V2 = Receptor
        scMlnet_results$LigRec$rec_n <- ifelse(is.na(idx), 0, rec_map[idx])
      } else {
        scMlnet_results$LigRec$rec_n <- 0
      }
      
      lig_df <- scMlnet_results$LigRec %>%
        dplyr::group_by(V1) %>%             # V1 = Ligand
        dplyr::summarise(score = get(main_fun)(.data$rec_n), .groups = "drop")
      lig_n <- setNames(lig_df$score, lig_df$V1)
      
      
      # ------- Ligand table + mapping to pathway -------
      lig_rank_all <- lig_n[order(lig_n, decreasing = TRUE)]
      lig_rank_all <- as.data.frame(lig_rank_all)
      colnames(lig_rank_all) <- "score"
      lig_rank_all$ligand <- rownames(lig_rank_all)
      
      lig_rank_all$pathway <- receptor_ligand$pathway_name[
        match(lig_rank_all$ligand, receptor_ligand$Ligand.ApprovedSymbol)
      ]
      lig_rank_all$score <- as.numeric(lig_rank_all$score)
      lig_rank_all <- lig_rank_all[!is.na(lig_rank_all$pathway), ]
      
      # Pathway summary
      pathway_n <- lig_rank_all %>%
        dplyr::group_by(pathway) %>%
        dplyr::summarise(score = as.numeric(get(path_fun)(as.numeric(score))), .groups = "drop") %>%
        dplyr::arrange(dplyr::desc(score))
      
      total_score <- sum(pathway_n$score)
      pathway_n$score.perc     <- (pathway_n$score / total_score) * 100
      pathway_n$score.perc.txt <- scales::percent(pathway_n$score / total_score, accuracy = 0.1)
      
      all_pathways <- names(palette)
      pathway_n$pathway    <- factor(pathway_n$pathway,    levels = all_pathways)
      lig_rank_all$pathway <- factor(lig_rank_all$pathway, levels = all_pathways)
      
      # Ligand percentages (for ligand-only treemap)
      lig_rank_all$score.perc     <- (lig_rank_all$score / sum(lig_rank_all$score)) * 100
      lig_rank_all$score.perc.txt <- scales::percent(lig_rank_all$score / sum(lig_rank_all$score), accuracy = 0.1)
      lig_rank_all$pathway        <- droplevels(lig_rank_all$pathway)
      lig_rank_all$score          <- as.numeric(lig_rank_all$score)
      lig_rank_all$score.perc     <- as.numeric(lig_rank_all$score.perc)
      
      pathway_n <- pathway_n[!is.na(pathway_n$pathway), ]
      pathway_n$pathway    <- droplevels(pathway_n$pathway)
      pathway_n$score      <- as.numeric(pathway_n$score)
      pathway_n$score.perc <- as.numeric(pathway_n$score.perc)
      
      
      # ---- Interaction (Ligand_Receptor) table ----
      # scMlnet_results$LigRec: V1=ligand, V2=receptor, rec_n column already filled above
      LR <- as.data.frame(scMlnet_results$LigRec, stringsAsFactors = FALSE)
      names(LR)[1:2] <- c("ligand","receptor")
      LR$interaction <- paste0(LR$ligand, "_", LR$receptor)
      
      # One score per pair (aggregate in case of dup edges)
      int_rank_all <- LR %>%
        dplyr::group_by(interaction) %>%
        dplyr::summarise(
          ligand   = dplyr::first(.data$ligand),
          receptor = dplyr::first(.data$receptor),
          score    = get(main_fun)(.data$rec_n),
          .groups  = "drop"
        ) %>%
        dplyr::arrange(dplyr::desc(.data$score))
      
      # Map pathway from ligand (matches your EGF_ERBB3 expectation)
      int_rank_all$pathway <- receptor_ligand$pathway_name[
        match(int_rank_all$ligand, receptor_ligand$Ligand.ApprovedSymbol)
      ]
      
      # Keep only pairs that mapped to a pathway (optional)
      int_rank_all <- int_rank_all[!is.na(int_rank_all$pathway), , drop = FALSE]
      
      # Percent-of-total for visualization
      int_rank_all$score.perc     <- (int_rank_all$score / sum(int_rank_all$score)) * 100
      int_rank_all$score.perc.txt <- sprintf("%.1f%%", int_rank_all$score.perc)
      
      # ------- NEW: Receptor table (+ pathway mapping) -------
      # ------- NEW: Receptor table (+ pathway mapping) -------
      # Build from names(rec_n) so we never rely on rownames()
      rec_rank_all <- data.frame(
        receptor = names(rec_n),
        score    = as.numeric(rec_n),
        stringsAsFactors = FALSE
      )
      rec_rank_all <- rec_rank_all[order(rec_rank_all$score, decreasing = TRUE), , drop = FALSE]
      rec_rank_all$pathway <- receptor_ligand$pathway_name[
        match(rec_rank_all$receptor, receptor_ligand$Receptor.ApprovedSymbol)
      ]
      rec_rank_all$score.perc     <- (rec_rank_all$score / sum(rec_rank_all$score)) * 100
      rec_rank_all$score.perc.txt <- sprintf("%.1f%%", rec_rank_all$score.perc)
      
      
      # ---- helpers (put once inside compute_pathway_network, before treemap calls) ----
      # ---- helpers: safe percent, stable treemap df, pre-aggregate, debug ----
      .safe_percent <- function(x) {
        x <- as.numeric(x)
        s <- sum(x, na.rm = TRUE)
        if (!is.finite(s) || s <= 0) return(rep(0, length(x)))
        (x / s) * 100
      }
      
      # Build a minimal, type-stable data.frame for treemap
      .make_tm_df <- function(df, index_cols, size_col, label_col = NULL) {
        df <- as.data.frame(df, stringsAsFactors = FALSE)
        
        # character-only indices
        for (nm in index_cols) df[[nm]] <- as.character(df[[nm]])
        
        # force vSize to "double" and sanitize
        vsize <- as.numeric(df[[size_col]]) * 1.0
        vsize[!is.finite(vsize)] <- 0
        df[[".vsize"]] <- vsize
        
        if (!is.null(label_col)) df[[label_col]] <- as.character(df[[label_col]])
        df[, c(index_cols, ".vsize", label_col), drop = FALSE]
      }
      
      # Pre-aggregate by the treemap hierarchy using base::aggregate (always returns double here)
      .preagg <- function(df, by_cols, size_col) {
        df <- as.data.frame(df, stringsAsFactors = FALSE)
        for (nm in by_cols) df[[nm]] <- as.character(df[[nm]])
        df[[size_col]] <- as.numeric(df[[size_col]]) * 1.0
        agg <- aggregate(df[[size_col]], by = df[by_cols], FUN = function(z) sum(as.numeric(z)))
        names(agg)[ncol(agg)] <- size_col
        agg
      }
      
      # Wrap treemap to print structure if anything throws
      .tm <- function(data, ...) {
        tryCatch(
          treemap::treemap(data, ...),
          error = function(e) {
            cat("\n---- treemap input that failed ----\n")
            print(utils::head(data, 6))
            str(data)
            stop(e)
          }
        )
      }
      
      
      
      # Recompute %s as clean doubles + labels
      pathway_n$score.perc     <- .safe_percent(pathway_n$score)
      pathway_n$score.perc.txt <- sprintf("%.1f%%", pathway_n$score.perc)
      
      lig_rank_all$score.perc     <- .safe_percent(lig_rank_all$score)
      lig_rank_all$score.perc.txt <- sprintf("%.1f%%", lig_rank_all$score.perc)
      
      rec_rank_all$score.perc     <- .safe_percent(rec_rank_all$score)
      rec_rank_all$score.perc.txt <- sprintf("%.1f%%", rec_rank_all$score.perc)
      
      # 1) Pathway treemap (2-level hierarchy)
      df_path_raw <- data.frame(
        pathway        = as.character(pathway_n$pathway),
        score.perc.txt = as.character(pathway_n$score.perc.txt),
        score.perc     = as.numeric(pathway_n$score.perc),
        stringsAsFactors = FALSE
      )
      # Pre-aggregate just in case there are accidental duplicates
      df_path_agg <- .preagg(df_path_raw, c("pathway","score.perc.txt"), "score.perc")
      df_path     <- .make_tm_df(df_path_agg, c("pathway","score.perc.txt"), "score.perc")
      
      # 2) Ligand-by-pathway (3-level)
      df_lbp_raw <- data.frame(
        pathway        = as.character(lig_rank_all$pathway),
        score.perc.txt = as.character(lig_rank_all$score.perc.txt),
        ligand         = as.character(lig_rank_all$ligand),
        score.perc     = as.numeric(lig_rank_all$score.perc),
        stringsAsFactors = FALSE
      )
      df_lbp_agg    <- .preagg(df_lbp_raw, c("pathway","score.perc.txt","ligand"), "score.perc")
      df_lig_by_path <- .make_tm_df(df_lbp_agg, c("pathway","score.perc.txt","ligand"), "score.perc")
      
      # 3) Ligand-only (single-level)
      df_lo_raw <- data.frame(
        ligand         = as.character(lig_rank_all$ligand),
        score.perc.txt = as.character(lig_rank_all$score.perc.txt),
        score.perc     = as.numeric(lig_rank_all$score.perc),
        stringsAsFactors = FALSE
      )
      df_lo_agg   <- .preagg(df_lo_raw, c("ligand","score.perc.txt"), "score.perc")
      df_lig_only <- .make_tm_df(df_lo_agg, c("ligand","score.perc.txt"), "score.perc")
      
      # 4) Receptor-only (single-level)
      df_ro_raw <- data.frame(
        receptor       = as.character(rec_rank_all$receptor),
        score.perc.txt = as.character(rec_rank_all$score.perc.txt),
        score.perc     = as.numeric(rec_rank_all$score.perc),
        stringsAsFactors = FALSE
      )
      df_ro_agg    <- .preagg(df_ro_raw, c("receptor","score.perc.txt"), "score.perc")
      df_rec_only  <- .make_tm_df(df_ro_agg, c("receptor","score.perc.txt"), "score.perc")
      
      
      
      # ------- Treemaps -------
      library(treemap)
      
      # Pathway treemap
      png(file.path(output_dir, sprintf("Tree_path_%s_to_%s_%s.png", receiver_label, target_label, normalization_label)), width=600, height=750)
      .tm(
        df_path,
        index   = c("pathway","score.perc.txt"),
        vSize   = ".vsize",
        type    = "categorical",
        vColor  = "pathway",
        palette = palette,
        title   = sprintf("%s --> %s", receiver_label, target_label), #"%s --> %s (%s norm)"
        border.col  = c("black","black"), border.lwds = c(3,1),
        mirror.y = mirror.y,
        fontsize.labels = c(8,8), fontcolor.labels = c("black","black"),
        fontface.labels = c(2,1), bg.labels = c("transparent"),
        align.labels = list(c("left","top"), c("right","bottom")),
        overlap.labels = 0, inflate.labels = TRUE, position.legend = "none"
      )
      dev.off()
      
      # Ligand-by-pathway treemap
      png(file.path(output_dir, sprintf("Tree_ligand_%s_to_%s_%s.png", receiver_label, target_label, normalization_label)), width=600, height=750)
      .tm(
        df_lig_by_path,
        index   = c("pathway","score.perc.txt","ligand"),
        vSize   = ".vsize",
        type    = "categorical",
        vColor  = "pathway",
        palette = palette,
        title   = sprintf("%s --> %s Ligand rank", receiver_label, target_label), #"%s --> %s Ligand rank (%s norm)"
        border.col = c("black","black","white"), border.lwds = c(3,1,3),
        mirror.y = mirror.y,
        fontsize.labels = c(8,0,5), fontcolor.labels = c("black","black","white"),
        fontface.labels = 2, bg.labels = c("transparent"),
        align.labels = list(c("left","top"), c("right","top"), c("center","bottom")),
        overlap.labels = 0, inflate.labels = TRUE, position.legend = "none"
      )
      dev.off()
      
      # Ligand-only treemap
      png(file.path(output_dir, sprintf("Tree_ligand_only_%s_to_%s_%s.png", receiver_label, target_label, normalization_label)), width=600, height=750)
      .tm(
        df_lig_only,
        index   = c("ligand","score.perc.txt"),
        vSize   = ".vsize",
        type    = "index",
        title   = sprintf("%s --> %s Ligand rank", receiver_label, target_label), #"%s --> %s Ligand rank (no pathway, %s norm)"
        border.col = c("black","black"), border.lwds = c(3,1),
        mirror.y = mirror.y,
        fontsize.labels = c(8,8), fontcolor.labels = c("black","black"),
        fontface.labels = c(2,1), bg.labels = c("transparent"),
        align.labels = list(c("left","top"), c("right","bottom")),
        overlap.labels = 0, inflate.labels = TRUE, position.legend = "none"
      )
      dev.off()
      
      # Receptor-only treemap
      png(file.path(output_dir, sprintf("Tree_receptor_only_%s_to_%s_%s.png", receiver_label, target_label, normalization_label)), width=600, height=750)
      .tm(
        df_rec_only,
        index   = c("receptor","score.perc.txt"),
        vSize   = ".vsize",
        type    = "index",
        title   = sprintf("%s --> %s Receptor rank", receiver_label, target_label), #"%s --> %s Receptor rank (no pathway, %s norm)"
        border.col = c("black","black"), border.lwds = c(3,1),
        mirror.y = mirror.y,
        fontsize.labels = c(8,8), fontcolor.labels = c("black","black"),
        fontface.labels = c(2,1), bg.labels = c("transparent"),
        align.labels = list(c("left","top"), c("right","bottom")),
        overlap.labels = 0, inflate.labels = TRUE, position.legend = "none"
      )
      dev.off()
      
      # ---- Interaction-by-pathway treemap ----
      df_ibp_raw <- data.frame(
        pathway        = as.character(int_rank_all$pathway),
        score.perc.txt = as.character(int_rank_all$score.perc.txt),
        interaction    = as.character(int_rank_all$interaction),
        score.perc     = as.numeric(int_rank_all$score.perc),
        stringsAsFactors = FALSE
      )
      df_ibp_agg     <- .preagg(df_ibp_raw, c("pathway","score.perc.txt","interaction"), "score.perc")
      df_inter_by_path <- .make_tm_df(df_ibp_agg, c("pathway","score.perc.txt","interaction"), "score.perc")
      
      png(file.path(output_dir, sprintf("Tree_interaction_%s_to_%s_%s.png", receiver_label, target_label, normalization_label)), width=600, height=750)
      .tm(
        df_inter_by_path,
        index   = c("pathway","score.perc.txt","interaction"),
        vSize   = ".vsize",
        type    = "categorical",
        vColor  = "pathway",
        palette = palette,
        title   = sprintf("%s → %s Interaction rank", receiver_label, target_label), #"%s → %s Interaction rank (%s norm)"
        border.col = c("black","black","white"), border.lwds = c(3,1,3),
        mirror.y = mirror.y,
        fontsize.labels = c(8,0,5), fontcolor.labels = c("black","black","white"),
        fontface.labels = 2, bg.labels = c("transparent"),
        align.labels = list(c("left","top"), c("right","top"), c("center","bottom")),
        overlap.labels = 0, inflate.labels = TRUE, position.legend = "none"
      )
      dev.off()
      
      # ---- Interaction-only treemap ----
      df_io_raw <- data.frame(
        interaction    = as.character(int_rank_all$interaction),
        score.perc.txt = as.character(int_rank_all$score.perc.txt),
        score.perc     = as.numeric(int_rank_all$score.perc),
        stringsAsFactors = FALSE
      )
      df_io_agg    <- .preagg(df_io_raw, c("interaction","score.perc.txt"), "score.perc")
      df_inter_only <- .make_tm_df(df_io_agg, c("interaction","score.perc.txt"), "score.perc")
      
      png(file.path(output_dir, sprintf("Tree_interaction_only_%s_to_%s_%s.png", receiver_label, target_label, normalization_label)), width=600, height=750)
      .tm(
        df_inter_only,
        index   = c("interaction","score.perc.txt"),
        vSize   = ".vsize",
        type    = "index",
        title   = sprintf("%s → %s Interaction rank", receiver_label, target_label),
        border.col = c("black","black"), border.lwds = c(3,1),
        mirror.y = mirror.y,
        fontsize.labels = c(8,8), fontcolor.labels = c("black","black"),
        fontface.labels = c(2,1), bg.labels = c("transparent"),
        align.labels = list(c("left","top"), c("right","bottom")),
        overlap.labels = 0, inflate.labels = TRUE, position.legend = "none"
      )
      dev.off()
      
      
      # ------- P-values (pathway, ligand, receptor) -------
      empirical_pvals <- function(value_vector, grouping, n_perm=1000) {
        observed <- tapply(value_vector, grouping, sum)
        n <- length(value_vector)
        null_mat <- matrix(NA, nrow=length(observed), ncol=n_perm)
        for (i in seq_len(n_perm)) {
          shuffled_group <- sample(grouping)
          null_mat[,i] <- tapply(value_vector, shuffled_group, sum)
        }
        null_means <- rowMeans(null_mat, na.rm=TRUE)
        pvals <- sapply(seq_along(observed), function(i) mean(null_mat[i, ] >= observed[i]))
        res <- data.frame(
          name = names(observed),
          observed = observed,
          null_mean = null_means,
          p.value = pvals,
          significant = ifelse(pvals < pval_cutoff, "yes", "no")
        )
        res[!is.na(res$observed), ]
      }
      pval_pathway <- empirical_pvals(lig_rank_all$score, lig_rank_all$pathway, n_perm=1000)
      pval_ligand  <- empirical_pvals(lig_rank_all$score, lig_rank_all$ligand,  n_perm=1000)
      # NEW: receptor p-values (one per receptor)
      pval_receptor <- empirical_pvals(rec_rank_all$score,  rec_rank_all$receptor, n_perm=1000)
      pval_interaction <- empirical_pvals(int_rank_all$score, int_rank_all$interaction, n_perm = 1000)
      
      CANON_COLS <- c("kind","entity","pathway","score","score.perc","p.value","significant")
      .unify <- function(kind, tbl_sum, key_col, pval_tbl, pval_key = "name", pathway_from = "pathway") {
        df <- as.data.frame(tbl_sum, stringsAsFactors = FALSE)
        df[[key_col]]   <- as.character(df[[key_col]])
        df[[pathway_from]] <- as.character(df[[pathway_from]])
        if (!is.null(pval_tbl) && nrow(pval_tbl)) {
          pv <- setNames(pval_tbl[c(pval_key,"p.value","significant")], c(key_col,"p.value","significant"))
          df <- merge(df, pv, by = key_col, all.x = TRUE)
        } else {
          df$p.value <- NA_real_; df$significant <- NA_character_
        }
        out <- data.frame(
          kind        = kind,
          entity      = df[[key_col]],
          pathway     = df[[pathway_from]],
          score       = as.numeric(df$score),
          score.perc  = as.numeric(df$score.perc),
          p.value     = suppressWarnings(as.numeric(df$p.value)),
          significant = as.character(df$significant),
          stringsAsFactors = FALSE
        )
        out <- out[CANON_COLS]
        out$score.perc <- ifelse(is.na(out$score.perc), NA, round(out$score.perc, 1))
        out
      }
      
      std <- list(
        pathway     = .unify("pathway",    pathway_n,     "pathway",    pval_pathway,    pval_key="name", pathway_from="pathway"),
        ligand      = .unify("ligand",     lig_rank_all,  "ligand",     pval_ligand,     pval_key="name", pathway_from="pathway"),
        receptor    = .unify("receptor",   rec_rank_all,  "receptor",   pval_receptor,   pval_key="name", pathway_from="pathway"),
        interaction = .unify("interaction",int_rank_all,  "interaction",pval_interaction,pval_key="name", pathway_from="pathway")
      )
      
      
      # ------- Save tables -------
      write.csv(pathway_n,      file.path(output_dir, paste0("pathway_summary_",   normalization_label, ".csv")), row.names=FALSE)
      write.csv(lig_rank_all,   file.path(output_dir, paste0("ligand_summary_",    normalization_label, ".csv")), row.names=FALSE)
      write.csv(rec_rank_all,   file.path(output_dir, paste0("receptor_summary_",  normalization_label, ".csv")), row.names=FALSE)
      write.csv(pval_pathway,   file.path(output_dir, paste0("pval_pathway_",      normalization_label, ".csv")), row.names=FALSE)
      write.csv(pval_ligand,    file.path(output_dir, paste0("pval_ligand_",       normalization_label, ".csv")), row.names=FALSE)
      write.csv(pval_receptor,  file.path(output_dir, paste0("pval_receptor_",     normalization_label, ".csv")), row.names=FALSE)
      write.csv(int_rank_all,   file.path(output_dir, paste0("interaction_summary_", normalization_label, ".csv")), row.names=FALSE)
      write.csv(pval_interaction, file.path(output_dir, paste0("pval_interaction_", normalization_label, ".csv")), row.names=FALSE)
      
      saveRDS(int_rank_all,     file.path(output_dir, paste0("interaction_summary_", normalization_label, ".rds")))
      saveRDS(pval_interaction, file.path(output_dir, paste0("pval_interaction_",    normalization_label, ".rds")))
      saveRDS(pathway_n,      file.path(output_dir, paste0("pathway_summary_",   normalization_label, ".rds")))
      saveRDS(lig_rank_all,   file.path(output_dir, paste0("ligand_summary_",    normalization_label, ".rds")))
      saveRDS(rec_rank_all,   file.path(output_dir, paste0("receptor_summary_",  normalization_label, ".rds")))
      saveRDS(pval_pathway,   file.path(output_dir, paste0("pval_pathway_",      normalization_label, ".rds")))
      saveRDS(pval_ligand,    file.path(output_dir, paste0("pval_ligand_",       normalization_label, ".rds")))
      saveRDS(pval_receptor,  file.path(output_dir, paste0("pval_receptor_",     normalization_label, ".rds")))
      
      diag_scenic <- NULL
      if (!is.null(tf_weight_vec) && length(tf_weight_vec)) {
        used_w <- tf_weight_vec[names(tf_n)]
        diag_scenic <- list(
          n_tf_weighted = sum(!is.na(used_w)),
          weight_summary = summary(used_w),
          # how much did it change TF scores?
          corr_before_after = suppressWarnings(cor(tf_target$aggregate / w, tf_target$aggregate, use="complete.obs"))
        )
      }
      
      # Return everything (now includes receptor)
      list(
        pathway_n       = pathway_n,
        lig_rank_all    = lig_rank_all,
        rec_rank_all    = rec_rank_all,
        int_rank_all    = int_rank_all,      # <-- NEW
        pval_pathway    = pval_pathway,
        pval_ligand     = pval_ligand,
        pval_receptor   = pval_receptor,
        pval_interaction= pval_interaction,  # <-- NEW
        std             = std                # <-- NEW (optional)
      )
      
    }
    
    
    # Database normalization
    result_db <- compute_pathway_network(
      tf_size_factor=tf_size_factor_db, rec_size_factor=rec_size_factor_db, lig_size_factor=lig_size_factor_db,
      normalization_label = "db"
    )
    # Detected normalization
    tf_size_factor_det  <- detected_size_factor(TFTar_split[,1])
    rec_size_factor_det <- detected_size_factor(RecTF_split[,1])
    lig_size_factor_det <- detected_size_factor(LigRec_split[,1])
    result_detected <- compute_pathway_network(
      tf_size_factor=tf_size_factor_det, rec_size_factor=rec_size_factor_det, lig_size_factor=lig_size_factor_det,
      normalization_label = "detected"
    )
    invisible(list(db=result_db, detected=result_detected))
  })
}
