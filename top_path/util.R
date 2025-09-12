normalize_cds_symbols <- function(cds, ortholog_csv = "database/mouse_to_human_orthologs.csv") {
  # Ensure gene_short_name exists
  if (is.null(rowData(cds)$gene_short_name)) {
    # try .var-like guesses: if rownames look like Ensembl, start from them
    rn <- rownames(cds)
    rowData(cds)$gene_short_name <- rn
  }
  syms <- as.character(rowData(cds)$gene_short_name)
  
  # Species guess from Ensembl pattern
  rn <- rownames(cds)
  is_mouse <- grepl("^ENSMUSG", rn[1])
  is_human <- grepl("^ENSG",    rn[1])
  
  # If symbols look empty, map Ensembl -> symbol (org.* if available)
  if (all(syms == rn)) {
    if (is_human && requireNamespace("org.Hs.eg.db", quietly = TRUE)) {
      syms <- AnnotationDbi::mapIds(org.Hs.eg.db::org.Hs.eg.db, keys = rn,
                                    keytype = "ENSEMBL", column = "SYMBOL", multiVals = "first")
      syms[is.na(syms) | syms == ""] <- rn
    } else if (is_mouse && requireNamespace("org.Mm.eg.db", quietly = TRUE)) {
      syms <- AnnotationDbi::mapIds(org.Mm.eg.db::org.Mm.eg.db, keys = rn,
                                    keytype = "ENSEMBL", column = "SYMBOL", multiVals = "first")
      syms[is.na(syms) | syms == ""] <- rn
    }
  }
  
  # If mouse → human conversion file exists, map; else just uppercase
  if (is_mouse && file.exists(ortholog_csv)) {
    map <- tryCatch(read.csv(ortholog_csv, stringsAsFactors = FALSE), error = function(e) NULL)
    if (!is.null(map) && nrow(map)) {
      # normalize map columns if present
      to_uc <- function(x) toupper(trimws(x))
      nm <- names(map)
      map[] <- lapply(map, function(col) if (is.character(col)) to_uc(col) else col)
      
      syms_uc <- to_uc(syms)         # current gene_short_name (symbols if available)
      rn_uc   <- to_uc(rownames(cds))# rownames (often Ensembl)
      
      # 1) try symbol→symbol
      msym_col <- grep("^mouse.*symbol$", nm, ignore.case = TRUE, value = TRUE)[1]
      hsym_col <- grep("^human.*symbol$", nm, ignore.case = TRUE, value = TRUE)[1]
      
      out <- rep(NA_character_, length(syms_uc))
      if (!is.na(msym_col) && !is.na(hsym_col)) {
        idx <- match(syms_uc, map[[msym_col]])
        out <- map[[hsym_col]][idx]
      }
      
      # 2) if still NA and Ensembl cols exist, try ensembl→symbol
      mens_col <- grep("^mouse.*ensembl", nm, ignore.case = TRUE, value = TRUE)[1]
      hens_col <- grep("^human.*ensembl", nm, ignore.case = TRUE, value = TRUE)[1]
      if (any(is.na(out)) && !is.na(mens_col) && !is.na(hens_col) && !is.na(hsym_col)) {
        idx2 <- match(rn_uc, map[[mens_col]])
        out2 <- map[[hsym_col]][idx2]
        fill <- is.na(out) | out == ""
        out[fill] <- out2[fill]
      }
      
      # fallback: uppercase whatever we have
      syms <- ifelse(is.na(out) | out == "", syms_uc, out)
    } else {
      syms <- toupper(syms)
    }
  } else {
    syms <- toupper(syms)
  }
  
  
  # Finalize
  rowData(cds)$gene_short_name <- syms
  # make unique just in case duplicates remain
  if (any(duplicated(rowData(cds)$gene_short_name))) {
    rowData(cds)$gene_short_name <- make.unique(rowData(cds)$gene_short_name, sep = "_")
  }
  cds
}




compute_scenic_tf_weights <- function(
    cds,
    scMlnet_results,
    cluster_col,
    receiver_cell,
    target_cell,
    beta    = 2,
    w_min   = 0.5,
    w_max   = 2,
    min_regulon_size = 5,
    nCores  = NULL,            # caller can pass a wish; we’ll sanitize
    cache_dir = NULL,          # optional AUCell rankings cache
    restrict_genes = FALSE     # optional: keep only TF targets to save RAM
) {
  suppressPackageStartupMessages({
    library(AUCell)
    library(Matrix)
  })
  
  # ---- cross-platform parallel plan (future + doFuture) ----
  suppressPackageStartupMessages({
    library(future)     # install.packages("future") if needed
    library(doFuture)   # install.packages("doFuture") if you use foreach %dopar%
  })
  
  # preserve previous plan and restore on exit (good hygiene inside functions)
  .old_plan <- future::plan()
  on.exit(future::plan(.old_plan), add = TRUE)
  
  # choose cores like before
  want_cores <- if (is.null(nCores)) max(1, parallel::detectCores() - 1) else as.integer(nCores)
  
  # foreach bridge (only needed if you use %dopar%)
  doFuture::registerDoFuture()
  
  # Select best backend per OS:
  # - Windows: multisession (processes)
  # - Linux (e.g., Posit Connect): multicore if supported, else multisession
  if (.Platform$OS.type == "windows") {
    future::plan(multisession, workers = want_cores)
  } else {
    if (future::supportsMulticore()) {
      future::plan(multicore, workers = want_cores)
    } else {
      future::plan(multisession, workers = want_cores)
    }
  }
  
  
  # --- subset to receiver+target and keep sparse ---
  expr <- counts(cds)
  if (!inherits(expr, "dgCMatrix")) expr <- Matrix(expr, sparse = TRUE)
  
  cl <- as.character(colData(cds)[[cluster_col]])
  cells_rec_idx <- which(cl == receiver_cell)
  cells_tar_idx <- which(cl == target_cell)
  if (length(cells_rec_idx) == 0L || length(cells_tar_idx) == 0L) return(NULL)
  
  keep_cols <- c(cells_rec_idx, cells_tar_idx)
  expr   <- expr[, keep_cols, drop = FALSE]
  cl_sub <- cl[keep_cols]
  colnames(expr) <- colnames(cds)[keep_cols]
  
  # --- gene symbols (uppercase), drop NA/dups to keep memory tight ---
  gs <- toupper(as.character(rowData(cds)$gene_short_name))
  names(gs) <- rownames(cds)
  rownames(expr) <- gs[rownames(expr)]
  ok <- !is.na(rownames(expr)) & nzchar(rownames(expr))
  expr <- expr[ok, , drop = FALSE]
  if (any(duplicated(rownames(expr)))) expr <- expr[!duplicated(rownames(expr)), , drop = FALSE]
  
  # --- regulons from your TF→Target edges ---
  TFt <- as.data.frame(stringr::str_split(scMlnet_results$TFTar, "_", simplify = TRUE),
                       stringsAsFactors = FALSE)
  if (ncol(TFt) != 2) stop("Unexpected scMlnet_results$TFTar format.")
  names(TFt) <- c("TF", "Target")
  TFt$TF     <- toupper(TFt$TF)
  TFt$Target <- toupper(TFt$Target)
  
  if (isTRUE(restrict_genes)) {
    keep_genes <- intersect(unique(TFt$Target), rownames(expr))
    if (length(keep_genes) >= 500L) {  # keep a reasonable floor
      expr <- expr[keep_genes, , drop = FALSE]
    }
  }
  
  regs <- split(TFt$Target, TFt$TF)
  regs <- lapply(regs, function(g) intersect(g, rownames(expr)))
  regs <- regs[sapply(regs, length) >= min_regulon_size]
  if (length(regs) == 0L) return(NULL)
  
  # --- AUCell rankings (with simple cache) ---
  rankings <- NULL
  if (!is.null(cache_dir)) {
    dir.create(cache_dir, recursive = TRUE, showWarnings = FALSE)
    cache_file <- file.path(cache_dir, "aucell_rankings_subset.rds")
    if (file.exists(cache_file)) {
      tmp <- try(readRDS(cache_file), silent = TRUE)
      if (!inherits(tmp, "try-error") && is.list(tmp) && !is.null(tmp$rankings)) {
        rn_ok <- identical(rownames(tmp$rankings), rownames(expr))
        cn_ok <- identical(colnames(tmp$rankings), colnames(expr))
        if (isTRUE(rn_ok) && isTRUE(cn_ok)) rankings <- tmp$rankings
      }
    }
  }
  
  if (is.null(rankings)) {
    # AUCell warns that nCores is deprecated; we pass 1 or a safe value and silence the warning.
    rankings <- suppressWarnings(
      AUCell_buildRankings(expr, nCores = run_cores, plotStats = FALSE, verbose = FALSE)
    )
    if (!is.null(cache_dir)) saveRDS(list(rankings = rankings), file = cache_file)
  }
  
  # AUC (use safe #cores; if run_cores==1 AUCell won't attempt doMC)
  auc <- suppressWarnings(
    AUCell_calcAUC(regs, rankings,
                   aucMaxRank = ceiling(0.05 * nrow(rankings)),
                   nCores = run_cores, verbose = FALSE)
  )
  aucMat <- as.matrix(getAUC(auc))   # rows=TF(regulon), cols=cells
  
  cells_rec_names <- colnames(expr)[cl_sub == receiver_cell]
  cells_tar_names <- colnames(expr)[cl_sub == target_cell]
  
  mean_rec <- rowMeans(aucMat[, cells_rec_names, drop = FALSE], na.rm = TRUE)
  mean_tar <- rowMeans(aucMat[, cells_tar_names, drop = FALSE], na.rm = TRUE)
  
  # z-score of target-vs-receiver delta → weight in [w_min, w_max]
  delta <- mean_tar - mean_rec
  mu  <- mean(delta, na.rm = TRUE); sdv <- stats::sd(delta, na.rm = TRUE)
  z <- if (!is.finite(sdv) || sdv == 0) rep(0, length(delta)) else (delta - mu) / sdv
  names(z) <- names(delta)
  
  w <- exp(beta * z)
  w[!is.finite(w)] <- 1
  w <- pmin(pmax(w, w_min), w_max)
  names(w) <- names(delta)
  w
}





abbreviate_clusters <- function(cluster_names,
                                overrides = NULL,
                                ignore_words = c("cell","cells","the","of","and")) {
  cluster_names <- as.character(cluster_names)
  
  # --- base abbreviation maker ---
  make_abbrev <- function(x) {
    s <- trimws(x)
    # split on underscores, hyphens, slashes, and CamelCase
    s <- gsub("_|\\-|/", " ", s)
    s <- gsub("([a-z])([A-Z])", "\\1 \\2", s, perl = TRUE)
    toks <- unlist(strsplit(s, "\\s+"))
    toks <- toks[nzchar(toks)]
    if (!length(toks)) return(toupper(gsub("\\s+", "", x)))
    
    # drop ignorable words
    keep <- !tolower(toks) %in% ignore_words
    toks <- toks[keep]
    if (!length(toks)) toks <- unlist(strsplit(s, "\\s+"))
    
    parts <- vapply(toks, function(tok) {
      # preserve tokens like CD4, SP6+, 09_11w as-is (uppercased)
      if (grepl("[0-9+_]", tok)) return(toupper(tok))
      # short tokens (<=3) keep whole; otherwise take initial
      if (nchar(tok) <= 3) toupper(tok) else toupper(substr(tok, 1, 1))
    }, character(1))
    
    ab <- paste0(parts, collapse = "")
    # ensure at least 2 chars
    if (nchar(ab) < 2) ab <- toupper(substr(gsub("\\s+", "", x), 1, 2))
    ab
  }
  
  # 1) base abbreviations
  abbr <- vapply(cluster_names, make_abbrev, character(1))
  
  # 2) apply overrides (named vector: names=full, values=abbr)
  if (!is.null(overrides)) {
    m <- match(cluster_names, names(overrides))
    hit <- !is.na(m)
    abbr[hit] <- overrides[m[hit]]
  }
  
  # 3) ensure uniqueness; prefer NOT to change overridden ones
  fixed <- if (!is.null(overrides)) cluster_names %in% names(overrides) else rep(FALSE, length(cluster_names))
  df <- data.frame(full = cluster_names, abbr = abbr, fixed = fixed, stringsAsFactors = FALSE)
  
  dup_keys <- unique(df$abbr[duplicated(df$abbr)])
  for (a in dup_keys) {
    idx <- which(df$abbr == a)
    # put fixed (overridden) first so we don't rename it
    ord <- order(!df$fixed[idx])  # fixed=TRUE comes first
    idx <- idx[ord]
    if (length(idx) > 1) {
      # keep first as-is; suffix the rest
      for (k in seq_along(idx)[-1]) {
        df$abbr[idx[k]] <- paste0(a, "_", k)
      }
    }
  }
  
  rownames(df) <- NULL
  df[, c("full","abbr")]
}


abbreviate_clusters_dict <- function(
    cluster_names,
    overrides = NULL,
    ignore_words = c("cell","cells","the","of","and"),
    # domain words: treated as "meaningful" (initials), extend as you like
    domain_words = c("epithelium","enamel","knot","outer","inner","stratum",
                     "intermedium","stellate","reticulum","pre","early","secretory",
                     "mesenchyme","ectomesenchyme","ecto","ameloblast","odontoblast",
                     "papilla","follicle")
) {
  cluster_names <- as.character(cluster_names)
  
  # --- dictionary backend (hunspell > qdapDictionaries > small built-in) ---
  is_word <- local({
    if (requireNamespace("hunspell", quietly = TRUE)) {
      dict <- hunspell::dictionary("en_US")
      function(tok) {
        w <- tolower(gsub("[^A-Za-z]", "", tok))
        if (!nzchar(w)) return(FALSE)
        w %in% tolower(domain_words) || hunspell::hunspell_check(w, dict = dict)
      }
    } else if (requireNamespace("qdapDictionaries", quietly = TRUE)) {
      lex <- tolower(qdapDictionaries::GradyAugmented)
      function(tok) {
        w <- tolower(gsub("[^A-Za-z]", "", tok))
        if (!nzchar(w)) return(FALSE)
        w %in% tolower(domain_words) || (w %in% lex)
      }
    } else {
      # small fallback lexicon (extend domain_words above to improve coverage)
      base_lex <- tolower(unique(c(
        "oral","dental","outer","inner","stratum","intermedium","stellate",
        "reticulum","pre","early","secretory","knot","epithelium","mesenchyme",
        domain_words
      )))
      function(tok) {
        w <- tolower(gsub("[^A-Za-z]", "", tok))
        if (!nzchar(w)) return(FALSE)
        w %in% base_lex
      }
    }
  })
  
  # split: spaces, _, -, /, and camelCase; special biomed rule for 'ectomesenchyme'
  split_tokens <- function(s) {
    s <- gsub("_|\\-|/", " ", s)
    s <- gsub("([a-z])([A-Z])", "\\1 \\2", s, perl = TRUE)
    toks <- unlist(strsplit(s, "\\s+"))
    toks <- toks[nzchar(toks)]
    # split 'ectomesenchyme' -> 'ecto', 'mesenchyme'
    toks <- unlist(lapply(toks, function(tk) {
      if (grepl("^ecto.*mesenchyme$", tolower(tk))) c("ecto","mesenchyme") else tk
    }))
    toks
  }
  
  make_one <- function(x) {
    s <- trimws(x)
    toks <- split_tokens(s)
    if (!length(toks)) return(toupper(gsub("\\s+", "", s)))
    
    # remove ignorable words
    keep <- !tolower(toks) %in% ignore_words
    toks <- toks[keep]; if (!length(toks)) toks <- split_tokens(s)
    
    parts <- vapply(toks, function(tok) {
      # tokens with digits/plus/colon likely IDs/genes -> keep
      if (grepl("[0-9+:]", tok)) return(toupper(tok))
      # tokens with '+' or '_' (markers) -> keep
      if (grepl("[+_]", tok)) return(toupper(tok))
      # short all-caps like OB, PA -> keep
      if (nchar(tok) <= 3 && grepl("^[A-Z]+$", tok)) return(tok)
      # dictionary (or domain) word -> initial
      if (is_word(tok)) return(toupper(substr(tok, 1, 1)))
      # otherwise keep token (likely an abbrev/gene)
      toupper(tok)
    }, character(1))
    
    ab <- paste0(parts, collapse = "")
    if (nchar(ab) < 2) ab <- toupper(substr(gsub("\\s+", "", s), 1, 2))
    ab
  }
  
  # base pass
  abbr <- vapply(cluster_names, make_one, character(1))
  
  overrides = c("Oral epithelium"="OE",
                "Dental epithelium"="DE",
                "Dental ectomesenchyme"="DEM",
                "Odontoblasts"="OB",
                "Odontoblast"="OB",
                "Ameloblasts"="AM",
                "Secretory Ameloblasts"="sAM",
                "Early Ameloblasts"="eAM",
                "pre-ameloblasts"="pAM")
  
  # overrides (case/space-insensitive)
  if (!is.null(overrides)) {
    norm_key <- function(z) tolower(gsub("\\s+", " ", trimws(z)))
    ov <- overrides; names(ov) <- norm_key(names(ov))
    key <- norm_key(cluster_names)
    hit <- key %in% names(ov)
    abbr[hit] <- ov[key[hit]]
  }
  
  # ensure uniqueness, do not alter overridden ones
  fixed <- rep(FALSE, length(cluster_names))
  if (!is.null(overrides)) {
    norm_key <- function(z) tolower(gsub("\\s+", " ", trimws(z)))
    fixed <- norm_key(cluster_names) %in% names(overrides)
  }
  df <- data.frame(full = cluster_names, abbr = abbr, fixed = fixed, stringsAsFactors = FALSE)
  dvals <- unique(df$abbr[duplicated(df$abbr)])
  for (a in dvals) {
    idx <- which(df$abbr == a)
    idx <- idx[order(!df$fixed[idx])]  # keep overridden one first
    if (length(idx) > 1) for (k in seq_along(idx)[-1]) df$abbr[idx[k]] <- paste0(a, "_", k)
  }
  df[, c("full","abbr")]
}



# top_path/util.R

# map must be a data.frame with columns: full, abbr
label_for_map <- function(full, map = NULL) {
  if (is.null(map)) return(full)
  i <- match(full, map$full)
  out <- ifelse(is.na(i), full, map$abbr[i])
  out
}

short_file_label_map <- function(full, map = NULL) {
  lbl <- label_for_map(full, map)
  if (is.null(lbl) || is.na(lbl) || !nzchar(lbl)) lbl <- full
  lbl <- gsub("\\s+", "", lbl)                         # drop spaces
  lbl <- gsub("[^A-Za-z0-9+_.-]", "_", lbl)           # safe for filenames
  lbl
}






.force_double <- function(x) {
  y <- suppressWarnings(as.numeric(x))
  y[!is.finite(y)] <- 0
  y * 1.0  # ensure "double" type
}

.preagg_double <- function(df, by_cols, size_col) {
  df <- as.data.frame(df, stringsAsFactors = FALSE)
  for (nm in by_cols) df[[nm]] <- as.character(df[[nm]])
  df[[size_col]] <- .force_double(df[[size_col]])
  agg <- stats::aggregate(df[[size_col]],
                          by = df[by_cols],
                          FUN = function(z) sum(as.numeric(z)))
  names(agg)[ncol(agg)] <- size_col
  agg
}



.open_png <- function(file, w = 7, h = 8.75, res = 220, pointsize = 10) {
  if (requireNamespace("ragg", quietly = TRUE)) {
    ragg::agg_png(filename = file, width = w, height = h, units = "in", res = res)
  } else {
    # cairo gives anti-aliased text/lines on all OSes
    png(filename = file, width = w, height = h, units = "in",
        res = res, pointsize = pointsize, type = "cairo")
  }
}





# --- helpers ---------------------------------------------------------------

.pairs_to_df <- function(x, left, right, sep = "_") {
  if (is.null(x)) {
    return(data.frame(setNames(list(character(0), character(0)), c(left, right)),
                      stringsAsFactors = FALSE))
  }
  if (is.data.frame(x)) {
    x <- x[, c(left, right), drop = FALSE]
    x[[left]]  <- as.character(x[[left]])
    x[[right]] <- as.character(x[[right]])
    return(x[nzchar(x[[left]]) & nzchar(x[[right]]), , drop = FALSE])
  }
  x <- as.character(x)
  parts <- strsplit(x, sep, fixed = TRUE)
  L <- vapply(parts, function(p) if (length(p)) p[1] else "", "", USE.NAMES = FALSE)
  R <- vapply(parts, function(p) if (length(p) >= 2) paste(p[-1], collapse = sep) else "", "", USE.NAMES = FALSE)
  out <- data.frame(setNames(list(trimws(L), trimws(R)), c(left, right)), stringsAsFactors = FALSE)
  out[nzchar(out[[left]]) & nzchar(out[[right]]), , drop = FALSE]
}

# Filter to a single pathway using your summary tables
# lig_tbl: result$db$lig_rank_all  (needs columns ligand, pathway[, score])
# rec_tbl: result$db$rec_rank_all  (needs columns receptor, pathway[, score])
# Use the receptor_ligand database to keep ONLY LR pairs that are annotated
# to the requested pathway. Then restrict RT and TT accordingly.
.filter_layers_for_pathway_pairs <- function(LR, RT, TT, pathway, rl_db) {
  toUC <- function(z) toupper(trimws(as.character(z)))
  
  db <- rl_db[rl_db$pathway_name == pathway,
              c("Ligand.ApprovedSymbol","Receptor.ApprovedSymbol")]
  names(db) <- c("ligand","receptor")
  db$ligand   <- toUC(db$ligand)
  db$receptor <- toUC(db$receptor)
  
  if (!nrow(LR)) return(list(LR=LR[0,], RT=RT[0,], TT=TT[0,]))
  
  LR$ligand   <- toUC(LR$ligand)
  LR$receptor <- toUC(LR$receptor)
  
  key_db <- paste(db$ligand, db$receptor)
  keep   <- paste(LR$ligand, LR$receptor) %in% key_db
  LR     <- unique(LR[keep, , drop = FALSE])
  
  if (!nrow(LR)) return(list(LR=LR, RT=RT[0,], TT=TT[0,]))
  
  RT <- RT[RT$receptor %in% unique(LR$receptor), , drop = FALSE]
  TT <- TT[TT$tf       %in% unique(RT$tf),       , drop = FALSE]
  
  list(LR = LR, RT = unique(RT), TT = unique(TT))
}



# Edge weights for drawing (use pathway-specific scores when available)
.decorate_weights <- function(LR, RT, TT, pathway, lig_tbl, rec_tbl) {
  if (nrow(LR)) {
    lt <- lig_tbl[lig_tbl$pathway == pathway, , drop = FALSE]
    sel <- intersect(c("ligand","score"), names(lt))
    if (!length(sel)) lt <- data.frame(ligand=character(0), score=numeric(0))
    wlr <- merge(LR, lt[, sel, drop = FALSE], by = "ligand", all.x = TRUE)
    mx  <- suppressWarnings(max(wlr$score, na.rm = TRUE)); if (!is.finite(mx)) mx <- 1
    LR$w <- ifelse(is.na(wlr$score), 0.5, pmax(0, wlr$score / mx))
  } else LR$w <- numeric(0)
  
  if (nrow(RT)) {
    rt <- rec_tbl[rec_tbl$pathway == pathway, , drop = FALSE]
    sel <- intersect(c("receptor","score"), names(rt))
    if (!length(sel)) rt <- data.frame(receptor=character(0), score=numeric(0))
    wrt <- merge(RT, rt[, sel, drop = FALSE], by = "receptor", all.x = TRUE)
    mx  <- suppressWarnings(max(wrt$score, na.rm = TRUE)); if (!is.finite(mx)) mx <- 1
    RT$w <- ifelse(is.na(wrt$score), 0.5, pmax(0, wrt$score / mx))
  } else RT$w <- numeric(0)
  
  TT$w <- if (nrow(TT)) rep(0.4, nrow(TT)) else numeric(0)
  list(LR = LR, RT = RT, TT = TT)
}


# Build an igraph from three layers with widths from 'w'
.build_explain_graph <- function(LR, RT, TT) {
  suppressPackageStartupMessages(library(igraph))
  edges <- rbind(
    if (nrow(LR)) data.frame(from = LR$ligand,   to = LR$receptor, w = LR$w,  kind = "LR") else NULL,
    if (nrow(RT)) data.frame(from = RT$receptor, to = RT$tf,       w = RT$w,  kind = "RT") else NULL,
    if (nrow(TT)) data.frame(from = TT$tf,       to = TT$target,   w = TT$w,  kind = "TT") else NULL
  )
  if (is.null(edges) || !nrow(edges)) return(make_empty_graph())
  g <- graph_from_data_frame(edges, directed = TRUE)
  E(g)$width <- 1 + 6 * pmin(1, pmax(0, edges$w))   # nicer dynamic range
  V(g)$layer <- ifelse(V(g)$name %in% LR$ligand, "ligand",
                       ifelse(V(g)$name %in% LR$receptor, "receptor",
                              ifelse(V(g)$name %in% RT$tf, "tf", "target")))
  g
}











# ---------- GNN-LIKE PATHWAY + EXPLANATION SUBGRAPHS ----------
plot_gnn_pathway_treemap <- function(gnn_df, out_png, title,
                                     palette, mirror.y = FALSE) {
  df <- as.data.frame(gnn_df, stringsAsFactors = FALSE)
  df$pathway <- as.character(df$pathway)
  df$score   <- suppressWarnings(as.numeric(df$score))
  df <- df[!is.na(df$pathway) & nzchar(df$pathway) & df$score > 0, , drop = FALSE]
  
  tm_df <- data.frame(
    pathway = df$pathway,
    score   = df$score,
    stringsAsFactors = FALSE
  )
  tot <- sum(tm_df$score, na.rm = TRUE); if (!is.finite(tot) || tot == 0) tot <- 1
  tm_df$pct <- sprintf("%.1f%%", 100 * tm_df$score / tot)
  
  dir.create(dirname(out_png), recursive = TRUE, showWarnings = FALSE)
  # old (breaks): png(out_png, width=600, height=520, dpi=300)
  png(out_png, width = 600, height = 750)
  
  # canonical palette (named)
  pal <- c(
    NOTCH="#9AA556", WNT="#E98FC5", ncWNT="#64B1DF", EGF="#66E1E1",
    IGF="#D5AB90", NRG="#8B9D63", BMP="#E88E50", ACTIVIN="#C2E579",
    TGFb="#DBC6DE", NT="#DE8D9F", VEGF="#AAD8E1", HH="#6EE07F",
    FGF="#8580DF", HGF="#D9E83D", PDGF="#62E4B7", GDNF="#E2E2A1"
  )
  # make an explicit color column (fall back to grey if unknown)
  tm_df$col <- unname(pal[tm_df$pathway]); tm_df$col[is.na(tm_df$col)] <- "#CCCCCC"
  
  treemap::treemap(
    tm_df,
    index   = c("pathway","pct"),
    vSize   = "score",
    type    = "color",
    vColor  = "col",                    # <-- use the explicit hex column
    title   = title,
    border.col  = c("black","black"),
    border.lwds = c(3,1),
    mirror.y    = mirror.y,
    fontsize.labels    = c(9,5),
    fontcolor.labels   = c("black","black"),
    fontface.labels    = c(2,1),
    align.labels       = list(c("left","top"), c("right","bottom")),
    bg.labels          = 0,
    position.legend    = "none",
    inflate.labels     = TRUE,
    lowerbound.cex.labels = 0.25
  )
  dev.off()
  
}







# --- helpers ---
.split_pairs <- function(x, nm=c("A","B")) {
  x <- as.character(x); x <- x[nzchar(x)]
  if (!length(x)) return(data.frame(setNames(list(character(0), character(0)), nm),
                                    stringsAsFactors = FALSE))
  sp <- strsplit(x, "_", fixed = TRUE)
  ok <- vapply(sp, function(v) length(v) == 2 && all(nzchar(v)), logical(1))
  sp <- sp[ok]
  if (!length(sp)) return(data.frame(setNames(list(character(0), character(0)), nm),
                                     stringsAsFactors = FALSE))
  A <- vapply(sp, `[[`, "", 1)
  B <- vapply(sp, `[[`, "", 2)
  out <- data.frame(A=A, B=B, stringsAsFactors = FALSE)
  names(out) <- nm
  out
}

.norm01 <- function(x) {
  x <- as.numeric(x)
  if (!length(x) || all(!is.finite(x))) return(rep(0.5, length(x)))
  r <- range(x, na.rm = TRUE)
  if (!is.finite(r[1]) || r[1] == r[2]) return(rep(0.5, length(x)))
  (x - r[1]) / (r[2] - r[1])
}


# -------- FULL 4-LAYER EXPLANATION SUBGRAPH --------
# Build magnitude for targets from DEG (abs(logFC) or invert FC<1)
.target_mag_from_deg <- function(deg) {
  v <- suppressWarnings(as.numeric(deg$foldChange))
  nm <- as.character(deg$Gene)
  names(v) <- nm
  v <- v[is.finite(v)]
  v <- abs(v)
  # collapse duplicates by max weight
  tapply(v, names(v), max)
}

# choose ≤ k_max reps per cluster to maximize weighted TF-target coverage
prune_by_within_cluster_cover <- function(emb_df, targets_by_pathway, w_tg, k_max = 2) {
  emb_df$cluster <- factor(emb_df$cluster)
  cov_score <- function(sel) {
    sel <- intersect(sel, names(targets_by_pathway))
    if (!length(sel)) return(0)
    tg <- unique(unlist(targets_by_pathway[sel], use.names = FALSE))
    if (!length(tg)) return(0)
    sum(w_tg[intersect(names(w_tg), tg)], na.rm = TRUE)
  }
  pieces <- lapply(split(emb_df$pathway, emb_df$cluster), function(pws) {
    pws <- intersect(pws, names(targets_by_pathway))
    if (!length(pws)) return(character(0))
    
    # best single
    s1 <- sapply(pws, cov_score)
    best1 <- names(which.max(s1))
    
    if (k_max <= 1 || length(pws) == 1) return(best1)
    
    # best pair (only if it improves over best1)
    if (length(pws) >= 2) {
      pairs <- combn(pws, 2, simplify = FALSE)
      s2 <- sapply(pairs, cov_score)
      best2 <- pairs[[which.max(s2)]]
      if (cov_score(best2) > cov_score(best1)) return(best2)
    }
    best1
  })
  
  sel <- unlist(pieces, use.names = FALSE)
  data.frame(
    cluster = emb_df$cluster[match(sel, emb_df$pathway)],
    pathway = sel,
    stringsAsFactors = FALSE
  )
}

# greedy top-K selector on sets (for LR pairs)
greedy_cover_k <- function(sets, w, k_max = 10) {
  chosen <- character(0); covered <- character(0)
  w <- w[is.finite(w)]
  for (i in seq_len(k_max)) {
    best <- NA_character_; best_gain <- 0
    for (nm in names(sets)) {
      if (nm %in% chosen) next
      tg <- setdiff(sets[[nm]], covered)
      gain <- sum(w[intersect(names(w), tg)], na.rm = TRUE)
      if (gain > best_gain) { best_gain <- gain; best <- nm }
    }
    if (is.na(best) || best_gain <= 0) break
    chosen  <- c(chosen, best)
    covered <- union(covered, sets[[best]])
  }
  chosen
}
# TF scores: sum of target magnitudes per TF (optionally * SCENIC weights)
.tf_scores_from_scmlnet <- function(scMlnet_results, targ_mag, tf_weight_vec = NULL) {
  # scMlnet_results$TFTar is "TF_Target"
  parts <- as.data.frame(stringr::str_split(scMlnet_results$TFTar, "_", simplify = TRUE),
                         stringsAsFactors = FALSE)
  if (ncol(parts) != 2L) return(setNames(numeric(0), character(0)))
  names(parts) <- c("TF","Target")
  parts$Target <- toupper(parts$Target); parts$TF <- toupper(parts$TF)
  parts$mag <- targ_mag[parts$Target]; parts$mag[is.na(parts$mag)] <- 0
  tf <- tapply(parts$mag, parts$TF, sum)
  tf <- as.numeric(tf); names(tf) <- names(tapply(parts$mag, parts$TF, sum))
  if (!is.null(tf_weight_vec) && length(tf_weight_vec)) {
    w <- tf_weight_vec[names(tf)]; w[is.na(w)] <- 1
    # re-center to mean 1 to keep scale stable
    w <- w / mean(w, na.rm = TRUE)
    tf <- tf * w
  }
  tf
}






# NEW full 4-layer explanation (below), so we also reshape deep_learning_run
# ========================= deep_learning_run ===============================
deep_learning_run <- function(
    pathway_n,
    int_rank_all,
    rec_rank_all,
    scMlnet_results,
    deg,
    tf_weight_vec = NULL,
    out_dir,
    receiver,
    target,
    norm_label = "db",
    topK = 3,
    edges_per_path = 20,
    k_lr = edges_per_path,      # NEW: cap LR edges per pathway
    k_rt = 40, k_tft = 80,
    q_lr = 0.80, q_rt = 0.80,   # NEW: quantile thresholds
    q_tt = 0.85,
    lig_rank_all = NULL
) 
{
  q_lr <- 0; k_lr <- Inf   # no ligand/receptor trimming
  q_rt <- 0; q_tt <- 0     # let .keep_top_tf_targets do the only trimming
  
  # ---- helpers (local, so this function is fully self-contained) ----------
  sanitize <- function(x) {
    x <- gsub("\\s+", "", x)
    gsub("[^A-Za-z0-9+_.-]", "_", x)
  }
  .pairs_to_df <- function(x, left, right, sep = "_") {
    if (is.null(x)) {
      return(data.frame(setNames(list(character(0), character(0)), c(left, right)),
                        stringsAsFactors = FALSE))
    }
    if (is.data.frame(x) && all(c(left, right) %in% names(x))) {
      out <- x[, c(left, right), drop = FALSE]
      out[[left]]  <- as.character(out[[left]])
      out[[right]] <- as.character(out[[right]])
      return(out[nzchar(out[[left]]) & nzchar(out[[right]]), , drop = FALSE])
    }
    x <- as.character(x)
    parts <- strsplit(x, sep, fixed = TRUE)
    L <- vapply(parts, function(p) if (length(p)) p[1] else "", "", USE.NAMES = FALSE)
    R <- vapply(parts, function(p) if (length(p) >= 2) paste(p[-1], collapse = sep) else "", "", USE.NAMES = FALSE)
    out <- data.frame(setNames(list(trimws(L), trimws(R)), c(left, right)), stringsAsFactors = FALSE)
    out[nzchar(out[[left]]) & nzchar(out[[right]]), , drop = FALSE]
  }
  .derive_lig_tbl <- function(lig_tbl, int_tbl) {
    if (!is.null(lig_tbl) && nrow(lig_tbl)) {
      want <- intersect(c("ligand","pathway","score"), names(lig_tbl))
      out  <- lig_tbl[, want, drop = FALSE]
      if (!"score" %in% names(out)) out$score <- 1
      return(out)
    }
    if (!is.null(int_tbl) && nrow(int_tbl)) {
      want <- intersect(c("ligand","pathway","score"), names(int_tbl))
      tmp  <- int_tbl[, want, drop = FALSE]
      if (!"score" %in% names(tmp)) tmp$score <- 1
      agg  <- aggregate(tmp$score, by = list(ligand = tmp$ligand, pathway = tmp$pathway), FUN = sum)
      names(agg)[3] <- "score"
      return(agg)
    }
    data.frame(ligand = character(0), pathway = character(0), score = numeric(0), stringsAsFactors = FALSE)
  }
  .derive_rec_tbl <- function(rec_tbl, int_tbl) {
    if (!is.null(rec_tbl) && nrow(rec_tbl)) {
      want <- intersect(c("receptor","pathway","score"), names(rec_tbl))
      out  <- rec_tbl[, want, drop = FALSE]
      if (!"score" %in% names(out)) out$score <- 1
      return(out)
    }
    if (!is.null(int_tbl) && nrow(int_tbl)) {
      want <- intersect(c("receptor","pathway","score"), names(int_tbl))
      tmp  <- int_tbl[, want, drop = FALSE]
      if (!"score" %in% names(tmp)) tmp$score <- 1
      agg  <- aggregate(tmp$score, by = list(receptor = tmp$receptor, pathway = tmp$pathway), FUN = sum)
      names(agg)[3] <- "score"
      return(agg)
    }
    data.frame(receptor = character(0), pathway = character(0), score = numeric(0), stringsAsFactors = FALSE)
  }
  .filter_layers_for_pathway <- function(LR, RT, TT, pathway, lig_tbl, rec_tbl) {
    keep_lig <- unique(as.character(lig_tbl$ligand   [lig_tbl$pathway == pathway]))
    keep_rec <- unique(as.character(rec_tbl$receptor [rec_tbl$pathway == pathway]))
    LR <- LR[LR$ligand %in% keep_lig & LR$receptor %in% keep_rec, , drop = FALSE]
    RT <- RT[RT$receptor %in% keep_rec, , drop = FALSE]
    TT <- TT[TT$tf %in% RT$tf, , drop = FALSE]
    list(LR = unique(LR), RT = unique(RT), TT = unique(TT))
  }
  .decorate_weights <- function(LR, RT, TT, pathway, lig_tbl, rec_tbl) {
    if (nrow(LR)) {
      wlr <- merge(LR, lig_tbl[lig_tbl$pathway == pathway, c("ligand","score")], by = "ligand", all.x = TRUE)
      mx  <- suppressWarnings(max(wlr$score, na.rm = TRUE)); if (!is.finite(mx)) mx <- 1
      LR$w <- ifelse(is.na(wlr$score), 0.5, pmax(0, wlr$score / mx))
    } else LR$w <- numeric(0)
    if (nrow(RT)) {
      wrt <- merge(RT, rec_tbl[rec_tbl$pathway == pathway, c("receptor","score")], by = "receptor", all.x = TRUE)
      mx  <- suppressWarnings(max(wrt$score, na.rm = TRUE)); if (!is.finite(mx)) mx <- 1
      RT$w <- ifelse(is.na(wrt$score), 0.5, pmax(0, wrt$score / mx))
    } else RT$w <- numeric(0)
    TT$w <- if (nrow(TT)) rep(0.4, nrow(TT)) else numeric(0)
    list(LR = LR, RT = RT, TT = TT)
  }
  .build_graph <- function(LR, RT, TT) {
    suppressPackageStartupMessages(library(igraph))
    .tag <- function(x, L) if (length(x)) paste0(x, "||", L) else character(0)
    
    e_lr <- if (nrow(LR)) data.frame(
      from = .tag(LR$ligand,   "Ligand"),
      to   = .tag(LR$receptor, "Receptor"),
      w=LR$w, kind="LR", stringsAsFactors=FALSE) else NULL
    
    e_rt <- if (nrow(RT)) data.frame(
      from = .tag(RT$receptor, "Receptor"),
      to   = .tag(RT$tf,       "TF"),
      w=RT$w, kind="RT", stringsAsFactors=FALSE) else NULL
    
    e_tt <- if (nrow(TT)) data.frame(
      from = .tag(TT$tf,       "TF"),
      to   = .tag(TT$target,   "Target"),
      w=TT$w, kind="TT", stringsAsFactors=FALSE) else NULL
    
    edges <- rbind(e_lr, e_rt, e_tt)
    if (is.null(edges) || !nrow(edges)) return(igraph::make_empty_graph())
    
    edges <- edges[nzchar(edges$from) & nzchar(edges$to), , drop = FALSE]
    
    vnames <- unique(c(edges$from, edges$to))
    parts  <- strsplit(vnames, "\\|\\|")
    nodes  <- data.frame(
      name  = vnames,
      layer = vapply(parts, function(p) tolower(p[2]), ""),
      label = vapply(parts, `[`, "", 1),
      stringsAsFactors = FALSE
    )
    
    g <- igraph::graph_from_data_frame(edges, vertices = nodes, directed = TRUE)
    E(g)$width <- 1 + 6 * pmin(1, pmax(0, edges$w))
    g
  }
  
  
  
  .plot_explain <- function(
    g, out_png, title,
    width_px = 1000, height_px = 595, dpi = 320,
    write_hires = TRUE, hires_px = 2200, hires_dpi = 320,
    show_legend = FALSE   # <- stays here in case you want it back later
  ) {
    suppressPackageStartupMessages(library(igraph))
    
    .open_px <- function(path, wpx, hpx, res) {
      if (missing(path) || is.null(path) || is.function(path) ||
          !is.character(path) || length(path) != 1L || !nzchar(path))
        stop("out_png must be a single filename string")
      dir.create(dirname(path), recursive = TRUE, showWarnings = FALSE)
      win <- wpx / res; hin <- hpx / res
      if (requireNamespace("ragg", quietly = TRUE)) {
        ragg::agg_png(path, width = win, height = hin, units = "in", res = res, background = "white")
      } else {
        png(path, width = win, height = hin, units = "in", res = res, type = "cairo", bg = "white")
      }
    }
    
    draw_once <- function() {
      if (vcount(g) == 0) { plot.new(); title(main = paste0(title, " — no edges")); return(invisible()) }
      
      layers <- c("ligand","receptor","tf","target")
      cols   <- setNames(c("#42a5f5","#ffb74d","#8d6e63","#66bb6a"), layers)
      vcol   <- cols[V(g)$layer]
      
      # fixed layered layout; occupy the full width
      y <- setNames(c(3, 2, 1, 0), layers)[V(g)$layer]
      x <- ave(seq_along(y), y, FUN = function(ix) seq(0.03, 0.97, length.out = length(ix)))
      lay <- cbind(x, y)
      
      ec <- ifelse(E(g)$kind == "LR", grDevices::adjustcolor("grey20", 0.55),
                   ifelse(E(g)$kind == "RT", grDevices::adjustcolor("grey40", 0.40),
                          grDevices::adjustcolor("grey20", 0.28)))
      
      op <- par(no.readonly = TRUE); on.exit(par(op), add = TRUE)
      # single panel; extra right margin to avoid clipping long labels
      par(mar = c(2.0, 1.4, 3.2, 1.4), xpd = NA, family = "sans")
      
      plot(g, layout = lay,
           rescale = FALSE, xlim = c(-0.02, 1.02), ylim = c(-0.05, 3.15), asp = 0,
           vertex.size = 8,
           vertex.color = vcol,
           vertex.frame.color = "grey25",
           vertex.label = V(g)$label,
           vertex.label.cex = 0.37,
           vertex.label.family = "sans",
           vertex.label.font = 2,
           edge.arrow.size = 0.15,
           edge.width = E(g)$width,
           edge.color = ec,
           main = title)
      
      if (isTRUE(show_legend)) {
        legend("topright", bty = "n",
               legend = c("Ligand","Receptor","TF","Target",
                          "LR score","RT score","TF→Target score"),
               pch    = c(16,16,16,16, NA, NA, NA),
               col    = c(cols, NA, NA, NA),
               pt.cex = .9, cex = .85, lwd = c(NA,NA,NA,NA, 3,3,3))
      }
    }
    
    .open_px(out_png, width_px, height_px, dpi); draw_once(); dev.off()
    
    if (isTRUE(write_hires)) {
      hi <- sub("\\.png$", "_hires.png", out_png)
      .open_px(hi, hires_px, round(hires_px * height_px / width_px), hires_dpi)
      draw_once(); dev.off()
    }
  }
  
  
  
  
  
  
  .prune_edges <- function(df, wcol = "w", q = 0.80, k = Inf) {
    if (!nrow(df)) return(df)
    df[[wcol]] <- as.numeric(df[[wcol]])
    if (any(df[[wcol]] > 0, na.rm = TRUE)) {
      thr <- stats::quantile(df[[wcol]], probs = q, na.rm = TRUE)
      df  <- df[df[[wcol]] >= thr, , drop = FALSE]
    }
    df <- df[order(-df[[wcol]]), , drop = FALSE]
    if (is.finite(k) && nrow(df) > k) df <- df[seq_len(k), , drop = FALSE]
    df
  }
  
  # keep top TFs and Targets (global, not per-node), drop any overhangs
  .keep_top_tf_targets <- function(RT, TT, max_tf = 12, max_targets = 12) {
    RT <- RT[, c("receptor","tf","w"), drop = FALSE]
    TT <- TT[, c("tf","target","w"),   drop = FALSE]
    RT$w <- as.numeric(RT$w); TT$w <- as.numeric(TT$w)
    
    # ---- choose TFs by summed RT weight ----
    if (!nrow(RT)) return(list(RT = RT[0,], TT = TT[0,], keep_receptors = character(0)))
    tf_rank <- aggregate(w ~ tf, data = RT, sum)
    keep_tf <- head(tf_rank[order(-tf_rank$w), "tf"], max_tf)
    RT <- RT[RT$tf %in% keep_tf, , drop = FALSE]
    
    # ---- choose Targets by summed TT weight from the kept TFs ----
    TT <- TT[TT$tf %in% keep_tf, , drop = FALSE]
    if (nrow(TT)) {
      tg_rank <- aggregate(w ~ target, data = TT, sum)
      keep_tg <- head(tg_rank[order(-tg_rank$w), "target"], max_targets)
      TT <- TT[TT$target %in% keep_tg, , drop = FALSE]
    }
    
    # ---- cascade clean-up: no hanging TFs or receptors ----
    if (nrow(TT)) {
      keep_tf2 <- unique(TT$tf)
      RT <- RT[RT$tf %in% keep_tf2, , drop = FALSE]
    } else {
      RT <- RT[0,]
    }
    keep_receptors <- unique(RT$receptor)
    
    list(RT = RT, TT = TT, keep_receptors = keep_receptors)
  }
  
  
  
  # -------------------- main ------------------------------------------------
  if (!dir.exists(out_dir)) dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)
  
  # prune once, then split into clean 2-col layers
  nl_pruned <- prune_net_to_full_paths(
    list(LigRec = scMlnet_results$LigRec,
         RecTF  = scMlnet_results$RecTF,
         TFTar  = scMlnet_results$TFTar),
    sep = "_"
  )
  LR_all <- .pairs_to_df(nl_pruned$LigRec, "ligand",   "receptor")
  RT_all <- .pairs_to_df(nl_pruned$RecTF,  "receptor", "tf")
  TT_all <- .pairs_to_df(nl_pruned$TFTar,  "tf",       "target")
  
  # build per-pathway gating tables
  lig_tbl <- .derive_lig_tbl(lig_rank_all, int_rank_all)
  rec_tbl <- .derive_rec_tbl(rec_rank_all, int_rank_all)
  
  # rank top pathways (fallbacks included)
  pw_tbl <- as.data.frame(pathway_n)
  if (!"pathway" %in% names(pw_tbl)) stop("pathway_n must contain 'pathway'")
  score_col <- if ("score" %in% names(pw_tbl)) "score" else if ("aggregate" %in% names(pw_tbl)) "aggregate" else names(pw_tbl)[sapply(pw_tbl, is.numeric)][1]
  pw_tbl <- pw_tbl[order(-as.numeric(pw_tbl[[score_col]])), , drop = FALSE]
  pw_top <- head(unique(as.character(pw_tbl$pathway)), topK)
  
  # 1) a compact “GNN-like” treemap (summary)
  treemap_path <- file.path(
    out_dir,
    sprintf("DeepGNN_Treemap_%s_to_%s_%s.png",
            sanitize(receiver), sanitize(target), sanitize(norm_label))
  )
  
  plot_gnn_pathway_treemap(
    gnn_df  = pw_tbl[, c("pathway","score")],                  # expects pathway, score
    out_png = treemap_path,
    title   = sprintf("Deep GNN Pathways — %s \u2192 %s", receiver, target),
    palette = c(
      NOTCH="#9AA556", WNT="#E98FC5", ncWNT="#64B1DF", EGF="#66E1E1",
      IGF="#D5AB90", NRG="#8B9D63", BMP="#E88E50", ACTIVIN="#C2E579",
      TGFb="#DBC6DE", NT="#DE8D9F", VEGF="#AAD8E1", HH="#6EE07F",
      FGF="#8580DF", HGF="#D9E83D", PDGF="#62E4B7", GDNF="#E2E2A1"
    ),
    mirror.y = FALSE
  )
  
  # 2) per-pathway explanation PNGs
  explain_paths <- character(0)
  rank_i <- 0L
  for (pw in pw_top) {
    rank_i <- rank_i + 1L
    edges_pw <- .filter_layers_for_pathway_pairs(LR_all, RT_all, TT_all, pw, receptor_ligand)
    edges_pw <- .decorate_weights(edges_pw$LR, edges_pw$RT, edges_pw$TT, pw, lig_tbl, rec_tbl)
    
    # keep all LR/RT/TT weights as-decorated; trim ONLY TFs/Targets to top 12
    sel <- .keep_top_tf_targets(edges_pw$RT, edges_pw$TT, max_tf = 12, max_targets = 12)
    RTp <- sel$RT
    TTp <- sel$TT
    
    # keep ligands/receptors only if they still have connections (prevents overhang)
    LRp <- edges_pw$LR[edges_pw$LR$receptor %in% sel$keep_receptors, , drop = FALSE]
    
    
    .limit_fanout <- function(df, from_col, wcol, k_per_from) {
      if (!nrow(df) || is.infinite(k_per_from)) return(df)
      df[[wcol]] <- as.numeric(df[[wcol]])
      df <- df[order(df[[from_col]], -df[[wcol]]), , drop = FALSE]
      do.call(rbind, lapply(split(df, df[[from_col]]), function(d) head(d, k_per_from)))
    }
    
    
    
    
    g_pw <- .build_graph(LRp, RTp, TTp)
    
    
    out_png <- file.path(out_dir,
                         sprintf("Explain_%d_%s_to_%s_%s_%s.png",
                                 rank_i, sanitize(receiver), sanitize(target),
                                 sanitize(norm_label), sanitize(pw)))
    .plot_explain(
      g_pw,
      out_png = out_png,
      title   = sprintf("%s — rank #%d", pw, rank_i),
      width_px = 860, height_px = 540, dpi = 200,
      write_hires = TRUE
    )
    
    
    
    explain_paths <- c(explain_paths, out_png)
  }
  
  list(paths = list(gnn = treemap_path, explain = explain_paths))
}
# ========================= end deep_learning_run ===========================







# ---------- SIMPLE LEARNING (no external DL pkgs) ----------
# L_R -> set of reachable targets via R->TF->Target
.lr_targets_by_interaction <- function(scMlnet_results) {
  toUC <- function(z) toupper(trimws(as.character(z)))
  LR <- .pairs_to_df(scMlnet_results$LigRec, "ligand","receptor");  LR[] <- lapply(LR, toUC)
  RT <- .pairs_to_df(scMlnet_results$RecTF,  "receptor","tf");      RT[] <- lapply(RT, toUC)
  TT <- .pairs_to_df(scMlnet_results$TFTar,  "tf","target");        TT[] <- lapply(TT, toUC)
  
  tf_by_rec <- split(RT$tf, RT$receptor)
  tg_by_tf  <- split(TT$target, TT$tf)
  
  inter <- paste(LR$ligand, LR$receptor, sep = "_")
  sets <- setNames(vector("list", length(inter)), inter)
  
  for (i in seq_len(nrow(LR))) {
    r  <- LR$receptor[i]
    fs <- tf_by_rec[[r]]
    if (is.null(fs) || !length(fs)) { sets[[i]] <- character(0); next }
    tg <- unique(unlist(tg_by_tf[fs], use.names = FALSE))
    sets[[i]] <- sort(unique(tg))
  }
  sets
}

# Active/Unique metrics for a *given set of interactions*
.lr_active_unique <- function(target_sets, w_tg) {
  totalW <- sum(w_tg, na.rm = TRUE); if (!is.finite(totalW) || totalW == 0) totalW <- 1
  inters <- names(target_sets)
  
  # precompute union of "others" per interaction
  all_union <- unique(unlist(target_sets, use.names = FALSE))
  out <- lapply(inters, function(k) {
    tg <- target_sets[[k]]
    active_frac <- sum(w_tg[intersect(names(w_tg), tg)], na.rm = TRUE) / totalW
    others <- setdiff(all_union, tg)
    uniq   <- setdiff(tg, others)
    uniq_w <- sum(w_tg[intersect(names(w_tg), uniq)], na.rm = TRUE)
    data.frame(interaction = k, active = active_frac, unique = uniq_w, stringsAsFactors = FALSE)
  })
  do.call(rbind, out)
}

.plot_treemap_onelevel <- function(df, out_png, title) {
  dir.create(dirname(out_png), recursive = TRUE, showWarnings = FALSE)
  png(out_png, width = 600, height = 750)
  treemap::treemap(
    df,
    index = "interaction",
    vSize = "score",
    type  = "color",
    vColor = "col",  
    title = title,
    palette = RColorBrewer::brewer.pal(8, "Set3"),
    border.col = "black",
    border.lwds = 2,
    fontsize.labels = 9,
    fontface.labels = 2,
    position.legend = "none"
  )
  dev.off()
}


drop_no_pathway <- function(df, col = "pathway") {
  if (is.null(df) || !nrow(df)) return(df[0, ])
  df[!is.na(df[[col]]) & df[[col]] != "no_pathway", , drop = FALSE]
}

.rescale01 <- function(x) {
  x <- as.numeric(x)
  if (!length(x) || all(!is.finite(x))) return(rep(0, length(x)))
  rng <- range(x, na.rm = TRUE)
  if (!is.finite(rng[1]) || rng[2] - rng[1] < 1e-12) return(rep(0, length(x)))
  (x - rng[1]) / (rng[2] - rng[1])
}

pick_one_per_cluster <- function(pw_meta, w_active = 0.5, w_unique = 0.5) {
  stopifnot(all(c("pathway","cluster","targ_coverage_frac","targ_unique_sum") %in% names(pw_meta)))
  rk <- w_active*.rescale01(pw_meta$targ_coverage_frac) +
    w_unique*.rescale01(pw_meta$targ_unique_sum)
  pw_meta$rk <- rk
  best <- lapply(split(seq_len(nrow(pw_meta)), pw_meta$cluster), function(ii) ii[which.max(pw_meta$rk[ii])])
  pw_meta[unlist(best), , drop = FALSE]
}


# --- Pareto front on (Active, Unique) ---
.pareto_front <- function(df, x = "targ_coverage_frac", y = "targ_unique_sum") {
  if (!nrow(df)) return(df)
  d <- df[, c(x, y)]
  # treat NA as worst possible (so they won’t dominate anything)
  d[, x][!is.finite(d[, x])] <- -Inf
  d[, y][!is.finite(d[, y])] <- -Inf
  keep <- rep(TRUE, nrow(d))
  for (i in seq_len(nrow(d))) {
    if (!keep[i]) next
    dom <- d[, x] >= d[i, x] & d[, y] >= d[i, y] &
      (d[, x] >  d[i, x] | d[, y] >  d[i, y])
    dom[i] <- FALSE
    if (any(dom, na.rm = TRUE)) keep[i] <- FALSE
  }
  df[keep, , drop = FALSE]
}


# --- tidy treemap wrapper with your aesthetic locked ---
.plot_treemap <- function(tm_df, out_png, title, mirror.y = FALSE) {
  dir.create(dirname(out_png), recursive = TRUE, showWarnings = FALSE)
  png(out_png, width = 600, height = 750)
  treemap::treemap(
    tm_df,
    index = c("cluster", "pathway"),
    vSize = "score",
    type = "color",
    vColor = "col",                  # already hex colors in tm_df$col
    title = title,
    border.col      = c("black","black"),
    border.lwds     = c(3,1),
    mirror.y        = mirror.y,
    fontsize.labels = c(0, 9),       # suppress huge cluster digits; pathway = 9pt
    fontcolor.labels= c("black","black"),
    fontface.labels = c(2,1),
    align.labels    = list(c("left","top"), c("center","center")),
    bg.labels       = 0,
    inflate.labels  = TRUE,
    position.legend = "none",
    lowerbound.cex.labels = 0.25
  )
  dev.off()
}

.robust_scale <- function(M) {
  M <- as.matrix(M)
  mu  <- colMeans(M)
  sdv <- apply(M, 2, sd)
  sdv[!is.finite(sdv) | sdv < 1e-8] <- 1
  sweep(sweep(M, 2, mu, "-"), 2, sdv, "/")
}


.plot_treemap_2 <- function(tm_df, out_png, title, mirror.y = FALSE, pal_by_cluster = NULL) {
  dir.create(dirname(out_png), recursive = TRUE, showWarnings = FALSE)
  tm_df <- as.data.frame(tm_df, stringsAsFactors = FALSE)
  
  need <- c("pathway","cluster","score")
  miss <- setdiff(need, names(tm_df))
  if (length(miss)) return(.placeholder_png(out_png, paste("Missing columns:", paste(miss, collapse=", "))))
  
  tm_df$pathway <- trimws(as.character(tm_df$pathway))
  tm_df$cluster <- trimws(as.character(tm_df$cluster))
  tm_df$score   <- suppressWarnings(as.numeric(tm_df$score))
  
  tm_df <- tm_df[is.finite(tm_df$score) & tm_df$score > 0, , drop = FALSE]
  tm_df <- tm_df[nzchar(tm_df$pathway), , drop = FALSE]
  tm_df$cluster[!nzchar(tm_df$cluster) | is.na(tm_df$cluster)] <- "Unclustered"
  if (!nrow(tm_df)) return(.placeholder_png(out_png, "No pathways with non-zero size"))
  
  clusters <- sort(unique(tm_df$cluster))
  if (is.null(pal_by_cluster)) {
    base12 <- RColorBrewer::brewer.pal(12, "Set3")
    pal_by_cluster <- setNames(colorRampPalette(base12)(length(clusters)), clusters)
  } else {
    pal_by_cluster <- pal_by_cluster[clusters]
    miss <- is.na(pal_by_cluster)
    if (any(miss)) {
      base12 <- RColorBrewer::brewer.pal(12, "Set3")
      pal_by_cluster[miss] <- colorRampPalette(base12)(sum(miss))
      names(pal_by_cluster) <- clusters
    }
  }
  
  tm_df$col <- unname(pal_by_cluster[tm_df$cluster])
  tm_df$col[is.na(tm_df$col) | !nzchar(tm_df$col)] <- "#CCCCCC"
  
  png(out_png, width = 600, height = 750)
  treemap::treemap(
    tm_df,
    index   = c("cluster","pathway"),
    vSize   = "score",
    type    = "color",
    vColor  = "col",
    title   = title,
    border.col      = c("black","black"),
    border.lwds     = c(3,1),
    mirror.y        = isTRUE(mirror.y),
    fontsize.labels = c(0, 9),
    fontcolor.labels= c("black","black"),
    fontface.labels = c(2,1),
    align.labels    = list(c("left","top"), c("center","center")),
    bg.labels       = 0,
    inflate.labels  = TRUE,
    position.legend = "none",
    lowerbound.cex.labels = 0.25
  )
  dev.off()
}





# --- tiny PNG with a text message (for “no LR data” instead of NULL) ---
.placeholder_png <- function(path, msg = "No data") {
  dir.create(dirname(path), recursive = TRUE, showWarnings = FALSE)
  png(path, width = 600, height = 750)
  op <- par(no.readonly = TRUE); on.exit(par(op), add = TRUE)
  par(mar = c(0,0,2,0))
  plot.new(); title(main = msg, cex.main = 1.2)
  dev.off()
  path
}


# Safe data.frame that auto-repeats scalars and fails with a readable message
safe_df <- function(..., .n = NULL) {
  cols <- list(...)
  nms  <- names(cols)
  lens <- vapply(cols, length, integer(1))
  if (length(lens) == 0) return(data.frame())
  target <- if (is.null(.n)) max(lens) else .n
  fix_len <- function(x, nm) {
    lx <- length(x)
    if (lx == target) return(x)
    if (lx == 1L)     return(rep(x, length.out = target))
    stop(sprintf("safe_df: column '%s' has length %d incompatible with target %d",
                 nm, lx, target), call. = FALSE)
  }
  cols <- Map(fix_len, cols, nms)
  data.frame(cols, stringsAsFactors = FALSE)
}

# Quick sanity checker that prints a compact summary of an object/data.frame
dbg <- function(label, x) {
  cat(sprintf("[DBG] %s: class=%s  nrow=%s  ncol=%s  length=%s\n",
              label, paste(class(x), collapse = "/"),
              ifelse(is.data.frame(x) || is.matrix(x), NROW(x), NA),
              ifelse(is.data.frame(x) || is.matrix(x), NCOL(x), NA),
              length(x)))
  invisible(x)
}

# Assert presence of columns; give a precise error with a prefix
assert_cols <- function(df, cols, where = "unknown") {
  miss <- setdiff(cols, names(df))
  if (length(miss))
    stop(sprintf("Missing columns in %s: %s", where, paste(miss, collapse = ", ")),
         call. = FALSE)
  invisible(TRUE)
}

`%||%` <- function(a,b) if (is.null(a) || (length(a)==1 && is.na(a))) b else a

.pick_col <- function(df, candidates) {
  hit <- intersect(candidates, names(df))
  if (length(hit)) hit[1] else NA_character_
}

# ---------- Feature builders ----------
# Replace the body of .pathway_features() for the TF-centered modes with this idea:

.pathway_features <- function(scMlnet_results, deg = NULL, tf_weight_vec = NULL,
                              mode = c("tf_centered","tftarget_centered",
                                       "tf_centered_unweighted","tftarget_centered_unweighted"),
                              rl_db = NULL) {
  mode <- match.arg(mode)
  
  
  # resolve receptor_ligand if not supplied
  if (is.null(rl_db)) rl_db <- get("receptor_ligand", inherits = TRUE)
  
  # 1) Precompute reachability from ligand-pathway annotation → (TFs, Targets)
  R <- build_reachability(scMlnet_results, rl_db)
  PWs <- names(R$tfs_by_pathway)
  
  if (mode %in% c("tf_centered","tf_centered_unweighted")) {
    TFs <- sort(unique(unlist(R$tfs_by_pathway, use.names = FALSE)))
    M <- matrix(0, nrow = length(PWs), ncol = length(TFs),
                dimnames = list(PWs, TFs))
    
    # default = unweighted
    wTF <- setNames(rep(1, length(TFs)), TFs)
    
    if (mode == "tf_centered" && !is.null(tf_weight_vec)) {
      # names already uppercased in the early normalizer
      w <- as.numeric(tf_weight_vec[TFs])
      match_rate <- mean(!is.na(w))
      if (!is.finite(match_rate)) match_rate <- 0
      # fallback to 1 for missing TFs
      w[is.na(w)] <- 1
      
      # if nearly flat, treat as unweighted to avoid fake differences
      if (sd(w) < 1e-9) {
        wTF <- setNames(rep(1, length(TFs)), TFs)
      } else {
        wTF <- w / mean(w, na.rm = TRUE)
      }
      
      # optional: message when names barely matched
      if (match_rate < 0.5)
        message(sprintf("tf_centered: only %.0f%% of TFs had a weight (check name casing).",
                        100*match_rate))
    }
    
    for (p in PWs) {
      tfs <- R$tfs_by_pathway[[p]]
      if (length(tfs)) M[p, tfs] <- wTF[tfs]
    }
    return(M)
  }
  
  
  TGs <- sort(unique(unlist(R$targets_by_pathway, use.names = FALSE)))
  M <- matrix(0, nrow = length(PWs), ncol = length(TGs),
              dimnames = list(PWs, TGs))
  wTG <- rep(1, length(TGs)); names(wTG) <- TGs
  if (mode == "tftarget_centered" && !is.null(deg)) {
    v <- abs(suppressWarnings(as.numeric(deg$foldChange)))
    names(v) <- toupper(as.character(deg$Gene))
    wTG <- v[TGs]; wTG[is.na(wTG)] <- 0
  }
  for (p in PWs) {
    tg <- R$targets_by_pathway[[p]]
    if (length(tg)) M[p, tg] <- wTG[tg]
  }
  M
}




# ---------- Clustering + 2D projection ----------
.cluster_from_features <- function(M, proj = c("mds","pca")) {
  proj <- match.arg(proj)
  if (!nrow(M) || !ncol(M)) return(list(hclust=NULL, coords=matrix(0,0,2)))
  # L2 row-normalize -> cosine space
  rs <- sqrt(rowSums(M^2)); rs[rs < 1e-12] <- 1e-12
  Mz <- M / rs
  
  # distance + clustering
  S <- tcrossprod(Mz)         # cosine similarity
  S[!is.finite(S)] <- 0
  D <- as.dist(1 - pmin(pmax(S, -1), 1))
  hc <- hclust(D, method = "average")
  
  # 2D coordinates
  if (proj == "mds") {
    coords <- cmdscale(D, k = 2)
  } else {
    pc <- prcomp(Mz, center = TRUE, scale. = FALSE)
    coords <- pc$x[, 1:2, drop = FALSE]
  }
  list(hclust = hc, coords = coords)
}

# ---- Receptor-pruned treemap (from rec_rank_all) ----
make_rec_pruned <- function(rec_rank_all, out_dir, receiver, target, norm_label, top_k = 8L) {
  df <- as.data.frame(rec_rank_all %||% data.frame())
  if (!nrow(df)) return(NULL)
  
  # standardize columns
  if (!"receptor" %in% names(df)) {
    rn <- rownames(df); if (!is.null(rn) && any(nzchar(rn))) df$receptor <- rn
  }
  if (!"pathway" %in% names(df)) df$pathway <- as.character(df$pathway %||% df$name)
  sc <- df$score %||% df$aggregate
  df$score <- suppressWarnings(as.numeric(sc))
  df <- df[!is.na(df$receptor) & !is.na(df$score), c("receptor","pathway","score")]
  
  # pooled score per receptor + dominant pathway
  rec_sum  <- aggregate(score ~ receptor, data = df, sum)
  rec_path <- aggregate(score ~ receptor + pathway, data = df, sum)
  rec_path <- rec_path[order(rec_path$receptor, -rec_path$score), ]
  rec_path <- rec_path[!duplicated(rec_path$receptor), c("receptor","pathway")]
  
  pf_rec <- merge(rec_sum, rec_path, by = "receptor", all.x = TRUE)
  pf_rec <- pf_rec[order(-pf_rec$score), ]
  pf_rec_top <- utils::head(pf_rec, top_k)
  
  if (!nrow(pf_rec_top)) return(NULL)
  
  rec_png <- file.path(out_dir, sprintf("Simple_rec_pruned_treemap_%s_to_%s_%s.png",
                                        receiver, target, norm_label))
  pal <- colorRampPalette(RColorBrewer::brewer.pal(8, "Set3"))
  png(rec_png, width = 720, height = 900)
  treemap::treemap(
    pf_rec_top,
    index   = "receptor",
    vSize   = "score",
    vColor  = "pathway",
    type    = "categorical",
    palette = pal(max(3, length(unique(pf_rec_top$pathway)))),
    title   = "Top receptors (size = pooled receptor score; color = dominant pathway)",
    fontsize.labels = 14,
    bg.labels       = "white",
    border.col      = "black"
  )
  dev.off()
  rec_png
}

.canon_interactions <- function(df) {
  df <- as.data.frame(df, stringsAsFactors = FALSE)
  
  # Locate likely columns
  sc <- .pick_col(df, c("score","aggregate","weight","value"))
  ic <- .pick_col(df, c("interaction","pair","lig_rec","ligand_receptor","Ligand_Receptor"))
  lc <- .pick_col(df, c("ligand","Ligand","Ligand.ApprovedSymbol","from"))
  rc <- .pick_col(df, c("receptor","Receptor","Receptor.ApprovedSymbol","to"))
  
  if (is.na(sc)) stop("int_rank_all has no score/aggregate column")
  
  # Derive missing pieces
  if (is.na(ic) && !is.na(lc) && !is.na(rc)) {
    df$interaction <- paste(df[[lc]], df[[rc]], sep = "_")
    ic <- "interaction"
  }
  if (!is.na(ic) && (is.na(lc) || is.na(rc))) {
    sp <- .split_pairs(df[[ic]], nm = c("ligand","receptor"))
    if (is.na(lc)) df$ligand   <- sp$ligand
    if (is.na(rc)) df$receptor <- sp$receptor
    lc <- "ligand"; rc <- "receptor"
  }
  
  # Standardize names
  names(df)[match(c(lc, rc, sc), names(df))] <- c("ligand","receptor","score")
  
  # Pathway if missing → infer from receptor_ligand by ligand (best-effort)
  if (!"pathway" %in% names(df)) {
    if (exists("receptor_ligand", inherits = TRUE)) {
      mp <- unique(receptor_ligand[, c("Ligand.ApprovedSymbol","pathway_name")])
      names(mp) <- c("ligand","pathway")
      df <- merge(df, mp, by = "ligand", all.x = TRUE)
    } else {
      df$pathway <- NA_character_
    }
  }
  
  # Hard requirements
  need <- c("interaction","ligand","receptor","score")
  miss <- setdiff(need, names(df))
  if (length(miss)) stop(sprintf("int_rank_all is missing required columns after canonization: %s",
                                 paste(miss, collapse=", ")))
  df
}

# ---- 2D embedding helper (robust, tiny, drop-in) ----
embed2D <- function(D) {
  # accept dist OR symmetric matrix
  if (inherits(D, "dist")) {
    dd   <- D
    labs <- attr(dd, "Labels"); n <- attr(dd, "Size")
  } else {
    M <- as.matrix(D)
    M[!is.finite(M)] <- 0
    M <- (M + t(M))/2
    diag(M) <- 0
    dd   <- stats::as.dist(pmax(0, M))
    labs <- rownames(M); n <- nrow(M)
  }
  if (is.null(n) || n == 0) {
    return(matrix(numeric(0), ncol = 2,
                  dimnames = list(character(0), c("x","y"))))
  }
  if (is.null(labs) || length(labs) != n) labs <- paste0("item", seq_len(n))
  if (n == 1L) {
    out <- matrix(c(0, 0), nrow = 1, dimnames = list(labs, c("x","y")))
    return(out)
  }
  if (n == 2L) {
    out <- matrix(c(-0.5, 0, 0.5, 0), ncol = 2, byrow = TRUE,
                  dimnames = list(labs, c("x","y")))
    return(out)
  }
  # prefer classical MDS; fall back to eigen trick if needed
  emb <- try(stats::cmdscale(dd, k = 2), silent = TRUE)
  if (inherits(emb, "try-error")) {
    DM <- as.matrix(dd)
    J  <- diag(n) - matrix(1/n, n, n)
    B  <- -0.5 * J %*% (DM^2) %*% J
    ev <- eigen(B, symmetric = TRUE)
    vals <- pmax(ev$values[1:2], 0)
    V    <- ev$vectors[, 1:2, drop = FALSE]
    emb  <- V %*% diag(sqrt(vals), 2)
  }
  colnames(emb) <- c("x","y"); rownames(emb) <- labs
  emb
}

# TF-IDF + cleanup for pathway×feature matrices
.prep_features <- function(M, do_tfidf = TRUE, drop_ubiq = 0.70, zscore_cols = TRUE) {
  M <- as.matrix(M); M[is.na(M)] <- 0
  if (do_tfidf) {
    df  <- colSums(M > 0)
    idf <- log(nrow(M) / pmax(1, df))
    M   <- sweep(M, 2, idf, `*`)
  }
  if (!is.null(drop_ubiq)) {
    frac <- colSums(M > 0) / nrow(M)
    keep <- frac <= drop_ubiq
    if (any(!keep)) M <- M[, keep, drop = FALSE]
  }
  if (zscore_cols && ncol(M) > 1) {
    sds <- apply(M, 2, sd); sds[!is.finite(sds) | sds < 1e-8] <- 1
    M <- scale(M, center = TRUE, scale = sds)
    M[is.na(M)] <- 0
  }
  M
}

.has_arg <- function(fn, arg) arg %in% names(formals(fn))

# drop-in: better defaults for small-N, uwot >= 0.23
# drop-in
# replace your .pick_umap_args() with this
.pick_umap_args <- function(M, seed = 1, user = list()) {
  n <- nrow(M)
  
  # very small N → just use MDS
  if (n < 8) return(list(.method = "mds"))
  
  # in .pick_umap_args()
  nn <- min(n - 1, max(6L, ceiling(0.7 * n)))   # a bit more global context
  args <- modifyList(list(
    X                  = M,
    n_components       = 2,
    metric             = "cosine",
    nn_method          = "hnsw",
    n_neighbors        = nn,
    min_dist           = 0.65,   # <- bigger = tighter clumps, less spread
    local_connectivity = 2,
    set_op_mix_ratio   = 0.98,   # <- closer to union = denser graph
    repulsion_strength = 0.7,    # <- softer repulsion
    init               = "spectral",
    n_sgd_threads      = 1L,
    verbose            = FALSE,
    ret_model          = FALSE,
    scale              = FALSE,
    seed               = seed
  ), user)
  
  
  args
}


.compact_coords <- function(X, target_sd = 0.45, clip = 3) {
  X <- scale(X, center = TRUE, scale = apply(X, 2, function(v) stats::mad(v) %||% sd(v)))
  X[X >  clip] <-  clip
  X[X < -clip] <- -clip
  X <- X / max(1e-8, sd(as.vector(X))) * target_sd
  X
}

auto_hdbscan_on_coords <- function(coords2d, minPts_grid = NULL, prob_min = 0.55, k_target = 4) {
  stopifnot(requireNamespace("dbscan", quietly = TRUE))
  X <- .robust_scale(as.matrix(coords2d))
  if (is.null(minPts_grid)) minPts_grid <- 2:max(6, floor(nrow(X)/3))
  
  best <- NULL; best_score <- -Inf
  for (m in minPts_grid) {
    hd <- dbscan::hdbscan(X, minPts = m)
    k  <- length(unique(hd$cluster[hd$cluster > 0]))
    if (k == 0) next
    mean_prob <- mean(hd$membership_prob[hd$cluster > 0])
    noise_fr  <- mean(hd$cluster == 0)
    score <- 1.0*mean_prob - 0.55*noise_fr - 0.05*(k - k_target)^2   # <-- prefer k≈4
    if (is.finite(score) && score > best_score) { best <- hd; best_score <- score }
  }
  if (is.null(best)) return(list(cluster = factor(rep("N", nrow(X))), model = NULL))
  lab <- ifelse(best$cluster > 0 & best$membership_prob >= prob_min, paste0("C", best$cluster), "N")
  
  # renumber to C1..Ck with no gaps
  pos <- sort(unique(lab[lab != "N"]))
  if (length(pos)) {
    remap <- setNames(paste0("C", seq_along(pos)), pos)
    lab[lab != "N"] <- remap[lab[lab != "N"]]
  }
  lab[!nzchar(lab) | is.na(lab)] <- "N"
  list(cluster = factor(lab, levels = c("N", sort(unique(lab[lab!="N"])))), model = best)
}





# drop-in replacement (robust)
.mds_cosine <- function(M) {
  M <- as.matrix(M); storage.mode(M) <- "double"; M[!is.finite(M)] <- 0
  n <- nrow(M); labs <- rownames(M); if (is.null(labs)) labs <- paste0("p", seq_len(n))
  
  if (n == 0L) return(matrix(numeric(0), ncol=2, dimnames=list(character(0), c("x","y"))))
  if (n == 1L) return(matrix(c(0,0), nrow=1, dimnames=list(labs, c("x","y"))))
  if (n == 2L) return(matrix(c(-0.5,0, 0.5,0), nrow=2, byrow=TRUE, dimnames=list(labs, c("x","y"))))
  
  # L2 row-normalize → cosine similarity
  rn <- sqrt(rowSums(M^2)); rn[rn < 1e-12] <- 1
  Xn <- M / rn
  S  <- tcrossprod(Xn)
  S  <- (S + t(S)) / 2
  diag(S) <- 1
  S[!is.finite(S)] <- 0
  
  # explicit square matrix for as.dist (prevents the “non-square” complaint)
  DM <- 1 - pmin(pmax(S, -1), 1)
  diag(DM) <- 0
  dd <- stats::as.dist(DM)
  
  coords <- stats::cmdscale(dd, k = 2)
  if (is.null(dim(coords))) coords <- cbind(coords, 0)
  colnames(coords) <- c("x","y"); rownames(coords) <- labs
  coords
}

auto_hdbscan_on_umap <- function(coords2d, minPts_grid = NULL, prob_min = 0.60) {
  stopifnot(requireNamespace("dbscan", quietly = TRUE))
  X <- scale(as.matrix(coords2d))
  
  if (is.null(minPts_grid))
    minPts_grid <- 2:max(6, floor(nrow(X)/3))
  
  best <- NULL; best_score <- -Inf
  for (m in minPts_grid) {
    hd <- dbscan::hdbscan(X, minPts = m)
    k  <- length(unique(hd$cluster[hd$cluster > 0]))
    if (k == 0) next
    mean_prob <- mean(hd$membership_prob[hd$cluster > 0])
    noise_fr  <- mean(hd$cluster == 0)
    # prefer high confidence, few noise points, and 2–4 clusters
    score <- 1.0*mean_prob - 0.5*noise_fr - 0.1*abs(k - 3)
    if (is.finite(score) && score > best_score) { best <- hd; best_score <- score }
  }
  if (is.null(best)) return(list(cluster = factor(rep("N", nrow(X))), model = NULL))
  
  lab <- ifelse(best$cluster > 0 & best$membership_prob >= prob_min,
                paste0("C", best$cluster),
                "N")
  list(cluster = factor(lab, levels = c("N", sort(unique(lab[lab != "N"])))),
       model   = best)
}


# UMAP (if available) -> else MDS
# ---- replace your .embed2D_plus with this version ----
.embed2D_plus <- function(M, method = c("auto","mds","umap","pca"),
                          umap_args = list(), seed = 1) {
  method <- match.arg(method)
  M <- as.matrix(M); storage.mode(M) <- "double"; M[!is.finite(M)] <- 0
  
  # Auto: prefer MDS to avoid uwot incompatibilities
  if (identical(method, "auto")) method <- "mds"
  
  if (method == "umap") {
    stopifnot(requireNamespace("uwot", quietly = TRUE))
    n <- nrow(M); if (n < 8) return(.mds_cosine(M))
    
    mds_init <- .mds_cosine(M)
    args <- modifyList(list(
      X                  = M,
      n_components       = 2,
      metric             = "cosine",
      nn_method          = "hnsw",
      n_neighbors        = min(n-1, max(6L, round(0.8*n))),
      min_dist           = 0.55,
      set_op_mix_ratio   = 0.98,
      local_connectivity = 2,
      repulsion_strength = 0.7,
      # densmap / dens_lambda are NOT universally available → don’t include
      init               = mds_init,
      n_sgd_threads      = 1L,
      verbose            = FALSE,
      ret_model          = FALSE,
      scale              = FALSE,
      seed               = seed
    ), umap_args)
    
    # Drop any args not supported by this uwot version
    allowed <- names(formals(uwot::umap))
    args <- args[intersect(names(args), allowed)]
    
    emb <- do.call(uwot::umap, args)
    emb <- .compact_coords(emb, target_sd = 0.45)
    rownames(emb) <- rownames(M); colnames(emb) <- c("x","y")
    attr(emb, "method_used") <- "UMAP/hnsw (MDS-anchor)"
    return(emb)
  }
  
  if (method == "pca") {
    pc <- prcomp(M, center = TRUE, scale. = FALSE)
    coords <- pc$x[, 1:2, drop = FALSE]; colnames(coords) <- c("x","y")
    attr(coords, "method_used") <- "PCA"; return(coords)
  }
  
  coords <- .mds_cosine(M)
  colnames(coords) <- c("x","y"); rownames(coords) <- rownames(M)
  attr(coords, "method_used") <- "MDS"
  coords
}



# Choose K by average silhouette on a hierarchical tree (cosine distance)
# Choose K by average silhouette (robust to 2D coords or high-D features)
# ---- replace your .choose_k_by_sil with this hardened version ----
# robust K chooser: works with high-D features OR 2-D coords
.choose_k_by_sil <- function(M, kmax = 8) {
  if (!requireNamespace("cluster", quietly = TRUE))
    return(max(2, min(6, round(sqrt(nrow(M))))))
  
  M <- as.matrix(M)
  storage.mode(M) <- "double"
  M[!is.finite(M)] <- 0
  n <- nrow(M)
  if (n < 3) return(2)
  
  # 1) Build a proper distance object
  if (ncol(M) == 2 && all(tolower(colnames(M)) %in% c("x","y"))) {
    # You passed 2-D coordinates → use Euclidean distances
    D <- stats::dist(M)
  } else {
    # You passed features → use cosine distance
    rn <- sqrt(rowSums(M^2)); rn[rn < 1e-12] <- 1
    Xn <- M / rn
    S  <- tcrossprod(Xn)
    S  <- (S + t(S))/2
    S[!is.finite(S)] <- 0
    diag(S) <- 1
    D  <- stats::as.dist(pmax(0, 1 - S))
  }
  
  # safety check so hclust never sees a bad length
  if (length(D) != n*(n-1)/2)
    stop(sprintf(".choose_k_by_sil: bad dist length %d for n=%d", length(D), n), call. = FALSE)
  
  hc <- stats::hclust(D, method = "average")
  ks <- 2:min(kmax, max(2, n - 1))
  best <- ks[1]; bestS <- -Inf
  for (k in ks) {
    cl <- stats::cutree(hc, k)
    si <- cluster::silhouette(cl, D)
    m  <- mean(si[,3])
    if (is.finite(m) && m > bestS) { bestS <- m; best <- k }
  }
  best
}






simple_learning_run <- function(
    lig_rank_all, pathway_n,
    int_rank_all = NULL,
    rec_rank_all = NULL,
    scMlnet_results = NULL,
    receptor_ligand = NULL,
    deg = NULL, tf_weight_vec = NULL,
    out_dir, receiver, target, norm_label = "db",
    mode = c("ligand_weighted","ligand_tf_weighted","interaction_weighted","ligand_set",
             "tf_centered","tftarget_centered","tf_centered_unweighted","tftarget_centered_unweighted"),
    mirror.y = FALSE,
    n_clusters = NULL,
    rank_method = c("pareto","weighted"),
    k_per_cluster = 1,
    w_active = 0.5, w_unique = 0.5,
    embed_method = c("auto","umap","mds","pca"), 
    umap_args = list(n_neighbors = NULL, min_dist = NULL),
    cluster_method = c("kmeans","hdbscan","dbscan"),
    dbscan_minPts = NULL, dbscan_eps = NULL,
    min_clusters       = 4,
    hdbscan_prob_min   = 0.50,             # was 0.60 → too strict
    hdbscan_snap_noise = TRUE,
    cluster_fallback   = c("kmeans","none"),
    seed = 1
) {
  embed_method <- match.arg(embed_method)
  # normalize TF weight vector names once (upper-case, trim)
  if (!is.null(tf_weight_vec)) {
    if (is.null(names(tf_weight_vec)))
      stop("tf_weight_vec must be a *named* numeric vector of TF weights.")
    nm <- toupper(trimws(names(tf_weight_vec)))
    # collapse duplicate TFs by mean (robust if names repeat)
    tf_weight_vec <- tapply(as.numeric(tf_weight_vec), nm, mean)
  }
  
  rank_method <- match.arg(rank_method)
  
  # normalize/alias BEFORE match.arg()
  if (is.character(mode) && length(mode) == 1L) {
    m <- tolower(mode)
    if (m %in% c("tftarget_weighted","tf_target_weighted","tf-target-weighted")) {
      warning("simple_learning_run: 'tftarget_weighted' is deprecated; using 'tftarget_centered'.", call. = FALSE)
      mode <- "tftarget_centered"
    }
    
  }
  
  mode <- match.arg(mode, c("ligand_weighted","ligand_tf_weighted","interaction_weighted","ligand_set",
                            "tf_centered","tftarget_centered","tf_centered_unweighted","tftarget_centered_unweighted"))
  if (identical(mode, "tftarget_weighted")) {
    warning("simple_learning_run: 'tftarget_weighted' not supported in SIMPLE; falling back to 'ligand_weighted'.")
    mode <- "ligand_weighted"
  }
  dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)
  
  suppressPackageStartupMessages({
    library(ggplot2); library(ggrepel); library(scales); library(reshape2)
  })
  
  
  paths <- list()
  emb_df <- NULL 
  
  
  
  # harden input tables
  numify <- function(x) suppressWarnings(as.numeric(x))
  
  need <- function(df, cols, where) {
    miss <- setdiff(cols, names(df))
    if (length(miss)) stop(sprintf("%s is missing: %s", where, paste(miss, collapse=", ")), call. = FALSE)
  }
  
  # lig_rank_all must have ligand, pathway, numeric score
  lig_rank_all <- as.data.frame(lig_rank_all, stringsAsFactors = FALSE)
  need(lig_rank_all, c("ligand","pathway"), "lig_rank_all")
  if (!"score" %in% names(lig_rank_all)) {
    alt <- lig_rank_all$aggregate
    if (is.null(alt)) alt <- lig_rank_all$value
    if (is.null(alt)) stop("lig_rank_all lacks 'score' (or 'aggregate'/'value').", call. = FALSE)
    lig_rank_all$score <- numify(alt)
  } else {
    lig_rank_all$score <- numify(lig_rank_all$score)
  }
  
  # int_rank_all must have interaction, ligand, receptor, pathway, numeric score
  if (!is.null(int_rank_all)) {
    int_rank_all <- as.data.frame(int_rank_all, stringsAsFactors = FALSE)
    need(int_rank_all, c("interaction","ligand","receptor","pathway"), "int_rank_all")
    if (!"score" %in% names(int_rank_all)) {
      alt <- int_rank_all$aggregate
      if (is.null(alt)) alt <- int_rank_all$value
      if (is.null(alt)) stop("int_rank_all lacks 'score' (or 'aggregate'/'value').", call. = FALSE)
      int_rank_all$score <- numify(alt)
    } else {
      int_rank_all$score <- numify(int_rank_all$score)
    }
  }
  
  # rec_rank_all needs receptor, pathway, numeric score when used later
  if (!is.null(rec_rank_all)) {
    rec_rank_all <- as.data.frame(rec_rank_all, stringsAsFactors = FALSE)
    need(rec_rank_all, c("receptor","pathway"), "rec_rank_all")
    if (!"score" %in% names(rec_rank_all)) {
      alt <- rec_rank_all$aggregate
      if (is.null(alt)) alt <- rec_rank_all$value
      if (is.null(alt)) stop("rec_rank_all lacks 'score' (or 'aggregate'/'value').", call. = FALSE)
      rec_rank_all$score <- numify(alt)
    } else {
      rec_rank_all$score <- numify(rec_rank_all$score)
    }
  }
  
  # ---------- sanitize & keep pathways ----------
  drop_no_pathway <- function(df, col = "pathway") {
    if (is.null(df) || !nrow(df)) return(df[0, ])
    df[!is.na(df[[col]]) & df[[col]] != "no_pathway", , drop = FALSE]
  }
  lig_rank_all <- drop_no_pathway(as.data.frame(lig_rank_all), "pathway")
  pathway_n    <- drop_no_pathway(as.data.frame(pathway_n),    "pathway")
  if (!nrow(lig_rank_all) || !nrow(pathway_n)) stop("Empty ligand/pathway tables.")
  
  # ---------- build distance for SIMPLE scatter ----------
  mk_cosine_dist <- function(M) {
    M <- as.matrix(M); if (is.null(dim(M))) M <- matrix(M, nrow = 1L)
    M[is.na(M)] <- 0
    rn <- sqrt(rowSums(M^2)); rn[rn == 0] <- 1
    Xn <- M / rn
    S <- Xn %*% t(Xn); S <- (S + t(S))/2; diag(S) <- 1
    D <- 1 - S; D[!is.finite(D)] <- 0; D[D < 0] <- 0
    stats::as.dist(D)
  }
  
  if (mode == "ligand_weighted") {
    # classic: pathway × ligand from lig_rank_all, trim to pathway_n, drop empties
    X <- reshape2::dcast(
      lig_rank_all[, c("pathway","ligand","score")],
      pathway ~ ligand, value.var = "score",
      fun.aggregate = function(z) sum(as.numeric(z)), fill = 0
    )
    rownames(X) <- X$pathway; X$pathway <- NULL
    
    
    
    # keep only pathways we are analyzing/plotting
    allowed <- unique(as.character(pathway_n$pathway))
    if (length(allowed)) X <- X[intersect(rownames(X), allowed), , drop = FALSE]
    
    # critical: drop all-zero rows (avoids MDS collapse to a single point)
    X <- X[rowSums(X) > 0, , drop = FALSE]
    stopifnot(nrow(X) >= 2)
    
    # embed here so we don't rely on embed2D()
    X <- .prep_features(X, do_tfidf = TRUE, drop_ubiq = NULL, zscore_cols = TRUE)
    coords <- .embed2D_plus(X, method = embed_method, seed = seed)   # UMAP if available, else MDS
    
    emb_df <- data.frame(pathway = rownames(coords),
                         x = coords[,1], y = coords[,2],
                         stringsAsFactors = FALSE)
  } else if (mode == "interaction_weighted") {
    stopifnot(!is.null(int_rank_all))
    Y <- dcast(int_rank_all[, c("pathway","interaction","score")],
               pathway ~ interaction, value.var = "score",
               fun.aggregate = function(z) sum(as.numeric(z)), fill = 0)
    rownames(Y) <- Y$pathway; Y$pathway <- NULL
    D <- mk_cosine_dist(Y); pths_use <- rownames(Y)
    
    Y <- .prep_features(Y, do_tfidf = TRUE, drop_ubiq = NULL, zscore_cols = TRUE)
    coords <- .embed2D_plus(Y, method = embed_method, seed = seed)
    
  } else if (mode == "ligand_set") {
    split_ligs <- split(as.character(lig_rank_all$ligand),
                        as.character(lig_rank_all$pathway))
    # keep only pathways we care about
    allowed <- unique(as.character(pathway_n$pathway))
    split_ligs <- split_ligs[names(split_ligs) %in% allowed]
    
    # drop empty sets
    split_ligs <- split_ligs[sapply(split_ligs, function(v) length(unique(v)) > 0)]
    stopifnot(length(split_ligs) >= 2)
    
    pths <- names(split_ligs)
    jacc <- function(a, b) {
      A <- unique(a); B <- unique(b); u <- length(unique(c(A,B))); if (!u) return(0)
      length(intersect(A,B)) / u
    }
    S <- outer(pths, pths, Vectorize(function(i,j) jacc(split_ligs[[i]], split_ligs[[j]])))
    diag(S) <- 1
    D <- as.dist(1 - S)
    
    coords <- cmdscale(D, k = 2)
    emb_df <- data.frame(pathway = rownames(coords),
                         x = coords[,1], y = coords[,2],
                         stringsAsFactors = FALSE)
    
  } else if (mode == "ligand_tf_weighted") {
    # Walk L->R->TF once; build pathway x TF matrix weighted by tf_weight_vec (AUCell) or 1
    stopifnot(!is.null(scMlnet_results), !is.null(receptor_ligand))
    Rch <- build_reachability(scMlnet_results, receptor_ligand)
    sets <- Rch$tfs_by_pathway
    all_tfs <- sort(unique(unlist(sets, use.names = FALSE)))
    stopifnot(length(all_tfs) > 0)
    
    # weights (AUCell) or 1, re-center to mean=1
    wTF <- rep(1, length(all_tfs)); names(wTF) <- all_tfs
    if (!is.null(tf_weight_vec) && length(tf_weight_vec)) {
      tfw <- tf_weight_vec
      names(tfw) <- toupper(names(tfw))        # case-insensitive match
      w <- tfw[toupper(all_tfs)]; w[is.na(w)] <- 1
      wTF <- w / mean(w, na.rm = TRUE)
    }
    
    X <- matrix(0, nrow = length(sets), ncol = length(all_tfs),
                dimnames = list(names(sets), all_tfs))
    for (p in names(sets)) if (length(sets[[p]])) X[p, sets[[p]]] <- wTF[sets[[p]]]
    
    # trim to pathways_n and drop empty rows → prevents all points at (0,0)
    allowed <- unique(as.character(pathway_n$pathway))
    X <- X[intersect(rownames(X), allowed), , drop = FALSE]
    X <- X[rowSums(X) > 0, , drop = FALSE]
    stopifnot(nrow(X) >= 2)
    
    D <- mk_cosine_dist(X)
    coords <- cmdscale(D, k = 2)
    emb_df <- data.frame(pathway = rownames(coords),
                         x = coords[,1], y = coords[,2],
                         stringsAsFactors = FALSE)
  } else if (mode %in% c("tf_centered","tftarget_centered",
                         "tf_centered_unweighted","tftarget_centered_unweighted")) {
    
    # Build features; rows = pathways, cols = TFs or Targets (weights as requested)
    M <- .pathway_features(scMlnet_results, deg = deg,
                           tf_weight_vec = tf_weight_vec,
                           mode = mode, rl_db = receptor_ligand)
    
    M <- .prep_features(M, do_tfidf = TRUE, drop_ubiq = 0.70, zscore_cols = TRUE)
    
    # Keep only pathways that appear in pathway_n or lig_rank_all
    allowed <- unique(c(
      as.character(pathway_n$pathway),
      if (!is.null(lig_rank_all)) as.character(lig_rank_all$pathway)
    ))
    allowed <- allowed[nzchar(allowed)]
    if (length(allowed)) {
      keep <- intersect(rownames(M), allowed)
      M <- M[keep, , drop = FALSE]
    }
    
    # coords <- .embed2D_plus(M, method = embed_method, seed = seed)
    
    # --- 2D embedding: MDS only ---
    coords <- .mds_cosine(M)                  # M is your pathway×feature matrix
    coords <- .compact_coords(coords)         # optional, helps HDBSCAN a bit
    
    
    message("First 3 rows:\n", capture.output(print(head(coords, 3))), collapse="\n")
    
    
    emb_df <- data.frame(
      pathway = rownames(coords),
      x = coords[,1], y = coords[,2],
      stringsAsFactors = FALSE
    )
    # From here on we fall through to the common embed+cluster+plot code
    
  } else {
    stop("Unknown mode.")
  }
  
  # helper
  numify <- function(x) {
    if (is.numeric(x)) return(x)
    x <- gsub(",", "", as.character(x), fixed = TRUE)
    x <- sub("%$", "", x)
    suppressWarnings(as.numeric(x))
  }
  
  # ---------- embed & cluster ----------
  # ---------- embed & cluster (single code path for ALL modes) ----------
  set.seed(seed)
  
  # If emb_df wasn't created by the mode above (ligand/interaction modes),
  # build it from the distance 'D' we computed there.
  if (is.null(emb_df)) {
    emb <- embed2D(D)  # your existing helper that handles n=1/2 nicely
    labs <- rownames(emb)
    emb_df <- data.frame(
      pathway = labs, x = emb[,1], y = emb[,2],
      stringsAsFactors = FALSE
    )
  }
  
  # add new args to simple_learning_run() signature (optional)
  # cluster_method = c("kmeans","hdbscan","dbscan"),
  # dbscan_minPts = NULL, dbscan_eps = NULL,
  
  # ... after emb_df has x,y from UMAP ...
  coords2d   <- as.matrix(emb_df[, c("x","y")])
  coords2d_s <- scale(coords2d)
  
  cluster_method <- match.arg(if (is.null(cluster_method)) "kmeans" else cluster_method,
                              c("kmeans","hdbscan","dbscan"))
  
  label_C <- function(z) paste0("C", as.integer(z))
  
  if (cluster_method == "kmeans") {
    if (is.null(n_clusters)) {
      n_clusters <- .choose_k_by_sil(coords2d, kmax = min(8, nrow(coords2d)-1))
    }
    km <- kmeans(.robust_scale(coords2d), centers = n_clusters, nstart = 200)
    lab <- paste0("C", as.integer(km$cluster))
    
  } else {
    if (!requireNamespace("dbscan", quietly = TRUE)) {
      warning("dbscan not installed; falling back to kmeans")
      if (is.null(n_clusters)) n_clusters <- max(2, round(sqrt(nrow(coords2d))))
      km <- kmeans(.robust_scale(coords2d), centers = n_clusters, nstart = 200)
      lab <- paste0("C", as.integer(km$cluster))
      
    } else if (cluster_method == "hdbscan") {
      
      fit <- auto_hdbscan_on_coords(
        coords2d,
        minPts_grid = 2:max(6, floor(nrow(coords2d)/3)),
        prob_min    = hdbscan_prob_min,        # <-- from args
        k_target    = min_clusters             # <-- bias toward ~4
      )
      
      lab <- as.character(fit$cluster)         # e.g. "N","C1","C3",...
      
      # (a) optionally snap noise to nearest positive cluster
      if (hdbscan_snap_noise && any(lab == "N")) {
        Z   <- .robust_scale(coords2d)
        pos <- lab != "N"
        if (sum(pos) >= 2) {
          ctr <- rowsum(Z[pos, , drop = FALSE], lab[pos]) / as.vector(table(lab[pos]))
          assign_one <- function(i) {
            d <- sqrt(rowSums((t(ctr) - Z[i, ])^2))
            names(which.min(d))
          }
          lab[lab == "N"] <- vapply(which(lab == "N"), assign_one, character(1))
        }
      }
      
      # (b) renumber to C1..Ck without gaps
      pos <- sort(unique(lab[lab != "N"]))
      if (length(pos)) {
        remap <- setNames(paste0("C", seq_along(pos)), pos)
        lab[lab != "N"] <- remap[lab[lab != "N"]]
      }
      
      # (c) if still < min_clusters → hard fallback to k-means(K = min_clusters)
      k_pos <- length(unique(lab[lab != "N"]))
      if (!is.null(min_clusters) && k_pos < min_clusters &&
          match.arg(cluster_fallback) == "kmeans") {
        km  <- kmeans(.robust_scale(coords2d), centers = max(2L, min_clusters), nstart = 200)
        lab <- paste0("C", as.integer(km$cluster))
      }
      
    }
  }
  # sanitize factor levels (put N first only if present)
  # lvl <- sort(unique(lab))
  # lvl <- if ("N" %in% lvl) c("N", setdiff(lvl, "N")) else lvl
  # emb_df$cluster <- factor(lab, levels = lvl)
  
  # ... inside simple_learning_run(), after you compute `lab` for ANY method
  lab <- as.character(lab)
  
  lvl <- sort(unique(lab))
  if ("N" %in% lvl) lvl <- c("N", setdiff(lvl, "N"))  # put N first only if it exists
  
  emb_df$cluster <- factor(lab, levels = lvl)         # <- always rebuild
  emb_df$cluster <- droplevels(emb_df$cluster)        # drop any stale levels
  
  
  
  
  # attach pathway activity for point sizes
  pt <- as.data.frame(pathway_n, stringsAsFactors = FALSE)
  if (!"pathway" %in% names(pt)) pt$pathway <- as.character(pt$entity %||% pt$name)
  sc_col <- intersect(c("score","aggregate","value"), names(pt))[1]
  pt$score <- suppressWarnings(as.numeric(pt[[sc_col]]))
  pt <- pt[!is.na(pt$pathway) & is.finite(pt$score), c("pathway","score")]
  emb_df <- merge(emb_df, pt, by = "pathway", all.x = TRUE)
  emb_df$score[!is.finite(emb_df$score)] <- 0
  
  # point sizes from activity
  rng <- range(emb_df$score, na.rm = TRUE)
  emb_df$pt_size <- if (!is.finite(rng[1]) || diff(rng) == 0) 2.2 else 1.2 + 2.0 * (emb_df$score - rng[1]) / diff(rng)
  
  # palette fixed by cluster labels
  ncl <- nlevels(emb_df$cluster)
  pal <- RColorBrewer::brewer.pal(max(3, min(12, ncl)), "Set3")
  names(pal) <- levels(emb_df$cluster)
  
  # ---------- scatter (ALWAYS the same code) ----------
  used <- attr(coords, "method_used") %||% toupper(embed_method)
  message(sprintf("[Simple] Embedding used: %s  (n_pathways=%d, n_features=%d)",
                  used, nrow(coords), if (exists("X")) ncol(X) else if (exists("Y")) ncol(Y) else ncol(M)))
  
  # M  <- .pathway_features(scMlnet_results, deg, tf_weight_vec,
  #                         mode="tftarget_centered", rl_db=receptor_ligand)
  # M  <- .prep_features(M, do_tfidf=TRUE, drop_ubiq=0.70, zscore_cols=TRUE)
  # um <- .embed2D_plus(M, method="umap", seed=1)
  # md <- .embed2D_plus(M, method="pca",  seed=1)
  # 
  # # Spearman correlation between pairwise distances of the two embeddings
  # print(cor(stats::dist(um), stats::dist(md), method="spearman"))
  # 
  scatter_png <- file.path(out_dir, sprintf("Simple_clusters_scatter_%s_to_%s_%s.png",
                                            receiver, target, norm_label))
  
  p <- ggplot(emb_df, aes(x, y)) +
    geom_point(aes(fill = cluster, size = pt_size),
               shape = 21, color = "black", stroke = 0.3) +
    ggrepel::geom_text_repel(aes(label = pathway), size = 3, max.overlaps = Inf) +
    scale_fill_manual(values = pal, drop = TRUE, name = "Cluster") +
    guides(size = "none") +
    labs(
      title = sprintf("Simple: %s — pathway clusters (embed = %s)",
                      gsub("_", " ", mode, fixed = TRUE), used),
      x = "Dim 1", y = "Dim 2"   # here can change the name of x and y axis
    ) +
    theme_minimal(base_size = 12) +
    theme(
      panel.background = element_rect(fill = "white", colour = NA),
      plot.background  = element_rect(fill = "white", colour = NA)
    )
  ggsave(scatter_png, p, width = 8.5, height = 6.2, dpi = 160, bg = "white")
  paths$scatter <- scatter_png
  
  
  
  
  # ---------- coverage metrics for ranking ----------
  # Use your pathway_coverage_tables to compute Active/Unique
  stopifnot(!is.null(scMlnet_results), !is.null(receptor_ligand), !is.null(deg))
  pw_cov <- pathway_coverage_tables(
    scMlnet_results = list(LigRec = scMlnet_results$LigRec,
                           RecTF  = scMlnet_results$RecTF,
                           TFTar  = scMlnet_results$TFTar),
    receptor_ligand = receptor_ligand,
    deg = deg, tf_weight_vec = tf_weight_vec,
    limit_pathways = emb_df$pathway
  )
  cov_df <- pw_cov[, c("pathway","targ_coverage_frac","targ_unique_sum")]
  emb_df <- merge(emb_df, cov_df, by = "pathway", all.x = TRUE)
  emb_df$targ_coverage_frac[!is.finite(emb_df$targ_coverage_frac)] <- 0
  emb_df$targ_unique_sum  [!is.finite(emb_df$targ_unique_sum)]   <- 0
  
  # attach a numeric 'score' to emb_df from pathway_n
  if (!"score" %in% names(emb_df)) {
    pt <- as.data.frame(pathway_n, stringsAsFactors = FALSE)
    if (!"pathway" %in% names(pt)) {
      # fallbacks if the column is named differently
      pt$pathway <- as.character(pt$entity %||% pt$name)
    }
    sc_col <- intersect(c("score","aggregate","value"), names(pt))[1]
    if (is.na(sc_col)) stop("pathway_n has no score/aggregate/value column")
    pt$score <- suppressWarnings(as.numeric(pt[[sc_col]]))
    pt <- pt[!is.na(pt$pathway) & !is.na(pt$score), c("pathway","score")]
    emb_df <- merge(emb_df, pt, by = "pathway", all.x = TRUE)
  }
  
  
  # cluster palette (fixed by cluster, not pathway)
  ncl <- nlevels(emb_df$cluster)
  pal <- RColorBrewer::brewer.pal(max(3, min(12, ncl)), "Set3")
  cluster_cols <- setNames(rep(pal, length.out = ncl), levels(emb_df$cluster))
  
  
  
  # ---------- ranking inside each cluster ----------
  .res01 <- function(z) {
    r <- range(z, na.rm = TRUE)
    if (!is.finite(r[1]) || r[2] - r[1] < 1e-12) return(rep(0, length(z)))
    (z - r[1]) / (r[2] - r[1])
  }
  
  chosen <- list()
  for (cl in levels(emb_df$cluster)) {
    sub <- emb_df[emb_df$cluster == cl, , drop = FALSE]
    if (!nrow(sub)) next
    
    if (identical(rank_method, "pareto")) {
      pf <- .pareto_front(sub, "targ_coverage_frac", "targ_unique_sum")
      # if more than K on front, pick highest activity
      if (nrow(pf) > k_per_cluster) {
        ord <- order(
          -.force_double(pf$score),
          -.force_double(pf$targ_coverage_frac),
          -.force_double(pf$targ_unique_sum)
        )
        
        pf <- pf[head(ord, k_per_cluster), , drop = FALSE]
      }
      chosen[[cl]] <- pf
    } else {
      rk <-  w_active*.res01(sub$targ_coverage_frac) + w_unique*.res01(sub$targ_unique_sum)
      sub$rk <- rk
      ord <- order(-sub$rk, -.force_double(sub$score))
      chosen[[cl]] <- sub[head(ord, k_per_cluster), , drop = FALSE]
    }
  }
  sel <- do.call(rbind, chosen)
  
  # ---------- treemap: CLUSTERED (activity size) ----------
  message("emb_df cols: ", paste(names(emb_df), collapse = ", "))
  message("head(emb_df):\n", paste(capture.output(print(utils::head(emb_df))), collapse = "\n"))
  
  
  #tm_all <- merge(emb_df[, c("pathway","cluster","score")], data.frame(pathway = unique(emb_df$pathway)), by="pathway", all.y=TRUE)
  tm_all <- unique(emb_df[, c("pathway","cluster"), drop = FALSE])
  tm_all <- merge(tm_all, unique(pathway_n[, c("pathway","score")]), by = "pathway", all.x = TRUE)
  tm_all$score[!is.finite(tm_all$score)] <- 0
  tm_all$cluster <- as.character(tm_all$cluster)
  #tm_all$cluster <- droplevels(tm_all$cluster)
  tm_all$col <- cluster_cols[tm_all$cluster]
  
  tm_all <- tm_all[!is.na(tm_all$pathway) & nzchar(tm_all$pathway), , drop = FALSE]
  tm_all$cluster[is.na(tm_all$cluster) | !nzchar(tm_all$cluster)] <- "Unclustered"
  
  map <- unique(tm_all[, c("cluster","col")])
  pal_by_cluster <- setNames(map$col, map$cluster)
  
  cluster_png <- file.path(out_dir, sprintf("Simple_clustered_treemap_%s_to_%s_%s.png",
                                            receiver, target, norm_label))
  .plot_treemap_2(tm_all, cluster_png, "Clustered pathways (size = activity score)",pal_by_cluster =pal_by_cluster, mirror.y=mirror.y)
  
  # ---------- treemap: PRUNED (ranked by Active/Unique/Pareto), size = activity ----------
  if (is.null(sel) || nrow(sel) == 0L) {
    pruned_png <- .placeholder_png(file.path(out_dir, sprintf(
      "Simple_pruned_treemap_%s_to_%s_%s.png", receiver, target, norm_label
    )), "No pathways selected by ranking")
  } else {
    tm_pruned <- sel[, c("pathway","cluster","score")]
    tm_pruned$cluster <- as.character(tm_pruned$cluster)
    tm_pruned$col <- cluster_cols[tm_pruned$cluster]
    map_p <- unique(tm_pruned[, c("cluster","col")])
    pal_p <- setNames(map_p$col, map_p$cluster)
    
    pruned_png <- file.path(out_dir, sprintf(
      "Simple_pruned_treemap_%s_to_%s_%s.png", receiver, target, norm_label
    ))
    .plot_treemap_2(
      tm_pruned, pruned_png,
      sprintf("Pruned (rank=%s, K=%d)", rank_method, k_per_cluster),
      mirror.y = mirror.y,                      # ← name it explicitly (see note)
      pal_by_cluster = pal_p
    )
  }
  
  # ---------- treemap: PRUNED LR (1–2 per receptor; color = pathway) ----------
  lr_png <- file.path(out_dir, sprintf("Simple_LR_pruned_treemap_%s_to_%s_%s.png",
                                       receiver, target, norm_label))
  
  if (!is.null(int_rank_all) && nrow(int_rank_all)) {
    sets_all <- .lr_targets_by_interaction(scMlnet_results)
    
    cand <- int_rank_all[order(-as.numeric(int_rank_all$score)),
                         c("interaction","ligand","receptor","score","pathway")]
    cand <- cand[!duplicated(cand$interaction), , drop = FALSE]
    
    # --- match reachability naming (UPPERCASE) ---
    cand$interaction <- toupper(as.character(cand$interaction))
    cand$ligand      <- toupper(as.character(cand$ligand))
    cand$receptor    <- toupper(as.character(cand$receptor))
    
    # scorables only + keep pathway for coloring
    cand <- cand[cand$interaction %in% names(sets_all), , drop = FALSE]
    cand <- cand[!is.na(cand$pathway) & nzchar(cand$pathway), , drop = FALSE]
    
    if (nrow(cand)) {
      w_tg <- .target_mag_from_deg(deg)
      au   <- .lr_active_unique(sets_all[cand$interaction], w_tg)
      
      lr_df <- merge(cand, au, by = "interaction", all.x = TRUE)
      lr_df$active[!is.finite(lr_df$active)] <- 0
      lr_df$unique[!is.finite(lr_df$unique)] <- 0
      lr_df$score  <- .force_double(lr_df$score)
      
      # ---- per-RECEPTOR selection: Pareto on (active, unique), cap to K ----
      K_PER_REC <- 1L
      pick_for_rec <- function(d) {
        if (!nrow(d)) return(d[0, ])
        # strong-first order as a tie-breaker
        d <- d[order(-d$unique, -d$active, -d$score), ]
        pf <- .pareto_front(d, x = "active", y = "unique")
        if (!nrow(pf)) pf <- d
        if (nrow(pf) > K_PER_REC) {
          ord <- order(-pf$unique, -pf$active, -pf$score)
          pf <- pf[head(ord, K_PER_REC), , drop = FALSE]
        }
        pf
      }
      keep <- do.call(rbind, lapply(split(lr_df, lr_df$receptor), pick_for_rec))
      
      if (!nrow(keep)) {
        lr_png <- .placeholder_png(lr_png, "No LR interactions after per-receptor pruning")
      } else {
        keep$size <- pmax(keep$score, 1e-8)
        
        # (optional) limit number of receptor groups for readability
        N_REC_SHOW <- 10L
        rec_sum <- aggregate(size ~ receptor, keep, sum)
        top_recs <- head(rec_sum[order(-rec_sum$size), "receptor"], N_REC_SHOW)
        keep <- keep[keep$receptor %in% top_recs, , drop = FALSE]
        
        dir.create(dirname(lr_png), recursive = TRUE, showWarnings = FALSE)
        png(lr_png, width = 720, height = 900)
        pal <- colorRampPalette(RColorBrewer::brewer.pal(12, "Set3"))(
          max(3, length(unique(keep$pathway)))
        )
        treemap::treemap(
          keep,
          index   = c("receptor", "interaction"),  # ← group by RECEPTOR
          vSize   = "size",
          type    = "categorical",
          vColor  = "pathway",                     # ← color by PATHWAY
          palette = pal,
          title   = "Top 10 interactions pruned (color = pathway)",
          border.col      = c("black","black"),
          border.lwds     = c(3,1),
          mirror.y        = mirror.y,
          fontsize.labels = c(0, 9),
          fontcolor.labels= c("black","black"),
          fontface.labels = c(2,1),
          align.labels    = list(c("left","top"), c("center","center")),
          bg.labels       = 0,
          inflate.labels  = TRUE,
          position.legend = "none",
          lowerbound.cex.labels = 0.25
        )
        dev.off()
      }
    } else {
      lr_png <- .placeholder_png(lr_png, "No LR interactions found in network")
    }
  } else {
    lr_png <- .placeholder_png(lr_png, "No LR scores available")
  }
  
  
  
  
  
  
  # ---------- AU bubble (merge the two AU plots into one) ----------
  # x=Active, y=Unique, size=Activity
  # ---------- AU bubble (labels, bigger cluster legend, hide size legend, thick axes) ----------
  au_png <- file.path(out_dir, sprintf("Simple_AU_bubble_%s_to_%s_%s.png",
                                       receiver, target, norm_label))
  bub <- emb_df[, c("pathway","targ_coverage_frac","targ_unique_sum","score","cluster")]
  names(bub) <- c("pathway","active","unique","activity","cluster")
  
  # Ensure cluster is character -> then factor on *observed* values
  bub$cluster <- as.character(bub$cluster)
  clusters <- sort(unique(bub$cluster))
  if (length(clusters) == 0L) clusters <- "C1"  # safety if data is empty
  bub$cluster <- factor(bub$cluster, levels = clusters)
  
  # Build a palette sized to the observed clusters
  base12 <- RColorBrewer::brewer.pal(12, "Set3")
  pal_fn <- colorRampPalette(base12)
  cluster_cols <- setNames(
    pal_fn(max(3, length(clusters)))[seq_along(clusters)],
    clusters
  )
  
  g2 <- ggplot(bub, aes(active, unique)) +
    geom_point(aes(size = activity, fill = cluster), shape = 21, color = "black", stroke = 0.3) +
    ggrepel::geom_text_repel(aes(label = pathway), size = 3.2, max.overlaps = Inf, box.padding = 0.3) +
    scale_fill_manual(values = cluster_cols, name = "Cluster", drop = FALSE) +  # <- robust
    scale_size_continuous(range = c(2.5, 11), name = "") +
    guides(size = "none") +
    labs(x = "Active (fraction of target signal covered)",
         y = "Unique target signal (sum)",
         title = "Active vs Unique (size = Activity)") +
    theme_minimal(base_size = 13) +
    theme(
      legend.title = element_text(size = 12),
      legend.text  = element_text(size = 11),
      axis.line = element_blank(),
      axis.line.x.bottom = element_line(linewidth  = 1.1, colour = "black"),
      axis.line.y.left   = element_line(linewidth  = 1.1, colour = "black"),
      panel.grid.minor = element_blank()
    )
  ggsave(au_png, g2, width = 9, height = 6.2, dpi = 160)
  
  
  return(list(paths = modifyList(paths, list(
    scatter            = scatter_png,
    cluster_treemap    = cluster_png,
    pruned_treemap     = pruned_png,
    lr_pruned_treemap  = lr_png,
    au_bubble          = au_png
  ))))
  
}














########### true GNN

# --- internal: load small PyTorch GNN into the Python session once ---
.ensure_python_gnn <- function() {
  if (!requireNamespace("reticulate", quietly = TRUE))
    stop("Package 'reticulate' is required for the Python GNN.")
  # If torch not present we won't error here (fallback happens in caller).
  if (isTRUE(getOption("true_gnn_loaded", FALSE))) return(invisible(TRUE))
  py <- reticulate::py_run_string
  py("
import math, random
import numpy as np
try:
    import torch
    import torch.nn as nn
    import torch.nn.functional as F
    TORCH_OK = True
except Exception as _e:
    TORCH_OK = False
")
  options(true_gnn_loaded = TRUE)
  invisible(TRUE)
}

# --- tiny utils we reuse ---
.split2 <- function(x) {
  x <- as.character(x); x <- x[nzchar(x)]
  sp <- strsplit(x, "_", fixed = TRUE)
  data.frame(A = vapply(sp, `[`, "", 1),
             B = vapply(sp, `[`, "", 2),
             stringsAsFactors = FALSE)
}
.mode_weighted_magnitude <- function(scMlnet_results, deg) {
  # same logic as in top_pathway(): assign Mode to each (TF, Target) edge,
  # then weight fold changes accordingly and convert to magnitudes
  tf_withmode <- read.delim("database/TFTargetGene_withmode.txt", stringsAsFactors = FALSE)
  # normalize
  to_uc <- function(z) toupper(trimws(z))
  tf_withmode$TF     <- to_uc(tf_withmode$TF)
  # tolerate Target(s)
  tcol <- intersect(names(tf_withmode), c("Target","Targets","target","targets"))
  names(tf_withmode)[match(tcol[1], names(tf_withmode))] <- "Target"
  tf_withmode$Target <- to_uc(tf_withmode$Target)
  tf_withmode$Mode   <- tolower(trimws(if (!"Mode" %in% names(tf_withmode)) "unknown" else tf_withmode$Mode))
  tf_withmode$Mode[!tf_withmode$Mode %in% c("activation","repression","unknown")] <- "unknown"
  
  pairs <- .split2(scMlnet_results$TFTar)
  names(pairs) <- c("TF","Target")
  pairs$TF     <- to_uc(pairs$TF)
  pairs$Target <- to_uc(pairs$Target)
  
  pm <- merge(pairs, tf_withmode[,c("TF","Target","Mode")], by=c("TF","Target"), all.x=TRUE)
  pm$Mode[is.na(pm$Mode)] <- "unknown"
  
  # majority mode by Target
  mm <- aggregate(Mode ~ Target, data = pm, FUN = function(v) {
    a <- sum(v=="activation"); r <- sum(v=="repression")
    if (a>r) "activation" else if (r>a) "repression" else "unknown"
  })
  mode_by_target <- setNames(mm$Mode, mm$Target)
  
  # fold changes → mode-aware magnitudes
  v <- as.numeric(deg$foldChange); names(v) <- toupper(as.character(deg$Gene))
  v <- v[!is.na(v) & !is.nan(v)]
  if (any(is.infinite(v))) {
    fv <- v[is.finite(v)]
    v[v== Inf] <- max(fv); v[v==-Inf] <- min(fv)
  }
  tg <- names(v)
  gene_modes <- mode_by_target[tg]; gene_modes[is.na(gene_modes)] <- "unknown"
  
  is_signed <- any(v < 0)
  w <- rep(1, length(v))
  if (is_signed) {
    pos <- v > 0; neg <- v < 0
    w[pos & gene_modes=="activation"] <- 2
    w[pos & gene_modes=="repression"] <- 0.5
    w[neg & gene_modes=="activation"] <- 0.5
    w[neg & gene_modes=="repression"] <- 2
    v <- abs(v * w)
  } else {
    up   <- v >= 1
    down <- v > 0 & v < 1
    w[up   & gene_modes=="activation"] <- 2
    w[up   & gene_modes=="repression"] <- 0.5
    w[down & gene_modes=="activation"] <- 0.5
    w[down & gene_modes=="repression"] <- 2
    v <- v * w
    idx <- v > 0 & v < 1
    v[idx] <- 1/v[idx]
  }
  v[!is.finite(v)] <- 0
  v
}

# --- robust diffusion fallback (no torch required) ---
# scMlnet_results: list(LigRec, RecTF, TFTar) as "A_B" strings
# deg: data frame with columns Gene and foldChange
# receptor_ligand: your database mapping ligands/receptors to pathways (already loaded globally)
.diffusion_scores_R <- function(scMlnet_results, deg, tf_weight_vec = NULL) {
  suppressPackageStartupMessages({
    library(Matrix)
  })
  
  # 0) Prepare signed -> magnitude vector for targets
  expr <- as.numeric(deg$foldChange); names(expr) <- toupper(as.character(deg$Gene))
  expr <- expr[!is.na(expr)]
  # make positive magnitudes (abs logFC or invert FC<1)
  if (any(expr < 0, na.rm = TRUE)) {
    expr <- abs(expr)
  } else {
    idx <- expr > 0 & expr < 1
    expr[idx] <- 1/expr[idx]
  }
  expr[!is.finite(expr)] <- 0
  
  # 1) Parse edges safely ------------------------------------------------------
  split_pairs <- function(x) {
    x <- as.character(x); x <- x[nzchar(x)]
    sp <- strsplit(x, "_", fixed = TRUE)
    ok <- vapply(sp, function(v) length(v) == 2 && all(nzchar(v)), logical(1))
    sp <- sp[ok]
    if (!length(sp)) return(data.frame(A=character(0), B=character(0)))
    data.frame(A = vapply(sp, `[[`, "", 1), B = vapply(sp, `[[`, "", 2), stringsAsFactors = FALSE)
  }
  
  LR <- split_pairs(scMlnet_results$LigRec)   # A=Ligand, B=Receptor
  RT <- split_pairs(scMlnet_results$RecTF)    # A=Receptor, B=TF
  TT <- split_pairs(scMlnet_results$TFTar)    # A=TF,      B=Target
  
  LR$A <- toupper(LR$A); LR$B <- toupper(LR$B)
  RT$A <- toupper(RT$A); RT$B <- toupper(RT$B)
  TT$A <- toupper(TT$A); TT$B <- toupper(TT$B)
  
  # Keep only targets we have expression magnitude for
  TT <- TT[TT$B %in% names(expr), , drop = FALSE]
  if (!nrow(LR) || !nrow(RT) || !nrow(TT)) {
    return(list(
      lig_rank_all = data.frame(), rec_rank_all = data.frame(),
      int_rank_all = data.frame(), pathway_n = data.frame()
    ))
  }
  
  # 2) Build consistent node order per layer ----------------------------------
  T.targets <- sort(unique(TT$B))
  F.tfs     <- sort(unique(intersect(TT$A, RT$B)))
  R.recs    <- sort(unique(intersect(RT$A, LR$B)))
  L.ligs    <- sort(unique(LR$A))
  
  # Re-restrict edges so all endpoints exist in both layers
  TT <- TT[TT$A %in% F.tfs & TT$B %in% T.targets, , drop = FALSE]
  RT <- RT[RT$A %in% R.recs & RT$B %in% F.tfs,    , drop = FALSE]
  LR <- LR[LR$A %in% L.ligs & LR$B %in% R.recs,   , drop = FALSE]
  
  # Index maps
  idx_T <- setNames(seq_along(T.targets), T.targets)
  idx_F <- setNames(seq_along(F.tfs),     F.tfs)
  idx_R <- setNames(seq_along(R.recs),    R.recs)
  idx_L <- setNames(seq_along(L.ligs),    L.ligs)
  
  # 3) Sparse adjacency (row-normalized) --------------------------------------
  # Shapes we want:  M_TF: (|F| x |T|),  M_RT: (|R| x |F|),  M_LR: (|L| x |R|)
  M_TF <- sparseMatrix(
    i = idx_F[TT$A], j = idx_T[TT$B],
    x = 1, dims = c(length(F.tfs), length(T.targets)), dimnames = list(F.tfs, T.targets)
  )
  M_RT <- sparseMatrix(
    i = idx_R[RT$A], j = idx_F[RT$B],
    x = 1, dims = c(length(R.recs), length(F.tfs)), dimnames = list(R.recs, F.tfs)
  )
  M_LR <- sparseMatrix(
    i = idx_L[LR$A], j = idx_R[LR$B],
    x = 1, dims = c(length(L.ligs), length(R.recs)), dimnames = list(L.ligs, R.recs)
  )
  
  # Row-normalize helper
  row_norm <- function(M) {
    rs <- Matrix::rowSums(M)
    rs[rs == 0] <- 1
    Diagonal(x = 1/rs) %*% M
  }
  M_TF <- row_norm(M_TF)
  M_RT <- row_norm(M_RT)
  M_LR <- row_norm(M_LR)
  
  # 4) Diffusion (forward propagation T -> F -> R -> L) -----------------------
  s_T <- Matrix(expr[T.targets], ncol = 1, dimnames = list(T.targets, NULL))  # |T| x 1
  
  s_F <- M_TF %*% s_T                                    # |F| x 1
  if (!is.null(tf_weight_vec) && length(tf_weight_vec)) {
    w <- as.numeric(tf_weight_vec[rownames(s_F)]); w[is.na(w)] <- 1
    w <- w / mean(w, na.rm = TRUE)
    s_F <- s_F * w
  }
  
  s_R <- M_RT %*% s_F                                    # |R| x 1
  s_L <- M_LR %*% s_R                                    # |L| x 1
  
  # 5) Build outputs (ligand, receptor, interaction, pathway) -----------------
  # Ligand scores
  lig_rank_all <- data.frame(
    ligand     = rownames(s_L),
    score      = as.numeric(s_L),
    stringsAsFactors = FALSE
  )
  lig_rank_all <- lig_rank_all[order(-lig_rank_all$score), , drop = FALSE]
  lig_rank_all$pathway <- receptor_ligand$pathway_name[
    match(lig_rank_all$ligand, receptor_ligand$Ligand.ApprovedSymbol)
  ]
  lig_rank_all$score.perc     <- 100 * lig_rank_all$score / pmax(sum(lig_rank_all$score), 1e-12)
  lig_rank_all$score.perc.txt <- sprintf("%.1f%%", lig_rank_all$score.perc)
  
  # Receptor scores
  rec_rank_all <- data.frame(
    receptor   = rownames(s_R),
    score      = as.numeric(s_R),
    stringsAsFactors = FALSE
  )
  rec_rank_all <- rec_rank_all[order(-rec_rank_all$score), , drop = FALSE]
  rec_rank_all$pathway <- receptor_ligand$pathway_name[
    match(rec_rank_all$receptor, receptor_ligand$Receptor.ApprovedSymbol)
  ]
  rec_rank_all$score.perc     <- 100 * rec_rank_all$score / pmax(sum(rec_rank_all$score), 1e-12)
  rec_rank_all$score.perc.txt <- sprintf("%.1f%%", rec_rank_all$score.perc)
  
  # Interaction scores (split receptor score across incoming ligands)
  deg_in <- as.numeric(Matrix::colSums(M_LR))           # per-receptor in-degree
  deg_in[deg_in == 0] <- 1
  pair_w <- 1/deg_in                                     # equal split among ligands to a receptor
  names(pair_w) <- colnames(M_LR)
  
  LR_use <- LR[LR$A %in% rownames(M_LR) & LR$B %in% colnames(M_LR), , drop = FALSE]
  if (nrow(LR_use)) {
    int_rank_all <- within(LR_use, {
      score <- as.numeric(s_R[LR_use$B, 1]) * pair_w[LR_use$B]
      interaction <- paste0(A, "_", B)
      pathway <- receptor_ligand$pathway_name[
        match(A, receptor_ligand$Ligand.ApprovedSymbol)
      ]
    })
    int_rank_all <- int_rank_all[, c("interaction","A","B","score","pathway")]
    names(int_rank_all)[2:3] <- c("ligand","receptor")
    int_rank_all <- int_rank_all[order(-int_rank_all$score), , drop = FALSE]
    int_rank_all$score.perc     <- 100 * int_rank_all$score / pmax(sum(int_rank_all$score), 1e-12)
    int_rank_all$score.perc.txt <- sprintf("%.1f%%", int_rank_all$score.perc)
  } else {
    int_rank_all <- data.frame()
  }
  
  # Pathway summary (sum ligand scores per pathway)
  pathway_n <- aggregate(score ~ pathway, data = na.omit(lig_rank_all[, c("pathway","score")]), sum)
  pathway_n <- pathway_n[order(-pathway_n$score), , drop = FALSE]
  total <- pmax(sum(pathway_n$score), 1e-12)
  pathway_n$score.perc     <- 100 * pathway_n$score / total
  pathway_n$score.perc.txt <- sprintf("%.1f%%", pathway_n$score.perc)
  
  list(
    lig_rank_all = lig_rank_all,
    rec_rank_all = rec_rank_all,
    int_rank_all = int_rank_all,
    pathway_n    = pathway_n
  )
}


# ======================================================================
# TRUE GNN (Python) + plotting wrapper
# ======================================================================
true_gnn_run <- function(scMlnet_results,
                         deg,
                         out_dir,
                         receiver, target, norm_label,
                         tf_weight_vec = NULL,     # optional AUCell weights (we use as feature)
                         epochs = 150, d_model = 32, learn_rate = 1e-2, seed = 1,
                         topK = 3, edges_per_path = 20,
                         k_rt = 40, k_tft = 80, q_rt = 0.80, q_tft = 0.85,
                         mirror.y = FALSE,
                         palette = c(
                           NOTCH="#9AA556", WNT="#E98FC5", ncWNT="#64B1DF", EGF="#66E1E1",
                           IGF="#D5AB90", NRG="#8B9D63", BMP="#E88E50", ACTIVIN="#C2E579",
                           TGFb="#DBC6DE", NT="#DE8D9F", VEGF="#AAD8E1", HH="#6EE07F",
                           FGF="#8580DF", HGF="#D9E83D", PDGF="#62E4B7", GDNF="#E2E2A1"
                         )) {
  dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)
  
  
  scMlnet_results
  
  # 1) Prepare graph (names only)
  LR <- .split2(scMlnet_results$LigRec); names(LR) <- c("L","R")
  RT <- .split2(scMlnet_results$RecTF); names(RT) <- c("R","T")
  TT <- .split2(scMlnet_results$TFTar); names(TT) <- c("T","G")
  
  ligs <- sort(unique(LR$L)); recs <- sort(unique(LR$R))
  tfs  <- sort(unique(RT$T)); tars <- sort(unique(TT$G))
  
  # 2) Mode-aware target magnitudes (same semantics as your code)
  targ_mag <- .mode_weighted_magnitude(scMlnet_results, deg)   # named by TARGET
  
  # TF feature = sum of target magnitudes for that TF (optionally * AUCell)
  TFt <- .split2(scMlnet_results$TFTar); names(TFt) <- c("TF","Target")
  TFt$Target <- toupper(TFt$Target); TFt$TF <- toupper(TFt$TF)
  mag <- targ_mag[TFt$Target]; mag[is.na(mag)] <- 0
  tf_feat <- tapply(mag, TFt$TF, sum); tf_feat <- tf_feat[sort(unique(names(tf_feat)))]
  # AUCell weights if present
  if (!is.null(tf_weight_vec) && length(tf_weight_vec)) {
    w <- tf_weight_vec[names(tf_feat)]; w[is.na(w)] <- 1
    w <- w / mean(w, na.rm = TRUE)
    tf_feat <- tf_feat * w
  }
  tf_feat <- tf_feat[match(tfs, names(tf_feat))]; tf_feat[is.na(tf_feat)] <- 0
  
  tar_feat <- targ_mag[tars]; tar_feat[is.na(tar_feat)] <- 0
  
  # 3) Try Python GNN; if torch missing, fall back to diffusion scoring
  .ensure_python_gnn()
  have_torch <- FALSE
  try({
    have_torch <- reticulate::py_eval("TORCH_OK")
  }, silent = TRUE)
  
  if (!isTRUE(have_torch)) {
    message("PyTorch not available → using fast diffusion fallback (still different from your scoring).")
    LR_sc <- .diffusion_scores_R(scMlnet_results, targ_mag)
  } else {
    # define / load the Python model code (only once per session)
    reticulate::py_run_string("
import numpy as np, math, random
import torch
import torch.nn as nn
import torch.nn.functional as F

def _to_tensor(x, dtype=torch.float32, device='cpu'):
    return torch.as_tensor(x, dtype=dtype, device=device)

class MiniHeteroGNN(nn.Module):
    def __init__(self, nL, nR, nT, nG, d=32):
        super().__init__()
        self.embL = nn.Embedding(nL, d)
        self.embR = nn.Embedding(nR, d)
        self.embT = nn.Embedding(nT, d)
        self.embG = nn.Embedding(nG, d)
        # feature projections (TF & Target have 1-d features)
        self.tf_proj  = nn.Linear(1, d)
        self.tar_proj = nn.Linear(1, d)
        # self transforms
        self.WL = nn.Linear(d, d); self.WR = nn.Linear(d, d)
        self.WT = nn.Linear(d, d); self.WG = nn.Linear(d, d)
        # relation transforms (source -> target)
        self.W_LR = nn.Linear(d, d, bias=False)
        self.W_RT = nn.Linear(d, d, bias=False)
        self.W_TG = nn.Linear(d, d, bias=False)
        # decoders for link probs
        self.D_LR = nn.Bilinear(d, d, 1)
        self.D_RT = nn.Bilinear(d, d, 1)
        self.D_TG = nn.Bilinear(d, d, 1)

    def forward(self, lr_idx, rt_idx, tg_idx, tf_feat, tar_feat):
        # initial states
        L = self.embL.weight
        R = self.embR.weight
        T = self.embT.weight + self.tf_proj(tf_feat)
        G = self.embG.weight + self.tar_proj(tar_feat)

        # one round message passing (L->R->T->G), with residual self transforms
        if lr_idx.numel() > 0:
            mR = torch.zeros_like(R)
            src, dst = lr_idx[:,0], lr_idx[:,1]
            m = self.W_LR(L[src])
            mR.index_add_(0, dst, m)
            deg = torch.bincount(dst, minlength=R.shape[0]).clamp(min=1).float().unsqueeze(1)
            R = torch.relu(self.WR(R) + mR / deg)

        if rt_idx.numel() > 0:
            mT = torch.zeros_like(T)
            src, dst = rt_idx[:,0], rt_idx[:,1]
            m = self.W_RT(R[src])
            mT.index_add_(0, dst, m)
            deg = torch.bincount(dst, minlength=T.shape[0]).clamp(min=1).float().unsqueeze(1)
            T = torch.relu(self.WT(T) + mT / deg)

        if tg_idx.numel() > 0:
            mG = torch.zeros_like(G)
            src, dst = tg_idx[:,0], tg_idx[:,1]
            m = self.W_TG(T[src])
            mG.index_add_(0, dst, m)
            deg = torch.bincount(dst, minlength=G.shape[0]).clamp(min=1).float().unsqueeze(1)
            G = torch.relu(self.WG(G) + mG / deg)

        return L, R, T, G

    def score_edges(self, L, R, T, G, lr_idx, rt_idx, tg_idx):
        def _score(Bi, X, Y, idx):
            if idx.numel() == 0: 
                return torch.zeros(0,1, device=X.device)
            s = torch.sigmoid(Bi(X[idx[:,0]], Y[idx[:,1]]))
            return s
        sLR = _score(self.D_LR, L, R, lr_idx)
        sRT = _score(self.D_RT, R, T, rt_idx)
        sTG = _score(self.D_TG, T, G, tg_idx)
        return sLR, sRT, sTG

def train_gnn_and_score(
        nL, nR, nT, nG,
        lr_pairs, rt_pairs, tg_pairs,
        tf_feat, tar_feat,
        epochs=150, d=32, lr=1e-2, seed=1, device='cpu'):
    torch.manual_seed(seed); np.random.seed(seed); random.seed(seed)
    model = MiniHeteroGNN(nL, nR, nT, nG, d=d).to(device)

    # tensors
    tf_feat = _to_tensor(tf_feat, device=device).view(-1,1)
    tar_feat = _to_tensor(tar_feat, device=device).view(-1,1)
    lr_idx = _to_tensor(lr_pairs, dtype=torch.long, device=device)
    rt_idx = _to_tensor(rt_pairs, dtype=torch.long, device=device)
    tg_idx = _to_tensor(tg_pairs, dtype=torch.long, device=device)

    opt = torch.optim.Adam(model.parameters(), lr=lr)
    bce = nn.BCELoss()

    def neg_sample(idx, n_src, n_dst):
        # corrupt destination uniformly
        if idx.numel() == 0:
            return torch.zeros(0,2, dtype=torch.long, device=device)
        src = idx[:,0]
        dst = torch.randint(0, n_dst, (idx.shape[0],), device=device)
        return torch.stack([src, dst], dim=1)

    for ep in range(epochs):
        model.train(); opt.zero_grad()
        L,R,T,G = model(lr_idx, rt_idx, tg_idx, tf_feat, tar_feat)

        # positives
        sLR_p, sRT_p, sTG_p = model.score_edges(L,R,T,G, lr_idx, rt_idx, tg_idx)
        # negatives
        lr_neg = neg_sample(lr_idx, nL, nR)
        rt_neg = neg_sample(rt_idx, nR, nT)
        tg_neg = neg_sample(tg_idx, nT, nG)
        sLR_n, sRT_n, sTG_n = model.score_edges(L,R,T,G, lr_neg, rt_neg, tg_neg)

        loss = 0
        if sLR_p.numel(): loss = loss + bce(sLR_p, torch.ones_like(sLR_p)) + bce(sLR_n, torch.zeros_like(sLR_n))
        if sRT_p.numel(): loss = loss + bce(sRT_p, torch.ones_like(sRT_p)) + bce(sRT_n, torch.zeros_like(sRT_n))
        if sTG_p.numel(): loss = loss + bce(sTG_p, torch.ones_like(sTG_p)) + bce(sTG_n, torch.zeros_like(sTG_n))

        loss.backward(); opt.step()

    # final scores
    model.eval()
    with torch.no_grad():
        L,R,T,G = model(lr_idx, rt_idx, tg_idx, tf_feat, tar_feat)
        sLR, sRT, sTG = model.score_edges(L,R,T,G, lr_idx, rt_idx, tg_idx)
    return {
        'lr_scores': sLR.squeeze(-1).cpu().numpy(),
        'rt_scores': sRT.squeeze(-1).cpu().numpy(),
        'tg_scores': sTG.squeeze(-1).cpu().numpy()
    }
")
    
    # Build index maps and numeric edge arrays for Python
    idx <- function(v) { i <- seq_along(v)-1L; names(i) <- v; i }  # 0-based for Python
    iL <- idx(ligs); iR <- idx(recs); iT <- idx(tfs); iG <- idx(tars)
    
    lr_pairs <- cbind(iL[LR$L], iR[LR$R])
    rt_pairs <- cbind(iR[RT$R], iT[RT$T])
    tg_pairs <- cbind(iT[TT$T], iG[TT$G])
    
    py <- reticulate::py_call
    out <- reticulate::py$train_gnn_and_score(
      as.integer(length(ligs)), as.integer(length(recs)),
      as.integer(length(tfs)),  as.integer(length(tars)),
      lr_pairs, rt_pairs, tg_pairs,
      as.numeric(tf_feat[match(tfs, names(tf_feat))]),
      as.numeric(tar_feat[match(tars, names(tar_feat))]),
      epochs = as.integer(epochs), d = as.integer(d_model),
      lr = learn_rate, seed = as.integer(seed)
    )
    
    LR_sc <- LR
    LR_sc$score <- as.numeric(out$lr_scores)
  }
  
  # 4) Build NEW ranking tables from learned scores
  LR_sc$interaction <- paste0(LR_sc$L, "_", LR_sc$R)
  LR_sc$ligand <- LR_sc$L; LR_sc$receptor <- LR_sc$R
  LR_sc$pathway <- receptor_ligand$pathway_name[
    match(LR_sc$ligand, receptor_ligand$Ligand.ApprovedSymbol)
  ]
  LR_sc <- LR_sc[!is.na(LR_sc$pathway), , drop = FALSE]
  
  # interaction table
  int_rank_all <- aggregate(score ~ interaction + ligand + receptor + pathway, LR_sc, sum)
  int_rank_all <- int_rank_all[order(-int_rank_all$score), , drop = FALSE]
  int_rank_all$score.perc <- 100 * int_rank_all$score / sum(int_rank_all$score)
  int_rank_all$score.perc.txt <- sprintf("%.1f%%", int_rank_all$score.perc)
  
  # receptor table (sum incoming LR scores)
  rec_rank_all <- aggregate(score ~ receptor, LR_sc, sum)
  rec_rank_all <- rec_rank_all[order(-rec_rank_all$score), , drop = FALSE]
  rec_rank_all$pathway <- receptor_ligand$pathway_name[
    match(rec_rank_all$receptor, receptor_ligand$Receptor.ApprovedSymbol)
  ]
  rec_rank_all$score.perc <- 100 * rec_rank_all$score / sum(rec_rank_all$score)
  rec_rank_all$score.perc.txt <- sprintf("%.1f%%", rec_rank_all$score.perc)
  
  # pathway table = sum of learned LR scores grouped by ligand pathway
  pathway_n <- aggregate(score ~ pathway, int_rank_all, sum)
  pathway_n <- pathway_n[order(-pathway_n$score), , drop = FALSE]
  pathway_n$score.perc <- 100 * pathway_n$score / sum(pathway_n$score)
  pathway_n$score.perc.txt <- sprintf("%.1f%%", pathway_n$score.perc)
  
  
  # ---- NEW: ligand table from learned LR scores ----
  lig_rank_all <- aggregate(score ~ ligand + pathway, data = int_rank_all, sum)
  lig_rank_all <- lig_rank_all[order(-lig_rank_all$score), , drop = FALSE]
  tot <- pmax(sum(lig_rank_all$score), 1e-12)
  lig_rank_all$score.perc     <- 100 * lig_rank_all$score / tot
  lig_rank_all$score.perc.txt <- sprintf("%.1f%%", lig_rank_all$score.perc)
  
  # 5) Plot the NEW treemap (GNN-only scores)
  gnn_png <- file.path(out_dir, sprintf("TrueGNN_pathways_%s_to_%s_%s.png", receiver, target, norm_label))
  plot_gnn_pathway_treemap(
    gnn_df = pathway_n[, c("pathway","score")],
    out_png = gnn_png,
    title = sprintf("GNN Pathways — %s \u2192 %s", receiver, target),
    palette = palette,
    mirror.y = mirror.y
  )
  
  # ------- helpers to mirror deep aesthetics + pruning (local to this fn) ----
  .keep_top_tf_targets <- function(RT, TT, max_tf = 12, max_targets = 12) {
    RT <- RT[, c("receptor","tf","w"), drop = FALSE]
    TT <- TT[, c("tf","target","w"),   drop = FALSE]
    RT$w <- as.numeric(RT$w); TT$w <- as.numeric(TT$w)
    if (!nrow(RT)) return(list(RT = RT[0,], TT = TT[0,], keep_receptors = character(0)))
    tf_rank <- aggregate(w ~ tf, data = RT, sum)
    keep_tf <- head(tf_rank[order(-tf_rank$w), "tf"], max_tf)
    RT <- RT[RT$tf %in% keep_tf, , drop = FALSE]
    TT <- TT[TT$tf %in% keep_tf, , drop = FALSE]
    if (nrow(TT)) {
      tg_rank <- aggregate(w ~ target, data = TT, sum)
      keep_tg <- head(tg_rank[order(-tg_rank$w), "target"], max_targets)
      TT <- TT[TT$target %in% keep_tg, , drop = FALSE]
    } else {
      RT <- RT[0,]
    }
    if (nrow(TT)) {
      keep_tf2 <- unique(TT$tf)
      RT <- RT[RT$tf %in% keep_tf2, , drop = FALSE]
    } else {
      RT <- RT[0,]
    }
    list(RT = RT, TT = TT, keep_receptors = unique(RT$receptor))
  }
  
  .build_graph <- function(LR, RT, TT) {
    suppressPackageStartupMessages(library(igraph))
    .tag <- function(sym, layer) paste(sym, layer, sep = "||")
    nodes <- rbind(
      data.frame(name = .tag(unique(LR$ligand),   "Ligand"),   layer = "ligand",  label = unique(LR$ligand)),
      data.frame(name = .tag(unique(LR$receptor), "Receptor"), layer = "receptor",label = unique(LR$receptor)),
      data.frame(name = .tag(unique(RT$tf),       "TF"),       layer = "tf",      label = unique(RT$tf)),
      data.frame(name = .tag(unique(TT$target),   "Target"),   layer = "target",  label = unique(TT$target))
    )
    nodes <- unique(nodes)
    e_lr <- if (nrow(LR)) data.frame(from = .tag(LR$ligand, "Ligand"),
                                     to   = .tag(LR$receptor,"Receptor"),
                                     w = LR$w, kind="LR") else NULL
    e_rt <- if (nrow(RT)) data.frame(from = .tag(RT$receptor,"Receptor"),
                                     to   = .tag(RT$tf,"TF"),
                                     w = RT$w, kind="RT") else NULL
    e_tt <- if (nrow(TT)) data.frame(from = .tag(TT$tf,"TF"),
                                     to   = .tag(TT$target,"Target"),
                                     w = TT$w, kind="TT") else NULL
    edges <- rbind(e_lr, e_rt, e_tt); if (is.null(edges) || !nrow(edges)) return(igraph::make_empty_graph())
    g <- igraph::graph_from_data_frame(edges, vertices = nodes, directed = TRUE)
    # widths normalized per layer
    norm01 <- function(x){ x <- as.numeric(x); r <- range(x, na.rm=TRUE); if(!is.finite(r[1])||r[1]==r[2]) return(rep(.5,length(x))); (x-r[1])/(r[2]-r[1]) }
    EW <- numeric(nrow(edges))
    for (kk in c("LR","RT","TT")) {
      ii <- which(edges$kind==kk)
      if (length(ii)) EW[ii] <- 1 + 7*norm01(edges$w[ii])
    }
    E(g)$width <- EW
    g
  }
  
  .plot_explain <- function(
    g, out_png, title,
    width_px = 860, height_px = 540, dpi = 200,
    write_hires = TRUE, hires_px = 2200, hires_dpi = 320,
    show_legend = FALSE
  ) {
    suppressPackageStartupMessages(library(igraph))
    .open_px <- function(path, wpx, hpx, res) {
      dir.create(dirname(path), recursive = TRUE, showWarnings = FALSE)
      win <- wpx / res; hin <- hpx / res
      if (requireNamespace("ragg", quietly = TRUE)) {
        ragg::agg_png(path, width = win, height = hin, units = "in", res = res, background = "white")
      } else {
        png(path, width = win, height = hin, units = "in", res = res, type = "cairo", bg = "white")
      }
    }
    draw_once <- function() {
      if (vcount(g) == 0) { plot.new(); title(main = paste0(title, " — no edges")); return(invisible()) }
      layers <- c("ligand","receptor","tf","target")
      cols   <- setNames(c("#42a5f5","#ffb74d","#8d6e63","#66bb6a"), layers)
      vcol   <- cols[V(g)$layer]
      y <- setNames(c(3,2,1,0), layers)[V(g)$layer]
      x <- ave(seq_along(y), y, FUN = function(ix) seq(0.03, 0.97, length.out = length(ix)))
      lay <- cbind(x, y)
      ec <- ifelse(E(g)$kind=="LR", grDevices::adjustcolor("grey20", 0.55),
                   ifelse(E(g)$kind=="RT", grDevices::adjustcolor("grey40", 0.40),
                          grDevices::adjustcolor("grey20", 0.28)))
      op <- par(no.readonly = TRUE); on.exit(par(op), add = TRUE)
      par(mar = c(2.0, 1.4, 3.2, 1.4), xpd = NA, family = "sans")
      plot(g, layout = lay, rescale = FALSE, xlim = c(-0.02,1.02), ylim = c(-0.05,3.15), asp = 0,
           vertex.size = 8, vertex.color = vcol, vertex.frame.color = "grey25",
           vertex.label = V(g)$label, vertex.label.cex = 0.37,
           vertex.label.family = "sans", vertex.label.font = 2,
           edge.arrow.size = 0.15, edge.width = E(g)$width, edge.color = ec, main = title)
      if (isTRUE(show_legend)) legend("topright", bty="n",
                                      legend=c("Ligand","Receptor","TF","Target","LR score","RT score","TF→Target score"),
                                      pch=c(16,16,16,16,NA,NA,NA), col=c(cols,NA,NA,NA), pt.cex=.9, cex=.85, lwd=c(NA,NA,NA,NA,3,3,3))
    }
    .open_px(out_png, width_px, height_px, dpi); draw_once(); dev.off()
    if (isTRUE(write_hires)) { hi <- sub("\\.png$", "_hires.png", out_png)
    .open_px(hi, hires_px, round(hires_px * height_px / width_px), hires_dpi); draw_once(); dev.off() }
  }
  # --------------------------------------------------------------------------
  
  # 6) Explanations using the learned scores — same pruning + look as deep ----
  ord <- order(-pathway_n$score)
  top_paths <- head(as.character(pathway_n$pathway[ord]), topK)
  
  # (i) prune network once to full L→R→TF→Target paths
  nl_pruned <- prune_net_to_full_paths(
    list(LigRec = scMlnet_results$LigRec,
         RecTF  = scMlnet_results$RecTF,
         TFTar  = scMlnet_results$TFTar),
    sep = "_"
  )
  LR_all <- .pairs_to_df(nl_pruned$LigRec, "ligand","receptor")
  RT_all <- .pairs_to_df(nl_pruned$RecTF,  "receptor","tf")
  TT_all <- .pairs_to_df(nl_pruned$TFTar,  "tf","target")
  
  # (ii) build per-pathway gating tables from learned scores
  lig_tbl <- aggregate(score ~ ligand + pathway,
                       data = int_rank_all[, c("ligand","pathway","score")], sum)
  rec_tbl <- rec_rank_all[, c("receptor","pathway","score")]
  
  explain_pngs <- character(0)
  for (i in seq_along(top_paths)) {
    pw <- top_paths[i]
    
    # filter pairs to the annotated pathway (receptor_ligand DB)
    edges_pw <- .filter_layers_for_pathway_pairs(LR_all, RT_all, TT_all, pw, receptor_ligand)
    # decorate weights from learned ligand/receptor tables (normalized per layer)
    edges_pw <- .decorate_weights(edges_pw$LR, edges_pw$RT, edges_pw$TT, pw, lig_tbl, rec_tbl)
    
    # keep only the strongest TFs/Targets globally; then drop overhang L/R
    sel <- .keep_top_tf_targets(edges_pw$RT, edges_pw$TT, max_tf = 12, max_targets = 12)
    RTp <- sel$RT
    TTp <- sel$TT
    LRp <- edges_pw$LR[edges_pw$LR$receptor %in% sel$keep_receptors, , drop = FALSE]
    
    g_pw <- .build_graph(LRp, RTp, TTp)
    
    f <- file.path(out_dir, sprintf("TrueGNN_Explain_%d_%s_to_%s_%s_%s.png",
                                    i, receiver, target, norm_label, gsub("[^A-Za-z0-9+_.-]","_", pw)))
    .plot_explain(
      g_pw,
      out_png = f,
      title   = sprintf("GNN . %s — rank #%d", pw, i),
      width_px = 860, height_px = 540, dpi = 200,
      write_hires = TRUE
    )
    explain_pngs <- c(explain_pngs, f)
  }
  
  invisible(list(
    tables = list(
      pathway_n    = pathway_n,
      int_rank_all = int_rank_all,
      rec_rank_all = rec_rank_all,
      lig_rank_all = lig_rank_all   # <- add this
    ),
    paths  = list(gnn = gnn_png, explain = explain_pngs),
    top_paths = top_paths
  ))
}













# ---- Robust pruning to full L→R→TF→Target paths (with safe fallback) ----
# --- ensure only full L->R->TF->Target paths are kept -----------------------
prune_net_to_full_paths <- function(nl, sep = "_", deg_genes = NULL, require_deg_targets = FALSE) {
  to_df <- function(x, a, b) {
    if (is.null(x)) return(data.frame(a=character(0), b=character(0)))
    if (is.data.frame(x) && all(c(a,b) %in% names(x))) {
      out <- x[, c(a,b), drop = FALSE]; names(out) <- c(a,b); out[] <- lapply(out, as.character); return(out)
    }
    x <- as.character(x); x <- x[nzchar(x)]
    if (!length(x)) return(data.frame(a=character(0), b=character(0)))
    parts <- strsplit(x, sep, fixed = TRUE)
    A <- vapply(parts, function(p) if (length(p)) p[1] else "", "")
    B <- vapply(parts, function(p) if (length(p) >= 2) paste(p[-1], collapse = sep) else "", "")
    out <- data.frame(A=trimws(A), B=trimws(B), stringsAsFactors = FALSE); names(out) <- c(a,b)
    out[nzchar(out[[a]]) & nzchar(out[[b]]), , drop = FALSE]
  }
  to_vec <- function(df, a, b) if (!nrow(df)) character(0) else paste(df[[a]], df[[b]], sep = sep)
  
  LR <- to_df(nl$LigRec, "ligand", "receptor")
  RT <- to_df(nl$RecTF,  "receptor", "tf")
  TT <- to_df(nl$TFTar,  "tf", "target")
  
  if (require_deg_targets && !is.null(deg_genes)) {
    deg_uc <- toupper(as.character(deg_genes))
    TT <- TT[toupper(TT$target) %in% deg_uc, , drop = FALSE]
  }
  
  keep_receptors <- intersect(unique(LR$receptor), unique(RT$receptor))
  RT <- RT[RT$receptor %in% keep_receptors, , drop = FALSE]
  LR <- LR[LR$receptor %in% keep_receptors, , drop = FALSE]
  
  keep_tf <- intersect(unique(RT$tf), unique(TT$tf))
  RT <- RT[RT$tf %in% keep_tf, , drop = FALSE]
  TT <- TT[TT$tf %in% keep_tf, , drop = FALSE]
  
  list(
    LigRec = to_vec(unique(LR), "ligand","receptor"),
    RecTF  = to_vec(unique(RT), "receptor","tf"),
    TFTar  = to_vec(unique(TT), "tf","target")
  )
}




empty_canon_df <- function() {
  data.frame(
    kind        = character(0),
    entity      = character(0),
    pathway     = character(0),
    score       = numeric(0),
    score.perc  = character(0),
    p.value     = numeric(0),
    significant = character(0),
    stringsAsFactors = FALSE
  )
}


`%||%` <- function(a,b) if (!is.null(a)) a else b

ensure_true_tables <- function(tr) {
  tbl <- tr$tables %||% tr
  # derive lig_rank_all if not present
  if (is.null(tbl$lig_rank_all) && !is.null(tbl$int_rank_all) && nrow(tbl$int_rank_all)) {
    lr <- aggregate(score ~ ligand + pathway, data = tbl$int_rank_all, sum)
    lr$score.perc     <- 100 * lr$score / sum(lr$score)
    lr$score.perc.txt <- sprintf("%.1f%%", lr$score.perc)
    tbl$lig_rank_all  <- lr
  }
  tbl
}


# returns lists of sets: which TFs/Targets each pathway or LR interaction can reach
# helpers (one-time; move to a utils file if you prefer)
# build pathway -> reachable targets (and TFs)
# --- build per-pathway target reachability from scMLnet + receptor_ligand ---
build_reachability <- function(scMlnet_results, receptor_ligand) {
  toUC <- function(z) toupper(trimws(as.character(z)))
  LR <- .pairs_to_df(scMlnet_results$LigRec, "ligand","receptor");  LR[] <- lapply(LR, toUC)
  RT <- .pairs_to_df(scMlnet_results$RecTF,  "receptor","tf");      RT[] <- lapply(RT, toUC)
  TT <- .pairs_to_df(scMlnet_results$TFTar,  "tf","target");        TT[] <- lapply(TT, toUC)
  
  rec_by_lig <- split(LR$receptor, LR$ligand)
  tf_by_rec  <- split(RT$tf,       RT$receptor)
  tg_by_tf   <- split(TT$target,   TT$tf)
  
  db <- receptor_ligand[, c("pathway_name","Ligand.ApprovedSymbol")]
  names(db) <- c("pathway","ligand"); db[] <- lapply(db, toUC)
  pths <- sort(unique(db$pathway))
  
  targets_by_pathway <- setNames(vector("list", length(pths)), pths)
  tfs_by_pathway     <- setNames(vector("list", length(pths)), pths)
  
  for (p in pths) {
    Lp <- unique(db$ligand[db$pathway == p])
    tf_set <- character(0)
    tg_set <- character(0)
    for (lg in Lp) {
      rc <- rec_by_lig[[lg]]; if (is.null(rc)) next
      for (r in unique(rc)) {
        tfs <- tf_by_rec[[r]]; if (is.null(tfs)) next
        tf_set <- union(tf_set, tfs)
        for (f in unique(tfs)) {
          tgs <- tg_by_tf[[f]]; if (is.null(tgs)) next
          tg_set <- union(tg_set, tgs)
        }
      }
    }
    tfs_by_pathway[[p]]     <- sort(unique(tf_set))
    targets_by_pathway[[p]] <- sort(unique(tg_set))
  }
  list(targets_by_pathway = targets_by_pathway,
       tfs_by_pathway     = tfs_by_pathway)
}




# table of Active vs Unique (targets)
# --- weighted coverage summary used in "Active vs Unique" plots -------------
pathway_coverage_tables <- function(scMlnet_results, receptor_ligand, deg, tf_weight_vec = NULL,
                                    limit_pathways = NULL) {
  R <- build_reachability(scMlnet_results, receptor_ligand)
  w_tg <- .target_mag_from_deg(deg)  # abs(logFC) with collapse
  
  pths <- names(R$targets_by_pathway)
  if (!is.null(limit_pathways)) pths <- intersect(pths, limit_pathways)
  if (!length(pths)) return(data.frame(pathway=character(), targ_coverage_frac=numeric(), targ_unique_sum=numeric()))
  
  totalW <- sum(w_tg, na.rm = TRUE); if (!is.finite(totalW) || totalW == 0) totalW <- 1
  union_others <- function(i) {
    others <- R$targets_by_pathway[setdiff(pths, pths[i])]
    sort(unique(unlist(others, use.names = FALSE)))
  }
  
  out <- lapply(seq_along(pths), function(i) {
    p  <- pths[i]
    tg <- R$targets_by_pathway[[p]]
    cov_w <- sum(w_tg[intersect(names(w_tg), tg)], na.rm = TRUE) / totalW
    uniq <- setdiff(tg, union_others(i))
    uniq_w <- sum(w_tg[intersect(names(w_tg), uniq)], na.rm = TRUE)
    data.frame(pathway = p, targ_coverage_frac = cov_w, targ_unique_sum = uniq_w, stringsAsFactors = FALSE)
  })
  do.call(rbind, out)
}


# --- greedy minimal cover up to a threshold of weighted targets -------------
greedy_cover <- function(sets, w, threshold = 0.8) {
  targetW <- sum(w, na.rm = TRUE); if (!is.finite(targetW) || targetW == 0) targetW <- 1
  need <- threshold * targetW
  chosen <- character(0); covered <- character(0); have <- 0
  while (have < need) {
    best <- NA_character_; best_gain <- -Inf
    for (nm in names(sets)) {
      if (nm %in% chosen) next
      gain_set <- setdiff(sets[[nm]], covered)
      gain <- sum(w[intersect(names(w), gain_set)], na.rm = TRUE)
      if (gain > best_gain + 1e-12) { best_gain <- gain; best <- nm }
    }
    if (is.na(best) || best_gain <= 0) break
    chosen  <- c(chosen, best)
    covered <- union(covered, sets[[best]])
    have    <- sum(w[intersect(names(w), covered)], na.rm = TRUE)
  }
  chosen
}