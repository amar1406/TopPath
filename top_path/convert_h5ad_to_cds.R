

# Robust .h5ad -> Monocle3 CDS with auto species detect, symbols, and mouse->human mapping

suppressPackageStartupMessages({
  library(zellkonverter)
  library(SingleCellExperiment)
  library(SummarizedExperiment)
  library(Matrix)
  library(monocle3)
  library(AnnotationDbi)
})

# --- config: where to cache/read the ortholog map
.ORTHO_CSV <- "database/mouse_to_human_orthologs.csv"

# --- helpers -----------------------------------------------------------------
.guess_species_from_ids <- function(ids) {
  x <- ids[which(!is.na(ids))[1]]
  if (is.na(x)) return(NA_character_)
  if (grepl("^ENSMUSG", x)) return("mouse")
  if (grepl("^ENSG",    x)) return("human")
  NA_character_
}

.map_ensembl_to_symbol <- function(ens_ids, species=c("auto","human","mouse")) {
  species <- match.arg(species)
  if (species == "auto") species <- .guess_species_from_ids(ens_ids)
  if (is.na(species)) return(ens_ids)
  
  if (species == "human" && requireNamespace("org.Hs.eg.db", quietly=TRUE)) {
    sym <- AnnotationDbi::mapIds(org.Hs.eg.db::org.Hs.eg.db, keys=ens_ids,
                                 keytype="ENSEMBL", column="SYMBOL", multiVals="first")
    return(ifelse(is.na(sym) | sym=="", ens_ids, unname(sym)))
  }
  if (species == "mouse" && requireNamespace("org.Mm.eg.db", quietly=TRUE)) {
    sym <- AnnotationDbi::mapIds(org.Mm.eg.db::org.Mm.eg.db, keys=ens_ids,
                                 keytype="ENSEMBL", column="SYMBOL", multiVals="first")
    return(ifelse(is.na(sym) | sym=="", ens_ids, unname(sym)))
  }
  ens_ids
}

.load_or_build_ortholog_map <- function(path=.ORTHO_CSV) {
  if (file.exists(path)) {
    df <- tryCatch(read.csv(path, stringsAsFactors=FALSE), error=function(e) NULL)
    if (!is.null(df)) return(df)
  }
  # Try to build once via biomaRt (best effort; if fails, we return NULL)
  if (!requireNamespace("biomaRt", quietly=TRUE)) return(NULL)
  df <- try({
    mm <- biomaRt::useEnsembl(biomart="genes", dataset="mmusculus_gene_ensembl")
    orth <- biomaRt::getBM(
      attributes=c("ensembl_gene_id","mgi_symbol",
                   "hsapiens_homolog_ensembl_gene",
                   "hsapiens_homolog_associated_gene_name",
                   "hsapiens_homolog_orthology_type"),
      mart=mm
    )
    names(orth) <- c("mouse_ensembl","mouse_symbol","human_ensembl","human_symbol","orthology_type")
    orth <- orth[orth$human_symbol!="" & orth$mouse_symbol!="", ]
    dir.create(dirname(path), recursive=TRUE, showWarnings=FALSE)
    write.csv(orth, path, row.names=FALSE)
    orth
  }, silent=TRUE)
  if (inherits(df, "try-error")) return(NULL)
  df
}

.mouse_symbols_to_human <- function(symbols, orth_map) {
  if (is.null(orth_map) || nrow(orth_map)==0) return(toupper(symbols))  # fallback
  df <- orth_map
  # prefer one2one when available
  if ("orthology_type" %in% names(df)) {
    ord <- match(df$orthology_type, c("ortholog_one2one","ortholog_one2many","ortholog_many2one"))
    ord[is.na(ord)] <- 99
    df <- df[order(ord), ]
    df <- df[!duplicated(toupper(df$mouse_symbol)), ]
  }
  mouse_uc <- toupper(df$mouse_symbol)
  human_uc <- toupper(df$human_symbol)
  idx <- match(toupper(symbols), mouse_uc)
  out <- human_uc[idx]
  ifelse(is.na(out) | out=="", toupper(symbols), out)
}

.to_dgc <- function(m) { if (inherits(m,"dgCMatrix")) m else suppressWarnings(as(m,"dgCMatrix")) }

.pick_best_counts <- function(sce) {
  an <- assayNames(sce)
  
  # 1) Prefer a real counts assay in the main experiment
  if ("counts" %in% an) {
    return(list(mat = assay(sce, "counts"), src = "counts"))
  }
  
  # 2) Check alternate experiments safely
  ae_names <- try(altExpNames(sce), silent = TRUE)
  if (!inherits(ae_names, "try-error") && length(ae_names) > 0) {
    for (nm in ae_names) {
      ae <- altExp(sce, nm)
      ae_assays <- assayNames(ae)
      if ("counts" %in% ae_assays) {
        return(list(mat = assay(ae, "counts"), src = paste0("altExp:", nm, "/counts")))
      }
      if ("X" %in% ae_assays) {
        return(list(mat = assay(ae, "X"), src = paste0("altExp:", nm, "/X")))
      }
    }
  }
  
  # 3) Fallback to X in the main experiment
  if ("X" %in% an) {
    return(list(mat = assay(sce, "X"), src = "X"))
  }
  
  stop("No suitable assay found (counts or X).")
}


.clean_to_integer_counts <- function(mat) {
  mat[is.na(mat)] <- 0
  mat[mat < 0] <- 0
  rounded <- FALSE
  if (any(abs(mat - round(mat)) > 1e-8)) { mat <- round(mat); rounded <- TRUE }
  list(mat=.to_dgc(mat), rounded=rounded)
}

.aggregate_duplicate_genes <- function(cnt, gene_ids, gene_symbols) {
  if (!any(duplicated(gene_ids))) return(list(counts=cnt, ids=gene_ids, symbols=gene_symbols))
  agg <- rowsum(cnt, group=gene_ids, reorder=FALSE)
  new_ids <- rownames(agg)
  sym_map <- tapply(gene_symbols, INDEX=gene_ids, FUN=function(v) v[which(!is.na(v) & v!="")[1]])
  sym_map[is.na(sym_map) | sym_map==""] <- new_ids[is.na(sym_map) | sym_map==""]
  sym_map <- sym_map[new_ids]
  list(counts=.to_dgc(agg), ids=new_ids, symbols=unname(sym_map))
}

# --- main -------------------------------------------------------------------
convert_h5ad_to_cds <- function(h5ad_path, use_disk_backed=TRUE) {
  sce <- zellkonverter::readH5AD(h5ad_path, use_hdf5=use_disk_backed)
  
  # pick counts-like matrix
  pick <- .pick_best_counts(sce)
  cnt  <- pick$mat
  i_src <- pick$src
  
  clean <- .clean_to_integer_counts(cnt)
  cnt   <- clean$mat
  i_rounded <- clean$rounded
  
  var <- as.data.frame(rowData(sce))
  ens_ids <- rownames(sce)
  sp_guess <- .guess_species_from_ids(ens_ids)
  
  # choose symbol column from .var if present; else map from Ensembl
  symbol_candidates <- c("gene_symbol","gene_symbols","gene","symbol","Gene","SYMBOL","gene_short_name")
  hit <- intersect(symbol_candidates, colnames(var))
  if (length(hit)) {
    gene_symbols <- as.character(var[[hit[1]]])
  } else {
    gene_symbols <- .map_ensembl_to_symbol(ens_ids, species = "auto")
  }
  gene_symbols[is.na(gene_symbols) | gene_symbols==""] <- ens_ids
  
  # collapse duplicate Ensembl rows before downstream
  agg <- .aggregate_duplicate_genes(cnt, ens_ids, gene_symbols)
  cnt          <- agg$counts
  ens_ids      <- agg$ids
  gene_symbols <- agg$symbols
  
  # auto humanize if mouse
  orth_used <- FALSE
  orth_rows <- 0L
  if (identical(sp_guess, "mouse")) {
    orth <- .load_or_build_ortholog_map(.ORTHO_CSV)
    old <- gene_symbols
    gene_symbols <- .mouse_symbols_to_human(gene_symbols, orth)
    orth_used <- !is.null(orth)
    orth_rows <- if (!is.null(orth)) nrow(orth) else 0L
  }
  # always uppercase final symbols to align with LR DB
  gene_symbols <- toupper(gene_symbols)
  
  # cell metadata
  cell_md <- as.data.frame(colData(sce))
  if (anyDuplicated(rownames(cell_md))) {
    rownames(cell_md) <- make.unique(rownames(cell_md), sep="_")
    colnames(cnt)     <- rownames(cell_md)
  }
  
  gene_md <- data.frame(
    id = ens_ids,
    gene_short_name = gene_symbols,
    row.names = ens_ids,
    stringsAsFactors = FALSE
  )
  
  cds <- monocle3::new_cell_data_set(
    expression_data = cnt,
    cell_metadata   = cell_md,
    gene_metadata   = gene_md
  )
  
  attr(cds, "import_info") <- list(
    source_assay   = i_src,
    rounded_counts = i_rounded,
    species_guess  = sp_guess,
    symbol_source  = if (length(hit)) paste0("var$", hit[1]) else "mapped_from_ensembl",
    ortholog_used  = orth_used,
    ortholog_rows  = orth_rows
  )
  cds
}


convert_seurat_to_cds <- function(seu, assay = NULL) {
  stopifnot(inherits(seu, "Seurat"))
  if (is.null(assay)) assay <- DefaultAssay(seu)
  
  # counts
  if (!assay %in% names(seu@assays)) stop("Assay '", assay, "' not found in Seurat object.")
  cnt <- tryCatch(seu@assays[[assay]]@counts, error = function(e) NULL)
  if (is.null(cnt) || length(cnt) == 0) {
    # fallback to data (lognorm) if counts missing; we’ll round to integer
    cnt <- seu@assays[[assay]]@data
  }
  cnt <- as(cnt, "dgCMatrix")
  cnt@x[cnt@x < 0] <- 0
  if (any(abs(cnt@x - round(cnt@x)) > 1e-8)) cnt@x <- round(cnt@x)
  
  # gene metadata
  ens_ids <- rownames(cnt)
  gene_symbols <- rownames(seu@assays[[assay]]@data)
  if (is.null(gene_symbols) || length(gene_symbols) != length(ens_ids)) {
    gene_symbols <- ens_ids
  }
  
  gene_md <- data.frame(
    id = ens_ids,
    gene_short_name = toupper(gene_symbols),
    row.names = ens_ids, stringsAsFactors = FALSE
  )
  
  # cell metadata
  cell_md <- seu@meta.data
  if (anyDuplicated(rownames(cell_md))) {
    rownames(cell_md) <- make.unique(rownames(cell_md), sep = "_")
  }
  # ensure columns match
  cnt <- cnt[, rownames(cell_md), drop = FALSE]
  
  cds <- monocle3::new_cell_data_set(
    expression_data = cnt,
    cell_metadata   = cell_md,
    gene_metadata   = gene_md
  )
  attr(cds, "import_info") <- list(
    source_assay   = paste0("Seurat:", assay, "/counts"),
    rounded_counts = TRUE,
    species_guess  = NA_character_,
    symbol_source  = "Seurat rownames",
    ortholog_used  = FALSE,
    ortholog_rows  = 0L
  )
  cds
}

