# Shared helper functions for ciToti paper scRNA-seq analyses

parse_cli_args <- function(defaults = list()) {
  args <- commandArgs(trailingOnly = TRUE)
  out <- defaults
  if (length(args) == 0) return(out)
  i <- 1
  while (i <= length(args)) {
    key <- args[[i]]
    if (!startsWith(key, '--')) stop('Arguments should be provided as --name value. Problematic argument: ', key)
    key <- sub('^--', '', key)
    if (i == length(args) || startsWith(args[[i + 1]], '--')) {
      out[[key]] <- TRUE
      i <- i + 1
    } else {
      out[[key]] <- args[[i + 1]]
      i <- i + 2
    }
  }
  out
}

arg_vec <- function(x, sep = ',') {
  if (is.null(x) || is.na(x) || x == '') return(character(0))
  trimws(unlist(strsplit(as.character(x), sep, fixed = TRUE)))
}

arg_int_vec <- function(x, sep = ',') {
  as.integer(arg_vec(x, sep = sep))
}

ensure_dir <- function(path) {
  if (!dir.exists(path)) dir.create(path, recursive = TRUE, showWarnings = FALSE)
  invisible(path)
}

message_section <- function(...) {
  msg <- paste0(...)
  cat('\n', paste(rep('=', nchar(msg) + 8), collapse = ''), '\n', sep = '')
  cat('=== ', msg, ' ===\n', sep = '')
  cat(paste(rep('=', nchar(msg) + 8), collapse = ''), '\n', sep = '')
}

require_packages <- function(packages) {
  missing <- packages[!vapply(packages, requireNamespace, logical(1), quietly = TRUE)]
  if (length(missing) > 0) {
    stop('Missing required R packages: ', paste(missing, collapse = ', '),
         '\nInstall them before running this script.')
  }
  invisible(lapply(packages, library, character.only = TRUE))
}

save_plot_pdf <- function(plot, filename, width = 6, height = 5, useDingbats = FALSE) {
  ensure_dir(dirname(filename))
  grDevices::pdf(filename, width = width, height = height, useDingbats = useDingbats)
  print(plot)
  grDevices::dev.off()
  invisible(filename)
}

load_seurat_rds <- function(path, project = NULL, default_assay = 'RNA') {
  obj <- readRDS(path)
  if (inherits(obj, 'Seurat')) return(obj)
  if (!is.list(obj) || is.null(obj$counts) || is.null(obj$meta.data)) {
    stop('Input RDS must be a Seurat object or a list with counts and meta.data: ', path)
  }
  assay_name <- default_assay
  if (!is.null(obj$assay) && is.character(obj$assay) && length(obj$assay) == 1) assay_name <- obj$assay
  Seurat::CreateSeuratObject(
    counts = obj$counts,
    meta.data = obj$meta.data,
    assay = assay_name,
    project = if (is.null(project)) tools::file_path_sans_ext(basename(path)) else project
  )
}

join_layers_if_needed <- function(obj, assay = 'RNA') {
  if ('JoinLayers' %in% getNamespaceExports('Seurat')) {
    obj <- tryCatch(
      Seurat::JoinLayers(obj, assay = assay),
      error = function(e) obj
    )
  }
  obj
}

get_assay_data_compat <- function(obj, assay = 'RNA', slot_or_layer = 'counts') {
  DefaultAssay(obj) <- assay
  mat <- tryCatch(
    Seurat::GetAssayData(obj, assay = assay, slot = slot_or_layer),
    error = function(e) NULL
  )
  if (is.null(mat)) {
    mat <- tryCatch(
      Seurat::GetAssayData(obj, assay = assay, layer = slot_or_layer),
      error = function(e) NULL
    )
  }
  if (is.null(mat)) stop('Could not retrieve assay data: assay=', assay, ', slot/layer=', slot_or_layer)
  mat
}

set_orig_ident <- function(obj, sample_name) {
  obj$orig.ident <- sample_name
  obj$dataset <- sample_name
  obj
}

assign_celltype_from_cluster_map <- function(obj,
                                             cluster_map,
                                             cluster_col = 'seurat_clusters',
                                             celltype_col = 'celltype') {
  if (!cluster_col %in% colnames(obj[[]])) stop('Missing cluster metadata column: ', cluster_col)
  Idents(obj) <- cluster_col
  map_use <- cluster_map[names(cluster_map) %in% levels(Idents(obj))]
  if (length(map_use) > 0) {
    obj <- do.call(Seurat::RenameIdents, c(list(object = obj), as.list(map_use)))
  }
  obj[[celltype_col]] <- as.character(Idents(obj))
  Idents(obj) <- celltype_col
  obj
}

run_harmony_pipeline <- function(obj,
                                 group.by.vars = 'orig.ident',
                                 dims = 1:20,
                                 resolution = 0.5,
                                 npcs = max(dims),
                                 assay = 'RNA') {
  DefaultAssay(obj) <- assay
  obj <- join_layers_if_needed(obj, assay = assay)
  obj <- Seurat::NormalizeData(obj)
  obj <- Seurat::FindVariableFeatures(obj)
  obj <- Seurat::ScaleData(obj)
  obj <- Seurat::RunPCA(obj, npcs = npcs, verbose = FALSE)
  obj <- harmony::RunHarmony(obj, group.by.vars = group.by.vars)
  obj <- Seurat::FindNeighbors(obj, reduction = 'harmony', dims = dims)
  obj <- Seurat::FindClusters(obj, resolution = resolution)
  obj$harmony_clusters <- as.character(obj$seurat_clusters)
  obj <- Seurat::RunUMAP(obj, reduction = 'harmony', dims = dims)
  obj
}

add_sample_cluster <- function(obj, sample_col = 'orig.ident', output_col = 'sample_cluster') {
  obj[[output_col]] <- paste(as.character(Idents(obj)), obj[[sample_col]][, 1], sep = '_')
  obj
}

export_umap_source <- function(obj,
                               out_csv,
                               reduction = 'umap',
                               sample_col = 'orig.ident',
                               celltype_col = NULL) {
  emb <- as.data.frame(Seurat::Embeddings(obj, reduction = reduction))
  emb$cell_barcode <- rownames(emb)
  emb$sample <- obj[[sample_col]][rownames(emb), 1]
  emb$identity <- as.character(Idents(obj))[match(rownames(emb), names(Idents(obj)))]
  if (!is.null(celltype_col) && celltype_col %in% colnames(obj[[]])) emb$celltype <- obj[[celltype_col]][rownames(emb), 1]
  utils::write.csv(emb, out_csv, row.names = FALSE)
  invisible(emb)
}

save_dimplot <- function(obj,
                         out_pdf,
                         reduction = 'umap',
                         split.by = NULL,
                         group.by = NULL,
                         colors = NULL,
                         label = FALSE,
                         pt.size = 0.25,
                         width = 6,
                         height = 5) {
  p <- Seurat::DimPlot(obj, reduction = reduction, split.by = split.by, group.by = group.by,
                       label = label, pt.size = pt.size)
  if (!is.null(colors)) p <- p + ggplot2::scale_color_manual(values = colors, na.value = 'grey80')
  p <- p + ggplot2::theme_classic(base_size = 12)
  save_plot_pdf(p, out_pdf, width = width, height = height)
  invisible(p)
}

save_dotplot <- function(obj,
                         features,
                         out_pdf,
                         out_csv,
                         colors = c('steelblue', 'white', 'darkred'),
                         limits = NULL,
                         width = 8,
                         height = 4.5) {
  features_present <- features[features %in% rownames(obj)]
  if (length(features_present) == 0) stop('None of the requested dot-plot genes were found in the object.')
  missing <- setdiff(features, features_present)
  if (length(missing) > 0) warning('Missing dot-plot genes: ', paste(missing, collapse = ', '))
  p <- Seurat::DotPlot(obj, features = features_present)
  if (is.null(limits)) {
    p <- p + ggplot2::scale_colour_gradientn(colors = colors)
  } else {
    p <- p + ggplot2::scale_colour_gradientn(colors = colors, limits = limits)
  }
  p <- p + ggplot2::theme_classic(base_size = 11) +
    ggplot2::theme(axis.text.x = ggplot2::element_text(angle = 45, hjust = 1)) +
    ggplot2::guides(color = ggplot2::guide_colourbar(title = 'Expression'))
  utils::write.csv(p$data, out_csv, row.names = FALSE)
  save_plot_pdf(p, out_pdf, width = width, height = height)
  invisible(p)
}

save_average_expression_heatmap <- function(obj,
                                            features,
                                            group.by,
                                            out_pdf,
                                            out_csv,
                                            assay = 'RNA',
                                            cluster_rows = FALSE,
                                            cluster_cols = FALSE,
                                            width = 5,
                                            height = 6) {
  obj <- join_layers_if_needed(obj, assay = assay)
  features_present <- features[features %in% rownames(obj)]
  if (length(features_present) == 0) stop('None of the requested genes were found for AverageExpression.')
  missing <- setdiff(features, features_present)
  if (length(missing) > 0) warning('Missing heatmap genes: ', paste(missing, collapse = ', '))
  avg <- Seurat::AverageExpression(obj, group.by = group.by, features = features_present, assays = assay)
  mat <- avg[[assay]][features_present, , drop = FALSE]
  utils::write.csv(mat, out_csv)
  grDevices::pdf(out_pdf, width = width, height = height, useDingbats = FALSE)
  pheatmap::pheatmap(
    mat,
    color = grDevices::colorRampPalette(c('steelblue', 'white', 'darkred'))(100),
    cluster_rows = cluster_rows,
    cluster_cols = cluster_cols,
    show_rownames = TRUE,
    show_colnames = TRUE,
    main = 'Average expression'
  )
  grDevices::dev.off()
  invisible(mat)
}

compute_fraction_matrix <- function(obj, sample_col = 'orig.ident') {
  tab <- table(as.character(Idents(obj)), as.character(obj[[sample_col]][, 1]))
  prop.table(tab, margin = 2)
}

save_fraction_correlation_heatmap <- function(obj,
                                              out_pdf,
                                              out_matrix_csv,
                                              out_fraction_csv,
                                              sample_col = 'orig.ident',
                                              width = 4,
                                              height = 4) {
  prop <- compute_fraction_matrix(obj, sample_col = sample_col)
  utils::write.csv(as.matrix(prop), out_fraction_csv)
  correlation <- stats::cor(as.data.frame.matrix(prop))
  utils::write.csv(correlation, out_matrix_csv)
  grDevices::pdf(out_pdf, width = width, height = height, useDingbats = FALSE)
  pheatmap::pheatmap(
    correlation,
    color = grDevices::colorRampPalette(c('#006695', 'white', '#950000'))(50),
    breaks = seq(-1, 1, length.out = 51),
    border_color = NA,
    clustering_method = 'complete',
    display_numbers = FALSE,
    main = 'Pearson correlation'
  )
  grDevices::dev.off()
  invisible(correlation)
}

save_fraction_area_plot <- function(obj,
                                    out_pdf,
                                    out_csv,
                                    sample_col = 'orig.ident',
                                    sample_levels = NULL,
                                    width = 7.5,
                                    height = 5) {
  prop <- as.data.frame(compute_fraction_matrix(obj, sample_col = sample_col))
  colnames(prop) <- c('celltype', 'sample', 'fraction')
  if (!is.null(sample_levels)) prop$sample <- factor(prop$sample, levels = sample_levels)
  prop$x <- as.integer(prop$sample)
  utils::write.csv(prop, out_csv, row.names = FALSE)
  p <- ggplot2::ggplot(prop, ggplot2::aes(x = x, y = fraction, fill = celltype, group = celltype)) +
    ggplot2::geom_area(position = 'fill', alpha = 0.95) +
    ggplot2::scale_x_continuous(breaks = seq_along(levels(prop$sample)), labels = levels(prop$sample)) +
    ggplot2::labs(x = NULL, y = 'Fraction', fill = NULL) +
    ggplot2::theme_classic(base_size = 12) +
    ggplot2::theme(axis.text.x = ggplot2::element_text(angle = 45, hjust = 1))
  save_plot_pdf(p, out_pdf, width = width, height = height)
  invisible(prop)
}

scale_rows <- function(mat) {
  out <- t(scale(t(as.matrix(mat))))
  out[is.na(out)] <- 0
  out
}

get_top_conserved_genes <- function(marker_df, top_n = 10) {
  if (is.null(marker_df) || nrow(marker_df) == 0) return(character(0))
  logfc_cols <- grep('avg_log2FC|avg_logFC', colnames(marker_df), value = TRUE)
  if (length(logfc_cols) == 0) stop('Cannot find avg_log2FC/avg_logFC columns in marker table.')
  marker_df$mean_avg_log2FC <- rowMeans(marker_df[, logfc_cols, drop = FALSE], na.rm = TRUE)
  rownames(marker_df)[order(marker_df$mean_avg_log2FC, decreasing = TRUE)][seq_len(min(top_n, nrow(marker_df)))]
}
