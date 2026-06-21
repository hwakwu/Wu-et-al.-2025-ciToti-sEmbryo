# Monocle2 trajectory helper functions

patch_monocle2_project2MST <- function() {
  ns <- asNamespace('monocle')
  f <- get('project2MST', envir = ns)
  txt <- paste(deparse(body(f)), collapse = '\n')
  if (grepl('class\\(projection\\)', txt)) {
    txt2 <- gsub(
      'if \\(class\\(projection\\) != "matrix"\\)\\s*projection <- as.matrix\\(projection\\)',
      'projection <- as.matrix(projection)',
      txt
    )
    if (!identical(txt, txt2)) {
      body(f) <- parse(text = txt2)[[1]]
      assignInNamespace('project2MST', f, ns = 'monocle')
      message('Patched monocle::project2MST() for recent Matrix/class behavior.')
    }
  }
  invisible(TRUE)
}

downsample_cells_by_group <- function(meta,
                                      group_col = 'orig.ident',
                                      total_n = 8000,
                                      seed = 2) {
  if (!group_col %in% colnames(meta)) stop('Missing downsampling group column: ', group_col)
  meta$cell_id <- rownames(meta)
  set.seed(seed)
  prop_take <- min(1, total_n / nrow(meta))
  cells <- meta |>
    dplyr::group_by(.data[[group_col]]) |>
    dplyr::group_modify(~ dplyr::slice_sample(.x, prop = prop_take)) |>
    dplyr::ungroup() |>
    dplyr::pull(cell_id)
  cells
}

build_monocle2_cds <- function(obj,
                               cells,
                               assay = 'RNA',
                               lower_detection_limit = 0.5) {
  obj <- join_layers_if_needed(obj, assay = assay)
  counts <- get_assay_data_compat(obj, assay = assay, slot_or_layer = 'counts')[, cells, drop = FALSE]
  meta <- obj[[]][cells, , drop = FALSE]
  genes <- data.frame(gene_short_name = rownames(counts), row.names = rownames(counts), stringsAsFactors = FALSE)
  pd <- new('AnnotatedDataFrame', data = meta)
  fd <- new('AnnotatedDataFrame', data = genes)
  cds <- monocle::newCellDataSet(
    counts,
    phenoData = pd,
    featureData = fd,
    lowerDetectionLimit = lower_detection_limit,
    expressionFamily = monocle::negbinomial.size()
  )
  cds <- monocle::estimateSizeFactors(cds)
  cds <- monocle::estimateDispersions(cds)
  cds
}

run_monocle2_ddrtree <- function(obj,
                                 group_col = 'orig.ident',
                                 assay = 'RNA',
                                 total_n = 8000,
                                 n_ordering_genes = 200,
                                 seed = 2,
                                 cores = 1,
                                 reverse_pseudotime = TRUE,
                                 lower_detection_limit = 0.5) {
  patch_monocle2_project2MST()
  meta <- obj[[]]
  cells <- downsample_cells_by_group(meta, group_col = group_col, total_n = total_n, seed = seed)
  message('Monocle2 downsampling: ', length(cells), ' cells selected from ', nrow(meta), ' total cells.')
  message('Cells per group after downsampling:')
  print(table(meta[cells, group_col, drop = TRUE]))

  cds <- build_monocle2_cds(
    obj = obj,
    cells = cells,
    assay = assay,
    lower_detection_limit = lower_detection_limit
  )

  message('Testing genes by ', group_col, ' for Monocle2 ordering...')
  deg <- monocle::differentialGeneTest(cds, fullModelFormulaStr = paste0('~', group_col), cores = cores)
  deg <- deg[order(deg$qval, decreasing = FALSE), , drop = FALSE]
  ordering_genes <- rownames(deg)[seq_len(min(n_ordering_genes, nrow(deg)))]

  cds <- monocle::setOrderingFilter(cds, ordering_genes = ordering_genes)
  cds <- monocle::reduceDimension(cds, max_components = 2, method = 'DDRTree')
  cds <- monocle::orderCells(cds)

  if (reverse_pseudotime) {
    Biobase::pData(cds)$Pseudotime <- max(Biobase::pData(cds)$Pseudotime, na.rm = TRUE) - Biobase::pData(cds)$Pseudotime
  }

  list(cds = cds, deg = deg, ordering_genes = ordering_genes, cells = cells)
}

save_monocle2_source_data <- function(result, outdir, prefix) {
  ensure_dir(outdir)
  cds <- result$cds
  pd <- Biobase::pData(cds)
  pd$cell_barcode <- rownames(pd)
  utils::write.csv(pd, file.path(outdir, paste0(prefix, '_pseudotime_metadata.csv')), row.names = FALSE)
  utils::write.csv(result$deg, file.path(outdir, paste0(prefix, '_differential_gene_test_all_genes.csv')))
  utils::write.csv(data.frame(ordering_gene = result$ordering_genes),
                   file.path(outdir, paste0(prefix, '_top_ordering_genes.csv')), row.names = FALSE)
  saveRDS(cds, file.path(outdir, paste0(prefix, '_monocle2_cds.rds')))
  invisible(pd)
}

save_monocle2_trajectory_plot <- function(cds,
                                           color_by,
                                           out_pdf,
                                           colors = NULL,
                                           width = 4.8,
                                           height = 3.1,
                                           cell_size = 0.5) {
  p <- monocle::plot_cell_trajectory(cds, color_by = color_by, cell_size = cell_size) +
    ggplot2::theme_classic(base_size = 11) +
    ggplot2::theme(legend.position = 'right') +
    ggplot2::guides(color = ggplot2::guide_legend(override.aes = list(size = 4)))
  if (!is.null(colors)) p <- p + ggplot2::scale_color_manual(values = colors, na.value = 'grey80')
  save_plot_pdf(p, out_pdf, width = width, height = height)
  invisible(p)
}

save_pseudotime_density_plot <- function(cds,
                                         color_col,
                                         out_pdf,
                                         colors = NULL,
                                         exclude = NULL,
                                         facet = FALSE,
                                         width = 7.3,
                                         height = 5.3,
                                         bw = 0.8) {
  df <- as.data.frame(Biobase::pData(cds))
  if (!is.null(exclude)) df <- df[!as.character(df[[color_col]]) %in% exclude, , drop = FALSE]
  p <- ggplot2::ggplot(df, ggplot2::aes(x = Pseudotime, colour = .data[[color_col]], fill = .data[[color_col]])) +
    ggplot2::geom_density(bw = bw, linewidth = 1, alpha = 0.6) +
    ggplot2::theme_classic(base_size = 14) +
    ggplot2::labs(color = NULL, fill = NULL)
  if (!is.null(colors)) {
    p <- p + ggplot2::scale_color_manual(values = colors, na.value = 'grey80') +
      ggplot2::scale_fill_manual(values = colors, na.value = 'grey80')
  }
  if (facet) p <- p + ggplot2::facet_wrap(stats::as.formula(paste('~', color_col)), ncol = 1, scales = 'free_y')
  save_plot_pdf(p, out_pdf, width = width, height = height)
  invisible(p)
}

save_ddrtree_stage_pseudotime_panels <- function(cds,
                                                 group_col,
                                                 out_pdf,
                                                 group_colors = NULL,
                                                 width = 9,
                                                 height = 4) {
  p1 <- monocle::plot_cell_trajectory(cds, color_by = group_col, cell_size = 0.5) +
    ggplot2::theme_classic(base_size = 11) + ggplot2::theme(legend.position = 'right')
  if (!is.null(group_colors)) p1 <- p1 + ggplot2::scale_color_manual(values = group_colors, na.value = 'grey80')
  p2 <- monocle::plot_cell_trajectory(cds, color_by = 'Pseudotime', cell_size = 0.5) +
    ggplot2::theme_classic(base_size = 11) + ggplot2::theme(legend.position = 'right')
  p <- p1 + p2
  save_plot_pdf(p, out_pdf, width = width, height = height)
  invisible(p)
}
