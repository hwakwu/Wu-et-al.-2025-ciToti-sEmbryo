#!/usr/bin/env Rscript

# Preprocess raw 10X scRNA-seq folders into filtered Seurat objects.
# This corresponds to the initial sample-processing section of the analysis code.

source('R/utils_common.R')
require_packages(c('Seurat', 'ggplot2'))

args <- parse_cli_args(list(
  manifest = 'data/sample_manifest.csv',
  outdir = 'data/processed_seurat',
  mt_pattern = '^mt-|^MT-'
))
ensure_dir(args$outdir)

manifest <- utils::read.csv(args$manifest, stringsAsFactors = FALSE, check.names = FALSE)
required_cols <- c('sample_name', 'data_path')
missing_cols <- setdiff(required_cols, colnames(manifest))
if (length(missing_cols) > 0) stop('Manifest missing required columns: ', paste(missing_cols, collapse = ', '))

process_sample <- function(sample_name,
                           data_path,
                           min_cells = 3,
                           nFeature_min = 500,
                           nFeature_max = 7000,
                           nCount_max = 40000,
                           percent_mt_max = 0.05,
                           mt_pattern = '^mt-|^MT-') {
  data <- Seurat::Read10X(data.dir = data_path)
  obj <- Seurat::CreateSeuratObject(counts = data, min.cells = min_cells, project = sample_name)
  obj$orig.ident <- sample_name
  obj[['percent.mt']] <- Seurat::PercentageFeatureSet(obj, pattern = mt_pattern)

  qc_pdf <- file.path(args$outdir, paste0(sample_name, '_QC_violin.pdf'))
  qc_plot <- Seurat::VlnPlot(obj, features = c('nFeature_RNA', 'nCount_RNA', 'percent.mt'), ncol = 3)
  save_plot_pdf(qc_plot, qc_pdf, width = 8, height = 3)

  obj <- subset(
    obj,
    subset = nFeature_RNA > nFeature_min &
      nFeature_RNA < nFeature_max &
      nCount_RNA < nCount_max &
      percent.mt < percent_mt_max
  )
  obj
}

processed <- list()
for (i in seq_len(nrow(manifest))) {
  sample_name <- manifest$sample_name[[i]]
  message_section('Processing ', sample_name)
  min_cells <- if ('min_cells' %in% colnames(manifest) && !is.na(manifest$min_cells[[i]])) manifest$min_cells[[i]] else 3
  nFeature_min <- if ('nFeature_min' %in% colnames(manifest) && !is.na(manifest$nFeature_min[[i]])) manifest$nFeature_min[[i]] else 500
  nFeature_max <- if ('nFeature_max' %in% colnames(manifest) && !is.na(manifest$nFeature_max[[i]])) manifest$nFeature_max[[i]] else 7000
  nCount_max <- if ('nCount_max' %in% colnames(manifest) && !is.na(manifest$nCount_max[[i]])) manifest$nCount_max[[i]] else 40000
  percent_mt_max <- if ('percent_mt_max' %in% colnames(manifest) && !is.na(manifest$percent_mt_max[[i]])) manifest$percent_mt_max[[i]] else 0.05

  processed[[sample_name]] <- process_sample(
    sample_name = sample_name,
    data_path = manifest$data_path[[i]],
    min_cells = min_cells,
    nFeature_min = nFeature_min,
    nFeature_max = nFeature_max,
    nCount_max = nCount_max,
    percent_mt_max = percent_mt_max,
    mt_pattern = args$mt_pattern
  )
  saveRDS(processed[[sample_name]], file.path(args$outdir, paste0(sample_name, '.rds')))
}

saveRDS(processed, file.path(args$outdir, 'processed_sample_list.rds'))
message('Done. Seurat objects saved in ', args$outdir)
