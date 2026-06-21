#!/usr/bin/env Rscript

# Fig. 3k,l: Integrated ciToti4-ciToti10 atlas and cell-type fraction plot.

source('R/utils_common.R')
source('R/celltype_annotations.R')
require_packages(c('Seurat', 'dplyr', 'ggplot2', 'patchwork', 'harmony'))

args <- parse_cli_args(list(
  ciToti4 = 'data/ciToti4_subset.rds',
  ciToti7 = 'data/ciToti7_subset.rds',
  ciToti8 = 'data/ciToti8_subset.rds',
  ciToti10 = 'data/ciToti10_subset.rds',
  outdir = 'results/Fig3kl',
  dims = '1:20',
  resolution = '0.5',
  assay = 'RNA'
))
ensure_dir(args$outdir)

parse_dims <- function(x) {
  if (grepl(':', x, fixed = TRUE)) {
    parts <- as.integer(strsplit(x, ':', fixed = TRUE)[[1]])
    seq(parts[1], parts[2])
  } else {
    as.integer(arg_vec(x))
  }
}
dims <- parse_dims(args$dims)

message_section('Loading ciToti4/7/8/10 subset objects')
ciToti4 <- load_seurat_rds(args$ciToti4, project = 'ciToti4') |> set_orig_ident('ciToti4')
ciToti7 <- load_seurat_rds(args$ciToti7, project = 'ciToti7') |> set_orig_ident('ciToti7')
ciToti8 <- load_seurat_rds(args$ciToti8, project = 'ciToti8') |> set_orig_ident('ciToti8')
ciToti10 <- load_seurat_rds(args$ciToti10, project = 'ciToti10') |> set_orig_ident('ciToti10')

# If the objects are not already annotated, apply the manual cluster maps used elsewhere.
if (!'celltype' %in% colnames(ciToti4[[]]) && 'seurat_clusters' %in% colnames(ciToti4[[]])) ciToti4 <- assign_celltype_from_cluster_map(ciToti4, ciToti4_cluster_map)
if (!'celltype' %in% colnames(ciToti7[[]]) && 'seurat_clusters' %in% colnames(ciToti7[[]])) ciToti7 <- assign_celltype_from_cluster_map(ciToti7, ciToti7_cluster_map)
if (!'celltype' %in% colnames(ciToti8[[]]) && 'seurat_clusters' %in% colnames(ciToti8[[]])) ciToti8 <- assign_celltype_from_cluster_map(ciToti8, ciToti8_cluster_map)
if (!'celltype' %in% colnames(ciToti10[[]]) && 'seurat_clusters' %in% colnames(ciToti10[[]])) ciToti10 <- assign_celltype_from_cluster_map(ciToti10, ciToti10_cluster_map)

# Use common genes before merging, matching the original analysis code.
common_genes <- Reduce(intersect, list(rownames(ciToti4), rownames(ciToti7), rownames(ciToti8), rownames(ciToti10)))
ciToti4_common <- ciToti4[common_genes, ]
ciToti7_common <- ciToti7[common_genes, ]
ciToti8_common <- ciToti8[common_genes, ]
ciToti10_common <- ciToti10[common_genes, ]

combined4 <- merge(
  ciToti4_common,
  y = list(ciToti7_common, ciToti8_common, ciToti10_common),
  add.cell.ids = c('ciToti4', 'ciToti7', 'ciToti8', 'ciToti10'),
  project = 'Fig3_ciToti4_10_atlas'
)
combined4$orig.ident <- factor(as.character(combined4$orig.ident), levels = c('ciToti4', 'ciToti7', 'ciToti8', 'ciToti10'))

message_section('Running Harmony integration for ciToti4-10 atlas')
combined4 <- run_harmony_pipeline(
  combined4,
  group.by.vars = 'orig.ident',
  dims = dims,
  resolution = as.numeric(args$resolution),
  assay = args$assay
)

# Prefer existing celltype metadata for Fig. 3k/l annotation. If absent after merge, use active identities.
if ('celltype' %in% colnames(combined4[[]])) {
  Idents(combined4) <- 'celltype'
}

message_section('Saving Fig. 3k atlas UMAP')
save_dimplot(
  combined4,
  out_pdf = file.path(args$outdir, 'Fig3k_ciToti4_10_integrated_UMAP_split_by_stage.pdf'),
  split.by = 'orig.ident',
  colors = fig3_atlas_colors,
  pt.size = 0.25,
  width = 12,
  height = 4
)
export_umap_source(combined4, file.path(args$outdir, 'Fig3k_ciToti4_10_integrated_UMAP_source_data.csv'), celltype_col = if ('celltype' %in% colnames(combined4[[]])) 'celltype' else NULL)

message_section('Saving Fig. 3l cell-type fraction plot')
save_fraction_area_plot(
  combined4,
  out_pdf = file.path(args$outdir, 'Fig3l_ciToti4_10_celltype_fraction_area.pdf'),
  out_csv = file.path(args$outdir, 'Fig3l_ciToti4_10_celltype_fraction_source_data.csv'),
  sample_col = 'orig.ident',
  sample_levels = c('ciToti4', 'ciToti7', 'ciToti8', 'ciToti10'),
  width = 7.5,
  height = 5
)

saveRDS(combined4, file.path(args$outdir, 'Fig3kl_ciToti4_10_integrated_annotated_seurat.rds'))
message('Done. Fig. 3k,l outputs saved to ', args$outdir)
