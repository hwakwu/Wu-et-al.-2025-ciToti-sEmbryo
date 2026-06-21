#!/usr/bin/env Rscript

# Fig. 2c: Monocle2 trajectory analysis for ciToti4, ciToti7, and ciToti8.

source('R/utils_common.R')
source('R/utils_monocle2.R')
source('R/celltype_annotations.R')
require_packages(c('Seurat', 'Matrix', 'dplyr', 'ggplot2', 'patchwork', 'monocle', 'Biobase'))

args <- parse_cli_args(list(
  ciToti4 = 'data/ciToti4_v4.rds',
  ciToti7 = 'data/ciToti7_v4.rds',
  ciToti8 = 'data/ciToti8_v4.rds',
  outdir = 'results/Fig2c',
  total_n = '8000',
  n_ordering_genes = '200',
  seed = '2',
  assay = 'RNA'
))
ensure_dir(args$outdir)

message_section('Loading and annotating ciToti4/7/8')
ciToti4 <- load_seurat_rds(args$ciToti4, project = 'ciToti4') |> set_orig_ident('ciToti4') |> assign_celltype_from_cluster_map(ciToti4_cluster_map)
ciToti7 <- load_seurat_rds(args$ciToti7, project = 'ciToti7') |> set_orig_ident('ciToti7') |> assign_celltype_from_cluster_map(ciToti7_cluster_map)
ciToti8 <- load_seurat_rds(args$ciToti8, project = 'ciToti8') |> set_orig_ident('ciToti8') |> assign_celltype_from_cluster_map(ciToti8_cluster_map)

merged_obj <- merge(
  x = ciToti4,
  y = list(ciToti7, ciToti8),
  add.cell.ids = c('ciToti4', 'ciToti7', 'ciToti8'),
  project = 'ciToti4_8'
)
merged_obj$orig.ident <- factor(as.character(merged_obj$orig.ident), levels = c('ciToti4', 'ciToti7', 'ciToti8'))
merged_obj$dataset <- merged_obj$orig.ident

stage_colors <- c('ciToti4' = '#B7A9A9', 'ciToti7' = '#B16E6E', 'ciToti8' = '#870000')

message_section('Running Monocle2 DDRTree trajectory for ciToti4-8')
traj <- run_monocle2_ddrtree(
  obj = merged_obj,
  group_col = 'orig.ident',
  assay = args$assay,
  total_n = as.integer(args$total_n),
  n_ordering_genes = as.integer(args$n_ordering_genes),
  seed = as.integer(args$seed),
  cores = 1,
  reverse_pseudotime = TRUE
)

save_monocle2_source_data(traj, args$outdir, 'Fig2c_ciToti4_8')
save_monocle2_trajectory_plot(traj$cds, 'orig.ident', file.path(args$outdir, 'Fig2c_ciToti4_8_trajectory_by_stage.pdf'), colors = stage_colors)
save_pseudotime_density_plot(traj$cds, 'orig.ident', file.path(args$outdir, 'Fig2c_ciToti4_8_density_by_stage.pdf'), colors = stage_colors)
save_pseudotime_density_plot(traj$cds, 'orig.ident', file.path(args$outdir, 'Fig2c_ciToti4_8_density_by_stage_facet.pdf'), colors = stage_colors, facet = TRUE, height = 5.5)
save_ddrtree_stage_pseudotime_panels(traj$cds, 'orig.ident', file.path(args$outdir, 'Fig2c_ciToti4_8_DDRTree_stage_pseudotime.pdf'), group_colors = stage_colors)
saveRDS(merged_obj, file.path(args$outdir, 'Fig2c_ciToti4_8_merged_annotated_seurat.rds'))

message('Done. Fig. 2c outputs saved to ', args$outdir)
