#!/usr/bin/env Rscript

# Fig. 2d: Monocle2 trajectory analysis for natural embryos from E4.5 to E7.5.
# The default input is the combined query object containing D_E45, D_E65, and D_E75 labels.

source('R/utils_common.R')
source('R/utils_monocle2.R')
require_packages(c('Seurat', 'Matrix', 'dplyr', 'ggplot2', 'patchwork', 'monocle', 'Biobase'))

args <- parse_cli_args(list(
  input = 'data/Combined_query_data.RDS',
  stages = 'D_E45,D_E65,D_E75',
  outdir = 'results/Fig2d',
  total_n = '8000',
  n_ordering_genes = '200',
  seed = '2',
  assay = 'RNA'
))
ensure_dir(args$outdir)

message_section('Loading natural embryo combined object')
obj <- load_seurat_rds(args$input, project = 'E45_E75')
obj <- join_layers_if_needed(obj, assay = args$assay)
stages <- arg_vec(args$stages)
missing <- setdiff(stages, unique(as.character(obj$orig.ident)))
if (length(missing) > 0) stop('Missing stages in orig.ident: ', paste(missing, collapse = ', '))
obj <- subset(obj, subset = orig.ident %in% stages)
obj$orig.ident <- factor(as.character(obj$orig.ident), levels = stages)
obj$dataset <- obj$orig.ident

stage_colors <- setNames(c('#B7A9A9', '#B16E6E', '#870000')[seq_along(stages)], stages)

message_section('Running Monocle2 DDRTree trajectory for natural E4.5-E7.5')
traj <- run_monocle2_ddrtree(
  obj = obj,
  group_col = 'orig.ident',
  assay = args$assay,
  total_n = as.integer(args$total_n),
  n_ordering_genes = as.integer(args$n_ordering_genes),
  seed = as.integer(args$seed),
  cores = 1,
  reverse_pseudotime = TRUE
)

save_monocle2_source_data(traj, args$outdir, 'Fig2d_E45_E75')
save_monocle2_trajectory_plot(traj$cds, 'orig.ident', file.path(args$outdir, 'Fig2d_E45_E75_trajectory_by_stage.pdf'), colors = stage_colors)
save_pseudotime_density_plot(traj$cds, 'orig.ident', file.path(args$outdir, 'Fig2d_E45_E75_density_by_stage.pdf'), colors = stage_colors)
save_pseudotime_density_plot(traj$cds, 'orig.ident', file.path(args$outdir, 'Fig2d_E45_E75_density_by_stage_facet.pdf'), colors = stage_colors, facet = TRUE, height = 5.5)
save_ddrtree_stage_pseudotime_panels(traj$cds, 'orig.ident', file.path(args$outdir, 'Fig2d_E45_E75_DDRTree_stage_pseudotime.pdf'), group_colors = stage_colors)
saveRDS(obj, file.path(args$outdir, 'Fig2d_E45_E75_subset_seurat.rds'))

message('Done. Fig. 2d outputs saved to ', args$outdir)
