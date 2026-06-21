#!/usr/bin/env Rscript

# Fig. 1m: Monocle2 trajectory analysis for ciToti1-ciToti4.
# The implementation follows the original workflow:
#   1. Load ciToti1-4 object.
#   2. Infer ciToti stage and cell type from sample_cluster if needed.
#   3. Proportionally downsample cells by orig.ident.
#   4. Build a Monocle2 CellDataSet using negbinomial.size().
#   5. Use differentialGeneTest(~orig.ident) and the top 200 q-value genes as ordering genes.
#   6. Run DDRTree, order cells, reverse pseudotime, and export source data.

source('R/utils_common.R')
source('R/utils_monocle2.R')
require_packages(c('Seurat', 'Matrix', 'dplyr', 'ggplot2', 'patchwork', 'monocle', 'Biobase'))

args <- parse_cli_args(list(
  input = 'data/ciToti1_4_v4.rds',
  outdir = 'results/Fig1m',
  total_n = '8000',
  n_ordering_genes = '200',
  seed = '2',
  assay = 'RNA'
))
ensure_dir(args$outdir)

message_section('Loading ciToti1-4 object')
obj <- load_seurat_rds(args$input, project = 'ciToti1_4')
DefaultAssay(obj) <- args$assay
obj <- join_layers_if_needed(obj, assay = args$assay)

# Some exported objects contain labels like "Intermediate_ciToti1" in sample_cluster.
# Recover orig.ident and celltype when these metadata columns are not already present.
if ('sample_cluster' %in% colnames(obj[[]])) {
  sc <- as.character(obj$sample_cluster)
  stage <- regmatches(sc, regexpr('ciToti[1-4]', sc))
  stage[stage == ''] <- NA_character_
  if (!'orig.ident' %in% colnames(obj[[]]) || length(unique(as.character(obj$orig.ident))) <= 1) {
    obj$orig.ident <- ifelse(is.na(stage), as.character(obj$orig.ident), stage)
  }
  if (!'celltype' %in% colnames(obj[[]])) {
    obj$celltype <- sub('_?ciToti[1-4]$', '', sc)
  }
}

if (!'orig.ident' %in% colnames(obj[[]])) stop('Metadata column orig.ident is required.')
if (!'celltype' %in% colnames(obj[[]])) stop('Metadata column celltype is required for Fig. 1m density plots.')
obj$orig.ident <- factor(as.character(obj$orig.ident), levels = c('ciToti1', 'ciToti2', 'ciToti3', 'ciToti4'))
obj$dataset <- obj$orig.ident

stage_colors <- c('ciToti1' = '#FFF3F3', 'ciToti2' = '#B7A9A9', 'ciToti3' = '#B16E6E', 'ciToti4' = '#870000')
celltype_colors <- c('Totipotent' = '#D9D9D9', 'Intermediate' = '#B7A9A9', 'TE' = '#7b1fa2', 'EPI' = '#B16E6E', 'Epiblast' = '#B16E6E', 'PrE' = '#870000')
celltype_colors <- celltype_colors[names(celltype_colors) %in% unique(as.character(obj$celltype))]

message_section('Running Monocle2 DDRTree trajectory')
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

save_monocle2_source_data(traj, args$outdir, 'Fig1m_ciToti1_4')
save_monocle2_trajectory_plot(traj$cds, 'orig.ident', file.path(args$outdir, 'Fig1m_ciToti1_4_trajectory_by_stage.pdf'), colors = stage_colors)
save_pseudotime_density_plot(traj$cds, 'celltype', file.path(args$outdir, 'Fig1m_ciToti1_4_density_by_celltype.pdf'), colors = celltype_colors, exclude = 'Totipotent')
save_pseudotime_density_plot(traj$cds, 'celltype', file.path(args$outdir, 'Fig1m_ciToti1_4_density_by_celltype_facet.pdf'), colors = celltype_colors, exclude = 'Totipotent', facet = TRUE, height = 6)
save_ddrtree_stage_pseudotime_panels(traj$cds, 'orig.ident', file.path(args$outdir, 'Fig1m_ciToti1_4_DDRTree_stage_pseudotime.pdf'), group_colors = stage_colors)

message('Done. Fig. 1m outputs saved to ', args$outdir)
