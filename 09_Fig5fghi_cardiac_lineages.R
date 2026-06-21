#!/usr/bin/env Rscript

# Fig. 5f-i: Cardiac-lineage analysis comparing ciToti12-H and E10.5-H natural embryo cardiac tissues.

source('R/utils_common.R')
source('R/celltype_annotations.R')
require_packages(c('Seurat', 'dplyr', 'tidyr', 'ggplot2', 'patchwork', 'harmony', 'pheatmap', 'ComplexHeatmap', 'circlize', 'grid'))

args <- parse_cli_args(list(
  ciToti12_H = 'data/ciToti12_H.rds',
  E105_H = 'data/E10.5_H_in_utero.rds',
  outdir = 'results/Fig5',
  dims = '1:20',
  resolution = '0.2',
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

message_section('Loading ciToti12-H and E10.5-H objects')
ciToti12_H <- load_seurat_rds(args$ciToti12_H, project = 'ciToti12_H') |> set_orig_ident('ciToti12_H')
E105_H <- load_seurat_rds(args$E105_H, project = 'E10.5_H') |> set_orig_ident('E10.5_H')

combined6 <- merge(
  ciToti12_H,
  y = list(E105_H),
  add.cell.ids = c('ciToti12_H', 'E10.5_H'),
  project = 'Fig5_cardiac'
)

message_section('Running Harmony/UMAP for cardiac lineages')
combined6 <- run_harmony_pipeline(
  combined6,
  group.by.vars = 'orig.ident',
  dims = dims,
  resolution = as.numeric(args$resolution),
  assay = args$assay
)
combined6$cardiac_subcluster <- as.character(combined6$seurat_clusters)
combined6 <- assign_celltype_from_cluster_map(combined6, cardiac_cluster_map, cluster_col = 'seurat_clusters', celltype_col = 'celltype')
Idents(combined6) <- 'celltype'

message_section('Fig. 5f: UMAP colored by sample')
sample_colors <- c('ciToti12_H' = '#8A2BE2', 'E10.5_H' = 'lightgray')
save_dimplot(
  combined6,
  out_pdf = file.path(args$outdir, 'Fig5f_cardiac_UMAP_by_sample.pdf'),
  group.by = 'orig.ident',
  colors = sample_colors,
  pt.size = 0.4,
  width = 5.5,
  height = 5
)

message_section('Fig. 5g: UMAP split by sample and colored by cardiac cell type')
save_dimplot(
  combined6,
  out_pdf = file.path(args$outdir, 'Fig5g_cardiac_UMAP_split_by_sample_celltype.pdf'),
  split.by = 'orig.ident',
  colors = cardiac_colors,
  pt.size = 0.4,
  width = 9,
  height = 4.5
)
export_umap_source(combined6, file.path(args$outdir, 'Fig5fg_cardiac_UMAP_source_data.csv'), celltype_col = 'celltype')

message_section('Fig. 5h: cardiac-lineage fraction correlation excluding blood lineages')
cardiac_keep <- c('V-CMs', 'OFT', 'Artery-EC', 'Epicardium', 'Endocardium', 'A-CMs', 'CNCCs')
keep_cells <- names(Idents(combined6))[as.character(Idents(combined6)) %in% cardiac_keep]
CT_without_Blood <- subset(combined6, cells = keep_cells)
Idents(CT_without_Blood) <- 'celltype'
save_fraction_correlation_heatmap(
  CT_without_Blood,
  out_pdf = file.path(args$outdir, 'Fig5h_cardiac_without_blood_fraction_correlation_heatmap.pdf'),
  out_matrix_csv = file.path(args$outdir, 'Fig5h_cardiac_without_blood_fraction_correlation_matrix.csv'),
  out_fraction_csv = file.path(args$outdir, 'Fig5h_cardiac_without_blood_fraction_matrix.csv'),
  sample_col = 'orig.ident',
  width = 4,
  height = 4
)

message_section('Fig. 5i: cardiac marker dot heatmap')
combined6 <- add_sample_cluster(combined6, sample_col = 'orig.ident', output_col = 'sample_cluster')
Idents(combined6) <- 'sample_cluster'
features_present <- fig5i_cardiac_genes[fig5i_cardiac_genes %in% rownames(combined6)]
if (length(features_present) == 0) stop('None of the Fig. 5i genes were found in the object.')
missing_features <- setdiff(fig5i_cardiac_genes, features_present)
if (length(missing_features) > 0) warning('Missing Fig. 5i genes: ', paste(missing_features, collapse = ', '))

P <- Seurat::DotPlot(combined6, features = features_present) +
  ggplot2::scale_colour_gradientn(colors = c('steelblue', 'white', 'darkred')) +
  ggplot2::theme_classic(base_size = 11) +
  ggplot2::theme(axis.text.x = ggplot2::element_text(angle = 45, hjust = 1)) +
  ggplot2::guides(color = ggplot2::guide_colourbar(title = 'Expression'))
utils::write.csv(P$data, file.path(args$outdir, 'Fig5i_cardiac_marker_dotplot_ggplot_source_data.csv'), row.names = FALSE)
save_plot_pdf(P, file.path(args$outdir, 'Fig5i_cardiac_marker_dotplot_Seurat.pdf'), width = 12, height = 5.5)

df <- P$data
exp_mat <- df |>
  dplyr::select(-pct.exp, -avg.exp) |>
  tidyr::pivot_wider(names_from = id, values_from = avg.exp.scaled) |>
  as.data.frame()
row.names(exp_mat) <- exp_mat$features.plot
exp_mat <- as.matrix(exp_mat[, -1, drop = FALSE])

percent_mat <- df |>
  dplyr::select(-avg.exp, -avg.exp.scaled) |>
  tidyr::pivot_wider(names_from = id, values_from = pct.exp) |>
  as.data.frame()
row.names(percent_mat) <- percent_mat$features.plot
percent_mat <- as.matrix(percent_mat[, -1, drop = FALSE])

utils::write.csv(exp_mat, file.path(args$outdir, 'Fig5i_cardiac_marker_dot_heatmap_scaled_expression_matrix.csv'))
utils::write.csv(percent_mat, file.path(args$outdir, 'Fig5i_cardiac_marker_dot_heatmap_percent_matrix.csv'))

col_fun <- circlize::colorRamp2(c(-1, 0, 2), c('steelblue', 'white', 'darkred'))
column_ha <- ComplexHeatmap::HeatmapAnnotation(cluster = colnames(exp_mat), na_col = 'Gray')

layer_fun <- function(j, i, x, y, w, h, fill) {
  grid::grid.rect(x = x, y = y, width = w, height = h, gp = grid::gpar(col = NA, fill = NA))
  grid::grid.circle(
    x = x,
    y = y,
    r = ComplexHeatmap::pindex(percent_mat, i, j) / 100 * grid::unit(2, 'mm'),
    gp = grid::gpar(fill = col_fun(ComplexHeatmap::pindex(exp_mat, i, j)), col = NA)
  )
}

lgd_list <- list(
  ComplexHeatmap::Legend(
    labels = c(0, 0.25, 0.5, 0.75, 1),
    title = 'Fraction of cells',
    graphics = list(
      function(x, y, w, h) grid::grid.circle(x = x, y = y, r = 0 * grid::unit(2, 'mm'), gp = grid::gpar(fill = 'black')),
      function(x, y, w, h) grid::grid.circle(x = x, y = y, r = 0.25 * grid::unit(2, 'mm'), gp = grid::gpar(fill = 'black')),
      function(x, y, w, h) grid::grid.circle(x = x, y = y, r = 0.5 * grid::unit(2, 'mm'), gp = grid::gpar(fill = 'black')),
      function(x, y, w, h) grid::grid.circle(x = x, y = y, r = 0.75 * grid::unit(2, 'mm'), gp = grid::gpar(fill = 'black')),
      function(x, y, w, h) grid::grid.circle(x = x, y = y, r = 1 * grid::unit(2, 'mm'), gp = grid::gpar(fill = 'black'))
    )
  )
)

set.seed(230809)
hp <- ComplexHeatmap::Heatmap(
  exp_mat,
  name = 'Expression',
  heatmap_legend_param = list(title = 'Average expression'),
  col = col_fun,
  rect_gp = grid::gpar(type = 'none'),
  layer_fun = layer_fun,
  row_names_gp = grid::gpar(fontsize = 8),
  row_km = 9,
  border = 'black',
  top_annotation = column_ha,
  show_column_names = FALSE,
  show_parent_dend_line = FALSE
)

grDevices::pdf(file.path(args$outdir, 'Fig5i_cardiac_marker_dot_heatmap_ComplexHeatmap.pdf'), width = 11, height = 8, useDingbats = FALSE)
ComplexHeatmap::draw(hp, annotation_legend_list = lgd_list)
grDevices::dev.off()

saveRDS(combined6, file.path(args$outdir, 'Fig5_cardiac_lineages_annotated_seurat.rds'))
message('Done. Fig. 5f-i outputs saved to ', args$outdir)
