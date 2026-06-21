#!/usr/bin/env Rscript

# Fig. 4c-f: Blood-lineage analysis comparing ciToti10-EMs and E8.5 natural embryos.

source('R/utils_common.R')
source('R/celltype_annotations.R')
require_packages(c('Seurat', 'dplyr', 'ggplot2', 'patchwork', 'harmony', 'pheatmap'))

args <- parse_cli_args(list(
  input = 'data/NS_harmony7.rds',
  outdir = 'results/Fig4',
  blood_clusters = '10,11,14,23',
  dims = '1:25',
  resolution = '0.1',
  top_marker_n = '10',
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

message_section('Loading integrated E8.5/ciToti10 object and subsetting blood clusters')
NS_harmony7 <- load_seurat_rds(args$input, project = 'Fig4_blood_source')
blood_clusters <- arg_vec(args$blood_clusters)
if (!'seurat_clusters' %in% colnames(NS_harmony7[[]])) stop('Input object must contain seurat_clusters metadata.')
cells_use <- rownames(NS_harmony7[[]])[as.character(NS_harmony7$seurat_clusters) %in% blood_clusters]
if (length(cells_use) == 0) stop('No cells found for blood clusters: ', paste(blood_clusters, collapse = ', '))
Blood <- subset(NS_harmony7, cells = cells_use)
Idents(Blood) <- 'seurat_clusters'

message_section('Re-running Harmony/UMAP for blood-lineage subset')
Blood <- run_harmony_pipeline(
  Blood,
  group.by.vars = 'orig.ident',
  dims = dims,
  resolution = as.numeric(args$resolution),
  assay = args$assay
)
Blood$blood_subcluster <- as.character(Blood$seurat_clusters)
Blood <- assign_celltype_from_cluster_map(Blood, blood_cluster_map, cluster_col = 'seurat_clusters', celltype_col = 'celltype')
Idents(Blood) <- 'celltype'

message_section('Fig. 4c: blood-lineage UMAP')
save_dimplot(
  Blood,
  out_pdf = file.path(args$outdir, 'Fig4c_blood_lineage_UMAP.pdf'),
  colors = blood_colors,
  pt.size = 1,
  width = 5,
  height = 4.5
)
export_umap_source(Blood, file.path(args$outdir, 'Fig4c_blood_lineage_UMAP_source_data.csv'), celltype_col = 'celltype')

message_section('Fig. 4d: cell-fraction correlation')
save_fraction_correlation_heatmap(
  Blood,
  out_pdf = file.path(args$outdir, 'Fig4d_blood_fraction_correlation_heatmap.pdf'),
  out_matrix_csv = file.path(args$outdir, 'Fig4d_blood_fraction_correlation_matrix.csv'),
  out_fraction_csv = file.path(args$outdir, 'Fig4d_blood_celltype_fraction_matrix.csv'),
  sample_col = 'orig.ident',
  width = 4,
  height = 4
)

message_section('Fig. 4e: conserved-marker heatmap')
cell_types <- c('Erythroblasts', 'Vascular Endothelia', 'Macrophage', 'HEP', 'Epithelia')
markers_list <- list()
for (ct in cell_types) {
  message('Finding conserved markers for ', ct)
  markers_list[[ct]] <- Seurat::FindConservedMarkers(
    Blood,
    ident.1 = ct,
    grouping.var = 'orig.ident',
    verbose = FALSE,
    only.pos = TRUE,
    min.diff.pct = 0.25,
    min.pct = 0.25,
    logfc.threshold = 0.25
  )
  utils::write.csv(markers_list[[ct]], file.path(args$outdir, paste0('Fig4e_', gsub('[^A-Za-z0-9]+', '_', ct), '_conserved_markers.csv')))
}

top_genes <- unique(unlist(lapply(markers_list, get_top_conserved_genes, top_n = as.integer(args$top_marker_n))))
genes_present <- top_genes[top_genes %in% rownames(Blood)]
if (length(genes_present) == 0) stop('No conserved marker genes found for Fig. 4e heatmap.')
missing_genes <- setdiff(top_genes, genes_present)
if (length(missing_genes) > 0) warning('Missing marker genes from object: ', paste(missing_genes, collapse = ', '))

expr <- get_assay_data_compat(Blood, assay = args$assay, slot_or_layer = 'data')[genes_present, , drop = FALSE]
expr_scaled <- scale_rows(expr)
annotation_col <- data.frame(
  Sample = Blood$orig.ident,
  Cluster = as.character(Idents(Blood)),
  row.names = colnames(expr_scaled)
)
annotation_col$Cluster <- factor(annotation_col$Cluster, levels = cell_types)
sample_levels <- unique(as.character(Blood$orig.ident))
annotation_col$Sample <- factor(annotation_col$Sample, levels = sample_levels)
new_col_order <- order(annotation_col$Cluster, annotation_col$Sample)
expr_ordered <- expr_scaled[, new_col_order, drop = FALSE]
annotation_col_ordered <- annotation_col[new_col_order, , drop = FALSE]

utils::write.csv(expr_ordered, file.path(args$outdir, 'Fig4e_conserved_marker_scaled_expression_matrix.csv'))
utils::write.csv(annotation_col_ordered, file.path(args$outdir, 'Fig4e_conserved_marker_column_annotation.csv'))

sample_palette <- c('#40B5AD', '#DC143C', '#999999', '#666666')
sample_colors <- setNames(sample_palette[seq_along(sample_levels)], sample_levels)
annotation_colors <- list(
  Cluster = blood_colors,
  Sample = sample_colors
)

grDevices::pdf(file.path(args$outdir, 'Fig4e_conserved_marker_heatmap.pdf'), width = 8, height = 7, useDingbats = FALSE)
pheatmap::pheatmap(
  expr_ordered,
  annotation_col = annotation_col_ordered,
  show_rownames = TRUE,
  show_colnames = FALSE,
  cluster_rows = TRUE,
  cluster_cols = FALSE,
  fontsize_row = 6,
  fontsize_col = 6,
  color = grDevices::colorRampPalette(c('navy', 'white', 'firebrick3'))(100),
  breaks = seq(-1, 1, length.out = 101),
  scale = 'none',
  annotation_colors = annotation_colors
)
grDevices::dev.off()

message_section('Fig. 4f: blood-lineage marker dot plot')
Blood <- add_sample_cluster(Blood, sample_col = 'orig.ident', output_col = 'sample_cluster')
Idents(Blood) <- 'sample_cluster'
save_dotplot(
  Blood,
  features = fig4f_blood_genes,
  out_pdf = file.path(args$outdir, 'Fig4f_blood_marker_dotplot.pdf'),
  out_csv = file.path(args$outdir, 'Fig4f_blood_marker_dotplot_source_data.csv'),
  colors = c('steelblue', 'white', 'darkred'),
  limits = c(-2.5, 2.5),
  width = 8,
  height = 5
)

saveRDS(Blood, file.path(args$outdir, 'Fig4_blood_lineages_annotated_seurat.rds'))
message('Done. Fig. 4c-f outputs saved to ', args$outdir)
