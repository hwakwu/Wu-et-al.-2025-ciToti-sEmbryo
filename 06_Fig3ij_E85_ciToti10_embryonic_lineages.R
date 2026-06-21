#!/usr/bin/env Rscript

# Fig. 3i,j: Integration of ciToti10-EMs and E8.5 natural embryos, followed by embryonic-lineage subclustering.

source('R/utils_common.R')
source('R/celltype_annotations.R')
require_packages(c('Seurat', 'dplyr', 'ggplot2', 'patchwork', 'harmony'))

args <- parse_cli_args(list(
  E85 = 'data/E8.5_in_utero.rds',
  ciToti10 = 'data/ciToti10.rds',
  outdir = 'results/Fig3ij',
  global_dims = '1:20',
  subset_dims = '1:15',
  global_resolution = '0.5',
  subset_resolution = '0.05',
  embryonic_clusters = '2,6,18,24',
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

global_dims <- parse_dims(args$global_dims)
subset_dims <- parse_dims(args$subset_dims)

message_section('Loading E8.5 and ciToti10 objects')
E85 <- load_seurat_rds(args$E85, project = 'E8.5') |> set_orig_ident('E8.5')
ciToti10 <- load_seurat_rds(args$ciToti10, project = 'ciToti10') |> set_orig_ident('ciToti10')

Combined3 <- merge(E85, y = list(ciToti10), add.cell.ids = c('E8.5', 'ciToti10'), project = 'Fig3_E85_ciToti10')

message_section('Global Harmony integration')
Combined3 <- run_harmony_pipeline(
  Combined3,
  group.by.vars = 'orig.ident',
  dims = global_dims,
  resolution = as.numeric(args$global_resolution),
  assay = args$assay
)
Combined3$global_harmony_cluster <- as.character(Combined3$seurat_clusters)

message_section('Subsetting embryonic clusters and reclustering')
embryonic_clusters <- arg_vec(args$embryonic_clusters)
cells_use <- rownames(Combined3[[]])[Combined3$global_harmony_cluster %in% embryonic_clusters]
if (length(cells_use) == 0) stop('No cells matched --embryonic_clusters: ', paste(embryonic_clusters, collapse = ', '))
Em_subset <- subset(Combined3, cells = cells_use)

# Rerun Harmony on the subset. The original exploratory script used integrated.cca at this point;
# this release uses Harmony consistently because the combined object above is Harmony-integrated.
Em_subset <- run_harmony_pipeline(
  Em_subset,
  group.by.vars = 'orig.ident',
  dims = subset_dims,
  resolution = as.numeric(args$subset_resolution),
  assay = args$assay
)
Em_subset$subcluster <- as.character(Em_subset$seurat_clusters)
Em_subset <- assign_celltype_from_cluster_map(Em_subset, fig3_embryonic_subset_map, cluster_col = 'seurat_clusters', celltype_col = 'celltype')
Idents(Em_subset) <- 'celltype'

message_section('Saving Fig. 3i UMAP')
save_dimplot(
  Em_subset,
  out_pdf = file.path(args$outdir, 'Fig3i_embryonic_lineages_UMAP_split_by_sample.pdf'),
  split.by = 'orig.ident',
  colors = fig3_embryonic_colors,
  pt.size = 0.4,
  width = 8,
  height = 4
)
export_umap_source(Em_subset, file.path(args$outdir, 'Fig3i_embryonic_lineages_UMAP_source_data.csv'), celltype_col = 'celltype')

message_section('Saving Fig. 3j marker dot plot')
Em_subset <- add_sample_cluster(Em_subset, sample_col = 'orig.ident', output_col = 'sample_cluster')
cluster_levels <- c('Somites', 'MidHindGut', 'Neuroectoderm', 'Notochord')
levels_order <- unlist(lapply(cluster_levels, function(cl) grep(paste0('^', cl), unique(Em_subset$sample_cluster), value = TRUE)))
Em_subset$sample_cluster <- factor(Em_subset$sample_cluster, levels = levels_order)
Idents(Em_subset) <- 'sample_cluster'

save_dotplot(
  Em_subset,
  features = fig3j_marker_genes,
  out_pdf = file.path(args$outdir, 'Fig3j_embryonic_lineage_marker_dotplot.pdf'),
  out_csv = file.path(args$outdir, 'Fig3j_embryonic_lineage_marker_dotplot_source_data.csv'),
  colors = c('#006695', 'white', '#950000'),
  width = 8.5,
  height = 5
)

saveRDS(Combined3, file.path(args$outdir, 'Fig3ij_global_integrated_seurat.rds'))
saveRDS(Em_subset, file.path(args$outdir, 'Fig3ij_embryonic_subset_annotated_seurat.rds'))
message('Done. Fig. 3i,j outputs saved to ', args$outdir)
