#!/usr/bin/env Rscript

# Fig. 1n: Four-dimensional lineage score analysis.
# This script deliberately keeps the core analysis code in-line rather than hiding it in a helper:
#   1. Define the four gene sets used for scoring: Pluripotency/EPI, PrE, TE, and totipotent/2C-like.
#   2. Add Seurat module scores with AddModuleScore(features = list(pluri, pre, te, twoC)).
#   3. For each sample, extract Pluri1, PrE2, TE3, TwoC4 scores.
#   4. Convert the four scores to probabilities by exp(score * softmax_scale) / rowSums(...).
#   5. Project the four probabilities into a regular tetrahedron.
#   6. Export source data and interactive Plotly HTML for each sample.

source('R/utils_common.R')
require_packages(c('Seurat', 'dplyr', 'plotly', 'htmlwidgets'))

args <- parse_cli_args(list(
  input = 'data/ciToti0_4_E35_45.rds',
  samples = 'ciToti0,ciToti1,ciToti2,ciToti3,ciToti4,E4.5',
  outdir = 'results/Fig1n',
  assay = 'RNA',
  softmax_scale = '2',
  add_module_scores = 'TRUE'
))
ensure_dir(args$outdir)

message_section('Loading object for 4D lineage scoring')
obj <- load_seurat_rds(args$input, project = 'ciToti0_4_E35_45')
DefaultAssay(obj) <- args$assay
obj <- join_layers_if_needed(obj, assay = args$assay)

# -------------------------------------------------------------------------
# 1. Gene sets used in the original four-dimensional scoring analysis
# -------------------------------------------------------------------------
toti_genes_pure <- c(
  'Zfp809', 'Zfp560', 'Zfp352', 'Tdpoz1', 'Tdpoz3', 'Tdpoz4', 'Gm2016',
  'Gm13078', 'Gm8300', 'Obox3', 'Obox1', 'Usp17la', 'Usp17ld', 'Usp17lc',
  'Usp17le', 'Zscan4c', 'Zscan4f', 'Zscan4d', 'Tcstv1', 'Tcstv3'
)

epi_massive <- c(
  'Pou5f1', 'Nanog', 'Sox2', 'Klf2', 'Esrrb', 'Zfp42', 'Fgf4', 'Tdgf1',
  'Otx2', 'Utf1', 'Dppa4', 'Dppa5a', 'Lin28a', 'Tet1', 'Nodal', 'Zscan10',
  'Etv5', 'Mphosph8', 'Gdf3', 'Ifitm3', 'Trh', 'Lck', 'Eif4ebp1'
)

te_massive <- c(
  'Tfap2c', 'Krt8', 'Krt18', 'Krt19', 'Gpx3', 'S100a11', 'Cldn6', 'Peg10',
  'Phlda2', 'Cdkn1c', 'H19', 'Igf2', 'Msx2', 'Ndrg1', 'Cldn7', 'Lgals1'
)

pre_massive <- c(
  'Gata6', 'Sox17', 'Pdgfra', 'Col4a1', 'Col4a2', 'Lama1', 'Cubn', 'Sox7',
  'Serpinh1', 'Hnf1b', 'Ihh', 'Apoc2', 'Ttr', 'Loxl2', 'Serpine2', 'Dpp4', 'Rbp4'
)

# The order below is intentional. Seurat appends the list index to each score name,
# producing Pluri1, PrE2, TE3, and TwoC4, matching the original analysis code.
pluri <- epi_massive
pre <- pre_massive
te <- te_massive
twoC <- toti_genes_pure

# Export the gene sets for source-data/reproducibility.
gene_set_df <- data.frame(
  gene = c(pluri, pre, te, twoC),
  score = c(rep('Pluri/EPI', length(pluri)), rep('PrE', length(pre)), rep('TE', length(te)), rep('TwoC/Totipotent-like', length(twoC)))
)
utils::write.csv(gene_set_df, file.path(args$outdir, 'Fig1n_4D_score_gene_sets.csv'), row.names = FALSE)

# -------------------------------------------------------------------------
# 2. Add module scores
# -------------------------------------------------------------------------
add_scores <- toupper(as.character(args$add_module_scores)) %in% c('TRUE', 'T', '1', 'YES')
if (add_scores) {
  message_section('Adding module scores with Seurat::AddModuleScore')
  obj <- Seurat::AddModuleScore(
    object = obj,
    features = list(pluri, pre, te, twoC),
    name = c('Pluri', 'PrE', 'TE', 'TwoC')
  )
}

score_cols <- c('Pluri1', 'PrE2', 'TE3', 'TwoC4')
missing_score_cols <- setdiff(score_cols, colnames(obj[[]]))
if (length(missing_score_cols) > 0) {
  stop('Missing expected module score columns: ', paste(missing_score_cols, collapse = ', '),
       '. Run with --add_module_scores TRUE or provide an object already containing these columns.')
}

saveRDS(obj, file.path(args$outdir, 'Fig1n_object_with_module_scores.rds'))

# -------------------------------------------------------------------------
# 3. Define softmax transformation and regular tetrahedron projection
# -------------------------------------------------------------------------
softmax_scores <- function(score_matrix, scale = 2) {
  M <- exp(as.matrix(score_matrix) * scale)
  M <- sweep(M, 1, rowSums(M), '/')
  colnames(M) <- c('Pluri', 'PrE', 'TE', 'TwoC')
  M
}

# Orthogonal regular tetrahedron vertices.
# Rows correspond to Pluri/EPI, PrE, TE, and TwoC/totipotent-like vertices.
tetrahedron_vertices <- matrix(
  c(
    sqrt(8 / 9), 0, -1 / 3,
    -sqrt(2 / 9), sqrt(2 / 3), -1 / 3,
    -sqrt(2 / 9), -sqrt(2 / 3), -1 / 3,
    0, 0, 1
  ),
  ncol = 3,
  byrow = TRUE
)
colnames(tetrahedron_vertices) <- c('x', 'y', 'z')
rownames(tetrahedron_vertices) <- c('Pluri', 'PrE', 'TE', 'TwoC')

project_sample_to_tetrahedron <- function(global_obj,
                                          target_sample,
                                          sample_col = 'orig.ident',
                                          score_cols = c('Pluri1', 'PrE2', 'TE3', 'TwoC4'),
                                          softmax_scale = 2) {
  cells <- rownames(global_obj[[]])[as.character(global_obj[[sample_col]][, 1]) == target_sample]
  if (length(cells) == 0) stop('No cells found for sample: ', target_sample)
  obj_sub <- subset(global_obj, cells = cells)

  raw_scores <- Seurat::FetchData(obj_sub, vars = score_cols)
  colnames(raw_scores) <- c('raw_Pluri1', 'raw_PrE2', 'raw_TE3', 'raw_TwoC4')

  probabilities <- softmax_scores(raw_scores, scale = softmax_scale)
  xyz <- probabilities %*% tetrahedron_vertices[c('Pluri', 'PrE', 'TE', 'TwoC'), ]
  colnames(xyz) <- c('x', 'y', 'z')

  df <- cbind(as.data.frame(xyz), as.data.frame(probabilities), as.data.frame(raw_scores))
  df$Cell_Barcode <- rownames(probabilities)
  df$sample <- target_sample
  df <- df[, c('sample', 'Cell_Barcode', 'x', 'y', 'z', 'Pluri', 'PrE', 'TE', 'TwoC',
               'raw_Pluri1', 'raw_PrE2', 'raw_TE3', 'raw_TwoC4')]
  df <- df[is.finite(df$x) & is.finite(df$y) & is.finite(df$z), , drop = FALSE]
  df
}

plot_tetrahedron <- function(df,
                             target_sample,
                             color_by = 'TwoC',
                             point_size = 10,
                             line_width = 8) {
  V <- tetrahedron_vertices
  verts <- data.frame(
    x = V[, 1], y = V[, 2], z = V[, 3],
    label = c('<b>Pluripotency</b>', '<b>PrE</b>', '<b>TE</b>', '<b>Totipotent-like</b>')
  )
  edge_pairs <- list(c(1, 2), c(1, 3), c(1, 4), c(2, 3), c(2, 4), c(3, 4))
  edge_df <- do.call(rbind, lapply(edge_pairs, function(e) {
    rbind(
      data.frame(x = verts$x[e], y = verts$y[e], z = verts$z[e]),
      data.frame(x = NA, y = NA, z = NA)
    )
  }))

  df$hover <- sprintf(
    'Pluri: %.1f%%<br>PrE: %.1f%%<br>TE: %.1f%%<br>2C/Toti: %.1f%%',
    df$Pluri * 100, df$PrE * 100, df$TE * 100, df$TwoC * 100
  )

  plotly::plot_ly() |>
    plotly::add_trace(
      type = 'mesh3d', x = verts$x, y = verts$y, z = verts$z,
      i = c(0, 0, 0, 1), j = c(1, 2, 3, 2), k = c(2, 3, 1, 3),
      opacity = 0, color = I('white'), hoverinfo = 'none', showscale = FALSE
    ) |>
    plotly::add_trace(
      type = 'scatter3d', mode = 'lines', data = edge_df,
      x = ~x, y = ~y, z = ~z,
      line = list(color = 'black', width = line_width),
      hoverinfo = 'none', showlegend = FALSE
    ) |>
    plotly::add_markers(
      data = df, x = ~x, y = ~y, z = ~z,
      text = ~hover, hoverinfo = 'text',
      marker = list(
        size = point_size,
        opacity = 0.8,
        color = df[[color_by]],
        colorscale = 'Viridis',
        cmin = 0,
        cmax = 1,
        showscale = TRUE,
        colorbar = list(title = 'Totipotency score', len = 0.6)
      ),
      name = target_sample
    ) |>
    plotly::add_text(
      data = transform(verts, x = 1.2 * x, y = 1.2 * y, z = 1.15 * z),
      x = ~x, y = ~y, z = ~z, text = ~label,
      textfont = list(size = 16, color = 'black'),
      showlegend = FALSE
    ) |>
    plotly::layout(
      title = list(text = paste(target_sample, '4D lineage score'), y = 0.95),
      scene = list(
        aspectmode = 'data',
        camera = list(eye = list(x = 1.3, y = -1.3, z = 0.6)),
        xaxis = list(visible = FALSE),
        yaxis = list(visible = FALSE),
        zaxis = list(visible = FALSE)
      ),
      margin = list(t = 50)
    ) |>
    plotly::config(
      toImageButtonOptions = list(
        format = 'png',
        filename = paste0(target_sample, '_4D_lineage_score'),
        width = 1000,
        height = 1000,
        scale = 4
      )
    )
}

# -------------------------------------------------------------------------
# 4. Run per-sample scoring, projection, and plotting
# -------------------------------------------------------------------------
target_samples <- arg_vec(args$samples)
if (length(target_samples) == 0) stop('No target samples supplied. Use --samples sample1,sample2,...')

all_source <- list()
for (sample in target_samples) {
  message_section('Projecting ', sample, ' into 4D lineage tetrahedron')
  df <- project_sample_to_tetrahedron(
    global_obj = obj,
    target_sample = sample,
    sample_col = 'orig.ident',
    score_cols = score_cols,
    softmax_scale = as.numeric(args$softmax_scale)
  )
  utils::write.csv(df, file.path(args$outdir, paste0('Fig1n_', sample, '_4D_lineage_score_source_data.csv')), row.names = FALSE)
  p <- plot_tetrahedron(df, target_sample = sample, color_by = 'TwoC')
  htmlwidgets::saveWidget(p, file.path(args$outdir, paste0('Fig1n_', sample, '_4D_lineage_score.html')), selfcontained = TRUE)
  all_source[[sample]] <- df
}

all_source_df <- dplyr::bind_rows(all_source)
utils::write.csv(all_source_df, file.path(args$outdir, 'Fig1n_all_samples_4D_lineage_score_source_data.csv'), row.names = FALSE)
saveRDS(all_source, file.path(args$outdir, 'Fig1n_all_samples_4D_lineage_score_source_data.rds'))
message('Done. Fig. 1n outputs saved to ', args$outdir)
