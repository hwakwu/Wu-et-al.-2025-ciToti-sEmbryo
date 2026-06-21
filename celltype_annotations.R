# Manual cluster-to-celltype annotations and marker lists used in figure scripts

ciToti4_cluster_map <- c(
  '0' = 'Intermediate', '1' = 'Epiblast', '2' = 'PrE', '3' = 'Intermediate',
  '4' = 'TE', '5' = 'Intermediate'
)

ciToti7_cluster_map <- c(
  '0' = 'Epiblast', '1' = 'Epiblast', '2' = 'PS', '3' = 'Intermediate',
  '4' = 'PS', '5' = 'Em VE', '6' = 'Epiblast', '7' = 'Epiblast',
  '8' = 'Epiblast', '9' = 'Epiblast', '10' = 'Intermediate', '11' = 'PS',
  '12' = 'ExE ecto', '13' = 'Primitive Node Cells', '14' = 'Epiblast',
  '15' = 'ExE ecto', '16' = 'Primitive Node Cells', '17' = 'ExE ecto',
  '18' = 'PGCs', '19' = 'Parietal endo', '20' = 'ExE VE'
)

ciToti8_cluster_map <- c(
  '0' = 'Amnion', '1' = 'Neuromesodermal progenitors', '2' = 'Epiblast',
  '3' = 'Mesenchyme', '4' = 'Visceral endo', '5' = 'ExE endo',
  '6' = 'Neuromesodermal progenitors', '7' = 'ExE ecto',
  '8' = 'Haematoendothelial progenitors', '9' = 'PS', '10' = 'ExE endo',
  '11' = 'Intermediate', '12' = 'ExE ecto', '13' = 'Gut progenitors',
  '14' = 'Amnion', '15' = 'Parietal endo', '16' = 'Intermediate'
)

ciToti10_cluster_map <- c(
  '0' = 'Mixed', '1' = 'Amnion', '2' = 'Neuroectoderm', '3' = 'Allantois',
  '4' = 'PS', '5' = 'Epithelium', '6' = 'Somites', '7' = 'Intermediate',
  '8' = 'PGCs', '9' = 'ExE endo', '10' = 'Erythroblast',
  '11' = 'Vascular endothelium', '12' = 'ExE endo', '13' = 'Surface ecto',
  '14' = 'Primitive erythroid cells', '15' = 'Parietal endo',
  '16' = 'Cardiacmyocyte', '17' = 'Intermediate', '18' = 'Gut progenitors',
  '19' = 'Somites', '20' = 'Endothelium', '21' = 'ExE ecto',
  '22' = 'Cardiac meso', '23' = 'Haematoendothelial progenitors', '24' = 'Notochord'
)

fig2_integrated_cluster_map <- c(
  '0' = 'Mesenchyme', '1' = 'Epiblast', '2' = 'Epiblast', '3' = 'Visceral endo',
  '4' = 'Nascent meso', '5' = 'Intermedia', '6' = 'Neuroectodermal progenitors',
  '7' = 'ExE ecto', '8' = 'Parietal endo', '9' = 'ExE endo', '10' = 'Trophoblast',
  '11' = 'EPC', '12' = 'Amnion', '13' = 'Haematoendothelial progenitors',
  '14' = 'ExE endo', '15' = 'PS', '16' = 'Definitive endo', '17' = 'Definitive endo',
  '18' = 'ExE ecto', '19' = 'ExE endo', '20' = 'Intermedia', '21' = 'Intermedia'
)

fig3_embryonic_subset_map <- c(
  '1' = 'Somites', '0' = 'Neuroectoderm', '2' = 'MidHindGut', '3' = 'Notochord'
)

blood_cluster_map <- c(
  '0' = 'Erythroblasts', '1' = 'Vascular Endothelia', '2' = 'Erythroblasts',
  '3' = 'Macrophage', '4' = 'HEP', '5' = 'Epithelia'
)

cardiac_cluster_map <- c(
  '0' = 'Erythrocytes', '1' = 'V-CMs', '2' = 'Undefined', '3' = 'OFT',
  '4' = 'Artery-EC', '5' = 'Epicardium', '6' = 'Endocardium', '7' = 'A-CMs',
  '8' = 'HSC', '9' = 'CNCCs', '10' = 'Erythrocytes'
)

fig2_colors <- c(
  'Mesenchyme' = 'mediumvioletred', 'Neuroectodermal progenitors' = 'lightblue',
  'Visceral endo' = '#e6b0aa', 'Nascent meso' = 'lightslategray',
  'Intermedia' = '#1565c0', 'Intermediate' = '#1565c0', 'ExE ecto' = '#7b1fa2',
  'Parietal endo' = '#ec407a', 'ExE endo' = '#880e4f', 'Trophoblast' = 'darkorange',
  'EPC' = 'orchid', 'Amnion' = '#ca6f1e', 'Aminion' = '#ca6f1e',
  'Haematoendothelial progenitors' = 'forestgreen', 'PS' = 'cadetblue',
  'Definitive endo' = '#ec7063', 'Epiblast' = 'darkseagreen', 'VE' = '#e6b0aa'
)

fig3_atlas_colors <- c(
  'PrE' = 'mediumvioletred', 'Intermedia' = 'gray', 'Intermediate' = 'gray',
  'Epiblast' = 'steelblue', 'Polar TE' = '#993399', 'Mural TE' = '#996699',
  'TE' = '#996699', 'PS' = '#1565c0', 'Em VE' = 'orchid', 'ExE ecto' = '#ec407a',
  'ExE VE' = '#7b1fa2', 'PGCs' = 'gold', 'Parietal endo' = '#d1c4e9',
  'Amnion' = '#ca6f1e', 'Aminion' = '#ca6f1e',
  'Neuromesodermal progenitors' = 'forestgreen', 'Mesenchyme' = 'cadetblue',
  'Visceral endo' = '#880e4f', 'Mixed' = 'slategray',
  'Haematoendothelial progenitors' = 'darkorange', 'Gut progenitors' = 'darkgoldenrod',
  'Neuroectoderm' = 'cornflowerblue', 'Allantois' = 'navajowhite',
  'Epithelium' = '#9ccc65', 'Somites' = 'mediumaquamarine', 'ExE endo' = '#e6b0aa',
  'Erythroblast' = '#DC143C', 'Vascular endothelium' = 'orangered',
  'Surface ecto' = 'darkseagreen', 'Primitive erythroid cells' = '#ec7063',
  'Cardiacmyocyte' = 'royalblue', 'Endothelium' = '#808000',
  'Cardiac meso' = 'lightslategray', 'Notochord' = 'lightgreen'
)

fig3_embryonic_colors <- c(
  'MidHindGut' = 'forestgreen', 'Neuroectoderm' = '#7b1fa2',
  'Somites' = 'darkgoldenrod', 'Notochord' = '#1565c0'
)

blood_colors <- c(
  'Erythroblasts' = 'Salmon', 'Vascular Endothelia' = '#457cba',
  'Macrophage' = '#DC143C', 'HEP' = '#a2c3e9', 'Epithelia' = '#14a61a'
)

cardiac_colors <- c(
  'V-CMs' = 'forestgreen', 'OFT' = '#1565c0', 'Artery-EC' = '#FFC000',
  'Epicardium' = 'steelblue', 'Endocardium' = 'lightgreen', 'A-CMs' = 'purple',
  'CNCCs' = 'lightslategray', 'HSC' = '#DC143C', 'Erythrocytes' = '#FF4433',
  'Undefined' = 'grey80'
)

fig2k_average_genes <- c(
  'Taldo1', 'Mrpl24', 'Mrpl44', 'Ruvbl2', 'Mrto4', 'Unc93b1', 'Dohh',
  'Fcf1', 'Ube2l3', 'Mrpl42', 'Col1a1', 'Dnajb13', 'Fli1', 'Noto',
  'Cfap44', 'Hoxd9', 'Pkd1l1', 'Unc5c', 'Nkx6-1', 'Gdf6'
)

fig2k_epiblast_genes <- c('Tbx3', 'Fbxo15', 'Esrrb', 'Zfp42', 'Nanog', 'Fgf5', 'Eomes', 'Pou3f1', 'Sox2', 'T')

fig3j_marker_genes <- c(
  'Meox1', 'Tcf15', 'Snai1', 'Foxc2', 'Apoa1', 'Cdh6', 'Sox17', 'Gata4',
  'Sox2', 'Otx2', 'Adgrv1', 'Sox1', 'Pax3', 'T', 'Shh', 'Sox9', 'Foxa2', 'Noto', 'Chrd'
)

fig4f_blood_genes <- c(
  'Hba-x', 'Hbb-bh1', 'Klf1', 'Smim1', 'Zfpm1', 'Gata1', 'Ank1',
  'Hba-a1', 'Hbb-bs', 'Epb42', 'Runx1', 'Spi1', 'Laptm5', 'Tyrobp',
  'Fcer1g', 'Plcg2', 'Tspo', 'Cebpb', 'Cd44', 'Lyz2'
)

fig5i_cardiac_genes <- c(
  'Nkx2-5', 'Cryab', 'Tnni1', 'Tnni3', 'Myl3', 'Myl7', 'Tnnt2', 'Myh6',
  'Myh7', 'Cacna1c', 'Ryr2', 'Ttn', 'Myom1',
  'Foxa1', 'Foxa2', 'Cldn4', 'Cldn6', 'Cldn7', 'Epcam', 'Crb3', 'Spint1',
  'Wnt4', 'Wnt6', 'Tfap2a', 'Gata3', 'Tacstd2', 'Tac2', 'Kcnk1', 'Dlx5',
  'Dll4', 'Robo4', 'Erg', 'Mfng', 'Cldn5', 'Cdh5', 'Sox18', 'Esam',
  'Wt1', 'Tbx18', 'Tcf21', 'Aldh1a1', 'Upk1b', 'Upk3b', 'Prkcb', 'Ildr2',
  'Nfatc1', 'Nrg1', 'Tgfbr3', 'Pde1c', 'Trpm3', 'Ednrb', 'Cped1', 'Itga9',
  'Rhag', 'Hbb-bh1', 'Hba-a2', 'Hba-x', 'Hbb-y', 'Hba-a1', 'Gypa', 'Alas2',
  'Csf1r', 'Ptprc', 'Tyrobp', 'Spi1', 'Fcer1g', 'Rac2', 'Coro1a', 'Laptm5'
)
