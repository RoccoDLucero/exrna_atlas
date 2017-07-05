#get_cell_type_references_from_excerpt_runs.R

input_pth <- '../cell_type_of_origin/input/'
bspn_sra_run_annots <- 'Brainspan+SRA_merged_results/run_annotations_full--fixed.txt'

bspn_sra_meta <- read.csv2(file = paste(input_pth, bspn_sra_run_annots, sep = ''),
                          sep = '\t')

liver_samples <- bspn_sra_meta[grep('liver', bspn_sra_meta$tissue), ]

kidney_cell_samples <- bspn_sra_meta[grep('kidney cells', bspn_sra_meta$tissue), ]
