#!/bin/bash

# uncorrected
./beat.R -g example_data/dataset_1/dataset1_gene_counts.csv -a example_data/dataset_1/dataset1_metadata.csv -u -o  output/dataset1/dataset1_uncorrected -f diagnosis -n dataset1_uncorrected

# limma
./beat.R -g example_data/dataset_1_limma/dataset1_limma_gene_counts.csv -a example_data/dataset_1_limma/dataset1_limma_metadata.csv -o output/dataset1/dataset1_limma -n dataset1_limma -f diagnosis

# combat
./beat.R -g example_data/dataset_1_combat/dataset1_combat_gene_counts.csv -a example_data/dataset_1_combat/dataset1_combat_metadata.csv -n dataset1_combat -o output/dataset1/dataset1_combat -f diagnosis

# z-score
./beat.R -g example_data/dataset_1_z_score/dataset1_z_score_gene_counts.csv -a example_data/dataset_1_z_score/dataset1_z_score_metadata.csv -n dataset1_z_score -o output/dataset1/dataset1_z_score

# mean-centering
./beat.R -g example_data/dataset_1_mean_centered/dataset_1_mean_centered_gene_counts.csv -a example_data/dataset_1_mean_centered/dataset_1_mean_centered_metadata.csv -n dataset1_mean_centered -o output/dataset1/dataset1_mean_centered -f diagnosis

