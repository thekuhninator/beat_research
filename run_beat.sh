#!/bin/bash

mkdir output
mkdir output/dataset1

# uncorrected
./beat.R -g example_data/dataset_1/dataset1_gene_counts.csv -a example_data/dataset_1/dataset1_metadata.csv -u -o  output/dataset1/dataset1_uncorrected -f diagnosis -n dataset1_uncorrected

# limma
./beat.R -g example_data/dataset_1_limma/dataset1_limma_gene_counts.csv -a example_data/dataset_1_limma/dataset1_limma_metadata.csv -o output/dataset1/dataset1_limma -n dataset1_limma -f diagnosis

# combat
./beat.R -g example_data/dataset_1_roman_combat/dataset1_roman_combat_gene_counts.csv -a example_data/dataset_1_roman_combat/dataset1_roman_combat_metadata.csv -n dataset1_roman_combat -o output/dataset1/dataset1_roman_combat -f diagnosis

# z-score
./beat.R -g example_data/dataset_1_z_score/dataset1_z_score_gene_counts.csv -a example_data/dataset_1_z_score/dataset1_z_score_metadata.csv -n dataset1_z_score -o output/dataset1/dataset1_z_score

# mean-centering
./beat.R -g example_data/dataset_1_mean_centered_roman/dataset_1_mean_centered_roman_gene_counts.csv -a example_data/dataset_1_mean_centered_roman/dataset_1_mean_centered_roman_metadata.csv -n dataset1_mean_centered_roman -o output/dataset1/dataset1_mean_centered_roman -f diagnosis

