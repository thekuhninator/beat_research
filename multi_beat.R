#!/usr/bin/env Rscript

#
# Install required libraries
#

# install neccessary packages
# install.packages("devtools", dependencies = TRUE)
# library(devtools)
# install_github('theislab/kBET')
# install.packages("bapred")
# install.packages("ggplot2")
# install.packages('dplyr')
# install.packages('gplots')
# install.packages('getopt')
# install.packages('sva')
# #install.packages('base64enc')
# install.packages('bapred')
# install.packages('readr')
# install.packages('matrixStats')
# # install kBEt
# install_github('theislab/kBET')
# # install biocmanager
# if (!requireNamespace("BiocManager", quietly = TRUE))
#   install.packages("BiocManager")
# # install m3c
# BiocManager::install("M3C")


#
# load in required libraries
#

#
# Install packages
#
# library(devtools)
# install.packages("png", repos = "http://cran.us.r-project.org")
# install.packages("raster", repos = "http://cran.us.r-project.org")
# install.packages('hash', repos = "http://cran.us.r-project.org")
# install.packages("grid", repos = "http://cran.us.r-project.org")
# install.packages("gridExtra", repos = "http://cran.us.r-project.org")
# install.packages("getopt", repos = "http://cran.us.r-project.org")
# install.packages("knitr", repos = "http://cran.us.r-project.org")
# install.packages("cowplot", repos = "http://cran.us.r-project.org")

# load in the neccessary libraries
library(ggplot2)
library(png)
library("raster")
library('hash')
library(gridExtra)
library(getopt)
library(knitr)


#
# necessary variables
#

# you can add in descirption for this
spec = matrix(c(
  'parent_dir', 'd', 1, "character",
  'output_dir','o', 1, "character",
  'output_name','n', 1, "character"
), byrow=TRUE, ncol=4)
opt = getopt(spec)


parent_dir <- opt$parent_dir
output_dir <- opt$output_dir
output_name <- opt$output_name

# create the direcotry if it does not already exist
dir.create(file.path(output_dir), showWarnings = FALSE)

# check if the parent_dir exists...
if (!dir.exists(parent_dir))
  stop("Cannot find parent directory. Are you sure the path is correct?")

#
# Get the kbet plot data
#
get_kbet_plot_data <- function(hvgs_hash, original_dataset_name, kbet_accept_hash)
{
  original_hvg_list <- hvgs_hash[[original_dataset_name]]
  #print('length of orignal hvgs list')
  num_original_hvgs <- length(original_hvg_list)
  #print(num_original_hvgs)
  
  # get the hvgs percent retained
  for (dataset_name in datasets)
  {
    #print('dataset name')
    #print(dataset_name)
    hvgs_retained <- Reduce(intersect, list(orig = original_hvg_list, other = hvgs_hash[[dataset_name]] ))
    #print('length of hvgs retained')
    #print(length(hvgs_retained))
    percent_retained <- length(hvgs_retained) / num_original_hvgs
    hvgs_retention[[dataset_name]] <- percent_retained
    #print('percent retained')
    #print(percent_retained)
  }
  
  hvgs_retention[["dataset1_combat"]]
  
  kbet_plot_data <- data.frame("kbet_acceptance" = numeric(), "dataset_name" = character(), "hvgs_retained" = numeric(), stringsAsFactors = FALSE)
  kbet_plot_data
  
  # get the kbet value for each dataset
  for (dataset_name in datasets)
  {
    kbet_val <- 1 - kbet_accept_hash[[dataset_name]]
    retained <- hvgs_retention[[dataset_name]]
    new_data <- list(kbet_val, dataset_name, retained)
    #print('new data')
    #print(new_data)
    #print(dataset_name)
    kbet_plot_data[nrow(kbet_plot_data) + 1,] = new_data
    #kbet_plot_data[nrow(kbet_plot_data) + 1, 1] = kbet_val
    #kbet_plot_data[nrow(kbet_plot_data) , 2] = dataset_name
    #kbet_plot_data[nrow(kbet_plot_data), 3] = retained
  }
  print(kbet_plot_data)
  return(kbet_plot_data)
}

kbet_hvg_scatterplot <- function(kbet_plot_data, output_dir, output_name)
{
  # create the output file
  file_name <- paste(output_name, '_kbet_hvg_scatterplot.png')
  output_file_path <- file.path(output_dir, file_name)
  png(output_file_path)
  
  print(kbet_plot_data)
  g <- ggplot(kbet_plot_data, aes(x=as.numeric(hvgs_retained), y = as.numeric(kbet_acceptance), color = dataset_name)) + 
    geom_point(size=4) + 
    labs(title = "KBET vs HVGs Plot", x = "Percent of HVGS Retained",
         y = "kBET Acceptance Rate", color="Dataset") + 
    scale_x_continuous(limits=c(0,1)) + 
    scale_y_continuous(limits=c(0,1))
  
  print(g)
  
  dev.off()
  return(output_file_path)
  
}

#
# Genereate grouped boxplot
#
grouped_boxplot <- function(boxplot_data, output_dir, output_name)
{
  # create the output file
  file_name <- paste(output_name, '_comparative_boxplot.png')
  output_file_path <- file.path(output_dir, file_name)
  png(output_file_path)
  
  # generate the plot
  g <- ggplot(comparative_boxplot_data, aes(x=dataset,
                                            y=mean,
                                            fill=as.factor(batch))) +
    geom_boxplot() +
    labs(title = "Comparative Grouped Boxplot", x ="Dataset", y = "Gene Mean Expression", fill = "Batch") +
    theme(axis.text.x = element_text(angle=30))
  
  
  # save the output file
  print(g)
  dev.off()
  
  return(output_file_path)
}

#
# Tile plots
#
tile_plots <- function(plots, dataset_names, plot_name, output_dir, output_name, combined_title)
{
  # create the output file
  file_name <- paste(output_name, plot_name, 'plot.png', sep = "_")
  output_file_path <- file.path(output_dir, file_name)
  png(output_file_path)
  plot_list <- list()
  for(name in dataset_names)
  {
    plot <- plots[[name]]
    plot <- plot + theme_bw()
    plot_list[[name]] <-plot
  }
  
  p_no_legend <- lapply(plot_list, function(x) x + theme(legend.position = "none"))
  legend <- cowplot::get_legend(plot_list[[1]] + theme(legend.position = "bottom"))
  title <- cowplot::ggdraw() + cowplot::draw_label(combined_title, fontface = "bold")
  p_grid <- cowplot::plot_grid(plotlist = p_no_legend, ncol = 2)
  g <- cowplot::plot_grid(title, p_grid, legend, ncol = 1, rel_heights = c(0.1, 1, 0.2))
  
  print(g)
  dev.off()
  return (output_file_path)
}

# load in necessary data from folders
beat_files <- list.files(path=parent_dir, recursive=TRUE, pattern="\\.beat$")

original_dataset_name <- NULL
hvgs_hash        <- hash()
hvgs_retention   <- hash()
kbet_accept_hash <- hash()
pca_plots        <- hash()
tsne_plots       <- hash()

retained_hvgs <- NULL
comparative_boxplot_data <- NULL
datasets <- c()

#print(beat_files)
#beat_file
for (beat_file in beat_files)
{
  load(file=file.path(parent_dir, beat_file))
  datasets <- c(datasets, dataset_name)
  
  if(original)
    original_dataset_name <- dataset_name
  
  # save the hvgs
  hvgs_hash[[dataset_name]] <- hvgs
  kbet_accept_hash[[dataset_name]] <- kbet_results['kBET.observed'][1,]
  boxplot_data['dataset'] <- dataset_name
  pca_plots[[dataset_name]] <- pca_plot
  tsne_plots[[dataset_name]] <- tsne_plot
  
  if(is.null(comparative_boxplot_data))
  {
    comparative_boxplot_data <- boxplot_data
  }
  else
  {
    comparative_boxplot_data <- rbind(comparative_boxplot_data, boxplot_data)
  }
  
}

generate_aggregate_report <- function(kbet_hvg_src, boxplot_src, pca_tile_src, tsne_tile_src, output_dir, output_name)
{
  html_string = 
    paste('
          <html>
          <head>
          <link rel="stylesheet" href="https://cdn.jsdelivr.net/gh/kognise/water.css@latest/dist/light.min.css">
          
          <style>body{ margin:0 100; background:whitesmoke; }</style>
          </head>
          <body>
          <h1>Batch Correction Aggregate Report for ', output_name, '</h1>
          
          
          <!-- *** Section 3 *** --->
          <h2>Principal Component Analysis</h2>
          <iframe width="1000" height="550" frameborder="0" seamless="seamless" scrolling="no" \
          src="', pca_tile_src, '" ></iframe> <p>Principal component analysis (PCA) is a technique used to emphasize variation and bring out strong patterns in a dataset. Its often used to make data easy to explore and visualize..</p>
          
          <h2>kBET - K-Nearest Neighbour Batch Effect Test Acceptance Rate vs HVGs Retention Rate</h2>
          <iframe width="1000" height="550" frameborder="0" seamless="seamless" scrolling="no" \
          src="', kbet_hvg_src, '"></iframe> <p>The K-Nearest Neighbor Batch Effect Test provides a test for batch effects in high-dimensional single-cell RNA sequencing data. It evaluates the accordance of replicates based on Pearsons chi^2 test. First, the algorithm creates k-nearest neighbour matrix and choses 10% of the samples to check the batch label distribution in its neighbourhood. If the local batch label distribution is sufficiently similar to the global batch label distribution, the chi^2-test does not reject the null hypothesis (that is "all batches are well-mixed"). The neighbourhood size k is fixed for all tests. Next, the test returns a binary result for each of the tested samples. Finally, the result of kBET is the average test rejection rate. The lower the test result, the less bias is introduced by the batch effect. kBET is very sensitive to any kind of bias. If kBET returns an average rejection rate of 1 for your batch-corrected data, you may also consider to compute the average silhouette width and PCA-based batch-effect measures to explore the degree of the batch effect. Learn more about kBET and batch effect correction in our publication.
          . The highly variable genes (HVGs) metric serves as a metric for biological preservation. The percent retained is calculated as followed (number of HVGs in both the corrrected dataset and uncorrected dataset)/ (number of hvgs in uncorrected dataset) </p>
          
          <!-- *** Section 1 *** --->
          <h2>T-Stochastic Neighbor Embedding (T-SNE) Tiled Plots</h2>
          <iframe width="1000" height="550" frameborder="0" seamless="seamless" scrolling="no" \
          src="' , tsne_tile_src, '"></iframe>
          <p>T-SNE is a t-distributed stochastic neighbor estimated.</p>
          
          <!-- *** Section 2 *** --->
          <h2>Comparative BoxPlot</h2>
          <iframe width="1000" height="550" frameborder="0" seamless="seamless" scrolling="no" \
          src="',  boxplot_src, '"></iframe>
          <p>The comparative boxplot shows the difference in the distributions between the different batches. If any of the boxes differ significanlty then there is a difference.</p>
          
          
          
          </body>
          </html>', sep="")
  
  file_name <- paste(output_name, '_batch_correction_report.html', sep="")
  output_file_path <- file.path(output_dir, file_name)
  html_report <- file(output_file_path)
  writeLines(c(html_string), html_report)
  close(html_report)
}

toBase64 <- function(image_file) {
  uri=image_uri(image_file)
}

# create the output
print('made it this far')
kbet_plot_data <- get_kbet_plot_data(hvgs_hash, original_dataset_name, kbet_accept_hash)  
print('made it this far2')
print(kbet_plot_data)
kbet_hvg_path  <- kbet_hvg_scatterplot(kbet_plot_data, output_dir, output_name)
print('made it this far2.5')
boxplot_path   <- grouped_boxplot(comparative_boxplot_data, output_dir, output_name)
print('made it this far3')
pca_tile_path  <- tile_plots(pca_plots, datasets, 'pca', output_dir, output_name, 'PCA Combined Plots')
print('made it this far4')
tsne_tile_path <- tile_plots(tsne_plots,  datasets,'tsne',   output_dir, output_name, 'T-SNE Combined Plots')

kbet_hvg_base64  <- toBase64(kbet_hvg_path)
boxplot_base64   <- toBase64(boxplot_path)
pca_tile_base64  <- toBase64(pca_tile_path)
tsne_tile_base64 <- toBase64(tsne_tile_path)

generate_aggregate_report(kbet_hvg_base64, boxplot_base64, pca_tile_base64, tsne_tile_base64, output_dir, output_name)