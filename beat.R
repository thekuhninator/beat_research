#!/usr/bin/env Rscript

#
# Install the libraries
#

#install.packages("devtools", dependencies = TRUE, repos = "http://cran.us.r-project.org")
# library(devtools)
# install_github('theislab/kBET')
# install.packages("bapred", repos = "http://cran.us.r-project.org")
# install.packages("ggplot2", repos = "http://cran.us.r-project.org" )
# install.packages('dplyr', repos = "http://cran.us.r-project.org")
# install.packages('gplots', repos = "http://cran.us.r-project.org")
# install.packages('getopt', repos = "http://cran.us.r-project.org")
# install.packages('sva', repos = "http://cran.us.r-project.org")
# install.packages('bapred', repos = "http://cran.us.r-project.org")
# install.packages('readr', repos = "http://cran.us.r-project.org")
# install.packages('matrixStats', repos = "http://cran.us.r-project.org")
# if (!requireNamespace("BiocManager", quietly = TRUE))
#   install.packages("BiocManager")
# BiocManager::install("M3C")


#
# load in required libraries
#


library(devtools, quietly = TRUE, warn.conflicts = FALSE) 
library(kBET, quietly = TRUE, warn.conflicts = FALSE) 
library(knitr, quietly = TRUE, warn.conflicts = FALSE) 
#library(base64enc)
#library(bapred)
library(ggplot2, quietly = TRUE, warn.conflicts = FALSE) 
library(dplyr, quietly = TRUE, warn.conflicts = FALSE) 
library(gplots, quietly = TRUE, warn.conflicts = FALSE) 
library(getopt, quietly = TRUE, warn.conflicts = FALSE) 
library(M3C, quietly = TRUE, warn.conflicts = FALSE)
library(readr, quietly = TRUE, warn.conflicts = FALSE) 
library(matrixStats, quietly = TRUE, warn.conflicts = FALSE) 
#library(Rtsne)


#
# necessary variables
#

# you can add in descirption for this
spec = matrix(c(
  'input_counts', 'g', 1, "character",
  'input_annot', 'a', 1, "character",
  'output_dir','o', 1, "character",
  'dataset_name','n',1,"character",
  'uncorrected','u',2,"logical",
  'factor_of_interest', 'f', '2', 'character'
), byrow=TRUE, ncol=4)
opt = getopt(spec)

if ( is.null(opt$uncorrected    ) ) { opt$uncorrected    = FALSE     }
if ( is.null(opt$factor_of_interest    ) ) { opt$factor_of_interest    = NULL     }


input_counts <- opt$input_counts
input_metadata <- opt$input_annot
output_dir <- opt$output_dir
dataset_name <- opt$dataset_name
original <- opt$uncorrected
factor_of_interest <- opt$factor_of_interest

#input_counts = "C:/Users/Roman/Documents/Work/Depression_and_Immunology/England_Research/Data_No_Filter_Final/dataset_1/toups_unfiltered_gene_counts.csv"
#input_metadata = "C:/Users/Roman/Documents/Work/Depression_and_Immunology/England_Research/Data_No_Filter_Final/dataset_1/toups_unfiltered_metadata.csv"
#output_dir = "C:/Users/Roman/Documents/Work/Depression_and_Immunology/beast"
#output_name = "combined_dataset"

# create the direcotry if it does not already exist
dir.create(file.path(output_dir), showWarnings = FALSE)

# check if file paths exist
if (!file.exists(input_counts))
  stop("Cannot find gene counts file. Are you sure the path is correct?")
if (!file.exists(input_metadata))
  stop("Cannot find metadata file. Are you sure the path is correct?")

# Read in data file
gene_counts <- t(read.csv(input_counts, header = TRUE, row.names = 1, check.names=FALSE))
#Read in annotation file
annot <- read.csv(input_metadata, header = TRUE, row.names = 1, check.names = FALSE)

# check if factor of interest was inputted and exists in the metadataw
if(!is.null(factor_of_interest) && !factor_of_interest %in% names(annot))
  stop("Factor of interest not in the columns of the metadata. Are you sure it exists and is spelled correctly?")

# kBET function
# Input: Takes in gene conts data table, metadata data table, a foi, output_directory, and output name
# Ouptut: returns a dataframe containing data regarding the kBET test

kbet <- function(gene_counts, annot, output_dir, output_name) 
{
  batch = annot$batch
  # get the factor of interest (not currently used)
  #foi = annot[factor_of_interest]
  
  # additionall quantitative metrics
  #avedistVal = avedist(new_data, as.factor(batch) )#, as.factor(diagnosis))
  #pvcamVal = pvcam(new_data, batch, foi)
  #skewdivVal = skewdiv(new_data, as.factor(batch))
  #kldistVal = kldist(new_data, as.factor(batch))
  #sepscoreVal = sepscore(new_data, as.factor(batch) )
  # diffexprmVal = 
  # cobraVal = 
  
  data = gene_counts
  batch.estimate <- kBET(data, batch, plot = FALSE)
  
  # Capitalize Function
  capitalize <- function(x) {
    substr(x, 1, 1) <- toupper(substr(x, 1, 1))
    x
  }
  
  # Generate the Plot title
  # plotTitle = paste(capitalize(unlist(strsplit(name, "_"))), collapse= " ")
  plotTitle = paste(output_name, " kBET Plot", sep="")
  
  file_path <- paste( file.path(output_dir, output_name), "_kbet_plot.png", sep="")
  png(file=file_path)
  
  batch.estimate <- kBET(data, batch, plot=FALSE)
  
  plot.data <- data.frame(class=rep(c('observed', 'expected'),
                                    each=length(batch.estimate$stats$kBET.observed)), 
                          data =  c(batch.estimate$stats$kBET.observed,
                                    batch.estimate$stats$kBET.expected))
  
  print(ggplot(plot.data, aes(class, data)) + geom_boxplot() + 
          labs(x='Test', y='Rejection rate',title=paste(plotTitle,'kBET test results',sep=' ')) +
          theme_bw() +  
          scale_y_continuous(limits=c(0,1)))
  
  dev.off()
  
  results <- batch.estimate$summary
  
  # add extra quantitative measures to results
  #results['avedist'] = avedistVal
  #results['pvcam'] = pvcamVal
  #results['skewdiv'] = skewdivVal
  #results['kldist'] = kldistVal
  
  return( list("path" = file_path, "results" = results))
  
}

simpleCap <- function(x) {
  s <- strsplit(x, " ")[[1]]
  paste(toupper(substring(s, 1,1)), substring(s, 2),
        sep="", collapse=" ")
}

# Capitalize Function
capitalize <- function(x) {
  substr(x, 1, 1) <- toupper(substr(x, 1, 1))
  x
}

# author shakeel jessa
pca <- function(gene_counts, annot, output_dir, output_name, foi=NULL) {
  
  
  #if ( is.null(opt$PCx ) ) { opt$PCx = 1}
  #if ( is.null(opt$PCy ) ) { opt$PCy = 2}
  #xaxis <- opt$PCx
  #yaxis<- opt$PCy
  xaxis <- 1
  yaxis <- 2
  
  # Read a txt file and convert gene names (first column) into indeces
  #my_data <- read.csv(input_counts, header = TRUE, row.names = 1, check.names = FALSE)
  
  # Read a text file and create indeces
  #annot <- read.csv(input_annot, header = TRUE, row.names = 1, check.names = FALSE )
  
  my_data <- gene_counts
  annot['batch'] <- lapply(annot['batch'], as.character)
  
  gene_df <- t(my_data)   # Transpose data table
  gene_df <- data.frame(gene_df)  # Convert transposed table into a data frame
  
  # Remove all non-expressed and 0-variance genes from data frame
  express_df <- gene_df[,sapply(gene_df,function(x){!all(x==0)})]
  express_final_df <- express_df[vapply(express_df, function(x) length(unique(x)) > 1, logical(1L))]
  
  # Scale data frame
  scale_df <- scale(express_final_df)
  
  pca <- prcomp(scale_df) # Perform PCA analysis on scaled data and store
  
  # Plot PCx vs PCy and convert pca data into a data frame
  pca_df <- data.frame(pca$x)
  
  # Count number of columns in pca_df and store it
  pre_ncol <- ncol(pca_df)
  
  # Add annotations to pca dataframe
  pca_df <- data.frame(pca$x, annot)
  
  # Count number of columns in pca_df after adding annot and store
  post_ncol <- ncol(pca_df)
  
  # Generate range of values to iterate through for plot labeling
  label_ncol <- pre_ncol + 1
  
  # Generate Importance of Components and store it
  components <- summary(pca)
  
  # Pull values of the two highest Proportion of Variance's
  Prop_ofPCx <- components$importance[2, xaxis] * 100
  Prop_ofPCy <- components$importance[2, yaxis] * 100
  
  # Convert Variance values into strings
  Prop_ofPCx <- toString(Prop_ofPCx)
  Prop_ofPCy <- toString(Prop_ofPCy)
  
  # Create x and y axis labels with variance %
  x_label <- paste("PC",xaxis,",",' Var %: ', Prop_ofPCx)
  y_label <- paste("PC",yaxis,",", ' Var %: ', Prop_ofPCy)
  
  # Capture the rotation matrix in a data frame
  rotation_data <- data.frame(pca$rotation, variable=row.names(pca$rotation))
  
  PCxstring <- paste("PC",xaxis)
  PCystring <- paste("PC", yaxis)
  
  # Sort PCx/PCy columns from G->L and store in new data frames
  sorted_PCx <- rotation_data[order(rotation_data[,xaxis], decreasing=TRUE),, drop=FALSE]
  sorted_PCy <- rotation_data[order(rotation_data[,yaxis], decreasing=TRUE),, drop=FALSE]
  
  # Store gene names and values for top 10 genes
  PCx_names <- rownames(sorted_PCx[1:10,])
  PCx_values <- sorted_PCx[1:10,xaxis]
  PCx_10list <- paste(PCx_names, PCx_values, sep=': ')
  PCy_names <- rownames(sorted_PCy[1:10,])
  PCy_values <- sorted_PCy[1:10,yaxis]
  PCy_10list <- paste(PCy_names, PCy_values, sep=': ')
  
  
  sorted_PCxTen <- sorted_PCx[1:10,]
  sorted_PCyTen <- sorted_PCy[1:10,]
  pc_vector<- c(xaxis, yaxis)
  
  PCx_Ten <- sorted_PCxTen[pc_vector]
  PCy_Ten <- sorted_PCyTen[pc_vector]
  topTens <- rbind(PCx_Ten,PCy_Ten)
  newdata <- pca$x
  
  # output PCA
  file_path <- paste( file.path(output_dir, output_name), "_pca_plot.png", sep="")
  png(file=file_path)
  
  #name = paste(sapply(paste(unlist(strsplit(output_name, "_")), sep=" "), simpleCap), collapse=" ")
  #plotTitle = paste(capitalize(unlist(strsplit(name, "_"))), collapse= " ")
  plotTitle <- paste( output_name, " PCA", sep=" ")
  
  # Generate range of values to iterate through for plot labeling
  label_ncol <- grep("batch", colnames(pca_df))
  
  plot.new()
  
  if(is.null(foi))
  {
    g <- ggplot(pca_df, aes(x=PC1, y=PC2, color= pca_df[, label_ncol])) + 
      geom_point(size = 2) + labs(x =x_label, y = y_label, color=colnames(pca_df[label_ncol])) +
      ggtitle(plotTitle) + 
      theme(plot.title   = element_text(size=19, hjust = .5),  
            axis.title.y = element_text(size = 15), 
            axis.title.x = element_text(size = 15), 
            legend.title = element_text(size=14), 
            axis.text.x =  element_text(size=13),
            axis.text.y =  element_text(size=13))
  }
  else
  {
    post_ncol <- grep(foi, colnames(pca_df))
    g <- ggplot(pca_df, aes(x=pca_df[,xaxis], y=pca_df[,yaxis], color= pca_df[, label_ncol], shape = factor(pca_df[,post_ncol]))) + 
      geom_point(size = 2) + 
      labs(x =x_label, y = y_label, color=colnames(pca_df[label_ncol]), shape=colnames(pca_df[post_ncol])) +
      ggtitle(plotTitle) + 
      theme(plot.title   = element_text(size=19, hjust = .5),  
            axis.title.y = element_text(size = 15), 
            axis.title.x = element_text(size = 15), 
            legend.title = element_text(size=14), 
            axis.text.x =  element_text(size=13),
            axis.text.y =  element_text(size=13))
    
  }
  print(g);
  
  
  dev.off()
  
  return( list("path" = file_path, "plot" = g))
  
}




tsne_batch <- function(gene_counts, annot, ouput_dir, output_name)
{
  # output PCA
  file_path <- paste( file.path(output_dir, output_name), "_tsne_plot.png", sep="")
  png(file=file_path)
  
  #name = paste(sapply(paste(unlist(strsplit(output_name, " ")), sep=" "), simpleCap), collapse=" ")
  #plotTitle = paste(capitalize(unlist(strsplit(name, "_"))), collapse= " ")
  plotTitle <- paste( output_name, "T-SNE", sep=" ")
  invisible(capture.output( 
    g <- tsne(t(gene_counts), labels=as.factor(annot$batch), legendtitle ="Batch", dotsize = 2) + 
      ggtitle(plotTitle) +
      theme(plot.title = element_text(size=20, hjust = .5))
  ))
  #print(typeof(g))
  print(g)
  dev.off()
  
  return (list("path" = file_path, "plot" = g))
}

pca_m3c <- function(gene_counts, annot, output_dir, output_name, foi=NULL)
{
  
  file_name <- paste(file.path(output_dir, output_name), '.png', sep="")
  #png(file_name)
  
  # output PCA
  file_path <- paste( file.path(output_dir, output_name), "_pca_picture.png", sep="")
  png(file=file_path)
  
  name = paste(sapply(paste(unlist(strsplit(output_name, "[_]")), sep=" "), simpleCap), collapse=" ")
  plotTitle = paste(capitalize(unlist(strsplit(name, "_"))), collapse= " ")
  plotTitle <- paste( plotTitle, " PCA", sep="")
  
  g <- M3C::pca(gene_counts, labels=as.factor(annot$batch), legendtitle ="Batch", dotsize = 2) + 
    ggtitle(plotTitle) +
    theme(plot.title = element_text(size=20, hjust = .5)) +
    scale_size_continuous(range = c(1, 2))
  print(g)
  dev.off()
  
  return (list("path" = file_path, "plot" = g))
  
  
}

grouped_boxplot <- function(gene_counts, annot, output_dir, dataset_name) 
{
  # get list of samples for all batchs
  # for each batch variable
  # get list of samples for this batch
  # pull out their gene expression values
  # get the means and put it in a list
  # add it to a list
  # add batch x to names
  
  boxplot_data = NULL
  
  
  for (x in unique(annot$batch))
  {
    batchX = annot[annot$batch == x,]
    #print(batchX)
    
    rownames(batchX)
    
    gene_counts_subset <- t(gene_counts[rownames(batchX), ])
    
    #print(dim(gene_counts_subset))
    
    means = rowMeans(gene_counts_subset)
    
    mean_values <- as.vector(data.frame(means = rowMeans(gene_counts_subset))$means)
    
    batch_name = paste(c("batch_", x), collapse="", sep="")
    
    if(is.null(boxplot_data))
    {
      #print('first time')
      boxplot_data = data.frame("mean" = mean_values, "batch" = rep(x, each=length(mean_values)))
    }
    else
    {
      #print('not first time')
      new_data = data.frame("mean" = mean_values, "batch" = rep(x, each=length(mean_values)))
      #print(dim(new_data))
      
      boxplot_data = rbind(boxplot_data, new_data)
    }
    #print(dim(boxplot_data))
    
  }
  
  
  #print(boxplot)data
  file_name <- paste(dataset_name, '_comparative_boxplot.png', sep ="")
  output_file_name <- file.path(output_dir, file_name)
  png(output_file_name)
  
  
  g <- ggplot(boxplot_data, aes(x = factor(batch), y = mean, fill=factor(batch))) +
    geom_boxplot() + 
    labs(title = "Comparative Boxplot", x = "Batch", y = "Mean Gene Expression", fill = "Batch") + 
    theme(plot.title   = element_text(size=19, hjust = .5),  
          axis.title.y = element_text(size = 15), 
          axis.title.x = element_text(size = 15), 
          legend.title = element_text(size=14), 
          axis.text.x =  element_text(size=16),
          axis.text.y =  element_text(size=16))
  
  
  print(g)
  
  dev.off()
  
  return (list( "path" = output_file_name, "data" = boxplot_data))
}

getHvgs <- function(gene_counts, annot)
{
  hvgs <- NULL
  
  # let's iterate through all the batches
  for (batch_i in unique(annot$batch))
  {
    batch_counts <- t(gene_counts[rownames(annot[annot$batch == batch_i,]),]);
    batch_counts <- as.data.frame(cbind(batch_counts, 'var'=rowVars(batch_counts)))
    top_percent <- .10
    n_genes <- dim(batch_counts)[1] * top_percent
    top_varying <- rownames(head(batch_counts[order(batch_counts$var,decreasing=T),],n_genes))
    
    if(is.null(hvgs))
    {
      hvgs <- top_varying
    }
    else
    {
      hvgs <- Reduce(intersect, list(hvgs, top_varying))
    }

  }
  return (hvgs)
}

generate_report <- function(dataset_name, kbet_src, pca_src, tsne_src, boxplot_src)
{
  
  #dataset_name = "combined_dataset1"
  #pca_path = 'C:/Users/Roman/Documents/Work/Depression_and_Immunology/beast/dataset_1_stuff_picture.png'
  #pca_path = 'C:/Users/Roman/Documents/Work/Depression_and_Immunology/beast/dataset_1_stuff_picture.png'
  #pca_path = 'C:/Users/Roman/Documents/Work/Depression_and_Immunology/beast/dataset_1_stuff_picture.png'
  html_string = 
    paste('
          <html>
          <head>
              <link rel="stylesheet" href="https://cdn.jsdelivr.net/gh/kognise/water.css@latest/dist/light.min.css">
  
              <style>body{ margin:0 100; background:whitesmoke; }</style>
          </head>
          <body>
              <h1>Batch Correction Report for ', dataset_name, '</h1>
  
  
              <!-- *** Section 3 *** --->
              <h2>Principal Component Analysis</h2>
               <iframe width="1000" height="550" frameborder="0" seamless="seamless" scrolling="no" \
  src="', pca_src, '" ></iframe> <p>Principal Component Analysis (PCA) is a dimensionality reduction technique that emphasizes the variation in the data and allows us to see patterns in the data. The X axis represents the first principal component and its contributor rate. The Y axis represents the second component and its contributor rate. Points represent each sample. Sample colors and shapes are according to a group the sample belongs to. If the plot shows many samples of the same color (same batch) clustering together, this means there is a strong batch effect presense in the data. If the plot shows colors well mixed the batch effect is not severe in the data.

          </p>
              
               <h2>kBET - K-Nearest Neighbour Batch Effect test</h2>
               <iframe width="1000" height="550" frameborder="0" seamless="seamless" scrolling="no" \
  src="', kbet_src, '"></iframe> <p>The K-Nearest Neighbor Batch Effect Test (kBET) is a test metric used for assessing the severity of a batch effect in the data. The algorithm creates k-nearest neighbour matrix and choses 10% of the samples to check the batch label distribution in its neighbourhood. If the local batch label distribution is sufficiently similar to the global batch label distribution, the chi-squared test does not reject the null hypothesis (that is "all batches are well-mixed"). The neighbourhood size k is fixed for all tests. Next, the test returns a binary result for each of the tested samples. Finally, the result of kBET is the average test rejection rate. The lower the test result, the less bias is introduced by the batch effect. kBET is very sensitive to any kind of bias. By comparing the distribution of the expected and observed plots, one can see how severe the batch effect is in the data. The closer the two boxes are, the less severe, the batch effect. For more information about kbet, see their their paper or github.

</p>
      
              <!-- *** Section 1 *** --->
              <h2>T-Stochastic Neighbor Embedding (T-SNE)</h2>
               <iframe width="1000" height="550" frameborder="0" seamless="seamless" scrolling="no" \
  src="' , tsne_src, '"></iframe>
              <p>T-distributed Stochastic Neighbor Embedding (t-sne) is a machine learning algorithm for visualization. It is also a dimensionality reduction technique like PCA and is also useful in determing the severity of the batch effect by examining how strongly the colors (batches) are clustering together.</p>
              
              <!-- *** Section 2 *** --->
              <h2>Comparative BoxPlot</h2>
               <iframe width="1000" height="550" frameborder="0" seamless="seamless" scrolling="no" \
  src="',  boxplot_src, '"></iframe>
              <p>The comparative boxplot is a useful way of visualizing how the batches vary in the distribution of each gene mean expression value. The mene gene expression value across all samples within a batch are used as data points in constructing the comparative boxplot. If the boxes appear to be similar in their distribution the batch effect is not as severe for the dataset.

.</p>
  
  
  
          </body>
      </html>', sep="")
  
  file_name <- paste(dataset_name, '_batch_correction_report.html', sep="")
  output_file_path <- file.path(output_dir, file_name)
  html_report <- file(output_file_path)
  writeLines(c(html_string), html_report)
  close(html_report)
}

generateLogFile <- function(dataset_name, original, kbet_results, tsne_plot, pca_plot, boxplot_results, hvgs)
{

  file_name <- paste(dataset_name, '_beat_log.beat',sep="")
  output_file_path <- file.path(output_dir, file_name)
  save(dataset_name, original, kbet_results, tsne_plot, boxplot_data, pca_plot, hvgs, file=output_file_path)
}

toBase64 <- function(image_file) {
  uri=image_uri(image_file)
}

# get all the results
print('Running kbet...')
kbet_results <- kbet(gene_counts, annot, output_dir, dataset_name)
print('Creating t-sne plot...')
tsne_results <- tsne_batch(gene_counts, annot, output_dir, dataset_name)
print('Creating pca plot...')
pca_results <- pca_m3c(t(gene_counts), annot, output_dir, dataset_name, factor_of_interest)
print('Creating boxplot...')
boxplot_results <- grouped_boxplot(gene_counts, annot, output_dir, dataset_name);
print('Getting HVGs...')
hvgs <- getHvgs(gene_counts, annot)

tsne_path <- tsne_results$path
tsne_plot <- tsne_results$plot
pca_path  <- pca_results$path
pca_plot  <- pca_results$plot
kbet_path <- kbet_results$path
kbet_data <- kbet_results$results
boxplot_path <- boxplot_results$path
boxplot_data <- boxplot_results$data

kbet_base64     <- toBase64(kbet_path)
pca_base64      <- toBase64(pca_path)
tsne_base64     <- toBase64(tsne_path)
boxplot_base64  <- toBase64(boxplot_path)

# get all the base64 data of images
#tsne_base64 <- toBase64(tsne_path)
#pca_base64  <- toBase64(pca_path)
#box_base64  <- toBase64(boxplot_path)

# write them to an html report
generate_report(dataset_name, kbet_base64, pca_base64, tsne_base64, boxplot_base64)
# write them to the log file
generateLogFile(dataset_name, original, kbet_data, tsne_plot, pca_plot, boxplot_results, hvgs)

print('The report and log file have succesfully been generated!')
