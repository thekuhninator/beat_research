# beat
BEAT - (Batch Effect Assessment Tool) is designed for researchers to be able to assess the severity of the batch effect present in their RNA-Seq data.

The batch effect assessment tool (BEAT) is a tool written to help researchers assess the severity of the batch effect in their dataset. The tool is written in R and run from the command line, taking in as input a gene counts matrix in the form of a csv file and metadata also in the form of a csv file and outputs a report with a PCA plot, tSNE plot, comparative boxplot, as well as results from kBET as well as a BEAT log file containing the information in the report.
After running several batch correction methods on a particular dataset and running BEAT on each dataset, the user can then input a folder of the outputted log files into multiBEAT to generate an aggregate report to compare the batch effect across the datasets using the metrics provided by BEAT. This can be used to help researchers choose which batch correction method works best for their particular dataset. 

<h1> Installation </h1>

To install the latest version 

`install.packages('beat.zip', repos = NULL, type = 'source')`

<h1> Usage of BEAT </h1>

The simplest way to run BEAT for a set of datasets that have been batch corrected is to place them in a directory with the following structure. The gene counts files must have gene_counts or genecounts in the name of the file or beat will not recognzie them. The metadata must have "metadata" "meta_data" or "annot" in the name of the file or else BEAT will not recognize it. This is done for simplicity and effecenciency. To see an example of how the gene counts csv file must be organized, see the example_data folder.

To run BEAT on every folder in a directory, one can construct the command by doing the following:

`for dataset_folder in ./*; do echo "beat -o <output_folder> -d $f" ;done`

Example of directory structure:
![Example of directory structure](https://github.com/thekuhninator/beat/blob/master/dir_structure.PNG?raw=true)

<h1> Dependencies and Docker </h1>

BEAT requires many dependencies to run. In addition to this many of the dependencies don't work with certain R versions. For this reason it is highly recommended that Docker is used to make the use of `BERP` stress free.

Docker documentation TBD.

<h1> Cite BEAT </h1>

If you find `BEAT` useful in your research please cite the related publication:

[Assessment of Batch Correction Methods for Whole Blood RNA-Seq Data](http://google.com)
<h1> BEAT Workflow Pipeline/Use Diagram </h1>

![BEAT Pipeline](https://github.com/thekuhninator/beat/blob/master/beat_pipeline.png?raw=true)
