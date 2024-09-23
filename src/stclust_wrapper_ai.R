#!/usr/bin/env Rscript

# Load necessary libraries
library(optparse)
library(spatialGE)
library(ggrepel)
library(tools)  # for file_path_sans_ext
library(utils)  # for untar

# Define command line options
option_list <- list(
  make_option(c("-e", "--exprmats"), type = "character", default = NULL, 
              help = "Comma-separated list of expression matrix file paths", metavar = "character"),
  make_option(c("-m", "--metas"), type = "character", default = NULL, 
              help = "Comma-separated list of metadata file paths", metavar = "character"),
  make_option(c("-a", "--samples"), type = "character", default = NULL, 
              help = "Comma-separated list of sample names", metavar = "character"),

  make_option(c("-w", "--ws"), type="double", default=0.025, 
              help="Weight to be applied to spatial distances (0-1). Default is 0.025", metavar="double"),
  make_option(c("-d", "--dist_metric"), type="character", default="euclidean", 
              help="Distance metric to be used. Default is 'euclidean'", metavar="character"),
  make_option(c("-l", "--linkage"), type="character", default="ward.D2", 
              help="Linkage method applied to hierarchical clustering. Default is 'ward.D2'", metavar="character"),
  make_option(c("-k", "--ks"), type="character", default="dtc", 
              help="Range of k values to assess. Default is 'dtc'", metavar="character"),
  make_option(c("-t", "--topgenes"), type="integer", default=2000, 
              help="Number of genes with highest spot-to-spot expression variation. Default is 2000", metavar="integer"),
  make_option(c("-s", "--deepSplit"), type="logical", default=FALSE, 
              help="Logical or integer (1-4) to control cluster resolution. Default is FALSE", metavar="logical"),
  make_option(c("-p", "--plot"), type="logical", default=FALSE, 
              help="Option to plot intermediate results. Default is FALSE", metavar="logical"),
  make_option(c("-o", "--output"), type="character", default="output.RData", 
              help="Path to save the output STlist object. Default is 'output.RData'", metavar="character"),
  make_option(c("-f", "--keep_fovs"), type = "character", default = NULL, 
              help = "Comma-separated list of FOVs to keep. If not specified, no filtering is done.", metavar = "character"),
  
  # Additional filter_data parameters
  make_option(c("--spot_minreads"), type = "integer", default = 20, 
              help = "Minimum number of reads per spot. Default is 20", metavar = "integer"),
  make_option(c("--rm_genes_expr"), type = "character", default = "^NegPrb", 
              help = "Regular expression to filter out genes. Default is '^NegPrb'", metavar = "character"),
  make_option(c("--rm_tissue"), type = "character", default = NULL, 
              help = "Comma-separated list of tissues to remove. Default is NULL (no filtering)", metavar = "character"),

  # Option for image archives
  make_option(c("-i", "--images"), type="character", default=NULL, 
              help="Path to a text file with tar archive file paths, one per row, for each sample.", metavar="character")
)

# Parse command line options
opt <- parse_args(OptionParser(option_list=option_list))

exprmats <- readLines(opt$exprmats)
metas <- readLines(opt$metas)
samples <- unlist(strsplit(opt$samples, ","))
image_archives <- readLines(opt$images)

# Load data
lung <- STlist(rnacounts = exprmats, spotcoords = metas, samples = samples)
lung <- filter_data(lung, spot_minreads=opt$spot_minreads, rm_genes_expr=opt$rm_genes_expr)
summ_df = summarize_STlist(lung)
lung <- transform_data(lung, method='sct')

# Filter FOVs if specified
if (!is.null(opt$keep_fovs)) {
  keep_fovs <- unlist(strsplit(opt$keep_fovs, ","))
  rm_fovs = summ_df$sample_name[!(summ_df$sample_name %in% keep_fovs)] 
  lung_subset <- filter_data(lung, rm_tissue=rm_fovs)
} else {
  lung_subset <- lung
}

# Extract images from tar archives and load them
for (tar_file in image_archives) {
  # Create a temporary directory for extraction
  temp_dir <- tempfile(pattern = "image_extract_")
  dir.create(temp_dir)
  
  # Extract the tar archive
  untar(tar_file, exdir = temp_dir)
  
  # Get all image file paths from the extracted directory
  image_files <- list.files(temp_dir, recursive = TRUE, full.names = TRUE)
  
  # Load images into the lung_subset
  lung_subset <- load_images(lung_subset, images=image_files)
  
  # Clean up temporary directory after use
  unlink(temp_dir, recursive = TRUE)
}

# Perform clustering using STclust
lung_subset <- STclust(x=lung_subset, ws=opt$ws, ks=opt$ks)
ti <- plot_image(lung_subset)
dom_p <- STplot(lung_subset, ks='dtc', ws=0.02, deepSplit=F, color_pal='discreterainbow')

# Optionally plot intermediate results
if (opt$plot) {
  cluster_p <- STplot(lung_subset, genes=c('KRT6A', 'HLA-B'), samples='Lung6_fov_4', color_pal='YlOrBr')
  print(cluster_p[[1]])
}

# Save the output STlist object
print(lung_subset)
save(lung_subset, file=opt$output)

# Print completion message
cat("Clustering completed and results saved to", opt$output, "\n")

traceback()

