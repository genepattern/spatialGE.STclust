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
  make_option(c("-f", "--keep_fovs"), type="character", default=NULL,
              help = "Comma-separated list of FOVs to keep", metavar = "character"),
  make_option(c("-i", "--image_paths"), type="character", default=NULL,
              help = "Comma-separated list of image paths", metavar = "character"),

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
              help="Path to save the output STlist object. Default is 'output.RData'", metavar="character")
)

# Parse command line options
opt <- parse_args(OptionParser(option_list=option_list))

# Split comma-separated file paths
exprmats <- unlist(strsplit(opt$exprmats, ","))
metas <- unlist(strsplit(opt$metas, ","))
samples <- unlist(strsplit(opt$samples, ","))
keep_fovs <- unlist(strsplit(opt$keep_fovs, ","))
image_paths <- unlist(strsplit(opt$image_paths, ","))

print(exprmats)

# Load data
lung <- STlist(rnacounts = exprmats, spotcoords = metas, samples = samples)
lung <- filter_data(lung, spot_minreads=20, rm_genes_expr='^NegPrb')
summ_df = summarize_STlist(lung)
lung <- transform_data(lung, method='sct')

keep_fovs
rm_fovs = summ_df$sample_name[!(summ_df$sample_name %in% keep_fovs)] 
lung_subset <- filter_data(lung, rm_tissue=rm_fovs)

# Perform clustering using STclust
lung_subset <- STclust(x=lung_subset, ws=opt$ws, ks = opt$ks)
for(image_path in image_paths){
  lung_subset <- load_images(lung_subset, images=image_path)
}
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

# Rscript --vanilla /genepattern/stclust_wrapper_ai.R \
# --exprmats /data/Lung5_Rep1/Lung5_Rep1-Flat_files_and_images/Lung5_Rep1_exprMat_file.csv,/data/Lung6/Lung6-Flat_files_and_images/Lung6_exprMat_file.csv \
# --metas /data/Lung5_Rep1/Lung5_Rep1-Flat_files_and_images/Lung5_Rep1_metadata_file.csv,/data/Lung6/Lung6-Flat_files_and_images/Lung6_metadata_file.csv \
# --samples Lung5_Rep1,Lung6 \
# --keep_fovs Lung5_Rep1_fov_6,Lung6_fov_6 \
# --image_paths /data/Lung5_Rep1/Lung5_Rep1-Flat_files_and_images/CellComposite,/data/Lung6/Lung6-Flat_files_and_images/CellComposite