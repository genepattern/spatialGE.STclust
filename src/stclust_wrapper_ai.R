#!/usr/bin/env Rscript

# Load necessary libraries
library(optparse)
library(spatialGE)
library(ggrepel)
library(tools)  # for file_path_sans_ext
library(utils)  # for untar
library(zip)    # for unzipping zip files
library(magrittr)
library(jsonlite)

# Define command line options
option_list <- list(
  make_option(c("-e", "--stlist"), type = "character", default = NULL, 
              help = "STList object in .RDS file format to cluster on.", metavar = "character"),

  make_option(c("-w", "--ws"), type="double", default=0.025, 
              help="Weight to be applied to spatial distances (0-1). Default is 0.025", metavar="double"),
  make_option(c("-d", "--dist_metric"), type="character", default="euclidean", 
              help="Distance metric to be used. Default is 'euclidean'", metavar="character"),
  make_option(c("-l", "--linkage"), type="character", default="ward.D2", 
              help="Linkage method applied to hierarchical clustering. Default is 'ward.D2'", metavar="character"),
  make_option(c("-k", "--ks"), type="character", default="dtc", 
              help="Use Dynamic Tree Cut?", metavar="character"),
  make_option(c("-m", "--mink"), type="integer", default="2", 
                help="Minimum K", metavar="integer"),
  make_option(c("-n", "--maxk"), type="integer", default="5", 
                help="Maximum K", metavar="integer"),
  make_option(c("-t", "--topgenes"), type="integer", default=2000, 
              help="Number of genes with highest spot-to-spot expression variation. Default is 2000", metavar="integer"),
  make_option(c("-s", "--deepSplit"), type="logical", default=FALSE, 
              help="Logical or integer (1-4) to control cluster resolution. Default is FALSE", metavar="logical"),
  make_option(c("-o", "--output"), type="character", default="spatialGE_STclustered.rds", 
              help="Path to save the output STlist object. Default is 'spatialGE_STclustered.rds'", metavar="character"),
  make_option(c("-f", "--keep_fovs"), type = "character", default = NULL, 
              help = "Comma-separated list of FOVs to keep. If not specified, no filtering is done.", metavar = "character"),
  
  # Additional filter_data parameters
  make_option(c("--spot_minreads"), type = "integer", default = 20, 
              help = "Minimum number of reads per spot. Default is 20", metavar = "integer"),
  make_option(c("--rm_genes_expr"), type = "character", default = "^NegPrb", 
              help = "Regular expression to filter out genes. Default is '^NegPrb'", metavar = "character"),
  make_option(c("--rm_tissue"), type = "character", default = NULL, 
              help = "Comma-separated list of tissues to remove. Default is NULL (no filtering)", metavar = "character"),

  # Option for image zip file
  make_option(c("-i", "--images"), type="character", default=NULL, 
              help="Path to a zip file containing the 'CellComposite' folder with images.", metavar="character"),
              
  # job number for plot URLs
  make_option(c("-j", "--job"), type="character", default=NULL, 
                help="Job number for creating URLs to resulting plots", metavar="character")
)

# Parse command line options
opt <- parse_args(OptionParser(option_list=option_list))

# exprmats <- readLines(opt$exprmats)
# metas <- readLines(opt$metas)
# samples <- unlist(strsplit(opt$samples, ","))
image_zip <- opt$images

# Load data OLD
# lung <- STlist(rnacounts = exprmats, spotcoords = metas, samples = samples)
# lung <- filter_data(lung, spot_minreads=opt$spot_minreads, rm_genes_expr=opt$rm_genes_expr)
# summ_df = summarize_STlist(lung)
# lung <- transform_data(lung, method='sct')

## Read RDS from input
lung <- readRDS(opt$stlist)

# Filter FOVs if specified this should be done in preprocessing
# if (!is.null(opt$keep_fovs)) {
#   keep_fovs <- unlist(strsplit(opt$keep_fovs, ","))
#   rm_fovs = summ_df$sample_name[!(summ_df$sample_name %in% keep_fovs)] 
#   lung_subset <- filter_data(lung, rm_tissue=rm_fovs)
# } else {
#   lung_subset <- lung
# }

# Extract images from the zip file and load them
if (!is.null(image_zip)) {
  temp_dir <- tempfile(pattern = "image_extract_")
  dir.create(temp_dir)
  
  unzip(image_zip, exdir = temp_dir)
  
  cell_composite_path <- file.path(temp_dir, "CellComposite")
  image_files <- list.files(cell_composite_path, recursive = TRUE, full.names = TRUE)
  
  for (image_file in image_files) {
    lung <- load_images(lung, images=image_file)
  }

  # Clean up temporary directory after use
  unlink(temp_dir, recursive = TRUE)
}

# Assemble number of clusters
ks <- 'dtc'
if (opt$ks == "FALSE") {
  ks <- opt$mink:opt$maxk
}

# Perform clustering using STclust
lung_subset <- STclust(x=lung, ws=opt$ws, ks=ks)
ti <- plot_image(lung_subset)
dom_p <- STplot(lung_subset, ks=ks, ws=0.02, deepSplit=F, color_pal='smoothrainbow', ptsize=2, txsize=14)

# Save the plots to disk
for (i in seq_along(ti)) {
  file_name <- paste0("image_", i, ".jpg")
  ggsave(file_name, plot=ti[i])
}

for (i in names(dom_p)) {
  file_name <- paste0(i, ".jpg")
  ggsave(file_name, plot=dom_p[[i]])
}


# Generate annotation_variables_clusters table for domains.json
cluster_values = tibble::tibble()
for(i in names(lung_subset@spatial_meta)) {
  for(cl in grep('spagcn_|stclust_', colnames(lung_subset@spatial_meta[[i]]), value=T)) {
    cluster_values = dplyr::bind_rows(cluster_values, tibble::tibble(cluster=as.character(unique(lung_subset@spatial_meta[[i]][[cl]]))) %>% tibble::add_column(annotation=as.character(cl)))
  }
}
cluster_values = dplyr::distinct(cluster_values)
# write.table(cluster_values, 'annotation_variables_clusters.csv', quote=F, row.names=F, col.names=F, sep=',')

# Generate annotation_clusters for domains.json

display_label <- function(value) {
  # Split the input string by underscores
  parts <- strsplit(value, "_")[[1]]
  
  if (parts[1] == "stclust") {
    # Handle stclust prefix
    label <- "STclust"
    spw <- "No spatial weight"
    domains <- NULL
    deepsplit <- NULL
    
    for (part in parts[-1]) {
      if (startsWith(part, "spw")) {
        if (part != "0") {
          spw <- sub("spw", "spatial weight:", part)
        }
      } else if (startsWith(part, "k")) {
        domains <- sprintf("Domains (k): %02d", as.numeric(sub("k", "", part)))
      } else if (startsWith(part, "dspl")) {
        deepsplit <- ifelse(grepl("True", part), "DeepSplit=True", "DeepSplit=False")
      }
    }
    
    # Combine components for the display value
    display <- label
    if (!is.null(domains)) display <- paste(display, domains, sep = "; ")
    if (!is.null(deepsplit)) display <- paste(display, deepsplit, "Automatic mode (DynamicTreeCut)", sep = "; ")
    if (!is.null(spw)) display <- paste(display, spw, sep = "; ")
    return(display)
    
  } else if (parts[1] == "spagcn") {
    # Handle spagcn prefix
    label <- "SpaGCN"
    domains <- NULL
    
    for (part in parts[-1]) {
      if (startsWith(part, "k")) {
        domains <- sprintf("Domains (k): %02d", as.numeric(sub("k", "", part)))
      }
    }
    
    # Combine components for the display value
    display <- label
    if (!is.null(domains)) display <- paste(display, domains, sep = "; ")
    return(display)
    
  } else if (parts[1] == "insitutype") {
    # Handle insitutype prefix
    return("InSituType")
  }
  
  # Default case (if input does not match any known prefix)
  return("Unknown format")
}

annotations <- function(t) {
  if (!"annotation" %in% colnames(t)) {
    stop("The tibble does not have an 'annotation' column.")
  }
  
  # Get unique values from the annotation column
  unique_values <- unique(t$annotation)
  
  # Generate display labels for each unique value
  result <- tibble::tibble(
    label = sapply(unique_values, display_label),
    value = unique_values
  )
  
  return(result)
}

annotation_values = annotations(cluster_values)

# Build cluster image URLs
cluster_images = c()
for (i in names(dom_p)) {
  file_path <- paste0("/gp/jobResults/", opt$job, "/", i, ".jpg")
  cluster_images = c(cluster_images, file_path)
}

# Write domains.json
domains_data <- list(annotation_variables = annotation_values, annotation_variables_clusters = cluster_values, image_paths = cluster_images)
write_json(domains_data, path = 'domains.json', pretty = TRUE)

# Save the output STlist object
print(lung_subset)
saveRDS(lung_subset, file=opt$output)

# Print completion message
cat("Clustering completed and results saved to", opt$output, "\n")

traceback()
