# Spatial Transcriptomics Clustering Script
=============================================

## Description
This script performs clustering analysis on gene expression data using the spatialGE package. It performs the STClust function from SpatialGE. The outputting .RDS object is a STList object for use in downstream Spatial GE analysis / visualization. 

## Module Details
- Authors: Edwin Huang
- Categories: Spatial transcriptomics
- Source Repo: Link: https://github.com/FridleyLab/spatialGE
- Contact: Link: https://groups.google.com/u/1/g/genepattern-help
- Programming Language: R


## Input Files (info + type)
| Info | Type |
| --- | --- |
| `input.RDS` | Required |
| `CellComposite.zip` | Optional |

## Output Files (info + type)
| Info | Type |
| --- | --- |
| `spatialGE_STclustered.rds` | Required |

## Parameters (formatted as a Markdown table with name, description, default value, and type)

| **Parameter** | **Name** | **Description** | **Default Value** | **Type** |
| --- | --- | --- | --- | --- |
| `--e (--stlist)` | `-e` | A required option that specifies a STList object in .RDS file format to cluster on. | `None` | `str` |
| `--w (--ws)` | `-w` | An optional parameter that sets the weight to be applied to spatial distances. The default value is 0.025. | `0.025` | `float` |
| `--d (--dist_metric)` | `-d` | An optional parameter that specifies the distance metric to be used. The default value is 'euclidean'. | `'euclidean'` | `str` |
| `--l (--linkage)` | `-l` | An optional parameter that specifies the linkage method to be used for hierarchical clustering. The default value is 'ward.D2'. | `'ward.D2'` | `str` |
| `--k (--ks)` | `-k` | An optional parameter that specifies the range of k values to assess. The default value is 'dtc'. | `2000` | `int` |
| `--t (--topgenes)` | `-t` | An optional integer parameter that sets the number of genes with highest spot-to-spot expression variation. The default value is 2000. | `2000` | `int` |
| `--s (--deepSplit)` | `-s` | An optional logical parameter that controls cluster resolution. The default value is FALSE. | `FALSE` | `bool` |
| `--p (--plot)` | `-p` | An optional logical parameter that enables or disables plotting intermediate results. The default value is FALSE. | `FALSE` | `bool` |
| `--o (--output)` | `-o` | An optional character parameter that specifies the path to save the output STList object. The default value is 'spatialGE_STclustered.rds'. | `'spatialGE_STclustered.rds'` | `str` |
| `--keep_fovs` | `--keep_fovs` | An optional character parameter that specifies a comma-separated list of FOVs (Focus Organizations) to keep. If not specified, no filtering is done. | `None` | `str` |
| `--spot_minreads` | `--spot_minreads` | An optional integer parameter that sets the minimum number of reads per spot. The default value is 20. | `20` | `int` |
| `--rm_genes_expr` | `--rm_genes_expr` | An optional character parameter that specifies a regular expression to filter out genes. The default value is '^NegPrb'. | `^NegPrb` | `str` |
| `--rm_tissue` | `--rm_tissue` | An optional character parameter that specifies a comma-separated list of tissues to remove. If not specified, no filtering is done. | `None` | `str` |
