# infercnvpip

<!-- badges: start -->
[![R-CMD-check](https://github.com/hungms/infercnvpip/actions/workflows/R-CMD-check/badge.svg)](https://github.com/hungms/infercnvpip/actions/workflows/R-CMD-check.yaml)
<!-- badges: end -->

## Overview

`infercnvpip` is a streamlined interface to the [InferCNV](https://github.com/broadinstitute/infercnv) package for analyzing copy number variations in single-cell RNA-seq data. It simplifies the workflow for running InferCNV on [Seurat](https://satijalab.org/seurat/) objects.

## Installation

You can install the development version of infercnvpip from GitHub with:

```r
# install.packages("devtools")
devtools::install_github("hungms/infercnvpip")
```

## Features

- Simple interface for running InferCNV on Seurat objects
- Support for both individual and batch processing of samples
- Human and mouse genome annotations included
- Integration with Seurat objects

## Usage

```r
library(infercnvpip)

# Run InferCNV on a single object
run_infercnv_individual(
  obj = seurat_obj,
  ref_obj = ref_obj,
  individual_name = "sample1",
  annot_column = "celltype",
  org = "human"
)

# Run InferCNV on multiple objects
run_infercnv_multi(
  obj = seurat_obj,
  split.by = "sample", 
  ref_obj = ref_obj,
  annot_column = "celltype",
  org = "human"
)
```

## Documentation

For more detailed documentation, see the package vignettes:

```r
browseVignettes("infercnvpip")
```


