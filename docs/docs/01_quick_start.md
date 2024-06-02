---
layout: article
title: Installation
sidebar:
  nav: document_nav
permalink: docs/01_quick_start
---

StarTail is implemented in R.

### Dependencies 
``` 
* R version >= 3.6.0.
* Dependent R packages: 
    - data.table
    - spNNGP
    - stats
```

### Installation

To install StarTrail, you can either install from GitHub or directly from source file. 

  <a class="button button--rounded button--sm" style= "background-color: #F285AD; color: white;" href="https://github.com/JiawenChenn/StarTrail"><i class="fas fa-download"></i> View on GitHub</a>
  <a class="button button--rounded button--sm" style= "background-color: #4B7BA6; color: white;" href="doc_data/StarTrail_1.0.0.tar.gz"><i class="fas fa-download"></i> Download source</a>

```r
# method1: install from GitHub
require('devtools')
devtools::install_github('JiawenChennn/StarTrail')

# method2: install from source file
install.packages('StarTrail_1.0.0.tar.gz',type="source",repos = NULL)
```

### Optional package for visualization
```r
library(dplyr)
library(data.table)
library(ggplot2)
library(cowplot)
library(ggsci)
library(egg)
library(paletteer)
```