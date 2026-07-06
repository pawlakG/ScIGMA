#' about_panel UI Function
#'
#' @description A shiny Module.
#'
#' @param id,input,output,session Internal parameters for {shiny}.
#'
#' @noRd
#'
#' @importFrom shiny NS tagList
mod_about_panel_ui <- function(id) {
    ns <- NS(id)
    tagList(
        bslib::card(
            bslib::card_header(
                class = "bg-primary text-white",
                shiny::tagList(shiny::icon("info-circle"), " About ScIGMA")
            ),
            bslib::card_body(
                shiny::markdown("
### ScIGMA: an open-source application for genotype and scADT-seq implementation for scDNA-seq data


Geoffrey Pawlak<sup>1,2,*</sup>, Benjamin Podvin<sup>2,3,*</sup>, Claude Preudhomme<sup>2,3</sup>, Salomon Manier<sup>2,4</sup>, Nicolas Duployez<sup>2,3</sup>, Augustin Boudry<sup>2,3</sup>

<sup>1</sup> **University of Lille**, ULR 2694 METRICS, Lille, France<br>
<sup>2</sup> **PERSTIM**, University of Lille, INSERM, CHU Lille, CNRS, Institut Pasteur of Lille, U1366-UMR9020 - CRCLille - Cancer Research Center of Lille, F-59000 Lille, France<br>
<sup>3</sup> **Laboratory of Hematology**, Lille University Hospital, Lille, France<br>
<sup>4</sup> **Department of Clinical Hematology**, Lille University Hospital, Lille, France<br>
<sup>*</sup> *These authors contributed equally.*

---

#### How to Cite
**Pawlak G., Podvin B., et al.** (202X). *ScIGMA: an open-source application for genotype and scADT-seq implementation for scDNA-seq data.* [Journal Name]. DOI: 10.xxxx/xxxx

---

#### Core Dependencies & Libraries
ScIGMA relies on a robust ecosystem of R packages for single-cell multi-omics integration and visualization.
                "),
                bslib::accordion(
                    id = ns("libraries_accordion"),
                    open = FALSE,
                    bslib::accordion_panel(
                        title = "View detailed package dependencies",
                        icon = shiny::icon("cubes"),
                        shiny::div(
                            style = "max-height: 400px; overflow-y: auto; font-size: 0.85em;",
                            shiny::markdown("
|Package Name|Version|Source|Category|Description|
|:---|:---|:---|:---|:---|
|GenomicRanges|1.60.0|Bioconductor|Single-cell & Genomics|Representation and manipulation of genomic intervals and annotations.|
|IRanges|2.42.0|Bioconductor|Single-cell & Genomics|Vector, list, and range-like data structures for biological sequences.|
|S4Vectors|0.46.0|Bioconductor|Single-cell & Genomics|Low-level S4 vectors and containers for bioconductor packages.|
|SingleCellExperiment|1.30.1|Bioconductor|Single-cell & Genomics|S4 class container for single-cell data storage and manipulation.|
|SummarizedExperiment|1.38.1|Bioconductor|Single-cell & Genomics|Container for matrix-like assays and genomic coordinates annotations.|
|MultiAssayExperiment|1.34.0|Bioconductor|Single-cell & Genomics|Integration of multi-omics data across multiple genomic assays.|
|Seurat|5.4.0|CRAN|Single-cell & Genomics|Comprehensive data analysis and visualization framework for single-cell genomics.|
|EnsDb.Hsapiens.v75|2.99.0|Bioconductor|Single-cell & Genomics|Ensembl-based genomic annotations for Homo sapiens (build GRCh37/hg19).|
|EnsDb.Hsapiens.v86|2.99.0|Bioconductor|Single-cell & Genomics|Ensembl-based genomic annotations for Homo sapiens (build GRCh38/hg38).|
|ensembldb|2.32.0|Bioconductor|Single-cell & Genomics|Utilities to build and query Ensembl-based annotation databases.|
|DelayedArray|0.34.1|Bioconductor|Data Formats & Large Datasets|Delayed operations on array-like datasets for memory-efficient handling.|
|DelayedMatrixStats|1.30.0|Bioconductor|Data Formats & Large Datasets|Fast matrix calculation methods for DelayedMatrix and DelayedArray objects.|
|HDF5Array|1.36.0|Bioconductor|Data Formats & Large Datasets|Read and write HDF5 datasets as DelayedArray representations.|
|rhdf5|2.52.1|Bioconductor|Data Formats & Large Datasets|R interface to the HDF5 high-performance data management library.|
|Matrix|1.7.5|CRAN|Data Formats & Large Datasets|Classes and methods for dense and sparse matrix computations.|
|dplyr|1.2.0|CRAN|Data Manipulation & Utilities|A grammar for data manipulation providing a consistent set of verbs.|
|forcats|1.0.1|CRAN|Data Manipulation & Utilities|Tools for working with categorical variables (factors) in R.|
|purrr|1.2.1|CRAN|Data Manipulation & Utilities|A functional programming toolkit for writing cleaner, robust code.|
|stringr|1.6.0|CRAN|Data Manipulation & Utilities|Consistent, simple wrappers for common string manipulation operations.|
|tibble|3.3.1|CRAN|Data Manipulation & Utilities|A modern reimagining of data frames with lazy evaluation.|
|tidyr|1.3.2|CRAN|Data Manipulation & Utilities|Tools to reshape and tidy messy datasets for analysis.|
|data.table|1.18.2.1|CRAN|Data Manipulation & Utilities|High-performance data manipulation and aggregation with fast syntax.|
|magrittr|2.0.4|CRAN|Data Manipulation & Utilities|Pipes to write cleaner, readable, and chaining-friendly R code.|
|rlang|1.1.7|CRAN|Data Manipulation & Utilities|Tidyverse engine for low-level R programming and metaprogramming.|
|utils|4.5.2|Base R|Data Manipulation & Utilities|Utility functions for system, files, and package management.|
|zip|2.3.3|CRAN|Data Manipulation & Utilities|Cross-platform compression utilities to zip and unzip files.|
|shiny|1.13.0|CRAN|Shiny & Interactive UI|Web application framework for creating interactive dashboards in R.|
|bslib|0.10.0|CRAN|Shiny & Interactive UI|Tools for customizing Bootstrap themes and HTML layouts.|
|DT|0.34.0|CRAN|Shiny & Interactive UI|R wrapper for the DataTables library for interactive tables.|
|gridlayout|0.2.1|CRAN|Shiny & Interactive UI|Layout manager for building structured, grid-based Shiny apps.|
|promises|1.5.0|CRAN|Shiny & Interactive UI|Abstractions for asynchronous programming and reactive systems.|
|shinybusy|0.3.3|CRAN|Shiny & Interactive UI|Busy indicator animations for R Shiny web applications.|
|shinycssloaders|1.1.0|CRAN|Shiny & Interactive UI|CSS loading animations for Shiny outputs during calculations.|
|shinyjs|2.1.1|CRAN|Shiny & Interactive UI|Easily perform common JavaScript operations in Shiny apps.|
|shinyWidgets|0.9.1|CRAN|Shiny & Interactive UI|Custom input widgets and UI controls for Shiny applications.|
|waiter|0.2.5.1|CRAN|Shiny & Interactive UI|Loading screens and progress bars for Shiny applications.|
|golem|0.5.1|CRAN|Shiny & Interactive UI|An opinionated framework for building production-grade Shiny applications.|
|gargoyle|0.0.1|CRAN|Shiny & Interactive UI|An event-driven programming framework for reactive Shiny apps.|
|ggplot2|4.0.2|CRAN|Data Visualization|Declarative graphics system based on the Grammar of Graphics.|
|ggprism|1.0.7|CRAN|Data Visualization|Add-on for ggplot2 to customize plots with Prism-like styling.|
|ggridges|0.5.7|CRAN|Data Visualization|Create ridgeline plots for visualizing density distributions.|
|plotly|4.12.0|CRAN|Data Visualization|Create interactive, web-based quality graphs from ggplot2 objects.|
|scales|1.4.0|CRAN|Data Visualization|Internal scaling infrastructure for mapping data to visual properties.|
|viridis|0.6.5|CRAN|Data Visualization|Colorblind-friendly, perceptually uniform color palettes for visualization.|
|viridisLite|0.4.3|CRAN|Data Visualization|Lightweight implementation of viridis color palettes without dependencies.|
|colorBlindness|0.1.9|CRAN|Data Visualization|Safe color palettes and simulation of color vision deficiencies.|
|ComplexHeatmap|2.24.1|Bioconductor|Data Visualization|Highly customizable, comprehensive heatmap plotting with annotation.|
|InteractiveComplexHeatmap|1.16.0|Bioconductor|Data Visualization|Render ComplexHeatmap plots as interactive widgets in Shiny.|
|circlize|0.4.17|CRAN|Data Visualization|Circular visualization layout engine for multi-omics data mapping.|
|DiagrammeR|1.0.11|CRAN|Data Visualization|Graph and network visualization using Graphviz and Mermaid syntax.|
|DiagrammeRsvg|0.1|CRAN|Data Visualization|Convert DiagrammeR SVG outputs into raw SVG strings.|
|grid|4.5.2|Base R|Data Visualization|Low-level graphics system for plotting primitives and viewports.|
|BiocParallel|1.42.2|Bioconductor|Parallel Computing & Async|Parallel execution framework designed for Bioconductor computations.|
|future|1.70.0|CRAN|Parallel Computing & Async|Unified parallel and distributed processing framework for R.|
|parallel|4.5.2|Base R|Parallel Computing & Async|Support for parallel computation via multicore or socket clusters.|
|config|0.3.2|CRAN|Development & Core|Manage environment-specific configuration files in YAML format.|
|httr2|1.2.2|CRAN|Development & Core|Modern HTTP client for sending requests and parsing API responses.|
|methods|4.5.2|Base R|Development & Core|Formal methods and classes for object-oriented R programming.|
|pkgload|1.5.0|CRAN|Development & Core|Simulate package installation and load packages during development.|
|R6|2.6.1|CRAN|Development & Core|Fast, lightweight object-oriented class system for R.|
|Rcpp|1.1.1|CRAN|Development & Core|Seamless R and C++ integration for high-performance computing.|
|stats|4.5.2|Base R|Development & Core|Standard statistical calculations, distributions, and models.|
|tools|4.5.2|Base R|Development & Core|Internal tools for package development, checks, and processing.|
|whereami|0.2.0|CRAN|Development & Core|Determine the source file and path of active code execution.|
|compositions|2.0.9|CRAN|Development & Core|Mathematical methods for compositional data analysis in R.|
|uwot|0.2.4|CRAN|Development & Core|Fast implementation of Uniform Manifold Approximation and Projection (UMAP).|
|testthat|3.3.2|CRAN|Development & Core|Unit testing framework to ensure package stability and correctness.|
                            ")
                        )
                    )
                ),
                shiny::markdown("
---

#### Contact & Support
For questions, bug reports, or feature requests, please contact:

**Dr. Geoffrey Pawlak**
METRICS, Lille University Hospital
1 Place de Verdun, 59000 Lille, France
**E-mail:** [geoffrey.pawlak@gmail.com](mailto:geoffrey.pawlak@gmail.com)

The source code is publicly available under the GPL-3 license. For bug reports or feature requests, please use our GitHub issue tracker.
**GitHub Repository:** <a href='https://github.com/geoffrey-pawlak/ScIGMA' target='_blank'>https://github.com/geoffrey-pawlak/ScIGMA</a>

---

    #### Application version

**Version:** 1.0.0

---

#### License
This software is released under the **GNU General Public License v3.0 (GPL-3)**. It is an open-source application intended for the scientific community. This 'copyleft' license guarantees your freedom to share and change all versions of the program, and ensures that any derived works must also remain free and open-source.


                ")
            )
        )
    )
}

#' about_panel Server Functions
#'
#' @noRd
mod_about_panel_server <- function(id, ScIGMA_data) {
    moduleServer(id, function(input, output, session) {
        ns <- session$ns
    })
}

## To be copied in the UI
# mod_about_panel_ui("about_panel_1")

## To be copied in the server
# mod_about_panel_server("about_panel_1", ScIGMA_data)
