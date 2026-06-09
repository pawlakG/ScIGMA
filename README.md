# ScIGMA: Single-cell Integrated Genomic & Multi-omics Analyzer

**An interactive computational framework for the exploration and integrative analysis of single-cell proteogenomic data.**

## Abstract

ScIGMA (Single-cell Integrated Genomic & Multi-omics Analyzer) is an advanced R-based application designed to process, integrate, and analyze multi-modal single-cell datasets. Built for precision medicine and systems biology, the framework specializes in the simultaneous evaluation of Single Nucleotide Variants (SNVs), Copy Number Variations (CNVs), and targeted cell-surface protein expression. By resolving the technical noise inherent to single-cell sequencing and coupling genomic phylogenies with phenotypic states, ScIGMA provides a comprehensive mapping of cellular heterogeneity and clonal evolution.

## Core Capabilities and Architecture

### 1. High-Fidelity Data Ingestion
- **Matrix Processing**: Native extraction and parsing of multi-assay hierarchical data formats (HDF5 / `.h5`), highly optimized for architectures such as the Mission Bio Tapestri platform.
- **Rigorous Filtering**: Dynamic, threshold-based exclusion of low-quality cells, doublets, and uninformative variants prior to downstream statistical modeling.

### 2. Probabilistic Clonal Inference
- **MCMC Genotype Imputation**: ScIGMA implements a native `Rcpp` interface to the COMPASS C++ backend. This integration allows for rigorous Markov Chain Monte Carlo (MCMC) inference to resolve allele dropouts and technical missingness.
- **Phylogenetic Reconstruction**: Joint probabilistic modeling of SNVs and CNVs to construct robust clonal architectures, returning objective, data-driven single-cell phylogenies.

### 3. Automated Variant Annotation
- **Clinical and Biological Mapping**: Integrated RESTful communication with the Ensembl Variant Effect Predictor (VEP) API. 
- **Consequence Stratification**: Strict extraction of the most severe biological consequence (e.g., missense, splice-site) and mapping of pathogenic status via ClinVar, resolving complex multi-allelic annotations into singular, clinically relevant metrics.

### 4. Proteogenomic Dimensionality Reduction
- **Protein Space Mapping**: Utilization of UMAP algorithms for the dimensionality reduction of targeted protein assays.
- **Integrative Projection**: Capability to cross-project inferred phylogenetic clones directly onto immunophenotypic spatial coordinates, quantifying the genotype-to-phenotype continuous transitions.

### 5. System Infrastructure
- **Object-Oriented Design**: Built upon an encapsulated `R6` object system, ensuring immutable raw data states and strict version control of imputed layers.
- **Memory Optimization**: Leverages `SummarizedExperiment` and `MultiAssayExperiment` classes from Bioconductor for memory-efficient matrix operations and metadata handling.

## Installation

ScIGMA is distributed as a self-contained R package containing a native C++ backend (COMPASS). Installation requires R version >= 4.2.0 and a standard C++ compiler appropriate for your operating system.

### 1. System Requirements (Compiler Toolchain)
Before installing the package, ensure your machine has the necessary compiler tools to build the C++ objects:
- **Windows**: Install [Rtools](https://cran.r-project.org/bin/windows/Rtools/) (ensure the version matches your R version).
- **macOS**: Open Terminal and run `xcode-select --install` to install the Command Line Tools.
- **Linux**: Install the build-essential tools via terminal: `sudo apt-get install build-essential r-base-dev`.

### 2. Package Installation
Install ScIGMA directly from the GitHub repository. The `remotes` package will automatically download the source code, compile the C++ binaries locally for your specific CPU architecture, and resolve all R dependencies.

```r
if (!require("remotes")) install.packages("remotes")
remotes::install_github("pawlakG/ScIGMA", dependencies = TRUE)
```

### 3. Launching the Application
Once the package is installed, the Shiny interface and the underlying analytic engines are fully encapsulated. Launch the environment using the native package entry point:

```r
library(ScIGMA)
ScIGMA::run_app()
```

## Protocol Overview

1. **Upload**: Load raw `.h5` assay matrices containing DNA variants and corresponding protein expression data.
2. **Pre-processing**: Define strict filtering parameters (e.g., minimum depth, VAF thresholds) to clean the mutational matrices.
3. **Inference**: Execute the COMPASS MCMC algorithms to impute missing genotypes and extract the clonal tree.
4. **Annotation**: Trigger the API module to functionally annotate selected variants against Ensembl and ClinVar.
5. **Integration**: Reduce protein dimensions, define phenotypic subpopulations, and compute the mathematical intersection between genomic clones and protein profiles.

## Target Domains

This computational tool is specifically tailored for:
- Systems Oncology and Hematology (e.g., Multiple Myeloma, Leukemias).
- Dissection of the tumor microenvironment and intra-tumoral heterogeneity.
- Tracking of sub-clonal resistance mechanisms during therapeutic pressure.
- High-resolution biomarker discovery via paired DNA/Protein analysis.

## Citation

If you utilize ScIGMA in your research or analytical pipelines, please cite:

> Pawlak G., *et al.* (2026). ScIGMA: A computational framework for the interactive exploration of single-cell proteogenomic architecture. *In preparation*.

## License

This project is licensed under the GNU General Public License v3.0 (GPLv3). See the `LICENSE` file for full terms and conditions.

## Contact

For academic inquiries, methodological questions, or system bug reports:

**Geoffrey Pawlak**
PharmD, PhD in Bioinformatics & Systems Oncology
pawlak.geo [at] pm.me
