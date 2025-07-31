# SingleCellProtExplorer

**Interactive R Shiny application for exploration and analysis of single-cell proteogenomic data**

## Overview

SingleCellProtExplorer is an R Shiny-based application designed to empower researchers and clinicians with an intuitive interface to explore, visualize, and analyze single-cell proteogenomic datasets. The app enables the integration of transcriptomic and proteomic single-cell data, offering a multi-modal perspective on cellular heterogeneity in complex biological systems, such as cancer.

This tool aims to bridge the gap between advanced bioinformatics pipelines and non-computational users, facilitating data-driven discoveries in precision medicine and systems biology.

---

## Features

- 🔬 **Multi-omic integration**: Joint exploration of scRNA-seq and single-cell proteomic data.
- 🧠 **Intuitive user interface**: Designed for ease of use by biologists, clinicians, and researchers.
- 📊 **Interactive visualizations**:
  - UMAP/t-SNE embedding with modality overlay
  - Feature plots for gene/protein expression
  - Cluster and cell-type annotations
- 📈 **Custom analysis modules**:
  - Differential expression analysis
  - Gene/protein co-expression exploration
  - Trajectory inference (optional)
- 🧩 **Modular design**: Easily extensible with new panels and functions.
- 💾 **Support for standard formats**: Compatible with `.h5ad`, `.csv`, and `.rds` input files.

---

## Installation

To install and launch the app locally:

```bash
git clone https://github.com/<your-username>/SingleCellProtExplorer.git
cd SingleCellProtExplorer
Rscript -e "shiny::runApp('.')"
```

**Dependencies**: R ≥ 4.2, and packages listed in `DESCRIPTION` or `renv.lock`. You can restore the environment using `renv`:

```r
renv::restore()
```

---

## Getting Started

1. Prepare your input files:
   - A normalized scRNA-seq matrix (`.rds`, `.csv`, or `.h5ad`)
   - Optional: matching proteomic data (same format, same cells)
   - Optional: metadata table with annotations

2. Launch the application and upload your datasets.

3. Explore clusters, expression profiles, correlations, and more via the app's interface.

---

## Use Cases

This application is particularly suited for:
- Single-cell studies in oncology, immunology, or developmental biology
- Integration of CITE-seq or REAP-seq datasets
- Analysis of tumor microenvironment heterogeneity
- Translational research and biomarker discovery

---

## Citation

If you use this tool in your research, please cite:

> Pawlak G., *et al.* (2025). SingleCellProtExplorer: a modular R Shiny app for interactive exploration of single-cell proteogenomic data. _In preparation_.

---

## License

This project is licensed under the MIT License. See the [LICENSE](./LICENSE) file for details.

---

## Contact

For questions, feature requests, or bug reports, please open an issue or contact:

**Geoffrey Pawlak**  
PharmD, PhD in Bioinformatics & Systems Oncology  
[LinkedIn](https://www.linkedin.com/in/geoffreypawlak) | pawlak.geo [at] pm.me

---

## Future Development

Planned features include:
- Integration with public datasets (e.g., Human Cell Atlas, CPTAC)
- Support for spatial transcriptomics and proteomics
- AI-powered biomarker recommendation
