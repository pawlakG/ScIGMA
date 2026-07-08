# Helper functions for generating mock data for tests

library(rhdf5)
library(Matrix)
library(SummarizedExperiment)
library(MultiAssayExperiment)

#' Create a mock Tapestri-like HDF5 file
#' @param filepath Path to save the mock H5 file
#' @param n_cells Number of cells
#' @param n_vars Number of DNA variants
#' @param n_cnv Number of CNV amplicons
#' @param n_ptn Number of proteins
#' @param include_protein Boolean to include protein assay
#' @return The filepath (invisible)
create_mock_h5 <- function(filepath, n_cells = 10, n_vars = 5, n_cnv = 3, n_ptn = 2, include_protein = TRUE) {
  if (file.exists(filepath)) {
    file.remove(filepath)
  }
  
  h5createFile(filepath)
  
  # Create cell metadata
  h5createGroup(filepath, "assays")
  h5createGroup(filepath, "assays/dna_variants")
  h5createGroup(filepath, "assays/dna_variants/ra")
  h5createGroup(filepath, "assays/dna_variants/ca")
  
  # Cell barcodes
  barcodes <- paste0("Cell_", seq_len(n_cells))
  h5write(barcodes, filepath, "assays/dna_variants/ra/barcode")
  h5write(rep(FALSE, n_cells), filepath, "assays/dna_variants/ra/filtered")
  h5write(rep("Sample1", n_cells), filepath, "assays/dna_variants/ra/sample_name")
  
  # Variants metadata
  var_ids <- paste0("chr1:", 1000 + seq_len(n_vars), ":A:T")
  h5write(var_ids, filepath, "assays/dna_variants/ca/id")
  h5write(rep("chr1", n_vars), filepath, "assays/dna_variants/ca/CHROM")
  h5write(as.character(1000 + seq_len(n_vars)), filepath, "assays/dna_variants/ca/POS")
  h5write(rep("A", n_vars), filepath, "assays/dna_variants/ca/REF")
  h5write(rep("T", n_vars), filepath, "assays/dna_variants/ca/ALT")
  
  # Variant layers
  h5createGroup(filepath, "assays/dna_variants/layers")
  h5write(matrix(runif(n_vars * n_cells, 0, 100), nrow = n_vars, ncol = n_cells), filepath, "assays/dna_variants/layers/AF")
  h5write(matrix(sample(0:3, n_vars * n_cells, replace = TRUE), nrow = n_vars, ncol = n_cells), filepath, "assays/dna_variants/layers/NGT")
  h5write(matrix(sample(10:100, n_vars * n_cells, replace = TRUE), nrow = n_vars, ncol = n_cells), filepath, "assays/dna_variants/layers/DP")
  h5write(matrix(sample(30:99, n_vars * n_cells, replace = TRUE), nrow = n_vars, ncol = n_cells), filepath, "assays/dna_variants/layers/GQ")
  h5write(matrix(0, nrow = n_vars, ncol = n_cells), filepath, "assays/dna_variants/layers/FILTER_MASK")
  
  # Global metadata
  h5createGroup(filepath, "assays/dna_variants/metadata")
  h5write("Mock Version", filepath, "assays/dna_variants/metadata/version")
  
  # CNV metadata
  h5createGroup(filepath, "assays/dna_read_counts")
  h5createGroup(filepath, "assays/dna_read_counts/ca")
  h5createGroup(filepath, "assays/dna_read_counts/ra")
  amp_ids <- paste0("AMPLICON_", seq_len(n_cnv))
  h5write(amp_ids, filepath, "assays/dna_read_counts/ca/id")
  h5write(rep("chr1", n_cnv), filepath, "assays/dna_read_counts/ca/CHROM")
  h5write(as.character(seq_len(n_cnv) * 1000), filepath, "assays/dna_read_counts/ca/start_pos")
  h5write(as.character(seq_len(n_cnv) * 1000 + 100), filepath, "assays/dna_read_counts/ca/end_pos")
  h5write(barcodes, filepath, "assays/dna_read_counts/ra/barcode")
  h5write(rep("Sample1", n_cells), filepath, "assays/dna_read_counts/ra/sample_name")
  
  # CNV layers
  h5createGroup(filepath, "assays/dna_read_counts/layers")
  h5write(matrix(sample(50:500, n_cnv * n_cells, replace = TRUE), nrow = n_cnv, ncol = n_cells), filepath, "assays/dna_read_counts/layers/read_counts")
  
  h5createGroup(filepath, "assays/dna_read_counts/metadata")
  h5write("Mock Version", filepath, "assays/dna_read_counts/metadata/version")
  
  if (include_protein) {
    h5createGroup(filepath, "assays/protein_read_counts")
    h5createGroup(filepath, "assays/protein_read_counts/ca")
    h5createGroup(filepath, "assays/protein_read_counts/ra")
    ptn_ids <- paste0("CD", seq_len(n_ptn))
    h5write(ptn_ids, filepath, "assays/protein_read_counts/ca/id")
    h5write(barcodes, filepath, "assays/protein_read_counts/ra/barcode")
    h5write(rep("Sample1", n_cells), filepath, "assays/protein_read_counts/ra/sample_name")
    
    h5createGroup(filepath, "assays/protein_read_counts/layers")
    h5write(matrix(sample(0:1000, n_ptn * n_cells, replace = TRUE), nrow = n_ptn, ncol = n_cells), filepath, "assays/protein_read_counts/layers/read_counts")
    
    h5createGroup(filepath, "assays/protein_read_counts/metadata")
    h5write("Mock Version", filepath, "assays/protein_read_counts/metadata/version")
  }
  
  invisible(filepath)
}
