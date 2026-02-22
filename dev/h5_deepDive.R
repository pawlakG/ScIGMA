devtools::load_all()

directory <- "../inputs/tapestriDatasets/4-cell-lines-AML-multiomics/4-cell-lines-AML-multiomics.dna+protein.h5"
directory <- "../inputs/bPodvinDatasets/Jas_Rec/Jas_rec.dna+protein.h5"
directory <- "../inputs/tapestriDatasets/2_PBMC_Mix_KG-1_Spike_In/Sample_Data_Set_2_PBMC_Mix_KG-1_Spike_In.labeled.dna+protein.h5"
directory <- "../inputs/tapestriDatasets/KG-1-Raji-50-50-Myeloid/KG-1-Raji-50-50-Myeloid.dna.h5"
directory <- "../inputs/bPodvinDatasets/Dut_Ev/Dut_Ev_Rec.dna+protein.h5"
directory <- "../inputs/bPodvinDatasets/dew_diag.h5"

rhdf5::h5ls(directory, all = TRUE, recursive = TRUE)

h5f <- H5Fopen(directory, flags = "H5F_ACC_RDONLY")

h5f$assays$dna_variants$layers$NGT |> table()

h5f$assays$dna_read_counts$ca

dim(h5f$assays$protein_read_counts$layers$read_counts)

h5f$assays$protein_read_counts$ra$sample_name


h5closeAll()


# --------------------------------------------------------------- #
#

obj <- loadH5_HDF5(directory, sample.name = "SampleA", omic.type = "DNA+protein")


obj <- loadH5_HDF5(directory, sample.name = "SampleA", omic.type = "DNA")

obj <- loadH5_dir_HDF5(directory, omic.type = "DNA+protein",feature_policy = "intersect")

bench::mark(loadH5_HDF5(directory, sample.name = "SampleA", omic.type = "DNA+protein"), iterations = 10)

file.exists("../inputs/bPodvinDatasets/")

# --------------------------------------------------------------- #
# More than one  h5

# obj <- loadH5_dir_HDF5("../inputs/bPodvinDatasets/all/", feature_policy = "intersect", omic.type = )
# obj$realize_all(dir = "store", file = "ScIGMA_merged.h5", chunkdim = c(1024,512), level = 6)

# --------------------------------------------------------------- #
# DEBUG

library(bench)

obj <- filter_variant_ScIGMA(obj = obj,
                             min.dp = 10,
                             min.gq = 30,
                             vaf.ref = 5,
                             vaf.hom = 95,
                             vaf.het = 30,
                             min.cell.pt = 10,
                             min.mut.cell.pt = 10)

obj$variant.annotation <- tryCatch(
    fetch_variants_batch_fields(obj$variants.filtered,
                                batch_size = 300,
                                paths = cfg$paths)
    , error = function(e){
        remove_modal_spinner()
        message(warning("Error during variant annotation: "),
                stop(e$message))
    })
# Add info about proportion of mutated cells per variants
obj$variant.annotation$probe <- gsub("^[^:]*:", "", obj$variant.annotation$variant_id)
obj$variant.annotation$cell_proportion <- apply(as.matrix(obj$vaf.mtx.filtered)[,obj$variant.annotation$probe], 2, \(x){
    sum(x != 0) / nrow(obj$vaf.mtx.filtered)
})


obj <- protein_run_pca(obj)


obj$seurat_object <- RunUMAP(obj$seurat_object,
                                     dims = 1:(nrow(obj$seurat_object)-2),
                                     min.dist = 0.15,
                                     n.neighbors = 30,
                                     future.seed=TRUE)



# [ NODE_ACCESS : TRY BITSC2 ]
# ----------------------------------------------------- _



library(hdf5r)

library(rhdf5)

# prepare_bitsc2_input: Extracts and transforms H5 layers for BiTSC2
prepare_bitsc2_input <- function(file_path) {
    # Define internal paths within the Tapestri H5 structure
    assay_path <- "/assays/dna_variants/"
    layer_path <- paste0(assay_path, "layers/")

    # Read AF (Allele Frequency) and DP (Total Depth)
    # rhdf5 reads dimensions in the order they are stored in H5
    af_matrix <- rhdf5::h5read(file_path, paste0(layer_path, "AF"))
    af_matrix_init <- af_matrix
    af_matrix <- af_matrix/100
    dp_matrix <- rhdf5::h5read(file_path, paste0(layer_path, "DP"))

    # Retrieve metadata for identification
    variant_ids <- rhdf5::h5read(file_path, paste0(assay_path, "ca/id"))
    cell_barcodes <- rhdf5::h5read(file_path, paste0(assay_path, "ra/barcode"))

    # Close all open H5 handles to prevent file locking
    rhdf5::h5closeAll()

    # Reconstruct alternate allele counts: A = AF * DP
    # We use round() to ensure integer values for the binomial likelihood
    alt_counts <- round(af_matrix * dp_matrix)

    # BiTSC2 requires the matrix shape: variants (rows) x cells (columns)
    # Note: verify dimensions as rhdf5 might transpose depending on H5 layout
    if (nrow(alt_counts) != length(variant_ids)) {
        alt_counts <- t(alt_counts)
        dp_matrix <- t(dp_matrix)
    }

    rownames(alt_counts) <- variant_ids
    colnames(alt_counts) <- cell_barcodes
    rownames(dp_matrix) <- variant_ids
    colnames(dp_matrix) <- cell_barcodes

    return(list(
        af_matrix_init = af_matrix_init,
        alternate_counts = alt_counts,
        total_depth = dp_matrix
    ))
}


library(rhdf5)

# extract_clean_bitsc2_data: Filters out low-quality genotypes before BiTSC2
extract_clean_bitsc2_data <- function(file_path, gq_threshold = 30) {
    # Load layers using rhdf5
    assay_path <- "/assays/dna_variants/"
    layer_path <- paste0(assay_path, "layers/")
    af_matrix <- rhdf5::h5read(file_path, "/assays/dna_variants/layers/AF")
    af_matrix <- af_matrix/100
    dp_matrix <- rhdf5::h5read(file_path, "/assays/dna_variants/layers/DP")
    gq_matrix <- rhdf5::h5read(file_path, "/assays/dna_variants/layers/GQ")

    # Retrieve metadata for identification
    variant_ids <- rhdf5::h5read(file_path, paste0(assay_path, "ca/id"))
    cell_barcodes <- rhdf5::h5read(file_path, paste0(assay_path, "ra/barcode"))

    # Close all open H5 handles to prevent file locking
    rhdf5::h5closeAll()

    # Reconstruct alternate allele counts: A = AF * DP
    # We use round() to ensure integer values for the binomial likelihood
    alt_counts <- round(af_matrix * dp_matrix)

    # BiTSC2 requires the matrix shape: variants (rows) x cells (columns)
    # Note: verify dimensions as rhdf5 might transpose depending on H5 layout
    if (nrow(alt_counts) != length(variant_ids)) {
        alt_counts <- t(alt_counts)
        dp_matrix <- t(dp_matrix)
    }

    rownames(alt_counts) <- variant_ids
    colnames(alt_counts) <- cell_barcodes
    rownames(dp_matrix) <- variant_ids
    colnames(dp_matrix) <- cell_barcodes

    # Identify low quality calls (GQ < threshold)
    # These should be treated as 'Missing' (NA) rather than Wild Type (0)
    low_quality_mask <- gq_matrix < gq_threshold

    # Apply mask: setting AF to NA for low quality entries
    # This forces BiTSC2 to rely on the clonal structure (imputation)
    af_matrix[low_quality_mask] <- NA

    # Final count reconstruction
    # A = round(AF * DP)
    alt_counts <- round(af_matrix * dp_matrix)

    rhdf5::h5closeAll()

    return(list(
        A = t(alt_counts),
        D = t(dp_matrix)
    ))
}

test <- extract_clean_bitsc2_data(directory)

saveRDS(test, "../../extract_clean_bitsc2_data.rds")
