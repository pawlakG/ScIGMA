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


# [ NODE_ACCESS : Test MAE pipeline ]
# ----------------------------------------------------- _
devtools::load_all()
ScIGMA_data <- loadH5_HDF5_biocond(
    filepath = directory,
    sample_name = "aml_4_lines",
    omic_type = "DNA+protein"
)
ScIGMA_data <- normalizeProtein(ScIGMA_data)
ScIGMA_data$seurat_object <- protein_run_pca(ScIGMA_data)

ScIGMA_data <- filter_and_annotate_variants(ScIGMA_data, paths = cfg$paths)

SummarizedExperiment::rowData(ScIGMA_data$mae[["dna_variants"]]) |>
    as.data.frame() |>
    dplyr::select(variant_id, gene, variant_type, gene_function, impact, clinvar, cell_proportion) |>
    arrange(desc(cell_proportion), desc(impact))




ScIGMA_data$seurat_object <- RunUMAP(ScIGMA_data$seurat_object,
                                     dims = 1:(nrow(ScIGMA_data$seurat_object)-2),
                                     min.dist = 0.15,
                                     n.neighbors = 30,
                                     future.seed=TRUE)
