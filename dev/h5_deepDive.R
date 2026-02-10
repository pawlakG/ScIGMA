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

obj <- loadH5_dir_HDF5("../inputs/bPodvinDatasets/all/", feature_policy = "intersect", omic.type = )
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



obj$seurat_object@meta.data <- cbind(obj$seurat_object@meta.data,
                                             obj$protein.mtx.filtered.normalized[rownames(obj$seurat_object@meta.data),])


umap_df <- cbind(obj$seurat_object@reductions$umap@cell.embeddings,
                 obj$seurat_object@meta.data[,colnames(obj$protein.mtx.filtered.normalized)])

umap_df_long <- umap_df |> pivot_longer(-c(umap_1, umap_2), values_to = "ptn_expression", names_to = "marker")

umap_df_long |> ggplot(aes(x=umap_1, y=umap_2, color=ptn_expression)) +
    geom_point(size = 1) +
    facet_wrap(~marker) +
    scale_color_viridis_c("inferno") +
    ggprism::theme_prism() +
    xlab("UMAP 1") +
    ylab("UMAP 2")





obj_marker_long <- obj$seurat_object@assays$RNA$data |>
    t() |>
    as.data.frame() |>
    rownames_to_column("barcode") |>
    pivot_longer(-barcode, names_to = "marker", values_to = "marker_expression")
obj_marker_long$seurat_cluster <- obj$seurat_object$seurat_clusters[obj_marker_long$barcode]

obj_marker_long |> group_by(marker, seurat_cluster) |> summarize(mean_expression = mean(marker_expression)) |> ggplot() +
    geom_bar( aes(x=marker, y=mean_expression), stat="identity",alpha=0.7) + facet_grid(.~seurat_cluster)

obj_marker_long |> ggplot() +
    geom_boxplot(aes(x=marker, y=marker_expression),alpha=0.7) +
    facet_grid(.~seurat_cluster)

obj_marker_long |> ggplot(aes(x=marker, y=marker_expression)) +
    ggbeeswarm::geom_beeswarm() +
    facet_grid(.~seurat_cluster)



umap_protein <- run_umap_protein(expression_matrix = obj$protein.mtx.filtered.normalized, n_neighbors = 30, min_dist = 0.1)
plot(umap_protein)


plot(obj$protein.mtx.filtered.normalized [,"CD19"],
     obj$protein.mtx.filtered.normalized [,"CD45"])

plan(sequential)
no_parral_norm <- bench::mark({
    NormalizeData(obj_seuratObject, normalization.method = "CLR")
})

avail_cores <- future::availableCores() - 2
plan(multisession, workers = avail_cores)
parral_norm <- bench::mark({
    NormalizeData(obj_seuratObject, normalization.method = "CLR")
})
