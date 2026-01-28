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



obj$protein.mtx.filtered.normalized <- normalize_linear_regression(as.matrix(obj$protein.mtx), jitter = 0.5)

# Initiate seurat object
library(Seurat)
obj_seuratObject <- CreateSeuratObject(counts = t(obj$protein.mtx.filtered.normalized) ,
                   project = "ScIGMA_data",
                   min.cells = 3,
                   min.features = floor(sqrt(ncol(obj$protein.mtx.filtered.normalized))))

# obj_seuratObject <- CreateSeuratObject(counts = as.matrix(t(obj$protein.mtx.filtered)),
#                                        project = "ScIGMA_data",
#                                        min.cells = 3,
#                                        min.features = floor(sqrt(ncol(obj$protein.mtx.filtered))))


# obj_seuratObject <- NormalizeData(obj_seuratObject, normalization.method = "CLR")
obj_seuratObject@assays$RNA$data <- normalize_linear_regression(as.matrix(t(obj$protein.mtx.filtered)), jitter = 0.1)

obj_seuratObject <- FindVariableFeatures(obj_seuratObject,
                                         selection.method = "vst",
                                         nfeatures = nrow(obj_seuratObject))

obj_seuratObject <- ScaleData(obj_seuratObject, features = rownames(obj_seuratObject))
obj_seuratObject <- RunPCA(obj_seuratObject,
                           features = VariableFeatures(object = obj_seuratObject),
                           npcs = nrow(obj_seuratObject)-2)

obj_seuratObject <- FindNeighbors(obj_seuratObject, dims = 1:(nrow(obj_seuratObject)-2))
obj_seuratObject <- FindClusters(obj_seuratObject, resolution = 0.15, future.rng.onMisuse = "ignore" )

obj_seuratObject <- RunUMAP(obj_seuratObject,
                            dims = 1:(nrow(obj_seuratObject)-2),
                            min.dist = 0.2,
                            n.neighbors = 30)
obj$seurat_object <- obj_seuratObject

UMAPPlot(obj_seuratObject)

obj_seuratObject_markers <- FindAllMarkers(obj_seuratObject, only.pos = TRUE)



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
