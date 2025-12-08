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
                             min.cell.pt = 50,
                             min.mut.cell.pt = 50)

dim(obj$dna.variant.filter.mask)
dim(obj$dna.variant.filter.mask.filtered)

variants_annotated <- fetch_variants_batch_fields(obj$variants.filtered, paths = paths)


# all_variants <- fetch_variants_batch_fields(obj$variants, paths = paths)

# get number cell filtered per variants

cols <- sub(x = variants_annotated$variant_id, pattern = "^([^:]+:)|^:", "")

apply(obj$vaf.mtx[,1:50], 2, \(x) sum(x == 0)) # Get number of unmutated cells

nrow(obj$dna.variant.filter.mask.filtered) - realize(obj$dna.variant.filter.mask.filtered[, cols] != as.raw(0)) |> colSums()

obj$dna.variant.filter.mask.filtered


# variant heatmap
## Select some variants
selected_variants <- cols[1:3]
# selected_variants <- cols[-9]


# test <- generate_dna_variant_heatmap(obj = obj, selected_variants = selected_variants, n_cluster = 3)
test <- generate_dna_variant_heatmap(obj = obj, selected_variants_df = variants_annotated[1:5,], n_cluster = 3)



# --------------------------------------------------------------- #
# CNV
## Compute ploidy

### After filtering
dna_variant_no_missing <- get_no_missing_dna_variant(obj, selected_variants_df = variants_annotated[1:2,])
sum(rowSums2(dna_variant_no_missing) == 0)




# --------------------------------------------------------------- #
# Protein analysis
library(ggplot2)
library(ggprism)
library(ggridges)
library(tidyr)
library(ADTnorm)
library(tibble)
# --------------------------------------------------------------- #
# Test ADTnorm

tmp_ptn_mtx <- as.matrix(obj$protein.mtx) # 1313 10
tmp_cell_mtx <- data.frame(cell_barcode = obj$cell.ids, sample = "4CL_AML")

test_ADT <- ADTnorm(cell_x_adt = tmp_ptn_mtx, cell_x_feature = tmp_cell_mtx, save_outpath = "test_adt/")
test_ADT_tmp <- as.data.frame(test_ADT) |> rownames_to_column("cell")
test_ADT_long <- test_ADT_tmp |> pivot_longer(names_to = "marker", values_to = "value", cols = -cell)

ggplot(test_ADT_long, aes(x = value, y = marker)) + geom_density_ridges()

plot(test_ADT[,13], test_ADT[,6])

# ---------------------------- #
# Normalization (CLR)
obj <- normalizeProtein(obj)
# ---------------------------- #
# Barplot

plot_protein_barplot(obj)

# ---------------------------- #
# Biplot
marker_1 <- "CD33"
marker_2 <- "CD45"

tmp_ptn_df <- obj$protein.mtx.filtered.normalized[,c(marker_1, marker_2)]
# |> as.data.frame()

tmp_ptn_df_norm <- apply(tmp_ptn_df, 2, scale) |> as.data.frame()

tmp_ptn_df_norm |>
    ggplot(aes_string(x = marker_1, y = marker_2)) +
    geom_point() + theme_prism()

plot(tmp_ptn_df[,1], tmp_ptn_df[,2])
plot(tmp_ptn_df_norm)


# --------------------------------------------------------------- #
# Test repreoduction of NSP algorithm
# counts : matrice (cells x proteins) de lectures brutes
tmp_mat <- obj$protein.mtx |> as.matrix()
res <- nsp_transform(
    counts_mat           = tmp_mat,
    jitter               = 2,        # doc NSP
    scale                = NULL,       # auto-scale si nécessaire
    sample_size          = 1000,       # ANSP-like si très gros n
    random_state         = 42,
    p_low                = 0.1,
    p_high               = 0.9,
    max_zero_read_cells  = 0.01
)

nsp_counts <- res$normalized
str(res$models[[1]])  # coefficients fond/signal de la 1re protéine
res$scaling_factor


marker_1 <- "CD19"
marker_2 <- "CD45"

tmp_df <- res$normalized[,c(marker_1, marker_2)]

tmp_df_norm <- apply(tmp_df, 2, log2)  |> as.data.frame()
# tmp_df_norm <- tmp_df

tmp_df_norm |>
    ggplot(aes_string(x = marker_1, y = marker_2)) +
    geom_point() + theme_prism()
