devtools::load_all()

directory <- "../inputs/tapestriDatasets/4-cell-lines-AML-multiomics/4-cell-lines-AML-multiomics.dna+protein.h5"
directory <- "../inputs/bPodvinDatasets/Jas_Rec/Jas_rec.dna+protein.h5"
directory <- "../inputs/tapestriDatasets/2_PBMC_Mix_KG-1_Spike_In/Sample_Data_Set_2_PBMC_Mix_KG-1_Spike_In.labeled.dna+protein.h5"
directory <- "../inputs/tapestriDatasets/KG-1-Raji-50-50-Myeloid/KG-1-Raji-50-50-Myeloid.dna.h5"
directory <- "../inputs/bPodvinDatasets/Dut_Ev/Dut_Ev_Rec.dna+protein.h5"

rhdf5::h5ls(directory, all = TRUE, recursive = TRUE)

h5f <- H5Fopen(directory, flags = "H5F_ACC_RDONLY")

h5f$assays$dna_variants$layers$NGT |> table()



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

# get number cell filtered per variants

cols <-sub(x = variants_annotated$variant_id, pattern = "^([^:]+:)|^:", "")

apply(obj$vaf.mtx[,1:50], 2, \(x) sum(x == 0)) # Get number of unmutated cells

nrow(obj$dna.variant.filter.mask.filtered) - realize(obj$dna.variant.filter.mask.filtered[, cols] != as.raw(0)) |> colSums()

obj$dna.variant.filter.mask.filtered


# variant heatmap
## Select some variants
# selected_variants <- cols[1:3]
selected_variants <- cols[-9]
tmp_heamtap_matrix <- obj$gt.mtx[,selected_variants] |> as.matrix()
ngt_matrix_selected <- obj$dna.variant.filter.mask.filtered[, selected_variants] |> as.matrix() != 0

## Apply NGT_MASK matrix
tmp_heamtap_matrix_filtered <- tmp_heamtap_matrix[intersect(rownames(tmp_heamtap_matrix),
                                                            rownames(obj$dna.variant.filter.mask.filtered)),]



tmp_heamtap_matrix_filtered[cbind(row(ngt_matrix_selected)[ngt_matrix_selected],
                                  col(ngt_matrix_selected)[ngt_matrix_selected])] <- 3


## Take a subset where no row has values equal to 3
tmp_heamtap_matrix_filtered_noMissing <- tmp_heamtap_matrix_filtered[rowSums(tmp_heamtap_matrix_filtered == 3) == 0,]

## Take rows with any values equal to 3
tmp_heamtap_matrix_filtered_withMissing <- tmp_heamtap_matrix_filtered[rowSums(tmp_heamtap_matrix_filtered == 3) > 0,]

## performe clustering on matrix with no value equal to 3
### Hamming distance
sample_dist_matrix <- proxy::dist(tmp_heamtap_matrix_filtered_noMissing, method = function(x,y) sum(x != y)) |>
    as.matrix()
variant_dist_matrix <- proxy::dist(t(tmp_heamtap_matrix_filtered_noMissing), method = function(x,y) sum(x != y)) |>
    as.matrix()

### reorder rows and columns
#### Get clusters
hc_sample_noMissing <- hclust(as.dist(sample_dist_matrix), method = "ward.D2") |> cutree(k = sqrt(nrow(sample_dist_matrix))) |> sort()
hc_variant_noMissing <- hclust(as.dist(variant_dist_matrix), method = "ward.D2")|> cutree(k = 6) |> sort()
#### Reassign too small samples clusters (compared to total samples) to a "small" group
res_table_clusters <- table(hc_sample_noMissing)
too_small_clusters <- names(res_table_clusters)[res_table_clusters/nrow(tmp_heamtap_matrix) < 0.01]
hc_sample_noMissing[hc_sample_noMissing %in% too_small_clusters] <- "small"
hc_sample_noMissing <- hc_sample_noMissing |> sort() |> as.factor()
#### Reorder
tmp_heamtap_matrix_filtered_noMissing_ordered <- tmp_heamtap_matrix_filtered_noMissing[names(hc_sample_noMissing),names(hc_variant_noMissing)]














## plot heatmap

# fig <- plotly::plot_ly(z = tmp_heamtap_matrix_filtered_noMissing_ordered, type = "heatmap")
fig <- heatmaply::heatmaply(tmp_heamtap_matrix_filtered_noMissing_ordered, scale = "none",
                            Rowv = FALSE, Colv = FALSE)
fig <- heatmaply::heatmaply(tmp_heamtap_matrix_filtered_noMissing_ordered)

fig
