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
tmp_heamtap_matrix_filtered_withMissing[tmp_heamtap_matrix_filtered_withMissing == 3] <- NA

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

## rbind with samples with missing info
tmp_heamtap_matrix_filtered_complete_ordered <- rbind(tmp_heamtap_matrix_filtered_noMissing_ordered,
                                                      tmp_heamtap_matrix_filtered_withMissing)

## Set labels / annotations
heatmap_split_vector <- c(as.character(hc_sample_noMissing),
                          rep("missing", nrow(tmp_heamtap_matrix_filtered_withMissing)))

## plot heatmap
library(ComplexHeatmap)
library(colorBlindness)
### Color palette
dna_variant_colorPalette <- setNames(c("#E9E8EC", "#BAB7D0", "#3C2692"), nm = c("0","1","2"))
### Legend
dna_variant_legend <- Legend(at = c("0","1","2"), labels = c("WT", "HET", "HOM"),
                             title = "Genotype", legend_gp = gpar(fill = c("#E9E8EC", "#BAB7D0", "#3C2692")))
### Annotation
heatmap_true_levels <- levels(hc_sample_noMissing)[levels(hc_sample_noMissing)!="small"]
annotationColor <- list(Cluster = setNames(c(colorBlindness::paletteMartin[-1][1:length(heatmap_true_levels)],"grey","#333333"),
                                           nm = c(heatmap_true_levels, "missing","small")))
dna_variant_annotation <- rowAnnotation(Cluster = heatmap_split_vector,
                                        col =annotationColor,
                                        show_annotation_name = FALSE,
                                        show_legend = FALSE)

### Heatmap
dna_variant_heatmap <- Heatmap(tmp_heamtap_matrix_filtered_complete_ordered,
                               column_names_rot = 70,
                               row_split = heatmap_split_vector, show_column_dend = FALSE,
                               cluster_rows = FALSE, show_row_names = FALSE,
                               na_col = "black", col = dna_variant_colorPalette,
                               show_heatmap_legend = FALSE, left_annotation = dna_variant_annotation)




