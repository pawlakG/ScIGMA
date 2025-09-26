devtools::load_all()

directory <- "../inputs/tapestriDatasets/4-cell-lines-AML-multiomics/4-cell-lines-AML-multiomics.dna+protein.h5"
directory <- "../inputs/bPodvinDatasets/Jas_Rec/Jas_rec.dna+protein.h5"
directory <- "../inputs/tapestriDatasets/2_PBMC_Mix_KG-1_Spike_In/Sample_Data_Set_2_PBMC_Mix_KG-1_Spike_In.labeled.dna+protein.h5"
directory <- "../inputs/tapestriDatasets/KG-1-Raji-50-50-Myeloid/KG-1-Raji-50-50-Myeloid.dna.h5"
directory <- "../inputs/bPodvinDatasets/Dut_Ev/Dut_Ev_Rec.dna+protein.h5"

rhdf5::h5ls(directory, all = TRUE, recursive = TRUE)

h5f <- H5Fopen(directory, flags = "H5F_ACC_RDONLY")

h5f$assays$dna_variants$layers$NGT |> table()

h5f$assays$dna_read_counts$ca

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
library(ggridges)
library(tidyr)
# ---------------------------- #
# Normalization (CLR)
obj <- normalizeProtein(obj)
obj$protein.normalize.method
dim(obj$protein.mtx.filtered.normalized)
View(obj$protein.mtx.filteed.normalized)
# ---------------------------- #
# Ridge plot
protein_plot <- obj$protein.mtx.filtered.normalized |> as_tibble() |> pivot_longer(everything()) |>
    ggplot(aes(x=value, y=name, fill = name)) +
    geom_density_ridges() +
    theme_ridges() +
    theme(legend.position = "none")
protein_plot |> plotly::ggplotly()

# ---------------------------- #
# barplot of relative percentage of counts
protein_rel_counts <- colSums(obj$protein.mtx.filtered) / sum(obj$protein.mtx.filtered |> colSums())

protein_rel_counts_tb <- tibble( protein = names(protein_rel_counts), percent = protein_rel_counts)

protein_rel_percent_barplot <- protein_rel_counts_tb |> ggplot(aes(x=reorder(protein, -percent), y=percent, fill=protein)) +
    geom_bar(stat="identity") +
    theme_minimal() +
    theme(legend.position = "none") +
    ylab("Relative percentage of counts") +
    xlab("Protein") +
    scale_y_continuous(labels = scales::percent_format(accuracy = 1))

protein_rel_percent_barplot |> plotly::ggplotly()

# ---------------------------- #
# Scatter plot of two proteins
protein_scatter <- obj$protein.mtx.filtered.normalized |> as_tibble() |>
    ggplot(aes(x=`CD19`, y=`CD45`)) +
    geom_point(alpha=0.5) +
    theme_minimal() +
    xlab("CD19 normalized counts") +
    ylab("CD45 normalized counts")



obj$protein.mtx.filtered.normalized |> as_tibble() |>
    ggplot(aes(x=`CD19`, y=`CD45`)) +
    geom_point(alpha=0.5) +
    theme_minimal() +
    xlab("CD19 normalized counts") +
    ylab("CD45 normalized counts")

