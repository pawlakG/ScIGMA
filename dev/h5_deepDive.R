library(rhdf5)

directory <- "../inputs/tapestriDatasets/4-cell-lines-AML-multiomics/4-cell-lines-AML-multiomics.dna+protein.h5"
directory <- "../inputs/bPodvinDatasets/Jas_Rec/Jas_rec.dna+protein.h5"
directory <- "../inputs/tapestriDatasets/2_PBMC_Mix_KG-1_Spike_In/Sample_Data_Set_2_PBMC_Mix_KG-1_Spike_In.labeled.dna+protein.h5"
directory <- "../inputs/tapestriDatasets/KG-1-Raji-50-50-Myeloid/KG-1-Raji-50-50-Myeloid.dna.h5"
directory <- "../inputs/bPodvinDatasets/Dut_Ev/Dut_Ev_Rec.dna+protein.h5"

rhdf5::h5ls(directory, all = TRUE, recursive = TRUE)

h5f <- H5Fopen(directory, flags = "H5F_ACC_RDONLY")
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

obj <- filter_variant_ScIGMA(obj = obj, min.dp = 10, min.gq = 30, vaf.ref = 5, vaf.hom = 95, 
vaf.het = 30, min.cell.pt = 50, min.mut.cell.pt = 50)


dim(obj$dna.variant.filter.mask)
dim(obj$dna.variant.filter.mask.filtered)

variants_annotated <- fetch_variants_batch_fields(obj$variants.filtered, paths = paths)

# get number cell filtered per variants

cols <-sub(x = variants_annotated$variant_id, pattern = "^([^:]+:)", "")

apply(obj$vaf.mtx[,1:50], 2, \(x) sum(x == 0)) # Get number of unmutated cells

nrow(obj$dna.variant.filter.mask.filtered) - realize(obj$dna.variant.filter.mask.filtered[, cols] != as.raw(0)) |> colSums()

obj$dna.variant.filter.mask.filtered
