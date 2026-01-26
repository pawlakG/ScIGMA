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

obj$protein.mtx.filtered.normalized <- normalize_linear_regression(as.matrix(obj$protein.mtx), jitter = 0.5)



umap_protein <- run_umap_protein(expression_matrix = obj$protein.mtx.filtered.normalized, n_neighbors = 30, min_dist = 0.1)
plot(umap_protein)


plot(obj$protein.mtx.filtered.normalized [,"CD19"],
     obj$protein.mtx.filtered.normalized [,"CD45"])

render_protein_ridge_plot(obj)



variants_annotated <- fetch_variants_batch_fields(obj$variants.filtered, paths = paths)

# all_variants <- fetch_variants_batch_fields(obj$variants, paths = paths)

# get number cell filtered per variants

cols <- sub(x = variants_annotated$variant_id, pattern = "^([^:]+:)|^:", "")

# variant heatmap
## Select some variants
selected_variants <- cols[1:4]
# selected_variants <- cols[-9]


dna_variant_clones <- generate_dna_variant_heatmap(obj = obj, selected_variants_df = variants_annotated[1:4,])

# --------------------------------------------------------------- #
# Compute ploidy

library(tibble)

## Filter
## Amplicon completeness
amp_completeness <- 10
## Amplicon read depth
amp_readDepth <- 10
## Mean cell read depth
amp_meanCellRead <- 10
# clone group set a reference
diploid_ref <- "1"

filtered_data <- filter_cnv_profile(obj,
                                    dna_variant_clones$clones,
                                    amp_completeness = amp_completeness,
                                    amp_readDepth = amp_readDepth,
                                    amp_meanCellRead = amp_meanCellRead)

ploidy_data <- process_cnv_to_clonal_profile(filtered_data,
                                             dna_variant_clones$clones,
                                             diploid_ref = diploid_ref,
                                             exclude_clone = "small")

obj$ploidy.mtx <- ploidy_data

render_annotation_table(obj = obj, ploidy_data = obj$ploidy.mtx)$symbol |> unique() |> sort()


cnv_heatmap <- plot_cnv_heatmap(obj = obj, ploidy_data = obj$ploidy.mtx)

cnv_heatmap <- plot_cnv_heatmap(obj = obj, ploidy_data = obj$ploidy.mtx, display_gene = TRUE)



## Annotate amplicons
test <- annotate_genomic_regions(region_data = obj$cnv_id_table, build = "hg38")



dumb_df <- data.frame(chrom = "1", start_pos = 438149280, end_pos = 438149280)
annotate_genomic_regions(region_data = dumb_df, build = "hg19")



# --------------------------------------------------------------- #
# Lineplot

display_gene = TRUE

mat_data <- t(ploidy_data)

# --- 1. Reorder according to chromosomal position ---
# Note: We can now drop 'dplyr::' because of the @importFrom above
tmp_var_table <- obj$cnv_id_table |>
    filter(dna_id %in% colnames(mat_data)) |>
    arrange(as.numeric(chrom), as.numeric(start_pos)) |>
    mutate(chr_lit = paste0("chr", chrom))

mat_data <- mat_data[, tmp_var_table$dna_id]

tmp_split_table <- tmp_var_table[match(colnames(mat_data), tmp_var_table$dna_id), ]
sorted_gen_levels <- sort_genomic_chromosomes(tmp_split_table$chrom)
# ---------------------------- #
# If used want gene names add a column with gene annotation
# tmp_annotation <- annotate_genomic_regions(region_data = obj$cnv_id_table, build = "hg19")
# tmp_var_table <- merge(tmp_var_table, tmp_annotation[,"dna_id","hgnc_symbol"])

tmp_split_vec <- annotate_genomic_regions(region_data = tmp_split_table,
                                          build = obj$cnv_metadata$genome_version)
split_vec <- factor(tmp_split_vec$symbol,
                    levels = unique(tmp_split_vec$symbol))


# --- 2. Color Definition ---
col_fun <- colorRamp2(
    breaks = c(
        quantile(mat_data, c(0, 0.25)),
        2,
        quantile(mat_data, c(0.75, 1))
    ),
    colors = c("black", "#4575B4", "#F0F0F0", "#D73027", "#67001F")
)

group_colors <- setNames(
    viridis(nrow(mat_data)),
    nm = rownames(mat_data)
)

clone_selected <- "1"

mat_data_restricted <-as.matrix(mat_data[clone_selected,])|>t()
rownames(mat_data_restricted) <- rownames(mat_data)[rownames(mat_data) == clone_selected]
# --- 3. Left Annotation ---
left_ann <- rowAnnotation(
    df = data.frame(Group = rownames(mat_data_restricted)),
    col = list(Group = group_colors),
    show_legend = FALSE,
    simple_anno_size = unit(1, "cm"),
    show_annotation_name = FALSE
)

gene_annotation = data.frame('Gene' = tmp_split_vec$symbol,
                             'Chromosome' = tmp_split_vec$chrom,
                             'Probe' = tmp_split_vec$dna_id,
                             'Chrom_pos' = factor(tmp_split_vec$chr_lit,
                                                  levels = unique(sort_genomic_chromosomes(tmp_split_vec$chrom))),
                             'Chrom_start' = tmp_split_vec$start_pos)


plot_cnv_genome(cnv_matrix = mat_data, sub_indices = "2",
                gene_annotation = gene_annotation, lineplot_type = "genes+amplicons")


# --------------------------------------------------------------- #
# Protein analysis

