# ==============================================================================
# MANUSCRIPT FIGURES GENERATION FUNCTIONS
# ==============================================================================

#' Helper to setup figure path
#' @noRd
setup_fig_path <- function(filename) {
    dir.create("results/figures", recursive = TRUE, showWarnings = FALSE)
    return(file.path("results/figures", filename))
}

# ==============================================================================
# 1. SNV - DNA VARIANT
# ==============================================================================

#' SNV - DNA variant heatmap (post-COMPASS genotype)
plot_fig_snv_heatmap <- function(ScIGMA_obj) {
    message("Generating SNV Heatmap...")
    fig_path <- setup_fig_path("Fig1_SNV_Genotype_Heatmap.pdf")

    ht_res <- ScIGMA:::generate_dna_variant_heatmap(
        obj = ScIGMA_obj,
        selected_variants_df = ScIGMA_obj$variants.filtered,
        heatmap_include_all_samples = FALSE,
        use_imputed = TRUE
    )
    pdf(fig_path)
    ComplexHeatmap::draw(ht_res$heatmap)

    dev.off()
    return(fig_path)
}

#' SNV - COMPASS Maximum likelihood phylogeny
plot_fig_snv_phylogeny <- function(ScIGMA_obj) {
    message("Generating SNV Phylogeny Tree...")
    fig_path <- setup_fig_path("Fig2_SNV_Phylogeny_Tree.pdf")
    tree_content <- S4Vectors::metadata(
        ScIGMA_obj$mae
    )$compass$tree_dot
    pdf(fig_path)
    DiagrammeR::grViz(tree_content)

    dev.off()
    return(fig_path)
}

# ==============================================================================
# 2. CNV
# ==============================================================================

#' CNV - Heatmap projection
plot_fig_cnv_heatmap <- function(ScIGMA_obj) {
    message("Generating CNV Heatmap...")

    # Get cnv x_axis
    cnv_id_table <- as.data.frame(SummarizedExperiment::rowData(ScIGMA_obj$mae[[
        "amplicons"
    ]]))
    features <- sort_genomic_chromosomes(
        cnv_id_table$chrom
    )

    fig_path <- setup_fig_path("Fig3_CNV_Heatmap.pdf")
    pdf(fig_path)
    ScIGMA:::generate_cnv_heatmap_filtered(
        obj = ScIGMA_obj,
        features = features,
        projection_type = "Position"
    )

    dev.off()
    return(fig_path)
}

#' CNV - Lineplot projection
plot_fig_cnv_lineplot <- function(ScIGMA_obj) {
    message("Generating CNV Lineplot...")
    fig_path <- setup_fig_path("Fig4_CNV_Lineplot.pdf")
    p <- ScIGMA:::generate_cnv_lineplot_filtered(
        obj = ScIGMA_obj,
        projection_type = "Position",
        clone = params$cnv_clonal_inference$ref_clone
    )
    pdf(fig_path)
    print(p)

    dev.off()
    return(fig_path)
}

# ==============================================================================
# 3. PROTEIN
# ==============================================================================

#' Protein - Description: Ridgeplot
plot_fig_protein_ridgeplot <- function(ScIGMA_obj) {
    message("Generating Protein Ridgeplot...")
    fig_path <- setup_fig_path("Fig5_Protein_Ridgeplot.pdf")
    p_ridge <- ScIGMA:::generate_protein_ridgeplot(ScIGMA_obj)
    pdf(fig_path)
    print(p_ridge)

    dev.off()
    return(fig_path)
}

#' Protein - Description: Barplot
plot_fig_protein_barplot <- function(ScIGMA_obj) {
    message("Generating Protein Barplot...")
    fig_path <- setup_fig_path("Fig6_Protein_Abundance_Barplot.pdf")
    p_bar <- ScIGMA:::generate_protein_barplot(
        ScIGMA_obj,
        title = "Protein Percentage Distribution"
    )
    pdf(fig_path)
    print(p_bar)

    dev.off()
    return(fig_path)
}

#' Protein - Immunophenotype gating: Biplot CD19 vs CD33 + DNA clone projection
plot_fig_protein_biplot_gating <- function(ScIGMA_obj) {
    message("Generating Protein Biplot Gating (CD19 vs CD33)...")
    fig_path <- setup_fig_path("Fig7_Protein_Biplot_DNA_Clones.pdf")
    p_biplot <- ScIGMA:::generate_protein_biplot(
        obj = ScIGMA_obj,
        xvar = "CD19",
        yvar = "CD33",
        logx = TRUE, # Applique log1p pour gérer les outliers de distribution
        logy = TRUE,
        color_genotype = c("1"), # Met en évidence le clone 1 (les autres sont en gris transparent au fond)
        title = "Biplot Immunophénotypique CD19 vs CD33"
    )
    pdf(fig_path)
    print(p_biplot)

    dev.off()
    return(fig_path)
}

#' Protein - UMAP: Global
plot_fig_protein_umap <- function(ScIGMA_obj) {
    message("Generating Protein Global UMAP...")
    fig_path <- setup_fig_path("Fig8_Protein_UMAP_Global.pdf")
    p <- ScIGMA:::generate_protein_umap(ScIGMA_obj)
    pdf(fig_path)
    print(p)
    dev.off()
    return(fig_path)
}

#' Protein - UMAP: Markers expression (Facetted)
plot_fig_protein_umap_markers <- function(ScIGMA_obj) {
    message("Generating Protein UMAP Markers (Facetted)...")
    fig_path <- setup_fig_path("Fig9_Protein_UMAP_Markers_Facetted.pdf")
    # We use a default set of markers (e.g. CD19, CD33) if available
    markers_to_plot <- intersect(
        c("CD19", "CD33", "CD34", "CD38"),
        rownames(ScIGMA_obj$mae[["proteins"]])
    )
    p <- ScIGMA:::generate_protein_umap_markers(ScIGMA_obj, markers_to_plot)
    pdf(fig_path)
    print(p)
    dev.off()
    return(fig_path)
}

#' Protein - UMAP: Immunophenotype gating projection
plot_fig_protein_umap_gating <- function(ScIGMA_obj) {
    message("Generating Protein UMAP Gating Projection...")
    fig_path <- setup_fig_path("Fig10_Protein_UMAP_Gating.pdf")
    p <- generate_protein_umap_gating(ScIGMA_obj)
    pdf(fig_path)
    print(p)
    dev.off()
    return(fig_path)
}

#' Protein - UMAP: Unsupervised clustering projection
plot_fig_protein_umap_clustering <- function(ScIGMA_obj) {
    message("Generating Protein UMAP Clustering Projection...")
    fig_path <- setup_fig_path("Fig11_Protein_UMAP_Clustering.pdf")
    p <- ScIGMA:::generate_protein_umap_clustering(ScIGMA_obj)
    pdf(fig_path)
    print(p)
    dev.off()
    return(fig_path)
}

# ==============================================================================
# 4. MULTIOMICS INTEGRATION
# ==============================================================================

#' Multiomics - UMAP with projection of DNA clones
plot_fig_multi_umap_dna_clones <- function(ScIGMA_obj) {
    message("Generating Multiomics UMAP (DNA Clones)...")
    fig_path <- setup_fig_path("Fig12_Multiomics_UMAP_DNA_Clones.pdf")
    p <- ScIGMA:::generate_multi_umap_dna_clones(ScIGMA_obj)
    pdf(fig_path)
    print(p)
    dev.off()
    return(fig_path)
}

#' Multiomics - UMAP with genotype of a specific DNA clone
plot_fig_multi_umap_dna_genotype <- function(ScIGMA_obj) {
    message("Generating Multiomics UMAP (Clone Genotype)...")
    fig_path <- setup_fig_path("Fig13_Multiomics_UMAP_Clone_Genotype.pdf")
    variant_to_plot <- rownames(ScIGMA_obj$mae[["dna_variants"]])[1]
    p <- ScIGMA:::generate_multi_umap_dna_genotype(ScIGMA_obj, variant_to_plot)
    pdf(fig_path)
    print(p)
    dev.off()
    return(fig_path)
}

#' Multiomics - Barplot of DNA clones across UMAP unsupervised clusters
plot_fig_multi_barplot_clones_per_cluster <- function(ScIGMA_obj) {
    message(
        "Generating Multiomics Barplot: DNA Clones vs Unsupervised Clusters..."
    )
    fig_path <- setup_fig_path(
        "Fig14_Multiomics_Barplot_Clones_vs_Clusters.pdf"
    )
    p <- ScIGMA:::generate_multi_barplot_clones_per_cluster(ScIGMA_obj)
    pdf(fig_path)
    print(p)
    dev.off()
    return(fig_path)
}

#' Multiomics - Barplot of DNA variants according to clusters and genes
plot_fig_multi_barplot_variant_genotype <- function(ScIGMA_obj) {
    message(
        "Generating Multiomics Barplot: DNA Variants vs Clusters & Genes..."
    )
    fig_path <- setup_fig_path(
        "Fig15_Multiomics_Barplot_Variants_vs_Clusters.pdf"
    )
    variants_to_plot <- rownames(ScIGMA_obj$mae[["dna_variants"]])[
        1:min(3, nrow(ScIGMA_obj$mae[["dna_variants"]]))
    ]
    p <- ScIGMA:::generate_multi_barplot_variant_genotype(
        ScIGMA_obj,
        variants_to_plot
    )
    pdf(fig_path)
    print(p)
    dev.off()
    return(fig_path)
}

#' Multiomics - Immunophenotype Gating: Barplot of gating subpopulation genotype vs DNA variants
plot_fig_multi_barplot_gating_genotype <- function(ScIGMA_obj) {
    message(
        "Generating Multiomics Barplot: Gating Subpop Genotype vs Variants..."
    )
    fig_path <- setup_fig_path(
        "Fig16_Multiomics_Barplot_Gating_vs_Variants.pdf"
    )
    variants_to_plot <- rownames(ScIGMA_obj$mae[["dna_variants"]])[
        1:min(3, nrow(ScIGMA_obj$mae[["dna_variants"]]))
    ]
    p <- ScIGMA:::generate_multi_barplot_gating_genotype(
        ScIGMA_obj,
        variants_to_plot
    )
    pdf(fig_path)
    print(p)
    dev.off()
    return(fig_path)
}
