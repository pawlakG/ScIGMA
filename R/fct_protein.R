#' Protein matrix normalization
#'
#' Normalizes the protein matrix within a `ScIGMA_object` using the
#' centered log-ratio (CLR) transform, inspired by the analogous function
#' in the optima package.
#'
#' @param obj A `ScIGMA_object` (R6).
#' @import compositions
#' @return The input `ScIGMA_object` with `protein.mtx` normalized and
#'  `protein.normalize.method` set to "normalized".
#' @keywords ScIGMA protein normalization
#' @noRd
#' @examples
#' \dontrun{
#' scigma <- normalize_protein_clr(scigma)
#' }
normalize_protein_clr <- function(obj) {
    inputMatrix <- SummarizedExperiment::assay(obj$mae[["proteins"]], "counts") |> as.matrix()

    tmp_clr <- compositions::clr(t(inputMatrix) + 1)

    tmp_mat <- as.matrix(unclass(tmp_clr))
    tmp_mat <- apply(tmp_mat, 2, \(x) x + abs(min(x))) # Translate to non-negative values

    ret <- t(tmp_mat)

    SummarizedExperiment::assay(obj$mae[["proteins"]], "normalized") <- ret

    S4Vectors::metadata(obj$mae)$protein_normalize_method <- "CLR normalized"

    obj$protein.filtered <- TRUE

    return(obj)
}

#' Normalize Protein Data using NSP
#'
#' @description
#' Normalizes protein expression data using the NSP method.
#'
#' @param obj A `ScIGMA_object` (R6).
#'
#' @return The input `ScIGMA_object` with `protein.mtx` normalized and
#'  `protein.normalize.method` set to "NSP normalized".
#' @export
normalize_protein_nsp <- function(obj) {
    counts <- SummarizedExperiment::assay(obj$mae[["proteins"]], "counts")
    counts_ram <- as.matrix(counts)
    
    # Transpose so that cells are in rows for NSP$transform
    counts_ram_t <- t(counts_ram)
    
    nsp <- NSP$new()
    norm_prot_t <- nsp$transform(counts_ram_t)
    
    # We retranspose to have proteins in rows
    norm_prot <- t(norm_prot_t)
    
    # Name restoration
    rownames(norm_prot) <- rownames(counts_ram)
    colnames(norm_prot) <- colnames(counts_ram)
    
    SummarizedExperiment::assay(obj$mae[["proteins"]], "normalized") <- norm_prot
    S4Vectors::metadata(obj$mae)$protein_normalize_method <- "NSP normalized"
    obj$protein.filtered <- TRUE
    
    return(obj)
}

#' Normalize Protein Data using asinh
#'
#' @description
#' Normalizes protein expression data using asinh transformation.
#'
#' @param obj A `ScIGMA_object` (R6).
#' @param cofactor Numeric. The cofactor for asinh transformation.
#'
#' @return The input `ScIGMA_object` with `protein.mtx` normalized and
#'  `protein.normalize.method` set to "asinh normalized".
#' @export
normalize_protein_asinh <- function(obj, cofactor = 5) {
    counts <- SummarizedExperiment::assay(obj$mae[["proteins"]], "counts")
    counts_ram <- as.matrix(counts)
    
    norm_prot <- asinh(counts_ram / cofactor)
    
    SummarizedExperiment::assay(obj$mae[["proteins"]], "normalized") <- norm_prot
    S4Vectors::metadata(obj$mae)$protein_normalize_method <- "asinh normalized"
    obj$protein.filtered <- TRUE
    
    return(obj)
}

#' Render ridge plot with ggplot2 and plotly
#' @noRd
render_protein_ridge_plot <- function(obj) {
    # We extract the CLR assay if it exists, otherwise fallback to counts
    # assay_to_use <- ifelse("clr" %in% SummarizedExperiment::assayNames(obj$mae[["proteins"]]), "clr", "counts")
    assay_to_use <- ifelse("normalized" %in% SummarizedExperiment::assayNames(obj$mae[["proteins"]]),
        "normalized", "counts"
    )

    tmp_data <- t(SummarizedExperiment::assay(obj$mae[["proteins"]], assay_to_use)) |>
        dplyr::as_tibble() |>
        tidyr::pivot_longer(dplyr::everything())

    tmp_plot <- tmp_data |>
        ggplot2::ggplot(ggplot2::aes(x = value, y = name, fill = name)) +
        ggridges::geom_density_ridges() +
        ggridges::theme_ridges() +
        ggplot2::theme(legend.position = "none")

    plotly::ggplotly(tmp_plot)
}

#' @noRd
plot_protein_barplot <- function(obj, title = "Protein Percentage Distribution") {
    ptn_mtx <- SummarizedExperiment::assay(obj$mae[["proteins"]], "counts")
    if (is.null(ptn_mtx)) {
        stop("Object does not have a protein count matrix.")
    }

    protein_barplot_df <- data.frame(
        protein = rownames(ptn_mtx),
        percent = round(DelayedMatrixStats::rowSums2(ptn_mtx) / sum(ptn_mtx), 5) * 100
    )

    protein_barplot <- protein_barplot_df |>
        ggplot2::ggplot(ggplot2::aes(x = protein, y = percent)) +
        ggplot2::geom_bar(stat = "identity", fill = "steelblue") +
        ggplot2::xlab("Antibody") +
        ggplot2::ylab("Percentage") +
        ggplot2::ggtitle(title) +
        ggplot2::theme_minimal() +
        ggplot2::theme(
            axis.text.x = ggplot2::element_text(angle = 45, vjust = 1, hjust = 1),
            plot.title = ggplot2::element_text(hjust = 0.5)
        )

    plotly::ggplotly(protein_barplot)
}

#' Normalize Protein Data using DSB (Denoised Scaled by Background)
#'
#' @description
#' Normalizes protein expression data using the DSB method, which corrects for background
#' noise using empty droplets.
#'
#' @param obj A `ScIGMA_object` (R6).
#' @param min_bg_counts Minimum number of background counts to include an empty droplet.
#' @param max_bg_counts Maximum number of background counts to include an empty droplet.
#'
#' @return The input `ScIGMA_object` with `protein.mtx` normalized and
#'  `protein.normalize.method` set to "DSB normalized".
#' @export
normalize_protein_dsb <- function(obj, min_bg_counts = 10, max_bg_counts = 500) {
    cells_mtx <- SummarizedExperiment::assay(obj$mae[["proteins"]], "counts")
    empty_mtx_delayed <- obj$get_empty_drops()
    
    if (is.null(empty_mtx_delayed)) {
        stop("Empty drops matrix is not available for DSB normalization.")
    }
    
    drop_sums <- DelayedMatrixStats::colSums2(empty_mtx_delayed)
    valid_empty_barcodes <- colnames(empty_mtx_delayed)[drop_sums > min_bg_counts & drop_sums < max_bg_counts]
    
    if (length(valid_empty_barcodes) == 0) {
        stop("No valid empty droplets found after filtering with the given background counts thresholds.")
    }
    
    empty_mtx_ram <- as.matrix(empty_mtx_delayed[, valid_empty_barcodes, drop = FALSE])
    cells_mtx_ram <- as.matrix(cells_mtx)
    
    # Absolute security: we force alignment of row names.
    # Both matrices share the same layout (same panel), but the main MAE object
    # may have undergone typographical cleaning (e.g. sanitize_protein_markers) that the empty drops did not.
    rownames(empty_mtx_ram) <- rownames(cells_mtx_ram)
    
    norm_prot <- dsb::DSBNormalizeProtein(
        cell_protein_matrix = cells_mtx_ram,
        empty_drop_matrix = empty_mtx_ram,
        denoise.counts = TRUE,
        use.isotype.control = FALSE
    )
    
    SummarizedExperiment::assay(obj$mae[["proteins"]], "normalized") <- norm_prot
    S4Vectors::metadata(obj$mae)$protein_normalize_method <- "DSB normalized"
    obj$protein.filtered <- TRUE
    
    return(obj)
}

#' Normalize Protein Data using Linear Regression (with Jitter & Scaling)
#'
#' @description
#' Corrects protein expression by regressing out the effect of library size.
#' It includes options for jittering (to handle discrete zero-inflation) and
#' scaling (to compress the dynamic range or standardize variance).
#'
#' @param raw_matrix Numeric matrix (Cells x Proteins).
#' @param use_log Logical. If TRUE (recommended), performs regression on log(x+1).
#' @param add_mean Logical. If TRUE, adds the original mean back to residuals
#' to preserve the biological baseline.
#' @param jitter Numeric. The standard deviation of Gaussian noise added to raw counts.
#' @param scale Numeric or NULL. The factor by which the corrected counts are divided.
#' If NULL, the algorithm estimates it using the global standard deviation of the
#' corrected data (robust estimation of distribution spread).
#'
#' @return A numeric matrix of normalized (and optionally scaled) values.
#' @noRd
normalize_linear_regression <- function(raw_matrix,
                                        use_log = TRUE,
                                        add_mean = TRUE,
                                        jitter = 0,
                                        scale = 1,
                                        seed = 42) {
    if (!is.matrix(raw_matrix)) raw_matrix <- as.matrix(raw_matrix)

    if (jitter > 0) {
        if (!is.null(seed)) {
            if (exists(".Random.seed", envir = .GlobalEnv, inherits = FALSE)) {
                old_seed <- get(".Random.seed", envir = .GlobalEnv, inherits = FALSE)
                on.exit(assign(".Random.seed", old_seed, envir = .GlobalEnv), add = TRUE)
            } else {
                on.exit(rm(list = ".Random.seed", envir = .GlobalEnv), add = TRUE)
            }
            set.seed(seed)
        }
        noise <- matrix(rnorm(length(raw_matrix), 0, jitter), nrow = nrow(raw_matrix))
        raw_matrix <- raw_matrix + noise
        raw_matrix[raw_matrix < 0] <- 0 # Clip negatives
    }

    # 3. Calculate Library Size (CRITIQUE: colSums car les cellules sont en colonnes)
    library_size <- colSums(raw_matrix)

    # 4. Transform for Regression
    if (use_log) {
        work_mat <- log1p(raw_matrix)
        work_lib_size <- log1p(library_size)
    } else {
        work_mat <- raw_matrix
        work_lib_size <- library_size
    }

    norm_mat <- apply(work_mat, 1, function(protein_counts) {
        if (var(protein_counts) == 0) {
            return(protein_counts)
        }

        fit <- lm(protein_counts ~ work_lib_size)
        resid <- residuals(fit)

        if (add_mean) {
            return(resid + mean(protein_counts, na.rm = TRUE))
        } else {
            return(resid)
        }
    })

    norm_mat <- t(norm_mat)
    rownames(norm_mat) <- rownames(raw_matrix)
    colnames(norm_mat) <- colnames(raw_matrix)

    # 6. Apply Scaling
    final_scale <- scale
    if (is.null(final_scale)) {
        message("Estimating scale factor from data distribution...")
        final_scale <- sd(norm_mat, na.rm = TRUE)
    }

    if (final_scale == 0) final_scale <- 1

    return(norm_mat / final_scale)
}


#' Run UMAP on Protein Expression Data
#'
#' @description
#' Performs Uniform Manifold Approximation and Projection (UMAP) on an expression matrix.
#' Optimized for single-cell protein data using the 'uwot' implementation (C++).
#' It is strongly recommended to use normalized data (e.g., CLR or NSP/Arcsinh) as input.
#'
#' @param expression_matrix Numeric matrix or data.frame (Cells x Features).
#' Rows should represent cells and columns should represent protein markers.
#' @param n_neighbors Integer. The size of the local neighborhood (default: 30).
#' Larger values (30-50) preserve more global structure (trajectories), while smaller
#' values (5-15) focus on local clusters.
#' @param min_dist Numeric. The effective minimum distance between embedded points (default: 0.3).
#' Controls how tightly points are packed.
#' @param n_components Integer. The dimension of the space to embed into (default: 2).
#' @param metric Character. The distance metric to use (default: "cosine").
#' "cosine" is often superior to "euclidean" for high-dimensional cytometry data.
#' @param seed Integer. Random seed for reproducibility.
#' @param n_threads Integer or NULL. Number of threads for parallel processing.
#' If NULL, uses all available cores minus one.
#'
#' @return A matrix of UMAP coordinates (Cells x n_components) with row names
#' matching the input matrix. Column names are "UMAP_1", "UMAP_2", etc.
#'
#' @noRd
#'
#' @examples
#' # 1. Normalize data first (Crucial!)
#' norm_mat <- normalize_protein_nsp(raw_counts, scale = 5)
#'
#' # 2. Run UMAP
#' umap_coords <- run_umap_protein(norm_mat, n_neighbors = 30, min_dist = 0.3)
#'
#' # 3. Merge with data for plotting
#' plot_data <- cbind(as.data.frame(umap_coords), as.data.frame(norm_mat))
#'
run_umap_protein <- function(expression_matrix,
                            n_neighbors = 30,
                            min_dist = 0.3,
                            n_components = 2,
                            metric = "cosine",
                            seed = 42,
                            n_threads = NULL) {
    # 1. Dependency Check
    if (!requireNamespace("uwot", quietly = TRUE)) {
        stop("Error: Package 'uwot' is required. Please install it with install.packages('uwot').")
    }

    # 2. Input Validation
    if (!is.matrix(expression_matrix) && !is.data.frame(expression_matrix)) {
        stop("Error: 'expression_matrix' must be a matrix or data.frame.")
    }

    mat <- as.matrix(expression_matrix)

    if (!is.numeric(mat)) {
        stop("Error: Input matrix must contain numeric values.")
    }

    # Scientific sanity check: Raw counts?
    # If max value > 100, it's likely raw counts. UMAP on raw counts is essentially
    # a "Library Size Map". We issue a warning.
    if (max(mat, na.rm = TRUE) > 100) {
        warning("Warning: Input values seem high (>100). Are you running UMAP on RAW counts?
            It is strongly recommended to use normalized data (CLR or NSP/Arcsinh)
            to avoid clustering based purely on sequencing depth.")
    }

    message(sprintf("Running UMAP on %d cells and %d features...", nrow(mat), ncol(mat)))

    # 3. Setup Parallelism
    if (is.null(n_threads)) {
        # Conservative default: All cores - 1
        n_threads <- max(1, parallel::detectCores() - 1)
    }

    # 4. Execution (uwot)

    umap_result <- uwot::umap(
        X = mat,
        n_neighbors = n_neighbors,
        min_dist = min_dist,
        n_components = n_components,
        metric = metric,
        n_threads = n_threads,
        ret_model = FALSE # We only want coordinates here
    )

    # 5. Formatting Output
    colnames(umap_result) <- paste0("UMAP_", seq_len(n_components))
    rownames(umap_result) <- rownames(mat)

    return(umap_result)
}


#' Comprehensive UMAP Evaluation Suite (True Ground Truth)
#' Computes Exact Trustworthiness, Continuity (Venna & Kaski) and Spearman topology
#' via stochastic sampling to ensure O(M^2) speed instead of O(N^2).
#' Fully vectorized to avoid R for-loops.
#' @param high_dim_mat Original normalized matrix (Cells x Proteins)
#' @param low_dim_mat UMAP coordinates matrix (Cells x 2)
#' @param k Number of neighbors for local structure evaluation
#' @param sample_size Number of cells to sample for exact calculation
#' @return A list of metric scores
#' @noRd
compute_umap_metrics <- function(high_dim_mat, low_dim_mat, k = 15, sample_size = 1000) {
    common_cells <- intersect(rownames(high_dim_mat), rownames(low_dim_mat))

    if (length(common_cells) < 3) {
        warning("Not enough common cells to compute UMAP metrics.")
        return(list(trustworthiness = 0, continuity = 0, global_spearman = 0))
    }

    high_dim_mat <- high_dim_mat[common_cells, , drop = FALSE]
    low_dim_mat <- low_dim_mat[common_cells, , drop = FALSE]
    N <- nrow(high_dim_mat)

    sample_n <- min(sample_size, N)
    k_safe <- min(k, sample_n - 2)

    idx <- sample(1:N, sample_n)
    high_sub <- high_dim_mat[idx, , drop = FALSE]
    low_sub <- low_dim_mat[idx, , drop = FALSE]

    d_high <- as.matrix(dist(high_sub))
    d_low <- as.matrix(dist(low_sub))

    if (sd(d_high, na.rm = TRUE) == 0 || sd(d_low, na.rm = TRUE) == 0) {
        spearman_score <- 0
    } else {
        spearman_score <- cor(
            d_high[lower.tri(d_high)],
            d_low[lower.tri(d_low)],
            method = "spearman"
        )
    }

    diag(d_high) <- Inf
    diag(d_low) <- Inf

    r_high <- t(apply(d_high, 1, rank, ties.method = "first"))
    r_low <- t(apply(d_low, 1, rank, ties.method = "first"))

    denom <- ifelse(k_safe < sample_n / 2,
        sample_n * k_safe * (2 * sample_n - 3 * k_safe - 1) / 2,
        sample_n * (sample_n - k_safe) * (sample_n - k_safe - 1) / 2
    )
    norm_factor <- 1 / denom

    # trust_penalty : On prend les rangs d'origine des K-voisins de l'UMAP
    trust_penalty <- sum(pmax(0, r_high[r_low <= k_safe] - k_safe))

    # cont_penalty : On prend les rangs UMAP des K-voisins d'origine
    cont_penalty <- sum(pmax(0, r_low[r_high <= k_safe] - k_safe))

    trustworthiness <- max(0, 1 - (trust_penalty * norm_factor))
    continuity <- max(0, 1 - (cont_penalty * norm_factor))

    return(list(
        trustworthiness = trustworthiness,
        continuity = continuity,
        global_spearman = spearman_score
    ))
}

#' Run Unsupervised Clustering on Protein Markers
#' @description Headless wrapper for Seurat-based unsupervised clustering on protein data.
#' @param obj A ScIGMA_object (R6) containing a valid seurat_object with PCA computed.
#' @param resolution Numeric. Resolution parameter for FindClusters (modularity optimization). Default is 0.8.
#' @param seed Integer. Seed for reproducibility. Default is 42.
#' @return The mutated ScIGMA_object with assigned seurat_clusters.
#' @export
run_protein_clustering <- function(obj, resolution = 0.8, seed = 42) {
    if (is.null(obj$seurat_object)) {
        stop("seurat_object is missing. Run PCA/UMAP steps first.")
    }
    
    message("1/2 - Computing nearest neighbor graph...")
    obj$seurat_object <- Seurat::FindNeighbors(
        object = obj$seurat_object,
        features = rownames(obj$seurat_object),
        dims = seq_len(ncol(obj$seurat_object@reductions$pca@cell.embeddings))
    )
    
    message("2/2 - Optimizing modularity (Louvain)...")
    obj$seurat_object <- Seurat::FindClusters(
        object = obj$seurat_object,
        resolution = resolution,
        random.seed = seed
    )
    
    return(obj)
}

#' Generate Protein Abundance Ridge Plot
#'
#' @description
#' API function to generate the protein ridge plot outside the Shiny context.
#'
#' @param obj A ScIGMA_data object.
#' @return An interactive Plotly object.
#' @export
generate_protein_ridgeplot <- function(obj) {
    render_protein_ridge_plot(obj)
}

#' Generate Protein Barplot
#'
#' @description
#' API function to generate the protein distribution barplot outside the Shiny context.
#'
#' @param obj A ScIGMA_data object.
#' @param title The plot title.
#' @return An interactive Plotly object.
#' @export
generate_protein_barplot <- function(obj, title = "Protein Percentage Distribution") {
    plot_protein_barplot(obj, title = title)
}

#' Generate Protein Biplot (e.g. CD19 vs CD33)
#'
#' @description
#' API function to generate an interactive biplot (scatter plot) between 
#' two protein markers outside the Shiny context, with the possibility to color 
#' by a specific genotype/clone.
#'
#' @param obj A ScIGMA_data object.
#' @param xvar The marker name for the X axis (e.g., "CD19").
#' @param yvar The marker name for the Y axis (e.g., "CD33").
#' @param logx Boolean. Apply log1p on the X axis. Default FALSE.
#' @param logy Boolean. Apply log1p on the Y axis. Default FALSE.
#' @param color_genotype Character vector indicating the clones to color (e.g., "1", "2"). If NULL or "None", uniform coloring.
#' @param title The plot title.
#' @return An interactive Plotly object.
#' @export
generate_protein_biplot <- function(
    obj, 
    xvar, 
    yvar, 
    logx = FALSE, 
    logy = FALSE, 
    color_genotype = NULL,
    title = "Protein Biplot"
) {
    if (!"proteins" %in% names(obj$mae)) {
        stop("The 'proteins' matrix is missing from the MAE object.")
    }

    # Prism style copied from the Shiny module
    prism_axis_style <- list(
        titlefont = list(size = 16, color = "black", family = "Arial"),
        tickfont = list(size = 14, color = "black", family = "Arial"),
        showline = TRUE,
        linewidth = 2,
        linecolor = "black",
        mirror = FALSE,
        ticks = "outside",
        tickwidth = 2,
        ticklen = 6,
        tickcolor = "black",
        showgrid = FALSE,
        zeroline = FALSE
    )

    assay_to_use <- ifelse(
        "normalized" %in% SummarizedExperiment::assayNames(obj$mae[["proteins"]]),
        "normalized",
        "counts"
    )

    ptn_matrix <- SummarizedExperiment::assay(obj$mae[["proteins"]], assay_to_use)

    if (!xvar %in% rownames(ptn_matrix)) stop(paste("The marker", xvar, "is not present."))
    if (!yvar %in% rownames(ptn_matrix)) stop(paste("The marker", yvar, "is not present."))

    plot_df <- data.frame(
        x = ptn_matrix[xvar, ],
        y = ptn_matrix[yvar, ],
        custom_id = colnames(ptn_matrix)
    )

    if (logx) plot_df$x <- log1p(plot_df$x)
    if (logy) plot_df$y <- log1p(plot_df$y)

    color_formula <- NULL
    color_palette <- NULL
    marker_config <- list(size = 5, opacity = 0.8, color = "#2c3e50")

    if (!is.null(color_genotype) && !identical(color_genotype, "None")) {
        if (is.null(obj$dna.clones)) stop("Missing clonal assignment (obj$dna.clones).")

        # Extraction of clones associated with cells
        plot_df$Genotype <- as.character(obj$dna.clones[rownames(plot_df)])

        # Grouping of non-targets
        plot_df$Genotype <- ifelse(
            plot_df$Genotype %in% color_genotype,
            plot_df$Genotype,
            "Other"
        )

        # Locking Z-Indexing ("Other" at the bottom)
        plot_df$Genotype <- factor(
            plot_df$Genotype,
            levels = c("Other", color_genotype)
        )
        plot_df <- plot_df[order(plot_df$Genotype), ]

        color_formula <- ~Genotype

        # Merging palettes
        if (!is.null(obj$dna_clone_colors)) {
            color_palette <- c(obj$dna_clone_colors, "Other" = "#e0e0e033")
        } else {
            warning("obj$dna_clone_colors not found, using the default palette.")
            color_palette <- "Set1" 
        }

        marker_config <- list(size = 5)
    }

    p <- plotly::plot_ly(
        data = plot_df,
        x = ~x,
        y = ~y,
        key = ~custom_id,
        color = color_formula,
        colors = color_palette,
        type = "scattergl",
        mode = "markers",
        marker = marker_config
    ) %>%
        plotly::layout(
            title = list(
                text = paste("<b>", title, "</b>"),
                font = list(family = "Arial", size = 18)
            ),
            plot_bgcolor = "white",
            paper_bgcolor = "white",
            xaxis = c(
                list(
                    title = paste("<b>", xvar, "</b>"),
                    range = list(0, max(plot_df$x, na.rm = TRUE) + 1)
                ),
                prism_axis_style
            ),
            yaxis = c(
                list(
                    title = paste("<b>", yvar, "</b>"),
                    range = list(0, max(plot_df$y, na.rm = TRUE) + 1)
                ),
                prism_axis_style
            ),
            legend = list(title = list(text = "<b>Genotype</b>")),
            margin = list(l = 60, r = 30, b = 60, t = 50)
        ) %>%
        plotly::config(displaylogo = FALSE)

    return(p)
}
#' Generate Global Protein UMAP
#' @export
generate_protein_umap <- function(obj) {
    if (is.null(obj$seurat_object) || is.null(obj$seurat_object@reductions$umap)) {
        stop("UMAP not computed.")
    }
    umap_df <- as.data.frame(obj$seurat_object@reductions$umap@cell.embeddings)
    p <- ggplot2::ggplot(umap_df, ggplot2::aes(x = umap_1, y = umap_2)) +
        ggplot2::geom_point(size = 1.5, color = "#2c3e50", alpha = 0.8) +
        ggplot2::xlab("UMAP 1") +
        ggplot2::ylab("UMAP 2") +
        ggprism::theme_prism(base_size = 14, base_fontface = "bold")
    return(p)
}

#' Generate Protein UMAP with Markers
#' @export
generate_protein_umap_markers <- function(obj, markers) {
    if (is.null(obj$seurat_object) || is.null(obj$seurat_object@reductions$umap)) stop("UMAP not computed.")
    umap_df <- as.data.frame(obj$seurat_object@reductions$umap@cell.embeddings)
    umap_df$barcode <- rownames(umap_df)
    
    assay_to_use <- ifelse("normalized" %in% SummarizedExperiment::assayNames(obj$mae[["proteins"]]), "normalized", "counts")
    expr_mat <- SummarizedExperiment::assay(obj$mae[["proteins"]], assay_to_use)[markers, , drop = FALSE]
    expr_df <- as.data.frame(t(expr_mat))
    expr_df$barcode <- rownames(expr_df)
    
    umap_df <- merge(umap_df, expr_df, by = "barcode")
    umap_df <- tidyr::pivot_longer(umap_df, cols = tidyr::all_of(markers), values_to = "ptn_expression", names_to = "marker")
    
    p <- ggplot2::ggplot(umap_df, ggplot2::aes(x = umap_1, y = umap_2, color = ptn_expression)) +
        ggplot2::geom_point(size = 1.2, alpha = 0.9) +
        ggplot2::facet_wrap(~marker) +
        ggplot2::scale_color_viridis_c(option = "inferno", name = "Expression") +
        ggplot2::xlab("UMAP 1") +
        ggplot2::ylab("UMAP 2") +
        ggprism::theme_prism(base_size = 14, base_fontface = "bold")
    return(p)
}

#' Generate Protein UMAP with Gating
#' @export
generate_protein_umap_gating <- function(obj) {
    if (is.null(obj$seurat_object) || is.null(obj$seurat_object@reductions$umap)) stop("UMAP not computed.")
    umap_df <- as.data.frame(obj$seurat_object@reductions$umap@cell.embeddings)
    umap_df$barcode <- rownames(umap_df)
    
    if (is.null(obj$protein_gates)) stop("No gates found.")
    
    # Assign gate to each cell (last gate wins if overlapping)
    umap_df$Gate <- "Ungated"
    for (gate_name in names(obj$protein_gates)) {
        cells <- obj$protein_gates[[gate_name]]
        umap_df$Gate[umap_df$barcode %in% cells] <- gate_name
    }
    
    p <- ggplot2::ggplot(umap_df, ggplot2::aes(x = umap_1, y = umap_2, color = Gate)) +
        ggplot2::geom_point(size = 1.5, alpha = 0.8) +
        ggplot2::xlab("UMAP 1") +
        ggplot2::ylab("UMAP 2") +
        ggprism::theme_prism(base_size = 14, base_fontface = "bold") +
        ggplot2::theme(legend.title = ggplot2::element_text(face = "bold"))
    return(p)
}

#' Generate Protein UMAP with Unsupervised Clustering
#' @export
generate_protein_umap_clustering <- function(obj) {
    if (is.null(obj$seurat_object) || is.null(obj$seurat_object@reductions$umap)) stop("UMAP not computed.")
    if (!"seurat_clusters" %in% colnames(obj$seurat_object@meta.data)) stop("Clustering not computed.")
    
    umap_df <- as.data.frame(obj$seurat_object@reductions$umap@cell.embeddings)
    umap_df$Cluster <- obj$seurat_object@meta.data$seurat_clusters
    
    p <- ggplot2::ggplot(umap_df, ggplot2::aes(x = umap_1, y = umap_2, color = Cluster)) +
        ggplot2::geom_point(size = 1.5, alpha = 0.8) +
        ggplot2::xlab("UMAP 1") +
        ggplot2::ylab("UMAP 2") +
        ggprism::theme_prism(base_size = 14, base_fontface = "bold")
    return(p)
}
