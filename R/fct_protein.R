#' Protein matrix normalization
#'
#' Normalizes the protein matrix within a `ScIGMA_object` using the
#' centered log-ratio (CLR) transform, inspired by the analogous function
#' in the optima package.
#'
#' @param ScIGMA_object A `ScIGMA_object` (R6).
#' @import compositions
#' @return The input `ScIGMA_object` with `protein.mtx` normalized and
#'  `protein.normalize.method` set to "normalized".
#' @keywords ScIGMA protein normalization
#' @noRd
#' @examples
#' \dontrun{
#' scigma <- normalizeProtein(scigma)
#' }
normalizeProtein <- function(ScIGMA_object) {
    inputMatrix <- SummarizedExperiment::assay(ScIGMA_object$mae[["proteins"]], "counts") |> as.matrix()

    tmp_clr <- compositions::clr(t(inputMatrix) + 1)

    tmp_mat <- as.matrix(unclass(tmp_clr))
    tmp_mat <- apply(tmp_mat, 2, \(x) x + abs(min(x))) # Translate to non-negative values

    ret <- t(tmp_mat)

    SummarizedExperiment::assay(ScIGMA_object$mae[["proteins"]], "clr") <- ret

    S4Vectors::metadata(ScIGMA_object$mae)$protein_normalize_method <- "CLR normalized"

    ScIGMA_object$protein.filtered <- TRUE

    return(ScIGMA_object)
}

#' Render ridge plot with ggplot2 and plotly
#' @noRd
render_protein_ridge_plot <- function(obj) {
    # On extrait l'assay CLR s'il existe, sinon on fallback sur counts
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
