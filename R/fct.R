
sanitize_filename <- function(filename, replacement = "_") {

    # 1. (Optional) Handle accents/special characters
    # Converts "Hélène" to "Helene".
    # Note: `iconv` results can vary slightly depending on OS locale.
    filename <- iconv(filename, to = "ASCII//TRANSLIT")

    # 2. Replace spaces with the replacement character
    filename <- gsub(" ", replacement, filename)

    # 3. Replace any character that is NOT alphanumeric, a dot, or a dash
    # Regex explanation:
    # [^...]      : Match any character NOT in this set
    # [:alnum:]   : Alphanumeric characters (letters and numbers)
    # \\.         : A literal dot (needs double escape)
    # -           : A literal dash
    filename <- gsub("[^[:alnum:]\\.-]", replacement, filename)

    # 4. Remove repeated instances of the replacement character
    # Example: "My___File" becomes "My_File"
    pattern_repeat <- paste0(replacement, "+")
    filename <- gsub(pattern_repeat, replacement, filename)

    # 5. Remove the replacement character from the start or end of the string
    # Example: "_file_" becomes "file"
    pattern_edge <- paste0("^", replacement, "|", replacement, "$")
    filename <- gsub(pattern_edge, "", filename)

    return(filename)
}

#' Protein matrix normalization
#'
#' Normalizes the protein matrix within a `ScIGMA_object` using the
#' centered log‑ratio (CLR) transform, inspired by the analogous function
#' in the optima package.
#'
#' @param ScIGMA_object A `ScIGMA_object` (R6).
#' @import compositions
#' @return The input `ScIGMA_object` with `protein.mtx` normalized and
#'  `protein.normalize.method` set to "normalized".
#' @keywords ScIGMA protein normalization
#' @export
#' @examples
#' \dontrun{
#' scigma <- normalizeProtein(scigma)
#' }

normalizeProtein <- function(ScIGMA_object) {
    # extract count matrix
    inputMatrix <- ScIGMA_object$protein.mtx.filtered |> as.matrix()
    # apply normalization CLR method
    ret <- (compositions::clr(inputMatrix + 1))

    ScIGMA_object$protein.mtx.filtered.normalized <- as.matrix(ret)
    ScIGMA_object$protein.normalize.method <- "CLR normalized"
    return(ScIGMA_object)
}



# --------------------------------------------------------------- #
#' Render ridge plot with ggplot2 and plotly
#'
#' @import ggplot2
#' @import ggridges
#' @import plotly
#' @import tidyr
render_protein_ridge_plot <- function(obj){
    tmp_data <- obj$protein.mtx.filtered.normalized |>
        as_tibble() |>
        pivot_longer(everything())
    tmp_plot <- tmp_data |>
        ggplot(aes(x=value, y=name, fill = name)) +
        geom_density_ridges() +
        theme_ridges() +
        theme(legend.position = "none")
    print("render plotlied plot")
    # print("tmp_plot")
    # print(tmp_plot)
    print(tmp_plot) |> plotly::ggplotly()
}

#' Génère un barplot de la proportion de protéines
#'
#' @import MatrixGenerics
#'
#' @param obj Un objet contenant la matrice protéique (obj$protein.mtx)
#' @param title Titre optionnel du graphique
#'
#' @return Un objet ggplot représentant la distribution relative des protéines
#' @export
#'
#' @examples
#' plot_protein_barplot(obj)
plot_protein_barplot <- function(obj, title = "Protein Percentage Distribution") {

    print("obj$protein.mtx")
    print(dim(obj$protein.mtx))
    print(head(obj$protein.mtx, 2))
    print(class(obj$protein.mtx))

    # Vérifications basiques
    if (is.null(obj$protein.mtx)) {
        stop("Object do not have a protein matrix (obj$protein.mtx).")
    }

    # Calcul des pourcentages
    protein_barplot_df <- data.frame(
        protein = colnames(obj$protein.mtx),
        percent = round(colSums2(obj$protein.mtx) / sum(obj$protein.mtx), 5) * 100
    )

    # Génération du plot
    protein_barplot <- protein_barplot_df |>
        ggplot(aes(x = protein, y = percent)) +
        geom_bar(stat = "identity", fill = "steelblue") +
        xlab("Antibody") +
        ylab("Percentage") +
        ggtitle(title) +
        theme_minimal() +
        theme(
            axis.text.x = element_text(angle = 45, vjust = 1, hjust = 1),
            plot.title = element_text(hjust = 0.5)
        )

    print(protein_barplot) |> plotly::ggplotly()
}

#----------------------------#
# Utilitaires internes
#----------------------------#

#' Estimation d'un facteur d'échelle global pour runs "oversequencés"
#' Heuristique : la doc NSP évoque 'scaling_factor(reads, jitter)' et une contrainte
#'   max_zero_read_cells : fraction max de cellules avec 0 reads après scaling.
#' Ici, on propose une stratégie robuste & simple, cohérente avec la doc.
estimate_scaling_factor <- function(total_reads,
                                    max_zero_read_cells = 0.05) {
    # total_reads : vecteur (longueur = nb cellules)
    # max_zero_read_cells : seuil (0–1) accepté de cellules à 0 après scaling
    #
    # 1) Si déjà peu de zéros, on n'impose pas de scaling.
    zero_frac <- mean(total_reads == 0)
    if (!is.finite(zero_frac) || zero_frac > max_zero_read_cells) {
        return(1.0)
    }
    # 2) Vise à réduire une longue traîne en utilisant un ratio médiane/p95
    #    (si p95 >> médiane, on rabaisse l'échelle)
    if (all(total_reads == 0)) return(1.0)
    tr_pos <- total_reads[total_reads > 0]
    if (length(tr_pos) < 10) return(1.0)

    med <- stats::median(tr_pos)
    p95 <- stats::quantile(total_reads, 0.95, names = FALSE)
    if (!is.finite(med) || !is.finite(p95) || p95 <= 0) return(1.0)

    # 3) Proposition : scale = med / p95 (<= 1 si p95 > med)
    proposal <- as.numeric(med / p95)
    proposal <- max(min(proposal, 1.0), 0.0)

    # 4) Vérifie que le scaling ne crée pas trop de zéros
    trial <- floor(total_reads * proposal)
    if (mean(trial == 0) <= max_zero_read_cells) {
        return(proposal)
    } else {
        return(1.0)
    }
}

#' Sélectionne indices "background" (bas) et "signal" (haut) pour un vecteur y
select_bg_sig_indices <- function(y, p_low = 0.10, p_high = 0.90,
                                  min_points = 20L) {
    # y : intensités d'une protéine sur un sous-ensemble de cellules
    # p_low/p_high : quantiles pour définir fond/signal
    # min_points : garde-fou si quantiles donnent trop peu de points
    q_low  <- stats::quantile(y, p_low,  names = FALSE, type = 7)
    q_high <- stats::quantile(y, p_high, names = FALSE, type = 7)

    bg_idx <- which(y <= q_low)
    sg_idx <- which(y >= q_high)

    # garde-fou : si trop peu de points, on prend extrêmes
    if (length(bg_idx) < min_points) {
        bg_idx <- order(y, decreasing = FALSE)[seq_len(min(min_points, length(y)))]
    }
    if (length(sg_idx) < min_points) {
        sg_idx <- order(y, decreasing = TRUE)[seq_len(min(min_points, length(y)))]
    }
    list(bg_idx = bg_idx, sg_idx = sg_idx, q_low = q_low, q_high = q_high)
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
#' @export
normalize_linear_regression <- function(raw_matrix,
                                        use_log = TRUE,
                                        add_mean = TRUE,
                                        jitter = 0,
                                        scale = 1) {

    # 1. Validation & Setup
    if (!is.matrix(raw_matrix)) raw_matrix <- as.matrix(raw_matrix)

    # 2. Apply Jitter (Noise Injection)
    if (jitter > 0) {
        n_elements <- length(raw_matrix)
        noise <- rnorm(n = n_elements, mean = 0, sd = jitter)
        raw_matrix <- raw_matrix + matrix(noise, nrow = nrow(raw_matrix), ncol = ncol(raw_matrix))
        raw_matrix[raw_matrix < 0] <- 0 # Clip negatives
    }

    # 3. Calculate Library Size (Technical Covariate)
    library_size <- rowSums(raw_matrix)

    # 4. Transform for Regression
    if (use_log) {
        work_mat <- log1p(raw_matrix)
        work_lib_size <- log1p(library_size)
    } else {
        work_mat <- raw_matrix
        work_lib_size <- library_size
    }

    # 5. Regression Loop (Regress Out Depth)
    norm_mat <- apply(work_mat, 2, function(protein_counts) {
        # Check for zero variance to avoid lm failure
        if (var(protein_counts) == 0) return(protein_counts)

        fit <- lm(protein_counts ~ work_lib_size)
        resid <- residuals(fit)

        if (add_mean) {
            return(resid + mean(protein_counts, na.rm = TRUE))
        } else {
            return(resid)
        }
    })

    # Restore dimensions/names lost during apply
    rownames(norm_mat) <- rownames(raw_matrix)
    colnames(norm_mat) <- colnames(raw_matrix)

    # 6. Apply Scaling (Dynamic Range Compression)
    # "The amount by which the read counts are scaled down."
    final_scale <- scale

    if (is.null(final_scale)) {
        # Estimation Strategy:
        # If no scale is provided, we use the global Standard Deviation of the
        # corrected matrix. This standardizes the data unit (similar to Z-score scaling
        # but maintaining the relative mean differences if add_mean=TRUE).
        message("Estimating scale factor from data distribution...")
        final_scale <- sd(norm_mat, na.rm = TRUE)
        message(sprintf("Estimated Scale Factor: %.4f", final_scale))
    }

    # Prevent division by zero
    if (final_scale == 0) final_scale <- 1

    # Apply scaling
    norm_mat <- norm_mat / final_scale

    return(norm_mat)
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
#' @export
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
    set.seed(seed)
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
