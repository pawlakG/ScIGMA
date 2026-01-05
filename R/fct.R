
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
