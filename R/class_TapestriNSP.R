# R6 and other dependencies are loaded via the package namespace

# =========================================================================
# CONSTANTES
# =========================================================================
COLORS <- c(
  '#1f77b4',
  '#ff7f0e',
  '#d62728',
  '#2ca02c',
  '#9467bd',
  '#8c564b',
  '#e377c2',
  '#17becf',
  '#bcbd22',
  '#bfbfbf',
  '#aec7e8',
  '#ffbb78',
  '#ff9896',
  '#98df8a',
  '#c5b0d5',
  '#c49c94',
  '#f7b6d2',
  '#9edae5',
  '#dbdb8d',
  '#9f9f9f',
  '#44708f',
  '#c3834a',
  '#aa5353',
  '#498349',
  '#937ca8',
  '#7c615b',
  '#c892b7',
  '#4598a1',
  '#969649',
  '#7f7f7f',
  '#bcc9da',
  '#ddbb9a',
  '#e5b1b0',
  '#a6ca9f',
  '#c4b9cc',
  '#b8a4a0',
  '#e7c6d4',
  '#b0ced3',
  '#c7c8a0',
  '#5f5f5f',
  '#576d7c',
  '#a58568',
  '#946969',
  '#587458',
  '#92879d',
  '#746663',
  '#baa0b2',
  '#5c868a',
  '#83835c',
  '#3f3f3f',
  '#c4cad2',
  '#ccbbab',
  '#d8bebd',
  '#adbfaa',
  '#c3bec7',
  '#b2a8a6',
  '#dfced5',
  '#b9c8ca',
  '#bebeaa',
  '#1f1f1f',
  '#000000'
)

# =========================================================================
# CLASS : ExpressionProfile
# =========================================================================
#' @export
ExpressionProfile <- R6Class(
  "ExpressionProfile",
  public = list(
    NUM_BINS = 1000,
    bandwidth = NULL,
    expression_ = NULL,
    liklihood_ = NULL,
    grid_ = NULL,
    peaks_ = NULL,
    valleys_ = NULL,

    initialize = function(bandwidth = 0.02) {
      self$bandwidth <- bandwidth
    },

    fit = function(expression) {
      self$expression_ <- expression
      expr_max <- max(expression)
      expr_norm <- if (expr_max == 0) expression else expression / expr_max

      # Strict equivalent of KDEUnivariate(cut=1)
      kde <- stats::density(
        expr_norm,
        bw = self$bandwidth,
        kernel = "gaussian",
        n = self$NUM_BINS,
        cut = 1
      )

      filt <- (kde$x >= min(expr_norm)) & (kde$x <= max(expr_norm))
      grid <- kde$x[filt]
      likelihood <- kde$y[filt]

      peaks <- data.frame(x = numeric(), p = numeric())
      valleys <- data.frame(x = numeric(), p = numeric())

      add_extrema <- function(index, type) {
        center <- grid[index]
        num_cells <- sum(
          expr_norm > (center - 2 * self$bandwidth) &
            expr_norm < (center + 2 * self$bandwidth)
        )
        if (num_cells > 5) {
          row <- data.frame(x = center * expr_max, p = likelihood[index])
          if (type == "peak") {
            peaks <<- rbind(peaks, row)
          } else {
            valleys <<- rbind(valleys, row)
          }
        }
      }

      half_window_size <- as.integer(self$bandwidth / (grid[2] - grid[1]) / 2)

      for (i in seq_along(grid)) {
        win_start <- max(1, i - half_window_size)
        win_end <- min(length(grid), i + half_window_size)
        window_lik <- likelihood[win_start:win_end]

        if (likelihood[i] == max(window_lik)) {
          add_extrema(i, "peak")
        }
        if (likelihood[i] == min(window_lik)) add_extrema(i, "valley")
      }

      self$liklihood_ <- likelihood
      self$grid_ <- grid
      self$peaks_ <- peaks
      self$valleys_ <- valleys
      invisible(self)
    },

    plot = function() {
      if (is.null(self$expression_)) {
        stop("Execute `fit` first.")
      }

      fig <- plot_ly() |>
        add_trace(
          x = ~ self$expression_,
          y = ~ rep(0, length(self$expression_)),
          type = "scatter",
          mode = "markers",
          name = "ticks",
          hoverinfo = "skip",
          marker = list(
            size = 5,
            opacity = 0.15,
            symbol = "line-ns",
            line = list(color = "black", width = 1)
          )
        ) |>
        add_trace(
          x = ~ self$grid_ * max(self$expression_),
          y = ~ self$liklihood_,
          type = "scatter",
          mode = "none",
          fill = "tozeroy",
          fillcolor = COLORS[1],
          opacity = 0.5,
          name = "histogram"
        )

      shapes <- list()
      if (nrow(self$peaks_) > 0) {
        for (x in self$peaks_$x) {
          shapes <- append(
            shapes,
            list(list(
              type = "line",
              x0 = x,
              x1 = x,
              y0 = 0,
              y1 = 1,
              yref = "paper",
              line = list(color = COLORS[3], width = 2, dash = "dash")
            ))
          )
        }
      }
      if (nrow(self$valleys_) > 0) {
        for (x in self$valleys_$x) {
          shapes <- append(
            shapes,
            list(list(
              type = "line",
              x0 = x,
              x1 = x,
              y0 = 0,
              y1 = 1,
              yref = "paper",
              line = list(color = COLORS[1], width = 2, dash = "dash")
            ))
          )
        }
      }

      fig |>
        layout(
          shapes = shapes,
          showlegend = FALSE,
          xaxis = list(title = "Normalized expression"),
          yaxis = list(title = "Distribution"),
          template = "plotly_white"
        )
    }
  )
)

# =========================================================================
# CLASS : NSP (Noise corrected and Scaled Protein counts)
# =========================================================================
#' @export
NSP <- R6Class(
  "NSP",
  public = list(
    jitter = NULL,
    random_state = NULL,
    sample_size = NULL,
    max_zero_read_cells = NULL,
    profiles_ = list(),
    subset_ = NULL,
    scale_ = NULL,
    scales_ = NULL,
    antibodies_ = NULL,
    reads_ = NULL,
    background_ = NULL,
    signal_ = NULL,
    f_sig_ = NULL,
    f_back_ = NULL,

    initialize = function(
      jitter = 0.5,
      random_state = NULL,
      sample_size = Inf,
      max_zero_read_cells = 0.05
    ) {
      self$jitter <- jitter
      self$random_state <- random_state
      self$sample_size <- sample_size
      self$max_zero_read_cells <- max_zero_read_cells
    },

    transform = function(reads, scale = NULL) {
      self$scale_ <- 1
      self$reads_ <- reads
      self$profiles_ <- list()

      no_reads <- rowSums(reads) == 0
      if (all(no_reads)) {
        return(matrix(0, nrow = nrow(reads), ncol = ncol(reads)))
      }
      if (any(no_reads)) {
        warning(sprintf(
          "NSP ignored %d cells with 0 protein reads.",
          sum(no_reads)
        ))
      }

      pos_reads <- reads[!no_reads, , drop = FALSE]
      if (nrow(pos_reads) == 1) {
        normal_counts <- (reads - min(reads)) / max(1, max(reads) - min(reads))
        return(private$expand_to_all_cells(reads, normal_counts))
      }

      if (is.null(scale)) {
        sf <- self$scaling_factor(pos_reads)
        scale <- sf$scale
        if (length(sf$ab) <= ncol(pos_reads) / 2) scale <- 1
      }

      normal_counts <- self$arcsinh(pos_reads / scale)
      total_reads <- log10(rowSums(pos_reads))

      num_cells <- min(self$sample_size, nrow(pos_reads))
      if (!is.null(self$random_state)) {
        # set.seed(self$random_state)
      }
      subset <- sample(seq_len(nrow(pos_reads)), num_cells, replace = FALSE)

      cat("Fitting GMMs...\n")
      gmm_res <- private$cell_signal_and_background(normal_counts[
        subset,
        ,
        drop = FALSE
      ])
      signal <- gmm_res$signal
      background <- gmm_res$background

      # Replacing np.polyfit with lm (simple linear model) with data.frame to avoid environment issues
      df_sig <- data.frame(signal = signal, tr = total_reads[subset])
      fit_sig <- lm(signal ~ tr, data = df_sig)

      df_back <- data.frame(background = background, tr = total_reads[subset])
      fit_back <- lm(background ~ tr, data = df_back)

      f_sig <- function(x) coef(fit_sig)[1] + coef(fit_sig)[2] * x
      f_back <- function(x) coef(fit_back)[1] + coef(fit_back)[2] * x

      cell_factor <- f_back(total_reads)
      true_signal <- f_sig(total_reads) - f_back(total_reads)
      ab_factor <- colSums(pos_reads > 0) / nrow(pos_reads)

      # numpy Broadcasting -> Outer product in R
      normal_counts <- normal_counts - outer(cell_factor, ab_factor, "*")
      normal_counts <- sweep(normal_counts, 1, pmax(true_signal, 1), `/`)

      self$scale_ <- scale
      self$subset_ <- subset
      self$reads_ <- pos_reads
      self$signal_ <- signal
      self$background_ <- background
      self$f_sig_ <- f_sig
      self$f_back_ <- f_back

      return(private$expand_to_all_cells(reads, normal_counts))
    },

    scaling_factor = function(reads, jitter = 0.5) {
      asinh_vals <- self$arcsinh(reads, jitter = jitter)
      bins <- seq(0, 19, by = 1)
      peaks <- private$get_ab_profile_peaks(asinh_vals)

      ps <- private$get_possible_scales(bins, peaks)
      self$scales_ <- ps$scales
      self$antibodies_ <- ps$antibodies

      hs <- private$get_highest_possible_scale(
        reads,
        self$scales_,
        self$antibodies_
      )
      list(scale = hs$scale, ab = hs$ab)
    },

    arcsinh = function(reads, jitter = NULL) {
      jitter <- if (is.null(jitter)) self$jitter else jitter
      if (!is.null(self$random_state)) {
        # set.seed(self$random_state)
      }
      noise <- matrix(
        rnorm(length(reads), 0, jitter),
        nrow = nrow(reads),
        ncol = ncol(reads)
      )
      asinh(reads + noise)
    }
  ),

  private = list(
    cell_signal_and_background = function(asinh_reads) {
      pbapply::pboptions(type = "txt")
      means_mat <- do.call(
        rbind,
        pbapply::pblapply(seq_len(nrow(asinh_reads)), function(i) {
          x <- asinh_reads[i, ]
          if (!is.null(self$random_state)) {
            # set.seed(self$random_state + i)
          }
          suppressWarnings({
            # We force K-means initialization (like scikit-learn)
            km <- kmeans(x, centers = 2)
            mc <- mclust::Mclust(
              x,
              G = 2,
              modelNames = "V",
              initialization = list(z = mclust::unmap(km$cluster)),
              verbose = FALSE
            )
            if (is.null(mc)) c(km$centers) else mc$parameters$mean
          })
        })
      )
      list(
        signal = apply(means_mat, 1, max),
        background = apply(means_mat, 1, min)
      )
    },

    get_ab_profile_peaks = function(asinh) {
      self$profiles_ <- list()
      lapply(seq_len(ncol(asinh)), function(i) {
        prof <- ExpressionProfile$new()
        prof$fit(asinh[, i])
        self$profiles_[[i]] <- prof
        if (nrow(prof$peaks_) > 0) {
          prof$peaks_$x[prof$peaks_$x > 1]
        } else {
          numeric(0)
        }
      })
    },

    get_possible_scales = function(bins, peaks) {
      if (all(vapply(peaks, length, integer(1)) == 0)) {
        return(list(scales = numeric(), antibodies = list()))
      }

      bin_ids <- private$get_bin_ids_with_most_peaks(bins, peaks)
      all_scales <- numeric()
      all_antibodies <- list()

      for (b in bin_ids) {
        low_bin <- bins[b]
        high_bin <- bins[b + 1] # Shifted R indexing (bins[1]=0, bins[2]=1)
        p_in_bin <- lapply(peaks, function(p) p[p >= low_bin & p < high_bin])
        n_peaks <- vapply(p_in_bin, length, integer(1))

        if (max(n_peaks) == 1) {
          sc <- 0
          ab <- integer()
          for (i in seq_along(p_in_bin)) {
            if (length(p_in_bin[[i]]) > 0) {
              sc <- sc + p_in_bin[[i]][1]
              ab <- c(ab, i)
            }
          }
          all_scales <- c(all_scales, sinh(sc / length(ab)))
          all_antibodies <- append(all_antibodies, list(ab))
        } else {
          sub_bins <- seq(low_bin, high_bin, length.out = 3)
          sub_res <- private$get_possible_scales(sub_bins, p_in_bin)
          all_scales <- c(all_scales, sub_res$scales)
          all_antibodies <- append(all_antibodies, sub_res$antibodies)
        }
      }
      list(scales = all_scales, antibodies = all_antibodies)
    },

    get_highest_possible_scale = function(reads, scales, abs) {
      if (length(scales) == 0) {
        return(list(scale = 1.0, ab = integer()))
      }
      valid <- integer()
      for (i in order(scales)) {
        if (
          mean(rowSums(round(reads / scales[i])) == 0) <=
            self$max_zero_read_cells
        ) {
          valid <- c(valid, i)
        }
      }
      if (length(valid) == 0) {
        return(list(scale = 1.0, ab = integer()))
      }

      n_abs <- vapply(abs[valid], length, integer(1))
      best_idx <- valid[n_abs == max(n_abs)]
      best_idx <- best_idx[length(best_idx)] # Takes the largest (last)
      list(scale = scales[best_idx], ab = abs[[best_idx]])
    },

    get_bin_ids_with_most_peaks = function(bins, peaks) {
      # findInterval is the strict and vectorized equivalent of np.digitize
      ids <- unlist(lapply(peaks, function(p) unique(findInterval(p, bins))))
      if (length(ids) == 0) {
        return(integer())
      }
      tbl <- table(ids)
      as.integer(names(tbl)[tbl == max(tbl)])
    },

    expand_to_all_cells = function(reads, norm) {
      no_reads <- rowSums(reads) == 0
      if (any(no_reads)) {
        res <- matrix(0, nrow(reads), ncol(reads))
        res[!no_reads, ] <- norm
        return(res)
      }
      norm
    }
  )
)

#' Normalize Protein Reads
#'
#' Wrapper function to normalize Tapestri protein read counts.
#'
#' @param reads A numeric matrix of raw reads (rows = cells, columns = antibodies).
#' @param method Character string: "CLR", "asinh", "NSP" or "ANSP".
#' @param jitter Standard deviation of the Gaussian noise to add (default = 0.5). Useful for "NSP" and "asinh".
#' @param scale Scaling factor (for "NSP" / "ANSP").
#' @param sample_size Number of cells to use to estimate model parameters ("ANSP").
#' @param random_state Integer to set the random seed (reproducibility).
#'
#' @return A list containing the `normalized_counts` matrix and the `nsp_model` object (if applicable).
#' @export
normalize_reads <- function(
  reads,
  method = "CLR",
  jitter = 0.5,
  scale = NULL,
  sample_size = 1000,
  random_state = NULL
) {
  if (method == "NSP") {
    nsp_obj <- NSP$new(jitter = jitter, random_state = random_state)
    normal_counts <- nsp_obj$transform(reads, scale = scale)
    return(list(normalized_counts = normal_counts, nsp_model = nsp_obj))
  } else if (method == "ANSP") {
    nsp_obj <- NSP$new(
      jitter = jitter,
      random_state = random_state,
      sample_size = sample_size
    )
    normal_counts <- nsp_obj$transform(reads, scale = scale)
    return(list(normalized_counts = normal_counts, nsp_model = nsp_obj))
  } else if (method == "CLR") {
    # log(x + 1)
    normal_counts <- log1p(reads)
    # Subtract the mean of each row (Translation of Numpy broadcasting: normal_counts.mean(axis=1)[:, None])
    normal_counts <- sweep(normal_counts, 1, rowMeans(normal_counts), "-")
    return(list(normalized_counts = normal_counts, nsp_model = NULL))
  } else if (method == "asinh") {
    if (!is.null(random_state)) {
      # set.seed(random_state)
    }
    # Gaussian noise creation
    noise <- matrix(
      rnorm(length(reads), mean = 0, sd = jitter),
      nrow = nrow(reads),
      ncol = ncol(reads)
    )
    normal_counts <- asinh(reads + noise)
    return(list(normalized_counts = normal_counts, nsp_model = NULL))
  } else {
    stop("Invalid value: Please provide one of {'CLR', 'NSP', 'ANSP', 'asinh'}")
  }
}
