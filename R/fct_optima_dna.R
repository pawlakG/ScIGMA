# Fonctions for optima object handling

#' loadH5
#'
#' @description Load h5 file with a modified version of optima readHdf5 fonction.
#' This is a modifies version of optima's readHdf5 to handle cell label.
#'
#' @param filePath str h5 file path
#'
#' @return The return value, if any, from executing the function.
#'
#' @noRd
#' @import optima
loadH5 <- function(directory, sample.name, omic.type = "DNA+protein") {
    # --------------------------------------------------------------- #
    # Load file
    h5f <- rhdf5::H5Fopen(directory, flags = "H5F_ACC_RDONLY")
    if (omic.type == "DNA+protein") {
        my.proteins = as.character(h5f$assays$protein_read_counts$ca$id)
        my.protein.normalize.method = "unnormalized"
        my.protein.mtx = t(h5f$assays$protein_read_counts$layers$read_counts)
        rownames(my.protein.mtx) <- h5f$assays$protein_read_counts$ra$sample_name
        colnames(my.protein.mtx) <- h5f$assays$protein_read_counts$ca$id

    }
    else if (omic.type == "DNA") {
        cat("DNA data only, skip reading protein data...\n")
        my.proteins = "non-protein"
        my.protein.normalize.method = "non-protein"
        my.protein.mtx = matrix(NA)
    }
    else {
        stop("illegale argument for omic.type")
    }

    # --------------------------------------------------------------- #
    # Modification on original function : Set cell labels
    if (!is.null(h5f$assays$dna_read_counts$ra$sample_name)){
        cell.labels <- as.character(h5f$assays$dna_read_counts$ra$sample_name)
    } else if (!is.null(h5f$assays$dna_variants$ra$sample_name)){
        cell.labels <- as.character(h5f$assays$dna_variants$ra$sample_name)
    } else if (!is.null(h5f$assays$protein_read_counts$ra$sample_name)){
        cell.labels <- as.character(h5f$assays$protein_read_counts$ra$sample_name)
    } else {
        cell.labels <-  rep("unassigned", length(h5f$assays$dna_read_counts$ra$barcode))
    }

    # --------------------------------------------------------------- #
    # Prepare DNA variant
    vaf.mtx.tmp <- t(h5f$assays$dna_variants$layers$AF)
    dimnames(vaf.mtx.tmp) <- list(h5f$assays$dna_variants$ra$sample_name, h5f$assays$dna_variants$ca$id)
    gt.mtx.tmp <- t(h5f$assays$dna_variants$layers$NGT)
    dimnames(gt.mtx.tmp) <- list(h5f$assays$dna_variants$ra$sample_name, h5f$assays$dna_variants$ca$id)
    dp.mtx.tmp <- t(h5f$assays$dna_variants$layers$DP)
    dimnames(dp.mtx.tmp) <- list(h5f$assays$dna_variants$ra$sample_name, h5f$assays$dna_variants$ca$id)
    gq.mtx.tmp <- t(h5f$assays$dna_variants$layers$GQ)
    dimnames(gq.mtx.tmp) <- list(h5f$assays$dna_variants$ra$sample_name, h5f$assays$dna_variants$ca$id)

    # --------------------------------------------------------------- #
    # Prepare DNA read counts
    amp.mtx.tmp <- t(h5f$assays$dna_read_counts$layers$read_counts)
    dimnames(gq.mtx.tmp) <- list(h5f$assays$dna_read_counts$ra$sample_name, h5f$assays$dna_variants$ca$id)

    # --------------------------------------------------------------- #
    # Filter dna_variants according
    # --------------------------------------------------------------- #
    # Setup optima object
    optima.obj <- new("optima",
                      meta.data = sample.name,
                      cell.ids = as.character(h5f$assays$dna_read_counts$ra$barcode),
                      cell.labels = cell.labels,
                      variants = as.character(h5f$assays$dna_variants$ca$id),
                      variant.filter = "unfiltered",
                      # vaf.mtx = t(h5f$assays$dna_variants$layers$AF),
                      vaf.mtx = vaf.mtx.tmp,
                      # gt.mtx = t(h5f$assays$dna_variants$layers$NGT),
                      gt.mtx = gt.mtx.tmp,
                      # dp.mtx = t(h5f$assays$dna_variants$layers$DP),
                      dp.mtx = dp.mtx.tmp,
                      # gq.mtx = t(h5f$assays$dna_variants$layers$GQ),
                      gq.mtx = gq.mtx.tmp,
                      amps = as.character(h5f$assays$dna_read_counts$ca$id),
                      amp.normalize.method = "unnormalized",
                      # amp.mtx = t(h5f$assays$dna_read_counts$layers$read_counts),
                      amp.mtx = gq.mtx.tmp,
                      ploidy.mtx = matrix(),
                      proteins = my.proteins,
                      protein.normalize.method = my.protein.normalize.method,
                      protein.mtx = my.protein.mtx)
    # --------------------------------------------------------------- #
    # CNV preprocess
    ## Normalization
    optima.obj_CNVnorm <-  normalizeCNV(optima.obj)

    if (omic.type == "DNA+protein"){
        # --------------------------------------------------------------- #
        # Protein preprocess
        opt.obj.preprocessed <- normalizeProtein(optima.obj_CNVnorm)
    } else {
        opt.obj.preprocessed <- optima.obj_CNVnorm
    }
    # --------------------------------------------------------------- #
    # Set result
    y <- opt.obj.preprocessed
    return(y)
}






#' --------------------------------------------------------------- #
#' Function to get DNA clones
