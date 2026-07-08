#' Fetch Clinical VEP Annotations
#'
#' @param custom_variant_vector Character vector of variants (e.g. "chr1-115256669-G-A")
#' @param genome_build Genome build, default "grch37"
#' @return A tibble with annotations
#' @noRd
fetch_clinical_vep_annotations <- function(
    custom_variant_vector,
    genome_build = "grch37"
) {
    # Dictionary translation (3-letter -> 1-letter clinical)
    aa_map <- c(
        "Ala" = "A",
        "Arg" = "R",
        "Asn" = "N",
        "Asp" = "D",
        "Cys" = "C",
        "Gln" = "Q",
        "Glu" = "E",
        "Gly" = "G",
        "His" = "H",
        "Ile" = "I",
        "Leu" = "L",
        "Lys" = "K",
        "Met" = "M",
        "Phe" = "F",
        "Pro" = "P",
        "Ser" = "S",
        "Thr" = "T",
        "Trp" = "W",
        "Tyr" = "Y",
        "Val" = "V",
        "Ter" = "*"
    )

    # 1. Parsing HGVS Input (Robust parsing for MAE native ids like 'chr1:115256669:G/A')
    hgvs_queries <- custom_variant_vector %>%
        as.character() %>%
        chartr(":/", "--", .) %>%
        stringr::str_remove("^chr") %>%
        stringr::str_split("-") %>%
        purrr::map_chr(~ paste0(.x[1], ":g.", .x[2], .x[3], ">", .x[4]))

    # 2. Configuration Endpoint
    server <- if (genome_build == "grch37") {
        "https://grch37.rest.ensembl.org"
    } else {
        "https://rest.ensembl.org"
    }
    endpoint <- "/vep/human/hgvs?refseq=1&hgvs=1&variant_class=1&canonical=1&phenotypes=1"

    # 3. Batch Request
    body_json <- jsonlite::toJSON(
        list(hgvs_notations = hgvs_queries),
        auto_unbox = TRUE
    )

    response <- httr::POST(
        url = paste0(server, endpoint),
        httr::content_type("application/json"),
        httr::accept("application/json"),
        body = body_json
    )

    if (httr::status_code(response) != 200) {
        stop(sprintf(
            " [FATAL] API VEP %s: %s",
            httr::status_code(response),
            httr::content(response, "text")
        ))
    }

    # 4. Flattening and initial filtering
    raw_data <- httr::content(response, "parsed", simplifyVector = TRUE)

    annotation_table <- raw_data %>%
        tibble::as_tibble() %>%
        tidyr::unnest(
            transcript_consequences,
            keep_empty = TRUE,
            names_repair = "unique"
        ) %>%
        dplyr::filter(canonical == 1)


    expected_cols <- c("hgvsc", "hgvsp")
    for (col in expected_cols) {
        if (!col %in% names(annotation_table)) {
            annotation_table[[col]] <- NA_character_
        }
    }

    if (!"colocated_variants" %in% names(annotation_table)) {
        annotation_table$colocated_variants <- list(NULL)
    }

    # 5. Data engineering and 1-letter mapping
    annotation_table <- annotation_table %>%
        dplyr::mutate(
            consequence_terms = purrr::map_chr(
                consequence_terms,
                ~ if (length(.x) > 0) .x[1] else NA_character_
            ),
            cDNA = dplyr::if_else(
                !is.na(hgvsc),
                stringr::str_extract(hgvsc, "c\\..+"),
                NA_character_
            ),
            extracted_p = stringr::str_extract(hgvsp, "p\\..+"),
            extracted_p = stringr::str_replace_all(extracted_p, "%3D", "="),
            extracted_p = stringr::str_replace_all(extracted_p, aa_map),
            PROTEIN = dplyr::case_when(
                is.na(extracted_p) ~ NA_character_,
                is.na(gene_symbol) ~ extracted_p,
                TRUE ~ paste0(gene_symbol, ":", extracted_p)
            ),
            original_variant = custom_variant_vector[match(
                input,
                hgvs_queries
            )],
            CLINVAR = purrr::map2_chr(
                colocated_variants,
                original_variant,
                ~ {
                    if (is.null(.x) || !is.data.frame(.x)) {
                        return(NA_character_)
                    }

                    alt_allele <- stringr::str_extract(.y, "[^/]+$")

                    if ("clin_sig_allele" %in% names(.x)) {
                        pattern <- paste0("\\b", alt_allele, ":([^;]+)")
                        extracted <- stringr::str_match(
                            unlist(.x$clin_sig_allele),
                            pattern
                        )[, 2]
                        extracted <- extracted[!is.na(extracted)]
                        if (length(extracted) > 0) {
                            return(paste(unique(extracted), collapse = ","))
                        }
                    }

                    if (
                        "allele_string" %in%
                            names(.x) &&
                            "clin_sig" %in% names(.x)
                    ) {
                        valid_rows <- stringr::str_detect(
                            .x$allele_string,
                            paste0("\\b", alt_allele, "\\b")
                        )
                        if (any(valid_rows, na.rm = TRUE)) {
                            sigs <- unlist(.x$clin_sig[valid_rows])
                            if (length(sigs) > 0) {
                                return(paste(
                                    unique(sigs[!is.na(sigs)]),
                                    collapse = ","
                                ))
                            }
                        }
                    }
                    return(NA_character_)
                }
            )
        ) %>%
        dplyr::filter(!is.na(cDNA)) %>%
        dplyr::mutate(
            consequence_terms = stringr::str_replace_all(
                consequence_terms,
                "_",
                " "
            ),
            CLINVAR = stringr::str_replace_all(CLINVAR, "_", " ")
        ) %>%
        dplyr::select(
            original_variant,
            gene = gene_symbol,
            transcript_id = transcript_id,
            protein = PROTEIN,
            cdna = cDNA,
            variant_type = dplyr::any_of("variant_class"),
            gene_function = consequence_terms,
            clinvar = CLINVAR
        )



    return(annotation_table)
}

#' Filter variants and fetch annotations
#'
#' @param obj ScIGMA_object instance.
#' @param paths List of API paths for fetch_variants_batch_fields.
#' @param min_dp Integer. Minimum depth.
#' @param min_gq Integer. Minimum genotype quality.
#' @param vaf_ref Numeric. Maximum VAF for reference call.
#' @param vaf_hom Numeric. Minimum VAF for homozygous alternate call.
#' @param vaf_het Numeric. Maximum VAF for heterozygous call.
#' @param min_cell_pt Numeric. Minimum percentage of cells covering a variant.
#' @param min_mut_cell_pt Numeric. Minimum percentage of mutated cells.
#' @param batch_size Integer. Batch size for API calls.
#' @return Filtered and annotated ScIGMA_object.
#' @noRd
filter_and_annotate_variants <- function(
    obj,
    paths,
    min_dp = 10,
    min_gq = 30,
    vaf_ref = 5,
    vaf_hom = 95,
    vaf_het = 30,
    min_cell_pt = 10,
    min_mut_cell_pt = 10,
    batch_size = 300
) {
    # 1. Pipeline execution: Filtering
    obj <- filter_variant_ScIGMA_mae(
        obj = obj,
        min.dp = min_dp,
        min.gq = min_gq,
        vaf.ref = vaf_ref,
        vaf.hom = vaf_hom,
        vaf.het = vaf_het,
        min.cell.pt = min_cell_pt,
        min.mut.cell.pt = min_mut_cell_pt
    )

    message("Fetching and injecting annotations...")

    # 2. Strict ID extraction
    variant_ids <- rownames(obj$mae[["dna_variants"]])

    # Gatekeeper: Do not call API if filtering removed everything
    if (length(variant_ids) == 0) {
        warning("No variants left after filtering. Skipping annotation.")
        return(obj)
    }

    # 3. Protected API call
    annot_df <- tryCatch(
        {
            fetch_clinical_vep_annotations(variant_ids)
        },
        error = function(e) {
            stop(sprintf("API variant annotation failed: %s", e$message))
        }
    )

    if (!is.null(annot_df) && nrow(annot_df) > 0) {
        # Format variant_id like GENE:chrX:POS:REF/ALT for legacy compatibility downstream
        annot_df$variant_id <- ifelse(
            !is.na(annot_df$gene),
            paste0(
                annot_df$gene,
                ":",
                sub(
                    "^([^-]+)-([^-]+)-([^-]+)-([^-]+)$",
                    "\\1:\\2:\\3/\\4",
                    chartr(":/", "--", annot_df$original_variant)
                )
            ),
            paste0("Unmapped:", annot_df$original_variant)
        )

        # 4. Dimensional synchronization
        current_rowdata <- as.data.frame(
            SummarizedExperiment::rowData(obj$mae[["dna_variants"]])
        )
        current_rowdata$query_id <- rownames(current_rowdata)

        merged_rowdata <- merge(
            x = current_rowdata,
            y = annot_df,
            by.x = "query_id",
            by.y = "original_variant",
            all.x = TRUE,
            sort = FALSE
        )

        # 5. MAE Topology Guardrail (CRITICAL)
        # Prevents 1-to-many API mappings from expanding the rowData and crashing the MAE
        if (nrow(merged_rowdata) > nrow(current_rowdata)) {
            warning(
                "API returned multiple annotations per variant. Deduplicating to preserve MAE strict dimensions."
            )
            merged_rowdata <- merged_rowdata[
                !duplicated(merged_rowdata$query_id),
            ]
        }

        # 6. Canonical order restoration
        rownames(merged_rowdata) <- merged_rowdata$query_id
        merged_rowdata <- merged_rowdata[variant_ids, ]
        merged_rowdata$query_id <- NULL

        merged_rowdata$variant_id[is.na(merged_rowdata$variant_id)] <- paste0(
            "Unmapped:",
            rownames(merged_rowdata)
        )[is.na(merged_rowdata$variant_id)]

        SummarizedExperiment::rowData(obj$mae[[
            "dna_variants"
        ]]) <- S4Vectors::DataFrame(merged_rowdata)
        message("Annotation matrix successfully integrated into MAE rowData.")
    } else {
        warning("API returned an empty object. rowData remains unannotated.")
    }

    message("Calculating variant cell proportions...")

    vaf_mtx <- SummarizedExperiment::assay(obj$mae[["dna_variants"]], "vaf")

    mutated_cells_count <- DelayedMatrixStats::rowSums2(vaf_mtx > 0)

    total_cells <- ncol(vaf_mtx)

    SummarizedExperiment::rowData(obj$mae[[
        "dna_variants"
    ]])$cell_proportion <- round(mutated_cells_count / total_cells, 2)

    message("Cell proportions successfully added to rowData.")
    S4Vectors::metadata(obj$mae)$variant_filter <- "filtered"

    invisible(obj)
}

#' Run COMPASS MCMC Inference on a ScIGMA Object
#' @description This headless wrapper allows running the COMPASS phylogeny algorithm outside the Shiny UI.
#' @param obj A ScIGMA_object (R6) containing the filtered MAE.
#' @param chains Integer. Number of Markov Chains to run.
#' @param chain_length Integer. Iterations per chain.
#' @param patient_sex Character. Patient biological sex ("male" or "female").
#' @return The mutated ScIGMA_object with inferred phylogeny and clonal assignments.
#' @export
run_compass_inference <- function(
    obj,
    chains = 4L,
    chain_length = 500L,
    patient_sex = "female"
) {
    message("1/3 - Extracting HDF5 arrays and preparing matrix bindings...")

    # 1. Target variants identification
    if (is.null(obj$variants.filtered)) {
        target_vars <- rownames(obj$mae[["dna_variants"]])
    } else {
        target_vars <- rownames(obj$variants.filtered)
    }

    # 2. Build Core Matrices
    compass_inputs <- build_compass_matrices(
        obj = obj,
        selected_variants = target_vars
    )

    mat_ref <- t(as.matrix(compass_inputs$M_ref))
    mat_alt <- t(as.matrix(compass_inputs$M_alt))
    mat_cna <- t(as.matrix(compass_inputs$C))

    gt_assay <- SummarizedExperiment::assay(obj$mae[["dna_variants"]], "gt")
    mat_gt <- as.matrix(gt_assay[target_vars, , drop = FALSE])
    mat_gt[mat_gt == 3L] <- NA

    storage.mode(mat_ref) <- "integer"
    storage.mode(mat_alt) <- "integer"
    storage.mode(mat_cna) <- "integer"
    storage.mode(mat_gt) <- "integer"

    variant_matrices <- list(REF = mat_ref, ALT = mat_alt, GT = mat_gt)

    # 3. Locus and Chromosome Mapping
    dna_se <- obj$mae[["dna_variants"]]
    snv_sub <- as.data.frame(SummarizedExperiment::rowData(dna_se))[
        target_vars,
    ]
    vec_locus_names <- snv_sub$gene
    vec_locus_chrom <- snv_sub$chrom

    amp_se <- obj$mae[["amplicons"]]
    cna_row_data <- as.data.frame(SummarizedExperiment::rowData(amp_se))
    vec_region_names <- unique(paste0(
        cna_row_data$chrom,
        "_",
        vapply(cna_row_data$dna_id, \(x) strsplit(x, "_")[[1]][3], character(1))
    ))
    vec_region_chrom <- sub(
        "^chr",
        "",
        vapply(
            vec_region_names,
            \(x) strsplit(x, "_")[[1]][1],
            character(1),
            USE.NAMES = FALSE
        ),
        ignore.case = TRUE
    )

    use_cna <- if (ncol(variant_matrices$REF) != ncol(mat_cna)) FALSE else TRUE
    vec_locus_regions <- as.integer(compass_inputs$locus_regions)

    prefix_out <- file.path(
        tempdir(),
        paste0("compass_headless_", as.integer(Sys.time()))
    )

    message("2/3 - Executing MCMC Inference (C++ Engine)...")
    res <- run_compass_mcmc(
        variant_matrices = variant_matrices,
        locus_regions = vec_locus_regions,
        region_matrix = mat_cna,
        output_prefix = prefix_out,
        locus_names = vec_locus_names,
        locus_chromosomes = vec_locus_chrom,
        region_names = vec_region_names,
        region_chromosomes = vec_region_chrom,
        chains = as.integer(chains),
        chain_length = as.integer(chain_length),
        patient_sex = patient_sex,
        use_cna = use_cna
    )

    message("3/3 - Updating object with imputed genotypes and clonal labels...")

    # mat_imputed <- get_imputed_genotypes(prefix_out = res$prefix)
    # message("res")
    # message(res)
    # message("prefix_out")
    # message(prefix_out)
    mat_imputed <- get_imputed_genotypes(prefix_out = prefix_out)
    cells_to_keep <- rownames(mat_imputed)

    if (length(cells_to_keep) == 0) {
        stop("COMPASS filtered all cells (pure doublets).")
    }

    # 3a. Subsetting MAE to surviving cells
    obj$mae <- obj$mae[, cells_to_keep, ]

    if (!is.null(obj$seurat_object)) {
        valid_seurat_cells <- intersect(
            cells_to_keep,
            colnames(obj$seurat_object)
        )
        obj$seurat_object <- subset(
            obj$seurat_object,
            cells = valid_seurat_cells
        )
    }

    # 3b. Reconstructing imputed genotype matrix
    mat_imputed_t <- t(mat_imputed)
    full_imputed <- matrix(
        3L,
        nrow = nrow(obj$mae[["dna_variants"]]),
        ncol = length(cells_to_keep)
    )
    rownames(full_imputed) <- rownames(obj$mae[["dna_variants"]])
    colnames(full_imputed) <- cells_to_keep
    # full_imputed[res$targets, ] <- mat_imputed_t
    full_imputed[target_vars, ] <- mat_imputed_t

    SummarizedExperiment::assay(
        obj$mae[["dna_variants"]],
        "compass_imputed"
    ) <- full_imputed

    # 3c. Clonal Assignment
    obj$dna_clones_renamed <- NULL

    # mat_imputed_for_clones <- t(full_imputed[res$targets, , drop = FALSE])
    mat_imputed_for_clones <- t(full_imputed[target_vars, , drop = FALSE])

    if (is.null(obj$variants.filtered)) {
        target_variants_df <- as.data.frame(
            SummarizedExperiment::rowData(obj$mae[["dna_variants"]])[target_vars, ]
        )
    } else {
        target_variants_df <- obj$variants.filtered
    }

    clustering_res <- generate_clonal_labels(
        ngt_matrix = mat_imputed_for_clones,
        target_variants_df = target_variants_df,
        ignore_missing = FALSE
    )

    obj$dna.clones <- setNames(
        as.factor(clustering_res$cell_metadata$clonal_cluster_id),
        clustering_res$cell_metadata$cell_barcode
    )

    obj$dna_clone_colors <- generate_clone_palette(obj$dna.clones)

    tree_gv_path <- res$tree_dot
    tree_dot_content <- if (file.exists(tree_gv_path)) {
        paste(readLines(tree_gv_path, warn = FALSE), collapse = "\n")
    } else NULL

    S4Vectors::metadata(obj$mae)$compass <- list(
        tree_dot = tree_dot_content
    )

    return(obj)
}
