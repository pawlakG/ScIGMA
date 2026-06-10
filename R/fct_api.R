#' Fetch Clinical VEP Annotations
#'
#' @param custom_variant_vector Character vector of variants (e.g. "chr1-115256669-G-A")
#' @param genome_build Genome build, default "grch37"
#' @return A tibble with annotations
#' @noRd
fetch_clinical_vep_annotations <- function(custom_variant_vector, genome_build = "grch37") {                                                       
  
  # DICTIONNAIRE TRANSLATION (3-lettres -> 1-lettre clinique)
  aa_map <- c(
    "Ala"="A", "Arg"="R", "Asn"="N", "Asp"="D", "Cys"="C", 
    "Gln"="Q", "Glu"="E", "Gly"="G", "His"="H", "Ile"="I", 
    "Leu"="L", "Lys"="K", "Met"="M", "Phe"="F", "Pro"="P", 
    "Ser"="S", "Thr"="T", "Trp"="W", "Tyr"="Y", "Val"="V", 
    "Ter"="*"
  )

  # 1. Parsing HGVS Input (Robust parsing for MAE native ids like 'chr1:115256669:G/A')                                                                                                                          
  hgvs_queries <- custom_variant_vector %>%
    as.character() %>%
    chartr(":/", "--", .) %>%
    stringr::str_remove("^chr") %>%                                                                                                                         
    stringr::str_split("-") %>%                                                                                                                             
    purrr::map_chr(~ paste0(.x[1], ":g.", .x[2], .x[3], ">", .x[4]))                                                                                      
  
  # 2. Configuration Endpoint
  server <- if (genome_build == "grch37") "https://grch37.rest.ensembl.org" else "https://rest.ensembl.org"                                        
  endpoint <- "/vep/human/hgvs?refseq=1&hgvs=1&variant_class=1&canonical=1&phenotypes=1"                                                                        
  
  # 3. Batch Request                                                                                                                               
  body_json <- jsonlite::toJSON(list(hgvs_notations = hgvs_queries), auto_unbox = TRUE)                                                                      
  
  response <- httr::POST(                                                                                                                                
    url = paste0(server, endpoint),                                                                                                                
    httr::content_type("application/json"),                                                                                                              
    httr::accept("application/json"),                                                                                                                    
    body = body_json                                                                                                                               
  )                                                                                                                                                
  
  if (httr::status_code(response) != 200) stop(sprintf(" [FATAL] API VEP Error %s: %s", httr::status_code(response), httr::content(response, "text")))             
  
  # 4. Flattening & Filtrage Initial                                                                                                               
  raw_data <- httr::content(response, "parsed", simplifyVector = TRUE)                                                                                   
  
  annotation_table <- raw_data %>%                                                                                                                 
    tibble::as_tibble() %>%                                                                                                                                
    tidyr::unnest(transcript_consequences, keep_empty = TRUE, names_repair = "unique") %>%                                                                
    dplyr::filter(canonical == 1)                                                                                                                         
  
  # ----------------------------------------------------------------------------                                                                   
  # ----------------------------------------------------------------------------                                                                   
  expected_cols <- c("hgvsc", "hgvsp")                                                                             
  for (col in expected_cols) {                                                                                                                     
    if (!col %in% names(annotation_table)) annotation_table[[col]] <- NA_character_                                                                                                     
  }
  
  if (!"colocated_variants" %in% names(annotation_table)) {
    annotation_table$colocated_variants <- list(NULL)
  }
  
  # 5. Data Engineering & Mapping 1-Lettre                                                                                                                
  annotation_table <- annotation_table %>%                                                                                                         
    dplyr::mutate(                                                                                                                                        
      consequence_terms = purrr::map_chr(consequence_terms, ~ if (length(.x) > 0) .x[1] else NA_character_),
      cDNA = dplyr::if_else(!is.na(hgvsc), stringr::str_extract(hgvsc, "c\\..+"), NA_character_),
      
      extracted_p = stringr::str_extract(hgvsp, "p\\..+"),
      extracted_p = stringr::str_replace_all(extracted_p, "%3D", "="),
      extracted_p = stringr::str_replace_all(extracted_p, aa_map),
      
      PROTEIN = dplyr::case_when(
        is.na(extracted_p) ~ NA_character_,
        is.na(gene_symbol) ~ extracted_p,
        TRUE ~ paste0(gene_symbol, ":", extracted_p)
      ),
      
      original_variant = custom_variant_vector[match(input, hgvs_queries)],
      
      CLINVAR = purrr::map2_chr(colocated_variants, original_variant, ~ {
        if (is.null(.x) || !is.data.frame(.x)) return(NA_character_)
        
        alt_allele <- stringr::str_extract(.y, "[^/]+$")
        
        if ("clin_sig_allele" %in% names(.x)) {
            pattern <- paste0("\\b", alt_allele, ":([^;]+)")
            extracted <- stringr::str_match(unlist(.x$clin_sig_allele), pattern)[, 2]
            extracted <- extracted[!is.na(extracted)]
            if (length(extracted) > 0) return(paste(unique(extracted), collapse = ","))
        }
        
        if ("allele_string" %in% names(.x) && "clin_sig" %in% names(.x)) {
            valid_rows <- stringr::str_detect(.x$allele_string, paste0("\\b", alt_allele, "\\b"))
            if (any(valid_rows, na.rm = TRUE)) {
                sigs <- unlist(.x$clin_sig[valid_rows])
                if (length(sigs) > 0) return(paste(unique(sigs[!is.na(sigs)]), collapse = ","))
            }
        }
        return(NA_character_)
      })
    ) %>%
    dplyr::filter(!is.na(cDNA)) %>%
    dplyr::mutate(
      consequence_terms = stringr::str_replace_all(consequence_terms, "_", " "),
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
  
  print("=== TABLE D'ANNOTATION ===")
  print(annotation_table)
  
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
filter_and_annotate_variants <- function(obj,
                                         paths,
                                         min_dp = 10,
                                         min_gq = 30,
                                         vaf_ref = 5,
                                         vaf_hom = 95,
                                         vaf_het = 30,
                                         min_cell_pt = 10,
                                         min_mut_cell_pt = 10,
                                         batch_size = 300) {

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
    annot_df <- tryCatch({
        fetch_clinical_vep_annotations(variant_ids)
    }, error = function(e) {
        stop(sprintf("API Error during variant annotation: %s", e$message))
    })

    if (!is.null(annot_df) && nrow(annot_df) > 0) {

        # Format variant_id like GENE:chrX:POS:REF/ALT for legacy compatibility downstream
        annot_df$variant_id <- ifelse(
            !is.na(annot_df$gene),
            paste0(annot_df$gene, ":", sub("^([^-]+)-([^-]+)-([^-]+)-([^-]+)$", "\\1:\\2:\\3/\\4", chartr(":/", "--", annot_df$original_variant))),
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
            warning("API returned multiple annotations per variant. Deduplicating to preserve MAE strict dimensions.")
            merged_rowdata <- merged_rowdata[!duplicated(merged_rowdata$query_id), ]
        }

        # 6. Canonical order restoration
        rownames(merged_rowdata) <- merged_rowdata$query_id
        merged_rowdata <- merged_rowdata[variant_ids, ]
        merged_rowdata$query_id <- NULL

        merged_rowdata$variant_id[is.na(merged_rowdata$variant_id)] <- paste0("Unmapped:",rownames(merged_rowdata))[is.na(merged_rowdata$variant_id)]

        SummarizedExperiment::rowData(obj$mae[["dna_variants"]]) <- S4Vectors::DataFrame(merged_rowdata)
        message("Annotation matrix successfully integrated into MAE rowData.")

    } else {
        warning("API returned an empty object. rowData remains unannotated.")
    }

    message("Calculating variant cell proportions...")

    vaf_mtx <- SummarizedExperiment::assay(obj$mae[["dna_variants"]], "vaf")

    mutated_cells_count <- DelayedMatrixStats::rowSums2(vaf_mtx > 0)

    total_cells <- ncol(vaf_mtx)

    SummarizedExperiment::rowData(obj$mae[["dna_variants"]])$cell_proportion <- round(mutated_cells_count / total_cells, 2)

    message("Cell proportions successfully added to rowData.")
    S4Vectors::metadata(obj$mae)$variant_filter <- "filtered"

    invisible(obj)
}
