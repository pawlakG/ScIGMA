#' analysis_right_DNA UI Function
#'
#' @description A shiny Module.
#'
#' @param id,input,output,session Internal parameters for {shiny}.
#'
#' @noRd
#'
#' @importFrom DT DTOutput datatable renderDT
mod_analysis_DNA_ui <- function(id) {
    ns <- NS(id)
    tagList(
        navset_card_underline(
            id = ns("dna_tabs"),
            nav_panel(
                "Variant selection",
                accordion(
                    id = ns("acc"),
                    open = FALSE,
                    accordion_panel(
                        "Select DNA variants",
                        DTOutput(ns("variant_selection")),
                        br(),
                        fluidRow(
                            actionButton(
                                ns("btn_filtrer"),
                                "Apply",
                                class = "btn-primary"
                            )
                        )
                    ),
                    accordion_panel(
                        "DNA variant heatmap",
                        fluidRow(
                            div(
                                plotOutput(
                                    ns("dna_variant_heatmap"),
                                    height = "600px",
                                    width = "900px"
                                ),
                                align = "center"
                            ),
                            br(),
                            fluidRow(
                                column(
                                    4,
                                    div(
                                        shinyWidgets::materialSwitch(
                                            inputId = ns(
                                                "heatmap_include_all_samples"
                                            ),
                                            label = "Show missings ?",
                                            value = TRUE,
                                            status = "success"
                                        ),
                                        align = "left"
                                    )
                                ),
                                column(
                                    4,
                                    div(
                                        shinyWidgets::materialSwitch(
                                            inputId = ns(
                                                "heatmap_use_compass_imputed"
                                            ),
                                            label = "Use COMPASS imputed matrix ?",
                                            value = FALSE,
                                            # status = "primary"
                                            status = "success"
                                        ),
                                        align = "left"
                                    )
                                ),
                                column(
                                    4,
                                    div(
                                        actionButton(
                                            ns("btn_dna_variant_download"),
                                            "Download plot",
                                            class = "btn-primary"
                                        ),
                                        align = "center"
                                    )
                                )
                            )
                        )
                    ),
                    accordion_panel(
                        "Rename clusters",
                        uiOutput(ns("rename_cluster_ui"))
                    ),
                )
            ),
            nav_panel(
                title = "COMPASS",
                value = "compass_tab",
                bslib::card(
                    fill = FALSE,
                    bslib::card_header(
                        shiny::icon("code-branch"),
                        "COMPASS: Clonal Architecture Inference",
                        class = "bg-dark text-white"
                    ),
                    bslib::card_body(
                        fillable = FALSE,
                        fill = FALSE,
                        shiny::p(
                            "COMPASS (COpy number and Mutation Phylogeny from
                            Amplicon Single-cell Sequencing) is a probabilistic
                            framework designed to infer high-resolution clonal
                            architecture from single-cell DNA sequencing."
                        ),
                        shiny::p(
                            "By jointly modeling somatic mutations (SNVs) and
                            copy number alterations (CNAs) through Markov Chain
                            Monte Carlo (MCMC) inference, the algorithm
                            rigorously resolves technical noise and performs
                            probabilistic imputation of allele dropouts
                            (missing data)."
                        ),
                        shiny::p(
                            shiny::strong("Utility in the current pipeline: "),
                            "This module computes the definitive genetic
                            backbone of the sample. By assigning each
                            individual cell to a mathematically validated
                            clonal lineage, it provides the strict ground
                            truth required to project mutational profiles onto
                            the downstream protein expression space (UMAP),
                            enabling robust genotype-to-phenotype mapping."
                        )
                    ),
                    br(),
                    helpText(
                        "Warning: Running COMPASS inference will reset the
                            clonal architecture and erase any custom clone
                            names. Please rename your clones only after your
                            final architecture is computed."
                    ),
                    shiny::br(),
                    fluidRow(
                        column(
                            4,
                            shiny::radioButtons(
                                inputId = ns("patient_sex"),
                                label = "Sex",
                                choices = c(
                                    "Female" = "female",
                                    "Male" = "male"
                                ),
                                selected = "female",
                                inline = TRUE
                            )
                        ),
                        column(
                            4,
                            sliderTextInput(
                                inputId = ns("run_compass_n_chains"),
                                label = "Number of Markov Chains",
                                choices = seq_len(8),
                                grid = TRUE,
                                selected = 4,
                                width = "100%"
                            )
                        ),
                        column(
                            4,
                            sliderTextInput(
                                inputId = ns("run_compass_length_chains"),
                                label = "Markov Chains length",
                                choices = seq(100, 5000, 50),
                                grid = TRUE,
                                selected = 800,
                                width = "100%"
                            )
                        ),
                        align = "center"
                    ),
                    shiny::br(),
                    shiny::actionButton(
                        inputId = ns("btn_run_compass"),
                        label = "Run MCMC Inference (Async)",
                        icon = shiny::icon("play"),
                        class = "btn-primary"
                    ),
                    shiny::br(),
                    shiny::br(),
                    shiny::uiOutput(ns("compass_tree_ui"))
                )
            )
        )
    )
}

#' analysis_right_DNA Server Functions
#'
#' @noRd
#'
#'
#' @import InteractiveComplexHeatmap
#' @importFrom ComplexHeatmap draw
#' @importFrom forcats fct_recode
mod_analysis_DNA_server <- function(id, ScIGMA_data) {
    moduleServer(id, function(input, output, session) {
        ns <- session$ns

        compass_tree_visible <- shiny::reactiveVal(FALSE)

        bslib::nav_hide(id = "dna_tabs", target = "compass_tab")

        observeEvent(
            input$btn_filtrer,
            {
                req(ScIGMA_data$mae)

                ScIGMA_data$variants.filtered <- NULL
                ScIGMA_data$dna_clones_renamed <- NULL

                output$rename_cluster_ui <- renderUI({
                    NULL
                })

                if (!is.null(S4Vectors::metadata(ScIGMA_data$mae)$compass)) {
                    S4Vectors::metadata(ScIGMA_data$mae)$compass <- NULL
                    if ("cnv.active.clones" %in% names(ScIGMA_data)) {
                        ScIGMA_data$cnv.active.clones <- NULL
                    }

                    shinyWidgets::updateMaterialSwitch(
                        session,
                        "heatmap_use_compass_imputed",
                        value = FALSE
                    )
                    compass_tree_visible(FALSE)

                    shiny::showNotification(
                        "Modified variants: Previous COMPASS model purged.",
                        type = "warning"
                    )
                }

                sel_indices <- input$variant_selection_rows_selected
                if (length(sel_indices) > 0) {
                    ScIGMA_profile("2. Variant filtering and annotation", {
                        sorted_annotation <- SummarizedExperiment::rowData(ScIGMA_data$mae[[
                            "dna_variants"
                        ]]) |>
                            as.data.frame() |>
                            dplyr::arrange(desc(cell_proportion), clinvar)

                        selected_df <- sorted_annotation[
                            sel_indices,
                            ,
                            drop = FALSE
                        ]
                        selected_df$label <- paste(
                            selected_df$protein,
                            selected_df$cdna,
                            sep = " / "
                        ) # Add label
                        ScIGMA_data$variants.filtered <- selected_df

                        short_vars <- rownames(selected_df)
                        gt_raw <- t(as.matrix(SummarizedExperiment::assay(
                            ScIGMA_data$mae[["dna_variants"]],
                            "gt"
                        ))[short_vars, , drop = FALSE])
                        msk_raw <- t(as.matrix(SummarizedExperiment::assay(
                            ScIGMA_data$mae[["dna_variants"]],
                            "variant_filter_mask"
                        ))[short_vars, , drop = FALSE]) !=
                            0

                        # Apply mask (3 = Missing/Dropout)
                        gt_raw[cbind(
                            row(msk_raw)[msk_raw],
                            col(msk_raw)[msk_raw]
                        )] <- 3L
                        gt_complete_cells <- gt_raw[
                            rowSums(gt_raw == 3L) == 0,
                            ,
                            drop = FALSE
                        ]

                        if (nrow(gt_complete_cells) > 0) {
                            res_raw <- generate_clonal_labels(
                                gt_complete_cells,
                                selected_df,
                                ignore_missing = TRUE
                            )
                            ScIGMA_data$dna.clones <- setNames(
                                as.factor(
                                    res_raw$cell_metadata$clonal_cluster_id
                                ),
                                res_raw$cell_metadata$cell_barcode
                            )
                            ScIGMA_data$dna.clones_pre_compass <- ScIGMA_data$dna.clones
                            ScIGMA_data$dna_clone_colors <- generate_clone_palette(
                                ScIGMA_data$dna.clones
                            )

                            bslib::nav_show(
                                id = "dna_tabs",
                                target = "compass_tab"
                            )
                        }
                    })
                } else {
                    bslib::nav_hide(id = "dna_tabs", target = "compass_tab")
                }
            },
            priority = 10,
            ignoreInit = TRUE
        )

        observeEvent(
            watch("dataLoaded"),
            {
                shinyWidgets::updateMaterialSwitch(
                    session,
                    "heatmap_use_compass_imputed",
                    value = FALSE
                )
                if (exists("compass_tree_visible")) {
                    compass_tree_visible(FALSE)
                }

                # Clear graphs
                output$dna_variant_heatmap <- renderPlot({
                    NULL
                })
                output$rename_cluster_ui <- renderUI({
                    NULL
                })

                bslib::nav_hide(id = "dna_tabs", target = "compass_tab")
            },
            priority = 20,
            ignoreInit = TRUE
        )

        # 1. Render DNA variants dataframe
        output$variant_selection <- renderDT({
            watch("dnaVariant_filtered")
            watch("dataLoaded")
            req(ScIGMA_data$mae)

            tmp_variant_annotation <- SummarizedExperiment::rowData(ScIGMA_data$mae[[
                "dna_variants"
            ]]) |>
                as.data.frame() |>
                dplyr::select(
                    gene,
                    transcript_id,
                    protein,
                    cdna,
                    variant_type,
                    gene_function,
                    clinvar,
                    cell_proportion
                ) |>
                # dplyr::filter(!is.na(variant_id)) |>
                dplyr::arrange(desc(cell_proportion), clinvar)

            datatable(
                tmp_variant_annotation,
                selection = "multiple",
                rownames = FALSE,
                colnames = stringr::str_replace_all(
                    names(tmp_variant_annotation),
                    "_",
                    " "
                ),
                options = list(
                    pageLength = 5,
                    lengthMenu = c(5, 10, 15)
                )
            )
        })

        observeEvent(
            {
                input$btn_filtrer
                input$heatmap_include_all_samples
                input$heatmap_use_compass_imputed
                watch("dna_clones_renamed")
                watch("compass_completed")
            },
            {
                sel_indices <- input$variant_selection_rows_selected
                heatmap_include_all_samples <- input$heatmap_include_all_samples
                use_imputed <- input$heatmap_use_compass_imputed

                if (
                    length(sel_indices) > 0 &&
                        max(sel_indices) <=
                            nrow(ScIGMA_data$mae[["dna_variants"]])
                ) {
                    if (
                        isTRUE(use_imputed) &&
                            is.null(
                                S4Vectors::metadata(ScIGMA_data$mae)$compass
                            )
                    ) {
                        shiny::showNotification(
                            "COMPASS inference missing. Please run COMPASS first.",
                            type = "warning"
                        )

                        shinyWidgets::updateMaterialSwitch(
                            session,
                            "heatmap_use_compass_imputed",
                            value = FALSE
                        )
                        use_imputed <- FALSE
                    }

                    sorted_annotation <- SummarizedExperiment::rowData(ScIGMA_data$mae[[
                        "dna_variants"
                    ]]) |>
                        as.data.frame() |>
                        dplyr::select(
                            dplyr::any_of("variant_id"),
                            gene,
                            transcript_id,
                            protein,
                            cdna,
                            variant_type,
                            gene_function,
                            clinvar,
                            cell_proportion
                        ) |>
                        dplyr::arrange(desc(cell_proportion), clinvar)

                    selected_df <- sorted_annotation[
                        sel_indices,
                        ,
                        drop = FALSE
                    ]
                    selected_df$label <- paste(
                        selected_df$protein,
                        selected_df$cdna,
                        sep = " / "
                    ) # Add label
                    ScIGMA_data$variants.filtered <- selected_df

                    # print("selected_df")
                    # print(selected_df)

                    ht_res <- generate_dna_variant_heatmap(
                        obj = ScIGMA_data,
                        selected_variants_df = selected_df,
                        heatmap_include_all_samples = heatmap_include_all_samples,
                        use_imputed = use_imputed
                    )

                    ht <- ComplexHeatmap::draw(ht_res$heatmap)

                    if (is.null(ScIGMA_data$dna_clones_renamed)) {
                        is_compass_locked <- !is.null(
                            S4Vectors::metadata(ScIGMA_data$mae)$compass
                        )

                        if (!is_compass_locked) {
                            new_clones <- ht_res$clones

                            ScIGMA_data$dna.clones <- new_clones

                            if (!isTRUE(use_imputed)) {
                                ScIGMA_data$dna.clones_pre_compass <- new_clones
                            }

                            ScIGMA_data$dna_clone_colors <- generate_clone_palette(
                                ScIGMA_data$dna.clones
                            )

                            ScIGMA_data$cnv_dp_filtered <- NULL
                        }
                    }

                    shiny::showNotification(
                        ui = paste0(
                            "Success: ",
                            length(sel_indices),
                            " DNA variants successfully selected and processed."
                        ),
                        type = "message",
                        duration = 5
                    )

                    trigger("dnaVariant_selected")

                    output$dna_variant_heatmap <- renderPlot({
                        ht
                    })
                }
            },
            ignoreInit = TRUE
        )

        observeEvent(watch("dnaVariant_selected"), {
            req(ScIGMA_data$dna.clones)
            output$rename_cluster_ui <- renderUI({
                tagList(
                    p(
                        "Here you can rename a cluster, select a cluster name on drop list on the left, write its new name in right box and click on 'Apply New Labels' button."
                    ),
                    fluidRow(
                        column(
                            6,
                            pickerInput(
                                inputId = ns("rename_cluster_ui_oldName"),
                                label = "Style : primary",
                                choices = levels(ScIGMA_data$dna.clones),
                                options = pickerOptions(
                                    container = "body",
                                    style = "btn-outline-primary"
                                ),
                                width = "100%"
                            )
                        ),
                        column(
                            6,
                            textInput(
                                ns("rename_cluster_ui_newName"),
                                "New name"
                            )
                        )
                    ),
                    div(
                        actionButton(
                            inputId = ns("btn_update_cluster_labels"),
                            label = "Apply New Labels",
                            icon = icon("check"),
                            class = "btn-primary w-100"
                        ),
                        style = "margin-top:10px;"
                    )
                )
            })
        })

        observeEvent(input$btn_update_cluster_labels, {
            req(
                ScIGMA_data$dna.clones,
                input$rename_cluster_ui_oldName,
                input$rename_cluster_ui_newName
            )

            ScIGMA_data$update_dna_clone_names(
                old_name = input$rename_cluster_ui_oldName,
                new_name = input$rename_cluster_ui_newName
            )

            output$rename_cluster_ui <- renderUI({
                tagList(
                    p(
                        "Here you can rename a cluster, select a cluster name on drop list on the left, write its new name in right box and click on 'Apply New Labels' button."
                    ),
                    fluidRow(
                        column(
                            6,
                            pickerInput(
                                inputId = ns("rename_cluster_ui_oldName"),
                                label = "Style : primary",
                                choices = levels(ScIGMA_data$dna.clones),
                                options = pickerOptions(
                                    container = "body",
                                    style = "btn-outline-primary"
                                ),
                                width = "100%"
                            )
                        ),
                        column(
                            6,
                            textInput(
                                ns("rename_cluster_ui_newName"),
                                "New name"
                            )
                        )
                    ),
                    div(
                        actionButton(
                            inputId = ns("btn_update_cluster_labels"),
                            label = "Apply New Labels",
                            icon = icon("check"),
                            class = "btn-primary w-100"
                        ),
                        style = "margin-top:10px;"
                    )
                )
            })

            trigger("dna_clones_renamed")
        })

        # Dynamic UI rendering
        output$compass_tree_ui <- shiny::renderUI({
            shiny::req(compass_tree_visible())

            bslib::card(
                fill = FALSE,
                bslib::card_header(
                    div(
                        class = "d-flex justify-content-between align-items-center",
                        span(
                            shiny::icon("project-diagram"),
                            " Maximum Likelihood Phylogeny"
                        )
                    ),
                    class = "bg-dark text-white"
                ),
                bslib::card_body(
                    fillable = FALSE,
                    fill = FALSE,
                    shinycssloaders::withSpinner(
                        DiagrammeR::grVizOutput(
                            ns("compass_tree_plot"),
                            width = "100%",
                            height = "500px"
                        ),
                        type = 4,
                        color = "#007bff"
                    )
                ),
                bslib::card_footer(
                    shiny::downloadButton(
                        outputId = ns("btn_download_tree"),
                        label = "Download Tree (.svg)",
                        class = "btn-sm btn-light"
                    )
                )
            )
        })

        compass_task <- shiny::ExtendedTask$new(
            function(
                compass_length_chains,
                compass_n_chains,
                patient_sex,
                variant_matrices,
                vec_locus_regions,
                mat_cna,
                vec_locus_names,
                vec_locus_chrom,
                vec_region_names,
                vec_region_chrom,
                use_cna,
                target_vars
            ) {
                # CRITICAL FIX: future_promise() instead of future()
                promises::future_promise(
                    {
                        if (!isNamespaceLoaded("ScIGMA")) {
                            pkgload::load_all(export_all = FALSE, quiet = TRUE)
                        }

                        prefix_out <- file.path(
                            tempdir(),
                            paste0("compass_async_", as.integer(Sys.time()))
                        )

                        ScIGMA_profile("3. COMPASS algorithm", {
                            run_compass_mcmc(
                                variant_matrices = variant_matrices,
                                locus_regions = vec_locus_regions,
                                region_matrix = mat_cna,
                                output_prefix = prefix_out,
                                locus_names = vec_locus_names,
                                locus_chromosomes = vec_locus_chrom,
                                region_names = vec_region_names,
                                region_chromosomes = vec_region_chrom,
                                chains = as.integer(compass_n_chains),
                                chain_length = compass_length_chains,
                                patient_sex = patient_sex,
                                use_cna = use_cna
                            )
                        })

                        return(list(prefix = prefix_out, targets = target_vars))
                    },
                    seed = TRUE
                )
            }
        )

        observeEvent(
            input$btn_run_compass,
            {
                req(ScIGMA_data$mae)

                shinyWidgets::confirmSweetAlert(
                    session = session,
                    inputId = ns("confirm_run_compass"),
                    type = "warning",
                    title = "MCMC Inference Initialization",
                    text = shiny::HTML(
                        "<b>Warning:</b> Running COMPASS will completely recalculate the clonal architecture.<br><br>Any custom clone names you have defined will be <b>permanently erased and reset</b>.<br><br>Do you want to proceed?"
                    ),
                    html = TRUE,
                    btn_labels = c("Cancel", "Run COMPASS"),
                    btn_colors = c("#d3d3d3", "#007bff")
                )
            },
            ignoreInit = TRUE
        )

        observeEvent(
            input$confirm_run_compass,
            {
                req(isTRUE(input$confirm_run_compass))

                req(ScIGMA_data$mae)
                req(input$run_compass_length_chains)

                compass_tree_visible(TRUE)

                shinyjs::disable("btn_run_compass")
                shiny::showNotification(
                    "1/2 - Extracting HDF5 and preparing matrices...",
                    id = "compass_notif",
                    duration = NULL,
                    type = "message"
                )

                if (is.null(ScIGMA_data$variants.filtered)) {
                    target_vars <- rownames(ScIGMA_data$mae[["dna_variants"]])
                } else {
                    target_vars <- rownames(ScIGMA_data$variants.filtered)
                }

                compass_length_chains <- input$run_compass_length_chains
                compass_n_chains <- input$run_compass_n_chains
                patient_sex <- input$patient_sex

                compass_inputs <- build_compass_matrices(
                    obj = ScIGMA_data,
                    selected_variants = target_vars
                )

                mat_ref <- t(as.matrix(compass_inputs$M_ref))
                mat_alt <- t(as.matrix(compass_inputs$M_alt))
                mat_cna <- t(as.matrix(compass_inputs$C))

                gt_assay <- SummarizedExperiment::assay(
                    ScIGMA_data$mae[["dna_variants"]],
                    "gt"
                )
                mat_gt <- as.matrix(gt_assay[target_vars, , drop = FALSE])
                mat_gt[mat_gt == 3L] <- NA

                storage.mode(mat_ref) <- "integer"
                storage.mode(mat_alt) <- "integer"
                storage.mode(mat_cna) <- "integer"
                storage.mode(mat_gt) <- "integer"

                variant_matrices <- list(
                    REF = mat_ref,
                    ALT = mat_alt,
                    GT = mat_gt
                )

                dna_se <- ScIGMA_data$mae[["dna_variants"]]
                snv_sub <- as.data.frame(SummarizedExperiment::rowData(dna_se))[
                    target_vars,
                ]
                vec_locus_names <- snv_sub$gene
                vec_locus_chrom <- snv_sub$chrom

                amp_se <- ScIGMA_data$mae[["amplicons"]]
                cna_row_data <- as.data.frame(SummarizedExperiment::rowData(
                    amp_se
                ))
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

                use_cna <- if (ncol(variant_matrices$REF) != ncol(mat_cna)) {
                    FALSE
                } else {
                    TRUE
                }

                vec_locus_regions <- as.integer(compass_inputs$locus_regions)

                shiny::showNotification(
                    "2/2 - Background MCMC. Session unlocked, you can run other analyses.",
                    id = "compass_notif",
                    type = "warning",
                    duration = NULL
                )

                compass_task$invoke(
                    compass_length_chains,
                    compass_n_chains,
                    patient_sex,
                    variant_matrices,
                    vec_locus_regions,
                    mat_cna,
                    vec_locus_names,
                    vec_locus_chrom,
                    vec_region_names,
                    vec_region_chrom,
                    use_cna,
                    target_vars
                )
            },
            ignoreInit = TRUE
        )

        observeEvent(
            compass_task$status(),
            {
                status <- compass_task$status()

                if (status == "success") {
                    res <- compass_task$result()

                    mat_imputed <- get_imputed_genotypes(
                        prefix_out = res$prefix
                    )
                    cells_to_keep <- rownames(mat_imputed)

                    if (length(cells_to_keep) == 0) {
                        stop("COMPASS filtered all cells (pure doublets).")
                    }

                    ScIGMA_data$mae <- ScIGMA_data$mae[, cells_to_keep, ]
                    if (!is.null(ScIGMA_data$seurat_object)) {
                        valid_seurat_cells <- intersect(
                            cells_to_keep,
                            colnames(ScIGMA_data$seurat_object)
                        )
                        ScIGMA_data$seurat_object <- subset(
                            ScIGMA_data$seurat_object,
                            cells = valid_seurat_cells
                        )
                    }

                    mat_imputed_t <- t(mat_imputed)

                    full_imputed <- matrix(
                        3L,
                        nrow = nrow(ScIGMA_data$mae[["dna_variants"]]),
                        ncol = length(cells_to_keep)
                    )
                    rownames(full_imputed) <- rownames(ScIGMA_data$mae[[
                        "dna_variants"
                    ]])
                    colnames(full_imputed) <- cells_to_keep

                    full_imputed[res$targets, ] <- mat_imputed_t

                    SummarizedExperiment::assay(
                        ScIGMA_data$mae[["dna_variants"]],
                        "compass_imputed"
                    ) <- full_imputed

                    ScIGMA_data$dna_clones_renamed <- NULL

                    mat_imputed_for_clones <- t(full_imputed[
                        res$targets,
                        ,
                        drop = FALSE
                    ])

                    clustering_res <- generate_clonal_labels(
                        ngt_matrix = mat_imputed_for_clones,
                        target_variants_df = ScIGMA_data$variants.filtered,
                        ignore_missing = FALSE
                    )

                    ScIGMA_data$dna.clones <- setNames(
                        as.factor(
                            clustering_res$cell_metadata$clonal_cluster_id
                        ),
                        clustering_res$cell_metadata$cell_barcode
                    )

                    ScIGMA_data$dna_clone_colors <- generate_clone_palette(
                        ScIGMA_data$dna.clones
                    )

                    shinyWidgets::updateMaterialSwitch(
                        session,
                        "heatmap_use_compass_imputed",
                        value = TRUE
                    )

                    tree_gv_path <- paste0(res$prefix, "_tree.gv")
                    tree_dot_content <- if (file.exists(tree_gv_path)) {
                        paste(readLines(tree_gv_path), collapse = "\n")
                    } else {
                        NULL
                    }

                    S4Vectors::metadata(ScIGMA_data$mae)$compass <- list(
                        singlet_barcodes = cells_to_keep,
                        target_variants = res$targets,
                        tree_dot = tree_dot_content
                    )

                    shiny::removeNotification(id = "compass_notif")
                    shiny::showNotification(
                        "COMPASS finished: Clonal architecture and Multi-Omics Matrix locked.",
                        duration = 10,
                        type = "message"
                    )
                    shinyjs::enable("btn_run_compass")

                    trigger("compass_completed")
                } else if (status == "error") {
                    shiny::removeNotification(id = "compass_notif")
                    shiny::showNotification(
                        paste(
                            "Critical C++ failure:",
                            compass_task$error()$message
                        ),
                        duration = NULL,
                        type = "error"
                    )
                    shinyjs::enable("btn_run_compass")
                }
            },
            ignoreInit = TRUE
        )

        output$compass_tree_plot <- DiagrammeR::renderGrViz({
            watch("compass_completed")
            req(ScIGMA_data$mae)

            tree_content <- S4Vectors::metadata(
                ScIGMA_data$mae
            )$compass$tree_dot
            req(tree_content)

            DiagrammeR::grViz(tree_content)
        })

        output$btn_download_tree <- shiny::downloadHandler(
            filename = function() {
                paste0(
                    "ScIGMA_compass_phylogeny_",
                    format(Sys.time(), "%Y%m%d_%H%M%S"),
                    ".svg"
                )
            },
            content = function(file) {
                req(ScIGMA_data$mae)
                tree_content <- S4Vectors::metadata(
                    ScIGMA_data$mae
                )$compass$tree_dot
                req(tree_content)

                # Native conversion of DOT to a pure Vector format
                svg_code <- DiagrammeRsvg::export_svg(DiagrammeR::grViz(
                    tree_content
                ))

                writeLines(svg_code, file)
            }
        )
    })
}
