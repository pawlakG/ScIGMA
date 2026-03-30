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
                            actionButton(ns("btn_filtrer"), "Apply",
                                         class = "btn-primary")
                        )
                    ),
                    accordion_panel(
                        "DNA variant heatmap",
                        fluidRow(
                            div(
                                plotOutput(ns("dna_variant_heatmap"),
                                           height = "600px", width = "900px"),
                                align = "center"),
                            fluidRow(
                                column(6,
                                       materialSwitch(
                                           inputId = ns("heatmap_include_all_samples"),
                                           label = "Show missings ?",
                                           value = TRUE,
                                           status = "success"
                                       )
                                ),
                                column(6,
                                       actionButton(ns("btn_dna_variant_download"), "Download plot",
                                                    class = "btn-primary")
                                )
                            )
                        )
                    ),
                    accordion_panel(
                        "Rename clusters",
                        uiOutput(ns("rename_cluster_ui")
                        )
                    ),
                )
            ),
            nav_panel(
                "COMPASS",
                bslib::card(
                    bslib::card_header(
                        shiny::icon("code-branch"),
                        "COMPASS: Clonal Architecture Inference",
                        class = "bg-dark text-white" # À adapter selon ton thème ggprism/application
                    ),
                    bslib::card_body(
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
                    shiny::br(),

                    div(
                        sliderTextInput(
                            inputId = ns("run_compass_length_chains"),
                            label = "Markov Chains length",
                            choices = seq(100, 5000, 50),
                            grid = TRUE,
                            selected = 800,
                            width = "150%"
                        ),
                        align = "center"),
                    shiny::br(),
                    shiny::actionButton(
                        inputId = ns("btn_run_compass"),
                        label = "Run MCMC Inference (Async)",
                        icon = shiny::icon("play"),
                        class = "btn-primary"
                    ),
                    shiny::br(),
                    shiny::br(),

                    # FIX CRITIQUE : Zone d'injection dynamique pour l'arbre
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
mod_analysis_DNA_server <- function(id, ScIGMA_data){
    moduleServer(id, function(input, output, session){
        # UPDATED
        # File: R/mod_analysis_overview.R (ou fichier contenant ce module)

        ns <- session$ns

        # 1. Render DNA variants dataframe
        output$variant_selection <- renderDT({
            watch("dnaVariant_filtered")
            req(ScIGMA_data$mae) # Sécurité

            # print("test_renderDT ")
            # print(SummarizedExperiment::rowData(ScIGMA_data$mae[["dna_variants"]]) )
            # Extraction et tri (La "Vue")
            tmp_variant_annotation <- SummarizedExperiment::rowData(ScIGMA_data$mae[["dna_variants"]]) |>
                as.data.frame() |>
                dplyr::select(variant_id, gene, variant_type, gene_function, impact, clinvar, cell_proportion) |>
                dplyr::arrange(desc(cell_proportion), desc(impact))

            datatable(tmp_variant_annotation,
                      selection = 'multiple',
                      rownames = FALSE, # Désactivé car variant_id est déjà présent
                      options = list(pageLength = 5,
                                     lengthMenu = c(5, 10, 15)))
        })

        # 2. Récupérer les lignes sélectionnées
        observeEvent({
            input$btn_filtrer
            input$heatmap_include_all_samples
            watch("dna_clones_renamed")
        }, {
            print("Rendering DNA heatmap")
            sel_indices <- input$variant_selection_rows_selected
            heatmap_include_all_samples <- input$heatmap_include_all_samples

            if (length(sel_indices) > 0) {

                # RECONSTRUCTION DE LA VUE : Indispensable pour mapper les index de l'UI (sel_indices)
                # avec les véritables identifiants biologiques, car arrange() a mélangé les lignes.
                sorted_annotation <- SummarizedExperiment::rowData(ScIGMA_data$mae[["dna_variants"]]) |>
                    as.data.frame() |>
                    dplyr::select(variant_id, gene, variant_type, gene_function, impact, clinvar, cell_proportion) |>
                    dplyr::arrange(desc(cell_proportion), desc(impact))

                # Extraction sécurisée des variants sélectionnés
                selected_df <- sorted_annotation[sel_indices, , drop = FALSE]

                # Mise à jour de l'objet global (au cas où d'autres modules l'utilisent)
                ScIGMA_data$variants.filtered <- selected_df

                # Génération de la Heatmap
                ht_res <- generate_dna_variant_heatmap(
                    obj = ScIGMA_data,
                    selected_variants_df = selected_df,
                    heatmap_include_all_samples = heatmap_include_all_samples
                )

                ht <- ComplexHeatmap::draw(ht_res$heatmap)
                print("New heatmap rendered")

                # Gestion des clones
                if (is.null(ScIGMA_data$dna_clones_renamed)) {
                    ScIGMA_data$dna.clones <- ht_res$clones
                }

                # Déclenchement des événements avals
                trigger("dnaVariant_selected")

                # Affichage
                output$dna_variant_heatmap <- renderPlot({
                    ht
                })
            }
        })


        observeEvent(watch("dnaVariant_selected"),
                     {
                         req(ScIGMA_data$dna.clones)
                         output$rename_cluster_ui <-  renderUI({
                             tagList(
                                 p("Here you can rename a cluster, select a cluster name on drop list on the left, write its new name in right box and click on 'Apply New Labels' button."),
                                 fluidRow(column(6,
                                                 pickerInput(
                                                     inputId = ns("rename_cluster_ui_oldName"),
                                                     label = "Style : primary",
                                                     choices = levels(ScIGMA_data$dna.clones),
                                                     options = pickerOptions(container = "body",
                                                                             style = "btn-outline-primary"),
                                                     width = "100%"
                                                 )
                                 ),
                                 column(6,
                                        textInput(ns("rename_cluster_ui_newName"),
                                                  "New name")
                                 )
                                 ),
                                 div(
                                     actionButton(
                                         inputId = ns("btn_update_cluster_labels"),
                                         label = "Apply New Labels",
                                         icon = icon("check"),
                                         class = "btn-primary w-100" # w-100 pour prendre toute la largeur
                                     ), style = "margin-top:10px;"
                                 )
                             )
                         })
                     })

        observeEvent(input$btn_update_cluster_labels,
                     {
                         req(ScIGMA_data$dna.clones)
                         # update dna.clones labels
                         oldName <- input$rename_cluster_ui_oldName
                         newName <- input$rename_cluster_ui_newName
                         levels <- oldName
                         names(levels) <- newName
                         # ScIGMA_data$dna.clones <- fct_recode(ScIGMA_data$dna.clones, newName = oldName)
                         ScIGMA_data$dna.clones <- fct_recode(ScIGMA_data$dna.clones, !!!levels)
                         ScIGMA_data$dna_clones_renamed <- ScIGMA_data$dna.clones
                         trigger("dna_clones_renamed")
                     })

        # [ NODE_ACCESS : COMPASS ]
        # ----------------------------------------------------- _

        compass_tree_visible <- shiny::reactiveVal(FALSE)

        # Rendu dynamique de la carte UI
        output$compass_tree_ui <- shiny::renderUI({
            shiny::req(compass_tree_visible()) # Ne s'affiche que si l'interrupteur est TRUE

            bslib::card(
                bslib::card_header(
                    div(class = "d-flex justify-content-between align-items-center",
                        span(shiny::icon("project-diagram"), " Maximum Likelihood Phylogeny")
                    ),
                    class = "bg-dark text-white"
                ),
                bslib::card_body(
                    # Le spinner tournera ici jusqu'à ce que renderGrViz renvoie l'image
                    shinycssloaders::withSpinner(
                        DiagrammeR::grVizOutput(ns("compass_tree_plot"), width = "100%", height = "500px"),
                        type = 4, color = "#007bff"
                    )
                ),
                br(),
                shiny::downloadButton(
                    outputId = ns("btn_download_tree"),
                    label = "Download Tree (.svg)",
                    class = "btn-sm btn-light"
                )
            )
        })

        compass_task <- shiny::ExtendedTask$new(
            function(compass_length_chains, variant_matrices, vec_locus_regions, mat_cna,
                     vec_locus_names, vec_locus_chrom,
                     vec_region_names, vec_region_chrom, use_cna, target_vars) {

                # FIX CRITIQUE : future_promise() au lieu de future()
                # Cela garantit la création d'une Promesse compatible avec l'Event Loop de Shiny
                promises::future_promise({

                    # ---------------------------------------------------------
                    # FIX CRITIQUE : Le "Vaccin" Devtools (Solution 2)
                    # Si le worker s'éveille et ne trouve pas le package ScIGMA
                    # (car lancé via load_all()), il le charge lui-même.
                    # ---------------------------------------------------------
                    if (!isNamespaceLoaded("ScIGMA")) {
                        suppressMessages(pkgload::load_all(export_all = FALSE))
                    }

                    prefix_out <- file.path(tempdir(), paste0("compass_async_", as.integer(Sys.time())))

                    run_compass_mcmc(
                        variant_matrices   = variant_matrices,
                        locus_regions      = vec_locus_regions,
                        region_matrix      = mat_cna,
                        output_prefix      = prefix_out,
                        locus_names        = vec_locus_names,
                        locus_chromosomes  = vec_locus_chrom,
                        region_names       = vec_region_names,
                        region_chromosomes = vec_region_chrom,
                        chains             = 4L,
                        # chain_length       = 500L,
                        chain_length       = compass_length_chains,
                        patient_sex        = "female",
                        use_cna            = use_cna
                    )

                    return(list(prefix = prefix_out, targets = target_vars))
                }, seed = TRUE)
            }
        )


        observeEvent(input$btn_run_compass, {
            req(ScIGMA_data$mae)
            req(input$run_compass_length_chains)

            compass_tree_visible(TRUE)

            shinyjs::disable("btn_run_compass")
            shiny::showNotification("1/2 - Extraction HDF5 et préparation des matrices...",
                                    id = "compass_notif",
                                    duration = NULL,
                                    type = "message")

            print("ScIGMA_data$variants.filtered")
            print(ScIGMA_data$variants.filtered)

            # 2. EXTRACTION SYNCHRONE (Thread Principal)
            # Opération vitale : on lit le HDF5 ici, on ne passe que des objets RAM au worker
            if (is.null(ScIGMA_data$variants.filtered)){
                target_vars <- rownames(ScIGMA_data$mae[["dna_variants"]])
            } else {
                target_vars <- rownames(ScIGMA_data$variants.filtered)
            }

            compass_length_chains <- input$run_compass_length_chains

            compass_inputs <- build_compass_matrices(obj = ScIGMA_data, selected_variants = target_vars)

            # Conversion en matrices denses + Transposition
            mat_ref <- t(as.matrix(compass_inputs$M_ref))
            mat_alt <- t(as.matrix(compass_inputs$M_alt))
            mat_cna <- t(as.matrix(compass_inputs$C))

            gt_assay <- SummarizedExperiment::assay(ScIGMA_data$mae[["dna_variants"]], "gt")
            mat_gt <- as.matrix(gt_assay[target_vars, , drop = FALSE])
            mat_gt[mat_gt == 3L] <- NA

            storage.mode(mat_ref) <- "integer"
            storage.mode(mat_alt) <- "integer"
            storage.mode(mat_cna) <- "integer"
            storage.mode(mat_gt)  <- "integer"


            print("test_4")

            variant_matrices <- list(REF = mat_ref, ALT = mat_alt, GT = mat_gt)

            # Métadonnées
            dna_se <- ScIGMA_data$mae[["dna_variants"]]
            snv_sub <- as.data.frame(SummarizedExperiment::rowData(dna_se))[target_vars, ]
            vec_locus_names <- snv_sub$gene
            vec_locus_chrom <- snv_sub$chrom

            amp_se <- ScIGMA_data$mae[["amplicons"]]
            cna_row_data <- as.data.frame(SummarizedExperiment::rowData(amp_se))
            vec_region_names <- unique(paste0(cna_row_data$chrom, "_", sapply(cna_row_data$dna_id, \(x) strsplit(x, "_")[[1]][3])))
            vec_region_chrom <- sub("^chr", "", sapply(vec_region_names, \(x) strsplit(x, "_")[[1]][1], USE.NAMES = FALSE), ignore.case = TRUE)

            use_cna <- if (ncol(variant_matrices$REF) != ncol(mat_cna)) FALSE else TRUE

            vec_locus_regions <- as.integer(compass_inputs$locus_regions)

            # La SEULE notification 2/2 qui doit rester (avec l'ID en dur)
            shiny::showNotification(
                "2/2 - MCMC en arrière-plan. La session est débloquée, vous pouvez lancer d'autres analyses.",
                id = "compass_notif",
                type = "warning",
                duration = NULL
            )

            # L'appel magique : Lance le calcul et LIBÈRE IMMÉDIATEMENT la session UI
            compass_task$invoke(
                compass_length_chains, variant_matrices, vec_locus_regions, mat_cna,
                vec_locus_names, vec_locus_chrom,
                vec_region_names, vec_region_chrom, use_cna, target_vars
            )
        })

        # 3. LE RÉCEPTEUR SILENCIEUX
        observeEvent(compass_task$status(), {
            status <- compass_task$status()

            if (status == "success") {
                res <- compass_task$result()

                mat_imputed <- get_imputed_genotypes(prefix_out = res$prefix)
                cells_to_keep <- rownames(mat_imputed)

                if (length(cells_to_keep) == 0) stop("COMPASS a filtré toutes les cellules (doublets purs).")

                # Purge et intégration S4 (Isolé, pas de boucle réactive)
                ScIGMA_data$mae <- ScIGMA_data$mae[, cells_to_keep, ]

                mat_imputed_t <- t(mat_imputed)
                rownames(mat_imputed_t) <- res$targets
                colnames(mat_imputed_t) <- cells_to_keep

                # Récupération de l'arbre phylogénétique généré par le C++
                tree_gv_path <- paste0(res$prefix, "_tree.gv")
                tree_dot_content <- if (file.exists(tree_gv_path)) paste(readLines(tree_gv_path), collapse = "\n") else NULL

                S4Vectors::metadata(ScIGMA_data$mae)$compass <- list(
                    imputed_gt = mat_imputed_t,
                    singlet_barcodes = cells_to_keep,
                    target_variants = res$targets,
                    tree_dot = tree_dot_content # Nécessaire pour DiagrammeR
                )

                shiny::removeNotification(id = "compass_notif")
                shiny::showNotification("COMPASS terminé : Architecture clonale verrouillée.", duration = 10, type = "message")
                shinyjs::enable("btn_run_compass")

                trigger("compass_completed")

            } else if (status == "error") {
                shiny::removeNotification(id = "compass_notif")
                shiny::showNotification(paste("Échec critique C++ :", compass_task$error()$message), duration = NULL, type = "error")
                shinyjs::enable("btn_run_compass")
            }
        }, ignoreInit = TRUE)

        output$compass_tree_plot <- DiagrammeR::renderGrViz({
            # On s'assure que le calcul est terminé et que l'objet contient l'arbre
            watch("compass_completed")
            req(ScIGMA_data$mae)

            tree_content <- S4Vectors::metadata(ScIGMA_data$mae)$compass$tree_dot
            req(tree_content) # Ne s'affiche que si l'arbre existe

            DiagrammeR::grViz(tree_content)
        })

        # ---------------------------------------------------------
        # Gestionnaire de téléchargement Vectoriel (SVG Q1 Standard)
        # ---------------------------------------------------------
        output$btn_download_tree <- shiny::downloadHandler(
            filename = function() {
                paste0("ScIGMA_compass_phylogeny_", format(Sys.time(), "%Y%m%d_%H%M%S"), ".svg")
            },
            content = function(file) {
                req(ScIGMA_data$mae)
                tree_content <- S4Vectors::metadata(ScIGMA_data$mae)$compass$tree_dot
                req(tree_content)

                # Conversion native du DOT vers un format Vectoriel pur
                # Protège contre les crashs de librairies graphiques UNIX
                svg_code <- DiagrammeRsvg::export_svg(DiagrammeR::grViz(tree_content))

                writeLines(svg_code, file)
            }
        )

    })
}

## To be copied in the UI
# mod_analysis_DNA_ui("analysis_right_DNA_1")

## To be copied in the server
# mod_analysis_DNA_server("analysis_right_DNA_1", ScIGMA_data)
