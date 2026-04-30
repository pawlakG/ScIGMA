#' download_panel UI Function
#'
#' @description A shiny Module.
#'
#' @param id,input,output,session Internal parameters for {shiny}.
#'
#' @noRd
#'
#' @importFrom shiny NS tagList uiOutput div h4
mod_download_panel_ui <- function(id) {
    ns <- NS(id)
    tagList(
        shiny::uiOutput(ns("download_main_ui"))
    )
}

#' download_panel Server Functions
#'
#' @noRd
#' @importFrom zip zip
#' @importFrom utils write.csv
mod_download_panel_server <- function(id, ScIGMA_data) {
    moduleServer(id, function(input, output, session) {
        ns <- session$ns

        # ---------------------------------------------------------------------
        # 1. MAIN UI RENDERER (Empty State Logic)
        # ---------------------------------------------------------------------
        output$download_main_ui <- shiny::renderUI({
            # Déclencheurs pour réévaluer l'UI au fil de l'analyse
            watch("dataLoaded")
            watch("dnaVariant_filtered")
            watch("CNV_filtered")

            # Condition 1 : Aucune donnée chargée
            if (is.null(ScIGMA_data$mae)) {
                return(
                    shiny::div(class = "text-center mt-5", style = "padding: 50px;",
                               shiny::icon("file-import", class = "fa-4x text-muted", style = "margin-bottom: 20px;"),
                               shiny::h3("No data available for export", class = "text-muted"),
                               shiny::p("Please upload an H5 file and process it in the 'Analysis' tab first.", class = "text-muted")
                    )
                )
            }

            # --- AIGUILLAGE DYNAMIQUE DES CHOIX D'EXPORT ---
            # De base : ADN et Métadonnées (qui s'auto-enrichiront des clones, gates, clusters)
            choices_list <- c("Cell Metadata" = "metadata",
                              "DNA Genotypes" = "dna")

            # Après filtrage des variants et normalisation des protéines
            if (!is.null(ScIGMA_data$variants.filtered)) {
                choices_list <- c(choices_list,
                                  "DNA Variants Info" = "variants",
                                  "Proteins (Normalized)" = "proteins")
            }

            # Après filtrage CNV (Bouton "Filter" dans l'onglet CNV)
            if (!is.null(ScIGMA_data$cnv_dp_filtered)) {
                choices_list <- c(choices_list, "Amplicons (CNV)" = "cnv")
            }
            # -----------------------------------------------

            # Condition 2 : Les données sont présentes, on injecte l'UI
            shiny::tagList(
                bslib::card(
                    bslib::card_header(
                        shiny::icon("download"),
                        "Export ScIGMA Project",
                        class = "bg-primary text-white"
                    ),
                    bslib::card_body(
                        shiny::p("Download your processed data and metadata in standard bioinformatics formats. Options will unlock automatically as you progress through the analysis."),
                        shiny::br(),

                        shiny::fluidRow(
                            # --- CONFIGURATION (Left Column) ---
                            shiny::column(4,
                                          bslib::card(
                                              shiny::h4("Export Configuration"),
                                              shiny::hr(),

                                              # 1. Sélection des Omiques (Dynamique)
                                              shiny::h5("Include Omics Data:"),
                                              shinyWidgets::checkboxGroupButtons(
                                                  inputId = ns("export_omics_selection"),
                                                  choices = choices_list,
                                                  selected = unname(choices_list), # Tout cocher par défaut
                                                  direction = "vertical",
                                                  status = "outline-primary",
                                                  width = "100%"
                                              ),
                                              shiny::br(),

                                              # 2. Options Clonales
                                              shiny::h5("Clonal Metadata Architecture:"),
                                              shiny::uiOutput(ns("export_clones_options_ui")),
                                              shiny::br()
                                          )
                            ),

                            # --- EXPORT ACTIONS (Right Column) ---
                            shiny::column(8,
                                          # Option 1: Seurat
                                          bslib::card(
                                              shiny::fluidRow(
                                                  shiny::column(8,
                                                                shiny::h4(shiny::icon("r-project"), " Seurat Object (.rds)"),
                                                                shiny::p("The gold standard for single-cell analysis in R. Includes PCA, UMAP, and all selected metadata integrated natively.", class = "text-muted")
                                                  ),
                                                  shiny::column(4,
                                                                shiny::div(
                                                                    shiny::uiOutput(ns("btn_download_seurat_ui")),
                                                                    align = "right"
                                                                )
                                                  )
                                              )
                                          ),

                                          # Option 2: MultiAssayExperiment
                                          bslib::card(
                                              shiny::fluidRow(
                                                  shiny::column(8,
                                                                shiny::h4(shiny::icon("database"), " MultiAssayExperiment (.rds)"),
                                                                shiny::p("The native BioConductor format. Contains the raw and processed matrices exactly as managed in the ScIGMA backend.", class = "text-muted")
                                                  ),
                                                  shiny::column(4,
                                                                shiny::div(
                                                                    shiny::downloadButton(ns("btn_download_mae"), "Download MAE", class = "btn-info w-100"),
                                                                    align = "right"
                                                                )
                                                  )
                                              )
                                          ),

                                          # Option 3: Flat CSV
                                          bslib::card(
                                              shiny::fluidRow(
                                                  shiny::column(8,
                                                                shiny::h4(shiny::icon("file-csv"), " Flat Tables (.zip)"),
                                                                shiny::p("Universal format. Generates a ZIP archive containing individual CSV files for matrices and unified cell metadata. Ideal for Python/Pandas.", class = "text-muted")
                                                  ),
                                                  shiny::column(4,
                                                                shiny::div(
                                                                    shiny::downloadButton(ns("btn_download_csv"), "Download ZIP", class = "btn-secondary w-100"),
                                                                    align = "right"
                                                                )
                                                  )
                                              )
                                          )
                            )
                        )
                    )
                )
            )
        })

        # ---------------------------------------------------------------------
        # 2. DYNAMIC UI SUB-CONTROLS
        # ---------------------------------------------------------------------

        # A. Contrôle des clones (Bruts vs COMPASS)
        output$export_clones_options_ui <- shiny::renderUI({
            watch("dataLoaded")
            watch("compass_completed")

            # Sécurité pour éviter un crash si l'UI est évaluée alors que mae est encore NULL
            if (is.null(ScIGMA_data$mae)) return(NULL)

            compass_exists <- !is.null(S4Vectors::metadata(ScIGMA_data$mae)$compass)

            if (compass_exists) {
                shinyWidgets::materialSwitch(
                    inputId = ns("export_use_compass_clones"),
                    label = "Use COMPASS imputed clones",
                    value = TRUE,
                    status = "success"
                )
            } else {
                shiny::div(
                    shinyWidgets::materialSwitch(
                        inputId = ns("export_use_compass_clones"),
                        label = "Run COMPASS first to unlock imputation",
                        value = FALSE,
                        status = "default"
                    ),
                    style = "pointer-events: none; opacity: 0.5;"
                )
            }
        })

        # B. Activation conditionnelle du bouton Seurat
        output$btn_download_seurat_ui <- shiny::renderUI({
            watch("dataLoaded")
            watch("umap_computed")

            if (is.null(ScIGMA_data$seurat_object)) {
                shiny::actionButton(ns("btn_seurat_disabled"), "Download Seurat", class = "btn-success w-100 disabled", icon = shiny::icon("lock"))
            } else {
                shiny::downloadButton(ns("btn_download_seurat"), "Download Seurat", class = "btn-success w-100")
            }
        })

        # ---------------------------------------------------------------------
        # 3. METADATA AGGREGATOR HELPER
        # ---------------------------------------------------------------------
        compile_global_metadata <- function() {
            req(ScIGMA_data$mae)

            # 1. Base stricte
            cells <- colnames(ScIGMA_data$mae)[["dna_variants"]]
            print("cells")
            print(head(cells))

            # FIX : Création sans rownames initiaux pour éviter l'erreur de longueur
            meta_df <- data.frame(Cell_Barcode = cells, stringsAsFactors = FALSE)

            # 2. Clones ADN
            if (isTRUE(input$export_use_compass_clones) && !is.null(ScIGMA_data$dna.clones)) {
                # Extraction par mapping explicite pour éviter les décalages
                print("compile_global_metadata test 0")
                meta_df$DNA_Clone <- as.character(ScIGMA_data$dna.clones[cells])
                print("compile_global_metadata test 1")
            } else if (!is.null(ScIGMA_data$dna.clones_pre_compass)) {
                print("compile_global_metadata test 2")
                meta_df$DNA_Clone <- as.character(ScIGMA_data$dna.clones_pre_compass[cells])
                print("compile_global_metadata test 3")
            } else {
                meta_df$DNA_Clone <- "Unassigned"
            }
            meta_df$DNA_Clone[is.na(meta_df$DNA_Clone)] <- "Missing"

            # 3. Seurat Clusters
            if (!is.null(ScIGMA_data$seurat_object) && "seurat_clusters" %in% colnames(ScIGMA_data$seurat_object@meta.data)) {
                # Conversion du facteur en character avant injection
                cluster_vec <- as.character(ScIGMA_data$seurat_object$seurat_clusters)
                names(cluster_vec) <- rownames(ScIGMA_data$seurat_object@meta.data)

                meta_df$Protein_UMAP_Cluster <- cluster_vec[cells]
                meta_df$Protein_UMAP_Cluster[is.na(meta_df$Protein_UMAP_Cluster)] <- "Unassigned"
            } else {
                meta_df$Protein_UMAP_Cluster <- "Unassigned"
            }

            # 4. Gating Manuel (Bi-plot)
            meta_df$Biplot_Gate <- "Background"
            if (!is.null(ScIGMA_data$protein_gating_tree) && length(ScIGMA_data$protein_gating_tree$gates_list) > 0) {
                gates <- ScIGMA_data$protein_gating_tree$gates_list
                meta <- ScIGMA_data$protein_gating_tree$meta_list

                valid_ids <- setdiff(names(gates), "root")
                if (length(valid_ids) > 0) {
                    depths <- sapply(valid_ids, function(id) meta[[id]]$depth)
                    ordered_ids <- valid_ids[order(depths)]

                    # Extraction stricte des barcodes depuis les protéines
                    all_cell_barcodes <- colnames(ScIGMA_data$mae[["proteins"]])

                    for (sid in ordered_ids) {
                        cell_indices <- gates[[sid]]
                        cell_barcodes <- all_cell_barcodes[cell_indices]

                        # Intersection sécurisée pour ne cibler que les cellules valides
                        valid_barcodes <- intersect(cell_barcodes, meta_df$Cell_Barcode)
                        if(length(valid_barcodes) > 0){
                            meta_df$Biplot_Gate[match(valid_barcodes, meta_df$Cell_Barcode)] <- meta[[sid]]$name
                        }
                    }
                }
            }

            # rename Cell_Barcode.value column
            colnames(meta_df)[colnames(meta_df) == "Cell_Barcode.value"] <- "Cell_Barcode"

            return(meta_df)
        }

        # ---------------------------------------------------------------------
        # 4. DOWNLOAD HANDLERS
        # ---------------------------------------------------------------------

        # --- A. Export SEURAT (.rds) ---
        output$btn_download_seurat <- shiny::downloadHandler(
            filename = function() {
                paste0("ScIGMA_SeuratProject_", format(Sys.time(), "%Y%m%d_%H%M"), ".rds")
            },
            content = function(file) {
                shiny::req(ScIGMA_data$seurat_object)
                shiny::showNotification("Building final Seurat object... This may take a moment.", type = "message", duration = 5)

                export_seurat <- ScIGMA_data$seurat_object

                # Injection des Métadonnées Globales
                global_meta <- compile_global_metadata()
                aligned_meta <- global_meta[colnames(export_seurat), , drop = FALSE]
                export_seurat <- Seurat::AddMetaData(export_seurat, metadata = aligned_meta)

                # Injection conditionnelle de la matrice ADN
                if ("dna" %in% input$export_omics_selection && !is.null(ScIGMA_data$mae[["dna_variants"]])) {
                    use_compass <- isTRUE(input$export_use_compass_clones)
                    if (use_compass && "compass_imputed" %in% SummarizedExperiment::assayNames(ScIGMA_data$mae[["dna_variants"]])) {
                        gt_mat <- SummarizedExperiment::assay(ScIGMA_data$mae[["dna_variants"]], "compass_imputed")
                        assay_name <- "DNA_Imputed"
                    } else {
                        gt_mat <- SummarizedExperiment::assay(ScIGMA_data$mae[["dna_variants"]], "gt")
                        assay_name <- "DNA_Raw"
                    }

                    gt_mat_aligned <- as.matrix(gt_mat[, colnames(export_seurat), drop = FALSE])
                    gt_mat_aligned[gt_mat_aligned == 3L] <- NA

                    if (!is.null(ScIGMA_data$variants.filtered)) {
                        valid_vars <- intersect(rownames(ScIGMA_data$variants.filtered), rownames(gt_mat_aligned))
                        gt_mat_aligned <- gt_mat_aligned[valid_vars, , drop = FALSE]
                        rownames(gt_mat_aligned) <- ScIGMA_data$variants.filtered[valid_vars, "label"]
                    }

                    dna_assay <- Seurat::CreateAssayObject(counts = gt_mat_aligned)
                    export_seurat[[assay_name]] <- dna_assay
                }

                saveRDS(export_seurat, file = file)
            }
        )

        # --- B. Export MultiAssayExperiment (.rds) ---
        output$btn_download_mae <- shiny::downloadHandler(
            filename = function() {
                paste0("ScIGMA_MAE_", format(Sys.time(), "%Y%m%d_%H%M"), ".rds")
            },
            content = function(file) {
                shiny::req(ScIGMA_data$mae)
                shiny::showNotification("Saving MultiAssayExperiment...", type = "message")

                export_mae <- ScIGMA_data$mae
                global_meta <- compile_global_metadata()
                current_colData <- as.data.frame(SummarizedExperiment::colData(export_mae))
                current_colData$Cell_Barcode <- rownames(current_colData)

                merged_colData <- merge(current_colData, global_meta, by = "Cell_Barcode", all.x = TRUE)
                print("test4")
                # rownames(merged_colData) <- merged_colData$Cell_Barcode

                print("rownames(SummarizedExperiment::colData(export_mae) ")
                print(head(rownames(SummarizedExperiment::colData(export_mae) )))

                print("new rownames")
                print(rownames(S4Vectors::DataFrame(merged_colData[match(colnames(export_mae)[["dna_variants"]] ,
                                                                         merged_colData$Cell_Barcode), ]) |> head()))

                print("duplicated colnames(export_mae)[[dna_variants]]")
                print(sum(duplicated(colnames(export_mae)[["dna_variants"]])))
                print("merged_colData")
                print(head(merged_colData))


                print(all(rownames(SummarizedExperiment::colData(export_mae) ) == rownames(S4Vectors::DataFrame(merged_colData[match(colnames(export_mae)[["dna_variants"]] ,
                                                                                                                                     merged_colData$Cell_Barcode), ]) |> head())))

                SummarizedExperiment::colData(export_mae) <- S4Vectors::DataFrame(merged_colData[match(colnames(export_mae)[["dna_variants"]] ,
                                                                                                       merged_colData$Cell_Barcode), ])

                saveRDS(export_mae, file = file)
            }
        )

        # --- C. Export CSV Tables (.zip) ---
        output$btn_download_csv <- shiny::downloadHandler(
            filename = function() {
                paste0("ScIGMA_Export_Tables_", format(Sys.time(), "%Y%m%d_%H%M"), ".zip")
            },
            content = function(file) {
                shiny::req(ScIGMA_data$mae)
                shiny::showNotification("Generating CSV files and compressing archive...", type = "message", duration = 5)

                temp_dir <- file.path(tempdir(), paste0("scigma_export_", as.numeric(Sys.time())))
                dir.create(temp_dir)
                files_to_zip <- c()

                # 1. Cell Metadata (Généré si la case est cochée)
                if ("metadata" %in% input$export_omics_selection) {
                    meta_df <- compile_global_metadata()
                    meta_file <- file.path(temp_dir, "metadata_cells.csv")
                    write.csv(meta_df, meta_file, row.names = FALSE)
                    files_to_zip <- c(files_to_zip, meta_file)
                }

                # 2. DNA Variants Info
                if ("variants" %in% input$export_omics_selection && !is.null(ScIGMA_data$mae[["dna_variants"]])) {
                    var_info <- as.data.frame(SummarizedExperiment::rowData(ScIGMA_data$mae[["dna_variants"]]))
                    var_info$Variant_ID <- rownames(var_info)
                    var_file <- file.path(temp_dir, "dna_variants_annotation.csv")
                    write.csv(var_info, var_file, row.names = FALSE)
                    files_to_zip <- c(files_to_zip, var_file)
                }

                # 3. DNA Genotypes Matrix
                if ("dna" %in% input$export_omics_selection && !is.null(ScIGMA_data$mae[["dna_variants"]])) {
                    use_compass <- isTRUE(input$export_use_compass_clones)
                    if (use_compass && "compass_imputed" %in% SummarizedExperiment::assayNames(ScIGMA_data$mae[["dna_variants"]])) {
                        gt_mat <- as.matrix(SummarizedExperiment::assay(ScIGMA_data$mae[["dna_variants"]], "compass_imputed"))
                        gt_file <- file.path(temp_dir, "matrix_dna_genotypes_imputed.csv")
                    } else {
                        gt_mat <- as.matrix(SummarizedExperiment::assay(ScIGMA_data$mae[["dna_variants"]], "gt"))
                        gt_file <- file.path(temp_dir, "matrix_dna_genotypes_raw.csv")
                    }
                    write.csv(gt_mat, gt_file, row.names = TRUE)
                    files_to_zip <- c(files_to_zip, gt_file)
                }

                # 4. Proteins Matrix
                if ("proteins" %in% input$export_omics_selection && !is.null(ScIGMA_data$mae[["proteins"]])) {
                    assay_to_use <- ifelse("normalized" %in% SummarizedExperiment::assayNames(ScIGMA_data$mae[["proteins"]]), "normalized", "counts")
                    ptn_mat <- as.matrix(SummarizedExperiment::assay(ScIGMA_data$mae[["proteins"]], assay_to_use))
                    ptn_file <- file.path(temp_dir, paste0("matrix_proteins_", assay_to_use, ".csv"))
                    write.csv(ptn_mat, ptn_file, row.names = TRUE)
                    files_to_zip <- c(files_to_zip, ptn_file)
                }

                # 5. CNV Matrix
                if ("cnv" %in% input$export_omics_selection && !is.null(ScIGMA_data$mae[["amplicons"]])) {
                    cnv_mat <- as.matrix(SummarizedExperiment::assay(ScIGMA_data$mae[["amplicons"]], "counts"))
                    cnv_file <- file.path(temp_dir, "matrix_cnv_amplicons_counts.csv")
                    write.csv(cnv_mat, cnv_file, row.names = TRUE)
                    files_to_zip <- c(files_to_zip, cnv_file)
                }

                current_wd <- getwd()
                setwd(temp_dir)
                on.exit(setwd(current_wd))

                zip::zip(zipfile = file, files = basename(files_to_zip))
            }
        )

    })
}
## To be copied in the UI
# mod_download_panel_ui("download_panel_1")

## To be copied in the server
# mod_download_panel_server("download_panel_1", ScIGMA_data)
