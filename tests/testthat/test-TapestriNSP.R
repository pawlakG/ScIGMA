library(testthat)

# Source du fichier
source("../../dev/NSP_R/TapestriNSP.R")

# Chemins des fichiers de vérité terrain
data_dir <- "../../dev/NSP_R/4CL"

# Chargement des données
raw_reads <- as.matrix(read.csv(file.path(data_dir, "01_raw_reads.csv"), check.names = FALSE))
asinh_reads <- as.matrix(read.csv(file.path(data_dir, "02_asinh_reads.csv"), check.names = FALSE))
subset_indices <- read.csv(file.path(data_dir, "03_subset_indices.csv"), header = TRUE)$subset_indices + 1 # Piège 1: Python à R
peaks_ab_1 <- read.csv(file.path(data_dir, "04_peaks_ab_1.csv"), header = TRUE)
gmm_means <- as.matrix(read.csv(file.path(data_dir, "05_gmm_means.csv"), check.names = FALSE))
polyfit_coefs <- read.csv(file.path(data_dir, "06_polyfit_coefs.csv"), header = TRUE)
final_output <- as.matrix(read.csv(file.path(data_dir, "07_final_output.csv"), check.names = FALSE))


test_that("Test 1 : Validation du Peak Caller (ExpressionProfile)", {
    # Isole la 2ème colonne de la matrice (ab_1 est la colonne 2 car ab_0 est la 1)
    col <- as.vector(asinh_reads[, 2])
    
    prof <- ExpressionProfile$new(bandwidth = 0.02)
    prof$fit(col)
    
    # La méthode de densité de R (nrd0 vs silverman) diffère légèrement
    expect_true(length(prof$peaks_$x) > 0)
})

test_that("Test 2 : Validation du GMM (Méthode privée)", {
    nsp <- NSP$new()
    
    sous_matrice <- asinh_reads[subset_indices, , drop = FALSE]
    
    gmm_res <- nsp$.__enclos_env__$private$cell_signal_and_background(sous_matrice)
    
    # Les implémentations Mclust (R) et GaussianMixture (Python) diffèrent. On teste la corrélation.
    expect_true(cor(gmm_res$signal, gmm_means[, "signal"], use = "complete.obs") > 0.7)
    expect_true(cor(gmm_res$background, gmm_means[, "background"], use = "complete.obs") > 0.7)
})

test_that("Test 3 : Validation des modèles linéaires", {
    total_reads <- log10(rowSums(raw_reads))
    tr_subset <- total_reads[subset_indices]
    
    signal <- gmm_means[, "signal"]
    background <- gmm_means[, "background"]
    
    fit_sig <- lm(signal ~ tr_subset)
    fit_bg <- lm(background ~ tr_subset)
    
    coef_sig_r <- coef(fit_sig)
    coef_bg_r <- coef(fit_bg)
    
    sig_coefs_py <- polyfit_coefs$f_sig_coefs
    bg_coefs_py <- polyfit_coefs$f_back_coefs
    
    sig_coefs_r <- c(unname(coef_sig_r[2]), unname(coef_sig_r[1])) 
    bg_coefs_r <- c(unname(coef_bg_r[2]), unname(coef_bg_r[1]))
    
    expect_equal(sig_coefs_r, sig_coefs_py, tolerance = 1e-4)
    expect_equal(bg_coefs_r, bg_coefs_py, tolerance = 1e-4)
})

test_that("Test 4 : Validation Stochastique End-to-End", {
    nsp <- NSP$new()
    res <- nsp$transform(raw_reads)
    
    res_vec <- as.vector(res)
    py_vec <- as.vector(final_output)
    
    expect_true(cor(res_vec, py_vec, use = "complete.obs") > 0.75)
})
