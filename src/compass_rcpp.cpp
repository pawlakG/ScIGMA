#include <Rcpp.h>
#include <exception>

#ifdef _OPENMP
#include <omp.h>
#endif

#include "Inference.h"
#include "Tree.h"
#include "Scores.h"
// Supprimer l'inclusion de "input.h" car nous l'avons supprimé

// Déclaration des variables globales requises par COMPASS
int n_cells;
int n_loci;
int n_regions;
std::vector<Cell> cells;
Data data;
Params parameters;

// NEW : Fonction d'ingestion In-Memory
// Remplace le vieux load_CSV de l'ETHZ
void ingest_r_data(
        Rcpp::IntegerMatrix ref_counts,
        Rcpp::IntegerMatrix alt_counts,
        Rcpp::IntegerMatrix genotypes,
        Rcpp::IntegerVector locus_region_mapping,
        Rcpp::IntegerMatrix region_counts, // Matrice vide si use_cna = false
        bool use_cna
) {
    n_loci = ref_counts.nrow();
    n_cells = ref_counts.ncol();

    // Validation dimensionnelle stricte
    if (alt_counts.nrow() != n_loci || alt_counts.ncol() != n_cells) {
        throw std::runtime_error("Dimensions mismatch between REF and ALT matrices.");
    }

    if (use_cna) {
        n_regions = region_counts.nrow();
        if (region_counts.ncol() != n_cells) {
            throw std::runtime_error("Dimensions mismatch between variant matrices and region matrices.");
        }
    } else {
        n_regions = 0;
    }

    // Remplissage de la structure Cell
    cells.clear();
    cells.reserve(n_cells);

    for (int j = 0; j < n_cells; j++) {
        Cell c;
        c.ref_counts.reserve(n_loci);
        c.alt_counts.reserve(n_loci);
        c.genotypes.reserve(n_loci);

        for (int i = 0; i < n_loci; i++) {
            c.ref_counts.push_back(ref_counts(i, j));
            c.alt_counts.push_back(alt_counts(i, j));
            c.genotypes.push_back(genotypes(i, j));
        }

        c.name = "Cell_" + std::to_string(j);

        int total_count = 0;
        if (use_cna) {
            c.region_counts.reserve(n_regions);
            for (int k = 0; k < n_regions; k++) {
                int count = region_counts(k, j);
                c.region_counts.push_back(count);
                total_count += count;
            }
            c.total_counts = total_count;
        }
        cells.push_back(c);
    }

    // Remplissage minimal de la structure Data
    // Dans une version complète, on passerait aussi les vecteurs de chr/pos depuis R
    data.locus_to_region.clear();
    for(int i = 0; i < n_loci; i++) {
        data.locus_to_region.push_back(locus_region_mapping[i]);
    }

    data.region_to_loci.clear();
    data.region_to_loci.resize(n_regions);
    for (int i = 0; i < n_loci; i++){
        data.region_to_loci[data.locus_to_region[i]].push_back(i);
    }

    // Initialisation des autres paramètres de Data (simplifié)
    data.variant_is_SNV = std::vector<bool>(n_loci, true);
    data.region_is_reliable = std::vector<bool>(n_regions, true);
}


// [[Rcpp::export]]
int run_compass_inference_cpp(
        Rcpp::IntegerMatrix ref_counts,
        Rcpp::IntegerMatrix alt_counts,
        Rcpp::IntegerMatrix genotypes,
        Rcpp::IntegerVector locus_region_mapping,
        Rcpp::IntegerMatrix region_counts,
        std::string output_prefix,
        int n_chains = 4,
        int chain_length = 5000,
        bool use_cna = true,
        std::string sex = "female"
) {
    try {
        // ... (init_params(), data.sex = sex, burn_in, etc.) ...

        // NEW : Appel de l'ingestion In-Memory au lieu de load_CSV
        ingest_r_data(
            ref_counts,
            alt_counts,
            genotypes,
            locus_region_mapping,
            region_counts,
            use_cna
        );

        // ... (Configuration OpenMP et boucle MCMC inchangée) ...

    } catch ( const std::exception& e ) {
        Rcpp::stop(e.what());
    } catch ( ... ) {
        Rcpp::stop("Fatal unknown C++ error during COMPASS execution.");
    }

    return 0;
}
