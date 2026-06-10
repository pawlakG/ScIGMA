// File: src/compass_rcpp.cpp
#include <Rcpp.h>
#include <exception>

#ifdef _OPENMP
#include <omp.h>
#endif

#include "Inference.h"
#include "Tree.h"
#include "Scores.h"

// Declaration of global variables required by COMPASS
int n_cells;
int n_loci;
int n_regions;
std::vector<Cell> cells;
Data data;
Params parameters;

// NEW: Configuration of hyperparameters
void init_params() {
    parameters.sequencing_error_rate = 0.02;
    parameters.omega_hom = 50.0;
    parameters.omega_het = 8.0;
    parameters.sequencing_error_rate_indel = 0.06;
    parameters.omega_hom_indel = 15.0;
    parameters.omega_het_indel = 4.0;

    parameters.prior_dropoutrate_mean = 0.05;
    parameters.prior_dropoutrate_omega = 100.0;

    parameters.theta = 6.0;
    parameters.doublet_rate = 0.08;

    parameters.use_doublets = true;
    parameters.filter_regions = true;
    parameters.filter_regions_CNLOH = true;
    parameters.verbose = true;

    parameters.node_cost = 1.0;
    parameters.CNA_cost = 85.0;
    parameters.LOH_cost = 85.0;
    parameters.mut_notAtRoot_cost = 10.0;
    parameters.mut_notAtRoot_freq_cost = 100000.0;
}

void ingest_r_data(
        Rcpp::IntegerMatrix ref_counts,
        Rcpp::IntegerMatrix alt_counts,
        Rcpp::IntegerMatrix genotypes,
        Rcpp::IntegerVector locus_region_mapping,
        Rcpp::IntegerMatrix region_counts,
        Rcpp::StringVector locus_names, // NEW
        bool use_cna
) {
    n_loci = ref_counts.nrow();
    n_cells = ref_counts.ncol();

    if ( alt_counts.nrow() != n_loci || alt_counts.ncol() != n_cells ) {
        throw std::runtime_error("Dimensions mismatch between REF and ALT matrices.");
    }

    int max_region = -1;
    for ( int i = 0; i < n_loci; i++ ) {
        if ( locus_region_mapping[i] > max_region ) {
            max_region = locus_region_mapping[i];
        }
    }

    if ( use_cna ) {
        n_regions = region_counts.nrow();
        if ( n_regions <= max_region ) {
            throw std::runtime_error("Region matrix has fewer rows than required by locus mapping.");
        }
        if ( region_counts.ncol() != n_cells ) {
            throw std::runtime_error("Dimensions mismatch between variants and regions.");
        }
    } else {
        n_regions = max_region + 1;
    }

    cells.clear();
    cells.reserve(n_cells);

    for ( int j = 0; j < n_cells; j++ ) {
        Cell c;
        c.ref_counts.reserve(n_loci);
        c.alt_counts.reserve(n_loci);
        c.genotypes.reserve(n_loci);

        for ( int i = 0; i < n_loci; i++ ) {
            c.ref_counts.push_back(ref_counts(i, j));
            c.alt_counts.push_back(alt_counts(i, j));
            c.genotypes.push_back(genotypes(i, j));
        }

        c.name = "Cell_" + std::to_string(j);

        if ( use_cna ) {
            c.region_counts.reserve(n_regions);
            int total_count = 0;
            for ( int k = 0; k < n_regions; k++ ) {
                int count = region_counts(k, j);
                c.region_counts.push_back(count);
                total_count += count;
            }
            c.total_counts = total_count;
        } else {
            c.region_counts = std::vector<int>(n_regions, 0);
            c.total_counts = 0;
        }
        cells.push_back(c);
    }

    data.locus_to_region.clear();
    for ( int i = 0; i < n_loci; i++ ) {
        data.locus_to_region.push_back(locus_region_mapping[i]);
    }

    data.region_to_loci.clear();
    data.region_to_loci.resize(n_regions);
    for ( int i = 0; i < n_loci; i++ ) {
        data.region_to_loci[data.locus_to_region[i]].push_back(i);
    }

    data.locus_to_name.clear();
    data.locus_to_chromosome.clear();
    data.locus_to_freq.clear();
    for ( int i = 0; i < n_loci; i++ ) {
        // UPDATED : Ingestion directe du nom de variant R vers le std::string C++
        data.locus_to_name.push_back(Rcpp::as<std::string>(locus_names[i]));
        data.locus_to_chromosome.push_back("1");
        data.locus_to_freq.push_back(0.0);
    }

    data.region_to_name.clear();
    data.region_to_chromosome.clear();
    for ( int k = 0; k < n_regions; k++ ) {
        data.region_to_name.push_back("Region_" + std::to_string(k));
        data.region_to_chromosome.push_back("1");
    }

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
        Rcpp::StringVector locus_names,
        Rcpp::StringVector locus_chromosomes,
        Rcpp::StringVector region_names,
        Rcpp::StringVector region_chromosomes,
        Rcpp::StringVector cell_names,          // <-- NEW : Les barcodes des cellules
        std::string output_prefix,
        int n_chains = 4,
        int chain_length = 5000,
        bool use_cna = true,
        std::string sex = "female"
) {
    try {
        init_params();
        parameters.verbose = false;
        data.sex = sex;
        int burn_in = chain_length / 2;

        ingest_r_data(
            ref_counts, alt_counts, genotypes,
            locus_region_mapping, region_counts,
            locus_names, use_cna
        );

        data.locus_to_name = Rcpp::as<std::vector<std::string>>(locus_names);
        data.locus_to_chromosome = Rcpp::as<std::vector<std::string>>(locus_chromosomes);
        data.region_to_name = Rcpp::as<std::vector<std::string>>(region_names);
        data.region_to_chromosome = Rcpp::as<std::vector<std::string>>(region_chromosomes);

        for (int i = 0; i < n_cells; i++) {
            cells[i].name = Rcpp::as<std::string>(cell_names[i]);
        }

        double betabin_overdisp = parameters.omega_het;

        parameters.omega_het = std::min(parameters.omega_het, betabin_overdisp);
        parameters.omega_het_indel = std::min(
            parameters.omega_het_indel,
            betabin_overdisp
        );
        parameters.omega_het = std::min(parameters.omega_het, betabin_overdisp);
        parameters.omega_het_indel = std::min(
            parameters.omega_het_indel,
            betabin_overdisp
        );

        std::vector<double> results{};
        results.resize(n_chains);
        std::vector<Tree> best_trees{};
        best_trees.resize(n_chains);

#ifdef _OPENMP
        if ( n_chains < omp_get_num_procs() ) {
            omp_set_num_threads(n_chains);
        } else {
            omp_set_num_threads(omp_get_num_procs());
        }
#endif

        Rcpp::Rcout << "Starting " << n_chains << " MCMC chains" << std::endl;

#ifdef _OPENMP
#pragma omp parallel for
#endif
        for ( int i = 0; i < n_chains; i++ ) {
            //std::srand(i);
            Inference infer{"", 10.0, i};
            best_trees[i] = infer.find_best_tree(use_cna, chain_length, burn_in);
            results[i] = best_trees[i].log_score;
        }

        double best_score = -DBL_MAX;
        int best_score_index = -1;
        for ( int i = 0; i < n_chains; i++ ) {
            if ( best_score < results[i] ) {
                best_score = results[i];
                best_score_index = i;
            }
        }

        best_trees[best_score_index].to_dot(output_prefix, true);
        Rcpp::Rcout << "Completed! Output written to " << output_prefix << std::endl;

    } catch ( const std::exception& e ) {
        Rcpp::stop(e.what());
    } catch ( ... ) {
        Rcpp::stop("Fatal unknown C++ error during COMPASS execution.");
    }

    return 0;
}
