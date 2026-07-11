# 1. Figer l'OS et la version de R/Bioconductor
FROM bioconductor/bioconductor_docker:RELEASE_3_23 AS builder

# 2. Figer les dépendances système C++ (COMPASS, HDF5, etc.)
RUN apt-get update && apt-get install -y \
    libhdf5-dev \
    libcurl4-openssl-dev \
    libxml2-dev \
    libssl-dev \
    librsvg2-dev \
    build-essential \
    && rm -rf /var/lib/apt/lists/*

# 3. Copier l'environnement de reproductibilité dans le conteneur
WORKDIR /reproducibility
COPY data/ data/
COPY renv.lock .
# COPY renv/ renv/ # Uncomment if renv folder is tracked
# COPY .Rprofile . # Uncomment if local .Rprofile is used
COPY _targets.R .
COPY config.yml .
COPY ScIGMA_0.99.3.tar.gz ScIGMA_0.99.3.tar.gz
# 4. Restauration bit-à-bit (avec patch dynamique pour Bioconductor Devel)

RUN Rscript -e "install.packages(c('jsonlite', 'BiocManager', 'renv'), repos='https://cloud.r-project.org'); \
    install.packages(c('targets', 'DiagrammeRsvg', 'rsvg', 'htmlwidgets')); \
    lock <- jsonlite::read_json('renv.lock'); \
    bioc_pkgs <- names(Filter(function(x) identical(x[['Source']], 'Bioconductor'), lock[['Packages']])); \
    bioc_pkgs <- c(bioc_pkgs, 'Seqinfo'); \
    options(renv.config.bioconductor.version = '3.21'); \
    BiocManager::install(bioc_pkgs, ask=FALSE, update=FALSE); \
    renv::record(paste0('bioc::', bioc_pkgs)); \
    renv::restore(prompt = FALSE); \
    if (!requireNamespace('remotes', quietly = TRUE)) install.packages('remotes'); \
    remotes::install_local('ScIGMA_0.99.3.tar.gz', dependencies = TRUE)"

COPY R/ R/

# 5. Exécution automatique de la science à la création de l'image (pour générer les figures)
RUN Rscript -e "targets::tar_make()"

# 6. Étape de copies des figures générées dans le dossier reproductibility
# Utiliser `docker build -o type=local,dest=. .` pour extraire les figures
FROM scratch AS export
COPY --from=builder /reproducibility/results/figures /figures
