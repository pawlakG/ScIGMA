# Worklow of ScIGMA application

## Input
h5 file from Tapestry pipeline
DNA only or DNA + Proteins

## Preprocessing
- Variant filtering
- Cells filtrering
- Protein normalization

## DNA
- Variant selection
- COMPASS imputation
- DNA variants heatmap
- DNA variants clones renaming

## CNV
- Amplicons filtering
- CNV Heatmap or lineplot figure by :
    - genes
    - chromosomes

## Proteins
- Descriptive: 
    - Ridge plot
    - Barplot
- Biplot gating
- UMAP modelisation:
    - protein expression projection
    - biplot gates projection
    - unsupervised clustering

## Muli-omics
- DNA x UMAP's protein model
- DNA x UMAP's protein clusters
- Biplot gates' DNA variants status

## Outputs
- Ready-to-publish figures
- Seurat object export (RDS)
- MultiAssayExperiment export (RDS)
- Flat tables (ZIP)
