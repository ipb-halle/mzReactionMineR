# R package for downstream analysis of mass spectrometry data after peak picking with external software.

mzReactionMineR is an R package for the analysis of metabolomics data that was pre-processed using external software outside of R (i.e. mzmine, MetaboScape). As its base it utilizes the default SummarizedExperiment object to allow for seamless integration to other R packages. Furthermore, it keeps the relatioship between sample and feature meta data with the measured data stable, minimizing the risk of processing errors.

# Installation

mzReactionMineR is hosted on GitHub and can be installed via devtools::install_github() or pak::pkg_install()

```
pak::pkg_install("ipb-halle/mzReactionMineR")
# devtools::install_github("ipb-halle/mzReactionMineR")
```

Alternatively, you can clone this repository and install from local

```
cd /path/to/local/clone
git clone git@github.com:ipb-halle/mzReactionMineR.git
install.packages("/path/to/local/clone/mzReactionMineR/", repos = NULL, type="source")

# or
# pak::local_install("/path/to/local/clone/mzReactionMineR/")
# devtools::install_local("/path/to/local/clone/mzReactionMineR/")

```

# Functionality overview

The package consist of many different functionalities, ranging from reading to manipulating metabolomics data:

## Reading data

- **mzmine_to_se**: Read output from mzmine and create a SummarizedExperiment.
- **feature_table_to_se**: Read generic feature tables and create SummarizeExperiment (i.e. MetaboScape export for SIRIUS).
- **join_se_sirius**: join output from SIRIUS to a SummarizedExperiment generated using mzReactionMineR.
- **join_align**: align multiple metabolomics experiments into a single SumamrizedExperiment.

## Filtering

- **filter_se**: Filter a metabolomics dataset stored in a SummarizedExperiment object by various measures.
- **blank_subtraction_se**: Remove features that were detected in blank samplesz.
- **remove_fature_mz**: Removes features based on specific m/z values within tolerance.
- **filter_spec**: filters a MS/MS data stored as a Spectra object.

## Imputation and normalization

- **impute_min_frac**: Impute mising values based on a fraction of the lowest value per feature.
- **normalize_pqn**: Use probabilistic quotient normalization to normalize an assay.
- **normalize_is**: Normalization based on one internal standard.

## Statistics

- **anova_limma**: execute an anova test using the limma package.
- **contrast_limma**: execute a pairwise contrast test using the limma package.
- **knn_clustering_samples**: cluster metabolomics samples based on community detection on a *k*-nn graph.
- **louvain_clustering_features**: cluster metabolomics features based on community detection on a correlation graph.

## Visualization

- **qc_plots**: generate qc plots (i.e. per sample m/z or rt deviation).
- **plot_is**: plot internal standard against the run index.


# Funding

Funded by the **German Federal Ministry of Education and Research** and the state of **Saxony-Anhalt** as part of the project **DiP-NA-WIR**.