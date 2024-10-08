# WMIDM
Data and code for paper "Study on Individual Differences in Visual Working Memory Tasks Based on Spatiotemporal Brain Functional Metrics and Biological Perspectives"

Bootstrap_gene.m:

This script performs a bootstrap analysis on gene expression data, resampling the data multiple times to evaluate the stability of statistical estimates for genes in a biological or neuroscience context, and finally calculate the PLS score for each gene.

Cell_types_permutationtest.m:

This script conducts permutation tests to evaluate the statistical significance of differences between different cell types based on gene expression data. The tests involve random shuffling (permuting) the data to compare actual results to a null distribution.

FractalDim.m:

This script is designed to calculate the fractal dimension of time series data or other signals, which is a measure used to describe the complexity or self-similarity of a structure or process, possibly in the context of brain signals or biological data.

HurstExponent.m:

This script computes the Hurst exponent, a measure used to evaluate the long-term memory of time series data. It indicates whether a time series is trending, mean-reverting, or random, and is commonly applied in neuroscience or financial data analysis.

PoincareIndex.m:
This script compute the Poincaré index, which is used in non-linear dynamic analysis to quantify the variability of data, often in physiological signals like heart rate variability or brain wave data.

Swaveletentropy.m:

This script calculates the wavelet entropy of a signal, which is a measure of the complexity and irregularity of a signal in the time-frequency domain, using wavelet transform techniques. This is commonly applied to analyze biological signals, such as brain or heart activity.

analyzeGeneCellTypes.m:

This script analyzes the relationship between gene expression and different cell types, by performing statistical tests to determine which genes are most associated with certain cell types in a given dataset.

trainLinearRegressionModel.m:

This script trains a linear regression model using a provided dataset, applying machine learning techniques to predict an outcome (response) based on a set of input features (predictors). It includes steps for training, validation, and testing of the model.

KEGGGO_Enrichment_CellType.R

This is a script for performing KEGG and GO gene enrichment analysis across different cell types. These categories include microglia (Micro), endothelial cells (Endo), oligodendrocyte precursor cells (OPC), oligodendrocytes (Oligo), astrocytes (Astro), excitatory neurons (Neuro-Ex), and inhibitory neurons (Neuro-In).

KEGGGO_Enrichment.R
This is a script for performing KEGG and GO gene enrichment analysis.
