# Standardized Generalized Hierarchical Factor Models (GenHFM)

This repository contains all the data, R code, Stan models, and supplementary materials associated with the paper: **"A Unified Framework for Psychometrics in Experimental Psychology: The Standardized Generalized Hierarchical Factor Model"**.

> **Note:** Due to the large size of certain output files (fitted models and large `.RData` objects), the heaviest results are hosted on our **OSF project [https://osf.io/5dutv/]**.

---

## 📂 Repository Structure

The directory is organized as follows:

### 1. Data
Contains the datasets used across the analyses.
* **1.1 Processed**: Cleaned data ready for model fitting.
* **1.2 Raw data**: Original raw datasets.

### 2. Figures
Contains all figures presented in the main text of the manuscript.

### 3. Manuscript
Drafts of the manuscript and Supplementary Material documents.

### 4. Model diagrams
Images and schematics illustrating the GenHFM architecture.

### 5. R scripts
The complete R codebase for reproducing the entire study.
* **5.1 Appendix**: Code for generating the unit vector visualization.
* **5.2 Empirical analysis**:
    * **5.2.1 Rey-Mermet**: Code to fit GenHFMs using criteria from Rey-Mermet et al. (2018) and posterior predictive check (PPC) figures.
        * *datasets*: Original data from Rey-Mermet et al.
        * *DDM_fitted_models*: Pre-fitted DDM models (note: these may be hosted on OSF due to size constraints).
    * **5.2.2 Viviani**: Scripts for data cleaning, model fitting, and result analysis.
    * **5.2.3 Whitehead**: Scripts for data cleaning, model fitting, and result analysis.
* **5.3 Meta-analysis**: Code for fitting ANOVA metamodels and estimating parameter recovery indices.
* **5.4 Model Diagrams**: R code used to generate the model schematics and export individual images.
* **5.5 R functions**: Custom functions used throughout the empirical analyses, meta-analysis, and simulations.
* **5.6 Simulation study**: Comprehensive code for running the simulation, generating plots, and analyzing recovery.
* **5.7 Supplementary Material**: Code used to generate results specifically for the supplement.

### 6. Results
Summary outputs and saved R objects.
* **6.1 Rdata**:
    * **6.1.1 ELPDs**: Results for model comparison using Expected Log-Pointwise Predictive Density.
    * **6.1.2 LMM**: Fitted Linear Mixed Models for the Viviani dataset.
    * **6.1.3 Meta-parameters**: Results and estimates from the meta-analysis/ANOVA.
    * **6.1.4 Simulation study**: Raw and processed results from the simulation study.
    * **6.1.5 Supplementary Material**: Result objects for the supplementary material.
    * **6.1.6 Tables**: Results tables included in the main manuscript.

### 7. Stan models
The source code for the Bayesian models.
* **7.1 BGHFM**: The primary Generalized Hierarchical Factor Model (GenHFM) implementations.
* **7.2 BGHFM (sim study)**: Versions of the GenHFM optimized specifically for the simulation study.
* **7.3 BLMM**: Bayesian Linear Mixed Models used in the meta-analysis.

---

## 🚀 How to Reproduce

To run the scripts in this repository, you will need **R**, **Stan**, and several R packages (primarily `cmdstanr`, `posterior` and `tidyverse`). 

To ensure all scripts run correctly without the need to manually set working directories, simply:
1. Download or clone the entire repository.
2. Open the `GenHFM.Rproj` file in RStudio.
3. Run the desired scripts; all file paths are relative to the project root.

## 🔗 External Links
* **OSF Repository:** [https://osf.io/5dutv/] (Contains large `.rds` model objects and high-volume datasets).

---
