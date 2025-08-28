# Statistics-Research-Project
This repository contains the data and R code for the MSc Statistics research project, "Static vs. Dynamic Approaches to Financial Modeling: A Performance Analysis of Lévy and BNS Frameworks on the S\&P 500 Index." All figures and results presented in the thesis are fully reproducible using the code provided.

## Abstract

Classical Gaussian models fail to capture the heavy tails and volatility clustering characteristic of financial returns. This project empirically compares models for S&P 500 daily log-returns, first evaluating a range of static Lévy process models (e.g., Generalized Hyperbolic, CGMY) before advancing to the Barndorff-Nielsen and Shephard (BNS) stochastic volatility framework. Model parameters are obtained via Maximum Likelihood Estimation, employing particle filtering (PF) and the Auxiliary Particle Filter (APF) to handle the latent volatility of the BNS models. The results demonstrate the superiority of the dynamic BNS framework, which successfully replicates the persistent autocorrelation in absolute returns—a key feature of volatility clustering entirely absent in static counterparts.

## Getting Started

Follow these instructions to set up the project environment and reproduce the analysis.

### Requirements

* **R** (version >= 4.5.1)
* The following R packages:
    * `stats4`
    * `here`
    * `tidyverse`
    * `stabledist`
    * `fBasics`
    * `ghyp`
    * `pracma`
    * `goftest`
    * `KSgeneral`
    * `twosamples`
    * `gridExtra`
    * `SuppDists`

### Installation

1.  **Clone this repository** to your local machine:
    ```bash
    git clone https://github.com/YiruiCui/Statistics-Research-Project.git
    cd Statistics-Research-Project
    ```

2.  **Install the required R packages.** Open an R console and run the following command to install all dependencies:
    ```R
    install.packages(c("stats4", "here", "tidyverse", "stabledist", "fBasics", "ghyp", "pracma", "goftest", "KSgeneral", "twosamples", "gridExtra", "SuppDists"))
    ```

---

## How to Reproduce the Analysis

All scripts required to run the analysis are located in the `analysis` directory. Before running a script, please check for any notes in the first few lines regarding dependencies or execution order.

All generated figures will be saved to the `outputs` directory. 

All other results (e.g., parameter estimates, p-values from GoF tests) will be printed directly to the console.

---

## Project Structure

* **`data`**: Contains the raw S&P 500 daily closing prices used in the analysis.
* **`src`**: Contains custom R functions for simulation (BNS-Gamma & BNS-IG), and evaluation (critical points for AD test) that are called by the main analysis scripts.
* **`analysis`**: Contains the primary R scripts used to perform the statistical analysis, fit the models, and generate all plots presented in the thesis.
* **`outputs`**: The default directory where all generated figures are saved.
