# 📌 Overview

This repository contains code, methods, and visualisations for a mini-project in **Uncertainty Quantification IV**, completed as part of the MATH module at Durham University. The project focuses on developing and applying **Bayes Linear emulators** to understand and approximate computational models with uncertain inputs and outputs.

The analysis follows a structured approach based on the course practicals, progressing from simple one-dimensional emulation to efficient multi-dimensional implementations, including **multi-wave history matching**, **design diagnostics**, and **prediction interval analysis**.

# 🎯 Objectives

- Explore the impact of uncertain inputs on physical models  
- Construct and test Bayes Linear emulators (1D and 2D)  
- Evaluate emulator performance through diagnostics and plots  
- Apply multi-wave history matching to refine input space  
- Compare emulator output to true function values for validation

# 📁 Repository Structure

| File               | Purpose |
|--------------------|---------|
| `Practicals`         | 1. Intro to 1D emulation with full BL formulation and prediction intervals.
                         2. Extension to 2D emulation using grid design and contour visualisation
                         3. History matching implementation and implausibility function diagnostics
                         4. Efficient emulation via matrix algebra and `pdist` for large-scale prediction |
| `New_Code.R`       | Final cleaned and integrated pipeline with all methods in one file |
| `Raul_Unnithan_Report.pdf` | Final mini-project report, with commentary and plots |
| `Mini-Project Report Guidelines.pdf` | Project instructions and marking criteria from the module |

# 🧠 Physical Model

The model used is the **rescaled Branin-Hoo function**, selected from the [Virtual Library of Simulation Experiments](https://www.sfu.ca/~ssurjano/index.html). It provides a realistic testbed for emulation due to its:

- Nonlinear behaviour
- Multiple global minima
- Interacting input dimensions

Inputs are scaled to the unit square \([0, 1]^2\) and evaluated using Latin Hypercube Sampling (LHS) for effective coverage.

# 🔬 Methodology

The main techniques demonstrated are:

- **Bayes Linear Emulation (BLE)**  
  - Construction of expectation and variance estimators  
  - Closed-form BL adjustments using covariance matrices  
  - Emulator training on sparse design points

- **1D and 2D Emulation**  
  - Fixing inputs and exploring response in 1D  
  - Full-grid contouring in 2D with uncertainty maps

- **Efficient Emulator Implementation**  
  - Optimised matrix inversion using Cholesky decomposition  
  - Pairwise distances computed via `pdist` for speed

- **Multi-Wave History Matching**  
  - Definition of observational discrepancy and error  
  - Implausibility scoring to exclude poor-fitting regions  
  - Sequential wave analysis to improve emulator focus

# 📦 Requirements

The following R packages are used across the scripts:

```r
install.packages(c("plot3D", "lhs", "pdist", "sensitivity", "viridisLite"))
