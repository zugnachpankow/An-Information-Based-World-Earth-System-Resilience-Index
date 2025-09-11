# An-Information-Based-World-Earth-System-Resilience-Index

# Project Title
Code base for the forthcoming publication on Information-Based World-Earth System Resilience

[![License: GPL v3](https://img.shields.io/badge/License-GPLv3-blue.svg)](LICENSE)
[![DOI](https://img.shields.io/badge/DOI-INSERT_DOI-blue)]()
[![Julia 1.10.3](https://img.shields.io/badge/Julia-1.10.3-blueviolet)](https://julialang.org/downloads/)

---

## Table of Contents
1. [Background](#background)
2. [Features](#features)
3. [Installation](#installation)
4. [Usage](#usage)
5. [License](#license)
6. [Citation](#citation)
7. [Contributing](#contributing)
8. [Contact](#contact)

---

## Background
Navigating the challenges of the anthropogenic World-Earth system requires robust indicators to quantify resilience. This project introduces an information-based resilience metric, defined as the probability of reaching a desired system regime given initial conditions and available information. The metric accounts for knowledge about system dynamics, boundaries, and perturbations, as well as intrinsic action capacities. Our approach enables operationalization and quantification of resilience in complex models, highlighting the interplay between information and action. The index is broadly applicable, supporting analysis and decision-making.

Link to the preprint (forthcoming): [TBD]().

---

## Features
- Implements an information-based resilience index for a World-Earth system
- Extends the [Anderies et al. 2023](https://iopscience.iop.org/article/10.1088/1748-9326/ace91d) model
- Supports ensemble simulations to explore resilience under varying information and action capacities
- Provides tools for deterministic and stochastic scenario analysis
- Includes scripts for generating publication-quality figures from simulation outputs
- Designed for reproducibility on local machines and HPC clusters

---

## Installation

The project uses the [julia](https://julialang.org/) language and some of its packages: 

### Requirements
- Julia = `1.10.3` (other versions may work, but exact reproducibility is guaranteed only with this one)

### Dependencies
- ModelingToolkit: https://docs.sciml.ai/ModelingToolkit/stable/
- DifferentialEquations: https://docs.sciml.ai/DiffEqDocs/stable/
- Distributions: https://juliastats.org/Distributions.jl/stable/
- DataFrames: https://dataframes.juliadata.org/stable/
- Arrow: https://arrow.apache.org/julia/stable/
- Plots: https://docs.juliaplots.org/stable/
- JLD2: https://github.com/JuliaIO/JLD2.jl
- Colors: https://juliagraphics.github.io/Colors.jl/stable/
- CSV: https://csv.juliadata.org/stable/
- Statistics: https://docs.julialang.org/en/v1/stdlib/Statistics/
- Roots: https://juliamath.github.io/Roots.jl/stable/
- NLsolve: https://docs.sciml.ai/NonlinearSolve/stable/api/nlsolve/
- ColorSchemes: https://github.com/JuliaGraphics/ColorSchemes.jl

### Clone the repository
git clone https://github.com/zugnachpankow/An-Information-Based-World-Earth-System-Resilience-Index.git
cd An-Information-Based-World-Earth-System-Resilience-Index

### Reproduce the Julia environment

#### 1. Start Julia with this project environment
julia --project=.

#### 2. Inside the Julia REPL, run:
using Pkg

Pkg.instantiate()

## Usage
How to use this and reproduce results?

### How to obtain data
- Download the dataset from [Zenodo](https://doi.org/10.5281/zenodo.17098743) and place the files in this folder.
- Alternatively, you can generate the data by running the scripts in `study_scripts/ensemble_simulations/`. If run on a hpc, remember to adjust paths and batch scripts. 

### How to reproduce publication figures
- Either download the dataset or run scripts in study_scripts/ensemble_simulations locally or on a cluster.
- Then run the scripts in study_scripts/results, this will create plots in /figures and CSV files with overall resiliences in /results.

### How to do Your own analysis
If you want to extend the model or do your own work with it please get in touch, we can provide access to the working repo that includes tutorials.

## License
This project is licensed under the GNU General Public License v3.0 (GPLv3). See the LICENSE file for details.

## Citation
TODO!

## Main Contributors
- Max Bechthold 
- John M. Anderies

## Main Contact
- Max Bechthold <maxbecht@pik-potsdam.de>
