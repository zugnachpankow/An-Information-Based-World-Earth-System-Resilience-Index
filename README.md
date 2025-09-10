# An-Information-Based-World-Earth-System-Resilience-Index

# Project Title
Code base for the forthcoming publication on Information-Based World-Earth System Resilience

[![License]()]()
[![DOI](https://img.shields.io/badge/DOI-INSERT_DOI-blue)]()

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
Brief context for the project. Explain the motivation, the problem addressed, and the main contribution.  
Link to the paper or preprint if available: [DOI or URL](INSERT_LINK).

---

## Features
- Key capability or feature #1
- Key capability or feature #2
- Key capability or feature #3

---

## Installation

The project uses the [julia](https://julialang.org/) language and some of its packages: 

### Requirements
- Julia = 

### Dependencies
(* optional)
- https://docs.sciml.ai/ModelingToolkit/stable/
- https://docs.sciml.ai/DiffEqDocs/stable/
- https://juliastats.org/Distributions.jl/stable/
- https://dataframes.juliadata.org/stable/
- https://arrow.apache.org/julia/stable/
- https://docs.juliaplots.org/stable/ (* for plotting)
- https://github.com/JuliaIO/JLD2.jl (* for saving)

### Clone the repository
git clone https://github.com/username/repo.git
cd repo

## Usage
How to use this and reproduce results?

### How to obtain data
- Download the dataset from [Zenodo](https://doi.org/10.5281/zenodo.17091722) and place the files in this folder.
- Alternatively, you can generate the data by running the scripts in `study_scripts/ensemble_simulations/`.

### Folder structure

### Locally:
1. add a folder in local_studies
2. use the provided example studies to build what you need
3. run locally

### On the cluster:
1. add a folder in cluster_studies
2. add a subfolder plots
3. add a study and a plot script (use the templates)
4. try to use explanatory naming and give some minimal explanation what your study does 
(do not add a docstring in the beginning of the file if you work on the cluster, it will misinterpret it)
5. use the submit_study.sh batch script to send it to the cluster

## License

## Citation

## Contributing

## Contact
