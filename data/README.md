# data

This folder contains the data files required for running analyses and reproducing results in this study.

## How to obtain data
- Download the dataset from [Zenodo](https://doi.org/10.5281/zenodo.17098743) and place the files in this folder.
- Alternatively, you can generate the data by running the scripts in `study_scripts/ensemble_simulations/`.

## Data formats
- Data files are in `.arrow` format. They contain post-processes resilience calculations, i.e. for the deterministic value each parameter combination is attributed 0 or 1 and for the stochastic version a mean resilience value for the Monte Carlo simulation between 0 and 1.

## Usage
- Ensure the required data files are present before running analysis or plotting scripts.