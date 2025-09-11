# study_scripts

This folder contains Julia scripts used for ensemble runs and results/plotting of figures in the publication. 

## Contents

### ensemble_simulations/

This folder contains Julia scripts used for large ensemble simulations in the study. These scripts were designed to be run on the PIK HPC cluster to generate the main results for the publication.

- **deterministic_many_worlds_compare.jl**  
  Compares deterministic simulations across different model configurations.

- **deterministic_many_worlds_gridded.jl**  
  Runs deterministic simulations over a gridded parameter space.

- **deterministic_many_worlds_threshold_distribution.jl**  
  Explores threshold distributions in deterministic many-worlds scenarios.

- **stochastic_many_worlds_mult_large_sigma.jl**  
  Performs stochastic simulations with multiple worlds and large noise (sigma) values.

- **stochastic_many_worlds_mult_small_sigma.jl**  
  Performs stochastic simulations with multiple worlds and small noise (sigma) values.

- **stochastic_many_worlds_shocks.jl**  
  Simulates stochastic many-worlds scenarios with shock events.

- **submit_study.sh**  
  Example SLURM submission script for running ensemble simulations on the PIK HPC cluster.  
  **Note:** This script must be adjusted to your specific needs and should be run for each Julia script individually. You may also modify it to submit all jobs at once.

### results/

Scripts for generating publication-quality figures from simulation outputs (i.e. the ones used in "Results" section), get saved in ../figures/ folder:

- **figure3.jl**  
  Generates Figure 3 for the manuscript.

- **figure4.jl**  
  Generates Figure 4 for the manuscript. Generates a CSV with the IBRI values.

- **figure5.jl**  
  Generates Figure 5 for the manuscript. 

- **figure6.jl**  
  Generates Figure 6 for the manuscript. Generates a CSV with the IBRI values.

### trajectories_for_figure5/

Scripts for creating and saving the trajectories for figure5. It is highly recommended to just download the provided files (see Usage for download location). You will need to start a new REPL for each script as they will have conflicting namespaces - this is a problem with similar naming between the different models which can lead to these conflict. In a future version this issue will be resolved.

## Usage

Either, you can 1. download the already produced data from [Zenodo](https://doi.org/10.5281/zenodo.17098743) into the data folder and directly run the plotting scripts to obtain our figures or 2. you can run the ensemble runs yourself, save the data and then plot. 

If 2.:

1. **Adjust the SLURM script:**  
   Edit `submit_study.sh` to specify the correct script, resources, and environment for your HPC setup.

2. **Submit jobs:**  
   Run the submission script for each Julia file as needed:
   ```bash
   sbatch submit_study.sh

For Plotting and results:

- **Plotting:**  
  After simulations are complete, use the plotting scripts to create figures for analysis and publication and to obtain the IBRI results that will be saved in ../results/ folder.
