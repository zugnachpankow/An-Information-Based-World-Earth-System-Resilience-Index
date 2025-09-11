# src

IMPORTANT NOTE: these files are not supposed to be run by themselves.

This folder contains the source code for the models and analyses used in the study. Each file implements a different model variant or analysis routine, primarily using Julia and ModelingToolkit. 

## Contents
- `ibri_analyses.jl`: Analysis functions for resilience and system dynamics.
- `ibri_drift_model.jl`: Implements the drift model variant.
- `ibri_shock_analyses.jl`: Analysis routines for shock scenarios.
- `ibri_shock_drift_model.jl`: Combines shock and drift dynamics.
- `ibri_shock_model.jl`: Implements the shock model variant.
- `ibri_std_model.jl`: Standard model implementation.

## Naming Convention
Files are named with the prefix `ibri_` followed by their main purpose or model type, for clarity and modularity.

## Usage
Import or include these files in your analysis scripts to construct and run different model scenarios.
