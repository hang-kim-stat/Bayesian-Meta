# **Supplementary Materials for ''Bayesian Random-Effects Meta-Analysis Integrating Individual Participant Data and Aggregate Data''**  

This repository follows the same structure as [the JASA reproducibility materials template on GitHub](https://github.com/jasa-acs/repro-template).

## `code` directory 

### Codes for simulation studies 

The directory has R codes for two simulation studies in Section 4 of the main text and our I-WIP meta-analysis study in Section 5. 

Each simulation study folder (`Simulation_1` and `Simulation_2`) contains codes for four methods: 
  - `1_Benchmark.R` Random-effects analysis using all 40 simulation IPD studies
  - `2_IPD-AD.R` Propsed meta-analysis with 10 IPD studies and 30 AD studies
  - `3_IPD-AD-pooled.R` Proposed meta-analysis with a mis-specification density ratio model
  - `4_IPD_only.R` Random-effects analysis using only 10 IPD studies

**Input**
  - The R codes use the simulation data stored in `data' folder.
  - A default setting is to import the 1st simulation dataset (out of 300 repeated simulations).
      - Other simulation dataset can be used by changing `rep_no` (in Line 3 of each R code) to a number from 1 to 300.
      - The authors used a batch script to run the 300 repeated simulations in multiple cores in parallel. 

**Output**
  - Posterior values from the MCMC (after burn-in) are stored in `output/[MethodName]/RData` directory.
      - If you run an R code in `code` directory, the output file for the 1st simulation data `rep_1.RData` will be produced.
      - For your convenience, the output files for all 300 simulation studies are already stored in `output/[MethodName]/RData` directory.
  - The current default setting `DrawDiagnostics = TRUE' on Line 14 will produce the diagnostic plot (while running MCMC) and store it in `output/[SimulationNumber]/[MethodName]/DiagnosticPlot` directory.
  - `W_Summary.R` generates two tables used for Figure 3 in the main text and Tables 5 and 6 in the Supplemenary Materials. 

### Code for I-WIP meta-analysis study

`I-WIP_Application_CodeOnly.R` is the code we used to apply our proposed method to the real data. 
  - The data are not released here for [  REASON ].  


## `data` directory 


## `output` directory 
