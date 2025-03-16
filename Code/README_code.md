## `Code` directory

The directory has R script files for simulation studies in Section 4 of the main text, the I-WIP meta-analysis study in Section 5, and generation of figures and tables. 

### `Simulation_1` directory

The directory contains R codes for the simulation study in Section 4.1 to implement four methods: 
  - `Simulation_1/1_Benchmark.R` Random-effects analysis using all 40 simulation IPD studies
  - `Simulation_1/2_IPD-AD.R` **Proposed integrated random-effects analysis** with 10 IPD studies and 30 AD studies
  - `Simulation_1/3_IPD-AD-pooled.R` Proposed integrated random-effects analysis with a mis-specification density ratio model
  - `Simulation_1/4_IPD_only.R` Random-effects analysis using only 10 IPD studies

**Input and output**
  - Each R script file uses the simulation data `Data/SimulationData_1.RData`
  - The default setting is to import the 1st simulation dataset (out of 300 repeated simulations) and run the MCMC. 
      - One can use other simulation dataset by changing `rep_no` (in Line 3 of each R code) to a number from 1 to 300.
      - The authors used a batch script to run the 300 repeated simulations in multiple cores in parallel.
  - Posterior values from the MCMC (after burn-in) are stored in `Output/Simulation_1/[Method]/RData` directory.
      - If you run an R code with the default setting, the output file for the 1st simulation data `rep_1.RData` will be produced.
      - For your convenience, the output files for all 300 simulation studies are already stored in the  `Output/Simulation_1/[Method]/RData` directory.
  - The current default setting `DrawDiagnostics = TRUE` on Line 14 of each R script file will produce the diagnostic plot of **μ** = E(**θ**<sub>l</sub>) (while running MCMC) and store it in `Output/Simulation_1/[Method]/DiagnosticPlot` directory.

### `Simulation_2` directory

The directory contains R codes for the simulation study in Section 4.2. The structure of this directory is the same with that of `Simulation_1` directory.

### `I-WIP_Application_CodeOnly.R` file

The R script was used to implement our proposed method on the real dataset for I-WIP meta-analysis study in Section 5. 
  - Due to public access restrictions, the data were not released here.
  - The complete dataset is available from the data custodian (OOO) at OOO@OOO, subject to the Terms and Conditions of Data Transfer.
