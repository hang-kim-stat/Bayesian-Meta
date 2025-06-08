## `Code` directory

The directory has R script files for simulation studies in Section 4 of the main text, the I-WIP meta-analysis study in Section 5, and the generation of figures and tables. 

### `Simulation_1` directory

The directory contains R codes for the simulation study in Section 4.1 to implement four methods: 
  - `1_Benchmark.R` Random-effects analysis using all 40 simulated IPD studies
  - `2_IPD-AD.R` **Proposed integrated random-effects analysis** with 10 IPD studies and 30 AD studies
  - `3_IPD-AD-pooled.R` Proposed integrated random-effects analysis with a misspecified density ratio model
  - `4_IPD_only.R` Random-effects analysis using only 10 IPD studies

**Input and output**
  - Each R script file uses the simulated data file `Data/SimulationData_1.RData`
  - With the default setting, each R code imports the **1st** simulated dataset (out of 300 repeated simulations) and apply each method. 
      - A reviewer can fit a method to another simulated dataset by changing `rep_no` (in Line 3 of each R code) to a number between 1 and 300.
      - The authors used a batch script to conduct the 300 repeated simulation studies in multiple cores in parallel.
  - Posterior values from the MCMC (after burn-in) are stored in `Output/Simulation_1/[Method]/RData` directory.
      - If you run the provided R code with the default setting, the output file for the 1st simulated data `rep_1.RData` will be produced.
      - For reviewer's convenience, the output files for all 300 simulation studies are already stored in the  `Output/Simulation_1/[Method]/RData` directory.
  - The current default setting `DrawDiagnostics = TRUE` on Line 14 of each R script file will produce the diagnostic plot of **μ** = E(**θ**<sub>l</sub>) (while running MCMC) and store it in `Output/Simulation_1/[Method]/DiagnosticPlot` directory.

### `Simulation_2` directory

The directory contains R codes for the simulation study in Section 4.2. Its structure is the same as that of `Simulation_1` directory.

### `Application_CodeOnly.R` file

The R script was used to implement our proposed method on the real dataset for I-WIP meta-analysis study in Section 5. 
  - Due to public access restrictions, the data were not released here.
  - The complete dataset is available from the data custodian (Queen Mary University of London) at smd-iwipdata@qmul.ac.uk, subject to the Terms and Conditions of Data Transfer.

### `Figure_Table` directory

The directory contains the R codes to generate tables in the Supplementary Material and figures in the main text. 
  - `Supplement_Table4.R` generates Table 4 of the Supplementary Material, storing to `Output/Figure_Table/Supplement_Table4.csv`
  - `Supplement_Table5.R` generates Table 5 of the Supplementary Material, storing to `Output/Figure_Table/Supplement_Table5.csv`
  - `Main_Figure3.R` generates Figure 3 of the main text, storing to `Output/Figure_Table/Maintext_Figure3.png`
  - `Main_Figures4and5.R` generates Figures 4 and 5 of the main text, storing to `Output/Figure_Table/Maintext_Figure4.png` and `Output/Figure_Table/Maintext_Figure5.png`
