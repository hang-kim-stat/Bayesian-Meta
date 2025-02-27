# **Supplementary Materials for ''Bayesian Random-Effects Meta-Analysis Integrating Individual Participant Data and Aggregate Data''**  

This repository follows the same structure as [the JASA reproducibility materials template on GitHub](https://github.com/jasa-acs/repro-template).

## **Folder Structure**  



This README file introduces the structure of the repository. 


[Click Here](https://example.com)

(). 
This repository contains **R script files** and **data/result files** to reproduce **Figure 3** from the main text.  

## **Folder Structure**  
The repository includes two main folders:  

- 📁 **`Simulation1/`** – Contains files for simulation studies in **Section 4.1**  
- 📁 **`Simulation2/`** – Contains files for simulation studies in **Section 4.2**  

Both folders share the **same structure**, containing the following files:  

### **1️⃣ `0_SimulData.R`**  
- Generates **300 replications** of simulation datasets.  
- Stores the generated datasets in the **`0_SimulData.RData`** file.  

### **2️⃣ `2_IPD-AD.R`**  
- Implements the **proposed integrated analysis** (labeled as **"IPD-AD"**).  

#### **Usage Instructions**  
- By default, the script runs on the **1st dataset** out of **300 replications**.  
- To use a different dataset, **modify Line 3** in the script:  
  ```r
  rep_no <- X  # Replace X with a number between 1 and 300
  ```
- The authors used a **batch script** to run all **300 simulations** in parallel using multiple cores.  

---

## **📊 Results**  

For each simulation dataset **(replication X)**, the following output files are generated:  

1️⃣ **Trace Plots**  
- **`W_rep_X_mu.png`** → Shows trace plots of the **mean of random regression coefficients (μ)** for dataset **X**.  

2️⃣ **MCMC Chain Results (`Rep_X.RData`)**  
- Contains results from the MCMC chain (**n_iter = 20,000**) with the following objects:  
  - **`draw_mu`** → `(n_iter × 4)` matrix, posterior draws of the mean of **random regression θ_l**.  
  - **`draw_Sigma_theta`** → `(n_iter × 4)` matrix, posterior draws of the **variance of random regression θ_l** (diagonal covariance matrix).  
  - **`draw_sig2`** → `(n_iter)` vector, posterior draws of the **variance of random errors (ε_li)**.  

---

This version improves clarity with **headings, bullet points, emojis, and consistent formatting** to make it **easier to read** on GitHub. 🚀 Let me know if you’d like any further refinements! 😊









Folder "Simulation2": files used for simulation studies in Section 4.1 


- `figure4_sensitivity_full_GMSS_reproduce.R`: This code first performs the GMSS and then performs sensitivity analysis using the intervention posterior. The GMSS sampling procedure takes approximately 3 days to run on Ohio Supercomputer Center cloud.

- `figure4_sensitivity_reproduce.R`: This code first downloads the results of the GMSS from a cloud storage and then conducts sensitivity analysis using the intervention posterior. The code detects the available cores on the machine and utilizes all cores except one. This code takes approximately 17 minutes to complete on a 16-core desktop (AMD 5950X), while utilizing 15 cores.

When the codes complete, they produce `overall_period_plot_combined.png` in the working directory, which corresponds to Figure 4 of the manuscript. Intermediate computation results are stored in `sensitivity_analysis` directory, which is created while running the code.

Included in the `experimental_dat` directory are four rds files:

1. `field_data_3_rep.rds`: Data for the three field experiments.
2. `init_20.rds`: Starting points for the 20 multisets.
3. `Max_param.rds`: Specified ranges for parameters. 
4. `rules.rds`: Binary (high/low prospects) rules for instrumental densities obtained from prognostic experiments.

`circadian_funs.r` has the auxiliary functions that are necessary to run GMSS and ODE models. 



Title: Supplement for "OOO"

In this supplement, we include R script files and data/result files to reproduce Figure 3 in the main text: 

Folders "Simulation1" and "Simulation2" have files used for simulation studies in Section 4.1 and 4.2, respectively. The folders have the exact same structure with following files:

- `0_SimulData.R': This code generates 300 replications of simulation datasets and stores them in the `0_SimulData.RData' file. 

- `2_IPD-AD.R': This code implements the proposed integrated analysis (labeled as ‘IPD-AD’).

    * The default is using the 1st simulation dataset out of 300 replications. To use other simulation datasets, a user needs to replace the number 1 for "rep_no" in Line 3 with another number between 1 and 300. (The authors used a batch script to run the 300 repeated simulations in multiple cores parallel.)
    * Results
      1. "W_rep_X_mu.png" shows trace plots of the mean of random regression coefficients, mu, from the Xth dataset
      2. "Rep_X.RData" contains results of the MCMC chain with (n_iter = 20,000) iterations in three objects:
         * draw_mu: (n_iter by 4)-matrix, posterior draws of the mean of random regression theta_l
         * draw_Sigma_theta: (n_iter by 4)-matrix, posterior draws of the variance of random regression theta_l (diagonal covariance matrix)
         * draw_sig2: (n_iter)-vector, posterior draws of the variance of random errors epsilon_li
 
      
Here’s your **GitHub README** with improved formatting for better readability while keeping the content unchanged:  

---
