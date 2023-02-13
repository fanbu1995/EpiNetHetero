### "EpiNetHetero"
Code repository for [Likelihood-based inference for partially observed epidemics with individual heterogeneity](https://arxiv.org/abs/2112.07892)

This contains the core code used for simulation experiments in the manuscript. 
As the real data are proprietary, we cannot share them publicly, and the real data analysis code isn't included either 
(but the analysis code is essentially the same as that for real data).

### 1. Simulate data
Run the last two code chunks in the `HeteroEpiNet.py` script to obtain example datasets for `N=200` and `N=500` sized populations.

Note that you have to **specify the `savepath` argument as your desired data directory**.

### 2. Run inference (the basics)
Run the `infer_sim_data.R` script to compile all necessary functions.
And then run the example code (with any desired modifications) attached at the bottom of `infer_sim_data.R`.

Note that you have to **set the `data_root` argument to where the example data folder is located**.

### 3. Run inference with variance estimation
Run the `hetero_coverage_study.R` script first, and then run necessary commands (refer to the example code at the end of the script).

Note that to get pooled estimates, multiple independent runs of the stEM have to be performed. 
(There is example code of handling multiple run results to produce variance estimates and 95% Wald-type confidence intervals.)
