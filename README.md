Microarray patch vaccines for typhoid conjugate vaccines: a global cost-effectiveness analysis
==============================================================================================

Administered by Marina Antillon, PhD.  

Collaborating centers:

Department of Epidemiology and Public Health   
Swiss Tropical and Public Health (Swiss TPH) Institute   
(An affiliated institute of the University of Basel)   
4123 Allschwil, Switzerland  

PATH
Seattle, WA, USA

License: GPL-3.0 

For questions, please contact Marina Antillon via LinkedIn.

---
# Project Objective 

In the current study, we leverage a previously developed model of typhoid transmission (LINKS) to evaluate the impact and cost-effectiveness of a potential micro-array patch typhoid conjugate vaccine not yet under development. The analysis assumes that traditional (needle-and-syringe) typhoid conjugate vaccines would be rolled out before 2032 in most or all low- and middle-income countries and that MAPs would be non-inferior in terms of the probability and duration of protection. The analysis is also novel in taking into account the 

The analysis encompasses two major geographical scopes:  
A. a global analysis where the cost-effectiveness of TCV-MAPS is considered on a national basis  
B. a sub-national analysis where the cost-effectiveness of TCV-MAPS is considered on a district-by-district basis  

---
# Brief description

The analysis is broadly defined in four parts, each of four parts is executed in various R code files:  
I. Re-fitting incidence parameters to accommodate for 5 wealth quintiles in all countries  
II. Simulating vaccine impact for national-level analysis   
    - Simulating the impact of N&S vaccines from deployment in the 2020's through the 2030's and 2040's as a comparator  
    - Simulating the impact of MAPS vaccines in the 2030's  
    - Cost-effectiveness analysis  
III. Simulating vaccine impact for sub-national (district) level analysis in five example countries  
    - Simulating the impact of N&S vaccines from deployment in the 2020's through the 2030's and 2040's as a comparator  
    - Simulating the impact of MAPS vaccines in the 2030's  
    - Cost-effectiveness analysis  
V. Collating the results  
    - National-level analyses  
    - Subnational analyses  

---
# File structure, output, and dependencies

## Directory `./A_incidence_code/`
These are the file to re-fit the parameters of the analysis for a model that now contains wealth quintiles. For the explanation of choices made here, see the supplement. The analysis is an expansion of the paper by Bilcke et al [here](https://doi.org/10.1016/s1473-3099(18)30804-1).  

The analysis is called by the file `run_code_r_A_incidence.sh` and is meant to run with a slurm scheduler. This portion of the analysis can take up to 2 hours or more.

## Directory `./B_vax_sims_code/`
These are the files to simulate the vaccine interventions. First to simulate the impact of TCV-N&S in the 2020's and then to simulate additional benefits of rolling out TCV-MAPs after 2032. One file is to simulate the impact for each country of the global analysis (`code_vax_sims_global.R`) and the other is toe simulate the impact for each region or district in the five countries for which we perform subnational analysis (`code_vax_sims_dd.R` in Burkina Faso, India, Kenya, Malawi, and Nepal).

The analysis is called by the files `run_code_r_Ba_vaxsims.sh`, `run_code_r_Ba_vaxsims_dd_reg.sh`, and `run_code_r_Ba_vaxsims_districts.sh` and is meant to run with a slurm scheduler. This portion of the analysis takes about 30 minutes to run.

## Directory `./C_cea_code/`
The cost-effectiveness analyses of TCV-MAPs compare to TCV_N&S and different sensitivity analyses. 

The analysis is called by the files `run_code_Ca_cea.sh`, `run_code_Cb_cea_dd_reg.sh`, `run_code_Ca_cea_dd_dist.sh` and is meant to run with a slurm scheduler.  The `.sh` files call the files `./C_cea_code/code00_master_cluster.R` and `./C_cea_code/code00_master_cluster_dd.sh` for the global and subnational analyses, respectively, and those call the rest of the files in the `./C_cea_code/` directory. This portion of the analysis takes about 30 minutes to run.

## Directory `./data/`
The set of data necessary to run all the analyses.

## Directory `./functions/`
Helpful functions to run the analyses. 

## Directory `./geographic_maps/`
The set of `.shp` file that makes the maps in the manuscript for the subnational analysis. No need to upload it to the cluster.

## Directories `./post_cluster_summaries_global/` and `./post_cluster_summaries_subnational/`
These directories contain the code that is used to call all the output from the cluster and make summaries as they appear in the manuscript. 

## Directory `./renv/`
This directory contains the necessary files to replicate the package environment that was used to run the analysis in RStudio. In order to load these, the `TCV_MAPS_CEA.RProj` must be called in order to run the analysis as a project in RStudio. In case of issues, see the section "Troubleshooting" below. 

## `.sh` files for the slurm scheduler (to be run in the cluster)

- `run_code_r_A_incidence` runs the incidence code, re-fitting parameters to be able to run a model with wealth quintiles.
- `run_code_r_Ba_vaxsims.sh`, `run_code_r_Bb_vaxsims_dd_reg.sh`, and `run_code_r_Bc_vaxsims_districts.sh` run the vaccine simulations at the global level, the subnational analysis at the region level, and the subnational analysis at the district level. 
- `run_code_Ca_cea.sh`, `run_code_Cb_cea_dd_reg.sh`, `run_code_Cc_cea_dd_dist.sh` run the cost-effectiveness analysis at the global level, the subnational analysis at the region level, and the subnational analysis at the district level. The `.sh` files call the files `./C_cea_code/code00_master_cluster.R` and `./C_cea_code/code00_master_cluster_dd.sh` for the global and subnational analyses, respectively, and those call the rest of the files in the `./C_cea_code/` directory.

## `TCV_MAPS_CEA.Rproj`
The `.RProj` file that makes the directory a proper project in RStudio. Only for use in RStudio. Not necessary if the analysis is to be run with 

---
# Computational considerations

## Hardware

**Hardware needs:** For optimal performance, the model was run for different places and different scenarios using a high-performance computing cluster at the University of Basel (Scicore, http://scicore.unibas.ch/). Scicore is run with a slurm scheduler, and the bash file is included in this repository. However, a user could use a parallel computing package in R or any other automated task management solution, but no such implementation is presented here.

**Memory needed:** Intermediate simulations and graphs need about 10 GB of storage for all outputs.

## Software
The tools for the analysis are coded in R. While we would highly recommend using the code within RStudio environment (in part because of its features to manage the project with .RProj and packrat) this is not strictly necessary and the benefits of git and renv are available from a classic R interface and a shell command line.

Some of the results tables for the project are produced automatically with the code in the project. This is done by employing the knitr and kableExtra packages. 

A single run of the code takes about 60 minutes in a MacBook Air (M1 chip, 2020 model) with 16 GB of RAM. 

---
# Installation to-do list (all free)
- [R](https://www.r-project.org) (required)
- [renv](https://rstudio.github.io/renv/index.html) (packages installed for R; required)
- [RStudio](https://www.rstudio.com/) (highly recommended though not required; it integrates well with the other management tools used in the project. The free version is enough.)
- Latex engine (not required, though some documentation created automatically won't run). See [here](https://support.rstudio.com/hc/en-us/articles/200532257-Customizing-LaTeX-Options) for the recommended engines.

---
# How to run this analysis 

## Part I: on the cluster

### 1. Upload the subdirectory to the cluster 
### 2. Re-parameterizing the model to accommodate wealth quintiles  
In the cluster, run _run_code_r_A_incidence.sh_ on terminal: `sbatch run_code_r_A_incidence.sh`. Adjust settings (for instance # of countries to run, if there are different numbers to run). This bash code calls `./A_incidence_code/code_R0_5WQ.R` which produces a directory `./out_fits/` with one subdirectory per country labeled with the IS03 label, and within each of those directories there is one file that is called `./fit_global.Rdata` which holds the new parameters for the simulations of vaccine impact.  
    – other helpful bash code: `squeue -u <username>` to see what jobs are running.  
    – other helpful bash code: `find . -type f -name "slurm*" -exec tail -1 {} \; -print` to see the last line of the slurm jobs to check which ones ended in an error.  
    – NOTE: change the `-1` to any number of last lines. I used this mostly to check which countries had failed. Then I would go to the output file to see why they had failed.  
    
### 3. Vaccine simulations 
This step is a little different if one wants to run the global analysis or the subnational analysis with regions or with districts. Note that in Burkina Faso, the subnational analysis will only run with regions.   

a. **Global Analysis:** In the cluster, run _run_code_r_Ba_vaxsims.sh_ on terminal: `sbatch run_code_r_Ba_vaxsims.sh`. Adjust settings (for instance # of countries to run, if there are different numbers to run). This bash code calls `./B_vax_sims_code/code_vax_sims_global.R` which produces a directory `./out_global/` with one subdirectory per sensitivity analysis (a combination of the MAPs coverage of the unvaccinated and the assumption about N&S coverage). Within each of those subdirectories, there will be one subdirectory per country labelled with the IS03 label, and within each of those directories, there are four files that hold the simulations for no MAPs and MAPs deployment with each of the three comparators.   

b. **Subnational analysis with regions:** In the cluster, run _run_code_r_Bb_vaxsims_dd_reg.sh_ on terminal: `sbatch run_code_r_Bb_vaxsims_dd_reg.sh`. Adjust settings (for instance # of regions to run, if there are different numbers to run). This bash code calls `./B_vax_sims_code/code_vax_sims_dd.R` which produces a directory `./out_deepdive_vax/reg/all` with one subdirectory per sensitivity analysis for the MAPs coverage of unvaccinated individuals (0.1, 0.2-the default, and 0.3). Within each subdirectory, there will be one subdirectory per country labelled with the IS03 label, and within each of those there will be one subdirectory per subnational region. That last subdirectory will have four files that hold the simulations for no MAPs and MAPs deployment with each of the three comparators.  

c. **Subnational analysis with districts (not possible for Burkina Faso):** In the cluster, run _run_code_r_Bc_vaxsims_dd_dist.sh_ on terminal: `sbatch run_code_r_Bc_vaxsims_dd_dist.sh`. Adjust settings (for instance # of regions to run, if there are different numbers to run). This bash code calls `./B_vax_sims_code/code_vax_sims_dd.R` which produces a directory `./out_deepdive_vax/dist/all` with one subdirectory per sensitivity analysis for the MAPs coverage of unvaccinated individuals (0.1, 0.2-the default, and 0.3). Within each subdirectory, there will be one subdirectory per country labeled with the IS03 label, and within each of those, there will be one subdirectory per subnational region. That last subdirectory will have four files that hold the simulations for no MAPs and MAPs deployment with each of the three comparators.  

### 4. CEA analyses  

a. **Global analysis:** In the cluster, run _run_code_r_Ca_cea.sh_ on terminal: `sbatch run_code_Ca_cea.sh`. Adjust settings (for instance # of countries to run, if there are different numbers to run). This bash code calls `./C_cea_code/code00_master_cluster.R` which calls the rest of the files in the `./C_cea_code/` directory. This portion of the analysis produces a directory `./out_global_cea/` with one subdirectory per sensitivity analysis (a combination of the MAPs coverage of the unvaccinated and the assumption about N&S coverage). Within each of those subdirectories, there will be a subdirectory for the cost of MAPs, and within those, there will be one subdirectory per country labelled with the IS03 label. Within each of the country-level subdirectories, there are intermediate and final output files for the cost-effectiveness analysis.  

b. **Subnational analysis by regions:** In the cluster, run _run_code_r_Cb_cea_dd_reg.sh_ on terminal: `sbatch run_code_r_Cb_cea_dd_reg.sh`. Adjust settings (for instance # of regions to run, if there are different numbers to run). This bash code calls `./C_cea_code/code00_master_cluster_dd.R` which calls the rest of the files in the `./C_cea_code/` directory. This portion of the analysis produces a directory `./out_deepdive_cea/reg/all` with one subdirectory per sensitivity analysis for the MAPs coverage of unvaccinated individuals (0.1, 0.2-the default, and 0.3). Within those, there will be one subdirectory for each of the prices of the MAPs that are examined in the analysis. Within each subdirectory, there will be one subdirectory per country labeled with the IS03 label, and within each of those there will be one subdirectory per subnational region.Within each region-level subdirectory, there are intermediate and final output files for the cost-effectiveness analysis. 

c. **Subnational analysis by district (not possible for Burkina Faso):** In the cluster, run _run_code_r_Cc_cea_dd_dist.sh_ on terminal: `sbatch run_code_r_Cc_cea_dd_dist.sh`. Adjust settings (for instance # of districts to run, if there are different numbers to run). This bash code calls `./C_cea_code/code00_master_cluster_dd.R` which calls the rest of the files in the `./C_cea_code/` directory. This portion of the analysis produces a directory `./out_deepdive_cea/dist/all` with one subdirectory per sensitivity analysis for the MAPs coverage of unvaccinated individuals (0.1, 0.2-the default, and 0.3). Within those, there will be one subdirectory each for each of the prices of the MAPs that are examined in the analysis. Within each of those subdirectories, there will be one subdirectory per country labeled with the IS03 label, and within each of those there will be one subdirectory per subnational district. Within each region-level subdirectory, there are intermediate and final output files for the cost-effectiveness analysis.
   
## Part II: on the local computer
### 5. Post-cluster output  

a. **Global analysis:** run these files in any order: `code06a_global_epibaselines.R` and `code06b_global_cea.R`. The two pieces of code do not depend on the outputs of the other. The first piece summarizes how the model behaves with wealth quintiles and what we expect the typhoid landscape to look like just before MAPs are deployed in 2033, and the second piece of code summarizes the impact, costs, and cost-effectiveness of MAPs introduction.  

b. **Subnational analysis:** run these files in this precise order: `code006a_deepdive_graphs_dist_strat_by_costs.R` and after that, `code06b_deepdive_graphs_dist_cost_sidebyside.R`. The second piece of code depends on the outputs of the first piece of code. These files call the district level cluster outputs for India, Kenya, Malawi, and Nepal, and the region-level outputs for Burkina Faso.   
    
---

# Troubleshooting

**Issues with the renv repository:** If there is an issue with the repository, try typing renv:: rinit() and when prompted type the number 1, for Restored the project from lockfile.

**Adding a package to the renv library:** just use the command install.packages(). Note that doing this won't install the package for use with other projects. Then update the lockfile by typing renv::snapshot().
