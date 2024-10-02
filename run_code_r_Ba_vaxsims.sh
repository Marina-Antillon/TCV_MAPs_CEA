#!/bin/bash
#SBATCH --job-name=maps_tcv_global_vaxsims
#SBATCH --cpus-per-task=1
#SBATCH --mem-per-cpu=1G
#SBATCH --time=00:30:00
#SBATCH --qos=30min
#SBATCH --mail-type=END,FAIL,TIME_LIMIT
#SBATCH --mail-user=marina.antillon@swisstph.ch

# Tell SLURM this is an array of jobs
######################
#SBATCH --array=1-133%133

#load your required modules below
#################################
module load R/4.2.2-foss-2021a
 #part of ls scicore/soft/modules/lang; or else Rscript means nothing.

#export your required environment variables below
#################################################
# export R_LIBS_USER=$HOME/R_package:$R_LIBS_USER

#add your command lines below
#############################
Rscript ./B_vax_sims_code/code_vax_sims_global.R maps_cov=0.2 mcv_sa="none" --default-packages=Hmisc, deSolve, MASS, matrixcalc, emdbook, RcolorBrewer, ggplot2, epitools, LaplacesDemon, tidyr, dplyr, abind, readxl, R.matlab
Rscript ./B_vax_sims_code/code_vax_sims_global.R maps_cov=0.1 mcv_sa="none" --default-packages=Hmisc, deSolve, MASS, matrixcalc, emdbook, RcolorBrewer, ggplot2, epitools, LaplacesDemon, tidyr, dplyr, abind, readxl, R.matlab
Rscript ./B_vax_sims_code/code_vax_sims_global.R maps_cov=0.3 mcv_sa="none" --default-packages=Hmisc, deSolve, MASS, matrixcalc, emdbook, RcolorBrewer, ggplot2, epitools, LaplacesDemon, tidyr, dplyr, abind, readxl, R.matlab

Rscript ./B_vax_sims_code/code_vax_sims_global.R maps_cov=0.2 mcv_sa="mcv1imp" --default-packages=Hmisc, deSolve, MASS, matrixcalc, emdbook, RcolorBrewer, ggplot2, epitools, LaplacesDemon, tidyr, dplyr, abind, readxl, R.matlab
Rscript ./B_vax_sims_code/code_vax_sims_global.R maps_cov=0.1 mcv_sa="mcv1imp" --default-packages=Hmisc, deSolve, MASS, matrixcalc, emdbook, RcolorBrewer, ggplot2, epitools, LaplacesDemon, tidyr, dplyr, abind, readxl, R.matlab
Rscript ./B_vax_sims_code/code_vax_sims_global.R maps_cov=0.3 mcv_sa="mcv1imp" --default-packages=Hmisc, deSolve, MASS, matrixcalc, emdbook, RcolorBrewer, ggplot2, epitools, LaplacesDemon, tidyr, dplyr, abind, readxl, R.matlab

Rscript ./B_vax_sims_code/code_vax_sims_global.R maps_cov=0.2 mcv_sa="mcv2" --default-packages=Hmisc, deSolve, MASS, matrixcalc, emdbook, RcolorBrewer, ggplot2, epitools, LaplacesDemon, tidyr, dplyr, abind, readxl, R.matlab
Rscript ./B_vax_sims_code/code_vax_sims_global.R maps_cov=0.1 mcv_sa="mcv2" --default-packages=Hmisc, deSolve, MASS, matrixcalc, emdbook, RcolorBrewer, ggplot2, epitools, LaplacesDemon, tidyr, dplyr, abind, readxl, R.matlab
Rscript ./B_vax_sims_code/code_vax_sims_global.R maps_cov=0.3 mcv_sa="mcv2" --default-packages=Hmisc, deSolve, MASS, matrixcalc, emdbook, RcolorBrewer, ggplot2, epitools, LaplacesDemon, tidyr, dplyr, abind, readxl, R.matlab
