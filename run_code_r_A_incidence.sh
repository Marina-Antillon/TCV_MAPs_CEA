#!/bin/bash
#SBATCH --job-name=maps_tcv_fits
#SBATCH --cpus-per-task=1
#SBATCH --mem-per-cpu=2G
#SBATCH --time=05:59:00
#SBATCH --qos=6hours
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
Rscript ./A_incidence_code/code_R0_5WQ.R --default-packages=Hmisc, deSolve, MASS, matrixcalc, emdbook, RcolorBrewer, ggplot2, epitools, LaplacesDemon, tidyr, dplyr, abind, readxl, R.matlab 
