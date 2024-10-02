#!/bin/bash
#SBATCH --job-name=maps_tcv_global2
#SBATCH --cpus-per-task=1
#SBATCH --mem-per-cpu=1G
#SBATCH --time=00:30:00
#SBATCH --qos=30min
#SBATCH --mail-type=END,FAIL,TIME_LIMIT
#SBATCH --mail-user=marina.antillon@swisstph.ch

# Tell SLURM this is an array of jobs, 1-133%133
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
Rscript ./C_cea_code/code00_master_cluster.R maps_cov=0.2 vax_sens=2.00 mcv_sa="none" --default-packages=zeallot, magrittr, RcolorBrewer, scales, ggplot2, grid, ggrepel, patchwork, tidyr, dplyr, cubelyr, abind, readxl, R.matlab
Rscript ./C_cea_code/code00_master_cluster.R maps_cov=0.2 vax_sens=2.25 mcv_sa="none" --default-packages=zeallot, magrittr, RcolorBrewer, scales, ggplot2, grid, ggrepel, patchwork, tidyr, dplyr, cubelyr, abind, readxl, R.matlab
Rscript ./C_cea_code/code00_master_cluster.R maps_cov=0.2 vax_sens=3.00 mcv_sa="none" --default-packages=zeallot, magrittr, RcolorBrewer, scales, ggplot2, grid, ggrepel, patchwork, tidyr, dplyr, cubelyr, abind, readxl, R.matlab
Rscript ./C_cea_code/code00_master_cluster.R maps_cov=0.2 vax_sens=4.50 mcv_sa="none" --default-packages=zeallot, magrittr, RcolorBrewer, scales, ggplot2, grid, ggrepel, patchwork, tidyr, dplyr, cubelyr, abind, readxl, R.matlab

Rscript ./C_cea_code/code00_master_cluster.R maps_cov=0.1 vax_sens=2.00 mcv_sa="none" --default-packages=zeallot, magrittr, RcolorBrewer, scales, ggplot2, grid, ggrepel, patchwork, tidyr, dplyr, cubelyr, abind, readxl, R.matlab
Rscript ./C_cea_code/code00_master_cluster.R maps_cov=0.1 vax_sens=2.25 mcv_sa="none" --default-packages=zeallot, magrittr, RcolorBrewer, scales, ggplot2, grid, ggrepel, patchwork, tidyr, dplyr, cubelyr, abind, readxl, R.matlab
Rscript ./C_cea_code/code00_master_cluster.R maps_cov=0.1 vax_sens=3.00 mcv_sa="none" --default-packages=zeallot, magrittr, RcolorBrewer, scales, ggplot2, grid, ggrepel, patchwork, tidyr, dplyr, cubelyr, abind, readxl, R.matlab
Rscript ./C_cea_code/code00_master_cluster.R maps_cov=0.1 vax_sens=4.50 mcv_sa="none" --default-packages=zeallot, magrittr, RcolorBrewer, scales, ggplot2, grid, ggrepel, patchwork, tidyr, dplyr, cubelyr, abind, readxl, R.matlab

Rscript ./C_cea_code/code00_master_cluster.R maps_cov=0.3 vax_sens=2.00 mcv_sa="none" --default-packages=zeallot, magrittr, RcolorBrewer, scales, ggplot2, grid, ggrepel, patchwork, tidyr, dplyr, cubelyr, abind, readxl, R.matlab
Rscript ./C_cea_code/code00_master_cluster.R maps_cov=0.3 vax_sens=2.25 mcv_sa="none" --default-packages=zeallot, magrittr, RcolorBrewer, scales, ggplot2, grid, ggrepel, patchwork, tidyr, dplyr, cubelyr, abind, readxl, R.matlab
Rscript ./C_cea_code/code00_master_cluster.R maps_cov=0.3 vax_sens=3.00 mcv_sa="none" --default-packages=zeallot, magrittr, RcolorBrewer, scales, ggplot2, grid, ggrepel, patchwork, tidyr, dplyr, cubelyr, abind, readxl, R.matlab
Rscript ./C_cea_code/code00_master_cluster.R maps_cov=0.3 vax_sens=4.50 mcv_sa="none" --default-packages=zeallot, magrittr, RcolorBrewer, scales, ggplot2, grid, ggrepel, patchwork, tidyr, dplyr, cubelyr, abind, readxl, R.matlab

Rscript ./C_cea_code/code00_master_cluster.R maps_cov=0.2 vax_sens=2.00 mcv_sa="mcv1imp" --default-packages=zeallot, magrittr, RcolorBrewer, scales, ggplot2, grid, ggrepel, patchwork, tidyr, dplyr, cubelyr, abind, readxl, R.matlab
Rscript ./C_cea_code/code00_master_cluster.R maps_cov=0.2 vax_sens=2.25 mcv_sa="mcv1imp" --default-packages=zeallot, magrittr, RcolorBrewer, scales, ggplot2, grid, ggrepel, patchwork, tidyr, dplyr, cubelyr, abind, readxl, R.matlab
Rscript ./C_cea_code/code00_master_cluster.R maps_cov=0.2 vax_sens=3.00 mcv_sa="mcv1imp" --default-packages=zeallot, magrittr, RcolorBrewer, scales, ggplot2, grid, ggrepel, patchwork, tidyr, dplyr, cubelyr, abind, readxl, R.matlab
Rscript ./C_cea_code/code00_master_cluster.R maps_cov=0.2 vax_sens=4.50 mcv_sa="mcv1imp" --default-packages=zeallot, magrittr, RcolorBrewer, scales, ggplot2, grid, ggrepel, patchwork, tidyr, dplyr, cubelyr, abind, readxl, R.matlab

Rscript ./C_cea_code/code00_master_cluster.R maps_cov=0.2 vax_sens=2.00 mcv_sa="mcv2" --default-packages=zeallot, magrittr, RcolorBrewer, scales, ggplot2, grid, ggrepel, patchwork, tidyr, dplyr, cubelyr, abind, readxl, R.matlab
Rscript ./C_cea_code/code00_master_cluster.R maps_cov=0.2 vax_sens=2.25 mcv_sa="mcv2" --default-packages=zeallot, magrittr, RcolorBrewer, scales, ggplot2, grid, ggrepel, patchwork, tidyr, dplyr, cubelyr, abind, readxl, R.matlab
Rscript ./C_cea_code/code00_master_cluster.R maps_cov=0.2 vax_sens=3.00 mcv_sa="mcv2" --default-packages=zeallot, magrittr, RcolorBrewer, scales, ggplot2, grid, ggrepel, patchwork, tidyr, dplyr, cubelyr, abind, readxl, R.matlab
Rscript ./C_cea_code/code00_master_cluster.R maps_cov=0.2 vax_sens=4.50 mcv_sa="mcv2" --default-packages=zeallot, magrittr, RcolorBrewer, scales, ggplot2, grid, ggrepel, patchwork, tidyr, dplyr, cubelyr, abind, readxl, R.matlab
