## Cost-effectiveness without uncertainty

# setwd("./maps_tcv_global")

source("./cea_code/code00b_load_packages.R")
# library(formattable)
# library(plyr)
# library(wbstats)

source("./cea_code/code00c_load_helper_functions.R")

# Had to update this:
# setRepositories(addURLs = c(CRANxtras = "http://lib.stat.cmu.edu/R/CRAN/bin/macosx/mavericks/contrib/3.3/"))

# Treatment model samples ----
## Run code that makes distribution for the treatment and cost of illness parameters
# iterations = 100
# keep samples common for countries within the same who region
CN = readMat("./data/burden_dynamic_link_Yale.mat")$cn
ISO = read.csv("./data/iso.csv")

# Bring in list of countries for which we do have vax simulations.
load("ISO_which_to_run_vaxsims.Rdata")

# bring in WB_group (with Countries.xlsx from MMGH) and Continent... and NS/MAPS intro dates
mmgh_data_ns_intro = read_xlsx("./data/Countries.xlsx", "Country Data w intro date")
colnames(mmgh_data_ns_intro) = gsub(" - ", "_",colnames(mmgh_data_ns_intro))
colnames(mmgh_data_ns_intro) = gsub(" ", "_",colnames(mmgh_data_ns_intro))
colnames(mmgh_data_ns_intro) = gsub("ISO_code", "ISO",colnames(mmgh_data_ns_intro))
colnames(mmgh_data_ns_intro) = gsub("-", "",colnames(mmgh_data_ns_intro))
colnames(mmgh_data_ns_intro) = gsub("\\(", "",colnames(mmgh_data_ns_intro))
colnames(mmgh_data_ns_intro) = gsub(")", "",colnames(mmgh_data_ns_intro))
colnames(mmgh_data_ns_intro) = gsub("WB_Group_June_2021", "WB_group",colnames(mmgh_data_ns_intro))
mmgh_data_ns_intro$NS_intro = sapply(1:183, function(i){max(mmgh_data_ns_intro[i, c("Low_intro_date", "Medium_Intro_date", "High_Intro_date")])})
mmgh_data_ns_intro$NS_intro[mmgh_data_ns_intro$NS_intro==0 & mmgh_data_ns_intro$WB_group %in% c("Low income", "Lower middle income")] = 2030
mmgh_data_ns_intro$NS_intro[mmgh_data_ns_intro$NS_intro==0 & mmgh_data_ns_intro$WB_group == "Upper middle income"] = 2030
#Continent TODO: change EMRO to respective continents
mmgh_data_ns_intro = mmgh_data_ns_intro %>% 
                        mutate(continent = ifelse(mmgh_data_ns_intro$WHO %in% c("WPRO", "SEARO"), "Asia", 
                                      ifelse(mmgh_data_ns_intro$WHO=="EMRO", "Mideast",
                                      ifelse(mmgh_data_ns_intro$WHO=="AFRO", "Africa", 
                                         ifelse(mmgh_data_ns_intro$WHO=="EURO", "Eurasia", "Americas")))))

tmp_z <- strtoi(Sys.getenv("SLURM_ARRAY_TASK_ID", unset=Sys.getenv("LINE_NUMBER", unset="1")))
z = which(which(ISO$countryiso==goodruns[tmp_z])==CN)

amr_data_unsd_region = read_xlsx("./data/AMR_data_Carey.xlsx", "Sheet1")

# the years of vaccination and of maps
intro_map = mmgh_data_ns_intro$TCVMAP_adoption_year[mmgh_data_ns_intro$ISO==ISO[CN[z],2]]
intro_sep = mmgh_data_ns_intro$TCVMAP_adoption_year[mmgh_data_ns_intro$ISO==ISO[CN[z],2]] - 
  mmgh_data_ns_intro$NS_intro[mmgh_data_ns_intro$ISO==ISO[CN[z],2]] # time between NS and MAPS
intro_sep[intro_sep>19] = 19

years=intro_sep:(intro_sep+19)
years_lbl= as.character(intro_map:(intro_map+19))
# scenario_ns = 2 # 1 for no vaccine, 2 for routine, 3 for routine + campaign for N&S
# scenario_switch = 1 # comparator

# args <- commandArgs()
# print(args)

# How it looks when it is running as a batch: 
# [1] "/scicore/soft/apps/R/4.2.2-foss-2021a/lib64/R/bin/exec/R"
# [2] "--no-echo"                                               
# [3] "--no-restore"                                            
# [4] "--file=./testing_command_args.R"                         
# [5] "--args"                                                  
# [6] "vax_sens=2.25"    

# pars_sens=gsub("--pars_sens=","",commandArgs()[grep("--pars_sens=", commandArgs())])
# eval(parse(text=commandArgs()[6])) # maps_cov
# eval(parse(text=commandArgs()[7])) # vax_sens <--- price of vaccine

maps_cov = as.numeric(gsub("maps_cov=","",commandArgs()[grep("maps_cov=", commandArgs())]))
vax_sens = as.numeric(gsub("vax_sens=","",commandArgs()[grep("vax_sens=", commandArgs())]))
mcv_sa = gsub("mcv_sa=","",commandArgs()[grep("mcv_sa=", commandArgs())])

outfolder = paste0("../out_global_cea/", round(maps_cov, 1), ifelse(mcv_sa!="none", paste0("_", mcv_sa), ""), "/", vax_sens, "dollars/")
inputfolder = paste0("../out_global/", round(maps_cov, 1), ifelse(mcv_sa!="none", paste0("_", mcv_sa), ""), "/")

subdir = paste0(outfolder, ISO[CN[z],2], "/figures/")
create_dir_if_necessary(subdir)

keep = c("CN", "ISO", "z", "keep", "years", "years_lbl", "mmgh_data_ns_intro", "vax_sens", 
         "outfolder", "inputfolder", "subdir", "amr_data_unsd_region")

# CEA parameter inputs ----

source('./cea_code/code01_CEA_inputs_no_unc.R')

fix$ns_price = 1.5
fix$maps_price = vax_sens

## save
save(z, unc, fix, threshold, lbl,  # pardata, pars, 
     years, years_lbl, # scenario_ns, scenario_switch, 
     file=paste0(outfolder, ISO[CN[z],2], '/trt_inputs.Rdata'))

# Dynamic model samples ----
rm(list=ls()[!(ls() %in% keep)]) # keep only those things in the important "keep" list from the top
gc()
load(paste0(outfolder, ISO[CN[z],2], '/trt_inputs.Rdata'))

source("./cea_code/00c_load_helper_functions.R")

load(paste0("../out/", ISO$countryiso[CN[z]], "/fit_global.Rdata"))

source("./cea_code/code02_dynamic_inputs_no_unc.R")
# somehow indicate when the maps time starts. That's year 0.
save(epipar, typhoid, doses, file=paste0(outfolder, ISO[CN[z],2], '/epi_inputs.Rdata'))

## Treatment ----
rm(list=ls()[!(ls() %in% keep)]) # keep only those things in the important "keep" list from the top
gc()

# begcea = Sys.time()

load(paste0(outfolder, ISO[CN[z],2], '/trt_inputs.Rdata'))
load(paste0(outfolder, ISO[CN[z],2], '/epi_inputs.Rdata'))
source("./cea_code/00c_load_helper_functions.R")

source('./cea_code/code03_treatment.R')
# endcea = Sys.time()

save(trt_model, trt_sims,  trt_sims_cond, trt_cost_sims_branch, yld_sims_branch,
     avg_cost, avg_deaths, avg_yld, 
     par_array_col_name,
     trt_costs, trt_costs_horizon, trt_costs_averted_horizon,
     yld_agegrp, yll_agegrp,
     daly, daly_agegrp, daly_averted_horizon, daly_horizon,
     deaths_agegrp, 
     file=paste0(outfolder, ISO[CN[z],2], '/code03_treatment.Rdata'))

## Intervention ----
rm(list=ls()[!(ls() %in% keep)]) # keep only those things in the important "keep" list from the top
gc()

# begcea = Sys.time()
load(paste0(outfolder, ISO[CN[z],2], '/trt_inputs.Rdata'))
load(paste0(outfolder, ISO[CN[z],2], '/epi_inputs.Rdata'))
source("./cea_code/00c_load_helper_functions.R")

source('./cea_code/code04_intervention.R')
# endcea = Sys.time()

save(int_costs, int_costs_dif_horizon, int_costs_horizon, doses_ns_maps,
     file=paste0(outfolder, ISO[CN[z],2], '/code04_intervention.Rdata'))

## CEA analysis ----

rm(list=ls()[!(ls() %in% keep)]) # keep only those things in the important "keep" list from the top
gc()
load(paste0(outfolder, ISO[CN[z],2], '/trt_inputs.Rdata'))
load(paste0(outfolder, ISO[CN[z],2], '/epi_inputs.Rdata'))
load(paste0(outfolder, ISO[CN[z],2], '/code03_treatment.Rdata'))
load(paste0(outfolder, ISO[CN[z],2], '/code04_intervention.Rdata'))
source("./cea_code/00c_load_helper_functions.R")

source('./cea_code/code05_cea_analysis.R')
# endcea = Sys.time()

save(icers_dalys_costs_summary, icers_dalys_costs_all_summary, 
     file=paste0(outfolder, ISO[CN[z],2], '/code05_cea.Rdata'))

## CEA make tables to report ----

# rm(list=ls()[!(ls() %in% keep)]) # keep only those things in the important "keep" list from the top
# gc()
# load(paste0('../out/', ISO[CN[z],2], '/trt_inputs.Rdata'))
# load(paste0('../out/', ISO[CN[z],2], '/epi_inputs.Rdata'))
# load(paste0('../out/', ISO[CN[z],2], '/code03_treatment.Rdata'))
# load(paste0('../out/', ISO[CN[z],2], '/code04_intervention.Rdata'))
# load(paste0('../out/', ISO[CN[z],2], '/code05_cea.Rdata'))
# source("./cea_code/00c_load_helper_functions.R")

# source('./cea_code/code05_make_simple_summaries.R')
