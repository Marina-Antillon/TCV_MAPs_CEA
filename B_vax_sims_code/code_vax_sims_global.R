library("Hmisc")
library('deSolve')
require("MASS")
require("matrixcalc")
library("emdbook")

library("RColorBrewer")
library("ggplot2")

library("epitools")
library("LaplacesDemon") # so it will give me a dirichlet distribution

library("tidyr")
library("dplyr")
library("abind")

library("readxl")
library("R.matlab")

# Analysis settings -----
maps_cov = as.numeric(gsub("maps_cov=","",commandArgs()[grep("maps_cov=", commandArgs())]))
mcv_sa = gsub("mcv_sa=","",commandArgs()[grep("mcv_sa=", commandArgs())])

# Read in functions -----
fls = list.files("./functions/", pattern="^[fcn]")
for (i in 1:length(fls)){source(paste0("./functions/", fls[i]))}

meanpi = function(x){c(mean=mean(x), lci=quantile(x, 0.025), hci=quantile(x, 0.975))}

# Set up high-level subdirs -----
subdir = paste0("../out_global/", round(maps_cov, 1), ifelse(mcv_sa!="none", paste0("_", mcv_sa), ""), "/")
create_dir_if_necessary(subdir)

#***************************************************************
## Step 0: Bring in data ----
#***************************************************************

TransPar = readMat("./data/model_samples_ALL_Aug2021.mat")
CN = readMat("./data/burden_dynamic_link_Yale.mat")$cn
agedist = readMat("./data/burden_dynamic_link_Yale.mat")$agedist.country

ISO = read.csv("./data/iso.csv")

# load("ISO_which_to_run.Rdata")

#Bring in coverage from DHS
indicators = read_xlsx("./data/DHS_MICS_GlobalAnalysis_updated.xlsx", "data_to_r")
indicators$ISO = countrycode::countrycode(indicators$Country, "country.name", "iso3c")
indicators$ISO[indicators$Country=="Kosovo"] = "KSV"
indicators$ISO[indicators$Country=="Dominica"] = "DMA" # there was an error with this before 22 May

# Country data from MMGH
Countries <- read_excel("./data/Countries.xlsx", sheet = "Country Data w intro date")
colnames(Countries) = gsub(" - ", "_",colnames(Countries))
colnames(Countries) = gsub(" ", "_",colnames(Countries))
colnames(Countries) = gsub("ISO_code", "ISO",colnames(Countries))
colnames(Countries) = gsub("-", "",colnames(Countries))
colnames(Countries) = gsub("\\(", "",colnames(Countries))
colnames(Countries) = gsub(")", "",colnames(Countries))
colnames(Countries) = gsub("WB_Group_June_2021", "WB_group",colnames(Countries))
Countries$NS_intro = sapply(1:183, function(i){max(Countries[i, c("Low_intro_date", "Medium_Intro_date", "High_Intro_date")])})
Countries$NS_intro[Countries$NS_intro==0 & Countries$WB_group %in% c("Low income", "Lower middle income")] = 2030
Countries$NS_intro[Countries$NS_intro==0 & Countries$WB_group == "Upper middle income"] = 2030
#Continent TODO: change EMRO to respective continents
Countries = Countries %>% 
                mutate(continent = ifelse(Countries$WHO %in% c("WPRO", "SEARO"), "Asia", 
                            ifelse(Countries$WHO=="EMRO", "Mideast",
                                   ifelse(Countries$WHO=="AFRO", "Africa", 
                                          ifelse(Countries$WHO=="EURO", "Eurasia", "Americas")))))
Countries$ISO[Countries$Country=="Dominica"] = "DMA" # there was an error with this before 22 May

# For those missing vaccination statistics, of mcv1 estimates from before 2010, let's replace it with 
# the estimate from WHO coverage estimates in 2021. Or whatever year MMMGH used (2020? 2019?)
# these countries we won't be able to do a ECEA, just a CEA, but that's ok.
# assume the difference between lowest to highest is 5%
mcv1_WHO = read.csv("./data/mcv1_cov_WHO.csv")

indicators = right_join(indicators, data.frame(ISO=Countries$ISO, WB_group=Countries$WB_group)) %>% 
  filter(WB_group!="High income") %>% 
  mutate(across(contains("sani_wq"), ~as.numeric(.x))) %>% 
  mutate(across(contains("water_wq"), ~as.numeric(.x))) %>% 
  mutate(across(contains("mcv1_wq"), ~as.numeric(.x))) %>% 
  mutate_at(c("mcv1_total", "sani_total", "water_total"), as.numeric)

for(k in which(is.na(indicators$mcv1_total))){
  indicators$mcv1_total[k] = ifelse(indicators$ISO[k] %in% mcv1_WHO$SpatialDimValueCode,
                                    mcv1_WHO$FactValueNumeric[mcv1_WHO$SpatialDimValueCode==indicators$ISO[k] & mcv1_WHO$Period==2021],
                                    78.3) # the average
}

indicators$mcv1_dif = indicators$mcv1_wq5-indicators$mcv1_wq1
indicators$mcv1_difa = indicators$mcv1_total-indicators$mcv1_wq1
indicators$mcv1_difb = indicators$mcv1_wq5-indicators$mcv1_total

# the other ones, for distribution
indicators = indicators %>% 
  mutate(mcv1_wq1 = ifelse(is.na(mcv1_wq1), pmax(mcv1_total-4, 0.1), mcv1_wq1),
         mcv1_wq2 = ifelse(is.na(mcv1_wq2), pmax(mcv1_total-2, 0.1), mcv1_wq2),
         mcv1_wq3 = ifelse(is.na(mcv1_wq3), pmin(mcv1_total, 99.9), mcv1_wq3),
         mcv1_wq4 = ifelse(is.na(mcv1_wq4), pmin(mcv1_total+2, 99.9), mcv1_wq4),
         mcv1_wq5 = ifelse(is.na(mcv1_wq5), pmin(mcv1_total+4, 99.9), mcv1_wq5)) %>% 
  mutate(across(contains("mcv1_wq"), ~pmax(pmin(.x,99.9),0.1)))

# if this is the sensitivity analysis for MCV2, then subtract 25 percentage points 
# from MCV1 coverage. Then make sure to bottom code.
# if this is the sensitivity analysis for MCV1 improvements, then decrease the 
# unvaccinated by 20% "naturally" in each WQ (not concentrating on WQ1+WQ2) 

## MMGH assumptions
# Forecasted coverage: MCV1 was used as a proxy in view of the likely time of 
# administration of TCV vaccine. 2019 MCV1 was projected to 2042 using the 
# standard MI4A growth rates assumption. If countries had MCV1 coverage of <70% 
# then MCV1 was increased by 3% per year; if countries had MCV1 coverage between 
# 70-85%, then a 1% per year growth was used; if countries had MCV1 coverage 
# > 85%, a 0.5% per year growth was used. [22] 
# Easier way: see DHS_mcv1_mcv2_national.xlsx 3rd sheet. 
# about 40% of the unvaccinated will be vaccinated by 2032. So make that the increase

if(mcv_sa=="mcv2"){
  indicators = indicators %>% 
    mutate(across(contains("mcv1_wq"), ~(.x-25))) %>% 
    mutate(across(contains("mcv1_wq"), ~pmax(pmin(.x,99.9),0.1)))
} else if(mcv_sa=="mcv1imp"){
  indicators = indicators %>% 
    mutate(across(contains("mcv1_wq"), ~(.x+(100-.x)*0.4))) %>% 
    mutate(across(contains("mcv1_wq"), ~pmax(pmin(.x,99.9),0.1)))
}

# Bring in UC proportions and MAPS intro from MMGH
mmgh_data_maps = read_xlsx("./data/2022_07_18_Base TCV-MAP.xlsx", "For PATH with MAPs", skip=3)
mmgh_data_ns = read_xlsx("./data/2022_07_18_Base TCV-MAP.xlsx", "For PATH without MAPs", skip=3)

mmgh_data_ns$ISO[mmgh_data_ns$CountryName=="Dominica"] = "DMA" # there was an error with this before 22 May
mmgh_data_maps$ISO[mmgh_data_maps$CountryName=="Dominica"] = "DMA" # there was an error with this before 22 May

wbdata = read_xlsx("./data/GDP_LE_pop.xlsx", sheet="Data")
wbdata = wbdata[wbdata$`Series Name`=="Population, total", c("ISO", "Yr_2019","Yr_2020")]

# Bring in military personnel. Using numbers from 2005 here...
militarypop = read_xlsx("./data/MilitaryPersonnelPerCapita.xlsx", "MilitaryPersonnelPerCapita")
militarypop$Amount = militarypop$Amount*100 # the data was persons per 1000, multiply by 100 to get per 100K, which is my population
militarypop$ISO = countrycode::countrycode(militarypop$Country, "country.name", "iso3c")
militarypop$ISO[militarypop$Country=="Kosovo"] = "KSV"
militarypop$ISO[militarypop$Country=="Dominica"] = "DMA" # there was an error with this before 22 May
militarypop = militarypop[!is.na(militarypop$ISO),]
# Alternative: 0.4% of the population 16-45 in LIC and LMIC, and 0.6% in UMIC and HIC [from MMGH]
# from MMGH
militarypop = right_join(militarypop, data.frame(ISO=Countries$ISO, WB_group=Countries$WB_group))

# some countries missing...
# so take MMGH's calculation:
# in MMGHs calculation 40% of people are in the 16-45 age group, and 0.006 of them are in the military itself (if UMIC or HIC) or 0.004 in LIC and LMIC.
for(k in which(is.na(militarypop$Amount))){
  militarypop$Amount[k] = ifelse(militarypop$WB_group %in% c("Upper middle income", "High income"), 
                                 1e5*0.4*0.006, 
                                 1e5*0.4*0.004) 
} # isn't 0.4*0.006*1e5 just better?

I2toCtf = F # option for secondary infections to go to be eligible for chronic carriage: no

tmp_z <- strtoi(Sys.getenv("SLURM_ARRAY_TASK_ID", unset=Sys.getenv("LINE_NUMBER", unset="1")))
z = which(which(ISO$countryiso==indicators$ISO[tmp_z])==CN)

create_dir_if_necessary(paste0(subdir,  ISO[CN[z],2]))

#***************************************************************
## Step 1: Simulations and graphs ----
#***************************************************************

# if(tmp_z!=133){
  load(paste0("../out/", ISO[CN[z],2], "/fit_global.Rdata"))
# }else{
  # load(paste0("../out/", ISO[CN[59],2], "/fit_global.Rdata")) # Use JAMAICA estimates for St. Vincent & Grenadines
# }
  
  pars_unc = SIR_fit$optimfit$par # mvrnorm(1000, SIR_fit$optimfit$par, ginv(SIR_fit$optimfit$hessian))
  vprob = 0.85 # rbeta(1000, 49.27, 8.595)
  omegav = 1/(10*52) # runif(1000, 1/(20*52), 1/(10*52))
  
  parameters = SIR_fit$parameters
  data = SIR_fit$data
  
  m1 = plogis(SIR_fit$optimfit$par[3])
  m2 = plogis(SIR_fit$optimfit$par[4])
  
  parameters$r = parameters$rC = plogis(SIR_fit$optimfit$par[5])
  
  # R0_WQ = exp(SIR_fit$optimfit$par[1] + (100-data$haves)/100*SIR_fit$optimfit$par[2])
  R0_WQ = exp(plogis(SIR_fit$optimfit$par[1])*2.609-1 + (100-data$haves)/100*(plogis(SIR_fit$optimfit$par[2])*2.303))
  
  coeff = c(m1, m1, m2, rep(1, parameters$al-3)) 
  
  tmp_newpop = rep(parameters$py/parameters$num_wq, parameters$num_wq)
  parameters$newpop_wq = tmp_newpop
  parameters$pop_wq = rep(1/parameters$num_wq, parameters$num_wq)
  
  # derive the effective contact rate from R0.
  tmp_theta = sum(tmp_newpop/sum(parameters$py)*parameters$theta) # should be around 0.01542 for IND
  parameters$betap_wq = c(t(R0_WQ %o% coeff*(parameters$mub+parameters$delta)/(1+parameters$delta*parameters$rC*tmp_theta/parameters$mub)))

  # ***************************************************************
  # SIMS: N&S years (2023-2032) ----
  # ***************************************************************
  
  t1=Sys.time()
  tmp_vax = as.numeric(indicators[indicators$ISO==ISO[CN[z],2], paste0("mcv1_wq", 1:5)])/100
  SIR_wq_vax = vax_sims_nowater_wq(parameters, vprob=vprob, omegav=omegav,
                                     rcov = tmp_vax,
                                     ccov = pmax(tmp_vax-0.1, 0.01))
  t2=Sys.time()
  
  intro_sep = Countries$TCVMAP_adoption_year[Countries$ISO==ISO[CN[z],2]] - 
            Countries$NS_intro[Countries$ISO==ISO[CN[z],2]] # time between NS and MAPS
  intro_sep[intro_sep>19] = 19
  
  # keep just what I need.
  SIR_wq_vax$outvax0 = SIR_wq_vax$outvax0[200*52+intro_sep*52,] # No NS
  SIR_wq_vax$outvax1 = SIR_wq_vax$outvax1[200*52+intro_sep*52,] # NS with R
  SIR_wq_vax$outvax2 = SIR_wq_vax$outvax2[200*52+intro_sep*52,] # NS with R + C
    
  SIR_wq_vax_unc = SIR_wq_vax # lapply(1:100, function(i){SIR_wq_vax}) 
  
  save(SIR_wq_vax_unc, pars_unc, vprob, omegav, file=paste0(subdir,  ISO[CN[z],2], "/SIR_wq_vax_NOunc_global.Rdata"))

  # ***************************************************************
  # SIMS: MAPS years (2033-2052) ----
  # ***************************************************************
# for (z in c(56, 97, 94, 62, 90)){
#   load(paste0("../out/", ISO[CN[z],2], "/fit_global.Rdata"))
#   load(paste0(subdir, ISO[CN[z],2], "/SIR_wq_vax_NOunc_global.Rdata"))
  
  for (scenario_switch in 1:2){
    # make lists so I can get things to play nice with the uncertainty cea code
    SIR_wq_maps_unc0 = list()
    SIR_wq_maps_unc1 = list()
    SIR_wq_maps_unc2 = list()
    
    SIR_wq_vax_unc0 = list() 
    SIR_wq_vax_unc1 = list()
    SIR_wq_vax_unc2 = list()
  
    tmp_vax = as.numeric(indicators[indicators$ISO==ISO[CN[z],2], paste0("mcv1_wq", 1:5)])/100
    extracov = mean(1-tmp_vax)*maps_cov/2/0.2 # div by 0.2 to get coverage relative to the WQ1 and WQ2 pop
    # uncov rate * whole pop (5) * cov_htr / two wqs where we will put these people
    
    tmp_uc1 = mmgh_data_ns$uc1_del_2037[mmgh_data_ns$ISO==ISO[CN[z],2]] # c(0.65, 0.62, 0.85, 0.30)
    # names(tmp_uc1) = c("IND", "KEN", "NGA", "PAK")
    tmp_uc2 = mmgh_data_ns$uc2_del_2037[mmgh_data_ns$ISO==ISO[CN[z],2]] # c(0.35, 0.27, 0.10, 0.65)
    # names(tmp_uc2) = c("IND", "KEN", "NGA", "PAK")
    tmp_uc3 = mmgh_data_ns$uc3_del_2037[mmgh_data_ns$ISO==ISO[CN[z],2]] # c(0.0, 0.12, 0.05, 0.05)
    # names(tmp_uc3) = c("IND", "KEN", "NGA", "PAK")
    
    # Not the new military, apparently 10% coverage of the WHOLE military each year...
    tmp_uc5_tmp = militarypop$Amount[militarypop$ISO==ISO[CN[z],2]]/(SIR_wq_vax_unc$parameters$py[5])*mmgh_data_ns$uc5_cov_2037[mmgh_data_ns$ISO==ISO[CN[z],2]]
    tmp_uc5 = -log(1-tmp_uc5_tmp)/52
    # the first part sees what part of the population 15-24, which are the people who enter military, are part of the military, 
    # which had been calculated as a proportion of the whole population (if there was data) or as a proportion of the pop 16-45 (if using the MMGH estimate.)
    # the second part sees out of what part of the military population is covered per year.
    
    tmp_uc6 = -log(1-0.01*mmgh_data_ns$uc6_cov_2037[mmgh_data_ns$ISO==ISO[CN[z],2]]/(sum(SIR_wq_vax_unc$parameters$py[5:6])/1e5))/52
    # traveler population is 1% of the whole population, but I only care about 
    # the part that is over 15 (the only part that will make a formidable difference)
    # and that's about 69% of people in the LMICs, and then derive the rate by log:
    # 0.01*0.05/0.69 ===> rate per week = -log(1-0.01*0.05/0.69)/52
    
    # Scenarios were there WAS no previous TCV (NS) -----
    # as a proportion of the people who age in to the 15-25 pop (I take a tenth of the population 5-14)
    init_ode_maps_wq0 = SIR_wq_vax_unc$outvax0
    SIR_wq_vax0=vax_maps_sims_nowater_wq(parameters, init_ode_maps_wq0[2:dim(init_ode_maps_wq0)[2]], vprob=vprob, omegav=omegav,
                                           uc1cov = rep(0, parameters$num_wq),
                                           uc2cov = rep(0, parameters$num_wq),
                                           uc3cov = rep(0, parameters$num_wq), 
                                           uc4cov = rep(0, parameters$num_wq),
                                           uc5cov = rep(0, parameters$num_wq), 
                                           uc6cov = rep(0, parameters$num_wq))
    
    if (scenario_switch %in% c(1, 3)){
      SIR_wq_maps0=vax_maps_sims_nowater_wq(parameters, init_ode_maps_wq0[2:dim(init_ode_maps_wq0)[2]], vprob=vprob, omegav=omegav,
                                            uc1cov = c(tmp_vax[1]+extracov, tmp_vax[2]+extracov, tmp_vax[3:5])*tmp_uc1,
                                            uc2cov = c(tmp_vax[1]+extracov, tmp_vax[2]+extracov, tmp_vax[3:5])*tmp_uc2,
                                            uc3cov = c(tmp_vax[1]+extracov, tmp_vax[2]+extracov, tmp_vax[3:5])*tmp_uc3, 
                                            uc4cov = c(tmp_vax[1]+extracov, tmp_vax[2]+extracov, tmp_vax[3:5])-0.1,
                                            # assume, like with Joke, that the C in RC for initiation is 10 points lower than R. BUT, shouldn't that mean MAPS has more potential?
                                            uc5cov = rep(tmp_uc5, parameters$num_wq), 
                                            uc6cov = rep(tmp_uc6, parameters$num_wq))
    } else {
      SIR_wq_maps0=vax_maps_sims_nowater_wq(parameters, init_ode_maps_wq0[2:dim(init_ode_maps_wq0)[2]], vprob=vprob, omegav=omegav,
                                            uc1cov = tmp_vax*tmp_uc1,
                                            uc2cov = tmp_vax*tmp_uc2 + c(rep(extracov,2), rep(0,3))*tmp_uc2/(tmp_uc2+tmp_uc3),
                                            uc3cov = tmp_vax*tmp_uc3 + c(rep(extracov,2), rep(0,3))*tmp_uc3/(tmp_uc2+tmp_uc3), 
                                            uc4cov = c(tmp_vax[1]+extracov, tmp_vax[2]+extracov, tmp_vax[3:5])-0.1,
                                            # assume, like with Joke, that the C in RC for initiation is 10 points lower than R. BUT, shouldn't that mean MAPS has more potential?
                                            uc5cov = rep(tmp_uc5, parameters$num_wq), 
                                            uc6cov = rep(tmp_uc6, parameters$num_wq))    
    }
    
    # Scenarios were there WAS previous campaign (started NS with R only, no catch-up) -----
    init_ode_maps_wq1 = SIR_wq_vax_unc$outvax1
    SIR_wq_vax1=vax_maps_sims_nowater_wq(parameters, init_ode_maps_wq1[2:dim(init_ode_maps_wq1)[2]], vprob=vprob, omegav=omegav,
                                           uc1cov = tmp_vax*tmp_uc1,
                                           uc2cov = tmp_vax*tmp_uc2,
                                           uc3cov = tmp_vax*tmp_uc3, 
                                           uc4cov = rep(0, parameters$num_wq), # if there was previous TCV, then no catch-up
                                           uc5cov = rep(tmp_uc5, parameters$num_wq), 
                                           uc6cov = rep(tmp_uc6, parameters$num_wq))
    
    if (scenario_switch %in% c(1, 3)){
      SIR_wq_maps1=vax_maps_sims_nowater_wq(parameters, init_ode_maps_wq1[2:dim(init_ode_maps_wq1)[2]], vprob=vprob, omegav=omegav,
                                            uc1cov = c(tmp_vax[1]+extracov, tmp_vax[2]+extracov, tmp_vax[3:5])*tmp_uc1,
                                            uc2cov = c(tmp_vax[1]+extracov, tmp_vax[2]+extracov, tmp_vax[3:5])*tmp_uc2,
                                            uc3cov = c(tmp_vax[1]+extracov, tmp_vax[2]+extracov, tmp_vax[3:5])*tmp_uc3, 
                                            uc4cov = rep(0, parameters$num_wq), # if there was previous TCV, then no catch-up
                                            uc5cov = rep(tmp_uc5, parameters$num_wq), 
                                            uc6cov = rep(tmp_uc6, parameters$num_wq))
    } else {
      SIR_wq_maps1=vax_maps_sims_nowater_wq(parameters, init_ode_maps_wq1[2:dim(init_ode_maps_wq1)[2]], vprob=vprob, omegav=omegav,
                                            uc1cov = tmp_vax*tmp_uc1,
                                            uc2cov = tmp_vax*tmp_uc2 + c(rep(extracov,2), rep(0,3))*tmp_uc2/(tmp_uc2+tmp_uc3),
                                            uc3cov = tmp_vax*tmp_uc3 + c(rep(extracov,2), rep(0,3))*tmp_uc3/(tmp_uc2+tmp_uc3), 
                                            uc4cov = rep(0, parameters$num_wq), # if there was previous TCV, then no catch-up
                                            uc5cov = rep(tmp_uc5, parameters$num_wq),
                                            uc6cov = rep(tmp_uc6, parameters$num_wq))
    }
    
    # Scenarios where there WAS previous campaign (started NS with RC: routine + catch-up) -----
    init_ode_maps_wq2 = SIR_wq_vax_unc$outvax2
      SIR_wq_vax2=vax_maps_sims_nowater_wq(parameters, init_ode_maps_wq2[2:dim(init_ode_maps_wq0)[2]], vprob=vprob, omegav=omegav,
                                           uc1cov = tmp_vax*tmp_uc1,
                                           uc2cov = tmp_vax*tmp_uc2,
                                           uc3cov = tmp_vax*tmp_uc3, 
                                           uc4cov = rep(0, parameters$num_wq), # if there was previous TCV, then no catch-up
                                           uc5cov = rep(tmp_uc5, parameters$num_wq), 
                                           uc6cov = rep(tmp_uc6, parameters$num_wq))
    
    if (scenario_switch %in% c(1, 3)){
      SIR_wq_maps2=vax_maps_sims_nowater_wq(parameters, init_ode_maps_wq2[2:dim(init_ode_maps_wq0)[2]], vprob=vprob, omegav=omegav,
                                            uc1cov = c(tmp_vax[1]+extracov, tmp_vax[2]+extracov, tmp_vax[3:5])*tmp_uc1,
                                            uc2cov = c(tmp_vax[1]+extracov, tmp_vax[2]+extracov, tmp_vax[3:5])*tmp_uc2,
                                            uc3cov = c(tmp_vax[1]+extracov, tmp_vax[2]+extracov, tmp_vax[3:5])*tmp_uc3,
                                            uc4cov = rep(0, parameters$num_wq), # if there was previous TCV, then no catch-up
                                            uc5cov = rep(tmp_uc5, parameters$num_wq), 
                                            uc6cov = rep(tmp_uc6, parameters$num_wq))
    } else {
      SIR_wq_maps2=vax_maps_sims_nowater_wq(parameters, init_ode_maps_wq2[2:dim(init_ode_maps_wq0)[2]], vprob=vprob, omegav=omegav,
                                            uc1cov = tmp_vax*tmp_uc1,
                                            uc2cov = tmp_vax*tmp_uc2 + c(rep(extracov,2), rep(0,3))*tmp_uc2/(tmp_uc2+tmp_uc3),
                                            uc3cov = tmp_vax*tmp_uc3 + c(rep(extracov,2), rep(0,3))*tmp_uc3/(tmp_uc2+tmp_uc3), 
                                            uc4cov = rep(0, parameters$num_wq), # if there was previous TCV, then no catch-up
                                            uc5cov = rep(tmp_uc5, parameters$num_wq), 
                                            uc6cov = rep(tmp_uc6, parameters$num_wq))
      
    }
    
    # keep just what I need.
    SIR_wq_vax0$outvax_maps = NULL
    SIR_wq_vax1$outvax_maps = NULL
    SIR_wq_vax2$outvax_maps = NULL
    SIR_wq_maps0$outvax_maps = NULL
    SIR_wq_maps1$outvax_maps = NULL
    SIR_wq_maps2$outvax_maps = NULL
    
    SIR_wq_vax_unc0 = SIR_wq_vax0 
    SIR_wq_vax_unc1 = SIR_wq_vax1 
    SIR_wq_vax_unc2 = SIR_wq_vax2  
    SIR_wq_maps_unc0 = SIR_wq_maps0 
    SIR_wq_maps_unc1 = SIR_wq_maps1 
    SIR_wq_maps_unc2 = SIR_wq_maps2 
  
    # the right comparisons are usually if the ones with the same number at the end.
    
    save(SIR_wq_vax_unc0, SIR_wq_vax_unc1, SIR_wq_vax_unc2, # these are if the MAPS never take place, what happens after the year of MAPS intro
        SIR_wq_maps_unc0, SIR_wq_maps_unc1, SIR_wq_maps_unc2, # this is if the MAPS are introduced. 
        pars_unc, vprob, omegav, 
        file=paste0(subdir, ISO[CN[z],2], "/wq_maps_NOunc_global_comp", scenario_switch, ".Rdata"))
  
    if(scenario_switch==1){
    save(SIR_wq_vax_unc0, SIR_wq_vax_unc1, SIR_wq_vax_unc2,
         SIR_wq_maps_unc0, SIR_wq_maps_unc1, SIR_wq_maps_unc2, 
         pars_unc, vprob, omegav, 
         file=paste0(subdir, ISO[CN[z],2], "/wq_maps_NOunc_global_comp", 3, ".Rdata"))
    }
  }

