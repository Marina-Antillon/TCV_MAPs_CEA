# Country in-depth analysis

# Packages ----

# options(bitmapType='cairo') # the cluster sometimes needs me to run this option.

library("Hmisc")
library('deSolve')
library("MASS")
library("matrixcalc")
library("epitools")

library("tidyr")
library("dplyr")
library("abind")

library("readxl")
library("R.matlab")

# Analysis settings -----
pars_sens=gsub("pars_sens=","",commandArgs()[grep("pars_sens=", commandArgs())]) # should be about pars_sens
maps_cov=gsub("maps_cov=","",commandArgs()[grep("maps_cov=", commandArgs())]) %>% as.numeric() # should be about maps_cov
reg_dist=gsub("reg_dist=","",commandArgs()[grep("reg_dist=", commandArgs())]) # should be about maps_cov

# reg_dist="dist"
# pars_sens="all" 
# maps_cov=0.2

# Read in functions -----
fls = list.files("./functions/", pattern="^[fcn]")
fls = fls[fls!="fcn_ggplot_settings.R"]
for (i in 1:length(fls)){source(paste0("./functions/", fls[i]))}

# Set up high-level subdirs -----
subdir = paste0("../out_deepdive_vax/", reg_dist, "/", pars_sens, "/", round(maps_cov, 1)) 
create_dir_if_necessary(subdir)

#***************************************************************
## Step 0: Bring in data ----
#***************************************************************

TransPar = readMat("./data/model_samples_ALL_Aug2021.mat")
CN = readMat("./data/burden_dynamic_link_Yale.mat")$cn
avgage_yale = readMat("./data/burden_dynamic_link_Yale.mat")$avgage.country
avgage_ihme = readMat("./data/burden_dynamic_link_IHME.mat")$avgage.country

agedist = readMat("./data/burden_dynamic_link_Yale.mat")$agedist.country

ISO = read.csv("./data/iso.csv")

## Real country populations. ----
# This helps calculate coverage for the military and travellers.
wbdata = read_xlsx("./data/GDP_LE_pop.xlsx", sheet="Data")
# tmp_wholepop= as.numeric(wbdata[wbdata$ISO==ISO[CN[z],2] & wbdata$`Series Name`=="Population, total", "Yr_2020"])

## Indicator data by subnational region ----
indicators = read_xlsx("./data/Deep_dive_indicators_2Aug.xlsx", paste0("data_to_r_", reg_dist))
# although is says Aug 2 is was updated in October.

if(reg_dist=="reg"){
  indicators_natl = indicators %>% filter(Region=="Total")
  indicators = indicators %>% filter(Region!="Total")
  indicators$RegionLbl = gsub(" ","", indicators$Region)
  indicators$RegionLbl = gsub("&","", indicators$RegionLbl)
} else {
  indicators_natl = indicators %>% filter(District=="Total")
  indicators = indicators %>% filter(District!="Total")
  indicators$RegionLbl = gsub(" ","", indicators$District)
  indicators$RegionLbl = gsub("&","", indicators$RegionLbl)
}

## Bring in UC proportions and MAPS intro from MMGH ----
mmgh_data_maps = read_xlsx("./data/2022_07_18_Base TCV-MAP.xlsx", "For PATH with MAPs", skip=3)
mmgh_data_ns = read_xlsx("./data/2022_07_18_Base TCV-MAP.xlsx", "For PATH without MAPs", skip=3)

## Bring in NS intro from MMGH ----
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

# Bring in military personnel. Using numbers from 2005 here...
militarypop = read_xlsx("./data/MilitaryPersonnelPerCapita.xlsx", "MilitaryPersonnelPerCapita")
militarypop$Prct = militarypop$Amount/100
militarypop$ISO = countrycode::countrycode(militarypop$Country, "country.name", "iso3c")
militarypop$ISO[militarypop$Country=="Kosovo"] = "KSV"
militarypop = militarypop[!is.na(militarypop$ISO),]
# Alternative: 0.4% of the population 16-45 in LIC and LMIC, and 0.6% in UMIC and HIC

#***************************************************************
# Step 1: Simulations and graphs ----
#***************************************************************

## Which thread is this? ----
zs <- strtoi(Sys.getenv("SLURM_ARRAY_TASK_ID", unset=Sys.getenv("LINE_NUMBER", unset="1")))
z=ISO$Var1[ISO$countryiso==indicators$ISO[zs]] #VERY necessary for the code that assigns costs according to the country.

## Low-level subdir -----
if(reg_dist=="reg"){
  tmp_subdir = paste0(subdir, "/", indicators$ISO[zs], "/", sprintf("%02d", zs), "_",indicators$RegionLbl[zs])
} else {
  tmp_subdir = paste0(subdir, "/", indicators$ISO[zs], "/", sprintf("%03d", zs), "_",indicators$RegionLbl[zs])
}

create_dir_if_necessary(tmp_subdir)

## Call in fits for that country ----
load(paste0("../out/", indicators$ISO[zs], "/fit_global.Rdata"))

pars_unc = SIR_fit$optimfit$par # mvrnorm(1000, SIR_fit$optimfit$par, ginv(SIR_fit$optimfit$hessian))
vprob = 0.85 # rbeta(1000, 49.27, 8.595)
omegav = 1/(10*52) # runif(1000, 1/(20*52), 1/(10*52))
parameters = SIR_fit$parameters
data = SIR_fit$data
if(pars_sens %in% c("all", "sanionly", "povsani")){
  data$haves = as.numeric(indicators[zs, paste0("sani_wq", 1:5)])
} else {
  data$haves = as.numeric(indicators_natl[indicators_natl$ISO==indicators$ISO[zs], paste0("sani_wq", 1:5)])
}
# no need to worry about mcv_cov yet because this is not about vax yet.
# mcv cov comes in in a few lines.
data$haves = pmax(pmin(data$haves, 0.99), 0.01)

## M1 and M2 parameters
m1 = plogis(SIR_fit$optimfit$par[3])
m2 = plogis(SIR_fit$optimfit$par[4])

## r and rC parameters
parameters$r = parameters$rC = plogis(SIR_fit$optimfit$par[5])

# R0_WQ = exp(SIR_fit$optimfit$par[1] + (100-data$haves*100)/100*SIR_fit$optimfit$par[2])
R0_WQ = exp(plogis(SIR_fit$optimfit$par[1])*2.609-1 + (100-data$haves*100)/100*(plogis(SIR_fit$optimfit$par[2])*2.303))

# from global: R0_WQ = exp(SIR_fit$optimfit$par[1] + (100-data$haves)/100*SIR_fit$optimfit$par[2])

# R0_WQ = pmax(pmin(R0_WQ, 25), 0.01)
coeff = c(m1, m1, m2, rep(1, parameters$al-3)) 

if(pars_sens %in% c("all", "povonly", "povsani", "povcov")){
  parameters$pop_wq = as.numeric(indicators[zs, paste0("pop_wq", 1:5)])
} else {
  parameters$pop_wq = rep(0.2, 5)
}

## Effective contact rate from R0. ----
# With DD countries, include poverty distribution
tmp_newpop = c(parameters$py %o% parameters$pop_wq)
parameters$newpop_wq = tmp_newpop
tmp_theta = sum(tmp_newpop/sum(parameters$py)*parameters$theta) # should be around 0.01542 for IND
parameters$betap_wq = c(
  t(R0_WQ %o% coeff*(parameters$mub+parameters$delta)/(1+parameters$delta*parameters$rC*tmp_theta/parameters$mub))
  )
parameters$R0_WQ = R0_WQ
parameters$R0 = sum(R0_WQ*parameters$pop_wq)
  
if(pars_sens %in% c("all", "covonly", "povcov")){
  tmp_vax = as.numeric(indicators[zs, paste0("mcv1_wq", 1:5)])
} else {
  tmp_vax = as.numeric(indicators_natl[indicators_natl$ISO==indicators$ISO[zs], paste0("mcv1_wq", 1:5)])
}
tmp_vax = pmax(pmin(tmp_vax, 0.99), 0.01)

# ***************************************************************
# SIMS: N&S years (2023-2032) ----
# ***************************************************************

SIR_wq_vax = vax_sims_nowater_wq(parameters, vprob=vprob, omegav=omegav,
                                 rcov = tmp_vax,
                                 ccov = pmax(tmp_vax-0.1, 0.01))

intro_sep = mmgh_data_ns_intro$TCVMAP_adoption_year[mmgh_data_ns_intro$ISO==ISO[CN[z],2]] - 
  mmgh_data_ns_intro$NS_intro[mmgh_data_ns_intro$ISO==ISO[CN[z],2]] # time between NS and MAPS
intro_sep[intro_sep>19] = 19

tmp = SIR_wq_vax$outvax0[,paste0("I1", 1:30)]
tmp = SIR_wq_vax$outvax0[,paste0("cumI1", 1:30)]

## Save N&S: keep just what I need & save----
SIR_wq_vax$outvax0 = SIR_wq_vax$outvax0[200*52+intro_sep*52,]
SIR_wq_vax$outvax1 = SIR_wq_vax$outvax1[200*52+intro_sep*52,]
SIR_wq_vax$outvax2 = SIR_wq_vax$outvax2[200*52+intro_sep*52,]

# colnames(SIR_wq_vax$outvax0)
# 
# tmp = SIR_wq_vax$outvax0[,2:181]
# tmp2 = array(as.numeric(tmp), dim=c(6,5,6)) 
# tmp2 %>% apply(2,"sum")

SIR_wq_vax_unc = SIR_wq_vax # lapply(1:100, function(i){SIR_wq_vax}) 

# SIR_wq_vax_unc$simcases_novac %>% apply(c(1,3), sum) %>% sweep(2, parameters$pop_wq, "/")

save(SIR_wq_vax_unc, pars_unc, vprob, omegav, parameters, 
     file=paste0(tmp_subdir, "/SIR_wq_vax_NOunc.Rdata"))

# ***************************************************************
# SIMS: MAPS years (2033-2052) ----
# ***************************************************************

  for (scenario_switch in 1:2){
    # make lists so I can get things to play nice with the uncertainty cea code
    SIR_wq_maps_unc0 = list()
    SIR_wq_maps_unc1 = list()
    SIR_wq_maps_unc2 = list()
    
    SIR_wq_vax_unc0 = list() 
    SIR_wq_vax_unc1 = list()
    SIR_wq_vax_unc2 = list()
    
    tmp_target_new_vax = sum((1-tmp_vax)*parameters$pop_wq)*maps_cov # /2/parameters$pop_wq[1]

    extracov=c()
    if(tmp_target_new_vax/2 < (1-tmp_vax[1])*parameters$pop_wq[1]){ # 
      extracov[1] = tmp_target_new_vax/2/max(parameters$pop_wq[1],0.0001)
    } else {
      extracov[1] = 1-tmp_vax[1]
    }
    
    if(tmp_target_new_vax/2 < (1-tmp_vax[2])*parameters$pop_wq[2]){ # 
      extracov[2] = tmp_target_new_vax/2/max(parameters$pop_wq[2],0.0001)
    } else {
      extracov[2] = 1-tmp_vax[2]
    }
    
    if(tmp_target_new_vax/2 > (1-tmp_vax[1])*parameters$pop_wq[1] & tmp_target_new_vax/2 > (1-tmp_vax[2])*parameters$pop_wq[2]){
      extracov[3:5] = 0
    } else {
      tmp1 = max(c(0, tmp_target_new_vax/2 - (1-tmp_vax[1])*parameters$pop_wq[1])) # max checks that it is not lower than zero
      tmp2 = max(c(0, tmp_target_new_vax/2 - (1-tmp_vax[2])*parameters$pop_wq[2])) # max checks that it is not lower than zero
      
      # then divide by the pop of 3 to see what the increase in coverage must be in terms of % of pop3
      extracov[3] = min(c((tmp1+tmp2)/max(parameters$pop_wq[3],0.0001), (1-tmp_vax[3]))) # if pop_wq[3] is 0, make sure it doesn't give nan
      extracov[4:5] = 0
    }
    
    extracov_uc4 = c()
    tmp_target_new_vax_uc4 = sum((1-tmp_vax)*parameters$pop_wq)*ifelse(maps_cov==0.2, 0.1, ifelse(maps_cov==0.1,0,0.2))
    if(tmp_target_new_vax_uc4/2 < (1-tmp_vax[1])*parameters$pop_wq[1]){
      extracov_uc4[1] = tmp_target_new_vax_uc4/2/parameters$pop_wq[1]
    } else {
      extracov_uc4[1] = 1-tmp_vax[1]
    }
    
    if(tmp_target_new_vax_uc4/2 < (1-tmp_vax[2])*parameters$pop_wq[2]){
      extracov_uc4[2] = tmp_target_new_vax_uc4/2/parameters$pop_wq[2]
    } else {
      extracov_uc4[2] = 1-tmp_vax[2]
    }
    
    if(tmp_target_new_vax_uc4/2 > (1-tmp_vax[1])*parameters$pop_wq[1] & tmp_target_new_vax_uc4/2 > (1-tmp_vax[2])*parameters$pop_wq[2]){
      extracov_uc4[3:5] = 0
    } else {
      tmp1 = max(c(0, tmp_target_new_vax_uc4/2 - (1-tmp_vax[1])*parameters$pop_wq[1])) # max checks that it is not lower than zero
      tmp2 = max(c(0, tmp_target_new_vax_uc4/2 - (1-tmp_vax[2])*parameters$pop_wq[2])) # max checks that it is not lower than zero
      
      # then divide by the pop of 3 to see what the increase in coverage must be in terms of % of pop3
      extracov_uc4[3] = min(c((tmp1+tmp2)/parameters$pop_wq[3], (1-tmp_vax[3])))
      extracov_uc4[4:5] = 0
    }
      
    tmp_uc1 = mmgh_data_ns$uc1_del_2037[mmgh_data_ns$ISO==indicators$ISO[zs]] 
    tmp_uc2 = mmgh_data_ns$uc2_del_2037[mmgh_data_ns$ISO==indicators$ISO[zs]] 
    tmp_uc3 = mmgh_data_ns$uc3_del_2037[mmgh_data_ns$ISO==indicators$ISO[zs]] 
    # according to MMGH, in LMICs, it is 0.4% of the population 16-45. 0.6% in UMIC and HIC.
    # 0.004*(SIR_wq_vax_unc$parameters$py[5] + 0.5*SIR_wq_vax_unc$parameters$py[6])
    # They assume 10% coverage of the WHOLE military every year, not just the new military... though perhaps 10% is an adequate est of the turnover.
    tmp_uc5 = militarypop$Prct[militarypop$ISO==indicators$ISO[zs]]*1e5/SIR_wq_vax_unc$parameters$py[5]*mmgh_data_ns$uc5_cov_2037[mmgh_data_ns$ISO==indicators$ISO[zs]]
    #as a proportion of the people who age in to the 15-25 pop (I take a tenth of the population 5-14)
    
    tmp_wholepop = as.numeric(wbdata[wbdata$ISO==indicators$ISO[zs] & wbdata$`Series Name`=="Population, total", "Yr_2020"])
    tmp_uc6 = mmgh_data_ns$uc6_tpop_2037[mmgh_data_ns$ISO==indicators$ISO[zs]]/tmp_wholepop*1e5/sum(SIR_wq_vax_unc$parameters$py[5:6])
    
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
                                            uc1cov = (tmp_vax+extracov)*tmp_uc1,
                                            uc2cov = (tmp_vax+extracov)*tmp_uc2,
                                            uc3cov = (tmp_vax+extracov)*tmp_uc3, 
                                            uc4cov = (tmp_vax+extracov_uc4)-0.1,
                                            uc5cov = rep(tmp_uc5, parameters$num_wq), 
                                            uc6cov = rep(0.0001466, parameters$num_wq))
    } else {
      SIR_wq_maps0=vax_maps_sims_nowater_wq(parameters, init_ode_maps_wq0[2:dim(init_ode_maps_wq0)[2]], vprob=vprob, omegav=omegav,
                                            uc1cov = tmp_vax*tmp_uc1,
                                            uc2cov = tmp_vax*tmp_uc2 + extracov*tmp_uc2/(tmp_uc2+tmp_uc3),
                                            uc3cov = tmp_vax*tmp_uc3 + extracov*tmp_uc3/(tmp_uc2+tmp_uc3), 
                                            uc4cov = (tmp_vax+extracov_uc4)-0.1,
                                            uc5cov = rep(tmp_uc5, parameters$num_wq), 
                                            uc6cov = rep(0.0001466, parameters$num_wq))    
    }
    
    init_ode_maps_wq1 = SIR_wq_vax_unc$outvax1
      SIR_wq_vax1=vax_maps_sims_nowater_wq(parameters, init_ode_maps_wq1[2:dim(init_ode_maps_wq1)[2]], vprob=vprob, omegav=omegav,
                                         uc1cov = tmp_vax*tmp_uc1,
                                         uc2cov = tmp_vax*tmp_uc2,
                                         uc3cov = tmp_vax*tmp_uc3, 
                                         uc4cov = rep(0, parameters$num_wq),
                                         uc5cov = rep(tmp_uc5, parameters$num_wq), 
                                         uc6cov = rep(0.0001466, parameters$num_wq))
    
    if (scenario_switch %in% c(1, 3)){
      SIR_wq_maps1=vax_maps_sims_nowater_wq(parameters, init_ode_maps_wq1[2:dim(init_ode_maps_wq1)[2]], vprob=vprob, omegav=omegav,
                                            uc1cov = (tmp_vax+extracov)*tmp_uc1,
                                            uc2cov = (tmp_vax+extracov)*tmp_uc2,
                                            uc3cov = (tmp_vax+extracov)*tmp_uc3, 
                                            uc4cov = rep(0, parameters$num_wq),
                                            uc5cov = rep(tmp_uc5, parameters$num_wq), 
                                            uc6cov = rep(0.0001466, parameters$num_wq))
    } else {
      SIR_wq_maps1=vax_maps_sims_nowater_wq(parameters, init_ode_maps_wq1[2:dim(init_ode_maps_wq1)[2]], vprob=vprob, omegav=omegav,
                                            uc1cov = tmp_vax*tmp_uc1,
                                            uc2cov = tmp_vax*tmp_uc2 + extracov*tmp_uc2/(tmp_uc2+tmp_uc3),
                                            uc3cov = tmp_vax*tmp_uc3 + extracov*tmp_uc3/(tmp_uc2+tmp_uc3), 
                                            uc4cov = rep(0, parameters$num_wq),
                                            uc5cov = rep(tmp_uc5, parameters$num_wq), 
                                            uc6cov = rep(0.0001466, parameters$num_wq))
    }
    
    init_ode_maps_wq2 = SIR_wq_vax_unc$outvax2
      SIR_wq_vax2=vax_maps_sims_nowater_wq(parameters, init_ode_maps_wq2[2:dim(init_ode_maps_wq0)[2]], vprob=vprob, omegav=omegav,
                                         uc1cov = tmp_vax*tmp_uc1,
                                         uc2cov = tmp_vax*tmp_uc2,
                                         uc3cov = tmp_vax*tmp_uc3, 
                                         uc4cov = rep(0, parameters$num_wq),
                                         uc5cov = rep(tmp_uc5, parameters$num_wq), 
                                         uc6cov = rep(0.0001466, parameters$num_wq))
    
    if (scenario_switch %in% c(1, 3)){
      SIR_wq_maps2=vax_maps_sims_nowater_wq(parameters, init_ode_maps_wq2[2:dim(init_ode_maps_wq0)[2]], vprob=vprob, omegav=omegav,
                                            uc1cov = (tmp_vax+extracov)*tmp_uc1,
                                            uc2cov = (tmp_vax+extracov)*tmp_uc2,
                                            uc3cov = (tmp_vax+extracov)*tmp_uc3, 
                                            uc4cov = rep(0, parameters$num_wq),
                                            uc5cov = rep(tmp_uc5, parameters$num_wq), 
                                            uc6cov = rep(0.0001466, parameters$num_wq))
    } else {
      SIR_wq_maps2=vax_maps_sims_nowater_wq(parameters, init_ode_maps_wq2[2:dim(init_ode_maps_wq0)[2]], vprob=vprob, omegav=omegav,
                                            uc1cov = tmp_vax*tmp_uc1,
                                            uc2cov = tmp_vax*tmp_uc2 + extracov*tmp_uc2/(tmp_uc2+tmp_uc3),
                                            uc3cov = tmp_vax*tmp_uc3 + extracov*tmp_uc3/(tmp_uc2+tmp_uc3), 
                                            uc4cov = rep(0, parameters$num_wq),
                                            uc5cov = rep(tmp_uc5, parameters$num_wq), 
                                            uc6cov = rep(0.0001466, parameters$num_wq))
    }
      # the previous way of doing comparator 2: 
      # uc1cov = tmp_vax*tmp_uc1,
      # uc2cov = c(tmp_vax[1]+extracov, tmp_vax[2]+extracov, tmp_vax[3:5])*tmp_uc2,
      # uc3cov = c(tmp_vax[1]+extracov, tmp_vax[2]+extracov, tmp_vax[3:5])*tmp_uc3,
    
# ***************************************************************
# Step 2: Save MAPS: keep just what I need & save----
# ***************************************************************
    
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
    
    save(SIR_wq_vax_unc0, SIR_wq_vax_unc1, SIR_wq_vax_unc2,
         SIR_wq_maps_unc0, SIR_wq_maps_unc1, SIR_wq_maps_unc2, 
         pars_unc, vprob, omegav, parameters, 
         file=paste0(tmp_subdir, "/wq_maps_NOunc_comp", scenario_switch, ".Rdata"))
    
    if(scenario_switch==1){
      save(SIR_wq_vax_unc0, SIR_wq_vax_unc1, SIR_wq_vax_unc2,
           SIR_wq_maps_unc0, SIR_wq_maps_unc1, SIR_wq_maps_unc2, 
           pars_unc, vprob, omegav, parameters, 
           file=paste0(tmp_subdir, "/wq_maps_NOunc_comp", 3, ".Rdata"))
    }
  }
  print(paste0("Done with ", indicators$ISO[zs], "/", sprintf("%02d", zs), "_",indicators$RegionLbl[zs]))
