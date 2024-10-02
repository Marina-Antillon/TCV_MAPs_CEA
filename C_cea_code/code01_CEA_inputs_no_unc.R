#################### Cost and DALY model results ####################

# CEA objects (lists)
# pardata = list() # data from literature, usually means and standard deviations. Informs pars.
# pars = list() # distribution parameters for uncertain parameters
unc = list() # uncertain parameters (sampled values) (cfr, hospitalization rate, duration of disease, costs)
fix = list() # fixed numbers (discount rate, daly weight)
threshold = list() # CEA thresholds

# Labels, to annotate arrays and data.
lbl = list()
lbl$age = c('0-9m','9m-2y','2-4y','5-14y','15-24y', '25y+') 
lbl$strat = c("N&S", "N&S + MAPS")
lbl$par <- c("R0", "R0m", "m1", "m2", "rC", "reporting", "vprob", "omegav") # c(dimnames(vax_par)[[3]][1:6], "reporting")
lbl$wq = paste0("WQ", 1:5)
lbl$horizon <- c('10y', '20y')
lbl$discounting = c("disc", "no_disc")

#################################
## DATA, PARS, UNC: Mortality, health care use + hospitalization, 
## antimicrobial resistance, duration of illness, & severity + weights
## Hospitalization parameterization still not correct.
#################################
# data: 0·57 (0·58) [0·42–0·77]
# beta.parms.from.quantiles(c(0.42, 0.77))
# Based on Antillon et al (2017); based on estimate of relative incidence for passive vs active surveillance
unc$seek = 0.57 # rbeta(iterations, 17.18, 11.37)
unc$sev_nohc = 0.5 # rbeta(iterations, 50, 50) # assume the prob of sev among those that don't seek care is about half the prob of hosp
  # sanity check later: make sure this all lines up with assumption that either
  # 1/3 deaths are outside of care (Joke's paper) or only 10% of deaths (Anita). Better Joke so that we have a citation

## Mortality given hospitalization
tmp = c(0.01, 0.01, 0.062, 0.01, 0.01)
names(tmp) = c("Asia", "Mideast", "Africa", "Eurasia", "Americas")
unc$ip_mort = as.numeric(tmp[mmgh_data_ns_intro$continent[mmgh_data_ns_intro$ISO==ISO[CN[z],2]]]) # rbeta(iterations, 18.61, 1842) 
# from Zoe: rbeta(iterations, 22.7, 851.2)

## Multiplier for mortality for those that do not seek care 
## compared to those that are hospitalized
## (parameters for a uniform distribution)
## lower. 0-0.05 gave us the final mortality we were looking for, I think
# pars$ip_vs_nohc_mort_low = 0.8
# pars$ip_vs_nohc_mort_high = 0.8

unc$ip_vs_nohc_mort_or = 1.5 # runif(iterations, pars$ip_vs_nohc_mort_low, pars$ip_vs_nohc_mort_high)

# make sev the same proportion...
# cfr_out_as = (1-0.03)*(1-0.57)*0.1*unc$ip_vs_nohc_mort*0.01
# cfr_out_af = (1-0.36)*(1-0.57)*0.1*unc$ip_vs_nohc_mort*0.062
# 
# cfr_out_as_AMR = 0.03*(1-0.70)*0.18*unc$ip_vs_nohc_mort*0.01*2
# cfr_out_af_AMR = 0.36*(1-0.70)*0.18*unc$ip_vs_nohc_mort*0.062*2
# 
# cfr_in_as = (1-0.03)*0.57*0.2*((1-0.007)*0.01 + 0.007*0.048)
# cfr_in_af = (1-0.36)*0.57*0.2*((1-0.076)*0.062 + 0.076*0.197)
# 
# cfr_in_as_no_AMR = 0.03*0.70*0.33*((1-0.007)*0.01 + 0.007*0.048)*2
# cfr_in_af_no_AMR = 0.36*0.70*0.33*((1-0.076)*0.062 + 0.076*0.197)*2
# 
# all_mort_nohs_as = (cfr_out_as + cfr_out_as_AMR)+(cfr_in_as+cfr_in_as_no_AMR)
# all_mort_nohs_af = (cfr_out_af + cfr_out_af_AMR)+(cfr_in_af+cfr_in_af_no_AMR)
# 
# val_mort_nohs_as = (cfr_out_as + cfr_out_as_AMR)/all_mort_nohs_as
# 0.2196566
# val_mort_nohs_af = (cfr_out_af + cfr_out_af_AMR)/all_mort_nohs_af
# 0.1600015

# Mortality among typhoid intestinal perforation: From Marchello
# 15.5% (6.7-24.1%)
# beta.parms.from.quantiles(c(0.067,0.241))
# In Africa: Hmisc::binconf(387, 1967-387)
# In Asia: Hmisc::binconf(46, 999-46)

# tmp = c(0.04605, 0.04605, 0.1967, 0.1967, 0.1967)
# names(tmp) = c("IND", "PAK", "NGA", "KEN", "MWI")
# unc$ipip_mort = tmp[ISO[CN[z],2]]

tmp = c(0.04827, 0.04827, 0.1967, 0.01, 0.01)
names(tmp) = c("Asia", "Mideast", "Africa", "Eurasia", "Americas")
unc$ipip_mort = as.numeric(tmp[mmgh_data_ns_intro$continent[mmgh_data_ns_intro$ISO==ISO[CN[z],2]]]) # rbeta(iterations, 18.61, 1842) 
# there is no documentation in the literature of IP in Eurasia or the Americas

# Probability of IP complications among the hospitalized.
# Africa: Hmisc::binconf(37, 486)
# Asia: Hmisc::binconf(34, 4622)

tmp = c(0.007356, 0.007356, 0.07613, 0.007356, 0.007356)
names(tmp) = c("Asia", "Mideast", "Africa", "Eurasia", "Americas")
unc$ipip_comp = as.numeric(tmp[mmgh_data_ns_intro$continent[mmgh_data_ns_intro$ISO==ISO[CN[z],2]]]) # rbeta(iterations, 18.61, 1842) 

## Hospitalization
# should be normal-binomial

# beta.parms.from.quantiles(c(0.03, 0.12))
# data: 0.06 (0.03, 0.12) # don't remember where this comes from
# Linda's estimate: 0·004-0·249 - is this a prediction interval?
# Active community surveillance only
# From Maile's in Malawi 0.04 (95% PI: 0.01-0.11), which used 
# Linda's estimate as a Bayesian prior and the 8/105 
# hospitalizations in STRATAA to update the prior
# The numbers from STRATAA: MWI: 8/105, NPL: 5/150
# BGD: 10/359 [See Meiring 2021 paper]
# tmpmeta=meta::metaprop(c(8, 5, 10), c(105, 150, 359))
# meta::forest(tmpmeta)

# Including studies that Linda had:
# tmpmeta=meta::metaprop(c(8, 5, 10, 11, 13, 8, 21),
# c(105, 150, 359, 80, 94, 11, 147))
# RE estimate: 0.11 (0.05-0.25)
# beta.parms.from.quantiles(c(0.05, 0.25))

# WHERE ARE THESE LAST FOUR STUDIES FROM??
# tmpmeta=meta::metaprop(c(11, 13, 8, 21), c(80, 94, 11, 147))
# meta::forest(tmpmeta)
# RE estimate: 0.22 (0.09-0.45)
# beta.parms.from.quantiles(c(0.09, 0.45))

# pars$hosp_alpha = 5
# pars$hosp_beta = 15

# Strataa: metaprop(c(8, 5, 10), c(105, 150, 359)) 
# SEFI, Tier 1: 46/299
# SEAP: c(1295, 455, 1054), c(4873, 1602, 2230)
# tmpmeta = meta::metaprop(c(8, 5, 10, 46, 1295, 455, 1054), 
#                          c(105, 150, 359, 299, 4873, 1602, 2230),
#                          c("MWI", "BGD", "NPL", "IND", "BGD", "NPL", "PAK"))
# meta::forest(tmpmeta)
# recover 0.2: sum all TyVAC but leave all others alone...
# Try stratifying by country and see then if two levels of REs makes a difference
# Make at least one figure about drivers that are based on assumptions or 
# CFR, OOHC CFR, HOSP, AMR impact. Does it need to be CEA or could it just be DALYs

unc$hosp = 0.14 # rbeta(iterations, pars$hosp_alpha, pars$hosp_beta)

## Length of stay (hospitalization, in days)
## independent of AMR
# pardata$hosp_stay_mean = 6
# pardata$hosp_stay_err = 2
# 
# pars$hosp_stay_alpha = (pardata$hosp_stay_mean)^2/(pardata$hosp_stay_err^2) #SHAPE
# pars$hosp_stay_beta = (pardata$hosp_stay_mean)/(pardata$hosp_stay_err^2) #RATE
# 
# unc$hosp_stay = rgamma(iterations, pars$hosp_stay_alpha, pars$hosp_stay_beta)

## AMR probability
# pars$amr_prob_alpha = 6.08	
# pars$amr_prob_beta = 233.75

# tmp1 = c(0.03, 0.07, 0.36, 0.01, 0.02) 
# Europe was made up, Africa midpoint eyeballed, but subregion within Africa is important
# names(tmp1) = c("Asia", "Mideast", "Africa", "Eurasia", "Americas")

tmp2 = c(0.03, 0.68, 0.39, 0.78, 0.93, 0.27, 0.02, 0.7, 0.0, 0.01, 
         0.6, 0.64, 0.36, 0.59, 0.13, 0.9, 0.61)
names(tmp2) = c("IND", "PAK", "NGA", "KEN", "MWI", "BGD", "NPL", "KHM", "IDN", "BFA", 
                "NGA", "GHA", "COD", "TZA", "UZB", "TJK", "ZAF")

if(ISO[CN[z],2] %in% names(tmp2)){
  unc$amr_prob = as.numeric(tmp2[ISO[CN[z],2]]) # rbeta(iterations, 18.61, 1842) 
}else{
  unc$amr_prob = as.numeric(amr_data_unsd_region$AMR_prev[amr_data_unsd_region$ISO==ISO[CN[z],2]]) # rbeta(iterations, 18.61, 1842) 
}
## effect of AMR along the care-seeking and treatment cascade (all relative to AMS case)
## on cost, on duration of illness, on prob of death, 
## on reporting odds (which affects # of moderate vs severe cases), 
## and on hospitalization odds (which affects mortality)
## (uniform distribution)
## Make different vectors for each, because they are not necessarily correlated
## must decide: would cost be higher because you stay longer in the hospital, 
## or because of the hosp AND more meds. (latter)
# pars$amr_rel_burden_low = 2
# pars$amr_rel_burden_high = 2

unc$amr_rel_burden_doi = 2 # runif(iterations, pars$amr_rel_burden_low, pars$amr_rel_burden_high)
unc$amr_rel_burden_cost = 2 # runif(iterations, pars$amr_rel_burden_low, pars$amr_rel_burden_high)

unc$amr_rel_death_or = 2 # Zoe calculated this: 1.7 (95% CI, .69–4.33)
# runif(iterations, pars$amr_rel_burden_low, pars$amr_rel_burden_high)
unc$amr_rel_reporting_or = 2 # runif(iterations, pars$amr_rel_burden_low, pars$amr_rel_burden_high)
unc$amr_rel_hosp_or = 2 # runif(iterations, pars$amr_rel_burden_low, pars$amr_rel_burden_high)

## DALY weights (mild, moderate, severe) [beta distributed]

# pars$daly_wt_mod_alpha = beta.parms.from.quantiles(c(0.032, 0.074))$a
# pars$daly_wt_mod_beta = beta.parms.from.quantiles(c(0.032, 0.074))$b
# pars$daly_wt_sev_alpha = beta.parms.from.quantiles(c(0.088, 0.19))$a
# pars$daly_wt_sev_beta = beta.parms.from.quantiles(c(0.088, 0.19))$b
# pars$daly_wt_sevip_alpha = beta.parms.from.quantiles(c(0.22, 0.442))$a
# pars$daly_wt_sevip_beta = beta.parms.from.quantiles(c(0.22, 0.442))$b

# unc$daly_wt_mild = rbeta(iterations, pars$daly_wt_mild_alpha, pars$daly_wt_mild_beta)
unc$daly_wt_mod = 0.052 # rbeta(iterations, pars$daly_wt_mod_alpha, pars$daly_wt_mod_beta)
unc$daly_wt_sev = 0.21 # rbeta(iterations, pars$daly_wt_sev_alpha, pars$daly_wt_sev_beta)
unc$daly_wt_sevip = 0.32 # rbeta(iterations, pars$daly_wt_sevip_alpha, pars$daly_wt_sevip_beta)

## Duration of illness, regardless of whether they are hospitalized all those days or not
# pardata$doi_hcs_mean = 16 
# pardata$doi_hcs_err = 2
# # same for ip and op
# pars$doi_hcs_alpha = pardata$doi_hcs_mean^2/(pardata$doi_hcs_err^2)
# pars$doi_hcs_beta = pardata$doi_hcs_mean/(pardata$doi_hcs_err^2)

unc$doi_hcs = 16 # rgamma(iterations, pars$doi_hcs_alpha, rate=pars$doi_hcs_beta)

# multiplier for duration of illness (doi) for people 
# who do not seek care
# pars$doi_nohcs_mult_low = 0.5
# pars$doi_nohcs_mult_high = 1

unc$doi_nohcs_mult = 0.5 # runif(iterations, pars$doi_nohcs_mult_low, pars$doi_nohcs_mult_high)

# multiplier for duration of illness (doi) for people 
# who have intestinal perforations
# pars$doi_hosp_vs_ip_mult_low = 1.5
# pars$doi_hosp_vs_ip_mult_high = 2 #must look up, but this is what I vaguely remember

unc$doi_hosp_vs_ip_mult = 2 # runif(iterations, pars$doi_hosp_vs_ip_mult_low, pars$doi_hosp_vs_ip_mult_high)

#################################
## DATA, PARS, UNC: Vaccine costs
## Not finished: campaign costs for countries with data
#################################

# Vaccine delivery

# convertcost(0.31, 2017, 2021, "IND", "IND", "usd", "usd")
# convertcost(3.68, 2017, 2021, "IND", "IND", "usd", "usd")
# find_gamma_pars(0.33, 3.91)

# pars$vax_delivery_alpha = 2.9698 # SHAPE
# pars$vax_delivery_beta = 0.5448 # RATE
# unc$vax_delivery_total = # rgamma(iterations, pars$vax_delivery_alpha, rate=pars$vax_delivery_beta)

# unc$vax_del_ns_uc_base = c(0.4137, 0.3691, 0.74, 0.5, 0.30, 0.29) # before I had country specific data
# unc$vax_del_maps_uc_base = c(0.7857, 0.8252, 1.04, 0.9, 0.50, 0.49) # before I had country specific data

vtia = read_xlsx("./data/vtia_vax_del_MA.xlsx", "toR")
vtia$ISO = countrycode::countrycode(vtia$Country, "country.name", "iso3c")
# vtia$ISO[vtia$Country=="Kosovo"] = "KSV"
# vtia$ISO[vtia$Country=="Dominica"] = "DMA" # there was an error with this before 22 May

unc$vax_del_ns_uc_base = c(vtia$uc1_ns[vtia$ISO==ISO[CN[z],2]], vtia$uc2_ns[vtia$ISO==ISO[CN[z],2]], vtia$uc3_ns[vtia$ISO==ISO[CN[z],2]], 0, 0.30, 0.29)
unc$vax_del_maps_uc_base = c(vtia$uc1_base[vtia$ISO==ISO[CN[z],2]], vtia$uc2_base[vtia$ISO==ISO[CN[z],2]], vtia$uc3_base[vtia$ISO==ISO[CN[z],2]], vtia$uc4_base[vtia$ISO==ISO[CN[z],2]], 0.50, 0.49)
unc$vax_del_maps_uc_pess = c(vtia$uc1_pess[vtia$ISO==ISO[CN[z],2]], vtia$uc2_pess[vtia$ISO==ISO[CN[z],2]], vtia$uc3_pess[vtia$ISO==ISO[CN[z],2]], vtia$uc4_pess[vtia$ISO==ISO[CN[z],2]], 0.50, 0.49)
unc$vax_del_maps_uc_opti = c(vtia$uc1_optim[vtia$ISO==ISO[CN[z],2]], vtia$uc2_optim[vtia$ISO==ISO[CN[z],2]], vtia$uc3_optim[vtia$ISO==ISO[CN[z],2]], vtia$uc4_optim[vtia$ISO==ISO[CN[z],2]], 0.50, 0.49)

# didn't use it...
# if(all(is.na(unc$vax_del_ns_uc_base))==4){
# unc$vax_del_ns_uc_base = c(max(vtia$uc1_ns,na.rm=T), max(vtia$uc2_ns,na.rm=T), max(vtia$uc3_ns,na.rm=T), 0, 0.30, 0.29)*2
# } 
# if(all(is.na(unc$vax_del_maps_uc_base))==4){
#   unc$vax_del_maps_uc_base = c(max(vtia$uc1_base,na.rm=T), max(vtia$uc2_base,na.rm=T), max(vtia$uc3_base,na.rm=T), 0, 0.30, 0.29)*2
# } 
# if(all(is.na(unc$vax_del_maps_uc_pess))==4){
#   unc$vax_del_maps_uc_pess = c(max(vtia$uc1_pess,na.rm=T), max(vtia$uc2_pess,na.rm=T), max(vtia$uc3_pess,na.rm=T), 0, 0.30, 0.29)*2
# }
# if(all(is.na(unc$vax_del_maps_uc_opti))==4){
#   unc$vax_del_maps_uc_opti = c(max(vtia$uc1_pess,na.rm=T), max(vtia$uc2_pess,na.rm=T), max(vtia$uc3_pess,na.rm=T), 0, 0.30, 0.29)*2
# }
  
# Vaccine SUPPLIES
# pardata$vax_sup_lci = 0.21
# pardata$vax_sup_uci = 0.24
# pars$vax_sup_mean = 0.5*(pardata$vax_sup_lci+pardata$vax_sup_uci)
# pars$vax_sup_err = 0.25*(pardata$vax_sup_uci-pardata$vax_sup_lci)

# Appendix B in Portnoy et al Vaccine 2015 says that in 2010 USD
# injection and safety supplies were 0.19-0.22. After inflation
# Joke assigned the min-max to be 0.21-0.24, and she left a note
# that we should assume that these are the 95% CI 
# she left another note to explore in sensitivity analyses whether 
# a wider CI would make much of a difference.

# pars$vax_sup_alpha = pars$vax_sup_mean^2/(pars$vax_sup_err^2) #SHAPE
# pars$vax_sup_beta = pars$vax_sup_mean/(pars$vax_sup_err^2) #RATE

unc$vax_sup = 0.25 # rgamma(iterations, pars$vax_sup_alpha, rate=pars$vax_sup_beta)

#################################
## DATA, PARS, UNC: RX costs (for people who do not come into healthcare)
#################################

## Drugs costs - outpatient

tmp1 =  c(13.38, 10.13, 0.35, 10.13, 10.13) # Europe was made up, Africa midpoint eyeballed, but subregion within Africa is important
names(tmp1) = c("Asia", "Mideast", "Africa", "Eurasia", "Americas")

tmp2 = c(0.35, 8.60, 17.45, 14.11)
names(tmp2) = c("MWI", "NPL", "BGD", "PAK")

if(ISO[CN[z],2] %in% names(tmp2)){
  unc$op_costs_rx = as.numeric(tmp2[ISO[CN[z],2]]) # rbeta(iterations, 18.61, 1842) 
}else{
  unc$op_costs_rx = as.numeric(tmp1[mmgh_data_ns_intro$continent[mmgh_data_ns_intro$ISO==ISO[CN[z],2]]]) # rbeta(iterations, 18.61, 1842) 
}

#################################
## DATA, PARS, UNC: Outpatient (OP) costs
#################################
## consult only 

tmp1 =  c(47.15, 42.19, 32.26, 42.19, 42.19) 
names(tmp1) = c("Asia", "Mideast", "Africa", "Eurasia", "Americas")

tmp2 = c(60.36, 4.15, 20.26, 24.97, 126.97, 16.41)
names(tmp2) = c("MWI", "UGA", "NPL", "BGD", "PAK", "IND")

if(ISO[CN[z],2] %in% names(tmp2)){
  unc$op_costs = as.numeric(tmp2[ISO[CN[z],2]]) # rbeta(iterations, 18.61, 1842) 
}else{
  unc$op_costs = as.numeric(tmp1[mmgh_data_ns_intro$continent[mmgh_data_ns_intro$ISO==ISO[CN[z],2]]]) # rbeta(iterations, 18.61, 1842) 
}

#################################
## DATA, PARS, UNC: Inpatient (IP) costs
#################################
# not including drugs or surgery

tmp1 =  c(143.81, 171.48, 226.81, 171.48, 171.48) 
names(tmp1) = c("Asia", "Mideast", "Africa", "Eurasia", "Americas")

tmp2 = c(397.48, 56.14, 98.99, 153.76, 232.47, 90.04)
names(tmp2) = c("MWI", "UGA", "NPL", "BGD", "PAK", "IND")

if(ISO[CN[z],2] %in% names(tmp2)){
  unc$ip_costs = as.numeric(tmp2[ISO[CN[z],2]]) # rbeta(iterations, 18.61, 1842) 
}else{
  unc$ip_costs = as.numeric(tmp1[mmgh_data_ns_intro$continent[mmgh_data_ns_intro$ISO==ISO[CN[z],2]]]) # rbeta(iterations, 18.61, 1842) 
}

#################################
## DATA, PARS, UNC: Drugs costs - inpatient
#################################

tmp1 =  c(64.70, 62.61, 56.35, 62.61, 62.61) 
names(tmp1) = c("Asia", "Mideast", "Africa", "Eurasia", "Americas")

tmp2 = c(56.35, 51.16, 54.77, 88.17)
names(tmp2) = c("MWI", "NPL", "BGD", "PAK")

if(ISO[CN[z],2] %in% names(tmp2)){
  unc$ip_costs_rx = as.numeric(tmp2[ISO[CN[z],2]]) # rbeta(iterations, 18.61, 1842) 
}else{
  unc$ip_costs_rx = as.numeric(tmp1[mmgh_data_ns_intro$continent[mmgh_data_ns_intro$ISO==ISO[CN[z],2]]]) # rbeta(iterations, 18.61, 1842) 
}

#################################
## DATA, PARS, UNC: Inpatient (IP) costs with intestinal perforation
## stratify more: data$ip_costs$bedday, data$ip_costs$labs, etc?
## for pars, indicate dist of parameters
#################################

tmp1 =  c(191.27, 183.08, 246.83, 183.08, 183.08) 
names(tmp1) = c("Asia", "Mideast", "Africa", "Eurasia", "Americas")

tmp2 = c(343.32, 150.34, 301.60, 115.06, 43.47, 304.95)
names(tmp2) = c("NER", "UGA", "NPL", "BGD", "PAK", "IND")

if(ISO[CN[z],2] %in% names(tmp2)){
  unc$ipip_costs_surgery = as.numeric(tmp2[ISO[CN[z],2]]) # rbeta(iterations, 18.61, 1842) 
}else{
  unc$ipip_costs_surgery = as.numeric(tmp1[mmgh_data_ns_intro$continent[mmgh_data_ns_intro$ISO==ISO[CN[z],2]]]) # rbeta(iterations, 18.61, 1842) 
}

#################################
## Fixed quantities 
#################################
# Weights and discounts
## Vaccine wastage (maybe should be a distribution later?)
fix$vax_waste_ns = 0.11
fix$vax_waste_maps = 0.01

## discount rate for YLD, YLS, and costs, all equal for the Gates reference case
fix$disc = 0.03

# From WDI as well: http://data.worldbank.org/indicator/SP.DYN.LE00.IN/countries/VN?display=default
wbdata = read_xlsx("./data/GDP_LE_pop.xlsx", sheet="Data")

fix$lifexp = as.numeric(wbdata[wbdata$ISO==ISO[CN[z],2] & wbdata$`Series Name`=="Life expectancy at birth, total (years)","Yr_2020"])
fix$lifexp = ifelse(!is.na(fix$lifexp), fix$lifexp, as.numeric(wbdata[wbdata$ISO==ISO[CN[z],2] & wbdata$`Series Name`=="Life expectancy at birth, total (years)","Yr_2019"]))
fix$lifexp = ifelse(!is.na(fix$lifexp), fix$lifexp, 65)
#years lived per age group (midpoint of that age group), age groups not the same for all settings
# (<9m, 9m-2y, 2-4y, 5-14y, 15-25y, 25+y)

tmplived = c(4.5/12, (24+9)/(2*12), 2+3/2, 10, 20, 0.5*(25+fix$lifexp))
fix$yll = fix$lifexp - tmplived
names(fix$yll) = lbl$age

fix$pop100k = as.numeric(wbdata[wbdata$ISO==ISO[CN[z],2] & wbdata$`Series Name`=="Population, total", "Yr_2020"])/1e5
fix$pop100k = ifelse(!is.na(fix$pop100k), fix$pop100k, as.numeric(wbdata[wbdata$ISO==ISO[CN[z],2] & wbdata$`Series Name`=="Population, total", "Yr_2019"])/1e5)

#################################
## Thresholds 
#################################
# GDP per capita in 2015
# Source: http://data.worldbank.org/indicator/NY.GDP.PCAP.PP.CD
# Variable name: GDP per capita, PPP (current international $)

threshold$vce = as.numeric(wbdata[wbdata$ISO==ISO[CN[z],2] & wbdata$`Series Name`=="GDP per capita (current US$)", "Yr_2021"])
threshold$vce = ifelse(!is.na(threshold$vce), threshold$vce, as.numeric(wbdata[wbdata$ISO==ISO[CN[z],2] & wbdata$`Series Name`=="GDP per capita (current US$)", "Yr_2020"]))
threshold$vce = ifelse(!is.na(threshold$vce), threshold$vce, as.numeric(wbdata[wbdata$ISO==ISO[CN[z],2] & wbdata$`Series Name`=="GDP per capita (current US$)", "Yr_2019"]))

tmp_vce = c(1071, 15975, 701.71, 640, 644)
names(tmp_vce) = c("SSD", "VEN", "YEM", "PRK", "ERI")

threshold$vce = ifelse(!is.na(threshold$vce), threshold$vce, as.numeric(tmp_vce[ISO[CN[z],2]]))
threshold$ce = 3*threshold$vce

