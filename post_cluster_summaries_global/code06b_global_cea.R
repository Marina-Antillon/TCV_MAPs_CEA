#############################
# Making graphs for the manuscript
# National analysis
#############################

library(readxl)
library(dplyr)
library(magrittr)
library(tidyr)
library(ggplot2)
library(ggrepel)
library(ggpubr)
library(cowplot)
library(patchwork)
library(scales)
library(countrycode)
library(wbstats)
library(flextable)
library(stringr)

# load("../ISO_which_to_run_vaxsims.Rdata")
source("./C_cea_code/00b_load_packages.R")
source("./C_cea_code/00c_load_helper_functions.R")

# useful functions
dollarmk = function(x){return(ifelse(x==0,"$0",ifelse(x<1e6, paste0("$", x/1e3, "K"), paste0("$", x/1e6, "M"))))} 
countmk = function(x){return(ifelse(x==0,"0",ifelse(x<1e6, paste0(x/1e3, "K"), paste0(x/1e6, "M"))))} 
countmk2 = function(x){return(ifelse(x==0,"0",ifelse(x<1e5, round(x,0), ifelse(x<1e6, paste0(round(x/1e3,0), "K"), paste0(round(x/1e6,0), "M")))))} 
countmk3 = function(x){return(ifelse(x==0,"",ifelse(x<1e6, paste0(x/1e3, "K"), ifelse(x<1e9, paste0(x/1e6, "M"), paste0(x/1e9, "B")))))} 

TransPar = readMat("./data/model_samples_ALL_Aug2021.mat")
CN = readMat("./data/burden_dynamic_link_Yale.mat")$cn
avgage_yale = readMat("./data/burden_dynamic_link_Yale.mat")$avgage.country
avgage_ihme = readMat("./data/burden_dynamic_link_IHME.mat")$avgage.country

agedist = readMat("./data/burden_dynamic_link_Yale.mat")$agedist.country

ISO = read.csv("./data/iso.csv")

# bring in continent information
library(readxl)
countries = read_excel("./data/Countries.xlsx", "Country Data w intro date") %>% 
  rename(ISO = `ISO code`, `WB_group` = `WB Group (June 2021)`, 
         `Gavi_eligibility` = `Gavi eligibility`, `TCV_MAPS_intro_date` = `TCV-MAP adoption year`,
         `TCV_NS_intro_date` = `TCV-NS adoption year`, `Typhoid_incidence` = `Typhoid incidence rating`,
         `TCVMAP_archetype`= `TCV-MAP archetype`) %>% 
  mutate(TCV_NS_intro_date = ifelse(TCV_NS_intro_date==0, 2030, TCV_NS_intro_date)) %>% 
  dplyr::select(ISO, Country, WB_group, Gavi_eligibility, Typhoid_incidence, 
                TCV_NS_intro_date, TCV_MAPS_intro_date, TCVMAP_archetype, WHO) %>% 
  filter(WB_group!="High income") %>% 
  mutate(continent = ifelse(WHO %in% c("WPRO", "SEARO"), "Asia", 
                            ifelse(WHO=="EMRO", "Mideast",
                                   ifelse(WHO=="AFRO", "Africa", 
                                          ifelse(WHO=="EURO", "Eurasia", "Americas")))))

#Bring in coverage from DHS
indicators = read_xlsx("./data/DHS_MICS_GlobalAnalysis_updated.xlsx", "data_to_r")
indicators$ISO = countrycode::countrycode(indicators$Country, "country.name", "iso3c")
indicators$ISO[indicators$Country=="Kosovo"] = "KSV"
indicators$ISO[indicators$Country=="Dominica"] = "DMA" # there was an error with this before 22 May

no_imp_water = read.csv("./data/share-without-improved-water_OWID_JMP.csv") %>% 
  filter(Year==2020)

no_imp_sanitation = read.csv("./data/sanitation_service_level_JMP.csv") %>% 
  filter(Year==2020, Service.level=="Unimproved")

indicators = right_join(indicators, data.frame(ISO=countries$ISO, WB_group=countries$WB_group, continent=countries$continent)) %>% 
  filter(WB_group!="High income") %>% 
  mutate(across(contains("sani_wq"), ~as.numeric(.x))) %>% 
  mutate(across(contains("water_wq"), ~as.numeric(.x))) %>% 
  mutate(across(contains("mcv1_wq"), ~as.numeric(.x))) %>% 
  mutate_at(c("mcv1_total", "sani_total", "water_total"), as.numeric)

for(k in which(is.na(indicators$sani_total))){
  indicators$sani_total[k] = ifelse(indicators$ISO[k] %in% no_imp_sanitation$ISO3,
                                    100-no_imp_sanitation$Coverage[no_imp_sanitation$ISO3==indicators$ISO[k]],
                                    99.9)
}
indicators$sani_total[indicators$sani_total>99.9] = 99.9

indicators$sani_dif = indicators$sani_wq5-indicators$sani_wq1
indicators$sani_difa = indicators$sani_total-indicators$sani_wq1
indicators$sani_difb = indicators$sani_wq5-indicators$sani_total

indicators = indicators %>% 
  mutate(sani_wq1 = ifelse(is.na(sani_wq1), pmax(sani_total-20, 0.1), sani_wq1),
         sani_wq2 = ifelse(is.na(sani_wq2), pmax(sani_total-10, 0.1), sani_wq2),
         sani_wq3 = ifelse(is.na(sani_wq3), pmin(sani_total, 99.9), sani_wq3),
         sani_wq4 = ifelse(is.na(sani_wq4), pmin(sani_total+10, 99.9), sani_wq4),
         sani_wq5 = ifelse(is.na(sani_wq5), pmin(sani_total+20, 99.9), sani_wq5)) %>% 
  mutate(across(contains("sani_wq"), ~pmax(pmin(.x,99.9),0.1))) # HERE: need to bottom code the lower ones to 0.1

# For those missing vaccination statistics, of mcv1 estimates from before 2010, let's replace it with 
# the estimate from WHO coverage estimates in 2021. Or whatever year MMMGH used (2020? 2019?)
# these countries we won't be able to do a ECEA, just a CEA, but that's ok.
# assume the difference between lowest to highest is 5%
mcv1_WHO = read.csv("./data/mcv1_cov_WHO.csv")

for(k in which(is.na(indicators$mcv1_total))){
  indicators$mcv1_total[k] = ifelse(indicators$ISO[k] %in% mcv1_WHO$SpatialDimValueCode,
                                    mcv1_WHO$FactValueNumeric[mcv1_WHO$SpatialDimValueCode==indicators$ISO[k] & mcv1_WHO$Period==2021],
                                    78.3) # the average
}

# indicators$mcv1_dif = indicators$mcv1_wq5-indicators$mcv1_wq1
# indicators$mcv1_difa = indicators$mcv1_total-indicators$mcv1_wq1
# indicators$mcv1_difb = indicators$mcv1_wq5-indicators$mcv1_total

indicators = indicators %>% 
  mutate(mcv1_wq1 = ifelse(is.na(mcv1_wq1), pmax(mcv1_total-4, 0.1), mcv1_wq1),
         mcv1_wq2 = ifelse(is.na(mcv1_wq2), pmax(mcv1_total-2, 0.1), mcv1_wq2),
         mcv1_wq3 = ifelse(is.na(mcv1_wq3), pmin(mcv1_total, 99.9), mcv1_wq3),
         mcv1_wq4 = ifelse(is.na(mcv1_wq4), pmin(mcv1_total+2, 99.9), mcv1_wq4),
         mcv1_wq5 = ifelse(is.na(mcv1_wq5), pmin(mcv1_total+4, 99.9), mcv1_wq5)) %>% 
  mutate(across(contains("mcv1_wq"), ~pmax(pmin(.x,99.9),0.1)))

wbdata = read_xlsx("./data/GDP_LE_pop_UN_DHS_combined.xlsx", sheet="Data")

gdp_cap = wb_data(indicator = "NY.GDP.PCAP.CD", 
                  country = countries$ISO, start_date = 2010, end_date=2023) %>% 
  pivot_wider(id_cols = iso3c, names_from = date, values_from = NY.GDP.PCAP.CD) %>%
  rename(ISO = iso3c) %>% 
  mutate(GDP_cap = `2022`) %>% 
  mutate(GDP_cap = ifelse(!is.na(GDP_cap), GDP_cap, `2021`)) %>% 
  mutate(GDP_cap = ifelse(!is.na(GDP_cap), GDP_cap, `2020`)) %>% 
  mutate(GDP_cap = ifelse(!is.na(GDP_cap), GDP_cap, `2015`)) %>% 
  mutate(GDP_cap = ifelse(ISO=="VEN", 1501, GDP_cap)) %>%
  mutate(GDP_cap = ifelse(!is.na(GDP_cap), GDP_cap, 650)) %>% 
  dplyr::select(ISO, GDP_cap)

# Bring in UC proportions and MAPS intro from MMGH
mmgh_data_maps = read_xlsx("./data/2022_07_18_Base TCV-MAP.xlsx", "For PATH with MAPs", skip=3)
mmgh_data_ns = read_xlsx("./data/2022_07_18_Base TCV-MAP.xlsx", "For PATH without MAPs", skip=3)

mmgh_data_ns$ISO[mmgh_data_ns$CountryName=="Dominica"] = "DMA" # there was an error with this before 22 May
mmgh_data_maps$ISO[mmgh_data_maps$CountryName=="Dominica"] = "DMA" # there was an error with this before 22 May

# Create necessary directories
create_dir_if_necessary("./maps_tcv_global/")

# Move to the right directory for outputs
setwd("./maps_tcv_global/")

# The Loop --
# for(i in c(2,2.25,3,4.5)){

# make table of options... then use that to call all results. 
# then bindrows.

sens = data.frame(vax_sens = rep(c(2, 2.25, 3, 4.5), times=5), 
          maps_sens = c(rep(c(0.2, 0.1, 0.3), each=4), rep(0.2, 8)), 
          mcv_sens = c(rep("mcv1current", 12), rep("mcv1imp", 4), rep("mcv2", 4)))
# MCV sens analysis is only done assumming 20% coverage, but with all prices.
  
icer_results_all = list()
icer_results_wq = list()
# Calling data ------

for(i in 1:20){
  
  tmp_icer_results_all = list()
  tmp_icer_results_wq = list()
  
  for(zz in 1:(length(indicators$ISO))){
    
    tmp_path = paste0('../out_global_cea-given in/', sens$maps_sens[i], ifelse(sens$mcv_sens[i]=="mcv1current", "", '_'), 
                      ifelse(sens$mcv_sens[i]=="mcv1current", "", sens$mcv_sens[i]), '/', sens$vax_sens[i], 'dollars/')
    
    load(paste0(tmp_path, indicators$ISO[zz], '/trt_inputs.Rdata'))
    load(paste0(tmp_path, indicators$ISO[zz], '/epi_inputs.Rdata'))
    load(paste0(tmp_path, indicators$ISO[zz], '/code03_treatment.Rdata')) # for the death rate per case
    load(paste0(tmp_path, indicators$ISO[zz], '/code05_cea.Rdata'))
    
    tmp = which(ISO[CN,2]==indicators$ISO[zz])
    # tmp_woods = wtp_woods$usd_2013_high[wtp_woods$ISO==goodruns[zz]]
    
    if(indicators$ISO[zz]=="JAM"){threshold$vce=5183.58}
    if(indicators$ISO[zz]=="MDA"){threshold$vce=5230.66}
    if(indicators$ISO[zz]=="MKD"){threshold$vce=6694.64}
    
    tmp_pop = as.numeric(wbdata[wbdata$ISO==ISO[CN[z],2] & wbdata$`Series Name`=="Population, total", "Yr_2021"])
    tmp_pop = ifelse(!is.na(tmp_pop), tmp_pop, as.numeric(wbdata[wbdata$ISO==ISO[CN[z],2] & wbdata$`Series Name`=="Population, total", "Yr_2020"]))
    tmp_pop = ifelse(!is.na(tmp_pop), tmp_pop, as.numeric(wbdata[wbdata$ISO==ISO[CN[z],2] & wbdata$`Series Name`=="Population, total", "Yr_2019"]))
    
    tmp_icer_results_all[[zz]] = icers_dalys_costs_all_summary %>% 
      mutate(gdp=gdp_cap$GDP_cap[gdp_cap$ISO==indicators$ISO[zz]], ISO=indicators$ISO[zz], 
             Country=countries$Country[countries$ISO==indicators$ISO[zz]], 
             WB_group=countries$WB_group[countries$ISO==indicators$ISO[zz]], 
             Country_archetype=countries$TCVMAP_archetype[countries$ISO==indicators$ISO[zz]], 
             Gavi_eligibility = countries$Gavi_eligibility[countries$ISO==indicators$ISO[zz]], 
             Continent=countries$continent[countries$ISO==indicators$ISO[zz]],
             TCV_NS_intro_date = countries$TCV_NS_intro_date[countries$ISO==indicators$ISO[zz]],
             TCV_MAPS_intro_date = countries$TCV_MAPS_intro_date[countries$ISO==indicators$ISO[zz]],
             lifexp=fix$lifexp, 
             pop=tmp_pop, # fix$pop100k*1e5,
             meaninc = mean(TransPar$incsamples[,tmp]), 
             meanage = 0.5*(avgage_yale[tmp]+avgage_ihme[tmp]),
             death_per_case=avg_deaths,
             mcv1_cov = indicators$mcv1_total[indicators$ISO==indicators$ISO[zz]],
             mcv1_cov_wq1 = indicators$mcv1_wq1[zz], 
             sani_cov = indicators$sani_total[zz],
             sani_cov_wq1 = indicators$sani_wq1[zz],
             infants_vax_fixedpost=mmgh_data_maps$uc1_del_2037[mmgh_data_maps$ISO==indicators$ISO[zz]],
             SIR_fit_rep = epipar["reporting"]) %>% 
      bind_cols(data.frame(t(unlist(unc)))) %>% 
      mutate(vax_cost_sens = sens$vax_sens[i],
             mcv_cov_sens = sens$mcv_sens[i], 
             maps_cov_sens = sens$maps_sens[i])
  
    tmp_icer_results_wq[[zz]] = icers_dalys_costs_summary %>% 
      mutate(gdp=gdp_cap$GDP_cap[gdp_cap$ISO==indicators$ISO[zz]], ISO=indicators$ISO[zz], 
             Country = countries$Country[countries$ISO==indicators$ISO[zz]], 
             WB_group = countries$WB_group[countries$ISO==indicators$ISO[zz]], 
             Country_archetype = countries$TCVMAP_archetype[countries$ISO==indicators$ISO[zz]], 
             Gavi_eligibility = countries$Gavi_eligibility[countries$ISO==indicators$ISO[zz]], 
             Continent=countries$continent[countries$ISO==indicators$ISO[zz]],
             TCV_NS_intro_date = countries$TCV_NS_intro_date[countries$ISO==indicators$ISO[zz]],
             TCV_MAPS_intro_date = countries$TCV_MAPS_intro_date[countries$ISO==indicators$ISO[zz]],
             lifexp = fix$lifexp, 
             pop=tmp_pop, 
             pop_wq=tmp_pop/5,
             meaninc = mean(TransPar$incsamples[,tmp]), 
             meanage = 0.5*(avgage_yale[tmp]+avgage_ihme[tmp]),
             death_per_case=avg_deaths,
             mcv1_cov = indicators$mcv1_total[indicators$ISO==indicators$ISO[zz]],
             mcv1_cov_wq1 = indicators$mcv1_wq1[zz], 
             sani_cov = indicators$sani_total[zz],
             sani_cov_wq1 = indicators$sani_wq1[zz],
             infants_vax_fixedpost=mmgh_data_maps$uc1_del_2037[mmgh_data_maps$ISO==indicators$ISO[zz]],
             SIR_fit_rep = epipar["reporting"]) %>% 
      bind_cols(data.frame(t(unlist(unc)))) %>% 
      mutate(vax_cost_sens = sens$vax_sens[i],
             mcv_cov_sens = sens$mcv_sens[i], 
             maps_cov_sens = sens$maps_sens[i])
  }
  tmp_icer_results_all = dplyr::bind_rows(tmp_icer_results_all) %>% 
    mutate(CEcat = ifelse(icers<0, "CS",
                          ifelse(icers/gdp<0.5, "HCE",
                                 ifelse(icers/gdp<1, "VCE", 
                                        ifelse(icers/gdp<3, "CE", "Not CE"))))) %>% 
    mutate(CEcat = factor(CEcat, levels=c("CS", "HCE", "VCE", "CE", "Not CE")))
  tmp_icer_results_wq = dplyr::bind_rows(tmp_icer_results_wq) %>% 
    mutate(CEcat = ifelse(icers<0, "CS", 
                          ifelse(icers/gdp<0.5, "HCE",
                                 ifelse(icers/gdp<1, "VCE", ifelse(icers/gdp<3, "CE", "Not CE"))))) %>% 
    mutate(CEcat = factor(CEcat, levels=c("CS", "HCE", "VCE", "CE", "Not CE")))
  
  icer_results_all[[i]] = tmp_icer_results_all
  icer_results_wq[[i]] = tmp_icer_results_wq
}

# Country TCV-MAP archetypes 
#   1: HIC/UMIC with low typhoid incidence and/or AMR; 
#   2: LMIC/LIC in rest of the world with high typhoid incidence and/or AMR; 
#   3: LMIC/LIC in the rest of the world with medium typhoid incidence and/or AMR; 
#   4: LMIC/LIC in the rest of the world with low typhoid incidence and/or AMR; 
#   5: LMIC/LIC with a significant private market (countries in SEAR or WPR) with medium or high typhoid incidence and/or AMR)

icer_results_all = dplyr::bind_rows(icer_results_all) %>% 
  mutate(Country_archetype=factor(Country_archetype, levels=as.character(1:5), 
                       labels=c("1) UMIC", "2) LMIC/LIC: high inc.", 
                                "3) LMIC/LIC: med inc", "4) LMIC/LIC: low inc",
                                "5) LMIC/LIC: sig pr. market"))) %>% 
  mutate(WB_group=ifelse(WB_group=="Lower middle income", "Lower-middle income", 
                         ifelse(WB_group=="Upper middle income", "Upper-middle income", WB_group))) %>% 
  mutate(Gavi_eligibility = factor(Gavi_eligibility, 
                                 levels=c("Initial self-financing", "Preparatory transition phase", 
                                          "Fully self-financing", "Not eligible"))) %>% 
  mutate(comp = factor(comp, levels=paste0("comparator", 1:3), 
                       labels=c("80% switch to MAPs", "Targeted MAPs introduction", "100% switch to MAPs"))) %>% 
  mutate(vax_cost_sens = factor(vax_cost_sens, levels=c(2, 2.25, 3, 4.5), 
                       labels=c("$2.00 per dose", "$2.25 per dose", 
                                "$3.00 per dose", "$4.50 per dose"))) 
  
icer_results_wq = dplyr::bind_rows(icer_results_wq) %>% 
  mutate(Country_archetype=factor(Country_archetype, levels=as.character(1:5), 
                                  labels=c("1) UMIC", "2) LMIC/LIC: high inc.", 
                                           "3) LMIC/LIC: med inc", "4) LMIC/LIC: low inc",
                                           "5) LMIC/LIC: sig pr. market"))) %>% 
  mutate(WB_group=ifelse(WB_group=="Lower middle income", "Lower-middle income", 
                         ifelse(WB_group=="Upper middle income", "Upper-middle income", WB_group))) %>% 
  mutate(Gavi_eligibility = factor(Gavi_eligibility, 
                                   levels=c("Initial self-financing", "Preparatory transition phase", 
                                            "Fully self-financing", "Not eligible"))) %>% 
  mutate(comp = factor(comp, levels=paste0("comparator", 1:3), 
                       labels=c("80% switch to MAPs", "Targeted MAPs introduction", "100% switch to MAPs"))) %>% 
  mutate(vax_cost_sens = factor(vax_cost_sens, levels=c(2, 2.25, 3, 4.5), 
                                labels=c("$2.00 per dose", "$2.25 per dose", 
                                         "$3.00 per dose", "$4.50 per dose"))) 

save(icer_results_all, icer_results_wq, file="./figures/EpiBaseline/icer_results_data.Rdata")

# Make Excel summary -----
# Should I make this a table... give it to them as an excel. Costs as columns rather than rows...
TotalsExcel = list()

TotalsExcel[["no_disc"]] = icer_results_all %>% filter(horizon=="20y", discounting=="no_disc") %>% 
  group_by(ns_strat, maps_profile, comp, vax_cost_sens, mcv_cov_sens, maps_cov_sens) %>% 
  dplyr::summarise(TotalCasesAverted = sum(`Cases_averted`/1e5*pop, na.rm=T),
                   TotalDeathsAverted = sum(`Deaths_averted`/1e5*pop, na.rm=T),
                   TotalDALYsAverted = sum(`DALYs_averted`/1e5*pop, na.rm=T),
                   TotalCostDif = sum(`Cost_dif`/1e5*pop, na.rm=T)) 

TotalsExcel[["disc"]] = icer_results_all %>% filter(horizon=="20y", discounting=="disc") %>% 
  group_by(ns_strat, maps_profile, comp, vax_cost_sens, mcv_cov_sens, maps_cov_sens) %>% 
  dplyr::summarise(TotalCasesAverted = sum(`Cases_averted`/1e5*pop, na.rm=T),
                   TotalDeathsAverted = sum(`Deaths_averted`/1e5*pop, na.rm=T),
                   TotalDALYsAverted = sum(`DALYs_averted`/1e5*pop, na.rm=T),
                   TotalCostDif = sum(`Cost_dif`/1e5*pop, na.rm=T)) 

TotalsExcel[["README"]] = data.frame(variable_note = c(colnames(TotalsExcel$no_disc), paste0("Note ", 1:3)), 
                    note = c("The strategy used to deploy the TCV needle and syringe between 2023-2032. This has substantial consequences for the projected incidence of typhoid when MAPS becomes available. None, deployment in Routine only at 9months old, and Routine + Campaign (default).", 
                             "The MAPS profiles are Base (default, 20 cm3, 70s admin time), Pessimistic (20 cm3, 5min admin time), Optimistic (5 cm3, 15s admin time). This has implications for delivery costs.",
                             "Comparators: 1) 80% market penetration and use in UC1-6 (default), 2) 100% market penetration but use in UC2-6 only , 3) 100% market penetration in UC1-6.",
                             "Sensitivity analysis of vaccine cost (per dose): $2, $2.25, $3 (default), $4.50.",
                             "N&S coverage assumptions, index to the MCV1 proxy (default), indexed according to improvements in MCV1 coverage, and indexed according to coverage in MCV2 (in case typhoid is administered at 15m instead of 9m).",
                             "MAPS coverage among the unvaccinated: 20% (default), 10%, 30%.",
                             "Cases averted compared to the scenario where only NS is available, unless ns_strat = None, in which case it is compared to no vaccination at all",
                             "Deaths averted compared to the scenario where only NS is available, unless ns_strat = None, in which case it is compared to no vaccination at all",
                             "DALYs averted compared to the scenario where only NS is available, unless ns_strat = None, in which case it is compared to no vaccination at all",
                             "Cost differences compared to the scenario where only NS is available, unless ns_strat = None, in which case it is compared to no vaccination at all",
                             "Notice one sheet is for results with disconuting and one is without discounting. Only DALYs and costs are discounted, not cases and deaths.",
                             "This all assumes a horizon of 20 years.",
                             "For the MCV coverage assumptions, I have run the simulations assuming maps_cov_sens of 20% only, to avoid too many outputs. It can be done for other values, however."))

TotalsExcel = TotalsExcel[c("README", "no_disc", "disc")]
openxlsx::write.xlsx(TotalsExcel, "./figures/EpiBaseline/TotalsExcel.xlsx", overwrite=T)

## IMPACT -----
# For section 3.1.1 of Milestone 4

grand_summaries_fcn = function(data){
  # consider making this more succint:
  # summarise_at # https://dplyr.tidyverse.org/reference/summarise_all.html
  # the one above may not work if because I also want to rename
  # AND:
  # mutate(across(v1:v2, ~ .x + n))
 tmp =  data %>%
    dplyr::summarise(Cases_NS = sum(`Status Quo_Cases`/1e5*pop, na.rm=T),
                     Deaths_NS = sum(`Status Quo_Deaths`/1e5*pop, na.rm=T),
                     DALYs_NS = sum(`Status Quo_DALYs`/1e5*pop, na.rm=T),
                     Cost_NS = sum(`Status Quo_Cost`/1e5*pop, na.rm=T), 
                     Cases_MAPS = sum(`MAPS add_Cases`/1e5*pop, na.rm=T),
                     Deaths_MAPS = sum(`MAPS add_Deaths`/1e5*pop, na.rm=T),
                     DALYs_MAPS = sum(`MAPS add_DALYs`/1e5*pop, na.rm=T),
                     Cost_MAPS = sum(`MAPS add_Cost`/1e5*pop, na.rm=T),
                     Cases_Dif = sum(`Cases_averted`/1e5*pop, na.rm=T),
                     Deaths_Dif = sum(`Deaths_averted`/1e5*pop, na.rm=T),
                     DALYs_Dif = sum(`DALYs_averted`/1e5*pop, na.rm=T),
                     Cost_Dif = sum(`Cost_dif`/1e5*pop, na.rm=T)) %>%
    mutate(Cases_NS = format(Cases_NS, big.mark=","),
           Deaths_NS = format(Deaths_NS, big.mark=","),
           DALYs_NS = format(DALYs_NS, big.mark=","),
           Cost_NS = format(Cost_NS, big.mark=","),
           Cases_MAPS = format(Cases_MAPS, big.mark=","),
           Deaths_MAPS = format(Deaths_MAPS, big.mark=","),
           DALYs_MAPS = format(DALYs_MAPS, big.mark=","),
           Cost_MAPS = format(Cost_MAPS, big.mark=","),
           Cases_Dif = format(Cases_Dif, big.mark=","),
           Deaths_Dif = format(Deaths_Dif, big.mark=","),
           DALYs_Dif = format(DALYs_Dif, big.mark=","),
           Cost_Dif = format(Cost_Dif, big.mark=",")) %>% 
    pivot_longer(cols=Cases_NS:Cost_Dif,
                 names_sep = "_",
                 names_to = c("Outcome", "Type"),
                 values_to = "value")  %>% 
   mutate(value = trimws(value))
  return(tmp)
    }

table2summaries_continent = icer_results_all %>% 
                    filter(ns_strat=="RoutineCampaign", horizon=="20y", discounting=="no_disc", 
                            maps_profile=="Base", comp=="80% switch to MAPs",
                            vax_cost_sens=="$3.00 per dose", 
                            mcv_cov_sens=="mcv1current", maps_cov_sens==0.2) %>% 
  group_by(Continent) %>%
  grand_summaries_fcn %>% 
  pivot_wider(id_cols = c(Continent, Outcome), 
              names_from = Type,
              values_from = value)

table2summaries_total = icer_results_all %>% 
                    filter(ns_strat=="RoutineCampaign", horizon=="20y", discounting=="no_disc", 
                            maps_profile=="Base", comp=="80% switch to MAPs",
                            vax_cost_sens=="$3.00 per dose", 
                            mcv_cov_sens=="mcv1current", maps_cov_sens==0.2) %>% 
  grand_summaries_fcn %>% 
  pivot_wider(id_cols = Outcome, 
              names_from = Type,
              values_from = value) %>%
  mutate(Continent="Total")

table2summaries = bind_rows(table2summaries_total, table2summaries_continent)
table2summaries = table2summaries[,c(5,1:4)] 

ft = flextable(table2summaries) %>%
  set_header_labels(`Continent` = "Continent", `NS` = "N&S", `MAPS` = "MAPs", `Dif` = "Difference") %>% 
  merge_v(j=c("Continent")) %>% 
  theme_vanilla() %>% 
  valign(j = "Continent", valign = "center") %>% 
  width(j="Continent", 1.25) %>% width(j="Outcome", 1.25) %>% width(j="NS", 1.25) %>% width(j="MAPS", 1.25) %>% width(j="Dif", 1.25) %>% 
  align(j=c("NS", "MAPS", "Dif"), align="right") # %>% 
  # append_chunks(i = c(1, 7), j = 1, as_chunk(paste0("\n(No TCV = ", Cases_totals$Totals[c(2,1)], " cases)"))) %>% 

save_as_docx(`Epidemiology outputs` = ft, path = "./figures/DefaultAssumptionsTotals.docx")


# delta Cases, DALYs, costs in a graph for the default options by WQ. Put global totals too.
# name place with gdp per capita...

icer_results_all %>% filter(ns_strat=="RoutineCampaign", horizon=="20y", discounting=="no_disc", 
                            maps_profile=="Base", comp=="80% switch to MAPs",
                            vax_cost_sens=="$3.00 per dose", 
                            mcv_cov_sens=="mcv1current", maps_cov_sens==0.2) %>% 
  dplyr::summarise(`Status Quo_IntCost` = sum(`Status Quo_IntCost`/1e5*pop, na.rm=T),
                   `Status Quo_TrtCost` = sum(`Status Quo_TrtCost`/1e5*pop, na.rm=T),
                   `MAPS add_IntCost` = sum(`MAPS add_IntCost`/1e5*pop, na.rm=T),
                   `MAPS add_TrtCost` = sum(`MAPS add_TrtCost`/1e5*pop, na.rm=T)) %>% c() %>% format(big.mark=",")

# Status Quo_IntCost Status Quo_TrtCost   MAPS add_IntCost   MAPS add_TrtCost 
# "6,893,858,657"   "11,966,851,139"   "10,717,541,043"   "11,684,760,830" 

tmp = icer_results_all %>% filter(ns_strat=="RoutineCampaign", horizon=="20y", discounting=="no_disc", 
                            maps_profile=="Base", comp=="80% switch to MAPs",
                            vax_cost_sens=="$3.00 per dose", 
                            mcv_cov_sens=="mcv1current", maps_cov_sens==0.2) %>% 
  dplyr::mutate(`Status Quo_Cost` = `Status Quo_Cost`/1e5*pop,
                `MAPS add_Cost` = `MAPS add_Cost`/1e5*pop,
                  `Status Quo_IntCost` = `Status Quo_IntCost`/1e5*pop,
                   `Status Quo_TrtCost` = `Status Quo_TrtCost`/1e5*pop,
                   `MAPS add_IntCost` = `MAPS add_IntCost`/1e5*pop,
                   `MAPS add_TrtCost` = `MAPS add_TrtCost`/1e5*pop) %>% 
  dplyr::select(ISO, `Status Quo_IntCost`, `Status Quo_TrtCost`, `MAPS add_IntCost`, `MAPS add_TrtCost`,
         `Status Quo_Cost`, `MAPS add_Cost`) %>% 
    dplyr::mutate(sanity_SQ = `Status Quo_IntCost`+`Status Quo_TrtCost`-`Status Quo_Cost`, 
                  ratio_StatusQuo = `Status Quo_IntCost`/`Status Quo_Cost`,
                  ratio_MAPS = `MAPS add_IntCost`/`MAPS add_Cost`,
                  diff_IntCost =  `MAPS add_IntCost` - `Status Quo_IntCost`,
                  diff_Cost_ratio = (`MAPS add_IntCost` - `Status Quo_IntCost`)/(`Status Quo_TrtCost`-`MAPS add_TrtCost`))
# The increase in Intervention costs is much bigger than the decrease in Treatment costs.

# Baseline outcomes by WQ ----

ByCountryExcel = list()

ByCountryExcel[["WQ_no_disc"]] = icer_results_wq %>% filter(ns_strat=="RoutineCampaign", horizon=="20y", discounting=="no_disc", 
                                           mcv_cov_sens=="mcv1current", maps_cov_sens==0.2) %>% 
  dplyr::select("ISO", "Country", "WB_group", "Gavi_eligibility", "Continent", "Country_archetype",
                "gdp", "pop", "wealth_quintile", "maps_profile", "vax_cost_sens", "comp",
                "Cases_averted", "Deaths_averted", "DALYs_averted", "Cost_dif", "icers", "CEcat") %>% 
  mutate(across(c("Cases_averted", "DALYs_averted", "Deaths_averted", "Cost_dif"), ~(.x/1e5*pop))) %>% 
  mutate(icer_scaled_gdp = icers/gdp) 

ByCountryExcel[["WQ_disc"]] = icer_results_wq %>% filter(ns_strat=="RoutineCampaign", horizon=="20y", discounting=="no_disc", 
                                                            mcv_cov_sens=="mcv1current", maps_cov_sens==0.2) %>% 
  dplyr::select("ISO", "Country", "WB_group", "Gavi_eligibility", "Continent", "Country_archetype",
                "gdp", "pop", "wealth_quintile", "maps_profile", "vax_cost_sens", "comp",
                "Cases_averted", "Deaths_averted", "DALYs_averted", "Cost_dif", "icers", "CEcat") %>% 
  mutate(across(c("Cases_averted", "DALYs_averted", "Deaths_averted", "Cost_dif"), ~(.x/1e5*pop))) %>% 
  mutate(icer_scaled_gdp = icers/gdp) 

ByCountryExcel[["Country_no_disc"]] = icer_results_all %>% 
  filter(ns_strat=="RoutineCampaign", horizon=="20y", discounting=="no_disc", 
         mcv_cov_sens=="mcv1current", maps_cov_sens==0.2) %>% 
dplyr::select("ISO", "Country", "WB_group", "Gavi_eligibility", "Continent", "Country_archetype",
              "gdp", "pop", "maps_profile", "vax_cost_sens", "comp",
              "Cases_averted", "Deaths_averted","DALYs_averted", "Cost_dif", "icers", "CEcat") %>% 
  mutate(across(c("Cases_averted", "DALYs_averted", "Deaths_averted", "Cost_dif"), ~(.x/1e5*pop))) %>% 
  mutate(icer_scaled_gdp = icers/gdp) 

ByCountryExcel[["Country_disc"]] = icer_results_all %>% 
  filter(ns_strat=="RoutineCampaign", horizon=="20y", discounting=="no_disc", 
         mcv_cov_sens=="mcv1current", maps_cov_sens==0.2) %>% 
  dplyr::select("ISO", "Country", "WB_group", "Gavi_eligibility", "Continent", "Country_archetype",
                "gdp", "pop", "maps_profile", "vax_cost_sens", "comp",
                "Cases_averted", "Deaths_averted", "DALYs_averted", "Cost_dif", "icers", "CEcat") %>% 
  mutate(across(c("Cases_averted", "DALYs_averted", "Deaths_averted", "Cost_dif"), ~(.x/1e5*pop))) %>% 
  mutate(icer_scaled_gdp = icers/gdp) 

ByCountryExcel[["README"]] = data.frame(variable_note = c(colnames(ByCountryExcel$WQ_no_disc), paste0("Note ", 1:6)), 
                                     note = c("ISO3 code for the country.",
                                              "Country name.",
                                              "World bank country group.",
                                              "Gavi eligibility as of 2021.",
                                              "Continent.",
                                              "Country archetype according to MMGH classification.",
                                              "GDP per capita for the country for 2021",
                                              "Population size of the country.",
                                              "Wealth quintiles. WQ1 is the poorest and WQ5 is the wealthiest.",
                                              "The MAPS profiles are Base (default, 20 cm3, 70s admin time), Pessimistic (20 cm3, 5min admin time), Optimistic (5 cm3, 15s admin time). This has implications for delivery costs.",
                                              "Vaccine price (per dose): $2, $2.25, $3 (default), $4.50.",
                                              "Comparators: 1) 80% market penetration and use in UC1-6 (default), 2) 100% market penetration but use in UC2-6 only , 3) 100% market penetration in UC1-6.",
                                              "Cases averted compared to the scenario where only NS is available.",
                                              "Deaths averted compared to the scenario where only NS is available.",
                                              "DALYs averted compared to the scenario where only NS is available.",
                                              "Cost differences compared to the scenario where only NS is available.",
                                              "Incremental cost-effectiveness ratios in USD per DALY averted.",
                                              "Cost-effectiveness category. CS=cost-effective, HCE=highly cost-effective, VCE=very cost-effective, CE=cost-effective, Not CE=not cost-effective.",
                                              "Incremental cost-effectiveness ratios scaled by GDP per capita.",
                                              "The first two sheets give results per country and wealth quintile -- both not discounted and discounted by 3% per year-- and the second two sheets give the results for the whole country -- also not discounted and discounted by 3% per year.",
                                              "Notice one sheet is for results with disconuting and one is without discounting. Only DALYs and costs are discounted, not cases or deaths.",
                                              "These results are under the assumptions of a horizon of 20 years.",
                                              "For the MCV coverage assumptions, I have run the simulations assuming maps_cov_sens of 20% only, to avoid too many outputs. It can be done for other values, however.",
                                              "NS coverage assumptions, index to the MCV1 proxy (default), indexed according to improvements in MCV1 coverage, and indexed according to coverage in MCV2 (in case typhoid is administered at 15m instead of 9m).",
                                              "MAPS coverage is fixed at 20% (default) in these results. Separately, results for 10%, 30% could be made available."
                                              ))

ByCountryExcel = ByCountryExcel[c("README", "WQ_no_disc", "WQ_disc", "Country_no_disc", "Country_disc")]
openxlsx::write.xlsx(ByCountryExcel, "./figures/EpiBaseline/ByCountryExcel.xlsx", overwrite=T)

totals_wq_def = icer_results_wq %>% filter(ns_strat=="RoutineCampaign", horizon=="20y", discounting=="no_disc", 
                           maps_profile=="Base", comp=="80% switch to MAPs",
                           vax_cost_sens=="$3.00 per dose", 
                           mcv_cov_sens=="mcv1current", maps_cov_sens==0.2) %>% 
  dplyr::select("ISO", "Country", "WB_group", "Gavi_eligibility", "Continent", "Country_archetype",
                "gdp", "pop", "wealth_quintile",
                "Cases_averted", "DALYs_averted", "Deaths_averted", "Cost_dif") %>% 
  mutate(across(c("Cases_averted", "DALYs_averted", "Deaths_averted", "Cost_dif"), ~(.x/1e5*pop))) %>% 
  pivot_longer(c("Cases_averted", "DALYs_averted", "Deaths_averted", "Cost_dif"), 
               names_to="Measure", values_to = "value") %>% 
  mutate(Measure = factor(Measure, levels=c("Cases_averted", "Deaths_averted", "DALYs_averted", "Cost_dif"),
                          labels=c("Cases Averted", "Deaths Averted", "DALYs Averted", "Cost Difference (USD)"))) %>% 
  mutate(name_place = paste0(Country, " - ", ISO, " ($", round(gdp,0), ")"))
  
iso_order = totals_wq_def$name_place[order(totals_wq_def$gdp)] %>% unique
totals_wq_def = totals_wq_def %>% mutate(name_place = factor(name_place, levels=iso_order)) 

ggplot((totals_wq_def), aes(x=name_place, y=value, fill=wealth_quintile)) + 
  geom_bar(stat="identity", position = position_stack(reverse = TRUE)) + 
  scale_fill_manual(values = c('#e66101','#fdb863','#ffffbf','#b2abd2','#5e3c99'),
                    labels = c("Lowest 20%", "Second 20%", "Third 20%", "Fourth 20%", "Highest 20%")) + 
  themebar + theme(legend.title = element_text(size = 16, face="bold"), 
                   legend.text = element_text(size = 16, face = "bold", lineheight=0.8, margin = margin(r = 0.3, unit = 'cm'))) + 
  coord_flip() + scale_x_discrete(limits=rev) + 
  scale_y_continuous(labels=countmk) +
  xlab("") + ylab("") + labs(fill="Wealth Quintile") + 
  facet_grid(WB_group~Measure, scales = "free", space="free_y", switch = "x") + 
  theme(strip.text.x.bottom = element_text(angle = 0, size=16, hjust=0.5), 
        strip.text.y = element_text(angle = 270, size=16, hjust=0.5), 
        axis.text.y = element_text(angle = 0, size=10), 
        axis.text.x = element_text(angle = 0, size=14), 
        strip.placement = "outside", 
        panel.spacing = unit(0.25,'lines'))

ggsave("./figures/EpiBaseline/BaselineCasesDALYsDeathsCostsWQ.pdf", width = 16, height = 20, units="in")

ggplot((totals_wq_def), aes(x=name_place, y=value, fill=wealth_quintile)) + 
  geom_bar(stat="identity", position = position_fill(reverse = TRUE)) + 
  scale_fill_manual(values = c('#e66101','#fdb863','#ffffbf','#b2abd2','#5e3c99'),
                    labels = c("Lowest 20%", "Second 20%", "Third 20%", "Fourth 20%", "Highest 20%")) + 
  themebar + theme(legend.title = element_text(size = 16, face="bold"), 
                   legend.text = element_text(size = 16, face = "bold", lineheight=0.8, margin = margin(r = 0.3, unit = 'cm'))) + 
  coord_flip() + scale_x_discrete(limits=rev) + 
  xlab("") + ylab("") + labs(fill="Wealth Quintile") + 
  facet_grid(WB_group~Measure, scales = "free", space="free_y", switch = "x") + 
  theme(strip.text.x.bottom = element_text(angle = 0, size=16, hjust=0.5), 
        strip.text.y = element_text(angle = 270, size=16, hjust=0.5), 
        axis.text.y = element_text(angle = 0, size=10), 
        axis.text.x = element_text(angle = 0, size=14), 
        strip.placement = "outside", 
        panel.spacing = unit(0.25,'lines'))

ggsave("./figures/EpiBaseline/BaselineCasesDALYsDeathsCostsWQ_Dist.pdf", width = 16, height = 20, units="in")

totals_wq_def_gavi = totals_wq_def %>% 
  group_by(Gavi_eligibility, Measure, wealth_quintile) %>% 
  dplyr::summarise(sumvalue = sum(value, na.rm=T)) %>% 
  rename(datcat=Gavi_eligibility)

totals_wq_def_WB_group = totals_wq_def %>% 
  group_by(WB_group, Measure, wealth_quintile) %>% 
  dplyr::summarise(sumvalue = sum(value, na.rm=T)) %>% 
  rename(datcat=WB_group)

totals_wq_def_archetype = totals_wq_def %>% 
  group_by(Country_archetype, Measure, wealth_quintile) %>% 
  dplyr::summarise(sumvalue = sum(value, na.rm=T)) %>% 
  rename(datcat=Country_archetype) 

totals_wq_continent = totals_wq_def %>% 
  group_by(Continent, Measure, wealth_quintile) %>% 
  dplyr::summarise(sumvalue = sum(value, na.rm=T)) %>% 
  rename(datcat=Continent)

summary_plot = function(dat){ggplot(dat, aes(x=datcat, y=sumvalue, fill=wealth_quintile)) + 
  geom_bar(stat="identity", position = position_stack(reverse = TRUE)) + 
  scale_fill_manual(values = c('#e66101','#fdb863','#ffffbf','#b2abd2','#5e3c99'),
                    labels = c("Lowest 20%", "Second 20%", "Third 20%", "Fourth 20%", "Highest 20%")) + 
  themebar + theme(legend.title = element_text(size = 14, face="bold"), 
                   legend.text = element_text(size = 14, face = "bold", lineheight=0.8, margin = margin(r = 0.3, unit = 'cm')),
                   plot.margin = unit(c(0, 15, 0, 2.5), "pt")) + 
    theme(strip.text.x= element_text(angle = 0, size=14, hjust=0.5),
          axis.text.y = element_text(angle = 0, size=14), 
          axis.text.x = element_text(angle = 0, size=12), 
          strip.placement = "outside", 
          panel.spacing = unit(0.5,'lines')) + 
  scale_x_discrete(limits=rev) + coord_flip() + 
  scale_y_continuous(labels=countmk3, expand = c(0, 0)) +
  xlab("") + ylab("") + labs(fill="Wealth Quintile") + 
  facet_grid(.~Measure, scale="free")}

ggsum1 = summary_plot(totals_wq_continent) + theme(legend.position = "none")
ggsum2 = summary_plot(totals_wq_def_WB_group) + theme(legend.position = "none") 
ggsum3 = summary_plot(totals_wq_def_gavi)

# ggsum4 = summary_plot(totals_wq_def_archetype)

### pull together --

ggarrange(ggsum1, ggsum2, ggsum3, 
          labels="AUTO",
          ncol = 1, heights = c(0.85, 0.65, 1), 
          common.legend = F, align="v") 
# list("Gavi eligibility", "World Bank income group", "Country archetype", "WHO region")
ggsave("./figures/EpiBaseline/BaselineCasesDALYsDeathsCostsWQ_Summaries.pdf", width = 12, height = 8, units="in")

totals_wq_def %>% filter(Measure!="DALYs Averted") %>%
  group_by(Continent, Measure) %>% 
  dplyr::summarise(sumvalue = sum(value, na.rm=T)) %>% 
  rename(datcat=Continent) %>% 
  ggplot(aes(x=datcat, y=sumvalue)) + 
  geom_bar(stat="identity", position = position_stack(reverse = TRUE), fill="royalblue2") + 
  themebar + theme(legend.title = element_text(size = 14, face="bold"), 
                   legend.text = element_text(size = 14, face = "bold", lineheight=0.8, margin = margin(r = 0.3, unit = 'cm')),
                   plot.margin = unit(c(0, 15, 0, 2.5), "pt")) + 
  theme(strip.text.x= element_text(angle = 0, size=14, hjust=0.5),
        axis.text.y = element_text(angle = 0, size=10), 
        axis.text.x = element_text(angle = 0, size=12), 
        strip.placement = "outside", 
        panel.spacing = unit(0.5,'lines')) + 
  scale_x_discrete(limits=rev) + coord_flip() + 
  scale_y_continuous(labels=countmk3, expand = c(0, 0)) +
  xlab("") + ylab("") + labs(fill="Wealth Quintile") + 
  facet_grid(.~Measure, scale="free")
ggsave("./figures/EpiBaseline/BaselineCasesDALYsDeathsCosts_Summaries.pdf", width = 9, height = 3, units="in")

totals_wq_def %>% 
  group_by(Continent, Measure) %>% 
  dplyr::summarise(sumvalue = sum(value, na.rm=T)) %>% 
  rename(datcat=Continent) %>% 
  pivot_wider(id_cols=datcat, names_from = Measure, values_from = sumvalue) 
  
# Refit, with 2021 pops
# datcat   `Cases Averted` `Deaths Averted` `DALYs Averted` `Cost Difference`
# 1 Africa          2431102.          41494.         2099090.        885121390.
# 2 Americas          89394.            114.            6998.        261986423.
# 3 Asia            1540765.           2504.          159873.       1843299471.
# 4 Eurasia           13915.             33.8           2020.        180204857.
# 5 EMRO            1097993.           2854.          163292.        370979936.

# CEs -----

# spot check a couple of countries to see what happens with Comp2. 
# maybe make Appendix table or something? For 5 countries?

SpotCheckComp2 = icer_results_all %>% 
  filter(ns_strat=="None", horizon=="20y", discounting=="disc", 
         maps_profile=="Base", # comp=="80% switch to MAPS",
         vax_cost_sens=="$3.00 per dose",  
         mcv_cov_sens=="mcv1current", maps_cov_sens==0.2, 
         ISO %in% c("BFA", "KEN", "IND", "MWI", "NPL", "NGA")) 

# for this, go fetch the Nigeria results by hand
  # tmp_doses = doses %>% apply(c(3:5), "sum")
  # tmp_doses = tmp_doses %>% cubelyr::as.tbl_cube() %>% as_tibble()
  # tmp_doses = doses %>% apply(c(3:5,6), "sum")
  # tmp_doses = tmp_doses %>% cubelyr::as.tbl_cube() %>% as_tibble()

basecea = function(xcol, vxprice){
  icer_results_all %>% rename("xcol"=xcol) %>% 
                      filter(ns_strat=="RoutineCampaign", horizon=="20y", discounting=="disc", 
                            maps_profile=="Base", comp=="80% switch to MAPs",
                            vax_cost_sens==vxprice,  
                            mcv_cov_sens=="mcv1current", maps_cov_sens==0.2) %>% 
  ggplot(aes(x=xcol, fill=CEcat)) + 
  geom_bar(position = position_fill(reverse = TRUE)) + 
  xlab(NULL) + ylab("Percent of countries\n(unweighted)") + 
  scale_y_continuous(labels = scales::percent) + 
  scale_fill_manual(values=c("#006837", "#66bd63", "#addd8e", "#d9ef8b",  "#d73027"),
                    label=c("Cost-saving", "Highly cost-effective", "Very cost-effective", "Cost-effective", "Not cost-effective"),
                    drop=F) + scale_x_discrete(labels = function(x) str_wrap(x, width = 12)) + 
  themebar + guides(fill=guide_legend(nrow=1)) + 
    theme(axis.text.y = element_text(angle = 0, size=10), 
          axis.text.x = element_text(angle = 0, size=10), 
          strip.placement = "outside", 
          panel.spacing = unit(0.5,'lines')) + 
    theme(legend.text = element_text(size=10), axis.title.y = element_text(size=10))}

basecea_weighted = function(xcol, vxprice){
  icer_results_all %>% rename("xcol"=xcol) %>% 
    filter(ns_strat=="RoutineCampaign", horizon=="20y", discounting=="disc", 
           maps_profile=="Base", comp=="80% switch to MAPs",
           vax_cost_sens==vxprice, 
           mcv_cov_sens=="mcv1current", maps_cov_sens==0.2, !is.na(icers)) %>% 
    group_by(xcol, CEcat) %>% 
    summarise(tot_pop = sum(pop)) %>%
    ggplot(aes(x=xcol, y=tot_pop, fill=CEcat)) + 
    geom_bar(position = position_fill(reverse = TRUE), stat="identity") + 
    xlab(NULL) + ylab("Percent of countries\n(weighted)") + 
    scale_y_continuous(labels = scales::percent) + 
    scale_fill_manual(values=c("#006837", "#66bd63", "#addd8e", "#d9ef8b",  "#d73027"),
                      label=c("Cost-saving", "Highly cost-effective", "Very cost-effective", "Cost-effective", "Not cost-effective"),
                      drop=F) + scale_x_discrete(labels = function(x) str_wrap(x, width = 12)) + 
    themebar + guides(fill=guide_legend(nrow=1)) + 
    theme(axis.text.y = element_text(angle = 0, size=10), 
          axis.text.x = element_text(angle = 0, size=10), 
          strip.placement = "outside", 
          panel.spacing = unit(0.5,'lines')) + 
    theme(legend.text = element_text(size=10), axis.title.y = element_text(size=10))}

# group together
ggarrange(basecea("Continent", "$3.00 per dose"), basecea_weighted("Continent", "$3.00 per dose"),
          basecea("WB_group", "$3.00 per dose"), basecea_weighted("WB_group", "$3.00 per dose"),
          # basecea("Country_archetype", "$3.00 per dose"), basecea_weighted("Country_archetype", "$3.00 per dose"), 
          basecea("Gavi_eligibility", "$3.00 per dose"), basecea_weighted("Gavi_eligibility", "$3.00 per dose"),
          labels = "AUTO", 
          ncol = 2, nrow = 3, heights = c(1, 1), 
          common.legend = T, legend = "bottom", align="v") 
# list("Gavi eligibility", "World Bank income group", "Country archetype", "WHO region")
ggsave("./figures/EpiBaseline/BaselineCEAs_3d.pdf", width = 11, height = 10, units="in")
ggsave("./figures/EpiBaseline/BaselineCEAs_3d.jpg", width = 11, height = 10, units="in")

ggarrange(basecea("Continent", "$2.25 per dose"), basecea_weighted("Continent", "$2.25 per dose"),
          basecea("WB_group", "$2.25 per dose"), basecea_weighted("WB_group", "$2.25 per dose"),
          # basecea("Country_archetype", "$2.25 per dose"), basecea_weighted("Country_archetype", "$2.25 per dose"), 
          basecea("Gavi_eligibility", "$2.25 per dose"), basecea_weighted("Gavi_eligibility", "$2.25 per dose"),
          labels = "AUTO", 
          ncol = 2, nrow = 3, heights = c(1, 1), 
          common.legend = T, legend = "bottom", align="v") 
# list("Gavi eligibility", "World Bank income group", "Country archetype", "WHO region")
ggsave("./figures/EpiBaseline/BaselineCEAs_2.25d.pdf", width = 11, height = 10, units="in")
ggsave("./figures/EpiBaseline/BaselineCEAs_2.25d.jpg", width = 11, height = 10, units="in")

# WEIGHTED, BASE, $3
icer_results_all %>%
  filter(ns_strat=="RoutineCampaign", horizon=="20y", discounting=="disc",
         maps_profile=="Base", comp=="80% switch to MAPs",
         vax_cost_sens=="$3.00 per dose",
         mcv_cov_sens=="mcv1current", maps_cov_sens==0.2, !is.na(icers)) %>%
  group_by(Continent, CEcat) %>%
  summarise(tot_pop = sum(pop)) %>%
  pivot_wider(id_cols=Continent, names_from = CEcat, values_from = tot_pop) %>%
  mutate(VCE = ifelse(is.na(VCE), 0, VCE),
         HCE = ifelse(is.na(HCE), 0, HCE),
         CE = ifelse(is.na(CE), 0, CE),
         CS = ifelse(is.na(CS), 0, CS)) %>%
  mutate(total=CE+`Not CE`+VCE+HCE+CS) %>%
  mutate(CS = CS/total,
        HCE = HCE/total,
         VCE = VCE/total,
         CE = CE/total,
         `Not CE` = `Not CE`/total)

# Continent     HCE    VCE      CE `Not CE`     CS       total
# 1 Africa    0.456   0.0919 0.339      0.113 0      1148593304 
# 2 Americas  0       0      0.0284     0.972 0       625859135 
# 3 Asia      0.00268 0      0.0435     0.954 0      3705635401.
# 4 Eurasia   0       0      0          1     0       418909060 
# 5 Mideast   0       0      0.00155    0.952 0.0462  713919117 

# %>% write.csv("global_by_continent.csv")

# UNWEIGHTED, BASE, $3
icer_results_all %>%
  filter(ns_strat=="RoutineCampaign", horizon=="20y", discounting=="disc",
         maps_profile=="Base", comp=="80% switch to MAPs",
         vax_cost_sens=="$3.00 per dose",
         mcv_cov_sens=="mcv1current", maps_cov_sens==0.2, !is.na(icers)) %>%
  group_by(Continent, CEcat) %>%
  summarise(N=n()) %>%
  pivot_wider(id_cols=Continent, names_from = CEcat, values_from = N) %>%
  mutate(VCE = ifelse(is.na(VCE), 0, VCE),
         HCE = ifelse(is.na(HCE), 0, HCE),
         CE = ifelse(is.na(CE), 0, CE),
         CS = ifelse(is.na(CS), 0, CS)) %>%
  mutate(total=CE+`Not CE`+VCE+HCE+CS) %>%
  mutate(CS = CS/total,
         HCE = HCE/total,
         VCE = VCE/total,
         CE = CE/total,
         `Not CE` = `Not CE`/total)

# UNWEIGHTED
# Continent    HCE   VCE     CE `Not CE`     CS total
# 1 Africa    0.222  0.133 0.422     0.222 0         45
# 2 Americas  0      0     0.0385    0.962 0         26
# 3 Asia      0.0385 0     0.192     0.769 0         26
# 4 Eurasia   0      0     0         1     0         20
# 5 Mideast   0      0     0.0625    0.875 0.0625    16

icer_results_all %>%
  filter(ns_strat=="RoutineCampaign", horizon=="20y", discounting=="disc",
         maps_profile=="Optimistic", comp=="Targeted MAPs introduction",
         vax_cost_sens=="$2.25 per dose",
         mcv_cov_sens=="mcv1current", maps_cov_sens==0.2) %>%
  group_by(CEcat) %>%
  summarise(N = n(), 
            Pop = sum(pop)) 

# CEcat      N         Pop
# <fct>  <int>       <dbl>
# 1 CS        11  467847101 
# 2 HCE       43 2947811203.
# 3 VCE        9  144228553 
# 4 CE        16  279638999 
# 5 Not CE    54 2773390161 

icer_results_all %>%
  filter(ns_strat=="RoutineCampaign", horizon=="20y", discounting=="disc",
         maps_profile=="Optimistic", comp=="Targeted MAPs introduction",
         vax_cost_sens=="$2.25 per dose",
         mcv_cov_sens=="mcv1current", maps_cov_sens==0.2) %>%
  summarise(Pop = sum(pop)) 

# Pop: 6612916017.

icer_results_all %>%
  filter(ns_strat=="RoutineCampaign", horizon=="20y", discounting=="disc",
         maps_profile %in% c("Base", "Optimistic"), 
         comp %in% c("80% switch to MAPs","Targeted MAPs introduction"),
         vax_cost_sens=="$4.50 per dose",
         mcv_cov_sens=="mcv1current", maps_cov_sens==0.2) %>%
  group_by(maps_profile, comp, CEcat) %>%
  summarise(Pop = sum(pop)) 

1-5656428139/6612916017 # Base, comp1
# [1] 0.1446
1-5382707830/6612916017 # Base, comp2
# [1] 0.186
1-5504667835/6612916017 # Optimistic, comp1
# [1] 0.1676
1-5296963496/6612916017 # Optimistic, comp2
# [1] 0.199

icer_results_all %>%
  filter(ns_strat=="RoutineCampaign", horizon=="20y", discounting=="disc",
         maps_profile %in% c("Base", "Optimistic"), 
         comp %in% c("80% switch to MAPs","Targeted MAPs introduction"),
         vax_cost_sens=="$2.00 per dose",
         mcv_cov_sens=="mcv1current", maps_cov_sens==0.2) %>%
  group_by(maps_profile, comp, CEcat) %>%
  summarise(Pop = sum(pop)) 

## Basic Sensitivity --- 

## Price, comparators, and MAPS profiles

punw = icer_results_all %>% filter(ns_strat=="RoutineCampaign", horizon=="20y", discounting=="disc", 
                            mcv_cov_sens=="mcv1current", maps_cov_sens==0.2, !is.na(icers)) %>% 
  ggplot(aes(x=comp, fill=CEcat)) + 
  geom_bar(position = position_fill(reverse = TRUE)) + 
  xlab(NULL) + ylab("Percent of countries (unweighted)") + 
  scale_y_continuous(labels = scales::percent) + 
  scale_fill_manual(values=c("#006837", "#66bd63", "#addd8e", "#d9ef8b",  "#d73027"),
                    label=c("Cost-saving", "Highly cost-effective", "Very cost-effective", "Cost-effective", "Not cost-effective"),
                    drop=F) + scale_x_discrete(labels = function(x) str_wrap(x, width = 10)) + 
  facet_grid(vax_cost_sens~maps_profile) + 
  themebar + guides(fill=guide_legend(nrow=1)) + theme(legend.text = element_text(size=9))

ggdraw(punw, xlim = c(0, 1), ylim = c(0, 1), clip = "off") + 
  draw_grob(
  grob = rectGrob(gp = gpar(col = "black", fill = NA, lwd = 4)),
  x = 0.102, y = 0.3625, width = 0.08, height = 0.175
)
  
ggsave(paste0("./figures/EpiBaseline/BasicSens_cea_unweighted.pdf"), width = 7.5, height = 7, units = "in")

pwei = icer_results_all %>% filter(ns_strat=="RoutineCampaign", horizon=="20y", discounting=="disc", 
                                   mcv_cov_sens=="mcv1current", maps_cov_sens==0.2, !is.na(icers))%>%
  group_by(comp, vax_cost_sens, maps_profile, CEcat) %>% 
  summarise(tot_pop = sum(pop)) %>%
  ggplot(aes(x=comp, y=tot_pop,fill=CEcat)) + 
  geom_bar(position = position_fill(reverse = TRUE), stat="identity") + 
  xlab(NULL) + ylab("Percent of countries (weighted)") + 
  scale_y_continuous(labels = scales::percent) + 
  scale_fill_manual(values=c("#006837", "#66bd63", "#addd8e", "#d9ef8b",  "#d73027"),
                    label=c("Cost-saving", "Highly cost-effective", "Very cost-effective", "Cost-effective", "Not cost-effective"),
                    drop=F) + scale_x_discrete(labels = function(x) str_wrap(x, width = 10)) + 
  facet_grid(vax_cost_sens~maps_profile) + 
  themebar + guides(fill=guide_legend(nrow=1)) + theme(legend.text = element_text(size=9))

ggdraw(pwei, xlim = c(0, 1), ylim = c(0, 1), clip = "off") + 
  draw_grob(
    grob = rectGrob(gp = gpar(col = "black", fill = NA, lwd = 4)),
    x = 0.102, y = 0.3625, width = 0.08, height = 0.175
  )

ggsave("./figures/EpiBaseline/BasicSens_cea_weighted.pdf", width = 7.5, height = 7, units = "in")

# Reordered
punw = icer_results_all %>% filter(ns_strat=="RoutineCampaign", horizon=="20y", discounting=="disc", 
                                   mcv_cov_sens=="mcv1current", maps_cov_sens==0.2, !is.na(icers), 
                                   comp=="80% switch to MAPs", vax_cost_sens %in% c("$2.25 per dose", "$3.00 per dose")) %>% 
  ggplot(aes(x=maps_profile, fill=CEcat)) + 
  geom_bar(position = position_fill(reverse = TRUE)) + 
  xlab(NULL) + ylab("Percent of countries (unweighted)") + 
  scale_y_continuous(labels = scales::percent) + 
  scale_fill_manual(values=c("#006837", "#66bd63", "#addd8e", "#d9ef8b",  "#d73027"),
                    label=c("Cost-saving", "Highly cost-effective", "Very cost-effective", "Cost-effective", "Not cost-effective"),
                    drop=F) + scale_x_discrete(labels = function(x) str_wrap(x, width = 10)) + 
  facet_grid(.~vax_cost_sens) + 
  themebar + guides(fill=guide_legend(nrow=1)) + theme(legend.text = element_text(size=9))

ggsave(paste0("./figures/EpiBaseline/BasicSens_cea_unweighted_reordered.pdf"), width = 8, height = 4, units = "in")

pwei = icer_results_all %>% filter(ns_strat=="RoutineCampaign", horizon=="20y", discounting=="disc", 
                                   mcv_cov_sens=="mcv1current", maps_cov_sens==0.2, !is.na(icers),
                                   comp=="80% switch to MAPs", vax_cost_sens %in% c("$2.25 per dose", "$3.00 per dose"))%>%
  group_by(comp, vax_cost_sens, maps_profile, CEcat) %>% 
  summarise(tot_pop = sum(pop)) %>%
  ggplot(aes(x=maps_profile, y=tot_pop,fill=CEcat)) + 
  geom_bar(position = position_fill(reverse = TRUE), stat="identity") + 
  xlab(NULL) + ylab("Percent of countries (weighted)") + 
  scale_y_continuous(labels = scales::percent) + 
  scale_fill_manual(values=c("#006837", "#66bd63", "#addd8e", "#d9ef8b",  "#d73027"),
                    label=c("Cost-saving", "Highly cost-effective", "Very cost-effective", "Cost-effective", "Not cost-effective"),
                    drop=F) + scale_x_discrete(labels = function(x) str_wrap(x, width = 10)) + 
  facet_grid(.~vax_cost_sens) + 
  themebar + guides(fill=guide_legend(nrow=1)) + theme(legend.text = element_text(size=9))
ggsave("./figures/EpiBaseline/BasicSens_cea_weighted_reordered.pdf", width = 8, height = 4, units = "in")


## Extended Sensitivity ------
# $2.25-3, Comp 1-2 only, Prof B and O only. Maybe even less.
### MAPS coverage ----
# OR switch to cost and profile, since comp is up to countries...

tmp_sens_mapcov = icer_results_all %>% filter(ns_strat=="RoutineCampaign", horizon=="20y", discounting=="disc", 
                                              comp=="80% switch to MAPs", # %in% c("80% switch to MAPS", "Targeted MAP introduction"),
                                              vax_cost_sens %in% c("$2.25 per dose", "$3.00 per dose"),
                                              maps_profile %in% c("Base", "Optimistic"), 
                                              mcv_cov_sens=="mcv1current", !is.na(icers)) %>% 
  mutate(maps_cov_sens=factor(maps_cov_sens, levels = c(0.1, 0.2, 0.3),
                              labels = c("10%", "20% (default)", "30%")))

punw = tmp_sens_mapcov %>% 
  ggplot(aes(x=maps_cov_sens, fill=CEcat)) + 
  geom_bar(position = position_fill(reverse = TRUE)) + 
  xlab(NULL) + ylab("Percent of countries (unweighted)") + 
  scale_y_continuous(labels = scales::percent) + 
  scale_fill_manual(values=c("#006837", "#66bd63", "#addd8e", "#d9ef8b",  "#d73027"),
                    label=c("Cost-saving", "Highly cost-effective", "Very cost-effective", "Cost-effective", "Not cost-effective"),
                    drop=F) + scale_x_discrete(labels = function(x) str_wrap(x, width = 10)) + 
  facet_grid(maps_profile+vax_cost_sens~.) + 
  themebar + guides(fill=guide_legend(nrow=1)) + theme(legend.text = element_text(size=9)) + 
  theme(strip.text.y = element_blank(), 
        # axis.title.y = element_text(vjust = -15),  # value by experiment
        # strip.placement = "outside", 
        panel.spacing = unit(0.25,'lines'), 
        plot.margin = unit(c(0, 0.25, 0, 0), "pt"))

pwei = tmp_sens_mapcov %>% 
  group_by(maps_cov_sens, vax_cost_sens, maps_profile, CEcat) %>% 
  summarise(tot_pop = sum(pop)) %>%
  ggplot(aes(x=maps_cov_sens, y=tot_pop, fill=CEcat)) + 
  geom_bar(position = position_fill(reverse = TRUE), stat="identity") + 
  xlab(NULL) + ylab("Percent of countries (weighted)") + 
  scale_y_continuous(labels = scales::percent) + 
  scale_fill_manual(values=c("#006837", "#66bd63", "#addd8e", "#d9ef8b",  "#d73027"),
                    label=c("Cost-saving", "Highly cost-effective", "Very cost-effective", "Cost-effective", "Not cost-effective"),
                    drop=F) + scale_x_discrete(labels = function(x) str_wrap(x, width = 10)) + 
  facet_grid(maps_profile+vax_cost_sens~.) + 
  themebar + guides(fill=guide_legend(nrow=1)) + theme(legend.text = element_text(size=9), 
                                                       plot.margin = unit(c(0, 0, 0, 0), "pt"))

punw_pwei = ggarrange(punw, pwei,
          # labels = "AUTO", 
          ncol = 2, nrow = 1, widths = c(1, 1), 
          common.legend = T, legend = "bottom", align="v") + 
  theme(plot.margin = margin(0.5,0.25,0.5,0.25, "cm")) 

ggdraw(punw_pwei, xlim = c(0, 1), ylim = c(0, 1), clip = "off") +
  draw_grob(
    grob = rectGrob(gp = gpar(col = "black", fill = NA, lwd = 4)),
    x = 0.21, y = 0.5475, width = 0.1, height = 0.2) + 
  draw_grob(
    grob = rectGrob(gp = gpar(col = "black", fill = NA, lwd = 4)),
    x = 0.695, y = 0.5475, width = 0.1, height = 0.2)

ggsave("./figures/EpiBaseline/ExtendedSens_MAPScov.pdf", width = 8, height = 10, units="in")

### MCV1 coverage ----

tmp_sens_mcvcov = icer_results_all %>% filter(ns_strat=="RoutineCampaign", horizon=="20y", discounting=="disc", 
                                              comp=="80% switch to MAPs", # %in% c("80% switch to MAPS", "Targeted MAP introduction"),
                                              maps_cov_sens==0.2,
                                              vax_cost_sens %in% c("$2.25 per dose", "$3.00 per dose"),
                                              maps_profile %in% c("Base", "Optimistic"), 
                                              !is.na(icers)) %>% 
  mutate(mcv_cov_sens=factor(mcv_cov_sens, levels = c("mcv1imp", "mcv1current", "mcv2"),
                             labels = c("MCV1 improved", "MCV1 current (default)", "MCV2")))

punw = tmp_sens_mcvcov %>% 
  ggplot(aes(x=mcv_cov_sens, fill=CEcat)) + 
  geom_bar(position = position_fill(reverse = TRUE)) + 
  xlab(NULL) + ylab("Percent of countries (unweighted)") + 
  scale_y_continuous(labels = scales::percent) + 
  scale_fill_manual(values=c("#006837", "#66bd63", "#addd8e", "#d9ef8b",  "#d73027"),
                    label=c("Cost-saving", "Highly cost-effective", "Very cost-effective", "Cost-effective", "Not cost-effective"),
                    drop=F) + scale_x_discrete(labels = function(x) str_wrap(x, width = 10)) + 
  facet_grid(maps_profile+vax_cost_sens~.) + 
  themebar + guides(fill=guide_legend(nrow=1)) + theme(legend.text = element_text(size=9)) + 
  theme(strip.text.y = element_blank(), 
        # axis.title.y = element_text(vjust = -15),  # value by experiment
        # strip.placement = "outside", 
        panel.spacing = unit(0.25,'lines'), 
        plot.margin = unit(c(0, 0.25, 0, 0), "pt"))

pwei = tmp_sens_mcvcov %>% 
  group_by(mcv_cov_sens, vax_cost_sens, maps_profile, CEcat) %>% 
  summarise(tot_pop = sum(pop)) %>%
  ggplot(aes(x=mcv_cov_sens, y=tot_pop, fill=CEcat)) + 
  geom_bar(position = position_fill(reverse = TRUE), stat="identity") + 
  xlab(NULL) + ylab("Percent of countries (weighted)") + 
  scale_y_continuous(labels = scales::percent) + 
  scale_fill_manual(values=c("#006837", "#66bd63", "#addd8e", "#d9ef8b",  "#d73027"),
                    label=c("Cost-saving", "Highly cost-effective", "Very cost-effective", "Cost-effective", "Not cost-effective"),
                    drop=F) + scale_x_discrete(labels = function(x) str_wrap(x, width = 10)) + 
  facet_grid(maps_profile+vax_cost_sens~.) + 
  themebar + guides(fill=guide_legend(nrow=1)) + 
  theme(legend.text = element_text(size=9), 
      # strip.text.y = element_blank(), 
      # axis.title.y = element_text(vjust = -15),  # value by experiment
      # strip.placement = "outside", 
      panel.spacing = unit(0.25,'lines'), 
      plot.margin = unit(c(0, 0, 0, 0), "pt"))

punw_pwei = ggarrange(punw, pwei,
          # labels = "AUTO", 
          ncol = 2, nrow = 1, widths = c(1, 1), 
          common.legend = T, legend = "bottom", align="v") + 
  theme(plot.margin = margin(0.5,0.25,0.5,0.25, "cm")) 

ggdraw(punw_pwei, xlim = c(0, 1), ylim = c(0, 1), clip = "off") +
  draw_grob(
    grob = rectGrob(gp = gpar(col = "black", fill = NA, lwd = 4)),
    x = 0.21, y = 0.5527, width = 0.1, height = 0.195) + 
  draw_grob(
    grob = rectGrob(gp = gpar(col = "black", fill = NA, lwd = 4)),
    x = 0.695, y = 0.5527, width = 0.1, height = 0.195)

ggsave("./figures/EpiBaseline/ExtendedSens_MVC1cov.pdf", width = 8, height = 10, units="in")

### ns_strat ----

tmp_sens_nstrat = icer_results_all %>% filter(horizon=="20y", discounting=="disc", 
                                              comp=="80% switch to MAPs", # %in% c("80% switch to MAPS", "Targeted MAP introduction"),
                                              maps_cov_sens==0.2,
                                              vax_cost_sens %in% c("$2.25 per dose", "$3.00 per dose"),
                                              maps_profile %in% c("Base", "Optimistic"), 
                                              mcv_cov_sens=="mcv1current", !is.na(icers)) %>% 
  mutate(ns_strat=factor(ns_strat, levels = c("RoutineCampaign", "Routine", "None"),
                         labels = c("Routine & Campaign (default)", "Routine", "None")))
  
punw = tmp_sens_nstrat %>% 
  ggplot(aes(x=ns_strat, fill=CEcat)) + 
  geom_bar(position = position_fill(reverse = TRUE)) + 
  xlab(NULL) + ylab("Percent of countries (unweighted)") + 
  scale_y_continuous(labels = scales::percent) + 
  scale_fill_manual(values=c("#006837", "#66bd63", "#addd8e", "#d9ef8b",  "#d73027"),
                    label=c("Cost-saving", "Highly cost-effective", "Very cost-effective", "Cost-effective", "Not cost-effective"),
                    drop=F) + scale_x_discrete(labels = function(x) str_wrap(x, width = 10)) + 
  facet_grid(maps_profile+vax_cost_sens~.) + 
  themebar + guides(fill=guide_legend(nrow=1)) + theme(legend.text = element_text(size=9)) + 
  theme(strip.text.y = element_blank(), 
        # axis.title.y = element_text(vjust = -15),  # value by experiment
        # strip.placement = "outside", 
        panel.spacing = unit(0.25,'lines'), 
        plot.margin = unit(c(0, 0.25, 0, 0), "pt"))

pwei = tmp_sens_nstrat %>% 
  group_by(ns_strat, vax_cost_sens, maps_profile, CEcat) %>% 
  summarise(tot_pop = sum(pop)) %>%
  ggplot(aes(x=ns_strat, y=tot_pop, fill=CEcat)) + 
  geom_bar(position = position_fill(reverse = TRUE), stat="identity") + 
  xlab(NULL) + ylab("Percent of countries (weighted)") + 
  scale_y_continuous(labels = scales::percent) + 
  scale_fill_manual(values=c("#006837", "#66bd63", "#addd8e", "#d9ef8b",  "#d73027"),
                    label=c("Cost-saving", "Highly cost-effective", "Very cost-effective", "Cost-effective", "Not cost-effective"),
                    drop=F) + scale_x_discrete(labels = function(x) str_wrap(x, width = 10)) + 
  facet_grid(maps_profile+vax_cost_sens~.) + 
  themebar + guides(fill=guide_legend(nrow=1)) + 
  theme(# strip.text.y = element_blank(), 
        # axis.title.y = element_text(vjust = -15),  # value by experiment
        # strip.placement = "outside", 
        panel.spacing = unit(0.25,'lines'), 
        plot.margin = unit(c(0, 0.25, 0, 0), "pt"))

punw_pwei = ggarrange(punw, pwei,
          # labels = "AUTO", 
          ncol = 2, nrow = 1, widths = c(1, 1), 
          common.legend = T, legend = "bottom", align="v") + 
  theme(plot.margin = margin(0.5,0.25,0.5,0.25, "cm")) 

ggdraw(punw_pwei, xlim = c(0, 1), ylim = c(0, 1), clip = "off") +
  draw_grob(
    grob = rectGrob(gp = gpar(col = "black", fill = NA, lwd = 4)),
    x = 0.1025, y = 0.5525, width = 0.1, height = 0.195) + 
  draw_grob(
    grob = rectGrob(gp = gpar(col = "black", fill = NA, lwd = 4)),
    x = 0.5875, y = 0.5525, width = 0.1, height = 0.195)

ggsave("./figures/EpiBaseline/ExtendedSens_NSstrat.pdf", width = 8, height = 10, units="in")

## CEA and disparities -----

ineq_analysis = icer_results_wq %>% dplyr::select(horizon, comp, ns_strat, wealth_quintile, discounting, maps_profile,
                                            maps_cov_sens, vax_cost_sens, mcv_cov_sens,
                                                  "Status Quo_DALYs", "MAPS add_DALYs", "Status Quo_Cases", "MAPS add_Cases", pop_wq, pop,
                                                  icers, gdp, ISO, CEcat, WB_group, Gavi_eligibility, Country_archetype, Continent) %>% 
  filter(wealth_quintile %in% c("WQ1", "WQ5"), ns_strat=="RoutineCampaign", horizon=="20y", discounting=="disc", maps_profile=="Base", 
         !is.na(icers), comp=="80% switch to MAPs", maps_cov_sens==0.2, vax_cost_sens=="$3.00 per dose", mcv_cov_sens=="mcv1current") %>% 
  dplyr::select(-c(ns_strat, horizon, discounting, comp, maps_cov_sens, vax_cost_sens, mcv_cov_sens)) %>%
  pivot_wider(id_cols=c(ISO, pop, gdp, WB_group, Gavi_eligibility, Country_archetype, Continent), names_from = wealth_quintile, 
              values_from = c(`Status Quo_DALYs`, `MAPS add_DALYs`, `Status Quo_Cases`, `MAPS add_Cases`, pop_wq, icers, CEcat)) 

tmp_overall_icer = icer_results_all %>% dplyr::select(horizon, comp, ns_strat, discounting, maps_profile,
                                                      maps_cov_sens, vax_cost_sens, mcv_cov_sens,
                                                      icers, ISO, CEcat) %>% 
  filter(ns_strat=="RoutineCampaign", horizon=="20y", discounting=="disc", !is.na(icers), maps_profile=="Base", 
         !is.na(icers), comp=="80% switch to MAPs", maps_cov_sens==0.2, vax_cost_sens=="$2.25 per dose", mcv_cov_sens=="mcv1current") %>% 
  rename(icers_all = icers, CEcat_all = CEcat)

ineq_analysis = ineq_analysis %>% 
  mutate(redORrat = 1-(`MAPS add_DALYs_WQ1`/`MAPS add_DALYs_WQ5`)/(`Status Quo_DALYs_WQ1`/`Status Quo_DALYs_WQ5`)) %>% 
  mutate(redORdif = (`Status Quo_DALYs_WQ1`/`Status Quo_DALYs_WQ5`-1)-(`MAPS add_DALYs_WQ1`/`MAPS add_DALYs_WQ5`-1)) %>% 
  mutate(redORbase = (`Status Quo_DALYs_WQ1`/`Status Quo_DALYs_WQ5`)-1) %>% 
  mutate(redORratC = 1-(`MAPS add_Cases_WQ1`/`MAPS add_Cases_WQ5`)/(`Status Quo_Cases_WQ1`/`Status Quo_Cases_WQ5`)) %>% 
  mutate(redORdifC = (`Status Quo_Cases_WQ1`/`Status Quo_Cases_WQ5`-1)-(`MAPS add_Cases_WQ1`/`MAPS add_Cases_WQ5`-1)) %>% 
  mutate(redORbaseC = (`Status Quo_Cases_WQ1`/`Status Quo_Cases_WQ5`)-1) %>% 
  left_join(tmp_overall_icer)

ineq_dalys = ineq_analysis %>%
  ggplot(aes(x=redORbase*100, y=redORrat*100)) + 
  geom_smooth(method = "lm", se = F, color="gray80") + 
  geom_text_repel(aes(label=ISO), size=3, max.time = 2, min.segment.length = 0.25) +
  xlab("Excess DALYs among the\nlowest 20% compared to the highest 20%\nbefore MAPs vaccination") + 
  ylab("Percent reduction in disparity of DALYs among\nthe lowest 20% compared to the highest 20%") + 
  labs(title="DALYs") +
  themebar +theme(plot.title=element_text(face="bold", size=14, hjust=0.5)) # + facet_grid(.~WB_group, scale="free_y")
# if I can put a correlation coefficient, that would be great
ggsave("./figures/EpiBaseline/Ineq_reduction_DALYS.pdf", width = 8, height = 10, units="in")

ineq_analysis %>%
  ggplot(aes(x=redORbaseC, y=redORratC)) + 
  geom_smooth(method = "lm", se = F, color="gray80") + 
  geom_text(aes(label=ISO)) +
  # geom_jitter(color="black", size=0.4, alpha=0.9) +
  xlab("Excess cases among the\nlowest 20% compared to the highest 20%\nbefore MAPS vaccination") + 
  ylab("Percent reduction in excess cases among\nthe lowest 20% compared to the highest 20%") + 
  themebar # + facet_grid(.~WB_group, scale="free_y")

# Do Dalys and cases give the same result (in case it's being driven by countries with high mortality)? because I never look at DALYs in the NS period.
cor.test(ineq_analysis$redORrat, ineq_analysis$redORbase)
cor.test(ineq_analysis$redORratC, ineq_analysis$redORbaseC)
# 

# Pearson's product-moment correlation
# 
# data:  ineq_analysis$redORrat and ineq_analysis$redORbase
# t = 2.7, df = 131, p-value = 0.007
# alternative hypothesis: true correlation is not equal to 0
# 95 percent confidence interval:
#  0.06566 0.38810
# sample estimates:
#    cor 
# 0.2333 
# 
# Pearson's product-moment correlation
# 
# data:  ineq_analysis$redORratC and ineq_analysis$redORbaseC
# t = 2.7, df = 131, p-value = 0.009
# alternative hypothesis: true correlation is not equal to 0
# 95 percent confidence interval:
#   0.05891 0.38232
# sample estimates:
#   cor 
# 0.2269 

# cor.test(ineq_analysis$redORdif, ineq_analysis$redORbase)
# cor.test(ineq_analysis$redORdifC, ineq_analysis$redORbaseC)

ineq_cases = ineq_analysis %>%
  ggplot(aes(x=redORbaseC*100, y=redORratC*100)) + 
  geom_smooth(method = "lm", se = F, color="gray80") + 
  geom_text_repel(aes(label=ISO), size=3, max.time = 2, min.segment.length = 0.25) +
  geom_point(color="black", size=0.4) +
  xlab("Excess cases among the\nlowest 20% compared to the highest 20%\nbefore MAPs vaccination") + 
  ylab("Percent reduction in disparity of cases among\nthe lowest 20% compared to the highest 20%") + 
  labs(title="Cases") +
  themebar +theme(plot.title=element_text(face="bold", size=14, hjust=0.5)) # + facet_grid(.~WB_group, scale="free_y")
ggsave("./figures/EpiBaseline/Ineq_reduction_cases.pdf", width = 8, height = 10, units="in")

# group cases and dalys

ggarrange(ineq_cases, ineq_dalys,
          labels = "AUTO", 
          ncol = 2, nrow = 1, widths = c(1, 1), align="v") + 
  theme(plot.margin = margin(0.5,0.25,0.5,0.25, "cm")) 

ggsave("./figures/EpiBaseline/Ineq_reduction_both.pdf", width = 9, height = 6, units="in")

ineq_analysis %>%
  ggplot(aes(x=redORbaseC, y=redORratC)) + 
  geom_smooth(method = "lm", se = F, color="gray80") + 
  geom_text_repel(aes(label=ISO), size=3, max.time = 2, min.segment.length = 0.25) +
  geom_point(color="black", size=0.4) +
  xlab("Excess DALYs among the\nlowest 20% compared to the highest 20%\nbefore MAPs vaccination") + 
  ylab("Percent reduction in disparity of cases among\nthe lowest 20% compared to the highest 20%") + 
  themebar + facet_grid(.~WB_group, scale="free_y")

ineq_analysis %>%
  ggplot(aes(x=redORbaseC, y=redORratC)) + 
  geom_smooth(method = "lm", se = F, color="gray80") + 
  geom_text_repel(aes(label=ISO), size=3) +
  geom_point(color="black", size=0.4) +
  xlab("Excess DALYs among the\nlowest 20% compared to the highest 20%\nbefore MAPs vaccination") + 
  ylab("Percent reduction in disparity of cases among\nthe lowest 20% compared to the highest 20%") + 
  themebar + facet_grid(.~Gavi_eligibility, scale="free_y")

ineq_analysis %>%
  ggplot(aes(x=redORbaseC, y=redORratC)) + 
  geom_smooth(method = "lm", se = F, color="gray80") + 
  geom_text_repel(aes(label=ISO), size=3) +
  geom_point(color="black", size=0.4) +
  xlab("Excess DALYs among the\nlowest 20% compared to the highest 20%\nbefore MAPs vaccination") + 
  ylab("Percent reduction in excess cases among\nthe lowest 20% compared to the highest 20%") + 
  themebar + facet_grid(.~Continent, scale="free_y")

ineq_analysis %>%
  ggplot(aes(x=redORbaseC, y=redORratC)) + 
  geom_smooth(method = "lm", se = F, color="gray80") + 
  geom_text_repel(aes(label=ISO), size=3) +
  geom_point(color="black", size=0.4) +
  xlab("Excess DALYs among the\nlowest 20% compared to the highest 20%\nbefore MAPs vaccination") + 
  ylab("Percent reduction in excess cases among\nthe lowest 20% compared to the highest 20%") + 
  themebar + facet_grid(.~Country_archetype, scale="free_y")

ineq_analysis %>%
  ggplot(aes(x=CEcat_all, y=redORrat)) + 
  geom_boxplot() +
  geom_jitter(color="black", size=0.4, alpha=0.9) +
  themebar + facet_grid(.~WB_group, scale="free_y")

# THIS ONE:

neglog10 <- scales::trans_new("signed_log",
                           transform=function(x) sign(x)*log(abs(x), 10),
                           inverse=function(x) sign(x)*10^(abs(x)))

ineq_icer_dalys = ineq_analysis %>% 
  ggplot(aes(x=redORrat, y=pmax(icers_all/gdp, 0.01))) + 
  geom_text_repel(aes(label=ISO), size=3, max.time = 2, min.segment.length = 0.25) +
  # geom_smooth(method = "lm", se = F, color="gray80") + 
  geom_point(color="gray50") + scale_y_log10(breaks=c(0.01, 0.1, 1, 10, 100),
                                             labels = c("<0", "0.10", "1.00", "10.00", "100.00")) + 
  geom_hline(yintercept=3, linetype = "dashed", color="darkred") + 
  xlab("Reduction in disparity of DALYs\nbetween richest and poorest") + 
  ylab("ICER scaled by GDP per capita (log-scale)") +
  labs(title="DALYs") +
  themebar +theme(plot.title=element_text(face="bold", size=14, hjust=0.5)) # + facet_grid(.~WB_group, scale="free_y")

ineq_icer_cases = ineq_analysis %>% 
  ggplot(aes(x=redORratC, y=pmax(icers_all/gdp, 0.01))) + 
  geom_text_repel(aes(label=ISO), size=3, max.time = 2, min.segment.length = 0.25) +
  # geom_smooth(method = "lm", se = F, color="gray80") + 
  geom_point(color="gray50") + scale_y_log10(breaks=c(0.01, 0.1, 1, 10, 100),
                                             labels = c("<0", "0.10", "1.00", "10.00", "100.00")) + 
  geom_hline(yintercept=3, linetype = "dashed", color="darkred") + 
  xlab("Reduction in disparity of cases\nbetween richest and poorest") + 
  ylab("ICER scaled by GDP per capita (log-scale)") +
  labs(title="Cases") +
  themebar +theme(plot.title=element_text(face="bold", size=14, hjust=0.5)) # + facet_grid(.~WB_group, scale="free_y")

ggarrange(ineq_icer_cases, ineq_icer_dalys, 
          labels = "AUTO", 
          ncol = 2, nrow = 1, widths = c(1, 1), align="v") + 
  theme(plot.margin = margin(0.5,0.25,0.5,0.25, "cm")) 

ggsave("./figures/EpiBaseline/Ineq_reduction_icers_both.pdf", width = 9, height = 6, units="in")
# corr anyways?
cor.test(ineq_analysis$redORrat, log10(pmax(ineq_analysis$icers_all/ineq_analysis$gdp, 0.01)))
cor.test(ineq_analysis$redORratC, log10(pmax(ineq_analysis$icers_all/ineq_analysis$gdp, 0.01)))
# 
# > cor.test(ineq_analysis$redORrat, log10(pmax(ineq_analysis$icers_all/ineq_analysis$gdp, 0.01)))
# 
# Pearson's product-moment correlation
# 
# data:  ineq_analysis$redORrat and log10(pmax(ineq_analysis$icers_all/ineq_analysis$gdp, 0.01))
# t = -7.3, df = 131, p-value = 2e-11
# alternative hypothesis: true correlation is not equal to 0
# 95 percent confidence interval:
#  -0.6488 -0.4048
# sample estimates:
#     cor 
# -0.5379 
# 
# > cor.test(ineq_analysis$redORratC, log10(pmax(ineq_analysis$icers_all/ineq_analysis$gdp, 0.01)))
# 
# 	Pearson's product-moment correlation
# 
# data:  ineq_analysis$redORratC and log10(pmax(ineq_analysis$icers_all/ineq_analysis$gdp, 0.01))
# t = -7.1, df = 131, p-value = 6e-11
# alternative hypothesis: true correlation is not equal to 0
# 95 percent confidence interval:
#   -0.6417 -0.3946
# sample estimates:
#   cor 
# -0.5293 

ineq_analysis %>% 
    ggplot(aes(x=redORrat, y=icers_all/gdp)) + 
    geom_smooth(method = "lm", se = F, color="gray80") + 
    geom_point(color="gray50") + scale_y_log10() + # scale_y_continuous(limits=c(0, 50)) + 
    geom_hline(yintercept=3, linetype = "dashed", color="darkred") + 
    xlab("Reduction in Inequality of DALYs\nbetween richest and poorest") + 
    ylab("ICER scaled by GDP per capita") +
  # scale_y_continuous(limits=c(-1,10)) + 
    themebar + facet_grid(.~WB_group, scale="free_y")

ineq_analysis %>% 
    ggplot(aes(x=redORrat, y=icers_all/gdp)) + 
    geom_smooth(method = "lm", se = F, color="gray80") + 
    geom_point(color="gray50") + scale_y_log10() + # scale_y_continuous(limits=c(0, 50)) + 
    geom_hline(yintercept=3, linetype = "dashed", color="darkred") + 
    xlab("Reduction in Inequality of DALYs\nbetween richest and poorest") + 
  # scale_y_continuous(limits=c(-1,10)) + 
    themebar + facet_grid(.~Gavi_eligibility, scale="free_y")
  
## CE by WQ histograms -----
tmp_wq_graph = bind_rows((icer_results_wq %>%  
                            filter(ns_strat=="RoutineCampaign", horizon=="20y", discounting=="disc", 
                                   comp=="80% switch to MAPs", # %in% c("80% switch to MAPS", "Targeted MAP introduction"),
                                   maps_cov_sens==0.2,
                                   vax_cost_sens %in% c("$2.25 per dose", "$3.00 per dose"),
                                   maps_profile %in% c("Base", "Optimistic"), 
                                   !is.na(icers), mcv_cov_sens=="mcv1current")),
                         (icer_results_all %>% filter(ns_strat=="RoutineCampaign", horizon=="20y", discounting=="disc", 
                                                      comp=="80% switch to MAPs", # %in% c("80% switch to MAPS", "Targeted MAP introduction"),
                                                      maps_cov_sens==0.2,
                                                      vax_cost_sens %in% c("$2.25 per dose", "$3.00 per dose"),
                                                      maps_profile %in% c("Base", "Optimistic"), 
                                                      !is.na(icers), mcv_cov_sens=="mcv1current") %>%
                            mutate(wealth_quintile="All")))
  
punw = tmp_wq_graph %>% 
    ggplot(aes(x=wealth_quintile, fill=CEcat)) + 
    geom_bar(position = position_fill(reverse = TRUE)) + 
  xlab(NULL) + ylab("Percent of countries (unweighted)") + 
  scale_y_continuous(labels = scales::percent) + 
  scale_fill_manual(values=c("#006837", "#66bd63", "#addd8e", "#d9ef8b", "#d73027"),
                      label=c("Cost-saving", "Highly cost-effective", "Very cost-effective", "Cost-effective", "Not cost-effective"),
                      drop=F) +
    themebar + guides(fill=guide_legend(nrow=1)) + theme(legend.text = element_text(size=10)) + 
    facet_grid(maps_profile+vax_cost_sens~.) + 
  theme(strip.text.y = element_blank(), 
        # axis.title.y = element_text(vjust = -15),  # value by experiment
        # strip.placement = "outside", 
        panel.spacing = unit(0.25,'lines'), 
        plot.margin = unit(c(0, 0.25, 0, 0), "pt"))

pwei = tmp_wq_graph %>%
  group_by(wealth_quintile, vax_cost_sens, maps_profile, CEcat) %>% 
  summarise(tot_pop = sum(pop)) %>%
  ggplot(aes(x=wealth_quintile, y=tot_pop, fill=CEcat)) + 
  geom_bar(position = position_fill(reverse = TRUE), stat="identity") + 
  xlab(NULL) + ylab("Percent of countries (weighted)") + 
  scale_y_continuous(labels = scales::percent) + 
  scale_fill_manual(values=c("#006837", "#66bd63", "#addd8e", "#d9ef8b", "#d73027"),
                    label=c("Cost-saving", "Highly cost-effective", "Very cost-effective", "Cost-effective", "Not cost-effective"),
                    drop=F) +
  themebar + guides(fill=guide_legend(nrow=1)) + theme(legend.text = element_text(size=10)) + 
  facet_grid(maps_profile+vax_cost_sens~.) 
  
punw_pwei = ggarrange(punw, pwei,
          # labels = "AUTO", 
          ncol = 2, nrow = 1, widths = c(1, 1), 
          common.legend = T, legend = "bottom", align="v") + 
  theme(plot.margin = margin(0.5,0.25,0.5,0.25, "cm")) 

ggdraw(punw_pwei, xlim = c(0, 1), ylim = c(0, 1), clip = "off") +
  draw_grob(
    grob = rectGrob(gp = gpar(col = "black", fill = NA, lwd = 4)),
    x = 0.1025, y = 0.54, width = 0.315, height = 0.2) + 
  draw_grob(
    grob = rectGrob(gp = gpar(col = "black", fill = NA, lwd = 4)),
    x = 0.59, y = 0.54, width = 0.315, height = 0.2)

ggsave("./figures/EpiBaseline/BasicSens_WQ.pdf", width = 8, height = 10, units="in")

# CE by WQ - weighted only, but showing both comp1 and comp2
tmp_wq_graph2 = bind_rows((icer_results_wq %>%  
                            filter(ns_strat=="RoutineCampaign", horizon=="20y", discounting=="disc", 
                                   comp %in% c("80% switch to MAPs", "Targeted MAPs introduction"),
                                   maps_cov_sens==0.2,
                                   vax_cost_sens %in% c("$3.00 per dose"), # "$2.25 per dose", 
                                   maps_profile %in% c("Base", "Optimistic"), 
                                   !is.na(icers), mcv_cov_sens=="mcv1current")),
                         (icer_results_all %>% filter(ns_strat=="RoutineCampaign", horizon=="20y", discounting=="disc", 
                                                      comp %in% c("80% switch to MAPS", "Targeted MAP introduction"),
                                                      maps_cov_sens==0.2,
                                                      vax_cost_sens %in% c("$3.00 per dose"),
                                                      maps_profile %in% c("Base", "Optimistic"), 
                                                      !is.na(icers), mcv_cov_sens=="mcv1current") %>%
                            mutate(wealth_quintile="All")))

tmp_wq_graph2 %>%
  group_by(wealth_quintile, maps_profile, comp, CEcat) %>% 
  summarise(tot_pop = sum(pop)) %>%
  ggplot(aes(x=wealth_quintile, y=tot_pop, fill=CEcat)) + 
  geom_bar(position = position_fill(reverse = TRUE), stat="identity") + 
  xlab(NULL) + ylab("Percent of countries (weighted)") + 
  scale_fill_manual(values=c("#006837", "#66bd63", "#addd8e", "#d9ef8b", "#d73027"),
                    label=c("Cost-saving", "Highly cost-effective", "Very cost-effective", "Cost-effective", "Not cost-effective"),
                    drop=F) +
  themebar + guides(fill=guide_legend(nrow=1)) + theme(legend.text = element_text(size=10)) + 
  facet_grid(maps_profile~comp) 

ggsave("./figures/EpiBaseline/BasicSens_WQ_3d_comp_profile.pdf", width = 8, height = 7, units="in")
# this didn't make it to the paper.

# Influential parameters ----
  # with the most default assumptions.
  
sensdata = icer_results_all %>% filter(horizon=="20y", discounting=="disc", 
                                       comp=="80% switch to MAPs", ns_strat=="RoutineCampaign",
                                       maps_cov_sens==0.2, vax_cost_sens=="$3.00 per dose",
                                       maps_profile=="Base", 
                                       mcv_cov_sens=="mcv1current", !is.na(icers)) %>% 
  dplyr::select(ISO, icers, gdp, lifexp, `Status Quo_Cases`, seek, SIR_fit_rep, `Status Quo_TrtCost`, mcv1_cov, mcv1_cov_wq1, death_per_case, vax_del_maps_uc_base1, vax_del_ns_uc_base1, 
                meaninc, meanage, sani_cov, sani_cov_wq1, amr_prob, ip_mort, ipip_mort, op_costs_rx, op_costs, 
                 ip_costs, ip_costs_rx, ipip_costs_surgery, TCV_NS_intro_date, TCV_MAPS_intro_date, pop) %>% 
  mutate(avgtrtcost = `Status Quo_TrtCost`/(`Status Quo_Cases`*seek),
        CEAindex=icers/gdp, add_del = vax_del_maps_uc_base1-vax_del_ns_uc_base1,
         years_between=TCV_MAPS_intro_date-TCV_NS_intro_date,
         meaninc = ifelse(meaninc>1000, 1000, meaninc))

## Plain corr -----
sa_graph = function(xcol, labcol){
  tmp_sensdata = sensdata %>% rename("xcol"=xcol)
  tmp_fabfive = tmp_sensdata %>% filter(ISO %in% c("IND", "NPL", "MWI", "KEN", "BFA")) 
  tmp_sa = cor.test(tmp_sensdata$xcol, tmp_sensdata$icers/tmp_sensdata$gdp)
  tmp_pval = ifelse(tmp_sa$p.value>0.001, paste0(" (p=",  round(tmp_sa$p.value, digits=3),")"), " (p<0.001)") 
  m = lm(CEAindex~xcol, data=tmp_sensdata)
  
  return(ggplot(data=tmp_sensdata, aes(x=xcol, y=CEAindex, group=1)) + 
    xlab(NULL) + ylab(NULL) + 
    labs(title=labcol) + 
    geom_point(color="gray80") + 
    stat_function(fun = ~predict(m, data.frame(xcol = .x))) +
    # geom_smooth(method = "lm", se = F, color="gray60") + 
    annotate("text", x=Inf,y=Inf,hjust=1.15,vjust=1.75,
             size=4, fontface="bold", 
             label=paste0("corr: ", round(tmp_sa$estimate, digits=2), tmp_pval), parse=F) + 
    geom_point(data=tmp_fabfive, aes(x=xcol, y=icers/gdp, color=ISO, shape=ISO), size=3) + 
      scale_y_continuous(limits = c(0, 100)) + 
    scale_shape_manual(values=15:19, labels=c("Burkina Faso", "India", "Kenya", "Malawi", "Nepal")) + 
    scale_color_manual(values=c(7,4,9:11), labels=c("Burkina Faso", "India", "Kenya", "Malawi", "Nepal")) +
    themebar + theme(plot.title=element_text(hjust=0.5, size=12))) 
  } 
  
sa_amr = sa_graph("amr_prob", "Prevalence of AMR")
sa_lifeexp = sa_graph("lifexp", "Life expectancy")
sa_mcv1cov = sa_graph("mcv1_cov", "Vaccine coverage - MCV1")
sa_death = sa_graph("death_per_case", "Case fatality rate") 
sa_incidence = sa_graph("meaninc", "Mean incidence per 100K\nwithout any vaccination") 
sa_meanage = sa_graph("meanage", "Mean age of infection\nwithout any vaccination") # p=0 should be p<0.001
sa_mcv1covwq1 = sa_graph("mcv1_cov_wq1", "Vaccine coverage -\nMCV1, WQ1")
sa_sanicov = sa_graph("sani_cov", "Improved sanitation\ncoverage") # p=0 should be p<0.001
sa_sanicovwq1 = sa_graph("sani_cov_wq1", "Improved sanitation\ncoverage, WQ1") # p=0 should be p<0.001
sa_add_del_cost = sa_graph("add_del", "Additional delivery costs\nof MAPs vs. N&S (USD)") 

sa_years_between = sa_graph("years_between", "Years between\nN&S & MAPs deployment")
ggsave("./figures/EpiBaseline/YearsBetween_3d.pdf", width = 4, height = 4, units="in")

sa_trt_cost = sa_graph("avgtrtcost", "Average treatment costs") 

# pllegend <- get_legend(
#   sa_trt_cost + 
#     guides(color=guide_legend(nrow=5), shape=guide_legend(nrow=5)) + 
#     theme(legend.text = element_text(size=9))
# )
# 
# # arrange with ggarrange...
# all_drivers = plot_grid(sa_add_del_cost+theme(legend.position = "none"), 
#                         sa_meanage+theme(legend.position = "none"), 
#                         sa_mcv1cov+theme(legend.position = "none"), 
#                         sa_amr+theme(legend.position = "none"),
#                         sa_incidence+theme(legend.position = "none"), 
#                         sa_death+theme(legend.position = "none"), 
#                         sa_lifeexp+theme(legend.position = "none"), 
#                         sa_years_between+theme(legend.position = "none"), 
#                         sa_trt_cost+theme(legend.position = "none"), 
#                         sa_sanicov+theme(legend.position = "none"), 
#                         sa_sanicovwq1+theme(legend.position = "none"), 
#                         sa_mcv1covwq1+theme(legend.position = "none"), 
#                         # pllegend,
#           ncol = 3, nrow = 4, rel_widths = c(1, 1, 1), 
#           align="hv") + 
#   theme(plot.margin = margin(0.25,0.25,0.5,0.25, "cm")) 
#   annotate_figure(all_drivers, left = text_grob("ICER scaled by GDP per capita", color = "black", rot = 90, face="bold")) 
  
  all_drivers = ggarrange(sa_add_del_cost,
            sa_meanage, 
            sa_mcv1cov,
            sa_incidence,
            sa_amr, 
            sa_death,
            sa_lifeexp,
            sa_years_between,
            sa_trt_cost,
            sa_sanicov,
            sa_sanicovwq1,
            sa_mcv1covwq1,
            nrow=4, ncol=3, labels=NULL,
            common.legend = T, legend = "bottom", align="hv") 
  annotate_figure(all_drivers, left = text_grob("ICER scaled by GDP per capita", color = "black", rot = 90, face="bold")) 
  
ggsave("./figures/EpiBaseline/InfluentialParameters_3d.pdf", width = 8, height = 10, units="in")

# National results for the 5 DD countries -----
tmp = icer_results_all %>% dplyr::filter(horizon=="20y", comp=="80% switch to MAPs", 
                                         maps_profile=="Base", ns_strat == "RoutineCampaign", 
                                         mcv_cov_sens=="mcv1current", maps_cov_sens==0.2,
                                         discounting=="disc", vax_cost_sens=="$3.00 per dose",
                                         ISO %in% c("BFA", "IND", "KEN", "MWI", "NPL")) %>% 
  dplyr::select("ISO", "maps_profile", "Cases_averted", 
                "Deaths_averted", "DALYs_averted", "Cost_dif", "icers", "gdp", "pop", "CEcat") %>% 
  mutate(Cases_averted = Cases_averted*pop/1e5,
         Deaths_averted = Deaths_averted*pop/1e5,
         DALYs_averted = DALYs_averted*pop/1e5,
         Cost_dif = Cost_dif*pop/1e5,
         CEindex=icers/gdp) 

tmp %>% t

# ISO            "BFA"        "IND"        "KEN"        "MWI"        "NPL"       
# maps_profile   "Base"       "Base"       "Base"       "Base"       "Base"      
# Cases_averted  " 21956"     "725691"     " 41011"     "  8406"     " 14022"    
# Deaths_averted "190.38"     "993.02"     "955.39"     "219.79"     " 18.73"    
# DALYs_averted  " 3681"      "22068"      "18389"      " 4282"      "  406"     
# Cost_dif       " 11876376"  "648520824"  " 40734965"  " 14909286"  " 14734343" 
# icers          " 3227"      "29387"      " 2215"      " 3482"      "36291"     
# gdp            " 830.0"     "2410.9"     "2099.3"     " 645.2"     "1336.5"    
# pop            "  16116845" "1358600466" "  55879174" "  19822929" "  30666598"
# CEcat          "Not CE"     "Not CE"     "CE"         "Not CE"     "Not CE"    
# CEindex        " 3.887"     "12.189"     " 1.055"     " 5.397"     "27.153"    

icer_results_all %>% dplyr::filter(horizon=="20y", comp=="80% switch to MAPs", 
                                   maps_profile %in% c("Base", "Optimistic"), ns_strat == "RoutineCampaign", 
                                   mcv_cov_sens=="mcv1current", maps_cov_sens==0.2,
                                   discounting=="disc", vax_cost_sens %in% c("$2.25 per dose", "$3.00 per dose"),
                                   ISO %in% c("BFA", "IND", "KEN", "MWI", "NPL")) %>% 
  dplyr::select("ISO", "maps_profile", "vax_cost_sens", "icers", "gdp") %>% 
  mutate(CEindex=icers/gdp) %>% 
  pivot_wider(id_cols=c(maps_profile, vax_cost_sens), names_from=ISO, values_from=CEindex)

# maps_profile vax_cost_sens    BFA    IND    KEN   MWI   NPL
# 1 Base         $2.25 per dose 2.08   6.08  0.545  2.89  13.5 
# 2 Optimistic   $2.25 per dose 0.443  0.688 0.0890 0.618  3.49
# 3 Base         $3.00 per dose 3.89  12.2   1.06   5.40  27.2 
# 4 Optimistic   $3.00 per dose 2.25   6.79  0.599  3.13  17.1 

# Provide as table in the appendix just to be able to show CE countries.
# Or provide the CE countries
# Or as a tile plot, so the CE cat is in colors. 8 columns as it is... CEindex in parenthesis.
tmp_base_icer = icer_results_all %>% dplyr::filter(horizon=="20y", comp=="80% switch to MAPs", 
                                   maps_profile=="Base", ns_strat == "RoutineCampaign", 
                                   mcv_cov_sens=="mcv1current", maps_cov_sens==0.2,
                                   discounting=="disc", vax_cost_sens=="$3.00 per dose") %>% 
  dplyr::select("Country", "ISO", "maps_profile", "Cases_averted", 
                "Deaths_averted", "DALYs_averted", "Cost_dif", "icers", "gdp", "pop", "CEcat") %>% 
  mutate(Cases_averted = Cases_averted*pop/1e5,
         Deaths_averted = Deaths_averted*pop/1e5,
         DALYs_averted = DALYs_averted*pop/1e5,
         Cost_dif = Cost_dif*pop/1e5,
         CEindex=icers/gdp) %>%
  arrange(CEindex) %>% dplyr::select("ISO", "Country", "icers", "CEindex", "gdp", "pop", "CEcat") %>% print(n=133)

# 44 countries are CE
# 18 countries are VCE or better (HCE, CS)
paste(tmp_base_icer$Country[tmp_base_icer$CEcat=="CS"], collapse=", ")
# "Yemen"
paste(tmp_base_icer$Country[tmp_base_icer$CEcat=="HCE"], collapse=", ")
# "Papua New Guinea, Mauritius, Gabon, Guinea, Benin, Nigeria, Botswana, Côte d'Ivoire, Angola, Ethiopia, Dem. Rep. Congo"
paste(tmp_base_icer$Country[tmp_base_icer$CEcat=="VCE"], collapse=", ")
# "Madagascar, Mauritania, Rep. of Congo, South Africa, Liberia, Equatorial Guinea"
paste(tmp_base_icer$Country[tmp_base_icer$CEcat=="CE"], collapse=", ")
# "Chad, Lao, Cameroon, Kenya, Philippines, Djibouti, Central Afr. Rep., Mali, Ghana, Comoros, Togo, Sierra Leone, Niger, Namibia, Tanzania, Uganda, Ecuador, Sri Lanka, Senegal, Eswatini, Cambodia, Zambia, Mozambique, Timor-Leste, Cape Verde, Eritrea"

# A tibble: 133 × 7
# ISO   Country                         icers   CEindex    gdp         pop CEcat 
# <chr> <chr>                           <dbl>     <dbl>  <dbl>       <dbl> <fct> 
#   1 YEM   Yemen                           -94.0   -0.135    699.   32981641  CS    
# 2 PNG   Papua New Guinea                122.     0.0393  3116.    9949437  HCE   
# 3 MUS   Mauritius                      1597.     0.156  10240.    1266060  HCE   
# 4 GAB   Gabon                          1658.     0.188   8820.    2341179  HCE   
# 5 GIN   Guinea                          314.     0.207   1515.   13531906  HCE   
# 6 BEN   Benin                           288.     0.220   1305.   12996895  HCE   
# 7 NGA   Nigeria                         652.     0.301   2163.  213401323  HCE   
# 8 CIV   Côte d'Ivoire                   805.     0.323   2492.   27478249  HCE   
#   9 BWA   Botswana                       2501.     0.324   7726.    2588423  HCE   
#  10 AGO   Angola                         1091.     0.372   2933.   34503774  HCE   
#  11 ETH   Ethiopia                        391.     0.380   1027.  120283026  HCE   
#  12 COD   Dem. Rep. Congo                 316.     0.475    665.   95894118  HCE   
#  13 MDG   Madagascar                      263.     0.510    517.   28915653  VCE   
#  14 MRT   Mauritania                     1454.     0.706   2057.    4614974  VCE   
#  15 COG   Rep. of Congo                  2065.     0.780   2649.    5835806  VCE   
#  16 GNQ   Equatorial Guinea              6566.     0.815   8052.    1634466  VCE   
#  17 ZAF   South Africa                   5856.     0.865   6766.   59392255  VCE   
#  18 LBR   Liberia                         662.     0.878    755.    5193416  VCE   
#  19 TCD   Chad                            719.     1.03     699.   17179740  CE    
#  20 LAO   Lao                            2114.     1.03    2054.    7425057  CE    
#  21 CMR   Cameroon                       1640.     1.05    1563.   27198628  CE    
#  22 KEN   Kenya                          2215.     1.06    2099.   55879174  CE    
#  23 PHL   Philippines                    3898.     1.11    3499.  113880328  CE    
#  24 DJI   Djibouti                       4227.     1.29    3278.    1105557  CE    
#  25 CAF   Central Afr. Rep.               596.     1.39     427.    5457154  CE    
#  26 MLI   Mali                           1208.     1.45     831.   21904983  CE    
#  27 GHA   Ghana                          3357.     1.51    2218.   32833031  CE    
#  28 COM   Comoros                        2373.     1.60    1485.     821625  CE    
#  29 TGO   Togo                           1527.     1.65     923.    8644829  CE    
#  30 SLE   Sierra Leone                    819.     1.72     476.    8420641  CE    
#  31 NER   Niger                          1050.     1.78     589.   25252722  CE    
#  32 TZA   Tanzania                       2198.     1.84    1193.   63588334  CE    
#  33 NAM   Namibia                        9116.     1.86    4896.    2530151  CE    
#  34 UGA   Uganda                         1835.     1.90     964.   45853778  CE    
#  35 ECU   Ecuador                       13204.     2.04    6477.   17797737  CE    
#  36 LKA   Sri Lanka                      7061.     2.11    3343.   22156000  CE    
#  37 SEN   Senegal                        3722.     2.33    1595.   16876720  CE    
#  38 SWZ   Eswatini                       9401.     2.36    3987.    1192271  CE    
#  39 KHM   Cambodia                       4586.     2.61    1760.   16589023  CE    
#  40 ZMB   Zambia                         3802.     2.61    1457.   19473125  CE    
#  41 MOZ   Mozambique                     1472.     2.64     558.   32077072  CE    
#  42 CPV   Cape Verde                    10625.     2.74    3884.     587925  CE    
#  43 TLS   Timor-Leste                    6751.     2.83    2389.    1320942  CE    
#  44 ERI   Eritrea                        1933.     2.97     650     3620312  CE    
#  45 MDV   Maldives                      36794.     3.12   11781.     521457  Not CE
#  46 SSD   South Sudan                    3687.     3.44    1072.   10748272  Not CE
#  47 SOM   Somalia                        2064.     3.49     592.   17065581  Not CE
#  48 BGD   Bangladesh                     9530.     3.55    2688.  169356251  Not CE
#  49 PAK   Pakistan                       5722.     3.60    1589.  231402117  Not CE
#  50 STP   Sao Tome & Principe            8819.     3.70    2387.     223107  Not CE
#  51 ZWE   Zimbabwe                       6208.     3.70    1677.   15993524  Not CE
#  52 GTM   Guatemala                     20759.     3.79    5473.   17109746  Not CE
#  53 BFA   Burkina Faso                   3227.     3.89     830.   16116845  Not CE
#  54 SLB   Solomon Islands                8738.     4.04    2163.     707851  Not CE
#  55 KIR   Kiribati                       8426.     4.09    2061.     128874  Not CE
#  56 VUT   Vanuatu                       13209.     4.22    3129.     319137  Not CE
#  57 IDN   Indonesia                     21045.     4.40    4788.  273753191  Not CE
#  58 PER   Peru                          34846.     4.81    7239.   33715471  Not CE
#  59 BDI   Burundi                        1301.     5.02     259.   12551213  Not CE
#  60 GNB   Guinea-Bissau                  4089.     5.02     814.    2060721  Not CE
#  61 FSM   Micronesia                    19143.     5.08    3767.     113131  Not CE
#  62 MWI   Malawi                         3482.     5.41     643.   19822929  Not CE
#  63 WSM   Samoa                         21894.     5.85    3746.     218764  Not CE
#  64 MYS   Malaysia                      75043.     6.26   11993.   33573874  Not CE
#  65 GUY   Guyana                       116338.     6.39   18199.     804567  Not CE
#  66 MMR   Myanmar                        7737.     6.73    1149.   53798084  Not CE
#  67 AFG   Afghanistan                    2389.     6.78     353.   40099462  Not CE
#  68 THA   Thailand                      48112.     6.96    6913.   71601103  Not CE
#  69 GMB   Gambia                         5707.     7.10     804.    2639916  Not CE
#  70 DZA   Algeria                       37449.     7.46    5023.   44177969  Not CE
#  71 LBN   Lebanon                       28580.     7.47    3824.    5592631  Not CE
#  72 VNM   Vietnam                       31334.     7.50    4179.   97468029  Not CE
#  73 LSO   Lesotho                        8058.     8.12     993.    2281454  Not CE
#  74 TUR   Turkey                        97071.     9.09   10675.   84775404  Not CE
#  75 BOL   Bolivia                       33493.     9.30    3600.   12079472  Not CE
#  76 ARG   Argentina                    135482.     9.92   13651.   45808747  Not CE
#  77 LCA   Saint Lucia                  145750.    11.2    13031.     179651  Not CE
#  78 IRQ   Iraq                          72207.    11.2     6442.   43533592  Not CE
#  79 IND   India                         29387.    12.4     2366. 1358600466. Not CE
#  80 MNE   Montenegro                   140857.    14.0    10093.     619211  Not CE
#  81 AZE   Azerbaijan                   111812.    14.4     7771.   10137750  Not CE
#  82 LBY   Libya                        118770.    14.5     8211.    6735277  Not CE
#  83 BRA   Brazil                       132204.    14.6     9065.  214326223  Not CE
#  84 MEX   Mexico                       174638.    15.2    11477.  126705138  Not CE
#  85 BTN   Bhutan                        66413.    17.9     3704.     777486  Not CE
#  86 MNG   Mongolia                      96899.    19.2     5046.    3347782  Not CE
#  87 SUR   Suriname                     112614.    19.2     5859.     612985  Not CE
#  88 TJK   Tajikistan                    21087.    19.6     1076.    9750064  Not CE
#  89 HTI   Haiti                         34295.    19.6     1748.   11447569  Not CE
#  90 RWA   Rwanda                        19803.    20.5      967.   13461888  Not CE
#  91 PAN   Panama                       363184.    20.9    17358.    4351267  Not CE
#  92 COL   Colombia                     140072.    21.0     6657.   51516562  Not CE
#  93 IRN   Iran                         100043.    21.4     4668.   87923432  Not CE
#  94 PRY   Paraguay                     140789.    22.8     6187.    6703799  Not CE
#  95 DOM   Dominican Republic           247312.    24.5    10111.   11117873  Not CE
#  96 BLZ   Belize                       186489.    26.7     6984.     400031  Not CE
#  97 NPL   Nepal                         36291.    26.9     1348.   30666598  Not CE
#  98 JOR   Jordan                       124907.    29.0     4311.   11148278  Not CE
#  99 GRD   Grenada                      287510.    29.5     9762.     124610  Not CE
# 100 ALB   Albania                      201622.    29.6     6810.    2811666  Not CE
# 101 TON   Tonga                        144079.    30.8     4682.     106017  Not CE
# 102 NIC   Nicaragua                     71749.    31.9     2252.    6850540  Not CE
# 103 KAZ   Kazakhstan                   367969.    32.0    11484.   19000988  Not CE
# 104 FJI   Fiji                         177151.    33.1     5356.     924610  Not CE
# 105 TKM   Turkmenistan                 307595.    35.0     8793.    6341855  Not CE
# 106 TUN   Tunisia                      154559.    42.8     3608.   12262946  Not CE
# 107 ROU   Romania                      766679.    48.9    15692.   19119880  Not CE
# 108 SDN   Sudan                         55562.    50.4     1102.   45657202  Not CE
# 109 CRI   Costa Rica                   697853.    52.2    13365.    5153957  Not CE
# 110 CUB   Cuba                         649471.    68.4     9500.   11256372  Not CE
# 111 HND   Honduras                     211400.    70.2     3012.   10278345  Not CE
# 112 EGY   Egypt                        343178.    79.9     4295.  109262178  Not CE
# 113 MAR   Morocco                      321827.    93.5     3442.   37076584  Not CE
# 114 VEN   Venezuela                    148629.    99.0     1501    28199867  Not CE
# 115 KGZ   Kyrgyzstan                   175195.   101.      1740.    6691800  Not CE
# 116 SYR   Syria                         44146.   105.       421.   21324367  Not CE
# 117 DMA   Dominica                    1002533.   120.      8347.      72412  Not CE
# 118 UZB   Uzbekistan                   293676.   129.      2276.   34915100  Not CE
# 119 SRB   Serbia                      1568694.   164.      9538.    6834326  Not CE
# 120 CHN   China                       2166486.   171.     12663. 1412360000  Not CE
# 121 UKR   Ukraine                      808183.   177.      4576.   43792855  Not CE
# 122 MKD   Macedonia                   1355824.   181.      7486.    2065092  Not CE
# 123 BGR   Bulgaria                    2762783.   198.     13974.    6877743  Not CE
# 124 GEO   Georgia                     1354253.   201.      6730.    3708610  Not CE
# 125 RUS   Russia                      3539859.   229.     15445.  143449286  Not CE
# 126 ARM   Armenia                     1614485.   230.      7018.    2790974  Not CE
# 127 SLV   El Salvador                 1300943.   258.      5048.    6314167  Not CE
# 128 JAM   Jamaica                     1715983.   284.      6047.    2827695  Not CE
# 129 BIH   Bosnia & Herzegovina        2306248.   304.      7588.    3270943  Not CE
# 130 VCT   Saint Vincent & Grenadines  4527309.   487.      9298.     104332  Not CE
# 131 PRK   North Korea                  351157.   540.       650    25971909  Not CE
# 132 BLR   Belarus                     8141212.  1018.      7995.    9340314  Not CE
# 133 MDA   Moldova                    15645771.  2738.      5715.    2615199  Not CE
