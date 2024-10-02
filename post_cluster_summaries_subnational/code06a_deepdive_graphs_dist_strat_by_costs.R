#############################
# Making graphs for the manuscript
# Subnational analysis, stratified by cost
#############################

# load("../out_cluster/ISO_which_to_run_cea.Rdata")
source("./C_cea_code/00b_load_packages.R")

library(readxl)
library(MASS)
library(dplyr)
library(magrittr)
library(tidyr)
library(ggplot2)
library(countrycode)
library(wbstats)
library(stringr)

source("./C_cea_code/00c_load_helper_functions.R")

TransPar = readMat("./data/model_samples_ALL_Aug2021.mat")
CN = readMat("./data/burden_dynamic_link_Yale.mat")$cn
avgage_yale = readMat("./data/burden_dynamic_link_Yale.mat")$avgage.country
avgage_ihme = readMat("./data/burden_dynamic_link_IHME.mat")$avgage.country

agedist = readMat("./data/burden_dynamic_link_Yale.mat")$agedist.country

ISO = read.csv("./data/iso.csv")

# bring in continent information
countries <- read_excel("./data/Countries.xlsx", 
                        sheet = "Country Data w intro date")
colnames(countries) = gsub(" - ", "_",colnames(countries))
colnames(countries) = gsub(" ", "_",colnames(countries))
colnames(countries) = gsub("ISO_code", "ISO",colnames(countries))
colnames(countries) = gsub("-", "",colnames(countries))
colnames(countries) = gsub("\\(", "",colnames(countries))
colnames(countries) = gsub(")", "",colnames(countries))
colnames(countries) = gsub("WB_Group_June_2021", "WB_group",colnames(countries))
countries$NS_intro = sapply(1:183, function(i){max(countries[i, c("Low_intro_date", "Medium_Intro_date", "High_Intro_date")])})
countries$NS_intro[countries$NS_intro==0 & countries$WB_group %in% c("Low income", "Lower middle income")] = 2030
countries$NS_intro[countries$NS_intro==0 & countries$WB_group == "Upper middle income"] = 2030
#Continent TODO: change EMRO to respective continents
countries = countries %>% 
  mutate(continent = ifelse(countries$WHO %in% c("WPRO", "SEARO"), "Asia", 
                            ifelse(countries$WHO=="EMRO", "Mideast",
                                   ifelse(countries$WHO=="AFRO", "Africa", 
                                          ifelse(countries$WHO=="EURO", "Eurasia", "Americas")))))

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

#Bring in coverage from DHS
indicators = read_xlsx("./data/Deep_dive_indicators_2Aug.xlsx", "data_to_r_dist") %>% 
  dplyr::select(-starts_with('care'))
indicators_bfa = read_xlsx("./data/Deep_dive_indicators_2Aug.xlsx", "data_to_r_reg") %>% 
          filter(ISO=="BFA") %>% rename("District" = "Region") %>% dplyr::select(-starts_with('care'))

indicators = indicators %>% bind_rows(indicators_bfa)

indicators_natl = indicators %>% filter(District=="Total")
indicators = indicators %>% filter(District!="Total")
indicators$RegionLbl = gsub(" ","", indicators$District)
indicators$RegionLbl = gsub("&","", indicators$RegionLbl)

indicators = indicators %>% mutate(mcv1_total=pmax(pmin(as.numeric(mcv1_total), 0.99), 0.01),
                                   mcv1_wq1=pmax(pmin(as.numeric(mcv1_wq1), 0.99), 0.01),
                                   mcv1_wq2=pmax(pmin(as.numeric(mcv1_wq2), 0.99), 0.01),
                                   mcv1_wq3=pmax(pmin(as.numeric(mcv1_wq3), 0.99), 0.01),
                                   mcv1_wq4=pmax(pmin(as.numeric(mcv1_wq4), 0.99), 0.01),
                                   mcv1_wq5=pmax(pmin(as.numeric(mcv1_wq5), 0.99), 0.01))

# Bring in UC proportions and MAPS intro from MMGH
mmgh_data_maps = read_xlsx("./data/2022_07_18_Base TCV-MAP.xlsx", "For PATH with MAPs", skip=3)
mmgh_data_ns = read_xlsx("./data/2022_07_18_Base TCV-MAP.xlsx", "For PATH without MAPs", skip=3)

# Create necessary directories
create_dir_if_necessary("./maps_tcv_country/")

# Move to the right directory for outputs
setwd("./maps_tcv_country/")

for(i in c(2.25,3)){

icer_results_all = list()
icer_results_wq = list()

create_dir_if_necessary(paste0("figures/dist/d", i, "/all"))
create_dir_if_necessary(paste0("figures/dist/d", i, "/IND"))
create_dir_if_necessary(paste0("figures/dist/d", i, "/BFA"))
create_dir_if_necessary(paste0("figures/dist/d", i, "/NPL"))
create_dir_if_necessary(paste0("figures/dist/d", i, "/KEN"))
create_dir_if_necessary(paste0("figures/dist/d", i, "/MWI"))

for(zs in 1:dim(indicators)[1]){ 
  
if(indicators$ISO[zs]!="BFA"){
  load(paste0('../out_deepdive_cea/dist/all/0.2/', i, 'dollars/', indicators$ISO[zs], "/", sprintf("%03d", zs), "_",indicators$RegionLbl[zs], '/trt_inputs.Rdata'))
  load(paste0('../out_deepdive_cea/dist/all/0.2/', i, 'dollars/', indicators$ISO[zs], "/", sprintf("%03d", zs), "_",indicators$RegionLbl[zs], '/epi_inputs.Rdata'))
  load(paste0('../out_deepdive_cea/dist/all/0.2/', i, 'dollars/', indicators$ISO[zs], "/", sprintf("%03d", zs), "_",indicators$RegionLbl[zs], '/code03_treatment.Rdata')) # for the death rate per case
  load(paste0('../out_deepdive_cea/dist/all/0.2/', i, 'dollars/', indicators$ISO[zs], "/", sprintf("%03d", zs), "_",indicators$RegionLbl[zs], '/code05_cea.Rdata'))
  # load(paste0('../out_deepdive_vax/dist/all/0.2/', indicators$ISO[zs], "/", sprintf("%03d", zs), "_",indicators$RegionLbl[zs], '/SIR_wq_vax_NOunc.Rdata'))
} else {
  tmpzs=zs-827 # 863+1-37
  load(paste0('../out_deepdive_cea/reg/all/0.2/', i, 'dollars/', indicators$ISO[zs], "/", sprintf("%02d", tmpzs), "_",indicators$RegionLbl[zs], '/trt_inputs.Rdata'))
  load(paste0('../out_deepdive_cea/reg/all/0.2/', i, 'dollars/', indicators$ISO[zs], "/", sprintf("%02d", tmpzs), "_",indicators$RegionLbl[zs], '/epi_inputs.Rdata'))
  load(paste0('../out_deepdive_cea/reg/all/0.2/', i, 'dollars/', indicators$ISO[zs], "/", sprintf("%02d", tmpzs), "_",indicators$RegionLbl[zs], '/code03_treatment.Rdata')) # for the death rate per case
  load(paste0('../out_deepdive_cea/reg/all/0.2/', i, 'dollars/', indicators$ISO[zs], "/", sprintf("%02d", tmpzs), "_",indicators$RegionLbl[zs], '/code05_cea.Rdata'))
  # load(paste0('../out_deepdive_vax/reg/all/0.2/', indicators$ISO[zs], "/", sprintf("%02d", tmpzs), "_",indicators$RegionLbl[zs], '/SIR_wq_vax_NOunc.Rdata'))
}
  # extract incidence in present for each province/district.
  tmp = ISO$Var1[ISO$countryiso==indicators$ISO[zs]] # indicators$CNnum[zs]
  
  tmp_pop_wq = as.numeric(indicators[zs, paste0("pop_wq", 1:5)])
  tmp_pop_region = as.numeric(indicators$pop_DHS_WB_calc[zs])
  
  # SIR_wq_vax_unc$simcases_novac
  tmp_baseline_inc_wq = as.numeric(apply(typhoid[1,,1,"None",1,], 2, sum))/tmp_pop_wq*epipar["reporting"] # /1e5 # *(tmp_pop_wq*tmp_pop_region)
  tmp_baseline_inc = sum(typhoid[1,,1,"None",1,])*epipar["reporting"] # /1e5 # *tmp_pop_region
  
  tmp_avg_age_wq = typhoid[1,,1,"None",1,] %>% sweep(2, colSums(typhoid[1,,1,"None",1,]), "/") %>%
    sweep(1, c(4.5/12, 2.75/2, 3.5, 10, 20, (fix$lifexp+25)/2), "*") %>% apply(2, sum) %>% as.numeric
  tmp_avg_age = sum((apply(typhoid[1,,1,"None",1,], 1, sum) %>% as.numeric)/sum(typhoid[1,,1,"None",1,])*c(4.5/12, 2.75/2, 3.5, 10, 20, (fix$lifexp+25)/2))
  
  tmp_rc_inc_wq = as.numeric(apply(typhoid[1,,1,"RoutineCampaign",1,], 2, sum))/tmp_pop_wq*epipar["reporting"] # /1e5 # *(tmp_pop_wq*tmp_pop_region)
  tmp_rc_inc = sum(typhoid[1,,1,"RoutineCampaign",1,])*epipar["reporting"] # /1e5 # *tmp_pop_region
  
  tmp_rc_avg_age_wq = typhoid[1,,1,"RoutineCampaign",1,] %>% sweep(2, colSums(typhoid[1,,1,"RoutineCampaign",1,]), "/") %>%
    sweep(1, c(4.5/12, 2.75/2, 3.5, 10, 20, (fix$lifexp+25)/2), "*") %>% apply(2, sum) %>% as.numeric
  tmp_rc_avg_age = sum((apply(typhoid[1,,1,"RoutineCampaign",1,], 1, sum) %>% as.numeric)/sum(typhoid[1,,1,"RoutineCampaign",1,])*c(4.5/12, 2.75/2, 3.5, 10, 20, (fix$lifexp+25)/2))
  
  # will have to be MERGED with icers_dalys_costs_all_summary (no WQ) a icers_dalys_costs_summary (w WQ)
  # need 
  tmp_doses_wq = apply(doses, c(1,3:6), sum) %>% cubelyr::as.tbl_cube(met_name="doses") %>% as_tibble() %>% 
    group_by(comp, ns_strat, strat, wealth_quintile) %>% 
    summarize(horizon_10y = sum(doses[year %in% unique(year)[1:10]]),
              horizon_20y = sum(doses)) %>% 
    pivot_longer(cols = starts_with("horizon_"), 
                 names_to = "horizon",
                 names_prefix = c("horizon_"),
                 values_to = "doses") %>% 
    pivot_wider(id_cols = c(horizon, comp, ns_strat, wealth_quintile), 
                names_from = strat, 
                names_glue = "{strat}_{.value}", 
                values_from = doses) 
  tmp_doses = apply(doses, c(1,3:5), sum) %>% cubelyr::as.tbl_cube(met_name="doses") %>% as_tibble() %>% 
    group_by(comp, ns_strat, strat) %>% 
    summarize(horizon_10y = sum(doses[year %in% unique(year)[1:10]]),
              horizon_20y = sum(doses)) %>% 
    pivot_longer(cols = starts_with("horizon_"), 
                 names_to = "horizon",
                 names_prefix = c("horizon_"),
                 values_to = "doses") %>% 
    pivot_wider(id_cols = c(horizon, comp, ns_strat), 
                names_from = strat, 
                names_glue = "{strat}_{.value}", 
                values_from = doses) 
  
  ### Compile together ---------
  icer_results_all[[zs]] = icers_dalys_costs_all_summary %>% left_join(tmp_doses) %>% 
    mutate(gdp=gdp_cap$GDP_cap[gdp_cap$ISO==indicators$ISO[zs]], ISO=indicators$ISO[zs], 
           # wtp_woods = ifelse(length(tmp_woods)==1, tmp_woods, threshold$vce/2), 
           # ochalek_wtp = wtp_ochalek$DALY1_USD2015[wtp_ochalek$ISO==goodruns[zz]], 
           Region=paste0(sprintf("%02d", zs), "_",indicators$RegionLbl[zs]),
           WB_group=countries$WB_group[countries$ISO==indicators$ISO[zs]], 
           Country_archetype=countries$TCVMAP_archetype[countries$ISO==indicators$ISO[zs]], 
           lifexp=fix$lifexp, 
           pop=indicators$pop_DHS_WB_calc[zs],
           meaninc = mean(TransPar$incsamples[,tmp]), # can we get this for the region. What did I calculate they would have?
           meanage = 0.5*(avgage_yale[tmp]+avgage_ihme[tmp]),
           inc_region = tmp_baseline_inc,
           mean_age_region = tmp_avg_age, 
           inc_region_rc = tmp_rc_inc, 
           mean_age_region_rc = tmp_rc_avg_age,
           mean_age_region_wq1 = tmp_avg_age_wq[1],
           death_per_case=avg_deaths,
           mcv1_cov = indicators$mcv1_total[zs],
           mcv1_cov_wq1 = indicators$mcv1_wq1[zs], 
           mcv1_cov_disp = (indicators$mcv1_wq5[zs] - indicators$mcv1_wq1[zs]),
           water_cov = indicators$water_total[zs],
           water_cov_wq1 = indicators$water_wq1[zs],
           sani_cov = indicators$sani_total[zs],
           sani_cov_wq1 = indicators$sani_wq1[zs],
           sani_cov_disp = (indicators$sani_wq5[zs] - indicators$sani_wq1[zs]),
           pop_perc_wq1 = indicators$pop_wq1[zs], 
           pop_perc_wq2 = indicators$pop_wq2[zs], 
           pop_perc_wq3 = indicators$pop_wq3[zs], 
           pov = sum(as.numeric(indicators[zs,paste0("pop_wq", 1:5)])*c(0.1, 0.3, 0.5, 0.7, 0.9)),
           infants_vax_fixedpost=mmgh_data_maps$uc1_del_2037[mmgh_data_maps$ISO==indicators$ISO[zs]]) %>% 
    bind_cols(data.frame(t(unlist(unc))))

  icer_results_wq[[zs]] = icers_dalys_costs_summary %>% left_join(tmp_doses_wq) %>%
    mutate(gdp=gdp_cap$GDP_cap[gdp_cap$ISO==indicators$ISO[zs]], ISO=indicators$ISO[zs], 
           # wtp_woods = ifelse(length(tmp_woods)==1, tmp_woods, threshold$vce/2), 
           # ochalek_wtp = wtp_ochalek$DALY1_USD2015[wtp_ochalek$ISO==goodruns[zz]], 
           Region=paste0(sprintf("%02d", zs), "_",indicators$RegionLbl[zs]),
           WB_group=countries$WB_group[countries$ISO==indicators$ISO[zs]], 
           Country_archetype=countries$TCVMAP_archetype[countries$ISO==indicators$ISO[zs]], 
           lifexp=fix$lifexp, 
           pop=indicators$pop_DHS_WB_calc[zs],
           meaninc = mean(TransPar$incsamples[,tmp]), # can we get this for the region. What did I calculate they would have?
           meanage = 0.5*(avgage_yale[tmp]+avgage_ihme[tmp]),
           pop_perc_wq1 = indicators$pop_wq1[zs], 
           pop_perc_wq2 = indicators$pop_wq2[zs], 
           death_per_case=avg_deaths,
           infants_vax_fixedpost=mmgh_data_maps$uc1_del_2037[mmgh_data_maps$ISO==indicators$ISO[zs]]) %>% 
    bind_cols(data.frame(t(unlist(unc)))) %>% 
    left_join(data.frame(wealth_quintile=paste0("WQ", 1:5), 
                         inc_region=tmp_baseline_inc_wq, 
                         mean_age_region=tmp_avg_age_wq, 
                         inc_region_rc = tmp_rc_inc_wq, 
                         mean_age_region_rc = tmp_rc_avg_age_wq,
                         # R0=parameters$R0,
                         pop_wq=as.numeric(indicators[zs,paste0("pop_wq", 1:5)]), 
                         mcv1=as.numeric(indicators[zs,paste0("mcv1_wq", 1:5)]),
                         water=as.numeric(indicators[zs,paste0("water_wq", 1:5)]),
                         sani=as.numeric(indicators[zs,paste0("sani_wq", 1:5)])))
}

icer_results_all = dplyr::bind_rows(icer_results_all) %>% 
  mutate(CEcat = ifelse(icers<0, "CS",
                        ifelse(icers/gdp<0.5, "HCE",
                               ifelse(icers/gdp<1, "VCE", 
                                      ifelse(icers/gdp<3, "CE", "Not CE"))))) %>% 
  mutate(CEcat = factor(CEcat, levels=c("CS", "HCE", "VCE", "CE", "Not CE"))) %>% 
  mutate(vax_price = paste0("$", format(i, digits=3, nsmall=2), " per dose")) 
  
icer_results_wq = dplyr::bind_rows(icer_results_wq) %>% 
  mutate(CEcat = ifelse(icers<0, "CS", 
                        ifelse(icers/gdp<0.5, "HCE",
                               ifelse(icers/gdp<1, "VCE", ifelse(icers/gdp<3, "CE", "Not CE"))))) %>% 
  mutate(CEcat = factor(CEcat, levels=c("CS", "HCE", "VCE", "CE", "Not CE"))) %>% 
  mutate(vax_price = paste0("$", format(i, digits=3, nsmall=2), " per dose")) 

save(icer_results_all, icer_results_wq, file=paste0("./figures/dist/d", i, "/all/icer_results_data.Rdata"))

# Influential parameters -----
tmp_ind_data = icer_results_all %>% filter(ns_strat=="RoutineCampaign", horizon=="20y", discounting=="disc",
                            maps_profile=="Base", comp=="comparator1") 

tmp_ind_data$poor_prev = tmp_ind_data$pop_perc_wq1+tmp_ind_data$pop_perc_wq2+tmp_ind_data$pop_perc_wq3

lm_bfa = lm(scale(icers) ~ scale(sani_cov) + scale(mcv1_cov) + scale(poor_prev), tmp_ind_data[tmp_ind_data$ISO=="BFA",]) %>% summary
lm_bfa_sani = lm(scale(icers) ~ scale(poor_prev), tmp_ind_data[tmp_ind_data$ISO=="BFA",]) %>% summary
lm_bfa_mcv = lm(scale(icers) ~ scale(sani_cov), tmp_ind_data[tmp_ind_data$ISO=="BFA",]) %>% summary
lm_bfa_pov = lm(scale(icers) ~ scale(mcv1_cov), tmp_ind_data[tmp_ind_data$ISO=="BFA",]) %>% summary

lm_ind = lm(scale(icers) ~ scale(sani_cov) + scale(mcv1_cov) + scale(poor_prev), tmp_ind_data[tmp_ind_data$ISO=="IND",]) %>% summary
lm_ind_sani = lm(scale(icers) ~ scale(sani_cov), tmp_ind_data[tmp_ind_data$ISO=="IND",]) %>% summary
lm_ind_mcv = lm(scale(icers) ~ scale(mcv1_cov), tmp_ind_data[tmp_ind_data$ISO=="IND",]) %>% summary
lm_ind_pov = lm(scale(icers) ~ scale(poor_prev), tmp_ind_data[tmp_ind_data$ISO=="IND",]) %>% summary

lm_ken = lm(scale(icers) ~ scale(sani_cov) + scale(mcv1_cov) + scale(poor_prev), tmp_ind_data[tmp_ind_data$ISO=="KEN",]) %>% summary
lm_ken_sani = lm(scale(icers) ~ scale(sani_cov), tmp_ind_data[tmp_ind_data$ISO=="KEN",]) %>% summary
lm_ken_mcv = lm(scale(icers) ~ scale(mcv1_cov), tmp_ind_data[tmp_ind_data$ISO=="KEN",]) %>% summary
lm_ken_pov = lm(scale(icers) ~ scale(poor_prev), tmp_ind_data[tmp_ind_data$ISO=="KEN",]) %>% summary

lm_mwi = lm(scale(icers) ~ scale(sani_cov) + scale(mcv1_cov) + scale(poor_prev), tmp_ind_data[tmp_ind_data$ISO=="MWI" & tmp_ind_data$icers<Inf,]) %>% summary
lm_mwi_sani = lm(scale(icers) ~ scale(sani_cov), tmp_ind_data[tmp_ind_data$ISO=="MWI" & tmp_ind_data$icers<Inf,]) %>% summary
lm_mwi_mcv = lm(scale(icers) ~ scale(mcv1_cov), tmp_ind_data[tmp_ind_data$ISO=="MWI" & tmp_ind_data$icers<Inf,]) %>% summary
lm_mwi_pov = lm(scale(icers) ~ scale(poor_prev), tmp_ind_data[tmp_ind_data$ISO=="MWI" & tmp_ind_data$icers<Inf,]) %>% summary

lm_npl = lm(scale(icers) ~ scale(sani_cov) + scale(mcv1_cov) + scale(poor_prev), tmp_ind_data[tmp_ind_data$ISO=="NPL",]) %>% summary
lm_npl_sani = lm(scale(icers) ~ scale(sani_cov), tmp_ind_data[tmp_ind_data$ISO=="NPL",]) %>% summary
lm_npl_mcv = lm(scale(icers) ~ scale(mcv1_cov), tmp_ind_data[tmp_ind_data$ISO=="NPL",]) %>% summary
lm_npl_pov = lm(scale(icers) ~ scale(poor_prev), tmp_ind_data[tmp_ind_data$ISO=="NPL",]) %>% summary

isolist = c("BFA", "IND", "KEN", "MWI", "NPL")
arche = c(2,5,2,3,5)
reg_graph_df = data.frame(ISO = rep(isolist, each=6), 
           Archetype = rep(arche, each=6),
           Measure = rep(c(rep("Unadjusted", 3), rep("Adjusted", each=3)), times=5), 
           Par = rep(c("Sanitation", "MCV1", "Poverty"), times=10), 
           Est = c(lm_bfa_sani$coefficients[2,"Estimate"], lm_bfa_mcv$coefficients[2,"Estimate"], lm_bfa_pov$coefficients[2,"Estimate"], lm_bfa$coefficients[2:4,"Estimate"], 
             lm_ind_sani$coefficients[2,"Estimate"], lm_ind_mcv$coefficients[2,"Estimate"], lm_ind_pov$coefficients[2,"Estimate"], lm_ind$coefficients[2:4,"Estimate"], 
             lm_ken_sani$coefficients[2,"Estimate"], lm_ken_mcv$coefficients[2,"Estimate"], lm_ken_pov$coefficients[2,"Estimate"], lm_ken$coefficients[2:4,"Estimate"], 
             lm_mwi_sani$coefficients[2,"Estimate"], lm_mwi_mcv$coefficients[2,"Estimate"], lm_mwi_pov$coefficients[2,"Estimate"], lm_mwi$coefficients[2:4,"Estimate"], 
             lm_npl_sani$coefficients[2,"Estimate"], lm_npl_mcv$coefficients[2,"Estimate"], lm_npl_pov$coefficients[2,"Estimate"], lm_npl$coefficients[2:4,"Estimate"]),
           StD = c(lm_bfa_sani$coefficients[2,"Std. Error"], lm_bfa_mcv$coefficients[2,"Std. Error"], lm_bfa_pov$coefficients[2,"Std. Error"], lm_bfa$coefficients[2:4,"Std. Error"], 
             lm_ind_sani$coefficients[2,"Std. Error"], lm_ind_mcv$coefficients[2,"Std. Error"], lm_ind_pov$coefficients[2,"Std. Error"], lm_ind$coefficients[2:4,"Std. Error"], 
             lm_ken_sani$coefficients[2,"Std. Error"], lm_ken_mcv$coefficients[2,"Std. Error"], lm_ken_pov$coefficients[2,"Std. Error"], lm_ken$coefficients[2:4,"Std. Error"], 
             lm_mwi_sani$coefficients[2,"Std. Error"], lm_mwi_mcv$coefficients[2,"Std. Error"], lm_mwi_pov$coefficients[2,"Std. Error"], lm_mwi$coefficients[2:4,"Std. Error"], 
             lm_npl_sani$coefficients[2,"Std. Error"], lm_npl_mcv$coefficients[2,"Std. Error"], lm_npl_pov$coefficients[2,"Std. Error"], lm_npl$coefficients[2:4,"Std. Error"]))
reg_graph_df$UCI = reg_graph_df$Est + 1.96*reg_graph_df$StD
reg_graph_df$LCI = reg_graph_df$Est - 1.96*reg_graph_df$StD
reg_graph_df$ISO = factor(reg_graph_df$ISO)

ggplot(reg_graph_df, aes(x=Par, y=Est, group=ISO, color=ISO)) + 
  geom_point(position=position_dodge(0.5))+
  geom_errorbar(aes(ymin=LCI, ymax=UCI), width=0.5,
                position=position_dodge(0.5)) + themebar + 
   coord_flip() + 
  geom_hline(yintercept =0) + 
  scale_x_discrete(limits=rev(c("Sanitation", "MCV1", "Poverty"))) + 
  ylab("Estimate (normalized)") + xlab("Parameter") + 
  facet_grid(.~factor(Measure, levels=c('Unadjusted','Adjusted')))

ggsave(paste0("../maps_tcv_country/figures/dist/d", i, "/all/pars_inf.pdf"), width = 6, height = 4, units = "in")

}
