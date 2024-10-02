#############################
# Making Epi output table going into the report ----
# (baseline/default settings only). 
#############################

library(readxl)
library(dplyr)
library(magrittr)
library(tidyr)
library(ggplot2)
library(ggpubr)
library(scales)
library(countrycode)
library(wbstats)
library(flextable)
library(stringr)

fls = list.files("./functions/", pattern="^[fcn]")
for (i in 1:length(fls)){source(paste0("./functions/", fls[i]))}
meanpi = function(x){c(mean=mean(x), lci=quantile(x, 0.025), hci=quantile(x, 0.975))}
strsplit_plus = function(x){strsplit(str_squish(gsub("\\)", "", gsub(",", "", gsub("\\(", " ", x)))), " ")[[1]]}

source("./C_cea_code/00c_load_helper_functions.R")

TransPar = readMat("./data/model_samples_ALL_Aug2021.mat")
CN = readMat("./data/burden_dynamic_link_Yale.mat")$cn
avgage_yale = readMat("./data/burden_dynamic_link_Yale.mat")$avgage.country
avgage_ihme = readMat("./data/burden_dynamic_link_IHME.mat")$avgage.country

agedist = readMat("./data/burden_dynamic_link_Yale.mat")$agedist.country

ISO = read.csv("./data/iso.csv")

# Graphic specifications for ggplot:
ggplot_theme <- theme(axis.text.x = element_text(face="bold", color="black", size=10, angle=0),
                      axis.title.x = element_text(size = 10, angle = 0, face="bold"),
                      axis.text.y = element_text(face="bold", color="black", size=8, angle=0), 
                      axis.title.y = element_text(size = 10, angle = 90, face="bold"),
                      plot.title = element_text(size = 14, angle = 0, face="bold"),
                      plot.subtitle = element_text(size = 14, angle = 0, face="bold"),
                      panel.border = element_rect(linetype = "solid", colour = "black", fill=NA),
                      legend.text = element_text(size = 8, face = "bold", lineheight=0.8, margin = margin(r = 0.3, unit = 'cm')),
                      legend.position = "bottom",
                      legend.box = "vertical",
                      legend.background = element_rect(fill=NA, linewidth =0.25, linetype="solid", colour ="black"),
                      legend.title = element_blank(),
                      legend.key = element_blank(),
                      legend.spacing.x = unit(0,"cm"), 
                      panel.grid.major = element_line(colour="gray", linetype = "dotted"),
                      panel.background = element_rect(fill = NA),
                      strip.background = element_rect(fill = NA),
                      strip.text = element_text(size=12, face="bold")) 

# bring in continent information ----

countries = read_excel("./data/Countries.xlsx", "Country Data w intro date") %>% 
  rename(ISO = `ISO code`, `WB_group` = `WB Group (June 2021)`, 
         `Gavi_eligibility` = `Gavi eligibility`, `TCV_MAPS_intro_date` = `TCV-MAP adoption year`,
         `TCV_NS_intro_date` = `TCV-NS adoption year`, `Typhoid_incidence` = `Typhoid incidence rating`,
         `TCVMAP_archetype`= `TCV-MAP archetype`) %>% 
  mutate(TCV_NS_intro_date = ifelse(TCV_NS_intro_date==0, 2030, TCV_NS_intro_date)) %>% 
  dplyr::select(ISO, Country, WB_group, Gavi_eligibility, Typhoid_incidence, TCV_NS_intro_date, TCV_MAPS_intro_date, TCVMAP_archetype) %>% 
  filter(WB_group!="High income")

# bring in per capita GDP ----

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

# bring in coverage from DHS ----
indicators = read_xlsx("./data/DHS_MICS_GlobalAnalysis_updated.xlsx", "data_to_r")
indicators$ISO = countrycode::countrycode(indicators$Country, "country.name", "iso3c")
indicators$ISO[indicators$Country=="Kosovo"] = "KSV"
indicators$ISO[indicators$Country=="Dominica"] = "DMA" # there was an error with this before 22 May

no_imp_water = read.csv("./data/share-without-improved-water_OWID_JMP.csv") %>% 
  filter(Year==2020)

no_imp_sanitation = read.csv("./data/sanitation_service_level_JMP.csv") %>% 
  filter(Year==2020, Service.level=="Unimproved")

indicators = right_join(indicators, data.frame(ISO=countries$ISO, WB_group=countries$WB_group)) %>% 
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

indicators$mcv1_dif = indicators$mcv1_wq5-indicators$mcv1_wq1
# indicators$mcv1_difa = indicators$mcv1_total-indicators$mcv1_wq1
# indicators$mcv1_difb = indicators$mcv1_wq5-indicators$mcv1_total

indicators = indicators %>% 
  mutate(mcv1_wq1 = ifelse(is.na(mcv1_wq1), pmax(mcv1_total-4, 0.1), mcv1_wq1),
         mcv1_wq2 = ifelse(is.na(mcv1_wq2), pmax(mcv1_total-2, 0.1), mcv1_wq2),
         mcv1_wq3 = ifelse(is.na(mcv1_wq3), pmin(mcv1_total, 99.9), mcv1_wq3),
         mcv1_wq4 = ifelse(is.na(mcv1_wq4), pmin(mcv1_total+2, 99.9), mcv1_wq4),
         mcv1_wq5 = ifelse(is.na(mcv1_wq5), pmin(mcv1_total+4, 99.9), mcv1_wq5)) %>% 
  mutate(across(contains("mcv1_wq"), ~pmax(pmin(.x,99.9),0.1)))

# Create necessary directories
create_dir_if_necessary("./maps_tcv_global/")

# Move to the right directory for outputs
setwd("./maps_tcv_global/")


# Call in all the outputs from the cluster ------
  icer_results_all = list()
  
  for(zz in countries$ISO){
    
    load(paste0("../out_fits/", zz, "/fit_global.Rdata"))
 
    load(paste0('../out_global/0.2/', zz, '/SIR_wq_vax_NOunc_global.Rdata'))
    
    distMAPS = (SIR_wq_vax_unc$simcases_unv_int2 + SIR_wq_vax_unc$simcases_vac_int2) %>% apply(c(1,3), sum) 
    distMAPS = distMAPS/(rowSums(distMAPS) %o% rep(1,5))
    tmpdate = countries$TCV_MAPS_intro_date[countries$ISO==zz] - countries$TCV_NS_intro_date[countries$ISO==zz]
    distMAPS = distMAPS[tmpdate,]
    
    z = which(which(ISO$countryiso==zz)==CN)
    
    input_R0sum = quantile(TransPar$R0samples[,z], c(0.5, 0.025, 0.975)) %>% ci_string_dec(2)
    cov = ginv(SIR_fit$optimfit$hessian)
    tmp = mvrnorm(1000, SIR_fit$optimfit$par[1:2], cov[1:2, 1:2])
    tmpR0 = exp((plogis(tmp[,1])*2.609-1) + (mean((100-SIR_fit$data$haves)/100) %o% (plogis(tmp[,2])*2.303-1)))
    output_R0sum = tmpR0 %>% apply(1, meanpi) %>% apply(2, ci_string_dec, 2)
    
    input_incsum = quantile(TransPar$incsamples[,z], c(0.5, 0.025, 0.975)) %>% ci_string_dec(2)
    output_incsum = sum(rowSums(SIR_fit$simcases)) # no fast way for CIs?
    
    tmp_haves_inc = TransPar$incsamples[,z]*1/2e4*1/rowSums((rep(1, length(SIR_fit$data$oddsrat)) %o% (SIR_fit$data$haves/100)) + (SIR_fit$data$oddsrat %o% (1-SIR_fit$data$haves/100)))
      # mean(TransPar$incsamples[,z])*1/2e4*1/sum(SIR_fit$data$haves/100 + (1-SIR_fit$data$haves/100)*mean(SIR_fit$data$oddsrat))
    
    input_wq_incsum = (((tmp_haves_inc %o% (SIR_fit$data$haves/100)) + (SIR_fit$data$oddsrat*tmp_haves_inc %o% (1-SIR_fit$data$haves/100)))*1e5) %>% apply(2, meanpi) %>% apply(2, ci_string_dec, 2)
    output_wq_incsum = rowSums(SIR_fit$simcases)*5 # no fast way for CIs?
    
    input_dist = SIR_fit$data$dirichletalphas/sum(SIR_fit$data$dirichletalphas) # 'input'
    output_dist = rowSums(SIR_fit$simcases)/sum(rowSums(SIR_fit$simcases)) # 'fitted'
    # chisq.test(rowSums(SIR_fit$simcases), SIR_fit$data$dirichletalphas)
    # not significantly different. So it's fine.
    
    icer_results_all[[zz]] = data.frame(ISO=zz, Country = countries$Country[countries$ISO==zz],
                                        GDP_cap = gdp_cap$GDP_cap[gdp_cap$ISO==zz],
             WB_group=countries$WB_group[countries$ISO==zz], 
             # put in GDP for the sake of graph ordering...
             Country_archetype = countries$TCVMAP_archetype[countries$ISO==zz], 
             Gavi_eligibility = countries$Gavi_eligibility[countries$ISO==zz], 
             
             mcv1_cov = indicators$mcv1_total[indicators$ISO==zz],
             mcv1_cov_wq1 = indicators$mcv1_wq1[indicators$ISO==zz], 
             sani_cov = indicators$sani_total[indicators$ISO==zz],
             sani_cov_wq1 = indicators$sani_wq1[indicators$ISO==zz], 
             
             # Input vs output outcomes
             input_incsum = input_incsum,
             output_incsum = output_incsum,
             input_wq1_incsum = input_wq_incsum[1],
             input_wq2_incsum = input_wq_incsum[2],
             input_wq3_incsum = input_wq_incsum[3],
             input_wq4_incsum = input_wq_incsum[4],
             input_wq5_incsum = input_wq_incsum[5],
             
             output_wq1_incsum = output_wq_incsum[1],
             output_wq2_incsum = output_wq_incsum[2],
             output_wq3_incsum = output_wq_incsum[3],
             output_wq4_incsum = output_wq_incsum[4],
             output_wq5_incsum = output_wq_incsum[5],
             
             input_dist1 = input_dist[1],
             input_dist2 = input_dist[2],
             input_dist3 = input_dist[3],
             input_dist4 = input_dist[4],
             input_dist5 = input_dist[5],
             
             output_dist1 = output_dist[1],
             output_dist2 = output_dist[2],
             output_dist3 = output_dist[3],
             output_dist4 = output_dist[4],
             output_dist5 = output_dist[5],
             
             intro_dist1 = distMAPS[1],
             intro_dist2 = distMAPS[2],
             intro_dist3 = distMAPS[3],
             intro_dist4 = distMAPS[4],
             intro_dist5 = distMAPS[5],
             
             # parameters that were re-calculated
             input_R0sum = input_R0sum, # this one is more important to graph than the others.
             output_R0sum = output_R0sum, # this one is more important to graph than the others.
             input_m1 = (quantile(TransPar$m1, c(0.5, 0.025, 0.975)) %>% ci_string_dec(2)),
             ouput_m1 = SIR_fit$results["m1"],
             input_m2 = (quantile(TransPar$m2, c(0.5, 0.025, 0.975)) %>% ci_string_dec(2)),
             output_m2 = SIR_fit$results["m2"],
             input_rep = (quantile(TransPar$repsamples[,z], c(0.5, 0.025, 0.975)) %>% ci_string_dec(2)),
             output_rep  = SIR_fit$results["rep"],
             input_rC = (quantile(TransPar$rC.sample, c(0.5, 0.025, 0.975)) %>% ci_string_dec(2)),
             output_rC = SIR_fit$results["rC"],

             # additional parameters that were not in the model before
             output_R0 = SIR_fit$results["R0"],
             output_R0m = SIR_fit$results["R0m"],
             output_R0wq1 = SIR_fit$R0_results[1],
             output_R0wq2 = SIR_fit$R0_results[2],
             output_R0wq3 = SIR_fit$R0_results[3],
             output_R0wq4 = SIR_fit$R0_results[4],
             output_R0wq5 = SIR_fit$R0_results[5],

             meanage = SIR_fit$data$mean_age)
    }
  
  icer_results_all = dplyr::bind_rows(icer_results_all) %>% 
        mutate(name_place = paste0(Country, " - ", ISO, " ($", round(GDP_cap,0), ")")) %>% 
        mutate(WB_group=ifelse(WB_group=="Lower middle income", "Lower-middle income", 
                          ifelse(WB_group=="Upper middle income", "Upper-middle income", WB_group)))
  
  write.csv(icer_results_all, file="./figures/EpiBaseline/EpiBaseline_IncidenceModel.csv", row.names = F)
    
  ## A) R0 obs and fit ----
  
  ## divided by background color
  # icer_results_all$input_R0sum, icer_results_all$output_R0sum
  R0_graph.df = data.frame(name_place = icer_results_all$name_place, WB_group = icer_results_all$WB_group, input_R0sum_mean=NA, input_R0sum_lci=NA, input_R0sum_hci=NA,
                           output_R0sum_mean=NA, output_R0sum_lci=NA, output_R0sum_hci=NA)
  
  for (i in 1:133){
  tmp1 = strsplit_plus(icer_results_all$input_R0sum[i])
  tmp2 = strsplit_plus(icer_results_all$output_R0sum[i])
  
  R0_graph.df$input_R0sum_mean[i] = as.numeric(tmp1[1])
  R0_graph.df$input_R0sum_lci[i] = as.numeric(tmp1[2])
  R0_graph.df$input_R0sum_hci[i] = as.numeric(tmp1[3])
  R0_graph.df$output_R0sum_mean[i] = as.numeric(tmp2[1])
  R0_graph.df$output_R0sum_lci[i] = as.numeric(tmp2[2])
  R0_graph.df$output_R0sum_hci[i] = as.numeric(tmp2[3])
  }
  
  # reshape
  R0_graph.df_rs = pivot_longer(R0_graph.df, cols=starts_with(c("input", "output")), names_to = c("type", "junk", "measure"), 
                             names_pattern = "(.*)_(.*)_(.*)") %>% 
    dplyr::select(-junk) %>% 
    pivot_wider(names_from = measure, values_from = value) 
  
  # reorganize by GDP 
  iso_order = icer_results_all$name_place[order(icer_results_all$GDP_cap)] %>% unique
  R0_graph.df_rs = R0_graph.df_rs %>% mutate(name_place = factor(name_place, levels=iso_order)) 

  ggplot() +
    # geom_rect(aes(xmin = 0, xmax = 100, ymin = seq(0.5, 20.5, by=2), ymax = seq(1.5, 22.5, 2)), 
    #           fill = 'gray85', color='gray85', alpha = 0.9) + 
    geom_point(data=R0_graph.df_rs, aes(x = mean, y = factor(name_place), color=type), 
               size = 2.5, position = position_dodge(0.8)) +
    geom_errorbarh(data=R0_graph.df_rs, aes(xmin = lci, xmax = hci, y = factor(name_place), color=type), 
                   linewidth = 1, position = position_dodge(0.8)) +
    # scale_x_continuous(limits = c(0,100), breaks = seq(0,100,10)) +
    ggplot_theme +
    guides(color=guide_legend()) + 
    xlab("") + ylab("") + scale_y_discrete(limits=rev, drop=T) + 
    facet_grid(WB_group~., scales = "free", space="free")
  ggsave("./figures/EpiBaseline/R0_input_output_percountry.pdf", width = 7, height = 20, units="in")
  # This graph didn't go into the paper or the supplement.
  
  ggplot() +
    geom_errorbarh(data=R0_graph.df, aes(xmin = input_R0sum_lci, xmax = input_R0sum_hci, y = output_R0sum_mean),
                   linewidth = 1, color="gray50", alpha=0.5) +
    geom_point(data=R0_graph.df, aes(x = input_R0sum_mean, y = output_R0sum_mean), 
               size = 2.5) +
    scale_x_continuous(limits = c(1,8), breaks = seq(0,12.5,1)) +
    scale_y_continuous(limits = c(1,8), breaks = seq(0,12.5,1)) +
    geom_abline(slope=1, intercept = 0) + 
    ggplot_theme + coord_fixed(ratio=1) + 
    guides(color=guide_legend()) + 
    xlab("Previous model (without wealth quintiles)") + ylab("Re-parameterized model\n(with wealth quintiles)") + 
    facet_grid(.~WB_group)
  
  ggsave("./figures/EpiBaseline/R0_input_output.pdf", width = 7, height = 2.75, units="in")
  # This didn't go into the paper or the supplement.
  
  ## B1) incidence obs and fit ----
  # icer_results_all$input_R0sum, icer_results_all$output_R0sum
  inc_graph.df = data.frame(name_place = icer_results_all$name_place, WB_group = icer_results_all$WB_group, input_incsum_mean=NA, input_incsum_lci=NA, input_incsum_hci=NA,
                            output_incsum_mean=NA) # , output_incsum_lci=NA, output_incsum_hci=NA)
  
  for (i in 1:133){
    tmp1 = strsplit_plus(icer_results_all$input_incsum[i])

    inc_graph.df$input_incsum_mean[i] = as.numeric(tmp1[1])
    inc_graph.df$input_incsum_lci[i] = as.numeric(tmp1[2])
    inc_graph.df$input_incsum_hci[i] = as.numeric(tmp1[3])
    inc_graph.df$output_incsum_mean[i] = icer_results_all$output_incsum[i]
  }
  
  # reshape
  inc_graph.df_rs = pivot_longer(inc_graph.df, cols=starts_with(c("input", "output")), 
                                 names_to = c("type", "junk", "measure"), 
                                names_pattern = "(.*)_(.*)_(.*)") %>% 
    dplyr::select(-junk) %>% 
    pivot_wider(names_from = measure, values_from = value) 
  
  # reorganize by GDP 
  iso_order = icer_results_all$name_place[order(icer_results_all$GDP_cap)] %>% unique
  inc_graph.df_rs = inc_graph.df_rs %>% mutate(name_place = factor(name_place, levels=iso_order)) 
  
  # The above problem looks a lot smaller here - I think that it's because rep is recalculated, 
  # it corrects for movement around previous R0....
  ggplot() +
    # geom_rect(aes(xmin = 0, xmax = 100, ymin = seq(0.5, 20.5, by=2), ymax = seq(1.5, 22.5, 2)), 
    #           fill = 'gray85', color='gray85', alpha = 0.9) + 
    geom_errorbarh(data=inc_graph.df_rs, aes(xmin = lci, xmax = hci, y = factor(name_place), color=type), 
                   size = 1, position = position_dodge(0.8)) +
    geom_point(data=inc_graph.df_rs, aes(x = mean, y = factor(name_place), color=type), 
               size = 2.5, position = position_dodge(0.8)) +
    scale_x_continuous(trans='log10') +
    ggplot_theme +
    guides(color=guide_legend()) + 
    xlab("") + ylab("") + scale_y_discrete(limits=rev, drop=T) + 
    facet_grid(WB_group~., scales = "free", space="free")
  ggsave("./figures/EpiBaseline/Cases_input_output_percountry.pdf", width = 7, height = 20, units="in")
  # This graph didn't go into the paper or the supplement.
  
  ggplot(data=inc_graph.df) +
    geom_errorbarh(aes(xmin = input_incsum_lci, xmax = input_incsum_hci, y = output_incsum_mean),
                   size = 1, color="gray50", alpha=0.5) +
    geom_point(aes(x = input_incsum_mean, y = output_incsum_mean), 
               size = 2.5) +
    scale_x_continuous(trans='log10', limits = c(0.1,1e4)) +
    scale_y_continuous(trans='log10', limits = c(0.1,1e4)) +
    geom_abline(slope=1, intercept = 0) + 
    ggplot_theme + coord_fixed(ratio=1) + 
    guides(color=guide_legend()) + 
    xlab("Previous model (without wealth quintiles)") + 
    ylab("Re-parameterized model\n(with wealth quintiles)") + 
    facet_grid(.~WB_group)
  
  ggsave("./figures/EpiBaseline/Cases_input_output.pdf", width = 7, height = 2.75, units="in")
  
  ## B2) incidence obs and fit, by WQ x WBG ----
    # icer_results_all$input_wq1-5_incsum, icer_results_all$output_wq1-5_incsum
  inc_wq_graph.df = data.frame(name_place = icer_results_all$name_place, WB_group = icer_results_all$WB_group, 
                               input_wq1_incsum_mean=NA, input_wq1_incsum_lci=NA, input_wq1_incsum_hci=NA,
                               input_wq2_incsum_mean=NA, input_wq2_incsum_lci=NA, input_wq2_incsum_hci=NA,
                               input_wq3_incsum_mean=NA, input_wq3_incsum_lci=NA, input_wq3_incsum_hci=NA,
                               input_wq4_incsum_mean=NA, input_wq4_incsum_lci=NA, input_wq4_incsum_hci=NA,
                               input_wq5_incsum_mean=NA, input_wq5_incsum_lci=NA, input_wq5_incsum_hci=NA,
                               output_wq1_incsum_mean=NA, output_wq2_incsum_mean=NA, output_wq3_incsum_mean=NA, 
                               output_wq4_incsum_mean=NA, output_wq5_incsum_mean=NA) # , output_incsum_lci=NA, output_incsum_hci=NA)
  
  for (i in 1:133){
    tmp1 = strsplit_plus(icer_results_all$input_wq1_incsum[i])
    tmp2 = strsplit_plus(icer_results_all$input_wq2_incsum[i])
    tmp3 = strsplit_plus(icer_results_all$input_wq3_incsum[i])
    tmp4 = strsplit_plus(icer_results_all$input_wq4_incsum[i])
    tmp5 = strsplit_plus(icer_results_all$input_wq5_incsum[i])
    
    inc_wq_graph.df$input_wq1_incsum_mean[i] = as.numeric(tmp1[1])
    inc_wq_graph.df$input_wq1_incsum_lci[i] = as.numeric(tmp1[2])
    inc_wq_graph.df$input_wq1_incsum_hci[i] = as.numeric(tmp1[3])
  
    inc_wq_graph.df$input_wq2_incsum_mean[i] = as.numeric(tmp2[1])
    inc_wq_graph.df$input_wq2_incsum_lci[i] = as.numeric(tmp2[2])
    inc_wq_graph.df$input_wq2_incsum_hci[i] = as.numeric(tmp2[3])
    
    inc_wq_graph.df$input_wq3_incsum_mean[i] = as.numeric(tmp3[1])
    inc_wq_graph.df$input_wq3_incsum_lci[i] = as.numeric(tmp3[2])
    inc_wq_graph.df$input_wq3_incsum_hci[i] = as.numeric(tmp3[3])
    
    inc_wq_graph.df$input_wq4_incsum_mean[i] = as.numeric(tmp4[1])
    inc_wq_graph.df$input_wq4_incsum_lci[i] = as.numeric(tmp4[2])
    inc_wq_graph.df$input_wq4_incsum_hci[i] = as.numeric(tmp4[3])
    
    inc_wq_graph.df$input_wq5_incsum_mean[i] = as.numeric(tmp5[1])
    inc_wq_graph.df$input_wq5_incsum_lci[i] = as.numeric(tmp5[2])
    inc_wq_graph.df$input_wq5_incsum_hci[i] = as.numeric(tmp5[3])
    
    inc_wq_graph.df$output_wq1_incsum_mean[i] = icer_results_all$output_wq1_incsum[i]
    inc_wq_graph.df$output_wq2_incsum_mean[i] = icer_results_all$output_wq2_incsum[i]
    inc_wq_graph.df$output_wq3_incsum_mean[i] = icer_results_all$output_wq3_incsum[i]
    inc_wq_graph.df$output_wq4_incsum_mean[i] = icer_results_all$output_wq4_incsum[i]
    inc_wq_graph.df$output_wq5_incsum_mean[i] = icer_results_all$output_wq5_incsum[i]
  }
  
  # reshape
  inc_wq_graph.df_rs = pivot_longer(inc_wq_graph.df, cols=starts_with(c("input", "output")), names_to = c("type", "wq", "junk", "measure"), 
                                 names_pattern = "(.*)_(wq.*)_(.*)_(.*)") %>% 
    filter(measure=="mean") %>% dplyr::select(-junk, -measure) 
  
  # reorganize by GDP 
  iso_order = icer_results_all$name_place[order(icer_results_all$GDP_cap)] %>% unique
  inc_wq_graph.df_rs = inc_wq_graph.df_rs %>% mutate(name_place = factor(name_place, levels=iso_order)) 
  
  # Don't do uncertainty here. Just ranges. 
  # mixed results; low income looks quite good, upper middle is so-so
  ggplot() +
     geom_point(data=inc_wq_graph.df_rs, aes(x = value, y = factor(name_place), color=type), 
               size = 2.5, position = position_dodge(0.8)) +
    scale_x_continuous(trans='log10') +
    ggplot_theme +
    guides(color=guide_legend()) + 
    xlab("") + ylab("") + scale_y_discrete(limits=rev, drop=T) + 
    facet_grid(WB_group~., scales = "free", space="free")
  ggsave("./figures/EpiBaseline/CasesWQ_input_output_percountry.pdf", width = 7, height = 20, units="in")
  
  # Far better..
  inc_wq_graph.df_rs2 = pivot_longer(inc_wq_graph.df, cols=starts_with(c("input", "output")), names_to = c("type", "wq", "junk", "measure"), 
                                    names_pattern = "(.*)_(wq.*)_(.*)_(.*)") %>% 
    dplyr::select(-junk) %>% 
    pivot_wider(names_from = c(type, measure), values_from = value) 
  
  wq_labeller <- c(wq1 = "Poorest", wq2 = "Second", wq3 = "Third", wq4 = "Fourth", wq5 = "Wealthiest")
  
  ggplot() +
    geom_errorbarh(data=inc_wq_graph.df_rs2, aes(xmin = input_lci, xmax = input_hci, y = output_mean),
                   size = 1, color="gray50", alpha=0.5) +
    geom_point(data=inc_wq_graph.df_rs2, aes(x = input_mean, y = output_mean), 
               size = 2.5) + coord_fixed(ratio=1) + 
    scale_x_continuous(trans='log10', limits = c(0.1,1e4)) +
    scale_y_continuous(trans='log10', limits = c(0.1,1e4)) +
    geom_abline(slope=1, intercept = 0) + 
    ggplot_theme +
    guides(color=guide_legend()) + 
    xlab("Target incidence (previous incidence & risk-factor adjustment)") + ylab("Re-parameterized model (with wealth quintiles)") + 
    facet_grid(wq~WB_group, labeller = labeller(wq = wq_labeller))
  
  ggsave("./figures/EpiBaseline/CasesWQ_input_output.pdf", width = 7, height = 10.5, units="in")
  
  ## C) burden distribution ----
    # icer_results_all$input_dist1-5, icer_results_all$output_dist1-5
  
  prop_wq_graph.df = data.frame(name_place = icer_results_all$name_place, WB_group = icer_results_all$WB_group, 
                               input_dist1=NA, input_dist2=NA, input_dist3=NA, input_dist4=NA, input_dist5=NA, 
                               output_dist1=NA, output_dist2=NA, output_dist3=NA, output_dist4=NA, output_dist5=NA,
                               intro_dist1=NA, intro_dist2=NA, intro_dist3=NA, intro_dist4=NA, intro_dist5=NA)
  
  for (i in 1:133){
    prop_wq_graph.df$input_dist1[i] = icer_results_all$input_dist1[i]
    prop_wq_graph.df$input_dist2[i] = icer_results_all$input_dist2[i]
    prop_wq_graph.df$input_dist3[i] = icer_results_all$input_dist3[i]
    prop_wq_graph.df$input_dist4[i] = icer_results_all$input_dist4[i]
    prop_wq_graph.df$input_dist5[i] = icer_results_all$input_dist5[i]

    prop_wq_graph.df$output_dist1[i] = icer_results_all$output_dist1[i]
    prop_wq_graph.df$output_dist2[i] = icer_results_all$output_dist2[i]
    prop_wq_graph.df$output_dist3[i] = icer_results_all$output_dist3[i]
    prop_wq_graph.df$output_dist4[i] = icer_results_all$output_dist4[i]
    prop_wq_graph.df$output_dist5[i] = icer_results_all$output_dist5[i]
    
    prop_wq_graph.df$intro_dist1[i] = icer_results_all$intro_dist1[i]
    prop_wq_graph.df$intro_dist2[i] = icer_results_all$intro_dist2[i]
    prop_wq_graph.df$intro_dist3[i] = icer_results_all$intro_dist3[i]
    prop_wq_graph.df$intro_dist4[i] = icer_results_all$intro_dist4[i]
    prop_wq_graph.df$intro_dist5[i] = icer_results_all$intro_dist5[i]
  }
  
  # reshape
  prop_wq_graph.df = pivot_longer(prop_wq_graph.df, cols=starts_with(c("input", "output", "intro")), names_to = c("type", "junk", "wq"), 
                                    names_pattern = "(.*)_(dist)(.*)") %>% dplyr::select(-junk) %>% 
    mutate(type = factor(type, levels=c("input", "output", "intro"), labels = c("T", "M", "I")))
  # c("Target Dist.", "Model Dist.", "MAP Intro Dist.")
  # reorganize by GDP 
  iso_order = icer_results_all$name_place[order(icer_results_all$GDP_cap)] %>% unique
  prop_wq_graph.df = prop_wq_graph.df %>% mutate(name_place = factor(name_place, levels=iso_order)) 


  disp_fig = function(i){
    prop_wq_graph.df %>%
        filter(WB_group == i) %>% 
        mutate(name_place=gsub(" - ","\n",name_place)) %>% 
        mutate(name_place=ifelse(name_place=="Saint Vincent & Grenadines\nVCT ($9298)", 
                                 "St Vct & Grenadines\nVCT ($9298)", name_place)) %>% 
    ggplot(aes(x=type, y=value, fill=wq)) + 
    geom_bar(stat="identity", position=position_fill(reverse = TRUE)) + 
    geom_abline(intercept = seq(0, 1, 0.20), slope=0) + 
    coord_flip() + 
    scale_fill_manual(values = c('#e66101','#fdb863','#ffffbf','#b2abd2','#5e3c99'),
                      labels = c("Lowest 20%", "Second 20%", "Third 20%", "Fourth 20%", "Highest 20%")) + 
    scale_x_discrete(limits=rev) + scale_y_continuous(breaks=seq(0, 1, 0.20), expand = c(0, 0)) + 
    ggplot_theme + theme(legend.title = element_text(size = 12, face="bold"), 
      legend.text = element_text(size = 12, face = "bold", lineheight=0.8, margin = margin(r = 0.3, unit = 'cm'))) + 
      xlab("") + ylab("Proportion of Cases") + labs(fill="Wealth Quintile") + 
    facet_wrap(name_place~., strip.position = "left", ncol=4) + # switch = "y"
    theme(strip.text.y.left = element_text(angle = 0, size=10, hjust=1), 
          # axis.title.x = element_text(vjust = -10),  # value by experiment
          strip.placement = "outside", 
          panel.spacing = unit(0.25,'lines'))}
  
  tmp_lic = disp_fig("Low income")
  tmp_lmic = disp_fig("Lower-middle income")
  tmp_umic = disp_fig("Upper-middle income")
  
  # arrange two plots into one column
  ggarrange(tmp_lic + ylab("") + ggtitle("Lower income"), 
            tmp_lmic + ylab("") + ggtitle("Lower-middle income"),
            tmp_umic + ggtitle("Upper-middle income"),
    # labels = list("LIC", "LMIC", "UMIC"), 
    ncol = 1, heights = c(0.5, 1, 1), 
    common.legend = T, legend = "bottom", align="v") # %>% 
    # annotate_figure(bottom = text_grob("Abbreviations: 
    #   T: the target distribution of cases among wealth quintiles according to calculations about IRRs of typhoid between people with and without improved sanitation. 
    #   M: the distribution of cases among wealth quintiles yielded by the re-calculated parameters of the typhoid dynamic transmission model. 
    #   I: the distribution of cases among wealth quintiles yielded by the model at the expected time of MAPS initiation.
    #   The vertical lines show the expected distribution under an assumption of perfect equity.", 
    #                                    color = "black", hjust = 0, x=0, face = "bold", size = 12))

  ggsave("./figures/EpiBaseline/BaselineCaseBurdenWQ.pdf", width = 14, height = 17, units="in")
  
  
  # Projections of 2023-2032 -----
  
  # collect cases when no NS been done?
  # collect cases *averted* while I am at it? under both assumptions of R and R+C? Assumptions of MCV coverage...

  wbdata = read_xlsx("~/Library/CloudStorage/Dropbox/My Drive in Dropbox/PATH Microneedle Patches/Model Assumptions/GDP_LE_pop_UN_DHS_combined.xlsx", sheet="Data")
  
  befMAPS_results_all = list()
  
  for(zz in countries$ISO){
    z = which(which(ISO$countryiso==zz)==CN)
    
    tmp_pop = as.numeric(wbdata[wbdata$ISO==ISO[CN[z],2] & wbdata$`Series Name`=="Population, total", "Yr_2021"])
    tmp_pop = ifelse(!is.na(tmp_pop), tmp_pop, as.numeric(wbdata[wbdata$ISO==ISO[CN[z],2] & wbdata$`Series Name`=="Population, total", "Yr_2020"]))
    tmp_pop = ifelse(!is.na(tmp_pop), tmp_pop, as.numeric(wbdata[wbdata$ISO==ISO[CN[z],2] & wbdata$`Series Name`=="Population, total", "Yr_2019"]))
    
    load(paste0("../out_fits/", zz, "/fit_global.Rdata"))
    
    output_incsum = sum(rowSums(SIR_fit$simcases)) # no fast way for CIs?
    rep_rate = as.numeric(strsplit_plus(SIR_fit$results["rep"])[1])
    # rep_rate = mean(TransPar$repsamples[,z])
    
    load(paste0('../out_global/0.2/', zz, '/SIR_wq_vax_NOunc_global.Rdata'))
    Cases_NSper_0_def = SIR_wq_vax_unc$simcases_novac %>% apply(c(1,3), sum) 
    Cases_NSper_R_def = (SIR_wq_vax_unc$simcases_unv_int1 + SIR_wq_vax_unc$simcases_vac_int1) %>% apply(c(1,3), sum) 
    Cases_NSper_RC_def = (SIR_wq_vax_unc$simcases_unv_int2 + SIR_wq_vax_unc$simcases_vac_int2) %>% apply(c(1,3), sum) 
    
    tmpdate = countries$TCV_MAPS_intro_date[countries$ISO==zz] - countries$TCV_NS_intro_date[countries$ISO==zz]
    Cases_NSper_0_def = Cases_NSper_0_def[1:tmpdate,] %>% sum
    Cases_NSper_R_def = Cases_NSper_R_def[1:tmpdate,] %>% sum
    Cases_NSper_RC_def = Cases_NSper_RC_def[1:tmpdate,] %>% sum
    
    ### when MCV is going to improve ----
    load(paste0('../out_global/0.2_mcv1imp/', zz, '/SIR_wq_vax_NOunc_global.Rdata'))
    Cases_NSper_0_mcv1imp = SIR_wq_vax_unc$simcases_novac %>% apply(c(1,3), sum) 
    Cases_NSper_R_mcv1imp = (SIR_wq_vax_unc$simcases_unv_int1 + SIR_wq_vax_unc$simcases_vac_int1) %>% apply(c(1,3), sum) 
    Cases_NSper_RC_mcv1imp = (SIR_wq_vax_unc$simcases_unv_int2 + SIR_wq_vax_unc$simcases_vac_int2) %>% apply(c(1,3), sum) 
    
    Cases_NSper_0_mcv1imp = Cases_NSper_0_mcv1imp[1:tmpdate,] %>% sum
    Cases_NSper_R_mcv1imp = Cases_NSper_R_mcv1imp[1:tmpdate,] %>% sum
    Cases_NSper_RC_mcv1imp = Cases_NSper_RC_mcv1imp[1:tmpdate,] %>% sum
    
    ### when MCV coverage is going to be MCV2 ----
    load(paste0('../out_global/0.2_mcv2/', zz, '/SIR_wq_vax_NOunc_global.Rdata'))
    Cases_NSper_0_mcv2 = SIR_wq_vax_unc$simcases_novac %>% apply(c(1,3), sum) 
    Cases_NSper_R_mcv2 = (SIR_wq_vax_unc$simcases_unv_int1 + SIR_wq_vax_unc$simcases_vac_int1) %>% apply(c(1,3), sum) 
    Cases_NSper_RC_mcv2 = (SIR_wq_vax_unc$simcases_unv_int2 + SIR_wq_vax_unc$simcases_vac_int2) %>% apply(c(1,3), sum) 
    
    Cases_NSper_0_mcv2 = Cases_NSper_0_mcv2[1:tmpdate,] %>% sum
    Cases_NSper_R_mcv2 = Cases_NSper_R_mcv2[1:tmpdate,] %>% sum
    Cases_NSper_RC_mcv2 = Cases_NSper_RC_mcv2[1:tmpdate,] %>% sum
    
    ## Impact of MAPS ----
    # only pulling comp1 because comp2 and 3 have almost the same coverage. 
    load(paste0('../out_global/0.2/', zz, '/wq_maps_NOunc_global_comp1.Rdata'))
    
    # SIR_wq_vax_unc0, SIR_wq_vax_unc1, SIR_wq_vax_unc2, # these are if the MAPS never take place, what happens after the year of MAPS intro
         # SIR_wq_maps_unc0, SIR_wq_maps_unc1, SIR_wq_maps_unc2, # this is if the MAPS are introduced. 
         # the right comparisons are usually if the ones with the same number at the end.
    CasesAv_MAPSper_MAPS_0_def20 = (SIR_wq_vax_unc0$simcases_maps-SIR_wq_maps_unc0$simcases_maps) %>% sum # this is for all 20 years
    CasesAv_MAPSper_MAPS_R_def20 = (SIR_wq_vax_unc1$simcases_maps-SIR_wq_maps_unc1$simcases_maps) %>% sum # this isfor all 20 years
    CasesAv_MAPSper_MAPS_RC_def20 = (SIR_wq_vax_unc2$simcases_maps-SIR_wq_maps_unc2$simcases_maps) %>% sum # this is for all 20 years
    CasesAv_MAPSper_NS_R_def20 = (SIR_wq_vax_unc0$simcases_maps - SIR_wq_vax_unc1$simcases_maps) %>% sum
    CasesAv_MAPSper_NS_RC_def20 = (SIR_wq_vax_unc0$simcases_maps - SIR_wq_vax_unc2$simcases_maps) %>% sum
    
    load(paste0('../out_global/0.1/', zz, '/wq_maps_NOunc_global_comp1.Rdata'))
    
    # SIR_wq_vax_unc0, SIR_wq_vax_unc1, SIR_wq_vax_unc2, # these are if the MAPS never take place, what happens after the year of MAPS intro
    # SIR_wq_maps_unc0, SIR_wq_maps_unc1, SIR_wq_maps_unc2, # this is if the MAPS are introduced. 
    # the right comparisons are usually if the ones with the same number at the end.
    CasesAv_MAPSper_MAPS_0_def10 = (SIR_wq_vax_unc0$simcases_maps-SIR_wq_maps_unc0$simcases_maps) %>% sum # this is for all 20 years
    CasesAv_MAPSper_MAPS_R_def10 = (SIR_wq_vax_unc1$simcases_maps-SIR_wq_maps_unc1$simcases_maps) %>% sum # this isfor all 20 years
    CasesAv_MAPSper_MAPS_RC_def10 = (SIR_wq_vax_unc2$simcases_maps-SIR_wq_maps_unc2$simcases_maps) %>% sum # this is for all 20 years
    CasesAv_MAPSper_NS_R_def10 = (SIR_wq_vax_unc0$simcases_maps - SIR_wq_vax_unc1$simcases_maps) %>% sum
    CasesAv_MAPSper_NS_RC_def10 = (SIR_wq_vax_unc0$simcases_maps - SIR_wq_vax_unc2$simcases_maps) %>% sum
    
    load(paste0('../out_global/0.3/', zz, '/wq_maps_NOunc_global_comp1.Rdata'))
    
    # SIR_wq_vax_unc0, SIR_wq_vax_unc1, SIR_wq_vax_unc2, # these are if the MAPS never take place, what happens after the year of MAPS intro
    # SIR_wq_maps_unc0, SIR_wq_maps_unc1, SIR_wq_maps_unc2, # this is if the MAPS are introduced. 
    # the right comparisons are usually if the ones with the same number at the end.
    CasesAv_MAPSper_MAPS_0_def30 = (SIR_wq_vax_unc0$simcases_maps-SIR_wq_maps_unc0$simcases_maps) %>% sum # this is for all 20 years
    CasesAv_MAPSper_MAPS_R_def30 = (SIR_wq_vax_unc1$simcases_maps-SIR_wq_maps_unc1$simcases_maps) %>% sum # this isfor all 20 years
    CasesAv_MAPSper_MAPS_RC_def30 = (SIR_wq_vax_unc2$simcases_maps-SIR_wq_maps_unc2$simcases_maps) %>% sum # this is for all 20 years
    CasesAv_MAPSper_NS_R_def30 = (SIR_wq_vax_unc0$simcases_maps - SIR_wq_vax_unc1$simcases_maps) %>% sum
    CasesAv_MAPSper_NS_RC_def30 = (SIR_wq_vax_unc0$simcases_maps - SIR_wq_vax_unc2$simcases_maps) %>% sum
    
    # when MCV is going to improve
    load(paste0('../out_global/0.2_mcv1imp/', zz, '/wq_maps_NOunc_global_comp1.Rdata'))
    CasesAv_MAPSper_MAPS_0_mcv1imp = (SIR_wq_vax_unc0$simcases_maps-SIR_wq_maps_unc0$simcases_maps) %>% sum # this is for all 20 years
    CasesAv_MAPSper_MAPS_R_mcv1imp = (SIR_wq_vax_unc1$simcases_maps-SIR_wq_maps_unc1$simcases_maps) %>% sum # this isfor all 20 years
    CasesAv_MAPSper_MAPS_RC_mcv1imp = (SIR_wq_vax_unc2$simcases_maps-SIR_wq_maps_unc2$simcases_maps) %>% sum # this is for all 20 years
    CasesAv_MAPSper_NS_R_mcv1imp = (SIR_wq_vax_unc0$simcases_maps - SIR_wq_vax_unc1$simcases_maps) %>% sum
    CasesAv_MAPSper_NS_RC_mcv1imp = (SIR_wq_vax_unc0$simcases_maps - SIR_wq_vax_unc2$simcases_maps) %>% sum
    
    # when MCV coverage is going to be MCV2
    load(paste0('../out_global/0.2_mcv2/', zz, '/wq_maps_NOunc_global_comp1.Rdata'))
    CasesAv_MAPSper_MAPS_0_mcv2 = (SIR_wq_vax_unc0$simcases_maps-SIR_wq_maps_unc0$simcases_maps) %>% sum # this is for all 20 years
    CasesAv_MAPSper_MAPS_R_mcv2 = (SIR_wq_vax_unc1$simcases_maps-SIR_wq_maps_unc1$simcases_maps) %>% sum # this isfor all 20 years
    CasesAv_MAPSper_MAPS_RC_mcv2 = (SIR_wq_vax_unc2$simcases_maps-SIR_wq_maps_unc2$simcases_maps) %>% sum # this is for all 20 years
    CasesAv_MAPSper_NS_R_mcv2 = (SIR_wq_vax_unc0$simcases_maps - SIR_wq_vax_unc1$simcases_maps) %>% sum
    CasesAv_MAPSper_NS_RC_mcv2 = (SIR_wq_vax_unc0$simcases_maps - SIR_wq_vax_unc2$simcases_maps) %>% sum
    
    befMAPS_results_all[[zz]] = data.frame(ISO=zz, Country = countries$Country[countries$ISO==zz],
                                        GDP_cap = gdp_cap$GDP_cap[gdp_cap$ISO==zz],
                                        WB_group=countries$WB_group[countries$ISO==zz], 
                                        # put in GDP for the sake of graph ordering...
                                        Country_archetype = countries$TCVMAP_archetype[countries$ISO==zz], 
                                        Gavi_eligibility = countries$Gavi_eligibility[countries$ISO==zz], 
                                        pop = tmp_pop,
                                        
                                        mcv1_cov = indicators$mcv1_total[indicators$ISO==zz],
                                        mcv1_cov_wq1 = indicators$mcv1_wq1[indicators$ISO==zz], 
                                        sani_cov = indicators$sani_total[indicators$ISO==zz],
                                        sani_cov_wq1 = indicators$sani_wq1[indicators$ISO==zz], 
                                        
                                        meanage = SIR_fit$data$mean_age,
                                        Cases_NSper_0_def_20 = Cases_NSper_0_def*tmp_pop*rep_rate/1e5,
                                        CasesAv_NSper_NS_R_def_20 = (Cases_NSper_0_def-Cases_NSper_R_def)*tmp_pop*rep_rate/1e5,
                                        CasesAv_NSper_NS_RC_def_20 = (Cases_NSper_0_def-Cases_NSper_RC_def)*tmp_pop*rep_rate/1e5,
                                        
                                        CasesAv_NSper_NS_R_mcv1imp_20 = (Cases_NSper_0_def-Cases_NSper_R_mcv1imp)*tmp_pop*rep_rate/1e5,
                                        CasesAv_NSper_NS_RC_mcv1imp_20 = (Cases_NSper_0_def-Cases_NSper_RC_mcv1imp)*tmp_pop*rep_rate/1e5,
                                        
                                        CasesAv_NSper_NS_R_mcv2_20 = (Cases_NSper_0_def-Cases_NSper_R_mcv2)*tmp_pop*rep_rate/1e5,
                                        CasesAv_NSper_NS_RC_mcv2_20 = (Cases_NSper_0_def-Cases_NSper_RC_mcv2)*tmp_pop*rep_rate/1e5,
                                        
                                        Cases_MAPSper_0_def = (SIR_wq_vax_unc0$simcases_maps %>% sum)*tmp_pop*rep_rate/1e5,
                                          
                                        CasesAv_MAPSper_MAPS_0_def_20 = CasesAv_MAPSper_MAPS_0_def20*tmp_pop*rep_rate/1e5,
                                        CasesAv_MAPSper_MAPS_R_def_20 = CasesAv_MAPSper_MAPS_R_def20*tmp_pop*rep_rate/1e5,
                                        CasesAv_MAPSper_MAPS_RC_def_20 = CasesAv_MAPSper_MAPS_RC_def20*tmp_pop*rep_rate/1e5,
                                        CasesAv_MAPSper_NS_R_def_20 = CasesAv_MAPSper_NS_R_def20*tmp_pop*rep_rate/1e5,
                                        CasesAv_MAPSper_NS_RC_def_20 = CasesAv_MAPSper_NS_RC_def20*tmp_pop*rep_rate/1e5,
                                        
                                        CasesAv_MAPSper_MAPS_0_def_10 = CasesAv_MAPSper_MAPS_0_def10*tmp_pop*rep_rate/1e5,
                                        CasesAv_MAPSper_MAPS_R_def_10 = CasesAv_MAPSper_MAPS_R_def10*tmp_pop*rep_rate/1e5,
                                        CasesAv_MAPSper_MAPS_RC_def_10 = CasesAv_MAPSper_MAPS_RC_def10*tmp_pop*rep_rate/1e5,
                                        CasesAv_MAPSper_NS_R_def_10 = CasesAv_MAPSper_NS_R_def10*tmp_pop*rep_rate/1e5,
                                        CasesAv_MAPSper_NS_RC_def_10 = CasesAv_MAPSper_NS_RC_def10*tmp_pop*rep_rate/1e5,
                                        
                                        CasesAv_MAPSper_MAPS_0_def_30 = CasesAv_MAPSper_MAPS_0_def30*tmp_pop*rep_rate/1e5,
                                        CasesAv_MAPSper_MAPS_R_def_30 = CasesAv_MAPSper_MAPS_R_def30*tmp_pop*rep_rate/1e5,
                                        CasesAv_MAPSper_MAPS_RC_def_30 = CasesAv_MAPSper_MAPS_RC_def30*tmp_pop*rep_rate/1e5,
                                        CasesAv_MAPSper_NS_R_def_30 = CasesAv_MAPSper_NS_R_def30*tmp_pop*rep_rate/1e5,
                                        CasesAv_MAPSper_NS_RC_def_30 = CasesAv_MAPSper_NS_RC_def30*tmp_pop*rep_rate/1e5,
                                        
                                        CasesAv_MAPSper_MAPS_0_mcv1imp_20 = CasesAv_MAPSper_MAPS_0_mcv1imp*tmp_pop*rep_rate/1e5,
                                        CasesAv_MAPSper_MAPS_R_mcv1imp_20 = CasesAv_MAPSper_MAPS_R_mcv1imp*tmp_pop*rep_rate/1e5,
                                        CasesAv_MAPSper_MAPS_RC_mcv1imp_20 = CasesAv_MAPSper_MAPS_RC_mcv1imp*tmp_pop*rep_rate/1e5,
                                        CasesAv_MAPSper_NS_R_mcv1imp_20 = CasesAv_MAPSper_NS_R_mcv1imp*tmp_pop*rep_rate/1e5,
                                        CasesAv_MAPSper_NS_RC_mcv1imp_20 = CasesAv_MAPSper_NS_RC_mcv1imp*tmp_pop*rep_rate/1e5,
                                        
                                        CasesAv_MAPSper_MAPS_0_mcv2_20 = CasesAv_MAPSper_MAPS_0_mcv2*tmp_pop*rep_rate/1e5,
                                        CasesAv_MAPSper_MAPS_R_mcv2_20 = CasesAv_MAPSper_MAPS_R_mcv2*tmp_pop*rep_rate/1e5,
                                        CasesAv_MAPSper_MAPS_RC_mcv2_20 = CasesAv_MAPSper_MAPS_RC_mcv2*tmp_pop*rep_rate/1e5,
                                        CasesAv_MAPSper_NS_R_mcv2_20 = CasesAv_MAPSper_NS_R_mcv2*tmp_pop*rep_rate/1e5,
                                        CasesAv_MAPSper_NS_RC_mcv2_20 = CasesAv_MAPSper_NS_RC_mcv2*tmp_pop*rep_rate/1e5)
  }
  
  befMAPS_results_all = dplyr::bind_rows(befMAPS_results_all) %>% as_tibble() %>% 
    mutate(name_place = paste0(Country, " - ", ISO, " ($", round(GDP_cap,0), ")"))
  
  Cases_totals = befMAPS_results_all %>% 
    dplyr::select(-starts_with("CasesAv")) %>% 
    pivot_longer(cols=starts_with("Cases_"), names_to = "Period", 
                 names_pattern = "Cases_(.*)") %>% 
    group_by(Period) %>% 
    dplyr::summarise(Totals = paste0(round(sum(value, na.rm=T)/1e6, 2), "M"))
    
  CasesAv_totals = befMAPS_results_all %>% 
    pivot_longer(cols=starts_with("CasesAv_"), names_to = c("Period", "Vax", "NSdepl", "MCVcov", "MAPScov"), 
                 names_pattern = "CasesAv_(.*)per_(.*)_(.*)_(.*)_(.*)") %>% 
    group_by(Period, Vax, NSdepl, MCVcov, MAPScov) %>% 
    dplyr::summarise(Totals = sum(value, na.rm=T)) %>% 
    mutate(Totals = paste0(round(Totals/1e6, 2), "M")) %>% 
    pivot_wider(names_from = Vax, values_from = Totals) %>% 
    relocate(Period, MCVcov, NSdepl, MAPScov, NS, MAPS) %>% 
    mutate(MAPScov = ifelse(Period=="NS", NA, MAPScov)) %>% 
    mutate(MCVcov = factor(MCVcov, levels=c("def", "mcv1imp", "mcv2"), 
                           labels=c("MCV1 (default)", "MCV1 improved (sensitivity)", "MCV2 (sensitivity)")),
           MAPScov = factor(MAPScov, levels=c("20", "10", "30"), labels=c("20% (default)", "10% (sensitivity)", "30% (sensitivity)")),
           NSdepl = factor(NSdepl, levels=c("RC", "R", "0"), labels=c("Routine & Campaign (default)", "Routine (sensitivity)", "None (sensitivity)"))) %>% 
    arrange(desc(Period), MCVcov, NSdepl, MAPScov)
    
  
    ft = flextable(CasesAv_totals) %>%
    set_header_labels(`Period` = "Period", `MCVcov` = "TCV NS coverage", 
                      `NSdepl` = "TCV NS deployment strategy", `MAPScov` = "TCV MAPS coverage of unvaccinated", 
                      `NS` = "NS",  `MAPS` = "MAPS") %>% 
      merge_v(j=c("Period", "MCVcov", "NSdepl", "MAPScov", "NS", "MAPS")) %>% 
      theme_vanilla() %>% 
          valign(j = "Period", valign = "top") %>% 
      width(j="Period", 1.25) %>% width(j="NSdepl", 2.25) %>% width(j="MAPScov", 1.5) %>% align(j=c("MAPS", "NS"), align="right") %>% 
      add_header_row(colwidths = c(4, 2), values = c("", "Cases averted")) %>% bold(bold = TRUE, part = "header") %>% 
      append_chunks(i = c(1, 7), j = 1, as_chunk(paste0("\n(No TCV = ", Cases_totals$Totals[c(2,1)], " cases)"))) %>% 
      add_footer_lines(values=paste("The Period NS refers to the period when only a needle & syringe",  
                       "presentation are available, roughly from NS deployment (2023-2032) until MAPS deployment (after 2032).",
                       "The Period MAPS refers to the period when both NS and MAPS presentations are available, after 2032.",
                       "The cases averted during the NS period are counted over a span of years that varies by country bounded by",  
                       "the years of NS deployment and MAPS deployment. The cases averted during the MAPS period are counted",  
                       "for 20 years after MAPS deployment")) 
    
    save_as_docx(`Epidemiology outputs` = ft, path = "./figures/EpiBaseline/EpiTable.docx")
  
