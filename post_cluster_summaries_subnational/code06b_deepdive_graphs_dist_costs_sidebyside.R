#############################
# Making graphs for the manuscript
# Subnational analysis, costs together 
# (run first code06a)
#############################

# load("../out_cluster/ISO_which_to_run_cea.Rdata")
source("./C_cea_code/00b_load_packages.R")
library(readxl)
library(dplyr)
library(magrittr)
library(tidyr)
library(ggplot2)
library(ggrepel)
library(cowplot)
library(patchwork)
library(scales)
library(flextable)
library(stringr)
library(sf)
source("./C_cea_code/00c_load_helper_functions.R")

setwd("./maps_tcv_country")
load("./figures/dist/d2.25/all/icer_results_data.Rdata")
icer_results_all_225 = icer_results_all
load("./figures/dist/d3/all/icer_results_data.Rdata")
icer_results_all_3 = icer_results_all

icer_results_all = bind_rows(icer_results_all_225, icer_results_all_3)

custom_labels <- c("BFA" = "Burkina Faso", "IND" = "India", "KEN" = "Kenya",
                   "MWI" = "Malawi", "NPL" = "Nepal")

pwei = icer_results_all %>% filter(ns_strat=="RoutineCampaign", horizon=="20y", discounting=="disc",
                                   comp=="comparator1", maps_profile!="Pessimistic") %>% # , maps_cov_sens==0.2, !is.na(icers)
  group_by(ISO, vax_price, maps_profile, CEcat) %>% 
  summarise(tot_pop = sum(pop)) %>%
  mutate(xax = paste0(maps_profile, "\n", vax_price)) %>%
  ggplot(aes(x=xax, y=tot_pop,fill=CEcat)) + 
  geom_bar(position = position_fill(reverse = TRUE), stat="identity") + 
  xlab(NULL) + ylab("Percent of countries (weighted)") + 
  scale_y_continuous(labels = scales::percent) + 
  scale_fill_manual(values=c("#006837", "#66bd63", "#addd8e", "#d9ef8b",  "#d73027"),
                    label=c("Cost-saving", "Highly cost-effective", "Very cost-effective", "Cost-effective", "Not cost-effective"),
                    drop=F) + scale_x_discrete(labels = function(x) str_wrap(x, width = 8)) + 
  facet_wrap(.~ISO, labeller = labeller(ISO = custom_labels)) + 
  themebar + guides(fill=guide_legend(nrow=1)) + theme(legend.text = element_text(size=9))

pllegend <- get_legend(
  pwei + 
    guides(fill=guide_legend(nrow=5)) + 
    theme(legend.text = element_text(size=9))
)

pwei + theme(legend.position = 'none') + 
  draw_grob(grob = rectGrob(gp = gpar(col = "black", fill = NA, lwd = 4)),
    x = 1.55, y = 0, width = 0.9, height = 1)  + 
  inset_element(pllegend, 0.725, 0.05, 0.95, 0.35) 
  
ggsave("./figures/dist/d3/BasicSens_cea_weighted.pdf", width = 10, height = 6, units = "in")

# MAPS ----- 
## IND ------
india = st_read("../geographic_maps/dhs_maps/IND_subnational_boundaries/shps/sdr_subnational_boundaries2simple.shp")
# ggplot(india) + geom_sf() + geom_sf_label(aes(label = ADM2_EN))

india$RegionLbl = tolower(gsub(" ", "", gsub(" & ", "", india$DHSREGEN)))

icer_results_all$RegionLbl = unname(sapply(icer_results_all$Region, "delete_prefix", "_"))
icer_results_all$RegionLbl[icer_results_all$RegionLbl=="leh(ladakh)"] = "leh"
icer_results_all$RegionLbl[icer_results_all$RegionLbl=="mahrajganj"] = "maharajganj"
icer_results_all$RegionLbl[icer_results_all$RegionLbl=="santravidasnagar(bhadohi)"] = "santravidasnagar"
icer_results_all$RegionLbl[icer_results_all$RegionLbl=="buxar"] = "buxer"

india = india %>% right_join((icer_results_all %>% filter(ISO=="IND")), by="RegionLbl")

india_map = india %>% filter(ns_strat=="RoutineCampaign", horizon=="20y", discounting=="disc", !is.na(icers),
                 comp=="comparator1", maps_profile!="Pessimistic") %>%
  ggplot() +
  geom_sf(aes(fill = CEcat)) +
  scale_fill_manual(values=c("#006837", "#66bd63", "#addd8e", "#d9ef8b",  "#d73027"), drop = FALSE) +
  facet_grid(vax_price~maps_profile) + themebar

ggdraw(india_map, xlim = c(0, 1), ylim = c(0, 1), clip = "off") +
  draw_grob(grob = rectGrob(gp = gpar(col = "black", fill = NA, lwd = 4)),
    x = 0.0725, y = 0.145, width = 0.425, height = 0.38) 

ggsave("./figures/dist/INDSens_225d_3d_maps.pdf", width = 6, height = 7, units = "in")

## KEN ------
kenya = st_read("../geographic_maps/dhs_maps/KEN_subnational_boundaries/shps/sdr_subnational_boundaries2.shp")
# ggplot(kenya) + geom_sf()

kenya = kenya %>% 
  mutate(RegionLbl = tolower(DHSREGEN)) %>%
  mutate(RegionLbl=ifelse(RegionLbl=="trans-nzoia", "trans nzoia", RegionLbl)) %>% 
  mutate(RegionLbl=ifelse(RegionLbl=="elgeyo marakwet", "elgeyo-marakwet", RegionLbl)) 
kenya$RegionLbl = gsub(" ", "", gsub(" & ", "", kenya$RegionLbl))

# ggplot(kenya) + geom_sf() +  geom_sf_label(aes(label = RegionLbl))

icer_results_all$RegionLbl = unname(sapply(icer_results_all$Region, "delete_prefix", "_"))

kenya = kenya %>% right_join((icer_results_all %>% filter(ISO=="KEN")), by="RegionLbl") 

kenya_map = kenya %>% filter(ns_strat=="RoutineCampaign", horizon=="20y", discounting=="disc", 
                 comp=="comparator1", maps_profile!="Pessimistic") %>%
  ggplot() +
  geom_sf(aes(fill = CEcat)) +
  scale_fill_manual(values=c("#006837", "#66bd63", "#addd8e", "#d9ef8b",  "#d73027"), drop = FALSE) +
  facet_grid(vax_price~maps_profile) + themebar

ggdraw(kenya_map, xlim = c(0, 1), ylim = c(0, 1), clip = "off") +
  draw_grob(grob = rectGrob(gp = gpar(col = "black", fill = NA, lwd = 4)),
    x = 0.0975, y = 0.115, width = 0.395, height = 0.4125) 

ggsave("./figures/dist/KENSens_225d_3d_maps.pdf", width = 6, height = 7, units = "in")

## MWI ------

malawi = st_read("../geographic_maps/mwi_adm_nso_hotosm_20230405_shp/mwi_admbnda_adm2_nso_hotosm_20230405.shp")
# ggplot(malawi) + geom_sf() + geom_sf_label(aes(label = ADM2_EN))

malawi$RegionLbl = gsub(" ", "", tolower(malawi$ADM2_EN))

icer_results_all$RegionLbl = unname(sapply(icer_results_all$Region, "delete_prefix", "_")) 
icer_results_all$RegionLbl = gsub("rural","",icer_results_all$RegionLbl)

malawi = malawi %>% right_join((icer_results_all %>% filter(ISO=="MWI")), by="RegionLbl") 

malawi_map = malawi %>% filter(ns_strat=="RoutineCampaign", horizon=="20y", discounting=="disc", 
                 comp=="comparator1", maps_profile!="Pessimistic") %>%
  ggplot() +
  geom_sf(aes(fill = CEcat)) +
  scale_x_continuous(breaks = seq(32, 36, by = 1)) +  # Custom longitude breaks
  scale_fill_manual(values=c("#006837", "#66bd63", "#addd8e", "#d9ef8b",  "#d73027"), drop = FALSE) +
  facet_grid(vax_price~maps_profile) + themebar

ggdraw(malawi_map, xlim = c(0, 1), ylim = c(0, 1), clip = "off") +
  draw_grob(grob = rectGrob(gp = gpar(col = "black", fill = NA, lwd = 4)),
    x = 0.2, y = 0.115, width = 0.3, height = 0.4125) 

ggsave("./figures/dist/MWISens_225d_3d_maps.pdf", width = 4, height = 7, units = "in")

## NPL ------

nepal = st_read("../geographic_maps/dhs_maps/Nepal_uscb_2015-2030_sdr/Nepal_adm2_uscb_2023.shp")

nepal$RegionLbl =  gsub(" ", "", tolower(nepal$AREA_NAME))

icer_results_all$RegionLbl = unname(sapply(icer_results_all$Region, "delete_prefix", "_")) 
icer_results_all$RegionLbl = gsub("dhanusha", "dhanusa", icer_results_all$RegionLbl)

tmp = icer_results_all %>% filter(ISO=="NPL") %>% 
  filter(ns_strat=="RoutineCampaign", horizon=="20y", discounting=="disc") # , !is.na(icers)

nepal = nepal %>% right_join(tmp, by="RegionLbl") 

nepal_map = nepal %>% filter(ns_strat=="RoutineCampaign", horizon=="20y", discounting=="disc", 
                  comp=="comparator1", maps_profile!="Pessimistic") %>%
  ggplot() +
  geom_sf(aes(fill = CEcat)) +
  scale_fill_manual(values=c("#006837", "#66bd63", "#addd8e", "#d9ef8b",  "#d73027"), drop = FALSE) +
  facet_grid(vax_price~maps_profile) + themebar

ggdraw(nepal_map, xlim = c(0, 1), ylim = c(0, 1), clip = "off") +
  draw_grob(grob = rectGrob(gp = gpar(col = "black", fill = NA, lwd = 4)),
            x = 0.1, y = 0.2, width = 0.4, height = 0.35) 

ggsave("./figures/dist/NPLSens_225d_3d_maps.pdf", width = 6, height = 4, units = "in")

## BFA -----
burkina = st_read("../geographic_maps/dhs_maps/BFA_subnational_boundaries/shps/sdr_subnational_boundaries.shp")
# ggplot(burkina) + geom_sf()

burkina$RegionLbl = gsub(" ", "", gsub(" & ", "", burkina$DHSREGEN))
burkina$RegionLbl[burkina$RegionLbl=="BoucleduMouhoun"] = "BoucledeMouhoun"
burkina$RegionLbl[burkina$RegionLbl=="Hauts-Bassins"] = "HautsBasins"

icer_results_all$RegionLbl = unname(sapply(icer_results_all$Region, "delete_prefix", "_"))

burkina = burkina %>% right_join((icer_results_all %>% filter(ISO=="BFA")), by="RegionLbl")

burkina_map = burkina %>% filter(ns_strat=="RoutineCampaign", horizon=="20y", discounting=="disc", 
                 comp=="comparator1", maps_profile!="Pessimistic") %>%
  ggplot() +
  geom_sf(aes(fill = CEcat)) +
  scale_fill_manual(values=c("#006837", "#66bd63", "#addd8e", "#d9ef8b",  "#d73027"), drop = FALSE) +
  facet_grid(vax_price~maps_profile) + themebar

ggdraw(burkina_map, xlim = c(0, 1), ylim = c(0, 1), clip = "off") +
  draw_grob(grob = rectGrob(gp = gpar(col = "black", fill = NA, lwd = 4)),
            x = 0.075, y = 0.165, width = 0.425, height = 0.375) 

ggsave("./figures/dist/BFASens_225d_3d_maps.pdf", width = 6, height = 5, units = "in")

# Basic results (sensitivity analysis below) ------
tmp = icer_results_all %>% dplyr::filter(horizon=="20y" & comp=="comparator1", 
                                         maps_profile=="Base" & ns_strat=="RoutineCampaign", 
                                         discounting=="disc", vax_price=="$3.00 per dose") %>% 
  dplyr::select("ISO", "RegionLbl", "comp", "maps_profile", "Cases_averted", 
                "Deaths_averted", "DALYs_averted", "Cost_dif", "icers", "gdp", "pop", "CEcat") %>% 
  mutate(Cases_averted = Cases_averted*pop/1e5,
         Deaths_averted = Deaths_averted*pop/1e5,
         DALYs_averted = DALYs_averted*pop/1e5,
         Cost_dif = Cost_dif*pop/1e5) %>% 
  group_by(ISO, CEcat) %>% 
  summarize(Cases_averted = sum(Cases_averted),
            Deaths_averted = sum(Deaths_averted),
            DALYs_averted = sum(DALYs_averted),
            Cost_dif = sum(Cost_dif),
            ICER = sum(Cost_dif)/sum(DALYs_averted),
            ICERscaled = sum(Cost_dif)/sum(DALYs_averted)/mean(gdp),
            gdp=mean(gdp),
            pop=sum(pop),
            N=n()) %>% 
  mutate(CS = ifelse(CEcat == "CS", "CS", "no"),
         HCE = ifelse(CEcat %in% c("CS", "HCE"), "HCE", "no"),
         VCE = ifelse(CEcat %in% c("CS", "HCE", "VCE"), "VCE", "no"),
         CE = ifelse(CEcat %in% c("CS", "HCE", "VCE", "CE"), "CE", "no"))

tmp_hce = tmp %>% 
  group_by(ISO, HCE) %>% 
  summarize(Cases_averted = sum(Cases_averted),
            Deaths_averted = sum(Deaths_averted),
            DALYs_averted = sum(DALYs_averted),
            Cost_dif = sum(Cost_dif),
            ICER = sum(Cost_dif)/sum(DALYs_averted),
            ICERscaled = sum(Cost_dif)/sum(DALYs_averted)/mean(gdp), 
            N = sum(N)) %>% 
  filter(HCE == "HCE") 

tmp_vce = tmp %>% 
  group_by(ISO, VCE) %>% 
  summarize(Cases_averted = sum(Cases_averted),
            Deaths_averted = sum(Deaths_averted),
            DALYs_averted = sum(DALYs_averted),
            Cost_dif = sum(Cost_dif),
            ICER = sum(Cost_dif)/sum(DALYs_averted),
            ICERscaled = sum(Cost_dif)/sum(DALYs_averted)/mean(gdp), 
            N = sum(N)) %>% 
  filter(VCE=="VCE")

tmp_ce = tmp %>% 
  group_by(ISO, CE) %>%
  summarize(Cases_averted = sum(Cases_averted),
            Deaths_averted = sum(Deaths_averted),
            DALYs_averted = sum(DALYs_averted),
            Cost_dif = sum(Cost_dif),
            ICER = sum(Cost_dif)/sum(DALYs_averted),
            ICERscaled = sum(Cost_dif)/sum(DALYs_averted)/mean(gdp),
            N = sum(N)) %>%
  filter(CE=="CE")

tmp %>% 
  group_by(ISO) %>%
  summarize(Cases_averted = sum(Cases_averted),
            Deaths_averted = sum(Deaths_averted),
            DALYs_averted = sum(DALYs_averted),
            Cost_dif = sum(Cost_dif),
            ICER = sum(Cost_dif)/sum(DALYs_averted),
            ICERscaled = sum(Cost_dif)/sum(DALYs_averted)/mean(gdp),
            N = sum(N))

tmp_all = tmp %>% 
  group_by(ISO) %>% 
  summarize(Cases_averted = sum(Cases_averted),
            Deaths_averted = sum(Deaths_averted),
            DALYs_averted = sum(DALYs_averted),
            Cost_dif = sum(Cost_dif),
            ICER = sum(Cost_dif)/sum(DALYs_averted),
            ICERscaled = sum(Cost_dif)/sum(DALYs_averted)/mean(gdp),
            pop=sum(pop), 
            N = sum(N)) %>% 
  dplyr::select(ISO, Cases_averted, Deaths_averted, DALYs_averted, Cost_dif, ICER, ICERscaled) %>% t

# [,1]        [,2]        [,3]        [,4]        [,5]       
# ISO            "BFA"       "IND"       "KEN"       "MWI"       "NPL"      
# Cases_averted  " 19743"    "571400"    " 37690"    "  6941"    "  9393"   
# Deaths_averted "171.19"    "781.89"    "878.04"    "181.48"    " 12.54"   
# DALYs_averted  " 3403.8"   "17713.5"   "17383.9"   " 3546.9"   "  287.1"  
# Cost_dif       " 14283349" "662227996" " 45792460" " 14532966" " 16452611"
# ICER           " 4196"     "37385"     " 2634"     " 4097"     "57300"    
# ICERscaled     " 5.056"    "15.507"    " 1.255"    " 6.351"    "42.871"   

# case_av_dif = 1- c(19/22, 571/726, 37/41, 69/84, 9.3/14)
# 0.34483 0.23356 0.02632 0.15854 0.30597 # unequal pops, 2020 pop for national
# 0.13636 0.21350 0.09756 0.17857 0.33571 # equal pops

tmp2 = icer_results_all %>% dplyr::filter(horizon=="20y" & comp=="comparator1", 
                                          maps_profile!="Pessimistic" & ns_strat == "RoutineCampaign", 
                                          discounting=="disc") %>% 
  dplyr::select("ISO", "RegionLbl", "comp", "maps_profile", "Cases_averted", "vax_price",
                "Deaths_averted", "DALYs_averted", "Cost_dif", "icers", "gdp", "pop", "CEcat") %>% 
  mutate(Cases_averted = Cases_averted*pop/1e5,
         Deaths_averted = Deaths_averted*pop/1e5,
         DALYs_averted = DALYs_averted*pop/1e5,
         Cost_dif = Cost_dif*pop/1e5) %>% 
  group_by(ISO, CEcat, comp, maps_profile, vax_price) %>% 
  summarize(Cases_averted = sum(Cases_averted),
            Deaths_averted = sum(Deaths_averted),
            DALYs_averted = sum(DALYs_averted),
            Cost_dif = sum(Cost_dif),
            ICER = sum(Cost_dif)/sum(DALYs_averted),
            ICERscaled = sum(Cost_dif)/sum(DALYs_averted)/mean(gdp),
            gdp=mean(gdp),
            pop=sum(pop), 
            N = n()) %>% 
  mutate(CS = ifelse(CEcat == "CS", "CS", "no"),
         HCE = ifelse(CEcat %in% c("CS", "HCE"), "HCE", "no"),
         VCE = ifelse(CEcat %in% c("CS", "HCE", "VCE"), "VCE", "no"),
         CE = ifelse(CEcat %in% c("CS", "HCE", "VCE", "CE"), "CE", "no"))

tmp_hce2 = tmp2 %>% 
  group_by(ISO, HCE, comp, maps_profile, vax_price) %>% 
  summarize(Cases_averted = sum(Cases_averted),
            Deaths_averted = sum(Deaths_averted),
            DALYs_averted = sum(DALYs_averted),
            Cost_dif = sum(Cost_dif),
            ICER = sum(Cost_dif)/sum(DALYs_averted),
            ICERscaled = sum(Cost_dif)/sum(DALYs_averted)/mean(gdp), 
            N = sum(N)) %>% 
  filter(HCE == "HCE") %>% rename(policy = HCE)

tmp_vce2 = tmp2 %>% 
  group_by(ISO, VCE, comp, maps_profile, vax_price) %>% 
  summarize(Cases_averted = sum(Cases_averted),
            Deaths_averted = sum(Deaths_averted),
            DALYs_averted = sum(DALYs_averted),
            Cost_dif = sum(Cost_dif),
            ICER = sum(Cost_dif)/sum(DALYs_averted),
            ICERscaled = sum(Cost_dif)/sum(DALYs_averted)/mean(gdp), 
            N = sum(N)) %>% 
  filter(VCE=="VCE") %>% rename(policy = VCE)

tmp_ce2 = tmp2 %>% 
  group_by(ISO, CE, comp, maps_profile, vax_price) %>% 
  summarize(Cases_averted = sum(Cases_averted),
            Deaths_averted = sum(Deaths_averted),
            DALYs_averted = sum(DALYs_averted),
            Cost_dif = sum(Cost_dif),
            ICER = sum(Cost_dif)/sum(DALYs_averted),
            ICERscaled = sum(Cost_dif)/sum(DALYs_averted)/mean(gdp), 
            N = sum(N)) %>% 
  filter(CE=="CE") %>% rename(policy = CE) 

tmp_condense = bind_rows(tmp_hce2, tmp_vce2, tmp_ce2) %>% 
  dplyr::select(ISO, policy, comp, maps_profile, vax_price, N) %>% 
  pivot_wider(id_cols = c(ISO, comp, maps_profile, vax_price), 
              names_from = policy,
              values_from = N) %>% 
  mutate(HCE = ifelse(is.na(HCE), 0, HCE),
         VCE = ifelse(is.na(VCE), 0, VCE),
         CE = ifelse(is.na(CE), 0, CE)) %>% 
  arrange(ISO,  maps_profile, vax_price) %>% 
  # filter(maps_profile=="Optimistic") %>% 
  ungroup %>% 
  dplyr::select(ISO, maps_profile, vax_price, HCE, VCE, CE) %>% 
  pivot_wider(id_cols=c(maps_profile, vax_price), names_from = ISO, values_from = CE)

# if any are NA the real value is 0
# maps_profile vax_price        BFA   IND   KEN   MWI   NPL
# <chr>        <chr>          <int> <int> <int> <int> <int>
# 1 Base         $2.25 per dose     6   138    45    14     3
# 2 Base         $3.00 per dose     1    63    38     2    NA
# 3 Optimistic   $2.25 per dose    13   488    47    29    10
# 4 Optimistic   $3.00 per dose     6   103    42    11     2

tmp_all2 = tmp2 %>% 
  group_by(ISO, maps_profile, vax_price) %>% 
  summarize(Cases_averted = sum(Cases_averted),
            Deaths_averted = sum(Deaths_averted),
            DALYs_averted = sum(DALYs_averted),
            Cost_dif = sum(Cost_dif),
            ICER = sum(Cost_dif)/sum(DALYs_averted),
            ICERscaled = sum(Cost_dif)/sum(DALYs_averted)/mean(gdp),
            pop=sum(pop), 
            N = sum(N)) %>% 
  dplyr::select(ISO, maps_profile, vax_price, ICERscaled) %>% 
  pivot_wider(id_cols=ISO, names_from = c(maps_profile, vax_price), values_from = ICERscaled) %>% t 

# ISO                       "BFA"     "IND"     "KEN"     "MWI"     "NPL"    
# Base_$2.25 per dose       " 2.7274" " 7.6514" " 0.6725" " 3.4001" "20.8392"
# Base_$3.00 per dose       " 5.056"  "15.507"  " 1.255"  " 6.351"  "42.871" 
# Optimistic_$2.25 per dose "0.6489"  "1.0944"  "0.1570"  "0.7462"  "6.0431" 
# Optimistic_$3.00 per dose " 2.9771" " 8.9500" " 0.7393" " 3.6971" "28.0753"

tmp_cfr = icer_results_all %>% dplyr::filter(horizon=="20y" & comp=="comparator1", 
                                             maps_profile=="Pessimistic" & ns_strat == "RoutineCampaign", 
                                             discounting=="disc") %>% 
  dplyr::select(ISO, `Status Quo_Cases`, `Status Quo_Deaths`) %>% 
  group_by(ISO) %>%
  summarize(`Status Quo_Cases` = sum(`Status Quo_Cases`),
            `Status Quo_Deaths` = sum(`Status Quo_Deaths`)) %>% 
  mutate(CFR_pract = `Status Quo_Deaths`/`Status Quo_Cases`*100)
# 
# ISO   `Status Quo_Cases` `Status Quo_Deaths` CFR_pract
# 1 BFA              212058.               1839.     0.867
# 2 IND             6558521.               8974.     0.137
# 3 KEN              523131.              12187.     2.33 
# 4 MWI              183371.               4794.     2.61 
# 5 NPL              705578.                942.     0.134
