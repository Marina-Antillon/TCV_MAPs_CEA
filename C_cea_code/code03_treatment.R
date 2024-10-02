#################### Calculating CEA metrics (dalys, costs) and performing calculations of NMB ####################

## See code03_objects.xlsx for notes on objects created here

###############################
## The model is essentially "coded" in an excel ----
###############################

trt_model = read_xlsx("./cea_code/code03_treatment_tree.xlsx", sheet="prob")

###############################
## Parameter arrays:treatment ----
###############################

# make these so the excel sheet will work

# because these are ORs, convert the probability to OR (1/(1-prob)) and then multiply by OR with AMR
# then convert the new OR to a probability (OR/(1+OR))
# this process makes sure that probabilities never get above 1
tmp_reporting_or = unc$seek/(1-unc$seek)*unc$amr_rel_reporting_or
tmp_hosp_or = unc$hosp/(1-unc$hosp)*unc$amr_rel_hosp_or
tmp_reporting_amr = tmp_reporting_or/(1+tmp_reporting_or)
tmp_hosp_amr = tmp_hosp_or/(1+tmp_hosp_or)

# severity classification (AMR will only affect death)
# not doing the OR thing because these probs are small
# tmp_mort_ams      = unc$ip_mort
tmp_mort_ams_nohc_or = unc$ip_mort/(1-unc$ip_mort)*unc$ip_vs_nohc_mort_or
tmp_mort_amr_or      = unc$ip_mort/(1-unc$ip_mort)*unc$amr_rel_death_or
tmp_mort_amr_nohc_or = unc$ip_mort/(1-unc$ip_mort)*unc$ip_vs_nohc_mort_or*unc$amr_rel_death_or

tmp_mort_ams_nohc = tmp_mort_ams_nohc_or/(1+tmp_mort_ams_nohc_or)
tmp_mort_amr      = tmp_mort_amr_or/(1+tmp_mort_amr_or)
tmp_mort_amr_nohc = tmp_mort_amr_nohc_or/(1+tmp_mort_amr_nohc_or)

par_array_col_name = c("AMS_nocare_mod", "AMS_nocare_sev", "AMS_nocare_death",
                       "AMS_op_mod", 
                       "AMS_ip_sev", "AMS_ip_death", "AMS_ipip_sev", "AMS_ipip_death",
                       "AMR_nocare_mod", "AMR_nocare_sev", "AMR_nocare_death",
                       "AMR_op_mod", 
                       "AMR_ip_sev", "AMR_ip_death", "AMR_ipip_sev", "AMR_ipip_death") 

trt_sims_cond = array(NA, dim=c(branches=dim(trt_model)[1], stages=4),
                      dimnames=list(branches=par_array_col_name, stages=1:4))

for(i in 1:dim(trt_model)[1]){
  trt_sims_cond[i,1] = eval(parse(text=trt_model$vars1[i]))
  trt_sims_cond[i,2] = eval(parse(text=trt_model$vars2[i]))
  trt_sims_cond[i,3] = eval(parse(text=trt_model$vars3[i]))
  trt_sims_cond[i,4] = eval(parse(text=trt_model$vars4[i]))
}

trt_sims = apply(trt_sims_cond, 1, prod) # for Sankey, need cumprod
names(trt_sims) = par_array_col_name

#sanity check
sum(trt_sims) # should all sum to exactly 1

###############################
## DALYs ---
###############################

yld_sims_branch = sapply(1:dim(trt_model)[1], function(i){eval(parse(text=trt_model$YLDs[i]))})/365.25
names(yld_sims_branch) = par_array_col_name

# pull out the probably of deaths vector to calculate YLL.
avg_deaths = sum(trt_sims[which(trt_model$stage4=="Death")])

###############################
## treatment cost ---
###############################

trt_cost_sims_branch = sapply(1:dim(trt_model)[1], function(i){eval(parse(text=trt_model$Costs[i]))})
names(trt_cost_sims_branch) = par_array_col_name

# averages
avg_yld = sum(trt_sims*yld_sims_branch)
avg_cost = sum(trt_sims*trt_cost_sims_branch)

# clean up time
rm(list=ls(pattern="tmp")[!(ls(pattern="tmp") %in% keep)])

# need to make absolute prob so far
# funky order bc of where I put UP
#might make things more conditional according to AMR, etc, so that it mirrors the tree better and the costs and deaths later too.

# current format would be good for DALY and costs at each step, but not marginalized.

bardata =  data.frame(branch=par_array_col_name,
              prob = trt_sims,
              cost = trt_sims*trt_cost_sims_branch,
              daly = trt_sims*yld_sims_branch + trt_sims*(trt_model$stage4=="Death")*fix$yll[4]) %>%
              as_tibble %>% pivot_longer(cols=c(prob, cost, daly), names_to = "type", values_to = "val") %>% 
            mutate(branch=factor(branch, levels=rev(par_array_col_name)),
                   type = factor(type, levels=c("prob", "cost", "daly")))
              
ggplot(bardata, aes(x=branch, y=val)) + 
  geom_col() + coord_flip() +  
  facet_wrap(.~type, scale="free_x") + themebar

if(exists("zs")){ # if this is the DD analysis, then zs exists in the environment
  ggsave(paste0(tmp_subdir_fig, "/treatment_tree.pdf"), 
         width = 7.5, height = 3, units="in", dpi=600)
} else {
  ggsave(paste0(outfolder, ISO[CN[z],2], "/figures/treatment_tree.pdf"), 
         width = 7.5, height = 3, units="in", dpi=600)
}

##############################
## Epi arrays: status quo ----
##############################

## Years lost to disability
tmp = dimnames(typhoid)
tmp$discounting = lbl$discounting
yld_agegrp = array(0, dim=c(dim(typhoid), 2), dimnames=tmp)

yld_agegrp[,,,,,,"no_disc"] = typhoid*epipar["reporting"]*avg_yld
yld_agegrp[,,,,,,"disc"] = sweep(yld_agegrp[,,,,,,"no_disc"], 1, 1/(1+fix$disc)^(1:length(years)), "*") # apply discounting

## Years of life lost (deaths)
deaths_agegrp = typhoid*epipar["reporting"]*avg_deaths
# dimnames(deaths_agegrp)

# apply years lost of life according to age group (not discount)
tmp = dimnames(typhoid)
tmp$discounting = lbl$discounting # add discounting labels
yll_agegrp = array(0, dim=c(dim(typhoid), 2), dimnames=tmp)
yll_agegrp[,,,,,,"no_disc"] = sweep(deaths_agegrp, 2, fix$yll, "*")

# double discounting for deaths
tmp$discountedyears=1/((1+fix$disc)^(0:100))
tmp$discountedlifeexpectancies=cumsum(tmp$discountedyears)
tmp$disc_yll = tmp$discountedlifeexpectancies[round(fix$yll)]
yll_agegrp[,,,,,,"disc"] = sweep(sweep(deaths_agegrp, 2, tmp$disc_yll, "*"), 1, 1/(1+fix$disc)^(1:length(years)), "*")
# HOW DO I DO A SANITY CHECK FOR THIS? comparing to no_disc is a must. All no-disc must be less than disc

## Dalys
daly_agegrp = yll_agegrp+yld_agegrp
# dimnames(daly_agegrp) = dimnames(yll_agegrp) nec?

daly = apply(daly_agegrp, c(1, 3:length(dim(daly_agegrp))), sum) # Sum over the age groups.

tmp = apply(daly, 2:6, cumsum)
daly_horizon = tmp[c(10,20),,,,,]
dimnames(daly_horizon) = list(horizon=lbl$horizon,
                              comp=paste("comparator", 1:3),
                              ns_strat=c("None", "Routine", "RoutineCampaign"),
                              strategy=lbl$strat, 
                              wealth_quintile=lbl$wq, 
                              discounting=lbl$discounting)

# clean up time
rm(list=ls(pattern="tmp")[!(ls(pattern="tmp") %in% keep)])

##############################
## Epi arrays: dalys averted ----
## SANITY CHECK necessary
##############################
daly_averted_horizon = daly_horizon[,,,1,,]-daly_horizon[,,,2,,]

##############################
## Cost arrays: treatment ----
## I am discounting at year 1, that year discount should = 0
## drug costs among non-HCS cases yet?
##############################
tmp = dimnames(typhoid)
tmp$discounting = lbl$discounting
trt_costs = array(0, dim=c(dim(typhoid), 2), dimnames=tmp)

trt_costs[,,,,,,"no_disc"] = typhoid*epipar["reporting"]*avg_cost
trt_costs[,,,,,,"disc"] = sweep(trt_costs[,,,,,,"no_disc"], 1, 1/(1+fix$disc)^(1:length(years)), "*") # apply discounting

# look at the costs summed over the horizons
tmp = apply(apply(trt_costs,c(1,3:7),sum), 2:6, cumsum)
trt_costs_horizon = tmp[c(10,20),,,,,]
dimnames(trt_costs_horizon) = list(horizon=lbl$horizon,
                                   comp=paste("comparator", 1:3),
                                   ns_strat=c("None", "Routine", "RoutineCampaign"),
                                   strategy=lbl$strat, 
                                   wealth_quintile=lbl$wq, 
                                   discounting=lbl$discounting)

# treatment costs averted
trt_costs_averted_horizon = trt_costs_horizon[,,,2,,]-trt_costs_horizon[,,,1,,]

# eyeball sanity check
tmp_trt_cost = trt_costs_horizon %>% 
  cubelyr::as.tbl_cube(met_name="trt_cost") %>% as_tibble %>% 
  pivot_wider(id_cols = c("horizon", "comp", "ns_strat",  
                          "wealth_quintile", "discounting"),
              names_glue = "{strategy}_TrtCost",
              names_from = "strategy", values_from = "trt_cost")

# clean up time
rm(list=ls(pattern="tmp")[!(ls(pattern="tmp") %in% keep)])

###############################
## Quick sanity checks ----
###############################

trt_df = trt_costs[,,,,,,"no_disc"] %>% apply(c(1,3:6), sum)
trt_df = trt_df %>% as.tbl_cube(met_name="cases") %>% as_tibble

ggplot(trt_df, aes(x=year, y=cases, color=strat)) + 
  geom_line() + 
  facet_grid(ns_strat+comp~wealth_quintile) + 
  theme(legend.position = "bottom")

if(exists("zs")){ # if this is the DD analysis, then zs exists in the environment
  ggsave(paste0(tmp_subdir_fig, "/trt_costs_sanity_byWQ.pdf"), 
         width = 7, height = 10, units="in", dpi=600)
} else {
  ggsave(paste0(outfolder, ISO[CN[z],2], "/figures/trt_costs_sanity_byWQ.pdf"), 
         width = 7, height = 10, units="in", dpi=600)
}

daly_df = daly[,,,,,"no_disc"]  %>% as.tbl_cube(met_name="cases") %>% as_tibble 

ggplot(daly_df, aes(x=year, y=cases, color=strat)) + 
  geom_line() + 
  facet_grid(ns_strat+comp~wealth_quintile) + 
  theme(legend.position = "bottom")

if(exists("zs")){ # if this is the DD analysis, then zs exists in the environment
  ggsave(paste0(tmp_subdir_fig, "/dalys_sanity_byWQ.pdf"), 
         width = 7, height = 10, units="in", dpi=600)
} else {
  ggsave(paste0(outfolder, ISO[CN[z],2], "/figures/dalys_sanity_byWQ.pdf"), 
         width = 7, height = 10, units="in", dpi=600)
}

deaths_df = deaths_agegrp %>% apply(c(1,3:6), sum) %>% as.tbl_cube(met_name="cases") %>% as_tibble 
  
ggplot(deaths_df, aes(x=year, y=cases, color=strat)) + 
  geom_line() + 
  facet_grid(ns_strat+comp~wealth_quintile) + 
  theme(legend.position = "bottom")

if(exists("zs")){ # if this is the DD analysis, then zs exists in the environment
  ggsave(paste0(tmp_subdir_fig, "/deaths_sanity_byWQ.pdf"), 
         width = 7, height = 10, units="in", dpi=600)
} else {
  ggsave(paste0(outfolder, ISO[CN[z],2], "/figures/deaths_sanity_byWQ.pdf"), 
         width = 7, height = 10, units="in", dpi=600)
}


