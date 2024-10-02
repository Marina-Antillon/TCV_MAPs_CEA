#################### Dynamic model results ####################

# How params vary (dynamic model only):
# By site: "beta"   "rep"  "m1"   "m2"   "rC" (although rC has the same distribution in all areas)
# Identical across site:  "veff" "omegav"

# Elements we get from sims: 
  # pop (population at the end of the 0 year marks, 10 year mark, 20 year mark, 30 year mark)
  # cumI1 (cum infections among unvaccinated people).
  # cumI1v (cum infections among vaccinated people).
  # cumdosesr (doses for routine - administered over time - no discounting)
  # cumdosesc (doses for campaigns - administered at the beginning - plenty of discounting)
# I pulled these two measures in case we do something with them later: chr, chrv (chronic carriers)

# Initialize arrays
  
# names = list(c("years" ,"ages", "strategies", "wealth quintiles", "use cases"))
# ADD: comparator, maps_profile, price (low, high, tiered), scenario_ns, vax_duration, 
# new pop coverage (0.2, 0.1, 0.3)

# others that could be added for SA: AMR burden multiplier, hosp,

lbl$age = c('0-9m','9m-2y','2-4y','5-14y','14-24y', '25y+') 
lbl$strat = c("Status Quo", "MAPS add")

typhoid <- array(NA, dim=c(length(years), # time
                           length(lbl$age), # Ages
                           3, # comparator
                           3, # scenario_ns
                           2, # strategies
                           5), # wealth quantiles
                     dimnames = list(year=years_lbl, age=lbl$age, 
                                     comp = c("comparator1", "comparator2", "comparator3"),
                                     ns_strat = c("None", "Routine", "RoutineCampaign"),
                                     strat=lbl$strat, wealth_quintile=lbl$wq))

doses <- array(NA, dim=c(length(years), # time
                            length(lbl$age), # Ages
                            3, # comparator
                            3, # scenario_ns
                            2, # strategies
                            5, # wealth quintiles
                            6), # use cases
                    dimnames = list(year=years_lbl, age=lbl$age, 
                                    comp = c("comparator1", "comparator2", "comparator3"),
                                    ns_strat = c("None", "Routine", "RoutineCampaign"),
                                    strat=lbl$strat, wealth_quintile=lbl$wq, use_case=1:6))

for(scenario_switch in 1:3){
  load(paste0(inputfolder, ISO[CN[z],2] ,'/wq_maps_NOunc_global_comp', scenario_switch, '.Rdata'))

  epipar = rep(NA,8)
  names(epipar) = lbl$par
  epipar[1:5] = c(exp(pars_unc[1]), exp(pars_unc[2]), # this is wrong, though it doesn't matter much. It's more than just exponentiated. 
                  plogis(pars_unc[3]), plogis(pars_unc[4]), 
                  plogis(pars_unc[5]))
  epipar['reporting'] = plogis(pars_unc[6]) # rbeta(iterations, 15, 85) # data$reporting
  epipar['vprob'] = vprob
  epipar['omegav'] = omegav
  
  for(scenario_ns in 1:3){
    
    if(scenario_ns==1){
      tmp_vax_out = SIR_wq_vax_unc0
      tmp_maps_out = SIR_wq_maps_unc0
    }else if(scenario_ns==2){
      tmp_vax_out = SIR_wq_vax_unc1
      tmp_maps_out = SIR_wq_maps_unc1
    }else if(scenario_ns==3){
      tmp_vax_out = SIR_wq_vax_unc2
      tmp_maps_out = SIR_wq_maps_unc2
    }
    
    typhoid[,,scenario_switch,scenario_ns,"Status Quo",] = tmp_vax_out$simcases_maps
    typhoid[,,scenario_switch,scenario_ns,"MAPS add",] = tmp_maps_out$simcases_maps
  
    doses[,,scenario_switch,scenario_ns,"Status Quo",,"1"] = array(t(tmp_vax_out$vaxuc1),dim=c(length(years), 6, 5))
    doses[,,scenario_switch,scenario_ns,"Status Quo",,"2"] = array(t(tmp_vax_out$vaxuc2),dim=c(length(years), 6, 5))
    doses[,,scenario_switch,scenario_ns,"Status Quo",,"3"] = array(t(tmp_vax_out$vaxuc3),dim=c(length(years), 6, 5))
    doses[,,scenario_switch,scenario_ns,"Status Quo",,"4"] = array(t(tmp_vax_out$vaxuc4),dim=c(length(years), 6, 5))
    doses[,,scenario_switch,scenario_ns,"Status Quo",,"5"] = array(t(tmp_vax_out$vaxuc5),dim=c(length(years), 6, 5))
    doses[,,scenario_switch,scenario_ns,"Status Quo",,"6"] = array(t(tmp_vax_out$vaxuc6),dim=c(length(years), 6, 5))
  
    doses[,,scenario_switch,scenario_ns,"MAPS add",,"1"] = array(t(tmp_maps_out$vaxuc1),dim=c(length(years), 6, 5))
    doses[,,scenario_switch,scenario_ns,"MAPS add",,"2"] = array(t(tmp_maps_out$vaxuc2),dim=c(length(years), 6, 5))
    doses[,,scenario_switch,scenario_ns,"MAPS add",,"3"] = array(t(tmp_maps_out$vaxuc3),dim=c(length(years), 6, 5))
    doses[,,scenario_switch,scenario_ns,"MAPS add",,"4"] = array(t(tmp_maps_out$vaxuc4),dim=c(length(years), 6, 5))
    doses[,,scenario_switch,scenario_ns,"MAPS add",,"5"] = array(t(tmp_maps_out$vaxuc5),dim=c(length(years), 6, 5))
    doses[,,scenario_switch,scenario_ns,"MAPS add",,"6"] = array(t(tmp_maps_out$vaxuc6),dim=c(length(years), 6, 5))
  }
}
  gc()
  
  typh_df = (typhoid*epipar["reporting"]) %>% apply(c(1,3:5), sum) %>% 
            as.tbl_cube(met_name="cases") %>% as_tibble 
  
  ggplot(typh_df, aes(x=year, y=cases, color=strat)) + 
    geom_line() + 
    facet_grid(comp~ns_strat) + 
    scale_y_continuous(limits=c(0,NA)) +
    theme(legend.position = "bottom")
  
  ggsave(paste0(outfolder, ISO[CN[z],2], "/figures/typh_cases_sanity.pdf"), width = 6, height = 6, units="in", dpi=600)
  
  typh_df = (typhoid*epipar["reporting"]) %>% apply(c(1,3:6), sum) %>% 
    as.tbl_cube(met_name="cases") %>% as_tibble 
  
  ggplot(typh_df, aes(x=year, y=cases, color=strat)) + 
    geom_line() + 
    facet_grid(ns_strat+comp~wealth_quintile) + 
    scale_y_continuous(limits=c(0,NA)) +
    theme(legend.position = "bottom")
  
  ggsave(paste0(outfolder, ISO[CN[z],2], "/figures/typh_cases_sanity_byWQ.pdf"), width = 10, height = 10, units="in", dpi=600)
  
  dose_df = doses %>% apply(c(1,3:7), sum) %>% 
   as.tbl_cube(met_name="doses") %>% as_tibble 
  
  ggplot(dose_df, aes(x=year, y=doses, color=strat)) + 
    geom_line() + 
    facet_grid(ns_strat+comp~use_case+wealth_quintile) + 
    scale_y_continuous(limits=c(0,NA)) +
    theme(legend.position = "bottom")
  
  ggsave(paste0(outfolder, ISO[CN[z],2], "/figures/doses_sanity_byWQ.pdf"), width = 20, height = 10, units="in", dpi=600)
