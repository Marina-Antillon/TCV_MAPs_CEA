#################### Calculating CEA metrics (dalys, costs) and performing calculations of NMB ####################

## See code03_objects.xlsx for notes on objects created here

##############################
# Cost of intervention ----
# SANITY CHECK NECESSARY
# Easy one: no intervention cost $0 
# must fix 2-dose scenario in campaigns. 2 doses only apply to children under 5
##############################

tmp = dimnames(doses)
tmp$vax_type = c("NS", "MAPS")
tmp$maps_profile = c("Base", "Pessimistic", "Optimistic")
doses_ns_maps = array(0, dim=sapply(1:length(tmp), function(i){length(tmp[[i]])}), 
                   dimnames=tmp)

doses_ns_maps[,,"comparator1",,"Status Quo",,,"NS",] = doses[,,"comparator1",,"Status Quo",,]
doses_ns_maps[,,"comparator1",,"MAPS add",,,"NS",] = doses[,,"comparator1",,"Status Quo",,]*0.2
doses_ns_maps[,,"comparator1",,"MAPS add",,,"MAPS",] = doses[,,"comparator1",,"MAPS add",,] - doses[,,"comparator1",,"Status Quo",,]*0.2

doses_ns_maps[,,"comparator2",,"Status Quo",,,"NS",] = doses[,,"comparator2",,"Status Quo",,]
doses_ns_maps[,,"comparator2",,"MAPS add",,as.character(c(1,3:6)),"MAPS",] = doses[,,"comparator2",,"MAPS add",,as.character(c(1,3:6))]
doses_ns_maps[,,"comparator2",,"MAPS add",,"2","NS",] = doses[,,"comparator2",,"MAPS add",,"2"]

doses_ns_maps[,,"comparator3",,"Status Quo",,,"NS",] = doses[,,"comparator3",,"Status Quo",,]
doses_ns_maps[,,"comparator3",,"MAPS add",,,"NS",] = 0
doses_ns_maps[,,"comparator3",,"MAPS add",,,"MAPS",] = doses[,,"comparator3",,"MAPS add",,]

# different profiles of the vaccine...
tmp_doses_wastage = sweep(doses_ns_maps, 8, 1+c(fix$vax_waste_ns,fix$vax_waste_maps), "*")

tmp_del_cost = abind(list(cbind(unc$vax_del_ns_uc, unc$vax_del_maps_uc_base), 
                    cbind(unc$vax_del_ns_uc, unc$vax_del_maps_uc_pess),
                    cbind(unc$vax_del_ns_uc, unc$vax_del_maps_uc_opti)), rev.along=0)
tmp_vax_del = sweep(tmp_doses_wastage, c(7,8,9), tmp_del_cost, "*")

tmp = dimnames(tmp_vax_del)
tmp_vax_sup = array(0, dim=sapply(1:length(tmp), function(i){length(tmp[[i]])}), 
                      dimnames=tmp) # remember: no supplies for the MAPS :)
# only need supplies when NS is actually administered, not when is it wasted
tmp_vax_sup[,,,,,,,"NS",] = doses_ns_maps[,,,,,,,"NS",]*unc$vax_sup

tmp_vax_price = sweep(tmp_doses_wastage, 8, c(fix$ns_price, fix$maps_price), "*")

tmp = dimnames(doses_ns_maps)
tmp$discounting = lbl$discounting
int_costs = array(0, dim=sapply(1:length(tmp), function(i){length(tmp[[i]])}), dimnames=tmp)

int_costs[,,,,,,,,,"no_disc"]  = tmp_vax_del + tmp_vax_sup + tmp_vax_price
int_costs[,,,,,,,,,"disc"] = sweep(int_costs[,,,,,,,,,"no_disc"], 1, 1/(1+fix$disc)^(1:length(years)),"*") 

# make int_costs_horizon
tmp = apply(apply(int_costs,c(1,3:10),sum), 2:9, cumsum)
tmp_horizons = tmp[c(10,20),,,,,,,,]
int_costs_horizon = tmp_horizons
dimnames(int_costs_horizon) = list(horizon=lbl$horizon,
                                   comp=paste("comparator", 1:3),
                                   ns_strat=c("None", "Routine", "RoutineCampaign"),
                                   strategy=lbl$strat, 
                                   wealth_quintile=lbl$wq, 
                                   use_case = as.character(1:6),
                                   vax_type = c("NS", "MAPS"),
                                   maps_profile = c("Base", "Pessimistic", "Optimistic"),
                                   discounting=lbl$discounting)

int_costs_dif_horizon = tmp_horizons[,,,2,,,,,] - tmp_horizons[,,,1,,,,,]
dimnames(int_costs_dif_horizon) = list(horizon=lbl$horizon,
                                       comp=paste("comparator", 1:3),
                                       ns_strat=c("None", "Routine", "RoutineCampaign"),
                                       wealth_quintile=lbl$wq, 
                                       use_case = as.character(1:6),
                                       vax_type = c("NS", "MAPS"),
                                       maps_profile = c("Base", "Pessimistic", "Optimistic"),
                                       discounting=lbl$discounting)

# clean up time
rm(list=ls(pattern="tmp")[!(ls(pattern="tmp") %in% keep)])

###############################
## Quick sanity checks ----
###############################

int_df = int_costs[,,"comparator1","Routine",,,,,,"no_disc"] %>% apply(c(1,3:5,7), sum) %>% 
  as.tbl_cube(met_name="costs") %>% as_tibble 

ggplot(int_df %>% filter(maps_profile=="Base"), aes(x=year, y=costs, color=strat)) + 
  geom_line() + 
  facet_grid(use_case~wealth_quintile) + 
  theme(legend.position = "bottom")

if(exists("zs")){ # if this is the DD analysis, then zs exists in the environment
  ggsave(paste0(tmp_subdir_fig, "/int_costs_sanity_byWQ_BaseProfile.pdf"), 
         width = 7.5, height = 8, units="in", dpi=600)
} else {
  ggsave(paste0(outfolder, ISO[CN[z],2], "/figures/int_costs_sanity_byWQ_BaseProfile.pdf"), 
         width = 7.5, height = 8, units="in", dpi=600)
}
