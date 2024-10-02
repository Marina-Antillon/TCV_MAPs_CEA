vax_maps_sims_nowater_wq = function(parameters, init_ode_maps_wq, vprob = 0.90, omegav = 1/(10*52), 
                               uc1cov = uc1cov, uc2cov = uc2cov, uc3cov = uc3cov,
                               uc4cov = uc4cov, uc5cov = uc5cov, uc6cov = uc6cov){
  
#***************************************************************
## expand init_ode ------
#***************************************************************

# now it comes from the previous run, but add trackers for vaccinations for each use case...
init_ode_maps_wq = c(init_ode_maps_wq,
  vaxuc1 = rep(0, parameters$al*parameters$num_wq),
  vaxuc2 = rep(0, parameters$al*parameters$num_wq),
  vaxuc3 = rep(0, parameters$al*parameters$num_wq),
  vaxuc4 = rep(0, parameters$al*parameters$num_wq),
  vaxuc5 = rep(0, parameters$al*parameters$num_wq),
  vaxuc6 = rep(0, parameters$al*parameters$num_wq))
  
# parameters$burntime = 200*52 # period after which there is routine vaccination
times = seq(0, 20*52, by = 1)

#***************************************************************
# Simulation of vaccine strategies ----------
#***************************************************************

# the ODE might be stiff if we do the campaign in one week, so 
# we simulate it in 4 weeks.
parameters$camp1 = 1  # time when the campaign starts
parameters$camp2 = 4 # time when the campaign ends
# only UC4 will be a campaign and 100% maps, other outreach will be spread out throughout the year

parameters$vprob = vprob # probability that vaccination immunizes recipient. sensitivity for this?
parameters$omegav = omegav # rate of waning. 

parameters$rcov=0 # turned off 
parameters$ccov=0 # turned off 

parameters$uc1cov = rep(uc1cov, each=parameters$al) # these need to be 30 elements long, just like each compartment
parameters$uc2cov = rep(uc2cov, each=parameters$al) # rep(-log(pmax(pmin(1-uc2cov/(1-uc1cov+0.001), 1),0))/52, each=parameters$al) # +0.001 is a correction for when uc2cov covers everyone - it would go to -Inf otherwise. pmin is to make sure it doesn't go over 1
parameters$uc3cov = rep(uc3cov, each=parameters$al) # rep(-log(pmax(pmin(1-uc3cov/(1-uc1cov+0.001), 1),0))/52, each=parameters$al) # +0.001 is a correction for when uc2cov covers everyone - it would go to -Inf otherwise. pmin is to make sure it doesn't go over 1
parameters$uc4cov = rep(1-(1-uc4cov)^0.25, each=parameters$al) # min(1-(1-ccov)^0.25, 0.62) # quad root because we perform the campaign over 4 weeks
parameters$uc4cov[parameters$uc4cov>0.62] = 0.62
parameters$uc5cov = rep(-log(1-uc5cov)/52, each=parameters$al) # inverse of p1 = 1-exp(-rate*time)
parameters$uc6cov = rep(-log(1-uc6cov)/52, each=parameters$al)

# Vaccine simulation with no campaign  
parameters$uinuc1 = rep(c(0,1,0,0,0,0), times = parameters$num_wq) # vector: age groups that go through routine vaccination
parameters$uinuc2 = rep(c(0,1,0,0,0,0), parameters$num_wq) 
parameters$uinuc3 = rep(c(0,1,0,0,0,0), parameters$num_wq) 
parameters$uinuc4 = rep(c(0,1,1,1,0,0), parameters$num_wq) # one-time campaign
parameters$uinuc5 = rep(c(0,0,0,0,1,0), parameters$num_wq) # military
parameters$uinuc6 = rep(c(0,0,0,0,1,1), parameters$num_wq) # travellers

outvax_maps = data.frame(ode(y = unlist(init_ode_maps_wq), times = times, func = typhoidsir_maps_wq, parms = parameters, method="euler"))

#***************************************************************
# Make the outputs more useful ----------
#***************************************************************

hor = seq(52, 20*52, 52)
times = dim(outvax_maps)[1] 

simcases_maps = outvax_maps[c(rev(times-hor)+52), paste0('cumI1', 1:(parameters$al*parameters$num_wq))] - outvax_maps[rev(times-hor), paste0('cumI1', 1:(parameters$al*parameters$num_wq))] + 
              outvax_maps[c(rev(times-hor)+52), paste0('cumIv', 1:(parameters$al*parameters$num_wq))] - outvax_maps[rev(times-hor), paste0('cumIv', 1:(parameters$al*parameters$num_wq))]
simcases_maps = abind(lapply(1:parameters$num_wq, function(i){simcases_maps[,((i-1)*parameters$al+1):(i*parameters$al)]}), along=3)
dimnames(simcases_maps) = list(years=1:20, ages=1:parameters$al, wq=1:parameters$num_wq)

simcases_unv_maps = outvax_maps[c(rev(times-hor)+52), paste0('cumI1', 1:(parameters$al*parameters$num_wq))] - outvax_maps[rev(times-hor), paste0('cumI1', 1:(parameters$al*parameters$num_wq))]
simcases_unv_maps = abind(lapply(1:parameters$num_wq, function(i){simcases_unv_maps[,((i-1)*parameters$al+1):(i*parameters$al)]}), along=3)
dimnames(simcases_unv_maps) = list(years=1:20, ages=1:parameters$al, wq=1:parameters$num_wq)
simcases_vac_maps = outvax_maps[c(rev(times-hor)+52), paste0('cumIv', 1:(parameters$al*parameters$num_wq))] - outvax_maps[rev(times-hor), paste0('cumIv', 1:(parameters$al*parameters$num_wq))]
simcases_vac_maps = abind(lapply(1:parameters$num_wq, function(i){simcases_vac_maps[,((i-1)*parameters$al+1):(i*parameters$al)]}), along=3)
dimnames(simcases_vac_maps) = list(years=1:20, ages=1:parameters$al, wq=1:parameters$num_wq)

# eyeballing it
# par(mfrow=c(2,3), mar = c(1, 2, 1, 0.5))
# for(i in 1:parameters$num_wq){matplot(simcases_maps[,,i], ylim=c(0,200))}
# 
# par(mfrow=c(2,3), mar = c(1, 2, 1, 0.5))
# for(i in 1:parameters$num_wq){matplot(simcases_unv_maps[,,i], ylim=c(0,200))}
# par(mfrow=c(2,3), mar = c(1, 2, 1, 0.5))
# for(i in 1:parameters$num_wq){matplot(simcases_vac_maps[,,i], ylim=c(0,200))}

# rout1 = t(outvax1[c(rev(times-hor)+52), paste0('vaxr', 1:(parameters$al*parameters$num_wq))] - 
#                          outvax1[rev(times-hor), paste0('vaxr', 1:(parameters$al*parameters$num_wq))])
# 
# rout2 = t(outvax2[c(rev(times-hor)+52), paste0('vaxr', 1:(parameters$al*parameters$num_wq))] - 
#                          outvax2[rev(times-hor), paste0('vaxr', 1:(parameters$al*parameters$num_wq))])
# 
# camp2 = as.numeric(outvax2[times, paste0('vaxc', 1:(parameters$al*parameters$num_wq))])

vaxuc1 = t(outvax_maps[c(rev(times-hor)+52), paste0('vaxuc1', 1:(parameters$al*parameters$num_wq))] -
            outvax_maps[rev(times-hor), paste0('vaxuc1', 1:(parameters$al*parameters$num_wq))])
vaxuc2 = t(outvax_maps[c(rev(times-hor)+52), paste0('vaxuc2', 1:(parameters$al*parameters$num_wq))] -
            outvax_maps[rev(times-hor), paste0('vaxuc2', 1:(parameters$al*parameters$num_wq))])
vaxuc3 = t(outvax_maps[c(rev(times-hor)+52), paste0('vaxuc3', 1:(parameters$al*parameters$num_wq))] -
             outvax_maps[rev(times-hor), paste0('vaxuc3', 1:(parameters$al*parameters$num_wq))])
vaxuc4 = t(outvax_maps[c(rev(times-hor)+52), paste0('vaxuc4', 1:(parameters$al*parameters$num_wq))] -
             outvax_maps[rev(times-hor), paste0('vaxuc4', 1:(parameters$al*parameters$num_wq))])
vaxuc5 = t(outvax_maps[c(rev(times-hor)+52), paste0('vaxuc5', 1:(parameters$al*parameters$num_wq))] -
             outvax_maps[rev(times-hor), paste0('vaxuc5', 1:(parameters$al*parameters$num_wq))])
vaxuc6 = t(outvax_maps[c(rev(times-hor)+52), paste0('vaxuc6', 1:(parameters$al*parameters$num_wq))] -
             outvax_maps[rev(times-hor), paste0('vaxuc6', 1:(parameters$al*parameters$num_wq))])

outvax_maps = list(outvax_maps=outvax_maps, parameters=parameters,
              simcases_maps=simcases_maps, simcases_unv_maps=simcases_unv_maps, simcases_vac_maps=simcases_vac_maps,
              vaxuc1=vaxuc1, vaxuc2=vaxuc2, vaxuc3=vaxuc3, vaxuc4=vaxuc4, vaxuc5=vaxuc5, vaxuc6=vaxuc6)

return(outvax_maps)

}
