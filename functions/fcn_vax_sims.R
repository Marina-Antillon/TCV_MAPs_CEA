vax_sims_nowater = function(parameters, vprob = 0.90, omegav = 1/(10*52), rcov = 0.94, ccov = 0.83){
  
#***************************************************************
## expand init_ode ------
#***************************************************************

init_ode_vax = c(S1 = parameters$py, 
                 I1 = rep(0, parameters$al), 
                 R = rep(0, parameters$al), 
                 C = rep(0, parameters$al), 
                 S2 = rep(0, parameters$al), 
                 I2 = rep(0, parameters$al), 
                 cumI1 = rep(0, parameters$al),
                 lambda = rep(0, parameters$al),
                 Sv = rep(0, parameters$al),
                 Iv = rep(0, parameters$al),
                 V1 = rep(0, parameters$al),
                 V2 = rep(0, parameters$al),
                 cumIv = rep(0, parameters$al),
                 vaxr = rep(0, parameters$al),
                 vaxc = rep(0, parameters$al)) # now = 0, 

init_ode_vax[['I13']] = 100
init_ode_vax[['S13']] = init_ode_vax[['S13']] - init_ode_vax[['I13']]

for (i in 1:parameters$al){
  init_ode_vax[[paste0('lambda', i)]] = init_ode_vax[['I13']]*parameters$betap[i]/parameters$population
}

# init_ode = c(S1 = parameters$population, 
#              I1 = rep(0, parameters$al), 
#              R = rep(0, parameters$al), 
#              C = rep(0, parameters$al), 
#              S2 = rep(0, parameters$al), 
#              I2 = rep(0, parameters$al),  
#              cumI1=rep(0, parameters$al),
#              lambda=rep(0, parameters$al))
# 
# init_ode[['I13']] = 100
# init_ode[['S13']] = init_ode[['S13']] - init_ode[['I13']]
# 
# for (i in 1:parameters$al){
#   init_ode[[paste0('lambda', i)]] = init_ode[['I13']]*parameters$betap[i]/sum(parameters$population)
# }

parameters$burntime = 200*52 # period after which there is routine vaccination
times = seq(0, parameters$burntime+20*52, by = 1)

#***************************************************************
# Simulation of vaccine strategies ----------
#***************************************************************

# the ODE might be stiff if we do the campaign in one week, so 
# we simulate it in 4 weeks.
parameters$camp1 = parameters$burntime-5 # time when the campaign starts
parameters$camp2 = parameters$burntime-1 # time when the campaign ends

parameters$vprob = vprob # probability that vaccination immunizes recipient. sensitivity for this?
parameters$omegav = omegav # rate of waning. 

parameters$uinr = rep(0, parameters$al)
parameters$uinc = rep(0, parameters$al)

parameters$rcov = rcov
parameters$ccov = min(1-(1-ccov)^0.25, 0.62) # quad root because we perform the campaign over 4 weeks

#z3=Sys.time()
outvax0 = data.frame(ode(y = init_ode_vax, times = times, func = typhoidsir_vax, parms = parameters, method="euler"))
#z4=Sys.time()

# plot(rowSums(outvax0[(length(times)-52*20):length(times), paste('cumI1', 1:10, sep="")] -
#                outvax0[(length(times)-52*21):(length(times)-52), paste('cumI1', 1:10, sep="")])*(1-exp(-optimfit$par[4]))*.75*.6, ylab="")

# Vaccine simulation with no campaign  
parameters$uinr = c(0,1,rep(0,parameters$al-2)) # vector: age groups that go through routine vaccination
parameters$uinc = rep(0, parameters$al)

outvax1 = data.frame(ode(y = init_ode_vax, times = times, func = typhoidsir_vax, parms = parameters, method="euler"))

# plot(rowSums(outvax1[(length(times)-52*20):length(times), paste('cumI1', 1:10, sep="")] -
#                outvax1[(length(times)-52*21):(length(times)-52), paste('cumI1', 1:10, sep="")] +
#                outvax1[(length(times)-52*20):length(times), paste('cumIv', 1:10, sep="")] -
#                outvax1[(length(times)-52*21):(length(times)-52), paste('cumIv', 1:10, sep="")])*logit_inv(optimfit$par[4])*.75*.6, ylab="")

# Vaccine simulation with campaign  
parameters$uinc = c(0,1,1,1,0,0) # vector: age groups that go through campaign vaccinations
outvax2 = data.frame(ode(y = init_ode_vax, times = times, func = typhoidsir_vax, parms = parameters, method="euler"))

#***************************************************************
# Make the outputs more useful ----------
#***************************************************************

hor = seq(52, 21*52, 52)
times = dim(outvax0)[1] 

simcases_novac = outvax0[c(rev(times-hor)+52), paste0('cumI1', 1:parameters$al)] - outvax0[rev(times-hor), paste0('cumI1', 1:parameters$al)]
simcases_unv_int1 = outvax1[c(rev(times-hor)+52), paste0('cumI1', 1:parameters$al)] - outvax1[rev(times-hor), paste0('cumI1', 1:parameters$al)]
simcases_vac_int1 = outvax1[c(rev(times-hor)+52), paste0('cumIv', 1:parameters$al)] - outvax1[rev(times-hor), paste0('cumIv', 1:parameters$al)]
simcases_unv_int2 = outvax2[c(rev(times-hor)+52), paste0('cumI1', 1:parameters$al)] - outvax2[rev(times-hor), paste0('cumI1', 1:parameters$al)]
simcases_vac_int2 = outvax2[c(rev(times-hor)+52), paste0('cumIv', 1:parameters$al)] - outvax2[rev(times-hor), paste0('cumIv', 1:parameters$al)]

# sanity checks
# matplot(simcases_novac/(rep(1, 21) %o% parameters$py)*1e5)
# matplot(simcases_unv_int1/(rep(1, 21) %o% parameters$py)*1e5)
# matplot(simcases_vac_int1/(rep(1, 21) %o% parameters$py)*1e5)
# matplot(simcases_unv_int2/(rep(1, 21) %o% parameters$py)*1e5)
# matplot(simcases_vac_int2/(rep(1, 21) %o% parameters$py)*1e5)

rout1 = t(outvax1[c(rev(times-hor)+52), paste0('vaxr', 1:parameters$al)] - 
                         outvax1[rev(times-hor), paste0('vaxr', 1:parameters$al)])
rout2 = t(outvax2[c(rev(times-hor)+52), paste0('vaxr', 1:parameters$al)] - 
                         outvax2[rev(times-hor), paste0('vaxr', 1:parameters$al)])
camp2 = as.numeric(outvax2[times, paste0('vaxc', 1:parameters$al)])

outvax = list(outvax0=outvax0, outvax1=outvax1, outvax2=outvax2, parameters=parameters,
              simcases_novac=simcases_novac, simcases_unv_int1=simcases_unv_int1, simcases_vac_int1=simcases_vac_int1,
              simcases_unv_int2=simcases_unv_int2, simcases_vac_int2=simcases_vac_int2, rout1=rout1, rout2=rout2, camp2=camp2)

return(outvax)

}
