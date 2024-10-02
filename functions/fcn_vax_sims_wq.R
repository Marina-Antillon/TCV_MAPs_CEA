vax_sims_nowater_wq = function(parameters, vprob = 0.90, omegav = 1/(10*52), rcov = rep(0.94, 5), ccov = rep(0.83, 5)){
  
#***************************************************************
## expand init_ode ------
#***************************************************************

init_ode_vax_wq = c(S1 = parameters$newpop_wq, # rep(parameters$py/parameters$num_wq, parameters$num_wq)
                 I1 = rep(0, parameters$al*parameters$num_wq), 
                 R = rep(0, parameters$al*parameters$num_wq), 
                 C = rep(0, parameters$al*parameters$num_wq), 
                 S2 = rep(0, parameters$al*parameters$num_wq), 
                 I2 = rep(0, parameters$al*parameters$num_wq), 
                 cumI1 = rep(0, parameters$al*parameters$num_wq),
                 lambda = rep(0, parameters$al*parameters$num_wq),
                 Sv = rep(0, parameters$al*parameters$num_wq),
                 Iv = rep(0, parameters$al*parameters$num_wq),
                 V1 = rep(0, parameters$al*parameters$num_wq),
                 V2 = rep(0, parameters$al*parameters$num_wq),
                 cumIv = rep(0, parameters$al*parameters$num_wq),
                 vaxr = rep(0, parameters$al*parameters$num_wq),
                 vaxc = rep(0, parameters$al*parameters$num_wq)) # now = 0, 

init_ode_vax_wq[['I13']] = 10 # put them in a specific quintile? spread them around?
init_ode_vax_wq[['S13']] = init_ode_vax_wq[['S13']] - init_ode_vax_wq[['I13']]

for (i in 1:(parameters$al*parameters$num_wq)){
  init_ode_vax_wq[[paste0('lambda', i)]] = init_ode_vax_wq[['I13']]*parameters$betap_wq[i]/(parameters$population/parameters$num_wq)
}

parameters$burntime = 200*52 # period after which there is routine vaccination
times = seq(0, parameters$burntime+30*52, by = 1)

#***************************************************************
# Simulation of vaccine strategies ----------
#***************************************************************

# the ODE might be stiff if we do the campaign in one week, so 
# we simulate it in 4 weeks.
parameters$camp1 = parameters$burntime-5 # time when the campaign starts
parameters$camp2 = parameters$burntime-1 # time when the campaign ends

parameters$vprob = vprob # probability that vaccination immunizes recipient. sensitivity for this?
parameters$omegav = omegav # rate of waning. 

parameters$uinr = rep(0, parameters$al*parameters$num_wq)
parameters$uinc = rep(0, parameters$al*parameters$num_wq)

parameters$rcov = rcov
parameters$ccov = 1-(1-ccov)^0.25 # min(1-(1-ccov)^0.25, 0.62) # quad root because we perform the campaign over 4 weeks
parameters$ccov[parameters$ccov>0.62] = 0.62

#z3=Sys.time()
outvax0 = data.frame(ode(y = init_ode_vax_wq, times = times, func = typhoidsir_vax_wq, parms = parameters, method = "euler"))
#z4=Sys.time()

# plot(rowSums(outvax0[(length(times)-52*20):length(times), paste('cumI1', 1:10, sep="")] -
#                outvax0[(length(times)-52*21):(length(times)-52), paste('cumI1', 1:10, sep="")])*(1-exp(-optimfit$par[4]))*.75*.6, ylab="")

# Vaccine simulation with no campaign  
parameters$uinr = rep(c(0,1,rep(0,parameters$al-2)), parameters$num_wq) # vector: age groups that go through routine vaccination
parameters$uinc = rep(0, parameters$al*parameters$num_wq)

outvax1 = data.frame(ode(y = init_ode_vax_wq, times = times, func = typhoidsir_vax_wq, parms = parameters, method="euler"))

# plot(rowSums(outvax1[(length(times)-52*20):length(times), paste('cumI1', 1:10, sep="")] -
#                outvax1[(length(times)-52*21):(length(times)-52), paste('cumI1', 1:10, sep="")] +
#                outvax1[(length(times)-52*20):length(times), paste('cumIv', 1:10, sep="")] -
#                outvax1[(length(times)-52*21):(length(times)-52), paste('cumIv', 1:10, sep="")])*logit_inv(optimfit$par[4])*.75*.6, ylab="")

# Vaccine simulation with campaign  
parameters$uinc = rep(c(0,1,1,1,0,0), parameters$num_wq) # vector: age groups that go through campaign vaccinations
outvax2 = data.frame(ode(y = init_ode_vax_wq, times = times, func = typhoidsir_vax_wq, parms = parameters, method="euler"))

#***************************************************************
# Make the outputs more useful ----------
#***************************************************************

hor = seq(52, 31*52, 52)
times = dim(outvax0)[1] 

simcases_novac = outvax0[c(rev(times-hor)+52), paste0('cumI1', 1:(parameters$al*parameters$num_wq))] - outvax0[rev(times-hor), paste0('cumI1', 1:(parameters$al*parameters$num_wq))]
simcases_novac = abind(lapply(1:parameters$num_wq, function(i){simcases_novac[,((i-1)*parameters$al+1):(i*parameters$al)]}), along=3)
dimnames(simcases_novac) = list(years=0:30, ages=1:parameters$al, wq=1:parameters$num_wq)

simcases_unv_int1 = outvax1[c(rev(times-hor)+52), paste0('cumI1', 1:(parameters$al*parameters$num_wq))] - outvax1[rev(times-hor), paste0('cumI1', 1:(parameters$al*parameters$num_wq))]
simcases_unv_int1 = abind(lapply(1:parameters$num_wq, function(i){simcases_unv_int1[,((i-1)*parameters$al+1):(i*parameters$al)]}), along=3)
dimnames(simcases_unv_int1) = list(years=0:30, ages=1:parameters$al, wq=1:parameters$num_wq)
simcases_vac_int1 = outvax1[c(rev(times-hor)+52), paste0('cumIv', 1:(parameters$al*parameters$num_wq))] - outvax1[rev(times-hor), paste0('cumIv', 1:(parameters$al*parameters$num_wq))]
simcases_vac_int1 = abind(lapply(1:parameters$num_wq, function(i){simcases_vac_int1[,((i-1)*parameters$al+1):(i*parameters$al)]}), along=3)
dimnames(simcases_vac_int1) = list(years=0:30, ages=1:parameters$al, wq=1:parameters$num_wq)

simcases_unv_int2 = outvax2[c(rev(times-hor)+52), paste0('cumI1', 1:(parameters$al*parameters$num_wq))] - outvax2[rev(times-hor), paste0('cumI1', 1:(parameters$al*parameters$num_wq))]
simcases_unv_int2 = abind(lapply(1:parameters$num_wq, function(i){simcases_unv_int2[,((i-1)*parameters$al+1):(i*parameters$al)]}), along=3)
dimnames(simcases_unv_int2) = list(years=0:30, ages=1:parameters$al, wq=1:parameters$num_wq)
simcases_vac_int2 = outvax2[c(rev(times-hor)+52), paste0('cumIv', 1:(parameters$al*parameters$num_wq))] - outvax2[rev(times-hor), paste0('cumIv', 1:(parameters$al*parameters$num_wq))]
simcases_vac_int2 = abind(lapply(1:parameters$num_wq, function(i){simcases_vac_int2[,((i-1)*parameters$al+1):(i*parameters$al)]}), along=3)
dimnames(simcases_vac_int2) = list(years=0:30, ages=1:parameters$al, wq=1:parameters$num_wq)

# eyeballing it
# par(mfrow=c(2,3), mar = c(1, 2, 1, 0.5))
# for(i in 1:parameters$num_wq){matplot(simcases_novac[,,i], ylim=c(0,200))}
# 
# par(mfrow=c(2,3), mar = c(1, 2, 1, 0.5))
# for(i in 1:parameters$num_wq){matplot(simcases_unv_int1[,,i], ylim=c(0,200))}
# par(mfrow=c(2,3), mar = c(1, 2, 1, 0.5))
# for(i in 1:parameters$num_wq){matplot(simcases_vac_int1[,,i], ylim=c(0,200))}
# 
# par(mfrow=c(2,3), mar = c(1, 2, 1, 0.5))
# for(i in 1:parameters$num_wq){matplot(simcases_unv_int2[,,i], ylim=c(0,200))}
# par(mfrow=c(2,3), mar = c(1, 2, 1, 0.5))
# for(i in 1:parameters$num_wq){matplot(simcases_vac_int2[,,i], ylim=c(0,200))}


rout1 = t(outvax1[c(rev(times-hor)+52), paste0('vaxr', 1:(parameters$al*parameters$num_wq))] - 
                         outvax1[rev(times-hor), paste0('vaxr', 1:(parameters$al*parameters$num_wq))])

rout2 = t(outvax2[c(rev(times-hor)+52), paste0('vaxr', 1:(parameters$al*parameters$num_wq))] - 
                         outvax2[rev(times-hor), paste0('vaxr', 1:(parameters$al*parameters$num_wq))])

camp2 = as.numeric(outvax2[times, paste0('vaxc', 1:(parameters$al*parameters$num_wq))])


outvax = list(outvax0=outvax0, outvax1=outvax1, outvax2=outvax2, parameters=parameters,
              simcases_novac=simcases_novac, simcases_unv_int1=simcases_unv_int1, simcases_vac_int1=simcases_vac_int1,
              simcases_unv_int2=simcases_unv_int2, simcases_vac_int2=simcases_vac_int2, rout1=rout1, rout2=rout2, camp2=camp2)

return(outvax)

}
