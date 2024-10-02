inc_ll_wq = function(pars, init_ode, data, parameters){
  # INPUTS: 
  # pars: parameters you fit
  # init_ode: vector of the population in the compartments at time 0 (usually the population
  # size by age in the S1 compartments with 100 people in the infected compartment for age group 3)
  # data: the site data
  # parameters: parameters you DON'T fit.
  # epsilonfit: whether you fit epsilon. REMOVED.
  # watertf: whether the water mechanism is turned on. REMOVED.
  
  # pars: R0, logm1, logrep, logepsilon, #logm2, 
  burn = 200*52 # 100 year burn-in enough? 
  # report = logit_inv(pars[8])
  
  # Now simulating SIR
  out = with(as.list(pars, init_ode=init_ode, parameters=parameters), {
  
    m1 = plogis(pars[3])
    m2 = plogis(pars[4])
    # rC = plogis(pars[5])
    parameters$rC = parameters$r = plogis(pars[5])
    coeff = c(m1, m1, m2, rep(1, parameters$al-3))
    # R0_WQ = plogis(pars[1])*10 + (100-data$haves)/100*plogis(pars[2])*10
    R0_WQ = exp(plogis(pars[1])*2.609-1 + (100-data$haves)/100*(plogis(pars[2])*2.303))
    
    # parameters$betap = exp(pars[1])*coeff*(parameters$mub+parameters$delta)/(1+parameters$delta*parameters$rC*sum(parameters$py/1e5*parameters$theta)/parameters$mub)
    tmp_theta = sum(rep(parameters$py/parameters$num_wq, parameters$num_wq)/sum(parameters$py)*parameters$theta) # should be around 0.01542 for IND.
    parameters$betap_wq = c(t(R0_WQ %o% coeff*(parameters$mub+parameters$delta)/(1+parameters$delta*parameters$rC*tmp_theta/parameters$mub)))
      
    times = 0:burn
    # bef = Sys.time()
    out_tmp = data.frame(ode(y = init_ode, times = times, func = typhoidsir_wq, parms = parameters, method="euler"))
    # aft= Sys.time()
    # data.frame(out)
  })
  
  # What was the incidence in the last year of the simulation?
  simcases = out[burn, paste('cumI1', 1:(parameters$al*parameters$num_wq), sep="")] - out[burn-52, paste('cumI1', 1:(parameters$al*parameters$num_wq), sep="")]
  simcases = abind(lapply(1:parameters$num_wq, function(i){simcases[,((i-1)*parameters$al+1):(i*parameters$al)]}), along=1)
  simcases = simcases * plogis(pars[6]) # taking into account reporting...
  disp_sim = apply(simcases, 1, sum)
  simavgage = sum(colSums(simcases) * data$mdpts_age/sum(simcases))
    
  # Now calculating the negative log-likelihood.
    nll = -(dnorm(sum(simcases), data$mean_inc, data$sd_inc, log=T) +
              ddirichlet(disp_sim/sum(disp_sim), data$dirichletalphas, log=T) +  
              dnorm(simavgage, data$mean_age, data$sd_age, log=T) + #avg age, to limit the sd with R0
              dlogis(pars[1], 0, 1, log=T) + # R0 if everyone had a toilet/improved water
              dlogis(pars[2], 0, 1, log=T) + # R0m multiplier times prevalence of have-nots
              dmvnorm(pars[3:4], c(data$m1_loc, data$m2_loc), data$m_Sigma, log=T) + 
              # dlogis(pars[3], location=data$m1_loc, scale=data$m1_sca, log=T) + # for transformed m1, data$m1_sca
              # dlogis(pars[4], location=data$m2_loc, scale=data$m2_sca, log=T) + # for transformed m2, data$m2_sca
              dlogis(pars[5], location=data$rC_loc, scale=data$rC_sca, log=T) + 
              dlogis(pars[6], location=data$rep_loc, scale=data$rep_sca, log=T)) #+

  return(nll) # min(max(nll, 0), 2000) # because of dirichlet, the nll is going under 0 and that's fine
}
