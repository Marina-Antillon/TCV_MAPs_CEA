fit_mod_wq = function(data, parameters, R0init=3, R0minit=1.2, repinit){
  # epsilonfit=F, watertf=T, 
  # INPUTS: 
  # data: the site data
  # epsilonfit: whether you fit epsilon. REMOVED
  # watertf: whether the water mechanism is turned on. REMOVED
  # the rest: parameters that are usually fixed. 
  
  # pars we DO fit: lnR0, R0 multiples for each WQ, logm1 (narrow sd), logm2 (narrow sd), logrep, logrC (narrow sd)
  # these will be the initial values given to the fitting command.
  logR0 = qlogis((log(R0init)+1)/2.609) # essentially, the limits are unif(0.37, 5) on the R0 scale (exponentiated), but this means ± ses will not make these go off the ranges
  logR0m = qlogis(log(R0minit)/2.303) 
  # essentially, the limits are unif(1, 10) on the R0 scale (exponentiated), but this means ± ses will not make these go off the ranges
  # in other words, R0 could be between 1-10 times as high if NO one has water and sanitation. 
  
  logm1 = data$m1_loc # qlogis(0.37)
  logm2 = data$m2_loc # qlogis(0.70)
  logrC = data$rC_loc # qlogis(0.25)
  
  #***************************************************************
  # parameters ----------
  #***************************************************************

  m1 = plogis(logm1)
  m2 = plogis(logm2)
  parameters$rC = parameters$r = plogis(logrC)
  # parameters$epsilon = ifelse(epsilontf==T, -log(0.5)/(exp(loginvepsilon)*52), epsilon)
  # report = logit_inv(logrep)
  # R0_WQ = plogis(pars[1])*10 + (100-data$haves)/100*plogis(pars[2])*10
  R0_WQ = exp(plogis(logR0)*2.609-1 + (100-data$haves)/100*(plogis(logR0m)*2.303))
  # c(1, R0m2, R0m3, R0m4, R0m5)
  
  coeff = c(m1, m1, m2, rep(1, parameters$al-3)) 
  
  # derive the effective contact rate from R0. 
  tmp_theta = sum(rep(parameters$py/parameters$num_wq, parameters$num_wq)/sum(parameters$py)*parameters$theta) # should be around 0.01542 for IND.
  parameters$betap_wq = c(t(R0_WQ %o% coeff*(parameters$mub+parameters$delta)/(1+parameters$delta*parameters$rC*tmp_theta/parameters$mub)))
  # parameters$rC = rC
  
  #***************************************************************
  # state of the world in disease-free equillibrium (the start of the simulation) ----------
  #***************************************************************
  init_ode = c(S1=rep(parameters$py/parameters$num_wq, parameters$num_wq), # 1e5*agedist
               I1=rep(0, parameters$al*parameters$num_wq), 
               R=rep(0, parameters$al*parameters$num_wq), 
               C=rep(0, parameters$al*parameters$num_wq), 
               S2=rep(0, parameters$al*parameters$num_wq), 
               I2=rep(0, parameters$al*parameters$num_wq), 
               cumI1=rep(0, parameters$al*parameters$num_wq),
               lambda=rep(0, parameters$al*parameters$num_wq))
  
  init_ode[['I13']] = 100 # starting with 100 primary infections in the third age group
  init_ode[['S13']] = init_ode[['S13']] - init_ode[['I13']] #take out the 100 people from the susceptibles
  
  for (i in 1:(parameters$al*parameters$num_wq)){
    init_ode[[paste0('lambda', i)]] = init_ode[['I13']]*parameters$betap_wq[i]/parameters$population
  }
  
  burn = 200*52 # Is a 100-year burn-in enough? 
  
  # pars = c(logR0, logR0m2, logR0m3, logR0m4, logR0m5, logm1, logm2, logrC) # 
  pars = c(logR0, logR0m, logm1, logm2, logrC, qlogis(mean(data$reporting))) #
  
  parameters$burntime = 200*52 # period after which there is routine vaccination
  times = seq(0, parameters$burntime+20*52, by = 1)
  
  #***************************************************************
  # the likelihood function ----
  #***************************************************************
  # to make optim work, we need to make wrap the kol_ll 
  # function in a function that only requires one parameter 
  # element (called p here)
  # (It seems redundant, I know, but it just works.)
  # ghn_ll comes from fcn_site_ll
  modelcost = function(p){
    nll_out = inc_ll_wq(p, init_ode, data, parameters)
    # epsilontf indicates whether we fit epsilon or not.
    # watertf indicates whether we are simulating water or not.
    return(nll_out)
  }
  
  # the code below just making sure ghn_ll works (don't need to run 
  # this all the time, only if you change something)
  # beg=Sys.time()
  # out2ll = inc_ll(p, init_ode, data, parameters)
  # end=Sys.time()
  
  #***************************************************************
  # fitting ----
  #***************************************************************
  # runs without Hessian but not with. Maybe boundary problem?
  # optimfit = optim(pars, modelcost, method="Nelder-Mead", control=c(maxit = 1000), hessian=F)
  
  t1 =Sys.time()
  optimfit = optim(pars, modelcost, method="Nelder-Mead", control=c(maxit = 3000), hessian=T)
  t2 = Sys.time()
  # t2-t1 tells you how long it takes
  
  # modelcost(c(log(3), log(1.5)))
  # modelcost(c(log(3.5), log(1.5)))
  
  print(c(Fit_time_min=as.integer(difftime(t2, t1, units="min")), 
          Convergence=as.character(optimfit$convergence==0),
          Positive_definite=as.character(is.positive.definite(optimfit$hessian)))) 
  # optimfit$convergence==0 means it reached convergence.
  # the positive definite command is from command from matrixcalc
  
  #***************************************************************
  # make summaries of the parameter values ----
  #***************************************************************
  
  results = c()

  cov = ginv(optimfit$hessian)
  
  ses = sqrt(c(cov[1,1], cov[2,2], cov[3,3], cov[4,4], cov[5,5], cov[6,6]))

  results[1] = paste(format(round(exp(plogis(optimfit$par[1])*2.609-1), 2), nsmall=2), "(", 
                     format(round(exp(plogis(optimfit$par[1] - 1.96*ses[1])*2.609-1), 2), nsmall=2), ", ",
                     format(round(exp(plogis(optimfit$par[1] + 1.96*ses[1])*2.609-1), 2), nsmall=2), ")", sep="")
  results[2] = paste(format(round(exp(1/10*(plogis(optimfit$par[2])*2.303)), 2), nsmall=2), "(", 
                     format(round(exp(1/10*(plogis(optimfit$par[2] - 1.96*ses[2])*2.303)), 2), nsmall=2), ", ",
                     format(round(exp(1/10*(plogis(optimfit$par[2] + 1.96*ses[2])*2.303)), 2), nsmall=2), ")", sep="")
  results[3] = paste(format(round(plogis(optimfit$par[3]), 2), nsmall=2), "(", 
                     format(round(plogis(optimfit$par[3] - 1.96*ses[3]), 2), nsmall=2), ", ",
                     format(round(plogis(optimfit$par[3] + 1.96*ses[3]), 2), nsmall=2), ")", sep="")
  results[4] = paste(format(round(plogis(optimfit$par[4]), 2), nsmall=2), "(", 
                      format(round(plogis(optimfit$par[4] - 1.96*ses[4]), 2), nsmall=2), ", ",
                      format(round(plogis(optimfit$par[4] + 1.96*ses[4]), 2), nsmall=2), ")", sep="")
  results[5] = paste(format(round(plogis(optimfit$par[5]), 2), nsmall=2), "(",
                     format(round(plogis(optimfit$par[5] - 1.96*ses[5]), 2), nsmall=2), ", ",
                     format(round(plogis(optimfit$par[5] + 1.96*ses[5]), 2), nsmall=2), ")", sep="")
  results[6] = paste(format(round(plogis(optimfit$par[6]), 2), nsmall=2), "(",
                     format(round(plogis(optimfit$par[6] - 1.96*ses[6]), 2), nsmall=2), ", ",
                     format(round(plogis(optimfit$par[6] + 1.96*ses[6]), 2), nsmall=2), ")", sep="")
  
  results[7] = ifelse(is.positive.definite(optimfit$hessian), "positive-definite",
                  ifelse(is.negative.definite(optimfit$hessian), "negative-definite","indefinite"))
  results[8] = ifelse(optimfit$convergence==0, "converged", "not converged")

  names(results) = c("R0", "R0m", "m1", "m2", "rC", "rep", "cov_matrix", "convergence") # "m2", 
  
  tmp = mvrnorm(10000, optimfit$par[1:2], cov[1:2, 1:2])
  tmpR0 = exp((rep(1, parameters$num_wq) %o% (plogis(tmp[,1])*2.609-1)) + (((100-data$haves)/100) %o% (plogis(tmp[,2])*2.303)))
  R0_results = tmpR0 %>% apply(1, meanpi) %>% apply(2, ci_string_dec, 2)
  
  #***************************************************************
  # Now simulating SIR # to see the predicted values ----
  #***************************************************************
  
    m1 = plogis(optimfit$par[3])
    m2 = plogis(optimfit$par[4])

    parameters$r=parameters$rC=plogis(optimfit$par[5])
    
    # R0_WQ = plogis(optimfit$par[1])*10 + (100-data$haves)/100*plogis(optimfit$pars[2])*10
    R0_WQ = exp(plogis(optimfit$par[1])*2.609-1 + (100-data$haves)/100*(plogis(optimfit$par[2])*2.303))
    coeff = c(m1, m1, m2, rep(1, parameters$al-3)) 
    
    # derive the effective contact rate from R0.
    tmp_theta = sum(parameters$py/sum(parameters$py)*parameters$theta[1:6])
    parameters$betap_wq = c(t(R0_WQ %o% coeff*(parameters$mub+parameters$delta)/(1+parameters$delta*parameters$rC*tmp_theta/parameters$mub)))

    init_ode = c(S1=rep(parameters$py/parameters$num_wq, parameters$num_wq), # 1e5*agedist
                 I1=rep(0, parameters$al*parameters$num_wq), 
                 R=rep(0, parameters$al*parameters$num_wq), 
                 C=rep(0, parameters$al*parameters$num_wq), 
                 S2=rep(0, parameters$al*parameters$num_wq), 
                 I2=rep(0, parameters$al*parameters$num_wq), 
                 cumI1=rep(0, parameters$al*parameters$num_wq),
                 lambda=rep(0, parameters$al*parameters$num_wq))
    
    init_ode[['I13']] = 100 # starting with 100 primary infections in the third age group
    init_ode[['S13']] = init_ode[['S13']] - init_ode[['I13']] #take out the 100 people from the susceptibles
    
    for (i in 1:(parameters$al*parameters$num_wq)){
      init_ode[[paste0('lambda', i)]] = init_ode[['I13']]*parameters$betap_wq[i]/parameters$population
    }
    
    times = seq(0, burn, by = 1)
    # bef = Sys.time()
    out_tmp = data.frame(ode(y = init_ode, times = times, func = typhoidsir_wq, parms = parameters, method="euler"))
    # aft= Sys.time()

  # What was the incidence in the last year of the simulation?
  simcases = out_tmp[burn, paste('cumI1', 1:(parameters$al*parameters$num_wq), sep="")] - 
    out_tmp[burn-52, paste('cumI1', 1:(parameters$al*parameters$num_wq), sep="")]
  simcases = abind(lapply(1:parameters$num_wq, function(i){simcases[,((i-1)*parameters$al+1):(i*parameters$al)]}), along=1)
  simcases = simcases*plogis(optimfit$par[6])
  # fit_results = list(optimfit, modelcost, ghn_ll, logit_inv, logit, data, parameters, pars, typhoidsir_water, simcases)
  
  fit_results = list(optimfit=optimfit, modelcost=modelcost, inc_ll_wq=inc_ll_wq, results=results,
                     data=data, parameters=parameters, pars=pars, simcases=simcases,
                     R0_results=R0_results)
  
  return(fit_results)
}
