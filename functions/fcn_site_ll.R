inc_ll = function(pars, init_ode, data, parameters, bloodsens, part_prob){
  # INPUTS: 
  # pars: parameters you fit
  # init_ode: vector of the population in the compartments at time 0 (usually the population
  # size by age in the S1 compartments with 100 people in the infected compartment for age group 3)
  # data: the site data
  # parameters: parameters you DON'T fit.
  # epsilonfit: whether you fit epsilon
  # watertf: whether the water mechanism is turned on.
  
  # pars: R0, logm1, logrep, logepsilon, #logm2, 
  burn = 200*52 # 100 year burn-in enough? 
  report = logit_inv(pars[4])
  
  # Now simulating SIR
  out = with(as.list(pars, init_ode=init_ode, parameters=parameters), {
    
    R0 = exp(pars[1])
    m1 = logit_inv(pars[2])
    m2 = logit_inv(pars[3])
    coeff = c(m1, m2, rep(1, dim(data)[1]-2))
    
    parameters$betap = exp(pars[1])*coeff*(parameters$mub+parameters$delta)/(1+parameters$delta*parameters$rC*sum(agedist*parameters$theta)/parameters$mub)
    parameters$epsilon = 0 # duration to rate

    times = seq(0, burn, by = 1)
    # bef = Sys.time()
    out_tmp = data.frame(ode(y = init_ode, times = times, func = typhoidsir, parms = parameters, method="euler"))
    # aft= Sys.time()
    # data.frame(out)
  })
  
  # What was the incidence in the last year of the simulation?
  simcases = out[burn, paste('cumI1', 1:(dim(data)[1]), sep="")] - out[burn-52, paste('cumI1', 1:(dim(data)[1]), sep="")]
  # plot(apply(as.matrix(outshort[53:burn, 50:57] - outshort[1:(burn-52), 50:57]), 1, "sum"))

  # Now calculating the negative log-likelihood.
  nll = -(sum(log(dpois(data$cases, as.numeric(simcases*report*bloodsens*part_prob))), na.rm=T) + 
              dunif(pars[1], 0, 2.996, log=T) + # R0
              dlogis(pars[2], location=0, scale=1, log=T) + # for transformed m1
              dlogis(pars[3], location=0, scale=1, log=T) + # for transformed m2
              dlogis(pars[4], location=0, scale=1, log=T)) # for transformed prob of symptomatic
  
  # introduce weights: data$py/sum(data$py)
  # nll = -(sum(log(dpois(data$cases, as.numeric(simcases*report*0.6*.75)))) + dunif(pars[1], 1, 3, log=T) + dexp(pars[2], rate=1, log=T)
  #         + dexp(pars[3], rate=1, log=T) + dexp(pars[4], rate=1, log=T) + dexp(pars[5], rate=0.00001, log=T))
  
  return(min(max(nll, 0), 2000))
}
