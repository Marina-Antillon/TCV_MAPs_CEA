fit_mod = function(data, birthrate=25.3, delta=1/4, 
                   alpha=0.01, omega=1/104, gamma=1, 
                   epsilon=0, xi=1/3, r=0.1, rC=0.1){
  # epsilonfit=F, watertf=T, 
  # INPUTS: 
  # data: the site data
  # epsilonfit: whether you fit epsilon
  # watertf: whether the water mechanism is turned on. REMOVED
  # the rest: parameters that are usually fixed. REMOVED
  
  # pars we DO fit: lnR0, logm1, logm2, logrep, 
  # these will be the initial values given to the fitting command.
  lnR0 = log(4)
  logm1 = qlogis(0.15)
  logm2 = qlogis(0.4)
  logrep = qlogis(0.25)
  # loginvepsilon = 5
  
  bloodsens = data$sens_BC
  # hcs = data$health_seek
  part_prob = data$part_prob
  
  #***************************************************************
  # parameters ----------
  #***************************************************************
  
  parameters = list(delta=delta, # rate of recovery
                    alpha=alpha, # 0.01, # prob dying because of typhoid
                    omega=omega, # duration of immunity. Alternative: -log(1-0.25)/52
                    # see article on typhoid outbreaks in British troops two years in a row.
                    theta = data$prob_chr, # prob I1 becomes chronic.
                    theta2 = if(I2toCtf == F){rep(0, dim(data)[1])}else{data$prob_chr},
                    # prob I2 becomes chronic (this is the focus of one kind of sensitivity analysis)
                    gamma=gamma,
                    xi=xi,
                    r=r, # contribution of secondary carriers to transmission
                    rC = rC, # contribution of chronic carriers to transmission
                    mub = log(1+birthrate/1000)/52,
                    # mu = -log(1-data$mort_per_py)/52,
                    u = data$aging_rate,
                    al = dim(data)[1], # number of age groups is equivalent to the number of rows in "data" table
                    population=sum(data$py))
  
  parameters$u[length(parameters$u)] = 0 # last age group doesn't age.
  
  # This is the way I derived the synthetic death rate for the 
  # demographic file. It is the death rate necessary to keep age 
  # distribution and population size constant
  
  # synthetic aging - see write up for more detail
  aging_in = c(parameters$mub*sum(data$py), parameters$u[1:(length(parameters$u)-1)]*data$py[1:(length(parameters$u)-1)])
  aging_out = data$aging_rate*data$py
  
  parameters$mu = (aging_in - aging_out)/data$py
  
  R0 = exp(lnR0)
  m1 = plogis(logm1)
  m2 = plogis(logm2)
  # parameters$epsilon = ifelse(epsilontf==T, -log(0.5)/(exp(loginvepsilon)*52), epsilon)
  report = plogis(logrep)
  
  coeff = c(m1, m2, rep(1, dim(data)[1]-2)) 
  
  # derive the effective contact rate from R0. Assume 50% of R0 is 
  # person-to-person transmission and 50% is water-borne transmission.
  parameters$betap = R0*coeff*(parameters$mub+parameters$delta)/(1+parameters$delta*parameters$rC*sum(agedist*parameters$theta)/parameters$mub)

  #***************************************************************
  # state of the world in disease-free equillibrium (the start of the simulation) ----------
  #***************************************************************
  init_ode = c(S1=data$py, # 1e5*agedist
               I1=rep(0, dim(data)[1]), 
               R=rep(0, dim(data)[1]), 
               C=rep(0, dim(data)[1]), 
               S2=rep(0, dim(data)[1]), 
               I2=rep(0, dim(data)[1]), 
               cumI1=rep(0, dim(data)[1]),
               lambda=rep(0, dim(data)[1]))
  
  init_ode[['I13']] = 100 # starting with 100 primary infections in the third age group
  init_ode[['S13']] = init_ode[['S13']] - init_ode[['I13']] #take out the 100 people from the susceptibles
  
  for (i in 1:parameters$al){
    init_ode[[paste('lambda', i, sep="")]] = init_ode[['I13']]*parameters$betap[i]/parameters$population
  }
  
  burn = 200*52 # Is a 100-year burn-in enough? 
  
  pars = c(lnR0, logm1, logm2, logrep) # 

  #***************************************************************
  # the likelihood function ----
  #***************************************************************
  # to make optim work, we need to make wrap the kol_ll 
  # function in a function that only requires one parameter 
  # element (called p here)
  # (It seems redundant, I know, but it just works.)
  # ghn_ll comes from fcn_site_ll
  modelcost = function(p){
    nll_out = inc_ll(p, init_ode, data, parameters, bloodsens, part_prob)
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
  
  t1 =Sys.time()
  optimfit = optim(pars, modelcost, method="Nelder-Mead", control=c(maxit = 1000), hessian=T)
  # took out lower and upper control statements
  t2 = Sys.time()
  # t2-t1 tells you how long it takes
  
  print(c(Fit_time=as.numeric((t2-t1)), Convergence=as.character(optimfit$convergence==0),
          Positive_definite=as.character(is.positive.definite(optimfit$hessian)))) 
  # optimfit$convergence==0 means it reached convergence.
  # the positive definite command is from command from matrixcalc
  
  #***************************************************************
  # make summaries of the parameter values ----
  #***************************************************************
  
  results = c()
  
  # print(optimfit$par)
  # print(ginv(optimfit$hessian)) # that's the VARIANCE
  # print(ifelse(is.positive.definite(optimfit$hessian), "positive-definite", 
  #              ifelse(is.negative.definite(optimfit$hessian), "negative-definite","indefinite")))
  cov = ginv(optimfit$hessian)
  
  ses = sqrt(c(cov[1,1], cov[2,2], cov[3,3], cov[4,4]))

  results[1] = paste(format(round(exp(optimfit$par[1]), 2), nsmall=2), "(", 
                     format(round(exp(optimfit$par[1] - 1.96*ses[1]), 2), nsmall=2), ", ",
                     format(round(exp(optimfit$par[1] + 1.96*ses[1]), 2), nsmall=2), ")", sep="")
  results[2] = paste(format(round(logit_inv(optimfit$par[2]), 2), nsmall=2), "(", 
                     format(round(logit_inv(optimfit$par[2] - 1.96*ses[2]), 2), nsmall=2), ", ",
                     format(round(logit_inv(optimfit$par[2] + 1.96*ses[2]), 2), nsmall=2), ")", sep="")
  results[3] = paste(format(round(logit_inv(optimfit$par[3]), 2), nsmall=2), "(", 
                      format(round(logit_inv(optimfit$par[3] - 1.96*ses[3]), 2), nsmall=2), ", ",
                      format(round(logit_inv(optimfit$par[3] + 1.96*ses[3]), 2), nsmall=2), ")", sep="")
  results[4] = paste(format(round(logit_inv(optimfit$par[4]), 2), nsmall=2), "(",
                     format(round(logit_inv(optimfit$par[4] - 1.96*ses[4]), 2), nsmall=2), ", ",
                     format(round(logit_inv(optimfit$par[4] + 1.96*ses[4]), 2), nsmall=2), ")", sep="")

  results[5] = ifelse(is.positive.definite(optimfit$hessian), "positive-definite",
                  ifelse(is.negative.definite(optimfit$hessian), "negative-definite","indefinite"))
  results[6] = ifelse(optimfit$convergence==0, "converged", "not converged")

  names(results) = c("R0", "m1", "m2", "reporting_prob", "cov_matrix", "convergence") # "m2", 
  
  #***************************************************************
  # Now simulating SIR # to see the predicted values ----
  #***************************************************************
  
    m1 = plogis(optimfit$par[2])
    m2 = plogis(optimfit$par[3])
    coeff = c(m1, m2, rep(1, dim(data)[1]-2))
    
    parameters$betap = exp(optimfit$par[1])*coeff*(parameters$mub+parameters$delta)/(1+parameters$delta*parameters$rC*sum(agedist*parameters$theta)/parameters$mub)
    parameters$epsilon = 0 # duration to rate

    times = seq(0, burn, by = 1)
    # bef = Sys.time()
    out_tmp = data.frame(ode(y = init_ode, times = times, func = typhoidsir, parms = parameters, method="euler"))
    # aft= Sys.time()

  # What was the incidence in the last year of the simulation?
  simcases = (out_tmp[burn, paste('cumI1', 1:length(data$pop), sep="")] - out_tmp[burn-52, paste('cumI1', 1:length(data$pop), sep="")])*
    logit_inv(optimfit$par[3])*data$sens_BC*data$part_prob

  # fit_results = list(optimfit, modelcost, ghn_ll, logit_inv, logit, data, parameters, pars, typhoidsir_water, simcases)
  
  fit_results = list(optimfit=optimfit, modelcost=modelcost, ghn_ll=ghn_ll, results=results,
                     data=data, parameters=parameters, pars=pars, simcases=simcases)
  
  return(fit_results)
}
