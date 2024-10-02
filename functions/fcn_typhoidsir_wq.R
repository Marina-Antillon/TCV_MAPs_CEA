typhoidsir_wq = function(time, state, parameters) {
  # state will have the compartment names (S1, I1, R, C, S2, I2, cumI1) and lambda
  # parameters will have the parameters (which I usually denote in Greek letters)
    # in this case, "births" is also a parameter, as will "population"
    # al, the number of age groups, and u, the rate of aging out, will also be a parameter
  # time will have the time points at which we want the function evaluated 
  
  S1 = as.matrix(unname(state[1:(parameters$al*parameters$num_wq)]))
  I1 = as.matrix(unname(state[(parameters$al*parameters$num_wq+1):(2*parameters$al*parameters$num_wq)]))
  R = as.matrix(unname(state[(2*parameters$al*parameters$num_wq+1):(3*parameters$al*parameters$num_wq)]))
  C = as.matrix(unname(state[(3*parameters$al*parameters$num_wq+1):(4*parameters$al*parameters$num_wq)]))
  S2 = as.matrix(unname(state[(4*parameters$al*parameters$num_wq+1):(5*parameters$al*parameters$num_wq)]))
  I2 = as.matrix(unname(state[(5*parameters$al*parameters$num_wq+1):(6*parameters$al*parameters$num_wq)]))
  cumI1 = as.matrix(unname(state[(6*parameters$al*parameters$num_wq+1):(7*parameters$al*parameters$num_wq)]))
  lambda = as.matrix(unname(state[(7*parameters$al*parameters$num_wq+1):(8*parameters$al*parameters$num_wq)]))
  
  births = as.matrix(rep(c(parameters$mub*parameters$population/parameters$num_wq, rep(0, parameters$al-1)), times=parameters$num_wq)) 
    #births don't come from the state matrix.
  
  S1[S1<0]=0
  S2[S2<0]=0
  
  I1[I1<0]=0
  I2[I2<0]=0

  with(as.list(S1=S1, I1=I1, R=R, C=C, S2=S2, I2=I2, cumI1=cumI1, lambda=lambda,
               births=births, parameters), {

      dS1 = births - lambda*S1 + epsilon*S2 - (u+mu)*S1 + 
        v*c(sapply(1:num_wq, function(i){c(0,S1[((i-1)*al+1):(i*al-1)])})) # aging IN, minus those who aged into being vaccinated.
      dI1 = lambda*S1 - delta*I1 - (u+mu)*I1 + 
        v*c(sapply(1:num_wq, function(i){c(0,I1[((i-1)*al+1):(i*al-1)])})) # aging in for I1
      dR = delta*(1-alpha-theta)*I1 + delta*(1-theta2)*I2 - omega*R - (u+mu)*R + 
        v*c(sapply(1:num_wq, function(i){c(0,R[((i-1)*al+1):(i*al-1)])}))
      dC = delta*(theta*I1+theta2*I2) - (u+mu)*C + 
        v*c(sapply(1:num_wq, function(i){c(0,C[((i-1)*al+1):(i*al-1)])})) # aging in for C
      dS2 = omega*R - (epsilon+lambda)*S2 - (u+mu)*S2 + 
        v*c(sapply(1:num_wq, function(i){c(0,S2[((i-1)*al+1):(i*al-1)])})) 
      dI2 = lambda*S2 - delta*I2 - (u+mu)*I2 + 
        v*c(sapply(1:num_wq, function(i){c(0,I2[((i-1)*al+1):(i*al-1)])})) # aging in for I2
      dcumI1 = lambda*S1 # Cumulative (symptomatic) incidence. No one ages out of here.
      dlambda = betap_wq*sum(I1+r*I2+rC*C)/population - lambda 

    list(c(dS1, dI1, dR, dC, dS2, dI2, dcumI1, dlambda)) 
  })
}
