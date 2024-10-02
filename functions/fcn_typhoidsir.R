typhoidsir = function(time, state, parameters) {
  # state will have the compartment names (S1, I1, R, C, S2, I2, cumI1) and lambda
  # parameters will have the parameters (which I usually denote in Greek letters)
    # in this case, "births" is also a parameter, as will "population"
    # al, the number of age groups, and u, the rate of aging out, will also be a parameter
  # time will have the time points at which we want the function evaluated 
  
  S1 = as.matrix(unname(state[1:parameters$al]))
  I1 = as.matrix(unname(state[(parameters$al+1):(2*parameters$al)]))
  R = as.matrix(unname(state[(2*parameters$al+1):(3*parameters$al)]))
  C = as.matrix(unname(state[(3*parameters$al+1):(4*parameters$al)]))
  S2 = as.matrix(unname(state[(4*parameters$al+1):(5*parameters$al)]))
  I2 = as.matrix(unname(state[(5*parameters$al+1):(6*parameters$al)]))
  cumI1 = as.matrix(unname(state[(6*parameters$al+1):(7*parameters$al)]))
  lambda = as.matrix(unname(state[(7*parameters$al+1):(8*parameters$al)]))

  births = as.matrix(c(parameters$mub*parameters$population, rep(0, parameters$al-1))) 
    #births don't come from the state matrix.
  
  S1[S1<0]=0
  S2[S2<0]=0
  
  I1[I1<0]=0
  I2[I2<0]=0
  
  with(as.list(S1=S1, I1=I1, R=R, C=C, S2=S2, I2=I2, cumI1=cumI1, lambda=lambda, births=births, parameters), {

      dS1 = births - lambda*S1 + epsilon*S2 - (u+mu)*S1 + c(0,u[1:(al-1)])*c(0,S1[1:(al-1)])
      dI1 = lambda*S1 - delta*I1 - (u+mu)*I1 + c(0,u[1:(al-1)])*c(0,I1[1:(al-1)])
      dR = delta*(1-alpha-theta)*I1 + delta*(1-theta2)*I2 - omega*R - (u+mu)*R + c(0,u[1:(al-1)])*c(0,R[1:(al-1)])
      dC = delta*(theta*I1+theta2*I2) - (u+mu)*C + c(0,u[1:(al-1)])*c(0,C[1:(al-1)])
      dS2 = omega*R - (epsilon+lambda)*S2 - (u+mu)*S2 + c(0,u[1:(al-1)])*c(0,S2[1:(al-1)])
      dI2 = lambda*S2 - delta*I2 - (u+mu)*I2 + c(0,u[1:(al-1)])*c(0,I2[1:(al-1)])
      dcumI1 = lambda*S1 # Cumulative (symptomatic) incidence. No one ages out of here.
      dlambda = betap*sum(I1+r*I2+rC*C)/population - lambda 

    list(c(dS1, dI1, dR, dC, dS2, dI2, dcumI1, dlambda))
  })
}
