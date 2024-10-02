typhoidsir_vax = function(time, state, parameters) {
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
  
  Sv = as.matrix(unname(state[(8*parameters$al+1):(9*parameters$al)]))
  Iv = as.matrix(unname(state[(9*parameters$al+1):(10*parameters$al)]))
  V1 = as.matrix(unname(state[(10*parameters$al+1):(11*parameters$al)]))
  V2 = as.matrix(unname(state[(11*parameters$al+1):(12*parameters$al)]))
  cumIv = as.matrix(unname(state[(12*parameters$al+1):(13*parameters$al)]))
  vaxr = as.matrix(unname(state[(13*parameters$al+1):(14*parameters$al)]))
  vaxc = as.matrix(unname(state[(14*parameters$al+1):(15*parameters$al)]))
  
  births = as.matrix(c(parameters$mub*parameters$population, rep(0, parameters$al-1))) 
    #births don't come from the state matrix.
  
  S1[S1<0]=0
  S2[S2<0]=0
  Sv[Sv<0]=0
  
  I1[I1<0]=0
  I2[I2<0]=0
  Iv[Iv<0]=0
  
  with(as.list(S1=S1, I1=I1, R=R, C=C, S2=S2, I2=I2, cumI1=cumI1, lambda=lambda,
               Sv=Sv, Iv=Iv, V1=V1, V2=V2, cumIv=cumIv, vaxr=vaxr, vaxc=vaxc, # now=now, 
               births=births, parameters), {

    # add terms for vaccination.
      dS1 = births - lambda*S1 + epsilon*S2 - (u+mu)*S1 + c(0,u[1:(al-1)])*c(0,S1[1:(al-1)]) -
        uinr*c(0,u[1:(al-1)])*c(0,S1[1:(al-1)])*ifelse(time>burntime, 1, 0)*rcov - # vprob???
        uinc*S1*(ifelse(time>camp1 & time<=camp2, 1, 0)*ccov)
      dI1 = lambda*S1 - delta*I1 - (u+mu)*I1 + c(0,u[1:(al-1)])*c(0,I1[1:(al-1)])
      dR = delta*(1-alpha-theta)*(I1+Iv) + delta*(1-theta2)*I2 - omega*R - (u+mu)*R + c(0,u[1:(al-1)])*c(0,R[1:(al-1)]) -
        uinr*c(0,u[1:(al-1)])*c(0,R[1:(al-1)])*ifelse(time>burntime, 1, 0)*rcov*vprob -  
        uinc*R*(ifelse(time>camp1 & time<=camp2, 1, 0)*ccov*vprob) 
      dC = delta*(theta*(I1+Iv)+theta2*I2) - (u+mu)*C + c(0,u[1:(al-1)])*c(0,C[1:(al-1)])
      dS2 = omega*R + omegav*V2 - (epsilon+lambda)*S2 - (u+mu)*S2 + c(0,u[1:(al-1)])*c(0,S2[1:(al-1)]) - 
        uinr*c(0,u[1:(al-1)])*c(0,S2[1:(al-1)])*ifelse(time>burntime, 1, 0)*rcov*vprob - # c(0,u[1:(al-1)]) is missing
        uinc*S2*ifelse(time>camp1 & time<=camp2, 1, 0)*ccov*vprob
      dI2 = lambda*S2 - delta*I2 - (u+mu)*I2 + c(0,u[1:(al-1)])*c(0,I2[1:(al-1)])
      dcumI1 = lambda*S1 # Cumulative (symptomatic) incidence. No one ages out of here.
      dlambda = betap*sum(I1+Iv+r*I2+rC*C)/population - lambda 

      # uin is a vector of 0's and 1's that indicates if an age groups is eligible for vaccination
      dSv = omegav*V1 - lambda*Sv - (u+mu)*Sv + c(0,u[1:(al-1)])*c(0,Sv[1:(al-1)]) + 
        uinr*c(0,u[1:(al-1)])*c(0,S1[1:(al-1)])*ifelse(time>burntime, 1, 0)*rcov*(1-vprob) + 
        uinc*S1*(ifelse(time>camp1 & time<=camp2, 1, 0)*ccov*(1-vprob))
      dIv = lambda*Sv - delta*Iv - (u+mu)*Iv + c(0,u[1:(al-1)])*c(0,Iv[1:(al-1)])
      dV1 = -omegav*V1 - (u+mu)*V1 + c(0,u[1:(al-1)])*c(0,V1[1:(al-1)]) + 
        uinr*c(0,u[1:(al-1)])*c(0,S1[1:(al-1)])*ifelse(time>burntime, 1, 0)*rcov*vprob + 
        uinc*S1*ifelse(time>camp1 & time<=camp2, 1, 0)*ccov*vprob
      dV2 = -omegav*V2 - (u+mu)*V2 + c(0,u[1:(al-1)])*c(0,V2[1:(al-1)]) + 
        uinr*c(0,u[1:(al-1)])*c(0,S2[1:(al-1)])*ifelse(time>burntime, 1, 0)*rcov*vprob + 
        uinc*S2*ifelse(time>camp1 & time<=camp2, 1, 0)*ccov*vprob +
        uinr*c(0,u[1:(al-1)])*c(0,R[1:(al-1)])*ifelse(time>burntime, 1, 0)*rcov*vprob +
        uinc*R*ifelse(time>camp1 & time<=camp2, 1, 0)*ccov*vprob
        
      dcumIv = lambda*Sv # Cumulative (symptomatic) incidence for vaccinated individuals. No one ages out of here.
      dvaxr = uinr*(c(0,S1[1:(al-1)])+c(0,I1[1:(al-1)])+
                    c(0,R[1:(al-1)])+c(0,C[1:(al-1)])+
                    c(0,S2[1:(al-1)])+c(0,I2[1:(al-1)]))*ifelse(time>burntime, 1, 0)*rcov*c(0,u[1:(al-1)])
            # include vax given to people who didn't change status (i.e. I1, R, C, etc)
      dvaxc = uinc*(S1+I1+R+C+S2+I2)*ifelse(time>camp1 & time<=camp2, 1, 0)*ccov
            # include vax given to people who didn't change status (i.e. I1, R, C, etc) 

    list(c(dS1, dI1, dR, dC, dS2, dI2, dcumI1, dlambda, dSv, dIv, dV1, dV2, dcumIv, dvaxr, dvaxc)) 
  })
}
