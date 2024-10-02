typhoidsir_maps_wq = function(time, state, parameters) {
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
  
  Sv = as.matrix(unname(state[(8*parameters$al*parameters$num_wq+1):(9*parameters$al*parameters$num_wq)]))
  Iv = as.matrix(unname(state[(9*parameters$al*parameters$num_wq+1):(10*parameters$al*parameters$num_wq)]))
  V1 = as.matrix(unname(state[(10*parameters$al*parameters$num_wq+1):(11*parameters$al*parameters$num_wq)]))
  V2 = as.matrix(unname(state[(11*parameters$al*parameters$num_wq+1):(12*parameters$al*parameters$num_wq)]))
  cumIv = as.matrix(unname(state[(12*parameters$al*parameters$num_wq+1):(13*parameters$al*parameters$num_wq)]))
  vaxr = as.matrix(unname(state[(13*parameters$al*parameters$num_wq+1):(14*parameters$al*parameters$num_wq)]))
  vaxc = as.matrix(unname(state[(14*parameters$al*parameters$num_wq+1):(15*parameters$al*parameters$num_wq)]))
  
  vaxuc1 = as.matrix(unname(state[(15*parameters$al*parameters$num_wq+1):(16*parameters$al*parameters$num_wq)]))
  vaxuc2 = as.matrix(unname(state[(16*parameters$al*parameters$num_wq+1):(17*parameters$al*parameters$num_wq)]))
  vaxuc3 = as.matrix(unname(state[(17*parameters$al*parameters$num_wq+1):(18*parameters$al*parameters$num_wq)]))
  vaxuc4 = as.matrix(unname(state[(18*parameters$al*parameters$num_wq+1):(19*parameters$al*parameters$num_wq)]))
  vaxuc5 = as.matrix(unname(state[(19*parameters$al*parameters$num_wq+1):(20*parameters$al*parameters$num_wq)]))
  vaxuc6 = as.matrix(unname(state[(20*parameters$al*parameters$num_wq+1):(21*parameters$al*parameters$num_wq)]))
  
  births = as.matrix(c(matrix(c(parameters$mub*parameters$population*parameters$pop_wq, rep(0, (parameters$al-1)*parameters$num_wq)), ncol=parameters$num_wq, byrow=T))) 
  #births don't come from the state matrix.
  
  S1[S1<0]=0
  S2[S2<0]=0
  Sv[Sv<0]=0
  
  I1[I1<0]=0
  I2[I2<0]=0
  Iv[Iv<0]=0
  
  with(as.list(S1=S1, I1=I1, R=R, C=C, S2=S2, I2=I2, cumI1=cumI1, lambda=lambda,
               Sv=Sv, Iv=Iv, V1=V1, V2=V2, cumIv=cumIv, vaxr=vaxr, vaxc=vaxc,  
               vaxuc1=vaxuc1, vaxuc2=vaxuc2, vaxuc3=vaxuc3, vaxuc4=vaxuc4, vaxuc5=vaxuc5, vaxuc6=vaxuc6,
               births=births, parameters), {
      dS1 = births - lambda*S1 + epsilon*S2 - (u+mu)*S1 + 
        v*(1-(uinuc1*uc1cov+uinuc2*uc2cov+uinuc3*uc3cov))*c(sapply(1:num_wq, function(i){c(0,S1[((i-1)*al+1):(i*al-1)])})) - # aging IN, minus those who aged into being vaccinated.
        uinuc4*S1*(ifelse(time>camp1 & time<=camp2, 1, 0)*uc4cov) -
        uinuc5*S1*uc5cov -
        uinuc6*S1*uc6cov
      dI1 = lambda*S1 - delta*I1 - (u+mu)*I1 + 
        v*c(sapply(1:num_wq, function(i){c(0,I1[((i-1)*al+1):(i*al-1)])})) # aging in for I1
      dR = delta*(1-alpha-theta)*(I1+Iv) + delta*(1-theta2)*I2 - omega*R - (u+mu)*R + 
        v*(1-(uinuc1*uc1cov+uinuc2*uc2cov+uinuc3*uc3cov)*vprob)*c(sapply(1:num_wq, function(i){c(0,R[((i-1)*al+1):(i*al-1)])})) -
        uinuc4*R*(ifelse(time>camp1 & time<=camp2, 1, 0)*uc4cov)*vprob -
        uinuc5*R*uc5cov*vprob -
        uinuc6*R*uc6cov*vprob 
      dC = delta*(theta*(I1+Iv)+theta2*I2) - (u+mu)*C + 
        v*c(sapply(1:num_wq, function(i){c(0,C[((i-1)*al+1):(i*al-1)])})) # aging in for C
      dS2 = omega*R + omegav*V2 - (epsilon+lambda)*S2 - (u+mu)*S2 + 
        v*(1-(uinuc1*uc1cov+uinuc2*uc2cov+uinuc3*uc3cov)*vprob)*c(sapply(1:num_wq, function(i){c(0,S2[((i-1)*al+1):(i*al-1)])})) - 
        uinuc4*S2*(ifelse(time>camp1 & time<=camp2, 1, 0)*uc4cov)*vprob -
        uinuc5*S2*uc5cov*vprob -
        uinuc6*S2*uc6cov*vprob 
      dI2 = lambda*S2 - delta*I2 - (u+mu)*I2 + 
        v*c(sapply(1:num_wq, function(i){c(0,I2[((i-1)*al+1):(i*al-1)])})) # aging in for I2
      
      dcumI1 = lambda*S1 # Cumulative (symptomatic) incidence. No one ages out of here.
      dlambda = betap_wq*sum(I1+Iv+r*I2+rC*C)/population - lambda 

      dSv = omegav*V1 - lambda*Sv - (u+mu)*Sv + 
        v*c(sapply(1:num_wq, function(i){c(0,Sv[((i-1)*al+1):(i*al-1)])})) + # aging in for Sv
        v*(uinuc1*uc1cov+uinuc2*uc2cov+uinuc3*uc3cov)*c(sapply(1:num_wq, function(i){c(0,S1[((i-1)*al+1):(i*al-1)])}))*(1-vprob) + # aging IN for S (not Sv), minus those who aged into being vaccinated.
        uinuc4*S1*(ifelse(time>camp1 & time<=camp2, 1, 0)*uc4cov)*(1-vprob) +
        uinuc5*S1*uc5cov*(1-vprob) +
        uinuc6*S1*uc6cov*(1-vprob) 
      dIv = lambda*Sv - delta*Iv - (u+mu)*Iv + 
        v*c(sapply(1:num_wq, function(i){c(0,Iv[((i-1)*al+1):(i*al-1)])})) # Iv aging in 
      dV1 = -omegav*V1 - (u+mu)*V1 + 
        v*c(sapply(1:num_wq, function(i){c(0,V1[((i-1)*al+1):(i*al-1)])})) + # V1 aging in
        v*(uinuc1*uc1cov+uinuc2*uc2cov+uinuc3*uc3cov)*c(sapply(1:num_wq, function(i){c(0,S1[((i-1)*al+1):(i*al-1)])}))*vprob + # aging INto vaccinated
        uinuc4*(ifelse(time>camp1 & time<=camp2, 1, 0)*uc4cov)*S1*vprob +
        uinuc5*uc5cov*S1*vprob +
        uinuc6*uc6cov*S1*vprob 
      dV2 = -omegav*V2 - (u+mu)*V2 + 
        v*c(sapply(1:num_wq, function(i){c(0,V2[((i-1)*al+1):(i*al-1)])})) + # V2 aging in
        v*(uinuc1*uc1cov+uinuc2*uc2cov+uinuc3*uc3cov)*c(sapply(1:num_wq, function(i){c(0,S2[((i-1)*al+1):(i*al-1)])}))*vprob + # successful vaccination
        uinuc4*(ifelse(time>camp1 & time<=camp2, 1, 0)*uc4cov)*S2*vprob +
        uinuc5*uc5cov*S2*vprob +
        uinuc6*uc6cov*S2*vprob +    
        v*(uinuc1*uc1cov+uinuc2*uc2cov+uinuc3*uc3cov)*c(sapply(1:num_wq, function(i){c(0,R[((i-1)*al+1):(i*al-1)])}))*vprob +
        uinuc4*(ifelse(time>camp1 & time<=camp2, 1, 0)*uc4cov)*R*vprob +
        uinuc5*uc5cov*R*vprob +
        uinuc6*uc6cov*R*vprob 
      
      dcumIv = lambda*Sv # Cumulative (symptomatic) incidence for vaccinated individuals. No one ages out of here.
      dvaxr = rep(0, num_wq*al)
      dvaxc = rep(0, num_wq*al)

      dvaxuc1 = uinuc1*(c(0,S1[1:(al-1)])+c(0,I1[1:(al-1)])+
                      c(0,R[1:(al-1)])+c(0,C[1:(al-1)])+
                      c(0,S2[1:(al-1)])+c(0,I2[1:(al-1)]))*uc1cov*v # v is the aging in parameter (counterpart to u)
      # include vax given to people who didn't change status (i.e. I1, R, C, etc). So, don't multiply by vprob
      dvaxuc2 = uinuc2*(c(0,S1[1:(al-1)])+c(0,I1[1:(al-1)])+
                          c(0,R[1:(al-1)])+c(0,C[1:(al-1)])+
                          c(0,S2[1:(al-1)])+c(0,I2[1:(al-1)]))*uc2cov*v 
      # include vax given to people who didn't change status (i.e. I1, R, C, etc) 
      dvaxuc3 = uinuc3*(c(0,S1[1:(al-1)])+c(0,I1[1:(al-1)])+
                          c(0,R[1:(al-1)])+c(0,C[1:(al-1)])+
                          c(0,S2[1:(al-1)])+c(0,I2[1:(al-1)]))*uc3cov*v 
      # include vax given to people who didn't change status (i.e. I1, R, C, etc) 
      dvaxuc4 = uinuc4*(S1+I1+R+C+S2+I2)*ifelse(time>camp1 & time<=camp2, 1, 0)*uc4cov
      dvaxuc5 = uinuc5*(S1+I1+R+C+S2+I2)*uc5cov
      dvaxuc6 = uinuc6*(S1+I1+R+C+S2+I2)*uc6cov
      
    list(c(dS1, dI1, dR, dC, dS2, dI2, dcumI1, dlambda, dSv, dIv, dV1, dV2, dcumIv, dvaxr, dvaxc,
           dvaxuc1, dvaxuc2,  dvaxuc3, dvaxuc4, dvaxuc5, dvaxuc6)) 
  })
}
