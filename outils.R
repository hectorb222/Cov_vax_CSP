compute_R0 = function(u, C){
  gamma <- 1/5 # recovery period (I -> R), ref: Davies
  
  # Davies NGM
  Du <- diag(u, 9)
  Dy <- diag(1/gamma, 9)
  
  NGM <- Du %*% C %*% Dy
  R0  <- abs(eigen(NGM)$values[1])
}

move_vaccinated = function(x, num_perday, vax_supply, vax_proportion, groups, v_e, sp, se) {
  # move those who are vaccinated in a given day
  num_compartment <- 13
  num_groups <- length(x)/num_compartment
  S    <- as.matrix(x[1:num_groups])
  Sv   <- as.matrix(x[(1*num_groups+1):(2*num_groups)])
  Sx   <- as.matrix(x[(2*num_groups+1):(3*num_groups)])
  E    <- as.matrix(x[(3*num_groups+1):(4*num_groups)])
  Ev   <- as.matrix(x[(4*num_groups+1):(5*num_groups)])
  Ex   <- as.matrix(x[(5*num_groups+1):(6*num_groups)])
  I    <- as.matrix(x[(6*num_groups+1):(7*num_groups)])
  Iv   <- as.matrix(x[(7*num_groups+1):(8*num_groups)])
  Ix   <- as.matrix(x[(8*num_groups+1):(9*num_groups)])
  R    <- as.matrix(x[(9*num_groups+1):(10*num_groups)])
  Rv   <- as.matrix(x[(10*num_groups+1):(11*num_groups)])
  Rx   <- as.matrix(x[(11*num_groups+1):(12*num_groups)])
  D    <- as.matrix(x[(12*num_groups+1):(13*num_groups)])
  
  if (vax_supply >= num_perday*pop_total){
    nvax <- num_perday*pop_total
  } else {
    nvax <- vax_supply
  }
  
  vax_distribution <- nvax*vax_proportion
  vax_eligible <- (S+E)*sp + R*(1-se)
  #vax_distribution[vax_distribution > (S+E+R)] <- S[vax_distribution > (S+E+R)] + E[vax_distribution > (S+E+R)] + R[vax_distribution > (S+E+R)]
  if (any(vax_distribution > vax_eligible)){
    # make sure everyone in the specificed age groups are vaccinated
    if (!all(vax_distribution[groups] > vax_eligible[groups])){
      temp <- vax_distribution
      temp[vax_distribution > vax_eligible] <- vax_eligible[vax_distribution > vax_eligible]
      leftover_vax <- sum(vax_distribution - temp)
      
      full_groups <- 1:9
      full_groups <- full_groups[vax_distribution > vax_eligible]
      leftover_groups <- groups[!groups %in% full_groups]
      people_to_vax <- sum(vax_eligible[leftover_groups])
      vax_proportion <- rep(0, num_groups)
      vax_proportion[leftover_groups] <- vax_eligible[leftover_groups]/people_to_vax
      
      vax_leftover_dist <- leftover_vax*vax_proportion
      vax_distribution <- temp + vax_leftover_dist
    }
    #vax_distribution[vax_distribution > (S+E+R)] <- S[vax_distribution > (S+E+R)] + E[vax_distribution > (S+E+R)] + R[vax_distribution > (S+E+R)]
    #distribute vaccines to all other age groups if there's doses left over after strategy specified
    #age groups
    if (any(vax_distribution > vax_eligible)){
      temp <- vax_distribution
      temp[vax_distribution > vax_eligible] <- vax_eligible[vax_distribution > vax_eligible]
      leftover_vax <- sum(vax_distribution - temp)
      
      leftover_groups <- 1:9
      leftover_groups <- leftover_groups[!leftover_groups %in% groups]
      people_to_vax <- sum(vax_eligible[leftover_groups])
      vax_proportion <- rep(0, num_groups)
      vax_proportion[leftover_groups] <- vax_eligible[leftover_groups]/people_to_vax
      
      vax_leftover_dist <- leftover_vax*vax_proportion
      vax_distribution <- temp + vax_leftover_dist
    }
    # don't go over the number eligible
    if (any(vax_distribution > vax_eligible)){
      vax_distribution[vax_distribution > vax_eligible] = vax_eligible
    }
  }
  
  alpha <- vax_distribution/vax_eligible
  alpha[alpha == Inf] <- 0 # no people in S,E,R
  alpha[is.nan(alpha)] <- 0 # no vax left and no people in S,E,R
  #print(paste0("t: ", t, " max alpha:", max(alpha))) ### FIXME
  if(any(alpha > 1)){print("ERROR: alpha > 1 in move_vaccinated")}
  

  # all-or-nothing
  dS  <- -as.matrix(S*alpha*sp) - as.matrix(S*alpha*(1-sp))
  dSv <- as.matrix(S*alpha*sp*v_e)
  dSx <- as.matrix(S*alpha*sp*(1-v_e)) + as.matrix(S*alpha*(1-sp)) 
    
  dE  <- - as.matrix(E*alpha*sp) - as.matrix(E*alpha*(1-sp))
  dEv <- as.matrix(E*alpha*sp*v_e)
  dEx <- as.matrix(E*alpha*sp*(1-v_e)) + as.matrix(E*alpha*(1-sp))
  }
  
  dR <- - as.matrix(R*alpha*(1-se)) - as.matrix(R*alpha*se)
  dRv <- as.matrix(R*alpha*(1-se))
  dRx <- as.matrix(R*alpha*se)
  
  # update compartments
  S    <- S + dS
  Sv   <- Sv + dSv
  Sx   <- Sx + dSx
  E    <- E + dE
  Ev   <- Ev + dEv
  Ex   <- Ex + dEx
  R    <- R + dR
  Rv   <- Rv + dRv
  Rx   <- Rx + dRx
  
  # output updated compartments
  out <- c(S,Sv,Sx,E,Ev,Ex,I,Iv,Ix,R,Rv,Rx,D)
}

calculate_derivatives_new=function(t, x, parameters){
  # x is a vector of length (# model compartment types)*(# age groups)
  # S, E, I, R etc. are vectors of length num_groups
  num_compartment <- 13
  num_groups <- length(x)/num_compartment
  S    <- as.matrix(x[1:num_groups])
  Sv   <- as.matrix(x[(1*num_groups+1):(2*num_groups)])
  Sx   <- as.matrix(x[(2*num_groups+1):(3*num_groups)])
  E    <- as.matrix(x[(3*num_groups+1):(4*num_groups)])
  Ev   <- as.matrix(x[(4*num_groups+1):(5*num_groups)])
  Ex   <- as.matrix(x[(5*num_groups+1):(6*num_groups)])
  I    <- as.matrix(x[(6*num_groups+1):(7*num_groups)])
  Iv   <- as.matrix(x[(7*num_groups+1):(8*num_groups)])
  Ix   <- as.matrix(x[(8*num_groups+1):(9*num_groups)])
  R    <- as.matrix(x[(9*num_groups+1):(10*num_groups)])
  Rv   <- as.matrix(x[(10*num_groups+1):(11*num_groups)])
  Rx   <- as.matrix(x[(11*num_groups+1):(12*num_groups)])
  D    <- as.matrix(x[(12*num_groups+1):(13*num_groups)])
  
  S[S < .Machine$double.eps] <- 0
  Sv[Sv < .Machine$double.eps] <- 0
  Sx[Sx < .Machine$double.eps] <- 0
  E[E < .Machine$double.eps] <- 0
  Ev[Ev < .Machine$double.eps] <- 0
  Ex[Ex < .Machine$double.eps] <- 0
  I[I < .Machine$double.eps] <- 0
  Iv[Iv < .Machine$double.eps] <- 0
  Ix[Ix < .Machine$double.eps] <- 0
  R[R < .Machine$double.eps] <- 0
  Rv[Rv < .Machine$double.eps] <- 0
  Rx[Rx < .Machine$double.eps] <- 0
  D[D < .Machine$double.eps] <- 0
  
  # I[I<0] = 0
  # Iv[Iv<0] = 0
  # Ix[Ix<0] = 0
  
  u <- parameters$u
  C <- parameters$C
  d_E <- parameters$d_E
  d_I <- parameters$d_I
  v_e <- parameters$v_e
  v_e_type <- parameters$v_e_type
  
  lambda <- as.matrix((C)%*%as.matrix((I+Iv+Ix)/N_i) * as.matrix(u) )
  

  # all-or-nothing
  dSv <- rep(0, num_groups)
  dEv <- -d_E*as.matrix(Ev)
  
  
  dS  <- -as.matrix(S*lambda)
  dSx <- -as.matrix(Sx*lambda)
  
  dE  <- as.matrix(S*lambda) - d_E*as.matrix(E)
  dEx <- as.matrix(Sx*lambda) - d_E*as.matrix(Ex)
  
  dI  <- as.matrix(E*d_E) - as.matrix(I*d_I)
  dIv <- as.matrix(Ev*d_E) - as.matrix(Iv*d_I)
  dIx <- as.matrix(Ex*d_E) - as.matrix(Ix*d_I)
  
  dR  <- as.matrix(I*d_I*(1-IFR)) 
  dRv <- as.matrix(Iv*d_I*(1-IFR)) 
  dRx <- as.matrix(Ix*d_I*(1-IFR)) 
  
  dD  <- as.matrix(I*d_I*IFR + Iv*d_I*IFR + Ix*d_I*IFR)
  
  out <- c(dS,dSv,dSx,dE,dEv,dEx,dI,dIv,dIx,dR,dRv,dRx,dD)
  list(out)
}