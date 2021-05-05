# Simulation

run_sim = function(C, percent_vax, strategy, num_perday, v_e = v_e_constant,
                       u = u_var, sp = 1, se = 0, syn_sero_compartments = NA){
  # New IC: start w/ 0.5% of each age group in I
  # Only run sim for 365 days post start of rollout
  # Vaccine rollout is continuous at 1% of total pop/day until all vaccines are distributed
  
  # Disease Tranmission
  d_E <- 1/3 # incubation period (E -> I), ref: Davies
  d_I <- 1/5 # recovery period (I -> R), ref: Davies
  
  E_0 <- Ev_0 <- Ex_0 <- Sv_0 <- Sx_0 <- Iv_0 <- Ix_0 <- Rv_0 <- Rx_0 <- D_0 <- rep(0, num_groups)
  R_0 <- 0
  
  if (num_perday == 1){
    I_0 <- rep(1, num_groups)
  } else {
    I_0 <- N_i*0.005 # 0.5% of each age group
  }
  
  S_0 <- N_i - I_0 - R_0
  
  # specify group to vaccinate according to allocation strategy
  if (strategy == "all"){ 
    groups <- 1:9
  } else if (strategy == "SPC1"){ 
    groups <- 1
  } else if (strategy == "SPC2") { 
    groups <- 2
  } else if (strategy == "SPC3") {
    groups <- 3
  } else if (strategy == "SPC4") {
    groups <- 4
  } else if (strategy == "SPC5"){ 
    groups <- 5
  } else if (strategy == "SPC6") { 
    groups <- 6
  } else if (strategy == "SPC7") {
    groups <- 7
  } else if (strategy == "SPC8") {
    groups <- 8
  } else if (strategy == "SPC9") {
    groups <- 9
  }
  people_to_vax <- sum(S_0[groups] + E_0[groups] + R_0[groups])
  
  vax_proportion <- rep(0, num_groups)
  vax_proportion[groups] <- (S_0[groups] + E_0[groups] + R_0[groups])/people_to_vax # Proportion of people to vaccinate by category
  
  vax_supply <- percent_vax*pop_total # Number of vaccines available
  
  # anticipatory rollout IC - vaccinate people at t=0
  if (num_perday == 1){
    nvax <- vax_supply
    vax_distribution <- nvax*vax_proportion # Number of vaccines available per category
    S <- S_0
    E <- E_0
    R <- R_0
    vax_eligible <- (S+E)*sp + R*(1-se) # Number of people eligible for vaccines per category
    
    if (any(vax_distribution > vax_eligible)){ # If there some categories where vaccines are available and the category isn't fully vaccinated
      # make sure everyone in the specificed age groups are vaccinated
      if (!all(vax_distribution[groups] > vax_eligible[groups])){ # If not all the people in the prioritized categories are not vaccinated
        temp <- vax_distribution
        temp[vax_distribution > vax_eligible] <- vax_eligible[vax_distribution > vax_eligible]
        leftover_vax <- sum(vax_distribution - temp)
        
        full_groups <- 1:9
        full_groups <- full_groups[vax_distribution > vax_eligible]
        leftover_groups <- groups[!groups %in% full_groups]
        people_to_vax <- sum(vax_eligible[leftover_groups])
        vax_proportion <- rep(0, num_groups)
        vax_proportion[leftover_groups] <- vax_eligible[leftover_groups]/people_to_vax
        # All available vaccines available are given to the prioritized categories
        vax_leftover_dist <- leftover_vax*vax_proportion
        vax_distribution <- temp + vax_leftover_dist 
      }
      
      # Distribute vaccines to all other age groups if there's doses left over after strategy specified
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
    }
    
    alpha <- as.matrix(vax_distribution)/(as.matrix(vax_eligible))
    alpha[alpha == Inf] <- 0
    alpha[is.nan(alpha)] <- 0
    
    S_0  <- S - as.matrix(S*alpha*sp) - as.matrix(S*alpha*(1-sp))
    Sv_0 <- as.matrix(S*alpha*sp*v_e)
    Sx_0 <- as.matrix(S*alpha*sp*(1-v_e)) + as.matrix(S*alpha*(1-sp)) 

    R_0  <- R - as.matrix(R*alpha*(1-se)) - as.matrix(R*alpha*se)
    Rv_0 <- as.matrix(R*alpha*(1-se))
    Rx_0 <- as.matrix(R*alpha*se)
    
    vax_supply <- 0
  }
  
  compartments_initial <- c(S_0,Sv_0,Sx_0,E_0,Ev_0,Ex_0,I_0,Iv_0,Ix_0,R_0,Rv_0,Rx_0,D_0)
  
 # if (length(syn_sero_compartments) > 1){
 #    compartments_initial <- syn_sero_compartments
 # }
  
  compartments_aftervax <- move_vaccinated(compartments_initial, num_perday, vax_supply, vax_proportion,
                                           groups, v_e, v_e_type, sp, se)
  
  parameters = list(u=u, C=C, d_E=d_E, d_I=d_I, v_e=v_e)
  
  running = TRUE
  t <- c(0:1)
  df <- as.data.frame(deSolve::lsoda(compartments_aftervax, t, calculate_derivatives_new, parameters))
  df[1,] <- c(0, compartments_initial)
  
  t <- t + 1
  
  vax_supply <- vax_supply - num_perday*pop_total
  vax_supply[vax_supply < 0] <- 0
  
  while(running == TRUE){
    compartments_initial <- as.numeric(df[t[2], -(1)])
    
    compartments_aftervax <- move_vaccinated(compartments_initial, num_perday, vax_supply, vax_proportion,
                                             groups, v_e, sp, se)
    
    parameters = list(u=u, C=C, d_E=d_E, d_I=d_I, v_e=v_e)
    temp <- as.data.frame(deSolve::lsoda(compartments_aftervax, t, calculate_derivatives_new, parameters))
    row.names(temp) <- t+1
    temp <- temp[-(1),]
    
    df <- rbind(df, temp)
    
    vax_supply <- vax_supply - num_perday*pop_total
    vax_supply[vax_supply < 0] <- 0
    
    if (vax_supply == 0){
      remaining_t = c(t[2]:365)
      inits <- as.numeric(df[t[2]+1, -(1)])
      parameters = list(u=u, C=C, d_E=d_E, d_I=d_I, v_e=v_e)
      temp <- as.data.frame(deSolve::lsoda(inits, remaining_t, calculate_derivatives_new, parameters))
      row.names(temp) <- remaining_t+1
      temp <- temp[-(1),]
      
      df <- rbind(df, temp)
      running = FALSE
    } else if (t[2] == 365){
      running = FALSE
    } else {
      t <- t + 1
    }
  }
  
  names(df) <- c("time", "S1", "S2", "S3", "S4", "S5", "S6", "S7", "S8", "S9",
                 "Sv1", "Sv2", "Sv3", "Sv4", "Sv5", "Sv6", "Sv7", "Sv8", "Sv9",
                 "Sx1", "Sx2", "Sx3", "Sx4", "Sx5", "Sx6", "Sx7", "Sx8", "Sx9",
                 "E1", "E2", "E3", "E4", "E5", "E6", "E7", "E8", "E9",
                 "Ev1", "Ev2", "Ev3", "Ev4", "Ev5", "Ev6", "Ev7", "Ev8", "Ev9",
                 "Ex1", "Ex2", "Ex3", "Ex4", "Ex5", "Ex6", "Ex7", "Ex8", "Ex9",
                 "I1", "I2", "I3", "I4", "I5", "I6", "I7", "I8", "I9",
                 "Iv1", "Iv2", "Iv3", "Iv4", "Iv5", "Iv6", "Iv7", "Iv8", "Iv9",
                 "Ix1", "Ix2", "Ix3", "Ix4", "Ix5", "Ix6", "Ix7", "Ix8", "Ix9",
                 "R1", "R2", "R3", "R4", "R5", "R6", "R7", "R8", "R9",
                 "Rv1", "Rv2", "Rv3", "Rv4", "Rv5", "Rv6", "Rv7", "Rv8", "Rv9",
                 "Rx1", "Rx2", "Rx3", "Rx4", "Rx5", "Rx6", "Rx7", "Rx8", "Rx9",
                 "D1", "D2", "D3", "D4", "D5", "D6", "D7", "D8", "D9")
  
  return(df)
}
