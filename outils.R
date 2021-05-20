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

plot_over_vax_avail_new = function(outcome, var_name, list_all, list_spc1, list_spc2, list_spc3, list_spc4, list_spc5, list_spc6, list_spc7, list_spc8, list_spc9){
  library(ggplot2)
  theme_set(theme_minimal(base_size = 12))
  
  x_spc1_switch        <- when_strat_switch(list_spc1,1)
  x_spc2_switch        <- when_strat_switch(list_spc2,2)
  x_spc3_switch        <- when_strat_switch(list_spc3,3)
  x_spc4_switch        <- when_strat_switch(list_spc4,4)
  x_spc5_switch        <- when_strat_switch(list_spc5,5)
  x_spc6_switch        <- when_strat_switch(list_spc6,6)
  x_spc7_switch        <- when_strat_switch(list_spc7,7)
  x_spc8_switch        <- when_strat_switch(list_spc8,8)
  x_spc9_switch        <- when_strat_switch(list_spc9,9)
  

  # get dataframe for specific outcome
  if (outcome == "cases"){
    df <- get_reduction_in_cases_df_new(list_all, list_spc1, list_spc2, list_spc3, list_spc4, list_spc5, list_spc6, list_spc7, list_spc8, list_spc9)
  } else if (outcome == "deaths"){
    df <- get_reduction_in_deaths_df_new(list_all, list_spc1, list_spc2, list_spc3, list_spc4, list_spc5, list_spc6, list_spc7, list_spc8, list_spc9)
  } 
  
  # reallocate to everyone after groups are completely vaccinated according to the strat
  points_x <- c(x_spc1_switch, x_spc2_switch, x_spc3_switch, x_spc4_switch, x_spc5_switch, x_spc6_switch, x_spc7_switch, x_spc8_switch, x_spc9_switch)
  
  points_y <- c(df[df$strat == "SPC1" & df$vax_avail == x_spc1_switch & df$variable == "constant", ]$reduction,
                df[df$strat == "SPC2" & df$vax_avail == x_spc2_switch & df$variable == "constant", ]$reduction,
                df[df$strat == "SPC3" & df$vax_avail == x_spc3_switch & df$variable == "constant", ]$reduction,
                df[df$strat == "SPC4" & df$vax_avail == x_spc4_switch & df$variable == "constant", ]$reduction,
                df[df$strat == "SPC5" & df$vax_avail == x_spc5_switch & df$variable == "constant", ]$reduction,
                df[df$strat == "SPC6" & df$vax_avail == x_spc6_switch & df$variable == "constant", ]$reduction,
                df[df$strat == "SPC7" & df$vax_avail == x_spc7_switch & df$variable == "constant", ]$reduction,
                df[df$strat == "SPC8" & df$vax_avail == x_spc8_switch & df$variable == "constant", ]$reduction,
                df[df$strat == "SPC9" & df$vax_avail == x_spc9_switch & df$variable == "constant", ]$reduction)

  
  # plot
  # for fig 2
  df_var <- df[df$variable == "var", ]
  df_const <- df[df$variable == "constant", ]
  p <- ggplot() +
    geom_line(aes(x = df_const$vax_avail, y = df_const$reduction, col = df_const$strat), 
              linetype = "solid", size = 1, alpha = 0.9) +
    geom_line(aes(x = df_var[df_var$strat == "optimal",]$vax_avail, 
                  y = df_var[df_var$strat == "optimal",]$reduction, 
                  col = df_var[df_var$strat == "optimal",]$strat),
              linetype = "solid", size = 1, alpha = 0.9) +
    xlab("Total vaccine supply (% of population)") +
    scale_color_brewer(palette = "Dark2", name = "Allocation Strategy",
                       labels =  c("All", "SPC1", "SPC2",  
                                   "SPC3", "SPC4", "SPC5", "SPC6", 
                                   "SPC7", "SPC8", "SPC9")) +
    scale_linetype_discrete(name = df$var_name,
                            labels = c("Constant", "Age-dependent")) +
    scale_y_continuous(expand = c(0,0), limit = c(0, 100)) +
    scale_x_continuous(expand = c(0,0), limit = c(0, 50)) +#, breaks = c(0,25,50)) +
    coord_fixed(0.5*4/5) +
    theme(legend.position = "none") +
    guides(colour = guide_legend(override.aes = list(size=3)))
  
  p <- p + geom_point(aes(x = points_x, y = points_y), size = 1.5) + 
    geom_point(aes(x = points_x_var, y = points_y_var), size = 1.5) 
  
  if (outcome == "cases"){ p <- p + ylab("Reduction in\ninfections (%)")}
  else if (outcome == "deaths"){ p <- p + ylab("Reduction in\ndeaths (%)")}
  else if (outcome == "YLL"){ p <- p + ylab("Reduction in\nYLL (%)")}
  return(p)
}

get_reduction_in_cases_df_new = function(list_all, list_spc1, list_spc2, list_spc3, list_spc4, list_spc5, list_spc6, list_spc7, list_spc8, list_spc9){
  total_cases <- rep(NA, 510)
  
  count <- 1
  for (i in list_all){
    total_cases[count] <- compute_total_cases_new(i)
    count <- count + 1
  }
  for (i in list_spc1){
    total_cases[count] <- compute_total_cases_new(i)
    count <- count + 1
  }
  for (i in list_spc2){
    total_cases[count] <- compute_total_cases_new(i)
    count <- count + 1
  }
  for (i in list_spc3){
    total_cases[count] <- compute_total_cases_new(i)
    count <- count + 1
  }
  for (i in list_spc4){
    total_cases[count] <- compute_total_cases_new(i)
    count <- count + 1
  }
  for (i in list_spc5){
    total_cases[count] <- compute_total_cases_new(i)
    count <- count + 1
  }
  for (i in list_spc6){
    total_cases[count] <- compute_total_cases_new(i)
    count <- count + 1
  }
  for (i in list_spc7){
    total_cases[count] <- compute_total_cases_new(i)
    count <- count + 1
  }
  for (i in list_spc8){
    total_cases[count] <- compute_total_cases_new(i)
    count <- count + 1
  }
  for (i in list_spc9){
    total_cases[count] <- compute_total_cases_new(i)
    count <- count + 1
  }
  
  num_strategies <- 10
  vax_avail <- c(rep(seq(0, 50, by = 1), num_strategies*2))
  num_per_list <- 51
  strat <- c(rep("all", num_per_list*2), rep("SPC1", num_per_list*2), 
             rep("SPC2", num_per_list*2), rep("SPC3", num_per_list*2), 
             rep("SPC4", num_per_list*2), rep("SPC5", num_per_list*2),
             rep("SPC6", num_per_list*2), rep("SPC7", num_per_list*2),
             rep("SPC8", num_per_list*2), rep("SPC9", num_per_list*2))
  temp <-  c(rep("constant", num_per_list), rep("var", num_per_list))
  variable <- c(rep(temp, num_strategies))
  
  baseline_cases <- compute_total_cases_new(list_all$`0`)
  baseline_cases_var <- compute_total_cases_new(list_all_var$`0`)
  
  temp <- c(rep(baseline_cases, num_per_list), rep(baseline_cases_var, num_per_list))
  baseline_cases <- c(rep(temp, num_strategies))
  
  reduction <- (1-(total_cases/baseline_cases))*100
  
  df <- data.frame(vax_avail, strat, reduction, variable)
}

get_reduction_in_deaths_df_new = function(list_all, list_spc1, list_spc2, list_spc3, list_spc4, list_spc5, list_spc6, list_spc7, list_spc8, list_spc9){
  total_deaths <- rep(NA, 510)
  
  count <- 1
  for (i in list_all){
    total_deaths[count] <- compute_total_deaths_new(i)
    count <- count + 1
  }
  for (i in list_spc1){
    total_deaths[count] <- compute_total_deaths_new(i)
    count <- count + 1
  }
  for (i in list_spc2){
    total_deaths[count] <- compute_total_deaths_new(i)
    count <- count + 1
  }
  for (i in list_spc3){
    total_deaths[count] <- compute_total_deaths_new(i)
    count <- count + 1
  }
  for (i in list_spc4){
    total_deaths[count] <- compute_total_deaths_new(i)
    count <- count + 1
  }
  for (i in list_spc5){
    total_deaths[count] <- compute_total_deaths_new(i)
    count <- count + 1
  }
  for (i in list_spc6){
    total_deaths[count] <- compute_total_deaths_new(i)
    count <- count + 1
  }
  for (i in list_spc7){
    total_deaths[count] <- compute_total_deaths_new(i)
    count <- count + 1
  }
  for (i in list_spc8){
    total_deaths[count] <- compute_total_deaths_new(i)
    count <- count + 1
  }
  for (i in list_spc9){
    total_deaths[count] <- compute_total_deaths_new(i)
    count <- count + 1
  }
  
  num_strategies <- 10
  vax_avail <- c(rep(seq(0, 50, by = 1), num_strategies*2))
  num_per_list <- 51
  strat <- c(rep("all", num_per_list*2), rep("SPC1", num_per_list*2), 
             rep("SPC2", num_per_list*2), rep("SPC3", num_per_list*2), 
             rep("SPC4", num_per_list*2), rep("SPC5", num_per_list*2),
             rep("SPC6", num_per_list*2), rep("SPC7", num_per_list*2),
             rep("SPC8", num_per_list*2), rep("SPC9", num_per_list*2))
  temp <-  c(rep("constant", num_per_list), rep("var", num_per_list))
  variable <- c(rep(temp, num_strategies))
  
  baseline_deaths <- compute_total_deaths_new(list_all$`0`)
  baseline_deaths_var <- compute_total_deaths_new(list_all_var$`0`)
  
  temp <- c(rep(baseline_deaths, num_per_list), rep(baseline_deaths_var, num_per_list))
  baseline_deaths <- c(rep(temp, num_strategies))
  
  reduction <- (1 - (total_deaths/baseline_deaths))*100
  
  df <- data.frame(vax_avail, strat, reduction, variable)
}

compute_total_deaths_new = function(df){
  deaths <- rep(0,num_groups)
  D_index <- 110

  for (i in 1:num_groups) {
    if (is.null(df[dim(df)[1], D_index])){
    }
    else {
      deaths[i] <- df[dim(df)[1], D_index]
      D_index <- D_index + 1
    }
  }
  tot_deaths <- sum(deaths)/pop_total * 100
}

compute_total_cases_new = function(df){
  infections <- rep(0,num_groups) # number of age groups
  
  R_index <- 83 # col number for R
  Rv_index <- 92 # col number for Rv
  Rx_index <- 101 # col number for Rx
  D_index <- 110
  
  for (i in 1:num_groups) {
    # infections = total # recovered - initial recovered (seropositive)
    if (is.null(df[dim(df)[1], R_index]) || is.null(df[1, R_index]) || is.null(df[dim(df)[1], Rv_index]) || is.null(df[1, Rv_index]) || is.null(df[dim(df)[1], Rx_index]) || is.null(df[1, Rx_index]) || is.null(df[dim(df)[1], D_index])){
    }
    else {
      infections[i] <- df[dim(df)[1], R_index] - df[1, R_index] +
        df[dim(df)[1], Rv_index] - df[1, Rv_index] + 
        df[dim(df)[1], Rx_index] - df[1, Rx_index] + 
        df[dim(df)[1], D_index]
    }
    R_index  <- R_index + 1
    Rv_index <- Rv_index + 1
    Rx_index <- Rx_index + 1
    D_index  <- D_index + 1
  }
  tot_infections <- sum(infections)/pop_total * 100
}

barplot_vax_strat = function(strategy){
  if (strategy == "no vax"){
    vaccinated <- rep(0,num_groups) * 100
    this_color <- "#808080"
    plot_title <- "No Vaccines"
  } else {
    if (strategy == "all"){
      groups <- 1:9
      this_color <- col_all
      plot_title <- "All ages"
    }
    else if (strategy == "SPC1"){
      groups <- 1
      this_color <- col_1
      plot_title <- "SPC1"
    }
    else if (strategy == "SPC2"){
      groups <- 2
      this_color <- col_2
      plot_title <- "SPC2"
    }
    else if (strategy == "SPC3"){
      groups <- 3
      this_color <- col_3
      plot_title <- "SPC3"
    }
    else if (strategy == "SPC4"){
      groups <- 4
      this_color <- col_4
      plot_title <- "SPC4"
    }
    else if (strategy == "SPC5"){
      groups <- 5
      this_color <- col_5
      plot_title <- "SPC5"
    }
    else if (strategy == "SPC6"){
      groups <- 6
      this_color <- col_6
      plot_title <- "SPC6"
    }
    else if (strategy == "SPC7"){
      groups <- 7
      this_color <- col_7
      plot_title <- "SPC7"
    }
    else if (strategy == "SPC8"){
      groups <- 8
      this_color <- col_8
      plot_title <- "SPC8"
    }
    else if (strategy == "SPC9"){
      groups <- 9
      this_color <- col_9
      plot_title <- "SPC9"
    }
    people_to_vax <- sum(N_i[groups])
    
    vax_proportion <- rep(0, num_groups)
    vax_proportion[groups] <- N_i[groups]/people_to_vax
  }
  spc_groups <- c(1, 2, 3, 4, 5, 6, 7, 8, 9)
  vax_proportion <- vax_proportion*100
  
  df <- data.frame(spc_groups, vax_proportion, spc_demo)
  
  # plot
  p <- ggplot(df, aes(x=spc_groups)) + 
    geom_bar(aes(y=vax_proportion), position="stack", stat="identity", fill = this_color) + 
    xlab("SPC") +
    ggtitle(plot_title) +
    scale_y_continuous(expand = c(0,0), limit = c(0, 52), breaks = c(0, 25, 50), labels = c(0, "", 50))
  
  if (strategy == "all"){
    p <- p + scale_x_continuous(limits = c(-5, NA), expand = c(0,0),
                                breaks=seq(-5, 75, by = 10),
                                labels=c("1","2", "3", "4", "5", "6", "7", "8", "9")) +
      theme(plot.title = element_text(color = this_color, size = 10),#, face = "bold"),
            axis.title.x = element_blank(),
            axis.text.x = element_text(size = 9))
  } else {
    p <- p + theme(plot.title = element_text(color = this_color,  size = 10),
                   axis.title.x =element_blank(),
                   axis.text.x = element_blank(),
                   axis.title.y = element_blank()) +
      scale_x_continuous(limits = c(-5, NA), expand = c(0,0),
                         breaks=seq(-5, 75, by = 10),
                         labels=c("", "", "", "", "", "", "", "", ""))
  }
  p <- p + theme(panel.grid.minor = element_blank(),
                 panel.grid.major = element_blank(),
                 axis.ticks = element_line(size = 0.35),
                 axis.line = element_line(size = 0.35),
                 plot.title = element_text(vjust = -1),
                 plot.margin = unit(c(-0.06,0,0,0), "cm"))
  
  return(p)
}

when_strat_switch = function(list_df, groups){
  # find the vaccine supply where you switch from strategy to vaccinating everyone else
  # input: list_df is a list with df for a strategy for vax supply 0-50
  #        groups: vector index of the strategy
  # output: vax supply where switch occurs
  
  other_groups <- 1:9
  other_groups <- other_groups[!other_groups %in% groups]
  x_switch <- -1
  
  for (i in 1:length(list_df)){
    temp <- list_df[[i]]
    if (any(temp[dim(temp)[1], other_groups + 10] > 0)){
      x_switch <- i-1 
      break
    }
  }
  return(x_switch)
}