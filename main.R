# Vaccination strategy

setwd("~/Desktop/Cov_vax_CSP")

library(deSolve)
library(doParallel)
library(ggplot2)
library(ggpubr)
library(grid)


source("sim.R")
source("outils.R")


#######################################################################################################
# PARAMETERS #
##############

# Counting matrices
C_work <- as.matrix(read.csv(file = "Matrices/Work.csv", header=FALSE))
C_home <- as.matrix(read.csv(file = "Matrices/Home.csv", header=FALSE))   
C_school <- as.matrix(read.csv(file = "Matrices/School.csv", header=FALSE)) 
C_other <- as.matrix(read.csv(file = "Matrices/Others.csv", header=FALSE)) 

# Scenarios
S1 <- as.matrix(read.csv(file = "Matrices/S1.csv", header=FALSE))
S2 <- as.matrix(read.csv(file = "Matrices/S2.csv", header=FALSE))
S3 <- as.matrix(read.csv(file = "Matrices/S3.csv", header=FALSE))
scenarii <- list(list(S1,1,0,0), list(S2,1,0.5,0), list(S3,1,0.7,0.4), list(diag(9),1,0,1), list(diag(9),1,1,1))


pop_total <- 67000000 # Country population (initial)
spc_demo <- c(0.007,0.034 , 0.103, 0.131, 0.13, 0.097, 0.0894, 0.3037, 0.1039) # Vector containing the size of all the different SCP categories (% of pop) (Tableau C)

N_i <- pop_total*spc_demo # Vector containing the population per SCP categories
num_groups <- length(spc_demo) # Number of SCPs chosen


IFR <- c(3.594256e-02, 3.594256e-02, 3.594256e-02, 3.594256e-02, 3.594256e-02, 
         3.594256e-02, 3.594256e-02, 4.545632e+00, 3.196070e-02) # IFR per CSP

#IFR <- c(9.530595e-01, 3.196070e-02, 1.071797e-01, 3.594256e-02, 1.205328e-01, 
         # 4.042049e-01, 1.355495e-01, 4.545632e+00, 1.524371e-2) # IFR per ages

IFR <- IFR/100 # en %

u <- c(0.86, 0.86, 0.86, 0.86, 0.86, 0.86, 0.86, 0.74, 0.68)/(39.80957/2) # Susceptibility per CSP

# u <- c(0.8, 0.68, 0.79, 0.86, 0.8, 0.82, 0.88, 0.74, 0.74)/(39.80957/2) # Susceptibility per ages

#R0 <- compute_R0(u_var, C) # Computing of R0

v_e <-  0.95 # Pfizer Covid-19 vaccine efficacy (ref. Israel observational study, 05/2021, Haas et al.)

#####################################################################################################
# SIMULATION #
##############

# Parallel core use to not overload computer
ptm <- proc.time()
cores=detectCores()
cl <- makeCluster(cores[1]-1)
registerDoParallel(cl)

# Different strategies
list_all <- list_strat1 <- list_strat2 <- list_strat3 <- list_strat4 <- list_strat5 <- vector(mode = "list")

# Vaccine rollout speed
num_per_day <- 500000/pop_total # (on average 500,000 people vaccinated per day in France in May 2021)

# Simulation by strategy

i = 0.85 # Number of vaccines to give (proportion of population)

for (k in 1:5){
  scenario <- scenarii[[k]] # Chosen scenario
  C = scenario[[1]] %*% C_work + scenario[[2]]*C_home + scenario[[3]]*C_school + scenario[[4]]*C_other 
  print(paste0("SIM ",k,"th scenario"))

  list_all[[paste0(k)]] <- run_sim(C, i, "All", num_per_day, v_e)
  list_strat1[[paste0(k)]] <- run_sim(C, i, "Elderly", num_per_day, v_e)
  list_strat2[[paste0(k)]] <- run_sim(C, i, "Economic players", num_per_day, v_e)
  list_strat3[[paste0(k)]] <- run_sim(C, i, "Essential workers", num_per_day, v_e)
  list_strat4[[paste0(k)]] <- run_sim(C, i, "No work-from-home", num_per_day, v_e)
  list_strat5[[paste0(k)]] <- run_sim(C, i, "Most contact", num_per_day, v_e)
}

df <- list(list_all, list_strat1, list_strat2, list_strat3, list_strat4, list_strat5)

#######################################################################################################
# GRAPH THEMES #
################
# 
# # ColorBrewer Dark2
# col_1 = "#1B9E77" # teal green
# col_all = "#D95F02" # orange
# col_2 = "#7570B3" # purple 
# col_3 = "#E7298A" # magenta
# col_4 = "#66A61E" # light green 
# col_5 = "#1B9E76"
# col_6 = "#D95F01" 
# col_7 = "#7570B8"
# col_8 = "#E7298A" 
# col_9 = "#61A61E" 
# 
# nolabels_theme <- theme(axis.title.x =element_blank(),
#                         axis.text.x = element_blank(),
#                         axis.title.y = element_blank(),
#                         axis.text.y = element_blank(),
#                         plot.title = element_text(size = 12, face = "plain"),
#                         legend.position = "none")
# onlyx_theme <- theme(axis.title.y = element_blank(),
#                      axis.text.y = element_blank(),
#                      plot.title = element_text(size = 12, face = "plain"),
#                      legend.position = "none")
# onlyy_theme <- theme(axis.title.x = element_blank(),
#                      axis.text.x = element_blank(),
#                      plot.title = element_text(size = 12, face = "plain"),
#                      legend.position = "none")
# 
#
#####################################################################################################
# STAT PANEL DISPLAY #
######################
#
#
# # Plotting of death and cases graphs
# p_mort <- plot_over_vax_avail_new("deaths", "None", list_all, list_spc1, list_spc2, list_spc3, list_spc4, list_spc5, list_spc6, list_spc7, list_spc8, list_spc9) 
# p_infect <- plot_over_vax_avail_new("cases", "None", list_all, list_spc1, list_spc2, list_spc3, list_spc4, list_spc5, list_spc6, list_spc7, list_spc8, list_spc9) 
# 
# print("DEATHS AND CASES PLOTED")
# 
# fig1_plots <- list(p_mort, p_infect)
# 
# stopCluster(cl) # End of parallel computing 
# proc.time() - ptm # Time the program spent running unitl now
# 
# # Graphic plotting
# mort_1 <- fig1_plots[[1]][[1]] + 
#   onlyy_theme + 
#   theme(plot.margin=unit(c(0.25,0.25,0.25,0.25), "cm"))
# mort_2 <- fig1_plots[[2]][[1]] + 
#   nolabels_theme
# 
# infect_1 <- fig1_plots[[1]][[2]]  + 
#   theme(axis.title.x = element_blank())
# infect_2 <- fig1_plots[[2]][[2]] + 
#   onlyx_theme + 
#   theme(axis.title.x = element_blank())
# 
# # Panel config
# print("SETTING UP PANEL...")
# panel <- ggarrange(mort_1, mort_2, infect_1, infect_2,
#                    nrow = 2,
#                    labels = c('B', 'C', 'D', 'E'),
#                    label.args = list(gp = grid::gpar(fontsize=12, fontface = "bold"),
#                                      hjust=0, vjust = 1.1),
#                    bottom = textGrob("Total vaccine supply (% of population)", gp = gpar(fontsize = 12), 
#                                      vjust = 0.3, hjust = 0.36),
#                    padding = unit(0.5, "line"))
# 
# print("CREATING BAR PLOT...")
# p1 <- barplot_vax_strat("SPC1") + 
#   theme(axis.title.y = element_blank())
# #ylab("Distribution\nof vaccines (%)")
# p2 <- barplot_vax_strat("SPC2") + 
#   theme(axis.title.y = element_blank())
# p3 <- barplot_vax_strat("SPC3") + 
#   theme(axis.title.y = element_blank())
# p4 <- barplot_vax_strat("SPC4") + 
#   theme(axis.title.y = element_blank())
# p5 <- barplot_vax_strat("SPC5") + 
#   theme(axis.title.y = element_blank())
# p6 <- barplot_vax_strat("SPC6") + 
#   theme(axis.title.y = element_blank())
# p7 <- barplot_vax_strat("SPC7") + 
#   theme(axis.title.y = element_blank())
# p8 <- barplot_vax_strat("SPC8") + 
#   theme(axis.title.y = element_blank())
# p9 <- barplot_vax_strat("SPC9") + 
#   theme(axis.title.y = element_blank())
# p_all <- barplot_vax_strat("all") + 
#   theme(axis.title.y = element_blank())
# 
# 
# strategy_panel <- ggarrange(p1, p2, p3, p4, p5, p6, p7, p8, p9, p_all,
#                             nrow = 5, 
#                             labels = c('A',  '', '', '', ''),
#                             label.args = list(gp = grid::gpar(fontsize=12, fontface = "bold"),
#                                               hjust=1.8, vjust = 1.1),
#                             left = textGrob("Distribution of vaccines (%)", rot = 90, hjust = 0.5),
#                             bottom = textGrob("SCP", vjust = 0))
# # export as pdf 9.5x4"
# grid.arrange(strategy_panel, panel,
#              ncol = 2, widths = c(2, 7.5),
#              padding = unit(1, "line"))
# 
# print("DONE!")