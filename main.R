# Vaccination strategy

setwd("~/Desktop/vaccine_MLC")
source("run_sim.R")
source("outils.R")

library(deSolve)
library(doParallel)

#######################################################################################################
# GRAPH THEMES #
################

nolabels_theme <- theme(axis.title.x =element_blank(),
                        axis.text.x = element_blank(),
                        axis.title.y = element_blank(),
                        axis.text.y = element_blank(),
                        plot.title = element_text(size = 12, face = "plain"),
                        legend.position = "none")
onlyx_theme <- theme(axis.title.y = element_blank(),
                     axis.text.y = element_blank(),
                     plot.title = element_text(size = 12, face = "plain"),
                     legend.position = "none")
onlyy_theme <- theme(axis.title.x = element_blank(),
                     axis.text.x = element_blank(),
                     plot.title = element_text(size = 12, face = "plain"),
                     legend.position = "none")


#######################################################################################################
# PARAMETERS #
##############

country <- 'FRA' # Choose country to study

C <- read.csv(file = paste0('C_'+country+'.csv')) # Counting matrix

spc_demo <- read.csv(file = paste0('SPC'+country+'.csv')) # Vector containing the size of all the different SCP categories (% of pop) and total population

pop_total <- spc_demo[10] # Country population (initial)
spc_demo <- spc_demo[1:9] # Vector containing the size of all the different SCP categories (% of pop)

N_i <- pop_total*spc_demo # Vector containing the population per SCP categories
num_groups <- length(spc_demo) # Number of SCPs chosen

IFR <- c() # IFR per CSP
IFR <- IFR/100 # en %

u_var <- c() # Susceptibility per CSP
R0 <- compute_R0(u_var, C) # Computing of R0

v_e <-  0.9 # Vaccine efficacy


#####################################################################################################
# SIMULATION #
##############

# Parallel core use to not overload computer
cores=detectCores()
cl <- makeCluster(cores[1]-1)
registerDoParallel(cl)

# Different SCP categories
list_all <- list_spc1 <- list_spc2 <- list_spc3 <- list_spc4 <- list_spc5 <- list_spc6 <- list_spc7 <- list_spc8 <- list_spc9 <- vector(mode = "list")
list_all_var <- list_1_var <- list_2_var <- list_3_var <- list_4_var <- list_5_var <- list_6_var <- list_7_var <- list_8_var <- list_9_var <- vector(mode = "list")

# Vaccine rollout speed
num_per_day <- 1

# Simulation   per SCP
for (i in seq(0, 50, by = 1)){
  j <- i/100
  list_all[[paste0(i)]] <- run_sim(C, j, "all", num_per_day, v_e)
  list_scp1[[paste0(i)]] <- run_sim(C, j, "SPC1", num_per_day, v_e)
  list_scp2[[paste0(i)]] <- run_sim(C, j, "SPC2", num_per_day, v_e)
  list_scp3[[paste0(i)]] <- run_sim(C, j, "SPC3", num_per_day, v_e)
  list_scp4[[paste0(i)]] <- run_sim(C, j, "SPC4", num_per_day, v_e)
  list_scp5[[paste0(i)]] <- run_sim(C, j, "SPC5", num_per_day, v_e)
  list_scp6[[paste0(i)]] <- run_sim(C, j, "SPC6", num_per_day, v_e)
  list_scp7[[paste0(i)]] <- run_sim(C, j, "SPC7", num_per_day, v_e)
  list_scp8[[paste0(i)]] <- run_sim(C, j, "SPC8", num_per_day, v_e)
  list_scp9[[paste0(i)]] <- run_sim(C, j, "SPC9", num_per_day, v_e)
}

#####################################################################################################
# STAT PANEL DISPLAY #
######################


# Plotting of death and cases graphs
p_mort <- plot_over_vax_avail_new("deaths", "None", list_all, list_spc1, list_spc2, list_spc3, list_spc4, list_spc5, list_spc6, list_spc7, list_spc8, list_spc9) 
p_infect <- plot_over_vax_avail_new("cases", "None", list_all, list_spc1, list_spc2, list_spc3, list_spc4, list_spc5, list_spc6, list_spc7, list_spc8, list_spc9) 
  
fig1_plots <- list(p_mort, p_infect)

stopCluster(cl) # End of parallel computing 
proc.time() - ptm # Time the program spent running unitl now

# Graphic plotting
mort_1 <- fig1_plots[[1]][[1]] + 
  onlyy_theme + 
  theme(plot.margin=unit(c(0.25,0.25,0.25,0.25), "cm"))
mort_2 <- fig1_plots[[2]][[1]] + 
  nolabels_theme
mort_3 <- fig1_plots[[3]][[1]] + 
  nolabels_theme

infect_1 <- fig1_plots[[1]][[2]]  + 
  theme(axis.title.x = element_blank())
infect_2 <- fig1_plots[[2]][[2]] + 
  onlyx_theme + 
  theme(axis.title.x = element_blank())
infect_3 <- fig1_plots[[3]][[2]] +
  onlyx_theme  + 
  theme(axis.title.x = element_blank())

# Panel config
panel <- ggarrange(mort_1, mort_2, mort_3, infect_1, infect_2, infect_3,
                   nrow = 2,
                   labels = c('B', 'C', 'D', 'E', 'F', 'G'),
                   label.args = list(gp = grid::gpar(fontsize=12, fontface = "bold"),
                                     hjust=0, vjust = 1.1),
                   bottom = textGrob("Total vaccine supply (% of population)", gp = gpar(fontsize = 12), 
                                     vjust = 0.3, hjust = 0.36),
                   padding = unit(0.5, "line"))

p1 <- barplot_vax_strat("SCP1") + 
  theme(axis.title.y = element_blank())
#ylab("Distribution\nof vaccines (%)")
p2 <- barplot_vax_strat("SCP2") + 
  theme(axis.title.y = element_blank())
p3 <- barplot_vax_strat("SCP3") + 
  theme(axis.title.y = element_blank())
p4 <- barplot_vax_strat("SCP4") + 
  theme(axis.title.y = element_blank())
p5 <- barplot_vax_strat("SCP5") + 
  theme(axis.title.y = element_blank())
p6 <- barplot_vax_strat("SCP6") + 
  theme(axis.title.y = element_blank())
p7 <- barplot_vax_strat("SCP7") + 
  theme(axis.title.y = element_blank())
p8 <- barplot_vax_strat("SCP8") + 
  theme(axis.title.y = element_blank())
p9 <- barplot_vax_strat("SCP9") + 
  theme(axis.title.y = element_blank())
p_all <- barplot_vax_strat("all") + 
  theme(axis.title.y = element_blank())

strategy_panel <- ggarrange(p1, p2, p3, p4, p5, p6, p7, p8, p9, p_all
                            nrow = 5, 
                            labels = c('A',  '', '', '', ''),
                            label.args = list(gp = grid::gpar(fontsize=12, fontface = "bold"),
                                              hjust=1.8, vjust = 1.1),
                            left = textGrob("Distribution of vaccines (%)", rot = 90, hjust = 0.5),
                            bottom = textGrob("SCP", vjust = 0))
# export as pdf 9.5x4"
grid.arrange(strategy_panel, panel,
             ncol = 2, widths = c(2, 7.5),
             padding = unit(1, "line"))


