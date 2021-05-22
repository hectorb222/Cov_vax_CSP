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
# GRAPH THEMES #
################

# ColorBrewer Dark2
col_1 = "#1B9E77" # teal green
col_all = "#D95F02" # orange
col_2 = "#7570B3" # purple 
col_3 = "#E7298A" # magenta
col_4 = "#66A61E" # light green 
col_5 = "#1B9E76"
col_6 = "#D95F01" 
col_7 = "#7570B8"
col_8 = "#E7298A" 
col_9 = "#61A61E" 

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

C <- as.matrix(read.csv(file = "C_FRA.csv")) # Counting matrix

pop_total <- 67000000 # Country population (initial)
spc_demo <- c(0.005, 0.027, 0.083, 0.107, 0.105, 0.078, 0.1, 0.1, 0.1) # Vector containing the size of all the different SCP categories (% of pop)

N_i <- pop_total*spc_demo # Vector containing the population per SCP categories
num_groups <- length(spc_demo) # Number of SCPs chosen

IFR <- c(9.530595e-01, 3.196070e-02, 1.071797e-01, 3.594256e-02, 1.205328e-01, 
         4.042049e-01, 1.355495e-01, 4.545632e+00, 1.524371e-2) # IFR per CSP
IFR <- IFR/100 # en %

u <- c(0.8, 0.68, 0.79, 0.86, 0.8, 0.82, 0.88, 0.74, 0.74)/(39.80957/2) # Susceptibility per CSP
R0 <- compute_R0(u_var, C) # Computing of R0

v_e <-  0.95 # Pfizer Covid-19 vaccine efficacy (ref. Israel observational study, 05/2021, Haas et al.)


#####################################################################################################
# SIMULATION #
##############

# Parallel core use to not overload computer
ptm <- proc.time()
cores=detectCores()
cl <- makeCluster(cores[1]-1)
registerDoParallel(cl)

# Different SCP categories
list_all <- list_spc1 <- list_spc2 <- list_spc3 <- list_spc4 <- list_spc5 <- list_spc6 <- list_spc7 <- list_spc8 <- list_spc9 <- vector(mode = "list")
list_all_var <- list_1_var <- list_2_var <- list_3_var <- list_4_var <- list_5_var <- list_6_var <- list_7_var <- list_8_var <- list_9_var <- vector(mode = "list")
# SPC1: "Agriculteurs exploitants"
# SPC2: "Artisans, commerçants, chefs d'entreprise"
# SPC3: "Cadres et professions intellectuelles supérieures"
# SPC4: "Professions intermédiaires"
# SPC5: "Employés"
# SPC6: "Ouvriers"
# SPC7: "Sans Emploi"
# SPC8: "Retraités"
# SPC9: "Etudiants"

# Vaccine rollout speed
num_per_day <- 0.01

# Simulation   per SCP
for (i in seq(0, 2, by = 1)){
  j <- i/100
  print(paste0("SIM ",i,"% of population vaccinated"))
  list_all[[paste0(i)]] <- run_sim(C, j, "all", num_per_day, v_e)
  print("SIM STRAT ALL")
  list_spc1[[paste0(i)]] <- run_sim(C, 1, "strategy A", num_per_day, v_e)
  print("SIM STRAT A")
  list_spc2[[paste0(i)]] <- run_sim(C, j, "SPC2", num_per_day, v_e)
  print("SIM STRAT 2")
  list_spc3[[paste0(i)]] <- run_sim(C, j, "SPC3", num_per_day, v_e)
  print("SIM STRAT 3")
  list_spc4[[paste0(i)]] <- run_sim(C, j, "SPC4", num_per_day, v_e)
  print("SIM STRAT 4")
  list_spc5[[paste0(i)]] <- run_sim(C, j, "SPC5", num_per_day, v_e)
  print("SIM STRAT 5")
  list_spc6[[paste0(i)]] <- run_sim(C, j, "SPC6", num_per_day, v_e)
  print("SIM STRAT 6")
  list_spc7[[paste0(i)]] <- run_sim(C, j, "SPC7", num_per_day, v_e)
  print("SIM STRAT 7")
  list_spc8[[paste0(i)]] <- run_sim(C, j, "SPC8", num_per_day, v_e)
  print("SIM STRAT 8")
  list_spc9[[paste0(i)]] <- run_sim(C, j, "SPC9", num_per_day, v_e)
  print("SIM STRAT 9")
}

print("SIM DONE")
#####################################################################################################
# STAT PANEL DISPLAY #
######################


# Plotting of death and cases graphs
p_mort <- plot_over_vax_avail_new("deaths", "None", list_all, list_spc1, list_spc2, list_spc3, list_spc4, list_spc5, list_spc6, list_spc7, list_spc8, list_spc9) 
p_infect <- plot_over_vax_avail_new("cases", "None", list_all, list_spc1, list_spc2, list_spc3, list_spc4, list_spc5, list_spc6, list_spc7, list_spc8, list_spc9) 

print("DEATHS AND CASES PLOTED")

fig1_plots <- list(p_mort, p_infect)

stopCluster(cl) # End of parallel computing 
proc.time() - ptm # Time the program spent running unitl now

# Graphic plotting
mort_1 <- fig1_plots[[1]][[1]] + 
  onlyy_theme + 
  theme(plot.margin=unit(c(0.25,0.25,0.25,0.25), "cm"))
mort_2 <- fig1_plots[[2]][[1]] + 
  nolabels_theme

infect_1 <- fig1_plots[[1]][[2]]  + 
  theme(axis.title.x = element_blank())
infect_2 <- fig1_plots[[2]][[2]] + 
  onlyx_theme + 
  theme(axis.title.x = element_blank())

# Panel config
print("SETTING UP PANEL...")
panel <- ggarrange(mort_1, mort_2, infect_1, infect_2,
                   nrow = 2,
                   labels = c('B', 'C', 'D', 'E'),
                   label.args = list(gp = grid::gpar(fontsize=12, fontface = "bold"),
                                     hjust=0, vjust = 1.1),
                   bottom = textGrob("Total vaccine supply (% of population)", gp = gpar(fontsize = 12), 
                                     vjust = 0.3, hjust = 0.36),
                   padding = unit(0.5, "line"))

print("CREATING BAR PLOT...")
p1 <- barplot_vax_strat("SPC1") + 
  theme(axis.title.y = element_blank())
#ylab("Distribution\nof vaccines (%)")
p2 <- barplot_vax_strat("SPC2") + 
  theme(axis.title.y = element_blank())
p3 <- barplot_vax_strat("SPC3") + 
  theme(axis.title.y = element_blank())
p4 <- barplot_vax_strat("SPC4") + 
  theme(axis.title.y = element_blank())
p5 <- barplot_vax_strat("SPC5") + 
  theme(axis.title.y = element_blank())
p6 <- barplot_vax_strat("SPC6") + 
  theme(axis.title.y = element_blank())
p7 <- barplot_vax_strat("SPC7") + 
  theme(axis.title.y = element_blank())
p8 <- barplot_vax_strat("SPC8") + 
  theme(axis.title.y = element_blank())
p9 <- barplot_vax_strat("SPC9") + 
  theme(axis.title.y = element_blank())
p_all <- barplot_vax_strat("all") + 
  theme(axis.title.y = element_blank())


strategy_panel <- ggarrange(p1, p2, p3, p4, p5, p6, p7, p8, p9, p_all,
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

print("DONE!")


