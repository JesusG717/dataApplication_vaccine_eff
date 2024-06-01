# Main script to perform the calculations
rm(list = ls())

setwd("C:/MSc Leiden/thesis/data application")
source('boundsLpSolve.R') # Adds the function calculate_bounds
# input: (q000,q100,q111,q001,q101,q110)
# output (lower_bound, upper_bound)
# Load the cumulative incidences for both trials control and treatment
load('cumulative_incidence.RData')

# Compute the vector of qs for each month (starting at month 2)

CHR_hiv_lb <- rep(NA, length(hiv_cuminc_control)-1)
CHR_hiv_ub <- rep(NA, length(hiv_cuminc_control)-1)

CHR_hiv_lb_mono <- rep(NA, length(hiv_cuminc_control)-1)
CHR_hiv_ub_mono <- rep(NA, length(hiv_cuminc_control)-1)

pi_LL_lb_hiv <- rep(NA, length(hiv_cuminc_control)-1)
pi_LL_ub_hiv <- rep(NA, length(hiv_cuminc_control)-1)

# First do the hiv
for (i in 1:(length(hiv_cuminc_control)-1)) {
  
  q000 <- hiv_cuminc_control[i]
  q001 <- hiv_cuminc_treatment[i]
  q110 <- 1 - hiv_cuminc_control[i+1]
  q111 <- 1 - hiv_cuminc_treatment[i+1]
  q100 <- 1 - q000 - q110
  q101 <- 1 - q001 - q111
  
  pi_LL_lb_hiv[i] <- max(0,1 - hiv_cuminc_treatment[i] - hiv_cuminc_control[i])
  pi_LL_ub_hiv[i] <- min(1 - hiv_cuminc_control[i], 1 - hiv_cuminc_treatment[i])
  
  # Bounds in the general case
  bounds <- calculate_bounds(q000 = q000, q100 = q100, q110 = q110,
                             q001 = q001, q101 = q101, q111 = q111)
  CHR_hiv_lb[i] <- bounds[[1]]
  CHR_hiv_ub[i] <- bounds[[2]]
  
  # Bounds under Monotonicity
  CHR_hiv_lb_mono[i] <- max(0, 1+(q110-q111)/q100)
  CHR_hiv_ub_mono[i] <- min(1, q101/q100)
  
}

# Plot results of the bounds (TEST: probably only good for monotone case)

df_CHR_hiv <- tibble(CHR_lb = CHR_malaria_lb_mono, CHR_ub = CHR_malaria_ub_mono, 
                     time = c(2:12))

# Plot with shaded areas between lower and upper bounds
ggplot(df_CHR_hiv, aes(x = time)) +
  geom_ribbon(aes(ymin = CHR_lb, ymax = CHR_ub), fill = "grey70", alpha = 0.5) +
  labs(title = "CHR Bounds Over Time", x = "Time (in months)", y = "CHR", color = "Group") +
  theme_minimal()

# Now malaria

CHR_malaria_lb <- rep(NA, length(malaria_cuminc_control)-1)
CHR_malaria_ub <- rep(NA, length(malaria_cuminc_control)-1)

CHR_malaria_lb_mono <- rep(NA, length(malaria_cuminc_control)-1)
CHR_malaria_ub_mono <- rep(NA, length(malaria_cuminc_control)-1)

pi_LL_ub_mal <- rep(NA, length(malaria_cuminc_control)-1)
pi_LL_lb_mal <- rep(NA, length(malaria_cuminc_control)-1)

for (i in 1:(length(malaria_cuminc_control)-1)) {
  q000 <- malaria_cuminc_control[i]
  q001 <- malaria_cuminc_treatment[i]
  q110 <- 1 - malaria_cuminc_control[i+1]
  q111 <- 1 - malaria_cuminc_treatment[i+1]
  q100 <- 1 - q000 - q110
  q101 <- 1 - q001 - q111
  
  pi_LL_lb_mal[i] <- max(0, 1 - malaria_cuminc_treatment[i] - malaria_cuminc_control[i])
  pi_LL_ub_mal[i] <- min(1 - malaria_cuminc_control[i], 1 - malaria_cuminc_treatment[i])
  
  # Bounds in the general case
  bounds <- calculate_bounds(q000 = q000, q100 = q100, q110 = q110,
                             q001 = q001, q101 = q101, q111 = q111)
  CHR_malaria_lb[i] <- bounds[[1]]
  CHR_malaria_ub[i] <- bounds[[2]]
  
  # Bounds under Monotonicity
  CHR_malaria_lb_mono[i] <- max(0, 1+(q110-q111)/q100)
  CHR_malaria_ub_mono[i] <- min(1, q101/q100)
}


