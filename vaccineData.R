# Script to get the data ready and estimate the qs
rm(list = ls())

setwd("C:/MSc Leiden/thesis/data application")
load('rtss.RData') # Malaria
# Need to select one df from the large list (they are all the same df)
malaria <- rtss[[1]]
load('rv144.RData') # HIV
hiv <- rv144

# Only first 3 columns are of interest in the VIH case
library(tidyverse)
hiv <- hiv |> select(ftime,ftype, vax)

# Same for Malaria
malaria <- malaria |> select(ftime, ftype, vaccine)
malaria <- malaria |> rename(vax = vaccine)
# All the type that are not 0 wil be grouped together for both trials
hiv <- hiv |> mutate(ftype = case_when(ftype %in% c(0) ~ 0,
                                       ftype %in% c(1,2) ~ 1))
malaria <- malaria |> mutate(ftype = case_when(ftype %in% c(0) ~ 0,
                                               ftype %in% c(1,2) ~ 1))


# Obtian cumulative incidence
calculate_incidence <- function(data){
  # Incidence for treatment
  time <- data |> filter(vax == 1) |> pull(ftime)
  status <- data |> filter(vax == 1) |> pull(ftype)
  
  ord <- order(time)
  time <- time[ord]
  status <- status[ord]
  
  tout <- time[status==1]
  tout <- unique(tout)
  
  nout <- length(tout)
  di <- Yi <- rep(NA, nout)
  for (i in 1:nout) {
    di[i] <- sum(time==tout[i] & status==1)
    Yi[i] <- sum(time>=tout[i])
  }
  inc_treatment <- di/Yi
  
  # Incidence for control
  time <- data |> filter(vax == 0) |> pull(ftime)
  status <- data |> filter(vax == 0) |> pull(ftype)
  
  ord <- order(time)
  time <- time[ord]
  status <- status[ord]
  
  tout <- time[status==1]
  tout <- unique(tout)
  
  nout <- length(tout)
  di <- Yi <- rep(NA, nout)
  for (i in 1:nout) {
    di[i] <- sum(time==tout[i] & status==1)
    Yi[i] <- sum(time>=tout[i])
  }
  inc_control <- di/Yi
  
  return(list(inc_treatment,inc_control))
}

malaria_inc <- calculate_incidence(malaria)
malaria_inc_treatment <- malaria_inc[[1]]
malaria_inc_control <- malaria_inc[[2]]

hiv_inc <- calculate_incidence(hiv)
hiv_inc_treatment <- hiv_inc[[1]]
hiv_inc_control <- hiv_inc[[2]]


# Compute cumulative incidence
malaria_cuminc_treatment <- cumsum(malaria_inc_treatment)
malaria_cuminc_control <- cumsum(malaria_inc_control)
hiv_cuminc_treatment <- cumsum(hiv_inc_treatment)
hiv_cuminc_control <- cumsum(hiv_inc_control)

# Visualizations

# Make cumulative incidence plots for both trials

# HIV

df_hiv <- tibble(vax_inc = hiv_inc_treatment, placebo_inc = hiv_inc_control,
                 vax_cuminc = hiv_cuminc_treatment, placebo_cuminc = hiv_cuminc_control,
                 time = c(1:6))

# Reshape the dataset for cumulative incidence
df_cuminc <- df_hiv %>%
  select(time, vax_cuminc, placebo_cuminc) %>%
  pivot_longer(cols = c(vax_cuminc, placebo_cuminc), names_to = "group", values_to = "cumulative_incidence") %>%
  mutate(group = recode(group, 
                        vax_cuminc = "Cumulative Incidence on Vaccinated", 
                        placebo_cuminc = "Cumulative Incidence on Placebo"))

# Plot for cumulative incidence
cuminc_plot_hiv <- ggplot(df_cuminc, aes(x = time, y = cumulative_incidence, color = group)) +
  geom_line() +
  geom_point() +
  labs(title = "Cumulative Incidence Over Time", x = "Time (in months)", y = "Cumulative Incidence", color = "Group") +
  theme_minimal()

# Reshape the dataset for incidence
df_inc <- df_hiv %>%
  select(time, vax_inc, placebo_inc) %>%
  pivot_longer(cols = c(vax_inc, placebo_inc), names_to = "group", values_to = "incidence") %>%
  mutate(group = recode(group, 
                        vax_inc = "Incidence on Vaccinated", 
                        placebo_inc = "Incidence on Placebo"))

# Plot for incidence
inc_plot_hiv <- ggplot(df_inc, aes(x = time, y = incidence, color = group)) +
  geom_line() +
  geom_point() +
  labs(title = "Incidence Over Time", x = "Time (in months)", y = "Incidence", color = "Group") +
  theme_minimal()

# Malaria

df_malaria <- tibble(vax_inc = malaria_inc_treatment, placebo_inc = malaria_inc_control,
                 vax_cuminc = malaria_cuminc_treatment, placebo_cuminc = malaria_cuminc_control,
                 time = c(1:12))

# Reshape the dataset for cumulative incidence
df_cuminc <- df_malaria %>%
  select(time, vax_cuminc, placebo_cuminc) %>%
  pivot_longer(cols = c(vax_cuminc, placebo_cuminc), names_to = "group", values_to = "cumulative_incidence") %>%
  mutate(group = recode(group, 
                        vax_cuminc = "Cumulative Incidence on Vaccinated", 
                        placebo_cuminc = "Cumulative Incidence on Placebo"))

# Plot for cumulative incidence
cuminc_plot_mal <- ggplot(df_cuminc, aes(x = time, y = cumulative_incidence, color = group)) +
  geom_line() +
  geom_point() +
  labs(title = "Cumulative Incidence Over Time", x = "Time (in months)", y = "Cumulative Incidence", color = "Group") +
  theme_minimal()

# Reshape the dataset for incidence
df_inc <- df_malaria %>%
  select(time, vax_inc, placebo_inc) %>%
  pivot_longer(cols = c(vax_inc, placebo_inc), names_to = "group", values_to = "incidence") %>%
  mutate(group = recode(group, 
                        vax_inc = "Incidence on Vaccinated", 
                        placebo_inc = "Incidence on Placebo"))

# Plot for incidence
inc_plot_mal <- ggplot(df_inc, aes(x = time, y = incidence, color = group)) +
  geom_line() +
  geom_point() +
  labs(title = "Incidence Over Time", x = "Time (in months)", y = "Incidence", color = "Group") +
  theme_minimal()



# Save the 4 cumulative incidences to use on the main script
save(malaria_cuminc_treatment, malaria_cuminc_control, 
     hiv_cuminc_treatment, hiv_cuminc_control, file = 'cumulative_incidence.RData')
